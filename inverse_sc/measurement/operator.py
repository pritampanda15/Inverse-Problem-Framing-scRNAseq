"""
Complete measurement operator: Composes all measurement components.

Implements the full forward model: Z_true → X_observed
"""

import torch
import torch.nn as nn
import numpy as np
from typing import Optional, Dict, Tuple

from .capture import CaptureModel
from .amplification import AmplificationModel
from .sequencing import SequencingModel


class MeasurementOperator(nn.Module):
    """
    Complete measurement model for scRNA-seq.

    Composes: Capture → Amplification → Sequencing

    This is the forward model M in the inverse problem:
        X_observed = M(Z_true) + noise

    Parameters
    ----------
    n_genes : int
        Number of genes
    n_cells : int
        Number of cells
    model_capture : bool
        Whether to model capture efficiency
    model_amplification : bool
        Whether to model amplification bias
    model_sequencing : bool
        Whether to model sequencing depth effects
    """

    def __init__(
        self,
        n_genes: int,
        n_cells: int,
        model_capture: bool = True,
        model_amplification: bool = True,
        model_sequencing: bool = True,
    ):
        super().__init__()
        self.n_genes = n_genes
        self.n_cells = n_cells

        # Sub-models
        self.capture_model = CaptureModel(
            n_genes=n_genes,
            n_cells=n_cells,
            use_gene_features=model_capture,
            use_cell_features=model_capture,
        )

        self.amplification_model = AmplificationModel(
            n_genes=n_genes,
            model_bias=model_amplification,
        )

        self.sequencing_model = SequencingModel(
            model_depth_effects=model_sequencing,
        )

    def forward(
        self,
        Z_true: torch.Tensor,
        library_sizes: Optional[torch.Tensor] = None,
        return_intermediates: bool = False,
        cell_indices: Optional[torch.Tensor] = None,
    ) -> torch.Tensor:
        """
        Forward model: True biology → Observed counts.

        Parameters
        ----------
        Z_true : torch.Tensor
            True transcriptional state (n_cells, n_genes)
        library_sizes : Optional[torch.Tensor]
            Sequencing depth per cell
        return_intermediates : bool
            If True, return intermediate steps
        cell_indices : Optional[torch.Tensor]
            Indices of cells in current batch (for mini-batching)

        Returns
        -------
        X_observed : torch.Tensor
            Observed count matrix
        intermediates : Dict[str, torch.Tensor] (if return_intermediates=True)
            Intermediate states
        """
        # Step 1: Capture
        Z_captured, capture_prob = self.capture_model(Z_true, cell_indices=cell_indices)

        # Step 2: Amplification
        Z_amplified = self.amplification_model(Z_captured)

        # Step 3: Sequencing
        if library_sizes is None:
            library_sizes = Z_amplified.sum(dim=1)

        X_observed = self.sequencing_model.deterministic_expectation(
            Z_amplified,
            library_sizes,
        )

        if return_intermediates:
            return X_observed, {
                'Z_true': Z_true,
                'Z_captured': Z_captured,
                'capture_prob': capture_prob,
                'Z_amplified': Z_amplified,
                'X_observed': X_observed,
            }
        else:
            return X_observed

    def sample_observation(
        self,
        Z_true: torch.Tensor,
        library_sizes: Optional[torch.Tensor] = None,
    ) -> torch.Tensor:
        """
        Stochastically sample observation (with noise).

        Parameters
        ----------
        Z_true : torch.Tensor
            True expression
        library_sizes : Optional[torch.Tensor]
            Sequencing depth

        Returns
        -------
        X_sampled : torch.Tensor
            Sampled counts (with stochasticity)
        """
        # Stochastic capture
        Z_captured, capture_prob = self.capture_model(Z_true)
        Z_captured_sampled = self.capture_model.sample_capture(
            Z_true,
            capture_prob,
        )

        # Amplification (deterministic bias, could add noise)
        Z_amplified = self.amplification_model(Z_captured_sampled)

        # Stochastic sequencing
        if library_sizes is None:
            library_sizes = Z_amplified.sum(dim=1)

        X_sampled = self.sequencing_model(Z_amplified, library_sizes)

        return X_sampled

    def calibrate(
        self,
        X_observed: np.ndarray,
        spike_ins: Optional[Dict[str, np.ndarray]] = None,
    ) -> Dict[str, Dict]:
        """
        Calibrate measurement model from observed data.

        Parameters
        ----------
        X_observed : np.ndarray
            Observed count matrix (n_cells, n_genes)
        spike_ins : Optional[Dict[str, np.ndarray]]
            Spike-in data for calibration

        Returns
        -------
        calibration_results : Dict[str, Dict]
            Calibration results for each sub-model
        """
        results = {}

        # Calibrate capture model from dropout patterns
        results['capture'] = self.capture_model.calibrate_from_data(X_observed)

        # Calibrate amplification from spike-ins if available
        if spike_ins is not None:
            results['amplification'] = self.amplification_model.calibrate_from_spike_ins(
                spike_in_counts=spike_ins['counts'],
                spike_in_true_amounts=spike_ins['true_amounts'],
            )

        return results

    def effective_sensitivity(
        self,
        Z_true: torch.Tensor,
    ) -> torch.Tensor:
        """
        Compute effective sensitivity: E[X | Z_true].

        This is the Jacobian diagonal: ∂E[X_i] / ∂Z_i

        For small perturbations: ΔX ≈ sensitivity * ΔZ

        Parameters
        ----------
        Z_true : torch.Tensor
            True expression

        Returns
        -------
        sensitivity : torch.Tensor
            How observed counts change with true expression
        """
        # Capture probability
        _, capture_prob = self.capture_model(Z_true)

        # Amplification factor
        amp_factor = torch.exp(self.amplification_model.log_amplification_factor)

        # Overall sensitivity (ignoring sequencing sampling)
        sensitivity = capture_prob * amp_factor

        return sensitivity


class InverseProblemDiagnostics:
    """
    Diagnostic tools for assessing inverse problem difficulty.
    """

    @staticmethod
    def condition_number(
        measurement_operator: MeasurementOperator,
        Z_test: torch.Tensor,
    ) -> float:
        """
        Estimate condition number of the inverse problem.

        High condition number → ill-posed problem,
        small changes in X lead to large changes in inferred Z.

        Parameters
        ----------
        measurement_operator : MeasurementOperator
            The forward model
        Z_test : torch.Tensor
            Test expression values

        Returns
        -------
        condition_number : float
            Estimated condition number
        """
        # Compute sensitivity (diagonal Jacobian)
        sensitivity = measurement_operator.effective_sensitivity(Z_test)

        # Condition number ≈ max(sensitivity) / min(sensitivity)
        cond = sensitivity.max() / (sensitivity.min() + 1e-8)

        return cond.item()

    @staticmethod
    def identifiability_analysis(
        X_observed: np.ndarray,
        measurement_operator: MeasurementOperator,
    ) -> Dict[str, np.ndarray]:
        """
        Analyze which aspects of Z are identifiable from X.

        Parameters
        ----------
        X_observed : np.ndarray
            Observed data
        measurement_operator : MeasurementOperator
            Forward model

        Returns
        -------
        identifiability : Dict[str, np.ndarray]
            Which genes/cells are well-identified
        """
        # Genes with high dropout are hard to identify
        dropout_rate = (X_observed == 0).mean(axis=0)

        # Genes with low total counts are uncertain
        total_counts = X_observed.sum(axis=0)

        # Heuristic identifiability score
        # High dropout or low counts → low identifiability
        identifiability_score = (1 - dropout_rate) * np.log1p(total_counts)

        return {
            'gene_identifiability': identifiability_score,
            'dropout_rate': dropout_rate,
            'total_counts': total_counts,
        }
