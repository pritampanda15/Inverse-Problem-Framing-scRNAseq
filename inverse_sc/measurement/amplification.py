"""
Amplification bias model: PCR and IVT amplification effects.

Models how different transcripts are amplified at different rates
during library preparation.
"""

import torch
import torch.nn as nn
import numpy as np
from typing import Optional, Dict


class AmplificationModel(nn.Module):
    """
    Models amplification bias in scRNA-seq library preparation.

    Amplification bias depends on:
    1. GC content (GC-rich transcripts amplify differently)
    2. Transcript length
    3. Secondary structure
    4. Primer binding efficiency

    Parameters
    ----------
    n_genes : int
        Number of genes
    model_bias : bool
        Whether to model gene-specific amplification bias
    """

    def __init__(
        self,
        n_genes: int,
        model_bias: bool = True,
    ):
        super().__init__()
        self.n_genes = n_genes
        self.model_bias = model_bias

        if model_bias:
            # Log-scale amplification factors
            # Positive = amplified more than average
            # Negative = amplified less than average
            self.log_amplification_factor = nn.Parameter(
                torch.zeros(n_genes)
            )
        else:
            self.register_buffer(
                'log_amplification_factor',
                torch.zeros(n_genes)
            )

    def forward(
        self,
        Z_captured: torch.Tensor,
        n_cycles: Optional[int] = None,
    ) -> torch.Tensor:
        """
        Apply amplification bias.

        Parameters
        ----------
        Z_captured : torch.Tensor
            Captured molecules (n_cells, n_genes)
        n_cycles : Optional[int]
            Number of PCR cycles (if modeling explicitly)

        Returns
        -------
        Z_amplified : torch.Tensor
            Amplified library
        """
        # Amplification factor per gene
        amp_factor = torch.exp(self.log_amplification_factor)

        # Apply gene-specific amplification
        Z_amplified = Z_captured * amp_factor

        return Z_amplified

    def calibrate_from_spike_ins(
        self,
        spike_in_counts: np.ndarray,
        spike_in_true_amounts: np.ndarray,
        spike_in_features: Optional[Dict[str, np.ndarray]] = None,
    ) -> Dict[str, np.ndarray]:
        """
        Calibrate amplification bias from ERCC spike-ins.

        Spike-ins have known input amounts, so we can measure
        how amplification distorts relative abundances.

        Parameters
        ----------
        spike_in_counts : np.ndarray
            Observed counts for spike-ins (n_cells, n_spike_ins)
        spike_in_true_amounts : np.ndarray
            Known input concentrations
        spike_in_features : Optional[Dict[str, np.ndarray]]
            Features like GC content, length

        Returns
        -------
        calibration_results : Dict[str, np.ndarray]
            Estimated amplification parameters
        """
        # Mean observed counts across cells
        mean_counts = spike_in_counts.mean(axis=0)

        # Expected proportionality: counts ∝ true_amounts
        # Deviation from proportionality = amplification bias

        # Normalize both to mean=1 for comparison
        norm_counts = mean_counts / mean_counts.mean()
        norm_amounts = spike_in_true_amounts / spike_in_true_amounts.mean()

        # Amplification factor
        amp_factor = norm_counts / (norm_amounts + 1e-8)

        # Log scale
        log_amp_factor = np.log(amp_factor + 1e-8)

        # If we have features (GC content, length), fit a model
        if spike_in_features is not None:
            # This would use regression to predict amp factor from features
            # Then apply to all genes
            # For now, placeholder
            pass

        return {
            'amplification_factor': amp_factor,
            'log_amplification_factor': log_amp_factor,
        }

    def estimate_from_technical_replicates(
        self,
        X_rep1: np.ndarray,
        X_rep2: np.ndarray,
    ) -> Dict[str, np.ndarray]:
        """
        Estimate amplification noise from technical replicates.

        If the same cell is sequenced twice, differences are due to
        technical variation (capture + amplification + sequencing).

        Parameters
        ----------
        X_rep1 : np.ndarray
            Counts from replicate 1 (n_cells, n_genes)
        X_rep2 : np.ndarray
            Counts from replicate 2 (same cells)

        Returns
        -------
        noise_estimates : Dict[str, np.ndarray]
            Technical variation parameters
        """
        # Coefficient of variation across replicates
        mean_expr = (X_rep1 + X_rep2) / 2
        diff = X_rep1 - X_rep2

        # CV = std / mean
        cv_squared = (diff ** 2).mean(axis=0) / (2 * mean_expr.mean(axis=0) ** 2 + 1e-8)

        return {
            'technical_cv': np.sqrt(cv_squared),
            'mean_expression': mean_expr.mean(axis=0),
        }


class GCBiasCorrector:
    """
    Corrects for GC content bias in amplification.

    GC-rich sequences often amplify poorly in PCR.
    """

    @staticmethod
    def estimate_gc_bias(
        observed_counts: np.ndarray,
        gc_content: np.ndarray,
    ) -> Dict[str, np.ndarray]:
        """
        Fit relationship between GC content and amplification.

        Parameters
        ----------
        observed_counts : np.ndarray
            Observed expression (n_cells, n_genes)
        gc_content : np.ndarray
            GC content per gene (0-1 scale)

        Returns
        -------
        gc_correction : Dict[str, np.ndarray]
            GC bias curve
        """
        # Normalize counts by median
        median_expr = np.median(observed_counts, axis=0)

        # Bin by GC content
        gc_bins = np.linspace(0, 1, 20)
        bin_indices = np.digitize(gc_content, gc_bins)

        # Median expression in each GC bin
        expr_by_gc = np.array([
            median_expr[bin_indices == i].mean()
            for i in range(len(gc_bins))
        ])

        # Normalize by overall median
        global_median = np.median(median_expr)
        gc_bias_curve = expr_by_gc / (global_median + 1e-8)

        return {
            'gc_bins': gc_bins,
            'gc_bias_curve': gc_bias_curve,
        }

    @staticmethod
    def correct_gc_bias(
        counts: np.ndarray,
        gc_content: np.ndarray,
        gc_bias_curve: np.ndarray,
        gc_bins: np.ndarray,
    ) -> np.ndarray:
        """
        Apply GC bias correction.

        Parameters
        ----------
        counts : np.ndarray
            Raw counts
        gc_content : np.ndarray
            GC content per gene
        gc_bias_curve : np.ndarray
            Bias values per GC bin
        gc_bins : np.ndarray
            GC bin edges

        Returns
        -------
        corrected_counts : np.ndarray
            GC-corrected counts
        """
        bin_indices = np.digitize(gc_content, gc_bins)
        correction_factors = gc_bias_curve[bin_indices]

        corrected = counts / (correction_factors + 1e-8)
        return corrected
