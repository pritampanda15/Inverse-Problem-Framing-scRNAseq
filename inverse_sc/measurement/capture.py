"""
Capture efficiency model: Cell and gene-specific dropout.

Models the probability that an mRNA molecule is captured during
cell lysis and reverse transcription.
"""

import torch
import torch.nn as nn
import numpy as np
from typing import Optional, Dict, Tuple


class CaptureModel(nn.Module):
    """
    Models capture efficiency in scRNA-seq.

    Capture probability depends on:
    1. Gene-specific factors (GC content, transcript length)
    2. Cell-specific factors (cell size, lysis efficiency)
    3. Technical factors (protocol, batch effects)

    Parameters
    ----------
    n_genes : int
        Number of genes
    n_cells : int
        Number of cells
    use_gene_features : bool
        Whether to model gene-specific capture rates
    use_cell_features : bool
        Whether to model cell-specific capture rates
    """

    def __init__(
        self,
        n_genes: int,
        n_cells: int,
        use_gene_features: bool = True,
        use_cell_features: bool = True,
    ):
        super().__init__()
        self.n_genes = n_genes
        self.n_cells = n_cells
        self.use_gene_features = use_gene_features
        self.use_cell_features = use_cell_features

        # Base capture rate (global)
        self.register_buffer('base_capture_logit', torch.tensor(0.0))

        # Gene-specific capture efficiency
        if use_gene_features:
            self.gene_capture_logit = nn.Parameter(
                torch.zeros(n_genes)
            )
        else:
            self.register_buffer('gene_capture_logit', torch.zeros(n_genes))

        # Cell-specific capture efficiency
        if use_cell_features:
            self.cell_capture_logit = nn.Parameter(
                torch.zeros(n_cells)
            )
        else:
            self.register_buffer('cell_capture_logit', torch.zeros(n_cells))

    def forward(
        self,
        Z_true: torch.Tensor,
        gene_features: Optional[torch.Tensor] = None,
        cell_features: Optional[torch.Tensor] = None,
        cell_indices: Optional[torch.Tensor] = None,
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        """
        Compute capture probabilities and simulate capture.

        Parameters
        ----------
        Z_true : torch.Tensor
            True expression matrix (n_cells, n_genes)
        gene_features : Optional[torch.Tensor]
            Gene-level features (e.g., GC content, length)
        cell_features : Optional[torch.Tensor]
            Cell-level features (e.g., size, quality metrics)
        cell_indices : Optional[torch.Tensor]
            Indices of cells in current batch (for mini-batching)

        Returns
        -------
        Z_captured : torch.Tensor
            Expected captured molecules
        capture_prob : torch.Tensor
            Capture probability matrix (n_cells, n_genes)
        """
        n_cells, n_genes = Z_true.shape

        # Compute capture probability
        logit = self.base_capture_logit + self.gene_capture_logit

        # Add cell-specific effects
        if self.use_cell_features:
            if cell_indices is not None:
                # Mini-batching: select cell-specific parameters for batch
                cell_logit = self.cell_capture_logit[cell_indices]
            else:
                cell_logit = self.cell_capture_logit[:n_cells]
            logit = logit + cell_logit.unsqueeze(1)

        capture_prob = torch.sigmoid(logit)

        # Expected captured molecules (continuous approximation of binomial)
        Z_captured = Z_true * capture_prob

        return Z_captured, capture_prob

    def sample_capture(
        self,
        Z_true: torch.Tensor,
        capture_prob: torch.Tensor,
    ) -> torch.Tensor:
        """
        Stochastically sample captured molecules.

        Uses Binomial sampling: for each molecule, independently
        decide if it's captured.

        Parameters
        ----------
        Z_true : torch.Tensor
            True expression (can be non-integer in expectation)
        capture_prob : torch.Tensor
            Capture probability per gene/cell

        Returns
        -------
        Z_captured : torch.Tensor
            Sampled captured molecules
        """
        # For large counts, use Normal approximation to Binomial
        # For small counts, use Binomial directly

        mean = Z_true * capture_prob
        var = Z_true * capture_prob * (1 - capture_prob)

        # Use Poisson approximation for efficiency
        # (exact Binomial is expensive for large counts)
        return torch.poisson(mean)

    def calibrate_from_data(
        self,
        X_obs: np.ndarray,
        initial_estimate: Optional[np.ndarray] = None,
    ) -> Dict[str, np.ndarray]:
        """
        Calibrate capture model from observed data.

        Uses zero-inflation patterns to estimate capture rates.
        Genes with high dropout likely have low capture efficiency.

        Parameters
        ----------
        X_obs : np.ndarray
            Observed count matrix (n_cells, n_genes)
        initial_estimate : Optional[np.ndarray]
            Initial estimate of true expression

        Returns
        -------
        calibration_results : Dict[str, np.ndarray]
            Estimated capture parameters
        """
        n_cells, n_genes = X_obs.shape

        # Estimate capture rate from dropout frequency
        # Higher dropout → lower capture rate
        dropout_rate = (X_obs == 0).mean(axis=0)

        # Convert to capture probability
        # Assumption: if true expression > 0, then dropout due to capture
        # This is a rough heuristic; Bayesian inference is more principled

        # Genes with 90% dropout → ~10% capture rate
        # Genes with 10% dropout → ~90% capture rate
        capture_rate_estimate = 1 - dropout_rate

        # Convert to logit scale
        capture_logit = np.log(capture_rate_estimate / (1 - capture_rate_estimate + 1e-8))

        with torch.no_grad():
            if self.use_gene_features:
                self.gene_capture_logit.copy_(
                    torch.from_numpy(capture_logit).float()
                )

        return {
            'capture_prob_estimate': capture_rate_estimate,
            'dropout_rate': dropout_rate,
        }


class DropoutModel:
    """
    Utility class for analyzing dropout patterns.

    Helps distinguish biological zeros from technical dropout.
    """

    @staticmethod
    def estimate_dropout_curve(
        X_obs: np.ndarray,
        mean_expression: Optional[np.ndarray] = None,
    ) -> Dict[str, np.ndarray]:
        """
        Estimate relationship between expression level and dropout.

        Parameters
        ----------
        X_obs : np.ndarray
            Observed counts (n_cells, n_genes)
        mean_expression : Optional[np.ndarray]
            Mean expression per gene

        Returns
        -------
        dropout_curve : Dict[str, np.ndarray]
            Dropout probability as function of expression
        """
        if mean_expression is None:
            mean_expression = X_obs.mean(axis=0)

        dropout_rate = (X_obs == 0).mean(axis=0)

        # Bin by expression level
        expr_bins = np.percentile(mean_expression, np.linspace(0, 100, 20))
        bin_indices = np.digitize(mean_expression, expr_bins)

        dropout_by_bin = np.array([
            dropout_rate[bin_indices == i].mean()
            for i in range(len(expr_bins))
        ])

        return {
            'expression_bins': expr_bins,
            'dropout_by_bin': dropout_by_bin,
            'gene_dropout': dropout_rate,
        }
