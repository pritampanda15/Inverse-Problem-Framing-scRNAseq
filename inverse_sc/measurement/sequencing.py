"""
Sequencing model: Library sampling and depth effects.

Models the final sequencing step where amplified libraries are
sampled to a fixed depth.
"""

import torch
import torch.nn as nn
import numpy as np
from typing import Optional, Dict
import pyro
import pyro.distributions as dist


class SequencingModel(nn.Module):
    """
    Models sequencing as multinomial sampling.

    Given amplified library, sequencing samples a fixed number
    of reads, creating depth variation across cells.

    Parameters
    ----------
    model_depth_effects : bool
        Whether to model cell-specific depth variation
    """

    def __init__(
        self,
        model_depth_effects: bool = True,
    ):
        super().__init__()
        self.model_depth_effects = model_depth_effects

    def forward(
        self,
        Z_amplified: torch.Tensor,
        library_sizes: Optional[torch.Tensor] = None,
    ) -> torch.Tensor:
        """
        Simulate sequencing as multinomial sampling.

        Parameters
        ----------
        Z_amplified : torch.Tensor
            Amplified library (n_cells, n_genes)
        library_sizes : Optional[torch.Tensor]
            Target sequencing depth per cell (n_cells,)
            If None, use observed library sizes

        Returns
        -------
        X_sequenced : torch.Tensor
            Observed counts after sequencing
        """
        n_cells, n_genes = Z_amplified.shape

        if library_sizes is None:
            # Use sum of amplified library as depth
            library_sizes = Z_amplified.sum(dim=1)

        # Normalize to probabilities
        probs = Z_amplified / (Z_amplified.sum(dim=1, keepdim=True) + 1e-8)

        # Sample from multinomial
        # Note: Pyro's Multinomial is more efficient than torch
        X_sequenced = pyro.sample(
            "sequencing",
            dist.Multinomial(
                total_count=library_sizes.unsqueeze(1),
                probs=probs,
            )
        )

        return X_sequenced

    def deterministic_expectation(
        self,
        Z_amplified: torch.Tensor,
        library_sizes: torch.Tensor,
    ) -> torch.Tensor:
        """
        Deterministic expectation (no sampling).

        Useful for gradient-based optimization.

        Parameters
        ----------
        Z_amplified : torch.Tensor
            Amplified library
        library_sizes : torch.Tensor
            Sequencing depth

        Returns
        -------
        E[X] : torch.Tensor
            Expected counts
        """
        # E[X_i] = N * p_i, where p_i = Z_i / sum(Z)
        probs = Z_amplified / (Z_amplified.sum(dim=1, keepdim=True) + 1e-8)
        expected_counts = library_sizes.unsqueeze(1) * probs

        return expected_counts


class LibrarySize:
    """
    Utilities for library size (sequencing depth) analysis.
    """

    @staticmethod
    def normalize_by_size_factors(
        X_obs: np.ndarray,
        method: str = "median",
    ) -> np.ndarray:
        """
        Classic library size normalization.

        Parameters
        ----------
        X_obs : np.ndarray
            Raw counts (n_cells, n_genes)
        method : str
            'median' (median-of-ratios) or 'total' (total counts)

        Returns
        -------
        X_normalized : np.ndarray
            Size-factor normalized counts
        """
        if method == "total":
            # Simple total count normalization
            lib_sizes = X_obs.sum(axis=1)
            median_lib = np.median(lib_sizes)
            size_factors = lib_sizes / median_lib

            X_normalized = X_obs / size_factors[:, np.newaxis]

        elif method == "median":
            # DESeq2-style median-of-ratios
            # More robust to highly expressed genes

            # Geometric mean per gene
            log_means = np.log(X_obs + 1).mean(axis=0)
            geo_means = np.exp(log_means)

            # Ratio to geometric mean
            ratios = X_obs / (geo_means + 1e-8)

            # Median ratio per cell = size factor
            size_factors = np.median(ratios, axis=1)

            X_normalized = X_obs / size_factors[:, np.newaxis]

        else:
            raise ValueError(f"Unknown normalization method: {method}")

        return X_normalized

    @staticmethod
    def estimate_size_factors_robust(
        X_obs: np.ndarray,
        highly_variable_genes: Optional[np.ndarray] = None,
    ) -> np.ndarray:
        """
        Robust size factor estimation.

        Excludes highly variable genes that might be biologically different.

        Parameters
        ----------
        X_obs : np.ndarray
            Raw counts
        highly_variable_genes : Optional[np.ndarray]
            Boolean mask of genes to exclude

        Returns
        -------
        size_factors : np.ndarray
            Size factors per cell
        """
        if highly_variable_genes is not None:
            # Use only stable genes
            X_stable = X_obs[:, ~highly_variable_genes]
        else:
            X_stable = X_obs

        # Median-of-ratios on stable genes
        log_means = np.log(X_stable + 1).mean(axis=0)
        geo_means = np.exp(log_means)

        ratios = X_stable / (geo_means + 1e-8)
        size_factors = np.median(ratios, axis=1)

        return size_factors


class DepthModel:
    """
    Model relationship between sequencing depth and detection.

    Deeper sequencing → detect more genes, especially low-abundance.
    """

    @staticmethod
    def saturation_curve(
        X_obs: np.ndarray,
        downsample_fractions: np.ndarray = np.linspace(0.1, 1.0, 10),
    ) -> Dict[str, np.ndarray]:
        """
        Estimate saturation curve: genes detected vs depth.

        Downsamples counts to different depths and measures
        number of detected genes.

        Parameters
        ----------
        X_obs : np.ndarray
            Observed counts (n_cells, n_genes)
        downsample_fractions : np.ndarray
            Fractions of total depth to sample

        Returns
        -------
        saturation_data : Dict[str, np.ndarray]
            Genes detected at each depth
        """
        n_cells, n_genes = X_obs.shape

        genes_detected = np.zeros((n_cells, len(downsample_fractions)))

        for i, frac in enumerate(downsample_fractions):
            # Downsample each cell
            for cell_idx in range(n_cells):
                cell_counts = X_obs[cell_idx, :]
                total = cell_counts.sum()

                # Multinomial downsample
                if total > 0:
                    target_depth = int(total * frac)
                    probs = cell_counts / total

                    downsampled = np.random.multinomial(
                        target_depth,
                        probs
                    )

                    genes_detected[cell_idx, i] = (downsampled > 0).sum()

        return {
            'downsample_fractions': downsample_fractions,
            'genes_detected': genes_detected,
            'mean_genes_detected': genes_detected.mean(axis=0),
        }

    @staticmethod
    def depth_adjustment_factor(
        library_size: np.ndarray,
        target_depth: float,
    ) -> np.ndarray:
        """
        Compute adjustment factors to normalize to target depth.

        Parameters
        ----------
        library_size : np.ndarray
            Observed library sizes (n_cells,)
        target_depth : float
            Target depth to normalize to

        Returns
        -------
        adjustment_factors : np.ndarray
            Multiplicative factors
        """
        return target_depth / (library_size + 1e-8)
