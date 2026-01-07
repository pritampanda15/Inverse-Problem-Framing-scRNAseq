"""
Synthetic data generation for validation.

Ground truth is known, so we can measure how well we recover it.
"""

import numpy as np
import torch
from anndata import AnnData
from typing import Optional, Dict, Tuple
import logging

from ..measurement.operator import MeasurementOperator

logger = logging.getLogger(__name__)


class SyntheticDataGenerator:
    """
    Generate synthetic scRNA-seq data with known ground truth.

    Process:
    1. Generate true biological state Z_true
    2. Apply realistic measurement model
    3. Observe distorted counts X
    4. Later: try to recover Z_true from X

    Parameters
    ----------
    n_cells : int
        Number of cells
    n_genes : int
        Number of genes
    n_programs : int
        Number of true transcriptional programs
    capture_efficiency : float
        Mean capture probability
    dropout_rate : float
        Target dropout rate
    sequencing_depth : int
        Mean library size
    """

    def __init__(
        self,
        n_cells: int = 1000,
        n_genes: int = 2000,
        n_programs: int = 5,
        capture_efficiency: float = 0.1,
        dropout_rate: float = 0.7,
        sequencing_depth: int = 10000,
    ):
        self.n_cells = n_cells
        self.n_genes = n_genes
        self.n_programs = n_programs
        self.capture_efficiency = capture_efficiency
        self.dropout_rate = dropout_rate
        self.sequencing_depth = sequencing_depth

    def generate(
        self,
        seed: Optional[int] = None,
    ) -> Tuple[AnnData, Dict[str, np.ndarray]]:
        """
        Generate synthetic dataset.

        Returns
        -------
        adata : AnnData
            Observed counts (what you'd get from sequencing)
        ground_truth : Dict[str, np.ndarray]
            True biological state (for validation)

        Examples
        --------
        >>> gen = SyntheticDataGenerator(n_cells=500, n_genes=1000)
        >>> adata, truth = gen.generate(seed=42)
        >>> # adata.X is observed (distorted)
        >>> # truth['Z_true'] is ground truth
        """
        if seed is not None:
            np.random.seed(seed)
            torch.manual_seed(seed)

        logger.info("Generating synthetic scRNA-seq data...")

        # 1. Generate true biological state
        Z_true, program_weights, program_signatures = self._generate_biology()

        # 2. Apply measurement model
        X_obs, intermediates = self._apply_measurement(Z_true)

        # 3. Create AnnData
        adata = AnnData(X=X_obs)

        # Add gene and cell names
        adata.obs_names = [f"Cell_{i}" for i in range(self.n_cells)]
        adata.var_names = [f"Gene_{i}" for i in range(self.n_genes)]

        # Add ground truth as obs/var annotations
        for prog_idx in range(self.n_programs):
            adata.obs[f'true_program_{prog_idx}'] = program_weights[:, prog_idx]

        # Store ground truth
        ground_truth = {
            'Z_true': Z_true,
            'program_weights': program_weights,
            'program_signatures': program_signatures,
            'Z_captured': intermediates['Z_captured'],
            'Z_amplified': intermediates['Z_amplified'],
            'capture_prob': intermediates['capture_prob'],
        }

        logger.info(f"Generated: {self.n_cells} cells × {self.n_genes} genes")
        logger.info(f"Dropout rate: {(X_obs == 0).mean():.2%}")
        logger.info(f"Mean library size: {X_obs.sum(axis=1).mean():.0f}")

        return adata, ground_truth

    def _generate_biology(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Generate true biological state.

        Returns
        -------
        Z_true : np.ndarray
            True expression (n_cells, n_genes)
        program_weights : np.ndarray
            Program composition per cell (n_cells, n_programs)
        program_signatures : np.ndarray
            Gene signatures per program (n_programs, n_genes)
        """
        # Create program signatures
        # Each program is a sparse gene expression pattern
        program_signatures = np.zeros((self.n_programs, self.n_genes))

        genes_per_program = self.n_genes // self.n_programs

        for prog_idx in range(self.n_programs):
            # Random subset of genes active in this program
            start_idx = prog_idx * genes_per_program
            end_idx = (prog_idx + 1) * genes_per_program

            # Log-normal expression for program genes
            program_signatures[prog_idx, start_idx:end_idx] = np.random.lognormal(
                mean=2.0,
                sigma=1.0,
                size=genes_per_program,
            )

        # Sample program weights per cell (Dirichlet)
        # Creates compositional data (sums to 1)
        program_weights = np.random.dirichlet(
            alpha=np.ones(self.n_programs) * 2.0,  # Concentration parameter
            size=self.n_cells,
        )

        # Compute true expression: weighted sum of programs
        Z_true = program_weights @ program_signatures

        # Add baseline expression
        baseline = np.random.lognormal(mean=0.5, sigma=0.5, size=self.n_genes)
        Z_true = Z_true + baseline

        return Z_true, program_weights, program_signatures

    def _apply_measurement(
        self,
        Z_true: np.ndarray,
    ) -> Tuple[np.ndarray, Dict[str, np.ndarray]]:
        """
        Apply realistic measurement model.

        Parameters
        ----------
        Z_true : np.ndarray
            True expression

        Returns
        -------
        X_obs : np.ndarray
            Observed counts
        intermediates : Dict[str, np.ndarray]
            Intermediate measurement steps
        """
        # Create measurement operator
        measurement_op = MeasurementOperator(
            n_genes=self.n_genes,
            n_cells=self.n_cells,
        )

        # Set capture efficiency
        with torch.no_grad():
            # Gene-specific capture rates (around target mean)
            capture_logits = np.random.normal(
                loc=np.log(self.capture_efficiency / (1 - self.capture_efficiency)),
                scale=0.5,
                size=self.n_genes,
            )
            measurement_op.capture_model.gene_capture_logit.copy_(
                torch.from_numpy(capture_logits).float()
            )

        # Convert to torch
        Z_torch = torch.from_numpy(Z_true).float()

        # Apply measurement
        X_obs_torch, intermediates_torch = measurement_op.forward(
            Z_torch,
            library_sizes=torch.ones(self.n_cells) * self.sequencing_depth,
            return_intermediates=True,
        )

        # Sample (add stochasticity)
        X_obs_torch = measurement_op.sample_observation(
            Z_torch,
            library_sizes=torch.ones(self.n_cells) * self.sequencing_depth,
        )

        # Convert back to numpy
        X_obs = X_obs_torch.cpu().numpy()

        intermediates = {
            k: v.cpu().numpy() if isinstance(v, torch.Tensor) else v
            for k, v in intermediates_torch.items()
        }

        return X_obs, intermediates


def generate_synthetic_data(
    n_cells: int = 1000,
    n_genes: int = 2000,
    n_programs: int = 5,
    **kwargs,
) -> Tuple[AnnData, Dict[str, np.ndarray]]:
    """
    Convenience function to generate synthetic data.

    Parameters
    ----------
    n_cells : int
        Number of cells
    n_genes : int
        Number of genes
    n_programs : int
        Number of programs
    **kwargs
        Additional arguments for SyntheticDataGenerator

    Returns
    -------
    adata : AnnData
        Observed data
    ground_truth : Dict[str, np.ndarray]
        True biological state

    Examples
    --------
    >>> adata, truth = generate_synthetic_data(n_cells=500, seed=42)
    >>> # Fit model
    >>> isc.pp.fit_inverse_model(adata)
    >>> # Compare to truth
    >>> Z_inferred = adata.obsm['Z_true_mean']
    >>> Z_true = truth['Z_true']
    >>> correlation = np.corrcoef(Z_inferred.flatten(), Z_true.flatten())[0, 1]
    >>> print(f"Recovery correlation: {correlation:.3f}")
    """
    generator = SyntheticDataGenerator(
        n_cells=n_cells,
        n_genes=n_genes,
        n_programs=n_programs,
        **kwargs,
    )

    return generator.generate()


def generate_realistic_benchmark(
    scenario: str = "simple",
) -> Tuple[AnnData, Dict[str, np.ndarray]]:
    """
    Generate benchmark datasets with different difficulty levels.

    Parameters
    ----------
    scenario : str
        'simple': Easy (low dropout, high depth)
        'moderate': Medium difficulty
        'hard': Hard (high dropout, low depth, batch effects)

    Returns
    -------
    adata : AnnData
        Data
    ground_truth : Dict
        Truth
    """
    if scenario == "simple":
        return generate_synthetic_data(
            n_cells=500,
            n_genes=1000,
            n_programs=3,
            capture_efficiency=0.3,
            dropout_rate=0.5,
            sequencing_depth=20000,
            seed=42,
        )

    elif scenario == "moderate":
        return generate_synthetic_data(
            n_cells=1000,
            n_genes=2000,
            n_programs=5,
            capture_efficiency=0.1,
            dropout_rate=0.7,
            sequencing_depth=10000,
            seed=42,
        )

    elif scenario == "hard":
        return generate_synthetic_data(
            n_cells=2000,
            n_genes=3000,
            n_programs=10,
            capture_efficiency=0.05,
            dropout_rate=0.85,
            sequencing_depth=5000,
            seed=42,
        )

    else:
        raise ValueError(f"Unknown scenario: {scenario}")
