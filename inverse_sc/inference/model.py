"""
Generative model for inverse problem.

Defines p(Z, X) = p(Z) * p(X | Z)

where:
- p(Z) is the biological prior
- p(X | Z) is the measurement model
"""

import torch
import torch.nn as nn
import pyro
import pyro.distributions as dist
from typing import Optional, Dict

from ..measurement.operator import MeasurementOperator


class InverseModel(nn.Module):
    """
    Full generative model for scRNA-seq inverse problem.

    Generative process:
    1. Sample biological state Z from prior p(Z)
    2. Apply measurement model: X = M(Z) + noise

    Parameters
    ----------
    n_genes : int
        Number of genes
    n_cells : int
        Number of cells
    n_latent : int
        Dimensionality of latent biological space
    n_programs : int
        Number of transcriptional programs
    measurement_operator : MeasurementOperator
        The forward measurement model
    """

    def __init__(
        self,
        n_genes: int,
        n_cells: int,
        n_latent: int = 30,
        n_programs: int = 20,
        measurement_operator: Optional[MeasurementOperator] = None,
    ):
        super().__init__()
        self.n_genes = n_genes
        self.n_cells = n_cells
        self.n_latent = n_latent
        self.n_programs = n_programs

        # Measurement operator
        if measurement_operator is None:
            self.measurement_operator = MeasurementOperator(
                n_genes=n_genes,
                n_cells=n_cells,
            )
        else:
            self.measurement_operator = measurement_operator

        # Biological prior: transcriptional programs
        # Each program is a gene expression signature
        self.program_signatures = nn.Parameter(
            torch.randn(n_programs, n_genes) * 0.1
        )

        # Baseline expression (independent of programs)
        self.baseline_expression = nn.Parameter(
            torch.zeros(n_genes)
        )

    def model(
        self,
        X_obs: Optional[torch.Tensor] = None,
        library_sizes: Optional[torch.Tensor] = None,
        cell_indices: Optional[torch.Tensor] = None,
    ):
        """
        Generative model p(Z, X).

        Parameters
        ----------
        X_obs : Optional[torch.Tensor]
            Observed counts (n_cells, n_genes). If provided, condition on it.
        library_sizes : Optional[torch.Tensor]
            Observed library sizes (n_cells,)
        cell_indices : Optional[torch.Tensor]
            Which cells to model (for mini-batching)
        """
        # Determine batch size
        if X_obs is not None:
            n_cells = X_obs.shape[0]
        elif cell_indices is not None:
            n_cells = len(cell_indices)
        else:
            n_cells = self.n_cells

        # Register modules with pyro
        pyro.module("inverse_model", self)

        with pyro.plate("cells", n_cells, subsample=cell_indices):
            # ===== BIOLOGICAL PRIOR p(Z) =====

            # 1. Sample program weights for each cell
            # Dirichlet ensures they sum to 1 (compositional)
            program_weights = pyro.sample(
                "program_weights",
                dist.Dirichlet(torch.ones(self.n_programs)),
            )

            # 2. Compute expected expression from programs
            # Z_mean = weighted sum of program signatures
            program_expression = torch.matmul(
                program_weights,
                self.program_signatures
            )

            # Add baseline
            Z_mean = torch.exp(self.baseline_expression + program_expression)

            # 3. Sample true expression with biological variability
            # LogNormal captures multiplicative noise
            Z_true = pyro.sample(
                "Z_true",
                dist.LogNormal(
                    loc=torch.log(Z_mean + 1e-8),
                    scale=0.5,  # Biological variability
                ).to_event(1)
            )

            # ===== MEASUREMENT MODEL p(X | Z) =====

            # Apply measurement operator
            X_expected = self.measurement_operator(
                Z_true,
                library_sizes=library_sizes,
                cell_indices=cell_indices,
            )

            # Observation model: X ~ Poisson(X_expected)
            # (or could use NegativeBinomial for overdispersion)
            pyro.sample(
                "X_obs",
                dist.Poisson(X_expected + 1e-8).to_event(1),
                obs=X_obs,
            )

    def forward(
        self,
        X_obs: torch.Tensor,
        library_sizes: Optional[torch.Tensor] = None,
    ):
        """
        Forward pass (just calls model).
        """
        return self.model(X_obs, library_sizes)

    def sample_prior(
        self,
        n_samples: int = 100,
    ) -> Dict[str, torch.Tensor]:
        """
        Sample from prior p(Z) to see what biology looks like
        before seeing data.

        Parameters
        ----------
        n_samples : int
            Number of cells to sample

        Returns
        -------
        samples : Dict[str, torch.Tensor]
            Samples from prior
        """
        samples = {}

        # Sample program weights
        program_weights = dist.Dirichlet(
            torch.ones(self.n_programs)
        ).sample((n_samples,))

        # Compute expression
        program_expression = torch.matmul(
            program_weights,
            self.program_signatures
        )
        Z_mean = torch.exp(self.baseline_expression + program_expression)

        # Sample with biological noise
        Z_true = dist.LogNormal(
            torch.log(Z_mean + 1e-8),
            0.5,
        ).sample()

        samples['program_weights'] = program_weights
        samples['Z_true'] = Z_true

        return samples


class BiologicalPrior(nn.Module):
    """
    Alternative prior models for biology.

    The default uses transcriptional programs, but you could use:
    - Factor models (PCA-like)
    - Gene regulatory networks
    - Pathway-based priors
    """

    def __init__(
        self,
        n_genes: int,
        prior_type: str = "programs",
        **kwargs,
    ):
        super().__init__()
        self.n_genes = n_genes
        self.prior_type = prior_type

        if prior_type == "programs":
            self.n_programs = kwargs.get('n_programs', 20)
            self.program_signatures = nn.Parameter(
                torch.randn(self.n_programs, n_genes) * 0.1
            )

        elif prior_type == "factor":
            self.n_factors = kwargs.get('n_factors', 30)
            self.factor_loadings = nn.Parameter(
                torch.randn(n_genes, self.n_factors) * 0.1
            )

        else:
            raise ValueError(f"Unknown prior type: {prior_type}")

    def sample(
        self,
        n_cells: int,
    ) -> torch.Tensor:
        """
        Sample Z from prior.

        Parameters
        ----------
        n_cells : int
            Number of cells

        Returns
        -------
        Z_true : torch.Tensor
            Sampled true expression
        """
        if self.prior_type == "programs":
            # Dirichlet composition
            program_weights = dist.Dirichlet(
                torch.ones(self.n_programs)
            ).sample((n_cells,))

            program_expr = torch.matmul(program_weights, self.program_signatures)
            Z_true = torch.exp(program_expr)

        elif self.prior_type == "factor":
            # Factor model: Z = W @ h, where h ~ Normal(0, 1)
            factors = dist.Normal(0, 1).sample((n_cells, self.n_factors))
            Z_true = torch.exp(torch.matmul(factors, self.factor_loadings.T))

        return Z_true


class MeasurementLikelihood:
    """
    Different likelihood models for p(X | Z, M).

    Options:
    - Poisson: Simple, assumes variance = mean
    - Negative Binomial: Overdispersion
    - Zero-Inflated: Explicit dropout modeling
    """

    @staticmethod
    def poisson_likelihood(
        X_obs: torch.Tensor,
        X_expected: torch.Tensor,
    ):
        """
        Poisson observation model.

        Parameters
        ----------
        X_obs : torch.Tensor
            Observed counts
        X_expected : torch.Tensor
            Expected counts from measurement model
        """
        return pyro.sample(
            "X_obs",
            dist.Poisson(X_expected + 1e-8).to_event(1),
            obs=X_obs,
        )

    @staticmethod
    def negative_binomial_likelihood(
        X_obs: torch.Tensor,
        X_expected: torch.Tensor,
        dispersion: torch.Tensor,
    ):
        """
        Negative Binomial observation model.

        Allows variance > mean (overdispersion).

        Parameters
        ----------
        X_obs : torch.Tensor
            Observed counts
        X_expected : torch.Tensor
            Expected mean
        dispersion : torch.Tensor
            Dispersion parameter (higher = more variance)
        """
        # NB parameterization: mean and concentration
        # var = mean + mean^2 / concentration
        return pyro.sample(
            "X_obs",
            dist.NegativeBinomial(
                total_count=dispersion,
                probs=dispersion / (dispersion + X_expected + 1e-8),
            ).to_event(1),
            obs=X_obs,
        )

    @staticmethod
    def zero_inflated_negative_binomial(
        X_obs: torch.Tensor,
        X_expected: torch.Tensor,
        dispersion: torch.Tensor,
        dropout_prob: torch.Tensor,
    ):
        """
        Zero-Inflated Negative Binomial (ZINB).

        Used by scVI. Mixture of:
        - Point mass at 0 (dropout)
        - Negative Binomial (biological + technical variation)

        Parameters
        ----------
        X_obs : torch.Tensor
            Observed counts
        X_expected : torch.Tensor
            Expected mean (when not dropout)
        dispersion : torch.Tensor
            Dispersion
        dropout_prob : torch.Tensor
            Probability of dropout zero
        """
        # This is more complex; simplified version:
        # In practice, use a mixture distribution

        # For now, return standard NB
        # (Full ZINB implementation would use dist.ZeroInflatedNegativeBinomial
        # if available, or manual mixture)
        return MeasurementLikelihood.negative_binomial_likelihood(
            X_obs,
            X_expected,
            dispersion,
        )
