"""
Variational guide (inference network) for posterior approximation.

The guide q(Z | X) approximates the true posterior p(Z | X).
"""

import torch
import torch.nn as nn
import pyro
import pyro.distributions as dist
from typing import Optional


class InferenceGuide(nn.Module):
    """
    Variational guide: encoder network q(Z | X).

    Implements amortized variational inference:
    Neural network maps X → parameters of q(Z).

    Parameters
    ----------
    n_genes : int
        Number of genes (input dimension)
    n_latent : int
        Latent dimension
    n_programs : int
        Number of transcriptional programs
    hidden_dims : list
        Hidden layer sizes
    """

    def __init__(
        self,
        n_genes: int,
        n_latent: int = 30,
        n_programs: int = 20,
        hidden_dims: list = [128, 64],
    ):
        super().__init__()
        self.n_genes = n_genes
        self.n_latent = n_latent
        self.n_programs = n_programs

        # Encoder network: X → latent representation
        layers = []
        in_dim = n_genes

        for hidden_dim in hidden_dims:
            layers.extend([
                nn.Linear(in_dim, hidden_dim),
                nn.BatchNorm1d(hidden_dim),
                nn.ReLU(),
                nn.Dropout(0.1),
            ])
            in_dim = hidden_dim

        self.encoder = nn.Sequential(*layers)

        # Output heads
        # 1. Program weights (Dirichlet parameters)
        self.program_head = nn.Sequential(
            nn.Linear(in_dim, n_programs),
            nn.Softplus(),  # Ensure positive (Dirichlet concentration)
        )

        # 2. True expression Z (LogNormal parameters)
        self.z_mean_head = nn.Linear(in_dim, n_genes)
        self.z_scale_head = nn.Sequential(
            nn.Linear(in_dim, n_genes),
            nn.Softplus(),  # Ensure positive scale
        )

    def guide(
        self,
        X_obs: Optional[torch.Tensor] = None,
        library_sizes: Optional[torch.Tensor] = None,
        cell_indices: Optional[torch.Tensor] = None,
    ):
        """
        Variational guide q(Z | X).

        Parameters
        ----------
        X_obs : torch.Tensor
            Observed counts (n_cells, n_genes)
        library_sizes : Optional[torch.Tensor]
            Library sizes (not used here, but must match model signature)
        cell_indices : Optional[torch.Tensor]
            Cell indices for mini-batching
        """
        # Register module
        pyro.module("guide", self)

        # Normalize input (log1p transform + library size normalization)
        if X_obs is None:
            raise ValueError("Guide requires observed data X_obs")

        X_normalized = self._normalize_input(X_obs, library_sizes)

        n_cells = X_obs.shape[0]

        with pyro.plate("cells", n_cells, subsample=cell_indices):
            # Encode X
            hidden = self.encoder(X_normalized)

            # ===== VARIATIONAL POSTERIOR q(program_weights | X) =====
            program_concentration = self.program_head(hidden)

            # Add small constant for numerical stability
            program_concentration = program_concentration + 0.1

            program_weights = pyro.sample(
                "program_weights",
                dist.Dirichlet(program_concentration),
            )

            # ===== VARIATIONAL POSTERIOR q(Z_true | X) =====
            z_loc = self.z_mean_head(hidden)
            z_scale = self.z_scale_head(hidden) + 0.01  # Minimum scale

            Z_true = pyro.sample(
                "Z_true",
                dist.LogNormal(z_loc, z_scale).to_event(1),
            )

    def _normalize_input(
        self,
        X_obs: torch.Tensor,
        library_sizes: Optional[torch.Tensor] = None,
    ) -> torch.Tensor:
        """
        Normalize input counts for encoder.

        Parameters
        ----------
        X_obs : torch.Tensor
            Raw counts
        library_sizes : Optional[torch.Tensor]
            Library sizes

        Returns
        -------
        X_normalized : torch.Tensor
            Normalized counts
        """
        if library_sizes is None:
            library_sizes = X_obs.sum(dim=1, keepdim=True)
        elif library_sizes.dim() == 1:
            library_sizes = library_sizes.unsqueeze(1)

        # Log1p + library size normalization
        X_normalized = torch.log1p(
            X_obs / (library_sizes + 1e-8) * 10000
        )

        return X_normalized

    def forward(
        self,
        X_obs: torch.Tensor,
        library_sizes: Optional[torch.Tensor] = None,
    ):
        """
        Forward pass (just calls guide).
        """
        return self.guide(X_obs, library_sizes)


class SimpleGuide(nn.Module):
    """
    Simpler mean-field guide without neural network.

    Uses per-cell, per-gene variational parameters.
    Less flexible but easier to interpret.

    Parameters
    ----------
    n_genes : int
        Number of genes
    n_cells : int
        Number of cells
    n_programs : int
        Number of programs
    """

    def __init__(
        self,
        n_genes: int,
        n_cells: int,
        n_programs: int = 20,
    ):
        super().__init__()
        self.n_genes = n_genes
        self.n_cells = n_cells
        self.n_programs = n_programs

        # Variational parameters (initialized randomly)
        # Program weights
        self.program_concentration = nn.Parameter(
            torch.ones(n_cells, n_programs)
        )

        # True expression
        self.z_loc = nn.Parameter(
            torch.randn(n_cells, n_genes)
        )
        self.z_scale_logit = nn.Parameter(
            torch.zeros(n_cells, n_genes)
        )

    def guide(
        self,
        X_obs: Optional[torch.Tensor] = None,
        library_sizes: Optional[torch.Tensor] = None,
        cell_indices: Optional[torch.Tensor] = None,
    ):
        """
        Mean-field guide.
        """
        pyro.module("simple_guide", self)

        if X_obs is not None:
            n_cells = X_obs.shape[0]
        else:
            n_cells = self.n_cells

        with pyro.plate("cells", n_cells, subsample=cell_indices):
            # Program weights
            if cell_indices is not None:
                conc = self.program_concentration[cell_indices]
                z_loc = self.z_loc[cell_indices]
                z_scale_logit = self.z_scale_logit[cell_indices]
            else:
                conc = self.program_concentration
                z_loc = self.z_loc
                z_scale_logit = self.z_scale_logit

            program_weights = pyro.sample(
                "program_weights",
                dist.Dirichlet(torch.softplus(conc) + 0.1),
            )

            # True expression
            z_scale = torch.sigmoid(z_scale_logit) * 2.0 + 0.01
            Z_true = pyro.sample(
                "Z_true",
                dist.LogNormal(z_loc, z_scale).to_event(1),
            )


class AutoGuideWrapper:
    """
    Wrapper for Pyro's automatic guides.

    Useful for quick prototyping without defining guide manually.
    """

    @staticmethod
    def create_auto_guide(
        model,
        guide_type: str = "normal",
    ):
        """
        Create automatic guide.

        Parameters
        ----------
        model : InverseModel
            The generative model
        guide_type : str
            'normal' (mean-field), 'mvn' (full covariance), or 'iaf' (flow)

        Returns
        -------
        guide : pyro.infer.autoguide
            Automatic guide
        """
        if guide_type == "normal":
            from pyro.infer.autoguide import AutoDiagonalNormal
            return AutoDiagonalNormal(model.model)

        elif guide_type == "mvn":
            from pyro.infer.autoguide import AutoMultivariateNormal
            return AutoMultivariateNormal(model.model)

        elif guide_type == "iaf":
            from pyro.infer.autoguide import AutoIAFNormal
            return AutoIAFNormal(model.model)

        else:
            raise ValueError(f"Unknown guide type: {guide_type}")
