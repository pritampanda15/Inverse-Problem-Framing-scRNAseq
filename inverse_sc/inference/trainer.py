"""
Training loop for inverse inference.

Optimizes the ELBO: L = E_q[log p(X, Z)] - E_q[log q(Z)]
"""

import torch
import numpy as np
import pyro
import pyro.optim as optim
from pyro.infer import SVI, Trace_ELBO, TraceMeanField_ELBO
from typing import Optional, Dict, List
from tqdm import tqdm

from .model import InverseModel
from .guide import InferenceGuide


class InverseTrainer:
    """
    Trainer for inverse problem inference.

    Handles:
    - ELBO optimization
    - Mini-batching
    - Learning rate scheduling
    - Convergence monitoring

    Parameters
    ----------
    model : InverseModel
        Generative model
    guide : InferenceGuide
        Variational guide
    learning_rate : float
        Initial learning rate
    use_cuda : bool
        Whether to use GPU
    """

    def __init__(
        self,
        model: InverseModel,
        guide: InferenceGuide,
        learning_rate: float = 1e-3,
        use_cuda: bool = False,
    ):
        self.model = model
        self.guide = guide
        self.learning_rate = learning_rate
        self.device = torch.device("cuda" if use_cuda and torch.cuda.is_available() else "cpu")

        # Move models to device
        self.model.to(self.device)
        self.guide.to(self.device)

        # Optimizer
        self.optimizer = optim.ClippedAdam({
            "lr": learning_rate,
            "clip_norm": 10.0,
        })

        # SVI
        self.svi = SVI(
            model=self.model.model,
            guide=self.guide.guide,
            optim=self.optimizer,
            loss=Trace_ELBO(),
        )

        # Training history
        self.loss_history = []

    def train(
        self,
        X_obs: np.ndarray,
        library_sizes: Optional[np.ndarray] = None,
        n_epochs: int = 100,
        batch_size: Optional[int] = None,
        verbose: bool = True,
    ) -> Dict[str, List[float]]:
        """
        Train the model via ELBO maximization.

        Parameters
        ----------
        X_obs : np.ndarray
            Observed count matrix (n_cells, n_genes)
        library_sizes : Optional[np.ndarray]
            Library sizes (n_cells,)
        n_epochs : int
            Number of training epochs
        batch_size : Optional[int]
            Mini-batch size (None = full batch)
        verbose : bool
            Whether to show progress bar

        Returns
        -------
        history : Dict[str, List[float]]
            Training history
        """
        n_cells, n_genes = X_obs.shape

        # Convert to tensors
        X_tensor = torch.from_numpy(X_obs).float().to(self.device)

        if library_sizes is not None:
            lib_tensor = torch.from_numpy(library_sizes).float().to(self.device)
        else:
            lib_tensor = None

        # Determine batch size
        if batch_size is None:
            batch_size = n_cells

        # Training loop
        iterator = range(n_epochs)
        if verbose:
            iterator = tqdm(iterator, desc="Training")

        for epoch in iterator:
            epoch_loss = 0.0
            n_batches = 0

            # Mini-batch training
            indices = torch.randperm(n_cells)

            for i in range(0, n_cells, batch_size):
                batch_indices = indices[i:i+batch_size]

                X_batch = X_tensor[batch_indices]
                lib_batch = lib_tensor[batch_indices] if lib_tensor is not None else None

                # SVI step - pass batch_indices for cell-specific parameters
                loss = self.svi.step(X_batch, lib_batch, batch_indices)

                epoch_loss += loss
                n_batches += 1

            # Average loss
            avg_loss = epoch_loss / n_batches
            self.loss_history.append(avg_loss)

            if verbose and epoch % 10 == 0:
                iterator.set_postfix({"ELBO": -avg_loss})

        return {"loss": self.loss_history}

    def predict(
        self,
        X_obs: np.ndarray,
        library_sizes: Optional[np.ndarray] = None,
        n_samples: int = 100,
    ) -> Dict[str, np.ndarray]:
        """
        Predict posterior over Z_true given X_obs.

        Parameters
        ----------
        X_obs : np.ndarray
            Observed counts
        library_sizes : Optional[np.ndarray]
            Library sizes
        n_samples : int
            Number of posterior samples

        Returns
        -------
        predictions : Dict[str, np.ndarray]
            Posterior mean and std for Z_true and programs
        """
        from pyro.infer import Predictive

        # Convert to tensors
        n_cells = X_obs.shape[0]
        X_tensor = torch.from_numpy(X_obs).float().to(self.device)
        lib_tensor = None
        if library_sizes is not None:
            lib_tensor = torch.from_numpy(library_sizes).float().to(self.device)

        # Sample from posterior
        predictive = Predictive(
            self.guide.guide,
            num_samples=n_samples,
        )

        # Create cell indices for full dataset
        cell_indices = torch.arange(n_cells).to(self.device)

        with torch.no_grad():
            samples = predictive(X_tensor, lib_tensor, cell_indices)

        # Extract samples
        Z_true_samples = samples['Z_true'].cpu().numpy()  # (n_samples, n_cells, n_genes)
        program_samples = samples['program_weights'].cpu().numpy()  # (n_samples, n_cells, n_programs)

        # Compute statistics
        predictions = {
            'Z_true_mean': Z_true_samples.mean(axis=0),
            'Z_true_std': Z_true_samples.std(axis=0),
            'Z_true_median': np.median(Z_true_samples, axis=0),
            'program_weights_mean': program_samples.mean(axis=0),
            'program_weights_std': program_samples.std(axis=0),
        }

        return predictions

    def evaluate_reconstruction(
        self,
        X_obs: np.ndarray,
        library_sizes: Optional[np.ndarray] = None,
    ) -> Dict[str, float]:
        """
        Evaluate reconstruction quality.

        Measures how well the model reconstructs X from inferred Z.

        Parameters
        ----------
        X_obs : np.ndarray
            Observed counts
        library_sizes : Optional[np.ndarray]
            Library sizes

        Returns
        -------
        metrics : Dict[str, float]
            Reconstruction metrics
        """
        # Get posterior predictions
        predictions = self.predict(X_obs, library_sizes, n_samples=100)
        Z_true_mean = predictions['Z_true_mean']

        # Forward through measurement model
        Z_tensor = torch.from_numpy(Z_true_mean).float().to(self.device)
        lib_tensor = None
        if library_sizes is not None:
            lib_tensor = torch.from_numpy(library_sizes).float().to(self.device)

        with torch.no_grad():
            X_reconstructed = self.model.measurement_operator(Z_tensor, lib_tensor)
            X_reconstructed = X_reconstructed.cpu().numpy()

        # Compute metrics
        # 1. Mean Squared Error (log scale)
        mse_log = np.mean((np.log1p(X_obs) - np.log1p(X_reconstructed)) ** 2)

        # 2. Pearson correlation (per cell)
        from scipy.stats import pearsonr
        correlations = [
            pearsonr(X_obs[i], X_reconstructed[i])[0]
            for i in range(X_obs.shape[0])
        ]
        mean_correlation = np.mean(correlations)

        # 3. Reconstruction error (counts)
        mae = np.mean(np.abs(X_obs - X_reconstructed))

        return {
            'mse_log': mse_log,
            'mean_correlation': mean_correlation,
            'mae': mae,
        }

    def save_checkpoint(
        self,
        path: str,
    ):
        """
        Save model checkpoint.

        Parameters
        ----------
        path : str
            Path to save checkpoint
        """
        checkpoint = {
            'model_state_dict': self.model.state_dict(),
            'guide_state_dict': self.guide.state_dict(),
            'optimizer_state_dict': self.optimizer.get_state(),
            'loss_history': self.loss_history,
        }
        torch.save(checkpoint, path)

    def load_checkpoint(
        self,
        path: str,
    ):
        """
        Load model checkpoint.

        Parameters
        ----------
        path : str
            Path to checkpoint
        """
        checkpoint = torch.load(path, map_location=self.device)
        self.model.load_state_dict(checkpoint['model_state_dict'])
        self.guide.load_state_dict(checkpoint['guide_state_dict'])
        self.optimizer.set_state(checkpoint['optimizer_state_dict'])
        self.loss_history = checkpoint['loss_history']


class EarlyStopping:
    """
    Early stopping to prevent overfitting.

    Monitors validation loss and stops training when it stops improving.

    Parameters
    ----------
    patience : int
        Number of epochs to wait before stopping
    min_delta : float
        Minimum change to qualify as improvement
    """

    def __init__(
        self,
        patience: int = 10,
        min_delta: float = 1e-4,
    ):
        self.patience = patience
        self.min_delta = min_delta
        self.best_loss = float('inf')
        self.counter = 0
        self.should_stop = False

    def __call__(
        self,
        val_loss: float,
    ) -> bool:
        """
        Check if training should stop.

        Parameters
        ----------
        val_loss : float
            Current validation loss

        Returns
        -------
        should_stop : bool
            Whether to stop training
        """
        if val_loss < self.best_loss - self.min_delta:
            self.best_loss = val_loss
            self.counter = 0
        else:
            self.counter += 1

        if self.counter >= self.patience:
            self.should_stop = True

        return self.should_stop
