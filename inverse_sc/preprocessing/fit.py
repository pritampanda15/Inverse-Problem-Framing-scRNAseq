"""
Main fitting function: Scanpy-compatible interface.

This is the primary user-facing function for inverse inference.
"""

import numpy as np
import torch
from anndata import AnnData
from typing import Optional, Dict
import logging

from ..measurement.operator import MeasurementOperator
from ..inference.model import InverseModel
from ..inference.guide import InferenceGuide
from ..inference.trainer import InverseTrainer

logger = logging.getLogger(__name__)


def fit_inverse_model(
    adata: AnnData,
    n_latent: int = 30,
    n_programs: int = 20,
    n_epochs: int = 100,
    batch_size: Optional[int] = 256,
    learning_rate: float = 1e-3,
    use_cuda: bool = False,
    layer: Optional[str] = None,
    copy: bool = False,
    **kwargs,
) -> Optional[AnnData]:
    """
    Fit inverse problem model to scRNA-seq data.

    This is the main entry point. After fitting, adata.obsm will contain:
    - 'Z_true_mean': Inferred true expression (posterior mean)
    - 'Z_true_std': Uncertainty (posterior std)
    - 'program_weights': Cell state composition

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix (Scanpy format)
    n_latent : int
        Latent dimension for inference network
    n_programs : int
        Number of transcriptional programs
    n_epochs : int
        Training epochs
    batch_size : Optional[int]
        Mini-batch size (None = full batch)
    learning_rate : float
        Learning rate for optimization
    use_cuda : bool
        Use GPU if available
    layer : Optional[str]
        Which layer to use (None = .X, 'raw' = .raw.X)
    copy : bool
        Return copy instead of modifying in-place
    **kwargs
        Additional arguments for measurement model

    Returns
    -------
    adata : Optional[AnnData]
        If copy=True, returns modified copy. Otherwise modifies in-place.

    Examples
    --------
    >>> import scanpy as sc
    >>> import inverse_sc as isc
    >>> adata = sc.read_h5ad("data.h5ad")
    >>> isc.pp.fit_inverse_model(adata, n_epochs=200)
    >>> # Now adata.obsm['Z_true_mean'] contains inferred expression
    >>> sc.pp.neighbors(adata, use_rep='Z_true_mean')
    >>> sc.tl.umap(adata)
    """
    logger.info("Fitting inverse problem model...")

    if copy:
        adata = adata.copy()

    # Extract count matrix
    if layer is None:
        X = adata.X
    elif layer == 'raw':
        X = adata.raw.X
    else:
        X = adata.layers[layer]

    # Convert to dense array if sparse
    if hasattr(X, 'toarray'):
        X = X.toarray()

    X = np.asarray(X)
    n_cells, n_genes = X.shape

    logger.info(f"Data shape: {n_cells} cells × {n_genes} genes")

    # Library sizes
    library_sizes = X.sum(axis=1)

    # Initialize measurement operator
    logger.info("Initializing measurement operator...")
    measurement_operator = MeasurementOperator(
        n_genes=n_genes,
        n_cells=n_cells,
        **kwargs,
    )

    # Calibrate from data
    logger.info("Calibrating measurement model from dropout patterns...")
    calibration = measurement_operator.calibrate(X)

    # Store calibration results
    adata.var['capture_prob_estimate'] = calibration['capture']['capture_prob_estimate']
    adata.var['dropout_rate'] = calibration['capture']['dropout_rate']

    # Initialize generative model
    logger.info("Initializing generative model...")
    model = InverseModel(
        n_genes=n_genes,
        n_cells=n_cells,
        n_latent=n_latent,
        n_programs=n_programs,
        measurement_operator=measurement_operator,
    )

    # Initialize guide (inference network)
    logger.info("Initializing inference network...")
    guide = InferenceGuide(
        n_genes=n_genes,
        n_latent=n_latent,
        n_programs=n_programs,
    )

    # Initialize trainer
    trainer = InverseTrainer(
        model=model,
        guide=guide,
        learning_rate=learning_rate,
        use_cuda=use_cuda,
    )

    # Train
    logger.info(f"Training for {n_epochs} epochs...")
    history = trainer.train(
        X_obs=X,
        library_sizes=library_sizes,
        n_epochs=n_epochs,
        batch_size=batch_size,
        verbose=True,
    )

    # Predict posterior
    logger.info("Computing posterior predictions...")
    predictions = trainer.predict(
        X_obs=X,
        library_sizes=library_sizes,
        n_samples=100,
    )

    # Store results in adata
    adata.obsm['Z_true_mean'] = predictions['Z_true_mean']
    adata.obsm['Z_true_std'] = predictions['Z_true_std']
    adata.obsm['Z_true_median'] = predictions['Z_true_median']
    adata.obsm['program_weights'] = predictions['program_weights_mean']
    adata.obsm['program_weights_std'] = predictions['program_weights_std']

    # Store training history
    adata.uns['inverse_model'] = {
        'loss_history': history['loss'],
        'n_epochs': n_epochs,
        'n_latent': n_latent,
        'n_programs': n_programs,
    }

    # Store model for later use
    adata.uns['inverse_model_obj'] = {
        'model': model,
        'guide': guide,
        'trainer': trainer,
    }

    logger.info("Done!")

    # Evaluate reconstruction
    metrics = trainer.evaluate_reconstruction(X, library_sizes)
    logger.info(f"Reconstruction metrics: {metrics}")
    adata.uns['inverse_model']['reconstruction_metrics'] = metrics

    if copy:
        return adata


def calibrate_measurement_model(
    adata: AnnData,
    spike_in_genes: Optional[list] = None,
    spike_in_amounts: Optional[np.ndarray] = None,
    layer: Optional[str] = None,
) -> Dict[str, np.ndarray]:
    """
    Calibrate measurement model from data.

    Uses dropout patterns and (optionally) spike-ins to estimate
    capture and amplification parameters.

    Parameters
    ----------
    adata : AnnData
        Annotated data
    spike_in_genes : Optional[list]
        Names of spike-in genes (e.g., ERCC)
    spike_in_amounts : Optional[np.ndarray]
        Known input amounts for spike-ins
    layer : Optional[str]
        Which layer to use

    Returns
    -------
    calibration_results : Dict[str, np.ndarray]
        Estimated measurement parameters

    Examples
    --------
    >>> calib = isc.pp.calibrate_measurement_model(adata)
    >>> adata.var['capture_efficiency'] = calib['capture_prob']
    """
    # Extract count matrix
    if layer is None:
        X = adata.X
    else:
        X = adata.layers[layer]

    if hasattr(X, 'toarray'):
        X = X.toarray()

    X = np.asarray(X)

    # Initialize measurement operator
    n_cells, n_genes = X.shape
    measurement_operator = MeasurementOperator(
        n_genes=n_genes,
        n_cells=n_cells,
    )

    # Calibrate
    calibration = measurement_operator.calibrate(X)

    # If spike-ins provided, calibrate amplification
    if spike_in_genes is not None and spike_in_amounts is not None:
        spike_in_mask = adata.var_names.isin(spike_in_genes)
        spike_in_counts = X[:, spike_in_mask]

        amp_calib = measurement_operator.amplification_model.calibrate_from_spike_ins(
            spike_in_counts=spike_in_counts,
            spike_in_true_amounts=spike_in_amounts,
        )

        calibration['amplification'] = amp_calib

    return calibration


def quick_fit(
    adata: AnnData,
    **kwargs,
) -> AnnData:
    """
    Quick fit with sensible defaults for exploration.

    Parameters
    ----------
    adata : AnnData
        Data
    **kwargs
        Override defaults

    Returns
    -------
    adata : AnnData
        Modified in-place
    """
    defaults = {
        'n_latent': 30,
        'n_programs': 15,
        'n_epochs': 50,
        'batch_size': 256,
        'learning_rate': 1e-3,
    }

    defaults.update(kwargs)

    return fit_inverse_model(adata, **defaults)
