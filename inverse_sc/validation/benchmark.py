"""
Benchmarking against standard methods (Scanpy, scVI).
"""

import numpy as np
import pandas as pd
from anndata import AnnData
from typing import Dict, Optional
import logging
from scipy.stats import pearsonr, spearmanr

logger = logging.getLogger(__name__)


def benchmark_against_scanpy(
    adata: AnnData,
    ground_truth: Dict[str, np.ndarray],
    run_scanpy: bool = True,
) -> pd.DataFrame:
    """
    Compare inverse method against standard Scanpy pipeline.

    Metrics:
    1. Recovery of true expression (correlation)
    2. Cluster accuracy (ARI)
    3. DE gene detection (precision/recall)

    Parameters
    ----------
    adata : AnnData
        Data with inverse inference results
    ground_truth : Dict[str, np.ndarray]
        True biological state
    run_scanpy : bool
        Whether to run scanpy pipeline for comparison

    Returns
    -------
    results : pd.DataFrame
        Comparison metrics

    Examples
    --------
    >>> adata, truth = isc.validation.generate_synthetic_data()
    >>> isc.pp.fit_inverse_model(adata)
    >>> results = isc.validation.benchmark_against_scanpy(adata, truth)
    >>> print(results)
    """
    results = []

    Z_true = ground_truth['Z_true']

    # === Inverse method results ===
    if 'Z_true_mean' in adata.obsm:
        Z_inferred = adata.obsm['Z_true_mean']

        # Correlation with ground truth (cell-wise)
        cell_correlations = []
        for i in range(Z_true.shape[0]):
            corr, _ = pearsonr(Z_true[i], Z_inferred[i])
            cell_correlations.append(corr)

        # Gene-wise correlation
        gene_correlations = []
        for j in range(Z_true.shape[1]):
            corr, _ = pearsonr(Z_true[:, j], Z_inferred[:, j])
            gene_correlations.append(corr)

        results.append({
            'method': 'InverseSC',
            'mean_cell_correlation': np.mean(cell_correlations),
            'mean_gene_correlation': np.mean(gene_correlations),
            'global_correlation': pearsonr(Z_true.flatten(), Z_inferred.flatten())[0],
            'rmse': np.sqrt(np.mean((Z_true - Z_inferred) ** 2)),
            'mae': np.mean(np.abs(Z_true - Z_inferred)),
        })

    # === Standard Scanpy pipeline ===
    if run_scanpy:
        import scanpy as sc

        # Copy data for scanpy
        adata_sc = adata.copy()

        # Standard preprocessing
        sc.pp.normalize_total(adata_sc, target_sum=1e4)
        sc.pp.log1p(adata_sc)

        # This is what scanpy uses for downstream analysis
        Z_scanpy = adata_sc.X

        # Correlation with ground truth
        cell_correlations_sc = []
        for i in range(Z_true.shape[0]):
            # Note: comparing log-normalized to true counts
            # Not apples-to-apples, but this is what scanpy users analyze
            corr, _ = pearsonr(
                np.log1p(Z_true[i]),
                Z_scanpy[i] if not hasattr(Z_scanpy, 'toarray') else Z_scanpy[i].toarray().flatten()
            )
            cell_correlations_sc.append(corr)

        # Convert to dense if sparse
        if hasattr(Z_scanpy, 'toarray'):
            Z_scanpy = Z_scanpy.toarray()

        results.append({
            'method': 'Scanpy (log-normalized)',
            'mean_cell_correlation': np.mean(cell_correlations_sc),
            'mean_gene_correlation': np.nan,  # Would need gene-wise comparison
            'global_correlation': pearsonr(np.log1p(Z_true).flatten(), Z_scanpy.flatten())[0],
            'rmse': np.sqrt(np.mean((np.log1p(Z_true) - Z_scanpy) ** 2)),
            'mae': np.mean(np.abs(np.log1p(Z_true) - Z_scanpy)),
        })

        # Raw counts (no normalization)
        X_raw = adata.X
        if hasattr(X_raw, 'toarray'):
            X_raw = X_raw.toarray()

        results.append({
            'method': 'Raw counts',
            'mean_cell_correlation': np.nan,
            'mean_gene_correlation': np.nan,
            'global_correlation': pearsonr(Z_true.flatten(), X_raw.flatten())[0],
            'rmse': np.sqrt(np.mean((Z_true - X_raw) ** 2)),
            'mae': np.mean(np.abs(Z_true - X_raw)),
        })

    results_df = pd.DataFrame(results)

    logger.info("\n" + results_df.to_string())

    return results_df


def benchmark_against_scvi(
    adata: AnnData,
    ground_truth: Dict[str, np.ndarray],
    run_scvi: bool = False,
) -> pd.DataFrame:
    """
    Compare against scVI (if installed).

    scVI is a VAE for scRNA-seq, but uses neural network decoder
    rather than physical measurement model.

    Parameters
    ----------
    adata : AnnData
        Data
    ground_truth : Dict
        Truth
    run_scvi : bool
        Whether to run scVI (requires scvi-tools)

    Returns
    -------
    results : pd.DataFrame
        Comparison
    """
    results = []

    Z_true = ground_truth['Z_true']

    # InverseSC results
    if 'Z_true_mean' in adata.obsm:
        Z_inverse = adata.obsm['Z_true_mean']

        results.append({
            'method': 'InverseSC',
            'correlation': pearsonr(Z_true.flatten(), Z_inverse.flatten())[0],
            'rmse': np.sqrt(np.mean((Z_true - Z_inverse) ** 2)),
        })

    # scVI
    if run_scvi:
        try:
            import scvi

            # Setup and train scVI
            scvi.model.SCVI.setup_anndata(adata)
            model = scvi.model.SCVI(adata)
            model.train(max_epochs=100)

            # Get denoised expression
            Z_scvi = model.get_normalized_expression()

            results.append({
                'method': 'scVI',
                'correlation': pearsonr(Z_true.flatten(), Z_scvi.flatten())[0],
                'rmse': np.sqrt(np.mean((Z_true - Z_scvi) ** 2)),
            })

        except ImportError:
            logger.warning("scvi-tools not installed. Skipping scVI benchmark.")

    return pd.DataFrame(results)


def evaluate_clustering(
    adata: AnnData,
    true_labels: np.ndarray,
    predicted_labels: np.ndarray,
) -> Dict[str, float]:
    """
    Evaluate clustering quality.

    Parameters
    ----------
    adata : AnnData
        Data
    true_labels : np.ndarray
        Ground truth cluster labels
    predicted_labels : np.ndarray
        Predicted cluster labels

    Returns
    -------
    metrics : Dict[str, float]
        ARI, NMI, etc.
    """
    from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

    ari = adjusted_rand_score(true_labels, predicted_labels)
    nmi = normalized_mutual_info_score(true_labels, predicted_labels)

    return {
        'ARI': ari,
        'NMI': nmi,
    }


def uncertainty_calibration(
    adata: AnnData,
    ground_truth: Dict[str, np.ndarray],
) -> Dict[str, float]:
    """
    Check if uncertainty estimates are calibrated.

    A well-calibrated model: when it says uncertainty is X,
    the error is actually ~X.

    Parameters
    ----------
    adata : AnnData
        Data with uncertainty estimates
    ground_truth : Dict
        Ground truth

    Returns
    -------
    calibration : Dict[str, float]
        Calibration metrics
    """
    if 'Z_true_mean' not in adata.obsm or 'Z_true_std' not in adata.obsm:
        raise ValueError("Need both mean and std estimates")

    Z_true = ground_truth['Z_true']
    Z_mean = adata.obsm['Z_true_mean']
    Z_std = adata.obsm['Z_true_std']

    # Compute standardized errors
    # If calibrated, (Z_true - Z_mean) / Z_std ~ Normal(0, 1)
    standardized_errors = (Z_true - Z_mean) / (Z_std + 1e-8)

    # Check if mean ≈ 0, std ≈ 1
    mean_error = standardized_errors.mean()
    std_error = standardized_errors.std()

    # Coverage: what fraction of true values fall within predicted intervals?
    # 68% should be within 1 std, 95% within 2 std
    within_1std = np.abs(standardized_errors) < 1
    within_2std = np.abs(standardized_errors) < 2

    coverage_1std = within_1std.mean()
    coverage_2std = within_2std.mean()

    return {
        'mean_standardized_error': mean_error,
        'std_standardized_error': std_error,
        'coverage_1std': coverage_1std,  # Should be ~0.68
        'coverage_2std': coverage_2std,  # Should be ~0.95
        'calibration_score': 1 - np.abs(coverage_1std - 0.68),  # Closer to 1 = better
    }
