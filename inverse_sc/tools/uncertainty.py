"""
Uncertainty quantification for downstream analysis.

Key insight: Unlike standard methods, we have uncertainty estimates
for every inferred value. Use them!
"""

import numpy as np
from anndata import AnnData
from typing import Optional, Dict, Tuple
from scipy import stats
import logging

logger = logging.getLogger(__name__)


def cluster_uncertainty(
    adata: AnnData,
    cluster_key: str = 'leiden',
    z_key: str = 'Z_true_mean',
    uncertainty_key: str = 'Z_true_std',
) -> None:
    """
    Quantify confidence in cluster assignments.

    For each cell, computes:
    1. How far is it from cluster centroid (in Z-space)?
    2. How uncertain is its inferred Z?
    3. Confidence score combining both

    Adds to adata.obs:
    - 'cluster_distance': Distance to assigned cluster centroid
    - 'cluster_confidence': Confidence in cluster assignment

    Parameters
    ----------
    adata : AnnData
        Data with cluster labels and inferred Z
    cluster_key : str
        Column in .obs with cluster labels
    z_key : str
        Key in .obsm with inferred expression
    uncertainty_key : str
        Key in .obsm with uncertainty

    Examples
    --------
    >>> isc.tl.cluster_uncertainty(adata)
    >>> # Low-confidence cells
    >>> uncertain = adata.obs['cluster_confidence'] < 0.5
    >>> adata.obs.loc[uncertain, 'cluster_uncertain'] = True
    """
    if cluster_key not in adata.obs.columns:
        raise ValueError(f"Cluster key '{cluster_key}' not found in adata.obs")

    if z_key not in adata.obsm.keys():
        raise ValueError(f"Z key '{z_key}' not found in adata.obsm")

    Z_mean = adata.obsm[z_key]
    Z_std = adata.obsm[uncertainty_key]

    clusters = adata.obs[cluster_key].values
    unique_clusters = np.unique(clusters)

    # Compute cluster centroids
    centroids = {}
    for cluster in unique_clusters:
        mask = clusters == cluster
        centroids[cluster] = Z_mean[mask].mean(axis=0)

    # Distance to assigned cluster centroid
    distances = np.zeros(adata.n_obs)
    for i, cluster in enumerate(clusters):
        centroid = centroids[cluster]
        distances[i] = np.linalg.norm(Z_mean[i] - centroid)

    # Average uncertainty per cell
    avg_uncertainty = Z_std.mean(axis=1)

    # Confidence score
    # High distance OR high uncertainty → low confidence
    # Normalize both to [0, 1] scale

    dist_normalized = distances / (distances.max() + 1e-8)
    uncertainty_normalized = avg_uncertainty / (avg_uncertainty.max() + 1e-8)

    # Confidence = 1 - combined measure
    confidence = 1 - np.sqrt(dist_normalized ** 2 + uncertainty_normalized ** 2) / np.sqrt(2)

    # Store results
    adata.obs['cluster_distance'] = distances
    adata.obs['cluster_uncertainty'] = avg_uncertainty
    adata.obs['cluster_confidence'] = confidence

    logger.info(f"Cluster confidence: mean={confidence.mean():.3f}, min={confidence.min():.3f}")


def differential_expression_robust(
    adata: AnnData,
    group_key: str,
    group1: str,
    group2: str,
    z_key: str = 'Z_true_mean',
    uncertainty_key: str = 'Z_true_std',
    min_fold_change: float = 1.5,
) -> Dict[str, np.ndarray]:
    """
    Differential expression accounting for uncertainty.

    Unlike standard DE (which compares observed counts), this compares
    inferred true expression and accounts for uncertainty.

    A gene is confidently DE if:
    1. Large difference in means
    2. Low uncertainty in both groups
    3. Difference is large relative to uncertainty

    Parameters
    ----------
    adata : AnnData
        Data
    group_key : str
        Column in .obs with group labels
    group1 : str
        First group
    group2 : str
        Second group
    z_key : str
        Key for inferred expression
    uncertainty_key : str
        Key for uncertainty
    min_fold_change : float
        Minimum fold change to consider

    Returns
    -------
    results : Dict[str, np.ndarray]
        DE results with uncertainty-aware statistics

    Examples
    --------
    >>> de = isc.tl.differential_expression_robust(
    ...     adata,
    ...     group_key='cell_type',
    ...     group1='T_cell',
    ...     group2='B_cell',
    ... )
    >>> # High-confidence DE genes
    >>> confident_de = de['gene_names'][de['confident']]
    """
    Z_mean = adata.obsm[z_key]
    Z_std = adata.obsm[uncertainty_key]

    # Get groups
    groups = adata.obs[group_key].values
    mask1 = groups == group1
    mask2 = groups == group2

    if mask1.sum() == 0 or mask2.sum() == 0:
        raise ValueError(f"One or both groups have no cells")

    logger.info(f"DE: {group1} (n={mask1.sum()}) vs {group2} (n={mask2.sum()})")

    # Expression in each group
    Z1_mean = Z_mean[mask1]
    Z1_std = Z_std[mask1]

    Z2_mean = Z_mean[mask2]
    Z2_std = Z_std[mask2]

    # Mean expression per gene per group
    mu1 = Z1_mean.mean(axis=0)
    mu2 = Z2_mean.mean(axis=0)

    # Uncertainty (standard error of the mean)
    # Combines biological variation + measurement uncertainty
    se1 = np.sqrt((Z1_std ** 2).mean(axis=0) + Z1_mean.var(axis=0)) / np.sqrt(mask1.sum())
    se2 = np.sqrt((Z2_std ** 2).mean(axis=0) + Z2_mean.var(axis=0)) / np.sqrt(mask2.sum())

    # Log fold change
    log_fc = np.log2((mu1 + 1) / (mu2 + 1))

    # Fold change
    fold_change = mu1 / (mu2 + 1e-8)

    # T-statistic accounting for uncertainty
    # t = (mu1 - mu2) / sqrt(se1^2 + se2^2)
    t_stat = (mu1 - mu2) / np.sqrt(se1 ** 2 + se2 ** 2 + 1e-8)

    # P-value (two-sided t-test)
    df = mask1.sum() + mask2.sum() - 2
    p_values = 2 * (1 - stats.t.cdf(np.abs(t_stat), df))

    # FDR correction (Benjamini-Hochberg)
    from statsmodels.stats.multitest import multipletests
    _, p_adj, _, _ = multipletests(p_values, method='fdr_bh')

    # Confident DE genes
    # Criteria:
    # 1. Significant after FDR correction
    # 2. Fold change > threshold
    # 3. Difference is large relative to combined uncertainty
    confident = (
        (p_adj < 0.05) &
        (np.abs(fold_change) > min_fold_change) &
        (np.abs(mu1 - mu2) > 2 * np.sqrt(se1 ** 2 + se2 ** 2))
    )

    results = {
        'gene_names': adata.var_names.values,
        'mean_group1': mu1,
        'mean_group2': mu2,
        'log_fold_change': log_fc,
        'fold_change': fold_change,
        'se_group1': se1,
        'se_group2': se2,
        't_statistic': t_stat,
        'p_value': p_values,
        'p_adj': p_adj,
        'confident': confident,
    }

    logger.info(f"Found {confident.sum()} confident DE genes")

    return results


def uncertainty_aware_umap(
    adata: AnnData,
    z_key: str = 'Z_true_mean',
    uncertainty_key: str = 'Z_true_std',
    n_samples: int = 10,
    **umap_kwargs,
) -> None:
    """
    UMAP that visualizes uncertainty.

    Instead of single embedding, samples from posterior and embeds each.
    Result: uncertainty cloud around each cell.

    Parameters
    ----------
    adata : AnnData
        Data
    z_key : str
        Mean expression
    uncertainty_key : str
        Std expression
    n_samples : int
        Number of posterior samples per cell
    **umap_kwargs
        Additional arguments for UMAP

    Adds to adata.obsm:
    - 'X_umap_samples': UMAP coordinates for all samples (n_cells * n_samples, 2)
    """
    try:
        import umap
    except ImportError:
        raise ImportError("Please install umap-learn: pip install umap-learn")

    Z_mean = adata.obsm[z_key]
    Z_std = adata.obsm[uncertainty_key]

    n_cells, n_genes = Z_mean.shape

    # Sample from posterior
    Z_samples = []
    for i in range(n_samples):
        # Sample: Z ~ Normal(mean, std)
        Z_sample = np.random.normal(
            loc=Z_mean,
            scale=Z_std,
        )
        Z_samples.append(Z_sample)

    Z_samples = np.concatenate(Z_samples, axis=0)  # (n_cells * n_samples, n_genes)

    # UMAP embedding
    logger.info(f"Computing UMAP on {n_samples} samples per cell...")
    reducer = umap.UMAP(**umap_kwargs)
    embedding = reducer.fit_transform(Z_samples)

    # Store
    adata.obsm['X_umap_samples'] = embedding

    # Also compute mean UMAP (standard approach)
    mean_embedding = reducer.fit_transform(Z_mean)
    adata.obsm['X_umap_mean'] = mean_embedding

    logger.info("Done. Use adata.obsm['X_umap_samples'] to visualize uncertainty.")


def estimate_identifiable_genes(
    adata: AnnData,
    z_key: str = 'Z_true_mean',
    uncertainty_key: str = 'Z_true_std',
) -> np.ndarray:
    """
    Identify which genes are well-identified from data.

    Genes with high uncertainty are poorly identified
    (e.g., due to extreme dropout).

    Parameters
    ----------
    adata : AnnData
        Data
    z_key : str
        Mean expression
    uncertainty_key : str
        Uncertainty

    Returns
    -------
    identifiable_mask : np.ndarray
        Boolean mask of identifiable genes
    """
    Z_mean = adata.obsm[z_key]
    Z_std = adata.obsm[uncertainty_key]

    # Coefficient of variation per gene
    # CV = std / mean
    # High CV → poorly identified

    mean_per_gene = Z_mean.mean(axis=0)
    std_per_gene = Z_std.mean(axis=0)

    cv = std_per_gene / (mean_per_gene + 1e-8)

    # Threshold: CV < 1 is reasonably well-identified
    identifiable = cv < 1.0

    adata.var['coefficient_variation'] = cv
    adata.var['identifiable'] = identifiable

    logger.info(f"{identifiable.sum()} / {len(identifiable)} genes are well-identified")

    return identifiable
