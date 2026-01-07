"""
Quality control functions compatible with Scanpy.
"""

import numpy as np
from anndata import AnnData
from typing import Optional, Tuple
import logging

logger = logging.getLogger(__name__)


def basic_qc(
    adata: AnnData,
    min_genes: int = 200,
    min_cells: int = 3,
    max_pct_mito: float = 20.0,
    mito_prefix: str = "MT-",
    copy: bool = False,
) -> Optional[AnnData]:
    """
    Basic quality control filtering.

    Parameters
    ----------
    adata : AnnData
        Data
    min_genes : int
        Minimum genes per cell
    min_cells : int
        Minimum cells per gene
    max_pct_mito : float
        Maximum mitochondrial percentage
    mito_prefix : str
        Prefix for mitochondrial genes
    copy : bool
        Return copy

    Returns
    -------
    adata : Optional[AnnData]
        Filtered data
    """
    if copy:
        adata = adata.copy()

    logger.info(f"Starting QC: {adata.n_obs} cells × {adata.n_vars} genes")

    # Compute QC metrics
    X = adata.X
    if hasattr(X, 'toarray'):
        X = X.toarray()

    # Genes per cell
    adata.obs['n_genes'] = (X > 0).sum(axis=1)

    # Cells per gene
    adata.var['n_cells'] = (X > 0).sum(axis=0)

    # Total counts
    adata.obs['total_counts'] = X.sum(axis=1)
    adata.var['total_counts'] = X.sum(axis=0)

    # Mitochondrial percentage
    mito_genes = adata.var_names.str.startswith(mito_prefix)
    adata.obs['pct_mito'] = (
        X[:, mito_genes].sum(axis=1) / (X.sum(axis=1) + 1e-8) * 100
    )

    # Filter cells
    cell_filter = (
        (adata.obs['n_genes'] >= min_genes) &
        (adata.obs['pct_mito'] <= max_pct_mito)
    )

    n_cells_removed = (~cell_filter).sum()
    logger.info(f"Removing {n_cells_removed} cells")

    adata = adata[cell_filter, :].copy()

    # Filter genes
    gene_filter = adata.var['n_cells'] >= min_cells
    n_genes_removed = (~gene_filter).sum()
    logger.info(f"Removing {n_genes_removed} genes")

    adata = adata[:, gene_filter].copy()

    logger.info(f"After QC: {adata.n_obs} cells × {adata.n_vars} genes")

    if copy:
        return adata


def filter_cells(
    adata: AnnData,
    min_genes: Optional[int] = None,
    max_genes: Optional[int] = None,
    min_counts: Optional[int] = None,
    max_counts: Optional[int] = None,
    copy: bool = False,
) -> Optional[AnnData]:
    """
    Filter cells by metrics.

    Parameters
    ----------
    adata : AnnData
        Data
    min_genes : Optional[int]
        Minimum genes detected
    max_genes : Optional[int]
        Maximum genes detected
    min_counts : Optional[int]
        Minimum total counts
    max_counts : Optional[int]
        Maximum total counts
    copy : bool
        Return copy

    Returns
    -------
    adata : Optional[AnnData]
        Filtered data
    """
    if copy:
        adata = adata.copy()

    # Compute metrics if not present
    if 'n_genes' not in adata.obs.columns:
        X = adata.X
        if hasattr(X, 'toarray'):
            X = X.toarray()
        adata.obs['n_genes'] = (X > 0).sum(axis=1)
        adata.obs['total_counts'] = X.sum(axis=1)

    # Build filter
    cell_filter = np.ones(adata.n_obs, dtype=bool)

    if min_genes is not None:
        cell_filter &= adata.obs['n_genes'] >= min_genes

    if max_genes is not None:
        cell_filter &= adata.obs['n_genes'] <= max_genes

    if min_counts is not None:
        cell_filter &= adata.obs['total_counts'] >= min_counts

    if max_counts is not None:
        cell_filter &= adata.obs['total_counts'] <= max_counts

    n_removed = (~cell_filter).sum()
    logger.info(f"Filtering {n_removed} cells")

    adata = adata[cell_filter, :].copy()

    if copy:
        return adata


def filter_genes(
    adata: AnnData,
    min_cells: Optional[int] = None,
    min_counts: Optional[int] = None,
    copy: bool = False,
) -> Optional[AnnData]:
    """
    Filter genes by detection.

    Parameters
    ----------
    adata : AnnData
        Data
    min_cells : Optional[int]
        Minimum cells expressing gene
    min_counts : Optional[int]
        Minimum total counts
    copy : bool
        Return copy

    Returns
    -------
    adata : Optional[AnnData]
        Filtered data
    """
    if copy:
        adata = adata.copy()

    # Compute metrics if not present
    if 'n_cells' not in adata.var.columns:
        X = adata.X
        if hasattr(X, 'toarray'):
            X = X.toarray()
        adata.var['n_cells'] = (X > 0).sum(axis=0)
        adata.var['total_counts'] = X.sum(axis=0)

    # Build filter
    gene_filter = np.ones(adata.n_vars, dtype=bool)

    if min_cells is not None:
        gene_filter &= adata.var['n_cells'] >= min_cells

    if min_counts is not None:
        gene_filter &= adata.var['total_counts'] >= min_counts

    n_removed = (~gene_filter).sum()
    logger.info(f"Filtering {n_removed} genes")

    adata = adata[:, gene_filter].copy()

    if copy:
        return adata
