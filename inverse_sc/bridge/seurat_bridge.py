"""
Seurat integration via reticulate/anndata conversion.
"""

import numpy as np
import pandas as pd
from anndata import AnnData
from typing import Optional, Dict, Any
import logging

logger = logging.getLogger(__name__)


def from_seurat(
    seurat_obj,
    assay: str = "RNA",
    use_raw: bool = False,
) -> AnnData:
    """
    Convert Seurat object to AnnData.

    Requires rpy2 for R-Python interop.

    Parameters
    ----------
    seurat_obj : rpy2.robjects
        Seurat object from R
    assay : str
        Which assay to extract
    use_raw : bool
        Use raw counts or normalized

    Returns
    -------
    adata : AnnData
        Converted AnnData object

    Examples
    --------
    In R:
    >>> library(Seurat)
    >>> library(reticulate)
    >>> seu <- readRDS("data.rds")

    In Python (via reticulate):
    >>> import inverse_sc as isc
    >>> adata = isc.bridge.from_seurat(seu)
    >>> isc.pp.fit_inverse_model(adata)
    """
    try:
        import rpy2.robjects as ro
        from rpy2.robjects import pandas2ri
        pandas2ri.activate()
    except ImportError:
        raise ImportError("rpy2 required for Seurat integration. Install: pip install rpy2")

    # Extract count matrix
    if use_raw:
        r_code = f"as.matrix(GetAssayData(seurat_obj, assay='{assay}', slot='counts'))"
    else:
        r_code = f"as.matrix(GetAssayData(seurat_obj, assay='{assay}', slot='data'))"

    X = np.array(ro.r(r_code)).T  # Transpose: Seurat is genes x cells

    # Extract metadata
    meta = ro.r("seurat_obj@meta.data")
    obs = pandas2ri.rpy2py(meta)

    # Extract gene names
    gene_names = ro.r(f"rownames(GetAssayData(seurat_obj, assay='{assay}'))")
    var = pd.DataFrame(index=gene_names)

    # Create AnnData
    adata = AnnData(
        X=X,
        obs=obs,
        var=var,
    )

    logger.info(f"Converted Seurat object: {adata.n_obs} cells × {adata.n_vars} genes")

    return adata


def to_seurat(
    adata: AnnData,
    z_key: str = 'Z_true_mean',
    program_key: str = 'program_weights',
) -> Any:
    """
    Convert AnnData back to Seurat.

    After inverse inference, creates a Seurat object with:
    - Raw counts in RNA assay
    - Inferred Z_true in a new assay
    - Program weights in metadata

    Parameters
    ----------
    adata : AnnData
        Data with inverse inference results
    z_key : str
        Key for inferred expression
    program_key : str
        Key for program weights

    Returns
    -------
    seurat_obj : rpy2.robjects
        Seurat object (accessible from R)

    Examples
    --------
    >>> seu = isc.bridge.to_seurat(adata)
    >>> # Now in R:
    >>> # seu has Z_true reduction, program weights in metadata
    """
    try:
        import rpy2.robjects as ro
        from rpy2.robjects import pandas2ri, numpy2ri
        pandas2ri.activate()
        numpy2ri.activate()
    except ImportError:
        raise ImportError("rpy2 required")

    # Transfer data to R
    # Count matrix
    ro.globalenv['counts_matrix'] = adata.X.T  # Transpose to genes x cells

    # Metadata
    ro.globalenv['metadata'] = adata.obs

    # Gene names
    ro.globalenv['gene_names'] = adata.var_names.values

    # Cell names
    ro.globalenv['cell_names'] = adata.obs_names.values

    # Inferred expression
    if z_key in adata.obsm:
        ro.globalenv['z_true'] = adata.obsm[z_key].T

    # Create Seurat object in R
    ro.r('''
        library(Seurat)

        # Set rownames and colnames
        rownames(counts_matrix) <- gene_names
        colnames(counts_matrix) <- cell_names

        # Create Seurat object
        seu <- CreateSeuratObject(
            counts = counts_matrix,
            meta.data = metadata
        )
    ''')

    # Add Z_true as a reduction
    if z_key in adata.obsm:
        ro.r('''
            # Add Z_true as a dimensional reduction
            colnames(z_true) <- cell_names
            rownames(z_true) <- gene_names

            # This goes in a custom assay
            seu[["Z_true"]] <- CreateAssayObject(data = z_true)
        ''')

    # Add program weights to metadata
    if program_key in adata.obsm:
        program_weights = adata.obsm[program_key]
        n_programs = program_weights.shape[1]

        for i in range(n_programs):
            col_name = f'Program_{i}'
            ro.globalenv[col_name] = program_weights[:, i]
            ro.r(f'seu@meta.data${col_name} <- {col_name}')

    seurat_obj = ro.r('seu')

    logger.info("Converted to Seurat object")

    return seurat_obj


def run_inverse_in_r(
    seurat_obj_name: str = "seu",
    **kwargs,
) -> str:
    """
    Generate R code to run inverse inference from Seurat.

    Parameters
    ----------
    seurat_obj_name : str
        Name of Seurat object in R environment
    **kwargs
        Arguments for fit_inverse_model

    Returns
    -------
    r_code : str
        R code to execute

    Examples
    --------
    >>> code = isc.bridge.run_inverse_in_r("pbmc_seurat", n_epochs=100)
    >>> print(code)
    >>> # Copy-paste into R console
    """
    r_code = f"""
# Load required libraries
library(Seurat)
library(reticulate)

# Import inverse_sc
isc <- import("inverse_sc")

# Convert Seurat to AnnData
adata <- isc$bridge$from_seurat({seurat_obj_name})

# Fit inverse model
isc$pp$fit_inverse_model(
    adata,
    n_epochs = {kwargs.get('n_epochs', 100)}L,
    n_latent = {kwargs.get('n_latent', 30)}L,
    n_programs = {kwargs.get('n_programs', 20)}L
)

# Convert back to Seurat
{seurat_obj_name} <- isc$bridge$to_seurat(adata)

# Now {seurat_obj_name} has:
# - Z_true assay with inferred expression
# - Program_0, Program_1, ... in metadata
# - Uncertainty estimates

# You can now do standard Seurat analysis on Z_true:
# DefaultAssay({seurat_obj_name}) <- "Z_true"
# {seurat_obj_name} <- FindNeighbors({seurat_obj_name})
# {seurat_obj_name} <- FindClusters({seurat_obj_name})
# {seurat_obj_name} <- RunUMAP({seurat_obj_name}, dims = 1:30)
"""

    return r_code
