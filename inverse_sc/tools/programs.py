"""
Transcriptional program interpretation and analysis.
"""

import numpy as np
from anndata import AnnData
from typing import Optional, Dict, List
import pandas as pd
import logging

logger = logging.getLogger(__name__)


def interpret_programs(
    adata: AnnData,
    program_key: str = 'program_weights',
    top_genes: int = 50,
) -> Dict[int, pd.DataFrame]:
    """
    Interpret learned transcriptional programs.

    For each program, identify:
    1. Top genes (highest loading)
    2. Cells with high program activity
    3. Program enrichment patterns

    Parameters
    ----------
    adata : AnnData
        Data with fitted model
    program_key : str
        Key in .obsm with program weights
    top_genes : int
        Number of top genes per program

    Returns
    -------
    program_info : Dict[int, pd.DataFrame]
        Information about each program

    Examples
    --------
    >>> program_info = isc.tl.interpret_programs(adata)
    >>> # See top genes for program 0
    >>> print(program_info[0].head(10))
    """
    if 'inverse_model_obj' not in adata.uns:
        raise ValueError("Model not found. Run fit_inverse_model first.")

    model = adata.uns['inverse_model_obj']['model']
    program_weights = adata.obsm[program_key]

    # Extract program signatures from model
    with torch.no_grad():
        program_signatures = model.program_signatures.cpu().numpy()

    n_programs, n_genes = program_signatures.shape

    program_info = {}

    for prog_idx in range(n_programs):
        signature = program_signatures[prog_idx]

        # Top genes (highest values in signature)
        top_indices = np.argsort(signature)[-top_genes:][::-1]

        # Program activity per cell
        activity = program_weights[:, prog_idx]

        # Cells with high activity
        high_activity_cells = np.argsort(activity)[-100:][::-1]

        # Create DataFrame
        df = pd.DataFrame({
            'gene': adata.var_names[top_indices],
            'loading': signature[top_indices],
        })

        program_info[prog_idx] = df

        logger.info(f"Program {prog_idx}: top gene = {df.iloc[0]['gene']}")

    return program_info


def program_enrichment(
    adata: AnnData,
    program_idx: int,
    gene_sets: Dict[str, List[str]],
    program_key: str = 'program_weights',
) -> pd.DataFrame:
    """
    Test for enrichment of gene sets in a program.

    Parameters
    ----------
    adata : AnnData
        Data
    program_idx : int
        Which program to analyze
    gene_sets : Dict[str, List[str]]
        Gene sets to test (e.g., from MSigDB)
    program_key : str
        Key for program weights

    Returns
    -------
    enrichment_results : pd.DataFrame
        Enrichment p-values for each gene set

    Examples
    --------
    >>> gene_sets = {
    ...     'T_cell_signature': ['CD3D', 'CD3E', 'CD8A'],
    ...     'B_cell_signature': ['CD79A', 'CD79B', 'MS4A1'],
    ... }
    >>> enrichment = isc.tl.program_enrichment(adata, program_idx=0, gene_sets=gene_sets)
    """
    if 'inverse_model_obj' not in adata.uns:
        raise ValueError("Model not found.")

    model = adata.uns['inverse_model_obj']['model']

    with torch.no_grad():
        program_signatures = model.program_signatures.cpu().numpy()

    signature = program_signatures[program_idx]

    results = []

    for gene_set_name, genes in gene_sets.items():
        # Find genes in data
        genes_in_data = [g for g in genes if g in adata.var_names]

        if len(genes_in_data) == 0:
            continue

        # Indices
        gene_indices = [np.where(adata.var_names == g)[0][0] for g in genes_in_data]

        # Mean loading in gene set vs background
        in_set = signature[gene_indices].mean()
        background = signature.mean()

        # Simple t-test
        from scipy.stats import ttest_ind

        t_stat, p_val = ttest_ind(
            signature[gene_indices],
            signature,
        )

        results.append({
            'gene_set': gene_set_name,
            'n_genes': len(genes_in_data),
            'mean_loading': in_set,
            'background_loading': background,
            'fold_enrichment': in_set / (background + 1e-8),
            'p_value': p_val,
        })

    results_df = pd.DataFrame(results)

    # FDR correction
    if len(results_df) > 0:
        from statsmodels.stats.multitest import multipletests
        _, p_adj, _, _ = multipletests(results_df['p_value'], method='fdr_bh')
        results_df['p_adj'] = p_adj

    return results_df


def assign_program_names(
    adata: AnnData,
    program_info: Dict[int, pd.DataFrame],
    auto_name: bool = True,
) -> Dict[int, str]:
    """
    Assign interpretable names to programs.

    Parameters
    ----------
    adata : AnnData
        Data
    program_info : Dict[int, pd.DataFrame]
        Program interpretation from interpret_programs
    auto_name : bool
        Automatically name based on top genes

    Returns
    -------
    program_names : Dict[int, str]
        Names for each program
    """
    program_names = {}

    for prog_idx, df in program_info.items():
        if auto_name:
            # Name after top 3 genes
            top_genes = df['gene'].values[:3]
            name = f"Program_{prog_idx}_{'/'.join(top_genes)}"
        else:
            name = f"Program_{prog_idx}"

        program_names[prog_idx] = name

    return program_names


def plot_program_activity(
    adata: AnnData,
    program_idx: int,
    program_key: str = 'program_weights',
    embedding: str = 'X_umap',
):
    """
    Plot program activity on UMAP.

    Parameters
    ----------
    adata : AnnData
        Data
    program_idx : int
        Program to plot
    program_key : str
        Key for program weights
    embedding : str
        Embedding to use for plotting
    """
    import matplotlib.pyplot as plt

    if embedding not in adata.obsm:
        raise ValueError(f"Embedding {embedding} not found. Run UMAP first.")

    coords = adata.obsm[embedding]
    activity = adata.obsm[program_key][:, program_idx]

    plt.figure(figsize=(8, 6))
    scatter = plt.scatter(
        coords[:, 0],
        coords[:, 1],
        c=activity,
        cmap='viridis',
        s=5,
        alpha=0.7,
    )
    plt.colorbar(scatter, label='Program Activity')
    plt.title(f'Program {program_idx} Activity')
    plt.xlabel('UMAP 1')
    plt.ylabel('UMAP 2')
    plt.tight_layout()
    plt.show()


# Add torch import at the top
import torch
