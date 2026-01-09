"""
Basic usage example of InverseSC.

This script demonstrates the simplest way to use InverseSC
for single-cell RNA-seq analysis.
"""

import scanpy as sc
import inverse_sc as isc

# 1. Load your data
adata = sc.read_h5ad("your_data.h5ad")

# 2. Basic QC (optional)
isc.pp.basic_qc(adata, min_genes=200, min_cells=3)

# 3. Fit inverse model
isc.pp.fit_inverse_model(
    adata,
    n_epochs=100,
    n_programs=20,
)

# 4. Downstream analysis on inferred expression
sc.pp.neighbors(adata, use_rep='Z_true_mean')
sc.tl.leiden(adata)
sc.tl.umap(adata)

# 5. Quantify uncertainty
isc.tl.cluster_uncertainty(adata)

# 6. Visualize
sc.pl.umap(adata, color=['leiden', 'cluster_confidence'])

# 7. Differential expression with uncertainty
de_results = isc.tl.differential_expression_robust(
    adata,
    group_key='leiden',
    group1='0',
    group2='1',
)

print(f"Confident DE genes: {de_results['confident'].sum()}")

# 8. Interpret programs
program_info = isc.tl.interpret_programs(adata)
for prog_idx, df in program_info.items():
    print(f"\nProgram {prog_idx} top genes:")
    print(df.head(10))

print("\nDone! Results stored in adata.obsm['Z_true_mean']")
