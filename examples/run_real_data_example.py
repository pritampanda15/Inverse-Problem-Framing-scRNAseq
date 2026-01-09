#!/usr/bin/env python
"""
Real Data Example: PBMC Dataset

Applying InverseSC to real single-cell data and comparing with standard approaches.

Usage:
    python scripts/run_real_data_example.py
"""

import numpy as np
import scanpy as sc
import inverse_sc as isc
import matplotlib.pyplot as plt
import argparse
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def main(args):
    # =========================================================================
    # 1. Load Data
    # =========================================================================
    logger.info("Loading PBMC3k dataset...")
    adata = sc.datasets.pbmc3k()
    logger.info(f"Raw data: {adata.shape}")

    # =========================================================================
    # 2. Basic Quality Control
    # =========================================================================
    logger.info("Running basic QC...")
    isc.pp.basic_qc(adata, min_genes=200, max_pct_mito=5)
    logger.info(f"After QC: {adata.shape}")

    # =========================================================================
    # 3. Standard Pipeline (Baseline)
    # =========================================================================
    logger.info("Running standard Scanpy pipeline...")

    # Save raw counts
    adata.layers['counts'] = adata.X.copy()

    # Standard Scanpy preprocessing
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)

    # PCA and clustering
    sc.tl.pca(adata, mask_var='highly_variable')
    sc.pp.neighbors(adata)
    sc.tl.leiden(adata, key_added='leiden_standard', resolution=args.resolution)
    sc.tl.umap(adata)

    # Save standard UMAP
    adata.obsm['X_umap_standard'] = adata.obsm['X_umap'].copy()
    logger.info(f"Standard pipeline: {adata.obs['leiden_standard'].nunique()} clusters")

    # =========================================================================
    # 4. Inverse Problem Approach
    # =========================================================================
    logger.info("Fitting inverse model on raw counts...")
    isc.pp.fit_inverse_model(
        adata,
        layer='counts',
        n_latent=args.n_latent,
        n_programs=args.n_programs,
        n_epochs=args.n_epochs,
        batch_size=args.batch_size,
        n_posterior_samples=args.n_samples,
    )

    # Clustering on inferred expression
    logger.info("Clustering on inferred expression...")
    # Run PCA on inferred expression for efficient neighbor computation
    from sklearn.decomposition import PCA
    pca = PCA(n_components=50)
    adata.obsm['Z_true_pca'] = pca.fit_transform(np.log1p(adata.obsm['Z_true_mean']))
    logger.info(f"PCA on Z_true_mean: {adata.obsm['Z_true_pca'].shape}")
    sc.pp.neighbors(adata, use_rep='Z_true_pca')
    sc.tl.leiden(adata, key_added='leiden_inverse', resolution=args.resolution)
    sc.tl.umap(adata)

    # Save inverse UMAP
    adata.obsm['X_umap_inverse'] = adata.obsm['X_umap'].copy()
    logger.info(f"Inverse pipeline: {adata.obs['leiden_inverse'].nunique()} clusters")

    # =========================================================================
    # 5. Compare Results
    # =========================================================================
    if not args.no_plots:
        logger.info("Generating comparison plots...")

        fig, axes = plt.subplots(1, 2, figsize=(14, 5))

        # Standard
        adata.obsm['X_umap'] = adata.obsm['X_umap_standard']
        sc.pl.umap(adata, color='leiden_standard', ax=axes[0], show=False, title='Standard Pipeline')

        # Inverse
        adata.obsm['X_umap'] = adata.obsm['X_umap_inverse']
        sc.pl.umap(adata, color='leiden_inverse', ax=axes[1], show=False, title='Inverse Problem')

        plt.tight_layout()
        plt.savefig(args.output_dir + '/comparison_umap.png', dpi=150, bbox_inches='tight')
        logger.info(f"Saved comparison plot to {args.output_dir}/comparison_umap.png")

    # =========================================================================
    # 6. Uncertainty Analysis
    # =========================================================================
    logger.info("Computing cluster uncertainty...")
    isc.tl.cluster_uncertainty(adata, cluster_key='leiden_inverse')

    # Use relative threshold (bottom 25% are "uncertain")
    confidence_threshold = np.percentile(adata.obs['cluster_confidence'], 25)
    uncertain_cells = adata.obs['cluster_confidence'] < confidence_threshold
    logger.info(f"Uncertain cells (bottom 25%): {uncertain_cells.sum()} / {adata.n_obs}")

    # Which clusters have most uncertainty?
    uncertainty_by_cluster = adata.obs.groupby('leiden_inverse')['cluster_confidence'].mean().sort_values()
    logger.info("Mean confidence by cluster:")
    print(uncertainty_by_cluster)

    if not args.no_plots:
        fig, axes = plt.subplots(1, 3, figsize=(18, 5))

        sc.pl.umap(adata, color='leiden_inverse', ax=axes[0], show=False)
        sc.pl.umap(adata, color='cluster_confidence', cmap='RdYlGn', ax=axes[1], show=False)
        sc.pl.umap(adata, color='cluster_distance', cmap='viridis', ax=axes[2], show=False)

        plt.tight_layout()
        plt.savefig(args.output_dir + '/uncertainty_analysis.png', dpi=150, bbox_inches='tight')
        logger.info(f"Saved uncertainty plot to {args.output_dir}/uncertainty_analysis.png")

    # =========================================================================
    # 7. Program Interpretation
    # =========================================================================
    logger.info("Interpreting learned programs...")
    program_info = isc.tl.interpret_programs(adata, top_genes=50)

    # Filter out ribosomal and mitochondrial genes for cleaner interpretation
    def filter_housekeeping(df):
        """Remove ribosomal (RPS/RPL), mitochondrial (MT-), and other housekeeping genes."""
        mask = ~df['gene'].str.match(r'^(RPS|RPL|MT-|MALAT1|ACTB|B2M|TMSB|EEF|FTL|FTH)')
        return df[mask].head(15)

    # Show top informative genes for first 5 programs
    print("\n" + "=" * 50)
    print("TOP GENES PER PROGRAM (filtered)")
    print("=" * 50)
    for prog_idx in range(min(5, len(program_info))):
        filtered = filter_housekeeping(program_info[prog_idx])
        if len(filtered) > 0:
            top_gene = filtered.iloc[0]['gene']
            print(f"\n=== Program {prog_idx} (top: {top_gene}) ===")
            print(filtered.head(10).to_string(index=False))
        else:
            print(f"\n=== Program {prog_idx} (mostly housekeeping genes) ===")

    if not args.no_plots:
        n_programs_to_plot = min(6, adata.obsm['program_weights'].shape[1])
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
        axes = axes.flatten()

        for prog_idx in range(n_programs_to_plot):
            # Add program weights as obs column for plotting
            adata.obs[f'program_{prog_idx}'] = adata.obsm['program_weights'][:, prog_idx]
            sc.pl.umap(
                adata,
                color=f'program_{prog_idx}',
                cmap='viridis',
                ax=axes[prog_idx],
                show=False,
                title=f'Program {prog_idx}',
                colorbar_loc='right',
            )

        plt.tight_layout()
        plt.savefig(args.output_dir + '/program_activity.png', dpi=150, bbox_inches='tight')
        logger.info(f"Saved program activity plot to {args.output_dir}/program_activity.png")

    # =========================================================================
    # 8. Cell Type Annotation & Comparison
    # =========================================================================
    logger.info("Annotating cell types using known markers...")

    # Known PBMC marker genes
    marker_genes = {
        'CD4+ T': ['IL7R', 'CD4', 'TCF7', 'LEF1'],
        'CD8+ T': ['CD8A', 'CD8B', 'GZMK'],
        'Naive T': ['CCR7', 'SELL', 'TCF7'],
        'Memory T': ['S100A4', 'IL7R'],
        'NK': ['NKG7', 'GNLY', 'KLRD1', 'KLRF1'],
        'B': ['CD79A', 'MS4A1', 'CD79B'],
        'CD14+ Mono': ['CD14', 'LYZ', 'S100A8', 'S100A9'],
        'FCGR3A+ Mono': ['FCGR3A', 'MS4A7'],
        'DC': ['FCER1A', 'CST3', 'CLEC10A'],
        'Platelet': ['PPBP', 'PF4'],
    }

    # Filter to available markers
    available_markers = {}
    for cell_type, genes in marker_genes.items():
        available = [g for g in genes if g in adata.var_names]
        if available:
            available_markers[cell_type] = available

    logger.info(f"Available marker sets: {list(available_markers.keys())}")

    def annotate_clusters(adata, cluster_key, markers_dict):
        """Annotate clusters based on marker gene expression."""
        annotations = {}
        cluster_scores = {}

        for cluster in adata.obs[cluster_key].unique():
            mask = (adata.obs[cluster_key] == cluster).values  # Convert to numpy array
            scores = {}

            for cell_type, genes in markers_dict.items():
                # Get mean expression of marker genes in this cluster
                gene_mask = adata.var_names.isin(genes)
                if gene_mask.sum() > 0:
                    # Use the log-normalized expression in adata.X
                    if hasattr(adata.X, 'toarray'):
                        expr = adata.X[mask][:, gene_mask].toarray()
                    else:
                        expr = np.asarray(adata.X[mask][:, gene_mask])
                    scores[cell_type] = np.mean(expr)

            # Assign to highest scoring cell type
            if scores:
                best_type = max(scores, key=scores.get)
                best_score = scores[best_type]
                annotations[cluster] = best_type if best_score > 0.5 else f'Unknown ({best_type}?)'
                cluster_scores[cluster] = scores
            else:
                annotations[cluster] = 'Unknown'
                cluster_scores[cluster] = {}

        return annotations, cluster_scores

    # Annotate both clustering results
    standard_annotations, standard_scores = annotate_clusters(adata, 'leiden_standard', available_markers)
    inverse_annotations, inverse_scores = annotate_clusters(adata, 'leiden_inverse', available_markers)

    # Add annotations to adata
    adata.obs['celltype_standard'] = adata.obs['leiden_standard'].map(standard_annotations)
    adata.obs['celltype_inverse'] = adata.obs['leiden_inverse'].map(inverse_annotations)

    # Print annotation summary
    print("\n" + "=" * 70)
    print("CELL TYPE ANNOTATION COMPARISON")
    print("=" * 70)

    print("\n--- Standard Pipeline Clusters ---")
    for cluster, celltype in sorted(standard_annotations.items(), key=lambda x: int(x[0])):
        n_cells = (adata.obs['leiden_standard'] == cluster).sum()
        print(f"  Cluster {cluster}: {celltype} ({n_cells} cells)")

    print("\n--- Inverse Problem Clusters ---")
    for cluster, celltype in sorted(inverse_annotations.items(), key=lambda x: int(x[0])):
        n_cells = (adata.obs['leiden_inverse'] == cluster).sum()
        print(f"  Cluster {cluster}: {celltype} ({n_cells} cells)")

    # Compare cell type distributions
    print("\n--- Cell Type Summary ---")
    standard_types = adata.obs['celltype_standard'].value_counts()
    inverse_types = adata.obs['celltype_inverse'].value_counts()

    all_types = set(standard_types.index) | set(inverse_types.index)

    print(f"\n{'Cell Type':<25} {'Standard':<15} {'Inverse':<15} {'Difference':<15}")
    print("-" * 70)
    for ct in sorted(all_types):
        std_n = standard_types.get(ct, 0)
        inv_n = inverse_types.get(ct, 0)
        diff = inv_n - std_n
        diff_str = f"+{diff}" if diff > 0 else str(diff)
        print(f"{ct:<25} {std_n:<15} {inv_n:<15} {diff_str:<15}")

    # Identify cell types found ONLY by inverse approach
    standard_unique_types = set(t for t in standard_types.index if 'Unknown' not in t)
    inverse_unique_types = set(t for t in inverse_types.index if 'Unknown' not in t)

    only_in_inverse = inverse_unique_types - standard_unique_types
    only_in_standard = standard_unique_types - inverse_unique_types

    print("\n--- Key Findings ---")
    if only_in_inverse:
        print(f"Cell types found ONLY by Inverse approach: {only_in_inverse}")
    if only_in_standard:
        print(f"Cell types found ONLY by Standard approach: {only_in_standard}")

    # Count how many more granular clusters per cell type
    print("\n--- Cluster Granularity by Cell Type ---")
    for ct in sorted(standard_unique_types & inverse_unique_types):
        std_clusters = sum(1 for v in standard_annotations.values() if v == ct)
        inv_clusters = sum(1 for v in inverse_annotations.values() if v == ct)
        if inv_clusters > std_clusters:
            print(f"  {ct}: {std_clusters} cluster(s) in Standard → {inv_clusters} cluster(s) in Inverse (+{inv_clusters-std_clusters} subtypes)")

    # Generate annotated comparison plot
    if not args.no_plots:
        logger.info("Generating annotated comparison plots...")

        fig, axes = plt.subplots(2, 2, figsize=(16, 14))

        # Standard - clusters
        adata.obsm['X_umap'] = adata.obsm['X_umap_standard']
        sc.pl.umap(adata, color='leiden_standard', ax=axes[0, 0], show=False,
                   title='Standard Pipeline (Clusters)', legend_loc='on data')

        # Standard - cell types
        sc.pl.umap(adata, color='celltype_standard', ax=axes[0, 1], show=False,
                   title='Standard Pipeline (Cell Types)')

        # Inverse - clusters
        adata.obsm['X_umap'] = adata.obsm['X_umap_inverse']
        sc.pl.umap(adata, color='leiden_inverse', ax=axes[1, 0], show=False,
                   title='Inverse Problem (Clusters)', legend_loc='on data')

        # Inverse - cell types
        sc.pl.umap(adata, color='celltype_inverse', ax=axes[1, 1], show=False,
                   title='Inverse Problem (Cell Types)')

        plt.tight_layout()
        plt.savefig(args.output_dir + '/annotated_comparison.png', dpi=150, bbox_inches='tight')
        logger.info(f"Saved annotated comparison to {args.output_dir}/annotated_comparison.png")

        # Marker gene dotplot for both
        fig, axes = plt.subplots(1, 2, figsize=(20, 8))

        # Flatten marker genes for dotplot
        flat_markers = []
        for genes in available_markers.values():
            flat_markers.extend(genes)
        flat_markers = list(dict.fromkeys(flat_markers))  # Remove duplicates, preserve order

        sc.pl.dotplot(adata, flat_markers, groupby='celltype_standard', ax=axes[0],
                      show=False, title='Standard Pipeline - Marker Expression')
        sc.pl.dotplot(adata, flat_markers, groupby='celltype_inverse', ax=axes[1],
                      show=False, title='Inverse Problem - Marker Expression')

        plt.tight_layout()
        plt.savefig(args.output_dir + '/marker_dotplot_comparison.png', dpi=150, bbox_inches='tight')
        logger.info(f"Saved marker dotplot to {args.output_dir}/marker_dotplot_comparison.png")

    # =========================================================================
    # 9. Summary
    # =========================================================================
    print("\n" + "=" * 50)
    print("ANALYSIS SUMMARY")
    print("=" * 50)
    print(f"Cells: {adata.n_obs}")
    print(f"Genes: {adata.n_vars}")
    print(f"Clusters (standard): {adata.obs['leiden_standard'].nunique()}")
    print(f"Clusters (inverse): {adata.obs['leiden_inverse'].nunique()}")
    print(f"Mean cluster confidence: {adata.obs['cluster_confidence'].mean():.3f}")
    print(f"Confidence range: {adata.obs['cluster_confidence'].min():.3f} - {adata.obs['cluster_confidence'].max():.3f}")
    print(f"Low confidence cells (bottom 25%): {uncertain_cells.sum()} (threshold: {confidence_threshold:.3f})")
    print("=" * 50)

    # Save results
    if args.save_adata:
        output_path = args.output_dir + '/pbmc3k_inverse.h5ad'
        adata.write(output_path)
        logger.info(f"Saved results to {output_path}")

    logger.info("Analysis complete!")
    return adata


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run inverse problem analysis on PBMC3k data')

    # Model parameters
    parser.add_argument('--n_latent', type=int, default=50, help='Latent dimension')
    parser.add_argument('--n_programs', type=int, default=15, help='Number of programs')
    parser.add_argument('--n_epochs', type=int, default=200, help='Training epochs')
    parser.add_argument('--batch_size', type=int, default=256, help='Batch size')
    parser.add_argument('--n_samples', type=int, default=50, help='Posterior samples (fewer=faster)')
    parser.add_argument('--resolution', type=float, default=1.0, help='Leiden clustering resolution (higher=more clusters)')

    # Output options
    parser.add_argument('--output_dir', type=str, default='results', help='Output directory')
    parser.add_argument('--no_plots', action='store_true', help='Skip generating plots')
    parser.add_argument('--save_adata', action='store_true', help='Save annotated data to h5ad')

    args = parser.parse_args()

    # Create output directory
    import os
    os.makedirs(args.output_dir, exist_ok=True)

    main(args)
