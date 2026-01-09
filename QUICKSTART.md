# InverseSC Quick Start Guide

## Installation

```bash
# Create environment
conda env create -f environment.yml
conda activate inverse-sc

# Install package
pip install -e .
```

## 5-Minute Example

```python
import scanpy as sc
import inverse_sc as isc

# Load data
adata = sc.datasets.pbmc3k()

# Basic QC
isc.pp.basic_qc(adata)

# Fit inverse model (this is the key step!)
isc.pp.fit_inverse_model(
    adata,
    n_epochs=100,
    n_programs=15,
)

# Downstream analysis on inferred expression
sc.pp.neighbors(adata, use_rep='Z_true_mean')
sc.tl.leiden(adata)
sc.tl.umap(adata)

# Quantify uncertainty
isc.tl.cluster_uncertainty(adata)

# Visualize
sc.pl.umap(adata, color=['leiden', 'cluster_confidence'])
```

## What Just Happened?

Instead of analyzing log-normalized counts (standard approach), you just:

1. **Modeled the measurement process**: Capture → Amplification → Sequencing
2. **Inverted it**: Inferred what biology must have existed before measurement
3. **Got uncertainty estimates**: Every value has a confidence score

## Key Outputs

After `fit_inverse_model`, your AnnData has:

- `adata.obsm['Z_true_mean']`: Inferred true expression (use this for downstream analysis)
- `adata.obsm['Z_true_std']`: Uncertainty estimates
- `adata.obsm['program_weights']`: Cell state composition
- `adata.var['capture_prob_estimate']`: Estimated capture efficiency per gene

## Next Steps

### Explore Uncertainty
```python
# Which cells are uncertain?
uncertain = adata.obs['cluster_confidence'] < 0.5
print(f"{uncertain.sum()} uncertain cells")

# Visualize uncertainty
sc.pl.umap(adata, color='cluster_confidence', cmap='RdYlGn')
```

### Robust Differential Expression
```python
# DE that accounts for uncertainty
de = isc.tl.differential_expression_robust(
    adata,
    group_key='leiden',
    group1='0',
    group2='1',
)

# Only confidently differential genes
confident_de = de['gene_names'][de['confident']]
print(f"Confident DE genes: {len(confident_de)}")
```

### Interpret Programs
```python
# What are the learned transcriptional programs?
programs = isc.tl.interpret_programs(adata)

# Top genes for each program
for prog_idx, df in programs.items():
    print(f"\nProgram {prog_idx}:")
    print(df.head(10))
```

## Validation

Want to test on data with known ground truth?

```python
# Generate synthetic data
adata, truth = isc.validation.generate_synthetic_data(seed=42)

# Fit model
isc.pp.fit_inverse_model(adata)

# Benchmark
results = isc.validation.benchmark_against_scanpy(adata, truth)
print(results)

# Check uncertainty calibration
calib = isc.validation.uncertainty_calibration(adata, truth)
print(f"Calibration score: {calib['calibration_score']:.3f}")
```

## Using with Seurat (R)

```r
library(Seurat)
library(reticulate)

# Import InverseSC
isc <- import("inverse_sc")

# Convert Seurat → AnnData
adata <- isc$bridge$from_seurat(seurat_obj)

# Fit inverse model
isc$pp$fit_inverse_model(adata, n_epochs=100L)

# Convert back to Seurat
seurat_obj <- isc$bridge$to_seurat(adata)

# Now use Z_true assay for downstream analysis
DefaultAssay(seurat_obj) <- "Z_true"
```

## Common Questions

### Q: How is this different from scVI?

**scVI**: Learns arbitrary reconstruction (neural network decoder)
**InverseSC**: Uses physical measurement model (fixed decoder based on real biology)

Result: InverseSC's latent space has physical meaning (true expression), not just good reconstruction.

### Q: When should I use this vs. standard normalization?

Use InverseSC when:
- High dropout (>70%)
- Need uncertainty estimates
- Comparing very different sequencing depths
- Want interpretable programs

Stick with standard when:
- Quick exploration
- Very deep data (>50k reads/cell)
- Computational constraints

### Q: How long does it take?

- 1k cells, 2k genes: ~2 minutes (CPU), ~30 seconds (GPU)
- 10k cells, 2k genes: ~10 minutes (CPU), ~2 minutes (GPU)
- Scales roughly linearly

### Q: Can I use this with spatial data?

Yes! Just treat spatial coordinates as additional cell features. Future versions will have explicit spatial priors.

## Troubleshooting

### Model not converging?

- Increase `n_epochs` (try 200-300)
- Decrease `learning_rate` (try 0.0001)
- Check for extreme outliers in data

### Out of memory?

- Decrease `batch_size` (try 128 or 64)
- Use fewer genes (filter to highly variable)
- Use GPU if available (`use_cuda=True`)

### Results don't make sense?

- Check QC (did you filter low-quality cells?)
- Verify data is raw counts (not already normalized)
- Try different `n_programs` (10-30 usually works)

## Learn More

- **Notebooks**: See `notebooks/` for detailed examples
- **Scripts**: See `scripts/` for run time examples on pbmc3k
- **API Reference**: See `docs/API_REFERENCE.md`
- **Methodology**: See `docs/METHODOLOGY.md` for mathematical details
- **GitHub**: [inverse-problem-scrna](https://github.com/pritampanda15/Inverse-Problem-Framing-scRNAseq)

## Citation

```bibtex
@article{inversesc2026,
  title={Inverse Problem Formulation for Single-Cell RNA Sequencing},
  author={Pritam Kumar Panda},
  journal={bioRxiv},
  year={2026}
}
```

---

**Happy analyzing! Questions? Open an issue on GitHub.**
