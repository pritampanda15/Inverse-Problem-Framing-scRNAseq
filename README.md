# InverseSC: Inverse Problem Framework for Single-Cell RNA-seq

A novel approach to single-cell analysis that inverts the measurement process rather than analyzing distorted observations directly.

## The Core Idea

Traditional pipelines (Scanpy/Seurat) analyze the observed count matrix **X** directly:
```
Raw Counts → Normalize → Reduce → Cluster
```

InverseSC asks: "Given observed **X**, what was the true biological state **Z** before measurement destroyed it?"
```
X = M(Z) + noise
Goal: Infer posterior p(Z | X)
```

## Key Differences from Standard Approaches

| Aspect | Scanpy/Seurat | InverseSC |
|--------|---------------|-----------|
| What you model | Observed matrix X | Pre-measurement state Z |
| Zeros | "Not expressed" | "Uncertain due to dropout" |
| Normalization | Heuristic scaling | Physical measurement model |
| Output | Point estimates | Posterior distributions with uncertainty |
| Batch effects | Empirical correction | Part of measurement operator |

## Installation

```bash
pip install inverse-sc
```

Or from source:
```bash
git clone https://github.com/yourusername/inverse-problem-scrna
cd inverse-problem-scrna
pip install -e .
```

## Quick Start

### Scanpy Integration

```python
import scanpy as sc
import inverse_sc as isc

# Load data as usual
adata = sc.read_h5ad("data.h5ad")

# Fit inverse model
isc.pp.fit_inverse_model(adata, n_latent=30, n_programs=20)

# Now adata.obsm contains:
# - 'Z_true_mean': Inferred true expression
# - 'Z_true_std': Uncertainty estimates
# - 'program_weights': Cell state composition

# Downstream analysis on inferred biology
sc.pp.neighbors(adata, use_rep='Z_true_mean')
sc.tl.leiden(adata)
sc.tl.umap(adata)

# Quantify cluster confidence
isc.tl.cluster_uncertainty(adata)
```

### Seurat Integration (via reticulate)

```r
library(Seurat)
library(reticulate)
isc <- import("inverse_sc")

# Convert Seurat to AnnData
adata <- isc$bridge$from_seurat(seurat_obj)

# Fit inverse model
isc$pp$fit_inverse_model(adata, n_latent=30L, n_programs=20L)

# Convert back
seurat_obj <- isc$bridge$to_seurat(adata)

# Z_true_mean is now in seurat_obj@reductions$Z_true
```

## What You Get

1. **True Expression Estimates**: Posterior over pre-measurement transcriptional state
2. **Uncertainty Quantification**: Confidence in every inferred value
3. **Program Decomposition**: Interpretable transcriptional programs
4. **Robust DE Testing**: Accounts for measurement uncertainty
5. **Dropout Correction**: Physical model, not imputation

## Architecture

```
inverse_sc/
├── measurement/         # Measurement operator M
│   ├── capture.py      # Cell/gene-specific capture efficiency
│   ├── amplification.py # PCR/IVT bias modeling
│   └── sequencing.py   # Depth and sampling effects
├── inference/          # Inverse problem solver
│   ├── model.py        # Pyro generative model
│   ├── guide.py        # Variational inference
│   └── optimizer.py    # Training loop
├── preprocessing/      # Scanpy-like interface
│   ├── fit.py          # Main fitting function
│   └── calibrate.py    # Measurement model calibration
├── tools/              # Analysis functions
│   ├── uncertainty.py  # Cluster/DE confidence
│   └── programs.py     # Program interpretation
└── bridge/             # Seurat integration
    └── seurat_bridge.py
```

## Validation

See [notebooks/validation.ipynb](notebooks/validation.ipynb) for:
- Synthetic data benchmarks
- Comparison with Scanpy/scVI
- Uncertainty calibration tests
- Real data case studies

## Citation

If you use InverseSC in your research, please cite:

```bibtex
@article{inversesc2026,
  title={Inverse Problem Formulation for Single-Cell RNA Sequencing},
  author={Pritam Kumar Panda},
  journal={bioRxiv},
  year={2026}
}
```

## License

MIT License
