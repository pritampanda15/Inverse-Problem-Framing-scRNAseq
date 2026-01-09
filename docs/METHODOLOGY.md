# Methodology: Inverse Problem Framework for scRNA-seq

## The Core Problem

### What Traditional Methods Do

Standard scRNA-seq analysis pipelines (Scanpy, Seurat) follow this pattern:

```
Raw Counts X → Normalize → Reduce Dimensions → Cluster → Interpret
```

The implicit assumption is that the observed count matrix **X** is a reasonable proxy for biology after simple transformations (log-normalization, scaling).

### The Hidden Reality

But **X** is not biology. It's the output of a destructive measurement process:

```
True Biology Z → Cell Lysis → Capture (dropout) → Amplification (bias) → Sequencing (sampling) → Observed X
```

Each step introduces distortion:
- **Capture**: Only ~10% of mRNA molecules are captured (gene- and cell-specific)
- **Amplification**: PCR introduces exponential bias (GC content, length effects)
- **Sequencing**: Multinomial sampling to fixed depth creates noise

Traditional methods essentially say: "Hope the distortion isn't too bad."

## The Inverse Problem Approach

Instead of analyzing **X** directly, we ask:

> **"Given that I observed X, what biological state Z must have existed before measurement?"**

This is an **inverse problem**: infer the input (Z) from the output (X) of a known process (M).

### Mathematical Formulation

**Forward Model (Measurement):**
```
X = M(Z) + noise
```

Where:
- **Z**: True transcriptional state (what we want)
- **M**: Measurement operator (capture, amplification, sequencing)
- **X**: Observed counts (what we have)

**Inverse Problem (Inference):**
```
p(Z | X) ∝ p(X | Z) × p(Z)
```

- **p(X | Z)**: Likelihood under measurement model (physics)
- **p(Z)**: Prior on biological states (biology)
- **p(Z | X)**: Posterior (what we infer)

## Key Differences from Existing Methods

### vs. Standard Normalization (Scanpy/Seurat)

| Aspect | Standard | InverseSC |
|--------|----------|-----------|
| What you analyze | log(X + 1) | Posterior over Z |
| Dropout handling | Ignore or impute | Explicit uncertainty |
| Output | Point estimates | Distributions |
| Biological assumptions | Implicit (in normalization) | Explicit (in prior) |

### vs. scVI (VAE approaches)

| Aspect | scVI | InverseSC |
|--------|------|-----------|
| Decoder | Neural network (learned) | Physical model (fixed) |
| What's learned | Arbitrary reconstruction | Biological state |
| Interpretability | Latent space is abstract | Latent = transcriptional programs |
| Zeros | ZINB parameter | Consequence of capture physics |

**Critical difference**: scVI asks "What latent representation reconstructs X well?"
InverseSC asks "What biology, passed through realistic measurement, produces X?"

## Implementation Details

### 1. Measurement Operator M

We explicitly model three stages:

**Capture Efficiency:**
```python
Z_captured = Z_true * p_capture
```
where `p_capture` is gene- and cell-specific, calibrated from dropout patterns.

**Amplification Bias:**
```python
Z_amplified = Z_captured * amplification_factor
```
Accounts for PCR bias (GC content, length).

**Sequencing:**
```python
X ~ Multinomial(library_size, Z_amplified / sum(Z_amplified))
```
Multinomial sampling to fixed depth.

### 2. Biological Prior p(Z)

We use a **transcriptional program** model:

```python
# Each cell is a mixture of programs
program_weights ~ Dirichlet(α)

# Each program has a gene signature
Z = Σ_k program_weights[k] * program_signature[k]

# With biological variability
Z_true ~ LogNormal(Z, σ_bio)
```

This encodes that:
- Cells are compositional mixtures of states
- Gene expression is multiplicative (log-normal)
- Programs are interpretable (gene sets)

### 3. Variational Inference

Since exact inference is intractable, we use **variational inference**:

**Encoder (Guide):**
```python
q(Z | X) = Neural_Network(X) → parameters of q
```

**ELBO Optimization:**
```
L = E_q[log p(X, Z)] - E_q[log q(Z)]
```

Train via stochastic gradient ascent on ELBO.

## What You Get

Unlike standard methods, InverseSC provides:

### 1. True Expression Estimates

`adata.obsm['Z_true_mean']`: Posterior mean of pre-measurement expression
`adata.obsm['Z_true_std']`: Uncertainty estimates

Use Z_true_mean for downstream analysis instead of log-normalized counts.

### 2. Uncertainty Quantification

Every inferred value has a confidence estimate:
- Cluster assignments have confidence scores
- DE genes are "confidently differential" or not
- UMAP can show uncertainty clouds

### 3. Interpretable Programs

`adata.obsm['program_weights']`: Each cell's composition of transcriptional programs
Programs are interpretable gene sets (unlike PCA components).

### 4. Calibrated Dropout Correction

Not imputation (filling zeros with guesses), but proper uncertainty:
- If a gene has high dropout, we admit we're uncertain
- Downstream analysis accounts for this

## Validation Strategy

### Synthetic Data

Generate data with known ground truth:
1. Sample true Z from program model
2. Apply realistic measurement model → observe X
3. Run inverse inference: X → inferred Z
4. Compare inferred Z to true Z

**Metrics:**
- Correlation between true and inferred Z
- Program recovery accuracy
- Uncertainty calibration (is uncertainty correct?)

### Real Data Benchmarks

Compare to standard methods:
1. Cluster stability
2. Marker gene detection
3. Batch correction quality
4. Agreement with known cell type annotations

## When to Use InverseSC

**Use InverseSC when:**
- High dropout is a concern (sparse data)
- You need uncertainty estimates (clinical decisions, rare cells)
- Comparing across very different sequencing depths
- You want interpretable cell states (programs)

**Stick with standard methods when:**
- Data is very deep (>50k reads/cell, <30% dropout)
- You just need quick exploratory analysis
- Computational resources are limited

## Computational Considerations

**Runtime:**
- Fit time: ~5-10 minutes for 10k cells, 2k genes (GPU)
- Scales roughly linearly with number of cells
- Recommend GPU for >5k cells

**Memory:**
- Model size: ~10MB per 1k genes
- Peak memory: ~2GB for 10k cells

**Parallelization:**
- Mini-batch training enables large datasets
- Multiple chains for uncertainty quantification

## Future Directions

### Methodological Extensions

1. **Spatial transcriptomics**: Add spatial prior on Z
2. **Multi-modal data**: Joint measurement model for CITE-seq
3. **Temporal dynamics**: Time-dependent measurement model
4. **Batch effects**: Explicit batch parameters in M

### Biological Priors

1. **Gene regulatory networks**: Replace programs with GRN prior
2. **Pathway-based**: Prior structured by GO/KEGG pathways
3. **Hierarchical**: Nested program structure (cell types → states)

### Computational Improvements

1. **Amortized inference**: Faster per-cell inference
2. **Importance sampling**: Better uncertainty estimates
3. **Sparse computation**: Efficient for very high-dimensional Z

## References & Related Work

**Inverse Problems in Biology:**
- Classical deconvolution problems
- Image reconstruction in microscopy
- Protein structure from NMR

**scRNA-seq Modeling:**
- scVI (Lopez et al. 2018): VAE approach
- SAVER (Huang et al. 2018): Gene-by-gene Bayesian
- DCA (Eraslan et al. 2019): Autoencoder with ZINB

**Key Distinction:**
Those methods learn arbitrary reconstructions. We fix the decoder to measurement physics.

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

## Contact

For questions, issues, or contributions:
- GitHub: [inverse-problem-scrna](https://github.com/yourusername/inverse-problem-scrna)
- Email: pritam@stanford.edu

---

*This is a novel methodological framework. We welcome feedback and collaboration!*
