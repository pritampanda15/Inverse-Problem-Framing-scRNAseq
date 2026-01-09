# API Reference

Complete reference for InverseSC functions and classes.

## Preprocessing (`inverse_sc.pp`)

### fit_inverse_model

```python
inverse_sc.pp.fit_inverse_model(
    adata,
    n_latent=30,
    n_programs=20,
    n_epochs=100,
    batch_size=256,
    learning_rate=0.001,
    use_cuda=False,
    layer=None,
    copy=False,
)
```

Fit inverse problem model to scRNA-seq data.

**Parameters:**
- `adata` (AnnData): Annotated data matrix
- `n_latent` (int): Latent dimension for inference network
- `n_programs` (int): Number of transcriptional programs
- `n_epochs` (int): Training epochs
- `batch_size` (int, optional): Mini-batch size
- `learning_rate` (float): Learning rate
- `use_cuda` (bool): Use GPU if available
- `layer` (str, optional): Which layer to use (None = .X)
- `copy` (bool): Return copy instead of modifying in-place

**Returns:**
- `adata` (AnnData, optional): If copy=True

**Adds to adata:**
- `adata.obsm['Z_true_mean']`: Inferred expression (posterior mean)
- `adata.obsm['Z_true_std']`: Uncertainty (posterior std)
- `adata.obsm['program_weights']`: Cell state composition
- `adata.uns['inverse_model']`: Training history and metadata

**Example:**
```python
import scanpy as sc
import inverse_sc as isc

adata = sc.read_h5ad("data.h5ad")
isc.pp.fit_inverse_model(adata, n_epochs=200)
```

---

### calibrate_measurement_model

```python
inverse_sc.pp.calibrate_measurement_model(
    adata,
    spike_in_genes=None,
    spike_in_amounts=None,
    layer=None,
)
```

Calibrate measurement model from data.

**Parameters:**
- `adata` (AnnData): Data
- `spike_in_genes` (list, optional): Names of spike-in genes (e.g., ERCC)
- `spike_in_amounts` (array, optional): Known input amounts
- `layer` (str, optional): Which layer to use

**Returns:**
- `calibration_results` (dict): Estimated measurement parameters

**Example:**
```python
calib = isc.pp.calibrate_measurement_model(adata)
adata.var['capture_efficiency'] = calib['capture']['capture_prob_estimate']
```

---

### basic_qc

```python
inverse_sc.pp.basic_qc(
    adata,
    min_genes=200,
    min_cells=3,
    max_pct_mito=20.0,
    mito_prefix="MT-",
    copy=False,
)
```

Basic quality control filtering.

**Parameters:**
- `adata` (AnnData): Data
- `min_genes` (int): Minimum genes per cell
- `min_cells` (int): Minimum cells per gene
- `max_pct_mito` (float): Maximum mitochondrial percentage
- `mito_prefix` (str): Prefix for mitochondrial genes
- `copy` (bool): Return copy

---

## Tools (`inverse_sc.tl`)

### cluster_uncertainty

```python
inverse_sc.tl.cluster_uncertainty(
    adata,
    cluster_key='leiden',
    z_key='Z_true_mean',
    uncertainty_key='Z_true_std',
)
```

Quantify confidence in cluster assignments.

**Parameters:**
- `adata` (AnnData): Data with cluster labels
- `cluster_key` (str): Column in .obs with cluster labels
- `z_key` (str): Key in .obsm with inferred expression
- `uncertainty_key` (str): Key in .obsm with uncertainty

**Adds to adata.obs:**
- `cluster_distance`: Distance to cluster centroid
- `cluster_confidence`: Confidence score (0-1)

**Example:**
```python
isc.tl.cluster_uncertainty(adata)
uncertain = adata.obs['cluster_confidence'] < 0.5
```

---

### differential_expression_robust

```python
inverse_sc.tl.differential_expression_robust(
    adata,
    group_key,
    group1,
    group2,
    z_key='Z_true_mean',
    uncertainty_key='Z_true_std',
    min_fold_change=1.5,
)
```

Differential expression accounting for uncertainty.

**Parameters:**
- `adata` (AnnData): Data
- `group_key` (str): Column in .obs with group labels
- `group1` (str): First group
- `group2` (str): Second group
- `z_key` (str): Key for inferred expression
- `uncertainty_key` (str): Key for uncertainty
- `min_fold_change` (float): Minimum fold change

**Returns:**
- `results` (dict): DE results with uncertainty-aware statistics
  - `gene_names`: Gene names
  - `log_fold_change`: Log2 fold change
  - `p_adj`: Adjusted p-values
  - `confident`: Boolean mask of confident DE genes

**Example:**
```python
de = isc.tl.differential_expression_robust(
    adata,
    group_key='cell_type',
    group1='T_cell',
    group2='B_cell'
)
confident_de = de['gene_names'][de['confident']]
```

---

### interpret_programs

```python
inverse_sc.tl.interpret_programs(
    adata,
    program_key='program_weights',
    top_genes=50,
)
```

Interpret learned transcriptional programs.

**Parameters:**
- `adata` (AnnData): Data with fitted model
- `program_key` (str): Key in .obsm with program weights
- `top_genes` (int): Number of top genes per program

**Returns:**
- `program_info` (dict): Dict mapping program index to DataFrame of top genes

**Example:**
```python
programs = isc.tl.interpret_programs(adata)
print(programs[0].head(10))  # Top genes for program 0
```

---

## Validation (`inverse_sc.validation`)

### generate_synthetic_data

```python
inverse_sc.validation.generate_synthetic_data(
    n_cells=1000,
    n_genes=2000,
    n_programs=5,
    capture_efficiency=0.1,
    dropout_rate=0.7,
    sequencing_depth=10000,
    seed=None,
)
```

Generate synthetic scRNA-seq data with known ground truth.

**Parameters:**
- `n_cells` (int): Number of cells
- `n_genes` (int): Number of genes
- `n_programs` (int): Number of transcriptional programs
- `capture_efficiency` (float): Mean capture probability
- `dropout_rate` (float): Target dropout rate
- `sequencing_depth` (int): Mean library size
- `seed` (int, optional): Random seed

**Returns:**
- `adata` (AnnData): Observed counts
- `ground_truth` (dict): True biological state
  - `Z_true`: True expression
  - `program_weights`: True program composition
  - `program_signatures`: Program gene signatures

**Example:**
```python
adata, truth = isc.validation.generate_synthetic_data(
    n_cells=500,
    seed=42
)
# Fit model
isc.pp.fit_inverse_model(adata)
# Compare
Z_inferred = adata.obsm['Z_true_mean']
correlation = np.corrcoef(Z_inferred.flatten(), truth['Z_true'].flatten())[0,1]
```

---

### benchmark_against_scanpy

```python
inverse_sc.validation.benchmark_against_scanpy(
    adata,
    ground_truth,
    run_scanpy=True,
)
```

Compare inverse method against standard Scanpy pipeline.

**Parameters:**
- `adata` (AnnData): Data with inverse inference results
- `ground_truth` (dict): True biological state
- `run_scanpy` (bool): Whether to run scanpy for comparison

**Returns:**
- `results` (DataFrame): Comparison metrics

**Example:**
```python
results = isc.validation.benchmark_against_scanpy(adata, truth)
print(results[['method', 'global_correlation', 'rmse']])
```

---

### uncertainty_calibration

```python
inverse_sc.validation.uncertainty_calibration(
    adata,
    ground_truth,
)
```

Check if uncertainty estimates are well-calibrated.

**Parameters:**
- `adata` (AnnData): Data with uncertainty estimates
- `ground_truth` (dict): Ground truth

**Returns:**
- `calibration` (dict): Calibration metrics
  - `coverage_1std`: Fraction within 1 std (should be ~0.68)
  - `coverage_2std`: Fraction within 2 std (should be ~0.95)
  - `calibration_score`: Overall calibration quality

---

## Bridge (`inverse_sc.bridge`)

### from_seurat

```python
inverse_sc.bridge.from_seurat(
    seurat_obj,
    assay="RNA",
    use_raw=False,
)
```

Convert Seurat object to AnnData.

**Parameters:**
- `seurat_obj` (rpy2.robjects): Seurat object from R
- `assay` (str): Which assay to extract
- `use_raw` (bool): Use raw counts or normalized

**Returns:**
- `adata` (AnnData): Converted data

**Example (in Python via reticulate):**
```python
import inverse_sc as isc
adata = isc.bridge.from_seurat(seu)
```

---

### to_seurat

```python
inverse_sc.bridge.to_seurat(
    adata,
    z_key='Z_true_mean',
    program_key='program_weights',
)
```

Convert AnnData back to Seurat.

**Parameters:**
- `adata` (AnnData): Data with inverse results
- `z_key` (str): Key for inferred expression
- `program_key` (str): Key for program weights

**Returns:**
- `seurat_obj` (rpy2.robjects): Seurat object

---

## Measurement Models

### MeasurementOperator

```python
from inverse_sc.measurement import MeasurementOperator

operator = MeasurementOperator(
    n_genes=2000,
    n_cells=1000,
    model_capture=True,
    model_amplification=True,
    model_sequencing=True,
)
```

Complete measurement model for scRNA-seq.

**Methods:**
- `forward(Z_true)`: Apply measurement model
- `calibrate(X_observed)`: Calibrate from data
- `effective_sensitivity(Z_true)`: Compute sensitivity

---

## Inference Models

### InverseModel

```python
from inverse_sc.inference import InverseModel

model = InverseModel(
    n_genes=2000,
    n_cells=1000,
    n_latent=30,
    n_programs=20,
)
```

Generative model for inverse problem.

**Methods:**
- `model(X_obs)`: Pyro generative model
- `sample_prior(n_samples)`: Sample from prior

---

### InferenceGuide

```python
from inverse_sc.inference import InferenceGuide

guide = InferenceGuide(
    n_genes=2000,
    n_latent=30,
    n_programs=20,
)
```

Variational guide (encoder network).

**Methods:**
- `guide(X_obs)`: Pyro guide function

---

### InverseTrainer

```python
from inverse_sc.inference import InverseTrainer

trainer = InverseTrainer(
    model=model,
    guide=guide,
    learning_rate=0.001,
)
```

Training loop for inference.

**Methods:**
- `train(X_obs, n_epochs)`: Train model
- `predict(X_obs)`: Get posterior predictions
- `save_checkpoint(path)`: Save model

---

## Common Workflows

### Workflow 1: Basic Analysis

```python
import scanpy as sc
import inverse_sc as isc

# Load and QC
adata = sc.read_h5ad("data.h5ad")
isc.pp.basic_qc(adata)

# Fit model
isc.pp.fit_inverse_model(adata, n_epochs=100)

# Downstream analysis
sc.pp.neighbors(adata, use_rep='Z_true_mean')
sc.tl.leiden(adata)
sc.tl.umap(adata)

# Uncertainty
isc.tl.cluster_uncertainty(adata)
```

### Workflow 2: Validation

```python
# Generate synthetic data
adata, truth = isc.validation.generate_synthetic_data()

# Fit
isc.pp.fit_inverse_model(adata)

# Benchmark
results = isc.validation.benchmark_against_scanpy(adata, truth)
calib = isc.validation.uncertainty_calibration(adata, truth)
```

### Workflow 3: Program Analysis

```python
# Fit model
isc.pp.fit_inverse_model(adata, n_programs=15)

# Interpret
programs = isc.tl.interpret_programs(adata)

# Visualize
for prog_idx in range(5):
    sc.pl.umap(adata, color=adata.obsm['program_weights'][:, prog_idx])
```
