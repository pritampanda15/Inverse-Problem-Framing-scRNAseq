"""
Basic tests for InverseSC package.
"""

import pytest
import numpy as np
import torch
from anndata import AnnData


def test_imports():
    """Test that all modules can be imported."""
    import inverse_sc
    import inverse_sc.measurement
    import inverse_sc.inference
    import inverse_sc.preprocessing as pp
    import inverse_sc.tools as tl
    import inverse_sc.bridge
    import inverse_sc.validation


def test_synthetic_data_generation():
    """Test synthetic data generation."""
    from inverse_sc.validation import generate_synthetic_data

    adata, truth = generate_synthetic_data(
        n_cells=100,
        n_genes=200,
        n_programs=3,
        seed=42,
    )

    assert adata.n_obs == 100
    assert adata.n_vars == 200
    assert 'Z_true' in truth
    assert truth['Z_true'].shape == (100, 200)


def test_measurement_operator():
    """Test measurement operator."""
    from inverse_sc.measurement import MeasurementOperator

    n_cells, n_genes = 50, 100

    op = MeasurementOperator(n_genes=n_genes, n_cells=n_cells)

    # Test forward pass
    Z_true = torch.randn(n_cells, n_genes).abs()
    X_obs = op.forward(Z_true)

    assert X_obs.shape == (n_cells, n_genes)
    assert torch.all(X_obs >= 0)  # Counts are non-negative


def test_inverse_model():
    """Test inverse model creation."""
    from inverse_sc.inference import InverseModel

    model = InverseModel(
        n_genes=100,
        n_cells=50,
        n_latent=10,
        n_programs=5,
    )

    assert model.n_genes == 100
    assert model.n_programs == 5


def test_inference_guide():
    """Test inference guide."""
    from inverse_sc.inference import InferenceGuide

    guide = InferenceGuide(
        n_genes=100,
        n_latent=10,
        n_programs=5,
    )

    assert guide.n_genes == 100


def test_fit_inverse_model_synthetic():
    """Test fitting on synthetic data."""
    from inverse_sc.validation import generate_synthetic_data
    from inverse_sc.preprocessing import fit_inverse_model

    # Generate small dataset
    adata, truth = generate_synthetic_data(
        n_cells=100,
        n_genes=200,
        n_programs=3,
        seed=42,
    )

    # Fit with minimal epochs for speed
    fit_inverse_model(
        adata,
        n_latent=10,
        n_programs=3,
        n_epochs=5,  # Minimal for testing
        batch_size=50,
    )

    # Check outputs exist
    assert 'Z_true_mean' in adata.obsm
    assert 'Z_true_std' in adata.obsm
    assert 'program_weights' in adata.obsm

    # Check shapes
    assert adata.obsm['Z_true_mean'].shape == (100, 200)
    assert adata.obsm['program_weights'].shape == (100, 3)


def test_cluster_uncertainty():
    """Test uncertainty quantification."""
    from inverse_sc.validation import generate_synthetic_data
    from inverse_sc.preprocessing import fit_inverse_model
    from inverse_sc.tools import cluster_uncertainty

    # Generate and fit
    adata, _ = generate_synthetic_data(n_cells=100, n_genes=200, seed=42)
    fit_inverse_model(adata, n_epochs=5)

    # Add fake cluster labels
    adata.obs['leiden'] = np.random.randint(0, 3, size=100).astype(str)

    # Compute uncertainty
    cluster_uncertainty(adata)

    # Check outputs
    assert 'cluster_confidence' in adata.obs.columns
    assert 'cluster_distance' in adata.obs.columns

    # Confidence should be between 0 and 1
    assert adata.obs['cluster_confidence'].min() >= 0
    assert adata.obs['cluster_confidence'].max() <= 1


def test_differential_expression():
    """Test differential expression with uncertainty."""
    from inverse_sc.validation import generate_synthetic_data
    from inverse_sc.preprocessing import fit_inverse_model
    from inverse_sc.tools import differential_expression_robust

    # Generate and fit
    adata, _ = generate_synthetic_data(n_cells=100, n_genes=200, seed=42)
    fit_inverse_model(adata, n_epochs=5)

    # Add fake group labels
    adata.obs['group'] = np.array(['A'] * 50 + ['B'] * 50)

    # Test DE
    de_results = differential_expression_robust(
        adata,
        group_key='group',
        group1='A',
        group2='B',
    )

    # Check outputs
    assert 'gene_names' in de_results
    assert 'log_fold_change' in de_results
    assert 'p_adj' in de_results
    assert 'confident' in de_results
    assert len(de_results['gene_names']) == 200


def test_benchmark():
    """Test benchmarking functions."""
    from inverse_sc.validation import (
        generate_synthetic_data,
        benchmark_against_scanpy,
        uncertainty_calibration,
    )
    from inverse_sc.preprocessing import fit_inverse_model

    # Generate and fit
    adata, truth = generate_synthetic_data(n_cells=50, n_genes=100, seed=42)
    fit_inverse_model(adata, n_epochs=5)

    # Benchmark
    results = benchmark_against_scanpy(adata, truth, run_scanpy=False)
    assert len(results) > 0
    assert 'method' in results.columns

    # Calibration
    calib = uncertainty_calibration(adata, truth)
    assert 'coverage_1std' in calib
    assert 'calibration_score' in calib


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
