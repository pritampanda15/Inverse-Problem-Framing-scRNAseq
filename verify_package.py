#!/usr/bin/env python3
"""
Verification script for InverseSC package.

Checks that all modules can be imported and basic functionality works.
"""

import sys


def test_imports():
    """Test that all modules import successfully."""
    print("Testing imports...")

    try:
        import inverse_sc
        print("✓ Main package imported")

        import inverse_sc.measurement
        print("✓ Measurement module imported")

        import inverse_sc.inference
        print("✓ Inference module imported")

        import inverse_sc.preprocessing as pp
        print("✓ Preprocessing module imported")

        import inverse_sc.tools as tl
        print("✓ Tools module imported")

        import inverse_sc.bridge
        print("✓ Bridge module imported")

        import inverse_sc.validation
        print("✓ Validation module imported")

        return True
    except ImportError as e:
        print(f"✗ Import failed: {e}")
        return False


def test_basic_functionality():
    """Test basic functionality."""
    print("\nTesting basic functionality...")

    try:
        from inverse_sc.validation import generate_synthetic_data

        # Generate small synthetic dataset
        print("  Generating synthetic data...")
        adata, truth = generate_synthetic_data(
            n_cells=50,
            n_genes=100,
            n_programs=3,
            seed=42
        )
        print(f"  ✓ Generated data: {adata.shape}")

        # Check ground truth
        assert 'Z_true' in truth
        assert truth['Z_true'].shape == (50, 100)
        print("  ✓ Ground truth verified")

        return True
    except Exception as e:
        print(f"  ✗ Functionality test failed: {e}")
        return False


def test_measurement_model():
    """Test measurement model."""
    print("\nTesting measurement model...")

    try:
        import torch
        from inverse_sc.measurement import MeasurementOperator

        op = MeasurementOperator(n_genes=100, n_cells=50)
        Z = torch.randn(50, 100).abs()
        X = op.forward(Z)

        assert X.shape == (50, 100)
        print("  ✓ Measurement operator works")

        return True
    except Exception as e:
        print(f"  ✗ Measurement model test failed: {e}")
        return False


def test_inference_model():
    """Test inference model."""
    print("\nTesting inference model...")

    try:
        from inverse_sc.inference import InverseModel, InferenceGuide

        model = InverseModel(
            n_genes=100,
            n_cells=50,
            n_latent=10,
            n_programs=5
        )
        print("  ✓ Inverse model created")

        guide = InferenceGuide(
            n_genes=100,
            n_latent=10,
            n_programs=5
        )
        print("  ✓ Inference guide created")

        return True
    except Exception as e:
        print(f"  ✗ Inference model test failed: {e}")
        return False


def check_file_structure():
    """Check that key files exist."""
    print("\nChecking file structure...")

    import os

    required_files = [
        'README.md',
        'setup.py',
        'requirements.txt',
        'LICENSE',
        'inverse_sc/__init__.py',
        'inverse_sc/measurement/__init__.py',
        'inverse_sc/inference/__init__.py',
        'inverse_sc/preprocessing/__init__.py',
        'inverse_sc/tools/__init__.py',
        'inverse_sc/bridge/__init__.py',
        'inverse_sc/validation/__init__.py',
        'tests/test_basic.py',
        'docs/METHODOLOGY.md',
        'docs/API_REFERENCE.md',
        'notebooks/01_quickstart.ipynb',
        'examples/basic_usage.py',
    ]

    all_exist = True
    for file in required_files:
        if os.path.exists(file):
            print(f"  ✓ {file}")
        else:
            print(f"  ✗ {file} missing")
            all_exist = False

    return all_exist


def main():
    """Run all verification tests."""
    print("="*60)
    print("InverseSC Package Verification")
    print("="*60)

    results = []

    # Test imports
    results.append(("Imports", test_imports()))

    # Test basic functionality
    results.append(("Basic Functionality", test_basic_functionality()))

    # Test measurement model
    results.append(("Measurement Model", test_measurement_model()))

    # Test inference model
    results.append(("Inference Model", test_inference_model()))

    # Check file structure
    results.append(("File Structure", check_file_structure()))

    # Summary
    print("\n" + "="*60)
    print("VERIFICATION SUMMARY")
    print("="*60)

    for test_name, passed in results:
        status = "PASS" if passed else "FAIL"
        symbol = "✓" if passed else "✗"
        print(f"{symbol} {test_name}: {status}")

    all_passed = all(result[1] for result in results)

    if all_passed:
        print("\n✓ All tests passed! Package is ready to use.")
        return 0
    else:
        print("\n✗ Some tests failed. Please check the output above.")
        return 1


if __name__ == "__main__":
    sys.exit(main())
