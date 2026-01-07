"""
Validation framework: Synthetic data generation and benchmarking.
"""

from .synthetic import generate_synthetic_data, SyntheticDataGenerator
from .benchmark import benchmark_against_scanpy, benchmark_against_scvi

__all__ = [
    "generate_synthetic_data",
    "SyntheticDataGenerator",
    "benchmark_against_scanpy",
    "benchmark_against_scvi",
]
