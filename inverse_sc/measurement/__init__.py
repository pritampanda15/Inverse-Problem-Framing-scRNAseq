"""
Measurement operator module: Models the physics of scRNA-seq measurement.

This module implements the forward model M: Z → X, where:
- Z is the true biological transcriptional state
- X is the observed count matrix
- M encapsulates capture, amplification, and sequencing
"""

from .operator import MeasurementOperator
from .capture import CaptureModel
from .amplification import AmplificationModel
from .sequencing import SequencingModel

__all__ = [
    "MeasurementOperator",
    "CaptureModel",
    "AmplificationModel",
    "SequencingModel",
]
