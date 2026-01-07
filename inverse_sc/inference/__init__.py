"""
Inverse inference module: Solves the inverse problem p(Z | X).

Given observed counts X, infer the posterior distribution over
true biological state Z.
"""

from .model import InverseModel
from .guide import InferenceGuide
from .trainer import InverseTrainer

__all__ = [
    "InverseModel",
    "InferenceGuide",
    "InverseTrainer",
]
