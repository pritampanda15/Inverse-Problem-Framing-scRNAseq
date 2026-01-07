"""
InverseSC: Inverse Problem Framework for Single-Cell RNA-seq

A novel approach that inverts the measurement process rather than
analyzing distorted observations directly.
"""

__version__ = "0.1.0"
__author__ = "Pritam Kumar Panda"

from . import measurement
from . import inference
from . import preprocessing as pp
from . import tools as tl
from . import bridge

__all__ = ["measurement", "inference", "pp", "tl", "bridge"]
