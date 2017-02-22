from pkg_resources import get_distribution

from .io import *
from .plots import *

__all__ = []
__version__ = get_distribution('smps').version
