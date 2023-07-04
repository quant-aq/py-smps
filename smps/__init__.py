from pkg_resources import get_distribution

from .io import *
from .plots import *
from .models import *
from .rcmod import *

__version__ = get_distribution("py-smps").version