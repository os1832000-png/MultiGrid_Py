"""MARFA: Molecular atmospheric Absorption with Rapid and Flexible Analysis."""
import logging

# Library-friendly: don't configure the root logger; just swallow our own
# records if the user hasn't set up logging.
logging.getLogger(__name__).addHandler(logging.NullHandler())

from .constants import Constants
from .molecules import MOLECULE_NAMES, MOLECULAR_MASSES
from .line_shapes import LineShapes
from .wing_corrections import WingCorrections
from .tips import TIPS
from .hitran import HITRANReader, SpectralLine
from .grids import GridState
from .grid_calculator import LineGridCalc
from .core import MARFA
from .profiles import AtmosphericProfile
from .pt_table import PTTableGenerator
from ._warnings import MarfaWarning, LineOutOfRangeWarning, InputFileWarning

__version__ = "3.0.2"
__all__ = [
    "Constants",
    "MOLECULE_NAMES", "MOLECULAR_MASSES",
    "LineShapes", "WingCorrections", "TIPS",
    "HITRANReader", "SpectralLine",
    "GridState", "LineGridCalc",
    "MARFA", "AtmosphericProfile", "PTTableGenerator",
    "MarfaWarning", "LineOutOfRangeWarning", "InputFileWarning",
]