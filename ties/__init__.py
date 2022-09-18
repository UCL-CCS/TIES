from pathlib import Path

from . import cli

from .ligand import Ligand
from .pair import Pair
from .config import Config
from .ligandmap import LigandMap
from .protein import Protein
from .md import MD

__version__ = open(Path(__file__).parent / 'version.txt').read().strip()

__all__ = [Ligand, Protein, Pair, Config, LigandMap, Protein, MD, __version__]