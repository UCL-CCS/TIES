"""
Superimpose the same molecule. This means the full superimposition should take place. Save it as another molecule.

This is because some of the molecules we have have a different atom orderthan what it should have.

Take two molecules as input. Match them fully, and then save a new molecule in the "right order".
"""
import os
import numpy as np
import shutil
import sys
import subprocess
from pathlib import Path, PurePosixPath
from generator import *
import topology_superimposer
import time

# load the .ac and .pdb file,
# find which atoms in .ac match which atoms in .pdb
# from .pdb we get the coords, and from .ac we get the charges,
# generate a new .mol2 file from the .ac with the right charges, and the order
# of atoms found in .pdb
suptop, mda_l1, mda_l2 = getSuptop('left_q.mol2', 'left_coor.pdb',
                                   ignore_charges_completely=True,
                                   ignore_bond_types=True,
                                   ignore_coords=True)
assert len(mda_l1.atoms) == len(mda_l2.atoms) == len(suptop)
# write_dual_top_pdb('morph.pdb', mda_l1, mda_l2, suptop)
# save the merged topologies as a .mol2 file
write_merged(suptop, 'morph.mol2', use_left_charges=True, use_left_coords=False)
