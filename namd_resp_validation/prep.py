"""
Superimpose the same molecule. This means the full superimposition should take place. Save it as another molecule.

This is because some of the molecules we have have a different atom orderthan what it should have.

Take two molecules as input. Match them fully, and then save a new molecule in the "right order".
"""
from generator import getSuptop, write_merged
import MDAnalysis as mda

# load the q from .mol2 and coords from .pdb file,
suptop, mda_l1, mda_l2 = getSuptop('left_q.mol2', 'left_coor.pdb',
                                   ignore_charges_completely=True,
                                   ignore_bond_types=True,
                                   left_coords_are_ref=False,
                                   align_molecules=False,
                                   use_only_gentype=True)
assert len(mda_l1.atoms) == len(mda_l2.atoms) == len(suptop)
write_merged(suptop, 'new_left.mol2', use_left_charges=True)

suptop, mda_l1, mda_l2 = getSuptop('right_q.mol2', 'right_coor.pdb',
                                   ignore_charges_completely=True,
                                   ignore_bond_types=True,
                                   left_coords_are_ref=False, # the .pdb is the coords
                                   align_molecules=True, # align at the end to the right molecule
                                   use_only_gentype=True)
assert len(mda_l1.atoms) == len(mda_l2.atoms) == len(suptop)
write_merged(suptop, 'new_right.mol2', use_left_charges=True) # use the charges from the left
