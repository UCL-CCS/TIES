#!/usr/bin/env python2
"""
Load the requested .pdb file
Load the .fep data about the atoms
Update the .pdb file
"""

import MDAnalysis as mda
import json

u = mda.Universe('merged_solvated.pdb')
fep_meta_str = open('l18_l39.fep').read()
fep_meta = json.loads(fep_meta_str)

left_matched = list(fep_meta['matching'].keys())
appearing_atoms = fep_meta['appearing']
disappearing_atoms = fep_meta['disappearing']

# update the Temp column
for atom in u.atoms:
    # if the atom was "matched", meaning present in both ligands (left and right)
    # then ignore
    # note: we only use the left ligand
    if atom.name in left_matched:
        continue
    elif atom.name in appearing_atoms:
        # appearing atoms should
        atom.tempfactor = 1
    elif atom.name in disappearing_atoms:
        atom.tempfactor = -1

u.atoms.write('merged_solvated_fep.pdb')
