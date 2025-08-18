"""
These tests focus on the Ligand
"""

import rdkit.Chem

from ties.ligand import Ligand


def test_correct_atom_names_HAY():
    # atom names do not follow the right format
    lig = Ligand("data/l_HAY.mol2")
    assert lig.are_atom_names_correct()


def test_correct_atom_names_HAY_renamed():
    # atom names do not follow the right format
    lig = Ligand("data/l_HAY.mol2")
    lig.correct_atom_names()
    assert lig.are_atom_names_correct()


def test_ligand_from_rdkit_mol():
    rdmol = rdkit.Chem.MolFromSmiles("COO")
    lig = Ligand(rdmol)
    # check if the file was saved into an SDF
    assert lig.current.exists()
