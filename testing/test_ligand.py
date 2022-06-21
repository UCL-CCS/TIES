"""
These tests focus on the Ligand
"""
import parmed

from ties.ligand import Ligand



def test_correct_atom_names_HAY():
    # atom names do not follow the right format
    lig = Ligand('data/l_HAY.mol2')
    assert lig.are_atom_names_correct()


def test_correct_atom_names_HAY_renamed():
    # atom names do not follow the right format
    lig = Ligand('data/l_HAY.mol2')
    lig.correct_atom_names()
    assert lig.are_atom_names_correct()
