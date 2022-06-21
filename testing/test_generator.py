"""
These tests focus on the generator (preprocessing of the input before applying superimpose_topologies
"""
import copy

import parmed

import ties.helpers
import ties.topology_superimposer
import ties.generator
from ties.ligand import Ligand
from ties import Pair
from ties import Config


def test_no_same_atom_names(dual_ring1):
    # make a copy of the atom list
    dual_ring_cmp = [copy.copy(a) for a in dual_ring1]
    ties.topology_superimposer.SuperimposedTopology.rename_ligands(dual_ring1, dual_ring_cmp)
    intersection = {a.name for a in dual_ring1}.intersection({a.name for a in dual_ring_cmp})
    # there should be no overlap
    assert len(intersection) == 0


def test_ligand_rename_atom_names_unique():
    lig = Ligand('data/p38_ligands_01.pdb', save=False)
    lig.correct_atom_names()
    assert len({a.name for a in lig.parmed.atoms}) == len([a.name for a in lig.parmed.atoms])


def test_are_correct_names():
    """
    Check if the pattern is correct, ie Element followed by number, C16
    :return:
    """
    lig = parmed.Structure()
    lig.add_atom(parmed.Atom(name='O1'), 'resname', 1)
    lig.add_atom(parmed.Atom(name='H1'), 'resname', 1)

    assert ties.helpers.are_correct_names([a.name for a in lig.atoms])

def test_atom_names_not_unique():
    """
    Check if the pattern is correct, ie Element followed by number, C16
    :return:
    """
    lig = parmed.Structure()
    lig.add_atom(parmed.Atom(name='O1'), 'resname', 1)
    lig.add_atom(parmed.Atom(name='H1'), 'resname', 1)
    lig.add_atom(parmed.Atom(name='O1'), 'resname', 1)

    assert ties.helpers.are_correct_names([a.name for a in lig.atoms])


# should be a test?
def test_input_prep():
    config = Config()
    config.ligand_net_charge = -1

    pair = Pair('../examples/mol2_2ligands_MCL1/l02.mol2', '../examples/mol2_2ligands_MCL1/l03.mol2', config=config)
    pair.make_atom_names_unique()
    hybrid = pair.superimpose()

    hybrid.prepare_inputs()
