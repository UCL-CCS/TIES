"""
These tests focus on the generator (preprocessing of the input before applying superimpose_topologies
"""
import copy

import MDAnalysis

import ties.helpers
import ties.topology_superimposer
import ties.generator
from ties.ligand import Ligand


def test_no_same_atom_names(dual_ring1):
    # make a copy of the atom list
    dual_ring_cmp = [copy.copy(a) for a in dual_ring1]
    ties.topology_superimposer.SuperimposedTopology.rename_ligands(dual_ring1, dual_ring_cmp)
    intersection = {a.name for a in dual_ring1}.intersection({a.name for a in dual_ring_cmp})
    # there should be no overlap
    assert len(intersection) == 0


def test_rename_ligand():
    lig = Ligand('data/p38_ligands_01.pdb', save=False)
    lig.make_atom_names_correct()
    assert len(set(lig.universe.atoms.names)) == len(lig.universe.atoms.names)


def test_are_correct_names():
    u = MDAnalysis.Universe.empty(n_atoms=2)
    # set the atom names
    u.add_TopologyAttr('name', ['O1', 'H1'])
    assert ties.helpers.are_correct_names(u.atoms.names)