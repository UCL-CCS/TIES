"""
Focuses on the _overlay function that actually traverses the molecule using the given starting points.
"""

from topology_superimposer import SuperimposedTopology, get_charges, \
    superimpose_topologies, _superimpose_topologies, assign_coords_from_pdb, \
    AtomNode, _overlay
import networkx as nx
from os import path
import numpy as np


def test_2diff_atoms_wrongStart():
    """
    create a simple molecule chain with an ester
     LIGAND 1        LIGAND 2
        C1              C11
        \                \
        N1              N11
    In this there is only one solution.
    """
    # construct the LIGAND 1
    # ignore the third dimension
    # c1
    c1 = AtomNode(1, 'C1', 'TEST', 1, 0, 'C')
    c1.set_coords(np.array([1, 1, 0], dtype='float32'))
    n1 = AtomNode(2, 'N1', 'TEST', 1, 0, 'N')
    n1.set_coords(np.array([1, 2, 0], dtype='float32'))
    c1.bindTo(n1)

    # construct the LIGAND 2
    # ignore the third dimension
    # c1
    c11 = AtomNode(11, 'C11', 'TEST', 1, 0, 'C')
    c11.set_coords(np.array([1, 1, 0], dtype='float32'))
    n11 = AtomNode(11, 'N11', 'TEST', 1, 0, 'N')
    n11.set_coords(np.array([1, 2, 0], dtype='float32'))
    c11.bindTo(n11)

    # should return a list with an empty sup_top
    suptops = _overlay(c1, n11, atol=9999)
    # it should return an empty suptop
    assert len(suptops) == 1
    assert len(suptops[0]) == 0


def test_2diff_atoms_rightStart():
    """
    Two different Atoms. Only one solution exists.

     LIGAND 1        LIGAND 2
        C1              C11
        \                \
        N1              N11
    """
    # construct the LIGAND 1
    # ignore the third dimension
    # c1
    c1 = AtomNode(1, 'C1', 'TEST', 1, 0, 'C')
    c1.set_coords(np.array([1, 1, 0], dtype='float32'))
    n1 = AtomNode(2, 'N1', 'TEST', 1, 0, 'N')
    n1.set_coords(np.array([1, 2, 0], dtype='float32'))
    c1.bindTo(n1)

    # construct the LIGAND 2
    # ignore the third dimension
    # c1
    c11 = AtomNode(11, 'C11', 'TEST', 1, 0, 'C')
    c11.set_coords(np.array([1, 1, 0], dtype='float32'))
    n11 = AtomNode(11, 'N11', 'TEST', 1, 0, 'N')
    n11.set_coords(np.array([1, 2, 0], dtype='float32'))
    c11.bindTo(n11)

    # should overlap 2 atoms
    suptops = _overlay(c1, c11, atol=9999)
    assert len(suptops) == 1
    suptop = suptops[0]

    # the number of overlapped atoms is two
    assert len(suptop) == 2
    correct_overlaps = [('C1', 'C11'), ('N1', 'N11')]
    for atomName1, atomName2 in correct_overlaps:
        assert suptop.contains_atomNamePair(atomName1, atomName2)


def test_3diff_atoms_rightStart():
    """
    Only one solution exists.

     LIGAND 1        LIGAND 2
        C1              C11
        \                \
        N1              N11
        /                /
        O1              O11
    """
    # construct the LIGAND 1
    # ignore the third dimension
    # c1
    c1 = AtomNode(1, 'C1', 'TEST', 1, 0, 'C')
    c1.set_coords(np.array([1, 1, 0], dtype='float32'))
    n1 = AtomNode(2, 'N1', 'TEST', 1, 0, 'N')
    n1.set_coords(np.array([1, 2, 0], dtype='float32'))
    c1.bindTo(n1)

    # construct the LIGAND 2
    # ignore the third dimension
    # c1
    c11 = AtomNode(11, 'C11', 'TEST', 1, 0, 'C')
    c11.set_coords(np.array([1, 1, 0], dtype='float32'))
    n11 = AtomNode(11, 'N11', 'TEST', 1, 0, 'N')
    n11.set_coords(np.array([1, 2, 0], dtype='float32'))
    c11.bindTo(n11)

    # should overlap 2 atoms
    suptops = _overlay(c1, c11, atol=9999)
    assert len(suptops) == 1
    suptop = suptops[0]

    # the number of overlapped atoms is two
    assert len(suptop) == 2
    correct_overlaps = [('C1', 'C11'), ('N1', 'N11')]
    for atomName1, atomName2 in correct_overlaps:
        assert suptop.contains_atomNamePair(atomName1, atomName2)


def test_esterSymmetry_rightStart():
    """
    A simple molecule chain with an ester.
    The ester allows for mapping (O1-O11, O2-O12) and (O1-O12, O2-O11)

     LIGAND 1        LIGAND 2
        C1              C11
        \                \
        N1              N11
        /\              / \
     O1    O2        O11   O12

    """
    # construct the LIGAND 1
    # ignore the third dimension
    # c1
    c1 = AtomNode(1, 'C1', 'TEST', 1, 0, 'C')
    c1.set_coords(np.array([1, 1, 0], dtype='float32'))
    n1 = AtomNode(2, 'N1', 'TEST', 1, 0, 'N')
    n1.set_coords(np.array([1, 2, 0], dtype='float32'))
    c1.bindTo(n1)

    # construct the LIGAND 2
    # ignore the third dimension
    # c1
    c11 = AtomNode(11, 'C11', 'TEST', 1, 0, 'C')
    c11.set_coords(np.array([1, 1, 0], dtype='float32'))
    n11 = AtomNode(11, 'N11', 'TEST', 1, 0, 'N')
    n11.set_coords(np.array([1, 2, 0], dtype='float32'))
    c11.bindTo(n11)

    # should overlap 2 atoms
    suptops = _overlay(c1, c11, atol=9999)
    assert len(suptops) == 1
    suptop = suptops[0]

    # the number of overlapped atoms is two
    assert len(suptop) == 2
    correct_overlaps = [('C1', 'C11'), ('N1', 'N11')]
    for atomName1, atomName2 in correct_overlaps:
        assert suptop.contains_atomNamePair(atomName1, atomName2)


def test_2atoms_2symmetries():
    """
    create a simple molecule chain with an ester
     LIGAND 1        LIGAND 2
        C1              C11
        \                \
        N1              N11
    In this there is only one solution.
    """
    # construct the LIGAND 1
    # ignore the third dimension
    # c1
    c1 = AtomNode(1, 'C1', 'TEST', 1, 0, 'C')
    c1.set_coords(np.array([1, 1, 0], dtype='float32'))
    n1 = AtomNode(2, 'N1', 'TEST', 1, 0, 'N')
    n1.set_coords(np.array([1, 2, 0], dtype='float32'))
    c1.bindTo(n1)

    # construct the LIGAND 2
    # ignore the third dimension
    # c1
    c11 = AtomNode(11, 'C11', 'TEST', 1, 0, 'C')
    c11.set_coords(np.array([1, 1, 0], dtype='float32'))
    n11 = AtomNode(11, 'N11', 'TEST', 1, 0, 'N')
    n11.set_coords(np.array([1, 2, 0], dtype='float32'))
    c11.bindTo(n11)

    # should return a list with an empty sup_top
    suptops = _overlay(c1, n11, atol=9999)
    assert suptops == []

    # should ovlap 2 atoms
    suptops = _overlay(c1, c11, atol=9999)
    assert len(suptops) == 1
    assert len(suptops[0]) == 2


def old_mcl1_l18l39_nocharges():
    """
    create a simple molecule chain with an ester
        C1              C11
        \                \
        C2              C12
        /                /
        C3               C13
        \                 \
        C4                C14
        /\                 /\
       O1  O2            O11 O12

    """
    #


    liglig_path = "asdf"
    # we are ignoring the charges by directly calling the superimposer
    suptops = _superimpose_topologies("asdf", "asdf", atol=9999)
    # in this case, there should be only one solution
    assert len(suptops) == 1

    suptop = suptops[0]
    assert len(suptop) == 43

    # check the core atoms of the superimposition for correctness
    # this ensures that the atoms which have no choice (ie they can only be superimposed in one way)
    # are superimposed that way
    core_test_pairs = [('C21', 'C43'), ('C15', 'C37'), ('O3', 'O6'), ('C9', 'C31'),
                  ('C1', 'C23'), ('N1', 'N2'), ('C6', 'C28'), ('C4', 'C26')]
    for atomName1, atomname2 in core_test_pairs:
        assert suptop.contains_atomNamePair(atomName1, atomname2)

    # resolve the multiple possbile matches
    avg_dst = suptop.correct_for_coordinates()

    # check if the mirrors were corrected
    corrected_symmetries = [('O2', 'O5'), ('O1', 'O4'), ('H7', 'H25'), ('H6', 'H24'),
                            ('H4', 'H22'), ('H5', 'H23'), ('H8', 'H26'), ('H9', 'H27')]
    for atomName1, atomname2 in corrected_symmetries:
        assert suptop.contains_atomNamePair(atomName1, atomname2)

    # refine against charges
    # ie remove the matches that change due to charge rather than spieces
    all_removed_pairs = suptop.refineAgainstCharges(atol=0.1)
    print(all_removed_pairs)
    removed_pairs = [('C5', 'C27'), ('C4', 'C26')]
    for atomName1, atomname2 in removed_pairs:
        assert not suptop.contains_atomNamePair(atomName1, atomname2)
