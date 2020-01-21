"""
Focuses on the _overlay function that actually traverses the molecule using the given starting points.
"""

from topology_superimposer import SuperimposedTopology, get_charges, \
    superimpose_topologies, _superimpose_topologies, assign_coords_from_pdb, \
    AtomNode, _overlay
import networkx as nx
from os import path
import numpy as np


def test_2diffAtoms_CN_wrongStart():
    """
    create a simple molecule chain with an ester
     LIGAND 1        LIGAND 2
        C1              C11
        \                \
        N1              N11
    In this there is only one solution.
    """
    # construct the LIGAND 1
    c1 = AtomNode(1, 'C1', 'TEST', 1, 0, 'C')
    c1.set_coords(np.array([1, 1, 0], dtype='float32'))
    n1 = AtomNode(2, 'N1', 'TEST', 1, 0, 'N')
    n1.set_coords(np.array([1, 2, 0], dtype='float32'))
    c1.bindTo(n1)

    # construct the LIGAND 2
    c11 = AtomNode(11, 'C11', 'TEST', 1, 0, 'C')
    c11.set_coords(np.array([1, 1, 0], dtype='float32'))
    n11 = AtomNode(11, 'N11', 'TEST', 1, 0, 'N')
    n11.set_coords(np.array([1, 2, 0], dtype='float32'))
    c11.bindTo(n11)

    # should return a list with an empty sup_top
    suptops = _overlay(c1, n11)
    # it should return an empty suptop
    assert len(suptops) == 1
    assert len(suptops[0]) == 0


def test_2diffAtoms_CN_rightStart():
    """
    Two different Atoms. Only one solution exists.

     LIGAND 1        LIGAND 2
        C1              C11
        \                \
        N1              N11
    """
    # construct the LIGAND 1
    c1 = AtomNode(1, 'C1', 'TEST', 1, 0, 'C')
    c1.set_coords(np.array([1, 1, 0], dtype='float32'))
    n1 = AtomNode(2, 'N1', 'TEST', 1, 0, 'N')
    n1.set_coords(np.array([1, 2, 0], dtype='float32'))
    c1.bindTo(n1)

    # construct the LIGAND 2
    c11 = AtomNode(11, 'C11', 'TEST', 1, 0, 'C')
    c11.set_coords(np.array([1, 1, 0], dtype='float32'))
    n11 = AtomNode(11, 'N11', 'TEST', 1, 0, 'N')
    n11.set_coords(np.array([1, 2, 0], dtype='float32'))
    c11.bindTo(n11)

    # should overlap 2 atoms
    suptops = _overlay(c1, c11)
    assert len(suptops) == 1
    suptop = suptops[0]

    # the number of overlapped atoms is two
    assert len(suptop) == 2
    correct_overlaps = [('C1', 'C11'), ('N1', 'N11')]
    for atomName1, atomName2 in correct_overlaps:
        assert suptop.contains_atomNamePair(atomName1, atomName2)

    # there is no other ways to traverse the molecule
    assert len(suptop.mirrors) == 0


def test_3diffAtoms_CNO_rightStart():
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
    c1 = AtomNode(1, 'C1', 'TEST', 1, 0, 'C')
    c1.set_coords(np.array([1, 1, 0], dtype='float32'))
    n1 = AtomNode(2, 'N1', 'TEST', 1, 0, 'N')
    n1.set_coords(np.array([1, 2, 0], dtype='float32'))
    c1.bindTo(n1)
    o1 = AtomNode(3, 'O1', 'TEST', 1, 0, 'O')
    o1.set_coords(np.array([1, 3, 0], dtype='float32'))
    o1.bindTo(n1)

    # construct the LIGAND 2
    c11 = AtomNode(11, 'C11', 'TEST', 1, 0, 'C')
    c11.set_coords(np.array([1, 1, 0], dtype='float32'))
    n11 = AtomNode(12, 'N11', 'TEST', 1, 0, 'N')
    n11.set_coords(np.array([1, 2, 0], dtype='float32'))
    c11.bindTo(n11)
    o11 = AtomNode(13, 'O11', 'TEST', 1, 0, 'O')
    o11.set_coords(np.array([1, 3, 0], dtype='float32'))
    o11.bindTo(n11)

    # should overlap 2 atoms
    suptops = _overlay(c1, c11)
    assert len(suptops) == 1
    suptop = suptops[0]

    # the number of overlapped atoms is two
    assert len(suptop) == 3
    correct_overlaps = [('C1', 'C11'), ('N1', 'N11'), ('O1', 'O11')]
    for atomName1, atomName2 in correct_overlaps:
        assert suptop.contains_atomNamePair(atomName1, atomName2)

    # no mirrors
    assert len(suptop.mirrors) == 0


def test_SimpleMultipleSolutions_rightStart():
    """
    A simple molecule chain with an ester.
    The ester allows for mapping (O1-O11, O2-O12) and (O1-O12, O2-O11)
    # fixme what if you start from the wrong O-O matching? should that be a case? how does
    comparies topologies would behave in that case?

     LIGAND 1        LIGAND 2
        C1              C11
        \                \
        N1              N11
        /\              / \
     O1    O2        O11   O12

    """
    # ignore the third coordinate dimension
    # construct the LIGAND 1
    c1 = AtomNode(1, 'C1', 'TEST', 1, 0, 'C')
    c1.set_coords(np.array([1, 1, 0], dtype='float32'))
    n1 = AtomNode(2, 'N1', 'TEST', 1, 0, 'N')
    n1.set_coords(np.array([1, 2, 0], dtype='float32'))
    c1.bindTo(n1)
    o1 = AtomNode(3, 'O1', 'TEST', 1, 0, 'O')
    o1.set_coords(np.array([1, 3, 0], dtype='float32'))
    o1.bindTo(n1)
    o2 = AtomNode(4, 'O2', 'TEST', 1, 0, 'O')
    o2.set_coords(np.array([2, 3, 0], dtype='float32'))
    o2.bindTo(n1)

    # construct the LIGAND 2
    c11 = AtomNode(11, 'C11', 'TEST', 1, 0, 'C')
    c11.set_coords(np.array([1, 1, 0], dtype='float32'))
    n11 = AtomNode(11, 'N11', 'TEST', 1, 0, 'N')
    n11.set_coords(np.array([1, 2, 0], dtype='float32'))
    c11.bindTo(n11)
    o11 = AtomNode(11, 'O11', 'TEST', 1, 0, 'O')
    o11.set_coords(np.array([1, 3, 0], dtype='float32'))
    o11.bindTo(n11)
    o12 = AtomNode(11, 'O12', 'TEST', 1, 0, 'O')
    o12.set_coords(np.array([2, 3, 0], dtype='float32'))
    o12.bindTo(n11)

    # should be two topologies
    suptops = _overlay(c1, c11)
    assert len(suptops) == 2

    # check if both representations were found
    # The ester allows for mapping (O1-O11, O2-O12) and (O1-O12, O2-O11)
    assert any(st.contains_atomNamePair('O1', 'O11') and st.contains_atomNamePair('O2', 'O12') for st in suptops)
    assert not all(st.contains_atomNamePair('O1', 'O11') and st.contains_atomNamePair('O2', 'O12') for st in suptops)
    assert any(st.contains_atomNamePair('O1', 'O12') and st.contains_atomNamePair('O2', 'O11') for st in suptops)
    assert not all(st.contains_atomNamePair('O1', 'O12') and st.contains_atomNamePair('O2', 'O11') for st in suptops)

    correct_overlaps = [('C1', 'C11'), ('N1', 'N11')]
    for st in suptops:
        for atomName1, atomName2 in correct_overlaps:
            assert st.contains_atomNamePair(atomName1, atomName2)

    # fixme - add a test case for the superimposer function that makes use of _overlay,
    # this is to resolve multiple solutions such as the one here


def test_2sameAtoms_2Cs_symmetry():
    """
    Two solutions with different starting points.

     LIGAND 1        LIGAND 2
        C1              C11
        \                \
        C2              C12
    """
    # construct the LIGAND 1
    c1 = AtomNode(1, 'C1', 'TEST', 1, 0, 'C')
    c1.set_coords(np.array([1, 1, 0], dtype='float32'))
    c2 = AtomNode(2, 'C2', 'TEST', 1, 0, 'C')
    c2.set_coords(np.array([1, 2, 0], dtype='float32'))
    c1.bindTo(c2)

    # construct the LIGAND 2
    c11 = AtomNode(11, 'C11', 'TEST', 1, 0, 'C')
    c11.set_coords(np.array([1, 1, 0], dtype='float32'))
    c12 = AtomNode(11, 'C12', 'TEST', 1, 0, 'C')
    c12.set_coords(np.array([1, 2, 0], dtype='float32'))
    c11.bindTo(c12)

    # should return a list with an empty sup_top
    suptops = _overlay(c1, c11)
    assert len(suptops) == 1
    assert suptops[0].contains_atomNamePair('C1', 'C11')
    assert suptops[0].contains_atomNamePair('C2', 'C12')

    suptops = _overlay(c1, c12)
    assert len(suptops) == 1
    assert suptops[0].contains_atomNamePair('C1', 'C12')
    assert suptops[0].contains_atomNamePair('C2', 'C11')


def test_mutation_separate_unique_match():
    """
    Two commponents separated by the mutation.

     LIGAND 1        LIGAND 2
        C1              C11
        |                |
        S1              O11
        |                |
        N1               N11
    """
    # construct the LIGAND 1
    c1 = AtomNode(1, 'C1', 'TEST', 1, 0, 'C')
    c1.set_coords(np.array([1, 1, 0], dtype='float32'))
    s1 = AtomNode(2, 'S2', 'TEST', 1, 0, 'S')
    s1.set_coords(np.array([2, 1, 0], dtype='float32'))
    c1.bindTo(s1)
    n1 = AtomNode(2, 'N1', 'TEST', 1, 0, 'N')
    n1.set_coords(np.array([3, 1, 0], dtype='float32'))
    n1.bindTo(s1)

    # construct the LIGAND 2
    c11 = AtomNode(11, 'C11', 'TEST', 1, 0, 'C')
    c11.set_coords(np.array([1, 1, 0], dtype='float32'))
    o11 = AtomNode(11, 'O11', 'TEST', 1, 0, 'O')
    o11.set_coords(np.array([2, 1, 0], dtype='float32'))
    c11.bindTo(o11)
    n11 = AtomNode(2, 'N11', 'TEST', 1, 0, 'N')
    n11.set_coords(np.array([3, 1, 0], dtype='float32'))
    n11.bindTo(o11)

    # should return a list with an empty sup_top
    suptops = _overlay(c1, c11)
    assert len(suptops) == 1
    assert suptops[0].contains_atomNamePair('C1', 'C11')

    suptops = _overlay(n1, n11)
    assert len(suptops) == 1
    assert suptops[0].contains_atomNamePair('N1', 'N11')


def test_mutation_separate_unique_match():
    """
    Two commponents separated by the mutation.

     LIGAND 1        LIGAND 2
        C1              C11
        |                |
        C2              C12
        |                |
        C3              O11
        |                |
        N1              N11
    """
    # construct the LIGAND 1
    c1 = AtomNode(1, 'C1', 'TEST', 1, 0, 'C')
    c1.set_coords(np.array([1, 1, 0], dtype='float32'))
    c2 = AtomNode(1, 'C2', 'TEST', 1, 0, 'C')
    c2.set_coords(np.array([2, 1, 0], dtype='float32'))
    c2.bindTo(c1)
    c3 = AtomNode(1, 'C3', 'TEST', 1, 0, 'C')
    c3.set_coords(np.array([2, 1, 0], dtype='float32'))
    c3.bindTo(c2)
    n1 = AtomNode(2, 'N1', 'TEST', 1, 0, 'N')
    n1.set_coords(np.array([3, 1, 0], dtype='float32'))
    n1.bindTo(c3)

    # construct the LIGAND 2
    c11 = AtomNode(11, 'C11', 'TEST', 1, 0, 'C')
    c11.set_coords(np.array([1, 1, 0], dtype='float32'))
    c12 = AtomNode(11, 'C12', 'TEST', 1, 0, 'C')
    c12.set_coords(np.array([2, 1, 0], dtype='float32'))
    c12.bindTo(c11)
    o11 = AtomNode(11, 'O11', 'TEST', 1, 0, 'O')
    o11.set_coords(np.array([3, 1, 0], dtype='float32'))
    o11.bindTo(c12)
    n11 = AtomNode(11, 'N11', 'TEST', 1, 0, 'N')
    n11.set_coords(np.array([3, 1, 0], dtype='float32'))
    n11.bindTo(o11)

    # should return a list with an empty sup_top
    suptops = _overlay(c1, c11)
    assert len(suptops) == 1
    assert suptops[0].contains_atomNamePair('C1', 'C11')

    suptops = _overlay(n1, n11)
    assert len(suptops) == 1
    assert suptops[0].contains_atomNamePair('N1', 'N11')


def test_3C_circle():
    """
    A circle should be detected.
    Many solutions (3 starting conditions, etc)

      LIGAND 1        LIGAND 2
         C1              C11
        /  \            /   \
      C2 - C3         C12 - C13
    """
    # construct the LIGAND 1
    c1 = AtomNode(1, 'C1', 'TEST', 1, 0, 'C')
    c1.set_coords(np.array([1, 1, 0], dtype='float32'))
    c2 = AtomNode(2, 'C2', 'TEST', 1, 0, 'C')
    c2.set_coords(np.array([1, 2, 0], dtype='float32'))
    c1.bindTo(c2)
    c3 = AtomNode(2, 'C3', 'TEST', 1, 0, 'C')
    c3.set_coords(np.array([2, 2, 0], dtype='float32'))
    c3.bindTo(c1)
    c3.bindTo(c2)

    # construct the LIGAND 2
    c11 = AtomNode(11, 'C11', 'TEST', 1, 0, 'C')
    c11.set_coords(np.array([1, 1, 0], dtype='float32'))
    c12 = AtomNode(11, 'C12', 'TEST', 1, 0, 'C')
    c12.set_coords(np.array([1, 2, 0], dtype='float32'))
    c11.bindTo(c12)
    c13 = AtomNode(11, 'C13', 'TEST', 1, 0, 'C')
    c13.set_coords(np.array([2, 2, 0], dtype='float32'))
    c13.bindTo(c11)
    c13.bindTo(c12)

    suptops = _overlay(c1, c11)
    # there are two solutions
    assert len(suptops) == 2
    assert any(st.contains_atomNamePair('C2', 'C12') and st.contains_atomNamePair('C3', 'C13') for st in suptops)
    assert any(st.contains_atomNamePair('C2', 'C13') and st.contains_atomNamePair('C3', 'C12') for st in suptops)
    # both solutions should have the same starting situation
    assert all(st.contains_atomNamePair('C1', 'C11') for st in suptops)
    # there should be one circle in each
    assert all(st.sameCircleNumber() for st in suptops)
    assert all(st.getCircleNumber() == (1, 1) for st in suptops)

    suptops = _overlay(c1, c12)
    assert len(suptops) == 2
    # there should be one circle in each
    assert all(st.sameCircleNumber() for st in suptops)
    assert all(st.getCircleNumber() == (1, 1) for st in suptops)

    suptops = _overlay(c1, c13)
    assert len(suptops) == 2
    # there should be one circle in each
    assert all(st.sameCircleNumber() for st in suptops)
    assert all(st.getCircleNumber() == (1, 1) for st in suptops)

    suptops = _overlay(c2, c11)
    assert len(suptops) == 2
    # there should be one circle in each
    assert all(st.sameCircleNumber() for st in suptops)
    assert all(st.getCircleNumber() == (1, 1) for st in suptops)

    suptops = _overlay(c2, c12)
    assert len(suptops) == 2
    # there should be one circle in each
    assert all(st.sameCircleNumber() for st in suptops)
    assert all(st.getCircleNumber() == (1, 1) for st in suptops)

    suptops = _overlay(c2, c13)
    assert len(suptops) == 2
    # there should be one circle in each
    assert all(st.sameCircleNumber() for st in suptops)
    assert all(st.getCircleNumber() == (1, 1) for st in suptops)

    suptops = _overlay(c3, c11)
    assert len(suptops) == 2
    # there should be one circle in each
    assert all(st.sameCircleNumber() for st in suptops)
    assert all(st.getCircleNumber() == (1, 1) for st in suptops)

    suptops = _overlay(c3, c12)
    assert len(suptops) == 2
    # there should be one circle in each
    assert all(st.sameCircleNumber() for st in suptops)
    assert all(st.getCircleNumber() == (1, 1) for st in suptops)

    suptops = _overlay(c3, c13)
    assert len(suptops) == 2
    # there should be one circle in each
    assert all(st.sameCircleNumber() for st in suptops)
    assert all(st.getCircleNumber() == (1, 1) for st in suptops)


def test_mcl1_l12l35():
    """
    Molecule inspired by Agastya's dataset (mcl1_l12l35).

    Ligand 1

         C1 - C2
         /      \
    Cl1-C3      C4
          \     /
          C5 - C6
          /     \
     C10-C7       N1
           \   /
             C8
             |
             C9


    Ligand 2
                 Cl11
                /
         C11 - C12
         /      \
        C13      C14
          \     /
          C15 - C16
          /     \
     C20-C17       N11
           \   /
             C18
             |
             C19
    """
    # construct LIGAND 1
    c1 = AtomNode(1, 'C1', 'TEST', 1, 0, 'C')
    c1.set_coords(np.array([1, 1, 0], dtype='float32'))
    c2 = AtomNode(2, 'C2', 'TEST', 1, 0, 'C')
    c2.set_coords(np.array([1, 2, 0], dtype='float32'))
    c1.bindTo(c2)
    c3 = AtomNode(2, 'C3', 'TEST', 1, 0, 'C')
    c3.set_coords(np.array([2, 2, 0], dtype='float32'))
    c3.bindTo(c1)
    cl1 = AtomNode(2, 'CL1', 'TEST', 1, 0, 'Cl')
    cl1.set_coords(np.array([2, 1, 0], dtype='float32'))
    cl1.bindTo(c3)
    c4 = AtomNode(2, 'C4', 'TEST', 1, 0, 'C')
    c4.set_coords(np.array([2, 3, 0], dtype='float32'))
    c4.bindTo(c2)
    c5 = AtomNode(2, 'C5', 'TEST', 1, 0, 'C')
    c5.set_coords(np.array([3, 1, 0], dtype='float32'))
    c5.bindTo(c3)
    c6 = AtomNode(2, 'C6', 'TEST', 1, 0, 'C')
    c6.set_coords(np.array([3, 2, 0], dtype='float32'))
    c6.bindTo(c5)
    c6.bindTo(c4)
    c7 = AtomNode(2, 'C7', 'TEST', 1, 0, 'C')
    c7.set_coords(np.array([4, 2, 0], dtype='float32'))
    c7.bindTo(c5)
    c10 = AtomNode(2, 'C10', 'TEST', 1, 0, 'C')
    c10.set_coords(np.array([4, 1, 0], dtype='float32'))
    c10.bindTo(c7)
    n1 = AtomNode(2, 'N1', 'TEST', 1, 0, 'N')
    n1.set_coords(np.array([4, 3, 0], dtype='float32'))
    n1.bindTo(c6)
    c8 = AtomNode(2, 'C8', 'TEST', 1, 0, 'C')
    c8.set_coords(np.array([5, 1, 0], dtype='float32'))
    c8.bindTo(c7)
    c8.bindTo(n1)
    c9 = AtomNode(2, 'C9', 'TEST', 1, 0, 'C')
    c9.set_coords(np.array([6, 1, 0], dtype='float32'))
    c9.bindTo(c8)

    # construct Ligand 2
    cl11 = AtomNode(2, 'Cl11', 'TEST', 1, 0, 'Cl')
    cl11.set_coords(np.array([1, 1, 0], dtype='float32'))
    c11 = AtomNode(1, 'C11', 'TEST', 1, 0, 'C')
    c11.set_coords(np.array([2, 1, 0], dtype='float32'))
    c12 = AtomNode(2, 'C12', 'TEST', 1, 0, 'C')
    c12.set_coords(np.array([2, 2, 0], dtype='float32'))
    c12.bindTo(c11)
    c12.bindTo(cl11)
    c13 = AtomNode(2, 'C13', 'TEST', 1, 0, 'C')
    c13.set_coords(np.array([3, 1, 0], dtype='float32'))
    c13.bindTo(c11)
    c14 = AtomNode(2, 'C14', 'TEST', 1, 0, 'C')
    c14.set_coords(np.array([3, 2, 0], dtype='float32'))
    c14.bindTo(c12)
    c15 = AtomNode(2, 'C15', 'TEST', 1, 0, 'C')
    c15.set_coords(np.array([4, 1, 0], dtype='float32'))
    c15.bindTo(c13)
    c16 = AtomNode(2, 'C16', 'TEST', 1, 0, 'C')
    c16.set_coords(np.array([4, 2, 0], dtype='float32'))
    c16.bindTo(c15)
    c16.bindTo(c14)
    c17 = AtomNode(2, 'C17', 'TEST', 1, 0, 'C')
    c17.set_coords(np.array([5, 2, 0], dtype='float32'))
    c17.bindTo(c15)
    c20 = AtomNode(2, 'C20', 'TEST', 1, 0, 'C')
    c20.set_coords(np.array([5, 1, 0], dtype='float32'))
    c20.bindTo(c17)
    n11 = AtomNode(2, 'N11', 'TEST', 1, 0, 'N')
    n11.set_coords(np.array([5, 3, 0], dtype='float32'))
    n11.bindTo(c16)
    c18 = AtomNode(2, 'C18', 'TEST', 1, 0, 'C')
    c18.set_coords(np.array([6, 1, 0], dtype='float32'))
    c18.bindTo(c17)
    c18.bindTo(n11)
    c19 = AtomNode(2, 'C19', 'TEST', 1, 0, 'C')
    c19.set_coords(np.array([7, 1, 0], dtype='float32'))
    c19.bindTo(c18)

    # the correct solution
    suptops = _overlay(c9, c19)
    assert len(suptops) == 1
    assert len(suptops[0]) == 11

    """
    fixme
    This is a rare case around which we'll have to find a work around.
    Basically, the best solution that obeyes our traversal allows for a superimposition
    that should not be allowed. So additional additional way has to be checked to discredit it 
    """
    suptops = _overlay(c5, c14)
    assert len(suptops) == 1
    assert len(suptops[0]) != 12
