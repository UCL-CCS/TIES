"""
Focuses on the _overlay function that actually traverses the molecule using the given starting points.

TODO:
- add more complex test case with multiple solutions (2 levels)  similar to test_SimpleMultipleSolutions_rightStart
- add a test case that shows that molecules cannot be fully superimposed when it is broken in half
- add more test cases that take into account the "differences"
"""

import pytest
from ties.topology_superimposer import Atom, _overlay, SuperimposedTopology


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
    c1 = Atom(name='C1', atom_type='C')
    c1.position = (1, 1, 0)
    n1 = Atom(name='N1', atom_type='N')
    n1.position = (1, 2, 0)
    c1.bind_to(n1, 'bondType1')
    left_atoms = [c1, n1]

    # construct the LIGAND 2
    c11 = Atom(name='C11', atom_type='C')
    c11.position = (1, 1, 0)
    n11 = Atom(name='N11', atom_type='N')
    n11.position = (1, 2, 0)
    c11.bind_to(n11, 'bondType1')
    right_atoms = [c11, n11]

    # should return a list with an empty sup_top
    #suptop = SuperimposedTopology(top1_nodes, top2_nodes, mda1_nodes, mda2_nodes)
    suptops = _overlay(c1, n11, parent_n1=None, parent_n2=None, bond_types=(None, None),
                       suptop=SuperimposedTopology(left_atoms, right_atoms))
    # it should return an empty suptop
    assert suptops == []


def test_2diffAtoms_CN_rightStart():
    """
    Two different Atoms. Only one solution exists.

     LIGAND 1        LIGAND 2
        C1              C11
        \                \
        N1              N11
    """
    # construct the LIGAND 1
    c1 = Atom(name='C1', atom_type='C')
    c1.position = (1, 1, 0)
    n1 = Atom(name='N1', atom_type='N')
    n1.position = (1, 2, 0)
    c1.bind_to(n1, 'bondType1')
    left_atoms=[c1, n1]

    # construct the LIGAND 2
    c11 = Atom(name='C11', atom_type='C')
    c11.position = (1, 1, 0)
    n11 = Atom(name='N11', atom_type='N')
    n11.position = (1, 2, 0)
    c11.bind_to(n11, 'bondType1')
    right_atoms = [c11, n11]

    #
    suptops = _overlay(c1, c11, parent_n1=None, parent_n2=None, bond_types=(None, None),
                      suptop=SuperimposedTopology(left_atoms, right_atoms))
    assert len(suptops) == 1
    suptop = suptops.pop()

    # should overlap 2 atoms
    assert len(suptop) == 2
    correct_overlaps = [('C1', 'C11'), ('N1', 'N11')]
    for atomName1, atomName2 in correct_overlaps:
        assert suptop.contains_atom_name_pair(atomName1, atomName2)

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
    c1 = Atom(name='C1', atom_type='C')
    c1.position = (1, 1, 0)
    n1 = Atom(name='N1', atom_type='N')
    n1.position = (1, 2, 0)
    c1.bind_to(n1, 'bondType1')
    o1 = Atom(name='O1', atom_type='O')
    o1.position = (1, 3, 0)
    o1.bind_to(n1, 'bondType1')
    left_atoms = [c1, n1, o1]

    # construct the LIGAND 2
    c11 = Atom(name='C11', atom_type='C')
    c11.position = (1, 1, 0)
    n11 = Atom(name='N11', atom_type='N')
    n11.position = (1, 2, 0)
    c11.bind_to(n11, 'bondType1')
    o11 = Atom(name='O11', atom_type='O')
    o11.position = (1, 3, 0)
    o11.bind_to(n11, 'bondType1')
    right_atoms = [c11, n11, o11]

    # should overlap 2 atoms
    suptops = _overlay(c1, c11, parent_n1=None, parent_n2=None, bond_types=(None, None),
                      suptop=SuperimposedTopology(left_atoms, right_atoms))
    assert len(suptops) == 1
    suptop = suptops.pop()

    # the number of overlapped atoms is two
    assert len(suptop) == 3
    correct_overlaps = [('C1', 'C11'), ('N1', 'N11'), ('O1', 'O11')]
    for atomName1, atomName2 in correct_overlaps:
        assert suptop.contains_atom_name_pair(atomName1, atomName2)

    # no mirrors
    assert len(suptop.mirrors) == 0


def test_SimpleMultipleSolutions_rightStart():
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
    # ignore the third coordinate dimension
    # construct the LIGAND 1
    c1 = Atom(name='C1', atom_type='C')
    c1.position = (1, 1, 0)
    n1 = Atom(name='N1', atom_type='N')
    n1.position = (1, 2, 0)
    c1.bind_to(n1, 'bondType1')
    o1 = Atom(name='O1', atom_type='O')
    o1.position = (1, 3, 0)
    o1.bind_to(n1, 'bondType1')
    o2 = Atom(name='O2', atom_type='O')
    o2.position = (2, 3, 0)
    o2.bind_to(n1, 'bondType1')
    left_atoms = [c1, n1, o1, o2]

    # construct the LIGAND 2
    c11 = Atom(name='C11', atom_type='C')
    c11.position = (1, 1, 0)
    n11 = Atom(name='N11', atom_type='N')
    n11.position = (1, 2, 0)
    c11.bind_to(n11, 'bondType1')
    o11 = Atom(name='O11', atom_type='O')
    o11.position = (1, 3, 0)
    o11.bind_to(n11, 'bondType1')
    o12 = Atom(name='O12', atom_type='O')
    o12.position = (2, 3, 0)
    o12.bind_to(n11, 'bondType1')
    right_atoms = [c11, n11, o11, o12]

    # the good solution is (O1-O11) and (O2-O12)

    # should generate one topology which has one "symmetry"
    suptops = _overlay(c1, c11, parent_n1=None, parent_n2=None, bond_types=(None, None),
                      suptop=SuperimposedTopology(left_atoms, right_atoms))
    assert len(suptops) == 1
    suptop = suptops.pop()
    assert len(suptop) == 4

    assert len(suptop.mirrors) == 1
    worse_st = suptop.mirrors[0]

    # check if both representations were found
    # The ester allows for mapping (O1-O11, O2-O12) and (O1-O12, O2-O11)
    assert suptop.contains_atom_name_pair('O1', 'O11') and suptop.contains_atom_name_pair('O2', 'O12')
    assert worse_st.contains_atom_name_pair('O1', 'O12') and worse_st.contains_atom_name_pair('O2', 'O11')

    correct_overlaps = [('C1', 'C11'), ('N1', 'N11')]
    for atomName1, atomName2 in correct_overlaps:
        assert suptop.contains_atom_name_pair(atomName1, atomName2)


def test_One2Many():
    """

     LIGAND 1        LIGAND 2
        N1              N11
        /              /   \
     O1              O11   O12

    """
    # ignore the third coordinate dimension
    # construct the LIGAND 1
    n1 = Atom(name='N1', atom_type='N')
    n1.position = (1, 1, 0)
    o1 = Atom(name='O1', atom_type='O')
    o1.position = (2, 1, 0)
    o1.bind_to(n1, 'bondType1')
    left_atoms = [n1, o1]

    # construct the LIGAND 2
    n11 = Atom(name='N11', atom_type='N')
    n11.position = (1, 1, 0)
    o11 = Atom(name='O11', atom_type='O')
    o11.position = (2, 1, 0)
    o11.bind_to(n11, 'bondType1')
    o12 = Atom(name='O12', atom_type='O')
    o12.position = (2, 2, 0)
    o12.bind_to(n11, 'bondType1')
    right_atoms = [n11, o11, o12]

    # the good solution is (O1-O11)

    # should generate one topology which has one "symmetry"
    suptops = _overlay(n1, n11, parent_n1=None, parent_n2=None, bond_types=(None, None),
                      suptop=SuperimposedTopology(left_atoms, right_atoms))
    assert len(suptops) == 2

    # extract the better suptop
    suptop = suptops[0]
    assert len(suptop) == 2

    # this is deduced by RMSD
    worse_suptop = suptops[1]

    # check if both representations were found
    # The ester allows for mapping (O1-O11, O2-O12) and (O1-O12, O2-O11)
    assert suptop.contains_atom_name_pair('O1', 'O11')
    assert worse_suptop.contains_atom_name_pair('O1', 'O12')

    assert suptop.contains_atom_name_pair('N1', 'N11')


def test_Many2One_part1():
    """

     LIGAND 1        LIGAND 2
        N1              N11
        / \              /
     O1    O2          O11

    """
    # ignore the third coordinate dimension
    # construct the LIGAND 1
    n1 = Atom(name='N1', atom_type='N')
    n1.position = (1, 1, 0)
    o1 = Atom(name='O1', atom_type='O')
    o1.position = (2, 1, 0)
    o1.bind_to(n1, 'bondType1')
    o2 = Atom(name='O2', atom_type='O')
    o2.position = (2, 2, 0)
    o2.bind_to(n1, 'bondType1')
    left_atoms = [n1, o1, o2]

    # construct the LIGAND 2
    n11 = Atom(name='N11', atom_type='N')
    n11.position = (1, 1, 0)
    o11 = Atom(name='O11', atom_type='O')
    o11.position = (2, 1, 0)
    o11.bind_to(n11, 'bondType1')
    right_atoms = [n11, o11]

    # the good solution is (O1-O11)

    # should generate one topology which has one "symmetry"
    suptops = _overlay(n1, n11, parent_n1=None, parent_n2=None, bond_types=(None, None),
                      suptop=SuperimposedTopology(left_atoms, right_atoms))
    assert len(suptops) == 2
    suptop = suptops[0]
    assert len(suptop) == 2

    worse_suptop = suptops[1]

    # check if both representations were found
    # The ester allows for mapping (O1-O11, O2-O12) and (O1-O12, O2-O11)
    assert suptop.contains_atom_name_pair('O1', 'O11')
    assert suptop.contains_atom_name_pair('N1', 'N11')
    assert worse_suptop.contains_atom_name_pair('O2', 'O11')


def test_Many2One_part2():
    """

     LIGAND 1        LIGAND 2
        N1              N11
        / \              \
     O1    O2             O12

    """
    # ignore the third coordinate dimension
    # construct the LIGAND 1
    n1 = Atom(name='N1', atom_type='N')
    n1.position = (1, 1, 0)
    o1 = Atom(name='O1', atom_type='O')
    o1.position = (2, 1, 0)
    o1.bind_to(n1, 'bondType1')
    o2 = Atom(name='O2', atom_type='O')
    o2.position = (2, 2, 0)
    o2.bind_to(n1, 'bondType1')
    left_atoms = [n1, o1, o2]

    # construct the LIGAND 2
    n11 = Atom(name='N11', atom_type='N')
    n11.position = (1, 1, 0)
    o12 = Atom(name='O12', atom_type='O')
    o12.position = (2, 2, 0)
    o12.bind_to(n11, 'bondType1')
    right_atoms = [n11, o12]

    # the good solution is (O1-O11)

    # should generate one topology which has one "symmetry"
    suptops = _overlay(n1, n11, parent_n1=None, parent_n2=None, bond_types=(None, None),
                      suptop=SuperimposedTopology(left_atoms, right_atoms))
    assert len(suptops) == 2
    suptop = suptops[0]
    assert len(suptop) == 2

    worse_suptop = suptops[1]

    # check if both representations were found
    # The ester allows for mapping (O1-O11, O2-O12) and (O1-O12, O2-O11)
    assert suptop.contains_atom_name_pair('O2', 'O12')
    assert suptop.contains_atom_name_pair('N1', 'N11')
    assert worse_suptop.contains_atom_name_pair('O1', 'O12')


def test_MultipleSolutions2Levels_rightStart():
    """
    A test with many different solutions.
    The ester allows for mapping (O1-O11, O2-O12) and (O1-O12, O2-O11)
    Then for each of the mapping there are two ways to map Nitrogens N
    e.g. for (O1-O11) there is (N1-N11)(N2-N12) or (N1-N12)(N2-N11)
    However, since we are eagerly accepting early matches, not all possible
    combinations will be explored.
    e.g. C1-O1- on the return will match (N1-N11)(N2-N12) correctly.

    So we have:

         LIGAND 1        LIGAND 2
            C1              C11
            /\              / \
         O1    O2        O11   O12
        / \    / \      /  \   /  \
       N1 N2  N3 N4   N11 N12 N13 N14
    """
    # ignore the third coordinate dimension
    # construct the LIGAND 1
    c1 = Atom(name='C1', atom_type='C')
    c1.position = (1, 1, 0)
    o1 = Atom(name='O1', atom_type='O')
    o1.position = (1, 3, 0)
    o1.bind_to(c1, 'bondType1')
    o2 = Atom(name='O2', atom_type='O')
    o2.position = (2, 3, 0)
    o2.bind_to(c1, 'bondType1')
    n1 = Atom(name='N1', atom_type='N')
    n1.position = (3, 1, 0)
    n1.bind_to(o1, 'bondType1')
    n2 = Atom(name='N2', atom_type='N')
    n2.position = (3, 2, 0)
    n2.bind_to(o1, 'bondType1')
    n3 = Atom(name='N3', atom_type='N')
    n3.position = (3, 3, 0)
    n3.bind_to(o2, 'bondType1')
    n4 = Atom(name='N4', atom_type='N')
    n4.position = (3, 4, 0)
    n4.bind_to(o2, 'bondType1')
    left_atoms = [c1, o1, o2, n1, n2, n3, n4]

    # construct the LIGAND 2
    c11 = Atom(name='C11', atom_type='C')
    c11.position = (1, 1, 0)
    o11 = Atom(name='O11', atom_type='O')
    o11.position = (1, 3, 0)
    o11.bind_to(c11, 'bondType1')
    o12 = Atom(name='O12', atom_type='O')
    o12.position = (2, 3, 0)
    o12.bind_to(c11, 'bondType1')
    n11 = Atom(name='N11', atom_type='N')
    n11.position = (3, 1, 0)
    n11.bind_to(o11, 'bondType1')
    n12 = Atom(name='N12', atom_type='N')
    n12.position = (3, 2, 0)
    n12.bind_to(o11, 'bondType1')
    n13 = Atom(name='N13', atom_type='N')
    n13.position = (3, 3, 0)
    n13.bind_to(o12, 'bondType1')
    n14 = Atom(name='N14', atom_type='N')
    n14.position = (3, 4, 0)
    n14.bind_to(o12, 'bondType1')
    right_atoms = [c11, o11, o12, n11, n12, n13, n14]

    # should generate one topology which has one "symmetry"
    suptops = _overlay(c1, c11, parent_n1=None, parent_n2=None, bond_types=(None, None),
                      suptop=SuperimposedTopology(left_atoms, right_atoms))
    assert len(suptops) == 1
    suptop = suptops.pop()

    # check if the main solution is correct
    correct_overlaps = [('C1', 'C11'),
                        ('O1', 'O11'), ('O2', 'O12'),
                        ('N1', 'N11'), ('N2', 'N12'), ('N3', 'N13'), ('N4', 'N14')]
    for atomName1, atomName2 in correct_overlaps:
        assert suptop.contains_atom_name_pair(atomName1, atomName2)

    # fixme - check the mirrors and esnure that they are properly stored.


def test_2sameAtoms_2Cs_symmetry():
    """
    Two solutions with different starting points.

     LIGAND 1        LIGAND 2
        C1              C11
        \                \
        C2              C12
    """
    # construct the LIGAND 1
    c1 = Atom(name='C1', atom_type='C')
    c1.position = (1, 1, 0)
    c2 = Atom(name='C2', atom_type='C')
    c2.position = (1, 2, 0)
    c1.bind_to(c2, 'bondType1')
    left_atoms = [c1, c2]

    # construct the LIGAND 2
    c11 = Atom(name='C11', atom_type='C')
    c11.position = (1, 1, 0)
    c12 = Atom(name='C12', atom_type='C')
    c12.position = (1, 2, 0)
    c11.bind_to(c12, 'bondType1')
    right_atoms = [c11, c12]

    # should return a list with an empty sup_top
    suptops = _overlay(c1, c11, parent_n1=None, parent_n2=None, bond_types=(None, None),
                      suptop=SuperimposedTopology(left_atoms, right_atoms))
    suptop = suptops.pop()
    assert len(suptop) == 2
    assert suptop.contains_atom_name_pair('C1', 'C11')
    assert suptop.contains_atom_name_pair('C2', 'C12')

    suptops = _overlay(c1, c12, parent_n1=None, parent_n2=None, bond_types=(None, None),
                      suptop=SuperimposedTopology(left_atoms, right_atoms))
    suptop = suptops.pop()
    assert len(suptop) == 2
    assert suptop.contains_atom_name_pair('C1', 'C12')
    assert suptop.contains_atom_name_pair('C2', 'C11')


def test_methyl():
    """
    Two solutions with different starting points.

     LIGAND 1        LIGAND 2
        C1              C11
       / | \            / | \
     H1 H2 H3        H11 H12 H13
    """
    # construct the LIGAND 1
    c1 = Atom(name='C1', atom_type='C')
    c1.position = (1, 1, 0)
    h1 = Atom(name='H1', atom_type='H')
    h1.position = (2, 1, 0)
    c1.bind_to(h1, 'bondType1')
    h2 = Atom(name='H2', atom_type='H')
    h2.position = (2, 2, 0)
    c1.bind_to(h2, 'bondType1')
    h3 = Atom(name='H3', atom_type='H')
    h3.position = (2, 3, 0)
    c1.bind_to(h3, 'bondType1')
    left_atoms = [c1, h1, h2, h3]

    # construct the LIGAND 2
    c11 = Atom(name='C11', atom_type='C')
    c11.position = (1, 1, 0)
    h11 = Atom(name='H11', atom_type='H')
    h11.position = (2, 1, 0)
    c11.bind_to(h11, 'bondType1')
    h12 = Atom(name='H12', atom_type='H')
    h12.position = (2, 2, 0)
    c11.bind_to(h12, 'bondType1')
    h13 = Atom(name='H13', atom_type='H')
    h13.position = (2, 3, 0)
    c11.bind_to(h13, 'bondType1')
    right_atoms = [c11, h11, h12, h13]

    # should return a list with an empty sup_top
    suptops = _overlay(c1, c11, parent_n1=None, parent_n2=None, bond_types=(None, None),
                      suptop=SuperimposedTopology(left_atoms, right_atoms))
    suptop = suptops.pop()
    assert suptop.contains_atom_name_pair('C1', 'C11')
    assert suptop.contains_atom_name_pair('H1', 'H11')
    assert suptop.contains_atom_name_pair('H2', 'H12')
    assert suptop.contains_atom_name_pair('H3', 'H13')


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
    c1 = Atom(name='C1', atom_type='C')
    c1.position = (1, 1, 0)
    s1 = Atom(name='S2', atom_type='S')
    s1.position = (2, 1, 0)
    c1.bind_to(s1, 'bondType1')
    n1 = Atom(name='N1', atom_type='N')
    n1.position = (3, 1, 0)
    n1.bind_to(s1, 'bondType1')

    # construct the LIGAND 2
    c11 = Atom(name='C11', atom_type='C')
    c11.position = (1, 1, 0)
    o11 = Atom(name='O11', atom_type='O')
    o11.position = (2, 1, 0)
    c11.bind_to(o11, 'bondType1')
    n11 = Atom(name='N11', atom_type='N')
    n11.position = (3, 1, 0)
    n11.bind_to(o11, 'bondType1')

    # should return a list with an empty sup_top
    suptops = _overlay(c1, c11, parent_n1=None, parent_n2=None, bond_types=(None, None))
    assert len(suptops) == 1
    assert suptops[0].contains_atom_name_pair('C1', 'C11')

    suptops = _overlay(n1, n11, parent_n1=None, parent_n2=None, bond_types=(None, None))
    assert len(suptops) == 1
    assert suptops[0].contains_atom_name_pair('N1', 'N11')


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
    c1 = Atom(name='C1', atom_type='C')
    c1.position = (1, 1, 0)
    c2 = Atom(name='C2', atom_type='C')
    c2.position = (2, 1, 0)
    c2.bind_to(c1, 'bondType1')
    c3 = Atom(name='C3', atom_type='C')
    c3.position = (2, 1, 0)
    c3.bind_to(c2, 'bondType1')
    n1 = Atom(name='N1', atom_type='N')
    n1.position = (3, 1, 0)
    n1.bind_to(c3, 'bondType1')
    left_atoms = [c1, c2, c3, n1]

    # construct the LIGAND 2
    c11 = Atom(name='C11', atom_type='C')
    c11.position = (1, 1, 0)
    c12 = Atom(name='C12', atom_type='C')
    c12.position = (2, 1, 0)
    c12.bind_to(c11, 'bondType1')
    o11 = Atom(name='O11', atom_type='O')
    o11.position = (3, 1, 0)
    o11.bind_to(c12, 'bondType1')
    n11 = Atom(name='N11', atom_type='N')
    n11.position = (3, 1, 0)
    n11.bind_to(o11, 'bondType1')
    right_atoms = [c11, c12, o11, n11]

    # should return a list with an empty sup_top
    suptops = _overlay(c1, c11, parent_n1=None, parent_n2=None, bond_types=(None, None),
                      suptop=SuperimposedTopology(left_atoms, right_atoms))
    suptop = suptops.pop()
    assert suptop.contains_atom_name_pair('C1', 'C11')

    suptops = _overlay(n1, n11, parent_n1=None, parent_n2=None, bond_types=(None, None),
                      suptop=SuperimposedTopology(left_atoms, right_atoms))
    suptop = suptops.pop()
    assert suptop.contains_atom_name_pair('N1', 'N11')


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
    c1 = Atom(name='C1', atom_type='C')
    c1.position = (1, 1, 0)
    c2 = Atom(name='C2', atom_type='C')
    c2.position = (1, 2, 0)
    c1.bind_to(c2, 'bondType1')
    c3 = Atom(name='C3', atom_type='C')
    c3.position = (2, 2, 0)
    c3.bind_to(c1, 'bondType1')
    c3.bind_to(c2, 'bondType1')
    left_atoms = [c1, c2, c3]

    # construct the LIGAND 2
    c11 = Atom(name='C11', atom_type='C')
    c11.position = (1, 1, 0)
    c12 = Atom(name='C12', atom_type='C')
    c12.position = (1, 2, 0)
    c11.bind_to(c12, 'bondType1')
    c13 = Atom(name='C13', atom_type='C')
    c13.position = (2, 2, 0)
    c13.bind_to(c11, 'bondType1')
    c13.bind_to(c12, 'bondType1')
    right_atoms = [c11, c12, c13]

    suptops = _overlay(c1, c11, parent_n1=None, parent_n2=None, bond_types=(None, None),
                      suptop=SuperimposedTopology(left_atoms, right_atoms))
    # there is one main component - that is the good solution
    suptop = suptops.pop()
    # there is one symmetrical way to traverse it
    assert len(suptop.mirrors) == 1
    wrong_st = suptop.mirrors[0]

    assert suptop.contains_atom_name_pair('C2', 'C12') and suptop.contains_atom_name_pair('C3', 'C13')
    assert wrong_st.contains_atom_name_pair('C2', 'C13') and wrong_st.contains_atom_name_pair('C3', 'C12')
    # both solutions should have the same starting situation
    assert all(st.contains_atom_name_pair('C1', 'C11') for st in [suptop, wrong_st])
    # there should be one circle in each
    assert all(st.same_circle_number() for st in [suptop, wrong_st])
    assert all(st.get_circle_number() == (1, 1) for st in [suptop, wrong_st])

    suptops = _overlay(c1, c12, parent_n1=None, parent_n2=None, bond_types=(None, None),
                      suptop=SuperimposedTopology(left_atoms, right_atoms))
    suptop = suptops.pop()
    # there should be one circle in each
    assert suptop.same_circle_number()
    assert suptop.get_circle_number() == (1, 1)

    suptops = _overlay(c1, c13, parent_n1=None, parent_n2=None, bond_types=(None, None),
                      suptop=SuperimposedTopology(left_atoms, right_atoms))
    suptop = suptops.pop()
    # there should be one circle in each
    assert suptop.same_circle_number()
    assert suptop.get_circle_number() == (1, 1)

    suptops = _overlay(c2, c11, parent_n1=None, parent_n2=None, bond_types=(None, None),
                      suptop=SuperimposedTopology(left_atoms, right_atoms))
    suptop = suptops.pop()
    # there should be one circle in each
    assert suptop.same_circle_number()
    assert suptop.get_circle_number() == (1, 1)

    suptops = _overlay(c2, c12, parent_n1=None, parent_n2=None, bond_types=(None, None),
                      suptop=SuperimposedTopology(left_atoms, right_atoms))
    suptop = suptops.pop()
    # there should be one circle in each
    assert suptop.same_circle_number()
    assert suptop.get_circle_number() == (1, 1)

    suptops = _overlay(c2, c13, parent_n1=None, parent_n2=None, bond_types=(None, None),
                      suptop=SuperimposedTopology(left_atoms, right_atoms))
    suptop = suptops.pop()
    # there should be one circle in each
    assert suptop.same_circle_number()
    assert suptop.get_circle_number() == (1, 1)

    suptops = _overlay(c3, c11, parent_n1=None, parent_n2=None, bond_types=(None, None),
                      suptop=SuperimposedTopology(left_atoms, right_atoms))
    suptop = suptops.pop()
    # there should be one circle in each
    assert suptop.same_circle_number()
    assert suptop.get_circle_number() == (1, 1)

    suptops = _overlay(c3, c12, parent_n1=None, parent_n2=None, bond_types=(None, None),
                      suptop=SuperimposedTopology(left_atoms, right_atoms))
    suptop = suptops.pop()
    # there should be one circle in each
    assert suptop.same_circle_number()
    assert suptop.get_circle_number() == (1, 1)

    suptops = _overlay(c3, c13, parent_n1=None, parent_n2=None, bond_types=(None, None),
                      suptop=SuperimposedTopology(left_atoms, right_atoms))
    suptop = suptops.pop()
    # there should be one circle in each
    assert suptop.same_circle_number()
    assert suptop.get_circle_number() == (1, 1)


def test_3C_circle_bonds():
    """
    Check if each pair is linked to two other pairs

      LIGAND 1        LIGAND 2
         C1              C11
        /  \            /   \
      C2 - C3         C12 - C13
    """
    # construct the LIGAND 1
    c1 = Atom(name='C1', atom_type='C')
    c1.position = (1, 1, 0)
    c2 = Atom(name='C2', atom_type='C')
    c2.position = (1, 2, 0)
    c1.bind_to(c2, 'bondType1')
    c3 = Atom(name='C3', atom_type='C')
    c3.position = (2, 2, 0)
    c3.bind_to(c1, 'bondType1')
    c3.bind_to(c2, 'bondType1')
    left_atoms = [c1, c2, c3]

    # construct the LIGAND 2
    c11 = Atom(name='C11', atom_type='C')
    c11.position = (1, 1, 0)
    c12 = Atom(name='C12', atom_type='C')
    c12.position = (1, 2, 0)
    c11.bind_to(c12, 'bondType1')
    c13 = Atom(name='C13', atom_type='C')
    c13.position = (2, 2, 0)
    c13.bind_to(c11, 'bondType1')
    c13.bind_to(c12, 'bondType1')
    right_atoms = [c11, c12, c13]

    suptop = _overlay(c1, c11, parent_n1=None, parent_n2=None, bond_types=(None, None),
                      suptop=SuperimposedTopology(left_atoms, right_atoms)).pop()

    # check that each pair is linked to two other pairs
    for pair, linked_pairs in suptop.matched_pairs_bonds.items():
        assert len(linked_pairs) == 2


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
    c1 = Atom(name='C1', atom_type='C')
    c1.position = (1, 1, 0)
    c2 = Atom(name='C2', atom_type='C')
    c2.position = (1, 2, 0)
    c1.bind_to(c2, 'bondType1')
    c3 = Atom(name='C3', atom_type='C')
    c3.position = (2, 2, 0)
    c3.bind_to(c1, 'bondType1')
    cl1 = Atom(name='CL1', atom_type='Cl')
    cl1.position = (2, 1, 0)
    cl1.bind_to(c3, 'bondType1')
    c4 = Atom(name='C4', atom_type='C')
    c4.position = (2, 3, 0)
    c4.bind_to(c2, 'bondType1')
    c5 = Atom(name='C5', atom_type='C')
    c5.position = (3, 1, 0)
    c5.bind_to(c3, 'bondType1')
    c6 = Atom(name='C6', atom_type='C')
    c6.position = (3, 2, 0)
    c6.bind_to(c5, 'bondType1')
    c6.bind_to(c4, 'bondType1')
    c7 = Atom(name='C7', atom_type='C')
    c7.position = (4, 2, 0)
    c7.bind_to(c5, 'bondType1')
    c10 = Atom(name='C10', atom_type='C')
    c10.position = (4, 1, 0)
    c10.bind_to(c7, 'bondType1')
    n1 = Atom(name='N1', atom_type='N')
    n1.position = (4, 3, 0)
    n1.bind_to(c6, 'bondType1')
    c8 = Atom(name='C8', atom_type='C')
    c8.position = (5, 1, 0)
    c8.bind_to(c7, 'bondType1')
    c8.bind_to(n1, 'bondType1')
    c9 = Atom(name='C9', atom_type='C')
    c9.position = (6, 1, 0)
    c9.bind_to(c8, 'bondType1')
    left_atoms = [c1, c2, c3, cl1, c4, c5, c6, c7, c8, c9, c10, n1]

    # construct Ligand 2
    cl11 = Atom(name='Cl11', atom_type='Cl')
    cl11.position = (1, 1, 0)
    c11 = Atom(name='C11', atom_type='C')
    c11.position = (2, 1, 0)
    c12 = Atom(name='C12', atom_type='C')
    c12.position = (2, 2, 0)
    c12.bind_to(c11, 'bondType1')
    c12.bind_to(cl11, 'bondType1')
    c13 = Atom(name='C13', atom_type='C')
    c13.position = (3, 1, 0)
    c13.bind_to(c11, 'bondType1')
    c14 = Atom(name='C14', atom_type='C')
    c14.position = (3, 2, 0)
    c14.bind_to(c12, 'bondType1')
    c15 = Atom(name='C15', atom_type='C')
    c15.position = (4, 1, 0)
    c15.bind_to(c13, 'bondType1')
    c16 = Atom(name='C16', atom_type='C')
    c16.position = (4, 2, 0)
    c16.bind_to(c15, 'bondType1')
    c16.bind_to(c14, 'bondType1')
    c17 = Atom(name='C17', atom_type='C')
    c17.position = (5, 2, 0)
    c17.bind_to(c15, 'bondType1')
    c20 = Atom(name='C20', atom_type='C')
    c20.position = (5, 1, 0)
    c20.bind_to(c17, 'bondType1')
    n11 = Atom(name='N11', atom_type='N')
    n11.position = (5, 3, 0)
    n11.bind_to(c16, 'bondType1')
    c18 = Atom(name='C18', atom_type='C')
    c18.position = (6, 1, 0)
    c18.bind_to(c17, 'bondType1')
    c18.bind_to(n11, 'bondType1')
    c19 = Atom(name='C19', atom_type='C')
    c19.position = (7, 1, 0)
    c19.bind_to(c18, 'bondType1')
    right_atoms = [cl11, c11, c12, c13, c14, c15, c16, c17, c18, c19, c20, n11]

    # the correct solution
    suptops = _overlay(c9, c19, parent_n1=None, parent_n2=None, bond_types=(None, None),
                      suptop=SuperimposedTopology(left_atoms, right_atoms))
    suptop = suptops.pop(0)
    assert len(suptop) == 11

    """
    This is a rare case which uses the rules of "consistent circles".
    A naive traversal of this is possible that reaches 12 nodes. However, 
    that means that there is atom match (L1-R1) such that L1 creates a 
    cricle in its own topology, and R1 does not. 
    """
    suptops = _overlay(c5, c14, parent_n1=None, parent_n2=None, bond_types=(None, None),
                      suptop=SuperimposedTopology(left_atoms, right_atoms))
    suptop = suptops.pop(0)
    assert len(suptop) != 12


def test_mcl1_l12l35_bonds():
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
    c1 = Atom(name='C1', atom_type='C')
    c1.position = (1, 1, 0)
    c2 = Atom(name='C2', atom_type='C')
    c2.position = (1, 2, 0)
    c1.bind_to(c2, 'bondType1')
    c3 = Atom(name='C3', atom_type='C')
    c3.position = (2, 2, 0)
    c3.bind_to(c1, 'bondType1')
    cl1 = Atom(name='CL1', atom_type='Cl')
    cl1.position = (2, 1, 0)
    cl1.bind_to(c3, 'bondType1')
    c4 = Atom(name='C4', atom_type='C')
    c4.position = (2, 3, 0)
    c4.bind_to(c2, 'bondType1')
    c5 = Atom(name='C5', atom_type='C')
    c5.position = (3, 1, 0)
    c5.bind_to(c3, 'bondType1')
    c6 = Atom(name='C6', atom_type='C')
    c6.position = (3, 2, 0)
    c6.bind_to(c5, 'bondType1')
    c6.bind_to(c4, 'bondType1')
    c7 = Atom(name='C7', atom_type='C')
    c7.position = (4, 2, 0)
    c7.bind_to(c5, 'bondType1')
    c10 = Atom(name='C10', atom_type='C')
    c10.position = (4, 1, 0)
    c10.bind_to(c7, 'bondType1')
    n1 = Atom(name='N1', atom_type='N')
    n1.position = (4, 3, 0)
    n1.bind_to(c6, 'bondType1')
    c8 = Atom(name='C8', atom_type='C')
    c8.position = (5, 1, 0)
    c8.bind_to(c7, 'bondType1')
    c8.bind_to(n1, 'bondType1')
    c9 = Atom(name='C9', atom_type='C')
    c9.position = (6, 1, 0)
    c9.bind_to(c8, 'bondType1')
    left_atoms = [c1, c2, c3, cl1, c4, c5, c6, c7, c8, c9, c10, n1]

    # construct Ligand 2
    cl11 = Atom(name='Cl11', atom_type='Cl')
    cl11.position = (1, 1, 0)
    c11 = Atom(name='C11', atom_type='C')
    c11.position = (2, 1, 0)
    c12 = Atom(name='C12', atom_type='C')
    c12.position = (2, 2, 0)
    c12.bind_to(c11, 'bondType1')
    c12.bind_to(cl11, 'bondType1')
    c13 = Atom(name='C13', atom_type='C')
    c13.position = (3, 1, 0)
    c13.bind_to(c11, 'bondType1')
    c14 = Atom(name='C14', atom_type='C')
    c14.position = (3, 2, 0)
    c14.bind_to(c12, 'bondType1')
    c15 = Atom(name='C15', atom_type='C')
    c15.position = (4, 1, 0)
    c15.bind_to(c13, 'bondType1')
    c16 = Atom(name='C16', atom_type='C')
    c16.position = (4, 2, 0)
    c16.bind_to(c15, 'bondType1')
    c16.bind_to(c14, 'bondType1')
    c17 = Atom(name='C17', atom_type='C')
    c17.position = (5, 2, 0)
    c17.bind_to(c15, 'bondType1')
    c20 = Atom(name='C20', atom_type='C')
    c20.position = (5, 1, 0)
    c20.bind_to(c17, 'bondType1')
    n11 = Atom(name='N11', atom_type='N')
    n11.position = (5, 3, 0)
    n11.bind_to(c16, 'bondType1')
    c18 = Atom(name='C18', atom_type='C')
    c18.position = (6, 1, 0)
    c18.bind_to(c17, 'bondType1')
    c18.bind_to(n11, 'bondType1')
    c19 = Atom(name='C19', atom_type='C')
    c19.position = (7, 1, 0)
    c19.bind_to(c18, 'bondType1')
    right_atoms = [cl11, c11, c12, c13, c14, c15, c16, c17, c18, c19, c20, n11]

    # the correct solution
    suptops = _overlay(c9, c19, parent_n1=None, parent_n2=None, bond_types=(None, None),
                      suptop=SuperimposedTopology(left_atoms, right_atoms))
    suptop = suptops.pop(0)

    # (c1, c11) should be linked to (c2, c12), (c3, c13)
    assert {x[0] for x in suptop.matched_pairs_bonds[(c1, c11)]} == {(c2, c12), (c3, c13)}
    # (c2, c12) should be linked to (c1, c11), (c4, c14)
    assert {x[0] for x in suptop.matched_pairs_bonds[(c2, c12)]} == {(c1, c11), (c4, c14)}
    # etc
    assert {x[0] for x in suptop.matched_pairs_bonds[(c3, c13)]} == {(c1, c11), (c5, c15)}
    assert {x[0] for x in suptop.matched_pairs_bonds[(c4, c14)]} == {(c2, c12), (c6, c16)}
    assert {x[0] for x in suptop.matched_pairs_bonds[(c5, c15)]} == {(c3, c13), (c6, c16), (c7, c17)}
    assert {x[0] for x in suptop.matched_pairs_bonds[(c6, c16)]} == {(c5, c15), (c4, c14), (n1, n11)}
    assert {x[0] for x in suptop.matched_pairs_bonds[(c7, c17)]} == {(c5, c15), (c10, c20), (c8, c18)}
    assert {x[0] for x in suptop.matched_pairs_bonds[(c8, c18)]} == {(c7, c17), (n1, n11), (c9, c19)}
    assert {x[0] for x in suptop.matched_pairs_bonds[(c9, c19)]} == {(c8, c18)}
    assert {x[0] for x in suptop.matched_pairs_bonds[(c10, c20)]} == {(c7, c17)}
    assert {x[0] for x in suptop.matched_pairs_bonds[(n1, n11)]} == {(c6, c16), (c8, c18)}


def test_tyk2_l11l14_part():
    """
    Molecule inspired by Agastya's dataset.
         n1             n11
         |               |
         C1-O1          C11-O11
         |                |
      H1-C2           H11-C12
        /  \            /    \
    H2-C3- C4-H3   H12-C13 - C14 - H13
       /    \           /      \
      F1    H4         CL1     H14
    """
    # construct Ligand 1
    n1 = Atom(name='N1', atom_type='N')
    n1.position = (1, 1, 0)
    c1 = Atom(name='C1', atom_type='C')
    c1.position = (2, 1, 0)
    c1.bind_to(n1, 'bondType1')
    o1 = Atom(name='O1', atom_type='O')
    o1.position = (2, 2, 0)
    o1.bind_to(c1, 'bondType1')
    h1 = Atom(name='H1', atom_type='H')
    h1.position = (3, 1, 0)
    c2 = Atom(name='C2', atom_type='C')
    c2.position = (3, 2, 0)
    c2.bind_to(h1, 'bondType1')
    c2.bind_to(c1, 'bondType1')
    h2 = Atom(name='H2', atom_type='H')
    h2.position = (4, 1, 0)
    c3 = Atom(name='C3', atom_type='C')
    c3.position = (4, 2, 0)
    c3.bind_to(c2, 'bondType1')
    c3.bind_to(h2, 'bondType1')
    c4 = Atom(name='C4', atom_type='C')
    c4.position = (4, 3, 0)
    c4.bind_to(c2, 'bondType1')
    c4.bind_to(c3, 'bondType1')
    h3 = Atom(name='H3', atom_type='H')
    h3.position = (4, 4, 0)
    h3.bind_to(c4, 'bondType1')
    f1 = Atom(name='F1', atom_type='F')
    f1.position = (5, 1, 0)
    f1.bind_to(c3, 'bondType1')
    h4 = Atom(name='H4', atom_type='H')
    h4.position = (5, 2, 0)
    h4.bind_to(c4, 'bondType1')
    left_atoms = [n1, c1, o1, h1, c2, h2, c3, c4, h3, f1, h4]

    # construct Ligand 2
    n11 = Atom(name='N11', atom_type='N')
    n11.position = (1, 1, 0)
    c11 = Atom(name='C11', atom_type='C')
    c11.position = (2, 1, 0)
    c11.bind_to(n11, 'bondType1')
    o11 = Atom(name='O11', atom_type='O')
    o11.position = (2, 2, 0)
    o11.bind_to(c11, 'bondType1')
    h11 = Atom(name='H11', atom_type='H')
    h11.position = (3, 1, 0)
    c12 = Atom(name='C12', atom_type='C')
    c12.position = (3, 2, 0)
    c12.bind_to(h11, 'bondType1')
    c12.bind_to(c11, 'bondType1')
    h12 = Atom(name='H12', atom_type='H')
    h12.position = (4, 1, 0)
    c13 = Atom(name='C13', atom_type='C')
    c13.position = (4, 2, 0)
    c13.bind_to(c12, 'bondType1')
    c13.bind_to(h12, 'bondType1')
    c14 = Atom(name='C14', atom_type='C')
    c14.position = (4, 3, 0)
    c14.bind_to(c12, 'bondType1')
    c14.bind_to(c13, 'bondType1')
    h13 = Atom(name='H13', atom_type='H')
    h13.position = (4, 4, 0)
    h13.bind_to(c14, 'bondType1')
    cl11 = Atom(name='CL11', atom_type='CL')
    cl11.position = (5, 1, 0)
    cl11.bind_to(c13, 'bondType1')
    h14 = Atom(name='H14', atom_type='H')
    h14.position = (5, 2, 0)
    h14.bind_to(c14, 'bondType1')
    right_atoms = [n11, c11, o11, h11, c12, h12, c13, c14, h13, cl11, h14]

    # the correct solution
    suptops = _overlay(n1, n11, parent_n1=None, parent_n2=None, bond_types=(None, None),
                      suptop=SuperimposedTopology(left_atoms, right_atoms))
    suptop = suptops.pop(0)
    assert len(suptop) == 10

    matching_pairs = [('N1', 'N11'), ('C1', 'C11'), ('O1', 'O11'),
            ('H1', 'H11'), ('C2', 'C12'), ('H2', 'H12'), ('C3', 'C13'),
            ('C4', 'C14'), ('H3', 'H13'), ('H4', 'H14')]
    for n1, n2 in matching_pairs[::-1]:
        if suptop.contains_atom_name_pair(n1, n2):
            matching_pairs.remove((n1, n2))
    assert len(matching_pairs) == 0, matching_pairs
