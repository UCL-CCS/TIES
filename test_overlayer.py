"""
Focuses on the _overlay function that actually traverses the molecule using the given starting points.

TODO:
- add more complex test case with multiple solutions (2 levels)  similar to test_SimpleMultipleSolutions_rightStart
- add a test case that shows that molecules cannot be fully superimposed when it is broken in half
- add more test cases that take into account the "differences"
"""

from topology_superimposer import AtomNode, _overlay


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
    c1 = AtomNode(name='C1', type='C')
    c1.set_position(x=1, y=1, z=0)
    n1 = AtomNode(name='N1', type='N')
    n1.set_position(x=1, y=2, z=0)
    c1.bindTo(n1, 'bondType1')

    # construct the LIGAND 2
    c11 = AtomNode(name='C11', type='C')
    c11.set_position(x=1, y=1, z=0)
    n11 = AtomNode(name='N11', type='N')
    n11.set_position(x=1, y=2, z=0)
    c11.bindTo(n11, 'bondType1')

    # should return a list with an empty sup_top
    suptops = _overlay(c1, n11, parent_n1=None, parent_n2=None, bond_types=(None, None))
    # it should return an empty suptop
    assert suptops is None


def test_2diffAtoms_CN_rightStart():
    """
    Two different Atoms. Only one solution exists.

     LIGAND 1        LIGAND 2
        C1              C11
        \                \
        N1              N11
    """
    # construct the LIGAND 1
    c1 = AtomNode(name='C1', type='C')
    c1.set_position(x=1, y=1, z=0)
    n1 = AtomNode(name='N1', type='N')
    n1.set_position(x=1, y=2, z=0)
    c1.bindTo(n1, 'bondType1')

    # construct the LIGAND 2
    c11 = AtomNode(name='C11', type='C')
    c11.set_position(x=1, y=1, z=0)
    n11 = AtomNode(name='N11', type='N')
    n11.set_position(x=1, y=2, z=0)
    c11.bindTo(n11, 'bondType1')

    # should overlap 2 atoms
    suptops = _overlay(c1, c11, parent_n1=None, parent_n2=None, bond_types=(None, None))
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
    c1 = AtomNode(name='C1', type='C')
    c1.set_position(x=1, y=1, z=0)
    n1 = AtomNode(name='N1', type='N')
    n1.set_position(x=1, y=2, z=0)
    c1.bindTo(n1, 'bondType1')
    o1 = AtomNode(name='O1', type='O')
    o1.set_position(x=1, y=3, z=0)
    o1.bindTo(n1, 'bondType1')

    # construct the LIGAND 2
    c11 = AtomNode(name='C11', type='C')
    c11.set_position(x=1, y=1, z=0)
    n11 = AtomNode(name='N11', type='N')
    n11.set_position(x=1, y=2, z=0)
    c11.bindTo(n11, 'bondType1')
    o11 = AtomNode(name='O11', type='O')
    o11.set_position(x=1, y=3, z=0)
    o11.bindTo(n11, 'bondType1')

    # should overlap 2 atoms
    suptop = _overlay(c1, c11, parent_n1=None, parent_n2=None, bond_types=(None, None))
    assert suptop is not None

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

     LIGAND 1        LIGAND 2
        C1              C11
        \                \
        N1              N11
        /\              / \
     O1    O2        O11   O12

    """
    # ignore the third coordinate dimension
    # construct the LIGAND 1
    c1 = AtomNode(name='C1', type='C')
    c1.set_position(x=1, y=1, z=0)
    n1 = AtomNode(name='N1', type='N')
    n1.set_position(x=1, y=2, z=0)
    c1.bindTo(n1, 'bondType1')
    o1 = AtomNode(name='O1', type='O')
    o1.set_position(x=1, y=3, z=0)
    o1.bindTo(n1, 'bondType1')
    o2 = AtomNode(name='O2', type='O')
    o2.set_position(x=2, y=3, z=0)
    o2.bindTo(n1, 'bondType1')

    # construct the LIGAND 2
    c11 = AtomNode(name='C11', type='C')
    c11.set_position(x=1, y=1, z=0)
    n11 = AtomNode(name='N11', type='N')
    n11.set_position(x=1, y=2, z=0)
    c11.bindTo(n11, 'bondType1')
    o11 = AtomNode(name='O11', type='O')
    o11.set_position(x=1, y=3, z=0)
    o11.bindTo(n11, 'bondType1')
    o12 = AtomNode(name='O12', type='O')
    o12.set_position(x=2, y=3, z=0)
    o12.bindTo(n11, 'bondType1')

    # the good solution is (O1-O11) and (O2-O12)

    # should generate one topology which has one "symmetry"
    suptop = _overlay(c1, c11, parent_n1=None, parent_n2=None, bond_types=(None, None))
    assert suptop is not None
    assert len(suptop) == 4

    assert len(suptop.mirrors) == 1
    worse_st = suptop.mirrors[0]

    # check if both representations were found
    # The ester allows for mapping (O1-O11, O2-O12) and (O1-O12, O2-O11)
    assert suptop.contains_atomNamePair('O1', 'O11') and suptop.contains_atomNamePair('O2', 'O12')
    assert worse_st.contains_atomNamePair('O1', 'O12') and worse_st.contains_atomNamePair('O2', 'O11')

    correct_overlaps = [('C1', 'C11'), ('N1', 'N11')]
    for atomName1, atomName2 in correct_overlaps:
        assert suptop.contains_atomNamePair(atomName1, atomName2)


def test_One2Many():
    """

     LIGAND 1        LIGAND 2
        N1              N11
        /              /   \
     O1              O11   O12

    """
    # ignore the third coordinate dimension
    # construct the LIGAND 1
    n1 = AtomNode(name='N1', type='N')
    n1.set_position(x=1, y=1, z=0)
    o1 = AtomNode(name='O1', type='O')
    o1.set_position(x=2, y=1, z=0)
    o1.bindTo(n1, 'bondType1')

    # construct the LIGAND 2
    n11 = AtomNode(name='N11', type='N')
    n11.set_position(x=1, y=1, z=0)
    o11 = AtomNode(name='O11', type='O')
    o11.set_position(x=2, y=1, z=0)
    o11.bindTo(n11, 'bondType1')
    o12 = AtomNode(name='O12', type='O')
    o12.set_position(x=2, y=2, z=0)
    o12.bindTo(n11, 'bondType1')

    # the good solution is (O1-O11)

    # should generate one topology which has one "symmetry"
    suptop = _overlay(n1, n11, parent_n1=None, parent_n2=None, bond_types=(None, None))
    assert suptop is not None
    assert len(suptop) == 2

    assert len(suptop.alternative_mappings) == 1
    worse_st = suptop.alternative_mappings[0]

    # check if both representations were found
    # The ester allows for mapping (O1-O11, O2-O12) and (O1-O12, O2-O11)
    assert suptop.contains_atomNamePair('O1', 'O11')
    assert worse_st.contains_atomNamePair('O1', 'O12')

    assert suptop.contains_atomNamePair('N1', 'N11')


def test_Many2One():
    """

     LIGAND 1        LIGAND 2
        N1              N11
        / \              /
     O1    O2          O11

    """
    # ignore the third coordinate dimension
    # construct the LIGAND 1
    n1 = AtomNode(name='N1', type='N')
    n1.set_position(x=1, y=1, z=0)
    o1 = AtomNode(name='O1', type='O')
    o1.set_position(x=2, y=1, z=0)
    o1.bindTo(n1, 'bondType1')
    o2 = AtomNode(name='O2', type='O')
    o2.set_position(x=2, y=2, z=0)
    o2.bindTo(n1, 'bondType1')

    # construct the LIGAND 2
    n11 = AtomNode(name='N11', type='N')
    n11.set_position(x=1, y=1, z=0)
    o11 = AtomNode(name='O11', type='O')
    o11.set_position(x=2, y=1, z=0)
    o11.bindTo(n11, 'bondType1')

    # the good solution is (O1-O11)

    # should generate one topology which has one "symmetry"
    suptop = _overlay(n1, n11, parent_n1=None, parent_n2=None, bond_types=(None, None))
    assert suptop is not None
    assert len(suptop) == 2

    assert len(suptop.alternative_mappings) == 1
    worse_st = suptop.alternative_mappings[0]

    # check if both representations were found
    # The ester allows for mapping (O1-O11, O2-O12) and (O1-O12, O2-O11)
    assert suptop.contains_atomNamePair('O1', 'O11')
    assert worse_st.contains_atomNamePair('O1', 'O12')

    assert suptop.contains_atomNamePair('N1', 'N11')


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
    c1 = AtomNode(name='C1', type='C')
    c1.set_position(x=1, y=1, z=0)
    o1 = AtomNode(name='O1', type='O')
    o1.set_position(x=1, y=3, z=0)
    o1.bindTo(c1, 'bondType1')
    o2 = AtomNode(name='O2', type='O')
    o2.set_position(x=2, y=3, z=0)
    o2.bindTo(c1, 'bondType1')
    n1 = AtomNode(name='N1', type='N')
    n1.set_position(3, 1, 0)
    n1.bindTo(o1, 'bondType1')
    n2 = AtomNode(name='N2', type='N')
    n2.set_position(3, 2, 0)
    n2.bindTo(o1, 'bondType1')
    n3 = AtomNode(name='N3', type='N')
    n3.set_position(3, 3, 0)
    n3.bindTo(o2, 'bondType1')
    n4 = AtomNode(name='N4', type='N')
    n4.set_position(3, 4, 0)
    n4.bindTo(o2, 'bondType1')

    # construct the LIGAND 2
    c11 = AtomNode(name='C11', type='C')
    c11.set_position(x=1, y=1, z=0)
    o11 = AtomNode(name='O11', type='O')
    o11.set_position(x=1, y=3, z=0)
    o11.bindTo(c11, 'bondType1')
    o12 = AtomNode(name='O12', type='O')
    o12.set_position(x=2, y=3, z=0)
    o12.bindTo(c11, 'bondType1')
    n11 = AtomNode(name='N11', type='N')
    n11.set_position(3, 1, 0)
    n11.bindTo(o11, 'bondType1')
    n12 = AtomNode(name='N12', type='N')
    n12.set_position(3, 2, 0)
    n12.bindTo(o11, 'bondType1')
    n13 = AtomNode(name='N13', type='N')
    n13.set_position(3, 3, 0)
    n13.bindTo(o12, 'bondType1')
    n14 = AtomNode(name='N14', type='N')
    n14.set_position(3, 4, 0)
    n14.bindTo(o12, 'bondType1')

    # the good solution is (O1-O11) and (O2-O12)

    # should generate one topology which has one "symmetry"
    suptops = _overlay(c1, c11, parent_n1=None, parent_n2=None, bond_types=(None, None))
    assert len(suptops) == 1
    suptop = suptops[0]

    # check if the main solution is correct
    correct_overlaps = [('C1', 'C11'),
                        ('O1', 'O11'), ('O2', 'O12'),
                        ('N1', 'N11'), ('N2', 'N12'), ('N3', 'N13'), ('N4', 'N14')]
    for atomName1, atomName2 in correct_overlaps:
        assert suptop.contains_atomNamePair(atomName1, atomName2)

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
    c1 = AtomNode(name='C1', type='C')
    c1.set_position(x=1, y=1, z=0)
    c2 = AtomNode(name='C2', type='C')
    c2.set_position(x=1, y=2, z=0)
    c1.bindTo(c2, 'bondType1')

    # construct the LIGAND 2
    c11 = AtomNode(name='C11', type='C')
    c11.set_position(x=1, y=1, z=0)
    c12 = AtomNode(name='C12', type='C')
    c12.set_position(x=1, y=2, z=0)
    c11.bindTo(c12, 'bondType1')

    # should return a list with an empty sup_top
    suptops = _overlay(c1, c11, parent_n1=None, parent_n2=None, bond_types=(None, None))
    assert len(suptops) == 1
    assert suptops[0].contains_atomNamePair('C1', 'C11')
    assert suptops[0].contains_atomNamePair('C2', 'C12')

    suptops = _overlay(c1, c12, parent_n1=None, parent_n2=None, bond_types=(None, None))
    assert len(suptops) == 1
    assert suptops[0].contains_atomNamePair('C1', 'C12')
    assert suptops[0].contains_atomNamePair('C2', 'C11')


def test_methyl():
    """
    Two solutions with different starting points.

     LIGAND 1        LIGAND 2
        C1              C11
       / | \            / | \
     H1 H2 H3        H11 H12 H13
    """
    # construct the LIGAND 1
    c1 = AtomNode(name='C1', type='C')
    c1.set_position(x=1, y=1, z=0)
    h1 = AtomNode(name='H1', type='H')
    h1.set_position(x=2, y=1, z=0)
    c1.bindTo(h1, 'bondType1')
    h2 = AtomNode(name='H2', type='H')
    h2.set_position(x=2, y=2, z=0)
    c1.bindTo(h2, 'bondType1')
    h3 = AtomNode(name='H3', type='H')
    h3.set_position(x=2, y=3, z=0)
    c1.bindTo(h3, 'bondType1')

    # construct the LIGAND 2
    c11 = AtomNode(name='C11', type='C')
    c11.set_position(x=1, y=1, z=0)
    h11 = AtomNode(name='H11', type='H')
    h11.set_position(x=2, y=1, z=0)
    c11.bindTo(h11, 'bondType1')
    h12 = AtomNode(name='H12', type='H')
    h12.set_position(x=2, y=2, z=0)
    c11.bindTo(h12, 'bondType1')
    h13 = AtomNode(name='H13', type='H')
    h13.set_position(x=2, y=3, z=0)
    c11.bindTo(h13, 'bondType1')

    # should return a list with an empty sup_top
    suptops = _overlay(c1, c11, parent_n1=None, parent_n2=None, bond_types=(None, None))
    assert len(suptops) == 1
    assert suptops[0].contains_atomNamePair('C1', 'C11')
    assert suptops[0].contains_atomNamePair('C2', 'C12')


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
    c1 = AtomNode(name='C1', type='C')
    c1.set_position(x=1, y=1, z=0)
    s1 = AtomNode(name='S2', type='S')
    s1.set_position(x=2, y=1, z=0)
    c1.bindTo(s1, 'bondType1')
    n1 = AtomNode(name='N1', type='N')
    n1.set_position(x=3, y=1, z=0)
    n1.bindTo(s1, 'bondType1')

    # construct the LIGAND 2
    c11 = AtomNode(name='C11', type='C')
    c11.set_position(x=1, y=1, z=0)
    o11 = AtomNode(name='O11', type='O')
    o11.set_position(x=2, y=1, z=0)
    c11.bindTo(o11, 'bondType1')
    n11 = AtomNode(name='N11', type='N')
    n11.set_position(x=3, y=1, z=0)
    n11.bindTo(o11, 'bondType1')

    # should return a list with an empty sup_top
    suptops = _overlay(c1, c11, parent_n1=None, parent_n2=None, bond_types=(None, None))
    assert len(suptops) == 1
    assert suptops[0].contains_atomNamePair('C1', 'C11')

    suptops = _overlay(n1, n11, parent_n1=None, parent_n2=None, bond_types=(None, None))
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
    c1 = AtomNode(name='C1', type='C')
    c1.set_position(x=1, y=1, z=0)
    c2 = AtomNode(name='C2', type='C')
    c2.set_position(x=2, y=1, z=0)
    c2.bindTo(c1, 'bondType1')
    c3 = AtomNode(name='C3', type='C')
    c3.set_position(x=2, y=1, z=0)
    c3.bindTo(c2, 'bondType1')
    n1 = AtomNode(name='N1', type='N')
    n1.set_position(x=3, y=1, z=0)
    n1.bindTo(c3, 'bondType1')

    # construct the LIGAND 2
    c11 = AtomNode(name='C11', type='C')
    c11.set_position(x=1, y=1, z=0)
    c12 = AtomNode(name='C12', type='C')
    c12.set_position(x=2, y=1, z=0)
    c12.bindTo(c11, 'bondType1')
    o11 = AtomNode(name='O11', type='O')
    o11.set_position(x=3, y=1, z=0)
    o11.bindTo(c12, 'bondType1')
    n11 = AtomNode(name='N11', type='N')
    n11.set_position(x=3, y=1, z=0)
    n11.bindTo(o11, 'bondType1')

    # should return a list with an empty sup_top
    suptops = _overlay(c1, c11, parent_n1=None, parent_n2=None, bond_types=(None, None))
    assert len(suptops) == 1
    assert suptops[0].contains_atomNamePair('C1', 'C11')

    suptops = _overlay(n1, n11, parent_n1=None, parent_n2=None, bond_types=(None, None))
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
    c1 = AtomNode(name='C1', type='C')
    c1.set_position(x=1, y=1, z=0)
    c2 = AtomNode(name='C2', type='C')
    c2.set_position(x=1, y=2, z=0)
    c1.bindTo(c2, 'bondType1')
    c3 = AtomNode(name='C3', type='C')
    c3.set_position(x=2, y=2, z=0)
    c3.bindTo(c1, 'bondType1')
    c3.bindTo(c2, 'bondType1')

    # construct the LIGAND 2
    c11 = AtomNode(name='C11', type='C')
    c11.set_position(x=1, y=1, z=0)
    c12 = AtomNode(name='C12', type='C')
    c12.set_position(x=1, y=2, z=0)
    c11.bindTo(c12, 'bondType1')
    c13 = AtomNode(name='C13', type='C')
    c13.set_position(x=2, y=2, z=0)
    c13.bindTo(c11, 'bondType1')
    c13.bindTo(c12, 'bondType1')

    suptops = _overlay(c1, c11, parent_n1=None, parent_n2=None, bond_types=(None, None))
    # there is one main component - that is the good solution
    assert len(suptops) == 1
    suptop = suptops[0]
    # there is one symmetrical way to traverse it
    assert len(suptop.mirrors) == 1
    wrong_st = suptop.mirrors[0]

    assert suptop.contains_atomNamePair('C2', 'C12') and suptop.contains_atomNamePair('C3', 'C13')
    assert wrong_st.contains_atomNamePair('C2', 'C13') and wrong_st.contains_atomNamePair('C3', 'C12')
    # both solutions should have the same starting situation
    assert all(st.contains_atomNamePair('C1', 'C11') for st in [suptop, wrong_st])
    # there should be one circle in each
    assert all(st.sameCircleNumber() for st in [suptop, wrong_st])
    assert all(st.getCircleNumber() == (1, 1) for st in [suptop, wrong_st])

    suptops = _overlay(c1, c12, parent_n1=None, parent_n2=None, bond_types=(None, None))
    assert len(suptops) == 1
    # there should be one circle in each
    assert suptops[0].sameCircleNumber()
    assert suptops[0].getCircleNumber() == (1, 1)

    suptops = _overlay(c1, c13, parent_n1=None, parent_n2=None, bond_types=(None, None))
    assert len(suptops) == 1
    # there should be one circle in each
    assert suptops[0].sameCircleNumber()
    assert suptops[0].getCircleNumber() == (1, 1)

    suptops = _overlay(c2, c11, parent_n1=None, parent_n2=None, bond_types=(None, None))
    assert len(suptops) == 1
    # there should be one circle in each
    assert suptops[0].sameCircleNumber()
    assert suptops[0].getCircleNumber() == (1, 1)

    suptops = _overlay(c2, c12, parent_n1=None, parent_n2=None, bond_types=(None, None))
    assert len(suptops) == 1
    # there should be one circle in each
    assert suptops[0].sameCircleNumber()
    assert suptops[0].getCircleNumber() == (1, 1)

    suptops = _overlay(c2, c13, parent_n1=None, parent_n2=None, bond_types=(None, None))
    assert len(suptops) == 1
    # there should be one circle in each
    assert suptops[0].sameCircleNumber()
    assert suptops[0].getCircleNumber() == (1, 1)

    suptops = _overlay(c3, c11, parent_n1=None, parent_n2=None, bond_types=(None, None))
    assert len(suptops) == 1
    # there should be one circle in each
    assert suptops[0].sameCircleNumber()
    assert suptops[0].getCircleNumber() == (1, 1)

    suptops = _overlay(c3, c12, parent_n1=None, parent_n2=None, bond_types=(None, None))
    assert len(suptops) == 1
    # there should be one circle in each
    assert suptops[0].sameCircleNumber()
    assert suptops[0].getCircleNumber() == (1, 1)

    suptops = _overlay(c3, c13, parent_n1=None, parent_n2=None, bond_types=(None, None))
    assert len(suptops) == 1
    # there should be one circle in each
    assert suptops[0].sameCircleNumber()
    assert suptops[0].getCircleNumber() == (1, 1)


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
    c1 = AtomNode(name='C1', type='C')
    c1.set_position(x=1, y=1, z=0)
    c2 = AtomNode(name='C2', type='C')
    c2.set_position(x=1, y=2, z=0)
    c1.bindTo(c2, 'bondType1')
    c3 = AtomNode(name='C3', type='C')
    c3.set_position(x=2, y=2, z=0)
    c3.bindTo(c1, 'bondType1')
    cl1 = AtomNode(name='CL1', type='Cl')
    cl1.set_position(x=2, y=1, z=0)
    cl1.bindTo(c3, 'bondType1')
    c4 = AtomNode(name='C4', type='C')
    c4.set_position(x=2, y=3, z=0)
    c4.bindTo(c2, 'bondType1')
    c5 = AtomNode(name='C5', type='C')
    c5.set_position(x=3, y=1, z=0)
    c5.bindTo(c3, 'bondType1')
    c6 = AtomNode(name='C6', type='C')
    c6.set_position(x=3, y=2, z=0)
    c6.bindTo(c5, 'bondType1')
    c6.bindTo(c4, 'bondType1')
    c7 = AtomNode(name='C7', type='C')
    c7.set_position(x=4, y=2, z=0)
    c7.bindTo(c5, 'bondType1')
    c10 = AtomNode(name='C10', type='C')
    c10.set_position(x=4, y=1, z=0)
    c10.bindTo(c7, 'bondType1')
    n1 = AtomNode(name='N1', type='N')
    n1.set_position(x=4, y=3, z=0)
    n1.bindTo(c6, 'bondType1')
    c8 = AtomNode(name='C8', type='C')
    c8.set_position(x=5, y=1, z=0)
    c8.bindTo(c7, 'bondType1')
    c8.bindTo(n1, 'bondType1')
    c9 = AtomNode(name='C9', type='C')
    c9.set_position(x=6, y=1, z=0)
    c9.bindTo(c8, 'bondType1')

    # construct Ligand 2
    cl11 = AtomNode(name='Cl11', type='Cl')
    cl11.set_position(x=1, y=1, z=0)
    c11 = AtomNode(name='C11', type='C')
    c11.set_position(x=2, y=1, z=0)
    c12 = AtomNode(name='C12', type='C')
    c12.set_position(x=2, y=2, z=0)
    c12.bindTo(c11, 'bondType1')
    c12.bindTo(cl11, 'bondType1')
    c13 = AtomNode(name='C13', type='C')
    c13.set_position(x=3, y=1, z=0)
    c13.bindTo(c11, 'bondType1')
    c14 = AtomNode(name='C14', type='C')
    c14.set_position(x=3, y=2, z=0)
    c14.bindTo(c12, 'bondType1')
    c15 = AtomNode(name='C15', type='C')
    c15.set_position(x=4, y=1, z=0)
    c15.bindTo(c13, 'bondType1')
    c16 = AtomNode(name='C16', type='C')
    c16.set_position(x=4, y=2, z=0)
    c16.bindTo(c15, 'bondType1')
    c16.bindTo(c14, 'bondType1')
    c17 = AtomNode(name='C17', type='C')
    c17.set_position(x=5, y=2, z=0)
    c17.bindTo(c15, 'bondType1')
    c20 = AtomNode(name='C20', type='C')
    c20.set_position(x=5, y=1, z=0)
    c20.bindTo(c17, 'bondType1')
    n11 = AtomNode(name='N11', type='N')
    n11.set_position(x=5, y=3, z=0)
    n11.bindTo(c16, 'bondType1')
    c18 = AtomNode(name='C18', type='C')
    c18.set_position(x=6, y=1, z=0)
    c18.bindTo(c17, 'bondType1')
    c18.bindTo(n11, 'bondType1')
    c19 = AtomNode(name='C19', type='C')
    c19.set_position(x=7, y=1, z=0)
    c19.bindTo(c18, 'bondType1')

    # the correct solution
    suptops = _overlay(c9, c19, parent_n1=None, parent_n2=None, bond_types=(None, None))
    assert len(suptops) == 1
    assert len(suptops[0]) == 11

    """
    This is a rare case around which we'll have to find a work around.
    Basically, the best solution that follows the basic traversal allows for a superimposition
    that should not be allowed. So additional additional way has to be checked to discredit it 
    """
    suptops = _overlay(c5, c14, parent_n1=None, parent_n2=None, bond_types=(None, None))
    assert len(suptops) == 1
    assert len(suptops[0]) != 12


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
    n1 = AtomNode(name='N1', type='N')
    n1.set_position(x=1, y=1, z=0)
    c1 = AtomNode(name='C1', type='C')
    c1.set_position(x=2, y=1, z=0)
    c1.bindTo(n1, 'bondType1')
    o1 = AtomNode(name='O1', type='O')
    o1.set_position(x=2, y=2, z=0)
    o1.bindTo(c1, 'bondType1')
    h1 = AtomNode(name='H1', type='H')
    h1.set_position(x=3, y=1, z=0)
    c2 = AtomNode(name='C2', type='C')
    c2.set_position(x=3, y=2, z=0)
    c2.bindTo(h1, 'bondType1')
    c2.bindTo(c1, 'bondType1')
    h2 = AtomNode(name='H2', type='H')
    h2.set_position(x=4, y=1, z=0)
    c3 = AtomNode(name='C3', type='C')
    c3.set_position(x=4, y=2, z=0)
    c3.bindTo(c2, 'bondType1')
    c3.bindTo(h2, 'bondType1')
    c4 = AtomNode(name='C4', type='C')
    c4.set_position(x=4, y=3, z=0)
    c4.bindTo(c2, 'bondType1')
    c4.bindTo(c3, 'bondType1')
    h3 = AtomNode(name='H3', type='H')
    h3.set_position(x=4, y=4, z=0)
    h3.bindTo(c4, 'bondType1')
    f1 = AtomNode(name='F1', type='F')
    f1.set_position(x=5, y=1, z=0)
    f1.bindTo(c3, 'bondType1')
    h4 = AtomNode(name='H4', type='H')
    h4.set_position(x=5, y=2, z=0)
    h4.bindTo(c4, 'bondType1')

    # construct Ligand 2
    n11 = AtomNode(name='N11', type='N')
    n11.set_position(x=1, y=1, z=0)
    c11 = AtomNode(name='C11', type='C')
    c11.set_position(x=2, y=1, z=0)
    c11.bindTo(n11, 'bondType1')
    o11 = AtomNode(name='O11', type='O')
    o11.set_position(x=2, y=2, z=0)
    o11.bindTo(c11, 'bondType1')
    h11 = AtomNode(name='H11', type='H')
    h11.set_position(x=3, y=1, z=0)
    c12 = AtomNode(name='C12', type='C')
    c12.set_position(x=3, y=2, z=0)
    c12.bindTo(h11, 'bondType1')
    c12.bindTo(c11, 'bondType1')
    h12 = AtomNode(name='H12', type='H')
    h12.set_position(x=4, y=1, z=0)
    c13 = AtomNode(name='C13', type='C')
    c13.set_position(x=4, y=2, z=0)
    c13.bindTo(c12, 'bondType1')
    c13.bindTo(h12, 'bondType1')
    c14 = AtomNode(name='C14', type='C')
    c14.set_position(x=4, y=3, z=0)
    c14.bindTo(c12, 'bondType1')
    c14.bindTo(c13, 'bondType1')
    h13 = AtomNode(name='H13', type='H')
    h13.set_position(x=4, y=4, z=0)
    h13.bindTo(c14, 'bondType1')
    cl11 = AtomNode(name='CL11', type='CL')
    cl11.set_position(x=5, y=1, z=0)
    cl11.bindTo(c13, 'bondType1')
    h14 = AtomNode(name='H14', type='H')
    h14.set_position(x=5, y=2, z=0)
    h14.bindTo(c14, 'bondType1')

    # the correct solution
    suptops = _overlay(n1, n11, parent_n1=None, parent_n2=None, bond_types=(None, None))
    assert len(suptops) == 1
    suptop = suptops[0]
    assert len(suptop) == 10

    matching_pairs = [('N1', 'N11'), ('C1', 'C11'), ('O1', 'O11'),
            ('H1', 'H11'), ('C2', 'C12'), ('H2', 'H12'), ('C3', 'C13'),
            ('C4', 'C14'), ('H3', 'H13'), ('H4', 'H14')]
    for n1, n2 in matching_pairs:
        assert suptop.contains_atomNamePair(n1, n2), (n1, n2)
