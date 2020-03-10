"""
Focuses on the _superimpose_topology function that
superimposes the molecule in many different and
then processes the outputs to ensure the best match is found.

TODO
- this testing module should focus on "separated" molecules
"""

from topology_superimposer import _superimpose_topologies, AtomNode


def test_2diffAtoms_CN():
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
    top1_list = [c1, n1]

    # construct the LIGAND 2
    c11 = AtomNode(name='C11', type='C')
    c11.set_position(x=1, y=1, z=0)
    n11 = AtomNode(name='N11', type='N')
    n11.set_position(x=1, y=2, z=0)
    c11.bindTo(n11, 'bondType1')
    top2_list = [c11, n11]

    # should return a list with an empty sup_top
    suptops = _superimpose_topologies(top1_list, top2_list)
    # it should return an empty suptop
    assert len(suptops) == 1
    assert len(suptops[0]) == 2


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
    top1_list = [c1, n1, o1]

    # construct the LIGAND 2
    c11 = AtomNode(name='C11', type='C')
    c11.set_position(x=1, y=1, z=0)
    n11 = AtomNode(name='N11', type='N')
    n11.set_position(x=1, y=2, z=0)
    c11.bindTo(n11, 'bondType1')
    o11 = AtomNode(name='O11', type='O')
    o11.set_position(x=1, y=3, z=0)
    o11.bindTo(n11, 'bondType1')
    top2_list = [c11, n11, o11]

    # should overlap 2 atoms
    suptops = _superimpose_topologies(top1_list, top2_list)
    assert len(suptops) == 1
    suptop = suptops[0]

    # the number of overlapped atoms is two
    assert len(suptop) == 3
    correct_overlaps = [('C1', 'C11'), ('N1', 'N11'), ('O1', 'O11')]
    for atomName1, atomName2 in correct_overlaps:
        assert suptop.contains_atomNamePair(atomName1, atomName2)

    # no mirrors
    assert len(suptop.mirrors) == 0


def test_SimpleMultipleSolutions():
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
    top1_list = [c1, n1, o1, o2]

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
    top2_list = [c11, n11, o11, o12]

    # should be two topologies
    suptops = _superimpose_topologies(top1_list, top2_list)
    # there is one solution
    assert len(suptops) == 1

    correct_overlaps = [('C1', 'C11'), ('N1', 'N11'), ('O1', 'O11'), ('O2', 'O12')]
    for st in suptops:
        for atomName1, atomName2 in correct_overlaps:
            assert st.contains_atomNamePair(atomName1, atomName2)

    # fixme - add a test case for the superimposer function that makes use of _overlay,
    # this is to resolve multiple solutions such as the one here


def test_SimpleMultipleSolutions_mirrors():
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
    top1_list = [c1, n1, o1, o2]

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
    top2_list = [c11, n11, o11, o12]

    # should be two topologies
    suptops = _superimpose_topologies(top1_list, top2_list)
    # there is one solution
    assert len(suptops) == 1
    #
    assert len(suptops[0].mirrors) == 1

    correct_overlaps = [('C1', 'C11'), ('N1', 'N11'), ('O1', 'O11'), ('O2', 'O12')]
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
    c1 = AtomNode(name='C1', type='C')
    c1.set_position(x=1, y=1, z=0)
    c2 = AtomNode(name='C2', type='C')
    c2.set_position(x=1, y=2, z=0)
    c1.bindTo(c2, 'bondType1')
    top1_list = [c1, c2]

    # construct the LIGAND 2
    c11 = AtomNode(name='C11', type='C')
    c11.set_position(x=1, y=1, z=0)
    c12 = AtomNode(name='C12', type='C')
    c12.set_position(x=1, y=2, z=0)
    c11.bindTo(c12, 'bondType1')
    top2_list = [c11, c12]

    # should return a list with an empty sup_top
    suptops = _superimpose_topologies(top1_list, top2_list)
    assert len(suptops) == 1
    assert len(suptops[0]) == 2


def test_2sameAtoms_2Cs_symmetry_mirrors():
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
    top1_list = [c1, c2]

    # construct the LIGAND 2
    c11 = AtomNode(name='C11', type='C')
    c11.set_position(x=1, y=1, z=0)
    c12 = AtomNode(name='C12', type='C')
    c12.set_position(x=1, y=2, z=0)
    c11.bindTo(c12, 'bondType1')
    top2_list = [c11, c12]

    # should return a list with an empty sup_top
    suptops = _superimpose_topologies(top1_list, top2_list)
    assert len(suptops) == 1
    assert len(suptops[0]) == 2
    assert len(suptops[0].mirrors) == 1
    assert len(suptops[0].mirrors[0]) == 2


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
    top1_list = [c1, c2, c3]

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
    top2_list = [c11, c12, c13]

    suptops = _superimpose_topologies(top1_list, top2_list)
    assert len(suptops) == 1
    assert len(suptops[0]) == 3
    assert len(suptops[0].mirrors) > 2


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
    top1_list = [c1, c2, c3, c4, cl1, c5, c6, c10, c7, n1, c8, c9]

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
    top2_list = [cl11, c11, c12, c13, c14, c15, c16, c20, c17, n11, c18, c19]

    # we have to discriminate against this case somehow
    suptops = _superimpose_topologies(top1_list, top2_list)
    assert len(suptops) == 1
    assert len(suptops[0]) == 11
    # two largest solutions are found?
    print('hi')


# def cprofile_mcl1():
# import cProfile
# cProfile.run("test_mcl1_l12l35()")