"""
Focuses on the _overlay function that actually traverses the molecule using the given starting points.

TODO:
- add more complex test case with multiple solutions (2 levels)  similar to test_SimpleMultipleSolutions_rightStart
- add a test case that shows that molecules cannot be fully superimposed when it is broken in half
- add more test cases that take into account the "differences"
"""

import copy

from ties.topology_superimposer import Atom, _overlay, SuperimposedTopology


def test_2diffAtoms_CN_wrongStart(CN):
    """
    CN mapped to CN with the wrong starting point does not create any solutions
    """
    CN2 = copy.deepcopy(CN)

    c1, n1 = CN
    c11, n11 = CN2

    suptop = _overlay(
        c1,
        n11,
        parent_n1=None,
        parent_n2=None,
        bond_types=(None, None),
        suptop=SuperimposedTopology(CN, CN2),
    )
    assert suptop is None


def test_2diffAtoms_CN_rightStart(CN):
    """
    CN mapped to CN when exploring from the right starting point.
    """
    CN2 = copy.deepcopy(CN)

    c1, _ = CN
    c11, _ = CN2

    suptop = _overlay(
        c1,
        c11,
        parent_n1=None,
        parent_n2=None,
        bond_types=(None, None),
        suptop=SuperimposedTopology(CN, CN2),
    )

    assert len(suptop) == 2
    suptop.contains_atom_name_pair("C1", "C1")
    suptop.contains_atom_name_pair("N1", "N1")

    # there is no other ways to traverse the molecule
    assert len(suptop.mirrors) == 0


def test_3diffAtoms_CNO_rightStart(CNO):
    """
    Only one solution exists.

     LIGAND 1
        C1
        \
        N1
        /
        O1
    """
    # construct the LIGAND 2
    CNO2 = copy.deepcopy(CNO)
    c1, n1, o1 = CNO
    c11, n11, o11 = CNO2

    # should overlap 2 atoms
    suptop = _overlay(
        c1,
        c11,
        parent_n1=None,
        parent_n2=None,
        bond_types=(None, None),
        suptop=SuperimposedTopology(CNO, CNO2),
    )

    # the number of overlapped atoms is two
    assert len(suptop) == 3
    correct_overlaps = [("C1", "C1"), ("N1", "N1"), ("O1", "O1")]
    for atomName1, atomName2 in correct_overlaps:
        assert suptop.contains_atom_name_pair(atomName1, atomName2)

    # no mirrors
    assert len(suptop.mirrors) == 0


def test_SimpleMultipleSolutions_rightStart(CNO_O):
    """
    The ester allows for mapping (O1-O11, O2-O12) and (O1-O12, O2-O11)
    """
    CNO_O2 = copy.deepcopy(CNO_O)
    [c1, n1, o1, o2] = CNO_O
    [c11, n11, o11, o12] = CNO_O2

    # should generate one topology which has one "symmetry"
    suptop = _overlay(
        c1,
        c11,
        parent_n1=None,
        parent_n2=None,
        bond_types=(None, None),
        suptop=SuperimposedTopology(CNO_O, CNO_O2),
    )
    assert len(suptop) == 4

    assert len(suptop.mirrors) == 1
    worse_st = suptop.mirrors[0]

    # check if both representations were found
    # The ester allows for mapping (O1-O11, O2-O12) and (O1-O12, O2-O11)
    assert suptop.contains_atom_name_pair(
        "O1", "O1"
    ) and suptop.contains_atom_name_pair("O2", "O2")
    assert worse_st.contains_atom_name_pair(
        "O1", "O2"
    ) and worse_st.contains_atom_name_pair("O2", "O1")

    assert suptop.contains_atom_name_pair("C1", "C1")
    assert suptop.contains_atom_name_pair("N1", "N1")


def test_One2Many(CN, CN_N):
    """
        C1              C11
        /              /   \
     N1              N11   N12

    """
    c1, n1 = CN
    c11, n11, n12 = CN_N

    n1.position = (2, 1, 0)
    n11.position = (2, 1, 0)
    n12.position = (2, 2, 0)

    # should generate one topology which has one "symmetry"
    suptop = _overlay(
        c1,
        c11,
        parent_n1=None,
        parent_n2=None,
        bond_types=(None, None),
        suptop=SuperimposedTopology(CN, CN_N),
    )
    assert len(suptop) == 2

    # the better solution uses coordinates to figure out which match is better
    assert suptop.contains_atom_name_pair("C1", "C1")
    assert suptop.contains_atom_name_pair("N1", "N1")


def test_Many2One_part1(CN_N, CN):
    """
    LIGAND 1        LIGAND 2
       C1              C1
       / \              /
    N1    N2          N1
    """
    c1, n1, n2 = CN_N
    c11, n11 = CN

    n1.position = (2, 1, 0)
    n11.position = (2, 1, 0)
    n2.position = (2, 2, 0)

    # should generate one topology which has one "symmetry"
    suptop = _overlay(
        n1,
        n11,
        parent_n1=None,
        parent_n2=None,
        bond_types=(None, None),
        suptop=SuperimposedTopology(CN_N, CN),
    )
    assert len(suptop) == 2

    assert suptop.contains_atom_name_pair("C1", "C1")
    assert suptop.contains_atom_name_pair("N1", "N1")


def test_Many2One_part2(CN_N, CN):
    """
     LIGAND 1        LIGAND 2
         C1             C1
        / \              \
     N1    N2             N1
    """
    c1, n1, n2 = CN_N
    c11, n12 = CN

    n1.position = (2, 1, 0)
    n2.position = (2, 2, 0)
    n12.position = (2, 2, 0)

    # should generate one topology which has one "symmetry"
    suptop = _overlay(
        c1,
        c11,
        parent_n1=None,
        parent_n2=None,
        bond_types=(None, None),
        suptop=SuperimposedTopology(CN_N, CN),
    )
    assert len(suptop) == 2

    assert suptop.contains_atom_name_pair("C1", "C1")
    assert suptop.contains_atom_name_pair("N2", "N1")


def test_MultipleSolutions2Levels_rightStart(CN_N):
    """
         LIGAND 1        LIGAND 2
            C1              C11
            /\              / \
         N1    N2        N11   N12
            / \    / \      /  \   /  \
       O1 O2  O3 O4   O11 O12 O13 O14
    """
    c1, n1, n2 = CN_N

    o1 = Atom(name="O1", atom_type="O")
    o1.position = (3, 1, 0)
    o1.bind_to(n1, "bondType1")
    o2 = Atom(name="O2", atom_type="O")
    o2.position = (3, 2, 0)
    o2.bind_to(n1, "bondType1")
    o3 = Atom(name="O3", atom_type="O")
    o3.position = (3, 3, 0)
    o3.bind_to(n2, "bondType1")
    o4 = Atom(name="O4", atom_type="O")
    o4.position = (3, 4, 0)
    o4.bind_to(n2, "bondType1")
    pyramid_left = [c1, n1, n2, o1, o2, o3, o4]

    pyramid_right = copy.deepcopy(pyramid_left)

    suptop = _overlay(
        c1,
        pyramid_right[0],
        parent_n1=None,
        parent_n2=None,
        bond_types=(None, None),
        suptop=SuperimposedTopology(pyramid_left, pyramid_right),
    )

    correct_overlaps = [
        ("C1", "C1"),
        ("N1", "N1"),
        ("N2", "N2"),
        ("O1", "O1"),
        ("O2", "O2"),
        ("O3", "O3"),
        ("O4", "O4"),
    ]
    for atomName1, atomName2 in correct_overlaps:
        assert suptop.contains_atom_name_pair(atomName1, atomName2)


def test_2sameAtoms_2Cs_symmetry(CC):
    """
    Two solutions with different starting points.

     LIGAND 1        LIGAND 2
        C1              C11
        \                \
        C2              C12
    """
    # construct the LIGAND 1
    c1, c2 = CC
    c11, c12 = copy.deepcopy(CC)

    # should return a list with an empty sup_top
    suptop = _overlay(
        c1,
        c11,
        parent_n1=None,
        parent_n2=None,
        bond_types=(None, None),
        suptop=SuperimposedTopology(CC, [c11, c12]),
    )
    assert len(suptop) == 2
    assert suptop.contains_atom_name_pair("C1", "C1")
    assert suptop.contains_atom_name_pair("C2", "C2")

    suptop = _overlay(
        c1,
        c12,
        parent_n1=None,
        parent_n2=None,
        bond_types=(None, None),
        suptop=SuperimposedTopology(CC, [c11, c12]),
    )
    assert len(suptop) == 2
    assert suptop.contains_atom_name_pair("C1", "C2")
    assert suptop.contains_atom_name_pair("C2", "C1")


def test_mutation_separate_unique_match(CNO):
    """
    Separated by the mutation.

     LIGAND 1        LIGAND 2
        C1              C11
        |                |
        N1              S11
        |                |
        O1               O11
    """
    CSO = copy.deepcopy(CNO)
    c1, n1, o1 = CNO
    c11, s11, o11 = CSO

    s11.element = "S"
    s11.name = "S11"

    # should return a list with an empty sup_top
    suptop = _overlay(
        c1,
        c11,
        parent_n1=None,
        parent_n2=None,
        bond_types=(None, None),
        suptop=SuperimposedTopology(CNO, CSO),
    )
    assert len(suptop) == 1
    assert suptop.contains_atom_name_pair("C1", "C1")

    suptop = _overlay(
        o1,
        o11,
        parent_n1=None,
        parent_n2=None,
        bond_types=(None, None),
        suptop=SuperimposedTopology(CNO, CSO),
    )
    assert len(suptop) == 1
    assert suptop.contains_atom_name_pair("O1", "O1")


def test_3C_circle(CCC):
    """
    A circle should be detected.
    Many solutions (3 starting conditions, etc)

      LIGAND 1        LIGAND 2
         C1              C11
        /  \            /   \
      C2 - C3         C12 - C13
    """
    c1, c2, c3 = CCC
    c2.bind_to(c3, "bondType1")
    CCC2 = copy.deepcopy([c1, c2, c3])
    c11, c12, c13 = CCC2

    suptop = _overlay(
        c1,
        c11,
        parent_n1=None,
        parent_n2=None,
        bond_types=(None, None),
        suptop=SuperimposedTopology(CCC, CCC2),
    )
    # there is one symmetrical way to traverse it
    assert len(suptop.mirrors) == 1

    assert suptop.contains_atom_name_pair("C1", "C1")
    assert suptop.contains_atom_name_pair("C2", "C2")
    assert suptop.contains_atom_name_pair("C3", "C3")

    # there should be one circle in each
    assert suptop.same_circle_number()
    assert suptop.get_circle_number() == (1, 1)

    suptop = _overlay(
        c1,
        c12,
        parent_n1=None,
        parent_n2=None,
        bond_types=(None, None),
        suptop=SuperimposedTopology(CCC, CCC2),
    )
    assert suptop.same_circle_number()
    assert suptop.get_circle_number() == (1, 1)

    suptop = _overlay(
        c1,
        c13,
        parent_n1=None,
        parent_n2=None,
        bond_types=(None, None),
        suptop=SuperimposedTopology(CCC, CCC2),
    )
    assert suptop.same_circle_number()
    assert suptop.get_circle_number() == (1, 1)

    suptop = _overlay(
        c2,
        c11,
        parent_n1=None,
        parent_n2=None,
        bond_types=(None, None),
        suptop=SuperimposedTopology(CCC, CCC2),
    )
    assert suptop.same_circle_number()
    assert suptop.get_circle_number() == (1, 1)

    suptop = _overlay(
        c2,
        c12,
        parent_n1=None,
        parent_n2=None,
        bond_types=(None, None),
        suptop=SuperimposedTopology(CCC, CCC2),
    )
    assert suptop.same_circle_number()
    assert suptop.get_circle_number() == (1, 1)

    suptop = _overlay(
        c2,
        c13,
        parent_n1=None,
        parent_n2=None,
        bond_types=(None, None),
        suptop=SuperimposedTopology(CCC, CCC2),
    )
    assert suptop.same_circle_number()
    assert suptop.get_circle_number() == (1, 1)

    suptop = _overlay(
        c3,
        c11,
        parent_n1=None,
        parent_n2=None,
        bond_types=(None, None),
        suptop=SuperimposedTopology(CCC, CCC2),
    )
    assert suptop.same_circle_number()
    assert suptop.get_circle_number() == (1, 1)

    suptop = _overlay(
        c3,
        c12,
        parent_n1=None,
        parent_n2=None,
        bond_types=(None, None),
        suptop=SuperimposedTopology(CCC, CCC2),
    )
    assert suptop.same_circle_number()
    assert suptop.get_circle_number() == (1, 1)

    suptop = _overlay(
        c3,
        c13,
        parent_n1=None,
        parent_n2=None,
        bond_types=(None, None),
        suptop=SuperimposedTopology(CCC, CCC2),
    )
    assert suptop.same_circle_number()
    assert suptop.get_circle_number() == (1, 1)


def test_3C_circle_bonds(CCC):
    """
    Check if each pair is linked to two other pairs
      LIGAND 1        LIGAND 2
         C1              C11
        /  \            /   \
      C2 - C3         C12 - C13
    """
    c1, c2, c3 = CCC
    c2.bind_to(c3, "bondType1")
    CCC2 = copy.deepcopy([c1, c2, c3])
    c11, c12, c13 = CCC2

    suptop = _overlay(
        c1,
        c11,
        parent_n1=None,
        parent_n2=None,
        bond_types=(None, None),
        suptop=SuperimposedTopology(CCC, CCC2),
    )

    # check that each pair is linked to two other pairs
    for pair, linked_pairs in suptop.matched_pairs_bonds.items():
        assert len(linked_pairs) == 2


def test_mcl1_l12l35(indole_cl1, indole_cl2):
    """
    Molecule inspired by mcl1_l12l35

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
    _, _, _, _, _, c5, _, _, _, _, _, c9, _, _ = indole_cl1
    _, _, _, c14, _, _, _, _, _, c19, _, _ = indole_cl2

    # the correct solution
    suptop = _overlay(
        c9,
        c19,
        parent_n1=None,
        parent_n2=None,
        bond_types=(None, None),
        suptop=SuperimposedTopology(indole_cl1, indole_cl2),
    )
    assert len(suptop) == 11

    """
    This is a rare case which uses the rules of "consistent circles".
    A naive traversal of this is possible that reaches 12 nodes. However, 
    that means that there is atom match (L1-R1) such that L1 creates a 
    cricle in its own topology, and R1 does not. 
    """
    suptop = _overlay(
        c5,
        c14,
        parent_n1=None,
        parent_n2=None,
        bond_types=(None, None),
        suptop=SuperimposedTopology(indole_cl1, indole_cl2),
    )
    assert len(suptop) != 12


def test_mcl1_l12l35_bonds(indole_cl1, indole_cl2):
    """
    Molecule inspired by mcl1_l12l35
    """
    c1, c2, c3, c4, cl1, c5, c6, c10, c7, n1, c8, c9, cl1, c10 = indole_cl1
    c11, c12, c13, c14, c15, c16, c17, n11, c18, c19, cl11, c20 = indole_cl2

    # the correct solution
    suptop = _overlay(
        c9,
        c19,
        parent_n1=None,
        parent_n2=None,
        bond_types=(None, None),
        suptop=SuperimposedTopology(indole_cl1, indole_cl2),
    )

    # (c1, c11) should be linked to (c2, c12), (c3, c13)
    assert {x[0] for x in suptop.matched_pairs_bonds[(c1, c11)]} == {
        (c2, c12),
        (c3, c13),
    }
    # (c2, c12) should be linked to (c1, c11), (c4, c14)
    assert {x[0] for x in suptop.matched_pairs_bonds[(c2, c12)]} == {
        (c1, c11),
        (c4, c14),
    }
    # etc
    assert {x[0] for x in suptop.matched_pairs_bonds[(c3, c13)]} == {
        (c1, c11),
        (c5, c15),
    }
    assert {x[0] for x in suptop.matched_pairs_bonds[(c4, c14)]} == {
        (c2, c12),
        (c6, c16),
    }
    assert {x[0] for x in suptop.matched_pairs_bonds[(c5, c15)]} == {
        (c3, c13),
        (c6, c16),
        (c7, c17),
    }
    assert {x[0] for x in suptop.matched_pairs_bonds[(c6, c16)]} == {
        (c5, c15),
        (c4, c14),
        (n1, n11),
    }
    assert {x[0] for x in suptop.matched_pairs_bonds[(c7, c17)]} == {
        (c5, c15),
        (c10, c20),
        (c8, c18),
    }
    assert {x[0] for x in suptop.matched_pairs_bonds[(c8, c18)]} == {
        (c7, c17),
        (n1, n11),
        (c9, c19),
    }
    assert {x[0] for x in suptop.matched_pairs_bonds[(c9, c19)]} == {(c8, c18)}
    assert {x[0] for x in suptop.matched_pairs_bonds[(c10, c20)]} == {(c7, c17)}
    assert {x[0] for x in suptop.matched_pairs_bonds[(n1, n11)]} == {
        (c6, c16),
        (c8, c18),
    }


def test_tyk2_l11l14_part():
    """

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
    n1 = Atom(name="N1", atom_type="N")
    n1.position = (1, 1, 0)
    c1 = Atom(name="C1", atom_type="C")
    c1.position = (2, 1, 0)
    c1.bind_to(n1, "bondType1")
    o1 = Atom(name="O1", atom_type="O")
    o1.position = (2, 2, 0)
    o1.bind_to(c1, "bondType1")
    h1 = Atom(name="H1", atom_type="H")
    h1.position = (3, 1, 0)
    c2 = Atom(name="C2", atom_type="C")
    c2.position = (3, 2, 0)
    c2.bind_to(h1, "bondType1")
    c2.bind_to(c1, "bondType1")
    h2 = Atom(name="H2", atom_type="H")
    h2.position = (4, 1, 0)
    c3 = Atom(name="C3", atom_type="C")
    c3.position = (4, 2, 0)
    c3.bind_to(c2, "bondType1")
    c3.bind_to(h2, "bondType1")
    c4 = Atom(name="C4", atom_type="C")
    c4.position = (4, 3, 0)
    c4.bind_to(c2, "bondType1")
    c4.bind_to(c3, "bondType1")
    h3 = Atom(name="H3", atom_type="H")
    h3.position = (4, 4, 0)
    h3.bind_to(c4, "bondType1")
    f1 = Atom(name="F1", atom_type="F")
    f1.position = (5, 1, 0)
    f1.bind_to(c3, "bondType1")
    h4 = Atom(name="H4", atom_type="H")
    h4.position = (5, 2, 0)
    h4.bind_to(c4, "bondType1")
    left_atoms = [n1, c1, o1, h1, c2, h2, c3, c4, h3, f1, h4]

    # construct Ligand 2
    n11 = Atom(name="N11", atom_type="N")
    n11.position = (1, 1, 0)
    c11 = Atom(name="C11", atom_type="C")
    c11.position = (2, 1, 0)
    c11.bind_to(n11, "bondType1")
    o11 = Atom(name="O11", atom_type="O")
    o11.position = (2, 2, 0)
    o11.bind_to(c11, "bondType1")
    h11 = Atom(name="H11", atom_type="H")
    h11.position = (3, 1, 0)
    c12 = Atom(name="C12", atom_type="C")
    c12.position = (3, 2, 0)
    c12.bind_to(h11, "bondType1")
    c12.bind_to(c11, "bondType1")
    h12 = Atom(name="H12", atom_type="H")
    h12.position = (4, 1, 0)
    c13 = Atom(name="C13", atom_type="C")
    c13.position = (4, 2, 0)
    c13.bind_to(c12, "bondType1")
    c13.bind_to(h12, "bondType1")
    c14 = Atom(name="C14", atom_type="C")
    c14.position = (4, 3, 0)
    c14.bind_to(c12, "bondType1")
    c14.bind_to(c13, "bondType1")
    h13 = Atom(name="H13", atom_type="H")
    h13.position = (4, 4, 0)
    h13.bind_to(c14, "bondType1")
    cl11 = Atom(name="CL11", atom_type="CL")
    cl11.position = (5, 1, 0)
    cl11.bind_to(c13, "bondType1")
    h14 = Atom(name="H14", atom_type="H")
    h14.position = (5, 2, 0)
    h14.bind_to(c14, "bondType1")
    right_atoms = [n11, c11, o11, h11, c12, h12, c13, c14, h13, cl11, h14]

    # the correct solution
    suptop = _overlay(
        n1,
        n11,
        parent_n1=None,
        parent_n2=None,
        bond_types=(None, None),
        suptop=SuperimposedTopology(left_atoms, right_atoms),
    )
    assert len(suptop) == 10

    matching_pairs = [
        ("N1", "N11"),
        ("C1", "C11"),
        ("O1", "O11"),
        ("H1", "H11"),
        ("C2", "C12"),
        ("H2", "H12"),
        ("C3", "C13"),
        ("C4", "C14"),
        ("H3", "H13"),
        ("H4", "H14"),
    ]
    for n1, n2 in matching_pairs[::-1]:
        if suptop.contains_atom_name_pair(n1, n2):
            matching_pairs.remove((n1, n2))
    assert len(matching_pairs) == 0, matching_pairs
