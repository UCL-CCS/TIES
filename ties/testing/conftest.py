import copy

import pytest
from ties.topology_superimposer import Atom


@pytest.fixture
def lr_atoms_branches():
    """
     LIGAND 1        LIGAND 2
        C1              C11
        \                \
        N1              N11
        /\              / \
     O1    O2        O11   O12
    """
    # ignore the third coordinate dimension
    # construct the LIGAND 1
    c1 = Atom(name="C1", atom_type="C")
    c1.position = (1, 1, 0)
    n1 = Atom(name="N1", atom_type="N")
    n1.position = (1, 2, 0)
    c1.bind_to(n1, "bondType1")
    o1 = Atom(name="O1", atom_type="O")
    o1.position = (1, 3, 0)
    o1.bind_to(n1, "bondType1")
    o2 = Atom(name="O2", atom_type="O")
    o2.position = (2, 3, 0)
    o2.bind_to(n1, "bondType1")
    top1_list = [c1, n1, o1, o2]

    # construct the LIGAND 2
    c11 = Atom(name="C11", atom_type="C")
    c11.position = (1, 1, 0)
    n11 = Atom(name="N11", atom_type="N")
    n11.position = (1, 2, 0)
    c11.bind_to(n11, "bondType1")
    o11 = Atom(name="O11", atom_type="O")
    o11.position = (1, 3, 0)
    o11.bind_to(n11, "bondType1")
    o12 = Atom(name="O12", atom_type="O")
    o12.position = (2, 3, 0)
    o12.bind_to(n11, "bondType1")
    top2_list = [c11, n11, o11, o12]

    return top1_list, top2_list


@pytest.fixture
def lr_atoms_chain_2c():
    """
      LIGAND 1        LIGAND 2
        C1              C11
        \                \
        C2              C12
    """
    # construct the LIGAND 1
    c1 = Atom(name="C1", atom_type="C")
    c1.position = (1, 1, 0)
    c2 = Atom(name="C2", atom_type="C")
    c2.position = (1, 2, 0)
    c1.bind_to(c2, "bondType1")
    top1_list = [c1, c2]

    # construct the LIGAND 2
    c11 = Atom(name="C11", atom_type="C")
    c11.position = (1, 1, 0)
    c12 = Atom(name="C12", atom_type="C")
    c12.position = (1, 2, 0)
    c11.bind_to(c12, "bondType1")
    top2_list = [c11, c12]

    return top1_list, top2_list


@pytest.fixture
def C():
    c1 = Atom(name="C1", atom_type="C")
    c1.position = (1, 1, 0)
    return [c1]


@pytest.fixture
def CC(C):
    c = copy.deepcopy(C)[0]
    c2 = Atom(name="C2", atom_type="C")
    c2.position = (2, 1, 0)
    c2.bind_to(c, "bondType1")
    return c, c2


@pytest.fixture
def CCC(CC):
    c1, c2 = copy.deepcopy(CC)
    c3 = Atom(name="C3", atom_type="C")
    c3.position = (2, 2, 0)
    c3.bind_to(c1, "bondType1")
    return c1, c2, c3


@pytest.fixture
def CN():
    """
    create a simple CN
     LIGAND 1
        C1
        \
        N1
    """
    c1 = Atom(name="C1", atom_type="C")
    c1.position = (1, 1, 0)
    n1 = Atom(name="N1", atom_type="N")
    n1.position = (1, 2, 0)
    c1.bind_to(n1, "bondType1")
    return [c1, n1]


@pytest.fixture
def CN_N(CN):
    c, n = copy.deepcopy(CN)
    n2 = Atom(name="N2", atom_type="N")
    n2.position = (2, 2, 0)
    n2.bind_to(c, "bondType1")
    return [c, n, n2]


@pytest.fixture
def CNO(CN):
    """
        C1
        \
        N1
        /
        O1
    """
    c, n = copy.deepcopy(CN)
    o = Atom(name="O1", atom_type="O")
    o.position = (1, 3, 0)
    o.bind_to(n, "bondType1")
    return [c, n, o]


@pytest.fixture
def CNO_O(CNO):
    """
     LIGAND 1
        C1
        \
        N1
        /\
     O1    O2
    """
    c, n, o = copy.deepcopy(CNO)
    o2 = Atom(name="O2", atom_type="O")
    o2.position = (2, 3, 0)
    o2.bind_to(n, "bondType1")
    return [c, n, o, o2]


@pytest.fixture
def indole_simple():
    """
         Ligand 1

         C1 - C2
         /      \
        C3      C4
          \     /
          C5 - C6
          /     \
         C7       N1
           \   /
             C8
             |
             C9
    """
    c1 = Atom(name="C1", atom_type="C")
    c1.position = (1, 1, 0)
    c2 = Atom(name="C2", atom_type="C")
    c2.position = (1, 2, 0)
    c1.bind_to(c2, "bondType1")
    c3 = Atom(name="C3", atom_type="C")
    c3.position = (2, 2, 0)
    c3.bind_to(c1, "bondType1")
    c4 = Atom(name="C4", atom_type="C")
    c4.position = (2, 3, 0)
    c4.bind_to(c2, "bondType1")
    c5 = Atom(name="C5", atom_type="C")
    c5.position = (3, 1, 0)
    c5.bind_to(c3, "bondType1")
    c6 = Atom(name="C6", atom_type="C")
    c6.position = (3, 2, 0)
    c6.bind_to(c5, "bondType1")
    c6.bind_to(c4, "bondType1")
    c7 = Atom(name="C7", atom_type="C")
    c7.position = (4, 2, 0)
    c7.bind_to(c5, "bondType1")
    n1 = Atom(name="N1", atom_type="N")
    n1.position = (4, 3, 0)
    n1.bind_to(c6, "bondType1")
    c8 = Atom(name="C8", atom_type="C")
    c8.position = (5, 1, 0)
    c8.bind_to(c7, "bondType1")
    c8.bind_to(n1, "bondType1")
    c9 = Atom(name="C9", atom_type="C")
    c9.position = (6, 1, 0)
    c9.bind_to(c8, "bondType1")

    return [c1, c2, c3, c4, c5, c6, c7, n1, c8, c9]


@pytest.fixture
def indole_cl1(indole_simple):
    """
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
    """
    c1, c2, c3, c4, c5, c6, c7, n1, c8, c9 = copy.deepcopy(indole_simple)
    cl1 = Atom(name="CL1", atom_type="Cl")
    cl1.position = (2, 1, 0)
    cl1.bind_to(c3, "bondType1")
    c10 = Atom(name="C10", atom_type="C")
    c10.position = (4, 1, 0)
    c10.bind_to(c7, "bondType1")

    return c1, c2, c3, c4, cl1, c5, c6, c10, c7, n1, c8, c9, cl1, c10


@pytest.fixture
def indole_cl2(indole_simple):
    """
        Ligand 2
                 Cl1
                /
         C1 - C2
         /      \
        C3      C4
          \     /
          C5 - C6
          /     \
     C10-C7      N1
           \   /
            C8
             |
             C9
    """
    c1, c2, c3, c4, c5, c6, c7, n1, c8, c9 = copy.deepcopy(indole_simple)

    cl1 = Atom(name="Cl11", atom_type="Cl")
    cl1.position = (1, 1, 0)
    c10 = Atom(name="C10", atom_type="C")
    c10.position = (5, 1, 0)
    c10.bind_to(c7, "bondType1")
    return c1, c2, c3, c4, c5, c6, c7, n1, c8, c9, cl1, c10
