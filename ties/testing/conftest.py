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
    top1_list = [c1, n1, o1, o2]

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
    c1 = Atom(name='C1', atom_type='C')
    c1.position = (1, 1, 0)
    c2 = Atom(name='C2', atom_type='C')
    c2.position = (1, 2, 0)
    c1.bind_to(c2, 'bondType1')
    top1_list = [c1, c2]

    # construct the LIGAND 2
    c11 = Atom(name='C11', atom_type='C')
    c11.position = (1, 1, 0)
    c12 = Atom(name='C12', atom_type='C')
    c12.position = (1, 2, 0)
    c11.bind_to(c12, 'bondType1')
    top2_list = [c11, c12]

    return top1_list, top2_list


@pytest.fixture
def dual_ring1():
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

    return [c1, c2, c3, c4, cl1, c5, c6, c10, c7, n1, c8, c9]


@pytest.fixture
def dual_ring2():
    """
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
    return [cl11, c11, c12, c13, c14, c15, c16, c20, c17, n11, c18, c19]
