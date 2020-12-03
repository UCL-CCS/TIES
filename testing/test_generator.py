"""
These tests focus on the generator (preprocessing of the input before applying superimpose_topologies
"""

import pytest
from ties.topology_superimposer import superimpose_topologies, AtomNode, SuperimposedTopology


def test_same_atomNames_renamed():
    """
    # fixme - this test should be moved to a spearate files for the generator

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
                 Cl1
                /
         C1 - C2
         /      \
        C3      C4
          \     /
          C5 - C6
          /     \
     C10-C7       N1
           \   /
             C18
             |
             C19
    """
    # construct LIGAND 1
    c1 = AtomNode(name='C1', atom_type='C')
    c1.set_position(x=1, y=1, z=0)
    c2 = AtomNode(name='C2', atom_type='C')
    c2.set_position(x=1, y=2, z=0)
    c1.bind_to(c2, 'bondType1')
    c3 = AtomNode(name='C3', atom_type='C')
    c3.set_position(x=2, y=2, z=0)
    c3.bind_to(c1, 'bondType1')
    cl1 = AtomNode(name='CL1', atom_type='Cl')
    cl1.set_position(x=2, y=1, z=0)
    cl1.bind_to(c3, 'bondType1')
    c4 = AtomNode(name='C4', atom_type='C')
    c4.set_position(x=2, y=3, z=0)
    c4.bind_to(c2, 'bondType1')
    c5 = AtomNode(name='C5', atom_type='C')
    c5.set_position(x=3, y=1, z=0)
    c5.bind_to(c3, 'bondType1')
    c6 = AtomNode(name='C6', atom_type='C')
    c6.set_position(x=3, y=2, z=0)
    c6.bind_to(c5, 'bondType1')
    c6.bind_to(c4, 'bondType1')
    c7 = AtomNode(name='C7', atom_type='C')
    c7.set_position(x=4, y=2, z=0)
    c7.bind_to(c5, 'bondType1')
    c10 = AtomNode(name='C10', atom_type='C')
    c10.set_position(x=4, y=1, z=0)
    c10.bind_to(c7, 'bondType1')
    n1 = AtomNode(name='N1', atom_type='N')
    n1.set_position(x=4, y=3, z=0)
    n1.bind_to(c6, 'bondType1')
    c8 = AtomNode(name='C8', atom_type='C')
    c8.set_position(x=5, y=1, z=0)
    c8.bind_to(c7, 'bondType1')
    c8.bind_to(n1, 'bondType1')
    c9 = AtomNode(name='C9', atom_type='C')
    c9.set_position(x=6, y=1, z=0)
    c9.bind_to(c8, 'bondType1')
    top1_list = [c1, c2, c3, c4, cl1, c5, c6, c10, c7, n1, c8, c9]

    # construct Ligand 2
    cl11 = AtomNode(name='Cl1', atom_type='Cl')
    cl11.set_position(x=1, y=1, z=0)
    c11 = AtomNode(name='C1', atom_type='C')
    c11.set_position(x=2, y=1, z=0)
    c12 = AtomNode(name='C2', atom_type='C')
    c12.set_position(x=2, y=2, z=0)
    c12.bind_to(c11, 'bondType1')
    c12.bind_to(cl11, 'bondType1')
    c13 = AtomNode(name='C3', atom_type='C')
    c13.set_position(x=3, y=1, z=0)
    c13.bind_to(c11, 'bondType1')
    c14 = AtomNode(name='C4', atom_type='C')
    c14.set_position(x=3, y=2, z=0)
    c14.bind_to(c12, 'bondType1')
    c15 = AtomNode(name='C5', atom_type='C')
    c15.set_position(x=4, y=1, z=0)
    c15.bind_to(c13, 'bondType1')
    c16 = AtomNode(name='C6', atom_type='C')
    c16.set_position(x=4, y=2, z=0)
    c16.bind_to(c15, 'bondType1')
    c16.bind_to(c14, 'bondType1')
    c17 = AtomNode(name='C7', atom_type='C')
    c17.set_position(x=5, y=2, z=0)
    c17.bind_to(c15, 'bondType1')
    c20 = AtomNode(name='C10', atom_type='C')
    c20.set_position(x=5, y=1, z=0)
    c20.bind_to(c17, 'bondType1')
    n11 = AtomNode(name='N1', atom_type='N')
    n11.set_position(x=5, y=3, z=0)
    n11.bind_to(c16, 'bondType1')
    c18 = AtomNode(name='C18', atom_type='C')
    c18.set_position(x=6, y=1, z=0)
    c18.bind_to(c17, 'bondType1')
    c18.bind_to(n11, 'bondType1')
    c19 = AtomNode(name='C19', atom_type='C')
    c19.set_position(x=7, y=1, z=0)
    c19.bind_to(c18, 'bondType1')
    top2_list = [cl11, c11, c12, c13, c14, c15, c16, c20, c17, n11, c18, c19]

    # we have to discriminate against this case somehow
    # suptops = superimpose_topologies(top1_list, top2_list)

    SuperimposedTopology.rename_ligands(top1_list, top2_list)
    same_names = []
    for n1 in top1_list:
        for n2 in top2_list:
            if n1.name == n2.name:
                same_names.append((n1, n2))
    assert len(same_names) == 0, same_names
