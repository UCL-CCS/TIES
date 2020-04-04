"""
Focuses on the _superimpose_topology function that
superimposes the molecule in many different and
then processes the outputs to ensure the best match is found.

TODO
- this testing module should focus on "separated" molecules
"""

from topology_superimposer import superimpose_topologies, AtomNode

def test_unconnected_component_removed():
    """
    The charges correction removed the C8-C18 pair, which leaves C9-C19 as a disconnected component.
    This should be removed

    Ligand 1

         C1 - C2
         /      \
    Cl1-C3      C4
          \     /
          C5 - C6
          /     \
     C10-C7       N1
           \   /
             C8 (CHARGE)
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
             C18 (CHARGE)
             |
             C19
    """
    # construct LIGAND 1
    c1 = AtomNode(name='C1', type='C', charge=0)
    c1.set_position(x=1, y=1, z=0)
    c2 = AtomNode(name='C2', type='C', charge=0)
    c2.set_position(x=1, y=2, z=0)
    c1.bindTo(c2, 'bondType1')
    c3 = AtomNode(name='C3', type='C', charge=0)
    c3.set_position(x=2, y=2, z=0)
    c3.bindTo(c1, 'bondType1')
    cl1 = AtomNode(name='CL1', type='Cl', charge=0)
    cl1.set_position(x=2, y=1, z=0)
    cl1.bindTo(c3, 'bondType1')
    c4 = AtomNode(name='C4', type='C', charge=0)
    c4.set_position(x=2, y=3, z=0)
    c4.bindTo(c2, 'bondType1')
    c5 = AtomNode(name='C5', type='C', charge=0)
    c5.set_position(x=3, y=1, z=0)
    c5.bindTo(c3, 'bondType1')
    c6 = AtomNode(name='C6', type='C', charge=0)
    c6.set_position(x=3, y=2, z=0)
    c6.bindTo(c5, 'bondType1')
    c6.bindTo(c4, 'bondType1')
    c7 = AtomNode(name='C7', type='C', charge=0)
    c7.set_position(x=4, y=2, z=0)
    c7.bindTo(c5, 'bondType1')
    c10 = AtomNode(name='C10', type='C', charge=0)
    c10.set_position(x=4, y=1, z=0)
    c10.bindTo(c7, 'bondType1')
    n1 = AtomNode(name='N1', type='N', charge=0)
    n1.set_position(x=4, y=3, z=0)
    n1.bindTo(c6, 'bondType1')
    c8 = AtomNode(name='C8', type='C', charge=0)
    c8.set_position(x=5, y=1, z=0)
    c8.bindTo(c7, 'bondType1')
    c8.bindTo(n1, 'bondType1')
    c9 = AtomNode(name='C9', type='C', charge=0)
    c9.set_position(x=6, y=1, z=0)
    c9.bindTo(c8, 'bondType1')
    top1_list = [c1, c2, c3, c4, cl1, c5, c6, c10, c7, n1, c8, c9]

    # construct Ligand 2
    cl11 = AtomNode(name='Cl11', type='Cl', charge=0)
    cl11.set_position(x=1, y=1, z=0)
    c11 = AtomNode(name='C11', type='C', charge=0)
    c11.set_position(x=2, y=1, z=0)
    c12 = AtomNode(name='C12', type='C', charge=0)
    c12.set_position(x=2, y=2, z=0)
    c12.bindTo(c11, 'bondType1')
    c12.bindTo(cl11, 'bondType1')
    c13 = AtomNode(name='C13', type='C', charge=0)
    c13.set_position(x=3, y=1, z=0)
    c13.bindTo(c11, 'bondType1')
    c14 = AtomNode(name='C14', type='C', charge=0)
    c14.set_position(x=3, y=2, z=0)
    c14.bindTo(c12, 'bondType1')
    c15 = AtomNode(name='C15', type='C', charge=0)
    c15.set_position(x=4, y=1, z=0)
    c15.bindTo(c13, 'bondType1')
    c16 = AtomNode(name='C16', type='C', charge=0)
    c16.set_position(x=4, y=2, z=0)
    c16.bindTo(c15, 'bondType1')
    c16.bindTo(c14, 'bondType1')
    c17 = AtomNode(name='C17', type='C', charge=0)
    c17.set_position(x=5, y=2, z=0)
    c17.bindTo(c15, 'bondType1')
    c20 = AtomNode(name='C20', type='C', charge=0)
    c20.set_position(x=5, y=1, z=0)
    c20.bindTo(c17, 'bondType1')
    n11 = AtomNode(name='N11', type='N', charge=0)
    n11.set_position(x=5, y=3, z=0)
    n11.bindTo(c16, 'bondType1')
    c18 = AtomNode(name='C18', type='C', charge=0.11)
    c18.set_position(x=6, y=1, z=0)
    c18.bindTo(c17, 'bondType1')
    c18.bindTo(n11, 'bondType1')
    c19 = AtomNode(name='C19', type='C', charge=0)
    c19.set_position(x=7, y=1, z=0)
    c19.bindTo(c18, 'bondType1')
    top2_list = [cl11, c11, c12, c13, c14, c15, c16, c20, c17, n11, c18, c19]

    # we have to discriminate against this case somehow
    suptops = superimpose_topologies(top1_list, top2_list)
    assert not suptops[0].contains_atomNamePair('C8', 'C18')
    assert not suptops[0].contains_atomNamePair('C9', 'C19')
