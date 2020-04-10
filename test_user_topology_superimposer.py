"""
Focuses on the _superimpose_topology function that
superimposes the molecule in many different and
then processes the outputs to ensure the best match is found.

TODO
- this testing module should focus on "separated" molecules
"""

from topology_superimposer import superimpose_topologies, AtomNode
import numpy as np

def test_unconnected_component_removed():
    """
    The charges correction removed the C8-C18 and N1-N11 pairs,
    which leaves C9-C19 as a disconnected component, and
    because of that, C9-C19 should be removed

    Ligand 1

         C1 - C2
         /      \
    Cl1-C3      C4
          \     /
          C5 - C6
          /     \
     C10-C7       N1 (1e)
           \   /
             C8 (-1e)
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
     C20-C17       N11 (-1e)
           \   /
             C18 (1e)
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
    n1 = AtomNode(name='N1', type='N', charge=1)
    n1.set_position(x=4, y=3, z=0)
    n1.bindTo(c6, 'bondType1')
    c8 = AtomNode(name='C8', type='C', charge=-1)
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
    n11 = AtomNode(name='N11', type='N', charge=-1)
    n11.set_position(x=5, y=3, z=0)
    n11.bindTo(c16, 'bondType1')
    c18 = AtomNode(name='C18', type='C', charge=1)
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
    assert not suptops[0].contains_atomNamePair('N1', 'N11')
    assert not suptops[0].contains_atomNamePair('C9', 'C19')


def test_filter_net_charge_too_large():
    """
    After removing the pair that does not match C2-C12 due to charge (0.14e difference),
    there is a net charge difference larger than 0.1:
    1.4 in total: 0.8 from N1-N11 0.6 from C8-C18.
    Therefore, remove the pair that has the highest
    charge difference, which is N1-N11

    Ligand 1

         C1 - C2 (+0.07e)
         /      \
    Cl1-C3      C4
          \     /
          C5 - C6
          /     \
     C10-C7      N1 (-0.04e)
           \   /
             C8 (-0.03e)
             |
             C9


    Ligand 2
                 Cl11
                /
         C11 - C12 (-0.07e)
         /      \
        C13      C14
          \     /
          C15 - C16
          /     \
     C20-C17       N11 (+0.04e)
           \   /
             C18 (+0.03e)
             |
             C19
    """
    # construct LIGAND 1
    c1 = AtomNode(name='C1', type='C', charge=0)
    c1.set_position(x=1, y=1, z=0)
    c2 = AtomNode(name='C2', type='C', charge=0.07)
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
    n1 = AtomNode(name='N1', type='N', charge=-0.04)
    n1.set_position(x=4, y=3, z=0)
    n1.bindTo(c6, 'bondType1')
    c8 = AtomNode(name='C8', type='C', charge=-0.03)
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
    c12 = AtomNode(name='C12', type='C', charge=-0.07)
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
    n11 = AtomNode(name='N11', type='N', charge=0.04)
    n11.set_position(x=5, y=3, z=0)
    n11.bindTo(c16, 'bondType1')
    c18 = AtomNode(name='C18', type='C', charge=0.03)
    c18.set_position(x=6, y=1, z=0)
    c18.bindTo(c17, 'bondType1')
    c18.bindTo(n11, 'bondType1')
    c19 = AtomNode(name='C19', type='C', charge=0)
    c19.set_position(x=7, y=1, z=0)
    c19.bindTo(c18, 'bondType1')
    top2_list = [cl11, c11, c12, c13, c14, c15, c16, c20, c17, n11, c18, c19]

    # we have to discriminate against this case somehow
    suptops = superimpose_topologies(top1_list, top2_list)
    # removed because the charge difference was too large
    assert not suptops[0].contains_atomNamePair('C2', 'C12')
    # removed to bring the net charge difference below 0.1e
    assert not suptops[0].contains_atomNamePair('N1', 'N11')
    # this pair should remain as the net charge difference is below 0.1
    assert suptops[0].contains_atomNamePair('C8', 'C18')
    assert len(suptops[0]) == 9


def test_averaging_charges():
    """
    Averages charges in the matched region.
    This simple case does not lead to any imbalances.
    Ie after the correction, both sides are neutral.

    Ligand 1

         C1 - C2
         /      \
    Cl1-C3      C4
          \     /
          C5 - C6
          /     \
     C10-C7      N1 (+0.02e)
           \   /
             C8 (-0.02e)
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
     C20-C17       N11 (-0.02e)
           \   /
             C18 (+0.02e)
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
    n1 = AtomNode(name='N1', type='N', charge=0.02)
    n1.set_position(x=4, y=3, z=0)
    n1.bindTo(c6, 'bondType1')
    c8 = AtomNode(name='C8', type='C', charge=-0.02)
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
    c12 = AtomNode(name='C12', type='C', charge=-0)
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
    n11 = AtomNode(name='N11', type='N', charge=-0.02)
    n11.set_position(x=5, y=3, z=0)
    n11.bindTo(c16, 'bondType1')
    c18 = AtomNode(name='C18', type='C', charge=0.02)
    c18.set_position(x=6, y=1, z=0)
    c18.bindTo(c17, 'bondType1')
    c18.bindTo(n11, 'bondType1')
    c19 = AtomNode(name='C19', type='C', charge=0)
    c19.set_position(x=7, y=1, z=0)
    c19.bindTo(c18, 'bondType1')
    top2_list = [cl11, c11, c12, c13, c14, c15, c16, c20, c17, n11, c18, c19]

    suptops = superimpose_topologies(top1_list, top2_list)
    # we are averaging the charges of the atoms in the left and right,
    assert n1.charge == n11.charge == 0
    assert c8.charge == c18.charge == 0


def test_averaging_charges_imbalance_distribution_single():
    """
    Averages charges in the matched region.
    The mismatched region C9-C19 had charges -0.04 in L and 0.02 in R,
    which balances the matched N1-N11 with 0.04 in L, and -0.02 in R.
    After averaging N1-N11 to 0.01, L has -0.03 imbalance, and R has 0.03 imbalance.
    Thus in L, we have 0.01 in matched, and -0.04 in unmatched, and imbalance of 0.03 across just one atom (unmatched area).
    So we spread 0.03 over the one atom with -0.04 from unmatched area, giving -0.01.
    In R, we have now 0.01 in matched, and 0.02 in unmatched, but we need -0.01 in unmatched. Similarly,
    we do 0.02 + -0.03 (imbalance) which leads to -0.01.

    Q: is it possible to have a charge misbalance if the unmatched regions have the same net charges? No.
    Therefore, it is really the unmatched regions that introduce the imbalance.

    fixme In L, We could distribute -0.03 over all atoms in common area. But why we treat the mismatching atoms
    in such a special way? Why unmatching atoms charges are more important than matched atoms?
    Is there a better way to redistribute charges? How does RESP does it?

    Ligand 1

         C1 - C2
         /      \
        C3      C4
          \     /
          C5 - C6
          /     \
     C10-C7      N1 (+0.04e)
           \   /
             C8
             |
             X9 (-0.04e)


    Ligand 2


         C11 - C12
         /      \
        C13      C14
          \     /
          C15 - C16
          /     \
     C20-C17       N11 (-0.02e)
           \   /
             C18
             |
             Y19 (+0.02e)
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
    n1 = AtomNode(name='N1', type='N', charge=0.04)
    n1.set_position(x=4, y=3, z=0)
    n1.bindTo(c6, 'bondType1')
    c8 = AtomNode(name='C8', type='C', charge=0)
    c8.set_position(x=5, y=1, z=0)
    c8.bindTo(c7, 'bondType1')
    c8.bindTo(n1, 'bondType1')
    c9 = AtomNode(name='C9', type='X', charge=-0.04)
    c9.set_position(x=6, y=1, z=0)
    c9.bindTo(c8, 'bondType1')
    top1_list = [c1, c2, c3, c4, c5, c6, c10, c7, n1, c8, c9]

    # construct Ligand 2
    c11 = AtomNode(name='C11', type='C', charge=0)
    c11.set_position(x=2, y=1, z=0)
    c12 = AtomNode(name='C12', type='C', charge=-0)
    c12.set_position(x=2, y=2, z=0)
    c12.bindTo(c11, 'bondType1')
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
    n11 = AtomNode(name='N11', type='N', charge=-0.02)
    n11.set_position(x=5, y=3, z=0)
    n11.bindTo(c16, 'bondType1')
    c18 = AtomNode(name='C18', type='C', charge=0)
    c18.set_position(x=6, y=1, z=0)
    c18.bindTo(c17, 'bondType1')
    c18.bindTo(n11, 'bondType1')
    c19 = AtomNode(name='C19', type='Y', charge=0.02)
    c19.set_position(x=7, y=1, z=0)
    c19.bindTo(c18, 'bondType1')
    top2_list = [c11, c12, c13, c14, c15, c16, c20, c17, n11, c18, c19]

    suptops = superimpose_topologies(top1_list, top2_list)
    suptop = suptops[0]

    assert not suptop.contains_atomNamePair('C9', 'C19')
    # we are averaging the charges of the atoms in the left and right,
    print('charges, ', n1.charge, n11.charge)
    np.testing.assert_array_almost_equal([n1.charge, n11.charge], 0.01)
    np.testing.assert_almost_equal(c9.charge, -0.01)
    np.testing.assert_almost_equal(c19.charge,-0.01)

def test_averaging_charges_imbalance_distribution_multiple():
    """
    Averages charges in the matched region.
    The mismatched region C9-C19 had charges -0.04 in L and 0.02 in R,
    which balances the matched N1-N11 with 0.04 in L, and -0.02 in R.
    After averaging N1-N11 to 0.01, L has -0.03 imbalance, and R has 0.03 imbalance.

    Q: is it possible to have a charge misbalance if the unmatched regions have the same net charges? No.
    Therefore, it is really the unmatched regions that introduce the imbalance.

    In L, We could distribute -0.03 over all atoms in common area. But why we treat the mismatching atoms
    in such a special way? Why unmatching atoms charges are more important than matched atoms?


    So in this case, in L we have to remove the imbalance by spreading -0.03 across just one atom (unmatched area),
    with a similar situation for R. In R, the 0.03 imbalance has to be spread over a single atom.


    Note: if we were able to map mismatching atoms, we would know about C9-C19
    representing the same area. Would that knowledge affect the charge distribtion?

    Ligand 1

         C1 - C2
         /      \
    Cl1-C3      C4
          \     /
          C5 - C6
          /     \
     C10-C7      N1 (+0.04e)
           \   /
             C8
             |
             X9 (-0.04e)


    Ligand 2
                 Cl11
                /
         C11 - C12
         /      \
        C13      C14
          \     /
          C15 - C16
          /     \
     C20-C17       N11 (-0.02e)
           \   /
             C18
             |
             Y19 (+0.02e)
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
    n1 = AtomNode(name='N1', type='N', charge=0.04)
    n1.set_position(x=4, y=3, z=0)
    n1.bindTo(c6, 'bondType1')
    c8 = AtomNode(name='C8', type='C', charge=0)
    c8.set_position(x=5, y=1, z=0)
    c8.bindTo(c7, 'bondType1')
    c8.bindTo(n1, 'bondType1')
    c9 = AtomNode(name='C9', type='X', charge=-0.04)
    c9.set_position(x=6, y=1, z=0)
    c9.bindTo(c8, 'bondType1')
    top1_list = [c1, c2, c3, c4, cl1, c5, c6, c10, c7, n1, c8, c9]

    # construct Ligand 2
    cl11 = AtomNode(name='Cl11', type='Cl', charge=0)
    cl11.set_position(x=1, y=1, z=0)
    c11 = AtomNode(name='C11', type='C', charge=0)
    c11.set_position(x=2, y=1, z=0)
    c12 = AtomNode(name='C12', type='C', charge=-0)
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
    n11 = AtomNode(name='N11', type='N', charge=-0.02)
    n11.set_position(x=5, y=3, z=0)
    n11.bindTo(c16, 'bondType1')
    c18 = AtomNode(name='C18', type='C', charge=0)
    c18.set_position(x=6, y=1, z=0)
    c18.bindTo(c17, 'bondType1')
    c18.bindTo(n11, 'bondType1')
    c19 = AtomNode(name='C19', type='Y', charge=0.02)
    c19.set_position(x=7, y=1, z=0)
    c19.bindTo(c18, 'bondType1')
    top2_list = [cl11, c11, c12, c13, c14, c15, c16, c20, c17, n11, c18, c19]

    suptops = superimpose_topologies(top1_list, top2_list)
    suptop = suptops[0]

    assert not suptop.contains_atomNamePair('C9', 'C19')
    # we are averaging the charges of the atoms in the left and right,
    assert n1.charge == n11.charge == 0.01
    assert c9.charge == -0.01
    assert c19.charge == 0.05