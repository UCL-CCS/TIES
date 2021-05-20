"""
Focuses on the _superimpose_topology function that
superimposes the molecule in many different and
then processes the outputs to ensure the best match is found.

TODO
- this testing module should focus on "separated" molecules
"""

import pytest
import numpy as np

from ties.topology_superimposer import superimpose_topologies, Atom


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
    n1 = Atom(name='N1', atom_type='N', charge=1)
    n1.position = (4, 3, 0)
    n1.bind_to(c6, 'bondType1')
    c8 = Atom(name='C8', atom_type='C', charge=-1)
    c8.position = (5, 1, 0)
    c8.bind_to(c7, 'bondType1')
    c8.bind_to(n1, 'bondType1')
    c9 = Atom(name='C9', atom_type='C')
    c9.position = (6, 1, 0)
    c9.bind_to(c8, 'bondType1')
    top1_list = [c1, c2, c3, c4, cl1, c5, c6, c10, c7, n1, c8, c9]

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
    n11 = Atom(name='N11', atom_type='N', charge=-1)
    n11.position = (5, 3, 0)
    n11.bind_to(c16, 'bondType1')
    c18 = Atom(name='C18', atom_type='C', charge=1)
    c18.position = (6, 1, 0)
    c18.bind_to(c17, 'bondType1')
    c18.bind_to(n11, 'bondType1')
    c19 = Atom(name='C19', atom_type='C')
    c19.position = (7, 1, 0)
    c19.bind_to(c18, 'bondType1')
    top2_list = [cl11, c11, c12, c13, c14, c15, c16, c20, c17, n11, c18, c19]

    # we have to discriminate against this case somehow
    suptops = superimpose_topologies(top1_list, top2_list, align_molecules=False)
    assert not suptops[0].contains_atom_name_pair('C8', 'C18')
    assert not suptops[0].contains_atom_name_pair('N1', 'N11')
    assert not suptops[0].contains_atom_name_pair('C9', 'C19')


def test_averaging_q():
    """
    Averages the charges across the matched pairs to 0

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
    n1 = Atom(name='N1', atom_type='N', charge=.02)
    n1.position = (4, 3, 0)
    n1.bind_to(c6, 'bondType1')
    c8 = Atom(name='C8', atom_type='C', charge=-0.02)
    c8.position = (5, 1, 0)
    c8.bind_to(c7, 'bondType1')
    c8.bind_to(n1, 'bondType1')
    c9 = Atom(name='C9', atom_type='C')
    c9.position = (6, 1, 0)
    c9.bind_to(c8, 'bondType1')
    top1_list = [c1, c2, c3, c4, cl1, c5, c6, c10, c7, n1, c8, c9]

    # construct Ligand 2
    cl11 = Atom(name='Cl11', atom_type='Cl')
    cl11.position = (1, 1, 0)
    c11 = Atom(name='C11', atom_type='C')
    c11.position = (2, 1, 0)
    c12 = Atom(name='C12', atom_type='C', charge=-0)
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
    n11 = Atom(name='N11', atom_type='N', charge=-0.02)
    n11.position = (5, 3, 0)
    n11.bind_to(c16, 'bondType1')
    c18 = Atom(name='C18', atom_type='C', charge=.02)
    c18.position = (6, 1, 0)
    c18.bind_to(c17, 'bondType1')
    c18.bind_to(n11, 'bondType1')
    c19 = Atom(name='C19', atom_type='C')
    c19.position = (7, 1, 0)
    c19.bind_to(c18, 'bondType1')
    top2_list = [cl11, c11, c12, c13, c14, c15, c16, c20, c17, n11, c18, c19]

    suptops = superimpose_topologies(top1_list, top2_list, align_molecules=False)
    # we are averaging the charges of the atoms in the left and right,
    assert n1.charge == n11.charge == 0
    assert c8.charge == c18.charge == 0


def test_averaging_q_case2():
    """
    Averages charges. All atoms are matched, and each pair is within the charge tolerance limit of 0.1.

    Averaging N1-N11 to 0.01, and CL1-CL11 to -0.01

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
             CL1 (-0.04e)


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
             CL11 (+0.02e)
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
    n1 = Atom(name='N1', atom_type='N', charge=.04)
    n1.position = (4, 3, 0)
    n1.bind_to(c6, 'bondType1')
    c8 = Atom(name='C8', atom_type='C')
    c8.position = (5, 1, 0)
    c8.bind_to(c7, 'bondType1')
    c8.bind_to(n1, 'bondType1')
    cl1 = Atom(name='CL1', atom_type='CL', charge=-0.04)
    cl1.position = (6, 1, 0)
    cl1.bind_to(c8, 'bondType1')
    top1_list = [c1, c2, c3, c4, c5, c6, c10, c7, n1, c8, cl1]

    # construct Ligand 2
    c11 = Atom(name='C11', atom_type='C')
    c11.position = (2, 1, 0)
    c12 = Atom(name='C12', atom_type='C', charge=-0)
    c12.position = (2, 2, 0)
    c12.bind_to(c11, 'bondType1')
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
    n11 = Atom(name='N11', atom_type='N', charge=-0.02)
    n11.position = (5, 3, 0)
    n11.bind_to(c16, 'bondType1')
    c18 = Atom(name='C18', atom_type='C')
    c18.position = (6, 1, 0)
    c18.bind_to(c17, 'bondType1')
    c18.bind_to(n11, 'bondType1')
    cl11 = Atom(name='CL11', atom_type='CL', charge=.02)
    cl11.position = (7, 1, 0)
    cl11.bind_to(c18, 'bondType1')
    top2_list = [c11, c12, c13, c14, c15, c16, c20, c17, n11, c18, cl11]

    suptops = superimpose_topologies(top1_list, top2_list, align_molecules=False)
    suptop = suptops[0]

    # 10 atoms are should be matched
    assert len(suptop) == 11

    # we are averaging the charges of the atoms in the left and right,
    print('charges, ', n1.charge, n11.charge)
    np.testing.assert_array_almost_equal([n1.charge, n11.charge], 0.01)
    np.testing.assert_array_almost_equal([cl1.charge, cl11.charge], -0.01)


def test_q_redistribution_no_affect_mutating_atoms():
    """
    Redistribution to mutating areas: off

    Mutating atoms should not be affected by charge redistribution.

    N1-N11 are matched cancelling each other out.
    CL and BR are mismatched leading to Q imbalance.
    ie L will have to spread -0.025 across the 5 matched atoms.

    Ligand L

          C2 - C1
          /     \
         C3      N1 (+0.025e)
           \   /
             C4
             |
             CL1 (-0.025e)


    Ligand R

          C12 - C11
          /     \
        C13       N11 (-0.025e)
           \   /
             C14
             |
             BR1 (+0.025e)
    """
    # construct Ligand 1
    c1 = Atom(name='C1', atom_type='C')
    c1.position = (1, 2, 0)
    c2 = Atom(name='C2', atom_type='C')
    c2.position = (1, 1, 0)
    c1.bind_to(c2, 'bondType1')
    c3 = Atom(name='C3', atom_type='C')
    c3.position = (2, 1, 0)
    c3.bind_to(c2, 'bondType1')
    n1 = Atom(name='N1', atom_type='N', charge=.025)
    n1.position = (2, 2, 0)
    n1.bind_to(c1, 'bondType1')
    c4 = Atom(name='C4', atom_type='C')
    c4.position = (3, 1, 0)
    c4.bind_to(c3, 'bondType1')
    c4.bind_to(n1, 'bondType1')
    cl1 = Atom(name='CL1', atom_type='CL', charge=-0.025)
    cl1.position = (4, 1, 0)
    cl1.bind_to(c4, 'bondType1')
    top1_list = [c1, c2, c3, c4, cl1, n1]

    # construct Ligand 2
    c11 = Atom(name='C11', atom_type='C')
    c11.position = (1, 2, 0)
    c12 = Atom(name='C12', atom_type='C')
    c12.position = (1, 1, 0)
    c11.bind_to(c12, 'bondType1')
    c13 = Atom(name='C13', atom_type='C')
    c13.position = (2, 1, 0)
    c13.bind_to(c12, 'bondType1')
    n11 = Atom(name='N11', atom_type='N', charge=-0.025)
    n11.position = (2, 2, 0)
    n11.bind_to(c11, 'bondType1')
    c14 = Atom(name='C14', atom_type='C')
    c14.position = (3, 1, 0)
    c14.bind_to(c13, 'bondType1')
    c14.bind_to(n11, 'bondType1')
    br1 = Atom(name='BR1', atom_type='BR', charge=.025)
    br1.position = (4, 1, 0)
    br1.bind_to(c14, 'bondType1')
    top2_list = [c11, c12, c13, c14, br1, n11]

    suptops = superimpose_topologies(top1_list, top2_list, align_molecules=False,
                                     redistribute_charges_over_unmatched=False)
    suptop = suptops[0]

    # 5 atoms are should be matched
    assert len(suptop) == 5

    # we are averaging the charges of the atoms in the left and right,
    print('charges, ', n1.charge, n11.charge)
    # np.testing.assert_array_almost_equal([n1.charge, n11.charge], 0.01)
    # np.testing.assert_array_almost_equal([cl1.charge, br1.charge], -0.01)

    # the mutating atoms retain their original q
    assert cl1.charge == -0.025
    assert br1.charge == 0.025



def test_redistribution_q2():
    """
    Averages charges. All atoms are matched, and each pair is within the charge tolerance limit of 0.1.

    The mismatched region CL1-CL11 has charges -0.04 in L and 0.02 in R,
    which balances the matched N1-N11 with 0.04 in L, and -0.02 in R.

    Averaging N1-N11 to 0.01, and CL1-CL11 to -0.01

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
             CL1 (-0.04e)


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
             CL11 (+0.02e)
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
    n1 = Atom(name='N1', atom_type='N', charge=.04)
    n1.position = (4, 3, 0)
    n1.bind_to(c6, 'bondType1')
    c8 = Atom(name='C8', atom_type='C')
    c8.position = (5, 1, 0)
    c8.bind_to(c7, 'bondType1')
    c8.bind_to(n1, 'bondType1')
    cl1 = Atom(name='CL1', atom_type='CL', charge=-0.04)
    cl1.position = (6, 1, 0)
    cl1.bind_to(c8, 'bondType1')
    top1_list = [c1, c2, c3, c4, c5, c6, c10, c7, n1, c8, cl1]

    # construct Ligand 2
    c11 = Atom(name='C11', atom_type='C')
    c11.position = (2, 1, 0)
    c12 = Atom(name='C12', atom_type='C', charge=-0)
    c12.position = (2, 2, 0)
    c12.bind_to(c11, 'bondType1')
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
    n11 = Atom(name='N11', atom_type='N', charge=-0.02)
    n11.position = (5, 3, 0)
    n11.bind_to(c16, 'bondType1')
    c18 = Atom(name='C18', atom_type='C')
    c18.position = (6, 1, 0)
    c18.bind_to(c17, 'bondType1')
    c18.bind_to(n11, 'bondType1')
    cl11 = Atom(name='CL11', atom_type='CL', charge=.02)
    cl11.position = (7, 1, 0)
    cl11.bind_to(c18, 'bondType1')
    top2_list = [c11, c12, c13, c14, c15, c16, c20, c17, n11, c18, cl11]

    suptops = superimpose_topologies(top1_list, top2_list, align_molecules=False)
    suptop = suptops[0]

    # 10 atoms are should be matched
    assert len(suptop) == 11

    # we are averaging the charges of the atoms in the left and right,
    print('charges, ', n1.charge, n11.charge)
    np.testing.assert_array_almost_equal([n1.charge, n11.charge], 0.01)
    np.testing.assert_array_almost_equal([cl1.charge, cl11.charge], -0.01)


def test_averaging_charges_imbalance_distribution_multiple():
    """
    Averages charges in the matched region.
    The mismatched region C9-C19 had charges -0.04 in L and 0.02 in R,
    which balances the matched N1-N11 with 0.04 in L, and -0.02 in R.
    After averaging N1-N11 to 0.01, L has -0.03 imbalance, and R has 0.03 imbalance.
    This imbalance is corrected in the unmatched area.

    In L, -0.03 is spread over X8 and X9 and Cl1, and
    in R, 0.03 is spread over X18 and X19 and Cl11


    Ligand 1

         C1 - C2
         /      \
    Cl1-C3      C4
          \     /
          C5 - C6
          /     \
     C10-C7      N1 (+0.04e)
           \   /
             CL8 (0e)
             |
             CL9 (-0.04e)


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
             BR18 (0e)
             |
             BR19 (+0.02e)
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
    n1 = Atom(name='N1', atom_type='N', charge=.04)
    n1.position = (4, 3, 0)
    n1.bind_to(c6, 'bondType1')
    c8 = Atom(name='CL8', atom_type='CL')
    c8.position = (5, 1, 0)
    c8.bind_to(c7, 'bondType1')
    c8.bind_to(n1, 'bondType1')
    c9 = Atom(name='CL9', atom_type='CL', charge=-0.04)
    c9.position = (6, 1, 0)
    c9.bind_to(c8, 'bondType1')
    top1_list = [c1, c2, c3, c4, cl1, c5, c6, c10, c7, n1, c8, c9]

    # construct Ligand 2
    cl11 = Atom(name='Cl11', atom_type='Cl')
    cl11.position = (1, 1, 0)
    c11 = Atom(name='C11', atom_type='C')
    c11.position = (2, 1, 0)
    c12 = Atom(name='C12', atom_type='C', charge=-0)
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
    n11 = Atom(name='N11', atom_type='N', charge=-0.02)
    n11.position = (5, 3, 0)
    n11.bind_to(c16, 'bondType1')
    c18 = Atom(name='BR18', atom_type='BR')
    c18.position = (6, 1, 0)
    c18.bind_to(c17, 'bondType1')
    c18.bind_to(n11, 'bondType1')
    c19 = Atom(name='BR19', atom_type='BR', charge=.02)
    c19.position = (7, 1, 0)
    c19.bind_to(c18, 'bondType1')
    top2_list = [cl11, c11, c12, c13, c14, c15, c16, c20, c17, n11, c18, c19]

    suptops = superimpose_topologies(top1_list, top2_list, align_molecules=False, partial_rings_allowed=True)
    suptop = suptops[0]

    assert not suptop.contains_atom_name_pair('C9', 'C19')
    # charge averaging
    np.testing.assert_array_almost_equal([n1.charge, n11.charge], 0.01)
    # charge imbalance distributed in L
    np.testing.assert_almost_equal(c9.charge, -0.03)
    np.testing.assert_almost_equal(c8.charge, 0.01)
    np.testing.assert_almost_equal(cl1.charge, 0.01)

    # charge imbalance distributed in R
    np.testing.assert_almost_equal(c19.charge, 0.01)
    np.testing.assert_almost_equal(c18.charge, -0.01)
    np.testing.assert_almost_equal(cl11.charge, -0.01)


def test_averaging_charges_imbalance_distribution_2to2():
    """
    Charge averaging across more atoms, and charge distribution across more atoms.

    Ligand 1
    Net Change in matched: 0.01 + -0.03 = -0.02
    This is matched by 0.02 distributed over Cl1, X8, X9 in L, so 0.00(6) for each
         C1 - C2 (-0.01e) after avg (0e)
         /      \
    Cl1-C3      C4
          \     /
          C5 - C6
          /     \
     C10-C7      N1 (+0.04e) after avg (0.01e)
           \   /
             X8 (0.01e)
             |
             X9 (-0.04e)


    Ligand 2
    Net Change in matched: -0.01 + 0.03 = 0.02,
    this is matched by -0.02 distributed -0.01 on Y18, and on Y19
         C11 - C12 (+0.01e) after averaging (0.00)
         /      \
        C13      C14
          \     /
          C15 - C16
          /     \
     C20-C17       N11 (-0.02e) after averaging (0.01)
           \   /
             Y18 (-0.01e) should be (-0.02e)
             |
             Y19 (+0.02e) should be (0.01e)
    """
    # construct LIGAND 1
    c1 = Atom(name='C1', atom_type='C')
    c1.position = (1, 1, 0)
    c2 = Atom(name='C2', atom_type='C', charge=-0.01)
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
    n1 = Atom(name='N1', atom_type='N', charge=.04)
    n1.position = (4, 3, 0)
    n1.bind_to(c6, 'bondType1')
    c8 = Atom(name='C8', atom_type='X', charge=.01)
    c8.position = (5, 1, 0)
    c8.bind_to(c7, 'bondType1')
    c8.bind_to(n1, 'bondType1')
    c9 = Atom(name='C9', atom_type='X', charge=-0.04)
    c9.position = (6, 1, 0)
    c9.bind_to(c8, 'bondType1')
    top1_list = [c1, c2, c3, c4, cl1, c5, c6, c10, c7, n1, c8, c9]

    # construct Ligand 2
    c11 = Atom(name='C11', atom_type='C')
    c11.position = (2, 1, 0)
    c12 = Atom(name='C12', atom_type='C', charge=.01)
    c12.position = (2, 2, 0)
    c12.bind_to(c11, 'bondType1')
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
    n11 = Atom(name='N11', atom_type='N', charge=-0.02)
    n11.position = (5, 3, 0)
    n11.bind_to(c16, 'bondType1')
    c18 = Atom(name='C18', atom_type='Y', charge=-0.01)
    c18.position = (6, 1, 0)
    c18.bind_to(c17, 'bondType1')
    c18.bind_to(n11, 'bondType1')
    c19 = Atom(name='C19', atom_type='Y', charge=.02)
    c19.position = (7, 1, 0)
    c19.bind_to(c18, 'bondType1')
    top2_list = [c11, c12, c13, c14, c15, c16, c20, c17, n11, c18, c19]

    suptops = superimpose_topologies(top1_list, top2_list)
    suptop = suptops[0]

    assert not suptop.contains_atom_name_pair('C9', 'C19')
    assert not suptop.contains_atom_name_pair('C8', 'C18')
    # charge averaging
    np.testing.assert_array_almost_equal([n1.charge, n11.charge], 0.01)
    np.testing.assert_array_almost_equal([c2.charge, c12.charge], 0)
    # charge imbalance distributed in L
    np.testing.assert_almost_equal(c9.charge, -0.04 - -0.02/3)
    np.testing.assert_almost_equal(c8.charge, 0.01 - -0.02/3)
    np.testing.assert_almost_equal(cl1.charge, 0 --0.02/3)

    # charge imbalance distributed in R
    np.testing.assert_almost_equal(c18.charge, -0.02)
    np.testing.assert_almost_equal(c19.charge, 0.01)


def test_averaging_charges_imbalance_distribution_2to2():
    """

         C1 - C2
         /      \
    Cl1-C3      C4
          \     /
          C5 - C6
          /     \
     C10-C7      N1
           \   /
             CL8
             |
             CL9

         C11 - C12
         /      \
        C13      C14
          \     /
          C15 - C16
          /     \
     C20-C17       N11
           \   /
             BR18
             |
             BR19
    """
    # construct LIGAND 1
    c1 = Atom(name='C1', atom_type='C')
    c1.position = (1, 1, 0)
    c2 = Atom(name='C2', atom_type='C', charge=-0.01)
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
    n1 = Atom(name='N1', atom_type='N', charge=.04)
    n1.position = (4, 3, 0)
    n1.bind_to(c6, 'bondType1')
    c8 = Atom(name='CL8', atom_type='CL', charge=.01)
    c8.position = (5, 1, 0)
    c8.bind_to(c7, 'bondType1')
    c8.bind_to(n1, 'bondType1')
    c9 = Atom(name='CL9', atom_type='CL', charge=-0.04)
    c9.position = (6, 1, 0)
    c9.bind_to(c8, 'bondType1')
    top1_list = [c1, c2, c3, c4, cl1, c5, c6, c10, c7, n1, c8, c9]

    # construct Ligand 2
    c11 = Atom(name='C11', atom_type='C')
    c11.position = (2, 1, 0)
    c12 = Atom(name='C12', atom_type='C', charge=.01)
    c12.position = (2, 2, 0)
    c12.bind_to(c11, 'bondType1')
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
    n11 = Atom(name='N11', atom_type='N', charge=-0.02)
    n11.position = (5, 3, 0)
    n11.bind_to(c16, 'bondType1')
    c18 = Atom(name='BR18', atom_type='BR', charge=-0.01)
    c18.position = (6, 1, 0)
    c18.bind_to(c17, 'bondType1')
    c18.bind_to(n11, 'bondType1')
    c19 = Atom(name='BR19', atom_type='BR', charge=.02)
    c19.position = (7, 1, 0)
    c19.bind_to(c18, 'bondType1')
    top2_list = [c11, c12, c13, c14, c15, c16, c20, c17, n11, c18, c19]

    suptops = superimpose_topologies(top1_list, top2_list, align_molecules=False, partial_rings_allowed=True)
    suptop = suptops[0]

    assert not suptop.contains_atom_name_pair('C9', 'C19')
    assert not suptop.contains_atom_name_pair('C8', 'C18')
    # charge averaging
    np.testing.assert_array_almost_equal([n1.charge, n11.charge], 0.01)
    np.testing.assert_array_almost_equal([c2.charge, c12.charge], 0)
    # charge imbalance distributed in L
    np.testing.assert_almost_equal(c9.charge, -0.04 - -0.02/3)
    np.testing.assert_almost_equal(c8.charge, 0.01 - -0.02/3)
    np.testing.assert_almost_equal(cl1.charge, 0 --0.02/3)

    # charge imbalance distributed in R
    np.testing.assert_almost_equal(c18.charge, -0.02)
    np.testing.assert_almost_equal(c19.charge, 0.01)


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
    c1 = Atom(name='C1', atom_type='C')
    c1.position = (1, 1, 0)
    c2 = Atom(name='C2', atom_type='C', charge=.07)
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
    n1 = Atom(name='N1', atom_type='N', charge=-0.04)
    n1.position = (4, 3, 0)
    n1.bind_to(c6, 'bondType1')
    c8 = Atom(name='C8', atom_type='C', charge=-0.03)
    c8.position = (5, 1, 0)
    c8.bind_to(c7, 'bondType1')
    c8.bind_to(n1, 'bondType1')
    c9 = Atom(name='C9', atom_type='C')
    c9.position = (6, 1, 0)
    c9.bind_to(c8, 'bondType1')
    top1_list = [c1, c2, c3, c4, cl1, c5, c6, c10, c7, n1, c8, c9]

    # construct Ligand 2
    cl11 = Atom(name='Cl11', atom_type='Cl')
    cl11.position = (1, 1, 0)
    c11 = Atom(name='C11', atom_type='C')
    c11.position = (2, 1, 0)
    c12 = Atom(name='C12', atom_type='C', charge=-0.07)
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
    n11 = Atom(name='N11', atom_type='N', charge=.04)
    n11.position = (5, 3, 0)
    n11.bind_to(c16, 'bondType1')
    c18 = Atom(name='C18', atom_type='C', charge=.03)
    c18.position = (6, 1, 0)
    c18.bind_to(c17, 'bondType1')
    c18.bind_to(n11, 'bondType1')
    c19 = Atom(name='C19', atom_type='C')
    c19.position = (7, 1, 0)
    c19.bind_to(c18, 'bondType1')
    top2_list = [cl11, c11, c12, c13, c14, c15, c16, c20, c17, n11, c18, c19]

    # we have to discriminate against this case somehow
    suptops = superimpose_topologies(top1_list, top2_list, align_molecules=False, partial_rings_allowed=True)
    # removed because the charge difference was too large
    assert not suptops[0].contains_atom_name_pair('C2', 'C12')
    # removed to bring the net charge difference below 0.1e
    assert not suptops[0].contains_atom_name_pair('N1', 'N11')
    # this pair should remain as the net charge difference is below 0.1
    assert suptops[0].contains_atom_name_pair('C8', 'C18')
    assert len(suptops[0]) == 9



def test_perfect_ring():
    """
    A simple case with a clear ring match.
    It should be matched and found.

    Ligand 1

         C1 - C2 (-0.02e)
         /      \
        C3      C4 (+0.02e)
          \     /
          C5 - C6
                /
                C7
                \
               N1
              /
             C8
             |
             C9
             \
             C10


    Ligand 2

         C11 - C12 (+0.02e)
         /      \
        C13      C14 (-0.02e)
          \     /
          C15 - C16
                /
                C17
                \
                N11
               /
             C18
             |
             C19
             \
             C20
    """
    # construct LIGAND 1
    c1 = Atom(name='C1', atom_type='C')
    c1.position = (1, 1, 0)
    c2 = Atom(name='C2', atom_type='C', charge=-0.02)
    c2.position = (1, 2, 0)
    c1.bind_to(c2, 'bondType1')
    c3 = Atom(name='C3', atom_type='C')
    c3.position = (2, 2, 0)
    c3.bind_to(c1, 'bondType1')
    c4 = Atom(name='C4', atom_type='C', charge=.02)
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
    c7.bind_to(c6, 'bondType1')
    n1 = Atom(name='N1', atom_type='N')
    n1.position = (4, 3, 0)
    n1.bind_to(c7, 'bondType1')
    c8 = Atom(name='C8', atom_type='C')
    c8.position = (5, 1, 0)
    c8.bind_to(n1, 'bondType1')
    c9 = Atom(name='C9', atom_type='C')
    c9.position = (6, 1, 0)
    c9.bind_to(c8, 'bondType1')
    c10 = Atom(name='C10', atom_type='C')
    c10.position = (4, 1, 0)
    c10.bind_to(c9, 'bondType1')
    top1_list = [c1, c2, c3, c4, c5, c6, c10, c7, n1, c8, c9]

    # construct Ligand 2
    c11 = Atom(name='C11', atom_type='C')
    c11.position = (2, 1, 0)
    c12 = Atom(name='C12', atom_type='C', charge=.02)
    c12.position = (2, 2, 0)
    c12.bind_to(c11, 'bondType1')
    c13 = Atom(name='C13', atom_type='C')
    c13.position = (3, 1, 0)
    c13.bind_to(c11, 'bondType1')
    c14 = Atom(name='C14', atom_type='C', charge=-0.02)
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
    c17.bind_to(c16, 'bondType1')
    n11 = Atom(name='N11', atom_type='N')
    n11.position = (5, 3, 0)
    n11.bind_to(c17, 'bondType1')
    c18 = Atom(name='C18', atom_type='C')
    c18.position = (6, 1, 0)
    c18.bind_to(n11, 'bondType1')
    c19 = Atom(name='C19', atom_type='C')
    c19.position = (7, 1, 0)
    c19.bind_to(c18, 'bondType1')
    c20 = Atom(name='C20', atom_type='C')
    c20.position = (5, 1, 0)
    c20.bind_to(c19, 'bondType1')
    top2_list = [c11, c12, c13, c14, c15, c16, c20, c17, n11, c18, c19]

    suptops = superimpose_topologies(top1_list, top2_list, align_molecules=False)
    suptop = suptops[0]

    ring_overlap = [('C1', 'C11'), ('C2', 'C12'), ('C3', 'C13'),
                    ('C4', 'C14'), ('C5', 'C15'), ('C6', 'C16')]
    for a1, a2 in ring_overlap[::-1]:
        if suptop.contains_atom_name_pair(a1, a2):
            ring_overlap.remove((a1, a2))
    assert len(ring_overlap) == 0



def test_partial_ring():
    """
    The ring in each molecule is not the same.
    Therefore, the entire ring should be "unmatched".
    fixme - charges should not define on the whether the ring is partial or not, but the atom atom_type only

    Ligand 1

         C1 - C2 (-0.2e)
         /      \
        C3      C4 (+0.2e)
          \     /
          C5 - C6
                /
                C7
                \
               N1
              /
             C8
             |
             C9
             \
             C10


    Ligand 2

         C11 - C12 (+0.2e)
         /      \
        C13      C14 (-0.2e)
          \     /
          C15 - C16
                /
                C17
                \
                N11
               /
             C18
             |
             C19
             \
             C20
    """
    # construct LIGAND 1
    c1 = Atom(name='C1', atom_type='C')
    c1.position = (1, 1, 0)
    c2 = Atom(name='C2', atom_type='C', charge=-0.2)
    c2.position = (1, 2, 0)
    c1.bind_to(c2, 'bondType1')
    c3 = Atom(name='C3', atom_type='C')
    c3.position = (2, 2, 0)
    c3.bind_to(c1, 'bondType1')
    c4 = Atom(name='C4', atom_type='C', charge=.2)
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
    c7.bind_to(c6, 'bondType1')
    n1 = Atom(name='N1', atom_type='N')
    n1.position = (4, 3, 0)
    n1.bind_to(c7, 'bondType1')
    c8 = Atom(name='C8', atom_type='C')
    c8.position = (5, 1, 0)
    c8.bind_to(n1, 'bondType1')
    c9 = Atom(name='C9', atom_type='C')
    c9.position = (6, 1, 0)
    c9.bind_to(c8, 'bondType1')
    c10 = Atom(name='C10', atom_type='C')
    c10.position = (4, 1, 0)
    c10.bind_to(c9, 'bondType1')
    top1_list = [c1, c2, c3, c4, c5, c6, c10, c7, n1, c8, c9]

    # construct Ligand 2
    c11 = Atom(name='C11', atom_type='C')
    c11.position = (2, 1, 0)
    c12 = Atom(name='C12', atom_type='C', charge=.2)
    c12.position = (2, 2, 0)
    c12.bind_to(c11, 'bondType1')
    c13 = Atom(name='C13', atom_type='C')
    c13.position = (3, 1, 0)
    c13.bind_to(c11, 'bondType1')
    c14 = Atom(name='C14', atom_type='C', charge=-0.2)
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
    c17.bind_to(c16, 'bondType1')
    n11 = Atom(name='N11', atom_type='N')
    n11.position = (5, 3, 0)
    n11.bind_to(c17, 'bondType1')
    c18 = Atom(name='C18', atom_type='C')
    c18.position = (6, 1, 0)
    c18.bind_to(n11, 'bondType1')
    c19 = Atom(name='C19', atom_type='C')
    c19.position = (7, 1, 0)
    c19.bind_to(c18, 'bondType1')
    c20 = Atom(name='C20', atom_type='C')
    c20.position = (5, 1, 0)
    c20.bind_to(c19, 'bondType1')
    top2_list = [c11, c12, c13, c14, c15, c16, c20, c17, n11, c18, c19]

    suptops = superimpose_topologies(top1_list, top2_list, align_molecules=False, partial_rings_allowed=False)
    suptop = suptops[0]
    # check if the mismatched charges removed the right atoms
    assert not suptop.contains_atom_name_pair('C2', 'C12')
    assert not suptop.contains_atom_name_pair('C4', 'C14')
    # partial ring should have been removed
    partial_ring = [('C1', 'C11'), ('C3', 'C13'), ('C5', 'C15'), ('C6', 'C16')]
    for a1, a2 in partial_ring[::-1]:
        if not suptop.contains_atom_name_pair(a1, a2):
            partial_ring.remove((a1, a2))
    assert len(partial_ring) == 0


def test_partial_ring_cascade():
    """
    http://www.alchemistry.org/wiki/Constructing_a_Pathway_of_Intermediate_States#Alchemically_Transforming_Rings
    When atoms are unmatched due to partial rings,
    they might disturb another ring.
    In such a situation, the second ring will not be removed (ie no cascading)

    Ligand 1

         C1 - C2 (0.1e)
         /      \
        C3      C4 (-0.1e)
          \     /
          C5 - C6
          /     \
         C7      N1
           \   /
             C8
             |
             C9
             |
             C10


    Ligand 2

         C11 - C12 (-0.1e)
         /      \
        C13      C14 (0.1e)
          \     /
          C15 - C16
          /     \
         C17       N11
           \   /
             C18
             |
             C19
             \
             C20
    """
    # construct LIGAND 1
    c1 = Atom(name='C1', atom_type='C')
    c1.position = (1, 1, 0)
    c2 = Atom(name='C2', atom_type='C', charge=.1)
    c2.position = (1, 2, 0)
    c1.bind_to(c2, 'bondType1')
    c3 = Atom(name='C3', atom_type='C')
    c3.position = (2, 2, 0)
    c3.bind_to(c1, 'bondType1')
    c4 = Atom(name='C4', atom_type='C', charge=-0.1)
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
    c10 = Atom(name='C10', atom_type='C')
    c10.position = (4, 1, 0)
    c10.bind_to(c9, 'bondType1')
    top1_list = [c1, c2, c3, c4, c5, c6, c10, c7, n1, c8, c9]

    # construct Ligand 2
    c11 = Atom(name='C11', atom_type='C')
    c11.position = (2, 1, 0)
    c12 = Atom(name='C12', atom_type='C', charge=-0.1)
    c12.position = (2, 2, 0)
    c12.bind_to(c11, 'bondType1')
    c13 = Atom(name='C13', atom_type='C')
    c13.position = (3, 1, 0)
    c13.bind_to(c11, 'bondType1')
    c14 = Atom(name='C14', atom_type='C', charge=.1)
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
    c20 = Atom(name='C20', atom_type='C')
    c20.position = (5, 1, 0)
    c20.bind_to(c19, 'bondType1')
    top2_list = [c11, c12, c13, c14, c15, c16, c20, c17, n11, c18, c19]

    suptops = superimpose_topologies(top1_list, top2_list, align_molecules=False,
                                     partial_rings_allowed=False)
    suptop = suptops[0]
    # atoms were removed due to the charge
    assert not suptop.contains_atom_name_pair('C2', 'C12')
    assert not suptop.contains_atom_name_pair('C4', 'C14')
    # the length of the entire match should be 2
    assert len(suptop) == 2
    # which are
    assert suptop.contains_atom_name_pair('C9', 'C19')
    assert suptop.contains_atom_name_pair('C10', 'C20')



def test_partial_rings_overlap():
    """
    http://www.alchemistry.org/wiki/Constructing_a_Pathway_of_Intermediate_States#Alchemically_Transforming_Rings
    The shared atom across the two rings
    has a different charge.
    Therefore, no partial rings are allowed.

    Ligand 1

         C1 - C2 (0.1e)
         /      \
        C3      C4
          \     /
          C5 - C6 (-0.1e)
          /     \
         C7      N1
           \   /
             C8
             |
             C9


    Ligand 2

         C11 - C12 (-0.1e)
         /      \
        C13      C14
          \     /
          C15 - C16 (0.1e)
          /     \
        C17       N11
           \   /
             C18
             |
             C19
    """
    # construct LIGAND 1
    c1 = Atom(name='C1', atom_type='C')
    c1.position = (1, 1, 0)
    c2 = Atom(name='C2', atom_type='C', charge=.1)
    c2.position = (1, 2, 0)
    c1.bind_to(c2, 'bondType1')
    c3 = Atom(name='C3', atom_type='C')
    c3.position = (2, 2, 0)
    c3.bind_to(c1, 'bondType1')
    c4 = Atom(name='C4', atom_type='C')
    c4.position = (2, 3, 0)
    c4.bind_to(c2, 'bondType1')
    c5 = Atom(name='C5', atom_type='C')
    c5.position = (3, 1, 0)
    c5.bind_to(c3, 'bondType1')
    c6 = Atom(name='C6', atom_type='C', charge=-0.1)
    c6.position = (3, 2, 0)
    c6.bind_to(c5, 'bondType1')
    c6.bind_to(c4, 'bondType1')
    c7 = Atom(name='C7', atom_type='C')
    c7.position = (4, 2, 0)
    c7.bind_to(c5, 'bondType1')
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
    top1_list = [c1, c2, c3, c4, c5, c6, c7, n1, c8, c9]

    # construct Ligand 2
    c11 = Atom(name='C11', atom_type='C')
    c11.position = (2, 1, 0)
    c12 = Atom(name='C12', atom_type='C', charge=-0.1)
    c12.position = (2, 2, 0)
    c12.bind_to(c11, 'bondType1')
    c13 = Atom(name='C13', atom_type='C')
    c13.position = (3, 1, 0)
    c13.bind_to(c11, 'bondType1')
    c14 = Atom(name='C14', atom_type='C')
    c14.position = (3, 2, 0)
    c14.bind_to(c12, 'bondType1')
    c15 = Atom(name='C15', atom_type='C')
    c15.position = (4, 1, 0)
    c15.bind_to(c13, 'bondType1')
    c16 = Atom(name='C16', atom_type='C', charge=.1)
    c16.position = (4, 2, 0)
    c16.bind_to(c15, 'bondType1')
    c16.bind_to(c14, 'bondType1')
    c17 = Atom(name='C17', atom_type='C')
    c17.position = (5, 2, 0)
    c17.bind_to(c15, 'bondType1')
    n11 = Atom(name='N11', atom_type='N', charge=-0)
    n11.position = (5, 3, 0)
    n11.bind_to(c16, 'bondType1')
    c18 = Atom(name='C18', atom_type='C')
    c18.position = (6, 1, 0)
    c18.bind_to(c17, 'bondType1')
    c18.bind_to(n11, 'bondType1')
    c19 = Atom(name='C19', atom_type='C')
    c19.position = (7, 1, 0)
    c19.bind_to(c18, 'bondType1')
    top2_list = [c11, c12, c13, c14, c15, c16, c17, n11, c18, c19]

    suptops = superimpose_topologies(top1_list, top2_list, align_molecules=False,
                                     partial_rings_allowed=False)
    suptop = suptops[0]
    # only one atom should be left,
    assert len(suptop) == 1
    assert suptop.contains_atom_name_pair('C9', 'C19')


def test_averaging_q():
    return
    """
    Two atoms cross the pair-tolerance so each should be removed.
    However, all other atoms are fine 

    Ligand 1
         C1
          \
         C2
          /
         C3
          \
          C4
          /
          C5
          \
          C6
          /
          C7
          \
          C8
          /
          C9 (-0.2e)
          \
          C10 (0.2e)


    Ligand 2
         C11
          \
         C12
          /
         C13
          \
          C14
          /
          C15
          \
          C16
          /
          C17
          \
          C18
          /
          C19 (0.2e)
          \
          C20 (-0.2e)
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
    n1 = Atom(name='N1', atom_type='N', charge=.02)
    n1.position = (4, 3, 0)
    n1.bind_to(c6, 'bondType1')
    c8 = Atom(name='C8', atom_type='C', charge=-0.02)
    c8.position = (5, 1, 0)
    c8.bind_to(c7, 'bondType1')
    c8.bind_to(n1, 'bondType1')
    c9 = Atom(name='C9', atom_type='C')
    c9.position = (6, 1, 0)
    c9.bind_to(c8, 'bondType1')
    top1_list = [c1, c2, c3, c4, cl1, c5, c6, c10, c7, n1, c8, c9]

    # construct Ligand 2
    cl11 = Atom(name='Cl11', atom_type='Cl')
    cl11.position = (1, 1, 0)
    c11 = Atom(name='C11', atom_type='C')
    c11.position = (2, 1, 0)
    c12 = Atom(name='C12', atom_type='C', charge=-0)
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
    n11 = Atom(name='N11', atom_type='N', charge=-0.02)
    n11.position = (5, 3, 0)
    n11.bind_to(c16, 'bondType1')
    c18 = Atom(name='C18', atom_type='C', charge=.02)
    c18.position = (6, 1, 0)
    c18.bind_to(c17, 'bondType1')
    c18.bind_to(n11, 'bondType1')
    c19 = Atom(name='C19', atom_type='C')
    c19.position = (7, 1, 0)
    c19.bind_to(c18, 'bondType1')
    top2_list = [cl11, c11, c12, c13, c14, c15, c16, c20, c17, n11, c18, c19]

    suptops = superimpose_topologies(top1_list, top2_list, align_molecules=False)
    # we are averaging the charges of the atoms in the left and right,
    assert n1.charge == n11.charge == 0
    assert c8.charge == c18.charge == 0
