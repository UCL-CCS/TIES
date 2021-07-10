"""
Focuses on the _superimpose_topology function that
superimposes the molecule in many different and
then processes the outputs to ensure the best match is found.

TODO
- this testing module should focus on "separated" molecules
"""

import pytest
import json
from ties import Pair
from ties import Config
from ties.topology_superimposer import _superimpose_topologies, Atom, get_starting_configurations


def test_2diff_atoms_cn():
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
    top1_list = [c1, n1]

    # construct the LIGAND 2
    c11 = Atom(name='C11', atom_type='C')
    c11.position = (1, 1, 0)
    n11 = Atom(name='N11', atom_type='N')
    n11.position = (1, 2, 0)
    c11.bind_to(n11, 'bondType1')
    top2_list = [c11, n11]

    # should return a list with an empty sup_top
    suptops = _superimpose_topologies(top1_list, top2_list, starting_pairs_heuristics=False)
    # it should return an empty suptop
    assert len(suptops) == 1
    assert len(suptops[0]) == 2


def test_3diff_atoms_cno_right_start():
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
    top1_list = [c1, n1, o1]

    # construct the LIGAND 2
    c11 = Atom(name='C11', atom_type='C')
    c11.position = (1, 1, 0)
    n11 = Atom(name='N11', atom_type='N')
    n11.position = (1, 2, 0)
    c11.bind_to(n11, 'bondType1')
    o11 = Atom(name='O11', atom_type='O')
    o11.position = (1, 3, 0)
    o11.bind_to(n11, 'bondType1')
    top2_list = [c11, n11, o11]

    # should overlap 2 atoms
    suptops = _superimpose_topologies(top1_list, top2_list, starting_pairs_heuristics=False)
    assert len(suptops) == 1
    suptop = suptops[0]

    # the number of overlapped atoms is two
    assert len(suptop) == 3
    correct_overlaps = [('C1', 'C11'), ('N1', 'N11'), ('O1', 'O11')]
    for atomName1, atomName2 in correct_overlaps:
        assert suptop.contains_atom_name_pair(atomName1, atomName2)

    # no mirrors
    assert len(suptop.mirrors) == 0


def test_simple_multiple_solutions(lr_atoms_branches):
    """
    A simple molecule chain with an ester.
    The ester allows for mapping (O1-O11, O2-O12) and (O1-O12, O2-O11)
    # fixme what if you start from the wrong O-O matching? should that be a case? how does
    comparies topologies would behave in that case?
    """
    top1_list, top2_list = lr_atoms_branches

    # should be two topologies
    suptops = _superimpose_topologies(top1_list, top2_list, starting_pairs_heuristics=False)
    # there is one solution
    assert len(suptops) == 1

    correct_overlaps = [('C1', 'C11'), ('N1', 'N11'), ('O1', 'O11'), ('O2', 'O12')]
    for st in suptops:
        for atomName1, atomName2 in correct_overlaps:
            assert st.contains_atom_name_pair(atomName1, atomName2)

    # fixme - add a test case for the superimposer function that makes use of _overlay,
    # this is to resolve multiple solutions such as the one here


def test_SimpleMultipleSolutions_mirrors(lr_atoms_branches):
    """
    A simple molecule chain with an ester.
    The ester allows for mapping (O1-O11, O2-O12) and (O1-O12, O2-O11)
    # fixme what if you start from the wrong O-O matching? should that be a case? how does
    comparies topologies would behave in that case?
    """
    top1_list, top2_list = lr_atoms_branches

    # should be two topologies
    suptops = _superimpose_topologies(top1_list, top2_list, starting_pairs_heuristics=False)
    # note that mirros have not been finished
    return
    # there is one solution
    assert len(suptops) == 1
    assert len(suptops[0].mirrors) == 1

    correct_overlaps = [('C1', 'C11'), ('N1', 'N11'), ('O1', 'O11'), ('O2', 'O12')]
    for st in suptops:
        for atomName1, atomName2 in correct_overlaps:
            assert st.contains_atom_name_pair(atomName1, atomName2)

    # fixme - add a test case for the superimposer function that makes use of _overlay,
    # this is to resolve multiple solutions such as the one here


def test_2same_atoms_2c_symmetry(lr_atoms_chain_2c):
    """
    Two solutions with different starting points.
    """
    top1_list, top2_list = lr_atoms_chain_2c
    # should return a list with an empty sup_top
    suptops = _superimpose_topologies(top1_list, top2_list, starting_pairs_heuristics=False)
    assert len(suptops) == 1
    assert len(suptops[0]) == 2


# def test_2same_atoms_2c_symmetry_mirrors(lr_atoms_chain_2c):
#     """
#     Two solutions with different starting points.
#     """
#     top1_list, top2_list = lr_atoms_chain_2c
#
#     # should return a list with an empty sup_top
#     suptops = _superimpose_topologies(top1_list, top2_list)
#     # fixme - note that this we have not explicitly solved mirrors
    # assert len(suptops) == 1
    # assert len(suptops[0]) == 2
    # assert len(suptops[0].mirrors) == 1
    # assert len(suptops[0].mirrors[0]) == 2


def test_3c_circle():
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
    top1_list = [c1, c2, c3]

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
    top2_list = [c11, c12, c13]

    suptops = _superimpose_topologies(top1_list, top2_list, starting_pairs_heuristics=False)
    assert len(suptops) == 1
    assert len(suptops[0]) == 3
    assert len(suptops[0].mirrors) > 2


def test_mcl1_l12l35_crossed_double_cycle(dual_ring1, dual_ring2):
    """
    Molecule inspired by Agastya's dataset (mcl1_l12l35).
    This test checks if cycles are properly tracked.
    There exists a pathway here that would overestimate the superimposition.
    The larger incorrect solution finds a way to match the Cl atoms.
    """

    # we have to discriminate against this case somehow
    suptops = _superimpose_topologies(dual_ring1, dual_ring2, starting_pairs_heuristics=False)
    assert len(suptops) == 1
    assert len(suptops[0]) == 11


def test_api_serializable_hybrid():
    # .toJSON metadata can be converted into JSON

    config = Config()
    config.ligand_net_charge = -1

    pair = Pair('../examples/mol2_2ligands_MCL1/l02.mol2', '../examples/mol2_2ligands_MCL1/l03.mol2', config=config)
    pair.make_atom_names_unique()

    # overwrite the previous config settings with relevant parameters
    hybrid = pair.superimpose(use_element_in_superimposition=True, redistribute_q_over_unmatched=True)

    to_json = json.dumps(hybrid.toJSON(), indent=4)
    json.loads(to_json)


def test_refine_against_charges_order_problem():
    return
    """
    When accounting for charge refinement (when a pair is impatible due to the charges), 
    we are removing the incompatible nodes together with their dangling hydrogens. 
    Note that in theory the dangling hydrogens shoould be removed by "disconnected components"

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
    top2_list = [cl11, c11, c12, c13, c14, c15, c16, c20, c17, n11, c18, c19]

    # we have to discriminate against this case somehow
    suptops = _superimpose_topologies(top1_list, top2_list, starting_pairs_heuristics=False)
    assert len(suptops) == 1
    assert len(suptops[0]) == 11


def test_get_starting_configuration(dual_ring1, dual_ring2):
    starting_configurations = get_starting_configurations(dual_ring1, dual_ring2, fraction=0.1, filter_ring_c=True)

    # we expect over 10% (0.1 fraction) of the possible atom match.
    # the maximum overlap in theory is 12 atoms
    # so that's 1.2 in this case, and the rare atoms should be N and CL
    assert len(starting_configurations) == 2