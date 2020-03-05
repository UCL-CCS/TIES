"""

TIES - test cases for how charges and deletions affect bond generation
The testing with the defined cases by Agastya.

TODO - add test cases which capture which atoms do not match together.
For example, if there is a linking atom C that mutates into O,
Then we might be able to detect that this exact mutation takes place. 1
"""

from topology_superimposer import SuperimposedTopology, get_atoms_bonds_from_ac, \
    superimpose_topologies, _superimpose_topologies, assign_coords_from_pdb
import networkx as nx
import MDAnalysis as mda
from os import path


def load_problem_from_dir(liglig_path):
    """
    Helper function to work with the Agastya's dataset.
    """
    ligand_from, ligand_to = path.basename(liglig_path).split('-')
    hybrid_pair_path = path.join(liglig_path, "hybrid_par")
    print("working now on: ", hybrid_pair_path)
    # fixme - make sure these two are superimposed etc, that data is used later
    mda_left_lig = mda.Universe(path.join(hybrid_pair_path, 'init_%s.pdb' % ligand_from))
    mda_right_lig = mda.Universe(path.join(hybrid_pair_path, 'final_%s.pdb' % ligand_to))

    # read the corresponding charge values for the l14
    leftlig_atoms, leftlig_bonds = get_atoms_bonds_from_ac(path.join(hybrid_pair_path, 'init_%s.ac' % ligand_from))
    rightlig_atoms, rightlig_bonds = get_atoms_bonds_from_ac(path.join(hybrid_pair_path, 'final_%s.ac' % ligand_to))

    # fixme - make sure these two are superimposed etc, that data is used later
    # get the atom location using the .pdb which are superimposed onto each other
    assign_coords_from_pdb(leftlig_atoms, mda_left_lig)
    assign_coords_from_pdb(rightlig_atoms, mda_right_lig)

    # create graphs
    # create the nodes and add edges for one ligand
    # create the nodes and add edges for the other ligand
    ligand1_nodes = {}
    for atomNode in leftlig_atoms:
        ligand1_nodes[atomNode.get_id()] = atomNode
    for nfrom, nto in leftlig_bonds:
        ligand1_nodes[nfrom].bindTo(ligand1_nodes[nto], 'bondType1')

    ligand2_nodes = {}
    for atomNode in rightlig_atoms:
        ligand2_nodes[atomNode.get_id()] = atomNode
    for nfrom, nto in rightlig_bonds:
        ligand2_nodes[nfrom].bindTo(ligand2_nodes[nto], 'bondType1')

    return ligand1_nodes, ligand2_nodes


def test_mcl1_l18l39():
    # Agastya's cases
    liglig_path = "agastya_dataset/mcl1/l18-l39"
    lig1_nodes, lig2_nodes = load_problem_from_dir(liglig_path)
    # we are ignoring the charges by directly calling the superimposer
    suptops = _superimpose_topologies(lig1_nodes.values(), lig2_nodes.values())
    # in this case, there should be only one solution
    assert len(suptops) == 1

    suptop = suptops[0]
    assert len(suptop) == 43

    # check the core atoms of the superimposition for correctness
    # this ensures that the atoms which have no choice (ie they can only be superimposed in one way)
    # are superimposed that way
    core_test_pairs = [('C21', 'C43'), ('C15', 'C37'), ('O3', 'O6'), ('C9', 'C31'),
                  ('C1', 'C23'), ('N1', 'N2'), ('C6', 'C28'), ('C4', 'C26')]
    for atomName1, atomname2 in core_test_pairs:
        assert suptop.contains_atomNamePair(atomName1, atomname2)

    # resolve the multiple possible matches
    avg_dst = suptop.correct_for_coordinates()

    # check if the mirrors were corrected
    corrected_symmetries = [('O2', 'O5'), ('O1', 'O4'), ('H7', 'H25'), ('H6', 'H24'),
                            ('H4', 'H22'), ('H5', 'H23'), ('H8', 'H26'), ('H9', 'H27')]
    for atomName1, atomname2 in corrected_symmetries:
        assert suptop.contains_atomNamePair(atomName1, atomname2)

    # refine against charges
    # ie remove the matches that change due to charge rather than spieces
    all_removed_pairs = suptop.refineAgainstCharges(atol=0.1)
    print(all_removed_pairs)
    removed_pairs = [('C5', 'C27'), ('C4', 'C26')]
    for atomName1, atomname2 in removed_pairs:
        assert not suptop.contains_atomNamePair(atomName1, atomname2)


def test_mcl1_l17l9():
    # Agastya's cases
    liglig_path = "agastya_dataset/mcl1/l17-l9"
    lig1_nodes, lig2_nodes = load_problem_from_dir(liglig_path)

    suptops = _superimpose_topologies(lig1_nodes.values(), lig2_nodes.values())
    assert len(suptops) == 2

    for suptop in suptops:
        # the core chain should always be the same
        core_test_pairs = [('O3', 'O6'), ('C9', 'C28'), ('C3', 'C22'), ('C6', 'C25')]
        for atomName1, atomname2 in core_test_pairs:
            assert suptop.contains_atomNamePair(atomName1, atomname2)

    # use coordinates to solve multiple matches
    [st.correct_for_coordinates() for st in suptops]
    # sort according to the rmsd
    suptops.sort(key=lambda suptop: suptop.rmsd())
    solution_suptop = suptops[0]

    # check the core atoms of the superimposition for correctness
    # this ensures that the atoms which have no choice (ie they can only be superimposed in one way)
    # are superimposed that way
    multchoice_test_pairs = [('C18', 'C37'), ('O2', 'O4'),
                             ('O1', 'O5'), ('H10', 'H28'),
                             ('H9', 'H27'), ('H8', 'H26'),
                             ('H7', 'H25'), ('H6', 'H24'),
                             ('H5', 'H23')]
    for atomName1, atomname2 in multchoice_test_pairs:
        assert solution_suptop.contains_atomNamePair(atomName1, atomname2)

    # refine against charges
    # ie remove the matches that change due to charge rather than spieces
    all_removed_pairs = suptop.refineAgainstCharges(atol=0.1)
    print(all_removed_pairs)
    removed_pairs = [('C14', 'C33'), ('C15', 'C34'), ('C16', 'C35'), ('C17', 'C36')]
    for atomName1, atomname2 in removed_pairs:
        assert not suptop.contains_atomNamePair(atomName1, atomname2)

    # check if the lonely hydrogens were removed together with charges
    removed_lonely_hydrogens = [('H14', 'H32'), ('H11', 'H29')]
    for atomName1, atomname2 in removed_lonely_hydrogens:
        assert not suptop.contains_atomNamePair(atomName1, atomname2)


def test_mcl1_l8l18():
    # Agastya's cases
    liglig_path = "agastya_dataset/mcl1/l8-l18"
    lig1_nodes, lig2_nodes = load_problem_from_dir(liglig_path)

    suptops = _superimpose_topologies(lig1_nodes.values(), lig2_nodes.values())
    assert len(suptops) == 1
    suptop = suptops[0]

    # two rings that not been mutated
    two_untouched_rings = [('C13', 'C35'), ('C16', 'C38'), ('H10', 'H26'), ('C15', 'C37'),
                             ('H15', 'H31'), ('C14', 'C36'), ('H16', 'H32'), ('C17', 'C39'),
                             ('C18', 'C40'), ('C22', 'C44'), ('H14', 'H30'), ('C21', 'C43'),
                             ('H13', 'H29'), ('C20', 'C42'), ('H12', 'H28'), ('C19', 'C41'),
                             ('H11', 'H27')]
    for atomName1, atomname2 in two_untouched_rings[::-1]:
        # the core chain should always be the same
        if suptop.contains_atomNamePair(atomName1, atomname2):
            two_untouched_rings.remove((atomName1, atomname2))
    # the pairs that remain were not found
    assert len(two_untouched_rings) == 0, two_untouched_rings

    # check the linker backbone
    linker_backbone = [('O3', 'O7'), ('C12', 'C34'), ('C11', 'C33'), ('C10', 'C32')]
    for atomName1, atomname2 in linker_backbone[::-1]:
        if suptop.contains_atomNamePair(atomName1, atomname2):
            linker_backbone.remove((atomName1, atomname2))
    assert len(linker_backbone) == 0, linker_backbone

    # check the rings that mutate
    mutating_area = [('O2', 'O6'), ('O1', 'O5'), ('C1', 'C23'), ('C2', 'C24'),
                     ('C9', 'C31'), ('C8', 'C30'), ('C3', 'C25'), ('C4', 'C26'),
                     ('H1', 'H17'), ('C5', 'C27'), ('C6', 'C28'), ('C7', 'C29'),
                     ('H3', 'H19')]
    for atomName1, atomname2 in mutating_area[::-1]:
        if suptop.contains_atomNamePair(atomName1, atomname2):
            mutating_area.remove((atomName1, atomname2))
    assert len(mutating_area) == 0, mutating_area

    # check the linker hydrogens
    linker_hydrogens = [('H5', 'H21'), ('H4', 'H20'),
                        ('H7', 'H23'), ('H6', 'H22'),
                        ('H9', 'H25'), ('H8', 'H24')]
    for atomName1, atomname2 in linker_hydrogens[::-1]:
        if suptop.contains_atomNamePair(atomName1, atomname2):
            linker_hydrogens.remove((atomName1, atomname2))
    assert len(linker_hydrogens) == 0, linker_hydrogens

    # refine against charges
    # ie remove the matches that change due to charge rather than spieces
    removed_pairs = suptop.refineAgainstCharges(atol=0.1)
    # extract the atom names
    removed_atom_names = [(left.atomName, right.atomName) for left, right in removed_pairs]
    should_be_removed = [('C7', 'C29'), ('C6', 'C28'), ('C5', 'C27'),
                           ('C4', 'C26'), ('C3', 'C25'), ('C2', 'C24')]
    assert should_be_removed == removed_atom_names

    # check if the dangling hydrogens were removed
    removed_dangling_hydrogens = [('H3', 'H19'), ('H1', 'H17')]
    for atomName1, atomname2 in removed_dangling_hydrogens:
        assert not suptop.contains_atomNamePair(atomName1, atomname2), (atomName1, atomname2)


def test_mcl1_l32_l42():
    # Agastya's cases
    liglig_path = "agastya_dataset/mcl1/l32-l42"
    lig1_nodes, lig2_nodes = load_problem_from_dir(liglig_path)

    # c17 = next(filter(lambda x: x.atomName == 'C17', lig1_nodes.values()))
    # c38 = next(filter(lambda x: x.atomName == 'C38', lig2_nodes.values()))
    # suptops = _superimpose_topologies(lig1_nodes.values(), lig2_nodes.values(),
    #                                   starting_node_pairs=[(c17, c38)])

    suptops = _superimpose_topologies(lig1_nodes.values(), lig2_nodes.values())
    assert len(suptops) == 1
    suptop = suptops[0]

    # two rings that not been mutated
    two_untouched_rings = [('C5', 'C26'), ('H2', 'H19'), ('C6', 'C27'), ('H3', 'H20'),
                           ('C4', 'C25'), ('H1', 'H18'), ('C7', 'C28'), ('H4', 'H21'),
                           ('C3', 'C24'), ('C8', 'C29'), ('N2', 'N3'), ('H17', 'H35'),
                           ('C2', 'C23'), ('C9', 'C30'), ('C1', 'C22'), ('O1', 'O5'),
                           ('O2', 'O4')]
    for atomName1, atomname2 in two_untouched_rings[::-1]:
        # the core chain should always be the same
        if suptop.contains_atomNamePair(atomName1, atomname2):
            two_untouched_rings.remove((atomName1, atomname2))
    # the pairs that remain were not found
    assert len(two_untouched_rings) == 0, two_untouched_rings

    # check the linker backbone
    linker_backbone = [('O3', 'O6'), ('C12', 'C33'), ('C11', 'C32'), ('C10', 'C31')]
    for atomName1, atomname2 in linker_backbone[::-1]:
        if suptop.contains_atomNamePair(atomName1, atomname2):
            linker_backbone.remove((atomName1, atomname2))
    assert len(linker_backbone) == 0, linker_backbone

    # check the rings that mutate
    mutating_area = [('C13', 'C34'), ('C15', 'C36'), ('H11', 'H28'), ('C14', 'C35'),
                     ('H12', 'H29'), ('C18', 'C39'), ('C17', 'C38'), ('C16', 'C37'),
                     ('H16', 'H34'), ('C21', 'C43'), ('H15', 'H33'), ('C20', 'C42'),
                     ('H14', 'H32'), ('C19', 'C41')]
    for atomName1, atomname2 in mutating_area[::-1]:
        if suptop.contains_atomNamePair(atomName1, atomname2):
            mutating_area.remove((atomName1, atomname2))
    assert len(mutating_area) == 0, mutating_area

    # this pair of hydrogens has actually a different type, so it is not found
    assert not suptop.contains_atomNamePair('H13', 'H31')

    # check the linker hydrogens
    # note that it is 1) not clear which are which, 2) should not matter
    linker_hydrogens = [('H6', 'H23'), ('H5', 'H22'), ('H7', 'H25'), ('H8', 'H24'), ('H10', 'H26')]
    for atomName1, atomname2 in linker_hydrogens[::-1]:
        if suptop.contains_atomNamePair(atomName1, atomname2):
            linker_hydrogens.remove((atomName1, atomname2))
    if len(linker_hydrogens) == 0:
        print(linker_hydrogens)

    # refine against charges
    # ie remove the matches that change due to charge rather than spieces
    removed_pairs, rm_h_pairs = suptop.refineAgainstCharges(atol=0.1)
    print('removed', removed_pairs)
    should_remove_pairs = [('O3', 'O6'), ('C9', 'C30'), ('C21', 'C43'),
        ('C20', 'C42'), ('C19', 'C41'), ('C18', 'C39'), ('C17', 'C38'),
        ('C14', 'C35'), ('C11', 'C32')]
    for n1, n2 in removed_pairs:
        should_remove_pairs.remove((n1.atomName, n2.atomName))
    assert len(should_remove_pairs) == 0, should_remove_pairs

    # check if the lonely hydrogens were removed together with charges
    removed_lonely_hydrogens = [('H15', 'H33'), ('H14', 'H32'),
                   ('H12', 'H29'), ('H8', 'H24'), ('H7', 'H25')]
    for n1, n2 in rm_h_pairs:
        removed_lonely_hydrogens.remove((n1.atomName, n2.atomName))
    assert len(removed_lonely_hydrogens) == 0, removed_lonely_hydrogens



def test_tyk2_l11l14():
    # Agastya's cases
    liglig_path = "agastya_dataset/tyk2/l11-l14"
    lig1_nodes, lig2_nodes = load_problem_from_dir(liglig_path)

    suptop = superimpose_topologies(lig1_nodes.values(), lig2_nodes.values(), atol=0.1)
    assert suptop is not None

    # the core chain should always be the same
    core_test_pairs = [('C4', 'C20'), ('C7', 'C23'), ('O1', 'O3'), ('N1', 'N4'), ('H6', 'H18'), ('C8', 'C24'), ('C9', 'C25'),
     ('H7', 'H19'), ('C10', 'C26'), ('H8', 'H20'), ('N2', 'N5'), ('C11', 'C27'), ('C12', 'C28'), ('H9', 'H21'),
     ('N3', 'N6'), ('H10', 'H22'), ('C13', 'C29'), ('O2', 'O4'), ('C14', 'C30'), ('H1', 'H13')]
    for atomName1, atomname2 in core_test_pairs:
        assert suptop.contains_atomNamePair(atomName1, atomname2)

    # check the core atoms of the superimposition for correctness
    # this ensures that the atoms which have no choice (ie they can only be superimposed in one way)
    # are superimposed that way
    # choice
    multchoice_test_pairs = [('C15', 'C31'), ('H12', 'H24'), ('H11', 'H23'),('CL1', 'CL4'),
                             ('CL3', 'CL5'), ('C3', 'C19'), ('C5', 'C21'),
                             ('C6', 'C22'), ('H5', 'H17'), ('C1', 'C17'),
                             ('C2', 'C18'), ('H4', 'H16'), ('CL3', 'CL5')]
    for atomName1, atomname2 in multchoice_test_pairs:
        assert suptop.contains_atomNamePair(atomName1, atomname2), (atomName1, atomname2)

    # removed due to charges
    assert not suptop.contains_atomNamePair('C16', 'C32')
    # hydrogens should not be dangling by themselves
    assert not suptop.contains_atomNamePair('H2', 'H14')


def test_tyk2_l13l12():
    # Agastya's cases
    liglig_path = "agastya_dataset/tyk2/l13-l12"
    lig1_nodes, lig2_nodes = load_problem_from_dir(liglig_path)

    suptop = superimpose_topologies(lig1_nodes.values(), lig2_nodes.values(), atol=99999)
    assert suptop is not None

    # the core chain should always be the same
    core_test_pairs = [('C4', 'C21'), ('C7', 'C24'), ('O1', 'O3'), ('N1', 'N4'),
                       ('H4', 'H19'), ('C8', 'C25'), ('C12', 'C29'), ('C11', 'C28'),
                       ('N3', 'N6'), ('C13', 'C30'), ('O2', 'O4'), ('C14', 'C31')]
    for atomName1, atomname2 in core_test_pairs[::-1]:
        if suptop.contains_atomNamePair(atomName1, atomname2):
            core_test_pairs.remove((atomName1, atomname2))
    assert len(core_test_pairs) == 0, core_test_pairs

    # there should be only one solution? here is a tricky situation with the ring overlap
    # assert len(suptops) == 1

    # check the core atoms of the superimposition for correctness
    # this ensures that the atoms which have no choice (ie they can only be superimposed in one way)
    # are superimposed that way
    # choice
    # multchoice_test_pairs = [('C15', 'C31'), ('H12', 'H24'), ('H11', 'H23'),('CL1', 'CL4'),
    #                          ('CL3', 'CL5'), ('C3', 'C19'), ('C5', 'C21'),
    #                          ('C6', 'C22'), ('H5', 'H17'), ('C1', 'C17'),
    #                          ('C2', 'C18'), ('H4', 'H16'), ('CL3', 'CL5')]
    # for atomName1, atomname2 in multchoice_test_pairs:
    #     assert suptop.contains_atomNamePair(atomName1, atomname2), (atomName1, atomname2)
    #
    # # removed due to charges
    # assert not suptop.contains_atomNamePair('C16', 'C32')
    # # hydrogens should not be dangling by themselves
    # assert not suptop.contains_atomNamePair('H2', 'H14')
