# we start the testing with the defined cases by Agastya,

from topology_superimposer import SuperimposedTopology, get_charges, superimpose_topologies, _superimpose_topologies
import networkx as nx
import MDAnalysis as mda
from os import path


def get_problem(liglig_path):
    ligand_from, ligand_to = path.basename(liglig_path).split('-')
    hybrid_pair_path = path.join(liglig_path, "hybrid_par")
    print("working now on: ", hybrid_pair_path)
    l11 = mda.Universe(path.join(hybrid_pair_path, 'init_%s.pdb' % ligand_from))
    l14 = mda.Universe(path.join(hybrid_pair_path, 'final_%s.pdb' % ligand_to))

    # read the corresponding charge values for the l14
    l11_atoms, l11_bonds = get_charges(path.join(hybrid_pair_path, 'init_%s.ac' % ligand_from))
    l14_atoms, l14_bonds = get_charges(path.join(hybrid_pair_path, 'final_%s.ac' % ligand_to))

    # create graphs
    # create the nodes and add edges for one ligand
    ligand2_nodes = {}
    for atomNode in l14_atoms:
        ligand2_nodes[atomNode.atomId] = atomNode
    for nfrom, nto in l14_bonds:
        ligand2_nodes[nfrom].bindTo(ligand2_nodes[nto])

    # create the nodes and add edges for the other ligand
    ligand1_nodes = {}
    for atomNode in l11_atoms:
        ligand1_nodes[atomNode.atomId] = atomNode
    for nfrom, nto in l11_bonds:
        ligand1_nodes[nfrom].bindTo(ligand1_nodes[nto])

    return ligand1_nodes, ligand2_nodes


def test_mcl1_l18l39_nocharges():
    # Agastya's cases
    liglig_path = "/home/dresio/ucl/dataset/agastya_extracted/mcl1/l18-l39"
    lig1_nodes, lig2_nodes = get_problem(liglig_path)
    # we are ignoring the charges by directly calling the superimposer
    suptops = _superimpose_topologies(lig1_nodes.values(), lig2_nodes.values(), atol=9999)
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

    # resolve the multiple possbile matches
    suptop.correct_for_coordinates()

    # check if the mirrors were corrected
    corrected_symmetries = [('O2', 'O5'), ('O1', 'O4'), ('H7', 'H25'), ('H6', 'H24'),
                            ('H4', 'H22'), ('H5', 'H23'), ('H8', 'H26'), ('H9', 'H27')]
    for atomName1, atomname2 in corrected_symmetries:
        assert suptop.contains_atomNamePair(atomName1, atomname2)
