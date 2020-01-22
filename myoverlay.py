#!/bin/usr/env python3
"""
Overlay two molecule graphs.

Consider / fixme
 - scan the different molecules and find what would be the appropriate absolute/relative tolerance.
 It seems that absolute tolerance would be more meaningful in this case. Relative tolerance has the issue
 that it means a different thing: taking 10% of a very small charge is a very small value, and a very
 small tolerance. Look up the units and figure out how much the actual absolute tolerance should be taken into account.

TODO
 - list all cases in a understanble fashion (
 - before doing any work, check if the topology of ligand/molecule is a strongly connected component. Use networkx with __hash__
 - consider the case where the superimposer finds two components, in which case there is a question about the linking area
 more specifically, whether the area is of the same dimensions. If there is a linking atom that was mutated, then
 it is the atoms that needs to be replaced. However, if there is an insertion of more atoms in the linking area,
 then one has to decide whether the smaller components should appear from scratch. This Shunzhou said could depend
 on whether the component is the crucial anchor for the docking. Another thought about this is that having MD
 should mean that that the simulation should be extended long enough so that the molecule is properly docked. Moreover,
 it is possible that one should reconsider and reapply docking with the new ligand to understand and see if there are other
 types of interactions that can affect the entire process.
 - update naming everywhere to call it a superimposer rather than "overlay"
 - dimension analysis should be carried out. When there is more components, there might be a component that has only
 a single atom. In that case, if the atom is spatially in the same position, that's fine, but if the atom is farther away,
 then it should be marked as a different atom
 - GUI should be used and allow for a manual intervention.
 - LATER: extension to the dimension analysis - one should be able to see if the "same" atom topologically is actually different
 due to its completely different forcefield. Imagine that someone made this atom to be restrained at a different angle,
 this would not be currently noticed. However, how important are these differences?
 - ?? should you consider, if there is a dilemma with how to pick topologies, to use more topologies and compare them?
 - ?? should you consider affinity of atoms

# fixme
"""

import networkx as nx
import MDAnalysis as mda
from os import path
import matplotlib.pyplot as plt
from topology_superimposer import *
import topology_superimposer
# from rdkit import Chem
import glob

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

"""
Good to check:
TYK2, trka, 

Not suitable for testing: 
 - thrombin - a different dir structure
 - TIES_PTP1B - a different dir structure
"""

"""
fixme cases
--mcl1/l18-l39/hybrid_par: bizzare mismatch between different parts, complete misunderstanding, PRIORITY

** Investigate the "symmetry cases" which might require structural superimposition
Complex case l17-l9 in mcl1 which shows a good problem with symmetry:
-  this requires a way of fixing symmetry, particularly when it is possible in both ways, 
is there a simple way to resolve it? no. The two cases actually are pretty much the same case. 
This is because even though they appear like one case: the problem is that they are free to rotate. 
So whatever path you take
- case mcl1/l18-l39/hybrid_par which is similar to the above l17-l9, but acutally this one
has only one solution that is correct, and therefore the symmetry can be corrected,
- case mcl1/l32-l38/hybrid_par another good symmetry case 
- l2-l32/hybrid_par : appears correct but an interesting test for symmetry


Another symmetrical case: 
- imagine that you have A-B-A and you are matching it to A-B, there are two ways in which you can match it. 
now this is actually important because the protein to which the ligand A-B-A -> A-B docks is not symmetrical. 
And therefore removing one face over other means that we might be creating a new problem. 


Tasks:
-TISE - mark the "deleted" items as "not matching due to charge"
-TIES - remove the hydrogens by themselves
"""

proteins_paths = [
                    #"/home/dresio/ucl/dataset/agastya_extracted/tyk2/",
                    #"/home/dresio/ucl/dataset/agastya_extracted/trka/",
                    "/home/dresio/ucl/dataset/agastya_extracted/mcl1/"   # why the errors?
                  ]
for protein_path in proteins_paths:
    for liglig_path in glob.glob(protein_path + 'l*-l*'):
        ligand_from, ligand_to = path.basename(liglig_path).split('-')
        if not (ligand_from == 'l12' and ligand_to == 'l35'):
            continue
        hybrid_pair_path = path.join(liglig_path, "hybrid_par")
        print("working now on: ", hybrid_pair_path)
        # continue
        l11 = mda.Universe(path.join(hybrid_pair_path, 'init_%s.pdb' % ligand_from))
        l14 = mda.Universe(path.join(hybrid_pair_path, 'final_%s.pdb' % ligand_to))

        # read the corresponding charge values for the l14
        l11_atoms, l11_bonds = get_atoms_bonds_from_ac(path.join(hybrid_pair_path, 'init_%s.ac' % ligand_from))
        l14_atoms, l14_bonds = get_atoms_bonds_from_ac(path.join(hybrid_pair_path, 'final_%s.ac' % ligand_to))

        assign_coords_from_pdb(l11_atoms, l11)
        assign_coords_from_pdb(l14_atoms, l14)

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

        # overlay
        print("About to overlay %d atoms with %d atoms" % (len(ligand2_nodes), len(ligand1_nodes)))
        # 0.1 e charge has been used by default: Paper "Rapid, accurate" by Agastya et al (doi: 10.1021/acs.jctc.6b00979)
        si_topologies = superimpose_topologies(ligand1_nodes.values(), ligand2_nodes.values(), atol=0.1, useCharges=False)

        topology_superimposer.verbose_log = False

        print("Found Superimposed Topologies ", len(si_topologies))
        for si_top in si_topologies:
            print("Topology Pairs ", len(si_top.matched_pairs), "Mirror Number", len(si_top.mirrors))
            rmsd, winner, isMirror = si_top.findLowestRmsdMirror()
            print('rmsd', rmsd, 'isMirror', isMirror)

        # print the match
        # for strongly_connected_component in overlays:
        #     print("Strongly Connected Component, length:", len(strongly_connected_component))
        # for atom_from, atom_to in strongly_connected_component:
        #     print('Bound', atom_from.atomName, atom_to.atomName)

        # extract all the unique nodes from the pairs
        all_matched_nodes = set()
        for si_top in si_topologies:
            si_top.print_summary()
            print(si_top.matched_pairs)
            unique_nodes = []
            for pair in si_top.matched_pairs:
                unique_nodes.extend(list(pair))
            all_matched_nodes = all_matched_nodes.union(unique_nodes)

        maintop = si_topologies[0]
        for i, si_top in enumerate(si_top.mirrors, start=1):
            # print only the mismatching pairs
            different = set(si_top.matched_pairs).difference(set(maintop.matched_pairs))
            print('Multiple Choice:', different)

        # extract the atoms that are appearing and disappearing
        # the atom that appears has to be in G2 and not in any of the overlaps
        appearing = [node for node in ligand2_nodes.values() if not node in all_matched_nodes]
        disappearing = [node for node in ligand1_nodes.values() if not node in all_matched_nodes]

        print("disappearing", 'name ' + ' '.join([n.atomName for n in disappearing]))
        print("appearing", 'name ' + ' '.join([n.atomName for n in appearing]))
        # fixme - make it clear which convert to which?

        # fixme - you should check if this is not empty, if you have not found anything, but the molecules
        # has a different number of atoms, etc, then it is wrong, if the molecule is the same, then this also needs
        # a proper message

        # check if you found the correct atoms. They should be a subset of the atoms picked by Agastya in his research
        # agastya_disapp_atom_list_path = path.join(hybrid_pair_path, "disappearing_atoms.txt")
        #
        # if not path.isfile(agastya_disapp_atom_list_path):
        #     print(bcolors.WARNING + "Test file does not exist %s" %agastya_disapp_atom_list_path + bcolors.ENDC)
        # else:
        #     agastya_disapp_atoms = open(agastya_disapp_atom_list_path).read().split()
        #     print(bcolors.OKGREEN + "Agastya's disapp list:" + bcolors.ENDC,
        #           agastya_disapp_atoms)
        #     for disapp_atom in disappearing:
        #         isin = disapp_atom.atomName.lower() in [dis_atom.lower() for dis_atom in agastya_disapp_atoms]
        #         if not isin:
        #             print(bcolors.FAIL + "A disappearing atom not found in Agastya list. Atom: " + bcolors.ENDC,
        #                   disapp_atom.atomName)
        #             #raise Exception("Fail")
        #
        # agastya_app_atom_list_path = path.join(hybrid_pair_path, "appearing_atoms.txt")
        # if not path.isfile(agastya_app_atom_list_path):
        #     print(bcolors.WARNING + "Test file does not exist %s" % agastya_app_atom_list_path + bcolors.ENDC)
        # else:
        #     agastya_app_atoms = open(agastya_app_atom_list_path).read().split()
        #     print(bcolors.OKGREEN+ "Agastya's app list:" + bcolors.ENDC, agastya_app_atoms)
        #     for app_atom in appearing:
        #         isin = app_atom.atomName.lower() in [atom.lower() for atom in agastya_app_atoms]
        #         if not isin:
        #             print(bcolors.FAIL + "An appearing atom not found in Agastya's list. Atom: " + bcolors.ENDC, app_atom.atomName)
        #             #raise Exception("Fail")

        # make the message more positive if they're exactly the same

        """
        In theory you have all the overlays. 
        If some atom is not in any of the overlays, that means that it is a new atom that appears in one structure, 
        but not in another. 
        """