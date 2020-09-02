
"""
To do:
 - initially, check if each graph is strongly connected
 - rarity rank: rank the atoms so that the atoms that are shared
 in both graphs and which are also rare (e.g. a minimum number of them)
 can be used as starting points for the comparison
 - if you reconstruct an entire graph of one molecule, finish
 this means that the second molecule
 - if you do not reconstruct either of the molecule, that means
 that some of the atoms have been mutated
 - if there is a node sequence A-B-C overlapped with A-X-C, then we
 will finish with subgraphs A and subgraphs C, so there is a possibility
 of finishing with multiple connected subgraphs. In this case,
 the different subgraphs cannot overlap
 - TEST: with multiple graphs at the end, and uncertainty, we should ensure that
 we finish with disjoints graphs: if two graphs are present and they overlap,
 the smaller one should be removed, because it is a subgraph of the larger
 - Low Priority: currently, if the ligands l1 and l2 are symmetric, then the "matching" 
 or superimposition of the topologies
 can be done in multiple ways and there is no difference between them, 
 however, the charges on the molecules can be used to help with the symmetry: ie even
 though they are similar enough (within the absolute tolerance atol), they 
 might on one side be more similar than on the other side. However, here it is ignored
 for now because it does not lead to any wrong results
 -consider an overall try-catch that attaches a message for the user to contact you in the case something
 does not work as expected
 - fixme - you should check if you can rely on atomName and __hash__ for uniqueness which you need
 - ensure that the mirrors are recorded in a better way (rather than as whole suptops?)

 -fixme - make sure that the "bonds" are obeyed, if an atom changes its biding, e.g. it has a double bond and then it has
 one more hydrogen, then that is a very different atom, and both need to correctly evolve into each other

TODO:
 - the _overlayer now always can return only one solution, is there a situation where it cannot?
 - we should demand that the particles have 3D coordinates, and base our work on that

 Improvements:
 - switch to mdanalysis for creating universes and bonds and etc and other stuff,
 - switch to mdanalysis for all kinds of atom traversal - ie try to use atom bonds etc instead of doing
 all these details myself
 
 Case:
 - in the case of a ring, if the ring stops existing, ie one of the bonds is removed, then it is a chain,
 but that crates a big issue, because they entire ring now is completely broken, so it should not be anything similar 
 anymore? that would be detected because the charges change
 
 
 fixme
 - ERROR: imagine that you have a olecule with 3 overlaid Strongly Connected Components that are all the same to
 each other. Currently, their matching will be arbitrary, which could potentially lead to problems
 - ERROR: some overlaps are clearly incorrect: in the case of a chiral molecule, clearly one group disappears,
 and another group appears, however, they are symmetric to each other and therefore the strongly connected component
 is found in that place. To solve this, one could superimpose only the same molecules to understand the global 
 superimposition and see what is the best way to force it, and use that as guidence, however,
 i have the impression that the same thing can be done by taking into the account the position of the gropus in question
 to other components on their own "topologies", but i have to think about this more
 
 - there is a lot of matches which are substandard - possibly. 
 For example, imagine some hydrogen wrongly matched with another hydrogen, 
 These have been used in larger strongly connected components, meaning
 that this pair should be removed, which is the approach I take,  but 
 what if there is no much for one of the hydrogens nowhere? I guess it is new 
 
 ??
 - how to spot errors and poor matches? we should be able to verify our approach and give it some measure.
 if we make this tool available online, then we will see what people give us, and therefore we can work
 on improving the tool. However, we might need hand adjustment - and if something scores poorly, 
 even though it should be perfect, then we might want to flag that. We could use e.g. rdkit 
 to see the "fintertips" which would help us with the maximum score, we could also find the best match,
 using just the atom type 
 - should we take into account the spatial information? what if they can be easily superimposed, 
 and such as structural superimposition would show which atoms are which, how different is that with topology? 
 in theory, we should be able to try to superimpose topologies while ignoring the type (how would we do that? )
 - should we match together the mismatching atoms? we already know which atoms disappear and which appear, but
 do we know which atoms "were mutated" into which others? This would be a useful information. 
 
 TODO - 
 - use the atom type only to create a match to see how different these things are
 - create test cases to see if everything behaves well
 - check if you should understand the number of bonds which is crucial for understanding the species, 
 it might be a better way to distinguish between the atoms than using any atom types etc
 - ensure the molecule is "connected" before doing any work
 
 
 Optimisation:
 - consider sorting your topologies to make it easier to compare them
 - you don't have check every n1 and n2 in the two topologies: 
 e.g. if you find a component which extended to NX-NY match,
 then NX-NY is not a necessary starting condition  
 - one major question that will have to be asked is to understand how after one traversal, we can take away
 certain node pairs and based on that say that these pairs do not need to be used/searched as the starting point
 - we should understand the limits of the search: how many max pairs can we actually find?
 - what you could do is to find "carbons" which do not allow for multiple paths. In other words, if we have A-X-B
 where X is a linear chain, we should be able to start for any molecule on X and arrive at the conclusion,
 for this reason X cannot be on a loop, as the loops are the ones that make the situation more tricky.
 So what I could do is to feed in the topology from both molecules, and find the atoms that are not on loops, and
 use them to traverse the molecule with the hope of finding the right one. However, this is still tricky because there
 will be a lot of atoms that we unnecessarily check.
 - instead of trying every starting AX-BX combination, employ a good gessuing algorithm based on the structural
 overlap, and ensure to terminate searching if you find the right answer
 
 
 Complex Testing - Generation Topologies: 
 - one way to automate testing and make it more robust is to generate molecules yourself knowing 
 which parts are mutating and therefore knowing what parts is the same and how the topology should
 be matched, and therefore 
 - check the number of cycles/benzene rings in both, we should ensure that in both cases this is okay
 
 
 
 - extension: apply the same superimposition but ignore the charges to understand what is happening: we will know
 how many of the atoms are changed due to charge, and what happens to the overlap. Furthermore, we can use it 
 as our template. For example, by having the fully superimposed structure, we would be able to assign to each
 node its "universal position" and use that to make some of the later decisions. Furthermore, this
 global match can be extended, and once having multiple components, we can see if we can "bridge" them over the
 areas that are not matched. For example, in the case of a mutated atom, we can jump over that one atom
 and in both cases we know how we are connected to the other component. In other words, connecting the components
 together. 
 

fixme
- approach to rare cases where "mirror" is not correct, but it is clearly a mirror. Mirror is checked naively,
which does not always work. The case where it does not work is test_mcl1_l17l9. This is because
one of the atoms used in one of the overlaps is not used in the "mirror". However,
we can score the two topolgies, and notice which one is better based on the RMSD.
- the above can extended to a further search, if we have coordinates, we can compute the distance for atoms and
check exactly which atoms do not fit with each other. If other topologies matchd the same area in a better way,
and have a lower RMSD or in general, a better match in this case, then we would know and be able to specifically
say that this is some kind of "unusual mirror".
 
fixme 
 - optimising too much - finding matches which should not be found? how come our topology superimposition
 ignores the chirality in the case of mcl1?

Suggestions:
 - should you keep track of which bonds are rotatable? just like FEP+ does?
 - growing RMSD: another way to grow the molecules and check local RMSD is to superimpose a smaller part of it,
 ie are these two atoms matching? Maybe the whole superimpositon should not decided about that,
 but rather, the local environment
"""
import hashlib
import copy
import itertools
import sys
import os
import warnings
from io import StringIO
from functools import reduce

import numpy as np
import networkx as nx

import MDAnalysis as mda
from MDAnalysis.analysis.distances import distance_array
from MDAnalysis.analysis.align import rotation_matrix

# fixme: change to CA:C dictionary
# fixme delete this for now
general_atom_types = {
    'C' : {'C', 'CA', 'CB', 'C3', 'CX'},
    'CL' : {'CL'},
    'H' : {'H', 'HA', 'HN', 'H4', 'HC', 'H1', 'HO'},
    'O' : {'O', 'OH'},
    'N' : {'N', 'NB', 'NS'},
}
# fixme - change this to a simpler version
general_atom_types2 = {
    # Source http://ambermd.org/antechamber/gaff.html#atomtype
    # However, some general atom types were added separately
    'C': 'C', 'CA': 'C', 'CB': 'C', 'C3': 'C', 'CX': 'C', 'C1': 'C', 'C2': 'C', 'CC': 'C',
    'CD': 'C', 'CE': 'C', 'CF': 'C', 'CP': 'C', 'CQ': 'C', 'CU': 'C', 'CV': 'C', 'CY': 'C',
    'H': 'H', 'HA':'H', 'HN': 'H', 'H4':'H', 'HC':'H', 'H1':'H', 'HX': 'H',
    'HO':'H', 'HS':'H', 'HP':'H',  'H2':'H', 'H3': 'H',  'H5':'H',
    'P2': 'P', 'P3': 'P', 'P4': 'P', 'P5': 'P', 'PB': 'P', 'PC': 'P',
    'PD': 'P', 'PE': 'P', 'PF': 'P', 'PX': 'P', 'PY': 'P',
    'O': 'O', 'OH':'O', 'OS':'O',
    'N': 'N', 'NB':'N', 'NS':'N', 'N1':'N', 'N2':'N', 'N3':'N',
    'N4':'N', 'NA':'N', 'NH':'N', 'NO':'N', 'NC': 'N',  'ND':'N',
    'S2': 'S', 'SH': 'S', 'SS': 'S', 'S4': 'S', 'S6': 'S', 'SX': 'S', 'SY': 'S',
    'CL': 'CL',
    'F': 'F',
    'BR': 'BR', 'B': 'BR',
    'I': 'I',
    'S': 'S',
}

class AtomNode:
    counter = 1
    def __init__(self, name, type, charge=None, use_general_type=False):
        self.atomId = None
        # this atom name might change
        self.atomName = name.upper()
        self.originalAtomName = self.atomName
        self.resname = None
        self.resId = None
        self.charge = charge
        self.type = type.upper()
        self.bonds = set()
        self.use_general_type = use_general_type

        # save the general type
        self.gentype = general_atom_types2[self.type]

        self.unique_counter = AtomNode.counter
        AtomNode.counter += 1

        self.hash_value = None

    def set_atomName(self, atomName):
        self.atomName = atomName

    def set_id(self, id):
        self.atomId = id

    def get_id(self):
        return self.atomId

    def set_resname(self, resname):
        self.resname = resname

    def set_resid(self, res_ID):
        self.resId = res_ID

    def set_charge(self, charge):
        self.charge = charge

    def set_original_charge(self, charge):
        # this charge should never change and should only be
        # used once.  This is for a test.
        # ie that appearing and disappearing atoms keep their charges intact.
        self.original_charge = charge

    def set_type(self, amber_type):
        self.amber_type = amber_type

    def set_position(self, x, y, z):
        corrected_type = np.array([x, y, z], dtype='float32')
        self.position = corrected_type

    def isHydrogen(self):
        if self.type == 'H':
            return True

        return False

    def isBoundTo(self, other):
        for atom, bond_type in self.bonds:
            if atom is other:
                return True

        return False

    def __hash__(self):
        # Compute the hash key once
        if self.hash_value is not None:
            return self.hash_value

        m = hashlib.md5()
        # fixme - ensure that each node is characterised by its chemical info,
        # fixme - the atomId might not be unique, so check before the input data
        m.update(str(self.charge).encode('utf-8'))
        m.update(str(self.unique_counter).encode('utf-8'))
        # so include the number of bonds which is basically an atom type
        m.update(str(len(self.bonds)).encode('utf-8'))
        self.hash_value = int(m.hexdigest(), 16)

        return self.hash_value

    def __str__(self):
        return self.atomName

    def __repr__(self):
        return self.atomName

    def bindTo(self, other, btype):
        self.bonds.add((other, btype))
        other.bonds.add((self, btype))

    def eq(self, atom, atol=0):
        """
        What does it mean that two atoms are the same? They are the same type and charge.
        5 % tolerance by default
        """
        if self.type == atom.type and \
                np.isclose(self.charge, atom.charge, atol=atol):
            return True

        return False

    def sameGenType(self, atom):
        # fixme - this one or other? some other configuration system?
        if self.use_general_type:
            if self.gentype == atom.gentype:
                return True
            return False

        if self.type == atom.type.upper():
            return True

        return False

    def sameExactType(self, atom):
        if self.type == atom.type.upper():
            return True

        return False

    def __deepcopy__(self, memodict={}):
        # https://stackoverflow.com/questions/1500718/how-to-override-the-copy-deepcopy-operations-for-a-python-object
        # it is a shallow copy, as this object is "immutable"
        return self

def are_same_general_type(type1, type2):
    # check if the two atom types are of the same general type
    # ie C and CB are both carbons
    for atom_type in general_atom_types:
        if type1 in atom_type:
            if type2 in atom_type:
                return True
            return False
    raise Exception

class AtomPair:
    """
    An atom pair for networkx.
    """

    def __init__(self, left_node, right_node):
        self.left_node = left_node
        self.right_node = right_node
        # generate the hash value for this match
        self.hash_value = self._gen_hash()

    def _gen_hash(self):
        m = hashlib.md5()
        m.update(str(hash(self.left_node)).encode('utf-8'))
        m.update(str(hash(self.right_node)).encode('utf-8'))
        return int(m.hexdigest(), 16)

    def __hash__(self):
        return self.hash_value

    def is_pair(self, old_pair_tuple):
        if self.left_node is old_pair_tuple[0] and self.right_node is old_pair_tuple[1]:
            return True

        return False


class SuperimposedTopology:
    """
    SuperimposedTopology contains in the minimal case two sets of nodes S1 and S2, which
    are paired together and represent a strongly connected component.

    However, it can also represent the symmetrical versions that were superimposed.
    """

    def __init__(self, topology1=None, topology2=None, mdaL=None, mdaR=None):
        self.mda_ligandL = mdaL
        self.mda_ligandR = mdaR

        self.can_use_mda = True
        if self.mda_ligandL is None or self.mda_ligandR is None:
            self.can_use_mda = False

        """
        @superimposed_nodes : a set of pairs of nodes that matched together
        """
        matched_pairs = []

        # TEST: with the list of matching nodes, check if each node was used only once,
        # the number of unique nodes should be equivalent to 2*len(common_pairs)
        all_matched_nodes = []
        [all_matched_nodes.extend(list(pair)) for pair in matched_pairs]
        assert len(matched_pairs) * 2 == len(all_matched_nodes)

        # fixme don't allow for initiating with matche pairs, it's not used anyway

        # todo convert to nx? some other graph theory package?
        matched_pairs.sort(key=lambda pair: pair[0].atomName)
        self.matched_pairs = matched_pairs
        self.top1 = topology1
        self.top2 = topology2
        # create graph representation for both in networkx library, initially to track the number of cycles
        #fixme

        self.mirrors = []
        self.alternative_mappings = []
        # this is a set of all nodes rather than their pairs
        self.nodes = set(all_matched_nodes)
        self.nodes_added_log = []

        self.internal_ids = None
        self.unique_atom_count = 0
        self.matched_pairs_bonds = {}

        # removed because
        # fixme - make this into a list
        self._removed_due_to_charge = []    # atom-atom charge decided by qtol
        self._removed_because_disjointed_cc = []    # disjointed segment
        self._removed_due_to_net_charge = []
        self._removed_because_unmatched_rings = []

        # save the cycles in the left and right molecules
        if self.top1 is not None and self.top2 is not None:
            self._init_nonoverlapping_cycles()

    def _init_nonoverlapping_cycles(self):
        """
        Compile the cycles separately for the left and right molecule.
        Then, across the cycles, remove the nodes that join rings (double rings).
        """
        lcycles, rcycles = self.getOriginalCircles()
        # remove any nodes that are shared between two cycles
        for c1, c2 in itertools.combinations(lcycles, r=2):
            common = c1.intersection(c2)
            for atom in common:
                c1.remove(atom)
                c2.remove(atom)

        # same for rcycles
        for c1, c2 in itertools.combinations(rcycles, r=2):
            common = c1.intersection(c2)
            for atom in common:
                c1.remove(atom)
                c2.remove(atom)

        self._nonoverlapping_lcycles = lcycles
        self._nonoverlapping_rcycles = rcycles

    def get_single_topology_region(self):
        # get the atoms which were unmatched for any reason
        # fixme: this should not work with disjointed cc and others?
        unmatched_due_to_reasons = self._removed_because_disjointed_cc + \
                          self._removed_due_to_net_charge + \
                          self._removed_due_to_charge

        # combine with the main matched area
        return self.matched_pairs + unmatched_due_to_reasons

    def get_single_topology_app(self):
        # get the apperaing and disappearing region in the hybrid single topology
        # use the single topology region and classify all other atoms not in it
        # ass either apperaing or disapearing
        single_top_area = self.get_single_topology_region()
        # turn it into a set
        single_top_set = set()
        for l, r in single_top_area:
            single_top_set.add(l)
            single_top_set.add(r)

        # these unmatched atoms could be due to charge etc.
        # so they historically refer to the dual-topology
        unmatched_app = self.get_appearing_atoms()
        app = {a for a in unmatched_app if a not in single_top_set}
        unmatched_dis = self.get_disappearing_atoms()
        dis = {a for a in unmatched_dis if a not in single_top_set}

        return app, dis

    def is_or_was_matched(self, atomName1, atomName2):
        """
        A helper function. For whatever reasons atoms get discarded.
        E.g. they had a different charge, or were part of the disjointed component, etc.
        This function simply checks if the most original match was made between the two atoms.
        It helps with verifying the original matching.
        """
        if self.contains_atomNamePair(atomName1, atomName2):
            return True

        # check if it was unmatched
        unmatched_lists = [self._removed_because_disjointed_cc,
                           self._removed_due_to_net_charge,
                           self._removed_due_to_charge]
        for unmatched_list in unmatched_lists:
            for atom1, atom2 in unmatched_list:
                if atom1.atomName == atomName1 and atom2.atomName == atomName2:
                    return True

        return False

    def find_pair_with_atom(self, atomName):
        for node1, node2 in self.matched_pairs:
            if node1.atomName == atomName or node2.atomName == atomName:
                return node1, node2

        return None

    def get_unmatched_atoms(self):
        """
        Find the atoms in both topologies which were unmatched and return them.
        These are both, appearing and disappearing.

        Note that some atoms were removed due to charges.
        """
        unmatched_atoms = []
        for node in self.top1:
            if not self.contains_node(node):
                unmatched_atoms.append(node)

        for node in self.top2:
            if not self.contains_node(node):
                unmatched_atoms.append(node)

        return unmatched_atoms

    def getUnqiueAtomCount(self):
        """
        Requires that the .assign_atoms_ids() was called.
        This should be rewritten. But basically, it needs to count each matched pair as one atom,
        and the apperaing and disappearing atoms separately.
        """
        return self.unique_atom_count

    def alignLigandsUsingMatched(self, overwrite_original=False):
        # return self.rmsd()
        """
        Align the two ligands using the matched area.
        Note: we assume that the left ligand is docked. The left ligand is the reference here.
        fixme -consider the use the aligning/RMSD could be using in scoring of symmetrical suptops.
        Q: is aligning a good idea to deal with the symmetry issue?

        Use MDAnalysis for this.
        #fixme - in the future, all the work should be carried out on MDAnalysis?
        """
        if not self.can_use_mda:
            # cannot use MDA for aligning the ligands.
            # Simply return the rmsd of the current setting instead
            print('Will not use mda positions for aligning. Simply taking .rmsd()')
            return self.rmsd()

        # extract the IDs and use them to pick the atoms in MDAnalysis
        # note that the order matters
        matched_l_ids = [l.atomId for l, r in self.matched_pairs]
        matched_r_ids = [r.atomId for l, r in self.matched_pairs]

        # save the original positions (deep copy)
        original_left_pos = np.empty_like(self.mda_ligandL.atoms.positions)
        np.copyto(dst=original_left_pos, src=self.mda_ligandL.atoms.positions)
        original_right_pos = np.empty_like(self.mda_ligandR.atoms.positions)
        np.copyto(dst=original_right_pos, src=self.mda_ligandR.atoms.positions)

        # select the same atoms in MDAnalysis,
        # select separately to keep the order correct
        selection_ids_l = ['bynum ' + str(id) for id in matched_l_ids]
        mda_l_matched = self.mda_ligandL.select_atoms(*selection_ids_l)
        selection_ids_r = ['bynum ' + str(id) for id in matched_r_ids]
        mda_r_matched = self.mda_ligandR.select_atoms(*selection_ids_r)

        # translate all atoms to the origin of the matched subcomponent
        # using Centre-of-geometry of the matched area
        left_ligand_original_cog = mda_l_matched.center_of_geometry()
        self.mda_ligandL.atoms.translate(-left_ligand_original_cog)
        right_ligand_original_cog = mda_r_matched.center_of_geometry()
        self.mda_ligandR.atoms.translate(-right_ligand_original_cog)

        # set the right ligand as reference/mobile
        if self.left_coords_are_ref:
            ref = mda_l_matched
            mob = mda_r_matched
            mob_ligand = self.mda_ligandR
            ref_cog = left_ligand_original_cog
        else:
            ref = mda_r_matched
            mob = mda_l_matched
            mob_ligand = self.mda_ligandL
            ref_cog = right_ligand_original_cog

        # apply the rotation matrix to all atoms to match the mobile ligand to the ref ligand
        # if they are already superimposed, set rmsd to 0
        if np.all(mob.positions == ref.positions):
            rmsd = 0
        else:
            R, rmsd = rotation_matrix(mob.positions, ref.positions)
            mob_ligand.atoms.rotate(R)

        # the new cog is that of the ref
        self.mda_ligandL.atoms.translate(ref_cog)
        self.mda_ligandR.atoms.translate(ref_cog)

        # update the atoms with the mapping done via IDs
        # for the left
        if overwrite_original:
            for mda_a in self.mda_ligandL.atoms:
                found = False
                for loaded_a in self.top1:
                    if mda_a.id == loaded_a.atomId:
                        loaded_a.set_position(mda_a.position[0], mda_a.position[1], mda_a.position[2])
                        found = True
                        break
                assert found
            # and for the right
            for mda_a in self.mda_ligandR.atoms:
                found = False
                for loaded_a in self.top2:
                    if mda_a.id == loaded_a.atomId:
                        loaded_a.set_position(mda_a.position[0], mda_a.position[1], mda_a.position[2])
                        found = True
                        break
                assert found

        # put back original pos
        self.mda_ligandL.atoms.positions = original_left_pos
        self.mda_ligandR.atoms.positions = original_right_pos

        if rmsd == None:
            #fixme ? why does it return None?
            return 9999999999

        return rmsd

    def getDualTopologyBonds(self):
        """
        Get the bonds between all the atoms.
        Use the atom IDs for the bonds.
        """
        assert self.top1 is not None and self.top2 is not None
        # fixme - check if the atoms IDs have been generated
        assert self.internal_ids is not None

        # extract the bonds for the matched molecules first
        bonds = set()
        for from_pair, bonded_pair_list in self.matched_pairs_bonds.items():
            from_pair_id = self.get_generated_atom_ID(from_pair)
            for bonded_pair, bond_type in bonded_pair_list:
                if not self.ignore_bond_types:
                    assert bond_type[0] == bond_type[1]
                # every bonded pair has to be in the topology
                assert bonded_pair in self.matched_pairs
                to_pair_id = self.get_generated_atom_ID(bonded_pair)
                # before adding them to bonds, check if they are not already there
                bond_sorted = sorted([from_pair_id, to_pair_id])
                bond_sorted.append(bond_type[0])
                bonds.add(tuple(bond_sorted))

        # extract the bond information from the unmatched
        unmatched_atoms = self.get_unmatched_atoms()
        # for every atom, check to which "pair" the bond connects,
        # and use that pair's ID to make the link

        # several iterations of walking through the atoms,
        # this is to ensure that we remove each atom one by one
        # e.g. imagine this PAIR-SingleA1-SingleA2-SingleA3
        # so only the first SingleA1 is connected to a pair,
        # so the first iteration would take care of that,
        # the next iteration would connect SingleA2 to SingleA1, etc
        # first, remove the atoms that are connected to pairs
        for atom in unmatched_atoms:
            for bonded_atom, bond_type in atom.bonds:
                unmatched_atom_id = self.get_generated_atom_ID(atom)
                # check if the unmatched atom is bonded to any pair
                pair = self.find_pair_with_atom(bonded_atom.atomName)
                if pair is not None:
                    # this atom is bound to a pair, so add the bond to the pair
                    pairId = self.get_generated_atom_ID(pair[0])
                    # add the bond between the atom and the pair
                    bond_sorted = sorted([unmatched_atom_id, pairId])
                    bond_sorted.append(bond_type[0])
                    bonds.add(tuple(bond_sorted))
                    linked_to_matched_pair = True
                else:
                    # it is not directly linked to a matched pair,
                    # simply add this missing bond to whatever atom it is bound
                    another_unmatched_atom_id = self.get_generated_atom_ID(bonded_atom)
                    bond_sorted = sorted([unmatched_atom_id, another_unmatched_atom_id])
                    bond_sorted.append(bond_type[0])
                    bonds.add(tuple(bond_sorted))

        # fixme - what about circles etc? these bonds
        # that form circles should probably be added while checking if the circles make sense etc
        # also, rather than checking if it is a circle, we could check if the new linked atom,
        # is in a pair to which the new pair refers (the same rule that is used currently)
        return bonds

    def only_largest_CC_survives(self):
        """
        CC - Connected Component.

        Removes any disjoint components. Only the largest CC will be left.
        In the case of of equal length CCs, an arbitrary is chosen.

        How:
        Generates the graph where each pair is a single node, connecting the nodes if the bonds exist.
        Uses then networkx to find CCs.
        """
        def lookup_up(pairs, tuple_pair):
            for pair in pairs:
                if pair.is_pair(tuple_pair):
                    return pair

            raise Exception('Did not find the AtomPair')

        G = nx.Graph()
        atom_pairs = []
        for pair in self.matched_pairs:
            ap = AtomPair(pair[0], pair[1])
            atom_pairs.append(ap)
            G.add_node(ap)

        # connect the atom pairs
        for pair_from, pair_list in self.matched_pairs_bonds.items():
            # lookup the corresponding atom pairs
            ap_from = lookup_up(atom_pairs, pair_from)
            for tuple_pair, bond_type in pair_list:
                ap_to = lookup_up(atom_pairs, tuple_pair)
                G.add_edge(ap_from, ap_to)

        # check for connected comoponents (CC)
        remove_ccs = []
        ccs = list(nx.connected_components(G))
        largest_cc = max([len(cc) for cc in ccs])
        if len(ccs) > 1:
            # there are disjoint fragments, remove the smaller one
            for cc in ccs[::-1]:
                # remove the cc if it smaller than the largest component
                if len(cc) < largest_cc:
                    remove_ccs.append(cc)
                    ccs.remove(cc)

            if len(ccs) > 1:
                # there are equally large CCs
                print("The Connected Componnets are equally large! Picking the first one")
                for cc in ccs[1:]:
                    remove_ccs.append(cc)
                    ccs.remove(cc)

            assert len(ccs) == 1, "At this point there should be left only one main component"

        # remove the smaller ccs
        for cc in remove_ccs:
            for atom_pair in cc:
                atom_tuple = (atom_pair.left_node, atom_pair.right_node)
                self.remove_node_pair(atom_tuple)
                self._removed_because_disjointed_cc.append(atom_tuple)

        return largest_cc, remove_ccs

    def get_node(self, atom_id):
        for node in self.nodes:
            if node.atomName == atom_id:
                return node

        return None

    # fixme - id_start should be a property of the suptop
    def assign_atoms_ids(self, id_start=1):
        """
        Assign an ID to each pair A1-B1. This means that if we request an atom ID
        for A1 or B1 it will be the same.

        Then assign different IDs for the other atoms
        """
        self.internal_ids = {}
        id_counter = id_start
        # for each pair assign an ID
        for left_atom, right_atom in self.matched_pairs:
            self.internal_ids[left_atom] = id_counter
            self.internal_ids[right_atom] = id_counter
            # make it possible to look up the atom ID with a pair
            self.internal_ids[(left_atom, right_atom)] = id_counter

            id_counter += 1
            self.unique_atom_count += 1

        # for each atom that was not mapped to any other atom,
        # but is still in the topology, generate an ID for it

        # find the not mapped atoms in the left topology and assign them an atom ID
        for node in self.top1:
            # check if this node was matched
            if not self.contains_node(node):
                self.internal_ids[node] = id_counter
                id_counter += 1
                self.unique_atom_count += 1

        # find the not mapped atoms in the right topology and assign them an atom ID
        for node in self.top2:
            # check if this node was matched
            if not self.contains_node(node):
                self.internal_ids[node] = id_counter
                id_counter += 1
                self.unique_atom_count += 1

        # return the last atom
        return id_counter

    def get_generated_atom_ID(self, atom):
        return self.internal_ids[atom]

    def get_appearing_atoms(self):
        """
        # fixme - should check first if atomName is unique
        Return a list of appearing atoms (atomName) which are the
        atoms that are
        """
        unmatched = []
        for top2_atom in self.top2:
            is_matched = False
            for _, matched_right_ligand_atom in self.matched_pairs:
                if top2_atom is matched_right_ligand_atom:
                    is_matched = True
                    break
            if not is_matched:
                unmatched.append(top2_atom)

        return unmatched

    def get_disappearing_atoms(self):
        """
        # fixme - should check first if atomName is unique
        # fixme - update to using the node set
        Return a list of appearing atoms (atomName) which are the
        atoms that are found in the topology, and that
        are not present in the matched_pairs
        """
        unmatched = []
        for top1_atom in self.top1:
            is_matched = False
            for matched_left_ligand_atom, _ in self.matched_pairs:
                if top1_atom is matched_left_ligand_atom:
                    is_matched = True
                    break
            if not is_matched:
                unmatched.append(top1_atom)

        return unmatched

    def remove_lonely_hydrogens(self):
        """
        You could also remove the hydrogens when you correct charges.
        """
        raise Exception('Not Implemented')
        # in order to see any hydrogens that are by themselves, we check for any connection
        removed_pairs = []
        for A1, B1 in self.matched_pairs:
            # fixme - assumes hydrogens start their names with H*
            if not A1.atomName.upper().startswith('H'):
                continue

            # check if any of the bonded atoms can be found in this sup top
            if not self.contains_any_node(A1.bonds) or not self.contains_node(B1.bonds):
                # we appear disconnected, remove us
                pass
            for bonded_atom in A1.bonds:
                assert not bonded_atom.atomName.upper().startswith('H')
                if self.contains_node(bonded_atom):
                    continue

        return removed_pairs

    def __len__(self):
        return len(self.matched_pairs)

    def __repr__(self):
        return str(len(self.matched_pairs)) + ":" + ', '.join([a.atomName + '-' + b.atomName for a,b in self.matched_pairs])

    def set_tops(self, top1, top2):
        self.top1 = top1
        self.top2 = top2

    def set_mdas(self, ligandLmda, ligandRmda):
        self.mda_ligandL = ligandLmda
        self.mda_ligandR = ligandRmda

    def print_summary(self):
        print("Topology Pairs ", len(self.matched_pairs), "Mirror Number", len(self.mirrors))

        # print the match
        # for strongly_connected_component in overlays:
        #     print("Strongly Connected Component, length:", len(strongly_connected_component))
        # for atom_from, atom_to in strongly_connected_component:
        #     print('Bound', atom_from.atomName, atom_to.atomName)

        # extract all the unique nodes from the pairs
        all_matched_nodes = set()
        print("VMD Superimposed topology: len %d :" % len(self.matched_pairs),
              'name ' + ' '.join([node1.atomName.upper() for node1, _ in self.matched_pairs]),
              '\nto\n',
              'name ' + ' '.join([node2.atomName.upper() for _, node2 in self.matched_pairs]))
        print("PYMOL Superimposed topology: len %d :" % len(self.matched_pairs),
              'sel left, name ' + '+'.join([node1.atomName.upper() for node1, _ in self.matched_pairs]),
              '\nto\n',
              'sel right, name ' + '+'.join([node2.atomName.upper() for _, node2 in self.matched_pairs]))
        print(', '.join([a.atomName + '-' + b.atomName for a,b in self.matched_pairs]))
        print("Creation Order: ", self.nodes_added_log)
        unique_nodes = []
        for pair in self.matched_pairs:
            unique_nodes.extend(list(pair))
        all_matched_nodes = all_matched_nodes.union(unique_nodes)

        for i, si_top in enumerate(self.mirrors, start=1):
            print('Mirror:', i)
            # print only the mismatching pairs
            different = set(si_top.matched_pairs).difference(set(self.matched_pairs))
            print(different)

    def matched_atom_types_are_the_same(self):
        # in order to get the best superimposition, the algorithm will rely on the
        # general atom type. Ie CA and CD might be matched together to maximise
        # the size of the superimposition.
        # This function removes atom types that are not exactly the same.
        for a1, a2 in self.matched_pairs[::-1]:
            if not a1.sameExactType(a2):
                # remove this pair now. It served its purpose to get the best superimposition.
                # but the atoms "might" have been mutated.
                self.remove_node_pair((a1, a2))
                print(f'Removed earlier matched general type:{a1}-{a2}')

    def get_net_charge(self):
        """
        Calculate the net charge difference across
        the matched pairs.
        """
        return sum(n1.charge - n2.charge for n1, n2 in self.matched_pairs)

    def get_worst_match_charge(self):
        """
        Returns the largest difference in charge found in the pairs.
        """
        return max(np.abs(n1.charge - n2.charge) for n1, n2 in self.matched_pairs)

    def remove_worst_charge_match(self):
        """
        Find a match for which the charge difference is the worst. Then remove it.
        If there is no charge differences, return 0.
        Otherwise, return the charge difference of the removed pair.
        """
        largest_difference = self.get_worst_match_charge()
        if largest_difference == 0:
            return 0

        # find the pair with the largest difference
        worst_match = [(n1, n2) for n1, n2 in self.matched_pairs if np.abs(n1.charge - n2.charge) == largest_difference][0]
        self.remove_node_pair(worst_match)
        # add to the list of removed because of the net charge
        self._removed_due_to_net_charge.append(worst_match)
        return np.abs(largest_difference)

    def remove_node_pair(self, node_pair):
        assert len(node_pair) == 2
        # remove the pair
        self.matched_pairs.remove(node_pair)
        # remove from the current set
        self.nodes.remove(node_pair[0])
        self.nodes.remove(node_pair[1])

        # update the log
        self.nodes_added_log.append(("Removed", node_pair))

        # get bonds to and from this pair
        bound_pairs = self.matched_pairs_bonds[node_pair]
        del self.matched_pairs_bonds[node_pair]

        # make sure any
        for bound_pair, bond_type in bound_pairs:
            if bond_type[0] != bond_type[1]:
                # fixme - this requires more attention
                log('While removing a pair noticed that it has a different bond type. ')
            # remove their binding to the removed pair
            bound_pair_bonds = self.matched_pairs_bonds[bound_pair]
            bound_pair_bonds.remove((node_pair, bond_type))

    def remove_attached_hydrogens(self, node_pair):
        """
        This means that the node_pair, to which these hydrogens are attached, was removed.
        Therefore these are dangling hydrogens.

        Check if these hydrogen are matched/superimposed. If that is the case. Remove the pairs.

        Note that if the hydrogens are paired and attached to node_pairA,
        they have to be attached to node_pairB, as a rule of being a match.
        """
        attached_pairs = self.matched_pairs_bonds[node_pair]

        removed_pairs = []
        for pair, bond_types in list(attached_pairs):
            # ignore non hydrogens
            if not pair[0].atomName.startswith('H'):
                continue

            self.remove_node_pair(pair)
            log('Removed attached hydrogen pair: ', pair)
            removed_pairs.append(pair)
        return removed_pairs

    def findLowestRmsdMirror(self):
        """
        Walk through the different mirrors and out of all options select the one
        that has the lowest RMSD. This way we increase the chance of getting a better match.
        However, long term it will be necessary to use the dihedrals to ensure that we match
        the atoms better.
        """
        # fixme - you have to also take into account the "weird / other symmetires" besdies mirrors
        winner = self
        lowest_rmsd = self.rmsd()
        for mirror in self.mirrors:
            mirror_rmsd = mirror.rmsd()
            if mirror_rmsd < lowest_rmsd:
                lowest_rmsd = mirror_rmsd
                winner = mirror

        if self is winner:
            # False here means that it is not a mirror
            return lowest_rmsd, self, False
        else:
            return lowest_rmsd, winner, True

    def is_subgraph_of_global_top(self):
        """
        Check if after superimposition, one graph is a subgraph of another
        :return:
        """
        # check if one topology is a subgraph of another topology
        if len(self.matched_pairs) == len(self.top1) or len(self.matched_pairs) == len(self.top2):
            log("This graph is a equivalent to the topology")
            return True

        return False

    def rmsd(self):
        """
        For each pair take the distance, and then get rmsd, so root(mean(square(devisation)))
        """

        assert len(self.matched_pairs) > 0

        sq_dsts = []
        for nodeA, nodeB in self.matched_pairs:
            dst = distance_array(np.array([nodeA.position, ]), np.array([nodeB.position, ]))[0]
            sq_dsts.append(dst**2)
        return np.sqrt(np.mean(sq_dsts))

    def add_node_pair(self, node_pair):
        # Argument: bonds are most often used to for parent, but it is a
        # set of "matched pairs"

        # fixme - use this function in the __init__ to initialise
        assert not node_pair in self.matched_pairs, 'already added'
        # check if a1 or a2 was used before
        for a1, a2 in self.matched_pairs:
            if node_pair[0] is a1 and node_pair[1] is a2:
                raise Exception('already exists')
        self.matched_pairs.append(node_pair)
        self.matched_pairs.sort(key=lambda pair: pair[0].atomName)
        # update the list of unique nodes
        n1, n2 = node_pair
        assert not n1 in self.nodes and not n2 in self.nodes, (n1, n2)
        self.nodes.add(n1)
        self.nodes.add(n2)
        assert len(self.matched_pairs) * 2 == len(self.nodes)

        # update the log to understand the order in which this sup top was created
        self.nodes_added_log.append(("Added", node_pair))

        # --------------------
        # update the bond information
        # create a list of bonds for this pair
        self.matched_pairs_bonds[node_pair] = set()

    def link_pairs(self, from_pair, pairs):
        """
        This helps take care of the bonds.
        """
        assert from_pair in self.matched_pairs_bonds
        for pair, bond_types in pairs:
            # the parent pair should have its list of pairs
            assert pair in self.matched_pairs_bonds, f'not found pair {pair}'

            # link X-Y
            self.matched_pairs_bonds[from_pair].add((pair, bond_types))
            # link Y-X
            self.matched_pairs_bonds[pair].add((from_pair, bond_types))

    def link_with_parent(self, pair, parent, type):
        assert len(pair) == 2
        assert len(parent) == 2

        # the parent pair should have its list of pairs
        assert pair in self.matched_pairs_bonds, f'not found pair {pair}'
        assert parent in self.matched_pairs_bonds

        # link X-Y
        self.matched_pairs_bonds[parent].add((pair, type))
        # link Y-X
        self.matched_pairs_bonds[pair].add((parent, type))

    def __copy__(self):
        # https://stackoverflow.com/questions/1500718/how-to-override-the-copy-deepcopy-operations-for-a-python-object
        newone = type(self)()
        newone.__dict__.update(self.__dict__)

        # make a shallow copy of the arrays
        newone.matched_pairs = copy.copy(self.matched_pairs)
        newone.nodes = copy.copy(self.nodes)
        newone.nodes_added_log = copy.copy(self.nodes_added_log)

        # copy the bond information
        # improve
        copied_bonds = {}
        for pair, bonded_pairs_set in self.matched_pairs_bonds.items():
            copied_bonds[pair] = copy.copy(bonded_pairs_set)
        newone.matched_pairs_bonds = copied_bonds

        # copy the mirrors
        newone.mirrors = copy.copy(self.mirrors)
        newone.alternative_mappings = copy.copy(self.alternative_mappings)

        # fixme - check any other lists that you keep track of
        return newone

    def findMirrorChoices(self):
        """
        For each pair (A1, B1) find all the other options in the mirrors where (A1, B2)
        # ie Ignore (X, B1) search, if we repair from A to B, then B to A should be repaired too

        # fixme - is this still necessary if we are traversing all paths?
        """
        choices = {}
        for A1, B1 in self.matched_pairs:
            options_for_A1 = []
            for mirror in self.mirrors:
                for A2, B2 in mirror.matched_pairs:
                    if A1 is A2 and not B1 is B2:
                        options_for_A1.append(B2)

            if options_for_A1:
                options_for_A1.insert(0, B1)
                choices[A1] = options_for_A1

        return choices

    def addAlternativeMapping(self, weird_symmetry):
        """
        This means that there is another way to traverse and overlap the two molecules,
        but that the self is better (e.g. lower rmsd) than the other one
        """
        self.alternative_mappings.append(weird_symmetry)

    def correct_for_coordinates(self):
        """
        Use the coordinates of the atoms, to figure out which symmetries are the correct ones.
        Rearrange so that the overall topology represents the one that has appropriate coordinates,
        whereas all the mirrors represent the other poor matches.

        # fixme - ensure that each node is used only once at the end
        """

        # check if you have coordinates
        # fixme - rn we have it, check

        # superimpose the coordinates, ensure a good match
        # fixme - this was done before, so let's leave this way for now

        # fixme - consider putting this conf as a mirror, and then modifying this

        # check which are preferable for each of the mirrors
        # we have to match mirrors to each other, ie say we have (O1=O3) and (O2=O4)
        # we should find the mirror matching (O1=O4) and (O2=O3)
        # so note that we have a closure here: All 4 atoms are used in both cases, and each time are paired differently.
        # So this is how we defined the mirror - and therefore we can reduce this issue to the minimal mirrors.
        # fixme - is this a cycle? O1-O3-O2-O4-O1
        # Let's try to define a chain: O1 =O3, and O1 =O4, and O2 is =O3 or =O4
        # So we have to define how to find O1 matching to different parts, and then decide
        choices_mapping = self.findMirrorChoices()

        # fixme - rewrite this method to elminiate one by one the hydrognes that fit in perfectly,
        # some of them will have a plural significant match, while others might be hazy,
        # so we have to eliminate them one by one, searching the best matches and then eliminating them

        removed_nodes = set()
        for A1, choices in choices_mapping.items():
            # remove the old tuple
            # fixme - not sure if this is the right way to go,
            # but we break all the rules when applying this simplistic strategy
            self.remove_node_pair((A1, choices[0]))
            removed_nodes.add(A1)
            removed_nodes.add(choices[0])

        shortest_dsts = []

        added_nodes = set()

        # better matches
        # for each atom that mismatches, scan all molecules and find the best match and eliminate it
        blacklisted_bxs = []
        for _ in range(len(choices_mapping)):
            # fixme - optimisation of this could be such that if they two atoms are within 0.2A or something
            # then they are straight away fixed
            closest_dst = 9999999
            closest_a1 = None
            closest_bx = None
            for A1, choices in choices_mapping.items():
                # so we have several choices for A1, and now naively we are taking the one that is closest, and
                # assuming the superimposition is easy, this would work

                # FIXME - you cannot use simply distances, if for A1 and A2 the best is BX, then BX there should be
                # rules for that
                for BX in choices:
                    if BX in blacklisted_bxs:
                        continue
                    # use the distance_array because of PBC correction and speed
                    A1_BX_dst = distance_array(np.array([A1.position, ]), np.array([BX.position, ]))[0]
                    if A1_BX_dst < closest_dst:
                        closest_dst = A1_BX_dst
                        closest_bx = BX
                        closest_a1 = A1

            # across all the possible choices, found the best match now:
            blacklisted_bxs.append(closest_bx)
            shortest_dsts.append(closest_dst)
            log(closest_a1.atomName, 'is matching best with', closest_bx.atomName)

            # remove the old tuple and insert the new one
            self.add_node_pair((closest_a1, closest_bx))
            added_nodes.add(closest_a1)
            added_nodes.add(closest_bx)
            # remove from consideration
            del choices_mapping[closest_a1]
            # blacklist

        # fixme - check that the added and the removed nodes are the same set
        assert removed_nodes == added_nodes

        # this is the corrected region score (there might not be any)
        if len(shortest_dsts) != 0:
            avg_dst = np.mean(shortest_dsts)
        else:
            # fixme
            avg_dst = 0

        return avg_dst

    def are_matched_sets(self, l_atoms, r_atoms):
        if len(l_atoms) != len(r_atoms):
            return False

        for atom in l_atoms:
            _, matched_r = self.get_pair_with_atom(atom)
            if matched_r not in r_atoms:
                return False

        return True

    def enforce_no_partial_rings(self):
        """
        http://www.alchemistry.org/wiki/Constructing_a_Pathway_of_Intermediate_States
        It is the openning or closing of the rings that is an issue.
        This means that if any atom on a ring disappears, it breaks the ring,
        and therefore the entire ring should be removed and appeared again.

        If any atom is removed, it should check if it affects other rings,
        therefore cascading removing furhter rings.
        """
        MAX_CIRCLE_SIZE = 7

        # get circles in the original ligands
        L_circles, R_circles = self.getOriginalCircles()
        L_matched_circles, R_matched_circles = self.getCircles()

        # right now we are filtering out circles that are larger than 7 atoms,
        L_circles = list(filter(lambda c: len(c) <= MAX_CIRCLE_SIZE, L_circles))
        R_circles = list(filter(lambda c: len(c) <= MAX_CIRCLE_SIZE, R_circles))
        L_matched_circles = list(filter(lambda c: len(c) <= MAX_CIRCLE_SIZE, L_matched_circles))
        R_matched_circles = list(filter(lambda c: len(c) <= MAX_CIRCLE_SIZE, R_matched_circles))

        # first, see which matched circles eliminate themselves (simply matched circles)
        correct_circles = []
        for l_matched_circle in L_matched_circles[::-1]:
            for r_matched_circle in R_matched_circles[::-1]:
                if self.are_matched_sets(l_matched_circle, r_matched_circle):
                    # These two circles fully overlap, so they are fine
                    L_matched_circles.remove(l_matched_circle)
                    R_matched_circles.remove(r_matched_circle)
                    # update the original circles
                    L_circles.remove(l_matched_circle)
                    R_circles.remove(r_matched_circle)
                    correct_circles.append((l_matched_circle, r_matched_circle))

        # at this point, we should not have any matched circles, in either R and L
        # this is because we do not allow one ligand to have a matched circle, while another ligand not
        assert len(L_matched_circles) == len(R_matched_circles) == 0

        while True:
            # so now we have to work with the original rings which have not been overlapped,
            # these most likely means that there are mutations preventing it from overlapping
            l_removed_pairs = self._remove_unmatched_ring_atoms(L_circles)
            r_removed_pairs = self._remove_unmatched_ring_atoms(R_circles)

            for l_circle, r_circle in correct_circles:
                # checked if any removed atom affected any of the correct circles
                affected_l_circle = any(l_atom in l_circle for l_atom, r_atom in l_removed_pairs)
                affected_r_circle = any(r_atom in r_circle for l_atom, r_atom in r_removed_pairs)
                # add the circle to be disassembled
                if affected_l_circle or affected_r_circle:
                    L_circles.append(l_circle)
                    R_circles.append(r_circle)

            if len(l_removed_pairs) == len(r_removed_pairs) == 0:
                break

    def _remove_unmatched_ring_atoms(self, circles):
        """
        A helper function. Removes pairs with the given atoms.

        The removed atoms are classified as unmatched_rings.

        Parameters
        ----------
        circles : list
            A list of iterables. Each atom in a circle, if matched, is removed together with
            the corresponding atom from the suptop.
            The user should ensure that the rings/circles are partial

        Returns
        -------
        removed : bool
            True if any atom was removed. False otherwise.
        """
        removed_pairs = []
        for circle in circles:
            for unmatched_ring_atom in circle:
                # find if the ring has a match
                if self.contains_node(unmatched_ring_atom):
                    # remove the pair from matched
                    pair = self.get_pair_with_atom(unmatched_ring_atom)
                    self.remove_node_pair(pair)
                    self._removed_because_unmatched_rings.append(pair)
                    removed_pairs.append(pair)
        return removed_pairs

    def get_pair_with_atom(self, atom):
        for a1, a2 in self.matched_pairs:
            if a1 is atom:
                return (a1, a2)
            elif a2 is atom:
                return (a1, a2)

        return None

    def get_toppology_similarity_score(self):
        """
        Having the superimposed A(Left) and B(Right), score the match.
        This is a rather naive approach. It compares A-B match by checking
        if any of the node X and X' in A and B have a bond to another node Y that is
        not present in A-B, but that is directly reachable from X and X' in a similar way.
        We ignore the charge of Y and focus here only on the topology.

        For every "external bond" from the component we try to see if topologically it scores well.
        So for any matched pair, we extend the topology and the score is equal to the number size of
        such outputed component. Then we do this for all other matching nodes and sum the score.

        # fixme - maybe you should use the entire graphs in order to see if this is good or not?
        so the simpler approach is to ignore charges for a second to only understand the relative place in the topology,
        in other words, the question is, how similar are two nodes A and B vs A and C? let's traverse A and B together,
        and then A and C together, and while doing that, ignore the charges. In this case, A and B could
        get together 20 parts, whereas A and C traverses together 22 parts, meaning that topologically,
        it is a more suitable one, because it closer corresponds to the actual atom.
        Note that this approach has problem:
        - you can imagine A and B traversing where B is in a completely wrong global place, but it
        happens to have a bigger part common to A, than C which globally is correct. Answer to this:
        at the same time, ideally B would be excluded, because it should have been already matched to another
        topology.

        Alternative approach: take into consideration other components and the distance from this component
        to them. Specifically, allows mismatches

        FIXME - allow flexible mismatches. Meaning if someone mutates one bonded atom, then it might be noticed
        that
        """
        overall_score = 0
        for node_a, node_b in self.matched_pairs:
            # for every neighbour in Left
            for bonded_atom in node_a.bonds:
                # if this bonded atom is present in this superimposed topology (or component), ignore
                # fixme - surely this can be done better, you could have "contains this atom or something"
                in_this_sup_top = False
                for other_a, _ in self.matched_pairs:
                    if bonded_atom == other_a:
                        in_this_sup_top = True
                        break
                if in_this_sup_top:
                    continue

                # a candidate is found that could make the node_a and node_b more similar,
                # so check if it is also present in node_b,
                # ignore the charges to focus only on the topology and put aside the parameterisation
                for bonded_atom_b in node_b.bonds:
                    # fixme - what if the atom is mutated into a different atom? we have to be able
                    # to relies on other measures than just this one, here the situation is that the topology
                    # is enough to answer the question (because only charges were modified),
                    # however, this gets more tricky
                    # fixme - hardcoded
                    score = len(_overlay(bonded_atom, bonded_atom_b, atol=99999))

                    # this is a purely topology based score, the bigger the overlap the better the match
                    overall_score += score

                # check if the neighbour points to any node X that is not used in Left,

                # if node_b leads to the same node X
        return overall_score

    def refineAgainstCharges(self, atol):
        """
        Removes the matched pairs which turn out to have charges more different
        than the given absolute tolerance (atol) [Electron units]

        After removing a pair, removes any hydrogens attached to the heavy atoms.
        """
        # walk through the superimposed topologies
        # and move the atom-atom pairs that suffer from being the same
        removed_pairs = []
        removed_attached_hydrogens = []
        for node1, node2 in self.matched_pairs[::-1]:
            if node1.eq(node2, atol=atol):
                continue

            # remove any dangling hydrogens from this pair
            removed_h_pairs = self.remove_attached_hydrogens((node1, node2))
            removed_attached_hydrogens.extend(removed_h_pairs)
            # remove this pair
            self.remove_node_pair((node1, node2))

            removed_pairs.append((node1, node2))
        # keep track of the removed atoms due to the charge
        # fixme
        # if not self._removed_due_to_charge is None:
        #     raise Exception('the charges have already been refined, should not be called twice')
        self._removed_due_to_charge = removed_pairs

        return removed_pairs, removed_attached_hydrogens

    def is_consistent_cycles(self, suptop):
        # check if each sup top has the same number of cycles
        # fixme - not sure?
        selfG1, selfG2 = self.getNxGraphs()
        self_cycles1, self_cycles2 = len(nx.cycle_basis(selfG1)), len(nx.cycle_basis(selfG2))
        if self_cycles1 != self_cycles2:
            raise Exception('left G has a different number of cycles than right G')

        otherG1, otherG2 = suptop.getNxGraphs()
        other_cycles1, other_cycles2 = len(nx.cycle_basis(otherG1)), len(nx.cycle_basis(otherG2))
        if other_cycles1 != other_cycles2:
            raise Exception('left G has a different number of cycles than right G')

        # check if merging the two is going to create issues
        # with the circle inequality
        # fixme - Optimise: reuse the above nx graphs rather than making an entire copy
        self_copy = copy.copy(self)
        # we ignore the parent here because it is only to check the number of circles
        self_copy.merge(suptop)
        if not self_copy.sameCircleNumber():
            return False

        return True

    def is_consistent_with(self, suptop):
        """
        Conditions:
            - There should be a minimal overlap of at least 1 node.
            - There is no no pair (A=B) in this sup top such that (A=C) or (B=C) exists in other.
            - The number of cycles in this suptop and the other suptop must be the same
            - merging cannot lead to new cycles?? (fixme). What is the reasoning behind this?
                I mean, I guess the assumption is that, if the cycles were comptabile,
                they would be created during the search, rather than now while merging. ??
        """

        # confirm that there is no mismatches, ie (A=B) in suptop1 and (A=C) in suptop2 where (C!=B)
        for node1, node2 in self.matched_pairs:
            for nodeA, nodeB in suptop.matched_pairs:
                if (node1 is nodeA) and not (node2 is nodeB):
                    return False
                elif (node2 is nodeB) and not (node1 is nodeA):
                    return False

        # ensure there is at least one common pair
        if self.count_common_node_pairs(suptop) == 0:
            return False

        if not self.is_consistent_cycles(suptop):
            return False

        return True

    @staticmethod
    def _rename_ligand(atoms, name_counter=None):
        """
        name_counter: a dictionary with atom as the key such as 'N', 'C', etc,
        the counter keeps track of the last used counter for each name.
        Empty means that the counting will start from 1.
        """
        if name_counter is None:
            name_counter = {}

        for atom in atoms:
            # get the first letters that is not a character
            afterLetters = [i for i, l in enumerate(atom.atomName) if l.isalpha()][-1] + 1

            atom_name = atom.atomName[:afterLetters]
            last_used_counter = name_counter.get(atom_name, 0)

            # rename
            last_used_counter += 1
            newAtomName = atom_name + str(last_used_counter)
            print(f'Renaming {atom.atomName} to {newAtomName}')
            atom.atomName = newAtomName

            # update the counter
            name_counter[atom_name] = last_used_counter

        return name_counter

    @staticmethod
    def _get_atom_names_counter(atoms):
        """
        name_counter: a dictionary with atom as the key such as 'N', 'C', etc,
        the counter keeps track of the last used counter for each name.
        Ie if there are C1, C2, C3, this will return {'C':3} as the last counter.
        """
        name_counter = {}

        for atom in atoms:
            # get the first letters that is not a character
            afterLetters = [i for i, l in enumerate(atom.atomName) if l.isalpha()][-1] + 1

            atom_name = atom.atomName[:afterLetters]
            atom_number = int(atom.atomName[afterLetters:])
            last_used_counter = name_counter.get(atom_name, 0)

            # update the counter
            name_counter[atom_name] = max(last_used_counter, atom_number)

        return name_counter

    @staticmethod
    def _is_correct_atom_name_format(atoms):
        # check if the atom format is C15, ie atom name followed by a number
        for atom in atoms:
            afterLetters = [i for i, l in enumerate(atom.atomName) if l.isalpha()][-1] + 1

            atom_name = atom.atomName[:afterLetters]
            if len(atom_name) == 0:
                return False

            atom_number = atom.atomName[afterLetters:]
            try:
                int(atom_number)
            except:
                return False

        return True

    @staticmethod
    def rename_ligands(L_nodes, R_nodes):
        # rename the ligand to ensure that no atom has the same name
        # name atoms using the first letter (C, N, ..) and count them
        # keep the names if possible (ie if they are already different)

        # first, ensure that all the atom names are unique
        L_atom_names = [a.atomName for a in L_nodes]
        L_names_unique = len(set(L_atom_names)) == len(L_atom_names)
        L_correct_format = SuperimposedTopology._is_correct_atom_name_format(L_nodes)

        if not L_names_unique or not L_correct_format:
            print('Renaming Left Molecule Atom Names (Because it is needed)')
            name_counter_L_nodes = SuperimposedTopology._rename_ligand(L_nodes)
            L_atom_names = [a.atomName for a in L_nodes]
        else:
            name_counter_L_nodes = SuperimposedTopology._get_atom_names_counter(L_nodes)

        R_atom_names = [a.atomName for a in R_nodes]
        R_names_unique = len(set(R_atom_names)) == len(R_atom_names)
        R_correct_format = SuperimposedTopology._is_correct_atom_name_format(R_nodes)
        L_R_overlap = len(set(R_atom_names).intersection(set(L_atom_names))) > 0

        if not R_names_unique or not R_correct_format or L_R_overlap:
            print('Renaming Right Molecule Atom Names (Because it is needed)')
            SuperimposedTopology._rename_ligand(R_nodes, name_counter=name_counter_L_nodes)
        # each atom name is unique, fixme - this check does not apply anymore
        # ie it is fine for a molecule to use general type
        #assert len(set(R_atom_names)) == len(R_atom_names)

        return

    def getNxGraphs(self):
        "maybe at some point this should be created and used internally more? "
        gl = nx.Graph()
        gr = nx.Graph()
        # add each node
        for nA, nB in self.matched_pairs:
            gl.add_node(nA)
            gr.add_node(nB)
        # add all the edges
        for nA, nB in self.matched_pairs:
            # add the edges from nA
            for bonded_to_nA, bondtype1 in nA.bonds:
                if bonded_to_nA in gl:
                    gl.add_edge(nA, bonded_to_nA)
            for bonded_to_nB, bondtype1 in nB.bonds:
                if bonded_to_nB in gr:
                    gr.add_edge(nB, bonded_to_nB)

        return gl, gr

    def getCircles(self):
        """
        REturn circles found in the matched pairs.
        """
        gl, gr = self.getNxGraphs()
        gl_circles = [set(circle) for circle in nx.cycle_basis(gl)]
        gr_circles = [set(circle) for circle in nx.cycle_basis(gr)]
        return gl_circles, gr_circles

    def getOriginalCircles(self):
        """
        Return the original circles present in the input topologies.
        """
        # create a circles
        L_original = self._getOriginalCircle(self.top1)
        R_original = self._getOriginalCircle(self.top2)

        L_circles = [set(circle) for circle in nx.cycle_basis(L_original)]
        R_circles = [set(circle) for circle in nx.cycle_basis(R_original)]
        return L_circles, R_circles

    def _getOriginalCircle(self, atom_list):
        """Create a networkx circle out of the list
        atom_list - list of AtomNode
        """
        G = nx.Graph()
        # add each node
        for atom in atom_list:
            G.add_node(atom)

        # add all the edges
        for atom in atom_list:
            # add the edges from nA
            for other in atom_list:
                if atom.isBoundTo(other):
                    G.add_edge(atom, other)

        return G

    def getCircleNumber(self):
        gl_circles, gr_circles = self.getCircles()
        return len(gl_circles), len(gr_circles)


    def sameCircleNumber(self):
        gl_num, gr_num = self.getCircleNumber()
        if gl_num == gr_num:
            return True

        return False

    def cycle_spans_multiple_cycles(self):
        # This filter checks whether a newly created suptop cycle spans multiple cycles
        # this is one of the filters (#106)
        # fixme - should this be applied whenever we work with more than 1 cycle?
        # it checks whether any cycles in the left molecule,
        # is paired with more than one cycle in the right molecule
        """
        What is the circle is shared?
        We are using cycles which excluded atoms that join different rings.
        fixme - could this lead to a special case?
        """

        for lcycle in self._nonoverlapping_lcycles:
            overlap_counter = 0
            for rcycle in self._nonoverlapping_rcycles:
                # check if the cycles overlap
                if self._cycles_overlap(lcycle, rcycle):
                    overlap_counter += 1

            if overlap_counter > 1:
                return True

        for rcycle in self._nonoverlapping_rcycles:
            overlap_counter = 0
            for lcycle in self._nonoverlapping_lcycles:
                # check if the cycles overlap
                if self._cycles_overlap(lcycle, rcycle):
                    overlap_counter += 1

            if overlap_counter > 1:
                return True

        return False

    def _cycles_overlap(self, lcycle, rcycle):
        # check if any nodes are paired across the two cycles
        # any to any pairing
        for l, r in itertools.product(lcycle, rcycle):
            if self.contains((l, r)):
                return True

        return False

    def merge(self, suptop):
        """
        Absorb the other suptop by adding all the node pairs that are not present
        in the current sup top.

        WARNING: ensure that the other suptop is consistent with this sup top.
        """
        # assert self.is_consistent_with(suptop)

        # print("About the merge two sup tops")
        # self.print_summary()
        # other_suptop.print_summary()

        merged_pairs = []
        for pair in suptop.matched_pairs:
            # check if this pair is present
            if not self.contains(pair):
                n1, n2 = pair
                if self.contains_node(n1) or self.contains_node(n2):
                    raise Exception('already uses that node')
                # pass the bonded pairs here
                self.add_node_pair(pair)
                merged_pairs.append(pair)
        # after adding all the nodes, now add the bonds
        for pair in merged_pairs:
            # add the connections
            bonded_pairs = suptop.matched_pairs_bonds[pair]
            assert len(bonded_pairs) > 0
            self.link_pairs(pair, bonded_pairs)

        # removed from the "merged" the ones that agree, so it contains only the new stuff
        # to make it easier to read
        self.nodes_added_log.append(("merged with", merged_pairs))

        # check if duplication occured, fixme - temporary
        return merged_pairs

    @staticmethod
    def Validate_Charges(atom_listL, atom_listR):
        """
        Check the original charges:
        - ensure that the total charge of L and R are integers
        - ensure that they are equal to the same integer
        """
        whole_left_charge = sum(a.charge for a in atom_listL)
        np.testing.assert_almost_equal(whole_left_charge, round(whole_left_charge), decimal=2,
                                       err_msg=f'left charges are not integral. Expected {round(whole_left_charge)}'
                                               f' but found {whole_left_charge}')

        whole_right_charge = sum(a.charge for a in atom_listR)
        np.testing.assert_almost_equal(whole_right_charge, round(whole_right_charge), decimal=2,
                                       err_msg=f'right charges are not integral. Expected {round(whole_right_charge)}'
                                               f' but found {whole_right_charge}'
                                       )
        # same integer
        np.testing.assert_almost_equal(whole_left_charge, whole_right_charge, decimal=2)

        return round(whole_left_charge)

    def redistribute_charges(self):
        """
        After the match is made and the user commits to the superimposed topology,
        the charges can be revised.
        We calculate the average charges between every match, and check how that affects
        the rest of the molecule (the unmatched atoms).
        Then, we distribute the charges to the unmatched atoms to get
        the net charge as a whole number/integer.

        This function should be called after removing the matches for whatever reason.
        ie at the end of anything that could modify the atom pairing.
        """

        SuperimposedTopology.Validate_Charges(self.top1, self.top2)

        # the total charge in the matched region before the changes
        matched_total_chargeL = sum(l.charge for l, r in self.matched_pairs)
        matched_total_chargeR = sum(r.charge for l, r in self.matched_pairs)

        # average charges between matched atoms
        l_delta_charge_total = 0
        r_delta_charge_total = 0
        for l, r in self.matched_pairs:
            if l.charge != r.charge:
                avg_charge = (l.charge + r.charge) / 2.0
                l_delta_charge_total += l.charge - avg_charge
                r_delta_charge_total += r.charge - avg_charge
                # this new charge is made to each molecule
                l.charge = r.charge = avg_charge

        # fixme should matched_total_chargeL be l_delta_charge_total?

        # get the unmatched nodes in L and R
        L_unmatched = self.get_disappearing_atoms()
        R_unmatched = self.get_appearing_atoms()


        if len(L_unmatched) == 0 and l_delta_charge_total != 0:
            print('----------------------------------------------------------------------------------------------')
            print('ERROR? AFTER AVERAGING CHARGES, THERE ARE NO UNMATCHED ATOMS TO ASSIGN THE CHARGE TO: left ligand.')
            print('----------------------------------------------------------------------------------------------')
        if len(R_unmatched) == 0 and r_delta_charge_total != 0:
            print('----------------------------------------------------------------------------------------------')
            print('ERROR? AFTER AVERAGING CHARGES, THERE ARE NO UNMATCHED ATOMS TO ASSIGN THE CHARGE TO: right ligand. ')
            print('----------------------------------------------------------------------------------------------')

        # distribute the charges over the unmatched regions
        if len(L_unmatched) != 0:
            l_delta_per_atom = float(l_delta_charge_total) / len(L_unmatched)
        else:
            l_delta_per_atom = 0

        if len(R_unmatched) != 0:
            r_delta_per_atom = float(r_delta_charge_total) / len(R_unmatched)
        else:
            r_delta_per_atom = 0

        # redistribute that delta q over the atoms in the left and right molecule
        for atom in L_unmatched:
            atom.charge += l_delta_per_atom
        for atom in R_unmatched:
            atom.charge += r_delta_per_atom

        # note that we are really modifying right now the original nodes.
        SuperimposedTopology.Validate_Charges(self.top1, self.top2)

    def contains_node(self, node):
        # checks if this node was used in this overlay
        if len(self.nodes.intersection(set([node,]))) == 1:
            return True

        return False

    def contains_any_node(self, node_list):
        if len(self.nodes.intersection(set(node_list))) > 0:
            return True

        return False

    def count_common_nodes(self, node_list):
        number_of_common_nodes = len(self.nodes.intersection(set(node_list)))
        return number_of_common_nodes

    def count_common_node_pairs(self, other_suptop):
        return len(set(self.matched_pairs).intersection(set(other_suptop.matched_pairs)))

    def contains_any_node_from(self, other_sup_top):
        if len(self.nodes.intersection(other_sup_top.nodes)) > 0:
            return True

        return False

    def contains(self, node_pair):
        for match_pair in self.matched_pairs:
            if match_pair == node_pair:
                return True

        return False

    def contains_atomNamePair(self, atomName1, atomName2):
        for m1, m2 in self.matched_pairs:
            if m1.atomName == atomName1 and m2.atomName == atomName2:
                return True

        return False

    def contains_left_atomName(self, atomName):
        for m1, _ in self.matched_pairs:
            if m1.atomName == atomName:
                return True

        return False

    def contains_right_atomName(self, atomName):
        for _ , m in self.matched_pairs:
            if m.atomName == atomName:
                return True

        return False

    def contains_all(self, other_sup_top):
        for pair in other_sup_top.matched_pairs:
            if not self.contains(pair):
                return False

        return True

    def contains_same_atoms_symmetric(self, other_sup_top):
        """
        The atoms can be paired differently, but they are the same.
        """
        if len(self.nodes.symmetric_difference(other_sup_top.nodes)) == 0:
            return True

        return False

    def has_in_contrast_to(self, sup_top):
        return set(self.matched_pairs).difference(set(sup_top.matched_pairs))

    def report_differences(self, suptop):
        selfHasNotSuptop = self.has_in_contrast_to(suptop)
        log("self has not suptop:", selfHasNotSuptop)
        suptopHasNotSelf = suptop.has_in_contrast_to(self)
        log('Suptop has not self', suptopHasNotSelf)
        return selfHasNotSuptop, suptopHasNotSelf

    def has_left_nodes_same_as(self, other):
        if len(self.matched_pairs) != len(other.matched_pairs):
            return False

        for node1, _ in self.matched_pairs:
            # check if each node exists in the other
            node_found = False
            for other_node, _ in other.matched_pairs:
                if node1 == other_node:
                    node_found = True

            if not node_found:
                return False

        return True

    def has_right_nodes_same_as(self, other):
        if len(self.matched_pairs) != len(other.matched_pairs):
            return False

        for _, right_node in self.matched_pairs:
            # check if each node exists in the other
            node_found = False
            for _, other_right in other.matched_pairs:
                if right_node == other_right:
                    node_found = True

            if not node_found:
                return False

        return True

    def is_subgraph_of(self, other_sup_top):
        """
        Checks if this superimposed topology is a subgraph of another superimposed topology.
        Or if any mirror topology is a subgraph.
        """
        # subgraph cannot be equivalent self.eq, it is only proper subgraph (ie proper subset)
        if len(self.matched_pairs) >= len(other_sup_top.matched_pairs):
            return False

        # self is smaller, so it might be a subgraph
        if other_sup_top.contains_all(self):
            return True

        # self is not a subgraph, but it could be a subgraph of one of the mirrors
        for mirror in self.mirrors:
            if other_sup_top.contains_all(mirror):
                return True

        # other is bigger than self, but not a subgraph of self
        return False

    def subgraph_relationship(self, other_sup_top):
        """
        Return
        1 if self is a supergraph of other,
        -1 if self is a subgraph of other
        0 if they have the same number of elements (regardless of what the nodes are)
        """
        if len(self.matched_pairs) == len(other_sup_top.matched_pairs):
            return 0

        if len(self.matched_pairs) > len(other_sup_top.matched_pairs):
            # self is bigger than other,
            # check if self contains all nodes in other
            if self.contains_all(other_sup_top):
                return 1

            # other is not a subgraph, but check the mirrors if any of them are
            for mirror in self.mirrors:
                if mirror.contains_all(other_sup_top):
                    return 1

            # other is smaller but not a subgraph of this graph or any of its mirrors
            return 0

        if len(self.matched_pairs) < len(other_sup_top.matched_pairs):
            # other is bigger, so self might be a subgraph
            # check if other contains all nodes in self
            if other_sup_top.contains_all(self):
                return -1

            # self is not a subgraph, but it could be a subgraph of one of the mirrors
            for mirror in self.mirrors:
                if other_sup_top.contains_all(mirror):
                    return -1

            # other is bigger than self, but it is not a subgraph
            return 0

    def is_mirror_of(self, other_sup_top):
        """
        this is a naive check
        fixme - check if the found superimposed topology is the same (ie the same matches), what then?

        some of the superimposed topologies represent symmetrical matches,
        for example, imagine T1A and T1B is a symmetrical version of T2A and T2B,
        this means that
         - the number of nodes in T1A, T1B, T2A, and T2B is the same
         - all the nodes in T1A are in T2A,
         - all the nodes in T1B are in T2B
        """

        if len(self.matched_pairs) != len(other_sup_top.matched_pairs):
            return False

        if self.contains_same_atoms_symmetric(other_sup_top):
            return True

        return False

    def add_mirror_suptop(self, suptop):
        assert len(self.matched_pairs) == len(suptop.matched_pairs)
        # check if this this mirror was already added
        for mirror in self.mirrors:
            if suptop.eq(mirror):
                # a mirror like that already exists
                return

        # when you "absorb" another suptop as a mirror, extract its mirrors too
        self.mirrors.extend(suptop.mirrors)
        suptop.mirrors = []

        # add the mirror
        self.mirrors.append(suptop)

    def eq(self, sup_top):
        """
        Check if the superimposed topology is "the same". This means that every pair has a corresponding pair in the
        other topology (but possibly in a different order)
        """
        # fixme - should replace this with networkx
        if len(self) != len(sup_top):
            return False

        for pair in self.matched_pairs:
            # find for every pair the matching pair
            if not sup_top.contains(pair):
                return False

        return True


# todo - move to logging rather than this
verbose_log = False
def log(*args):
    if verbose_log:
        print(*args)


def get_largest(lists):
    """
    return a list of largest solutions
    """
    solution_sizes = [len(st) for st in lists]
    largest_sol_size = max(solution_sizes)
    return list(filter(lambda st: len(st) == largest_sol_size, lists))


def long_merge(suptop1, suptop2):
    """
    Carry out a merge and apply all checks.
    Merge suptop2 into suptop1.

    """
    if suptop1 is suptop2:
        return suptop1

    if suptop1.eq(suptop2):
        log("Merge: the two are the equal. Ignoring")
        return suptop1

    if suptop2.is_subgraph_of(suptop1):
        log("Merge: this is already a superset. Ignoring")
        return suptop1

    # check if the two are consistent
    # ie there is no clashes
    if not suptop1.is_consistent_with(suptop2):
        log("Merge: cannot merge - not consistent")
        return -1

    # fixme - this can be removed because it is now taken care of in the other functions?
    # g1, g2 = suptop1.getNxGraphs()
    # assert len(nx.cycle_basis(g1)) == len(nx.cycle_basis(g2))
    # g3, g4 = suptop2.getNxGraphs()
    # assert len(nx.cycle_basis(g3)) == len(nx.cycle_basis(g4))
    #
    # print("Will merge", suptop1, 'and', suptop2)
    # assert suptop1.sameCircleNumber()
    newly_added_pairs = suptop1.merge(suptop2)

    # if not suptop1.sameCircleNumber():
    #     raise Exception('something off')
    # # remove sol2 from the solutions:
    # all_solutions.remove(sol2)
    return newly_added_pairs


def solve_one_combination(one_atom_spieces, ignore_coords):
    atoms = one_atom_spieces
    if len(atoms) == 1:
        # simple case,  one in the left ligands, many in the right ligand, pick the best
        atom, candidates = list(atoms.items())[0]
        largest_candidates = get_largest(list(candidates.values()))
        best = extractBestSuptop(largest_candidates, ignore_coords)
        return best
    elif len(atoms) > 1:
        # many in the left ligand, unknown in the right ligand
        # case: only 1 in the right ligand, so Many2One case
        # check if all the atoms in the left ligand map to the same atom in the right ligand
        unique_atoms = set()
        for keys in atoms.values():
            for name in keys.keys():
                unique_atoms.add(name)
        if len(unique_atoms) == 1:
            log('Many (left) to one (right)')
            # just pick the best match fro the right ligand
            candidates = [list(v.values())[0] for v in atoms.values()]
            largest_candidates = get_largest(candidates)
            best = extractBestSuptop(largest_candidates, ignore_coords)
            return best

        # ie one type
        # many to many, ie
        # L1 L2 and R1 R2, so we have to create to create mappings,
        # and then merge, and then evaluate
        left_ligand_keys = list(atoms.keys())
        right_ligand_keys = list(unique_atoms)

        # generate all pairs first
        # fixme - this should be moved up?
        # fixme - these functions should be factored out and tested separately
        all_pairs = list(itertools.product(left_ligand_keys, right_ligand_keys))
        # if the left ligand can only match 2 atoms but the right 3 atoms,
        # then the best possible match is that of 2 pairs
        longest_match = min(len(left_ligand_keys), len(right_ligand_keys))
        # generate all possible combinations
        all_combinations = list(itertools.combinations(all_pairs, longest_match))
        # filter out the ones that are impossible
        chosen = list(filter(lambda s: len({pair[0] for pair in s}) == longest_match
                                       and len({pair[1] for pair in s}) == longest_match,
                             all_combinations))

        # fixme - use itertools instead?

        def combine(all_pairs):
            combined = set()
            # used pairs keep track of which pairs have been combined
            used_pairs = set()
            for i, pair in enumerate(all_pairs):
                # take one item, and try to combine it with as many as possible
                basis = [pair, ]
                # for j, pair2 in enumerate(all_pairs):
                # if the pair does not have any common elements, combine
                # if not any(a1 == a or a2 == b for a,b in basis):
                #     basis.append((a1, a2))
                #     used_pairs.add(basis[0])
                #     used_pairs.add((a1, a2))
                # if the basis was not combined with anything
                if len(basis) == 1 and basis[0] in used_pairs:
                    continue
                combined.add(tuple(basis))
            # remove the subsets of larger sets
            return combined

        generated_combinations = chosen
        # now that we have the pairs combined, we have to attempt to marge them this way,
        alternatives = []
        for combined_pairs in generated_combinations:
            merged = None
            for lk, rk in combined_pairs:
                # note that if the key is not present, this means it is not an option,
                # for example, it might have been tried, but because of the circle consistency rule,
                # it might have been deleted, for that reason that generated theoretical combination
                # will not work
                if lk not in atoms:
                    continue
                lk_map = atoms[lk]
                if rk not in lk_map:
                    continue

                top = lk_map[rk]
                if merged is None:
                    merged = copy.copy(top)
                    continue

                # check the output,
                # fixme - what to do when this is wrong?
                long_merge(merged, top)
            if merged is None:
                # this probably means that this combination is not possible e.g. due to cycles
                continue
            alternatives.append(merged)
        assert len(alternatives) >= 1
        # now that all alternatives have been computed,
        # decide which is best
        largest_candidates = get_largest(alternatives)
        return extractBestSuptop(largest_candidates, ignore_coords)

    raise Exception('not implemented')


def _overlay(n1, n2, parent_n1, parent_n2, bond_types, suptop, ignore_coords=False, useGenType=True):
    """
    Jointly and recurvisely traverse the molecule while building up the suptop.

    If n1 and n2 are the same, we will be traversing through both graphs, marking the jointly travelled areas.

    Return the topology of the common substructure between the two molecules.

    *n1 from the left molecule,
    *n2 from the right molecule
    """

    # if either of the nodes has already been matched, ignore
    if suptop.contains_any_node([n1, n2]):
        return None

    # if the two nodes are "the same"
    if useGenType and not n1.sameGenType(n2):
        # these two atoms have a different type, so return None
        return None
    elif not useGenType and not n1.sameExactType(n2):
        # these two atoms have a different type, so return None
        return None
   
    # Check for cycles
    # if a new cycle is created by adding this node,
    # then the cycle should be present in both, left and right ligand
    safe = True
    # if n1 is linked with node in suptop other than parent
    a_new_cycle_is_present = False
    for b1 in n1.bonds:
        if b1[0] != parent_n1 and suptop.contains_node(b1[0]):
            safe = False # n1 forms cycle, now need to check n2
            a_new_cycle_is_present = True
            for b2 in n2.bonds:
                if b2[0] != parent_n2 and suptop.contains_node(b2[0]):
                    # b2 forms cycle, now need to check it's the same in both
                    if suptop.contains((b1[0],b2[0])):
                        safe = True
                        break
            if not safe: # only n1 forms a cycle
                break
    if not safe: # either only n1 forms cycle or both do but different cycles
        return None
    
    # now the same for any remaining unchecked bonds in n2
    safe = True
    for b2 in n2.bonds:
        if b2[0] != parent_n2 and suptop.contains_node(b2[0]):
            safe = False
            a_new_cycle_is_present = True
            for b1 in n1.bonds:
                if b1[0] != parent_n1 and suptop.contains_node(b1[0]):
                    if suptop.contains((b1[0], b2[0])):
                        safe = True
                        break
            if not safe:
                break
    if not safe:
        return None

    # check if the cycle spans multiple cycles present in the left and right molecule,
    if suptop.cycle_spans_multiple_cycles():
        log('Found a cycle spanning multiple cycles')
        return None

    log("Adding ", (n1, n2), "in", suptop.matched_pairs)
    # append both nodes as a pair to ensure that we keep track of the mapping
    # having both nodes appended also ensure that we do not revisit/readd neither n1 and n2
    suptop.add_node_pair((n1, n2))
    if not (parent_n1 is parent_n2 is None):
        suptop.link_with_parent((n1, n2), (parent_n1, parent_n2), bond_types)

    # the extra bonds are legitimate
    # so let's make sure they are added
    for n1bonded_node, btype1 in n1.bonds:
        # ignore left parent
        if n1bonded_node is parent_n1:
            continue
        for n2bonded_node, btype2 in n2.bonds:
            # ignore right parent
            if n2bonded_node is parent_n2:
                continue

            # if the pair exists, add a bond between the two pairs
            if suptop.contains((n1bonded_node, n2bonded_node)):
                suptop.link_pairs((n1, n2),
                    [((n1bonded_node, n2bonded_node), (btype1, btype2)), ])

    # filter out parents
    n1_bonds_no_parent = list(filter(lambda bond: not bond[0] is parent_n1, n1.bonds))
    n2_bonds_no_parent = list(filter(lambda bond: not bond[0] is parent_n2, n2.bonds))

    # sort for the itertools.groupby
    n1_bonds_no_parent_srt = sorted(n1_bonds_no_parent, key=lambda b: b[0].gentype)
    n2_bonds_no_parent_srt = sorted(n2_bonds_no_parent, key=lambda b: b[0].gentype)

    # so first we have to group with groupby
    combinations = []
    for n1_type, n1_bonds in itertools.groupby(n1_bonds_no_parent_srt,
                                               key=lambda bonded: bonded[0].gentype):
        for n2_type, n2_bonds in itertools.groupby(n2_bonds_no_parent_srt,
                                                   key=lambda bonded: bonded[0].gentype):
            # the types will only match once
            # fixme - note that it always uses a general type here,
            # also, the group is not done correctly atm
            if n1_type != n2_type:
                continue

            # convert into a list
            n1_bonds = list(n1_bonds)
            n2_bonds = list(n2_bonds)

            # these two groups are of the same type, so we have to do each to each combination,

            # if one atom list of length 1, then we have L1-R1 and L1-R2 and
            # these two are in contradiction to each other

            # For 2 C in left and right, we have 2*2=4 combinations,
            # ie LC1-RC1 LC2-RC2 and LC1-RC2 LC2-RC1 which shuld be split into two groups

            # For 3 C ... fixme
            atom_type_solutions = {}
            for n1bonded_node, btype1 in n1_bonds:
                # so we first try each combintion with atom1
                solutions_for_this_left_atom = {}
                # so for every atom we have a list of choices,
                for n2bonded_node, btype2 in n2_bonds:
                    # a copy of the sup_top is needed because the traversal can take place
                    # using different pathways
                    bond_solutions = _overlay(n1bonded_node, n2bonded_node,
                                              parent_n1=n1, parent_n2=n2,
                                              bond_types=(btype1, btype2),
                                              suptop=copy.copy(suptop),
                                              ignore_coords=ignore_coords,
                                              useGenType=useGenType)
                    # if they were different type, etc, ignore
                    if bond_solutions is None:
                        continue

                    assert type(bond_solutions) is not list
                    solutions_for_this_left_atom[n2bonded_node] = bond_solutions

                # record all possible solution for this one atom
                if len(solutions_for_this_left_atom) > 0:
                    atom_type_solutions[n1bonded_node] = solutions_for_this_left_atom
            if len(atom_type_solutions) != 0:
                combinations.append(atom_type_solutions)

    # fixme
    # so we have all combinations for each atom,
    # and we have to put them together sensibly,

    # if there is no solution, return the current suptop
    if len(combinations) == 0:
        return suptop
    # simplest case: 1 atom type
    elif len(combinations) == 1:
        return solve_one_combination(combinations[0], ignore_coords)

    all_solutions = []
    # there is multiple atom types,
    # and each which have its own solution, which in turn have to be merged
    for atom_type in combinations:
        all_solutions.append(solve_one_combination(atom_type, ignore_coords))

    assert len(all_solutions) > 1
    log('combons done')

    # for n1bonded_node, btype1 in n1.bonds:
    #     for n2bonded_node, btype2 in n2.bonds:
    #         # if either of the nodes has already been matched, ignore
    #         if n1bonded_node is parent_n1 or n2bonded_node is parent_n2:
    #             continue
    #
    #         # a copy of the sup_top is needed because the traversal can take place
    #         # using different pathways
    #         bond_solutions = _overlay(n1bonded_node, n2bonded_node,
    #                                   parent_n1=n1, parent_n2=n2,
    #                                   bond_types=(btype1, btype2),
    #                                   suptop=copy.copy(suptop))
    #         # if they were different type, etc, ignore
    #         if bond_solutions is None:
    #             continue
    #
    #         all_solutions.extend(bond_solutions)
            # fixme - when you have a mirror like ester O1-O2 matches, then you could store them differently here

    # fixme - there should never the "equal" suptops in the solutions?
    # in relation to #15
    # combine the different walks,
    # Take two of -C-CH2-C- and focus on the middle part.
    # some solutions will be shorter: for example
    # the two H can be matched in two ways.
    # so we have to understand which solutions are "mirror"
    # images of themselves.
    # ie only one out of the two ways to match the hydrogens is corrects
    # here, we find which walks are alternatives to each other,
    # and we pick the best one
    # we could try to merge each with each which would create
    # lots of different combinations, but ideally
    # we would be able to say which solutions are in conflict with which
    # other solutions,
    # Say that there are 5 solutions and 4 of them are in conflict,

    # try every possible pathway

    # sort in the descending order
    all_solutions.sort(key=lambda st: len(st), reverse=True)
    for sol1 in all_solutions:
        for sol2 in all_solutions[::-1]:
            if sol1 is sol2:
                continue

            if sol1.eq(sol2):
                log("Found the same solution and removing, solution", sol1.matched_pairs)
                all_solutions.remove(sol2)
                continue

            if sol2.is_subgraph_of(sol1):
                continue

            # TODO: Remove?
            if sol1.is_consistent_with(sol2):
                # print("merging, current pair", (n1, n2))
                # join sol2 and sol1 because they're consistent
                log("Will merge", sol1, 'and', sol2)
                newly_added_pairs = sol1.merge(sol2)

                # remove sol2 from the solutions:
                all_solutions.remove(sol2)

    # after all the merges, return the best matches only,
    # the other mergers are wrong (and they are smaller)
    solution_sizes = [len(s2) for s2 in all_solutions]
    largest_sol_size = max(solution_sizes)
    # if there is more than one solution at this stage, reduce it by checking the rmsd
    # fixme - add dihedral angles etc and make sure you always return with just one best solution
    # maybe we can note the other solutions as part of this solution to see understand the differences
    # ie merge the other solutions here with this suptop, if it is a mirror then add it to the suptop,
    # for later analysis, rather than returning which would have to be compared with a lot of other parts
    bestSuptop = extractBestSuptop(list(filter(lambda st: len(st) == largest_sol_size, all_solutions)), ignore_coords)
    return bestSuptop


class Topology:
    """
    A helper class to organise a topology and its associated functions
    """
    def __init__(self, nodes):
        self.nodes = nodes

        # generate a networkx graph
        graph = nx.Graph()
        [graph.add_node(node) for node in nodes]
        for node in nodes:
            for bonded in node.bonds:
                graph.add_edge(node, bonded)
        self.nxgraph = graph

        self.cycles()

    def cycles(self):
        # find the cycles
        cycles = [frozenset(c) for c in nx.cycle_basis(self.nxgraph)]

        self.joined_cycles = {}
        for cycle in cycles:
            self.joined_cycles[cycle] = set()

        # see if any cycles overlap with 2 atoms
        # this would mean that they are in the same 2D plane
        for i, cycle1 in enumerate(cycles):
            for cycle2 in cycles[i + 1:]:
                if len(set(cycle1).intersection(set(cycle2))) >= 2:
                    # the two cycles overlap and are in the same plane
                    self.joined_cycles[cycle1].add(cycle2)
                    self.joined_cycles[cycle2].add(cycle1)


    def inSamePlane(self, cycle1, cycle2):
        if cycle1 in self.joined_cycles[cycle1]:
            return True

        return False

def superimpose_topologies(top1_nodes, top2_nodes, pair_charge_atol=0.1, use_charges=True,
                           use_coords=True, starting_node_pairs=None,
                           force_mismatch=None, no_disjoint_components=True,
                           net_charge_filter=True, net_charge_threshold=0.1,
                           redistribute_charges_over_unmatched=True,
                           ligandLmda=None, ligandRmda=None,
                           align_molecules=True,
                           partial_rings_allowed=True,
                           ignore_charges_completely=False,
                           ignore_bond_types=False,
                           ignore_coords=False,
                           left_coords_are_ref=True,
                           use_general_type=True,
                           use_only_gentype=False,
                           check_atom_names_unique=True):
    """
    This is a helper function that managed the entire process.

    TODO:
    - check if each molecule topology is connected
    - run the superimpose while ignoring the charges
    - run the superimpose with charges
    - check if any charges components are subcomponent of a larger charge-ignoring component,
    this would be useful with solving some dillemas, assign them to each other

    Other to think about:
    - what would happen if you have mutation that separates the molecule? what happens when you multiple of them?
    how do you match them together?

    """

    if not ignore_charges_completely:
        whole_charge = SuperimposedTopology.Validate_Charges(top1_nodes, top2_nodes)

    # ensure that none of the atom names across the two molecules are the different
    if check_atom_names_unique:
        sameAtomNames = {a.atomName for a in top1_nodes}.intersection({a.atomName for a in top2_nodes})
        assert len(sameAtomNames) == 0, f"The molecules have the same atom names. This is now allowed. They are: {sameAtomNames}"

    # prealign the 3D coordinates before applaying further changes
    # todo
    # if align_molecules:
    #     take_largest = lambda x, y: x if len(x) > len(y) else y
    #     reduce(take_largest, suptops).alignLigandsUsingMatched()

    suptops = _superimpose_topologies(top1_nodes, top2_nodes, ligandLmda, ligandRmda, starting_node_pairs=starting_node_pairs,
                                      ignore_coords=ignore_coords, left_coords_are_ref=left_coords_are_ref,
                                      use_general_type=use_general_type)

    if not suptops:
        raise Exception('Did not find a single superimposition state. ')

    for st in suptops:
        print('Original suptop len %d' % len(st))
        # st.print_summary()

    # ignore bond types
    for st in suptops:
        # fixme - do proper
        st.ignore_bond_types = ignore_bond_types

    # connect the suptops to their original molecules
    for suptop in suptops:
        # fixme this should have been done in the constructor?
        suptop.set_tops(top1_nodes, top2_nodes)
        suptop.set_mdas(ligandLmda, ligandRmda)

    # align the 3D coordinates before applaying further changes
    # use the largest suptop to align the molecules
    if align_molecules:
        take_largest = lambda x,y: x if len(x) > len(y) else y
        reduce(take_largest, suptops).alignLigandsUsingMatched()

    # fixme - you might not need because we are now doing this on the way back
    # if useCoords:
    #     for sup_top in sup_tops:
    #         sup_top.correct_for_coordinates()

    # ensure the actual atom types is correct (general atom type can be used to match atoms)
    # fixme - this is going to be another stage
    if not use_only_gentype:
        for st in suptops:
            st.matched_atom_types_are_the_same()

    # note that charges need to be checked before assigning IDs.
    # ie if charges are different, the matched pair
    # becomes two different atoms with different IDs
    if use_charges and not ignore_charges_completely:
        for sup_top in suptops:
            sup_top.refineAgainstCharges(atol=pair_charge_atol)
            if sup_top._removed_due_to_charge:
                print(f'Removed due to charge incompatibility: {sup_top._removed_due_to_charge}')

    # apply the force mismatch at the end
    # this is an interactive feature
    if force_mismatch is not None:
        for suptop in suptops:
            for an1, an2 in force_mismatch:
                if suptop.contains_atomNamePair(an1, an2):
                    n1, n2 = suptop.get_node(an1), suptop.get_node(an2)
                    suptop.remove_node_pair((n1, n2))

    if net_charge_filter and not ignore_charges_completely:
        # ensure that each found component has net charge < 0.1
        for suptop in suptops[::-1]:
            # fixme this should be function within suptop
            while np.abs(suptop.get_net_charge()) > net_charge_threshold:
                largest_difference = suptop.remove_worst_charge_match()
                if largest_difference == 0:
                    raise Exception("How can there be net charge but no pair can be found that actually is different?")

                # check if there are any atoms left
                if len(suptop) == 0:
                    # remove the suptop from the list
                    suptops.remove(suptop)
                    break
            print (f'Removed pairs due to net charge: {suptop._removed_due_to_net_charge}')

    if not partial_rings_allowed:
        # remove partial rings, note this is a cascade problem if there are double rings
        for suptop in suptops:
            suptop.enforce_no_partial_rings()
            print(f'Removed pairs because partial rings are not allowed {suptop._removed_because_unmatched_rings}')

    # remove the suptops that are empty
    for st in suptops[::-1]:
        if len(st) == 0:
            suptops.remove(st)

    if no_disjoint_components:
        print(f'Checking for disjoint components in the {len(suptops)} suptops')
        # ensure that each suptop represents one CC
        # check if the graph was divided after removing any pairs (e.g. due to charge mismatch)
        [st.only_largest_CC_survives() for st in suptops]

        for st in suptops:
            print('Removed disjoint components: ', st._removed_because_disjointed_cc)

        # remove the smaller suptop, or one arbitrary if they are equivalent
        if len(suptops) > 1:
            max_len = max([len(suptop) for suptop in suptops])
            for suptop in suptops[::-1]:
                if len(suptop) < max_len:
                    suptops.remove(suptop)

            # if there are equal length suptops left, take only the first one
            if len(suptops) > 1:
                suptops = [suptops[0]]

        assert len(suptops) == 1, suptops

    #
    if redistribute_charges_over_unmatched and not ignore_charges_completely:
        if len(suptops) > 1:
            raise NotImplementedError(
                'Currently distributing charges works only if there is no disjointed components')
        suptops[0].redistribute_charges()

    # atom ID assignment has to come after any removal of atoms due to their mismatching charges
    start_atom_id = 1
    for suptop in suptops:
        start_atom_id = suptop.assign_atoms_ids(start_atom_id)
        # increase the start ID by 1 so that the next sup top assigns them
        start_atom_id += 1

    # there might be several best solutions, order them according the RMSDs
    suptops.sort(key=lambda st: st.rmsd())

    # fixme - remove the hydrogens without attached heavy atoms

    # resolve_sup_top_multiple_match(sup_tops_charges)
    # sup_top_correct_chirality(sup_tops_charges, sup_tops_no_charges, atol=atol)

    # carry out a check. Each
    if align_molecules:
        for st in suptops:
            main_rmsd = st.alignLigandsUsingMatched()
            for mirror in st.mirrors:
                mirror_rmsd = mirror.alignLigandsUsingMatched()
                if mirror_rmsd < main_rmsd:
                    print('THE MIRROR RMSD IS LOWER THAN THE MAIN RMSD')
            st.alignLigandsUsingMatched(overwrite_original=True)

    return suptops


def calculateRmsd(atomPairs):
    deviations = []
    for atom1, atom2 in atomPairs:
        # check how far the atoms are to each other
        deviations.append(atom1.position - atom2.position)
    return np.sqrt(np.mean(((np.array(deviations)) ** 2)))


def extractBestSuptop(suptops, ignore_coords):
    # fixme - ignore coords currently does not work
    if ignore_coords:
        raise NotImplementedError('Ignoring coords during superimposition is currently not possible')
    # multiple different paths to traverse the topologies were found
    # this means some kind of symmetry in the topologies
    # For example, in the below drawn case (starting from C1-C11) there are two
    # solutions: (O1-O11, O2-O12) and (O1-O12, O2-O11).
    #     LIGAND 1        LIGAND 2
    #        C1              C11
    #        \                \
    #        N1              N11
    #        /\              / \
    #     O1    O2        O11   O12
    # Here we decide which of the mappings is better.
    # fixme - uses coordinates to decide which mapping is better.
    #  - Improve: use dihedral angles to decide which mapping is better too
    if len(suptops) == 0:
        raise Exception('Cannot generate the best mapping without any suptops')

    if len(suptops) == 1:
        # there is only one solution
        return suptops[0]

    candidates = copy.copy(suptops)

    # align the subcomponent and check
    rmsds = []
    for suptop in candidates:
        rmsd = suptop.alignLigandsUsingMatched()
        rmsds.append(rmsd)

    # get the best rmsd
    best = rmsds.index(min(rmsds))
    candidate_superimposed_top = candidates[best]

    for suptop in suptops:
        if suptop is candidate_superimposed_top:
            continue

        if suptop.is_mirror_of(candidate_superimposed_top):
            candidate_superimposed_top.add_mirror_suptop(suptop)
            continue

        candidate_superimposed_top.addAlternativeMapping(suptop)

    return candidate_superimposed_top


def getBestRmsdMatch(suptops):
    # multiple different paths to traverse the topologies were found
    # this means some kind of symmetry in the topologies
    # For example, in the below drawn case (starting from C1-C11) there are two
    # solutions: (O1-O11, O2-O12) and (O1-O12, O2-O11).
    #     LIGAND 1        LIGAND 2
    #        C1              C11
    #        \                \
    #        N1              N11
    #        /\              / \
    #     O1    O2        O11   O12
    # Here we decide which of the mappings is better.
    # fixme - uses coordinates to decide which mapping is better.
    #  - Improve: use dihedral angles to decide which mapping is better too
    if len(suptops) == 0:
        raise Exception('Cannot generate the best mapping without any suptops')

    if len(suptops) == 1:
        # there is only one solution
        return suptops[0]

    best_suptop = None
    best_rmsd = np.finfo('float32').max
    rmsd = None
    for suptop in suptops:
        # use the avg dst after the correction to understand which one is better,
        # and assign the worse
        # fixme - avg dst would ideally not be np.NaN
        avg_dst = suptop.correct_for_coordinates()
        rmsd = suptop.rmsd()
        # get a global match (RMSD)
        if rmsd < best_rmsd:
            best_suptop = suptop
            best_rmsd = rmsd

    candidate_superimposed_top = best_suptop

    for suptop in suptops:
        if suptop is candidate_superimposed_top:
            continue

        candidate_superimposed_top.addWeirdSymmetry(suptop)

    return candidate_superimposed_top


def alreadySeen(suptop, suptops):
    for next_suptop in suptops:
        if next_suptop.eq(suptop):
            return True

    return False


def isMirrorOfOne(suptop, suptops, ignore_coords):
    """
    "Mirror" in the sense that it is an alternative topological way to traverse the molecule.

    Depending on the "better" fit between the two mirrors, we pick the one that is better.
    """
    for next_suptop in suptops:
        if next_suptop.is_mirror_of(suptop):
            # the suptop saved as the mirror should be the suptop
            # that is judged to be of a lower quality
            bestSuptop = extractBestSuptop([suptop, next_suptop], ignore_coords)
            if bestSuptop is next_suptop:
                next_suptop.add_mirror_suptop(suptop)

            # the new suptop is better than the previous one
            suptops.remove(next_suptop)
            bestSuptop.add_mirror_suptop(next_suptop)
            suptops.append(bestSuptop)

            return True

    return False


def solve_partial_overlaps(candidate_suptop, suptops):
    # check if this candidate suptop uses a node that is used by a larger sup top
    # fixme: optimisation: this whole processing should happen probably at the end?
    for suptop in suptops[::-1]:  # reverse traversal in case deleting is necessary
        # there is an overlap, some of the ndoes are reused
        if suptop.contains_any_node_from(candidate_suptop):
            if len(suptop.matched_pairs) > len(candidate_suptop.matched_pairs):
                # the overlapping existing suptop is larger
                # ignore this sup top
                return True
            elif len(suptop.matched_pairs) < len(candidate_suptop.matched_pairs):
                # candidate sup top is larger, so that is the one we want to keep
                # delete the smaller sup top
                suptops.remove(suptop)
            else:
                # the two sup tops could be of the same length and yet traverse different atoms
                # e.g. case: mcl1_l17l9, where due to the symmetry, two different traversals
                # of the same length are found, which use different set of atoms

                # note that we already checked if it is a mirror - this is not a mirror
                # fixme: these two would ideally be chained together closer
                # to avoid any future bugs

                # there is a partial overlap, so two different ways to score these
                # you could call them "symmetries". Here we have to pick
                # which is the "worse symmetry",
                # let us use atom coordinates to score them
                if suptop.rmsd() < candidate_suptop.rmsd():
                    suptop.addAlternativeMapping(candidate_suptop)
                    ignore_cand_sup_top = True
                else:
                    candidate_suptop.addAlternativeMapping(suptop)
                    suptops.remove(suptop)


def isSubGraphOf(candidate_suptop, suptops):
    # check if the newly found subgraph is a subgraph of any other sup top
    # fixme - is this even possible?
    # fixme can any subgraph be a subgraph of another?
    for suptop in suptops:
        if candidate_suptop.is_subgraph_of(suptop):
            return True

    return False


def removeCandidatesSubgraphs(candidate_suptop, suptops):
    # check if this new sup top should replace other sup_tops, because they are its subgraphs
    # fixme - I am not sure if this can happen twice
    removed_subgraphs = False
    for suptop in suptops[::-1]:
        if suptop.is_subgraph_of(candidate_suptop):
            log('Removing candidate\'s subgraphs')
            suptops.remove(suptop)
            removed_subgraphs = True
    return removed_subgraphs


def _superimpose_topologies(top1_nodes, top2_nodes, mda1_nodes=None, mda2_nodes=None,
                            starting_node_pairs=None,
                            ignore_coords=False,
                            left_coords_are_ref=True,
                            use_general_type=True):
    """
    Superimpose two molecules.
    """
    # generate the graph for the top1 and top2
    top1 = Topology(top1_nodes)
    top2 = Topology(top2_nodes)

    # superimposed topologies
    suptops = []
    # grow the topologies using every combination node1-node2 as the starting point
    # fixme - Test/Optimisation: create a theoretical maximum of a match between two molecules
    # - Proposal 1: find junctions and use them to start the search
    # - Analyse components of the graph (ie rotatable due to a single bond connection) and
    #   pick a starting point from each component
    if not starting_node_pairs:
        # generate each to each nodes
        starting_node_pairs = itertools.product(top1_nodes, top2_nodes)

    for node1, node2 in starting_node_pairs:
        # fixme - optimisation - reduce the number of starting nX and nY pairs

        # if node1.get_id() != 10 or node2.get_id() != 25:
        #     continue

        # with the given starting two nodes, generate the maximum common component
        suptop = SuperimposedTopology(top1_nodes, top2_nodes, mda1_nodes, mda2_nodes)
        # fixme turn into a property
        suptop.left_coords_are_ref = left_coords_are_ref
        candidate_suptop = _overlay(node1, node2, parent_n1=None, parent_n2=None, bond_types=(None, None),
                                    suptop=suptop,
                                    ignore_coords=ignore_coords,
                                    useGenType=use_general_type)
        if candidate_suptop is None or len(candidate_suptop) == 0:
            # there is no overlap, ignore this case
            continue

        # check if the maximal possible solution was found
        # Optimise - can you at this point finish the superimposition if the molecules are fully superimposed?
        candidate_suptop.is_subgraph_of_global_top()

        # ignore if this topology was found before
        if alreadySeen(candidate_suptop, suptops):
            continue

        # ignore if it is a subgraph of another solution
        if isSubGraphOf(candidate_suptop, suptops):
            continue

        # check if this superimposed topology is a mirror of one that already exists
        # fixme the order matters in this place
        # fixme - what if the mirror has a lower rmsd match? in that case, pick that mirror here
        if isMirrorOfOne(candidate_suptop, suptops, ignore_coords):
            continue

        removed_subgraphs = removeCandidatesSubgraphs(candidate_suptop, suptops)
        if removed_subgraphs:
            suptops.append(candidate_suptop)
            continue

        # while comparing partial overlaps, suptops can be modified
        andIgnore = solve_partial_overlaps(candidate_suptop, suptops)
        if andIgnore:
            continue

        # fixme - what to do when about the odd pairs randomH-randomH etc? they won't be found in other graphs
        # follow a rule: if this node was used before in a larger superimposed topology, than it should
        # not be in the final list (we guarantee that each node is used only once)
        suptops.append(candidate_suptop)

    # if there are only hydrogens superimposed without a connection to any heavy atoms, ignore these too
    for suptop in suptops[::-1]:
        all_hydrogens = True
        for node1, _ in suptop.matched_pairs:
            if not node1.type == 'H':
                all_hydrogens = False
                break
        if all_hydrogens:
            log("Removing sup top because only hydrogens found", suptop.matched_pairs)
            suptops.remove(suptop)

    # TEST: check that each node was used only once
    all_nodes = []
    pair_count = 0
    for suptop in suptops:
        [all_nodes.extend([node1, node2]) for node1, node2 in suptop.matched_pairs]
        pair_count += len(suptop.matched_pairs)
    # fixme
    # assert len(set(all_nodes)) == 2 * pair_count

    # TEST: check that the nodes on the left are always from topology 1 and the nodes on the right are always from top2
    for suptop in suptops:
        for node1, node2 in suptop.matched_pairs:
            assert node1 in list(top1_nodes)
            assert node2 in list(top2_nodes)

    # clean the overlays by removing sub_overlays.
    # ie if all atoms in an overlay are found to be a bigger part of another overlay,
    # then that overlay is better
    log("Found altogether overlays", len(suptops))

    # fixme - return other info
    return suptops


def resolve_sup_top_multiple_match(sup_tops):
    # initially, see if you can resolve the problem of multiple match by using the larger superimposed element
    # TODO - check if some components can be mapped to multiple other components
    # find first the repeating component that can be matched to other components,
    same_left_sup_tops = []
    same_right_sup_tops = []
    for i, sup_top1 in enumerate(sup_tops):
        # fixme - add special messages when 3-1 or some other combination is found!
        # fixme - is 2-2 possible?
        for sup_top2 in sup_tops[i + 1:]:
            if sup_top1 is sup_top2:
                continue

            if sup_top1.has_left_nodes_same_as(sup_top2):
                # if [A,B] and [B,C] then combine to [A,B,C]
                # check if either A or B was was found before
                added_to_previous = False
                for same_left_sup_top in same_left_sup_tops:
                    if same_left_sup_top == sup_top1:
                        same_left_sup_top.append(sup_top2)
                        added_to_previous = True
                        break
                    elif same_left_sup_top == sup_top2:
                        same_left_sup_top.append(sup_top1)
                        added_to_previous = True
                        break
                if not added_to_previous:
                    same_left_sup_tops.append([sup_top1, sup_top2])
                log('found same left', sup_top1.matched_pairs, 'with', sup_top1.matched_pairs)

            if sup_top1.has_right_nodes_same_as(sup_top2):
                added_to_previous = False
                for same_right_sup_top in same_right_sup_tops:
                    if same_right_sup_top == sup_top1:
                        same_right_sup_top.append(sup_top2)
                        added_to_previous = True
                        break
                    elif same_right_sup_top == sup_top2:
                        same_right_sup_top.append(sup_top1)
                        added_to_previous = True
                        break
                if not added_to_previous:
                    same_right_sup_tops.append([sup_top1, sup_top2])
                log('found same right', sup_top1.matched_pairs)

    # first, attempt to see if you can resolve the conflict by going back to the sup_top without charges,
    for same_left_sup_top_list in same_left_sup_tops:
        multiple_match_that_have_superset = []
        multiple_match_no_superset = []
        for same_left_sup_top in same_left_sup_top_list:
            # check if this sup_top is correct according to the global sup_top without charges
            if same_left_sup_top.has_uncharged_superset_sup_top():
                multiple_match_that_have_superset.append(same_left_sup_top)
            else:
                multiple_match_no_superset.append(same_left_sup_top)

        # remove the sup tops that have no super set # fixme - is this correct?
        for sup_top in multiple_match_no_superset:
            sup_tops.remove(sup_top)

        assert len(multiple_match_that_have_superset) == 1

        # mark this list as resolved by emptying it
        [same_left_sup_top_list.remove(l) for l in same_left_sup_top_list[::-1]]

    # every one should be solved
    assert all([len(l) == 0 for l in same_left_sup_tops])

    assert len(same_right_sup_tops) == 0, 'not implemented yet'

    return

    # if [A,B] and [B,C] then combine to [A,B,C]
    # now that you extracted sup tops that are the same, check what is the mapping,
    # easiest example is 2-to-1, so 2 on the left map to the same one to the right,
    # to do this, we kind of need to construct a graph
    # fixme - think of other mappings 2-to-2 etc - create an error for now
    # check if we are working with n-to-1
    for same_left_sup_top_list in same_left_sup_tops:
        # so we know that the left top is the same, so we need to figure out which of the right top is the right match
        # this means basically ranking them on this list
        # if, we example, 0-A-B-C is with 0-A'-B'-C where A==A' and A==B', then A==A' because they both connect to 0,
        # in other words, check for the highest number of common connections between A and A' and B'

        multiple_match_that_have_superset = []
        for same_left_sup_top in same_left_sup_top_list:
            # for each of the bonded atoms in Left, check if the Right matched atom has also the same bonded atoms
            # for each such atom add one point
            # FIXME - this could be moved from here, and since this score requires left and right,
            # it should be computed inside of the SuperimposedTopology class
            score = same_left_sup_top.get_toppology_similarity_score()
            # keep track of the score for this match
            multiple_match_that_have_superset.append(score)

        # make sure that they all are not equal to 0
        assert all([0 != score for score in multiple_match_that_have_superset])
        # make sure that the top scoring find is not duplicated, ie that we have a clear winner
        assert multiple_match_that_have_superset.count(max(multiple_match_that_have_superset)) == 1

        # choose the best scoring match based on the previous work
        winner_index = multiple_match_that_have_superset.index(max(multiple_match_that_have_superset))
        log("multiple match winner is", same_left_sup_top_list[winner_index].matched_pairs)

        # remove the losers
        for index, worse_match in enumerate(same_left_sup_top_list):
            if index == winner_index:
                pass
            else:
                log("Removing a worse match", worse_match.matched_pairs)
                sup_tops.remove(worse_match)

        # remove the deleted not chosen topologies but keep track of them and return them as well,
    # fixme - what about the right side?


def sup_top_correct_chirality(sup_tops, sup_tops_no_charge, atol):
    # fixme chirality
    # some can elements can be chiral, e.g. imagine that you have three components A-B-C and A-B'-C
    # where B' is reversed B in the other direction, this means that all components together with B will
    # be found as separate components, meaning that nothing will mutate,
    # but actually if you look globally at it, it is clear that B' and B have different ends, and they should
    # be mutated. For that reason, it is important to check the connection from B and B' and see what that means
    # so first we want to test if the component has neighbour components, meaning if you can take a step
    # from B that would directly take you to C (but that would not fully solve it)

    # there is a more minimal case for chirality/assymetry, A-B and A-B'. This is not equivalent to a mutation of a node
    # separates A and B. However, there might be cases where it is very close. The mutation would separate A and B
    # such that even if you reversed B, it would not match (in most cases). However, in the case of chirality,
    # the reversal can match better. Say we have to sequences X=a-b-c-d-e and Y=a-b-e-d-c, such that A=a-b and
    # B=c-d-e (and therefore =e-d-c). You might notcie that d is always in the same place, and therefore should
    # be considered to be the same atom, whereas e and c swapped their places, and therefore should be
    # appearing/disappearing.

    # first, we need to detect chirality/assymetry. The condition in the last example is that
    # that B is connected via b in both cases - by the same atom (.eq which considers charges), even though
    # there is a superimposed component, which is supposed to maximise its space.
    # For that reason, B should be flagged as a "symmetric component" which should not be.


    # in the case of MCL1 you would think that we can check against the super set sup_top without charges.
    # however, this could be more tricky, because the superset sup_top without charges happens not to be a superset
    # but it provides some information: basically, there is no structure like that,
    # sup top without charges is our template now, because we know we only check the atom types, and therefore get
    # the match in a better way (based on the atop type)
    # fixme? finish this paragraph
    for sup_top in sup_tops:
        # check if any of the nodes are present in any of the discharged
        for sup_top_no_charge in sup_tops_no_charge:
            # check if the sup_top has misassigned nl-nr according to the sup top without charges
            # ie identify the sup top without charges that overlaps with this one

            pass

    # For each component, check if the matching topologies are reversed. If B=B' but they have a node x,
    # which is connected to y that is not in B, and x connects to y in B but not in B', then we know we have
    # the reverse relationship
    for sup_top in sup_tops:
        for node1, _ in sup_top.matched_pairs:
            for _, node2 in sup_top.matched_pairs:
                # fixme - i think you want the matched ones

                # check if they are topologically linked to the same molecule,
                # because it is impossible for two different nodes to be linked to exactly the same atom (same .eq)
                # because in that case that atom would belong to this component (ie it is the same,
                # and it is reachable).
                for bond1 in list(node1.bonds):
                    # ignore the bonds that are part of the component,
                    # could build "external bonds" method in sup top
                    if sup_top.contains_node(set([bond1, ])):
                        continue
                    for bond2 in list(node2.bonds):
                        if sup_top.contains_node(set([bond2, ])):
                            continue

                        if bond1.eq(bond2, atol=atol):
                            # there might be at any time any two nodes that are similar enough (eq), which means
                            # this is not a universal approach in itself, however, do we gain anything more
                            # knowing that one of the nodes is a part of another component? fixme
                            log("found assymetry", node1.atomName, node2.atomName,
                                  "due to", bond1.atomName, bond2.atomName)
            pass
        pass

    # add a test against the overall match (global match that ignores the charges)
    # FIXME - if two superimposed components come from two different places in the global map, then something's off
    # particularly, it could help with chirality - if with respect to the global match,
    # a local sup top travels in the wrong direction, then we have a clear issue


def get_atoms_bonds_from_ac(ac_file):
    # returns
    # 1) a dictionary with charges, e.g. Item: "C17" : -0.222903
    # 2) a list of bonds

    ac_lines = open(ac_file).readlines()

    # fixme - hide hydrogens
    # ac_lines = filter(lambda l:not('h' in l or 'H' in l), ac_lines)

    # extract the atoms
    # ATOM      1  C17 MOL     1      -5.179  -2.213   0.426 -0.222903        ca
    atom_lines = filter(lambda l:l.startswith('ATOM'), ac_lines)

    atoms = []
    for line in atom_lines:
        ATOM_phrase, atomID, atomName, resName, resID, x, y, z, charge, atom_colloq = line.split()
        x, y, z = float(x), float(y), float(z)
        charge = float(charge)
        resID = int(resID)
        atomID = int(atomID)
        atom = AtomNode(name=atomName, type=atom_colloq)
        atom.set_charge(charge)
        atom.set_id(atomID)
        atom.set_position(x, y, z)
        atom.set_resname(resName)
        atoms.append(atom)

    # fixme - add a check that all the charges come to 0 as declared in the header

    # extract the bonds, e.g.
    #     bondID atomFrom atomTo ????
    # BOND    1    1    2    7    C17  C18
    bond_lines = filter(lambda l: l.startswith('BOND'), ac_lines)
    bonds = [(int(bondFrom), int(bondTo)) for _, bondID, bondFrom, bondTo, something, atomNameFrom, atomNameTo in
             [l.split() for l in bond_lines]]

    return atoms, bonds


def load_mol2_wrapper(filename):
    # Load a .mol2 file
    # fixme - fuse with the .pdb loading, no distinction needed
    # ignore the .mol2 warnings about the mass
    warnings.filterwarnings(action='ignore', category=UserWarning,
                            message='Failed to guess the mass for the following atom types: '  # warning to ignore
                            )
    # squash the internal warning about parsing .mol2 within MDAnalysis
    warnings.filterwarnings(action='ignore', category=UserWarning,
                            message='Creating an ndarray from ragged nested sequences '  # warning to ignore
                            )
    u = mda.Universe(filename)
    return u


def get_atoms_bonds_from_mol2(ref_filename, mob_filename, use_general_type=True):
    """
    Use MDAnalysis to load the .mol2 files.

    Use MDAnalysis to superimpose the second structure onto the first structure.
    Examples: https://www.mdanalysis.org/MDAnalysisTutorial/analysismodule.html
    """
    # returns
    # 1) a dictionary with charges, e.g. Item: "C17" : -0.222903
    # 2) a list of bonds

    uref = load_mol2_wrapper(ref_filename)
    umob = load_mol2_wrapper(mob_filename)

    # this RMSD superimposition requires the same number of atoms to be superimposed
    # find out the RMSD between them and the rotation matrix
    # uref0 = uref.atoms.positions - uref.atoms.center_of_geometry()
    # umob0 = umob.atoms.positions - umob.atoms.center_of_geometry()
    #
    # # get the rotation matrix and rmsd
    # # fixme - make use of rmsd
    # R, rmsd = MDAnalysis.analysis.align.rotation_matrix(umob0, uref0)
    #
    # # update the umob atoms, the new coordinates is what we want to rely on
    # umob.atoms.translate(-umob.atoms.center_of_geometry())
    # umob.atoms.rotate(R)
    # umob.atoms.translate(uref.atoms.center_of_geometry())

    # create the atoms for left ligands
    def createAtoms(mda_atoms):
        atoms = []
        for mda_atom in mda_atoms:
            atom = AtomNode(name=mda_atom.name, type=mda_atom.type, use_general_type=use_general_type)
            try:
                atom.set_charge(mda_atom.charge)
                atom.set_original_charge(mda_atom.charge)
            except AttributeError as Att:
                # fixme - expand on the message
                print('Missing charge attribute, setting to N/A')
                atom.set_charge('N/A')
            atom.set_id(mda_atom.id)
            atom.set_position(mda_atom.position[0], mda_atom.position[1], mda_atom.position[2])
            atom.set_resname(mda_atom.resname)
            atoms.append(atom)
        return atoms

    uref_atoms = createAtoms(uref.atoms)
    # note that these coordinate should be superimposed
    umob_atoms = createAtoms(umob.atoms)

    # fixme - add a check that all the charges come to 0 as declared in the header
    uref_bonds = [(bond[0].id, bond[1].id, bond.order) for bond in uref.bonds]
    umob_bonds = [(bond[0].id, bond[1].id, bond.order) for bond in umob.bonds]

    return uref_atoms, uref_bonds, umob_atoms, umob_bonds, uref, umob



def assign_coords_from_pdb(atoms, pdb_atoms):
    """
    Match the atoms from the MDAnalysis object based on a .pdb file
    and copy the coordinates from the MDAnalysis atoms to the
    corresponding atoms.
    """
    for atom in atoms:
        # find the corresponding atom
        found_match = False
        for pdb_atom in pdb_atoms.atoms:
            if pdb_atom.name.upper() == atom.atomName.upper():
                # assign the charges
                pos = pdb_atom.position
                atom.set_position(pos[0], pos[1], pos[2])
                found_match = True
                break
        if not found_match:
            log("Did not find atom?", atom.atomName)
            raise Exception("wait a minute")
