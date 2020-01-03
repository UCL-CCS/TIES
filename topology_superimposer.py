import hashlib
import numpy as np
import networkx

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
 
 
 fixme
 - ERROR: imagine that you have a olecule with 3 overlaid Strongly Connected Components that are all the same to
 each other. Currently, their matching will be arbitrary, which could potentially lead to problems
 
 - there is a lot of matches which are substandard - possibly. 
 For example, imagine some hydrogen wrongly matched with another hydrogen, 
 These have been used in larger strongly connected components, meaning
 that this pair should be removed, which is the approach I take,  but 
 what if there is no much for one of the hydrogens nowhere? I guess it is new 
"""

class AtomNode:
    def __init__(self, atomId, atomName, resName, resId, charge, atom_colloq):
        self.atomId = atomId
        self.atomName = atomName
        self.resName = resName
        self.resId = resId
        self.charge = charge
        self.atom_colloq = atom_colloq
        self.bonds = set()

    def __hash__(self):
        m = hashlib.md5()
        m.update(str(self.atomId).encode('utf-8'))
        m.update(str(self.charge).encode('utf-8'))
        return int(m.hexdigest(), 16)

    def __str__(self):
        return self.atom_colloq

    def __repr__(self):
        return "%s_%s" % (self.atom_colloq, self.atomName)

    def bindTo(self, other):
        self.bonds.add(other)
        other.bonds.add(self)

    def eq(self, other, atol=0):
        """
        What does it mean that two atoms are the same? They are the same type and charge.
        5 % tolerance by default
        """
        if self.atom_colloq == other.atom_colloq and \
                np.isclose(self.charge, other.charge, atol=atol):
            return True

        return False


class SuperimposedTopology:
    """
    SuperimposedTopology contains in the minimal case two sets of nodes S1 and S2, which
    are paired together and represent a strongly connected component.

    However, it can also represent the symmetrical versions that were superimposed.
    """

    def __init__(self, matched_pairs, topology1, topology2):
        """
        @superimposed_nodes : a set of pairs of nodes that matched together
        """

        # TEST: with the list of matching nodes, check if each node was used only once,
        # the number of unique nodes should be equivalent to 2*len(common_pairs)
        all_matched_nodes = []
        [all_matched_nodes.extend(list(pair)) for pair in matched_pairs]
        assert len(matched_pairs) * 2 == len(all_matched_nodes)

        # todo convert to nx

        self.matched_pairs = matched_pairs
        self.top1 = topology1
        self.top2 = topology2
        self.mirrors = []
        self.nodes = set(all_matched_nodes)


    def is_subgraph_of_global_top(self):
        """
        Check if after superimposition, one graph is a subgraph of another
        :return:
        """
        # check if one topology is a subgraph of another topology
        if len(self.matched_pairs) == len(self.top1) or len(self.matched_pairs) == len(self.top2):
            print("This graph is a equivalent to the topology")
            return True

        return False


    def contains_node(self, node):
        # checks if this node was used in this overlay
        if len(self.nodes.intersection(node)) == 1:
            return True

        return False


    def contains_any_node_from(self, other_sup_top):
        if len(self.nodes.intersection(other_sup_top.nodes)) > 0:
            return True

        return False


    def contains(self, node_pair):
        for match_pair in self.matched_pairs:
            if match_pair == node_pair:
                return True

        return False


    def contains_all(self, other_sup_top):
        for pair in other_sup_top.matched_pairs:
            if not self.contains(pair):
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

        # this should not be applied on the same topology match
        if self.eq(other_sup_top):
            raise Exception("They are already the same, cannot be a mirror")

        if len(self.matched_pairs) != len(other_sup_top.matched_pairs):
            return False

        for nodeA, nodeB in self.matched_pairs:
            # find each node in the other sup top
            nodeA_found = False
            nodeB_found = False
            for node1, node2 in other_sup_top.matched_pairs:
                if nodeA == node1:
                    nodeA_found = True
                if nodeB == node2:
                    nodeB_found = True

            if not (nodeA_found and nodeB_found):
                return False

        return True


    def add_mirror_sup_top(self, mirror_sup_top):
        assert len(self.matched_pairs) == len(mirror_sup_top.matched_pairs)
        # print("mirror added")
        self.mirrors.append(mirror_sup_top)


    def eq(self, other_sup_top):
        """
        Check if the superimposed topology is "the same". This means that every pair has a corresponding pair in the
        other topology (but possibly in a different order)
        """
        # fixme - should replace this with networkx
        if len(self.matched_pairs) != len(other_sup_top.matched_pairs):
            return False

        for pair1 in self.matched_pairs:
            # find for every pair the matching pair
            exists_in_other = False
            for pair2 in other_sup_top.matched_pairs:
                if set(pair1) == set(pair2):
                    exists_in_other = True

            if not exists_in_other:
                return False

        return True


def _overlay(n1, n2, matched_nodes=None, atol=0):
    """
    n1 should be from one graph, and n2 should be from another.

    If n1 and n2 are the same, we will be traversing through both graphs, marking the jointly travelled areas.
    RETURN: the maximum common overlap with the n1 and n2 as the starting conditions.

    Here, recursively the graphs of n1 and of n2 will be explored as they are the same.
    """

    if matched_nodes == None:
        matched_nodes = []

    # if either of the nodes has already been matched, ignore this potential match
    for pair in matched_nodes:
        if n1 in pair or n2 in pair:
            return matched_nodes

    # if the two nodes are "the same", append them to the list of matched nodes
    if n1.eq(n2, atol=atol):
        # append both nodes as a pair to ensure that we keep track of the mapping
        # having both nodes appended also ensure that we do not revisit/readd neither n1 and n2
        matched_nodes.append([n1, n2])
        # continue traversing
        for n1bonded_node in n1.bonds:
            for n2bonded_node in n2.bonds:
                _overlay(n1bonded_node, n2bonded_node, matched_nodes, atol=atol)

    return matched_nodes


def superimpose_topologies(top1, top2, atol):
    """
    atol 1 means the charge of the atom can be up to 1 electron different,
    as long as the atom has the same type
    """
    sup_tops = []
    # fixme - TEST: create a theoretical maximum that could actually match if you ignored the topology completely
    # grow the topologies using every combination node1-node2 as the starting point
    for node1 in top1:
        for node2 in top2:
            # grow the topologies to see if they overlap
            matched_pairs = _overlay(node1, node2, atol=atol)
            # ignore if no atoms were superimposed
            if len(matched_pairs) == 0:
                continue

            # construct topology class
            candidate_superimposed_top = SuperimposedTopology(matched_pairs, top1, top2)
            candidate_superimposed_top.is_subgraph_of_global_top()

            # check if this superimposed topology was found before, if so, ignore
            sup_top_seen = False
            for other_sup_top in sup_tops:
                if other_sup_top.eq(candidate_superimposed_top):
                    sup_top_seen = True
                    break
            if sup_top_seen:
                continue

            # check if this new sup top should replace other sup_tops, because they are its subgraphs
            # fixme - I am not sure if this can happen twice
            for sup_top in sup_tops[::-1]:
                if sup_top.is_subgraph_of(candidate_superimposed_top):
                    print('removing previous smaller tops')
                    sup_tops.remove(sup_top)

            # check if this candidate sup top uses a node that is used by a larger sup top
            ignore_cand_sup_top = False
            for sup_top in sup_tops[::-1]:  # reverse traversal in case deleting is necessary
                # if there is a sup top with more elements
                if len(sup_top.matched_pairs) > len(candidate_superimposed_top.matched_pairs):
                    # and that sup top contains some nodes from the candidate sup top
                    # ignore this sup top
                    if sup_top.contains_any_node_from(candidate_superimposed_top):
                        ignore_cand_sup_top = True
                        break
                elif len(sup_top.matched_pairs) < len(candidate_superimposed_top.matched_pairs):
                    # this sup_top has fewer elements than candidate sup top
                    if sup_top.contains_any_node_from(candidate_superimposed_top):
                        # and there is an overlap, delete the smaller sup top
                        sup_tops.remove(sup_top)
            if ignore_cand_sup_top:
                continue


            # check if the newly found subgraph is a subgraph of any other sup top
            # fixme - is this even possible?
            # fixme can any subgraph be a subgraph of another?
            cand_is_subgraph = False
            for sup_top in sup_tops:
                if candidate_superimposed_top.is_subgraph_of(sup_top):
                    cand_is_subgraph = True
                    break
            if cand_is_subgraph:
                raise Exception('can this hapen?')
                continue

            # check if this superimposed topology is a mirror of one that already exists
            # fixme the order matters in this place
            is_mirror = False
            for other_sup_top in sup_tops:
                if other_sup_top.is_mirror_of(candidate_superimposed_top):
                    other_sup_top.add_mirror_sup_top(candidate_superimposed_top)
                    is_mirror = True
            if is_mirror:
                # print('mirror found, skipping')
                continue

            # fixme - what to do when about the odd pairs randomH-randomH etc? they won't be found in other graphs
            # follow a rule: if this node was used before in a larger superimposed topology, than it should
            # not be in the final list (we guarantee that each node is used only once)

            sup_tops.append(candidate_superimposed_top)

    # TEST: check that each node was used only once
    all_nodes = []
    pair_count = 0
    for sup_top in sup_tops:
        [all_nodes.extend([node1, node2]) for node1, node2 in sup_top.matched_pairs]
        pair_count += len(sup_top.matched_pairs)
    # assert len(set(all_nodes)) == 2 * pair_count

    # clean the overlays by removing sub_overlays.
    # ie if all atoms in an overlay are found to be a bigger part of another overlay,
    # then that overlay is better
    print("Found altogether overlays", len(sup_tops))
    return sup_tops


def get_charges(ac_file):
    # returns
    # 1) a dictionary with charges, e.g. Item: "C17" : -0.222903
    # 2) a list of bonds

    ac_lines = open(ac_file).readlines()

    # fixme - hide hydrogens
    # ac_lines = filter(lambda l:not('h' in l or 'H' in l), ac_lines)

    # extract the atoms
    # ATOM      1  C17 MOL     1      -5.179  -2.213   0.426 -0.222903        ca
    atom_lines = filter(lambda l:l.startswith('ATOM'), ac_lines)
    atoms = [AtomNode(atomID, atomName, resName, resID, float(charge), atom_colloq) for
             _, atomID, atomName, resName, resID, x, y, z, charge, atom_colloq in
             [l.split() for l in atom_lines]]
    # fixme - add a check that all the charges come to 0 as declared in the header

    # extract the bonds, e.g.
    #     bondID atomFrom atomTo ????
    # BOND    1    1    2    7    C17  C18
    bond_lines = filter(lambda l: l.startswith('BOND'), ac_lines)
    bonds = [(bondFrom, bondTo) for _, bondID, bondFrom, bondTo, something, atomNameFrom, atomNameTo in
             [l.split() for l in bond_lines]]

    return atoms, bonds