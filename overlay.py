import hashlib
import numpy as np

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

    # def __eq__(self, other):
    #     return self.eq(other)

def rankCommonAtomsRarity(g1, g2):
    pass


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
        matched_nodes.append(frozenset([n1, n2]))
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
    sccs = []
    # fixme - TEST: create a theoretical maximum that could actually match if you ignored the topology completely
    # grow the topologies using every combination node1-node2 as the starting point
    for node1 in top1:
        for node2 in top2:
            # if not (node1.atomName == 'H9' and node2.atomName == 'H27'):
            #     continue

            # grow the topologies to see if they overlap
            matched_pairs = _overlay(node1, node2, atol=atol)
            # check if no atoms were superimposed
            if len(matched_pairs) == 0:
                continue

            # check if one topology is a subgraph of another topology
            if len(matched_pairs) == len(top1) or len(matched_pairs) == len(top2):
                print("One graph is a subgraph of another")
                return [matched_pairs, ]

            # TEST: with the list of matching nodes, check if each node was used only once,
            # the number of unique nodes should be equivalent to 2*len(common_pairs)
            all_matched_nodes = []
            for pair in matched_pairs:
                all_matched_nodes.extend(list(pair))
            assert len(matched_pairs) * 2 == len(all_matched_nodes)

            # fixme - reduce the number by removing "subgraphs"
            smatch = set(matched_pairs)
            # before adding, check if this combination already exists
            relevant_scc = True
            for smatch2 in sccs:
                # this match was already found before, or it is a subset of another match, so ignore
                if not smatch.difference(smatch2):
                    relevant_scc = False
                    break

            if relevant_scc:
                # print('len of the next sccs', len(smatch), "starting", node1.atomName, node2.atomName)
                sccs.append(smatch)

    # remove sub_overlays which were missed because of the order of generating overlays
    for scc1 in sccs[::-1]:
        for scc2 in sccs[::-1]:
            if scc1 is scc2:
                continue

            # check if overlap1 is a subset of overlap2
            if not scc1.difference(scc2):
                sccs.remove(scc1)

    # clean the overlays by removing sub_overlays.
    # ie if all atoms in an overlay are found to be a bigger part of another overlay,
    # then that overlay is better
    print("Found altogether overlays", len(sccs))
    return sccs


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