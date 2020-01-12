import hashlib
import numpy as np
import networkx as nx
import copy

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
 - optimising too much - finding matches which should not be found? how come our topology superimposition
 ignores the chirality in the case of mcl1? 
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
        # fixme - ensure that each node is characterised by its chemical info,
        # fixme - the atomId might not be unique, so check before the input data
        m.update(str(self.atomId).encode('utf-8'))
        m.update(str(self.charge).encode('utf-8'))
        # so include the number of bonds which is basically an atom type
        m.update(str(len(self.bonds)).encode('utf-8'))
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

    def __init__(self, matched_pairs=None, topology1=None, topology2=None):
        """
        @superimposed_nodes : a set of pairs of nodes that matched together
        """

        if matched_pairs is None:
            matched_pairs = []

        # TEST: with the list of matching nodes, check if each node was used only once,
        # the number of unique nodes should be equivalent to 2*len(common_pairs)
        all_matched_nodes = []
        [all_matched_nodes.extend(list(pair)) for pair in matched_pairs]
        assert len(matched_pairs) * 2 == len(all_matched_nodes)

        # todo convert to nx? some other graph theory package?
        matched_pairs.sort(key=lambda pair: pair[0].atomName)
        self.matched_pairs = matched_pairs
        self.top1 = topology1
        self.top2 = topology2
        self.mirrors = []
        # this is a set of all nodes rather than their pairs
        self.nodes = set(all_matched_nodes)
        # the uncharged sup top is a sup top that was generated when the charges were ignored
        # if the found sup_top was larger and contains the current sup top, then the two can be linked
        self.uncharged_sup_top = None


    def __len__(self):
        return len(self.matched_pairs)


    def set_tops(self, top1, top2):
        self.top1 = top1
        self.top2 = top2


    def add_node_pair(self, node_pair):
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
        assert not n1 in self.nodes and not n2 in self.nodes
        self.nodes.add(n1)
        self.nodes.add(n2)
        assert len(self.matched_pairs) * 2 == len(self.nodes)


    # fixme - switch to this
    # def __copy__(self):
    #     newone = type(self)()
    #     newone.__dict__.update(self.__dict__)
    #     return newone


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
        # walk through the superimposed topologies
        # and move the atom-atom pairs that suffer from being the same
        for node1, node2 in self.matched_pairs[::-1]:
            if not node1.eq(node2, atol=atol):
                # remove this pair
                self.matched_pairs.remove([node1, node2])
                print('removed a pair due to the not-matching charges', node1, node2)


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


    def is_consistent_with(self, other_suptop):
        """
        There should be a minimal overlap of at least 1 node.

        There is no noe pair (A=B) in this sup top such that (A=C) or (B=C) exists in other.
        """

        # fixme - check for minimal partial overlap - this is more or a test than anything

        for node1, node2 in self.matched_pairs:
            for nodeA, nodeB in other_suptop.matched_pairs:
                if (node1 is nodeA) and not (node2 is nodeB):
                    return False
                elif (node2 is nodeB) and not (node1 is nodeA):
                    return False

        return True


    def merge(self, other_suptop):
        """
        Absorb the other suptop by adding all the node pairs that are not present
        in the current sup top.

        WARNING: ensure that the other suptop is consistent with this sup top.
        """
        assert self.is_consistent_with(other_suptop)

        for pair in other_suptop.matched_pairs:
            # check if this pair is present
            if not self.contains(pair):
                n1, n2 = pair
                if self.contains_node(n1) or self.contains_node(n2):
                    raise Exception('already uses that node')
                self.add_node_pair(pair)

        # check if duplication occured, fixme - temporary



    def contains_node(self, node):
        # checks if this node was used in this overlay
        if len(self.nodes.intersection(set([node,]))) == 1:
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
        print("a mirror sup top added")
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


def _overlay(n1, n2, n1_parent=None, n2_parent=None, sup_top=None, nxgl=None, nxgr=None, atol=0):
    """
    n1 should be from one graph, and n2 should be from another.

    If n1 and n2 are the same, we will be traversing through both graphs, marking the jointly travelled areas.
    RETURN: the maximum common overlap with the n1 and n2 as the starting conditions.

    Here, recursively the graphs of n1 and of n2 will be explored as they are the same.
    """

    if sup_top == None:
        sup_top = SuperimposedTopology()
        nxgl = nx.Graph()
        nxgr = nx.Graph()

    # if either of the nodes has already been matched, ignore this potential match
    for pair in sup_top.matched_pairs:
        if n1 in pair or n2 in pair:
            return sup_top

    # if n1.atomName == 'C13' and n2.atomName == 'C35':
    #     print('found')

    # if the two nodes are "the same"
    if n1.eq(n2, atol=atol):
        # update the graphs (which allows us to access finding circles)
        # nxgl.add_node(n1)
        # nxgr.add_node(n2)

        # if the parents are also present, connect them
        # nxgl.add_edge(n1_parent, n1)
        # nxgl.add_edge(n2_parent, n2)
        # connect the other nodes that are present in the bonds that could close circles
        # ie if n1 closes a circle with some nX in the the graph, but n2 does not have that bond
        # then they represent a different thing,
        # Summary: connect n1 and n2 to any nodes in the graph that they are bound to
        # for bonded_to_n1 in n1.bonds:
        #     if bonded_to_n1 in nxgl:
        #         nxgl.add_edge(n1, bonded_to_n1)
        # for bonded_to_n2 in n2.bonds:
        #     if bonded_to_n2 in nxgr:
        #         nxgr.add_edge(n2, bonded_to_n2)

        # fixme - it is possible (see mcl1 case) to find a situation where a larger match that is found
        # is actually the wrong match, add one more criteria:
        # if the new atom completes a ring, but the other atom does not, then the two are heading in the wrong dimension
        # fixme - add tests which check of features like rings/double rings etc and see if they are found
        # in both superimpositions and yet "not present" in the sup-top, meaning that the sup-top is wrong
        # if len(nx.cycle_basis(nxgl)) != len(nx.cycle_basis(nxgr)):
        #     # clearly the newly added edge changes the circles, and so it is not equivalent
        #     # even though it might appear like it (mcl1 case)
        #     nxgl.remove_node(n1)
        #     nxgr.remove_node(n2)
        #     return sup_top

        # append both nodes as a pair to ensure that we keep track of the mapping
        # having both nodes appended also ensure that we do not revisit/readd neither n1 and n2
        sup_top.add_node_pair((n1, n2))

        # continue traversing
        # try every possible pathway for creating an overlap,
        # and as the winner, pick the one that is the largest
        # fixme - this should find mirrors which can be tackled right here
        solutions = []
        for n1bonded_node in n1.bonds:
            for n2bonded_node in n2.bonds:
                # copy the sup_top and update the list
                # fixme - overload the __copy__ method for this
                copy_sup_top = copy.copy(sup_top)
                copy_sup_top.matched_pairs = copy.copy(sup_top.matched_pairs)
                copy_sup_top.nodes = copy.copy(sup_top.nodes)

                solution = _overlay(n1bonded_node, n2bonded_node,
                         n1_parent=n1, n2_parent=n2,
                         nxgl=copy.copy(nxgl), nxgr=copy.copy(nxgr),  # fixme - do I need a deep copy for the graph? test
                         sup_top=copy_sup_top, # a shallow copy here should be just fine
                            atol=atol)
                solutions.append(solution)

        # fixme - you should merge all solutions that are consistent? you should not take the largest solution
        # you should take the different pathways that you discovered that are consistent,
        # ie how to combine the multiple pathways?
        # as long as the two superimpositions are consistent, then we should be fine with the solution,
        # so now I have to combine the different walks, and then see the emerging consistent combinations,
        # return largest_solution

        # merge all the paths that are consistent with each other
        # fixme - optimise,
        for sol1 in solutions:
            for sol2 in solutions[::-1]:
                if sol1 is sol2:
                    continue

                if sol1.is_consistent_with(sol2):
                    # print("merging")
                    # join sol2 and sol1 because they're consistent
                    sol1.merge(sol2)
                    # remove sol2 from the solutions:
                    solutions.remove(sol2)

        # take the bigger out of the two,
        # print('solutions lengths', [len(s) for s in solutions])
        # fixme - check if there is more than two - ie the mirror case
        # so you hve the right solution, so merge it with the current browsed solution
        for s in solutions:
            if len(s) == max([len(s2) for s2 in solutions]):
                # sup top should always be compatible with the largest solution sup_top
                # fixme - how is it possible that the next solution is not consistent with the starting sup top?
                if sup_top.is_consistent_with(s):
                    sup_top.merge(s)
                else:
                    print('found an inconsistent largest result')

        if len(solutions) == 0:
            raise Exception('no solution? ')

    return sup_top


def superimpose_topologies(top1, top2, atol):
    # fixme replace with theoretical max
    large_value = 99999
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

    sup_tops_no_charges = _superimpose_topologies(top1, top2, atol=large_value)

    # apply_charges(sup_tops_no_charges, atol=atol)

    # fixme - remove the hydrogens without attached heavy atoms

    # sup_tops_charges = _superimpose_topologies(top1, top2, atol=atol)

    # link_components_to_supercomponents(sup_tops_charges, sup_tops_no_charges)

    # fixme - currently you hold on to "mirror superimpositions" but there is no need,
    # if they can be solved with the global uncharged sup top

    # resolve_sup_top_multiple_match(sup_tops_charges)
    # find_symmetric_components(sup_tops_charges, sup_tops_no_charges)
    # sup_top_correct_chirality(sup_tops_charges, sup_tops_no_charges, atol=atol)

    return sup_tops_no_charges


def _superimpose_topologies(top1, top2, atol):
    """
    atol 1 means the charge of the atom can be up to 1 electron different,
    as long as the atom has the same type
    """
    sup_tops = []
    # fixme - TEST: create a theoretical maximum that could actually match if you ignored the topology completely
    # grow the topologies using every combination node1-node2 as the starting point
    for node1 in top1:
        for node2 in top2:
            # stick to a case that should work
            if not (node1.atomName == 'C3' and node2.atomName == 'C25'):
                continue

            # fixme - optimisation
            # you don't have to start from matches that have already been found,
            # ie if something discovers that match, it should have explored all other ways in which
            # the two topologies can be reassembled,

            # grow the topologies to see if they overlap
            # fixme - do you still need to set up top1 and top2?
            candidate_superimposed_top = _overlay(node1, node2, atol=atol)
            # ignore if no atoms were superimposed
            if len(candidate_superimposed_top) == 0:
                continue

            candidate_superimposed_top.set_tops(top1, top2)
            candidate_superimposed_top.is_subgraph_of_global_top()

            # check if there is a pair C5==C27 and check where it is heading
            for n1, n2 in candidate_superimposed_top.matched_pairs:
                if n1.atomName == 'C7' and n2.atomName == 'C29' and \
                        len(candidate_superimposed_top.matched_pairs) > 30:
                    print(candidate_superimposed_top.matched_pairs)
                    print("FOUND CORRECT MATCH: len %d :" % len(candidate_superimposed_top.matched_pairs),
                          "Started with", node1.atomName, 'and', node2.atomName,'\n',
                          'name ' + ' '.join([node1.atomName.upper() for node1, _ in candidate_superimposed_top.matched_pairs]),
                          '\nto\n',
                          'name ' + ' '.join([node2.atomName.upper() for _, node2 in candidate_superimposed_top.matched_pairs]),
                          '\n\n')
                    break

                    # check if there is a pair C5==C27 and check where it is heading
            # for n1, n2 in candidate_superimposed_top.matched_pairs:
            #     if n1.atomName == 'C5' and n2.atomName == 'C27' and \
            #             len(candidate_superimposed_top.matched_pairs) == 42:
            #         print("found a 42 with the right atoms? ")
            #         print("Superimposed topology: len %d :" % len(candidate_superimposed_top.matched_pairs),
            #               'name ' + ' '.join([node1.atomName.upper() for node1, _ in
            #                                   candidate_superimposed_top.matched_pairs]),
            #               '\nto\n',
            #               'name ' + ' '.join([node2.atomName.upper() for _, node2 in
            #                                   candidate_superimposed_top.matched_pairs]))
            #         break

            # check if this superimposed topology was found before, if so, ignore
            sup_top_seen = False
            for other_sup_top in sup_tops:
                if other_sup_top.eq(candidate_superimposed_top):
                    sup_top_seen = True
                    print("The next candidate superimposed topology was seen before, len:", len(candidate_superimposed_top.matched_pairs))
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
                        #fixme
                        if len(candidate_superimposed_top.matched_pairs) > 37:
                            print("A larger sup top was found that already uses some of these nodes,"
                                  " ignoring this sup top, len", len(candidate_superimposed_top.matched_pairs))
                        ignore_cand_sup_top = True
                        break
                elif len(sup_top.matched_pairs) < len(candidate_superimposed_top.matched_pairs):
                    # this sup_top has fewer elements than candidate sup top
                    if sup_top.contains_any_node_from(candidate_superimposed_top):
                        # and there is an overlap, delete the smaller sup top
                        if len(candidate_superimposed_top.matched_pairs) > 37:
                            print('This new candidate sup top removes the previous sup top'
                                  'that is smaller, len',
                                  len(candidate_superimposed_top.matched_pairs))
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
                print('mirror found, skipping, sup top len', len(candidate_superimposed_top.matched_pairs))
                continue

            # fixme - what to do when about the odd pairs randomH-randomH etc? they won't be found in other graphs
            # follow a rule: if this node was used before in a larger superimposed topology, than it should
            # not be in the final list (we guarantee that each node is used only once)

            sup_tops.append(candidate_superimposed_top)

    # if there are only hydrogens superimposed without a connection to any heavy atoms, ignore these too
    for sup_top in sup_tops[::-1]:
        all_hydrogens = True
        for node1, _ in sup_top.matched_pairs:
            if not node1.atom_colloq.upper().startswith('H'):
                all_hydrogens = False
                break
        if all_hydrogens:
            print("Removing sup top because only hydrogens found", sup_top.matched_pairs)
            sup_tops.remove(sup_top)


    # TEST: check that each node was used only once
    all_nodes = []
    pair_count = 0
    for sup_top in sup_tops:
        [all_nodes.extend([node1, node2]) for node1, node2 in sup_top.matched_pairs]
        pair_count += len(sup_top.matched_pairs)
    # fixme
    # assert len(set(all_nodes)) == 2 * pair_count

    # TEST: check that the nodes on the left are always from topology 1 and the nodes on the right are always from top2
    for sup_top in sup_tops:
        for node1, node2 in sup_top.matched_pairs:
            assert node1 in list(top1)
            assert node2 in list(top2)

    # clean the overlays by removing sub_overlays.
    # ie if all atoms in an overlay are found to be a bigger part of another overlay,
    # then that overlay is better
    print("Found altogether overlays", len(sup_tops))
    # fixme - return other info
    return sup_tops


def apply_charges(sup_tops_no_charges, atol=0):
    for sup_top in sup_tops_no_charges:
        sup_top.refineAgainstCharges(atol=atol)


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
                print('found same left', sup_top1.matched_pairs, 'with', sup_top1.matched_pairs)

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
                print('found same right', sup_top1.matched_pairs)

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
        print("multiple match winner is", same_left_sup_top_list[winner_index].matched_pairs)

        # remove the losers
        for index, worse_match in enumerate(same_left_sup_top_list):
            if index == winner_index:
                pass
            else:
                print("Removing a worse match", worse_match.matched_pairs)
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
                            print("found assymetry", node1.atomName, node2.atomName,
                                  "due to", bond1.atomName, bond2.atomName)
            pass
        pass

    # add a test against the overall match (global match that ignores the charges)
    # FIXME - if two superimposed components come from two different places in the global map, then something's off
    # particularly, it could help with chirality - if with respect to the global match,
    # a local sup top travels in the wrong direction, then we have a clear issue


def find_symmetric_components(sup_top_charges, sup_top_no_charges):
    # in the case of mcl1
    pass


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