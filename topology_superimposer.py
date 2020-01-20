
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
from MDAnalysis.analysis.distances import distance_array
import hashlib
import numpy as np
import networkx as nx
import copy


class AtomNode:
    def __init__(self, atomId, atomName, resName, resId, charge, atom_type):
        self.atomId = atomId
        self.atomName = atomName
        self.resName = resName
        self.resId = resId
        self.charge = charge
        self.type = atom_type.upper()
        self.bonds = set()


    def set_coords(self, coords):
        assert len(coords) == 3
        assert coords.dtype == np.dtype('float32')
        self.coords = coords


    def __hash__(self):
        m = hashlib.md5()
        # fixme - ensure that each node is characterised by its chemical info,
        # fixme - the atomId might not be unique, so check before the input data
        m.update(self.atomName.encode('utf-8'))
        m.update(str(self.charge).encode('utf-8'))
        # so include the number of bonds which is basically an atom type
        m.update(str(len(self.bonds)).encode('utf-8'))
        return int(m.hexdigest(), 16)

    def __str__(self):
        return self.atomName

    def __repr__(self):
        return self.atomName

    def bindTo(self, other):
        self.bonds.add(other)
        other.bonds.add(self)

    def eq(self, atom, atol=0):
        """
        What does it mean that two atoms are the same? They are the same type and charge.
        5 % tolerance by default
        """
        if self.type == atom.type and \
                np.isclose(self.charge, atom.charge, atol=atol):
            return True

        return False


    def sameType(self, atom):
        if self.type == atom.type:
            return True

        return False


    def __deepcopy__(self, memodict={}):
        # https://stackoverflow.com/questions/1500718/how-to-override-the-copy-deepcopy-operations-for-a-python-object
        # it is a shallow copy, as this object is "immutable"
        return self


class SuperimposedTopology:
    """
    SuperimposedTopology contains in the minimal case two sets of nodes S1 and S2, which
    are paired together and represent a strongly connected component.

    However, it can also represent the symmetrical versions that were superimposed.
    """

    def __init__(self, topology1=None, topology2=None):
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
        self.weird_symmetries = []
        # this is a set of all nodes rather than their pairs
        self.nodes = set(all_matched_nodes)
        # the uncharged sup top is a sup top that was generated when the charges were ignored
        # if the found sup_top was larger and contains the current sup top, then the two can be linked
        self.uncharged_sup_top = None
        self.nodes_added_log = []


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


    def print_summary(self):
        print("Topology Pairs ", len(self.matched_pairs), "Mirror Number", len(self.mirrors))

        # print the match
        # for strongly_connected_component in overlays:
        #     print("Strongly Connected Component, length:", len(strongly_connected_component))
        # for atom_from, atom_to in strongly_connected_component:
        #     print('Bound', atom_from.atomName, atom_to.atomName)

        # extract all the unique nodes from the pairs
        all_matched_nodes = set()
        print("Superimposed topology: len %d :" % len(self.matched_pairs),
              'name ' + ' '.join([node1.atomName.upper() for node1, _ in self.matched_pairs]),
              '\nto\n',
              'name ' + ' '.join([node2.atomName.upper() for _, node2 in self.matched_pairs]))
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


    def remove_node_pair(self, node_pair):
        assert len(node_pair) == 2
        self.matched_pairs.remove(node_pair)
        # remove from the current set
        self.nodes.remove(node_pair[0])
        self.nodes.remove(node_pair[1])

        # update the log to understand the order in which it was created
        self.nodes_added_log.append(("Removed", node_pair))


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
            print("This graph is a equivalent to the topology")
            return True

        return False


    def rmsd(self):
        """
        For each pair take the distance, and then get rmsd, so root(mean(square(devisation)))
        """

        assert len(self.matched_pairs) > 0

        sq_dsts = []
        for nodeA, nodeB in self.matched_pairs:
            dst = distance_array(np.array([nodeA.coords, ]), np.array([nodeB.coords, ]))[0]
            sq_dsts.append(dst**2)
        return np.sqrt(np.mean(sq_dsts))


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
        assert not n1 in self.nodes and not n2 in self.nodes, (n1, n2)
        self.nodes.add(n1)
        self.nodes.add(n2)
        assert len(self.matched_pairs) * 2 == len(self.nodes)

        # update the log to understand the order in which this sup top was created
        self.nodes_added_log.append(("Added", node_pair))

        # fixme - ideally sup top would internally manage/use the networkx Graph
        # update the networkx graphs
        #self.nxlg.add_node(n1)
        # self.nxrg.add_node(n2)
        #


    def __copy__(self):
        # https://stackoverflow.com/questions/1500718/how-to-override-the-copy-deepcopy-operations-for-a-python-object
        newone = type(self)()
        newone.__dict__.update(self.__dict__)
        # make a shallow copy of the arrays
        newone.matched_pairs = copy.copy(self.matched_pairs)
        newone.nodes = copy.copy(self.nodes)
        newone.nodes_added_log = copy.copy(self.nodes_added_log)
        return newone


    def findMirrorChoices(self):
        """
        For each pair (A1, B1) find all the other options in the mirrors where (A1, B2)
        # ie Ignore (X, B1) search, if we repair from A to B, then B to A should be repaired too
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


    def addWeirdSymmetry(self, weird_symmetry):
        """
        This means that there is another way to traverse and overlap the two molecules,
        but that the self is better (e.g. lower rmsd) than the other one
        """
        self.weird_symmetries.append(weird_symmetry)


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
                    A1_BX_dst = distance_array(np.array([A1.coords, ]), np.array([BX.coords, ]))[0]
                    if A1_BX_dst < closest_dst:
                        closest_dst = A1_BX_dst
                        closest_bx = BX
                        closest_a1 = A1

            # across all the possible choices, found the best match now:
            blacklisted_bxs.append(closest_bx)
            shortest_dsts.append(closest_dst)
            print(closest_a1.atomName, 'is matching best with', closest_bx.atomName)

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
        for node1, node2 in self.matched_pairs[::-1]:
            if not node1.eq(node2, atol=atol):
                # remove this pair
                self.matched_pairs.remove((node1, node2))
                # get hydrogens attached to the heavy atoms
                node1_hydrogens = filter(lambda a: a.atomName.upper().startswith('H'), node1.bonds)
                node2_hydrogens = filter(lambda a: a.atomName.upper().startswith('H'), node2.bonds)
                # fixme - ideally the hydrogens would match each other before being removed
                for n1_hyd in node1_hydrogens:
                    # find the matching n2 hydrogen and remove them
                    n1_pair_removed = False
                    for n2_hyd in node2_hydrogens:
                        try:
                            self.matched_pairs.remove((n1_hyd, n2_hyd))
                            print('Removed lonely hydrogen pair', (n1_hyd.atomName, n2_hyd.atomName))
                            n1_pair_removed = True
                            removed_pairs.append((n1_hyd, n2_hyd))
                            break
                        except:
                            pass
                removed_pairs.append((node1, node2))
                print('removed a pair due to the not-matching charges', node1.atomName, node2.atomName)
        return removed_pairs


    def is_consistent_with(self, other_suptop):
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
            for nodeA, nodeB in other_suptop.matched_pairs:
                if (node1 is nodeA) and not (node2 is nodeB):
                    return False
                elif (node2 is nodeB) and not (node1 is nodeA):
                    return False

        # ensure there is at least one common pair
        if self.count_common_node_pairs(other_suptop) == 0:
            return False

        # check if each sup top has the same number of cycles
        # fixme - not sure?
        selfG1, selfG2 = self.getNxGraphs()
        self_cycles1, self_cycles2 = len(nx.cycle_basis(selfG1)), len(nx.cycle_basis(selfG2))
        if self_cycles1 != self_cycles2:
            raise Exception('left G has a different number of cycles than right G')

        otherG1, otherG2 = other_suptop.getNxGraphs()
        other_cycles1, other_cycles2 = len(nx.cycle_basis(otherG1)), len(nx.cycle_basis(otherG2))
        if other_cycles1 != other_cycles2:
            raise Exception('left G has a different number of cycles than right G')

        # if self_cycles1 != other_cycles1:
        #     # rings should be created during the traversal
        #     # ie new cycles are important for the right topology superimposition
        #     # so merging different graphs should not lead to new rings
        #     print('merging graphs should not take place')
        #     # clean up
        #     return False

        # check if they have a bigger number of cycles after merging
        # todo
        mergedG1 = nx.compose(selfG1, otherG1)
        mergedG2 = nx.compose(selfG2, otherG2)
        mergedG1_cycle_num = len(nx.cycle_basis(mergedG1))
        mergedG2_cycle_num = len(nx.cycle_basis(mergedG2))
        if mergedG1_cycle_num != mergedG2_cycle_num:
            return False

        # if mergedG1_cycle_num != self_cycles1:
        #     return False

        return True


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
            for bonded_to_nA in nA.bonds:
                if bonded_to_nA in gl:
                    gl.add_edge(nA, bonded_to_nA)
            for bonded_to_nB in nB.bonds:
                if bonded_to_nB in gr:
                    gr.add_edge(nB, bonded_to_nB)

        return gl, gr


    def getCircles(self):
        gl, gr = self.getNxGraphs()
        gl_circles = nx.cycle_basis(gl)
        gr_circles = nx.cycle_basis(gr)
        return gl_circles, gr_circles


    def getCircleNumber(self):
        gl_circles, gr_circles = self.getCircles()
        return len(gl_circles), len(gr_circles)


    def sameCircleNumber(self):
        gl_num, gr_num = self.getCircleNumber()
        if gl_num == gr_num:
            return True

        return False


    def merge(self, other_suptop):
        """
        Absorb the other suptop by adding all the node pairs that are not present
        in the current sup top.

        WARNING: ensure that the other suptop is consistent with this sup top.
        """
        assert self.is_consistent_with(other_suptop)

        # print("About the merge two sup tops")
        # self.print_summary()
        # other_suptop.print_summary()

        for pair in other_suptop.matched_pairs:
            # check if this pair is present
            if not self.contains(pair):
                n1, n2 = pair
                if self.contains_node(n1) or self.contains_node(n2):
                    raise Exception('already uses that node')
                self.add_node_pair(pair)

        self.nodes_added_log.append(("merged with", copy.deepcopy(other_suptop.nodes_added_log)))

        # check if duplication occured, fixme - temporary



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


    def contains_all(self, other_sup_top):
        for pair in other_sup_top.matched_pairs:
            if not self.contains(pair):
                return False

        return True


    def difference(self, sup_top):
        return set(self.matched_pairs).difference(set(sup_top.matched_pairs))


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

        # this should not be applied on the same topology match
        if self.eq(other_sup_top):
            raise Exception("They are already the same, cannot be a mirror")

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
        # print("a mirror sup top added")
        self.mirrors.append(mirror_sup_top)


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


verbose_log = True
def log(*args):
    if verbose_log:
        print(*args)


def _overlay(n1, n2, sup_top=None):
    """
    n1 should be from one graph, and n2 should be from another.

    If n1 and n2 are the same, we will be traversing through both graphs, marking the jointly travelled areas.
    RETURN: Returns a list of topologies (ie solutions)

    Recursively traverse the graphs of n1 and of n2 at the same time.

    fixme: return should keep track of symmetries
    - while returning the current system has to make sense of the different returning journies, this means
    that symmetry might show up here. For example, the same ring can be traversed in two different ways,
    therefore, we should be possible to continue returning and forming the different "variants".
    Currently, only one way is chosen, despite several different candidates.
    One way to store them it to make the entire different topologies, in a way that makes sense,
    right now we return with only one of the "symmetries", while others are being discovered
    when searching through the other node pairs.
    """
    if sup_top == None:
        sup_top = SuperimposedTopology()

    # if either of the nodes has already been matched, ignore
    if sup_top.contains_any_node([n1, n2]):
        return [sup_top, ]

    # if the two nodes are "the same"
    # fixme - remove the charges from here, you're not using them anyway for now
    if n1.sameType(n2):
        # Alternative to checking for cycles:
        # check if both n1 and n2 contain a bonded atom that already is in the topology
        # this way we know both create some connection back to the molecule
        # fixme - did not seem to work?
        # n1_common_nodes_tot = sup_top.count_common_nodes(n1.bonds)
        # n2_common_nodes_tot = sup_top.count_common_nodes(n2.bonds)
        # if n1_common_nodes_tot != n2_common_nodes_tot:
        #     return [sup_top, ]

        nxgl, nxgr = sup_top.getNxGraphs()
        nxgl_cycles_num, nxgr_cycles_num = len(nx.cycle_basis(nxgl)), len(nx.cycle_basis(nxgr))
        assert nxgl_cycles_num == nxgr_cycles_num

        log("Adding ", (n1, n2), "in", sup_top.matched_pairs)

        # append both nodes as a pair to ensure that we keep track of the mapping
        # having both nodes appended also ensure that we do not revisit/readd neither n1 and n2
        sup_top.add_node_pair((n1, n2))

        # fixme - it is possible (see mcl1 case) to find a situation where a larger match that is found
        # is actually the wrong match, add one more criteria:
        # if the new atom completes a ring, but the other atom does not, then the two are heading in the wrong dimension
        # fixme - add tests which check of features like rings/double rings etc and see if they are found
        # in both superimpositions and yet "not present" in the sup-top, meaning that the sup-top is wrong
        nxgl_after, nxgr_after = sup_top.getNxGraphs()
        nxgl_cycles_num_after, nxgr_cycles_num_after = len(nx.cycle_basis(nxgl_after)), len(nx.cycle_basis(nxgr_after))
        if nxgl_cycles_num_after != nxgr_cycles_num_after:
            # clearly the newly added edge changes the circles, and so it is not equivalent
            # even though it might appear like it (mcl1 case)
            sup_top.remove_node_pair((n1, n2))
            log("Removing pair because one creates a cycle and the other does not", (n1, n2))
            return None

        if nxgl_cycles_num_after > nxgl_cycles_num:
            log("Added a new cycle to both nxgr and nxgl by adding the pair", (n1, n2))

        # try every possible pathway
        all_solutions = []
        for n1bonded_node in n1.bonds:
            for n2bonded_node in n2.bonds:
                # if either of the nodes has already been matched, ignore
                if sup_top.contains_any_node([n1bonded_node, n2bonded_node]):
                    # the advantage of this approach is that we do not have to evaluate extra returned sup_tops
                    continue

                # a copy of the sup_top is needed because the traversal can take place
                # using different pathways
                bond_solutions = _overlay(n1bonded_node, n2bonded_node, sup_top=copy.copy(sup_top))
                if bond_solutions is None:
                    continue
                all_solutions.extend(bond_solutions)
                # fixme - when you have a mirror like ester O1-O2 matches, then you could store them differently here

        # if this was the last atom traversed together, return the current sup top
        if len(all_solutions) == 0:
            return [sup_top, ]

        # sort in the descending order
        all_solutions.sort(key=lambda st:len(st), reverse=True)

        # combine the different walks, the walks that are smaller can be thrown away?
        # we try as many mergings as possible?
        for sol1 in all_solutions:
            for sol2 in all_solutions[::-1]:
                if sol1 is sol2:
                    continue

                if sol1.eq(sol2):
                    log("Found the same solution and removing, solution", sol1.matched_pairs)
                    all_solutions.remove(sol2)
                    continue

                if sol1.is_consistent_with(sol2):
                    # print("merging, current pair", (n1, n2))
                    # join sol2 and sol1 because they're consistent
                    g1, g2 = sol1.getNxGraphs()
                    assert len(nx.cycle_basis(g1)) == len(nx.cycle_basis(g2))
                    g3, g4 = sol2.getNxGraphs()
                    assert len(nx.cycle_basis(g3)) == len(nx.cycle_basis(g4))

                    print("Will merge", sol1, 'and', sol2)
                    sol1.merge(sol2)
                    # remove sol2 from the solutions:
                    all_solutions.remove(sol2)

        # after all the merges, return the best matches only,
        # the other mergers are wrong (and there are smaller)
        largest_sol_size = max([len(s2) for s2 in all_solutions])
        return list(filter(lambda st:len(st) == largest_sol_size, all_solutions))

    return [sup_top, ]


def superimpose_topologies(top1, top2, atol, useCharges=True, useCoords=True):
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
    if useCharges:
        ensure_charges_match(sup_tops_no_charges, atol=atol)

    if useCoords:
        for sup_top in sup_tops_no_charges:
            sup_top.correct_for_coordinates()

    # fixme - remove the hydrogens without attached heavy atoms

    # resolve_sup_top_multiple_match(sup_tops_charges)
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
            # fixme - optimisation
            # you don't have to start from matches that have already been found,
            # ie if something discovers that match, it should have explored all other ways in which
            # the two topologies can be reassembled,

            # for testing speed C11=C33
            # (ca_C4, ca_C28) (ca_C8, ca_C24)
            # if (node1.atomName == 'C4' and node2.atomName == 'C28') or \
            #     (node1.atomName == 'C8' and node2.atomName == 'C24'):
            #     continue

            # if not (node1.atomName == 'C4' and node2.atomName == 'C28'):
            #     continue

            # grow the topologies to see if they overlap
            # fixme - do you still need to set up top1 and top2?
            candidate_superimposed_tops = _overlay(node1, node2)
            if candidate_superimposed_top is None:
                continue

            # _overlay returns a list of solutions, which can be traversed with that specific initial
            # starting two nodes on the two topologies.
            # multiple paths could have been taken
            if len(candidate_superimposed_tops) == 1:
                candidate_superimposed_top = candidate_superimposed_tops[0]
                # ignore if no atoms were superimposed
                if len(candidate_superimposed_top) == 0:
                    continue
            elif len(candidate_superimposed_tops) > 1:
                # so multiple different representation of the same molecule was found
                # this means that some kind of symmetry was found in the topology,
                # for example, cyclehexane connected via one atom can be traversed both ways
                # here now we have to decide which of them is better
                # fixme - simple approach - uses coordinates to decide, which part is better
                best_suptop = None
                min_rmsd = 9999999
                rmsd = None
                for suptop in candidate_superimposed_tops:
                    # use the avg dst after the correction to understand which one is better,
                    # and assign the worse
                    # fixme - avg dst would ideally not be np.NaN
                    avg_dst = suptop.correct_for_coordinates()
                    rmsd = suptop.rmsd()
                    # get a global match (RMSD)
                    if rmsd < min_rmsd:
                        best_suptop = suptop
                        min_rmsd = rmsd

                candidate_superimposed_top = best_suptop

                for suptop in candidate_superimposed_tops:
                    if suptop is candidate_superimposed_top:
                        continue

                    candidate_superimposed_top.addWeirdSymmetry(suptop)

            candidate_superimposed_top.set_tops(top1, top2)
            candidate_superimposed_top.is_subgraph_of_global_top()

            # check if this superimposed topology was found before, if so, ignore
            sup_top_seen = False
            for other_sup_top in sup_tops:
                if other_sup_top.eq(candidate_superimposed_top):
                    sup_top_seen = True
                    # print("The next candidate superimposed topology was seen before, len:", len(candidate_superimposed_top.matched_pairs))
                    break
            if sup_top_seen:
                continue

            # check if this superimposed topology is a mirror of one that already exists
            # fixme the order matters in this place
            # fixme - what if the mirror has a lower rmsd match? in that case, pick that mirror here
            is_mirror = False
            for other_sup_top in sup_tops:
                if other_sup_top.is_mirror_of(candidate_superimposed_top):
                    other_sup_top.add_mirror_sup_top(candidate_superimposed_top)
                    is_mirror = True
            if is_mirror:
                # print('mirror found, skipping, sup top len', len(candidate_superimposed_top.matched_pairs))
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
                        if len(candidate_superimposed_top.matched_pairs) > 37:
                            print('This new candidate sup top removes the previous sup top'
                                  'that is smaller, new len: ',
                                  len(candidate_superimposed_top.matched_pairs),
                                  'old len:', len(sup_top))
                        sup_tops.remove(sup_top)
                else:
                    # the two sup tops could be of the same length and yet traverse different atoms
                    # e.g. case: mcl1_l17l9, where due to the symmetry, two different traversals
                    # of the same length are found

                    # fixme - note that we already checked if it is a mirror,
                    # however, these two would ideally be chained together closer
                    # to avoid any future bugs

                    # check if there is an overlap
                    if sup_top.contains_any_node_from(candidate_superimposed_top):
                        # there is a partial overlap, so two different ways to score these
                        # you could call them "symmetries". Here we have to pick
                        # which is the "worse symmetry",
                        # let us use atom coordinates to score them
                        if sup_top.rmsd() < candidate_superimposed_top.rmsd():
                            # the
                            sup_top.addWeirdSymmetry(candidate_superimposed_top)
                            ignore_cand_sup_top = True
                        else:
                            candidate_superimposed_top.addWeirdSymmetry(sup_top)
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


def ensure_charges_match(sup_tops_no_charges, atol=0):
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
                atom.set_coords(pdb_atom.position)
                found_match = True
                break
        if not found_match:
            print("Did not find atom?", atom.atomName)
            raise Exception("wait a minute")