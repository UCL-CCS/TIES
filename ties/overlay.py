import warnings

import itertools

import copy

import logging
import numpy as np

logger = logging.getLogger(__name__)


def _overlay(
    n1,
    n2,
    parent_n1,
    parent_n2,
    bond_types,
    suptop,
    ignore_coords=False,
    use_element_type=True,
    exact_coords_cue=False,
    weights=(1, 0.01),
):
    """
    Jointly and recursively traverse the molecule while building up the suptop.

    If n1 and n2 are the same, we will be traversing through both graphs, marking the jointly travelled areas.

    Return the topology of the common substructure between the two molecules.

    *n1 from the left molecule,
    *n2 from the right molecule
    """

    # ignore if either of the nodes is part of the suptop
    if suptop.contains_node(n1) or suptop.contains_node(n2):
        return None

    if use_element_type and not n1.same_element(n2):
        return None

    # make more specific, ie if "use_specific_type"
    if not use_element_type and not n1.same_type(n2):
        return None

    # Check for cycles
    # if a new cycle is created by adding this node,
    # then the cycle should be present in both, left and right ligand
    safe = True
    # if n1 is linked with node in suptop other than parent
    for b1 in n1.bonds:
        # if this bound atom is not a parent and is already a suptop
        if b1.atom != parent_n1 and suptop.contains_node(b1.atom):
            safe = False  # n1 forms cycle, now need to check n2
            for b2 in n2.bonds:
                if b2.atom != parent_n2 and suptop.contains_node(b2.atom):
                    # b2 forms cycle, now need to check it's the same in both
                    if suptop.contains((b1.atom, b2.atom)):
                        safe = True
                        break
            if not safe:  # only n1 forms a cycle
                break
    if not safe:  # either only n1 forms cycle or both do but different cycles
        return None

    # now the same for any remaining unchecked bonds in n2
    safe = True
    for b2 in n2.bonds:
        if b2.atom != parent_n2 and suptop.contains_node(b2.atom):
            safe = False
            for b1 in n1.bonds:
                if b1.atom != parent_n1 and suptop.contains_node(b1.atom):
                    if suptop.contains((b1.atom, b2.atom)):
                        safe = True
                        break
            if not safe:
                break
    if not safe:
        return None

    # check if the cycle spans multiple cycles present in the left and right molecule,
    if suptop.cycle_spans_multiple_cycles():
        logger.debug("Found a cycle spanning multiple cycles")
        return None

    # logger.debug(f"Adding {(n1, n2)} to suptop.matched_pairs")

    # all looks good, create a new copy for this suptop
    suptop = copy.copy(suptop)

    # append both nodes as a pair to ensure that we keep track of the mapping
    # having both nodes appended also ensure that we do not revisit/read neither n1 and n2
    suptop.add_node_pair((n1, n2))
    if not (parent_n1 is parent_n2 is None):
        # fixme - adding a node pair should automatically take care of the bond, maybe using inner data?
        # fixme why is this link different than a normal link?
        suptop.link_with_parent((n1, n2), (parent_n1, parent_n2), bond_types)

    # the extra bonds are legitimate
    # so let's make sure they are added
    # fixme: add function get_bonds_without_parent? or maybe make them "subtractable" even without the type
    # for this it would be enough that the bonds is an object too, it will make it more managable
    # bookkeeping? Ideally adding "add_node_pair" would take care of this
    for n1_bonded in n1.bonds:
        # ignore left parent
        if n1_bonded.atom is parent_n1:
            continue
        for n2_bonded in n2.bonds:
            # ignore right parent
            if n2_bonded.atom is parent_n2:
                continue

            # if the pair exists, add a bond between the two pairs
            if suptop.contains((n1_bonded.atom, n2_bonded.atom)):
                # fixme: this linking of pairs should also be corrected
                # 1) add "pair" as an object rather than a tuple (n1, n2)
                # 2) this always has to happen, ie it is impossible to find (n1, n2)
                # ie make it into a more sensible method,
                # fixme: this does not link pairs?
                suptop.link_pairs(
                    (n1, n2),
                    [
                        (
                            (n1_bonded.atom, n2_bonded.atom),
                            (n1_bonded.type, n2_bonded.type),
                        ),
                    ],
                )

    # fixme: sort so that heavy atoms go first
    p1_bonds = n1.bonds.without(parent_n1)
    p2_bonds = n2.bonds.without(parent_n2)
    candidate_pairings = list(itertools.product(p1_bonds, p2_bonds))

    # check if any of the pairs have exactly the same location, use that as a hidden signal
    # it is possible at this stage to use predetermine the distances
    # and trick it to use the ones that have exactly the same distances,
    # and treat that as a signal
    # now the issue here is that someone might "predetermine" one part, ia CA1 mapping to CB1 rathern than CB2
    # but if CA1 and CA2 is present, and CA2 is not matched to CB2 in a predetermined manner, than CB2 should not be deleted
    # so we have to delete only the offers where CA1 = CB2 which would not be correct to pursue
    if exact_coords_cue:
        predetermined = {
            a1: a2
            for a1, a2 in candidate_pairings
            if np.array_equal(a1.atom.position, a2.atom.position)
        }
        predetermined.update(
            zip(list(predetermined.values()), list(predetermined.keys()))
        )

        # skip atom pairings that have been predetermined for other atoms
        for n1_bond, n2_bond in candidate_pairings:
            if n1_bond in predetermined or n2 in predetermined:
                if (
                    predetermined[n1_bond] != n2_bond
                    or predetermined[n2_bond] != n1_bond
                ):
                    candidate_pairings.remove((n1_bond, n2_bond))

    # but they will be considered as a group
    larger_suptops = []
    pairing_and_suptop = {}
    for n1_bond, n2_bond in candidate_pairings:
        # fixme - ideally we would allow other typing than just the chemical element
        if n1_bond.atom.element is not n2_bond.atom.element:
            continue

        # create a copy of the sup_top to allow for different traversals
        # fixme: note that you could just send bonds, and that would have both parent etc with a bit of work
        larger_suptop = _overlay(
            n1_bond.atom,
            n2_bond.atom,
            parent_n1=n1,
            parent_n2=n2,
            bond_types=(n1_bond.type, n2_bond.type),
            suptop=suptop,
            ignore_coords=ignore_coords,
            use_element_type=use_element_type,
            exact_coords_cue=exact_coords_cue,
            weights=weights,
        )

        if larger_suptop is not None:
            larger_suptops.append(larger_suptop)
            pairing_and_suptop[(n1_bond, n2_bond)] = larger_suptop

    # todo
    # check for "predetermined" atoms. Ie if they have the same coordinates,
    # then that's the path to take, rather than a competing path??

    # nothing further grown out of this suptop, so it is final
    if not larger_suptops:
        return suptop

    # fixme: compare every two pairs of returned suptops, if they are compatible, join them
    # fixme - note that we are repeating this partly below
    # it also removes subgraph suptops
    # all_solutions = merge_compatible_suptops(larger_suptops)
    all_solutions = merge_compatible_suptops_faster(
        pairing_and_suptop, min(len(p1_bonds), len(p2_bonds))
    )

    # if you couldn't merge any solutions, return the largest one
    if not all_solutions:
        all_solutions = list(pairing_and_suptop.values())

    # sort in the descending order
    all_solutions.sort(key=lambda st: len(st), reverse=True)
    for sol1, sol2 in itertools.combinations(all_solutions, r=2):
        if sol1.eq(sol2):
            logger.debug(f"Removing duplicate {sol1.matched_pairs}")
            if sol2 in all_solutions:
                all_solutions.remove(sol2)

    best_suptop = extract_best_suptop(all_solutions, ignore_coords, weights=weights)
    return best_suptop


def merge_compatible_suptops_faster(pairing_suptop: dict, min_bonds: int):
    """

    :param pairing_suptop:
    :param min_bonds: if the End molecule at this point has only two bonds, they can be mapped to two other bonds
        in the start molecule.
    :return:
    """

    if len(pairing_suptop) == 1:
        return [pairing_suptop.popitem()[1]]

    # any to any
    all_pairings = list(itertools.combinations(pairing_suptop.keys(), r=min_bonds))

    if min_bonds == 3:
        all_pairings += list(itertools.combinations(pairing_suptop.keys(), r=2))
    selected_pairings = all_pairings

    # selected_pairings = []
    # for pairings in all_pairings:
    #     n = set()
    #     for pairing in pairings:
    #         n.add(pairing[0])
    #         n.add(pairing[1])
    #     #
    #     if 2 * len(pairings) == len(n):
    #         selected_pairings.append(pairings)

    # start with all the suptops as starting points
    # this is because it might be impossible to merge
    # any of the paths
    # in which case the default paths will be the best
    built_topologies = list(pairing_suptop.values())

    # attempt to combine the different traversals
    for mapping in selected_pairings:
        # mapping the different bonds to different bonds

        # check if the suptops are consistent with each other
        if not are_consistent_topologies([pairing_suptop[key] for key in mapping]):
            continue

        # merge them!
        large_suptop = copy.copy(pairing_suptop[mapping[0]])
        for next_map in mapping[1:]:
            next_suptop = pairing_suptop[next_map]

            # add both the pairs and the bonds that are not present in st1
            large_suptop.merge(next_suptop)

        built_topologies.append(large_suptop)

    return built_topologies


def extract_best_suptop(suptops, ignore_coords, weights, get_list=False):
    """
    Assumes that any merging possible already took place.
    We now have a set of solutions and have to select the best ones.

    :param suptops:
    :param ignore_coords:
    :return:
    """

    # fixme - ignore coords currently does not work
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
    def item_or_list(suptops):
        if get_list:
            return suptops
        else:
            return suptops[0]

    if len(suptops) == 0:
        warnings.warn("Cannot decide on the best mapping without any suptops...")
        return None

    elif len(suptops) == 1:
        return item_or_list(suptops)

    # candidates = copy.copy(suptops)

    # sort from largest to smallest
    suptops.sort(key=lambda st: len(st), reverse=True)

    if ignore_coords:
        return item_or_list(suptops)

    # when length is the same, take the smaller RMSD
    # most likely this is about hydrogens
    different_length_suptops = []
    for key, same_length_suptops in itertools.groupby(suptops, key=lambda st: len(st)):
        # order by RMSD
        sorted_by_rmsd = sorted(
            same_length_suptops, key=lambda st: st.align_ligands_using_mcs()
        )
        # these have the same lengths and the same RMSD, so they must be mirrors
        for suptop in sorted_by_rmsd[1:]:
            if suptop.is_mirror_of(sorted_by_rmsd[0]):
                sorted_by_rmsd[0].add_mirror_suptop(suptop)
            else:
                # add it as a different solution
                different_length_suptops.append(suptop)
        different_length_suptops.append(sorted_by_rmsd[0])

    # sort using weights
    # score = mcs_score * weight - rmsd * weight ;
    def score(st):
        # inverse for 0 to be optimal
        mcs_score = (1 - st.mcs_score()) * weights[0]

        # rmsd 0 is best as well
        rmsd_score = st.align_ligands_using_mcs() * weights[1]

        return (mcs_score + rmsd_score) / len(weights)

    different_length_suptops.sort(key=score)
    # if they have a different length, there must be a reason why it is better.
    # todo

    return item_or_list(different_length_suptops)


def are_consistent_topologies(suptops: list["SuperimposedTopology"]):
    # each to each topology has to be check if they are all consistent
    for p1, p2 in itertools.combinations(suptops, r=2):
        if not p1.is_consistent_with(p2):
            return False

    return True
