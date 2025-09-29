import argparse
import copy
from argparse import ArgumentParser
from collections import OrderedDict, defaultdict
from itertools import islice
from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from networkx.exception import NetworkXNoPath


class Mapper:

    dst_prop = "weight"
    active = "0"

    dst_cutoff = 0.40
    do_not_make_me_repeat_myself_limit = 0.30
    soft_limit = 0.20

    # global data
    rev_sol_diG = nx.DiGraph()
    sol_diG = nx.DiGraph()
    sol_G = nx.Graph()
    selected_paths = []

    def __init__(self, mcs_data):
        self.full_G = nx.read_graphml(mcs_data)

    # extract plausible shortest paths
    def get_shortest_paths(self, num_shortest_paths=10, edges_threshold=0.4):

        G = mapper.remove_edges(mapper.full_G, threshold=edges_threshold)

        sols = defaultdict(list)
        for node in G.nodes:
            if node == self.active:
                continue

            try:
                for path in islice(
                    nx.shortest_simple_paths(
                        G, node, self.active, weight=self.dst_prop
                    ),
                    num_shortest_paths,
                ):
                    sols[node].append(tuple(path))
                    # print("-".join(path), [f"{d:.2f}" for d in get_dsts(path)])
            except NetworkXNoPath:
                print(f"Could not find the shortest path for {node}")

        return sols

    # drop edges that will never be used
    # remove edges that do not meet the threshold
    def remove_edges(self, G, threshold=0.4):
        G = copy.deepcopy(G)
        for edge_str in G.edges:
            edge = G.edges[edge_str]

            # weight is dst
            if edge[self.dst_prop] > threshold:
                G.remove_edge(*edge_str)

        return G

    def bfs_printer(self, G):
        # traverse BFS
        print("Simple - for docking and bash scripting")
        for source, node in nx.bfs_edges(G, self.active):
            print(node, source, f"{G.get_edge_data(source, node)[self.dst_prop]:.2f}")

    def bfs_depth_nodes(self, G):
        # extract nodes by distance from origin
        nodes = []
        for source, node in nx.bfs_edges(G, self.active):
            if node not in nodes:
                nodes.append(node)

        return nodes

    def node_to_dst_from_origin(self, paths):
        """
        Get the
        :param paths:
        :return:
        """

        dsts = defaultdict(lambda: float("inf"))
        for path in paths:
            dsts[path[0]] = min(len(path), dsts[path[0]])

        return dsts

    def bfs_edge_printer(self):
        # traverse BFS but print all edges for each node
        seen = set()
        print("Simple - for docking and bash scripting")
        for source, _ in nx.bfs_edges(self.sol_G, self.active):
            for edge in self.sol_G.edges(source):
                fs = frozenset(edge)
                if fs in seen:
                    continue

                u, v = edge
                print(u, v, f"{self.full_G.get_edge_data(u, v)[self.dst_prop]:.2f}")

                seen.add(fs)

    def edge_printer(self, G):
        # traverse BFS
        print("Simple - for docking and bash scripting")
        for source, node in nx.bfs_edges(G, self.active):
            print(node, source, f"{G.get_edge_data(source, node)[self.dst_prop]:.2f}")

    def get_dsts(self, path):
        dsts = []
        for a, b in zip(path, path[1:]):
            # print(f"\t {G.edges[a, b]['weight']:.2f}")
            dsts.append(self.full_G.edges[a, b][self.dst_prop])

        return dsts

    def get_cycle_dsts(self, path):
        dsts = []
        cycle_path = list(path) + [path[0]]
        for a, b in zip(path, cycle_path[1:]):
            # print(f"\t {G.edges[a, b]['weight']:.2f}")
            dsts.append(self.full_G.edges[a, b][self.dst_prop])

        return dsts

    def get_dst(self, path):
        return sum(self.get_dsts(path))

    def get_dst_cycle(self, path):
        return sum(self.get_dsts(list(path) + [path[0]]))

    def mult_dst_cycle(self, path):
        """
        So 0.1*0.1*01 > 0.2*0.05 * 0.05
        So a larger number is better

        :param path:
        :return:
        """
        return np.prod(self.get_dsts(list(path) + [path[0]]))

    def norm_mult_dst_cycle(self, path):
        """
        So 0.1*0.1*01 > 0.2*0.05 * 0.05
        So a larger number is better

        :param path:
        :return:
        """
        return np.prod(self.get_dsts(list(path) + [path[0]])) * self.get_dst_cycle(path)

    def get_max_step(self, path):
        return max(self.get_dsts(path))

    def select_path(self, paths, paths_to_ignore=None, limit=0.35):
        """
        Always use path with better perturbations.

        :param paths:
        :return:
        """

        if paths_to_ignore is None:
            paths_to_ignore = []

        perturbation_quality_limits = np.linspace(0.10, limit, round(limit / 0.05 - 1))
        for pert_quality in perturbation_quality_limits:
            for path in paths:
                if path in paths_to_ignore:
                    continue

                # below soft, number of steps is the best
                if self.get_max_step(path) < pert_quality:
                    return path

        # worst case scenario
        return paths[0]

    def add_path(self, path):
        if path in self.selected_paths:
            # no need to add a path that already exists
            return

        self.selected_paths.append(path)

        for node1, node2 in zip(path, path[1:]):
            # fixme - use the dst_prop
            self.sol_G.add_edge(
                node1,
                node2,
                weight=self.full_G.get_edge_data(node1, node2)[self.dst_prop],
            )
            self.sol_diG.add_edge(
                node1,
                node2,
                weight=self.full_G.get_edge_data(node1, node2)[self.dst_prop],
            )

        # reverse
        for node1, node2 in zip(path[::-1], path[-2::-1]):
            # fixme - use the dst_prop
            self.rev_sol_diG.add_edge(
                node1,
                node2,
                weight=self.full_G.get_edge_data(node1, node2)[self.dst_prop],
            )

    def add_shortest_paths_hub(self):
        """
        HUB approach.

        Take a list of nodes that have more than "1 edge".
        This means that they have other edges going via them.

        Hub size:
         - 1 extra node: just add one edge (cumulative)
         - 2 or more nodes: add 2 more edges

        Hub distance:
         - 1 (direct): add one edge
         - 2: add two edges
         - 3: 4 edges
         - 4: 5 edges

        Hub size/distance -> edges:
         - 0/*  -> 0
         - 1/1  -> 0
         - 2/1  -> 1
         - =>3/1 -> 2

         - 1/2  -> 1
         - 1/3  -> 1

         - 2/2  -> 2


         - 2/3  -> 3 ?
         - >3/2 -> 3
         - 2/>3 -> 3

         - 3/3  -> 4
         - >3/>3 -> 5

        It will be important that these new paths will have cycles.
        """
        G = mapper.rev_sol_diG
        # start with 0
        original_G = copy.deepcopy(G)
        for source, node in nx.bfs_edges(original_G, self.active):

            # ignore if not a hub
            hub_size = original_G.out_degree[node]
            # the original (best) distance
            dst = [len(p) - 1 for p in self.selected_paths if p[0] == node][0]

            add_paths = 0
            if hub_size == 0:
                continue
            if hub_size == 1 and dst == 1:
                continue
            elif {hub_size, dst} == {2, 1}:
                add_paths = 1
            elif hub_size == 1 and dst in {2, 3}:
                add_paths = 2
            elif {hub_size, dst} == {2} or hub_size > 3 and dst == 1:
                add_paths = 2
            elif hub_size == 2 and dst == 3:
                add_paths = 3
            elif hub_size >= 3 and dst == 2:
                add_paths = 3
            elif hub_size == 2 and dst > 3:
                add_paths = 3
            elif {hub_size, dst} == {3}:
                add_paths = 4
            elif hub_size > 3 and dst > 3:
                add_paths = 5
            else:
                print("there is some other option?", hub_size, dst)

            for i in range(add_paths):
                path = self.select_path(
                    shortest_paths[node], self.selected_paths, limit=0.30
                )
                self.add_path(path)

            # print(node, source, f"{G.get_edge_data(source, node)[dst_prop]:.2f}")

    def cycle_closure_hub(self):
        """
        Cycle Closure

        For hubs, we need cycle closure. Ideally, it would not an "edge" cycle closure.

        Find all simple cycles
        """
        rev_sol_diG = mapper.rev_sol_diG
        G03 = self.remove_edges(self.full_G, threshold=0.3)
        all_cycles = nx.cycle_basis(G03)
        dst_to_0 = self.node_to_dst_from_origin(self.selected_paths)

        def get_min_node_to0(path):
            # is this MST basically?
            return min(dst_to_0[p] for p in path)

        # get for each node, the distance it is from the active

        # start with 0
        original_G = copy.deepcopy(rev_sol_diG)
        cycles_added = 0
        for source, node in nx.bfs_edges(original_G, self.active):
            hub_size = original_G.out_degree[node]
            if hub_size < 1:
                continue

            # print(node)

            cycles = [c for c in all_cycles if node in c and len(c) == 3]

            # add and order by the cycle score (path length)
            cycle_score = OrderedDict(
                sorted(
                    [(tuple(c), self.get_dst_cycle(c)) for c in cycles],
                    key=lambda x: x[1],
                )[:5]
            )

            # for cycle, score in cycle_score.items():
            #     if len(cycle) > 3:
            #         continue
            #     print(
            #         "close",
            #         [dst_to_0[n] for n in cycle],
            #         "norm mult",
            #         norm_mult_dst_cycle(cycle),
            #         "mult",
            #         mult_dst_cycle(cycle),
            #         "dst",
            #         get_dst_cycle(cycle),
            #         "all dsts",
            #         get_cycle_dsts(cycle),
            #     )

            # pick the cycle that ideally takes you closer
            # pick the one that either is closer, or the first item

            if len(cycle_score) == 0:
                print("No cycles found for", node)
                continue

            shortest_cycle, dst = sorted(
                cycle_score.items(), key=lambda x: get_min_node_to0(x[0])
            )[0]
            self.add_path(list(shortest_cycle) + [shortest_cycle[0]])
            cycles_added += 1

        return cycles_added


parser = argparse.ArgumentParser(
    description="",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)
parser.add_argument(
    "-mcs_data",
    metavar="filename",
    dest="mcs_matrix",
    type=Path,
    required=False,
    default="mcs_data.graphml",
    help="A networkx file .graphml with the MCS data. ",
)
parser.add_argument(
    "-out_network",
    metavar="filename",
    dest="output_network",
    type=Path,
    required=False,
    default="perturbation_network.graphml",
    help="An output perturbation network to be computed (networkx file). ",
)

if __name__ == "__main__":
    args = parser.parse_args()

    mapper = Mapper(args.mcs_matrix)

    shortest_paths = mapper.get_shortest_paths()

    # add the shortest paths
    for node, paths in shortest_paths.items():
        path = mapper.select_path(paths)
        mapper.add_path(path)
        # path2 = select_path(paths, [path])
        # add_path(path2)

    mapper.add_shortest_paths_hub()

    cycles_added = mapper.cycle_closure_hub()

    print("---------------------")
    print("Data available for: ", len(shortest_paths))
    print("Selected: ", len(mapper.sol_G))
    print("Final paths: ", len(mapper.selected_paths))
    print("Cycles added", cycles_added)
    # bfs_printer(rev_sol_diG)
    print("Number of unique steps/perturbations: ", len(mapper.sol_G.edges))
    # nodes = mapper.bfs_depth_nodes(mapper.sol_G)

    mapper.bfs_edge_printer()

    nx.write_graphml(mapper.sol_G, args.output_network)

    # nx.draw(mapper.sol_G, with_labels=True)
    # plt.savefig("bps.png")
    # plt.show()
