import matplotlib.pyplot as plt
import numpy
import networkx
import dwave_networkx.algorithms
import dimod
import tabulate


class LigandMap():
    """
    Work on a list of morphs and use their information to generate a each to each map.
    This class then uses the information for
     * clustering,
     * path finding (traveling salesman, minimum spanning tree)
     * visualisation, etc.
    """

    def __init__(self, ligands, morphs):
        self.morphs = morphs
        self.ligands = ligands
        # similarity map
        self.map = None
        self.map_weights = None
        self.graph = None

    def generate_map(self):
        """
        Use the underlying morphs to extract the each to each cases.
        """
        # a simple 2D map of the ligands
        self.map = [list(range(len(self.ligands))) for l1 in range(len(self.ligands))]

        # weights based on the size of the superimposed topology
        self.map_weights = numpy.zeros([len(self.ligands), len(self.ligands)])
        for morph in self.morphs:
            self.map[morph.ligA.index][morph.ligZ.index] = morph
            self.map[morph.ligZ.index][morph.ligA.index] = morph

            matched_left, matched_right, disappearing_atoms, appearing_atoms = morph.overlap_fractions()
            # use the average number of matched fractions in both ligands
            weight = 1 - (matched_left + matched_right) / 2.0
            self.map_weights[morph.ligA.index][morph.ligZ.index] = weight
            self.map_weights[morph.ligZ.index][morph.ligA.index] = weight

            # update also the morph
            morph.set_distance(weight)

    def generate_graph(self):
        # create the nodes first
        graph = networkx.Graph()
        for ligand in self.ligands:
            graph.add_node(ligand)
        for morph in self.morphs:
            graph.add_edge(morph.ligA, morph.ligZ, weight=self.map_weights[morph.ligA.index][morph.ligZ.index])

        self.graph = graph

    def load_map(self, filename):
        self.map = numpy.loadtxt(filename)

    def traveling_salesmen(self):
        print('Traveling Salesmen (QUBO approximation): ')
        ts = dwave_networkx.traveling_salesperson(self.graph, dimod.ExactSolver())

        distances = [self.map[ligA.index][ligZ.index].distance for ligA, ligZ in zip(ts, ts[1:])]
        # one extra between the first and last
        distances.append(self.map[ts[0].index][ts[-1].index].distance)
        print(f'Sum: {sum(distances):.2f}')
        print(f'Average: {numpy.mean(distances):.2f}')

        print(ts)
        morphs = [self.map[ligA.index][ligZ.index] for ligA, ligZ in zip(ts, ts[1:])]
        return morphs

    def kruskal(self):
        print('Minimum spanning trees (Kruskal): ')
        mst = networkx.minimum_spanning_tree(self.graph)

        #
        pos = networkx.spring_layout(mst)
        networkx.draw(mst, pos=pos, with_labels=True, font_weight='bold')
        # draw the weights,
        edge_weights = {(u, v,): f"{d['weight']:.2f}" for u, v, d in mst.edges(data=True)}
        networkx.draw_networkx_edge_labels(mst, pos, edge_labels=edge_weights)
        plt.savefig(self.ligands[0].config.workdir / 'ties_map.png', dpi=300)

        # sum the distances
        mst_dsts = [item[1]['weight'] for item in mst.edges.items()]
        print(f'MST Average: {numpy.mean(mst_dsts):.2f}')
        print(f'MST Max: {numpy.max(mst_dsts):.2f}')
        print(f'Worst cases: {sorted(mst_dsts)[-10:]}')
        print("MST Pairs: \n" + '\n'.join([f'{e1} {e2} dst {w["weight"]}' for (e1, e2), w in mst.edges.items()]))

        # look up the selected MST morphs
        chosen_transformations = []
        for (ligA, ligZ), _ in mst.edges.items():
            for morph in self.morphs:
                if morph.ligA == ligA and morph.ligZ == ligZ:
                    chosen_transformations.append(morph)
                    break
        return chosen_transformations

    def print_map(self):
        # combine the maps
        print('Ligand Map')
        tabulate_map = []
        for ri, row in enumerate(self.map):
            line = []
            for ci, col in enumerate(row):
                if ci >= ri:
                    break
                line.append(repr(self.map[ri][ci]))
            if not line:
                continue
            tabulate_map.append(line)
        print(tabulate.tabulate(tabulate_map, tablefmt="grid"))
        print('LOMAP weights/similarities')
        numpy.set_printoptions(precision=2)
        print(self.map_weights)

        # save the map
        print('Saving the map as a 2D array')
        numpy.savetxt(self.ligands[0].config.workdir / 'map_weights.dat', self.map_weights, fmt='%10.5f')