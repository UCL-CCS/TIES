import networkx as nx
import MDAnalysis as mda
from os import path
import hashlib
import matplotlib.pyplot as plt
import numpy as np

"""
Use the networkx .compose G1 G2 which finds the overlap
nx.compose(G1,G2)           - combine graphs identifying nodes common to both
Composition is the simple union of the node sets and edge sets. 
The node sets of G and H do not need to be disjoint.
https://networkx.github.io/documentation/latest/reference/algorithms/generated/networkx.algorithms.operators.binary.compose.html?highlight=compose#networkx.algorithms.operators.binary.compose

nx.compose_all(list_of_graphs)
"""

# @dataclass
class AtomNode:
    def __init__(self, atomId, atomName, resName, resId, charge, atom_colloq):
        self.atomId = atomId
        self.atomName = atomName
        self.resName = resName
        self.resId = resId
        self.charge = charge
        self.atom_colloq = atom_colloq

    def __hash__(self):
        m = hashlib.md5()
        m.update(str(self.atomId).encode('utf-8'))
        m.update(str(self.charge).encode('utf-8'))
        return int(m.hexdigest(), 16)

    def __str__(self):
        return self.atom_colloq

    def __eq__(self, other):
        print('true')
        rtol = 0.05  # 5 % tolerance
        if np.isclose(self.charge, other.charge, rtol):
            return True

        return False


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
    atoms = [AtomNode(atomID, atomName, resName, resID, charge, atom_colloq) for
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


prefix = "/home/dresio/ucl/dataset/agastya_extracted/tyk2/l11-l14"
l11 = mda.Universe(path.join(prefix, 'init_l11.pdb'))
l14 = mda.Universe(path.join(prefix, 'final_l14.pdb'))

# read the corresponding charge values for the l14
l14_atoms, l14_bonds = get_charges(path.join(prefix, 'final_l14.ac'))
l11_atoms, l11_bonds = get_charges(path.join(prefix, 'init_l11.ac'))

# now that you have your atoms and their charges, enter the data into your networkx
# and use network overlay nx.compose(G1,G2)

GL1 = nx.Graph()
# create the nodes
l14_nodes = {}
for atomNode in l14_atoms:
    l14_nodes[atomNode.atomId] = atomNode
    GL1.add_node(atomNode)
# verify the number of added edges
assert len(GL1.nodes()) == len(l14_atoms)

# add the edges
for nfrom, nto in l14_bonds:
    GL1.add_edge(l14_nodes[nfrom], l14_nodes[nto])
assert len(GL1.edges()) == len(l14_bonds)

# the second graph
GL2 = nx.Graph()
l11_nodes = {}
for atomNode in l11_atoms:
    l11_nodes[atomNode.atomId] = atomNode
    GL2.add_node(atomNode)
# verify the number of added edges
assert len(GL2.nodes()) == len(l11_atoms)

# add the edges
for nfrom, nto in l11_bonds:
    GL2.add_edge(l11_nodes[nfrom], l11_nodes[nto])
assert len(GL2.edges()) == len(l11_bonds)

# import pylustrator
# pylustrator.start()

plt.subplot(131)
nx.draw(GL1, with_labels=True, font_weight='bold')
plt.subplot(132)
nx.draw(GL2, with_labels=True, font_weight='bold')
plt.subplot(133)
nx.draw(nx.compose(GL1, GL2), with_labels=True, font_weight='bold')
plt.show()
print('hi')
