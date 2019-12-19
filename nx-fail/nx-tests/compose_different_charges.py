# -*- coding: utf-8 -*-
import networkx as nx
import matplotlib.pyplot as plt
import hashlib

"""
Cases:
 - different charges are recognized as not the same.
 - what if there are multiple matches on the two graphs?  

"""

class Node:
    def __init__(self, label, charge):
        self.label = label
        self.charge = charge

    def __hash__(self):
        m = hashlib.md5()
        m.update(str(self.label).encode('utf-8'))
        m.update(str(self.charge).encode('utf-8'))
        return int(m.hexdigest(), 16)

    def __str__(self):
        return self.label + '_' + str(self.charge)

    def __eq__(self, other):
        if self.charge == other.charge:
            return True

        return False

G1 = nx.Graph()
n1 = Node('L1', 0.7)
G1.add_node(n1)
n2 = Node('L2', -0.3)
G1.add_node(n2)
n3 = Node('L3', -0.1)
G1.add_node(n3)
n4 = Node('L1', 0.9)
G1.add_node(n4)
G1.add_edge(n1, n2)
G1.add_edge(n2, n3)
G1.add_edge(n3, n1)
G1.add_edge(n4, n3)

G2 = nx.Graph()
n1 = Node('L1', 0.6)
G2.add_node(n1)
n2 = Node('L2', -0.3)
G2.add_node(n2)
n3 = Node('L3', -0.1)
G2.add_node(n3)
n4 = Node('L1', 0.9)
G2.add_node(n4)
G2.add_edge(n1, n2)
G2.add_edge(n2, n3)
G2.add_edge(n3, n1)
G2.add_edge(n4, n3)

plt.subplot(131)
nx.draw(G1, with_labels=True, font_weight='bold')
plt.subplot(132)
nx.draw(G2, with_labels=True, font_weight='bold')
plt.subplot(133)
plt.title('nx.compose')
GC = nx.compose(G1, G2)
nx.draw(GC, with_labels=True, font_weight='bold')
plt.show()