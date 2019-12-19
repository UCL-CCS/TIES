# -*- coding: utf-8 -*-
import networkx as nx
import matplotlib.pyplot as plt
import hashlib

class Node:
    def __init__(self, label, value):
        self.label = label
        self.value = value

    def __hash__(self):
        m = hashlib.md5()
        m.update(str(self.label).encode('utf-8'))
        m.update(str(self.value).encode('utf-8'))
        return int(m.hexdigest(), 16)

    def __str__(self):
        return str(self.value)

    def __eq__(self, other):
        if self.value == other.value:
            return True

        return False

G1 = nx.Graph()
n1 = Node('Arbitrary1', 1)
G1.add_node(n1)
n2 = Node('Arbitrary2', 2)
G1.add_node(n2)
n3 = Node('Arbitrary3', 3)
G1.add_node(n3)
n4 = Node('Arbitrary4', 1)
G1.add_node(n4)
G1.add_edge(n1, n2)
G1.add_edge(n2, n3)
G1.add_edge(n3, n1)
G1.add_edge(n4, n3)

G2 = nx.Graph()
n1 = Node('Arbitrary5', 1)
G2.add_node(n1)
n2 = Node('Arbitrary6', 2)
G2.add_node(n2)
n3 = Node('Arbitrary7', 3)
G2.add_node(n3)
n4 = Node('Arbitrary8', 1)
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
plt.title('nx.compose G1 G2')
GC = nx.compose(G1, G2)
nx.draw(GC, with_labels=True, font_weight='bold')
plt.show()