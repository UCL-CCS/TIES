# -*- coding: utf-8 -*-
import networkx as nx
import matplotlib.pyplot as plt
import hashlib

"""
Case: what if there are multiple matches on the two graphs?

In G1 there is one node with charge 0.7, and in G2 there are two such nodes. 
From these two nodes, one is disconnected from the graph - and most likely for that reason, 
it is never chosen as a match for the 0.7 node in G1. In other words, 
the edges and over topology match appears to be taken into account.   

Notes:
Having the merged version we can traverse through it. There will be nodes from both. This means that as we
 traverse through it. We will run into nodes which are missing in the other place. So this way we know which atoms
 can go missing. 
 
Fixme:
 - I should put more effort into making it a real life example in order to get some real-life cases.  

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
        print('match', self, other)
        if self.charge == other.charge:
            return True

        print('no match', self, other)
        return False

    def __lt__(self, other):
        print('lt')
    def __le__(self, other):
        print('le')
    def __ne__(self, other):
        print('ne')
    def __gt__(self, other):
        print('gt')
    def __ge__(self, other):
        print('ge')

G1 = nx.Graph()
n1 = Node('L1', 0.7)
G1.add_node(n1)
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
n1 = Node('L1', 0.7)
G2.add_node(n1)
n2 = Node('L2', -0.3)
G2.add_node(n2)
n3 = Node('L5', -0.1)
G2.add_node(n3)
n4 = Node('L1', 0.9)
G2.add_node(n4)
G2.add_edge(n1, n2)
G2.add_edge(n2, n3)
G2.add_edge(n3, n1)
G2.add_edge(n4, n3)
# extend with the same charge atom,
n = Node('L1 Match', 0.7)
G2.add_node(n)
G2.add_edge(n4, n)

plt.subplot(141)
nx.draw(G1, with_labels=True, font_weight='bold')
plt.subplot(142)
nx.draw(G2, with_labels=True, font_weight='bold')
plt.subplot(143)
plt.title('nx.compose G1 G2')
GC = nx.compose(G1, G2)
nx.draw(GC, with_labels=True, font_weight='bold')
plt.subplot(144)
plt.title('nx.compose G2 G1')
GC = nx.compose(G2, G1)
nx.draw(GC, with_labels=True, font_weight='bold')
plt.show()