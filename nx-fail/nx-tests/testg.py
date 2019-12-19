import networkx as nx
import matplotlib.pyplot as plt

G1 = nx.Graph()
[G1.add_node(x) for x in range(6)]
G1.add_edge(0, 1)
G1.add_edge(1, 2)
G1.add_edge(2, 3)
G1.add_edge(3, 4)
G1.add_edge(4, 0)
G1.add_edge(4, 5)

G2 = nx.Graph()
G2.add_node(1, vdw=0.7)
[G2.add_node(x) for x in range(2, 7)]

G2.add_edge(1, 2)
G2.add_edge(2, 3)
G2.add_edge(3, 4)
G2.add_edge(4, 5)
G2.add_edge(5, 1)
G2.add_edge(5, 6)

GC = nx.compose(G1, G2)

plt.subplot(131)
nx.draw(G1, with_labels=True, font_weight='bold')
plt.subplot(132)
nx.draw(G2, with_labels=True, font_weight='bold')
plt.subplot(133)
nx.draw(GC, with_labels=True, font_weight='bold')
plt.show()