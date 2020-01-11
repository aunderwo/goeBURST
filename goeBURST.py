#!/usr/bin/env python
import os
import sys
import csv
from functools import cmp_to_key
import networkx as nx
import networkx.drawing.nx_pydot as pydot
import matplotlib.pyplot as plt

#Class UF from 	https://www.ics.uci.edu/~eppstein/PADS/UnionFind.py
class UF:
    """An implementation of union find data structure.
    It uses weighted quick union by rank with path compression.
    """

    def __init__(self, N):
        """Initialize an empty union find object with N items.
        Args:
            N: Number of items in the union find object.
        """

        self._id = list(range(N))
        self._count = N
        self._rank = [0] * N

    def find(self, p):
        """Find the set identifier for the item p."""

        id = self._id
        while p != id[p]:
            id[p] = id[id[p]]   # Path compression using halving.
            p = id[p]
        return p

    def count(self):
        """Return the number of items."""

        return self._count

    def connected(self, p, q):
        """Check if the items p and q are on the same set or not."""

        return self.find(p) == self.find(q)

    def union(self, p, q):
        """Combine sets containing p and q into a single set."""

        id = self._id
        rank = self._rank

        i = self.find(p)
        j = self.find(q)
        if i == j:
            return

        self._count -= 1
        if rank[i] < rank[j]:
            id[i] = j
        elif rank[i] > rank[j]:
            id[j] = i
        else:
            id[j] = i
            rank[i] += 1

    def __str__(self):
        """String representation of the union find object."""
        return " ".join([str(x) for x in self._id])

    def __repr__(self):
        """Representation of the union find object."""
        return "UF(" + str(self) + ")"


def load_profiles(profiles_in):
    sts = []
    profiles=[]
    with open(profiles_in) as tsvfile:
        reader = csv.reader(tsvfile, delimiter='\t')
        next(reader) # skip header
        for row in reader:
            sts.append(int(row[0]))
            profiles.append(row[1:])

    #TODO Return Unique profiles only 
    return sts, profiles

def hamm_vect(v1,v2):
    ndif=sum(1 for i, j in zip(v1, v2) if i != j)
    return ndif

def edge_comp(e1,e2):
    u,v=e1
    x,y=e2
    leveluv = hamm_vect(profiles[u],profiles[v])  
    levelxy = hamm_vect(profiles[x],profiles[y])

    if leveluv != levelxy:
        return leveluv - levelxy

    for k in range (maxlen):				
        maxuv = max(lvs[u][k], lvs[v][k])	
        maxxy = max(lvs[x][k], lvs[y][k])
        if maxuv != maxxy:
            return maxxy - maxuv

        minuv = min(lvs[u][k], lvs[v][k])
        minxy = min(lvs[x][k], lvs[y][k])

        if minuv != minxy:
            return minxy - minuv 

        maxuv = max(u,v) 
        maxxy = max(x,y)

        if maxuv != maxxy:
            return maxxy - maxuv

        minuv = min(u,v)
        minxy = min(x,y)

        if minuv != minxy:
            return minxy - minuv

def kruskal():
    edges=[] 
    nprof=len(profiles)
    for i in range(nprof):
        for j in range(i +1, nprof):
            edges.append([i,j])
    edges.sort(key=cmp_to_key(edge_comp)) 

    # var uf = new UnionFind(n)
    uf = UF(nprof)

    tree = []
    i=0
    while i<len(edges) and len(tree)<nprof-1:
         
        if uf.find(edges[i][0]) != uf.find(edges[i][1]): 
            tree.append(edges[i])		
            uf.union(edges[i][0], edges[i][1])
        
        i+=1
    
    return tree

def calc_lvs(profiles):
    maxlen= len(profiles[1])
    nprof=len(profiles)
    hamming_distances = {}
    lvs=[ [0]*maxlen for i in range(nprof)]

    for i in range(nprof):
        hamming_distances[i] = {}
        for j in range(i+1,nprof):
            diff=hamm_vect(profiles[i],profiles[j])
            hamming_distances[i][j] = diff
            if j not in hamming_distances:
                hamming_distances[j] = {}
            hamming_distances[j][i] = diff
            lvs[i][diff-1]+=1
            lvs[j][diff-1]+=1

    return lvs,maxlen,hamming_distances

def get_edge_weight(graph, sts, st1, st2):
    st1_index = sts.index(st1)
    st2_index = sts.index(st2)
    weight = graph.get_edge_data(st1_index, st2_index)['weight']
 
def main():
    try:
        profiles_in = sys.argv[1]
    except IndexError:
        print("Please supply a text file with profiles")
    if len(sys.argv) > 2:
        threshold = int(sys.argv[2])
    else:
        threshold = 1
    global profiles
    sts, profiles = load_profiles(profiles_in)

    global lvs, maxlen
    lvs,maxlen,hamming_distances=calc_lvs(profiles)
    tree=kruskal()

    labelled_tree_with_distances = [ (sts[edge[0]], sts[edge[1]], hamming_distances[edge[0]][edge[1]]) for edge in tree]
    # make graph
    G=nx.Graph()
    for edge in labelled_tree_with_distances:
        G.add_edge(edge[0],edge[1], weight=edge[2])
    
    # split Graph 
    splitG= nx.Graph()
    # add all nodes first since some maybe singletons
    splitG.add_nodes_from(G.nodes)
    for edge in G.edges():
        weight = (G.get_edge_data(edge[0], edge[1]))['weight']
        if weight <= threshold:
            splitG.add_edge(edge[0], edge[1], weight=weight)
    # plot figure
    plt.figure(figsize=(10,10))
    nx.draw_kamada_kawai(G,with_labels=True, font_weight='bold', font_size=8)
    plt.savefig("goeBURST.png")
    plt.clf()
    nx.draw_spring(splitG,with_labels=True, font_weight='bold', font_size=8)
    plt.savefig(f"goeBURST_{threshold}LV_threshold.png")
    # make dot file
    pydot.write_dot(G,'goeBURST.dot')

    # investigate subgraphs
    most_central_nodes = []
    for i, connected in enumerate(nx.connected_components(splitG)):
        sg = splitG.subgraph(connected)
        # print("subgraph {} has {} nodes".format(i, sg.number_of_nodes()))
        # print("\tNodes:", sg.nodes(data=True))
        # print("\tEdges:", sg.edges())
        sg_centrality = nx.degree_centrality(sg)
        sorted_centrality = sorted(sg_centrality.items(), key=lambda x: x[1], reverse=True)
        # print(sorted_centrality)
        max_degree_centrality = sorted_centrality[0][1]
        nodes_with_max_degree_centrality = [node  for node, degree_centrality in sorted_centrality if degree_centrality == max_degree_centrality]
        print(nodes_with_max_degree_centrality)
        most_central_nodes.extend(nodes_with_max_degree_centrality)
    print(sorted(most_central_nodes))


if __name__ == "__main__":
    main()
