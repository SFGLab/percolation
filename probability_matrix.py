import networkx as nx
import sys
import matplotlib.pyplot as plt
import random as rd
import time
import re


def dis_graph(graph):
    # print the adjacency list
    for line in nx.generate_adjlist(graph):
        print(line)

    # write edgelist to grid.edgelist
    nx.write_edgelist(graph, path="grid.edgelist", delimiter=" ")

    # read edgelist from grid.edgelist
    H = nx.read_edgelist(path="grid.edgelist", delimiter=" ")

    nx.draw(H)
    plt.show()
    nx.write_edgelist(graph, "connected_componets.csv", delimiter=",")


def sorting(l):
    """

    :type l: object
    """
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key=alphanum_key)


def creat_graph(nodes, edges):
    graph = nx.Graph()
    for i in nodes:
        graph.add_node(i)

    for edge in edges:
        a = edge.split(",")[1]
        b = edge.split(",")[2]
        graph.add_edge(a, b)

    lcc = max(nx.connected_component_subgraphs(graph), key=len)
    dis_graph(lcc)


def load_data():
    nodes = []
    file_whole_genome = open("GM12878.CTCF.clusters_PET4+_nodes.csv", "r")

    for k, link in enumerate(file_whole_genome):
        if k == 0:
            continue
        link = link.split("\n")[0]
        nodes.append(link)

    file_whole_genome.close()

    nodes = list(dict.fromkeys(nodes))
    x = sorting(nodes)
    nodes = x
    nodes = list(dict.fromkeys(nodes))

    edges = []
    file_whole_genome = open("GM12878.CTCF.clusters_PET4+_edges.csv", "r")
    for k, link in enumerate(file_whole_genome):
        if k == 0:
            continue
        # if k > 500:
        #     break
        link = link.split("\n")[0]
        edges.append(link)

    file_whole_genome.close()

    print(len(nodes))
    print(len(edges))

    return nodes, edges


if __name__ == '__main__':
    nodes, edges = load_data()
    creat_graph(nodes, edges)
