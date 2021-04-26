import networkx as nx
import sys
import matplotlib.pyplot as plt
import random as rd
import time
import re
from multiprocessing import Pool


def draw_plot(chr_name, x_axis, y_axis):
    print("Drawing graph")
    print('x_axis values : ')
    print(x_axis)
    print('y_axis values : ')
    print(y_axis)
    fig, ax = plt.subplots()
    ax.plot(x_axis, y_axis)
    ax.set(xlabel='time (s)', ylabel='C1/n',
           title=chr_name + ' CTCF 2+ Random Edges Plot')
    ax.grid()
    fig.savefig(chr_name + " CTCF 2+ Random Edges Plot.png")
    plt.show()


def sorting(l):
    """

    :type l: object
    """
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key=alphanum_key)


def creat_graph(chr_mane, n, e):
    print("No of Nodes after Graph Creation in " + chr_mane + " : " + str(len(n)))
    print("No of Edges after Graph Creation in " + chr_mane + " : " + str(len(e)))
    x = []
    y = []
    k = 0.5

    graph = nx.Graph()
    for i in n:
        graph.add_node(i)

    while 1:
        print(len(e))
        edge = rd.choice(e)
        e.remove(edge)
        a = edge.split(",")[0]
        b = edge.split(",")[1]
        graph.add_edge(a, b)
        lcc = max(nx.connected_component_subgraphs(graph), key=len)
        y.append(len(lcc.nodes) / len(n))
        x.append(k)
        k += 0.5
        if len(e) == 0:
            break

    return x, y


def erdos_renyi(graph, p, n):
    k = 0
    for i in graph.nodes:
        print(i)
        for j in graph.nodes:
            if i < j:
                rand = rd.random()
                if rand <= p: # and check_edge(i, j)
                    graph.add_edge(i, j)
                    # time.sleep(10)
                    # dis_graph(graph)
                    lcc = max(nx.connected_component_subgraphs(graph), key=len)
                    y_axis.append(len(lcc.nodes) / n)
                    x_axis.append(k)
                    k += 1
                    # print(x_axis)
                else:
                    continue


# def load_data():
#     nodes = []
#     file_whole_genome = open("data/GM12878_RNAPII_PET2+_nodes.csv", "r")
#
#     for k, link in enumerate(file_whole_genome):
#         if k == 0:
#             continue
#         link = link.split("\n")[0]
#         nodes.append(link)
#
#     file_whole_genome.close()
#
#     nodes = list(dict.fromkeys(nodes))
#     x = sorting(nodes)
#     nodes = x
#     nodes = list(dict.fromkeys(nodes))
#
#     edges = []
#     file_whole_genome = open("data/GM12878_RNAPII_PET2+_edges.csv", "r")
#     for k, link in enumerate(file_whole_genome):
#         if k == 0:
#             continue
#         # if k > 500:
#         #     break
#         link = link.split("\n")[0]
#         edges.append(link)
#
#     file_whole_genome.close()
#
#     print(len(nodes))
#     print(len(edges))
#
#     return nodes, edges

def load_data(edge):
    edg = []
    nod = []

    for k, link in enumerate(edge):
        if k == 0:
            continue
        link = link.split("\n")[0]
        edg.append(link)
        nod.append(link.split(",")[0])
        nod.append(link.split(",")[1])

    nod = list(dict.fromkeys(nod))
    x = sorting(nod)
    nod = x
    nod = list(dict.fromkeys(nod))

    # print(len(nod))
    # print(len(edg))

    return nod, edg


def extract_data(chr_n):
    file_whole_genome = open("data/GM12878.CTCF.clusters_all_edges.csv", "r")
    edges = []
    for line in file_whole_genome:
        line = line.split("\n")[0]
        chr_name = line.split(",")[0].split(":")[0]
        if chr_name == chr_n:
            # print(line.split(",")[0] + "," + line.split(",")[1] + "," + line.split(",")[2])
            edges.append(line.split(",")[0] + "," + line.split(",")[1] + "," + line.split(",")[2])

    file_whole_genome.close()
    # print(len(edges))
    return edges


def add_linear_edge(edges):
    ccd_edge = []
    linear_edge = []
    for k, line2 in enumerate(edges):
        line2 = line2.split("\n")[0]
        ccd_edge.append(line2.split(",")[0])
        ccd_edge.append(line2.split(",")[1])
        # print(line2)
        linear_edge.append(line2)

    ccd_edge = list(dict.fromkeys(ccd_edge))
    x = sorting(ccd_edge)
    ccd_edge = x
    # print(len(ccd_edge))
    ccd_edge = list(dict.fromkeys(ccd_edge))
    # print(len(ccd_edge))

    for i in range(len(ccd_edge) - 1):
        # print(ccd_edge[i])
        # print(line.split("\t")[0] + "," + ccd_edge[i] + "," + ccd_edge[i+1] + ",1,Linear_edge,linear_edge")
        linear_edge.append(ccd_edge[i] + "," + ccd_edge[i + 1] + ",1")
    edge = []
    edge.append("anchor_id_A,anchor_id_B,link_score\n")
    for line_write in linear_edge:
        # print(line_write)
        edge.append(line_write + "\n")
    return edge


def filter_edges(edges):
    filter_edge = []
    for edge in edges:
        node_a = edge.split(",")[0]
        node_start_a = int(node_a.split("-")[0].split(":")[1])
        node_end_a = int(node_a.split("-")[1])
        dist_a = node_end_a - node_start_a
        node_b = edge.split(",")[1]
        node_start_b = int(node_b.split("-")[0].split(":")[1])
        node_end_b = int(node_b.split("-")[1])
        dist_b = node_end_b - node_start_b
        if dist_a > 100 and dist_b > 100:
            # print(edge)
            link_wt = int(edge.split(",")[2])
            if link_wt > 2:
                filter_edge.append(edge)
    # print(len(filter_edge))
    return filter_edge


def multi_driver(chr_mane):
    # chr_mane = "chr17"  # input("Enter chromosome for which plot should be shown : ")
    edges = extract_data(chr_mane)
    print("No of Edges in " + chr_mane + " : " + str(len(edges)))
    ed = filter_edges(edges)
    print("No of Edges after filtering in " + chr_mane + " : " + str(len(ed)))
    edge = add_linear_edge(ed)
    print("No of Edges after adding liner edges in " + chr_mane + " : " + str(len(edge)))
    n, e = load_data(edge)
    print("No of Nodes after Loading Data in " + chr_mane + " : " + str(len(n)))
    print("No of Edges after Loading Data in " + chr_mane + " : " + str(len(e)))
    x_axis, y_axis = creat_graph(chr_mane, n, e)
    # return x_axis, y_axis
    draw_plot(chr_mane, x_axis, y_axis)


if __name__ == '__main__':
    chr_mat = ["chr1", "chr2", "chr3"]
    p = Pool(5)
    print(p.map(multi_driver, chr_mat))

    # x = [[], [], []]
    # for i, chr_mane in enumerate(chr_mat):
    #     # chr_mane = "chr17"  # input("Enter chromosome for which plot should be shown : ")
    #     edges = extract_data(chr_mane)
    #     print("No of Edges in " + chr_mane + " : " + str(len(edges)))
    #     ed = filter_edges(edges)
    #     print("No of Edges after filtering in " + chr_mane + " : " + str(len(ed)))
    #     edge = add_linear_edge(ed)
    #     print("No of Edges after adding liner edges in " + chr_mane + " : " + str(len(edge)))
    #     n, e = load_data(edge)
    #     print("No of Nodes after Loading Data in " + chr_mane + " : " + str(len(n)))
    #     print("No of Edges after Loading Data in " + chr_mane + " : " + str(len(e)))
    #     x[i], y_axis = creat_graph(n, e)
    # draw_plot(chr_mane, x, y_axis)

if __name__ == '__main__':
    for k in range(1):
        del x_axis[:]
        del y_axis[:]
        print("Takking Values")
        n = 100  # int(input("Enter value of n : "))
        p = 0.5  # float(input("Enter value of p : "))
        print(n, p)
        print("Creat Matrix")
        creat_matrix(n)
        # edge_matrix = [[0, 1, 0, 1, 1], [1, 0, 1, 0, 0], [0, 1, 0, 1, 0], [1, 0, 1, 0, 0], [1, 0, 0, 0, 0]]
        # print(edge_matrix)
        graph = nx.Graph()
        for i in range(n):
            graph.add_node(i)
        # dis_graph(graph)
        print("Creat Graph")
        erdos_renyi(graph, p, n)
        print("Draw Plot")
        draw_plot()
        # print(type(graph))
