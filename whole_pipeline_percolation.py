import argparse
import random as rd
import time
from heapq import nlargest
from multiprocessing import Pool
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import collections
from sklearn import metrics
from networkx.readwrite import json_graph
import json as js
from tqdm import tqdm
import math
import os
import py2cytoscape
import seaborn as sns
from sklearn import preprocessing
from matplotlib import pyplot as plt
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis


def critical_point(graph, chr_name, k_values, cc1_sizes, cc2_sizes, moments_2):
    indx = np.argmax(cc2_sizes)
    print(f"{chr_name}", end=",")
    print(f"{len(list(graph.nodes()))}", end=",")
    print(f"{len(list(graph.edges()))}", end=",")
    print(f"{k_values[indx]}", end=",")
    print(f"{np.sum(np.gradient(cc1_sizes), axis=0)}", end=",")
    print(f"{metrics.mutual_info_score(cc1_sizes, moments_2)}", end=",")
    print(f"{metrics.mutual_info_score(cc1_sizes, cc2_sizes)}", end=",")
    return indx


def delta_fun(graph, chr_mane, k_values, cc1_sizes, gama=(1 / 2), a=0.1):
    nodes = len(list(graph.nodes()))
    edges = len(list(graph.edges()))

    g = int(pow(nodes, gama))
    a = int(nodes * a)
    print("g and a : ", g, a)

    v1 = 0
    v2 = 0
    for k, line in enumerate(cc1_sizes):
        p = float(line)
        if p > g:
            v1 = k
            break

    for k, line in enumerate(cc1_sizes):
        p = int(line)
        if p >= a:
            v2 = k
            break

    no_of_edge = v2 - v1
    print(no_of_edge)
    # pos1 = np.where(cc1_sizes == v1)
    # pos2 = np.where(cc1_sizes == v2)
    # print(pos1, pos2)


def draw_plot(chr_name, n, k_values, cc1_sizes, cc2_sizes, moments_2, type):
    fig, ax = plt.subplots(figsize=(6, 4), dpi=600)
    ax.set_facecolor((1.0, 1.0, 1.0))
    ax.se
    ax.plot(k_values, cc1_sizes, c='C0', label='CC1/n')
    ax.set(xlabel='Proportion of Edges added', ylabel='Largest cluster size' , title=f"{type}")
    ax.plot(k_values, cc2_sizes, c='C1', label='CC2/n')
    ax.plot(k_values, moments_2, c='C2', label='m2 of |CC| (norm.)')
    ax.legend(loc='upper left')
    ax.grid()
    fig.savefig(f"plots/{chr_name}/CTCF 2+ {type} Plot.png")
    plt.show()


def draw_plot_part(start, end, chr_name, n, k_values, cc1_sizes, cc2_sizes, moments_2, type):
    k_values = k_values[start:end]
    cc1_sizes = cc1_sizes[start:end]
    cc2_sizes = cc2_sizes[start:end]
    moments_2 = moments_2[start:end]
    fig, ax = plt.subplots(figsize=(6, 4), dpi=150)
    ax.plot(k_values, cc1_sizes, c='C0', label='C1/n')
    ax.set(xlabel='no.edges', title=f"{chr_name} CTCF graph {type}")
    ax.plot(k_values, cc2_sizes, c='C1', label='C2/n')
    ax.plot(k_values, moments_2, c='C2', label='m2 of |CC| (norm.)')
    ax.legend(loc='best')
    ax.grid()
    fig.savefig(f"plots/{chr_name} CTCF 2+ {type} {start} {end} Plot.png")
    plt.show()


def draw_plot_cmp(chr_name, n, k_values_er, cc1_sizes_er, k_values_ae, cc1_sizes_ae, k_values_te, cc1_sizes_te,
                  k_values_fe, cc1_sizes_fe, k_values_bm, cc1_sizes_bm, k_values_em, cc1_sizes_em, k_values_pm,
                  cc1_sizes_pm, type):
    fig, ax = plt.subplots(figsize=(6, 4), dpi=600)
    ax.set(xlabel='Proportion of Edges added', ylabel='Largest cluster size', title=f'Comparison between {type}')
    ax.plot(k_values_er, cc1_sizes_er, c='C0', label='ER Model')
    ax.plot(k_values_ae, cc1_sizes_ae, c='C1', label='AE Model')
    ax.plot(k_values_te, cc1_sizes_te, c='C2', label='TE Model')
    ax.plot(k_values_fe, cc1_sizes_fe, c='C3', label='LF Model')
    ax.plot(k_values_bm, cc1_sizes_bm, c='C4', label='CC Model')
    ax.plot(k_values_em, cc1_sizes_em, c='C5', label='CL Model')
    ax.plot(k_values_pm, cc1_sizes_pm, c='C6', label='CA Model')
    ax.legend(loc='upper left')
    ax.grid()
    fig.savefig(f"plots/{chr_name}/CTCF_2+_Comparison_Plot_{type}.png")
    plt.show()


def draw_degree_histogram(graph, chr_name, type):
    G = graph
    chr_no = chr_name[3:]
    chr_namex = "Chromosome " + chr_no
    print(chr_name)
    degree_sequence = sorted([d for n, d in G.degree()], reverse=True)  # degree sequence
    # print "Degree sequence", degree_sequence
    degreeCount = collections.Counter(degree_sequence)
    deg, cnt = zip(*degreeCount.items())

    fig, ax = plt.subplots()
    plt.bar(deg, cnt, width=0.80, color='b')
    plt.yscale('log')
    plt.title(f"Degree Histogram {chr_namex}")
    plt.ylabel("Count")
    plt.xlabel("Degree")
    ax.set_xticks([d + 0.4 for d in deg])
    ax.set_xticklabels(deg)
    fig.savefig(f"plots/{chr_name}/CTCF 2+ Degree Histogram {type} Plot.png")
    plt.show()


def find_connected_components(k, graph, prev_cc):
    # print(len(list(graph.edges())))
    if prev_cc != nx.number_connected_components(graph):
        print(nx.number_connected_components(graph), "->", k)
    return nx.number_connected_components(graph)


def get_prob_edge(graph, a, b):
    k = graph[a][b]["kind"]
    w = graph[a][b]["weight"]
    weight = int(w)
    weight = float(weight / graph.number_of_edges())
    print(w)
    return weight


def run_percolation_er(graph, chr_name):
    graph_new = nx.Graph()
    edges_shuffled = list(graph.edges())
    rd.shuffle(edges_shuffled)

    roots = list(graph.nodes())
    uf = nx.utils.union_find.UnionFind(roots)

    data_size = len(edges_shuffled)
    k_values = np.empty(data_size)
    cc1_sizes = np.empty(data_size)
    density_cc1 = np.empty(data_size)
    cc2_sizes = np.empty(data_size)
    moments_2 = np.empty(data_size)
    # flag = np.zeros(len(sub))
    # ccd_percolation = np.zeros(len(sub))

    for k, (a, b) in enumerate(edges_shuffled):

        graph.add_edge(a, b)
        w = graph[a][b]["weight"]
        graph_new.add_edge(a, b, weight=w)
        # print(k + 1)
        file_graph_write = open(
            f"/mnt/raid/kaustavs/Desktop/Project_data/perculation/percolation/data/graph_test/{chr_name}/graph_test_er/er_ctcf_2+_" + str(
                k + 1) + ".csv", "w")
        # file_graph_write.write("Source,Target,Weight\n")
        for (v1, v2) in list(graph_new.edges()):
            na = v1[0] + ":" + str(v1[1]) + "-" + str(v1[2])
            nb = v2[0] + ":" + str(v2[1]) + "-" + str(v2[2])
            x = na + "," + nb + "," + str(graph_new[v1][v2]["weight"])
            file_graph_write.write(x + "\n")
        file_graph_write.close()
        root_a = uf[a]
        root_b = uf[b]
        if root_a != root_b:
            uf.union(a, b)
            root_ab = uf[a]
            roots.remove(root_a)
            roots.remove(root_b)
            roots.append(root_ab)
        sizes = [uf.weights[v] for v in roots]
        unique_sizes, counts = np.unique(sizes, return_counts=True)
        moment_2 = ((unique_sizes ** 2) * counts).sum()
        if len(roots) > 1:
            cc1_size, cc2_size = nlargest(2, sizes)
        else:
            cc1_size, cc2_size = uf.weights[roots[0]], 0.0
        cc1_sizes[k] = cc1_size / len(list(graph.nodes()))
        density_cc1[k] = cc1_size / len(list(graph.nodes()))
        cc2_sizes[k] = cc2_size / len(list(graph.nodes()))
        k_values[k] = k + 1
        moments_2[k] = moment_2
        # if k == 1:
        #     break

    m = k_values.max()
    for i, values in enumerate(k_values):
        k_values[i] = values / m

    m = moments_2.max()
    for i, values in enumerate(moments_2):
        moments_2[i] = values / m

    # nx.draw(graph, pos=nx.spring_layout(graph))
    return graph_new, k_values, cc1_sizes, cc2_sizes, moments_2, density_cc1  # int(np.amin(ccd_percolation)), int(np.amax(ccd_percolation))


def run_percolation_ae(graph, chr_name):
    graph_new = nx.Graph()
    edges_shuffled = list(graph.edges())
    rd.shuffle(edges_shuffled)

    roots = list(graph.nodes())
    uf = nx.utils.union_find.UnionFind(roots)

    data_size = len(edges_shuffled)
    k_values = np.empty(data_size)
    cc1_sizes = np.empty(data_size)
    density_cc1 = np.empty(data_size)
    cc2_sizes = np.empty(data_size)
    moments_2 = np.empty(data_size)
    # flag = np.zeros(len(sub))
    # ccd_percolation = np.zeros(len(sub))

    for k, (a, b) in enumerate(edges_shuffled):

        neb = list(graph.edges([a, b]))
        v1 = a, b
        neb.remove(v1)
        v2 = rd.choice(neb)
        r = rd.randint(0, 500)
        if r % 2 == 0:
            a, b = v1
        else:
            a, b = v2
        graph.add_edge(a, b)
        w = graph[a][b]["weight"]
        graph_new.add_edge(a, b, weight=w)

        file_graph_write = open(
            f"/mnt/raid/kaustavs/Desktop/Project_data/perculation/percolation/data/graph_test/{chr_name}/graph_test_ae/ae_ctcf_2+_" + str(
                k + 1) + ".csv", "w")
        # file_graph_write.write("Source,Target,Weight\n")
        for (a1, b1) in list(graph_new.edges()):
            na = a1[0] + ":" + str(a1[1]) + "-" + str(a1[2])
            nb = b1[0] + ":" + str(b1[1]) + "-" + str(b1[2])
            x = na + "," + nb + "," + str(graph_new[a1][b1]["weight"])
            file_graph_write.write(x + "\n")
        file_graph_write.close()
        root_a = uf[a]
        root_b = uf[b]
        if root_a != root_b:
            uf.union(a, b)
            root_ab = uf[a]
            roots.remove(root_a)
            roots.remove(root_b)
            roots.append(root_ab)
        sizes = [uf.weights[v] for v in roots]
        unique_sizes, counts = np.unique(sizes, return_counts=True)
        moment_2 = ((unique_sizes ** 2) * counts).sum()
        if len(roots) > 1:
            cc1_size, cc2_size = nlargest(2, sizes)
        else:
            cc1_size, cc2_size = uf.weights[roots[0]], 0.0
        cc1_sizes[k] = cc1_size / len(list(graph.nodes()))
        density_cc1[k] = cc1_size / len(list(graph.nodes()))
        cc2_sizes[k] = cc2_size / len(list(graph.nodes()))
        k_values[k] = k + 1
        moments_2[k] = moment_2

    m = k_values.max()
    for i, values in enumerate(k_values):
        k_values[i] = values / m

    m = moments_2.max()
    for i, values in enumerate(moments_2):
        moments_2[i] = values / m

    return graph_new, k_values, cc1_sizes, cc2_sizes, moments_2, density_cc1


def run_percolation_te(graph, chr_name):
    graph_new = nx.Graph()
    edges_shuffled = list(graph.edges())
    rd.shuffle(edges_shuffled)

    roots = list(graph.nodes())
    uf = nx.utils.union_find.UnionFind(roots)

    data_size = len(edges_shuffled)
    k_values = np.empty(data_size)
    cc1_sizes = np.empty(data_size)
    density_cc1 = np.empty(data_size)
    cc2_sizes = np.empty(data_size)
    moments_2 = np.empty(data_size)
    # flag = np.zeros(len(sub))
    # ccd_percolation = np.zeros(len(sub))

    for k, (a, b) in enumerate(edges_shuffled):

        neb = list(graph.edges([a, b]))
        v1 = a, b
        neb.remove(v1)
        v2 = rd.choice(neb)
        c, d = v2
        neb2 = list(graph.edges([a, b, c, d]))
        neb2.remove(v1)
        neb2.remove(v2)
        neb2 = list(dict.fromkeys(neb2))
        v3 = rd.choice(neb2)

        r = rd.randint(0, 1200)
        if r % 3 == 0:
            v = v1
        elif r % 3 == 1:
            v = v2
        else:
            v = v3

        a, b = v
        graph.add_edge(a, b)
        w = graph[a][b]["weight"]
        graph_new.add_edge(a, b, weight=w)
        # print(k + 1)
        file_graph_write = open(
            f"/mnt/raid/kaustavs/Desktop/Project_data/perculation/percolation/data/graph_test/{chr_name}/graph_test_te/te_ctcf_2+_" + str(
                k + 1) + ".csv", "w")
        # file_graph_write.write("Source,Target,Weight\n")
        for (a1, b1) in list(graph_new.edges()):
            na = a1[0] + ":" + str(a1[1]) + "-" + str(a1[2])
            nb = b1[0] + ":" + str(b1[1]) + "-" + str(b1[2])
            x = na + "," + nb + "," + str(graph_new[a1][b1]["weight"])
            file_graph_write.write(x + "\n")
        file_graph_write.close()
        root_a = uf[a]
        root_b = uf[b]
        if root_a != root_b:
            uf.union(a, b)
            root_ab = uf[a]
            roots.remove(root_a)
            roots.remove(root_b)
            roots.append(root_ab)
        sizes = [uf.weights[v] for v in roots]
        unique_sizes, counts = np.unique(sizes, return_counts=True)
        moment_2 = ((unique_sizes ** 2) * counts).sum()
        if len(roots) > 1:
            cc1_size, cc2_size = nlargest(2, sizes)
        else:
            cc1_size, cc2_size = uf.weights[roots[0]], 0.0
        cc1_sizes[k] = cc1_size / len(list(graph.nodes()))
        density_cc1[k] = cc1_size / len(list(graph.nodes()))
        cc2_sizes[k] = cc2_size / len(list(graph.nodes()))
        k_values[k] = k + 1
        moments_2[k] = moment_2

    m = k_values.max()
    for i, values in enumerate(k_values):
        k_values[i] = values / m

    m = moments_2.max()
    for i, values in enumerate(moments_2):
        moments_2[i] = values / m

    return graph_new, k_values, cc1_sizes, cc2_sizes, moments_2, density_cc1


def run_percolation_fe(graph, chr_name):
    # for line in nx.generate_edgelist(graph):
    #     print(line)
    edges_shuffled = sorted(graph.edges(data=True), key=lambda t: t[2].get('weight'))
    edges_shuffled.reverse()
    roots = list(graph.nodes())
    uf = nx.utils.union_find.UnionFind(roots)

    graph_new = nx.Graph()
    data_size = len(edges_shuffled)
    k_values = np.empty(data_size)
    cc1_sizes = np.empty(data_size)
    density_cc1 = np.empty(data_size)
    cc2_sizes = np.empty(data_size)
    moments_2 = np.empty(data_size)
    count_ed = 0

    for k, (a, b, wt) in enumerate(edges_shuffled):

        graph.add_edge(a, b, weight=wt['weight'])
        # print(a, b, wt)
        count_ed += 1
        graph_new.add_edge(a, b, weight=wt['weight'])
        # print(k + 1)
        if count_ed % 1 == 0:
            file_graph_write = open(
                f"/mnt/raid/kaustavs/Desktop/Project_data/perculation/percolation/data/graph_test/{chr_name}/graph_test_fe/fe_ctcf_2+_" + str(
                    k + 1) + ".csv", "w")
            for (a_new, b_new) in list(graph_new.edges()):
                na = a_new[0] + ":" + str(a_new[1]) + "-" + str(a_new[2])
                nb = b_new[0] + ":" + str(b_new[1]) + "-" + str(b_new[2])
                wval = na + "," + nb + "," + str(graph_new[a_new][b_new]["weight"])
                # print(graph_new[a_new][b_new]["weight"])
                file_graph_write.write(wval + "\n")
            file_graph_write.close()
        root_a = uf[a]
        root_b = uf[b]
        if root_a != root_b:
            uf.union(a, b)
            root_ab = uf[a]
            roots.remove(root_a)
            roots.remove(root_b)
            roots.append(root_ab)
        sizes = [uf.weights[v] for v in roots]
        unique_sizes, counts = np.unique(sizes, return_counts=True)
        moment_2 = ((unique_sizes ** 2) * counts).sum()

        if len(roots) > 1:
            cc1_size, cc2_size = nlargest(2, sizes)
        else:
            cc1_size, cc2_size = uf.weights[roots[0]], 0.0

        cc1_sizes[k] = cc1_size / len(list(graph.nodes()))
        density_cc1[k] = cc1_size / len(list(graph.nodes()))
        cc2_sizes[k] = cc2_size / len(list(graph.nodes()))
        k_values[k] = k + 1
        moments_2[k] = moment_2

    m = k_values.max()
    for i, values in enumerate(k_values):
        k_values[i] = values / m

    m = moments_2.max()
    for i, values in enumerate(moments_2):
        moments_2[i] = values / m

    # print(k_values, cc1_sizes, cc2_sizes, moments_2, density_cc1)
    return graph_new, k_values, cc1_sizes, cc2_sizes, moments_2, density_cc1


def run_percolation_bio_model(graph, chr_name):
    graph_new = nx.Graph()
    edges_compartment = []
    edges_linear = []
    nodes = list(graph.nodes())
    edges_shuffled = list(graph.edges())
    for (a, b) in edges_shuffled:
        kind = graph[a][b]["kind"]
        if kind == 'A' or kind == 'B' or kind == "linear_A" or kind == "linear_B":
            buff = a + "," + b
            edges_compartment.append(buff)
        else:
            buff = a + "," + b
            edges_linear.append(buff)
    rd.shuffle(edges_compartment)
    print(edges_compartment[0])
    roots = nodes
    uf = nx.utils.union_find.UnionFind(roots)

    data_size = len(edges_shuffled)
    k_values = np.empty(data_size)
    cc1_sizes = np.empty(data_size)
    density_cc1 = np.empty(data_size)
    cc2_sizes = np.empty(data_size)
    moments_2 = np.empty(data_size)
    # flag = np.zeros(len(sub))
    # ccd_percolation = np.zeros(len(sub))

    file_write_edge = open("edge_order_compartment" + chr_name + ".txt", "w")
    for k, edge in enumerate(edges_compartment):
        a, b = edge.split(",")
        graph.add_edge(a, b)
        file_write_edge.write(edge + '\n')
        w = graph[a][b]["weight"]
        kind = graph[a][b]["kind"]
        graph_new.add_edge(a, b, weight=w, kind=kind)
        print(k)

        # file_graph_write = open(
        #     f"/mnt/raid/kaustavs/Desktop/Project_data/perculation/percolation/data/graph_test/{chr_name}/graph_test_compertment/compartment_" + chr_name + "_ctcf_2+_" + str(
        #         k) + ".csv", "w")
        # # file_graph_write.write("Source,Target,Weight\n")
        # for (x, y) in list(graph_new.edges()):
        #     z = x + "," + y + "," + str(graph_new[x][y]["weight"]) + "," + str(graph_new[x][y]['kind'])
        #     file_graph_write.write(z + "\n")
        # file_graph_write.close()

        root_a = uf[a]
        root_b = uf[b]
        if root_a != root_b:
            uf.union(a, b)
            root_ab = uf[a]
            roots.remove(root_a)
            roots.remove(root_b)
            roots.append(root_ab)
        sizes = [uf.weights[v] for v in roots]
        unique_sizes, counts = np.unique(sizes, return_counts=True)
        moment_2 = ((unique_sizes ** 2) * counts).sum()
        if len(roots) > 1:
            cc1_size, cc2_size = nlargest(2, sizes)
        else:
            cc1_size, cc2_size = uf.weights[roots[0]], 0.0
        cc1_sizes[k] = cc1_size / len(list(graph.nodes))
        density_cc1[k] = cc1_size / len(list(graph.nodes))
        cc2_sizes[k] = cc2_size / len(list(graph.nodes))
        k_values[k] = k + 1
        moments_2[k] = moment_2

    for edge in edges_linear:
        a, b = edge.split(",")
        graph.add_edge(a, b)
        file_write_edge.write(edge + '\n')
        w = graph[a][b]["weight"]
        kind = graph[a][b]["kind"]
        graph_new.add_edge(a, b, weight=w, kind=kind)
        k += 1
        # print(k)
        # file_graph_write1 = open(
        #     f"/mnt/raid/kaustavs/Desktop/Project_data/perculation/percolation/data/graph_test/{chr_name}/graph_test_compertment/compartment_" + chr_name + "_ctcf_2+_" + str(
        #         k) + ".csv", "w")
        # # file_graph_write.write("Source,Target,Weight\n")
        # for (x, y) in list(graph_new.edges()):
        #     z = x + "," + y + "," + str(graph_new[x][y]["weight"]) + "," + str(graph_new[x][y]['kind'])
        #     file_graph_write1.write(z + "\n")
        # file_graph_write1.close()

        root_a = uf[a]
        root_b = uf[b]
        if root_a != root_b:
            uf.union(a, b)
            root_ab = uf[a]
            roots.remove(root_a)
            roots.remove(root_b)
            roots.append(root_ab)
        sizes = [uf.weights[v] for v in roots]
        unique_sizes, counts = np.unique(sizes, return_counts=True)
        moment_2 = ((unique_sizes ** 2) * counts).sum()
        if len(roots) > 1:
            cc1_size, cc2_size = nlargest(2, sizes)
        else:
            cc1_size, cc2_size = uf.weights[roots[0]], 0.0
        cc1_sizes[k] = cc1_size / len(list(graph.nodes))
        density_cc1[k] = cc1_size / len(list(graph.nodes))
        cc2_sizes[k] = cc2_size / len(list(graph.nodes))
        k_values[k] = k + 1
        moments_2[k] = moment_2

    m = k_values.max()
    for i, values in enumerate(k_values):
        k_values[i] = values / m

    m = moments_2.max()
    for i, values in enumerate(moments_2):
        moments_2[i] = values / m

    print(cc1_sizes, cc2_sizes)
    return graph, k_values, cc1_sizes, cc2_sizes, moments_2, density_cc1


def run_percolation_epig_model(graph, chr_name):
    print("mapping Complete")
    graph_new = nx.Graph()
    edges_compartment_A = []
    edges_compartment_B = []
    edges_linear = []
    nodes = list(graph.nodes())
    edges_shuffled = list(graph.edges())
    for (a, b) in edges_shuffled:
        kind = graph[a][b]["kind"]
        # print(kind)
        if float(kind) >= 0:
            buff = a + "," + b
            edges_compartment_A.append(buff)
        elif float(kind) < -1:
            buff = a + "," + b
            edges_compartment_B.append(buff)
        else:
            buff = a + "," + b
            edges_linear.append(buff)

    edges_compartment_A = sorted(edges_compartment_A, key=lambda edges_compartment_B: edges_compartment_B[4])
    edges_compartment_B = sorted(edges_compartment_B, key=lambda edges_compartment_B: edges_compartment_B[4])
    print(len(edges_compartment_A))
    print(len(edges_compartment_B))
    # print(edges_compartment_A[0])
    roots = nodes
    uf = nx.utils.union_find.UnionFind(roots)

    data_size = len(edges_shuffled)
    k_values = np.empty(data_size)
    cc1_sizes = np.empty(data_size)
    density_cc1 = np.empty(data_size)
    cc2_sizes = np.empty(data_size)
    moments_2 = np.empty(data_size)
    # flag = np.zeros(len(sub))
    # ccd_percolation = np.zeros(len(sub))

    # min_len = min(len(edges_compartment_A), len(edges_compartment_B))
    # posa = 0
    # posb = 0
    # edge_comp = []

    # for i in range(0, min_len * 2):
    #     if i % 2 == 0:
    #         edge_comp.append(edges_compartment_A[posa])
    #         posa += 1
    #     else:
    #         edge_comp.append(edges_compartment_A[posb])
    #         posb += 1
    # sub = edges_compartment_A[min_len:]
    # for line in sub:
    #     edge_comp.append(line)

    graph_A = nx.Graph()
    file_write_edge = open("edge_order_compartment" + chr_name + ".csv", "w")
    for k, edge in enumerate(edges_compartment_A):
        a, b = edge.split(",")
        graph.add_edge(a, b)
        file_write_edge.write(edge + '\n')
        w = graph[a][b]["weight"]
        kind = graph[a][b]["kind"]

        graph_new.add_edge(a, b, weight=w, kind=kind)
        graph_A.add_edge(a, b, weight=w, kind=kind)
        # graph_viz.add_edge(nodes_int_map[a], nodes_int_map[b])
        # print(k)
        # file_graph_write = open(
        #     f"/mnt/raid/kaustavs/Desktop/Project_data/perculation/percolation/data/graph_test/{chr_name}/graph_test_epigenomic/epigenomic_" + chr_name + "_ctcf_2+_" + str(
        #         k) + ".csv", "w")
        # # file_graph_write.write("Source,Target,Weight\n")
        # for (x, y) in list(graph_new.edges()):
        #     z = x + "," + y + "," + str(graph_new[x][y]["weight"]) + "," + str(graph_new[x][y]['kind'])
        #     file_graph_write.write(z + "\n")
        # file_graph_write.close()

        root_a = uf[a]
        root_b = uf[b]
        if root_a != root_b:
            uf.union(a, b)
            root_ab = uf[a]
            roots.remove(root_a)
            roots.remove(root_b)
            roots.append(root_ab)
        sizes = [uf.weights[v] for v in roots]
        unique_sizes, counts = np.unique(sizes, return_counts=True)
        moment_2 = ((unique_sizes ** 2) * counts).sum()
        if len(roots) > 1:
            cc1_size, cc2_size = nlargest(2, sizes)
        else:
            cc1_size, cc2_size = uf.weights[roots[0]], 0.0
        cc1_sizes[k] = cc1_size / len(list(graph.nodes))
        density_cc1[k] = cc1_size / len(list(graph.nodes))
        cc2_sizes[k] = cc2_size / len(list(graph.nodes))
        k_values[k] = k + 1
        moments_2[k] = moment_2

    print("Edge Add Complete Compartment A")
    save_graph(graph_A, f'Epigenomic_method_compA', chr_name)

    graph_B = nx.Graph()
    for edge in edges_compartment_B:
        a, b = edge.split(",")
        graph.add_edge(a, b)
        file_write_edge.write(edge + '\n')
        w = graph[a][b]["weight"]
        kind = graph[a][b]["kind"]
        graph_new.add_edge(a, b, weight=w, kind=kind)
        graph_B.add_edge(a, b, weight=w, kind=kind)
        k += 1
        # print(k)
        # file_graph_write = open(
        #     f"/mnt/raid/kaustavs/Desktop/Project_data/perculation/percolation/data/graph_test/{chr_name}/graph_test_epigenomic/epigenomic_" + chr_name + "_ctcf_2+_" + str(
        #         k) + ".csv", "w")
        # # file_graph_write.write("Source,Target,Weight\n")
        # for (x, y) in list(graph_new.edges()):
        #     z = x + "," + y + "," + str(graph_new[x][y]["weight"]) + "," + str(graph_new[x][y]['kind'])
        #     file_graph_write.write(z + "\n")
        # file_graph_write.close()

        root_a = uf[a]
        root_b = uf[b]
        if root_a != root_b:
            uf.union(a, b)
            root_ab = uf[a]
            roots.remove(root_a)
            roots.remove(root_b)
            roots.append(root_ab)
        sizes = [uf.weights[v] for v in roots]
        unique_sizes, counts = np.unique(sizes, return_counts=True)
        moment_2 = ((unique_sizes ** 2) * counts).sum()
        if len(roots) > 1:
            cc1_size, cc2_size = nlargest(2, sizes)
        else:
            cc1_size, cc2_size = uf.weights[roots[0]], 0.0
        cc1_sizes[k] = cc1_size / len(list(graph.nodes))
        density_cc1[k] = cc1_size / len(list(graph.nodes))
        cc2_sizes[k] = cc2_size / len(list(graph.nodes))
        k_values[k] = k + 1
        moments_2[k] = moment_2

    print("Edge Add Complete Compartment B")
    print(len(list(graph_B.nodes())))
    save_graph(graph_B, f'Epigenomic_method_compB', chr_name)

    for edge in edges_linear:
        a, b = edge.split(",")
        graph.add_edge(a, b)
        file_write_edge.write(edge + '\n')
        w = graph[a][b]["weight"]
        kind = graph[a][b]["kind"]
        graph_new.add_edge(a, b, weight=w, kind=kind)
        k += 1
        # print(k)
        # file_graph_write1 = open(
        #     f"/mnt/raid/kaustavs/Desktop/Project_data/perculation/percolation/data/graph_test/{chr_name}/graph_test_epigenomic/epigenomic" + chr_name + "_ctcf_2+_" + str(
        #         k) + ".csv", "w")
        # # file_graph_write.write("Source,Target,Weight\n")
        # for (x, y) in list(graph_new.edges()):
        #     z = x + "," + y + "," + str(graph_new[x][y]["weight"]) + "," + str(graph_new[x][y]['kind'])
        #     file_graph_write1.write(z + "\n")
        # file_graph_write1.close()

        root_a = uf[a]
        root_b = uf[b]
        if root_a != root_b:
            uf.union(a, b)
            root_ab = uf[a]
            roots.remove(root_a)
            roots.remove(root_b)
            roots.append(root_ab)
        sizes = [uf.weights[v] for v in roots]
        unique_sizes, counts = np.unique(sizes, return_counts=True)
        moment_2 = ((unique_sizes ** 2) * counts).sum()
        if len(roots) > 1:
            cc1_size, cc2_size = nlargest(2, sizes)
        else:
            cc1_size, cc2_size = uf.weights[roots[0]], 0.0
        cc1_sizes[k] = cc1_size / len(list(graph.nodes))
        density_cc1[k] = cc1_size / len(list(graph.nodes))
        cc2_sizes[k] = cc2_size / len(list(graph.nodes))
        k_values[k] = k + 1
        moments_2[k] = moment_2

    m = k_values.max()
    for i, values in enumerate(k_values):
        k_values[i] = values / m

    m = moments_2.max()
    for i, values in enumerate(moments_2):
        moments_2[i] = values / m

    # print(cc1_sizes, cc2_sizes)
    return graph, k_values, cc1_sizes, cc2_sizes, moments_2, density_cc1


def peedictive_model(chromosome):
    file_order_parameter_loop_node = open(
        f"data/order_parameter_loop_node/order_parameter_loop_node_{chromosome}.csv", "r")

    nodes = []
    edges = []
    for line in file_order_parameter_loop_node:
        line = line.split("\n")[0]
        sub = line.split(",")
        edges.append([sub[0], sub[1], sub[2], sub[3], sub[4], sub[5]])
        nodes.append(sub[0])
        nodes.append(sub[1])

    nodes = list(dict.fromkeys(nodes))

    graph = nx.Graph()
    for node in nodes:
        graph.add_node(node)

    for edge in edges:
        graph.add_edge(edge[0], edge[1])
        graph.nodes[edge[0]]['op'] = float(edge[2])
        graph.nodes[edge[1]]['op'] = float(edge[3])
        graph.edges[edge[0], edge[1]]['weight'] = edge[4]
        graph.edges[edge[0], edge[1]]['op_edge'] = edge[5]
        node_a = np.array((float(edge[2]), float(edge[4]), float(edge[5])))
        node_b = np.array((1, 1, 1))
        graph.edges[edge[0], edge[1]]['diff'] = abs(float(edge[2]) - float(edge[3]))

    nodes = sorted(nodes)
    for node1, node2 in zip(nodes[:-1], nodes[1:]):
        if not graph.has_edge(node1, node2):
            graph.add_edge(node1, node2)
            graph.edges[node1, node2]['weight'] = 1

            # if graph.nodes[node1]['op'] >= 1 and graph.nodes[node2]['op'] >= 1:
            #     graph.edges[node1, node2]['op_edge'] = 1
            # elif graph.nodes[node1]['op'] < -1 and graph.nodes[node2]['op'] < -1:
            #     graph.edges[node1, node2]['op_edge'] = -1
            # else:
            #     graph.edges[node1, node2]['op_edge'] = 0
            graph.edges[node1, node2]['op_edge'] = 0
            graph.edges[node1, node2]['diff'] = abs(
                float(graph.nodes[node1]['op']) - float(graph.nodes[node2]['op']))

    print(len(list(graph.nodes)))
    print(len(list(graph.edges)))

    inter_compartmental = []
    compartment_a = []
    compartment_b = []

    for edge in list(graph.edges):
        a, b = edge
        op_edge = float(graph.edges[a, b]["op_edge"])
        op_node_a = float(graph.nodes[a]['op'])
        op_node_b = float(graph.nodes[b]['op'])
        # diff = float(graph.edges[a, b]['diff'])
        diff = abs(op_edge)
        if op_node_a >= 0 and op_node_b >= 0:
            compartment_a.append([a, b, op_edge, op_node_a, op_node_b, diff, graph.edges[a, b]["weight"], "A"])
        elif op_node_a < 0 and op_node_b < 0:
            compartment_a.append([a, b, op_edge, op_node_a, op_node_b, diff, graph.edges[a, b]["weight"], "A"])
            compartment_b.append([a, b, op_edge, op_node_a, op_node_b, diff, graph.edges[a, b]["weight"], "B"])
        else:
            inter_compartmental.append(
                [a, b, op_edge, op_node_a, op_node_b, diff, graph.edges[a, b]["weight"], "inter"])

    compartment_a = sorted(compartment_a, key=lambda x: x[5], reverse=True)
    print(compartment_a[0])
    print(len(compartment_a))
    print(len(compartment_b))
    # print(edges_compartment_A[0])
    graph_new = nx.Graph()

    roots = nodes
    uf = nx.utils.union_find.UnionFind(roots)

    data_size = len(list(graph.edges))
    k_values = np.empty(data_size)
    cc1_sizes = np.empty(data_size)
    density_cc1 = np.empty(data_size)
    cc2_sizes = np.empty(data_size)
    moments_2 = np.empty(data_size)
    # flag = np.zeros(len(sub))
    # ccd_percolation = np.zeros(len(sub))

    # min_len = min(len(edges_compartment_A), len(edges_compartment_B))
    # posa = 0
    # posb = 0
    # edge_comp = []

    # for i in range(0, min_len * 2):
    #     if i % 2 == 0:
    #         edge_comp.append(edges_compartment_A[posa])
    #         posa += 1
    #     else:
    #         edge_comp.append(edges_compartment_A[posb])
    #         posb += 1
    # sub = edges_compartment_A[min_len:]
    # for line in sub:
    #     edge_comp.append(line)

    graph_A = nx.Graph()
    # file_write_edge = open("edge_order_compartment" + str(chromosome) + ".csv", "w")
    for k, edge in enumerate(compartment_a):
        a = edge[0]
        b = edge[1]
        graph.add_edge(a, b)
        # file_write_edge.write(edge + '\n')
        w = graph[a][b]["weight"]
        kind = graph[a][b]["op_edge"]

        graph_new.add_edge(a, b, weight=w, kind=kind)
        graph_A.add_edge(a, b, weight=w, kind=kind)
        # graph_viz.add_edge(nodes_int_map[a], nodes_int_map[b])
        print(k)
        file_graph_write = open(
            f"/mnt/raid/kaustavs/Desktop/Project_data/perculation/percolation/data/graph_test/{chromosome}/graph_test_predictive/predictive_{chromosome}_ctcf_2+_" + str(
                k) + ".csv", "w")
        # file_graph_write.write("Source,Target,Weight\n")
        for (x, y) in list(graph_new.edges()):
            z = x + "," + y + "," + str(graph_new[x][y]["weight"]) + "," + str(graph_new[x][y]['kind'])
            file_graph_write.write(z + "\n")
        file_graph_write.close()

        root_a = uf[a]
        root_b = uf[b]
        if root_a != root_b:
            uf.union(a, b)
            root_ab = uf[a]
            roots.remove(root_a)
            roots.remove(root_b)
            roots.append(root_ab)
        sizes = [uf.weights[v] for v in roots]
        unique_sizes, counts = np.unique(sizes, return_counts=True)
        moment_2 = ((unique_sizes ** 2) * counts).sum()
        if len(roots) > 1:
            cc1_size, cc2_size = nlargest(2, sizes)
        else:
            cc1_size, cc2_size = uf.weights[roots[0]], 0.0
        cc1_sizes[k] = cc1_size / len(list(graph.nodes))
        density_cc1[k] = cc1_size / len(list(graph.nodes))
        cc2_sizes[k] = cc2_size / len(list(graph.nodes))
        k_values[k] = k + 1
        moments_2[k] = moment_2

    print("Edge Add Complete Compartment A")
    save_graph(graph_A, f'Epigenomic_method_nodes_compA', chromosome)

    # draw_plot(chromosome, 2, k_values, cc1_sizes, cc2_sizes, moments_2, "Epigenomic Nodes")

    # compartment_b = sorted(compartment_b, key=lambda x: x[5], reverse=True)
    # graph_B = nx.Graph()
    # for edge in compartment_b:
    #     a = edge[0]
    #     b = edge[1]
    #     graph.add_edge(a, b)
    #     # file_write_edge.write(edge + '\n')
    #     w = graph[a][b]["weight"]
    #     kind = graph[a][b]["op_edge"]
    #     graph_new.add_edge(a, b, weight=w, kind=kind)
    #     graph_B.add_edge(a, b, weight=w, kind=kind)
    #     k += 1
    #     # print(k)
    #     # file_graph_write = open(
    #     #     f"/mnt/raid/kaustavs/Desktop/Project_data/perculation/percolation/data/graph_test/{chr_name}/graph_test_epigenomic/epigenomic_" + chr_name + "_ctcf_2+_" + str(
    #     #         k) + ".csv", "w")
    #     # # file_graph_write.write("Source,Target,Weight\n")
    #     # for (x, y) in list(graph_new.edges()):
    #     #     z = x + "," + y + "," + str(graph_new[x][y]["weight"]) + "," + str(graph_new[x][y]['kind'])
    #     #     file_graph_write.write(z + "\n")
    #     # file_graph_write.close()
    #
    #     root_a = uf[a]
    #     root_b = uf[b]
    #     if root_a != root_b:
    #         uf.union(a, b)
    #         root_ab = uf[a]
    #         roots.remove(root_a)
    #         roots.remove(root_b)
    #         roots.append(root_ab)
    #     sizes = [uf.weights[v] for v in roots]
    #     unique_sizes, counts = np.unique(sizes, return_counts=True)
    #     moment_2 = ((unique_sizes ** 2) * counts).sum()
    #     if len(roots) > 1:
    #         cc1_size, cc2_size = nlargest(2, sizes)
    #     else:
    #         cc1_size, cc2_size = uf.weights[roots[0]], 0.0
    #     cc1_sizes[k] = cc1_size / len(list(graph.nodes))
    #     density_cc1[k] = cc1_size / len(list(graph.nodes))
    #     cc2_sizes[k] = cc2_size / len(list(graph.nodes))
    #     k_values[k] = k + 1
    #     moments_2[k] = moment_2
    #
    # print("Edge Add Complete Compartment B")
    # print(len(list(graph_B.nodes())))
    # save_graph(graph_B, f'Epigenomic_method_nodes_compB', chromosome)

    inter_compartmental = sorted(inter_compartmental, key=lambda x: x[5], reverse=False)
    for edge in inter_compartmental:
        a = edge[0]
        b = edge[1]
        graph.add_edge(a, b)
        # file_write_edge.write(edge + '\n')
        w = graph[a][b]["weight"]
        kind = graph[a][b]["op_edge"]
        graph_new.add_edge(a, b, weight=w, kind=kind)
        k += 1
        print(k)
        file_graph_write1 = open(
            f"/mnt/raid/kaustavs/Desktop/Project_data/perculation/percolation/data/graph_test/{chromosome}/graph_test_predictive/predictive_chr{chromosome}_ctcf_2+_" + str(
                k) + ".csv", "w")
        # file_graph_write.write("Source,Target,Weight\n")
        for (x, y) in list(graph_new.edges()):
            z = x + "," + y + "," + str(graph_new[x][y]["weight"]) + "," + str(graph_new[x][y]['kind'])
            file_graph_write1.write(z + "\n")
        file_graph_write1.close()

        root_a = uf[a]
        root_b = uf[b]
        if root_a != root_b:
            uf.union(a, b)
            root_ab = uf[a]
            roots.remove(root_a)
            roots.remove(root_b)
            roots.append(root_ab)
        sizes = [uf.weights[v] for v in roots]
        unique_sizes, counts = np.unique(sizes, return_counts=True)
        moment_2 = ((unique_sizes ** 2) * counts).sum()
        if len(roots) > 1:
            cc1_size, cc2_size = nlargest(2, sizes)
        else:
            cc1_size, cc2_size = uf.weights[roots[0]], 0.0
        cc1_sizes[k] = cc1_size / len(list(graph.nodes))
        density_cc1[k] = cc1_size / len(list(graph.nodes))
        cc2_sizes[k] = cc2_size / len(list(graph.nodes))
        k_values[k] = k + 1
        moments_2[k] = moment_2

    m = k_values.max()
    for i, values in enumerate(k_values):
        k_values[i] = values / m

    m = moments_2.max()
    for i, values in enumerate(moments_2):
        moments_2[i] = values / m

    # print(cc1_sizes, cc2_sizes)
    return graph, k_values, cc1_sizes, cc2_sizes, moments_2, density_cc1


def create_bio_graph(csv_file, chromosome, min_weight=2, linear_edge_weight=1):
    # Chromosome_no,Anchor_ID_A,Anchor_ID_B,Link_Weight
    t0 = time.time()
    # print("Input file : ", csv_file)
    # print("Chromosome : ", chromosome)
    # print("min_weight : ", min_weight)
    raw_df = pd.read_csv(
        csv_file,
        dtype={
            # 'Chromosome_no': 'str',  # for some files
            'Anchor_ID_A': 'str',
            'Anchor_ID_B': 'str',
            'Link_Weight': 'int'
        }
    )
    compertments = pd.read_csv(
        "data/GM12878_subcompartments.csv",
        dtype={
            # 'Chromosome_no': 'str',  # for some files
            'chr_name': 'str',
            'comp_start': 'int',
            'comp_end': 'int',
            'compertment': 'str',
            'a': 'int',
            'b': 'str',
            'd': 'int',
            'e': 'int',
            'f': 'int',
            'g': 'int'
        }
    )

    compertments_chr = []
    for _, (c_name, c_start, c_end, cmp, a, b, c, d, e, f, g) in compertments.iterrows():
        if chromosome == c_name:
            buff_comp_chr = c_name + "," + str(c_start) + "," + str(c_end) + "," + str(cmp)
            compertments_chr.append(buff_comp_chr)

    def _parse_anchor(anchor):
        chr_, rest = anchor.split(':')
        start, end = rest.split('-')
        return chr_, int(start), int(end)

    count_loops = 0
    loop_length = 0
    node_length = 0

    nodes = []
    loops = []
    comp = []
    node_compartment = {}
    inter_node_a = []
    inter_node_b = []
    start_anchor = 99999999
    end_anchor = 0

    for _, (a1, a2, w) in raw_df.iterrows():

        node1 = _parse_anchor(a1)
        node2 = _parse_anchor(a2)

        if w < min_weight:
            continue
        if not (node1[0] == chromosome and node2[0] == chromosome):
            continue

        if start_anchor > node1[1]:
            start_anchor = node1[1]
        if end_anchor < node2[2]:
            end_anchor = node2[2]

        count_loops += 1
        loop_length += (((node2[2] - node2[1]) / 2) - ((node1[2] - node1[1]) / 2))
        node_length += node2[2] - node2[1]
        node_length += node1[2] - node1[1]
        loop_buf = a1 + "," + a2 + "," + str(w)
        loops.append(loop_buf)
        nodes.append(a1)
        nodes.append(a2)

    nodes = list(dict.fromkeys(nodes))
    nodes.sort()
    print(f"Chromosome starts {chromosome}: ", start_anchor)
    print(f"Chromosome Ends {chromosome}: ", end_anchor)
    avg_node_length = node_length / (2 * count_loops)
    avg_loop_length = loop_length / count_loops

    for node in nodes:
        flag = 0
        anchor = _parse_anchor(node)
        chromosome_name = anchor[0]
        mid = int((int(anchor[1]) + int(anchor[2])) / 2)
        for comp_line in compertments_chr:
            c_name, c_start, c_end, cmp = comp_line.split(",")
            if int(c_start) <= mid <= int(c_end) and chromosome_name == c_name:
                buf = str(anchor[0]) + ":" + str(anchor[1]) + "-" + str(anchor[2]) + "," + str(cmp)
                comp.append(buf)
                node_compartment[anchor[0] + ":" + str(anchor[1]) + "-" + str(anchor[2])] = str(cmp)
                flag += 1
        # if flag == 0:
        #     # print(node)

    comp_write = open("./data/compartment_nodes/" + chromosome + ".csv", 'w')
    for line in comp:
        comp_write.write(line + "\n")
    comp_write.close()

    print(f"Total no of nodes in {chromosome}: ", len(nodes))
    print(f"Average node length in {chromosome}", avg_node_length)
    print(f"Total no of loops in {chromosome}: ", count_loops)
    print(f"Average Length of loops in {chromosome}:", avg_loop_length)
    print(f"Total no of mapped nodes in {chromosome}: ", len(comp))
    print(f"Total no of unmapped nodes in {chromosome}: ", len(nodes) - len(comp))

    # for loop in loops:
    #     sub_loop = loop.split(",")
    #     comp_a1 = comp_a2 = " "
    #     for node in comp:
    #         sub_node = node.split(",")
    #         if sub_loop[0] == sub_node[0]:

    for node in comp:
        sub_node = node.split(",")
        if sub_node[1][0] == "A":
            inter_node_a.append(sub_node[0])
        if sub_node[1][0] == "B":
            inter_node_b.append(sub_node[0])
    print(f"Total no of nodes in compartment A in {chromosome}: ", len(inter_node_a))
    print(f"Total no of nodes in compartment B in {chromosome}: ", len(inter_node_b))

    count_A_compertment_loops = 0
    count_B_compertment_loops = 0
    count_intercompartmental_loops = 0
    loops_A = []
    loops_B = []
    loops_intercompermental = []
    compartmental_loops = []
    for loop in loops:
        sub_loop = loop.split(",")
        if inter_node_a.__contains__(sub_loop[0]) and inter_node_a.__contains__(sub_loop[1]):
            count_A_compertment_loops += 1
            buf = loop + "," "A"
            compartmental_loops.append(buf)
            loops_A.append(loop)
        elif inter_node_b.__contains__(sub_loop[0]) and inter_node_b.__contains__(sub_loop[1]):
            count_B_compertment_loops += 1
            buf = loop + "," "B"
            compartmental_loops.append(buf)
            loops_B.append(loop)
        else:
            loops_intercompermental.append(loop)
            buf = loop + "," "inter"
            compartmental_loops.append(buf)
            count_intercompartmental_loops += 1

    count_same_compertment = count_A_compertment_loops + count_B_compertment_loops

    print(f"Total no of loops belong in same compertment in {chromosome}: ", count_same_compertment)
    print(f"Total inter-compartmental loops in {chromosome}: ", count_intercompartmental_loops)
    print(f"Total no of loops in compartment A in {chromosome}: ", count_A_compertment_loops)
    print(f"Total no of loops in compartment B in {chromosome}: ", count_B_compertment_loops)

    graph = nx.Graph()

    edges = []
    nodes = set()
    for edg in compartmental_loops:
        sub1 = edg.split(",")
        node1 = _parse_anchor(sub1[0])  # (chr, start, end)
        node2 = _parse_anchor(sub1[1])
        wei = sub1[2]
        comp = sub1[3]
        edges.append(node1[0] + ":" + str(node1[1]) + "-" + str(node1[2]) + "," +
                     node2[0] + ":" + str(node2[1]) + "-" + str(node2[2]) + "," + wei + "," + comp)
        nodes.add(node1[0] + ":" + str(node1[1]) + "-" + str(node1[2]))
        nodes.add(node2[0] + ":" + str(node2[1]) + "-" + str(node2[2]))

    graph.add_nodes_from(nodes)
    for loop in edges:
        # print(loop)
        sub = loop.split(",")
        graph.add_edge(sub[0], sub[1], weight=sub[2], kind=sub[3])

    edge_linear = []
    if linear_edge_weight is not None:
        nodes = sorted(nodes)
        for node1, node2 in zip(nodes[:-1], nodes[1:]):
            if not graph.has_edge(node1, node2):
                if node_compartment[node1] == node_compartment[node2]:
                    edge_linear.append(
                        node1 + "," + node2 + "," + str(linear_edge_weight) + ",linear_" + node_compartment[node1][0])
                    # print(node1 + "," + node2 + "," + str(linear_edge_weight) + ",linear_" + node_compartment[node1][0])
                else:
                    edge_linear.append(
                        node1 + "," + node2 + "," + str(linear_edge_weight) + ",linear_inter")

    graph.add_nodes_from(nodes)
    for loop in edge_linear:
        # print(loop)
        sub = loop.split(",")
        graph.add_edge(sub[0], sub[1], weight=sub[2], kind=sub[3])
    print(nx.number_connected_components(graph))

    print(
        f'Created network for {chromosome} with |V|={graph.number_of_nodes()}, |E|= {graph.number_of_edges()} ({time.time() - t0:.2f}s)')
    # print(graph.edges.data())
    file_wr = open("./data/compartment_graph/rao_graph_without_linear_edges_" + chromosome + ".csv", "w")
    for line in edges:
        # print(line)
        file_wr.write(line + "\n")
    file_wr.close()
    cyto = json_graph.node_link_data(graph)
    jn = js.dumps(cyto)
    f = open("./data/compartment_graph/rao_graph_without_linear_edges_" + chromosome + ".json", "w")
    f.write(jn)
    f.close()
    return graph


def create_epig_graph(csv_file, chromosome, min_weight=2, linear_edge_weight=1):
    # Chromosome_no,Anchor_ID_A,Anchor_ID_B,Link_Weight
    t0 = time.time()
    print("Input file : ", csv_file)
    print("Chromosome : ", chromosome)
    print("min_weight : ", min_weight)
    raw_df = pd.read_csv(
        csv_file,
        dtype={
            # 'Chromosome_no': 'str',  # for some files
            'Anchor_ID_A': 'str',
            'Anchor_ID_B': 'str',
            'Link_Weight': 'int'
        }
    )
    compertments = pd.read_csv(
        "data/GM12878_subcompartments.csv",
        dtype={
            # 'Chromosome_no': 'str',  # for some files
            'chr_name': 'str',
            'comp_start': 'int',
            'comp_end': 'int',
            'compertment': 'str',
            'a': 'int',
            'b': 'str',
            'd': 'int',
            'e': 'int',
            'f': 'int',
            'g': 'int'
        }
    )

    compertments_chr = []
    for _, (c_name, c_start, c_end, cmp, a, b, c, d, e, f, g) in compertments.iterrows():
        if chromosome == c_name:
            buff_comp_chr = c_name + "," + str(c_start) + "," + str(c_end) + "," + str(cmp)
            compertments_chr.append(buff_comp_chr)

    def _parse_anchor(anchor):
        chr_, rest = anchor.split(':')
        start, end = rest.split('-')
        return chr_, int(start), int(end)

    count_loops = 0
    loop_length = 0
    node_length = 0

    nodes = []
    loops = []
    comp = []
    node_compartment = {}
    inter_node_a = []
    inter_node_b = []
    start_anchor = 99999999
    end_anchor = 0

    for _, (a1, a2, w) in raw_df.iterrows():

        node1 = _parse_anchor(a1)
        node2 = _parse_anchor(a2)

        if w < min_weight:
            continue
        if not (node1[0] == chromosome and node2[0] == chromosome):
            continue

        if start_anchor > node1[1]:
            start_anchor = node1[1]
        if end_anchor < node2[2]:
            end_anchor = node2[2]

        count_loops += 1
        loop_length += (((node2[2] - node2[1]) / 2) - ((node1[2] - node1[1]) / 2))
        node_length += node2[2] - node2[1]
        node_length += node1[2] - node1[1]
        loop_buf = a1 + "," + a2 + "," + str(w)
        loops.append(loop_buf)
        nodes.append(a1)
        nodes.append(a2)

    nodes = list(dict.fromkeys(nodes))
    nodes.sort()
    print(f"Chromosome starts {chromosome}: ", start_anchor)
    print(f"Chromosome Ends {chromosome}: ", end_anchor)
    avg_node_length = node_length / (2 * count_loops)
    avg_loop_length = loop_length / count_loops

    for node in nodes:
        flag = 0
        anchor = _parse_anchor(node)
        chromosome_name = anchor[0]
        mid = int((int(anchor[1]) + int(anchor[2])) / 2)
        for comp_line in compertments_chr:
            c_name, c_start, c_end, cmp = comp_line.split(",")
            if int(c_start) <= mid <= int(c_end) and chromosome_name == c_name:
                buf = str(anchor[0]) + ":" + str(anchor[1]) + "-" + str(anchor[2]) + "," + str(cmp)
                comp.append(buf)
                node_compartment[anchor[0] + ":" + str(anchor[1]) + "-" + str(anchor[2])] = str(cmp)
                flag += 1
        # if flag == 0:
        #     # print(node)

    # comp_write = open("./data/compartment_nodes/" + chromosome + ".csv", 'w')
    # for line in comp:
    #     comp_write.write(line + "\n")
    # comp_write.close()

    print(f"Total no of nodes in {chromosome}: ", len(nodes))
    print(f"Average node length in {chromosome}", avg_node_length)
    print(f"Total no of loops in {chromosome}: ", count_loops)
    print(f"Average Length of loops in {chromosome}:", avg_loop_length)
    print(f"Total no of mapped nodes in {chromosome}: ", len(comp))
    print(f"Total no of unmapped nodes in {chromosome}: ", len(nodes) - len(comp))

    # for loop in loops:
    #     sub_loop = loop.split(",")
    #     comp_a1 = comp_a2 = " "
    #     for node in comp:
    #         sub_node = node.split(",")
    #         if sub_loop[0] == sub_node[0]:

    for node in comp:
        sub_node = node.split(",")
        if sub_node[1][0] == "A":
            inter_node_a.append(sub_node[0])
        if sub_node[1][0] == "B":
            inter_node_b.append(sub_node[0])
    print(f"Total no of nodes in compartment A in {chromosome}: ", len(inter_node_a))
    print(f"Total no of nodes in compartment B in {chromosome}: ", len(inter_node_b))

    count_A_compertment_loops = 0
    count_B_compertment_loops = 0

    loops_A = []
    loops_B = []
    loops_intercompermental = []
    compartmental_loops = []

    file_epigenomic = open("data/Order.parameter/epigenomic_driven_edges_" + chromosome + ".csv", "r")

    for line in file_epigenomic:
        line = line.split("\n")[0]
        sub = line.split(",")
        if float(sub[3]) >= 0:
            count_A_compertment_loops += 1
            loops_A.append(line)
            compartmental_loops.append(line + "," + "A")
        else:
            count_B_compertment_loops += 1
            loops_B.append(line)
            compartmental_loops.append(line + "," + "B")

    # for loop in loops:
    #     sub_loop = loop.split(",")
    #     if inter_node_a.__contains__(sub_loop[0]) and inter_node_a.__contains__(sub_loop[1]):
    #         count_A_compertment_loops += 1
    #         buf = loop + "," "A"
    #         compartmental_loops.append(buf)
    #         loops_A.append(loop)
    #     elif inter_node_b.__contains__(sub_loop[0]) and inter_node_b.__contains__(sub_loop[1]):
    #         count_B_compertment_loops += 1
    #         buf = loop + "," "B"
    #         compartmental_loops.append(buf)
    #         loops_B.append(loop)
    #     else:
    #         loops_intercompermental.append(loop)
    #         buf = loop + "," "inter"
    #         compartmental_loops.append(buf)
    #         count_intercompartmental_loops += 1

    count_same_compertment = count_A_compertment_loops + count_B_compertment_loops

    print(f"Total no of loops belong in same compertment in {chromosome}: ", count_same_compertment)
    print(f"Total no of loops in compartment A in {chromosome}: ", count_A_compertment_loops)
    print(f"Total no of loops in compartment B in {chromosome}: ", count_B_compertment_loops)

    graph = nx.Graph()

    edges = []
    nodes = set()
    for edg in compartmental_loops:
        sub1 = edg.split(",")
        node1 = _parse_anchor(sub1[0])  # (chr, start, end)
        node2 = _parse_anchor(sub1[1])
        wei = sub1[2]
        confidence_score = sub1[3]
        comp = sub1[4]
        edges.append(node1[0] + ":" + str(node1[1]) + "-" + str(node1[2]) + "," +
                     node2[0] + ":" + str(node2[1]) + "-" + str(
            node2[2]) + "," + wei + "," + "," + confidence_score + "," + comp + "," + comp + "," + comp)
        nodes.add(node1[0] + ":" + str(node1[1]) + "-" + str(node1[2]))
        nodes.add(node2[0] + ":" + str(node2[1]) + "-" + str(node2[2]))

    graph.add_nodes_from(nodes)
    for loop in edges:
        # print(loop)
        sub = loop.split(",")
        graph.add_edge(sub[0], sub[1], weight=sub[3], kind=sub[4])

    edge_linear = []
    if linear_edge_weight is not None:
        nodes = sorted(nodes)
        for node1, node2 in zip(nodes[:-1], nodes[1:]):
            if not graph.has_edge(node1, node2):
                if node_compartment[node1] == node_compartment[node2] == "A":
                    edge_linear.append(
                        node1 + "," + node2 + "," + str(linear_edge_weight) + ",1")
                    edges.append(
                        node1 + "," + node2 + "," + str(linear_edge_weight) + ",1")
                elif node_compartment[node1] == node_compartment[node2] == "B":
                    edge_linear.append(
                        node1 + "," + node2 + "," + str(linear_edge_weight) + ",-1")
                    edges.append(
                        node1 + "," + node2 + "," + str(linear_edge_weight) + ",-1")

                else:
                    edge_linear.append(
                        node1 + "," + node2 + "," + str(linear_edge_weight) + ",0")
                    edges.append(node1 + "," + node2 + "," + str(
                        linear_edge_weight) + ",0")

    graph.add_nodes_from(nodes)
    for loop in edge_linear:
        # print(loop)
        sub = loop.split(",")
        graph.add_edge(sub[0], sub[1], weight=sub[2], kind=sub[3])
    print(nx.number_connected_components(graph))

    print(
        f'Created network for {chromosome} with |V|={graph.number_of_nodes()}, |E|= {graph.number_of_edges()} ({time.time() - t0:.2f}s)')
    # print(graph.edges.data())
    file_wr = open("./data/epigenome_graph_with_linear_edges_" + chromosome + ".csv", "w")
    for line in edges:
        # print(line)
        file_wr.write(line + "\n")
    file_wr.close()
    # cyto = json_graph.node_link_data(graph)
    # jn = js.dumps(cyto)
    # f = open("./data/compartment_graph/rao_graph_without_linear_edges_" + chromosome + ".json", "w")
    # f.write(jn)
    # f.close()
    return graph


def save_results(file_path, *results):
    np.savetxt(file_path, np.concatenate(results).T)


def save_graph(graph, type, chr_name):
    nx.write_adjlist(graph, f'results/{chr_name}/graph_{type}.adjlist')


def scalar_value_calculate_loops(chrN):
    file = 'data/hg19_PET2+_featuretable/'
    file = file + 'CTCFplus2.' + chrN + '.featuretable.txt'

    df = pd.read_csv(file, sep="\t")
    X = df.iloc[:, 11:].values
    y_A = df["A.merged_fraction"].values
    y_B = df["B.merged_fraction"].values

    scaler = preprocessing.StandardScaler()
    X = scaler.fit_transform(X)

    # compute and plot distribution of scalar parameter
    ind_all = np.arange(len(y_A))
    # ind_trn, ind_tst = (y_A==1.0) | (y_B==1.0), ind_all #case1
    ind_trn, ind_tst = (y_A + y_B >= 0.90), ind_all  # case2
    # y_bin = y_A > 0.5
    y = y_A > y_B
    X_trn, y_trn = X[ind_trn, :], y[ind_trn]
    X_tst, y_tst = X[ind_tst, :], y[ind_tst]

    # build the model
    lda = LinearDiscriminantAnalysis(n_components=1)
    lda = lda.fit(X_trn, y_trn)

    sns.set()
    plt.figure(figsize=(8, 6))

    # plot for Compartment A
    # X_A = X_trn[y_trn>0.5, :]
    X_A = X[(y_A == 1.0) & (y_B == 0.0), :]
    X_ldaA = lda.transform(X_A)
    # for i in range(len(X_ldaA)):
    #    X_ldaA[i] = -10 if X_ldaA[i] < -10 else X_ldaA[i]
    #    X_ldaA[i] = +10 if X_ldaA[i] > +10 else X_ldaA[i]
    data = pd.DataFrame(X_ldaA, columns=['Compartment A'])
    sns.kdeplot(data['Compartment A'], shade=True, legend=True)

    # plot for Compartment B
    # X_B = X_trn[y_trn<=0.5, :]
    X_B = X[(y_A == 0.0) & (y_B == 1.0), :]
    X_ldaB = lda.transform(X_B)
    # for i in range(len(X_ldaB)):
    #    X_ldaB[i] = -10 if X_ldaB[i] < -10 else X_ldaB[i]
    #    X_ldaB[i] = +10 if X_ldaB[i] > +10 else X_ldaB[i]
    data = pd.DataFrame(X_ldaB, columns=['Compartment B'])
    sns.kdeplot(data['Compartment B'], shade=True, legend=True)

    # plot for mixed-mode loops
    X_M = X[(y_A != 1.0) & (y_B != 1.0), :]
    X_ldaM = lda.transform(X_M)
    for i in range(len(X_ldaM)):
        X_ldaM[i] = -10 if X_ldaM[i] < -10 else X_ldaM[i]
        X_ldaM[i] = +10 if X_ldaM[i] > +10 else X_ldaM[i]
    data = pd.DataFrame(X_ldaM, columns=['Mixed-mode'])
    sns.kdeplot(data['Mixed-mode'], shade=True, legend=True)

    # plot for all loops/regions
    X_lda = lda.transform(X_tst)
    X_cln = X_lda
    for i in range(len(X_cln)):
        X_cln[i] = -10 if X_cln[i] < -10 else X_cln[i]
        X_cln[i] = +10 if X_cln[i] > +10 else X_cln[i]
    data = pd.DataFrame(X_cln, columns=['All'])
    # sns.kdeplot(data['All'], shade=True, legend=True)

    plt.xlabel('Scalar value')
    plt.xlim([-7, 7])
    plt.ylabel('Frequency')
    title = 'Scalar value distribution Plot for loops'
    plt.title(title, fontsize=14)
    file = 'outdir/loops/CTCFplus2.' + chrN + '.dist.png'
    plt.savefig(file)
    plt.show()
    plt.close()

    f1 = df.iloc[:, :5]
    is_trn = np.zeros(len(y_A))
    is_trn[ind_trn] = 1
    is_trn = is_trn.reshape(len(y_A), 1)
    f2 = pd.DataFrame(np.concatenate((X_lda, is_trn), axis=1),
                      columns=['Order.parameter', 'Training.sample'])
    f3 = pd.concat([f1, f2], axis=1)
    file = 'outdir/loops/CTCFplus2.' + chrN + '.score.csv'
    f3.to_csv(file, index=False)


def scalar_value_calculate_anchors(chrN):
    file = 'data/hg19_PET2+_featuretable/'
    file = file + 'CTCFplus2.' + chrN + '.featuretable.txt'

    if chrN == "chr1-23":
        chrN = "chr1-22"

    df = pd.read_csv(file, sep="\t")
    X = df.iloc[:, 11:].values
    y_A = df["A.merged_fraction"].values
    y_B = df["B.merged_fraction"].values

    scaler = preprocessing.StandardScaler()
    X = scaler.fit_transform(X)

    # compute and plot distribution of scalar parameter
    ind_all = np.arange(len(y_A))
    # ind_trn, ind_tst = (y_A==1.0) | (y_B==1.0), ind_all #case1
    ind_trn, ind_tst = (y_A + y_B >= 0.90), ind_all  # case2
    # y_bin = y_A > 0.5
    y = y_A > y_B
    X_trn, y_trn = X[ind_trn, :], y[ind_trn]
    X_tst, y_tst = X[ind_tst, :], y[ind_tst]

    # build the model
    lda = LinearDiscriminantAnalysis(n_components=1)
    lda = lda.fit(X_trn, y_trn)

    sns.set()
    plt.figure(figsize=(8, 6))

    # plot for Compartment A
    # X_A = X_trn[y_trn>0.5, :]
    X_A = X[(y_A == 1.0) & (y_B == 0.0), :]
    X_ldaA = lda.transform(X_A)
    # for i in range(len(X_ldaA)):
    #    X_ldaA[i] = -10 if X_ldaA[i] < -10 else X_ldaA[i]
    #    X_ldaA[i] = +10 if X_ldaA[i] > +10 else X_ldaA[i]
    data = pd.DataFrame(X_ldaA, columns=['Compartment A'])
    sns.kdeplot(data['Compartment A'], shade=True, legend=True)

    # plot for Compartment B
    # X_B = X_trn[y_trn<=0.5, :]
    X_B = X[(y_A == 0.0) & (y_B == 1.0), :]
    X_ldaB = lda.transform(X_B)
    # for i in range(len(X_ldaB)):
    #    X_ldaB[i] = -10 if X_ldaB[i] < -10 else X_ldaB[i]
    #    X_ldaB[i] = +10 if X_ldaB[i] > +10 else X_ldaB[i]
    data = pd.DataFrame(X_ldaB, columns=['Compartment B'])
    sns.kdeplot(data['Compartment B'], shade=True, legend=True)

    # plot for mixed-mode loops
    X_M = X[(y_A != 1.0) & (y_B != 1.0), :]
    X_ldaM = lda.transform(X_M)
    for i in range(len(X_ldaM)):
        X_ldaM[i] = -10 if X_ldaM[i] < -10 else X_ldaM[i]
        X_ldaM[i] = +10 if X_ldaM[i] > +10 else X_ldaM[i]
    data = pd.DataFrame(X_ldaM, columns=['Mixed-mode'])
    sns.kdeplot(data['Mixed-mode'], shade=True, legend=True)

    # plot for all loops/regions
    X_lda = lda.transform(X_tst)
    X_cln = X_lda
    for i in range(len(X_cln)):
        X_cln[i] = -10 if X_cln[i] < -10 else X_cln[i]
        X_cln[i] = +10 if X_cln[i] > +10 else X_cln[i]
    data = pd.DataFrame(X_cln, columns=['All'])
    # sns.kdeplot(data['All'], shade=True, legend=True)

    plt.xlabel('Scalar value')
    plt.xlim([-7, 7])
    plt.ylabel('Frequency')
    title = 'Scalar value distribution Plot for loops'
    plt.title(title, fontsize=14)
    file = 'outdir/anchors/CTCFplus2.' + chrN + '.dist.png'
    plt.savefig(file)
    plt.show()
    plt.close()

    f1 = df.iloc[:, :5]
    is_trn = np.zeros(len(y_A))
    is_trn[ind_trn] = 1
    is_trn = is_trn.reshape(len(y_A), 1)
    f2 = pd.DataFrame(np.concatenate((X_lda, is_trn), axis=1),
                      columns=['Order.parameter', 'Training.sample'])
    f3 = pd.concat([f1, f2], axis=1)
    file = 'outdir/anchors/CTCFplus2.' + chrN + '.score.csv'
    f3.to_csv(file, index=False)


def run_percolation(chr_name, graph, graph_bio, graph_epigi):
    global graph_per
    random_run = 1
    t0 = time.time()
    data_size = len(list(graph.edges))
    k_values = np.zeros(data_size)
    cc1_sizes = np.zeros(data_size)
    density_cc1 = np.zeros(data_size)
    cc2_sizes = np.zeros(data_size)
    moments_2 = np.zeros(data_size)

    for no_of_iteration in tqdm(range(0, random_run)):
        graph_per, k_values_per, cc1_sizes_per, cc2_sizes_per, moments_2_per, density_cc1_per = run_percolation_er(
            graph, chr_name)
        k_values = [a + b for a, b in zip(k_values, k_values_per)]
        cc1_sizes = [a + b for a, b in zip(cc1_sizes, cc1_sizes_per)]
        cc2_sizes = [a + b for a, b in zip(cc2_sizes, cc2_sizes_per)]
        moments_2 = [a + b for a, b in zip(moments_2, moments_2_per)]
        density_cc1 = [a + b for a, b in zip(density_cc1, density_cc1_per)]

    k_values = np.divide(k_values, random_run)
    cc1_sizes = np.divide(cc1_sizes, random_run)
    cc2_sizes = np.divide(cc2_sizes, random_run)
    moments_2 = np.divide(moments_2, random_run)
    density_cc1 = np.divide(density_cc1, random_run)
    results_er = k_values, cc1_sizes, cc2_sizes, moments_2, density_cc1
    print(f'Simulation completed for {chr_name} Erdos Renyi ({time.time() - t0:.2f}s)')
    results_filename = f'results/{chr_name}/results_Erdos_Renyi.txt'
    save_results(results_filename, results_er)
    save_graph(graph_per, f'Erdos_Renyi', chr_name)
    print("No of nodes in ER :", len(list(graph_per.nodes())))
    print("No of edges in ER :", len(list(graph_per.edges())))
    #
    t1 = time.time()
    graph_fe, k_values, cc1_sizes, cc2_sizes, moments_2, density_cc1 = run_percolation_fe(graph, chr_name)
    results_fe = k_values, cc1_sizes, cc2_sizes, moments_2, density_cc1
    print(f'Simulation completed for {chr_name} Frequency Edges ({time.time() - t1:.2f}s)')
    results_filename = f'results/{chr_name}/results_Frequency_Edges.txt'
    save_results(results_filename, results_fe)
    save_graph(graph_fe, f'Frequency_Edges', chr_name)
    print("No of nodes in FE :", len(list(graph_fe.nodes())))
    print("No of edges in FE :", len(list(graph_fe.edges())))

    t2 = time.time()
    data_size = len(list(graph.edges))
    k_values = np.zeros(data_size)
    cc1_sizes = np.zeros(data_size)
    density_cc1 = np.zeros(data_size)
    cc2_sizes = np.zeros(data_size)
    moments_2 = np.zeros(data_size)

    for no_of_iteration in tqdm(range(0, random_run)):
        graph_per, k_values_per, cc1_sizes_per, cc2_sizes_per, moments_2_per, density_cc1_per = run_percolation_ae(
            graph, chr_name)
        k_values = [a + b for a, b in zip(k_values, k_values_per)]
        cc1_sizes = [a + b for a, b in zip(cc1_sizes, cc1_sizes_per)]
        cc2_sizes = [a + b for a, b in zip(cc2_sizes, cc2_sizes_per)]
        moments_2 = [a + b for a, b in zip(moments_2, moments_2_per)]
        density_cc1 = [a + b for a, b in zip(density_cc1, density_cc1_per)]

    k_values = np.divide(k_values, random_run)
    cc1_sizes = np.divide(cc1_sizes, random_run)
    cc2_sizes = np.divide(cc2_sizes, random_run)
    moments_2 = np.divide(moments_2, random_run)
    density_cc1 = np.divide(density_cc1, random_run)
    results_ae = k_values, cc1_sizes, cc2_sizes, moments_2, density_cc1
    print(f'Simulation completed for {chr_name} Adjacent Edges ({time.time() - t2:.2f}s)')
    results_filename = f'results/{chr_name}/results_Adjacent_Edges.txt'
    save_results(results_filename, results_ae)
    save_graph(graph_per, f'Adjacent_Edges', chr_name)
    print("No of nodes in AE :", len(list(graph_per.nodes())))
    print("No of edges in AE :", len(list(graph_per.edges())))

    t3 = time.time()
    data_size = len(list(graph.edges))
    k_values = np.zeros(data_size)
    cc1_sizes = np.zeros(data_size)
    density_cc1 = np.zeros(data_size)
    cc2_sizes = np.zeros(data_size)
    moments_2 = np.zeros(data_size)

    for no_of_iteration in tqdm(range(0, random_run)):
        graph_per, k_values_per, cc1_sizes_per, cc2_sizes_per, moments_2_per, density_cc1_per = run_percolation_te(
            graph, chr_name)
        k_values = [a + b for a, b in zip(k_values, k_values_per)]
        cc1_sizes = [a + b for a, b in zip(cc1_sizes, cc1_sizes_per)]
        cc2_sizes = [a + b for a, b in zip(cc2_sizes, cc2_sizes_per)]
        moments_2 = [a + b for a, b in zip(moments_2, moments_2_per)]
        density_cc1 = [a + b for a, b in zip(density_cc1, density_cc1_per)]

    k_values = np.divide(k_values, random_run)
    cc1_sizes = np.divide(cc1_sizes, random_run)
    cc2_sizes = np.divide(cc2_sizes, random_run)
    moments_2 = np.divide(moments_2, random_run)
    density_cc1 = np.divide(density_cc1, random_run)
    results_te = k_values, cc1_sizes, cc2_sizes, moments_2, density_cc1
    print(f'Simulation completed for {chr_name} Triangular Edges ({time.time() - t3:.2f}s)')
    results_filename = f'results/{chr_name}/results_Triangular_Edges.txt'
    save_results(results_filename, results_te)
    save_graph(graph_per, f'Triangular_Edges', chr_name)
    print("No of nodes in TR :", len(list(graph_per.nodes())))
    print("No of edges in TR :", len(list(graph_per.edges())))

    t4 = time.time()
    graph_bio_per, k_values_per, cc1_sizes_per, cc2_sizes_per, moments_2_per, density_cc1_per \
        = run_percolation_bio_model(graph_bio, chr_name)
    results_per = k_values_per, cc1_sizes_per, cc2_sizes_per, moments_2_per, density_cc1_per
    print(f'Simulation completed for {chr_name} Biological Model ({time.time() - t4:.2f}s)')
    results_filename = f'results/{chr_name}/results_Biological_Model.txt'
    save_results(results_filename, results_per)
    save_graph(graph_bio_per, f'Biological_Model', chr_name)
    print("No of nodes in BM :", len(list(graph_bio_per.nodes())))
    print("No of edges in BM :", len(list(graph_bio_per.edges())))

    t5 = time.time()
    graph_epige_per, k_values_per, cc1_sizes_per, cc2_sizes_per, moments_2_per, density_cc1_per \
        = run_percolation_epig_model(graph_epigi, chr_name)
    results_per = k_values_per, cc1_sizes_per, cc2_sizes_per, moments_2_per, density_cc1_per
    print(f'Simulation completed for {chr_name} Epigenomic Model ({time.time() - t5:.2f}s)')
    results_filename = f'results/{chr_name}/results_Epigenomic_Model.txt'
    save_results(results_filename, results_per)
    save_graph(graph_epige_per, f'Epigenomic_Model', chr_name)
    print("No of nodes in EM :", len(list(graph_epige_per.nodes())))
    print("No of edges in EM :", len(list(graph_epige_per.edges())))

    t5 = time.time()
    graph_predictive_per, k_values_per, cc1_sizes_per, cc2_sizes_per, moments_2_per, density_cc1_per \
        = peedictive_model(chr_name)
    results_per = k_values_per, cc1_sizes_per, cc2_sizes_per, moments_2_per, density_cc1_per
    print(f'Simulation completed for {chr_name} Predictive Model ({time.time() - t5:.2f}s)')
    results_filename = f'results/{chr_name}/results_Predictive_Model.txt'
    save_results(results_filename, results_per)
    save_graph(graph_predictive_per, f'Predictive_Model', chr_name)
    print("No of nodes in PM :", len(list(graph_predictive_per.nodes())))
    print("No of edges in PM :", len(list(graph_predictive_per.edges())))

    # return graph_er, graph_ae, graph_te, graph_fe


def load_data(chr_name, file_path):
    X = np.loadtxt(f'results/{chr_name}/results_{file_path}.txt')
    results = [
        X[:, i]
        for i in range(X.shape[1])
    ]

    return results


def run_draw_plot(chr_name, er, ae, te, fe, bm, em, pm, graph_er, graph_ae, graph_te, graph_fe, graph_bm, graph_epigi,
                  graph_pm, start=70000, end=85000):
    k_value_er, cc1_er, cc2_er, m_er, density_er = load_data(chr_name, er)
    draw_plot(chr_name, graph_er.number_of_nodes(), k_value_er, cc1_er, cc2_er, m_er, "Erdös-Rényi (ER) model")
    # start = 0
    # end = critical_point(graph, chr_name, k_value, cc1, cc2, m) - 6000
    # draw_plot_part(start, end, chr_name, graph_er.number_of_nodes(), k_value_er, cc1_er, cc2_er, m_er, er)

    k_value_ae, cc1_ae, cc2_ae, m_ae, density_ae = load_data(chr_name, ae)
    draw_plot(chr_name, graph_ae.number_of_nodes(), k_value_ae, cc1_ae, cc2_ae, m_ae, "Adjacent Edge (AE) model")
    # start = 0
    # end = critical_point(graph, chr_name, k_value, cc1, cc2, m) - 6000
    # draw_plot_part(start, end, chr_name, graph_ae.number_of_nodes(), k_value_ae, cc1_ae, cc2_ae, m_ae, ae)

    k_value_te, cc1_te, cc2_te, m_te, density_te = load_data(chr_name, te)
    draw_plot(chr_name, graph_te.number_of_nodes(), k_value_te, cc1_te, cc2_te, m_te, "Triple Edge (TE) model")
    # start = 0
    # end = critical_point(graph, chr_name, k_value, cc1, cc2, m) - 6000
    # draw_plot_part(start, end, chr_name, graph_te.number_of_nodes(), k_value_te, cc1_te, cc2_te, m_te, te)

    k_value_fe, cc1_fe, cc2_fe, m_fe, density_fe = load_data(chr_name, fe)
    draw_plot(chr_name, graph_fe.number_of_nodes(), k_value_fe, cc1_fe, cc2_fe, m_fe, "Loop Frequency (LF) model")
    # start = 0
    # end = critical_point(graph, chr_name, k_value, cc1, cc2, m) - 6000
    # draw_plot_part(start, end, chr_name, graph_fe.number_of_nodes(), k_value_fe, cc1_fe, cc2_fe, m_fe, fe)

    k_value_bm, cc1_bm, cc2_bm, m_bm, density_bm = load_data(chr_name, bm)
    draw_plot(chr_name, graph_bm.number_of_nodes(), k_value_bm, cc1_bm, cc2_bm, m_bm,
              "Chromatin Compartment (CC) model")
    # start = 0
    # end = critical_point(graph, chr_name, k_value, cc1, cc2, m) - 6000
    # draw_plot_part(start, end, chr_name, graph_fe.number_of_nodes(), k_value_fe, cc1_fe, cc2_fe, m_fe, fe)
    k_value_em, cc1_em, cc2_em, m_em, density_em = load_data(chr_name, em)
    draw_plot(chr_name, graph_epigi.number_of_nodes(), k_value_em, cc1_em, cc2_em, m_em, "Chromatin Loop (CL) model")

    k_value_pm, cc1_pm, cc2_pm, m_pm, density_pm = load_data(chr_name, pm)
    draw_plot(chr_name, graph_pm.number_of_nodes(), k_value_pm, cc1_pm, cc2_pm, m_pm, "Chromatin Anchor (CA) model")

    draw_plot_cmp(chr_name, graph_er.number_of_nodes(), k_value_er, cc1_er, k_value_ae, cc1_ae, k_value_te, cc1_te,
                  k_value_fe, cc1_fe, k_value_bm, cc1_bm, k_value_em, cc1_em, k_value_pm, cc1_pm, "largest cluster")
    draw_plot_cmp(chr_name, graph_ae.number_of_nodes(), k_value_er, cc2_er, k_value_ae, cc2_ae, k_value_te, cc2_te,
                  k_value_fe, cc2_fe, k_value_bm, cc2_bm, k_value_em, cc2_em, k_value_pm, cc2_pm,
                  "Second largest cluster")
    draw_plot_cmp(chr_name, graph_te.number_of_nodes(), k_value_er, m_er, k_value_ae, m_ae, k_value_te, m_te,
                  k_value_fe, m_fe, k_value_bm, m_bm, k_value_em, m_em, k_value_pm, m_pm, "moments")


# def run_compare_plot_part():
#     # start = 0
#     # end = int(len(list(graph.edges)) / 2) - 2000
#     # draw_plot_cmp_part(start, end, chr_name, graph.number_of_nodes(), k_value_er, cc1_er, k_value_ae, cc1_ae,
#     #                    k_value_te, cc1_te,
#     #                    "largest cluster")
#     # draw_plot_cmp_part(start, end, chr_name, graph.number_of_nodes(), k_value_er, cc2_er, k_value_ae, cc2_ae,
#     #                    k_value_te, cc2_te,
#     #                    "Second largest cluster")
#     # draw_plot_cmp_part(start, end, chr_name, graph.number_of_nodes(), k_value_er, m_er, k_value_ae, m_ae, k_value_te,
#     #                    m_te, "moments")


def creat_cyto_graph(chr_name):
    nodes = pd.read_csv(
        f'./data/compartment_nodes/chr1.csv',
        dtype={
            # 'Chromosome_no': 'str',  # for some files
            'Anchor_ID': 'str',
            'Compartment': 'str'
        }
    )
    edges = pd.read_csv(
        f'data/compartment_graph/rao_graph_without_linear_edges_chr1.csv',
        dtype={
            # 'Chromosome_no': 'str',  # for some files
            'Anchor_ID_A': 'str',
            'Anchor_ID_B': 'str',
            'Link_Weight': 'int',
            'Compartment': 'str'
        }
    )
    # file_cyto_out = open(f"results/{chr_name}/cytoscape_graph.csv", 'w')
    # for _, (a1, a2, w, comp) in edges.iterrows():
    #     comp_a1 = ''
    #     comp_a2 = ''
    #     flag = 0
    #     for _, (anchor, c) in nodes.iterrows():
    #         if a1 == anchor:
    #             comp_a1 = c
    #             flag += 1
    #         if a2 == anchor:
    #             comp_a2 = c
    #             flag += 1
    #         if flag == 2:
    #             break
    #     buffer = a1 + "," + comp_a1 + "," + a2 + "," + comp_a2 + "," + str(w) + "," + comp + "\n"
    #     file_cyto_out.write(buffer)
    #     # print(buffer)
    # file_cyto_out.close()


def read_json_file(filename):
    graph_read = nx.read_adjlist(filename)
    return graph_read


def create_graph_from_csv(csv_file, chromosome, min_weight=2, linear_edge_weight=25, min_distance=1):
    t0 = time.time()
    # Chromosome_no,Anchor_ID_A,Anchor_ID_B,Link_Weight
    # print("Input file : ", csv_file)
    # print("Chromosome : ", chromosome)
    # print("min_weight : ", min_weight)
    # print("min_distance : ", min_distance)
    raw_df = pd.read_csv(
        csv_file,
        dtype={
            # 'Chromosome_no': 'str',  # for some files
            'Anchor_ID_A': 'str',
            'Anchor_ID_B': 'str',
            'Link_Weight': 'int'
        }
    )

    def _parse_anchor(anchor):
        chr_, rest = anchor.split(':')
        start, end = rest.split('-')
        return chr_, int(start), int(end)

    edges = []
    nodes = set()
    for _, (a1, a2, w) in raw_df.iterrows():
        if w < min_weight:
            continue
        node1 = _parse_anchor(a1)  # (chr, start, end)
        node2 = _parse_anchor(a2)
        if not (node1[0] == chromosome and node2[0] == chromosome):
            continue
        if node1 > node2:
            node1, node2 = node2, node1
        if node2[1] - node1[1] >= min_distance:
            # print(node1, node2, w)
            edges.append((node1, node2, w))
            nodes.add(node1)
            nodes.add(node2)

    graph = nx.Graph()
    graph.add_weighted_edges_from(edges, kind='chiapet')

    if linear_edge_weight is not None:
        nodes = sorted(nodes)
        for v1, v2 in zip(nodes[:-1], nodes[1:]):
            if graph.has_edge(v1, v2):
                graph[v1][v2]['kind'] = 'both'
            else:
                graph.add_edge(v1, v2, kind='linear', weight=linear_edge_weight)

    print(
        f'Created network for {chromosome} with |V|={graph.number_of_nodes()}, |E|= {graph.number_of_edges()} ({time.time() - t0:.2f}s)')
    return graph


def multi_driver(chr_name):
    print(f"Creating output directories for {chr_name}")
    my_dir = "/mnt/raid/kaustavs/Desktop/Project_data/perculation/percolation/data/graph_test"
    # file_list = [f for f in os.listdir(my_dir)]
    # for f in file_list:
    #     os.remove(os.path.join(my_dir, f))
    # output_path_graphs = os.path.join(my_dir, chr_name)
    # os.mkdir(output_path_graphs)
    # os.mkdir(os.path.join(output_path_graphs, "graph_test_er"))
    # os.mkdir(os.path.join(output_path_graphs, "graph_test_ae"))
    # os.mkdir(os.path.join(output_path_graphs, "graph_test_te"))
    # os.mkdir(os.path.join(output_path_graphs, "graph_test_fe"))
    # os.mkdir(os.path.join(output_path_graphs, "graph_test_compertment"))
    # os.mkdir(os.path.join("./plots", chr_name))
    # os.mkdir(os.path.join("./results", chr_name))
    data_file = 'data/GM12878.CTCF.clusters_all_edges.csv'

    print("Creating Graph")
    graph = create_graph_from_csv(data_file, chr_name)
    graph_bio = create_bio_graph(data_file, chr_name)
    graph_epigi = create_epig_graph(data_file, chr_name)
    # graph = 0
    # graph_bio = 0
    # graph_epigi = 0

    # for (a, b) in graph.edges():
    #     get_prob_edge(graph, a, b)
    #     break

    print("Saving Graph")
    print("No of nodes in original graph :", len(list(graph.nodes())))
    print("No of edges in original graph :", len(list(graph.edges())))
    save_graph(graph, f"original", chr_name)
    save_graph(graph_bio, f"original_bio_graph", chr_name)

    print("Lodding Original Graphs")
    # graph = read_json_file(f"results/{chr_name}/graph_original.adjlist")
    # graph_bio = read_json_file(f"results/{chr_name}/graph_original_bio_graph_{chr_name}.adjlist")

    print("Drawing Degree Histogram")
    draw_degree_histogram(graph, chr_name, "Random")
    draw_degree_histogram(graph_bio, chr_name, "Compartment")
    draw_degree_histogram(graph_bio, chr_name, "Epigenomic")
    # creat_cyto_graph(chr_name)

    print("Calculate scalar values")
    scalar_value_calculate_anchors(chr_name)
    scalar_value_calculate_loops(chr_name)


    print("Running Percolation on Graph")
    # run_percolation(chr_name, graph, graph_bio, graph_epigi)
    #
    print("Drawing Simple Plots and Comparative Plots")

    graph_original = read_json_file(f"results/{chr_name}/graph_original.adjlist")
    print("No of nodes in original graph :", len(list(graph_original.nodes())))
    print("No of edges in original graph :", len(list(graph_original.edges())))

    graph_per_er = read_json_file(f"results/{chr_name}/graph_Erdos_Renyi.adjlist")
    print("No of nodes in ER :", len(list(graph_per_er.nodes())))
    print("No of edges in ER :", len(list(graph_per_er.edges())))

    graph_per_ae = read_json_file(f"results/{chr_name}/graph_Adjacent_Edges.adjlist")
    print("No of nodes in AE :", len(list(graph_per_ae.nodes())))
    print("No of edges in AE :", len(list(graph_per_ae.edges())))

    graph_per_te = read_json_file(f"results/{chr_name}/graph_Triangular_Edges.adjlist")
    print("No of nodes in TR :", len(list(graph_per_te.nodes())))
    print("No of edges in TR :", len(list(graph_per_te.edges())))

    graph_per_fe = read_json_file(f"results/{chr_name}/graph_Frequency_Edges.adjlist")
    print("No of nodes in FE :", len(list(graph_per_fe.nodes())))
    print("No of edges in FE :", len(list(graph_per_fe.edges())))

    graph_per_bio = read_json_file(f"results/{chr_name}/graph_Biological_Model.adjlist")
    print("No of nodes in CM :", len(list(graph_per_fe.nodes())))
    print("No of edges in CM :", len(list(graph_per_fe.edges())))

    graph_per_epigi = read_json_file(f"results/{chr_name}/graph_Epigenomic_Model.adjlist")
    print("No of nodes in EM :", len(list(graph_per_epigi.nodes())))
    print("No of edges in EM :", len(list(graph_per_epigi.edges())))

    graph_per_predictive = read_json_file(f"results/{chr_name}/graph_Predictive_Model.adjlist")
    print("No of nodes in PM :", len(list(graph_per_predictive.nodes())))
    print("No of edges in PM :", len(list(graph_per_predictive.edges())))

    run_draw_plot(chr_name, 'Erdos_Renyi', 'Adjacent_Edges', 'Triangular_Edges', 'Frequency_Edges', 'Biological_Model',
                  'Epigenomic_Model', 'Predictive_Model',
                  graph_per_er, graph_per_ae, graph_per_te, graph_per_fe, graph_per_bio, graph_per_epigi,
                  graph_per_predictive)


GENOME = [
    f'chr{str(s)}' for s in list(range(1, 22 + 1))
]  # we omit Y

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Percolation experiments")
    parser.add_argument('chromosomes', nargs='+', help="genome file")
    args = parser.parse_args()

    if 'all' in args.chromosomes:
        to_do = list(GENOME)
    else:
        to_do = [s for s in set(args.chromosomes) if s in GENOME]

    print('Jobs: ' + ' '.join(to_do))

    if len(to_do) == 1:
        multi_driver(to_do[0])
    else:
        p = Pool()
        res = p.map(multi_driver, to_do)

    print('ALL DONE')
