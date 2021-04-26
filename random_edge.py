import argparse
import random as rd
import re
import time
from heapq import nlargest
from multiprocessing import Pool

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd


def draw_plot(chr_name, k_values, cc1_sizes, cc2_sizes, moments_2, plots_dir=None):
    n = int(cc1_sizes.max())
    fig, ax = plt.subplots(figsize=(6, 4), dpi=150)
    k_norm = k_values / k_values.max()
    ax.plot(k_norm, cc1_sizes / n, c='C0', label='C1/n')
    ax.plot(k_norm, cc2_sizes / n, c='C1', label='C2/n')
    ax.plot(k_norm, moments_2 / moments_2.max(), c='C2', label='m2 of |CC| (norm.)')
    ax.set(xlabel='ratio of edges added', title=f'{chr_name} CTCF 2+ graph (|V|={n}, |E|={len(k_values)})')
    ax.legend(loc='upper left')
    ax.grid()
    if plots_dir is not None:
        fig.savefig(f"{plots_dir}/{chr_name} CTCF 2+ Random Edges Plot.png")
        # fig.savefig(f"{plots_dir}/{chr_name} CTCF 2+ Random Edges Plot.pdf")
    plt.show()

    fig, ax = plt.subplots(figsize=(6, 4), dpi=150)
    ax.plot(k_norm, cc1_sizes / n, c='C0', label='C1/n')
    ax.plot(k_norm, cc2_sizes / n, c='C1', label='C2/n')
    ax.plot(k_norm, moments_2 / moments_2.max(), c='C2', label='m2 of |CC| (norm.)')
    ax.set(xlabel='ratio of edges added', title=f'{chr_name} CTCF 2+ graph (|V|={n}, |E|={len(k_values)})')
    ax.set_xlim(0.4, 0.6)
    ax.set_ylim(0.0, 0.2)
    ax.legend(loc='upper left')
    ax.grid()
    if plots_dir is not None:
        fig.savefig(f"{plots_dir}/{chr_name} CTCF 2+ Random Edges Plot zoomed.png")
        # fig.savefig(f"{plots_dir}/{chr_name} CTCF 2+ Random Edges Plot zoomed.pdf")
    plt.show()


def sorting(l):
    """

    :type l: object
    """
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key=alphanum_key)


def run_percolation(graph):
    edges_shuffled = list(graph.edges())
    rd.shuffle(edges_shuffled)

    roots = list(graph.nodes())
    uf = nx.utils.union_find.UnionFind(roots)

    data_size = len(edges_shuffled)
    k_values = np.empty(data_size)
    cc1_sizes = np.empty(data_size)
    cc2_sizes = np.empty(data_size)
    moments_2 = np.empty(data_size)
    # cluster_sizes = np.zeros((data_size, data_size), dtype='int32')

    for k, (a, b) in enumerate(edges_shuffled):
        graph.add_edge(a, b)
        root_a = uf[a]
        root_b = uf[b]
        if root_a != root_b:
            uf.union(a, b)
            root_ab = uf[a]
            roots.remove(root_a)
            roots.remove(root_b)
            roots.append(root_ab)
        sizes = [uf.weights[v] for v in roots]
        # cluster_sizes[k, :len(sizes)] = sizes
        unique_sizes, counts = np.unique(sizes, return_counts=True)
        moment_2 = ((unique_sizes ** 2) * counts).sum()
        if len(roots) > 1:
            cc1_size, cc2_size = nlargest(2, sizes)
        else:
            cc1_size, cc2_size = uf.weights[roots[0]], 0.0
        cc1_sizes[k] = cc1_size
        cc2_sizes[k] = cc2_size
        k_values[k] = k + 1
        moments_2[k] = moment_2

    # return cluster_sizes, k_values, cc1_sizes, cc2_sizes, moments_2
    return k_values, cc1_sizes, cc2_sizes, moments_2


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
    file_whole_genome = open("data/GM12878.CTCF.clusters_all_edges.csv")
    edges = []
    for line in file_whole_genome:
        line = line.split("\n")[0]
        chr_name = line.split(",")[0].split(":")[0]
        if chr_name == chr_n:
            # print(line.split(",")[0] + "," + line.split(",")[1] + "," + line.split(",")[2])
            edges.append(line.split(",")[0] + "," + line.split(",")[1] + "," + line.split(",")[2])
    file_whole_genome.close()
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


def filter_edges(edges, thereshold=100):
    filter_edge = []    
    
    def split_node(node):
        chr_, rest = node.split(':')
        node_start, node_end = rest.split('-')
        return chr_, int(node_start), int(node_end)
    
    for edge in edges:
        _, node_a, node_b = edge.split(",")
        # link_wt = int(link_wt)
        # if link_wt < 2:
        #    continue
        chr_a, node_start_a, node_end_a = split_node(node_a)
        chr_b, node_start_b, node_end_b = split_node(node_b)
        if node_end_a - node_start_a > thereshold and node_end_b - node_start_b > thereshold:
            filter_edge.append(edge)
    return filter_edge


def multi_driver(chr_mane):
    # edges = extract_data(chr_mane)
    # print("No of Edges in " + chr_mane + " : " + str(len(edges)))
    # ed = filter_edges(edges)
    # print("No of Edges after filtering in " + chr_mane + " : " + str(len(ed)))
    # edge = add_linear_edge(ed)
    # print("No of Edges after adding liner edges in " + chr_mane + " : " + str(len(edge)))
    # n, e = load_data(edge)
    # print("No of Nodes after Loading Data in " + chr_mane + " : " + str(len(n)))
    # print("No of Edges after Loading Data in " + chr_mane + " : " + str(len(e)))

    t0 = time.time()
    data_file = 'data/GM12878.CTCF.clusters_all_edges.csv'
    # data_file = 'data/GM12878.CTCF.chr8.long_interactions.csv'
    # graph = create_graph_from_csv(data_file, chr_mane)
    # print(f'Created network for {chr_mane} with |V|={graph.number_of_nodes()}, |E|= {graph.number_of_edges()} ({time.time() - t0:.2f}s)')
    # t0 = time.time()
    # sizes, *results = run_percolation(graph)
    # results = run_percolation(graph)
    # print(f'Simualtion completed for {chr_mane}  ({time.time() - t0:.2f}s)')
    # save_cluster_sizes(f'results/{chr_mane}_cluster_sizes.hdf5', sizes)
    results_filename = f'results/{chr_mane}_results.txt'
    # save_results(results_filename, results)
    X = np.loadtxt(results_filename)
    results = [
        X[:, i]
        for i in range(X.shape[1])
    ]
    # chr_name, n, k_values, cc1_sizes, cc2_sizes, moments_2, plots_dir = None
    draw_plot(chr_mane, *results, plots_dir='plots')


def create_graph_from_csv(csv_file, chromosome, min_weight=2, linear_edge_weight=1, min_distance=100):
    # Chromosome_no,Anchor_ID_A,Anchor_ID_B,Link_Weight
    raw_df = pd.read_csv(
        csv_file,
        dtype={
            # 'Chromosome_no': 'str', for some files
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
            edges.append((node1, node2, w))
            nodes.add(node1)
            nodes.add(node2)

    graph = nx.Graph()
    # graph.add_nodes_from(nodes)  # do not add singleton nodes
    graph.add_weighted_edges_from(edges, kind='chiapet')

    if linear_edge_weight is not None:
        nodes = sorted(nodes)
        for v1, v2 in zip(nodes[:-1], nodes[1:]):
            if graph.has_edge(v1, v2):
                graph[v1][v2]['kind'] = 'both'
            else:
                graph.add_edge(v1, v2, weight=linear_edge_weight, kind='linear')

    return graph


# def save_cluster_sizes(file_path, sizes):
#     with h5py.File(file_path, 'w') as hf:
#         hf.create_dataset('cluster_sizes', data=sizes)

def save_results(file_path, *results):
    np.savetxt(file_path, np.concatenate(results).T)


GENOME = [
    f'chr{str(s)}' for s in list(range(1, 22 + 1)) + ['X']
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
