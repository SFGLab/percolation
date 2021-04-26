import argparse
import time
from collections import OrderedDict
from multiprocessing import Pool
from typing import Tuple, Type, Dict, List
from graph_datatypes import TimingMessage

import networkx as nx
import numpy as np
import pandas as pd
import dask.dataframe

from numpy.random import default_rng
from graph_datatypes import DisjointSetPlus


timing_info = TimingMessage


def graph_edges_from_networkx(graph: nx.Graph, weights='weight'):
    node_to_idx = {u: i for i, u in enumerate(graph.nodes())}
    dtype = [('u', int), ('v', int), ('w', int), ('kind', bool)]
    edges = np.empty(graph.number_of_edges(), dtype=dtype)
    for k, (u, v) in enumerate(graph.edges()):
        i, j = node_to_idx[u], node_to_idx[v]
        if i > j:
            i, j = j, i
        w = graph[u][v]['weight']
        kind = graph[u][v]['kind'] == 'linear'
        edges[k] = i, j, w, kind
    return edges, node_to_idx


class PercolationMethod(object):
    name = 'XX'

    def __init__(self, graph: nx.Graph, edges: np.ndarray, node_to_idx: Dict[str, int], rng):
        self._graph = graph
        self._node_to_idx = node_to_idx
        self._rng = rng
        self.edges = edges # TODO: field exposed for performance reasons

    @property
    def rng(self):
        return self._rng

    @property
    def graph(self):
        return self._graph

    def prepare(self):
        pass

    def get_next_edge(self, i_step: int) -> Tuple[int, int, int, int]:
        raise NotImplemented()


class ErdosRenyiMethod(PercolationMethod):
    name = 'ER'

    def __init__(self, graph: nx.Graph, edges: np.ndarray, node_to_idx: Dict[str, int], rng):
        super(ErdosRenyiMethod, self).__init__(graph, edges.copy(), node_to_idx, rng)

    def prepare(self):
        self.rng.shuffle(self.edges, axis=0)

    def get_next_edge(self, i_step: int) -> Tuple[int, int]:
        return self.edges[i_step]


class FrequencyEdgeMethod(PercolationMethod):
    name = 'FE'

    def __init__(self, graph: nx.Graph, edges: np.ndarray, node_to_idx: Dict[str, int], rng):
        edges_sorted = edges.copy()
        edges_sorted.sort(axis=0, order='w')
        edges_sorted = edges_sorted[::-1]
        super(FrequencyEdgeMethod, self).__init__(graph, edges_sorted, node_to_idx, rng)

    def prepare(self):
        pass

    def get_next_edge(self, i_step: int) -> Tuple[int, int]:
        return self.edges[i_step]


# class FrequencyEdgeSoftMethod(PercolationMethod):
#     name = 'FES'
#
#     def __init__(self, graph: nx.Graph, edges: np.ndarray, node_to_idx: Dict[str, int], rng):
#         super(FrequencyEdgeSoftMethod, self).__init__(graph, edges, node_to_idx, rng)
#         edges_shuffled = edges.copy()
#         weights = edges_shuffled['w']
#         prob = weights / weights.sum()
#         idx = rng.choice(len(edges_shuffled), replace=False, p=prob)
#         self._edges_shuffled = edges_shuffled[idx]
#
#     def get_next_edge(self, i_step: int) -> Tuple[int, int]:
#         u, v, *_ = self._edges_shuffled[i_step]
#         return u, v


def run_percolation(
        clusters: DisjointSetPlus,
        cc1_sizes: np.ndarray, cc2_sizes: np.ndarray, sizes_sq_sums: np.ndarray, n_clusters: np.ndarray, is_strand: np.ndarray,
        method: PercolationMethod
):
    data_size = len(cc1_sizes)
    n_nodes = clusters.n_clusters
    sizes_sq_sums_current = n_nodes  # n_nodes of size 1

    for i_step in range(data_size):
        a, b, _, kind = method.get_next_edge(i_step)
        size_a = clusters.get_size(a)
        size_b = clusters.get_size(b)
        if clusters.union(a, b):
            sizes_sq_sums_current += 2 * size_a * size_b  # simple school algebra:)
        sizes_sq_sums[i_step] = sizes_sq_sums_current
        cc1_sizes[i_step] = clusters.cc1_size
        cc2_sizes[i_step] = clusters.cc2_size
        n_clusters[i_step] = clusters.n_clusters
        is_strand[i_step] = kind


def run_percolation_batch(graph, edges, node_to_idx, method_cls, n_replicates):
    data_size = graph.number_of_edges()
    n_nodes = graph.number_of_nodes()
    cc1_sizes = np.empty((n_replicates, data_size), dtype=int)
    cc2_sizes = np.empty((n_replicates, data_size), dtype=int)
    sizes_sq_sums = np.empty((n_replicates, data_size), dtype=int)
    n_clusters = np.empty((n_replicates, data_size), dtype=int)
    is_strand = np.empty((n_replicates, data_size), dtype=bool)

    rng = default_rng()
    track_cc = 2
    clusters = DisjointSetPlus(n_nodes, track_cc=track_cc)
    method = method_cls(graph, edges, node_to_idx, rng)

    for i in range(n_replicates):
        method.prepare()
        clusters.reset(track_cc)
        run_percolation(clusters, cc1_sizes[i], cc2_sizes[i], sizes_sq_sums[i], n_clusters[i], is_strand[i], method)

    return {
        'cc1_sizes': cc1_sizes,
        'cc2_sizes': cc2_sizes,
        'sizes_sq_sums': sizes_sq_sums,
        'n_clusters': n_clusters,
        'is_strand': is_strand
    }


def get_raw_data_from_csv(csv_file):
    # Chromosome_no,Anchor_ID_A,Anchor_ID_B,Link_Weight
    raw_df = pd.read_csv(
        csv_file,
        dtype={
            # 'Chromosome_no': 'str',  # for some files
            'Anchor_ID_A': 'str',
            'Anchor_ID_B': 'str',
            'Link_Weight': 'int'
        }
    )
    return raw_df


def create_chromosome_graph(raw_df, chromosome, min_weight=2, linear_edge_weight=1, min_distance=1):
    t0 = time.time()

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

    return graph


def save_percolation_results(output_file, results, i_repl_start, i_repl_end, n_steps, method_name, chromosome):
    n_replicates = i_repl_end - i_repl_start
    data_size = n_replicates * n_steps
    cols = [
           ('chromosome', np.repeat(np.array(chromosome), data_size)),
           ('method', np.repeat(np.array(method_name), data_size)),
           ('i_repl', np.repeat(np.arange(i_repl_start, i_repl_end), n_steps)),
           ('i_step', np.tile(np.arange(n_steps), n_replicates))
    ] + [
        (name, data.reshape(data_size)) for name, data in results.items()
    ]
    df = pd.DataFrame.from_dict(OrderedDict(cols))
    ddf = dask.dataframe.from_pandas(df, npartitions=1).compute()
    ddf.to_hdf(output_file, 'raw', mode='w', format='table', complevel=5)  # compression='gzip', compression_opts=5


METHODS = [
    ErdosRenyiMethod
    # FrequencyEdgeSoftMethod
]

MULTIPROCESSING_SHARED_DATA = {}


def setup_shared_input_data(chiapet_file, chromosomes):
    for chromosome in chromosomes:
        with TimingMessage(f'Loading {chiapet_file}') as ti:
            raw_data = get_raw_data_from_csv(chiapet_file)
            ti.finished(f'Loaded {len(raw_data)} rows.')
        with TimingMessage(f'Creating graph for {chromosome}') as ti:
            graph = create_chromosome_graph(raw_data, chromosome)
            edges, node_to_idx = graph_edges_from_networkx(graph)
            ti.finished(f'Graph size: |V|={graph.number_of_nodes()}, |E|={graph.number_of_edges()}')
        MULTIPROCESSING_SHARED_DATA[chromosome + '_graph'] = graph
        MULTIPROCESSING_SHARED_DATA[chromosome + '_edges'] = edges
        MULTIPROCESSING_SHARED_DATA[chromosome + '_node_to_idx'] = node_to_idx


def run_job(job_id, chromosome, method, i_repl_start, i_repl_end, output_file, save_trajectories):
    graph = MULTIPROCESSING_SHARED_DATA[chromosome + '_graph']
    edges = MULTIPROCESSING_SHARED_DATA[chromosome + '_edges']
    node_to_idx = MULTIPROCESSING_SHARED_DATA[chromosome + '_node_to_idx']
    n_replicates = i_repl_end - i_repl_start
    with TimingMessage('Simulation ' + job_id):
        results = run_percolation_batch(graph, edges, node_to_idx, method, n_replicates)
    if save_trajectories:
        with TimingMessage('Saving ' + job_id):
            save_percolation_results(output_file, results, i_repl_start, i_repl_end, len(edges), method.name, chromosome)


def make_jobs(args) -> List[Tuple[int, str, PercolationMethod, int, int]]:
    n_replicates = args.nreplicates
    batch_size = args.batch_size if args.batch_size > 0 else n_replicates

    n_full_batches = args.nreplicates // args.batch_size
    spare_batch = args.nreplicates % args.batch_size
    # TODO: warning, ignores spare batch!

    setup_shared_input_data(args.chiapet_file, args.chromosomes)

    def _mk_job(i_batch, chromosome, method, start, end, save_trajectories):
        return (
            f'{chromosome}-{method.name}-{i_batch}',
            chromosome, method, start, end,
            f'{args.output_dir}/perc_res_{chromosome}_{method.name}_{i_batch:03d}.h5', save_trajectories
        )

    jobs = [
        _mk_job(i, chrom, met, rep, rep + batch_size, args.trajectories)
        for chrom in args.chromosomes
        for met in METHODS
        for i, rep in enumerate(range(0, n_replicates, batch_size))
    ]
    jobs += [
        _mk_job(1, chrom, met, 0, 1, args.trajectories)
        for chrom in args.chromosomes
        for met in [FrequencyEdgeMethod]
    ]
    return jobs


def main(args):
    n_cores = None if args.ncores == 0 else args.ncores
    if not args.trajectories:
        print('NOT saving trajectories.')
    else:
        print(f'Saving trajectories to {args.output_dir}')
    jobs = make_jobs(args)

    with timing_info('All jobs') as ti:
        ti.print(f'Running {len(jobs)} jobs.')
        if n_cores == 1 or len(jobs) == 1:
            ti.print('Not using parallel execution.')
            for job in jobs:
                run_job(*job)
        else:
            ti.print(f'Using {n_cores} cores.')
            with Pool(n_cores) as pool:
                pool.starmap(run_job, jobs)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Percolation experiments")
    parser.add_argument('chiapet_file', help="Input chiapet file")
    parser.add_argument('output_dir', help="Output dir")
    parser.add_argument('chromosomes', nargs='+', help="Chromosomes to run")
    parser.add_argument('-n', '--nreplicates', type=int, default=10)
    parser.add_argument('-b', '--batch_size', type=int, default=0)
    parser.add_argument('--trajectories', action='store_true')
    parser.add_argument('--ncores', type=int, default=0)
    main(parser.parse_args())
