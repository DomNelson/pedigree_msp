import sys
sys.path.insert(0, "../msprime")
import os
import msprime
import numpy as np
import scipy.sparse
from collections import defaultdict
from itertools import combinations
import pandas as pd
from tqdm import tqdm
import argparse
from datetime import datetime
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # NOQA


class TreeIBD:
    def __init__(self, max_time, ts=None, ts_file=None, split_chroms=True):
        if ts is None and ts_file is None:
            print("One of ts or ts_file must be specified")
            raise ValueError

        if ts is None:
            ts = msprime.load(ts_file)
        self.ts = ts

        self.max_time = max_time
        self.bps = list(ts.breakpoints())
        self.node_times = self.get_node_times()
        self.ca_times = scipy.sparse.lil_matrix(
                (ts.num_samples, ts.num_samples))
        self.ca_last = scipy.sparse.lil_matrix(
                (ts.num_samples, ts.num_samples))
        self.ca_count = scipy.sparse.lil_matrix(
                (ts.num_samples, ts.num_samples))
        self.ibd_list = []

        self.recombination_map = None
        if split_chroms:
            self.recombination_map = get_WG_recombination_map()

    def get_node_times(self):
        node_times = {}
        for node in self.ts.nodes():
            node_times[node.id] = node.time

        node_times[-1] = np.inf

        return node_times

    def get_node_diff_candidates(self, diff):
        """
        Returns nodes which may have been added or removed
        in the diff - will be involved in the most edge diffs
        """
        counts = defaultdict(int)
        for edge in diff:
            counts[edge.parent] += 1
            counts[edge.child] += 1

        max_count = max(counts.values())
        out_candidates = [node for node, count in counts.items()
                          if count == max_count]

        return out_candidates

    def get_first_ancestor_below_max_time(self, tree, node):
        if self.node_times[node] > self.max_time:
            return -1

        parent = tree.get_parent(node)
        while parent >= 0:
            parent_time = self.node_times[parent]
            if parent_time > self.max_time:
                break
            node = parent
            parent = tree.get_parent(node)

        assert self.node_times[node] <= self.max_time
        return node

    def get_samples_below_edge_diff(self, tree, edge_diff):
        nodes = set()
        samples = set()
        segment, edges_in, edges_out = edge_diff
        assert type(edges_in) == list

        for edge in edges_in + edges_out:
            oldest_relevant_anc = self.get_first_ancestor_below_max_time(
                    tree,
                    edge.child
                    )
            if oldest_relevant_anc >= 0:
                nodes.add(oldest_relevant_anc)

        for node in nodes:
            samples.update(tree.get_leaves(node))

        return list(samples)

    def get_tree_min_common_ancestor_times(self, tree, edge_diff):
        nodes_in_diff = self.get_samples_below_edge_diff(tree, edge_diff)

        for a, b in combinations(nodes_in_diff, 2):
            if a == b:
                continue
            a, b = sorted([a, b])
            ca = tree.get_mrca(a, b)
            ca_time = self.node_times[ca]
            if ca_time > self.max_time:
                continue
            stored_time = self.ca_times[a, b]
            if stored_time == 0 or ca_time < stored_time:
                # If stored_time == 0 we're recording the first common ancestor
                self.ca_times[a, b] = ca_time
                self.ca_last[a, b] = ca
                self.ca_count[a, b] = 1

            elif ca_time == stored_time and self.ca_last[a, b] != ca:
                # Here we've found a second common ancestor from the same
                # generation as the first
                self.ca_count[a, b] = 2

    def get_all_min_common_ancestor_times(self):
        with tqdm(total=self.ts.num_trees) as pbar:
            for tree, edge_diff in zip(self.ts.trees(), self.ts.edge_diffs()):
                self.get_tree_min_common_ancestor_times(tree, edge_diff)
                pbar.update(1)

    def write_ibd_pair(self, ibd_start, ibd_end, pair):
        if ibd_start[pair] != 0:
            assert ibd_end[pair] != 0
            ind1, ind2 = pair

            # Here we undo the offset we used when building segments
            start_idx = ibd_start[pair] - 1
            end_idx = ibd_end[pair] - 1

            start = self.bps[start_idx]
            end = self.bps[end_idx]

            if self.recombination_map:
                # Splits IBD segment on chromosome boundaries, updating
                # start position so the last split is added as normal
                for p in self.recombination_map.get_positions():
                    if start < p < end:
                        record = [ind1, ind2, start, p]
                        self.ibd_list.append(record)
                        start = p

            record = [ind1, ind2, start, end]
            self.ibd_list.append(record)

            # Reset the pair for building the next segment
            ibd_start[pair] = 0
            ibd_end[pair] = 0

    def write_ibd_all(self, ibd_start, ibd_end):
        coo = ibd_start.tocoo()
        for ind1, ind2, start in zip(coo.row, coo.col, coo.data):
            pair = (ind1, ind2)
            self.write_ibd_pair(ibd_start, ibd_end, pair)

    def get_children_of_edges_crossing_max_time(self, locus=None):
        children = set()
        for edge in self.ts.edges():
            if (locus is not None and
                    (edge.left > locus or edge.right <= locus)):
                continue
            parent_time = self.node_times[edge.parent]
            child_time = self.node_times[edge.child]
            if parent_time > self.max_time and child_time <= self.max_time:
                # print("Adding", edge)
                children.add(edge.child)

        return children

    def get_initial_clusters(self, first_tree):
        ancs = self.get_children_of_edges_crossing_max_time(locus=0)
        clusters = defaultdict(set)

        # If no clusters cross max_time, all nodes must be below max_time,
        # so samples cluster by root node (everyone is IBD unless the
        # tree is not fully coalesced and has multiple roots)
        if len(ancs) == 0:
            ancs = first_tree.roots

        for a in ancs:
            clusters[a] = set(first_tree.get_leaves(a))

        return clusters

    def update_clusters(self, tree, diff, clusters):
        new_clusters = defaultdict(set)
        changed_samples = self.get_samples_below_edge_diff(tree, diff)
        # print("Changed samples:", len(changed_samples))
        # print("Num clusters:", len(clusters))

        # First remove changed samples from existing clusters and update
        # oldest anc if necessary
        for anc, cluster in clusters.items():
            if len(cluster) == 0:
                continue
            new_anc = self.get_first_ancestor_below_max_time(tree, anc)
            new_clusters[new_anc].update(cluster.difference(changed_samples))

        # Now add to new clusters
        for sample in changed_samples:
            anc = self.get_first_ancestor_below_max_time(tree, sample)
            new_clusters[anc].add(sample)

        # Sanity check - no clusters should share samples
        # for c1, c2 in combinations(new_clusters.values(), 2):
        #     assert len(c1.intersection(c2)) == 0

        return new_clusters

    def get_ibd(self):
        n = self.ts.num_samples
        ibd_start = scipy.sparse.lil_matrix((n, n), dtype=int)
        ibd_end = scipy.sparse.lil_matrix((n, n), dtype=int)

        # Initialize with the first tree
        trees = self.ts.trees()
        diffs = self.ts.edge_diffs()
        first_tree = next(trees)
        # first_diff = next(diffs)
        clusters = self.get_initial_clusters(first_tree)

        # NOTE 1: We compute the upper-triangular portion of the IBD
        # matrix only. If a < b, IBD[b, a] == 0 does NOT imply no IBD
        # shared, only IBD[a, b] is meaningful.
        #
        # NOTE 2: Indices are offset +1 so that sparse matrix default fill
        # of zero implies no IBD, not IBD starting at 0 in genomic coords.
        # IBD starting at 0 is denoted by index 1.
        with tqdm(total=self.ts.num_trees, desc="Writing IBD pairs") as pbar:
            for i, (tree, diff) in enumerate(zip(trees, diffs)):
                # if i % 1000 == 0:
                #     print("Num clusters:", len(clusters))
                for cluster in clusters.values():
                    ibd_pairs = combinations(sorted(cluster), 2)
                    for pair in ibd_pairs:
                        # Check if we are starting a new IBD segment
                        # or continuing one
                        if ibd_end[pair] == i + 1:
                            ibd_end[pair] += 1
                        else:
                            # If we start a new segment, write the old
                            # one first
                            self.write_ibd_pair(ibd_start, ibd_end, pair)
                            ibd_start[pair] = i + 1
                            ibd_end[pair] = i + 2
                clusters = self.update_clusters(tree, diff, clusters)
                pbar.update(1)

        # TODO: Add last tree - integrate into main loop above
            if self.ts.num_trees > 1:
                i += 1
                for cluster in clusters.values():
                    ibd_pairs = combinations(sorted(cluster), 2)
                    for pair in ibd_pairs:
                        # Check if we are starting a new IBD segment
                        # or continuing one
                        if ibd_end[pair] == i + 1:
                            ibd_end[pair] += 1
                        else:
                            # If we start a new segment, write the
                            # old one first
                            self.write_ibd_pair(ibd_start, ibd_end, pair)
                            ibd_start[pair] = i + 1
                            ibd_end[pair] = i + 2
                pbar.update(1)

        # Write out all remaining segments, which reached the end of the
        # simulated region
        self.write_ibd_all(ibd_start, ibd_end)


def get_positions_rates(chrom_lengths, rho):
    """
    Takes a list of chromosome lengths and returns lists of positions and
    rates to pass to msprime.RecombinationMap
    """
    positions = []
    rates = []
    total_length = 0
    for length in chrom_lengths:
        positions.extend([int(total_length),
                          int(total_length) + int(length) - 1])
        rates.extend([rho, 0.5])
        total_length += length

    rates[-1] = 0

    return positions, rates


def plot_ibd(ibd_list, ca_times=None, min_length=0, out_file=None):
    print("Plotting!")
    fig, ax = plt.subplots()
    if out_file is None:
        out_file = '~/temp/ibd_plot.png'
    out_file = os.path.expanduser(out_file)

    cols = ["ind1", "ind2", "start", "end"]

    ibd_df = pd.DataFrame(ibd_list, columns=cols)
    ibd_df['len'] = ibd_df['end'] - ibd_df['start']
    ibd_df = ibd_df[ibd_df['len'] >= min_length]

    ibd_df['count'] = 1
    inds = set(ibd_df['ind1'].values)
    for ind in inds:
        ind_df = ibd_df[ibd_df['ind1'] == ind]
        ind_pairwise_ibd_df = ind_df.groupby('ind2').sum()
        ind_pairwise_ibd_df = ind_pairwise_ibd_df.reset_index()

        total_IBD = ind_pairwise_ibd_df['len'].values
        num_segments = ind_pairwise_ibd_df['count'].values
        ind2 = ind_pairwise_ibd_df['ind2'].values
        sizes = np.ones(total_IBD.shape[0]) * 8

        colours = None
        cmap = None
        if ca_times is not None:
            colours = []
            for i in range(len(ind2)):
                x, y = sorted([ind, ind2[i]])
                assert(x <= y)
                colours.append(ca_times[x, y])
                if colours[-1] == 0:
                    print(x, y)

            colours = np.array(colours)
            cmap = 'viridis'

        ax.scatter(total_IBD, num_segments, s=sizes, c=colours, cmap=cmap)

    ax.set_ylabel('Number of IBD segments')
    ax.set_xlabel('Total IBD')
    ax.set_xscale('log')
    print("Saving...")
    fig.savefig(out_file)
    print("Done!")


def get_WG_recombination_map():
    rho = 1e-8
    chrom_lengths_morgans = [2.77693825, 2.633496065, 2.24483368,
                             2.12778391, 2.03765845, 1.929517394,
                             1.867959329, 1.701765192, 1.68073935,
                             1.789473882, 1.594854258, 1.72777271,
                             1.26940475, 1.16331251, 1.2554709,
                             1.348911043, 1.29292106, 1.18978483,
                             1.077960694, 1.079243479, 0.61526812,
                             0.72706815]
    chrom_lengths = [l * 1e8 for l in chrom_lengths_morgans]
    num_loci = int(chrom_lengths[-1] + 1)

    positions, rates = get_positions_rates(chrom_lengths, rho)
    recombination_map = msprime.RecombinationMap(
            positions, rates, num_loci=num_loci
    )

    return recombination_map


def simulate(Ne, sample_size, model, max_time):
    max_time = max_time + 1e-9  # To include 'max_time' in node times
    recombination_map = get_WG_recombination_map()
    population_configuration = msprime.PopulationConfiguration(
            sample_size=sample_size, initial_size=Ne)

    ts = msprime.simulate(
            population_configurations=[population_configuration],
            model=model,
            recombination_map=recombination_map,
            end_time=max_time
            )

    return ts


def draw_trees(ts):
    for t in ts.trees():
        print(t.draw(format='unicode'))


def build_fname(label, ext, timestamp, args):
    fname = timestamp + '_'

    if args.Ne and args.sample_size and args.model:
        fname += 'Ne' + str(args.Ne) + '_'
        fname += 'samplesize' + str(args.sample_size) + '_'
        fname += 'maxtime' + str(args.max_time) + '_'
        fname += 'model' + str(args.model).upper() + '_'

    fname += label + '.' + ext
    dirname = os.path.expanduser(args.output_dir)

    return os.path.join(dirname, fname)


def check_args(args):
    if args.ts_file is None and (args.Ne is None or args.sample_size is None
                                 or args.model is None):
        raise ValueError("Must specify either tree sequence file or " +
                         "Ne, sample size, and model type " +
                         "('dtwf' or 'hudson').")

    if args.ts_file and (args.Ne or args.sample_size or args.model):
        raise ValueError("Cannot specify both tree sequence file and " +
                         "simulation parameters")

    if args.plot is None and args.output_dir is None:
        raise ValueError("Must specify at least one type of output!")


def main(args):
    check_args(args)
    timestamp = datetime.strftime(datetime.now(), '%Y-%m-%d_%H-%M-%S')

    ts_file = None
    if args.ts_file:
        ts_file = os.path.expanduser(args.ts_file)

    ts = None
    if args.Ne and args.sample_size and args.model:
        ts = simulate(args.Ne, args.sample_size, args.model,
                      max_time=args.max_time)
        if args.output_dir:
            ts_out_file = build_fname('ts', 'h5', timestamp, args)
            ts.dump(ts_out_file)

    t_ibd = TreeIBD(args.max_time, ts_file=ts_file, ts=ts,
                    split_chroms=args.split_chroms)
    t_ibd.get_ibd()
    t_ibd.get_all_min_common_ancestor_times()

    if args.output_dir:
        ibd_out_file = build_fname('ibd', 'npz', timestamp, args)
        ibd_array = np.array(t_ibd.ibd_list)
        np.savez_compressed(ibd_out_file, ibd_array=ibd_array)

        ca_out_file = build_fname('ca_times', 'npz', timestamp, args)
        scipy.sparse.save_npz(ca_out_file, t_ibd.ca_times.tocoo())

    if args.plot:
        plot_out_file = build_fname('plot_ibd', 'png', timestamp, args)
        plot_ibd(t_ibd.ibd_list, t_ibd.ca_times, min_length=args.min_length,
                 out_file=plot_out_file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--ts_file')
    parser.add_argument('--Ne', type=int)
    parser.add_argument('--sample_size', type=int)
    parser.add_argument('--model', choices=['dtwf', 'hudson'])
    parser.add_argument('--max_time', type=int, default=5)
    parser.add_argument('--min_length', type=float, default=5e6)
    parser.add_argument('--output_dir')
    parser.add_argument('--split_chroms', action='store_true')
    parser.add_argument('--plot', action='store_true')

    args = parser.parse_args()

    main(args)
