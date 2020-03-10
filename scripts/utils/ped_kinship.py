import sys
sys.path.append('/home/dnelson/project/msprime')
import msprime
import networkx as nx
import itertools
import numpy as np
import matplotlib.pyplot as plt
import os
import pytest

from tqdm import tqdm
from collections import defaultdict
        
    
# def compute_ind_kinship(self, ind1, ind2):
#     if len(list(self.G.predecessors(ind1))) == 0:
#         if ind1 == ind2:
#             # Founder kinship with self
#             return 0.5
#         if len(list(self.G.predecessors(ind2))) == 0:
#             # Founder kinship with other founder
#             return 0

#     if ind1 == ind2:
#         parents = self.G.predecessors(ind1)
#         return 0.5 * (1 + self.compute_ind_kinship(*parents))
#     else:
#         parents = self.G.predecessors(ind1)
#         return 0.5 * sum(
#             [self.compute_ind_kinship(ind2, p) for p in parents])

# def compute_kinship_recursive(self, samples):
#     n = len(samples)
#     kinship = np.zeros((n, n))
#     sample_idx_dict = dict(zip(samples, range(n)))

#     for s1, s2 in itertools.combinations_with_replacement(samples, 2):
#         sample_idx1 = sample_idx_dict[s1]
#         sample_idx2 = sample_idx_dict[s2]
#         sorted_idx = sorted([sample_idx1, sample_idx2])
#         sorted_idx = tuple(sorted_idx)
#         kinship[sorted_idx] = self.compute_ind_kinship(s1, s2)

#     return kinship

def digraph_kinship(G, samples, kinship, max_time=-1):
    ancestors = defaultdict(dict)
    climbers = {s: set([s]) for s in samples}
    
    i = 0
    while len(climbers) > 0:
        # Adding one only works if branch lengths are also == 1
        # TODO: Make more flexible ^^
        i += 1
        next_climbers = defaultdict(set)
        for sample, climbers in climbers.items():
            for climber in climbers:
                for anc in G.predecessors(climber):
                    next_climbers[sample].add(anc)
                    ancestors[anc][sample] = i

                    # Add kinship with other descendants who
                    # have reached this ancestor
                    for descendant in ancestors[anc]:
                        # print("Adding kinship between", sample, "and", descendant)
                        j = ancestors[anc][descendant]
                        idx = tuple(sorted([sample, descendant]))
                        kinship[idx] += 1 / (2 ** (i + j))
                    
        climbers = next_climbers
                    
    return kinship 


def tree_to_dict_of_dicts(tree):
    dod = {}
    for parent in tree.nodes():
        dod[parent] = {}
        for child in tree.children(parent):
            dod[parent][child] = {
                'branch_length': tree.branch_length(child)
            }
            
    return dod


class PedKinship:
    def __init__(self, G):
        self.G = G
        self.samples = None
        self.ancestor_dists = None
        self.kinship = None
        
    @staticmethod
    def from_array(pedarray):
        G = PedKinship.ped_to_digraph(pedarray)
        
        return PedKinship(G)   
        
    @staticmethod
    def from_txt(pedfile):
        pedarray = np.genfromtxt(
            os.path.expanduser(pedfile),
            skip_header=True,
            dtype=int)
                
        return PedKinship.from_array(pedarray)
        
    @staticmethod
    def ped_to_digraph(pedarray):
        G = nx.DiGraph()

        for child, father, mother in pedarray[:, :3]:
            if father > 0:
                G.add_edge(father, child)
            if mother > 0:
                G.add_edge(mother, child)

        return G
    
    def get_probands(self):
        return [v for v, d in self.G.out_degree() if d == 0]
    
    def _set_samples(self, samples):
        self.samples = samples
        n = len(samples)
        self.kinship = np.zeros((n, n))
    
    def _build_ancestor_dists(self, max_depth=7):
        """
        Builds a dict of dicts. Outer dict links ancestors to probands,
        inner dict holds list of path lengths between the two.
        """
        if self.samples is None:
            raise ValueError("Samples not set")
            
        # Each ancestor has a dict holding a list of path lengths
        # keyed by child ID
        self.ancestor_dists = defaultdict(lambda: defaultdict(list))
        climbers = {x: set([s]) for x, s in enumerate(self.samples)}
        
        i = 0
        while len(climbers) > 0:
            if i > max_depth:
                break

            next_climbers = defaultdict(set)
            for sample_ix, climbers in climbers.items():
                for climber in climbers:
                    self.ancestor_dists[climber][sample_ix].append(i)
                    for anc in self.G.predecessors(climber):
                        next_climbers[sample_ix].add(anc)

            climbers = next_climbers
            # Adding one only works if branch lengths are also == 1
            # TODO: Make more flexible ^^
            i += 1

    def compute_kinship(self, samples, max_depth=7):
        self._set_samples(samples)
        self._build_ancestor_dists(max_depth=max_depth)

        for anc, dists in tqdm(self.ancestor_dists.items(),
                               total=len(self.ancestor_dists),
                               smoothing=0.1):
            children = list(self.G.successors(anc))
            children_dists = [self.ancestor_dists[c] for c in children
                              if c in self.ancestor_dists]

            for paths1, paths2 in itertools.combinations(children_dists, 2):
                for sample_idx1, lengths1 in paths1.items():
                    for sample_idx2, lengths2 in paths2.items():
                        for l1, l2 in itertools.product(lengths1, lengths2):
                            # Since we're looking at the path lengths to 
                            # children of the common ancestor, we add one
                            # to each length
                            kinship = 1 / (2 ** (1 + l1 + 1 + l2 + 1))

                            # Make sure we're only filling upper-triangular
                            # portion of matrix
                            sorted_idx = sorted([sample_idx1, sample_idx2])
                            sorted_idx = tuple(sorted_idx)

                            self.kinship[sorted_idx] += kinship

        return self.kinship


    @staticmethod
    def test_kinship(G=None, pedfile=None):
        if ((G is None and pedfile is None) or 
            (G is not None and pedfile is not None)):
            raise ValueError("Must provide either Graph or pedigree text file")
        
        if G is not None:
            PK = PedKinship(G)
        else:
            PK = PedKinship.from_txt(pedfile)
            
        probands = PK.get_probands()
        kinship = PK.compute_kinship(probands[:10])
        
        return kinship
    

class TSKinship:
    # TODO: Pull out individuals based on ID stored in
    #       metadata
    
    @staticmethod
    def get_sample_inds(ts):
        # NOTE: This only checks if the first associated
        #       node is a sample, which I think is enough
        sample_inds = []
        sample_nodes = set(ts.samples())
        for ind in ts.individuals():
            if len(ind.nodes) > 0:
                if ind.nodes[0] in sample_nodes:
                    sample_inds.append(ind)
                    
        return sample_inds
    
    @staticmethod
    def build_sample_ind_idx_dict(ts, individuals):
        ind_node_dict = {}
        for i, ind in enumerate(individuals):
            if ind.id in set(ts.samples()):
                ind_node_dict.update({n: i for n in ind.nodes})
                
        return ind_node_dict
        
    @staticmethod
    def compute_node_kinship(ts, sample_nodes, max_depth=5):
        n = len(sample_nodes)

        kinship = np.zeros((n, n))
        sample_idx_dict = dict(zip(sample_nodes, range(n)))
        
        for tree in ts.trees():
            L = tree.get_length()
            for s1, s2 in itertools.combinations(sample_nodes, 2):
                if tree.tmrca(s1, s2) > max_depth:
                    continue
                    
                # Make sure we're only filling upper-triangular
                # portion of matrix
                s1_idx = sample_idx_dict[s1]
                s2_idx = sample_idx_dict[s2]
                sorted_idx = sorted([s1_idx, s2_idx])
                sorted_idx = tuple(sorted_idx)
                
                kinship[sorted_idx] += L / ts.get_sequence_length()
                
        return kinship
    
    @staticmethod
    def compute_ind_kinship(ts, individuals, max_depth=5):
        n = len(individuals)
        kinship = np.zeros((n, n))
        min_tmrca = 1000
        closest_pair = None
        
        sample_ind_idx_dict = TSKinship.build_sample_ind_idx_dict(ts, individuals)
        sample_nodes = sample_ind_idx_dict.keys()
        for k, v in sample_ind_idx_dict.items():
            print(k, v)
        
        for tree in tqdm(ts.trees(), total=ts.num_trees):
            L = tree.get_length()
            for s1, s2 in itertools.combinations(sample_nodes, 2):

                if tree.tmrca(s1, s2) > max_depth:
                    continue
                    
#                 if np.abs(s2 - s1) > 1:
#                     print(s1, s2)
#                     print(tree.tmrca(s1, s2))
#                     assert(1 == 2)
                # Make sure we're only filling upper-triangular
                # portion of matrix
                s1_idx = sample_ind_idx_dict[s1]
                s2_idx = sample_ind_idx_dict[s2]
                sorted_idx = sorted([s1_idx, s2_idx])
                sorted_idx = tuple(sorted_idx)
                
                if min_tmrca > tree.tmrca(s1, s2):
                    min_tmrca = tree.tmrca(s1, s2)
                    closest_pair = (s1_idx, s2_idx)
                    
                kinship[sorted_idx] += L / (4 * ts.get_sequence_length())
#                 kinship[sorted_idx] += 1 / 2 ** (2 * tree.tmrca(s1, s2))
    
        print(min_tmrca, closest_pair)
                
        return kinship


def sibs_plus_cousin_data():
    script_root = os.path.join(os.path.dirname(__file__), '..')
    pedfile = os.path.join(
            script_root, 'test', 'test_data', 'sibs_plus_cousin.txt')
    pedigree = msprime.Pedigree.read_txt(pedfile, time_col=3)

    length = 10
    n = 3
    ts = msprime.simulate(n, pedigree=pedigree, length=length, model="wf_ped")

    return (ts, pedfile)


class TestKinship:
    def test_sibs_plus_counsin(self):
        ts, pedfile = sibs_plus_cousin_data()

        PK = PedKinship.from_txt(pedfile)
        probands = PK.get_probands()
        print(PK.compute_kinship(probands))

        ts_sample_inds = TSKinship.get_sample_inds(ts)
        print(ts_sample_inds)
        print(TSKinship.compute_ind_kinship(ts, ts_sample_inds))
