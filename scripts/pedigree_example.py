import sys, os
project_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(os.path.join(project_dir, 'msprime'))
import msprime
import argparse
from IPython import embed
import numpy as np


def main(args):
    if args.pedfile and args.pedarray:
        raise ValueError("Cannot specify both pedfile and pedarray")

    samples = 10
    if args.samples:
        samples = args.samples

    pedigree = None
    if args.pedfile:
        pedigree = msprime.Pedigree.read_txt(args.pedfile)
    elif args.pedarray:
        pedigree = msprime.Pedigree.read_npy(os.path.expanduser(args.pedarray))
        pedigree.set_samples(args.samples)

    ts = msprime.simulate(samples, Ne=args.popsize, pedigree=pedigree,
            model='wf_ped', mutation_rate=args.mu, length=args.length,
            recombination_rate=args.rho, end_time=args.end_time)

    ## Check that all IDs in the tree sequence are contained in the pedigree
    ids = [int(ind.metadata) for ind in ts.individuals()]
    id_diff = list(set(ids).difference(pedigree.inds))
    if len(id_diff) > 0:
        print("Invalid inds in tree sequence:", id_diff)

    if args.outfile:
        outfile = os.path.expanduser(args.outfile)
        ts.dump(outfile)
    else:
        embed()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--pedfile')
    parser.add_argument('-a', '--pedarray')
    parser.add_argument('-o', '--outfile')
    parser.add_argument('-n', '--samples', type=int)
    parser.add_argument('-e', '--end_time', type=int)
    parser.add_argument('-N', '--popsize', type=int, default=100)
    parser.add_argument('-m', '--mu', type=float, default=1e-8)
    parser.add_argument('-l', '--length', type=float, default=1e6)
    parser.add_argument('-r', '--rho', type=float, default=1e-8)

    args = parser.parse_args()

    main(args)
