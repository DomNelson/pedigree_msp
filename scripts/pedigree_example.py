import sys, os
sys.path.append("msprime")
import msprime
import argparse
from IPython import embed
import numpy as np


def check_inds(ts, pedigree):
    try:
        ids = [int(ind.metadata) for ind in ts.individuals()]
        id_diff = list(set(ids).difference(pedigree.inds))
        if len(id_diff) > 0:
            print("Invalid inds in tree sequence:", id_diff)
    except:
        print("Unexpected error!")
        embed()
        raise


def main(args):
    if args.pedfile and args.pedarray:
        raise ValueError("Cannot specify both pedfile and pedarray")

    samples = args.samples

    pedigree = None
    if args.pedfile:
        time_col = None
        if args.load_time:
            time_col = 3
        pedigree = msprime.Pedigree.read_txt(args.pedfile, time_col=time_col)
    elif args.pedarray:
        pedigree = msprime.Pedigree.read_npy(os.path.expanduser(args.pedarray))
        pedigree.set_samples(args.samples)

    ## Build demographic events for model changes after pedigree sims
    des = []
    if args.model is not None:
        pedigree_end_time = max(pedigree.times)
        des.append(msprime.SimulationModelChange(pedigree_end_time, args.model))

    ## For testing, sometimes need to specify num_loci directly
    rm = msprime.RecombinationMap(
            [0, int(args.length)],
            [args.rho, 0],
            int(args.length))

    replicates = msprime.simulate(samples, Ne=args.popsize, pedigree=pedigree,
            model='wf_ped', mutation_rate=args.mu,
            # length=args.length, recombination_rate=args.rho,
            recombination_map=rm,
            end_time=args.end_time,
            demographic_events=des, num_replicates=args.replicates)

    ## Check that all IDs in the tree sequence are contained in the pedigree
    if args.replicates is None and args.outfile is not None:
        outfile = os.path.expanduser(args.outfile)
        ts = replicates
        if args.check_inds:
            check_inds(ts, pedigree)
        ts.dump(outfile)
    else:
        if args.check_inds:
            for ts in replicates:
                check_inds(ts, pedigree)

        embed()



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--pedfile')
    parser.add_argument('-a', '--pedarray')
    parser.add_argument('-o', '--outfile')
    parser.add_argument('-n', '--samples', type=int, default=10)
    parser.add_argument('-e', '--end_time', type=int)
    parser.add_argument('-d', '--model')
    parser.add_argument('-N', '--popsize', type=int, default=100)
    parser.add_argument('-m', '--mu', type=float, default=1e-8)
    parser.add_argument('-l', '--length', type=float, default=1e6)
    parser.add_argument('-r', '--rho', type=float, default=1e-8)
    parser.add_argument('-R', '--replicates', type=int)
    parser.add_argument('-c', '--check_inds', action="store_true")
    parser.add_argument('-t', '--load_time', action="store_true")

    args = parser.parse_args()

    main(args)
