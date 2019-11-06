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

    if args.pedarray and args.samples:
        raise ValueError("Cannot specify samples - already set in pedarray")

    samples = 10
    if args.samples:
        samples = args.samples

    pedigree = None
    if args.pedfile:
        pedigree = os.path.expanduser(args.pedfile)
    elif args.pedarray:
        pedigree = np.load(os.path.expanduser(args.pedarray))
        sample_col = 4
        samples = np.sum(pedigree[:, sample_col])

    ts = msprime.simulate(samples, Ne=args.popsize, pedigree=pedigree,
            model='dtwf', mutation_rate=args.mu, length=args.length,
            recombination_rate=args.rho, end_time=args.end_time)

    outfile = os.path.expanduser(args.outfile)
    ts.dump(outfile)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--pedfile')
    parser.add_argument('-a', '--pedarray')
    parser.add_argument('-o', '--outfile', required=True)
    parser.add_argument('-n', '--samples', type=int)
    parser.add_argument('-e', '--end_time', type=int)
    parser.add_argument('-N', '--popsize', type=int, default=100)
    parser.add_argument('-m', '--mu', type=float, default=1e-8)
    parser.add_argument('-l', '--length', type=float, default=1e6)
    parser.add_argument('-r', '--rho', type=float, default=1e-8)

    args = parser.parse_args()

    main(args)
