import sys, os
sys.path.append(os.path.expanduser('~/project/msprime'))
import msprime
import argparse
from IPython import embed


def main(args):
    pedfile = os.path.expanduser('~/project/pedigree_msp/data/pedEx.txt')
    if args.pedfile:
        pedfile = os.path.expanduser(args.pedfile)

    ts = msprime.simulate(args.samples, Ne=args.popsize, pedigree=pedfile,
            model='dtwf', mutation_rate=args.mu)
    embed()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--pedfile')
    parser.add_argument('-n', '--samples', type=int, default=10)
    parser.add_argument('-N', '--popsize', type=int, default=100)
    parser.add_argument('-m', '--mu', type=float, default=1e-7)

    args = parser.parse_args()

    main(args)
