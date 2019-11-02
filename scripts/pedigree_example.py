import sys, os
sys.path.append(os.path.expanduser('~/project/msprime'))
import msprime
import argparse
from IPython import embed


def main(args):
    pedfile = None
    if args.pedfile:
        pedfile = os.path.expanduser(args.pedfile)

    ts = msprime.simulate(args.samples, Ne=args.popsize, pedigree=pedfile,
            model='dtwf', mutation_rate=args.mu, length=args.length,
            recombination_rate=args.rho)

    outfile = os.path.expanduser(args.outfile)
    ts.dump(outfile)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--pedfile')
    parser.add_argument('-o', '--outfile', required=True)
    parser.add_argument('-n', '--samples', type=int, default=10)
    parser.add_argument('-N', '--popsize', type=int, default=100)
    parser.add_argument('-m', '--mu', type=float, default=1e-7)
    parser.add_argument('-l', '--length', type=float, default=1e6)
    parser.add_argument('-r', '--rho', type=float, default=1e-8)

    args = parser.parse_args()

    main(args)
