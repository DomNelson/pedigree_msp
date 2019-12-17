import sys, os
sys.path.insert(0, "msprime")  # NOQA
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
    except:  # NOQA
        print("Unexpected error!")
        embed()
        raise


def main(args):
    if args.pedarray is not None and args.time_col:
        raise ValueError(
                "Times are always loaded from pedigrees specified in " +\
                "numpy .npy format - cannot also specify `time_col`")

    if args.num_samples and args.samples_file:
        raise ValueError("Cannot specify both num_samples and samples_file")

    if args.convert:
        if args.outfile is None:
            raise ValueError(
                    "Must provide output file for .txt -> .npy conversion")
        if args.pedarray is not None:
            raise ValueError(
                    "Can only specify pedigree in text format for conversion.")

    num_samples = args.num_samples
    sample_IDs = None
    if args.samples_file:
        with open(os.path.expanduser(args.samples_file), 'r') as f:
            sample_IDs = [int(ID) for ID in f]
        num_samples = len(sample_IDs)

    pedigree = None
    if args.pedfile:
        time_col = None
        if args.load_times:
            time_col = 3
        pedigree = msprime.Pedigree.read_txt(args.pedfile, time_col=time_col)
    elif args.pedarray:
        pedigree = msprime.Pedigree.read_npy(os.path.expanduser(args.pedarray))

    if args.convert:
        print("Saving", args.pedfile, "as .npy file in:", args.outfile)
        pedigree_outfile = os.path.expanduser(args.outfile)
        pedigree.save_npy(pedigree_outfile)
        sys.exit()

    pedigree.set_samples(num_samples, sample_IDs)

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

    replicates = msprime.simulate(
            num_samples, Ne=args.popsize, pedigree=pedigree,
            model='wf_ped', mutation_rate=args.mu,
            # length=args.length, recombination_rate=args.rho,
            recombination_map=rm,
            end_time=args.end_time,
            demographic_events=des, num_replicates=args.replicates)

    ## Check that all IDs in the tree sequence are contained in the pedigree
    if args.replicates is None:
        ts = replicates
        if args.check_inds:
            check_inds(ts, pedigree)
        if args.outfile is not None:
            outfile = os.path.expanduser(args.outfile)
            ts.dump(outfile)
    else:
        if args.check_inds:
            for ts in replicates:
                check_inds(ts, pedigree)

    if args.outfile is None:
        embed()


if __name__ == "__main__":
    example_text = \
"""Pedigree text column format:
    Ind ID, Father ID, Mother ID, Time (generations in past)

Simulation example:
    python pedigree_example.py --pedfile path/to/ped_file.txt

Text-to-numpy conversion example:
    python pedigree_example.py --pedfile path/to/ped_file.txt \\
            --outfile /path/to/new/ped_file.npy --convert"""

    parser = argparse.ArgumentParser(epilog=example_text,
            formatter_class=argparse.RawDescriptionHelpFormatter)

    ped_group = parser.add_argument_group("Require one of")
    pedigree_input = ped_group.add_mutually_exclusive_group(required=True)
    pedigree_input.add_argument('-p', '--pedfile', metavar='',
            help="""Path to pedigree text file""")
    pedigree_input.add_argument('-a', '--pedarray', metavar='',
            help="""Path to pedigree in numpy .npy format. Loading time is much
            faster than text format.""")

    utils_group = parser.add_argument_group("Non-simulation options")
    utils_group.add_argument("-v", "--convert", action='store_true',
            help="""Flag to convert pedigree text file to numpy .npy format,
            and save to `outfile`. Simulations will not be performed and all
            other arguments are ignored""")

    add_argument_no_metavar = \
            lambda *a, **kwa: parser.add_argument(*a, **kwa, metavar='')
    add_argument_no_metavar('-o', '--outfile',
            help="""Path to store simulation output in tree sequence hdf5
            format. If not given, drops into IPython interpreter where the
            output tree sequence `ts` can be examined, and saved using
            `ts.dump("path/to/output/file")`""")
    add_argument_no_metavar('-n', '--num_samples', type=int, default=10,
            help="""default: %(default)d""")
    add_argument_no_metavar('-s', '--samples_file',
            help="""Path to file containing pedigree IDs of individuals to
            simulate, one ID per line. Overrides `num_samples`.""")
    add_argument_no_metavar('-e', '--end_time', type=int)
    add_argument_no_metavar('-d', '--model', choices=['dtwf', 'hudson'],
            help="""Simulation model to switch to after reaching the top of the
            pedigree. Leave unspecified to stop when pedigree founders are
            reached.""")
    add_argument_no_metavar('-N', '--popsize', type=int, default=100,
            help="""Effective population size to use when simulating beyond the
            top of the pedigree. Default: %(default)d""")
    add_argument_no_metavar('-m', '--mu', type=float, default=1e-8,
            help="""Mutation rate. Default: 1e-8""")
    add_argument_no_metavar('-l', '--length', type=float, default=1e6,
            help="""Length of simulated region. Default: 1e6""")
    add_argument_no_metavar('-r', '--rho', type=float, default=1e-8,
            help="""Recombination rate. Default: 1e-8""")
    add_argument_no_metavar('-R', '--replicates', type=int,
            help="""Number of replicate simulations to perform. Default: 1""")
    parser.add_argument('-c', '--check_inds', action="store_true",
            help="""Debug option - verify all individuals IDs in the simulated
            tree sequence are contained within the pedigree""")
    parser.add_argument('-t', '--load_times', action="store_true",
            help="""Flag to load times as 4th column of the pedigree. Without
            this flag, times are (roughly) inferred as a partial ordering of
            ancestors. Note that pedigrees in numpy .npy format always have
            times baked in.""")

    args = parser.parse_args()

    main(args)
