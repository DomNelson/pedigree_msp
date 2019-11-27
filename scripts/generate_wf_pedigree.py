import numpy as np
import sys
import os
import argparse


class PedWriter:
    def __init__(self, ninds, ngens=None):
        self.ninds = np.array(ninds).astype(int)
        if hasattr(ninds, '__iter__'):
            self.ngens = len(ninds)
        else:
            assert ngens is not None
            self.ngens = ngens
            self.ninds = [ninds] * ngens
            
        max_inds_per_gen = max([n + 1 for n in self.ninds])
        self.gen_id_digit = 10 ** (np.floor(np.log10(max_inds_per_gen)) + 1)
        self.ped_list = []
        
    def build_ped_list(self, monogamous=False):
        for i in range(self.ngens - 1):
            inds = np.arange(1, self.ninds[i] + 1) + self.gen_id_digit * i
            parents = np.arange(1, self.ninds[i + 1] + 1) + self.gen_id_digit * (i + 1)
            
            if monogamous is True:
                couples = np.random.choice(parents, size=(self.ninds[i+1] // 2, 2), replace=False)
                parent_choices_idx = np.random.choice(range(couples.shape[0]), size=self.ninds[i])
                parent_choices = couples[parent_choices_idx]
            else:
                parent_choices = np.random.choice(parents, size=(self.ninds[i], 2))
            
            for j in range(len(inds)):
                self.ped_list.append([inds[j], parent_choices[j][0], parent_choices[j][1], i])
                
        # Founders denoted by 0 for parents
        for p in set(parent_choices.ravel()):
            self.ped_list.append([p, 0, 0, self.ngens - 1])
                
    def write_ped(self, outfile):
        f = sys.stdout
        if outfile:
            f = open(outfile, 'w')

        for row in self.ped_list:
            line = '\t'.join([str(int(x)) for x in row])
            f.write(line + '\n')

        if f != sys.stdout:
            f.close()

    def print_inds_per_gen(self):
        print("Generation (in the past)\tPopulation size")
        for i, n in enumerate(self.ninds):
            print("{}\t{}".format(i, n))



def main(args):
    ngens = args.ngens
    ## FIXME: Clean up messy logic
    if args.ninds:
        ninds = [int(x) for x in args.ninds.split(',')]
        if args.ngens and len(ninds) > 1 and args.ngens != len(ninds):
            raise ValueError("Invalid number of generations specified")
        if len(ninds) == 1:
            if ngens is None:
                raise ValueError(
                        "Must specify ngens if a single value is passed to ninds")
            ninds = ninds[0]
        else:
            ngens = len(ninds)
    else:
        ninds = [int(np.exp(args.growth_rate * -t) * args.initial_size) for t in range(args.ngens)]

    outfile = None
    if args.outfile:
        outfile = os.path.expanduser(args.outfile)

    PW = PedWriter(ninds, ngens)

    if args.test:
        PW.print_inds_per_gen()
        sys.exit()

    PW.build_ped_list(args.monogamous)
    PW.write_ped(args.outfile)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", "--ninds")
    parser.add_argument("-g", "--ngens", type=int)
    parser.add_argument("-i", "--initial_size", type=float)
    parser.add_argument("-r", "--growth_rate", type=float)
    parser.add_argument("-t", "--test", action="store_true")
    parser.add_argument("-m", "--monogamous", action="store_true")
    parser.add_argument("-o", "--outfile")

    args = parser.parse_args()

    if args.initial_size is not None or args.growth_rate is not None:
        if args.initial_size is None or args.growth_rate is None or args.ngens is None:
            raise ValueError(
                    "Must specify all of initial_size, growth_rate, and "
                    "ngens if not explicitly passing ninds")

    main(args)
