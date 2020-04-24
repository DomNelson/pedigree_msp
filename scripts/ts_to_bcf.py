import sys
import os
import argparse
import subprocess
import numpy as np
import tempfile

sys.path.insert(0, 'msprime')
import msprime # NOQA


class Runner:
    def __init__(self, args):
        self.test = args.test
        self.verbose = (args.verbose or args.test)

    def run(self, cmd):
        prefix = "Running:"
        if self.test:
            prefix = "Testing:"

        if self.verbose:
            print(prefix, cmd)

        if not self.test:
            subprocess.run(cmd, shell=True, check=True)


def ts_clean_inds(ts):
    ll_tables = ts.dump_tables().asdict()
    sample_nodes = set(ts.samples())

    # Collect sampled individuals, checking they have two sample nodes
    sample_individuals = []
    for ind in ts.individuals():
        if len(ind.nodes) == 0:
            continue

        if ind.nodes[0] in sample_nodes:
            assert len(ind.nodes) == 2
            assert ind.nodes[1] in sample_nodes
            sample_individuals.append(ind)

    # Since we depend on the sampled individuals being the first $n$
    # individuals, check that this is in fact the case, and pull out the
    # last individual
    sample_individual_ids = [ind.id for ind in sample_individuals]
    assert np.min(sample_individual_ids) == 0
    assert np.max(np.diff(sorted(sample_individual_ids))) == 1
    last_sample_idx = np.max(sample_individual_ids)

    # Trim non-sample individuals from the tables
    ind_tables = ll_tables['individuals']
    ind_tables['flags'] = ind_tables['flags'][:last_sample_idx + 1]

    # We add one to offset indices. Since their first value is always 0,
    # our data is shifted right by 1
    ind_tables['location_offset'] = \
        ind_tables['location_offset'][:last_sample_idx + 2]
    idx = ind_tables['location_offset'][-1]
    ind_tables['location'] = ind_tables['location'][:idx]

    ind_tables['metadata_offset'] = \
        ind_tables['metadata_offset'][:last_sample_idx + 2]
    idx = ind_tables['metadata_offset'][-1]
    ind_tables['metadata'] = ind_tables['metadata'][:idx]

    # This is redundant (assigned them equal above) but is more explicit if
    # if wasn't clear that we in fact updated the original dict in-place
    ll_tables['individuals'] = ind_tables

    # Clear out references to individuals which we've removed
    for i in range(ll_tables['nodes']['individual'].shape[0]):
        ind_idx = ll_tables['nodes']['individual'][i]
        if ind_idx not in set(sample_individual_ids):
            ll_tables['nodes']['individual'][i] = -1

    # Create new tree sequence with only sample individuals and save to
    # tempfile to be converted to VCF
    ts_only_sample_inds = msprime.tskit.tables.TableCollection.fromdict(
                                ll_tables
                            ).tree_sequence()

    return ts_only_sample_inds


def ts_to_bcf_single(ts_file, out_file, runner):
    # Need to remove non-sample individuals from ts or else tskit gets confused
    ts = msprime.load(ts_file)
    if ts.num_individuals == 0:
        bcf_cmd = "tskit vcf --ploidy 2 {} | bcftools view -O b > {}".format(
                           ts_file, out_file)
        runner.run(bcf_cmd)
    else:
        ts_only_sample_inds = ts_clean_inds(ts)
        # TODO: This will probably fail if metadata isn't present
        ids = [ind.metadata.decode('utf8')
               for ind in ts_only_sample_inds.individuals()]

        read_fd, write_fd = os.pipe()
        write_pipe = os.fdopen(write_fd, "w")
        with open(out_file, "w") as f:
            proc = subprocess.Popen(
                ["bcftools", "view", "-O", "b"], stdin=read_fd, stdout=f
            )
            ts_only_sample_inds.write_vcf(
                write_pipe, individual_names=ids
            )
            write_pipe.close()
            os.close(read_fd)
            proc.wait()
            if proc.returncode != 0:
                raise RuntimeError(
                    "bcftools failed with status:", proc.returncode
                )


def bcf_convert_chrom(bcf_file, chrom_num, runner):
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as f:
        f.write('1\t' + str(chrom_num) + '\n')

        convert_chrom_cmd = (
                'bcftools annotate --rename-chrs {} {} | '
                'bcftools view -O b > tmp && mv tmp {} && '
                'rm {}').format(
                f.name, bcf_file, bcf_file, f.name)
        runner.run(convert_chrom_cmd)


def concat_bcf(bcf_files, out_file, runner):
    concat_chrom_cmd = 'bcftools concat -o {} -O b --threads 2'

    for f in bcf_files:
        concat_chrom_cmd += ' ' + f

    concat_chrom_cmd = concat_chrom_cmd.format(out_file)

    runner.run(concat_chrom_cmd)


def main(args):
    out_file = os.path.expanduser(args.out_file)
    runner = Runner(args)

    tmp_bcf_files = []
    with tempfile.TemporaryDirectory() as tmpdirname:
        for i, tsf in enumerate(args.ts_file):
            chrom_num = i + 1
            tmp_bcf_file = os.path.join(tmpdirname, '.tmp' + str(i) + '.bcf')
            tmp_bcf_files.append(tmp_bcf_file)

            ts_to_bcf_single(tsf, tmp_bcf_file, runner)
            bcf_convert_chrom(tmp_bcf_file, chrom_num, runner)

        concat_bcf(tmp_bcf_files, out_file, runner)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--ts-file', nargs="*", required=True)
    parser.add_argument('-o', '--out-file', required=True)
    parser.add_argument('-v', '--verbose', action='store_true')
    parser.add_argument('-T', '--test', action='store_true')

    args = parser.parse_args()

    main(args)
