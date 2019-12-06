import sys
import os
import argparse
import subprocess

sys.path.insert(0, 'msprime')
import msprime


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


def ts_to_bcf_single(ts_file, out_file, runner):
    make_bcf_cmd = "tskit vcf --ploidy 2 {} | bcftools view -O b > {}".format(
            ts_file, out_file)
    runner.run(make_bcf_cmd)

    ts = msprime.load(ts_file)
    assert (ts.num_samples % 2 == 0)
    num_samples = int(ts.num_samples / 2)
    new_sample_ids_file = ".tmp_sample_ids.txt"
    with open(new_sample_ids_file, 'w') as f:
        for i in range(1, num_samples+1):
            f.write('tsk_' + str(i) + '\n')
    reheader_cmd = ("bcftools reheader --samples {} {} > tmp && "
            "mv tmp {} && rm {}").format(
            new_sample_ids_file, out_file, out_file, new_sample_ids_file)
    runner.run(reheader_cmd)


def bcf_convert_chrom(bcf_file, chrom_num, runner):
    tmp_conv_file = '.tmp_chr_conv.txt'

    with open(tmp_conv_file, 'w') as f:
        f.write('1\t' + str(chrom_num) + '\n')

    convert_chrom_cmd = (
            'bcftools annotate --rename-chrs {} {} | '
            'bcftools view -O b > tmp && mv tmp {} && '
            'rm {}').format(
            tmp_conv_file, bcf_file, bcf_file, tmp_conv_file)
    runner.run(convert_chrom_cmd)


def concat_bcf(bcf_files, out_file, runner, delete_originals=True):
    concat_chrom_cmd = 'bcftools concat -o {} -O b --threads 2'

    for f in bcf_files:
        concat_chrom_cmd += ' ' + f

    concat_chrom_cmd = concat_chrom_cmd.format(out_file)

    runner.run(concat_chrom_cmd)

    if delete_originals:
        rm_cmd = 'rm {}'
        for f in bcf_files:
            f = os.path.expanduser(f)
            if not os.path.isabs(f):
                f = os.path.abspath(f)

            runner.run(rm_cmd.format(f))


def main(args):
    out_file = os.path.expanduser(args.out_file)
    runner = Runner(args)

    tmp_bcf_files = []
    for i, tsf in enumerate(args.ts_file):
        chrom_num = i + 1
        tmp_bcf_file = '.tmp' + str(i) + '.bcf'
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
