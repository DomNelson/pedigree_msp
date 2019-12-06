import os
import numpy as np
import umap
import subprocess
import argparse

import sys
sys.path.append('msprime')
import msprime

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from IPython import embed


def prototype():
    pass


def get_prefix(args):
    if args.prefix:
        prefix = os.path.expanduser(args.prefix)
        if not os.path.isabs(prefix):
            prefix = os.path.abspath(args.prefix)
    else:
        f = os.path.expanduser(args.ts_file)
        base, ext = os.path.splitext(f)
        _, basename = os.path.split(base)
        out_dir = os.path.expanduser(args.output_dir)
        if not os.path.isabs(out_dir):
            out_dir = os.path.abspath(out_dir)
        prefix = os.path.join(os.path.abspath(out_dir), basename)

    return prefix


def main(args):
    out_dir = os.path.expanduser(args.output_dir)
    ts_file = os.path.expanduser(args.ts_file)
    ts_time = os.path.getmtime(args.ts_file)

    prefix = get_prefix(args)
    bcf_file = prefix + '.bcf'

    if (not os.path.exists(bcf_file) or
            os.path.getmtime(bcf_file) < ts_time):
        print("Writing BCF from tree sequence file...")
        make_bcf_cmd = "tskit vcf --ploidy {} {} | bcftools view -O b > {}".format(
                args.ploidy, ts_file, bcf_file)
        subprocess.run(make_bcf_cmd, shell=True, check=True)

        ts = msprime.load(ts_file)
        assert (ts.num_samples % args.ploidy == 0)
        num_samples = int(ts.num_samples / args.ploidy)
        new_sample_ids_file = "tempfile"
        with open(new_sample_ids_file, 'w') as f:
            for i in range(1, num_samples+1):
                f.write('tsk_' + str(i) + '\n')

        print("Adjusting sample IDs for PCA...")
        reheader_cmd = ("bcftools reheader --samples {} {} > tmp && "
                "mv tmp {} && rm {}").format(
                new_sample_ids_file, bcf_file, bcf_file, new_sample_ids_file)
        subprocess.run(reheader_cmd, shell=True, check=True)
    else:
        print("Using existing BCF file:", bcf_file)
    bcf_time = os.path.getmtime(bcf_file)

    evecs_file = prefix + '.eigenvec'
    if (not os.path.exists(evecs_file) or
            os.path.getmtime(evecs_file) < bcf_time):

        print("Running PCA")
        plink_pca_cmd = "plink --pca --bcf {} --memory 4000 --out {}".format(
                bcf_file, prefix)
        subprocess.run(plink_pca_cmd, shell=True, check=True)
    else:
        print("Using existing PCA file:", evecs_file)
    pca_time = os.path.getmtime(bcf_file)

    print("Running UMAP")
    umap_file = prefix + '_umap.npy'
    data = np.genfromtxt(evecs_file)
    IDs, evecs = data[:, :2], data[:, 2:]

    embedding = umap.UMAP(n_neighbors=5,
                          min_dist=0.1).fit_transform(evecs)
    np.save(umap_file, embedding)

    print("Plotting")
    plot_file = prefix + '_umap.png'
    fig, ax = plt.subplots()
    ax.scatter(embedding[:, 0], embedding[:, 1], s=1)
    fig.savefig(plot_file)

    if args.ipython:
        embed()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--ts-file')
    parser.add_argument('-p', '--prefix')
    parser.add_argument('-o', '--output-dir', required=True)
    parser.add_argument('-d', '--ploidy', type=int, default=2)
    parser.add_argument('-I', '--ipython', action='store_true')

    args = parser.parse_args()

    if not (args.ts_file or args.prefix):
        raise ValueError("Must provide either a tree sequence or a prefix")

    main(args)
