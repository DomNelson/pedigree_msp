import sys
import os
import numpy as np
import umap
import subprocess
import argparse

sys.path.append('msprime')
import msprime

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from IPython import embed


def get_prefix(args):
    f = os.path.expanduser(args.bcf_file)
    base, ext = os.path.splitext(f)
    _, basename = os.path.split(base)

    out_dir = os.path.expanduser(args.output_dir)
    if not os.path.isabs(out_dir):
        out_dir = os.path.abspath(out_dir)

    prefix = os.path.join(os.path.abspath(out_dir), basename)

    return prefix


def main(args):
    out_dir = os.path.expanduser(args.output_dir)
    bcf_file = os.path.expanduser(args.bcf_file)
    bcf_time = os.path.getmtime(bcf_file)
    prefix = get_prefix(args)

    evecs_file = prefix + '.eigenvec'
    if (not os.path.exists(evecs_file) or
            os.path.getmtime(evecs_file) < bcf_time):

        print("Running PCA")
        plink_pca_cmd = "plink --pca --bcf {} --memory 4000 --out {}".format(
                bcf_file, prefix)

        if args.threads:
            plink_pca_cmd += " --threads {}".format(args.threads)

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

    plot_file = prefix + '_umap.png'
    print("Plotting to", plot_file)
    fig, ax = plt.subplots()
    ax.scatter(embedding[:, 0], embedding[:, 1], s=1)
    fig.savefig(plot_file)

    if args.ipython:
        embed()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--bcf-file', required=True)
    parser.add_argument('-o', '--output-dir', required=True)
    parser.add_argument('-t', '--threads', type=int)
    parser.add_argument('-d', '--ploidy', type=int, default=2)
    parser.add_argument('-I', '--ipython', action='store_true')

    args = parser.parse_args()

    main(args)
