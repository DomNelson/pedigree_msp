import os
import numpy as np
import umap
import subprocess
import argparse

import sys
sys.path.append(os.path.expanduser('~/project/msprime'))
import msprime

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from IPython import embed


def prototype():
    pass

def main(args):
    ts_file = os.path.expanduser(args.ts_file)
    out_dir = os.path.expanduser(args.output_dir)

    base, ext = os.path.splitext(ts_file)
    _, basename = os.path.split(base)
    prefix = os.path.join(os.path.abspath(out_dir), basename)

    bcf_file = prefix + '.bcf'

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

    print("Running PCA")
    plink_pca_cmd = "plink --pca --bcf {} --memory 4000 --out {}".format(
            bcf_file, prefix)
    subprocess.run(plink_pca_cmd, shell=True, check=True)

    print("Running UMAP")
    evecs_file = prefix + '.eigenvec'
    data = np.genfromtxt(evecs_file)
    IDs, evecs = data[:, :2], data[:, 2:]

    embedding = umap.UMAP(n_neighbors=5,
                          min_dist=0.3,
                          metric='correlation').fit_transform(evecs)

    fig, ax = plt.subplots()
    ax.scatter(embedding[:, 0], embedding[:, 1], s=1)
    fig.savefig(prefix + '.png')

    embed()



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--ts-file', required=True)
    parser.add_argument('-o', '--output-dir', required=True)
    parser.add_argument('-p', '--ploidy', type=int, default=2)

    args = parser.parse_args()

    main(args)
