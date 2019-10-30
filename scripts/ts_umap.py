import sys
import os
sys.path.append('~/project/msprime')
import msprime
import numpy as np
import umap
import argparse

make_bcf_cmd = "tskit vcf --ploidy 2 {} | bcftools view -O b > {}".format(
        ts_file, bcf_file)

new_sample_ids_file = "tempfile"

num_samples = ts.num_samples
with open(new_sample_ids_file, 'w') as f:
    for i in range(1, num_samples+1):
        f.write('tsk_' + str(i) + '\n')

reheader_cmd = "bcftools reheader --samples {} {} > {}".format(
        new_sample_ids_file, bcf_file, new_bcf_file)

plink_pca_cmd = "plink --pca --bcf {} --memory 4000".format(new_bcf_file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argumant('-t', '--ts-file', required=True)
    parser.add_argumant('-o', '--output-dir', required=True)
    parser.add_argumant('-t', '--ts-file', required=True)
    parser.add_argumant('-t', '--ts-file', required=True)
