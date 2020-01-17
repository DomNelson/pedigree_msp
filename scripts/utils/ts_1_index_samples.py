import os
import sys
sys.path.insert(0, 'msprime')
import msprime
import subprocess

def get_bcf_num_samples(vcf_file):
    cmd = "awk '{{print NF}}' <(bcftools view {} | tail) | sort -nu | tail -n 1".format(vcf_file)
    output = subprocess.run(cmd, shell=True, capture_output=True, check=True, executable='/bin/bash')
    num_samples = int(output.stdout) - 9

    return num_samples
    

if len(sys.argv) != 3:
    print("Usage: python ts_1_index_samples.py vcf_file out_file")
    sys.exit()

vcf_file = os.path.expanduser(sys.argv[1])
out_file = os.path.expanduser(sys.argv[2])

num_samples = get_bcf_num_samples(vcf_file)
assert (num_samples % 2 == 0)
new_sample_ids_file = ".tmp_sample_ids.txt"
with open(new_sample_ids_file, 'w') as f:
    for i in range(1, num_samples+1):
        f.write('tsk_' + str(i) + '\n')
reheader_cmd = ("bcftools reheader --samples {} {} > tmp && "
        "mv tmp {} && rm {}").format(
        new_sample_ids_file, vcf_file, out_file, new_sample_ids_file)
subprocess.run(reheader_cmd, shell=True, check=True)

# Make sure we have the same number of samples as before
assert(get_bcf_num_samples(out_file) == num_samples)
