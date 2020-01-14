import sys, os
import numpy as np
sys.path.insert(0, "../msprime")
import msprime
import argparse
import subprocess

def example_sim():
    pedfile = '/home/dnelson/project/pedigree_msp/data/pedEx.txt'
    pedigree = msprime.Pedigree.read_txt(pedfile)
    max_ped_time = max(pedigree.time)

    population_configurations = [msprime.PopulationConfiguration(100)]
    demographic_events = [msprime.SimulationModelChange(max_ped_time, 'dtwf')]

    ts = msprime.simulate(100, Ne=100, pedigree=pedigree, model='wf_ped',
                          demographic_events=demographic_events, length=1e6,
                          recombination_rate=1e8)
    
    return ts


def ts_ibd_to_plink(ibd_array, chr_num=1, LOD=100):
    """
    plink format uses the following columns:
    
        ID_1 copy_1 ID_2 copy_2 chr_num start_pos end_pos LOD
    """
    
    for row in ibd_array:
        ID_1, ID_2, start, end = row
        
        if not (isinstance(start, np.integer) and isinstance(end, np.integer)):
        # if type(start) != int or type(end) != int:
            raise ValueError("Start and end positions must be integers.")
            
        copy_1 = 1
        copy_2 = 1
        plink_vals = [ID_1, copy_1, ID_2, copy_2, chr_num, start, end, LOD]
        plink_row = '\t'.join([str(x) for x in plink_vals])
    
        yield plink_row + '\n'
        
        
def ibd_npz_to_plink(npz_file, plink_out_file, max_rows=None):
    ibd_plink = ts_ibd_to_plink(npz_file)
    
    with open(plink_out_file, 'w') as f:
        for i, row in enumerate(ibd_plink):
            if max_rows and i > max_rows:
                break
                
            f.write(row)
    

def ibd_write_map_file(ibd_array, outfile, chrom_num=1):
    """
    Recombination map format:
    
        Chrom ID, variant ID, cM position, base-pair position
        
    where cM position can have a dummy value of 0
    """
    boundary_positions = set(ibd_array[:, 2]).union(ibd_array[:, 3])
    boundary_positions = sorted(boundary_positions)

    with open(outfile, 'w') as f:
        for pos in boundary_positions:
            cM_pos = pos / 1e6
            row = '{} . {} {}\n'.format(chrom_num, cM_pos, pos)
            f.write(row)
            
            
def start_end_to_int(ibd_array):
    start_col = 2
    end_col = 3
    
    ibd_array[:, start_col] = ibd_array[:, start_col].astype(int)
    ibd_array[:, end_col] = ibd_array[:, end_col].astype(int)
    
    return ibd_array


def write_plink_files(ibd_file, plink_ibd_out, plink_map_out):
    ## IBD Ne expects integer-valued positions, so we convert everything to
    ## ints
    ibd_array = np.load(ibd_file)['ibd_array'].astype(int)

    ibd_npz_to_plink(ibd_array, plink_ibd_out)
    ibd_write_map_file(ibd_array, plink_map_out)


def run_ibd_ne(jar_file, plink_ibd_out, **kwargs):
    cmd = 'cat {} | java -jar {}'.format(
            plink_ibd_out, jar_file)

    for key, value in kwargs.items():
        cmd += ' {}={}'.format(key, value)

    subprocess.run(cmd, shell=True, check=True)


def main(args):
    plink_ibd_out = args.prefix + '.ibd'
    plink_map_out = args.prefix + '.map'

    write_plink_files(args.ibd_file, plink_ibd_out, plink_map_out)
    run_ibd_ne(args.jar_file, plink_ibd_out,
            map=plink_map_out,
            out=args.prefix,
            nboots=0)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--ibd_file', metavar='', required=True,
                        help="""Path to IBD file in .npz format""")
    parser.add_argument('-o', '--prefix', metavar='', required=True,
                        help="""Prefix for all output files""")
    parser.add_argument('-j', '--jar_file', metavar='', default='ibdne.19Sep19.268.jar',
                        help="""Path to jar file for IBD Ne""")

    args = parser.parse_args()
    main(args)
