#!/usr/bin/python$

import argparse
import subprocess
import os
import numpy
import sys
import multiprocessing as mp
#from gen_tree_overlappingCNA import gen_tree
from gen_readcount import gen_readcount
from gen_readcount import get_beta_dist
#from Gen_Ref_Fa import make_fa, make_fa_wABs 
from Gen_Ref_Fa import get_ref_chr_len_name
from joblib import Parallel, delayed
from tqdm import tqdm
import logging
## set logger
logging.basicConfig(
format='%(asctime)s %(levelname)-8s %(message)s',
level=logging.INFO,
datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(" ")
#from Gen_Ref_Fa import getlen_ref

# given normal genome size (diploid)
#genome_size = 6.32 * 10^9

def run_one_seg(dir,fa_f, chr_name, start, end,this_readcount,l,out_fq1,out_fq2):
    # chr_name_array[j]

    tmp_fa_file = "_".join([dir+'/'+fa_f.split('/')[-1], str(chr_name), str(start), str(end)]) + ".fa" 
    # use samtools faidx to get this sequence that is to be sequenced into reads
    args = "samtools faidx " + fa_f + " " + chr_name + ":" + str(start) + "-" + str(end) + " > " + tmp_fa_file
# TODO check how many N's 
# conclusion: Popen cannot exit. check_call can. But check_call does not work on all Ns as wgsim does not work on it. 
    print(args)
    try:
        subprocess.check_call(args, shell=True)
    except subprocess.CalledProcessError:
        print("Cannot work on " + args) 

    # check N's
    N_true = check_Ns(tmp_fa_file)

    if N_true == 1:
        args = wgsim_dir + "wgsim -h -N " + str(this_readcount) + " -1 " + str(l) + " -2 " + str(l) + " " + tmp_fa_file + " " + out_fq1 + " " + out_fq2
        print(args)
        try:
            subprocess.check_call(args, shell=True)
        except subprocess.CalledProcessError:
            print("Cannot work on " + args) 

    args = "rm " + tmp_fa_file 
    try:
        subprocess.check_call(args, shell=True)
    except subprocess.CalledProcessError:
        print("Cannot remove " + tmp_fa_file)
    return 

# for both bulk and single-cell sampling purposes
def gen_reads(dir, chr_len, chr_name, fq_prefix, cell_i, template_fa, Alpha, Beta, x0, y0, cov, l, window_size, u):
    logger.info("Now sequencing cell " + str(cell_i))

    out_fq1 = dir + "/" + fq_prefix + "_1.fq"
    out_fq2 = dir + "/" + fq_prefix + "_2.fq"
    os.system('rm '+ out_fq1+' '+out_fq2)

    # each chromosome
    for j in range(len(chr_len)):
        
        this_chrlen = chr_len[j]
        readcounts = gen_readcount(dir, Alpha, Beta, x0, y0, cov, l, window_size, u, this_chrlen)
        # each bin
        start = 0
        end = 0
        for k in readcounts:
            this_readcount = k
            # extract this segment of fa
            start = end + 1
            end = start + window_size - 1
            if end > this_chrlen:
                end = this_chrlen
            if start > this_chrlen:
                break

            run_one_seg(dir, template_fa, chr_name[j], start, end, k, l, out_fq1, out_fq2)

def check_Ns(file):
    N_line = 0
    total_line = 0
    with open(file, "r") as f:
        for line in f:
            total_line = total_line + 1
            if line.find('N') != -1:
                N_line = N_line + 1
    f.close()
    if total_line != 0 and N_line/float(total_line) > 0.2:
        return 0
    return 1

if len(sys.argv) <= 1:
    print("""
    Simulate a set number of normal cells given a genome in fasta format, the coverage of which fluctuates and mimics the real single cell data. 
    Usage: python sim_normcells.py -t [ref.fa] -n [number_cells] -S [wgsim_dir]
        -p (--processors)   Numbers of processors available.
        -r (--directory)    Location of simulated data. The program will remove the whole directory if it already exists. Otherwise it will create one. (default: test)
        -S (--wgsim-dir)    The directory of the binary of wgsim. It is in the same folder of this main.py. (need to specify) 
        -n (--normal-num)     Number of the cells on a level of interest. Always greater than -F treewidth. Treewidth controls the total number of clones whereas cell-num controls the total number of cells sequenced at a certain tree depth. (default: 8)
        -t (--template-ref) The reference file to sequence the reads. 
        -o (--outfile)      The standard output file, will be saved in output folder, just give the file name. (default: std.out)
        -f (--fq-prefix)    The prefix of the fastq names. (default: normal) This is the only line that is different from main.par.overlapping.py.
        -x (--Lorenz-x)     The value on the x-axis of the point furthest from the diagonal on the Lorenz curve imitating the real coverage uneveness. (default: 0.5) 
        -y (--Lorenz-y)     The value on the y-axis of the Lorenz curve imitating the real coverage unevenness. x > y. The closer (x, y) to the diagonal, the better the coverage evenness. (default: 0.4) 
        -v (--coverage)     The average coverage of the sequence. (default: 0.02)
        -l (--readlen)      Read length for each read sequenced. (default: 35bp)
        -w (--window-size)  Within a window, the coverage is according to a Gaussian distribution. Neighboring windows' read coverage is according to a Metropolis Hasting process. (default: 200000bp)
        -u (--acceptance-rate)  The probability to accept a proposal in Metropolis Hasting. (default: 0.5)
        """)
    sys.exit(0)

parser = argparse.ArgumentParser(description='Simulate the normal cells whose coverage fluctuates to mimic single cell sequencing . ')
parser.add_argument('-p', '--processors', default=20)
parser.add_argument('-r', '--directory', default="test")
parser.add_argument('-S', '--wgsim-dir', default="")
parser.add_argument('-n', '--normal-num', default=8)
parser.add_argument('-t', '--template-ref', default="hg19.fa")
parser.add_argument('-o', '--outfile', default="std.out")
parser.add_argument('-f', '--fq-prefix', default="normal")
parser.add_argument('-x', '--Lorenz-x', default=0.5)
parser.add_argument('-y', '--Lorenz-y', default=0.4)
parser.add_argument('-v', '--coverage', default=0.02)
parser.add_argument('-l', '--readlen', default=35)
parser.add_argument('-w', '--window-size', default=200000)
parser.add_argument('-u', '--acceptance-rate', default=0.5)


args = parser.parse_args()
NUM_OF_PROCESSES = int(args.processors)
dir = args.directory
wgsim_dir = args.wgsim_dir
n = int(args.leaf_num)
template_ref = args.template_ref
outfile = dir + "/" + args.outfile
fq_prefix = args.fq_prefix
x0 = float(args.Lorenz_x)
y0 = float(args.Lorenz_y)
cov = float(args.coverage)
l = int(args.readlen)
window_size = int(args.window_size)
u = float(args.acceptance_rate)

[chr_name, chr_len] = get_ref_chr_len_name(template)

if not os.path.exists(dir):
    subprocess.check_call("mkdir " + dir, shell=True)

[Alpha, Beta] = get_beta_dist(x0, y0)
print("Number of processes: " + str(NUM_OF_PROCESSES))
print("Number of normal cells: " + str(n))

# start the sequencing
for cell_i in range(n):
    # sequence each cell  
    os.system("mkdir -p "+dir+'/cell%d/'%cell_i)

    # n_thread = 15
sequences = Parallel(n_jobs=NUM_OF_PROCESSES )(delayed(gen_reads)(dir+'/cell%d/'%cell_i, 
                chr_len, 
                chr_name,
                fq_prefix,
                cell_i,
                template_fa, 
                Alpha, 
                Beta, 
                x0, 
                y0, 
                cov, 
                l, 
                window_size, 
                u) for cell_i in tqdm(range(n)))


   
