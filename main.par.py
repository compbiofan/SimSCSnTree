#!/usr/bin/python$

import argparse
import subprocess
import os
import numpy
import sys
import multiprocessing as mp
from gen_tree import gen_tree
from gen_readcount import gen_readcount
from gen_readcount import get_beta_dist
from Gen_Ref_Fa import make_fa 
from Gen_Ref_Fa import init_ref

# see how many processors I have
#NUM_OF_PROCESSES = mp.cpu_count()

def gen_reads(dir, index, leaf_index, all_chrlen, fa_prefix, Alpha, Beta, x0, y0, cov, l, window_size, u, chr_name_array):
    this_leaf_index = leaf_index[index]
    print(all_chrlen)
    print("####" + str(index))
    # each allele
    for i in range(len(all_chrlen)):
        fa_f = fa_prefix + str(this_leaf_index) + "_" + str(i + 1) + ".fa"
        # sequence the reads using wgsim
        out_fq_prefix = "leaf" + str(this_leaf_index) + "_allele" + str(i)
        out_fq1 = dir + "/" + out_fq_prefix + "_1.fq"
        out_fq2 = dir + "/" + out_fq_prefix + "_2.fq"
        args = "samtools faidx " + fa_f
        popen = subprocess.Popen(args, stdout=subprocess.PIPE, shell=True)
        popen.wait()
        #popen.terminate()
        #output = popen.stdout.read()
        #print output

        # each chromosome
        for j in range(len(all_chrlen[i])):
            print(j)
            this_chrlen = all_chrlen[i][j]
            readcounts = gen_readcount(dir, Alpha, Beta, x0, y0, cov, l, window_size, u, all_chrlen[i][j])
            # each bin
            start = 0
            end = 0
            for k in readcounts:
                this_readcount = int(k/2)
                # extract this segment of fa
                start = end + 1
                end = start + window_size - 1
                print(end)
                if end > this_chrlen:
                    end = this_chrlen
                if start > this_chrlen:
                    break
                tmp_fa_file = "_".join([fa_f, str(chr_name_array[j]), str(start), str(end)]) + ".fa" 
                # use samtools faidx to get this sequence that is to be sequenced into reads
                args = "samtools faidx " + fa_f + " " + chr_name_array[j] + ":" + str(start) + "-" + str(end) + " > " + tmp_fa_file
# TODO check how many N's 
# conclusion: Popen cannot exit. check_call can. But check_call does not work on all Ns as wgsim does not work on it. 
                print(args)
                try:
                    subprocess.check_call(args, shell=True)
                except subprocess.CalledProcessError:
                    print(e)
                    #error_out("Cannot work on '%s'" % args) 

                # check N's
                N_true = check_Ns(tmp_fa_file)

                #popen = subprocess.Popen(args, shell=False)
                # stdout=subprocess.PIPE, shell=True)
                #popen.wait()
                #popen.terminate()
                #output = popen.stdout.read()
                #print output

                # sequence the reads using wgsim
                #out_fq_prefix = "leaf" + str(this_leaf_index) + "_allele" + str(i) + "_chr" + str(j) + "_" + str(start) + "_" + str(end)
                #out_fq1 = out_fq_prefix + "_1.fq"
                #out_fq2 = out_fq_prefix + "_2.fq"
                if N_true == 1:
                    args = wgsim_dir + "wgsim -h -N " + str(this_readcount) + " -1 " + str(l) + " -2 " + str(l) + " " + tmp_fa_file + " " + out_fq1 + " " + out_fq2
                    print(args)
                    try:
                        subprocess.check_call(args, shell=True)
                    except subprocess.CalledProcessError as e:
                        print(e)
                        #error_out("Cannot work on '%s'" % args) 

                args = "rm " + tmp_fa_file 
                try:
                    subprocess.check_call(args, shell=True)
                except subprocess.CalledProcessError as e:
                    print(e)
                    #error_out("Cannot remove '%s'" % tmp_fa_file)
                #popen = subprocess.Popen(args)
                #, stdout=subprocess.PIPE, shell=True)
                #popen.wait()
                #popen.communicate()
                #output = popen.stdout.read()
                #print output
                #popen.terminate()
    print("Done with node " + str(this_leaf_index) + ", will remove it. ")
    
    ref_files = fa_prefix + str(this_leaf_index) + "_*.fa"
    args = "rm " + ref_files 
    try:
        subprocess.check_call(args, shell=True)
    except subprocess.CalledProcessError:
        print(e)
        #error_out("Cannot remove '%s'" % ref_files)


if len(sys.argv) <= 1:
    print("""
    A single cell simulator generating low coverage data. The program automatically generates a phylogenetic tree with copy number variations on the branches. On each leave of the tree, it generates the reads whose error profile, such as uneven coverage, mimics the real single cell data. 
    Usage: python main.py -t [ref.fa] -n [number_leafs] -S [wgsim_dir]
        -p (--processors)   Numbers of processors available.
        -r (--directory)    Location of simulated data. The program will remove the whole directory if it already exists. Otherwise it will create one. (default: test)
        -S (--wgsim-dir)    The directory of the binary of wgsim. It is in the same folder of this main.py. (need to specify) 
        -n (--leaf-num)     Number of the leafs. 
        -B (--Beta)         The program uses the Beta-splitting model to generate the phylogenetic tree. Specify a value between [0, 1]. (default: 0.5)
        -A (--Alpha)        The Alpha in Beta-splitting model. Specify a value between [0, 1]. The closer Alpha and Beta, the more balanced the tree. (default: 0.5).
        -D (--Delta)        The rate of a node to disappear. Specify a value between [0, 1]. If all nodes have daughter nodes, take 0. (default: 0)
        -c (--cn-num)       The average number of copy number variations to be added on a branch. (default: 1)
        -d (--del-rate)     The rate of deletion as compared to amplification. (default: 0.5)
        -m (--min-cn-size)  Minimum copy number size. (default: 200,000bp)
        -e (--exp-theta)    The parameter for the Exponential distribution for copy number size, beyond the minimum one. (default: 0.000001)
        -a (--amp-p)        The parameter for the Genometric distribution for the number of copies amplified. (default: 0.5)
        -t (--template-ref) The reference file to sequence the reads. 
        -o (--outfile)      The standard output file, will be saved in output folder, just give the file name. (default: std.out)
        -f (--fa-prefix)    The prefix of the alleles and read names. (default: ref)
        -x (--Lorenz-x)     The value on the x-axis of the point furthest from the diagonal on the Lorenz curve imitating the real coverage uneveness. (default: 0.5) 
        -y (--Lorenz-y)     The value on the y-axis of the Lorenz curve imitating the real coverage unevenness. x > y. The closer (x, y) to the diagonal, the better the coverage evenness. (default: 0.4) 
        -v (--coverage)     The average coverage of the sequence. (default: 0.02)
        -l (--readlen)      Read length for each read sequenced. (default: 35bp)
        -w (--window-size)  Within a window, the coverage is according to a Gaussian distribution. Neighboring windows' read coverage is according to a Metropolis Hasting process. (default: 200000bp)
        -u (--acceptance-rate)  The probability to accept a proposal in Metropolis Hasting. (default: 0.5)
        -k (--skip-first-step)  If the alleles for all nodes have been made, the step can be skipped. Make it 1 then. (default: 0)
        -R (--snv-rate)     The rate of the snv. snv-rate * branch-length = # snvs. (default: 1)
        -X (--multi-root)   The multiplier of the mean CNV on root. (default: 4)
        -W (--whole-amp)    If there is whole chromosome amplification, 1 as yes. (default: 1) 
        -C (--whole-amp-rate)   Whole amplification rate: rate of an allele chosen to be amplified (default: 0.2)
        -E (--whole-amp-num)    Whole amplification copy number addition, which occurs to one allele at a time. (default: 1)
        -J (--amp-num-geo-par)  Whole amplification copy number distribution (geometric distribution parameter: the smaller, the more evenly distributed). (default: 1)
        -Y (--leaf-index-range) For parallele job submission. >= min, < max leaf index will be processed. min.max (default: -1)
        -I (--leaf-ID-range) For parallele job submission. >= min, < max leaf ID will be processed. min.max (default: -1). When both -Y and -I are -1, all leaves will be processed.
        """)
    sys.exit(0)

parser = argparse.ArgumentParser(description='A copy number variation simulator for single cells.')
parser.add_argument('-p', '--processors', default=8)
parser.add_argument('-r', '--directory', default="test")
parser.add_argument('-S', '--wgsim-dir', default="")
parser.add_argument('-n', '--leaf-num', default=4)
parser.add_argument('-B', '--Beta', default=0.5)
parser.add_argument('-A', '--Alpha', default=0.5)
parser.add_argument('-D', '--Delta', default=0)
#parser.add_argument('-O', '--Output', default="test")
parser.add_argument('-c', '--cn-num', default=1)
parser.add_argument('-d', '--del-rate', default=0.5)
parser.add_argument('-m', '--min-cn-size', default=200000)
#parser.add_argument('-m', '--min-cn-size', default=10)
parser.add_argument('-e', '--exp-theta', default=0.000001)
parser.add_argument('-a', '--amp-p', default=0.5)
parser.add_argument('-t', '--template-ref', default="hg19.fa")
#parser.add_argument('-t', '--template-ref', default="chr1.fa")
parser.add_argument('-o', '--outfile', default="std.out")
parser.add_argument('-f', '--fa-prefix', default="ref")
parser.add_argument('-x', '--Lorenz-x', default=0.5)
parser.add_argument('-y', '--Lorenz-y', default=0.4)
parser.add_argument('-v', '--coverage', default=0.02)
#parser.add_argument('-v', '--coverage', default=10)
parser.add_argument('-l', '--readlen', default=35)
parser.add_argument('-w', '--window-size', default=200000)
#parser.add_argument('-w', '--window-size', default=20)
parser.add_argument('-u', '--acceptance-rate', default=0.5)
parser.add_argument('-k', '--skip-first-step', default=0)
parser.add_argument('-R', '--snv_rate', default=0)
parser.add_argument('-X', '--multi-root', default=4)
parser.add_argument('-W', '--whole-amp', default=0)
parser.add_argument('-C', '--whole-amp-rate', default=0.2)
parser.add_argument('-E', '--whole-amp-num', default=1)
parser.add_argument('-J', '--amp-num-geo-par', default=1)
parser.add_argument('-Y', '--leaf-index-range', default="-1")
parser.add_argument('-I', '--leaf-ID-range', default="-1")


args = parser.parse_args()
NUM_OF_PROCESSES = int(args.processors)
skip = int(args.skip_first_step)
dir = args.directory
save_prefix = dir + "/" + "from_first_step" 
wgsim_dir = args.wgsim_dir
n = int(args.leaf_num)
Beta = float(args.Beta)
Alpha = float(args.Alpha)
Delta = float(args.Delta)
#Output = args.Output
cn_num = int(args.cn_num)
del_rate = float(args.del_rate)
min_cn_size = int(args.min_cn_size)
exp_theta = float(args.exp_theta)
amp_p = float(args.amp_p)
template_ref = args.template_ref
outfile = dir + "/" + args.outfile
fa_prefix = dir + "/" + args.fa_prefix
x0 = float(args.Lorenz_x)
y0 = float(args.Lorenz_y)
cov = float(args.coverage)
l = int(args.readlen)
window_size = int(args.window_size)
u = float(args.acceptance_rate)
#sigma = args.sigma_for_Gaussian
snv_rate = float(args.snv_rate)
root_mult = int(args.multi_root)
whole_amp = int(args.whole_amp)
whole_amp_rate = float(args.whole_amp_rate)
whole_amp_num = int(args.whole_amp_num)
amp_num_geo_par = float(args.amp_num_geo_par)
leaf_index_range = args.leaf_index_range
index_min = 0
index_max = 10000000
if leaf_index_range != "-1":
    index_min, index_max = leaf_index_range.split(".")
    index_min = int(index_min)
    index_max = int(index_max)
leaf_ID_range = args.leaf_ID_range
leaf_index_min = 0
leaf_index_max = 10000000
if leaf_ID_range != "-1":
    leaf_index_min, leaf_index_max = leaf_ID_range.split(".")
    leaf_index_min = int(leaf_index_min)
    leaf_index_max = int(leaf_index_max)

if skip == 0: 
    if not os.path.exists(dir):
        subprocess.check_call("mkdir " + dir, shell=True)
    else:
        print(dir + " exists. Exit. ")
        sys.exit(0)

leaf_chrlen = []
leaf_index = []
chr_name_array = []

def check_Ns(file):
    N_line = 0
    total_line = 0
    with open(file, "r") as f:
        for line in f:
            total_line = total_line + 1
            if line.find('N') != -1:
                N_line = N_line + 1
    f.close()
    if N_line/float(total_line) > 0.2:
        return 0
    return 1



# Step 1. generate a phylogentic tree that has copy number alterations on branches.
# Now all CNs are in tree, not generating fa or remember it in the tree nodes.
if skip == 0:
    [leaf_chrlen, leaf_index, chr_name_array, tree] = gen_tree(n, Beta, Alpha, Delta, dir, cn_num, del_rate, min_cn_size, exp_theta, amp_p, template_ref, outfile, fa_prefix, snv_rate, root_mult, whole_amp, whole_amp_rate, whole_amp_num, amp_num_geo_par)
    numpy.save(save_prefix + ".leaf_chrlen.npy", leaf_chrlen)
    numpy.save(save_prefix + ".leaf_index.npy", leaf_index)
    numpy.save(save_prefix + ".chr_name_array.npy", chr_name_array)
    # save the tree for parallele job submission afterwards
    numpy.save(save_prefix + ".tree.npy", tree)
    #numpy.save(save_prefix + ".ref.npy", ref)
    print("Done with generating the tree. Save to npy. ")
#print leaf_index

# Step 2. Consider even coverage, use metropolis hasting to sample read count in each bin based on a given point on Lorenz curve. 

# if first step was skipped, read the previous stored file
if skip == 1:
    print("Skip the first step. Reading ")
    #abs_path = os.getcwd()
    leaf_chrlen_f = save_prefix + ".leaf_chrlen.npy"
    leaf_index_f = save_prefix + ".leaf_index.npy"
    chr_name_array_f = save_prefix + ".chr_name_array.npy"
    tree_f = save_prefix + ".tree.npy"
    #ref_f = save_prefix + ".ref.npy"
    leaf_chrlen = numpy.load(leaf_chrlen_f)
    leaf_index = numpy.load(leaf_index_f)
    chr_name_array = numpy.load(chr_name_array_f)
    tree = numpy.load(tree_f, allow_pickle=True)
    [ref, tmp_chr_name, tmp_len_chr] = init_ref(template_ref)
    #ref = numpy.load(ref_f)
    #print(leaf_chrlen_f)
    #print(leaf_index_f)
    #print(chr_name_array_f)

    # make it either making a tree, or generating the leaves
    [Alpha, Beta] = get_beta_dist(x0, y0)
    index = 0
    leaf_index_ = 0
    print("Alpha = %.2f, Beta = %.2f", Alpha, Beta)
    print("Number of processes: " + str(NUM_OF_PROCESSES))
    print("Number of leaves: " + str(n))
    # each leaf
    processes = []
    for all_chrlen in leaf_chrlen:
    #TODO add a function here to retrieve the CN and apply to reference, and write to a fa file for the two alleles.
        leaf_index_ = leaf_index[index]
        if leaf_index_range != "-1" and index < index_max and index >= index_min or leaf_ID_range != "-1" and leaf_index_ < leaf_index_max and leaf_index_ >= leaf_index_min or leaf_index_range == "-1" and leaf_ID_range == "-1": 
            make_fa(leaf_index[index], tree, ref, chr_name_array, fa_prefix)
            
            #if index % NUM_OF_PROCESSES == 0 and index != 0:
            
                # Run processes
                #for p in processes:
                #    p.start()
                
                # Exit the completed processes
                #for p in processes:
                #    p.join()
                
                # clean for the next batch
                #processes = []
        
                # clean the temp fa files
                #args = "rm " + fa_prefix + "*_*_*_*_*.fa"
                #popen = subprocess.Popen(args, stdout=subprocess.PIPE, shell=True)
                #popen.wait()
        
            processes.append(mp.Process(target=gen_reads, args=(dir, index, leaf_index, all_chrlen, fa_prefix, Alpha, Beta, x0, y0, cov, l, window_size, u, chr_name_array)))
    
        index = index + 1

    # clean the rest
    for p in processes:
        p.start()
    
    for p in processes:
        p.join()
    
    # clean up the temp fa files
    #args = "rm " + fa_prefix + "*_*_*_*_*.fa"
    #popen = subprocess.Popen(args, stdout=subprocess.PIPE, shell=True)
    #popen.wait()
