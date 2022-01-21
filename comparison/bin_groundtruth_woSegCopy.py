import re
import sys, argparse

if len(sys.argv) <= 1:
    print("""
    Given the ground truth file (gt.all.csv) and segcopy file, generate the ground truth binned file in the SegCopy file format (chr, start, end, copynumber for cell1, 2, 3 ...). This will be the input (from ground truth or from inferred results) to a phylogenetic tree inference algorithm for comparison.
    Usage: python bin_groundtruth_colNotSpecified.py -a [segcopy_file] -b [ground_truth_all] --leafonly -c [chr] -l [leaf_list_file, mandatory if leafonly]
    """)
    sys.exit(0)

parser = argparse.ArgumentParser(description='Make segcopy formatted file from ground truth file, with the columns from the segcopy file.')
parser.add_argument('-a', '--segcopyfile', default="NA")
parser.add_argument('-b', '--gtfile', default="NA")
parser.add_argument('--leafonly', default=False, action='store_true')
parser.add_argument('-c', '--CHR', default="NA")
parser.add_argument('-l', '--leaves', default="NA")

args = parser.parse_args()

gt = args.gtfile
segcopy = args.segcopyfile 
leafonly = args.leafonly
CHR = args.CHR
leaf_list_f = args.leaves

#segcopy = "/storage/hpc/work/nakhleh/xf2/benchmark/large_dataset/rep_2/ginkgo/SegCopy.p02"
#gt = "/storage/hpc/work/nakhleh/xf2/benchmark/large_dataset/rep_2/misc/gt.all.csv"

file_a = open(segcopy, "r")
# read all the coordinates and the leaves from file a
tag = 0
leaves = {}
pos = {}
leaves_order = []
chrs_order = []
for line in file_a:
    arrays = re.split(r'\t+', line.rstrip()) 
    chr_ = arrays[0]
    if arrays[0] == "CHR":
        continue
    if CHR != "NA":
        # chromosome specified, record the position only on this chromosome
        if CHR == chr_:
            if chr_ not in pos.keys():
                pos[chr_] = []
                chrs_order.append(chr_)
            pos[chr_].append(arrays[1] + "." + arrays[2])
    else:
        # chromosome not specified, use all chromosomes
        if chr_ not in pos.keys():
            pos[chr_] = []
            chrs_order.append(chr_)
        pos[chr_].append(arrays[1] + "." + arrays[2])
file_a.close()

if leafonly: 
    # read the leaf list file
    file_l = open(leaf_list_f, "r")
    for line in file_l:
        leaf_str = "leaf" + line.rstrip()
        leaves[leaf_str] = 1
        leaves_order.append(leaf_str)
    file_l.close()
else:
    # read all leaves in the groundtruth file
    file_gt = open(gt, "r")
    for line in file_gt:
        arrays = re.split(r'\t+', line.rstrip())
        [chr, s, e, cn, leaf_id] = arrays
        leaf_str = "leaf" + leaf_id
        if leaf_str not in leaves:
            leaves[leaf_str] = 1
            leaves_order.append(leaf_str)
    file_gt.close()

        
# convert file b to the format of a
file_b = open(gt, "r")
# normalize the whole matrix
mat = {}
for chr in pos.keys():
    if chr not in mat.keys():
        mat[chr] = {}
    for bin in range(len(pos[chr])):
        if bin not in mat[chr].keys():
            mat[chr][bin] = {}
        for leaf in leaves.keys():
            mat[chr][bin][leaf] = 2

for line in file_b:
    arrays = re.split(r'\t+', line.rstrip())
    [chr, s, e, cn, leaf_id] = arrays
    s = int(s)
    e = int(e)
    cn = int(cn)
    num = 0
    bin_s = 0
    bin_e = 0
    sel_s = 0
    sel_e = 0
    if "leaf" + leaf_id in leaves:
        for i in pos[chr]:
            (bin_s, bin_e) = re.split(r'\.', i)
            bin_s = int(bin_s)
            bin_e = int(bin_e)
            if s >= bin_s  and s < bin_e:
                #if abs(s - bin_s) < abs(s - bin_e):
                #    sel_s = num
                #else:
                #    sel_s = num + 1
                sel_s = num
            if e >= bin_s  and e < bin_e:
                #if abs(e - bin_s) < abs(e - bin_e):
                #    sel_e = num - 1
                #else:
                #    sel_e = num
                sel_e = num - 1
                for bin in range(sel_s, sel_e + 1):
                    #print chr + "," + str(bin) + "," + leaf_id
                    #if bin not in mat[chr].keys():
                    #    bin -= 1
                    if bin in mat[chr].keys():
                        mat[chr][bin]["leaf" + leaf_id] = cn

            num += 1

        if e >= bin_e and s <= bin_e:
            for bin in range(sel_s, num + 1):
                if bin in mat[chr].keys():
                    mat[chr][bin]["leaf" + leaf_id] = cn
                    
file_b.close()

# print result
print "CHR\tSTART\tEND",
for leaf in leaves_order:
    print "\t" + leaf,
print ""

for chr in chrs_order:
    num = 0
    for bin in pos[chr]:
        [s, e] = re.split(r'\.', bin)
        print chr + "\t" + str(s) + "\t" + str(e),
        for leaf in leaves_order:
            print "\t" + str(mat[chr][num][leaf]),
        num += 1
        print ""


