import re
import sys, argparse

if len(sys.argv) <= 1:
    print("""
    Given the binning file (format like Ginkgo's SegCopy) file and the ground truth file (gt.all.csv), generate the ground truth binned file in the same format as the first file (a matrix). This will be the input (from ground truth or from inferred results) to a phylogenetic tree inference algorithm for comparison.
    Usage: python bin_groundtruth.py -a [segcopy_f] -b [ground_truth_all] -c [leaf_only]
    """)
    sys.exit(0)

parser = argparse.ArgumentParser(description='Make segcopy formatted file from ground truth.')
parser.add_argument('-a', '--segcopy', default="NA")
parser.add_argument('-b', '--gtfile', default="NA")
parser.add_argument('--leafonly', default=False, action='store_true')

args = parser.parse_args()
segcopy = args.segcopy
gt = args.gtfile
leafonly = args.leafonly

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
    if tag == 0:
        for i in range(3, len(arrays)):
            leaves[arrays[i]] = 1
            leaves_order.append(arrays[i])
        tag = 1
    else:
        chr = arrays[0]
        if chr not in pos.keys():
            pos[chr] = []
            chrs_order.append(chr)
        pos[chr].append(arrays[1] + "." + arrays[2])
file_a.close()

if not leafonly: 
    file_b = open(gt, "r")
    leaves = {}
    leaves_order = []
    for line in file_b:
        arrays = re.split(r'\t+', line.rstrip())
        [chr, s, e, cn, leaf_id] = arrays
        leaf_str = "leaf" + leaf_id
        if leaf_str not in leaves.keys():
            leaves_order.append(leaf_str)
            leaves[leaf_str] = 1
    file_b.close()
        
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
                if abs(s - bin_s) < abs(s - bin_e):
                    sel_s = num
                else:
                    sel_s = num + 1
            if e >= bin_s  and e < bin_e:
                if abs(e - bin_s) < abs(e - bin_e):
                    sel_e = num
                else:
                    sel_e = num + 1
                for bin in range(sel_s, sel_e + 1):
                    #print chr + "," + str(bin) + "," + leaf_id
                    if bin not in mat[chr].keys():
                        bin -= 1
                    mat[chr][bin]["leaf" + leaf_id] = cn

            num += 1

file_b.close()

# print result
print("CHR\tSTART\tEND"),
for leaf in leaves_order:
    print("\t" + leaf),
print("")

for chr in chrs_order:
    num = 0
    for bin in pos[chr]:
        [s, e] = re.split(r'\.', bin)
        print(chr + "\t" + str(s) + "\t" + str(e)),
        for leaf in leaves_order:
            print("\t" + str(mat[chr][num][leaf])),
        num += 1
        print("")


