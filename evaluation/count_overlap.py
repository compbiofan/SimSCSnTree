#!/usr/bin/python
import sys
import re
import os
import argparse
import numpy as np
sys.path.append(os.getenv("HOME") + "/python_modules")
from tree import add_children_from_MyNode, calculate_overlap, add_ancestor_edges, edge_statistics, edge_statistics_chro_focal
if len(sys.argv) <= 1:
    print("""
        Calculate the number of overlapping CNAs and base pairs given a tree.npy file. If given the leaf.id file, then evaluate only the leaves. 
        Usage: python count_overlap.py -i [leaf_id] -t [tree] -f [all.gt.cnp]
            -i (--leaf-id) leaf.id from single cell simulator
            -t (--tree) tree.npy 
            -f (--gt-cnp) all.gt.cnp file (chr, s, e, change_cn, node_id, parent_id), in which each node has only its CNAs listed
            -a (--fai-f) the reference fai file to get the chromosome length
            -T (--chromosome-t) threshold above which percentage the event is counted as chromosomal 
    """)
    sys.exit(0)

parser = argparse.ArgumentParser(description="Evaluate overlapping extent of CNAs. ")
parser.add_argument('-i', '--leaf-id', default="NA")
parser.add_argument('-t', '--tree', default="NA")
parser.add_argument('-f', '--gt-cnp', default="NA")
parser.add_argument('-a', '--fai-f', default="NA")
parser.add_argument('-T', '--chromosome-t', default="NA")

args = parser.parse_args()

l = args.leaf_id
t = args.tree
f = args.gt_cnp
ref_fa_f = args.fai_f
T = float(args.chromosome_t)


tree = add_children_from_MyNode(t)

leaves = []
if l != "NA":
    fh = open(l, "r")
    line = fh.readline().rstrip("\n")
    while(line != ""):
        leaves.append(int(line))
        line = fh.readline().rstrip("\n")
    fh.close()
else:
    for i in range(len(tree)):
        leaves.append(i)

fh = open(f, "r")
line = fh.readline().rstrip("\n")
while(line != ""):
    a = re.split(r'\s+', line)
    e = a[0] + ":" + a[1] + "-" + a[2]
    tree[int(a[-2])].edge.append(e)
    line = fh.readline().rstrip("\n")

fh.close()

# add all ancestor's edges to this node
add_ancestor_edges(tree, 0)

ocna_number, ocna_length, number_perc = calculate_overlap(tree, leaves)
cna_number, cna_length = edge_statistics(tree, leaves) 

chro_n, nonchro_n, chro_l, nonchro_l, ov_nonchro_n, ov_nonchro_l, ov_nonchro_p, ov_l = edge_statistics_chro_focal(tree, leaves, ref_fa_f, T)

#print "Median of number of overlapping CNAs: " + str(np.median(np.array(ocna_number)))
#print "Median of overlapping bases: " + str(np.median(np.array(ocna_length)))
print "Average number of overlapping CNAs per cell: " + str(np.average(np.array(ocna_number)))
print "Average overlapping bases per cell: " + str(np.average(np.array(ocna_length)))
print "Average percentage of overlapping CNAs per cell: " + str(np.average(np.array(number_perc)))
print "Average number of CNAs per cell: " + str(np.average(np.array(cna_number)))
print "Average length of CNAs per cell: " + str(np.average(np.array(cna_length)))
print "Average chromosomal CNA number per cell: " + str(np.average(np.array(chro_n)))
print "Average nonchromosomal CNA number per cell: " + str(np.average(np.array(nonchro_n)))
print "Average chromosomal CNA length per cell: " + str(np.average(np.array(chro_l)))
print "Average nonchromosomal CNA length per cell: " + str(np.average(np.array(nonchro_l)))
print "Average overlapping nonchromosomal CNA number per cell: " + str(np.average(np.array(ov_nonchro_n)))
print "Average overlapping nonchromosomal CNA length per cell: " + str(np.average(np.array(ov_nonchro_l)))
print "Average overlapping nonchromosomal CNA percentage per cell: " + str(np.average(np.array(ov_nonchro_p)))
print "Average overlapping CNA bases (nonchromosomal with nonchromosomal, nonchromosomal with chromosomal) per cell: " + str(np.average(np.array(ov_l)))

