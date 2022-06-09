#!/usr/bin/python$
import numpy
import copy
from gen_tree import gen_tree 
from Gen_Ref_Fa import getlen_ref
from gen_tree import get_cn_from_corres
import CN
import argparse
import sys
import re

class mynode:
    def __init__(self, ID, parent=None):
        self.ID = ID
        # ID
        self.parent = parent
        self.children = []

# print which level has how many and which nodes 
def print_levels(tree_element_f):
    tree_elements_arr = numpy.load(tree_elements_f, allow_pickle=True)
    tree_elements = tree_elements_arr[0]
    level_indice = tree_elements.level_indice
    # level_indice is a dictionary pointing from level to the node indices in an array
    for i in range(len(level_indice)):
        nodes = level_indice[i]
        if i == 0:
            print("Skip level " + str(i) + " which represents the normal cell. ")
        else:
            print("On level " + str(i) + ", there are " + str(len(nodes)) + " nodes, and they are " + ",".join([str(k) for k in nodes]))

# print all snvs if no node is specified
def print_snvs(tree, cell):
    print("\t".join(["chr", "pos", "ori_nuc", "new_nuc", "allele", "cell"]))
    cells = []
    if cell == "NA":
        for i in range(len(tree)):
            cells.append(i)
    else:
        array = re.split(r';', cell)
        for i in array:
            cells.append(int(i))
    
    for i in cells:
        snv = tree[i].snvs
        for s in snv:
            print("\t".join([str(s.chr), str(s.ref_pos), s.ori_nuc, s.new_nuc, str(s.ale), str(i)]))
    

# a new function added in Feb. 2020
# output a matrix, in which the entry of (i, j) represents the number of events between these two profiles (here assume i and j are both leaves)
def get_pairwise_diff_matrix(tree):
    leaves = get_leafid_array(tree)
    for i in range(len(leaves)):
        for j in range(len(leaves)):
            x = get_pairwise_diff(tree, leaves[i], leaves[j])
            print(str(x) + "\t", end="")
        print("")
    

# a new function added in Feb. 2020
# for two leaf cells, output their difference (different event #)
def get_pairwise_diff(tree, i, j):
    p = LCA(tree, i, j)
    #print "LCA of " + str(i) + " and " + str(j) + " is " + str(p)
    e_num_i = get_event_num_path(tree, i, p)
    e_num_j = get_event_num_path(tree, j, p)
    return e_num_i + e_num_j
    
# augmentary function
# given two nodes, get the least common ancestor
def LCA(tree, i, j): 
    P_i = get_path(tree, i)
    P_j = get_path(tree, j)
    for x in P_i:
        for y in P_j:
            if x == y:
                return x
     
# augmentary function
# given the tree and a node, find the path bottom up till it reaches the root
def get_path(tree, i):
    p = tree[i].parent.id
    P = [i]
    while p != -1:
        P.append(p)
        #print p
        p = tree[p].parent.id
    return P 
 
# augmentary funciton
# given a tree, two nodes, one is the ancestor of the other, calculate the total event number from the ancestor to the node
def get_event_num_path(tree, i, p):
    cn_array_num = 0
    while i != p:
        cn_array_num += len(tree[i].cn)
        i = tree[i].parent.id
    return cn_array_num

# a new function added in Jan. 2020
# for each daughter cell compared with a parent cell, output its new CNAs
def get_event_num(tree):
    for t in tree:
        ID = t.id
        cn_array = t.cn
        print(str(ID) + "\t" + str(len(cn_array)) + "\t" + str(t.if_leaf))


def get_leaf(tree):
    for t in tree:
        if t.if_leaf:
            print(str(t.id))

def get_leafid_array(tree):
    a = []
    for t in tree:
        if t.if_leaf:
            a.append(t.id)
    return a

def retrieve_new_CNAs(tree):
    s = {}
    S = {}
    for t in tree:
        ID = t.id
        s[ID] = str(ID)
        S[ID] = str(ID) 
    # now count
    for t in tree:
        ID = t.id
        # s contains only new CNAs to this node
        s[ID] = {}
        # S contains all
        S[ID] = {}
        summary = t.cn_summary
        p = t.parentID        
        for chr in summary:
            for pos in summary[chr]:
                chr_ = chr + 1
                if chr_ == 23:
                    chr_ = "X"
                if chr_ == 24:
                    chr_ = "Y"
                str_ = str(chr_) + ":" + str(pos) + "_" + str(summary[chr][pos])
                if p == -1:
                    s[ID][str_] = 1
                elif str_ not in S[p]:
                    s[ID][str_] = 1
                S[ID][str_] = 1
    return s

# given a segcopy file with all nodes, a tree with the parent-children relationship, output the new CNAs of a child compared to the parent, including internal nodes
def retrieve_new_overlappingCNAs(segcopy_f, Tree):
    f = open(segcopy_f, "r")
    line = f.readline().rstrip("\n")
    names = []
    # from leaf id to the column index
    names_h = {}
    # h is from the position to the leaf index (0-based) to CN
    h = {}
    first = True
    while(line != ""):
        array = re.split(r'\s+', line)
        if first:
            names = array[3:]
            for i in range(len(names)):
                if names[i] == "":
                    break
                names_h[int(names[i][4:])] = i
            first = False
        else:
            cnas = array[3:]
            loc = str(array[0]) + ":" + str(array[1]) + "-" + str(array[2])
            h[loc] = {}
            #locs.append(loc)
            for i in range(len(cnas)):
                if cnas[i] == "":
                    break
                h[loc][i] = int(cnas[i])
        line = f.readline().rstrip("\n")
    f.close()
    # now retrieve the new CNAs comparing child's with parent's CNA profile
    # new_cna_h is from the leaf id to chromosome to positions to CNA (+: increase, -: decrease) with the absolute copy number change
    new_cna_h = {}
    for i in range(len(Tree)):
        p = Tree[i].parentID
        new_cna_h = {}
        i_col_id = names_h[i]
        if p == -1:
            for k in list(h.keys()):
                chr, pos = re.split(r':', k)
                s, e = re.split(r'-', pos)
                s = int(s)
                if chr not in list(new_cna_h.keys()):
                    new_cna_h[chr] = {}
                if h[k][i_col_id] != 2:
                    new_cna_h[chr][s] = e + ";" + str(h[k][i_col_id] - 2)
        else:
            p_col_id = names_h[p]
            for k in list(h.keys()):
                diff = h[k][i_col_id] - h[k][p_col_id] 
                if diff != 0:
                    chr, pos = re.split(r':', k)
                    s, e = re.split(r'-', pos)
                    s = int(s)
                    if chr not in list(new_cna_h.keys()):
                        new_cna_h[chr] = {}
                    new_cna_h[chr][s] = e + ";" + str(diff)
        combined = combine_cnas(new_cna_h)
        for l in combined:
            a = copy.deepcopy(l)
            a.append(str(i))
            a.append(str(p))
            print("\t".join(a))

# connect the CNA bins together, return an array with chr, start, end
def combine_cnas(h):
    h_ret = []
    for chr in list(h.keys()):
        start = "NA"
        end = "NA"
        prev_e = "NA"
        prev_cn = "NA"
        interval = 0
        for ss in sorted(list(h[chr].keys())):
            s = str(ss)
            e, cn = re.split(r';', h[chr][ss])
            if start == "NA":
                start = s
                end = e
                prev_e = e
                prev_cn = cn
                interval = 0
            elif int(s) == int(prev_e) + 1 and prev_cn == cn:
                interval += 1
                end = e
                prev_e = e
            else:
                # temporarily not look at the actual cn 
                if interval > 1:
                    h_ret.append([chr, start, end, prev_cn])
                start = s
                end = e
                prev_e = e
                prev_cn = cn
                interval = 0

        if interval > 1:
            h_ret.append([chr, start, end, cn])
    return h_ret
   
# processed are those that have already to poped up from the bucket 
def convert2newick_general(Tree, str_, bucket):
    if len(bucket) == 0:
        return str_
    splitted = re.split(r'(\d+)', str_) 
    new_bucket = []
    for i in bucket:
        for j in range(len(splitted)):
            if str(i) == splitted[j]:
                # replace it with children
                splitted[j] = "(" + ",".join([str(x) for x in Tree[i].children]) + ")" + str(i)
                for k in Tree[i].children:
                    if len(Tree[k].children) != 0:
                        new_bucket.append(k)
    
    return convert2newick_general(Tree, "".join(splitted), new_bucket)
 
# print the tree in newick format, level is to control when to stop to dig into the tree (in which case only leaves in that cluster and the head node will be printed)
def convert2newick(Tree, str_, i):
    # stop
    if i >= len(Tree):
        return str_
    node = Tree[i]
    parent = node.parent
    children = node.children
    # get rid of the node itself
    c_arr = []
    for c in children:
        if c != str(i):
            c_arr.append(c)
    children = c_arr

    if len(children) >= 1:
        #print children
        # prepare for where to replace
        splitted = re.split(r'(\d+)', str_)
        i_ = -1
        for i_ in range(len(splitted)):
            if str(splitted[i_]) == str(i):
                break
        # now make a new string to replace
        new_substr = "(" + "," . join(children) + ")"
        #print str(splitted[i_]) + " is to be replaced by " + new_substr
        # concatenate
        splitted[i_] = new_substr + splitted[i_]
        str_ = "".join(splitted)
        str_ = convert2newick(Tree, str_, i + 1)
    else:
        str_ = convert2newick(Tree, str_, i + 1)
    return str_
        
        
def print_newick(tree):

    # a new tree with a simplified structure
    Tree = []
    # in case a certain parent has not been added, initialize Tree
    for i in range(len(tree)):
        Tree.append(mynode(i, -1)) 
    
    # read the npy file and put them in MyNode
    for i in range(len(tree)):
        Tree[i].parent = tree[i].parent.id
        if tree[i].parent.id != -1:
            Tree[tree[i].parent.id].children.append(i)

    str_ = convert2newick_general(Tree, "0", [0]) 
    #str_ = convert2newick(Tree, "(0)", 0)
    return str_

def print_tree(tree):
    leaf_only = True
    dec = get_descendants(tree, leaf_only, False, "NA") 
    s = retrieve_new_CNAs(tree)
    children = get_children(tree)
    for t in tree:
        ID = t.id
        #if ID == 0:
        #    continue
        #else:
        decs = dec[ID]
        num = len(decs)
        print(str(ID) + "\t" + ";".join(children[ID]) + "\t" + str(num) + "\t" + ";".join(list(s[ID].keys())) + "\t" + ";".join(decs))

def get_summary(tree, select_leaf):
    for t in tree:
        ID = t.id
        if select_leaf:
            if not t.if_leaf:
                continue
        summary = t.cn_summary
        for chr in summary:
            for pos in summary[chr]:
                pos_s, pos_e = pos.split(".")
                cn = summary[chr][pos]
                chr_ = chr + 1
                chr_ = "chr" + str(chr_)
                if chr_ == "chr23":
                    chr_ = "chrX"
                if chr_ == "chr24":
                    chr_ = "chrY"
                print("\t".join([chr_, str(pos_s), str(pos_e), str(cn), str(ID)])) 

def print_large_clusters(tree, leaf_only, size):
    d = get_descendants(tree, leaf_only, True, size)
    for i in d:
        if len(d[i]) != 0:
            y = [str(x) for x in d[i]]
            print(" ".join(y))

def retrieve_new_CNAs(tree):
    for i in range(len(tree)):
        t = tree[i]
        tCNs = t.true_CNs
        for tCN in tCNs: 
            for k in list(tCN.keys()):
                chr, interval = k.split(":")
                s, e = interval.split(".")
                print("\t".join([chr, s, e, str(tCN[k]), str(i), str(t.parentID)])) 
    return
        
# for each node in the tree, find the descendants of it
# in a reverse order, traverse the tree. For every node, add the descendants of its daughter cells, and its daughter cells into its descendant list. 
# if leaf_only is on, output only leaf as the descendants; 
# if cut_by_size is on, then remove the descendants if the parent's descendants do not fall into the size range
def get_descendants(tree, leaf_only, cut_by_size, size):
    d = {}
    # intiialize the dictionary
    for t in tree[::-1]:
        ID = t.id
        if leaf_only:
            if t.if_leaf:
                d[ID] = [str(ID)]
            else:
                d[ID] = []
        else:
            d[ID] = [str(ID)]
    # now count
    for t in tree[::-1]:
        ID = t.id
        p = t.parentID
        if p != -1:
            for i in d[ID]:
                d[p].append(i) 

    # cut by size
    if cut_by_size:
        min_, max_ = size.split(",")
        for t in tree:
            if len(d[t.id]) < int(min_) or len(d[t.id]) > int(max_):
                d[t.id] = []
                
    return d

# return a dictionary with the keys the ids, the values an array of children IDs 
def get_children(tree):
    d = {}
    for t in tree:
        ID = t.id
        d[ID] = []
    for t in tree:
        ID = t.id
        p = t.parentID
        if p != -1:
            d[p].append(str(ID))
    return d

# in case the tree does not have the cn summary 
def make_summary_func(tree, ref):
    tmp_name_array, chr_sz = getlen_ref(ref)
    for i in tree:
        i.cn_detail, i.cn_summary = get_cn_from_corres(i.corres, chr_sz)
        ID = i.id
        for chr in i.cn_summary:
            for pos in i.cn_summary[chr]:
                pos_s, pos_e = pos.split(".")
                cn = i.cn_summary[chr][pos]
                chr_ = chr + 1
                chr_ = "chr" + str(chr_)
                if chr_ == "chr23":
                    chr_ = "chrX"
                if chr_ == "chr24":
                    chr_ = "chrY"
                print("\t".join([chr_, str(pos_s), str(pos_e), str(cn), str(ID)])) 


parser = argparse.ArgumentParser(description='Read a tree and output specific items of it. ')
parser.add_argument('-l', '--leaf', action='store_true')  
parser.add_argument('-s', '--summary', action='store_true')  
parser.add_argument('-L', '--select-leaf', action='store_true')  
parser.add_argument('-m', '--make-summary', action='store_true')
parser.add_argument('-r', '--ref', default="")
parser.add_argument('-f', '--file', default="") 
parser.add_argument('-e', '--event', action='store_true')
parser.add_argument('-d', '--pairdist', action='store_true')
parser.add_argument('-P', '--printtree', action='store_true')
parser.add_argument('-N', '--printnewick', action='store_true')
parser.add_argument('-C', '--printlargecluster', action='store_true')
parser.add_argument('-S', '--largeclustersize', default="8,10")
parser.add_argument('-O', '--retrieveoverlapping', action='store_true')
parser.add_argument('-F', '--segcopyf', default="NA")
parser.add_argument('-o', '--retrievealloverlaps', action='store_true')
parser.add_argument('-n', '--printsnvs', action='store_true')
parser.add_argument('-x', '--printsnvsforcell', default="NA")
parser.add_argument('-v', '--printlevels', action='store_true')
parser.add_argument('-E', '--treeelementsf', default="NA")


args = parser.parse_args()
if_leaf = args.leaf
event_num = args.event
pairdist = args.pairdist
if_summary = args.summary
select_leaf = args.select_leaf
make_summary = args.make_summary
printtree = args.printtree
ref_f = args.ref
npy_f = args.file
printnewick = args.printnewick
printlargecluster = args.printlargecluster
cluster_size = args.largeclustersize
retrieveoverlapping = args.retrieveoverlapping
segcopy_f = args.segcopyf
retrievealloverlaps = args.retrievealloverlaps
printsnvs = args.printsnvs
printsnvsforcell = args.printsnvsforcell
printlevels = args.printlevels
tree_elements_f = args.treeelementsf

# main starts here
if npy_f == "": 
    print("""
    Given a tree in npy format, output its leaf index or the CNV summary. 
    Usage: python read_tree.py -l -s -e -d -o -f [tree.npy] -E [tree_elements.npy]
        -l  (--leaf)    Print leaf index, one per line. (default: off)
        -s  (--summary) Print CNV summary, one per line (chr, start, end, CN). (default: off)
        -L  (--selectleaf)  select to print leaf, in conjunction with if_summary. (default: off)
        -m  (--make-summary)    When the tree does not have CNV summary info but the correspondence, make it and print the CNV info, like -s. Need the -r info to retrieve the info back. (default: off)
        -r  (--ref)     Reference file in .fa. Necessary only when -m is turned on. 
        -f  (--file)    The npy file storing the tree. (mandatory)
        -e  (--event)   Get the number of events (col 2) for each node (col 1) and indicate whether it is a leaf node or not (col 3). 
        -d  (--pairdist)    Get the pairwise distance between all leaf cells. Output the matrix to the stdout. (default: off) 
        -P  (--printtree)   Print the tree with each line representing a new CNA, with the columns node, the new CNA to this node and its descendants including itself.
        -N  (--printnewick) Print the newick formatted tree.
        -C  (--printlargecluster)   Print large clusters.
        -S  (--largeclustersize)   Size range separated by comma.
        -O  (--retrieveoverlapping) Retrieve new overlapping CNAs for each cell (including internal nodes). 
        -F  (--segcopyf)    Segcopy file for retrieving new overlapping CNA. 
        -o  (--retrievealloverlaps)  Retrieve new overlapping CNAs, even for those occurring on the same edge from true_CNs in the tree. Output the new CNAs for each node compared with its parent node. 
        -n  (--printsnvs)   Print SNVs for all cells.
        -x  (--printsnvsforcell)    Print SNVs for a particular set of nodes separated by colon specified here. 
        -v  (--printlevels) Print which level has which nodes.
        -E  (--tree_elements_f) Tree elements file is used when -v is on. 
        """)
    sys.exit(0)


tree = numpy.load(npy_f, allow_pickle=True)
if if_leaf:
    get_leaf(tree)

if if_summary:
    get_summary(tree, select_leaf) 

if make_summary:
    make_summary_func(tree, ref_f)

if event_num:
    get_event_num(tree)

if pairdist:
    get_pairwise_diff_matrix(tree)
    
if printtree:
    print_tree(tree)

if printnewick:
    print(print_newick(tree))

if printlargecluster:
    print_large_clusters(tree, True, cluster_size)

if retrieveoverlapping:
    retrieve_new_overlappingCNAs(segcopy_f, tree)

if retrievealloverlaps:
    retrieve_new_CNAs(tree)

if printsnvs:
    print_snvs(tree, "NA")

if printsnvsforcell != "NA":
    print_snvs(tree, printsnvsforcell) 

if printlevels:
    print_levels(tree_elements_f)


