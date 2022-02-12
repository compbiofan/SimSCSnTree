# rewrite it and wrap it in gen_tree.py
#!/usr/bin/env python2
# -*- coding: utf-8 -*-
# add copy number, each branch has a copy number variation
# assume the genome is diploid, and not go aneuploid
"""
Created on Sat Mar 24 13:59:59 2018

@authors: Xian Fan and Mohammad Edrisi
Contacting email: xf2@rice.edu
"""

# -*- coding: utf-8 -*- 
from anytree import Node, RenderTree
import numpy as np
from anytree.dotexport import RenderTreeGraph
from CN import CN
from Gen_Ref_Fa import gen_ref, init_ref, write_ref

n = 4
Beta = 0.5
Alpha = 0.5
Delta = 0
Output = "test"
cn_num = 1
del_rate = 0.5
#del_rate = 0
min_cn_size = 200000
#min_cn_size = 2000000
# exponential distribution
# smaller exp_theta means larger chance to get larger CNV 
exp_theta = 0.000001
# geometric distribution
# like simulated annealing, lower amp_p means larger chance to get large CN amp
amp_p = 0.5
#template_ref = "ref.fasta"
template_ref = "/home1/03626/xfan/reference/hg19.fa"
outfile = "/work/03626/xfan/lonestar/std.out"
fa_prefix = "/work/03626/xfan/lonestar/ref"

ref_array = []
chr_name_array = []
chr_sz = []
#n = int(raw_input("n:"))
#Beta = float(raw_input("beta:"))
#Alpha = float(raw_input("alpha:"))
#Delta = float(raw_input("delta:"))
#Output = raw_input("output file:")
#cn_num = int(raw_input("mean copy number:"))
#del_rate = float(raw_input("deletion rate [0, 1]:"))
#min_cn_size = int(raw_input("minimum copy number size, recommend > 2000000:"))
#exp_theta = float(raw_input("parameter for copy number size:"))
#amp_p = float(raw_input("parameter for amplification allele #:"))
#template_ref = raw_input("template fasta file:")
#outfile = raw_input("Output file name:")
#fa_f_prefix = raw_input("fasta prefix:")

f = open(outfile, "w")




class MyNode(Node):
    def __init__(self, name, parent=None):
        Node.__init__(self, name, parent)
        self.id = 0
        self.name = name
        self.parent=parent
        self.tuple=[]
        self.is_dead=False
        self.edge_length = 0
        # alelle length for each chromosome, root has the same as reference
        self.cn=[]
        self.chrlen=[]
        self.ref=[]
    def getTuple(self):
        return self.tuple
    def setDead(self):
        self.is_dead=True
    def getID(self):
        return self.id


def get_range(chr_len, min_cn_size, exp_theta):
    # get pos1 and pos2 for the copy number. 
    # skip the first and last "skip" nucleotides (CN does not happen in that range) 
# TODO
    skip = 0
    #skip = 5000000 
    # not try any more if cannot find a good one, just use the end
    trials = 0
    max_n = 10
    cn_size = np.random.exponential(exp_theta) + min_cn_size
    p1 = np.random.uniform(skip, chr_len - skip)
    # make sure p2 not out of boundary
    while chr_len - p1 - skip - cn_size < 0 and trials < max_n:
        p1 = np.random.uniform(skip, chr_len - skip)
        trials = trials + 1
    p2 = p1 + cn_size
    return int(p1), int(p2)

def binary_search(x, array, l, r):
    # find the index i where x < array[i] and x > array[i-1]
    t = l + int((r - l)/2)
    if r > l:
        if x < array[t]:
            r = t
            return binary_search(x, array, l, r)
        else:
            l = t + 1
            return binary_search(x, array, l, r)
    return r

def get_chr(chrlen_array):
    # draw a chromosome with the chance linear to the chromosome length
    acum = []
    for i in range(len(chrlen_array)):
        if i == 0:
            acum.append(chrlen_array[i])
        else:
            acum.append(acum[i-1] + chrlen_array[i])
    # randomly draw a position from the whole genome
    x = np.random.uniform(1, acum[len(chrlen_array)-1])
    return binary_search(x, acum, 0, len(chrlen_array) - 1)

def add_CN(chrlen, cn_num, del_rate, min_cn_size, exp_theta, amp_p):
    # for each branch, the copy number change happens on only one allele
    CN_array = []
    new_chrlen = [row[:] for row in chrlen]
    CN_Tot = np.random.poisson(cn_num, 1)
    for i in range(CN_Tot):
        # allele
        CN_Ale = np.random.binomial(1, 0.5)
        # deletion versus amplification 
        CN_Del = np.random.binomial(1, del_rate)
        CN_chromosome = get_chr(new_chrlen[CN_Ale])
        #print CN_Ale, CN_chromosome
        CN_p1, CN_p2 = get_range(new_chrlen[CN_Ale][CN_chromosome], min_cn_size, exp_theta)
        CN_amp_num = 0
        #print new_chrlen
        #print "before changing:"
        #print chrlen, new_chrlen
        if CN_Del == 0:
            # get amplification copy number
            # starting from 0
            CN_amp_num = int(np.random.geometric(amp_p, 1) - 1)
            #print CN_amp_num, CN_p2, CN_p1
            new_chrlen[CN_Ale][CN_chromosome] = new_chrlen[CN_Ale][CN_chromosome] + CN_amp_num * (CN_p2 - CN_p1)
            #print new_chrlen
        else:
            new_chrlen[CN_Ale][CN_chromosome] = new_chrlen[CN_Ale][CN_chromosome] - (CN_p2 - CN_p1)
            #print new_chrlen
        #print "after changing:"
        #print chrlen, new_chrlen
        CN_array.append(CN(CN_Ale, CN_Del, CN_chromosome, CN_p1, CN_p2, CN_amp_num))
    return CN_array, new_chrlen

def is_in(a, mytuple):
    if float(a)>float(mytuple[0]) and float(a)<=float(mytuple[1]):
        return True
    return False

def print_chr_len(chrlen_array):
    print("chr len:")
    print(chrlen_array)
    return ""




#n= int(n)
#Alpha = float(Alpha)
#Beta = float(Beta)
#Delta = float(Delta)
# add a root (node 0) to the tree
# edge length (there are at most 2*n - 1))
#           root
#            | CN0
#          node 0
#        / CN1   \ CN2
#    node 1    node 2
ti = np.random.exponential(1,2*n-1)
#print len(ti)
Ui = np.random.uniform(0.0,1.0,n-1)
Vi = np.random.uniform(0.0,1.0,n-1)
Di = np.random.uniform(0.0,1.0,n-1)
Bi = np.random.beta(float(Alpha+1),float(Beta+1),n-1)

#Normalizing the branch lengths
summation = 0
for t in ti:
    summation += t

for T in range(0,len(ti)):
    ti[T]=float(ti[T])/float(summation)

#print ti


#Contructing the phylogeny
# by default chromosome size
# from hg19, Navin's 2012 paper
#chr_sz = [249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566, 155270560, 59373566]

# root is the node before node 0 in tree

root=MyNode("0: [0,1]")
root.tuple=[0,1]
ref_array, chr_name_array, chr_sz = init_ref(template_ref)
chr_sz1 = []
# copy so that the two arrays of allele length are independent
for i in chr_sz:
    chr_sz1.append(i)
root.chrlen=[chr_sz, chr_sz1]
#print chr_sz
root.id = -1

Tree = []
Tree.append(MyNode("0: [0,1]"))
Tree[0].tuple=[0,1]
Tree[0].id = 0
# assume most of the CN happens on the root branch
Tree[0].cn, Tree[0].chrlen = add_CN(root.chrlen, (cn_num * 4), del_rate, min_cn_size, exp_theta, amp_p)
#print "Node 0:"
#print Tree[0].chrlen
Tree[0].parent=root
Tree[0].edge_length = np.random.exponential(1,1)

# update the reference on the node
Tree[0].ref = gen_ref(ref_array, Tree[0].cn)

Tree.append(MyNode(str(1)+":[0,"+"{0:.2f}".format(Bi[0])+"]"+","+"{0:.4f}".format(ti[0])))
Tree.append(MyNode(str(2)+":["+"{0:.2f}".format(Bi[0])+",1]"+","+"{0:.4f}".format(ti[1])))
# add copy number
Tree[1].cn, Tree[1].chrlen = add_CN(Tree[0].chrlen, cn_num, del_rate, min_cn_size, exp_theta, amp_p)
#print "Node 1:"
#print Tree[1].chrlen
Tree[2].cn, Tree[2].chrlen = add_CN(Tree[0].chrlen, cn_num, del_rate, min_cn_size, exp_theta, amp_p)
#print "Node 2:"
#print Tree[2].chrlen

# update the reference
Tree[1].ref = gen_ref(Tree[0].ref, Tree[1].cn)
Tree[2].ref = gen_ref(Tree[0].ref, Tree[2].cn)

Tree[1].parent=Tree[0]
Tree[2].parent=Tree[0]
Tree[1].id = 1
Tree[2].id = 2
Tree[1].tuple=[0,Bi[0]]
Tree[2].tuple=[Bi[0],1]
Tree[1].edge_length = ti[0]
Tree[2].edge_length = ti[1]

node_number=2
j=1

while j<n-1:
    if Vi[j] < Delta :
        for tr in Tree:
            if tr.is_leaf and is_in(Di[j], tr.getTuple()):
                if (not tr.is_dead):
                    tr.name = tr.name+"*"
                tr.setDead()
                break
    else:
        for tree in Tree:
            if tree.is_leaf and is_in(Ui[j], tree.getTuple()) and (not tree.is_dead) :
                #print "the node from " + str(node_number + 1) + " to " + str(node_number+2) + "s' parent id: " + str(tree.getID())
                a,b = tree.getTuple()
                node_number+=2
                #Two new children are born here
                middle = float(Bi[j])*float((float(b)-float(a)))+float(a)
                Tree.append(MyNode(str(node_number-1)+":["+"{0:.4f}".format(a)+","+"{0:.4f}".format(middle)+"]"+","+"{0:.4f}".format(ti[node_number-1]), parent=tree))
                Tree.append(MyNode(str(node_number)+":["+"{0:.4f}".format(middle)+","+"{0:.4f}".format(b)+"]"+","+"{0:.4f}".format(ti[node_number]), parent=tree))

                #The new intervals are assigned here
                Tree[node_number-1].tuple=[a,middle]
                Tree[node_number].tuple=[middle,b]
                Tree[node_number-1].edge_length = ti[node_number-1]
                Tree[node_number].edge_length = ti[node_number]
                #add copy number
                this_chrlen = tree.chrlen[:]
                #print this_chrlen
                #print node_number, tree.getID()
                #print "node " + str(node_number - 1)
                Tree[node_number-1].cn, Tree[node_number-1].chrlen = add_CN(this_chrlen, cn_num, del_rate, min_cn_size, exp_theta, amp_p)
                this_chrlen = tree.chrlen[:]
                #print this_chrlen
                #print node_number, tree.getID()
                #print "node " + str(node_number)
                Tree[node_number].cn, Tree[node_number].chrlen = add_CN(this_chrlen, cn_num, del_rate, min_cn_size, exp_theta, amp_p)
                this_chrlen = tree.chrlen[:]
                #print this_chrlen
                #print node_number, tree.getID()

                # add reference
                Tree[node_number-1].ref = gen_ref(tree.ref, Tree[node_number-1].cn) 
                Tree[node_number].ref = gen_ref(tree.ref, Tree[node_number].cn) 

                # set id
                Tree[node_number-1].id = node_number - 1
                Tree[node_number].id = node_number

                break

    j+=1

#Changing names of the leaves
#leaf_name=0
#for nd in Tree:
#    if nd.is_leaf:
#        nd.name = leaf_name
#        leaf_name+=1

#for pre, fill, node in RenderTree(Tree[0]):
#    print("%s%s" % (pre, node.name))

f.write("Before the tree, chromosomomal length is " + str(root.chrlen) + "\n")
for i in range(len(Tree)):
    f.write("node %d: \n" % i)
    f.write("    parent = %d\n" % Tree[i].parent.getID())
    f.write("    name = " + str(Tree[i].name) + "\n")
for i in range(len(Tree)):
    cn = Tree[i].cn
    f.write("node %d from %d: total CN # = %d\n" % (i, Tree[i].parent.getID(), len(cn)))
    for j in range(len(cn)):
        f.write("    copy number %d: allele: %d, is del: %d, chromosome: %d, position: [%d, %d], amplification #: %d\n" % (j, cn[j].CN_Ale, cn[j].CN_Del, cn[j].CN_chromosome, cn[j].CN_p1, cn[j].CN_p2, cn[j].CN_amp_num))
    f.write("    " + str(Tree[i].chrlen) + "\n")
        #print_chr_len(Tree[i].chrlen)
#RenderTreeGraph(Tree[0]).to_picture(str(Output))


# generate reference for each leaf
for i in range(len(Tree)):
    fa_f_prefix = fa_prefix + str(i) + "_"
    write_ref(Tree[i].ref, chr_name_array, fa_f_prefix)

f.close()
