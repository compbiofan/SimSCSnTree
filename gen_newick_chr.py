import gen_tree
import numpy
import sys
import re

if len(sys.argv) <= 1:
    print("""
    Read the npy tree file, convert to newick file without branch length, considering only one chromosome's contribution to the tree. If two genomes have the same CN profile for this chromosome, they are collapsed into one node.  
    Usage: python
    """ + sys.argv[0] + """ [tree.npy] [chr]""")
    sys.exit(0)

class mynode:
    def __init__(self, ID, parent=None):
        self.ID = ID
        # ID
        self.parent = parent
        self.children = []

# in par, given a node ID, find the very ancestor of ID
def seek_parent(ID, par):
    if ID not in par.keys():
        return ID
    else:
        while ID in par.keys():
            ID = par[ID]
        return ID

# now get rid of the non-leaves in the string
def look_for_leaves(str_, non_leaves):
    if str_.isdigit():
        return str_
    s = str_.find("(")
    e = str_.rfind(")")
    input_str = str_
    str_ = str_[s+1:e]
    nums = []
    if "," in str_:
        nums = str_.split(",")
    else:
        nums.append(str_)
    new_nums = []
    i = 0
    while i < len(nums):
        num = nums[i]
        while num.count('(') != num.count(')') and i < len(nums):
            # concatenate the strings when the ) is not the same as the left (
            i += 1
            num += "," + nums[i]
        new_nums.append(num)
        i += 1
    arr = []
    for i in new_nums:
        if i not in non_leaves:
            #print "Now look at " + i
            out = look_for_leaves(i, non_leaves)
            arr.append(out)
        #else:
            #print i + " is out"

    ret = "(" + ",".join(arr) + ")" + input_str[e+1:]
    #print "This is the return " + ret
    return ret
 
# step 2. convert to newick format
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
        
        
    
tree_file = sys.argv[1]
tree = numpy.load(tree_file, allow_pickle=True)
selected_chr = sys.argv[2]

# a new tree with a simplified structure
Tree = []

# all the keys in this dictionary are those that will be skipped. Their kids will be pointed to their parents (iteratively until no more such parent can be found in the key list of this dictionary)
par = {}

# non_leaves
non_leaves = {}
for i in range(len(tree)):
    node = tree[i]
    parent = node.parentID
    non_leaves[str(parent)] = 1
 

# step 1, read the npy file and put them in MyNode, if their chromosome of CN happens to be at the selected chromosome
for i in range(len(tree)):
    node = tree[i]
    parent = node.parentID
    new_parent = seek_parent(parent, par)
    CNs = node.cn
    tag = 0
    for j in range(len(CNs)):
        if str(CNs[j].CN_chromosome) == selected_chr:
            tag = 1
    if tag == 0:
        # connect all its children to its parent (which should have CNAs on selected chromosome)
        par[i] = new_parent
        Tree.append(mynode(-1, -1))
        continue
    ID = i
    Tree.append(mynode(ID, new_parent))
    # update the children of its parent
    Tree[new_parent].children.append(str(ID))

str_ = convert2newick(Tree, "(0)", 0)
str_ = "(" + look_for_leaves(str_[1:-1], non_leaves) + ")"
print(str_)

           

#print str_[1:-1]
#print str_


