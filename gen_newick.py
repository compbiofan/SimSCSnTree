import gen_tree
import numpy
import sys
import re

if len(sys.argv) <= 1:
    print("""
    Read the npy tree file, convert to newick file without branch length. 
    Usage: python
    """ + sys.argv[0] + """ [tree.npy]""")
    sys.exit(0)

class mynode:
    def __init__(self, ID, parent=None):
        self.ID = ID
        # ID
        self.parent = parent
        self.children = []

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

# a new tree with a simplified structure
Tree = []

# step 1, read the npy file and put them in MyNode
for i in range(len(tree)):
    node = tree[i]
    parent = node.parentID
    ID = i
    Tree.append(mynode(ID, parent))
    # update the children of its parent
    Tree[parent].children.append(str(ID))

str_ = convert2newick(Tree, "(0)", 0)
#print str_[1:-1]
print(str_)



