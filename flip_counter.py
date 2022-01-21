import ete3
from ete3 import Tree
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
import sys

def construct_Tree(dictionary, tree):
	flag = False
	for leaf in tree:
		if leaf.name in dictionary:
			flag = True
	if flag:
		for leaf in tree:
			if leaf.name in dictionary:
				children_names = dictionary[leaf.name]
				for child_name in children_names:
					leaf.add_child(name=child_name)
		tree = construct_Tree(dictionary,tree)
		return tree
	else:
		return tree


def parse(address):
	node_children_dict = {}
	node_assignment_dict = {}
	root_name = None
	with open(address,"r") as f:
		for line in f:
			line = line.strip().split('\t')
			# print line[0],line[2]
			#### line[0] is the name of the child 
			#### line[2] is the name of the parent
			#### line[3:] is the array of the bin states 
			node_assignment_dict[line[0]] = line[3:]
			if line[2]=="root":
				root_name = line[0]
				if line[0] not in node_children_dict:
					node_children_dict[line[0]]=[]
			else:
				if line[2] in node_children_dict:
					node_children_dict[line[2]].append(line[0])
				else:
					node_children_dict[line[2]] = [line[0]]
	###### Construct the tree according to the dictionary 
	t = Tree(name=root_name)
	t = construct_Tree(dictionary=node_children_dict,tree=t)
	for node in t.traverse('postorder'):
		node.add_features(states=node_assignment_dict[node.name])
	return t


def node_flip_counter(node):
	node_flip_array = [0 for i in range(len(node.states))]
	if len(node.get_ancestors())>0:
		ancestor = node.get_ancestors()[0]
		recent_changes = [0 if i==j else 1 for i,j in zip(node.states,ancestor.states)]
		past_changes = node_flip_counter(ancestor)
		node_flip_array = [x+y for x,y in zip(past_changes,recent_changes)]
	return node_flip_array

def tree_flip_counter(tree):
	flip_dict = {}
	for leaf in tree.get_leaves():
		flip_dict[leaf.name] = node_flip_counter(leaf)
	return flip_dict


def main():
	address = sys.argv[1]
	outfile_address = sys.argv[2]
	color_4_samples = ['royalblue','darkviolet','green','darkgoldenrod','red','sienna','magenta']
	
	fig, ax = plt.subplots(1,1)
	tree = parse(address=address)
	flips_dictionary = tree_flip_counter(tree=tree)
	flips_flat_list = []
	for key in flips_dictionary:
		flips_flat_list.extend(flips_dictionary[key])
	count_dict = Counter(flips_flat_list)
	y_pos = []
	count_arr = []
	for key in sorted(count_dict):
		if key <=9:
			y_pos.append(key)
			count_arr.append(count_dict[key])
	ax.set_ylabel("Percentage of the bins",fontsize=22,fontweight='bold')
	ax.set_xlabel("Number of copy number changes",fontsize=22,fontweight='bold')
	weights = [float(cnt)/float(sum(count_arr)) for cnt in count_arr]
	ax.bar(y_pos,weights,align='center',width=0.5,color=color_4_samples[0])
	ax.set_xticks([i for i in range(10)])
	ax.set_xticklabels([i+1 for i in range(10)],fontsize=22,fontweight='bold')
	# axis.set_ylim([0,1])
	for item in ax.get_yticklabels():
		item.set_fontsize(22)
		item.set_fontweight('bold')

	fig.set_size_inches(7,7)
	fig.tight_layout()
	plt.savefig(outfile_address,dpi=400)
	

if __name__ == "__main__":
	main()