#!/usr/bin/env python

import Header
from Header import * 

##-----------------------------------------------------
# this function reads the input tree list file
# parameters: ROOTED_TREE - whether the treelist to be read as rooted format
# PRESERVE_UNDERSCORE: whether underscores of the taxa name will be preserved or not
# INPUT_FILE_FORMAT: data is read from the file according to NEWICK or NEXUS format
# INPUT_FILENAME: file containing the input treelist

def Read_Input_Treelist(ROOTED_TREE, PRESERVE_UNDERSCORE, INPUT_FILE_FORMAT, INPUT_FILENAME):
	Inp_TreeList = dendropy.TreeList.get_from_path(INPUT_FILENAME, schema=INPUT_FILE_FORMAT, \
		preserve_underscores=PRESERVE_UNDERSCORE, default_as_rooted=ROOTED_TREE)
	return Inp_TreeList

###-----------------------------------------------------
## this function finds the MRCA of this two input taxa labels
## this is a custom function
## without using standard dendropy routine
#def Find_MRCA(Inp_Tree, spec_list):
	#node1 = Inp_Tree.find_node_with_taxon_label(spec_list[0])
	#pn = node1.parent_node
	#while (pn is not None):
		#leaf_labels = []
		#for n in pn.leaf_nodes():
			#leaf_labels.append(n.taxon.label)
		#if set(spec_list).issubset(set(leaf_labels)):
			#return pn
		#pn = pn.parent_node
			
	#return None

#--------------------------------------------------
# this function returns the label of an internal or a leaf node 
# in terms of newick representation
def Node_Label(inp_node):
	return str(inp_node.as_newick_string(suppress_edge_lengths=True))

#----------------------------------------
#def Complementary_Reln(inp_reln):
  #if (inp_reln == RELATION_R3) or (inp_reln == RELATION_R4):
    #return inp_reln
  #elif (inp_reln == RELATION_R1):
    #return RELATION_R2
  #else:
    #return RELATION_R1
