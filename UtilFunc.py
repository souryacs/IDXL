#!/usr/bin/env python

import Header
from Header import * 

#-----------------------------------------------------
# this function reads the input tree list file
# parameters: ROOTED_TREE - whether the treelist to be read as rooted format
# PRESERVE_UNDERSCORE: whether underscores of the taxa name will be preserved or not
# INPUT_FILE_FORMAT: data is read from the file according to NEWICK or NEXUS format
# INPUT_FILENAME: file containing the input treelist

def Read_Input_Treelist(ROOTED_TREE, PRESERVE_UNDERSCORE, INPUT_FILE_FORMAT, INPUT_FILENAME):
	Inp_TreeList = dendropy.TreeList.get_from_path(INPUT_FILENAME, schema=INPUT_FILE_FORMAT, \
		preserve_underscores=PRESERVE_UNDERSCORE, default_as_rooted=ROOTED_TREE)
	return Inp_TreeList

#--------------------------------------------------
# this function returns the label of an internal or a leaf node 
# in terms of newick representation
def Node_Label(inp_node):
	return str(inp_node.as_newick_string(suppress_edge_lengths=True))

#---------------------------------------------
"""
the function is used for checking equality of two floating point numbers
"""
def FlEq( a, b, eps=0.000001):
	#return (abs(math.log( a ) - math.log(b)) <= eps)
	return (abs(a - b) <= eps)

