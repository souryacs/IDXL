#!/usr/bin/env python

##---------------------------------------------
''' 
this program is used to generate a spcies tree from a set of gene trees
the gene trees generally associate conflicts among the constituent genes (representing individual taxa)
in terms of the topological relationships
species tree follows simple approach to produce species tree with high performance index

Author: Sourya Bhattacharyya
Dept of CSE, IIT Kharagpur
V1.0 - 12.01.2015 - basic code
V1.1 - 20.02.2015 - written to handle deep coalescence and level sum, and accumulated rank statistics information
V1.2 - 27.02.2015 - handles accumulated coalescence rank sum, and accumulated branch length count
V1.3 - 30.03.2015 - code clean and necessary comments
V1.4 - 10.04.2015 - frequency + mode based couplet feature
V1.5 - 11.04.2015 - couplet feature extraction - timing optimization
V1.6 - 19.06.2015 - added LCA rank and accumulated coalescence rank based feature for species tree estimation
V1.7 - 20.06.2015 - added mode of LCA rank and corresponding species tree estimation technique
V1.8 - 17.07.2015 - added mode of LCA rank and mode of accumulated branch count as a new species tree estimation technique
V2.0 - 18.07.2015 - separated the branch count based computation and corresponding methods
''' 

## Copyright 2015 Sourya Bhattacharyya and Jayanta Mukherjee.
## All rights reserved.
##
## See "LICENSE.txt" for terms and conditions of usage.
##
##---------------------------------------------

import Header
from Header import *
import NJ_SpecTree
from NJ_SpecTree import *
import RankInfo
from RankInfo import *
import UtilFunc
from UtilFunc import *

##-----------------------------------------------------
# this function is useful to parse various options for input data processing
def parse_options():  
	parser = OptionParser()
		
	parser.add_option("-I", "--INPFILE", \
				type="string", \
				action="store", \
				dest="INP_FILENAME", \
				default="", \
				help="name of the input file containing gene trees")

	parser.add_option("-O", "--OUTFILE", \
				type="string", \
				action="store", \
				dest="OUT_FILENAME", \
				default="", \
				help="name of the output file to contain target species tree")  

	parser.add_option("-p", "--inpform", \
				type="int", \
				action="store", \
				dest="inp_file_format", \
				default=1, \
				help="1 - input file format is NEWICK (default) \
				2 - input file format is NEXUS")

	parser.add_option("-m", "--method", \
				type="int", \
				action="store", \
				dest="method_type", \
				default=1, \
				help="1 - NJSTXL (Average of internode count and average of excess gene count) \
				2 - MedNJSTXL (Average of internode count with min(average, median) of excess gene count) \
				3 - ModeNJSTXL (Average of internode count with min(average, median, mode) of excess gene count)")
			
	#parser.add_option("-n", "--njrule", \
				#type="int", \
				#action="store", \
				#dest="NJ_type", \
				#default=1, \
				#help="1 - classical NJ method (Default) \
				#2 - Normalized couplet statistic for agglomeration")     
	
	#parser.add_option("-d", "--distmat", \
				#type="int", \
				#action="store", \
				#dest="dist_mat_type", \
				#default=1, \
				#help="1 - Average of XL values \
				#2 - Minimum of median and average of XL values \
				#3 - Minimum of median, average, and mode of XL values")     
	
	parser.add_option("-r", "--ROOT", \
			type="string", \
			action="store", \
			dest="outgroup_taxon_name", \
			default="", \
			help="name of the taxon by which the tree can be rooted (outgroup based rooting)")
			
	opts, args = parser.parse_args()
	return opts, args
  
##-----------------------------------------------------
# main function
def main():  
	opts, args = parse_options()

	ROOTED_TREE = True
	PRESERVE_UNDERSCORE = True
	if (opts.inp_file_format == 1):
		INPUT_FILE_FORMAT = 'newick'
	else:
		INPUT_FILE_FORMAT = 'nexus'
	INPUT_FILENAME = opts.INP_FILENAME
	OUTPUT_FILENAME = opts.OUT_FILENAME
	METHOD_USED = opts.method_type
	
	OUTGROUP_TAXON_NAME = opts.outgroup_taxon_name

	NJ_RULE_USED = TRADITIONAL_NJ	#opts.NJ_type
	
	## add - sourya
	## we select the NJ option based on the couplet feature employed
	## if it is average accumulated rank or average branch count then we use traditional NJ method
	## otherwise we use the weighted Normalized agglomerative clustering
	#if (METHOD_USED == NJ_ST_XL):
		#NJ_RULE_USED = TRADITIONAL_NJ
	#else:
		#NJ_RULE_USED = AGGLO_CLUST
	## end add - sourya
		
	if (INPUT_FILENAME == ""):
		print '******** THERE IS NO INPUT FILE (CONTAINING GENE TREE LIST) SPECIFIED - RETURN **********'
		return
	else:
		print 'input filename: ', INPUT_FILENAME

	# according to the location of input filename
	# adjust the locations of the output files as well
	k = INPUT_FILENAME.rfind("/")
	if (k == -1):
		dir_of_inp_file = './'
	else:
		dir_of_inp_file = INPUT_FILENAME[:(k+1)]
	inp_filename = INPUT_FILENAME[(k+1):]
	print 'dir_of_inp_file: ', dir_of_inp_file  

	if (OUTPUT_FILENAME == ""):
		# derive the output directory which will contain different output text results
		if (METHOD_USED == NJSTXL):
			dir_of_curr_exec = dir_of_inp_file + 'NJSTXL'
		elif (METHOD_USED == MedNJSTXL):
			dir_of_curr_exec = dir_of_inp_file + 'MedNJSTXL'
		elif (METHOD_USED == ModeNJSTXL):
			dir_of_curr_exec = dir_of_inp_file + 'ModeNJSTXL'
			
		# append the current output directory in the text file
		Output_Text_File = dir_of_curr_exec + '/' + 'Complete_Desription.txt'
		# create the directory
		if (os.path.isdir(dir_of_curr_exec) == False):
			mkdr_cmd = 'mkdir ' + dir_of_curr_exec
			os.system(mkdr_cmd)               
		print 'Output_Text_File: ', Output_Text_File      
	else:
		# when we have specified the output file then corresponding directory becomes the dir_of_curr_exec
		k1 = OUTPUT_FILENAME.rfind("/")
		if (k1 == -1):
			dir_of_curr_exec = './'
		else:
			dir_of_curr_exec = OUTPUT_FILENAME[:(k1+1)]
		Output_Text_File = OUTPUT_FILENAME + '_complete_text_description'   
		print 'Output_Text_File: ', Output_Text_File
		# create the directory
		if (os.path.isdir(dir_of_curr_exec) == False):
			mkdr_cmd = 'mkdir ' + dir_of_curr_exec
			os.system(mkdr_cmd)               

	# note the program beginning time 
	start_timestamp = time.time()

	#-------------------------------------------------------------
	""" 
	read the source trees collection and store it in a tree collection structure
	individual elements of this collection is thus a source tree 
	"""
	Gene_TreeList = Read_Input_Treelist(ROOTED_TREE, PRESERVE_UNDERSCORE, INPUT_FILE_FORMAT, INPUT_FILENAME)
	#-------------------------------------------------------------
	# from the input source trees, note the number of taxa (total)
	# and also define the class instances corresponding to single taxa
	for tr_idx in range(len(Gene_TreeList)):
		taxa_labels_curr_tree = Gene_TreeList[tr_idx].infer_taxa().labels()
		for i in range(len(taxa_labels_curr_tree)):
			if taxa_labels_curr_tree[i] not in COMPLETE_INPUT_TAXA_LIST:
				COMPLETE_INPUT_TAXA_LIST.append(taxa_labels_curr_tree[i])
		
	# now process individual trees to find the couplet relations of those trees
	for tr_idx in range(len(Gene_TreeList)):
		DeriveCoupletRelations(Gene_TreeList[tr_idx], METHOD_USED)

	# note the time
	end_timestamp = time.time() 
		
	# note the time taken for reading the dataset
	time_read_dataset = (end_timestamp - start_timestamp)

	fp = open(Output_Text_File, 'w')    
	fp.write('\n  total no of taxa: ' + str(len(COMPLETE_INPUT_TAXA_LIST)))
	if (DEBUG_LEVEL >= 2):
		fp.write('\n COMPLETE_INPUT_TAXA_LIST ' + str(COMPLETE_INPUT_TAXA_LIST))
	fp.write('\n \n ===============>>>>>>>>>>>>>>> TIME TAKEN FOR READING THE INPUT GENE TREES ' + str(time_read_dataset))
	fp.close()
	
	"""
	this to to print individual couplet information
	"""
	if (DEBUG_LEVEL >= 2):
		for coup in TaxaPair_Reln_Dict:
			TaxaPair_Reln_Dict[coup]._PrintTaxaPairRelnInfo(coup, Output_Text_File, METHOD_USED)
			
	# generate a star network from the input taxa labels
	# form a newick formatted string containing the tree
	star_net_str = ""
	for i in range(len(COMPLETE_INPUT_TAXA_LIST)):
		star_net_str = star_net_str + str(COMPLETE_INPUT_TAXA_LIST[i])
		if (i < (len(COMPLETE_INPUT_TAXA_LIST) - 1)):
			star_net_str = star_net_str + ","
	star_net_str = "(" + star_net_str + ")"

	# this tree denotes the initial star configuration (rooted at the central hub)
	Output_Tree = dendropy.Tree.get_from_string(star_net_str, schema="newick", \
								preserve_underscores=PRESERVE_UNDERSCORE, \
								default_as_rooted=True)          

	if (DEBUG_LEVEL > 2):
		fp = open(Output_Text_File, 'a')    
		fp.write('\n star_net_str: ' + str(star_net_str))
		fp.write('\n from tree ---: ' + Output_Tree.as_newick_string())
		fp.close()

	# now perform the agglomerative clustering technique based on the extra lineages
	Form_Species_Tree_NJ_Cluster(Output_Tree, METHOD_USED, NJ_RULE_USED, Output_Text_File)

	fp = open(Output_Text_File, 'a')
	fp.write('\n --- output species tree (in newick format): ' + Output_Tree.as_newick_string())    
	fp.close()

	out_treefilename = dir_of_curr_exec + '/' + 'outtree_Newick.tre'

	# we write the unweighted supertree
	outfile = open(out_treefilename, 'w')
	outfile.write(Output_Tree.as_newick_string())	
	outfile.close()

	# add  -sourya
	"""
	if user specifies one outgroup taxon name for rooting the tree, 
	we should reroot the tree and save it
	"""
	if (OUTGROUP_TAXON_NAME != ""):
		# first we root the tree using the specified outgroup node
		outgroup_node = Output_Tree.find_node_with_taxon_label(OUTGROUP_TAXON_NAME)
		if (outgroup_node is not None):
			Output_Tree.to_outgroup_position(outgroup_node, update_splits=False)
		
		out_treefilename = dir_of_curr_exec + '/' + 'outtree_Newick.tre.rooted'

		# we write the unweighted supertree
		outfile = open(out_treefilename, 'w')
		outfile.write(Output_Tree.as_newick_string())	
		outfile.close()
	# end add  -sourya

	# note the time
	end_timestamp = time.time()      

	# we write the time associated with the execution of this method
	fp = open(Output_Text_File, 'a')
	fp.write('\n \n\n ===============>>>>>>>>>>>>>>> TIME COMPLEXITY OF THE METHOD (in seconds) ')
	fp.write('\n \n Total time taken (in seconds) : ' + str(end_timestamp - start_timestamp))
	fp.close()

	# clear the storage  variables
	TaxaPair_Reln_Dict.clear()
	if (len(COMPLETE_INPUT_TAXA_LIST) > 0):
		COMPLETE_INPUT_TAXA_LIST[:] = []  
  
#-----------------------------------------------------
if __name__ == "__main__":
	main() 




