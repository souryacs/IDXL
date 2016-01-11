import Header
from Header import *
import UtilFunc
from UtilFunc import *

#--------------------------------------------------------
# this function is a shortcut to obtain the normalized expression 
# used in the agglomerative clustering proposed in this code
# as various methods are experimented, corresponding various forms of 
# agglomerative clustering is tried
#--------------------------------------------------------
def ObtainNormalizedVal(num, denom1, denom2):
	if ((denom1 + denom2) > 0):
		return (num * 1.0) / (denom1 + denom2)
	else:
		return 0

#---------------------------------------------
""" 
function to print the matrix content
N is the matrix dimension
"""
def PrintMatrixContent(N, TaxaList, inp_data, inp_str, textfile):
	fp = open(textfile, 'a')
	fp.write('\n printing contents of ' + str(inp_str) + ' ---- ')
	for i in range(N):
		fp.write('\n ' + str(i) + '--' + str(TaxaList[i]) + '--->>')
		if (i > 0):
			for j in range(i):
				fp.write(' ' + str(inp_data[i][j]))
	fp.close()

#----------------------------------------------------
"""
this function fills the distance matrix using accumulated internode count
"""
def Fill_DistMat_BranchInfo(DistMat, METHOD_USED):
	for l in TaxaPair_Reln_Dict:
		spec1 = l[0]
		spec2 = l[1]
		spec1_idx = COMPLETE_INPUT_TAXA_LIST.index(spec1)
		spec2_idx = COMPLETE_INPUT_TAXA_LIST.index(spec2)
		DistMat[spec2_idx][spec1_idx] = TaxaPair_Reln_Dict[l]._GetAvgSumLevel()
		DistMat[spec1_idx][spec2_idx] = DistMat[spec2_idx][spec1_idx]
	
	return

#----------------------------------------------------
"""
if excess gene count information is used, this 
function fills the distance matrix using normalized average excess gene count
"""
def Fill_DistMat_ExcessGeneCount(DistMat, METHOD_USED, DIST_MAT_TYPE):
	for l in TaxaPair_Reln_Dict:
		spec1 = l[0]
		spec2 = l[1]
		spec1_idx = COMPLETE_INPUT_TAXA_LIST.index(spec1)
		spec2_idx = COMPLETE_INPUT_TAXA_LIST.index(spec2)
		
		# code used earlier  - commented - sourya
		if (METHOD_USED == NJSTXL):
			DistMat[spec2_idx][spec1_idx] = TaxaPair_Reln_Dict[l]._GetAvgXLVal()
		elif (METHOD_USED == MedNJSTXL):
			# minimum of average and median of XL values
			DistMat[spec2_idx][spec1_idx] = min(TaxaPair_Reln_Dict[l]._GetAvgXLVal(), \
				TaxaPair_Reln_Dict[l]._MedianXLVal())
		elif (METHOD_USED == ProdNJSTXL):
			# minimum of average and median of XL values
			DistMat[spec2_idx][spec1_idx] = min(TaxaPair_Reln_Dict[l]._GetAvgXLVal(), \
				TaxaPair_Reln_Dict[l]._MedianXLVal())
			#DistMat[spec2_idx][spec1_idx] = TaxaPair_Reln_Dict[l]._MedianXLVal()
		# end comment - sourya

		## add - sourya
		#if (DIST_MAT_TYPE == 1):
			#DistMat[spec2_idx][spec1_idx] = TaxaPair_Reln_Dict[l]._GetAvgXLVal()
		#elif (DIST_MAT_TYPE == 2):
			#DistMat[spec2_idx][spec1_idx] = TaxaPair_Reln_Dict[l]._MedianXLVal()
		#elif (DIST_MAT_TYPE == 3):
			#DistMat[spec2_idx][spec1_idx] = TaxaPair_Reln_Dict[l]._GetMultiModeXLVal()
		#elif (DIST_MAT_TYPE == 4):
			## minimum of average and median of XL values
			#DistMat[spec2_idx][spec1_idx] = min(TaxaPair_Reln_Dict[l]._GetAvgXLVal(), TaxaPair_Reln_Dict[l]._MedianXLVal())
		#elif (DIST_MAT_TYPE == 5):
			## minimum of average, median, and mode of XL values
			#DistMat[spec2_idx][spec1_idx] = min(TaxaPair_Reln_Dict[l]._GetAvgXLVal(), \
				#TaxaPair_Reln_Dict[l]._MedianXLVal(), TaxaPair_Reln_Dict[l]._GetMultiModeXLVal())
		#elif (DIST_MAT_TYPE == 6):
			## minimum of median, and mode of XL values
			#DistMat[spec2_idx][spec1_idx] = min(TaxaPair_Reln_Dict[l]._MedianXLVal(), TaxaPair_Reln_Dict[l]._GetMultiModeXLVal())
		## end add - sourya
		
		# symmetric property of the distance matrix
		DistMat[spec1_idx][spec2_idx] = DistMat[spec2_idx][spec1_idx]

	return

#---------------------------------------------
"""
computing the row wise sum for individual taxca clusters
"""
def ComputeSumRowsDistMat(sum_list, nclust, DistMat, Output_Text_File, inpstr):
	for i in range(nclust):
		t1 = 0
		for j in range(nclust):
			t1 = t1 + DistMat[i][j]
		sum_list.append(t1)
		
	if (DEBUG_LEVEL > 2):
		fp = open(Output_Text_File, 'a')
		fp.write('\n content of ' + str(inpstr) + ' : ' + str(sum_list))
		fp.close()
	
	return

#----------------------------------------------------------
"""
fill the normalized matrix entries for agglomerative clustering based method
"""
def FillAggloClustNormalizeMatrix(NormMat, DistMat, sum_list, nclust):
	for i in range(nclust - 1):
		for j in range(i+1, nclust):
			# comment - sourya
			#NormMat[i][j] = ObtainNormalizedVal(DistMat[i][j], sum_list[i], sum_list[j])
			# add - sourya
			NormMat[i][j] = ObtainNormalizedVal(DistMat[i][j], (sum_list[i] - DistMat[i][j]), \
				(sum_list[j] - DistMat[i][j]))
			# end add - sourya
			NormMat[j][i] = NormMat[i][j]

	return

#----------------------------------------------------------
"""
fill the normalized matrix entries for NJ based method
"""
def FillNJNormalizeMatrix(NormMat, DistMat, sum_list, nclust):
	for i in range(nclust - 1):
		for j in range(i+1, nclust):
			ri = sum_list[i] / (nclust - 2)
			rj = sum_list[j] / (nclust - 2)
			NormMat[i][j] = (DistMat[i][j] - ri - rj)
			NormMat[j][i] = NormMat[i][j]

	return

#---------------------------------------------
"""
finds the minimum of the distance matrix
when the product of both accumulated coalescence rank and excess gene count based 
measures are used
"""
def Find_Unique_Min_XL(DistMat_CoalRank, Norm_DistMat_CoalRank, DistMat_XL, \
	Norm_DistMat_XL, nclust, clust_species_list, NJ_RULE_USED):
	
	target_val = (Norm_DistMat_CoalRank[0][1] * Norm_DistMat_XL[0][1])
	min_idx_i = 0
	min_idx_j = 1
	for i in range(nclust - 1):
		for j in range(i+1, nclust):
			if (i == j):
				continue
			if (NJ_RULE_USED == AGGLO_CLUST):
				if ((Norm_DistMat_CoalRank[i][j] * Norm_DistMat_XL[i][j]) < target_val):
					target_val = (Norm_DistMat_CoalRank[i][j] * Norm_DistMat_XL[i][j])
					min_idx_i = i
					min_idx_j = j
			else:
				if ((Norm_DistMat_CoalRank[i][j] * Norm_DistMat_XL[i][j]) > target_val):
					target_val = (Norm_DistMat_CoalRank[i][j] * Norm_DistMat_XL[i][j])
					min_idx_i = i
					min_idx_j = j
				
	return min_idx_i, min_idx_j

##---------------------------------------------
# sourya - keep this function but as a commented entry
##----------------------------------------------
#"""
#finds the minimum of the distance matrix
#when only branch count based measure is used
#"""
#def Find_Unique_Min(Norm_DistMat_Branch, nclust):
	#target_val = Norm_DistMat_Branch[0][1]
	#min_idx_i = 0
	#min_idx_j = 1
	#for i in range(nclust - 1):
		#for j in range(i+1, nclust):
			#if (i == j):
				#continue
			#if (Norm_DistMat_Branch[i][j] < target_val):
				#target_val = Norm_DistMat_Branch[i][j]
				#min_idx_i = i
				#min_idx_j = j
	
	#return min_idx_i, min_idx_j

#-------------------------------------------
"""
this function is a new version, to select the couplet with minimum aggregated rank
ranks are computed for both XL and branch count
minimum sum of ranks are considered
"""
def Find_Unique_Min_RankBased(DistMat_Branch, Norm_DistMat_Branch, DistMat_XL, \
	Norm_DistMat_XL, nclust, clust_species_list, outfile, RankMergeMethod):
	
	if (DEBUG_LEVEL >= 2):
		fp = open(outfile, 'a')
	
	Norm_DistMat_List = []
	DistMat_XL_List = []
	for i in range(nclust - 1):
		for j in range(i+1, nclust):
			if (i == j):
				continue
			Norm_DistMat_List.append(Norm_DistMat_Branch[i][j])
			DistMat_XL_List.append((DistMat_XL[i][j]))
	
	Norm_DistMat_List.sort()
	DistMat_XL_List.sort()
	
	if (RankMergeMethod == SIMPLE_SUM_RANK):
		min_idx_i = 0
		min_idx_j = 1
		min_rank_Norm_DistMat_List = Norm_DistMat_List.index(Norm_DistMat_Branch[0][1])
		min_rank_DistMat_XL_List = DistMat_XL_List.index(DistMat_XL[0][1])
		min_rank = min_rank_Norm_DistMat_List + min_rank_DistMat_XL_List
	elif (RankMergeMethod == MEAN_RECIPROCAL_RANK):
		max_idx_i = 0
		max_idx_j = 1
		# we allow rank values from 1, instead of 0
		max_rank_Norm_DistMat_List = Norm_DistMat_List.index(Norm_DistMat_Branch[0][1]) + 1
		max_rank_DistMat_XL_List = DistMat_XL_List.index(DistMat_XL[0][1]) + 1
		# compute the mean reciprocal rank value
		max_rank = 0.5 * ((1.0 / max_rank_Norm_DistMat_List) + (1.0 / max_rank_DistMat_XL_List))
		
	for i in range(nclust - 1):
		for j in range(i+1, nclust):
			if (i == j):
				continue
			if (RankMergeMethod == SIMPLE_SUM_RANK):
				Norm_DistMat_Rank = Norm_DistMat_List.index(Norm_DistMat_Branch[i][j])
				DistMat_XL_rank = DistMat_XL_List.index(DistMat_XL[i][j])
				total_rank = Norm_DistMat_Rank + DistMat_XL_rank
			elif (RankMergeMethod == MEAN_RECIPROCAL_RANK):
				# we allow rank values from 1, instead of 0
				Norm_DistMat_Rank = Norm_DistMat_List.index(Norm_DistMat_Branch[i][j]) + 1
				DistMat_XL_rank = DistMat_XL_List.index(DistMat_XL[i][j]) + 1
				total_rank = 0.5 * ((1.0 / Norm_DistMat_Rank) + (1.0 / DistMat_XL_rank))
				
			if (RankMergeMethod == SIMPLE_SUM_RANK):
				if (total_rank < min_rank):
					min_idx_i = i
					min_idx_j = j
					min_rank = total_rank
					min_rank_Norm_DistMat_List = Norm_DistMat_Rank
					min_rank_DistMat_XL_List = DistMat_XL_rank
					if (DEBUG_LEVEL >= 2):
						fp.write('\n Within iteration --- min_idx_i ' + str(min_idx_i) + ' min_idx_j : ' + str(min_idx_j) + \
							' min_rank_Norm_DistMat_List: ' + str(min_rank_Norm_DistMat_List) \
								+ '  min_rank_DistMat_XL_List: ' + str(min_rank_DistMat_XL_List))
				elif (total_rank == min_rank):
					if (Norm_DistMat_Rank < min_rank_Norm_DistMat_List) and (min_rank_DistMat_XL_List > 0):	# add - sourya
						min_idx_i = i
						min_idx_j = j
						min_rank = total_rank
						min_rank_Norm_DistMat_List = Norm_DistMat_Rank
						min_rank_DistMat_XL_List = DistMat_XL_rank
						if (DEBUG_LEVEL >= 2):
							fp.write('\n Within iteration --- min_idx_i ' + str(min_idx_i) + ' min_idx_j : ' + str(min_idx_j) + \
								' min_rank_Norm_DistMat_List: ' + str(min_rank_Norm_DistMat_List) \
									+ '  min_rank_DistMat_XL_List: ' + str(min_rank_DistMat_XL_List))

					elif (min_rank_Norm_DistMat_List > 0) and (min_rank_DistMat_XL_List > 0) \
						and ((Norm_DistMat_Rank == 0) or (DistMat_XL_rank == 0)):	# add - sourya
						min_idx_i = i
						min_idx_j = j
						min_rank = total_rank
						min_rank_Norm_DistMat_List = Norm_DistMat_Rank
						min_rank_DistMat_XL_List = DistMat_XL_rank
						if (DEBUG_LEVEL >= 2):
							fp.write('\n Within iteration --- min_idx_i ' + str(min_idx_i) + ' min_idx_j : ' + str(min_idx_j) + \
								' min_rank_Norm_DistMat_List: ' + str(min_rank_Norm_DistMat_List) + '  min_rank_DistMat_XL_List: ' + \
									str(min_rank_DistMat_XL_List))
							
			elif (RankMergeMethod == MEAN_RECIPROCAL_RANK):
				if (total_rank > max_rank):
					max_idx_i = i
					max_idx_j = j
					max_rank = total_rank
					max_rank_Norm_DistMat_List = Norm_DistMat_Rank
					max_rank_DistMat_XL_List = DistMat_XL_rank
					if (DEBUG_LEVEL >= 2):
						fp.write('\n Within iteration --- max_idx_i ' + str(max_idx_i) \
							+ ' max_idx_j : ' + str(max_idx_j) + \
							' max_rank_Norm_DistMat_List: ' + str(max_rank_Norm_DistMat_List) \
								+ '  max_rank_DistMat_XL_List: ' + str(max_rank_DistMat_XL_List))

	if (DEBUG_LEVEL >= 2):
		fp.close()
	
	if (RankMergeMethod == SIMPLE_SUM_RANK):
		return min_idx_i, min_idx_j
	elif (RankMergeMethod == MEAN_RECIPROCAL_RANK):
		return max_idx_i, max_idx_j

#-------------------------------------------
"""
checks whether a taxa cluster specified by the input index is a leaf
"""
def IsLeafCluster(clust_species_list, idx):
	if (len(clust_species_list[idx]) == 1):
		return True
	return False

#-------------------------------------------
"""
this function has following parameters:
1) first_cluster_mrca_node: root of 1st subtree 
2) second_cluster_mrca_node: root of 2nd subtree 
3) all_taxa_mrca_node: root of all these trees
4) Curr_tree: Tree containing all these subtrees

It creates one new internal node as a child of all_taxa_mrca_node
and places above mentioned subtrees as its children
"""
def MergeSubtrees(Curr_tree, first_cluster_mrca_node, second_cluster_mrca_node, \
	all_taxa_mrca_node, taxa_list, Output_Text_File):
	if (DEBUG_LEVEL >= 2):
		fp = open(Output_Text_File, 'a')
		fp.write('\n label of first_cluster_mrca_node: ' + str(Node_Label(first_cluster_mrca_node)))      
		fp.write('\n label of second_cluster_mrca_node: ' + str(Node_Label(second_cluster_mrca_node)))
		fp.write('\n label of all_taxa_mrca_node: ' + str(Node_Label(all_taxa_mrca_node)))
		fp.close()

	# create new internal node 
	newnode = dendropy.Node()
	
	# its parent node will be the previous MRCA node of all the taxa in two clusters
	all_taxa_mrca_node.add_child(newnode)
	newnode.parent_node = all_taxa_mrca_node
	all_taxa_mrca_node.remove_child(first_cluster_mrca_node)
	first_cluster_mrca_node.parent_node = None
	all_taxa_mrca_node.remove_child(second_cluster_mrca_node)
	second_cluster_mrca_node.parent_node = None
	
	# add these individual clusters' MRCA node as its children
	newnode.add_child(first_cluster_mrca_node)
	first_cluster_mrca_node.parent_node = newnode
	newnode.add_child(second_cluster_mrca_node)
	second_cluster_mrca_node.parent_node = newnode
	
	# update splits of the resulting tree
	Curr_tree.update_splits(delete_outdegree_one=False)

	if (DEBUG_LEVEL >= 2):
		fp = open(Output_Text_File, 'a')
		fp.write('\n label of newnode: ' + str(Node_Label(newnode)))
		fp.write('\n label of all taxa mrca node (recomputed): ' + \
			str(Node_Label(Curr_tree.mrca(taxon_labels=taxa_list))))  
		fp.close()

	return Curr_tree

#-------------------------------------------
"""
this function merges a pair of clusters whose indices are pointed by the min_idx_i and min_idx_j entries
this is part of the proposed agglomerative clustering
taxa_list is the union of these two clusters (species contents)
"""
def Merge_Cluster_Pair(Curr_tree, clust_species_list, min_idx_i, min_idx_j, taxa_list, Output_Text_File):
	isleaf_clust1 = IsLeafCluster(clust_species_list, min_idx_i)
	isleaf_clust2 = IsLeafCluster(clust_species_list, min_idx_j)
	
	if (isleaf_clust1):
		first_cluster_mrca_node = Curr_tree.find_node_with_taxon_label(clust_species_list[min_idx_i][0])
	else:
		first_cluster_mrca_node = Curr_tree.mrca(taxon_labels=clust_species_list[min_idx_i])
	
	if (isleaf_clust2):
		second_cluster_mrca_node = Curr_tree.find_node_with_taxon_label(clust_species_list[min_idx_j][0])
	else:
		second_cluster_mrca_node = Curr_tree.mrca(taxon_labels=clust_species_list[min_idx_j])
	
	all_taxa_mrca_node = Curr_tree.mrca(taxon_labels=taxa_list)
	
	Curr_tree = MergeSubtrees(Curr_tree, first_cluster_mrca_node, second_cluster_mrca_node, \
		all_taxa_mrca_node, taxa_list, Output_Text_File)
	
	return Curr_tree


#---------------------------------------------
"""
this function refines the initial species tree (in terms of a star network) to 
find the true species tree
it does using agglomerative clustering (NJ principle)
the distance metric employed for NJ algorithm can vary depending on experimentation 
"""
def Form_Species_Tree_NJ_Cluster(Curr_tree, METHOD_USED, NJ_RULE_USED, \
	Output_Text_File, RankMergeMethod, DIST_MAT_TYPE):	#, DIST_MAT_UPDATE):

	# initially we have N of clusters for N taxa, where individual clusters are isolated
	# agglomerating technique introduces a bipartition (speciation) which contains two taxa as its children
	no_of_taxa_clust = len(COMPLETE_INPUT_TAXA_LIST)

	# initialize the taxa clusters
	# copying the taxa list is done since initial clusters contain single species  
	# comment - sourya - we do not just copy ordinarily
	#clust_species_list = COMPLETE_INPUT_TAXA_LIST[:]
	# add - sourya - we enclose individual elements within a list and then copy these single element lists
	clust_species_list = []
	for i in range(len(COMPLETE_INPUT_TAXA_LIST)):
		subl = []
		subl.append(COMPLETE_INPUT_TAXA_LIST[i])
		clust_species_list.append(subl)

	# allocate a 2D square matrix of no_of_taxa_clust dimension
	# for a pair of taxa clusters Cx and Cy, it contains the employed main distance metric for the cluster pairs
	Mean_DistMat_ClustPair_NJ = numpy.zeros((no_of_taxa_clust, no_of_taxa_clust), dtype=numpy.float)

	# allocate one new square matrix which will contain the NJ based modified distance matrix (used for minimum finding routine)
	Norm_Mean_DistMat_ClustPair_NJ = numpy.zeros((no_of_taxa_clust, no_of_taxa_clust), dtype=numpy.float)

	# we allocate matrices containing excess gene measure computed for individual couplets
	XLVal_DistMat_ClustPair_NJ = numpy.zeros((no_of_taxa_clust, no_of_taxa_clust), dtype=numpy.float)
	
	# allocate one new square matrix which will contain the normalized excess gene measure
	# used in NJ iterations
	Norm_XLVal_DistMat_ClustPair_NJ = numpy.zeros((no_of_taxa_clust, no_of_taxa_clust), dtype=numpy.float)

	"""
	fill input distance matrices using accumulated branch count or internode count information
	"""
	Fill_DistMat_BranchInfo(Mean_DistMat_ClustPair_NJ, METHOD_USED)

	"""
	this function fills input distance matrices with excess gene information
	"""
	Fill_DistMat_ExcessGeneCount(XLVal_DistMat_ClustPair_NJ, METHOD_USED, DIST_MAT_TYPE)

	#--------------------------------------------------------
	# loop to execute the agglomerative clustering
	while(no_of_taxa_clust > 2): 
		if (DEBUG_LEVEL >= 2):
			fp = open(Output_Text_File, 'a')
			fp.write('\n iteration start --- number of clusters: ' + str(no_of_taxa_clust))
			fp.write('\n clust_species_list : ' + str(clust_species_list))
			fp.close()
			PrintMatrixContent(no_of_taxa_clust, clust_species_list, Mean_DistMat_ClustPair_NJ, \
				'Mean_DistMat_ClustPair_NJ', Output_Text_File)
			PrintMatrixContent(no_of_taxa_clust, clust_species_list, XLVal_DistMat_ClustPair_NJ, \
				'XLVal_DistMat_ClustPair_NJ', Output_Text_File)

		# for individual cluster Cx, it contains XL(Cx, :) - sum of extra lineages considering the cluster pair 
		# (Cx, Cy) for all other clusters Cy
		sum_DistMat_Clust = []
		ComputeSumRowsDistMat(sum_DistMat_Clust, no_of_taxa_clust, \
			Mean_DistMat_ClustPair_NJ, Output_Text_File, 'sum_DistMat_Clust')
		sum_XLVal_Clust = []
		ComputeSumRowsDistMat(sum_XLVal_Clust, no_of_taxa_clust, \
			XLVal_DistMat_ClustPair_NJ, Output_Text_File, 'sum_XLVal_Clust')
		
		"""
		fill the normalized matrix entries
		depending on the NJ type method employed
		"""
		if (NJ_RULE_USED == AGGLO_CLUST):
			FillAggloClustNormalizeMatrix(Norm_Mean_DistMat_ClustPair_NJ, \
				Mean_DistMat_ClustPair_NJ, sum_DistMat_Clust, no_of_taxa_clust)
			FillAggloClustNormalizeMatrix(Norm_XLVal_DistMat_ClustPair_NJ, \
				XLVal_DistMat_ClustPair_NJ, sum_XLVal_Clust, no_of_taxa_clust)
		else:
			FillNJNormalizeMatrix(Norm_Mean_DistMat_ClustPair_NJ, \
				Mean_DistMat_ClustPair_NJ, sum_DistMat_Clust, no_of_taxa_clust)
			FillNJNormalizeMatrix(Norm_XLVal_DistMat_ClustPair_NJ, \
				XLVal_DistMat_ClustPair_NJ, sum_XLVal_Clust, no_of_taxa_clust)

		if (DEBUG_LEVEL >= 2):
			PrintMatrixContent(no_of_taxa_clust, clust_species_list, Norm_Mean_DistMat_ClustPair_NJ, \
				'Norm_Mean_DistMat_ClustPair_NJ', Output_Text_File)
			if 0:	#(METHOD_USED == NJSTXL) or (METHOD_USED == MedNJSTXL) or (METHOD_USED == ModeNJSTXL):
				PrintMatrixContent(no_of_taxa_clust, clust_species_list, Norm_XLVal_DistMat_ClustPair_NJ, \
					'Norm_XLVal_DistMat_ClustPair_NJ', Output_Text_File)
		
		"""
		find the cluster pairs having minimum distance values
		
		here we have used agglomerative clustering using normalized matrix entries
		
		here both Norm_Mean_DistMat_ClustPair_NJ and 
		Norm_XLVal_DistMat_ClustPair_NJ matrices 
		contain positive but small fraction entries
		so we have to find the minimum among these elements 
		(product of these two matrices in position specific way)
		and find out the corresponding indices
		
		the minimum is stored in the element target_val
		we have used "<" operator to find out the minimum of the elements
		"""
		if (METHOD_USED == ProdNJSTXL):
			min_idx_i, min_idx_j = Find_Unique_Min_XL(Mean_DistMat_ClustPair_NJ, Norm_Mean_DistMat_ClustPair_NJ, \
				XLVal_DistMat_ClustPair_NJ, Norm_XLVal_DistMat_ClustPair_NJ, no_of_taxa_clust, clust_species_list, NJ_RULE_USED)
		else:
			min_idx_i, min_idx_j = Find_Unique_Min_RankBased(Mean_DistMat_ClustPair_NJ, \
				Norm_Mean_DistMat_ClustPair_NJ, XLVal_DistMat_ClustPair_NJ, \
					Norm_XLVal_DistMat_ClustPair_NJ, no_of_taxa_clust, clust_species_list, \
						Output_Text_File, RankMergeMethod)
		
		
		#---------------------------------------------------------------      
		# note down the taxa list in these two indices (min_idx_i and min_idx_j) 
		# of the clust_species_list
		taxa_list = []
		for x in clust_species_list[min_idx_i]:
			taxa_list.append(x)
		for x in clust_species_list[min_idx_j]:
			taxa_list.append(x)

		if (DEBUG_LEVEL >= 2):
			fp = open(Output_Text_File, 'a')
			fp.write('\n min_idx_i ' + str(min_idx_i) + ' min_idx_j : ' + str(min_idx_j))
			fp.write('\n min_idx_i species list ' + str(clust_species_list[min_idx_i]))
			fp.write('\n min_idx_j species list ' + str(clust_species_list[min_idx_j]))
			fp.write('\n complete taxa list (union) ' + str(taxa_list))
			fp.close()
		#----------------------------------------------------
		"""
		now we merge the pair of clusters pointed by these indices
		"""
		Curr_tree = Merge_Cluster_Pair(Curr_tree, clust_species_list, min_idx_i, min_idx_j, taxa_list, Output_Text_File)
		#----------------------------------------------------
		# adjust the Mean_DistMat_ClustPair_NJ by inserting one new row and column corresponding to the new cluster
		# and then deleting the information of earlier two clusters

		# first append one row
		Mean_DistMat_ClustPair_NJ = numpy.vstack((Mean_DistMat_ClustPair_NJ, \
			numpy.zeros((1, no_of_taxa_clust), dtype=numpy.float)))
		# then append one column
		Mean_DistMat_ClustPair_NJ = numpy.hstack((Mean_DistMat_ClustPair_NJ, \
			numpy.zeros((no_of_taxa_clust + 1, 1), dtype=numpy.float)))
		
		# now reshape the distance matrix
		Mean_DistMat_ClustPair_NJ = numpy.reshape(Mean_DistMat_ClustPair_NJ, \
			((no_of_taxa_clust + 1), (no_of_taxa_clust + 1)), order='C')
		
		# apply these operations on the excess gene count matrix as well
		XLVal_DistMat_ClustPair_NJ = numpy.vstack((XLVal_DistMat_ClustPair_NJ, \
			numpy.zeros((1, no_of_taxa_clust), dtype=numpy.float)))
		XLVal_DistMat_ClustPair_NJ = numpy.hstack((XLVal_DistMat_ClustPair_NJ, \
			numpy.zeros((no_of_taxa_clust + 1, 1), dtype=numpy.float)))
		XLVal_DistMat_ClustPair_NJ = numpy.reshape(XLVal_DistMat_ClustPair_NJ, \
			((no_of_taxa_clust + 1), (no_of_taxa_clust + 1)), order='C')
		
		# add taxa_list as a new element of clust_species_list
		clust_species_list.append(taxa_list)          
		
		# now recompute the entries of this new row and column (which is indexed by no_of_taxa_clust), according to the NJ principle
		# compute Mean_DistMat_ClustPair_NJ[no_of_taxa_clust][m] entries where m != min_idx_i and m != min_idx_j
		for m in range(no_of_taxa_clust):
			if (m == min_idx_i) or (m == min_idx_j):
				continue
			Mean_DistMat_ClustPair_NJ[no_of_taxa_clust][m] = (Mean_DistMat_ClustPair_NJ[min_idx_i][m] + \
				Mean_DistMat_ClustPair_NJ[min_idx_j][m] - Mean_DistMat_ClustPair_NJ[min_idx_i][min_idx_j]) / 2.0
			Mean_DistMat_ClustPair_NJ[m][no_of_taxa_clust] = Mean_DistMat_ClustPair_NJ[no_of_taxa_clust][m]

			if 1:	#(DIST_MAT_UPDATE == 1):
				XLVal_DistMat_ClustPair_NJ[no_of_taxa_clust][m] = \
					(XLVal_DistMat_ClustPair_NJ[min_idx_i][m] + XLVal_DistMat_ClustPair_NJ[min_idx_j][m]) / 2.0
			else:
				XLVal_DistMat_ClustPair_NJ[no_of_taxa_clust][m] = max(XLVal_DistMat_ClustPair_NJ[min_idx_i][m], \
					XLVal_DistMat_ClustPair_NJ[min_idx_j][m])
				#XLVal_DistMat_ClustPair_NJ[no_of_taxa_clust][m] = max(XLVal_DistMat_ClustPair_NJ[min_idx_i][m], \
					#XLVal_DistMat_ClustPair_NJ[min_idx_j][m], XLVal_DistMat_ClustPair_NJ[min_idx_i][min_idx_j])
				
			# symmetric property
			XLVal_DistMat_ClustPair_NJ[m][no_of_taxa_clust] = XLVal_DistMat_ClustPair_NJ[no_of_taxa_clust][m]

		# now remove the rows and columns corresponding to min_idx_i and min_idx_j
		Mean_DistMat_ClustPair_NJ = numpy.delete(Mean_DistMat_ClustPair_NJ, (min_idx_i), axis=0)	# delete the row
		Mean_DistMat_ClustPair_NJ = numpy.delete(Mean_DistMat_ClustPair_NJ, (min_idx_i), axis=1)	# delete the column
		Mean_DistMat_ClustPair_NJ = numpy.delete(Mean_DistMat_ClustPair_NJ, (min_idx_j - 1), axis=0)	# delete the row
		Mean_DistMat_ClustPair_NJ = numpy.delete(Mean_DistMat_ClustPair_NJ, (min_idx_j - 1), axis=1)	# delete the column

		# update the LCA rank matrix as well
		XLVal_DistMat_ClustPair_NJ = numpy.delete(XLVal_DistMat_ClustPair_NJ, (min_idx_i), axis=0)	# delete the row
		XLVal_DistMat_ClustPair_NJ = numpy.delete(XLVal_DistMat_ClustPair_NJ, (min_idx_i), axis=1)	# delete the column
		XLVal_DistMat_ClustPair_NJ = numpy.delete(XLVal_DistMat_ClustPair_NJ, (min_idx_j - 1), axis=0)	# delete the row
		XLVal_DistMat_ClustPair_NJ = numpy.delete(XLVal_DistMat_ClustPair_NJ, (min_idx_j - 1), axis=1)	# delete the column

		# clear Norm_Mean_DistMat_ClustPair_NJ
		Norm_Mean_DistMat_ClustPair_NJ = numpy.delete(Norm_Mean_DistMat_ClustPair_NJ, (min_idx_i), axis=0)	# delete the row
		Norm_Mean_DistMat_ClustPair_NJ = numpy.delete(Norm_Mean_DistMat_ClustPair_NJ, (min_idx_i), axis=1)	# delete the column
		Norm_Mean_DistMat_ClustPair_NJ.fill(0)
		
		# clear norm LCA matrix
		Norm_XLVal_DistMat_ClustPair_NJ = numpy.delete(Norm_XLVal_DistMat_ClustPair_NJ, (min_idx_i), axis=0)	# delete the row
		Norm_XLVal_DistMat_ClustPair_NJ = numpy.delete(Norm_XLVal_DistMat_ClustPair_NJ, (min_idx_i), axis=1)	# delete the column
		Norm_XLVal_DistMat_ClustPair_NJ.fill(0)    
		
		# remove individual clusters' taxa information from the clust_species_list
		clust_species_list.pop(min_idx_i)
		clust_species_list.pop(min_idx_j - 1)
		
		# decrement the number of clusters considered
		no_of_taxa_clust = no_of_taxa_clust - 1

	return
