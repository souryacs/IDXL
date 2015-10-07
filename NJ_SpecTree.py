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

##---------------------------------------------
""" 
function to print the matrix content
N is the matrix dimension
"""
def PrintMatrixContent(N, TaxaList, inp_data, inp_str, textfile):
	fp = open(textfile, 'a')
	fp.write('\n printing contents of ' + str(inp_str) + ' ---- ')
	for i in range(N):
		fp.write('\n ' + str(i) + '--' + str(TaxaList[i]) + '--->>')
		for j in range(i+1):
			fp.write(' ' + str(inp_data[i][j]))
	fp.close()

##---------------------------------------------
"""
this function refines the initial species tree (in terms of a star network) to 
find the true species tree
it does using agglomerative clustering (NJ principle)
the distance metric employed for NJ algorithm can vary depending on experimentation 
"""
def Form_Species_Tree_NJ_Cluster(Star_Tree_Initial, METHOD_USED, NJ_RULE_USED, Output_Text_File):

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

	# for individual cluster pairs, we compute the sum of extra lineages
	if (DEBUG_LEVEL >= 2):
		fp = open(Output_Text_File, 'a')
		fp.write('\n COMPLETE_INPUT_TAXA_LIST ' + str(COMPLETE_INPUT_TAXA_LIST))
		fp.write('\n Initial formed clust_species_list ' + str(clust_species_list))
		fp.close()        

	# allocate a 2D square matrix of no_of_taxa_clust dimension
	# for a pair of taxa clusters Cx and Cy, it contains the employed main distance metric for the cluster pairs
	if (METHOD_USED == NJ_ST) or (METHOD_USED == M_NJ_ST):
		Mean_DistMat_ClustPair_NJ = numpy.zeros((no_of_taxa_clust, no_of_taxa_clust), dtype=numpy.float)
	else:
		Mean_DistMat_ClustPair_NJ = numpy.zeros((no_of_taxa_clust, no_of_taxa_clust), dtype=numpy.int)

	# allocate one new square matrix which will contain the NJ based modified distance matrix (used for minimum finding routine)
	Norm_Mean_DistMat_ClustPair_NJ = numpy.zeros((no_of_taxa_clust, no_of_taxa_clust), dtype=numpy.float)

	if (METHOD_USED == M_NJ_ST_XL):
		# we allocate matrices containing LCA rank among individual couplets
		XLVal_DistMat_ClustPair_NJ = numpy.zeros((no_of_taxa_clust, no_of_taxa_clust), dtype=numpy.int)
		# allocate one new square matrix which will contain the normalized mean LCA rank
		# used in NJ iterations
		Norm_XLVal_DistMat_ClustPair_NJ = numpy.zeros((no_of_taxa_clust, no_of_taxa_clust), dtype=numpy.float)

	# now fill the Mean_DistMat_ClustPair_NJ according to the specified metric used for NJ like clustering
	for l in TaxaPair_Reln_Dict:
		spec1 = l[0]
		spec2 = l[1]
		spec1_idx = COMPLETE_INPUT_TAXA_LIST.index(spec1)
		spec2_idx = COMPLETE_INPUT_TAXA_LIST.index(spec2)
		# fill the distance matrix entries according to the distance metric and the statistics used (input parameters) 
		if (METHOD_USED == NJ_ST):
			Mean_DistMat_ClustPair_NJ[spec1_idx][spec2_idx] = Mean_DistMat_ClustPair_NJ[spec2_idx][spec1_idx] = TaxaPair_Reln_Dict[l]._GetAvgSumLevel()
		elif (METHOD_USED == M_NJ_ST) or (METHOD_USED == M_NJ_ST_XL):
			Mean_DistMat_ClustPair_NJ[spec1_idx][spec2_idx] = Mean_DistMat_ClustPair_NJ[spec2_idx][spec1_idx] = TaxaPair_Reln_Dict[l]._GetMultiModeSumLevel()
		
		if (METHOD_USED == M_NJ_ST_XL):
			# mean of XL value
			XLVal_DistMat_ClustPair_NJ[spec2_idx][spec1_idx] = TaxaPair_Reln_Dict[l]._GetAvgXLVal()
			XLVal_DistMat_ClustPair_NJ[spec1_idx][spec2_idx] = XLVal_DistMat_ClustPair_NJ[spec2_idx][spec1_idx]
					
	#--------------------------------------------------------
	# loop to execute the agglomerative clustering
	while(no_of_taxa_clust > 2): 
		if (DEBUG_LEVEL >= 2):
			fp = open(Output_Text_File, 'a')
			fp.write('\n iteration start --- number of clusters: ' + str(no_of_taxa_clust))
			fp.write('\n clust_species_list : ' + str(clust_species_list))
			fp.close()
			PrintMatrixContent(no_of_taxa_clust, clust_species_list, Mean_DistMat_ClustPair_NJ, 'Mean_DistMat_ClustPair_NJ', Output_Text_File)
			if (METHOD_USED == M_NJ_ST_XL):
				PrintMatrixContent(no_of_taxa_clust, clust_species_list, XLVal_DistMat_ClustPair_NJ, 'XLVal_DistMat_ClustPair_NJ', Output_Text_File)
					
		# for individual cluster Cx, it contains XL(Cx, :) - sum of extra lineages considering the cluster pair 
		# (Cx, Cy) for all other clusters Cy
		sum_DistMat_Clust = []
		if (METHOD_USED == M_NJ_ST_XL):
			sum_XLVal_Clust = []
		
		for i in range(no_of_taxa_clust):
			t1 = 0
			t2 = 0
			for j in range(no_of_taxa_clust):
				t1 = t1 + Mean_DistMat_ClustPair_NJ[i][j]
				if (METHOD_USED == M_NJ_ST_XL):
					t2 = t2 + XLVal_DistMat_ClustPair_NJ[i][j]
			sum_DistMat_Clust.append(t1)
			if (METHOD_USED == M_NJ_ST_XL):
				sum_XLVal_Clust.append(t2)
			
		if (DEBUG_LEVEL > 2):
			fp = open(Output_Text_File, 'a')
			fp.write('\n content of sum_DistMat_Clust : ' + str(sum_DistMat_Clust))
			if (METHOD_USED == M_NJ_ST_XL):
				fp.write('\n content of sum_XLVal_Clust : ' + str(sum_XLVal_Clust))
			fp.close()
			
		# fill the normalized matrix, which will be used for clustering
		for i in range(no_of_taxa_clust - 1):
			for j in range(i+1, no_of_taxa_clust):
				if (NJ_RULE_USED == AGGLO_CLUST):
					Norm_Mean_DistMat_ClustPair_NJ[i][j] = ObtainNormalizedVal(Mean_DistMat_ClustPair_NJ[i][j], sum_DistMat_Clust[i], sum_DistMat_Clust[j])
					Norm_Mean_DistMat_ClustPair_NJ[j][i] = Norm_Mean_DistMat_ClustPair_NJ[i][j]
					if (METHOD_USED == M_NJ_ST_XL):
						# update the normalized LCA rank matrix
						Norm_XLVal_DistMat_ClustPair_NJ[i][j] = ObtainNormalizedVal(XLVal_DistMat_ClustPair_NJ[i][j], sum_XLVal_Clust[i], sum_XLVal_Clust[j])
						Norm_XLVal_DistMat_ClustPair_NJ[j][i] = Norm_XLVal_DistMat_ClustPair_NJ[i][j]
				else:	
					# here ri , rj are the sum of all distances
					# standard NJ method implementation
					ri = sum_DistMat_Clust[i] / (no_of_taxa_clust - 2)
					rj = sum_DistMat_Clust[j] / (no_of_taxa_clust - 2)
					Norm_Mean_DistMat_ClustPair_NJ[i][j] = (Mean_DistMat_ClustPair_NJ[i][j] - ri - rj)
					Norm_Mean_DistMat_ClustPair_NJ[j][i] = Norm_Mean_DistMat_ClustPair_NJ[i][j]
					if (METHOD_USED == M_NJ_ST_XL):
						# update the normalized LCA rank matrix
						ri1 = sum_XLVal_Clust[i] / (no_of_taxa_clust - 2)
						rj1 = sum_XLVal_Clust[j] / (no_of_taxa_clust - 2)
						Norm_XLVal_DistMat_ClustPair_NJ[i][j] = (XLVal_DistMat_ClustPair_NJ[i][j] - ri1 - rj1)
						Norm_XLVal_DistMat_ClustPair_NJ[j][i] = Norm_XLVal_DistMat_ClustPair_NJ[i][j]

		if (DEBUG_LEVEL >= 2):
			PrintMatrixContent(no_of_taxa_clust, clust_species_list, Norm_Mean_DistMat_ClustPair_NJ, 'Norm_Mean_DistMat_ClustPair_NJ', Output_Text_File)
			if (METHOD_USED == M_NJ_ST_XL):
				PrintMatrixContent(no_of_taxa_clust, clust_species_list, Norm_XLVal_DistMat_ClustPair_NJ, 'Norm_XLVal_DistMat_ClustPair_NJ', Output_Text_File)
		
		#---------------------------------------------------------------  
		# add - sourya
		if (NJ_RULE_USED == AGGLO_CLUST):
			""" 
			here we have used agglomerative clustering using normalized matrix entries
			
			here both Norm_Mean_DistMat_ClustPair_NJ and 
			Norm_VarianceAccRank_DistMat_ClustPair_NJ /  Norm_XLVal_DistMat_ClustPair_NJ matrices 
			contain positive but small fraction entries
			so we have to find the minimum among these elements 
			(product of these two matrices in position specific way)
			and find out the corresponding indices
			
			the minimum is stored in the element target_val
			we have used "<" operator to find out the minimum of the elements
			"""
			if (METHOD_USED == M_NJ_ST_XL):
				target_val = (Norm_Mean_DistMat_ClustPair_NJ[0][1] * Norm_XLVal_DistMat_ClustPair_NJ[0][1])
				min_idx_i = 0
				min_idx_j = 1
				for i in range(no_of_taxa_clust - 1):
					for j in range(i+1, no_of_taxa_clust):
						if (i == j):
							continue
						if ((Norm_Mean_DistMat_ClustPair_NJ[i][j] * Norm_XLVal_DistMat_ClustPair_NJ[i][j]) < target_val):
							target_val = (Norm_Mean_DistMat_ClustPair_NJ[i][j] * Norm_XLVal_DistMat_ClustPair_NJ[i][j])
							min_idx_i = i
							min_idx_j = j
						elif ((Norm_Mean_DistMat_ClustPair_NJ[i][j] * Norm_XLVal_DistMat_ClustPair_NJ[i][j]) == target_val):
							# here we prioritize the cluster pair having minimum number of species
							if (len(clust_species_list[i]) + len(clust_species_list[j])) < (len(clust_species_list[min_idx_i]) + len(clust_species_list[min_idx_j])):
								min_idx_i = i
								min_idx_j = j
			
			else:
				target_val = Norm_Mean_DistMat_ClustPair_NJ[0][1]
				min_idx_i = 0
				min_idx_j = 1
				for i in range(no_of_taxa_clust - 1):
					for j in range(i+1, no_of_taxa_clust):
						if (i == j):
							continue
						if (Norm_Mean_DistMat_ClustPair_NJ[i][j] < target_val):
							target_val = Norm_Mean_DistMat_ClustPair_NJ[i][j]
							min_idx_i = i
							min_idx_j = j
						elif (Norm_Mean_DistMat_ClustPair_NJ[i][j] == target_val):
							# here we prioritize the cluster pair having minimum number of species
							if (len(clust_species_list[i]) + len(clust_species_list[j])) < (len(clust_species_list[min_idx_i]) + len(clust_species_list[min_idx_j])):
								min_idx_i = i
								min_idx_j = j
		
		else:	# NJ_RULE_USED = TRADITIONAL_NJ 
			""" 
			here we have used agglomerative clustering using standard NJ rules
			
			here both Norm_Mean_DistMat_ClustPair_NJ and 
			Norm_VarianceAccRank_DistMat_ClustPair_NJ / Norm_XLVal_DistMat_ClustPair_NJ matrices 
			contain negative entries (absolute value will be very high, so as to obtain very high negative values)
			so we have to find the maximum among the product matrix elements (individual high negative entries)
			(product of these two matrices in position specific way)
			and find out the corresponding indices
			
			the maximum is stored in the element target_val
			we have used ">" operator to find out the maximum of the elements
			"""
			if (METHOD_USED == M_NJ_ST_XL):
				target_val = (Norm_Mean_DistMat_ClustPair_NJ[0][1] * Norm_XLVal_DistMat_ClustPair_NJ[0][1])
				min_idx_i = 0
				min_idx_j = 1
				for i in range(no_of_taxa_clust - 1):
					for j in range(i+1, no_of_taxa_clust):
						if (i == j):
							continue
						if ((Norm_Mean_DistMat_ClustPair_NJ[i][j] * Norm_XLVal_DistMat_ClustPair_NJ[i][j]) > target_val):
							target_val = (Norm_Mean_DistMat_ClustPair_NJ[i][j] * Norm_XLVal_DistMat_ClustPair_NJ[i][j])
							min_idx_i = i
							min_idx_j = j
			
			else:
				target_val = Norm_Mean_DistMat_ClustPair_NJ[0][1]
				min_idx_i = 0
				min_idx_j = 1
				for i in range(no_of_taxa_clust - 1):
					for j in range(i+1, no_of_taxa_clust):
						if (i == j):
							continue
						if (Norm_Mean_DistMat_ClustPair_NJ[i][j] < target_val):
							target_val = Norm_Mean_DistMat_ClustPair_NJ[i][j]
							min_idx_i = i
							min_idx_j = j
		
		# end add - sourya
		#---------------------------------------------------------------      
		
		# note down the taxa list in these two indices (min_idx_i and min_idx_j) of the clust_species_list
		taxa_list = []
		for x in clust_species_list[min_idx_i]:
			taxa_list.append(x)
		for x in clust_species_list[min_idx_j]:
			taxa_list.append(x)

		if (DEBUG_LEVEL >= 2):
			fp = open(Output_Text_File, 'a')
			fp.write('\n min_idx_i ' + str(min_idx_i) + ' min_idx_j : ' + str(min_idx_j) + 'min val : ' + str(target_val))
			fp.write('\n min index - mean branch count: ' + str(Mean_DistMat_ClustPair_NJ[min_idx_i][min_idx_j]))
			if (METHOD_USED == M_NJ_ST_XL):
				fp.write('\n target index - LCA rank: ' + str(XLVal_DistMat_ClustPair_NJ[min_idx_i][min_idx_j]))
			fp.write('\n min_idx_i species list ' + str(clust_species_list[min_idx_i]))
			fp.write('\n min_idx_j species list ' + str(clust_species_list[min_idx_j]))
			fp.write('\n complete taxa list (union) ' + str(taxa_list))
			fp.close()

		#---------------------------------------------------------      
		# for individual clusters, we check if the cluster contains one or more species
		# case 1 - both the clusters have > 1 species
		# and the clusters are represented by an internal node which is the MRCA of the constituent species set
		if (len(clust_species_list[min_idx_i]) > 1) and (len(clust_species_list[min_idx_j]) > 1):
			# comment - sourya
			first_cluster_mrca_node = Star_Tree_Initial.mrca(taxon_labels=clust_species_list[min_idx_i])
			second_cluster_mrca_node = Star_Tree_Initial.mrca(taxon_labels=clust_species_list[min_idx_j])
			all_taxa_mrca_node = Star_Tree_Initial.mrca(taxon_labels=taxa_list)
			## add - sourya
			#first_cluster_mrca_node = Find_MRCA(Star_Tree_Initial, clust_species_list[min_idx_i])
			#second_cluster_mrca_node = Find_MRCA(Star_Tree_Initial, clust_species_list[min_idx_j])
			#all_taxa_mrca_node = Find_MRCA(Star_Tree_Initial, taxa_list)
			## end add - sourya
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
			Star_Tree_Initial.update_splits(delete_outdegree_one=False)
			
		# case 2 and 3 - one cluster has at least 2 species, while other is a leaf
		elif (len(clust_species_list[min_idx_i]) == 1) and (len(clust_species_list[min_idx_j]) > 1):
			first_cluster_leaf_node = Star_Tree_Initial.find_node_with_taxon_label(clust_species_list[min_idx_i][0])
			# comment - sourya
			second_cluster_mrca_node = Star_Tree_Initial.mrca(taxon_labels=clust_species_list[min_idx_j])
			all_taxa_mrca_node = Star_Tree_Initial.mrca(taxon_labels=taxa_list)
			# add - sourya
			#second_cluster_mrca_node = Find_MRCA(Star_Tree_Initial, clust_species_list[min_idx_j])
			#all_taxa_mrca_node = Find_MRCA(Star_Tree_Initial, taxa_list)
			# end add - sourya
			if (DEBUG_LEVEL >= 2):
				fp = open(Output_Text_File, 'a')
				fp.write('\n first cluster is a leaf - its label: ' + str(Node_Label(first_cluster_leaf_node)))      
				fp.write('\n label of second_cluster_mrca_node: ' + str(Node_Label(second_cluster_mrca_node)))
				fp.write('\n label of all_taxa_mrca_node: ' + str(Node_Label(all_taxa_mrca_node)))
				fp.close()
			
			# create new internal node 
			newnode = dendropy.Node()
			# its parent node will be the previous MRCA node of all the taxa in two clusters
			all_taxa_mrca_node.add_child(newnode)
			newnode.parent_node = all_taxa_mrca_node
			all_taxa_mrca_node.remove_child(first_cluster_leaf_node)
			first_cluster_leaf_node.parent_node = None
			all_taxa_mrca_node.remove_child(second_cluster_mrca_node)
			second_cluster_mrca_node.parent_node = None
			# add these individual clusters' MRCA node as its children
			newnode.add_child(first_cluster_leaf_node)
			first_cluster_leaf_node.parent_node = newnode
			newnode.add_child(second_cluster_mrca_node)
			second_cluster_mrca_node.parent_node = newnode
			# update splits of the resulting tree
			Star_Tree_Initial.update_splits(delete_outdegree_one=False)
			
		elif (len(clust_species_list[min_idx_i]) > 1) and (len(clust_species_list[min_idx_j]) == 1):
			# comment - sourya
			first_cluster_mrca_node = Star_Tree_Initial.mrca(taxon_labels=clust_species_list[min_idx_i])
			# add - sourya
			#first_cluster_mrca_node = Find_MRCA(Star_Tree_Initial, clust_species_list[min_idx_i])
			# end add - sourya
			second_cluster_leaf_node = Star_Tree_Initial.find_node_with_taxon_label(clust_species_list[min_idx_j][0])
			# comment - sourya
			all_taxa_mrca_node = Star_Tree_Initial.mrca(taxon_labels=taxa_list)
			# add - sourya
			#all_taxa_mrca_node = Find_MRCA(Star_Tree_Initial, taxa_list)
			# end add - sourya
			if (DEBUG_LEVEL >= 2):
				fp = open(Output_Text_File, 'a')
				fp.write('\n label of first_cluster_mrca_node: ' + str(Node_Label(first_cluster_mrca_node)))      
				fp.write('\n label of second_cluster_mrca_node: ' + str(Node_Label(second_cluster_leaf_node)))
				fp.write('\n label of all_taxa_mrca_node: ' + str(Node_Label(all_taxa_mrca_node)))
				fp.close()
			
			# create new internal node 
			newnode = dendropy.Node()
			# its parent node will be the previous MRCA node of all the taxa in two clusters
			all_taxa_mrca_node.add_child(newnode)
			newnode.parent_node = all_taxa_mrca_node
			all_taxa_mrca_node.remove_child(first_cluster_mrca_node)
			first_cluster_mrca_node.parent_node = None
			all_taxa_mrca_node.remove_child(second_cluster_leaf_node)
			second_cluster_leaf_node.parent_node = None
			# add these individual clusters' MRCA node as its children
			newnode.add_child(first_cluster_mrca_node)
			first_cluster_mrca_node.parent_node = newnode
			newnode.add_child(second_cluster_leaf_node)
			second_cluster_leaf_node.parent_node = newnode
			# update splits of the resulting tree
			Star_Tree_Initial.update_splits(delete_outdegree_one=False)
			
		# case 4 - when both child clusters are leaf nodes 
		else:
			first_cluster_leaf_node = Star_Tree_Initial.find_node_with_taxon_label(clust_species_list[min_idx_i][0])
			second_cluster_leaf_node = Star_Tree_Initial.find_node_with_taxon_label(clust_species_list[min_idx_j][0])
			# comment - sourya
			all_taxa_mrca_node = Star_Tree_Initial.mrca(taxon_labels=taxa_list)
			# add - sourya
			#all_taxa_mrca_node = Find_MRCA(Star_Tree_Initial, taxa_list)
			# end add - sourya      
			if (DEBUG_LEVEL >= 2):
				fp = open(Output_Text_File, 'a')
				fp.write('\n first cluster is a leaf - its label: ' + str(Node_Label(first_cluster_leaf_node)))      
				fp.write('\n second cluster is a leaf - its label: ' + str(Node_Label(second_cluster_leaf_node)))
				fp.write('\n label of all_taxa_mrca_node: ' + str(Node_Label(all_taxa_mrca_node)))
				fp.close()
			
			# create new internal node 
			newnode = dendropy.Node()
			# its parent node will be the previous MRCA node of all the taxa in two clusters
			all_taxa_mrca_node.add_child(newnode)
			newnode.parent_node = all_taxa_mrca_node
			all_taxa_mrca_node.remove_child(first_cluster_leaf_node)
			first_cluster_leaf_node.parent_node = None
			all_taxa_mrca_node.remove_child(second_cluster_leaf_node)
			second_cluster_leaf_node.parent_node = None
			# add these individual clusters' MRCA node as its children
			newnode.add_child(first_cluster_leaf_node)
			first_cluster_leaf_node.parent_node = newnode
			newnode.add_child(second_cluster_leaf_node)
			second_cluster_leaf_node.parent_node = newnode
			# update splits of the resulting tree
			Star_Tree_Initial.update_splits(delete_outdegree_one=False)

		#---------------------------------------------------------------------
		if (DEBUG_LEVEL >= 2):
			fp = open(Output_Text_File, 'a')
			fp.write('\n label of newnode: ' + str(Node_Label(newnode)))
			fp.write('\n label of all taxa mrca node (recomputed): ' + str(Node_Label(all_taxa_mrca_node)))      
			#fp.write('\n before inserting row col, Mean_DistMat_ClustPair_NJ dimension: ' + str(Mean_DistMat_ClustPair_NJ.size))
			fp.close()

		# adjust the Mean_DistMat_ClustPair_NJ by inserting one new row and column corresponding to the new cluster
		# and then deleting the information of earlier two clusters
		if (METHOD_USED == NJ_ST) or (METHOD_USED == M_NJ_ST):
			# first append one row
			Mean_DistMat_ClustPair_NJ = numpy.vstack((Mean_DistMat_ClustPair_NJ, numpy.zeros((1, no_of_taxa_clust), dtype=numpy.float)))
			# then append one column
			Mean_DistMat_ClustPair_NJ = numpy.hstack((Mean_DistMat_ClustPair_NJ, numpy.zeros((no_of_taxa_clust + 1, 1), dtype=numpy.float)))
		else:
			# first append one row
			Mean_DistMat_ClustPair_NJ = numpy.vstack((Mean_DistMat_ClustPair_NJ, numpy.zeros((1, no_of_taxa_clust), dtype=numpy.int)))
			# then append one column
			Mean_DistMat_ClustPair_NJ = numpy.hstack((Mean_DistMat_ClustPair_NJ, numpy.zeros((no_of_taxa_clust + 1, 1), dtype=numpy.int)))
		
		# now reshape the distance matrix
		Mean_DistMat_ClustPair_NJ = numpy.reshape(Mean_DistMat_ClustPair_NJ, ((no_of_taxa_clust + 1), (no_of_taxa_clust + 1)), order='C')
		
		#if (DEBUG_LEVEL > 2):
			#fp = open(Output_Text_File, 'a')
			#fp.write('\n after inserting row col, Mean_DistMat_ClustPair_NJ dimension: ' + str(Mean_DistMat_ClustPair_NJ.size))
			#fp.close()
		
		if (METHOD_USED == M_NJ_ST_XL):
			# apply these operations on the LCA matrix as well
			XLVal_DistMat_ClustPair_NJ = numpy.vstack((XLVal_DistMat_ClustPair_NJ, numpy.zeros((1, no_of_taxa_clust), dtype=numpy.int)))
			XLVal_DistMat_ClustPair_NJ = numpy.hstack((XLVal_DistMat_ClustPair_NJ, numpy.zeros((no_of_taxa_clust + 1, 1), dtype=numpy.int)))
			XLVal_DistMat_ClustPair_NJ = numpy.reshape(XLVal_DistMat_ClustPair_NJ, ((no_of_taxa_clust + 1), (no_of_taxa_clust + 1)), order='C')
		
		# add taxa_list as a new element of clust_species_list
		clust_species_list.append(taxa_list)          
		
		# now recompute the entries of this new row and column (which is indexed by no_of_taxa_clust), according to the NJ principle
		# compute Mean_DistMat_ClustPair_NJ[no_of_taxa_clust][m] entries where m != min_idx_i and m != min_idx_j
		for m in range(no_of_taxa_clust):
			if (m == min_idx_i) or (m == min_idx_j):
				continue
			# comment - sourya
			Mean_DistMat_ClustPair_NJ[no_of_taxa_clust][m] = (Mean_DistMat_ClustPair_NJ[min_idx_i][m] + Mean_DistMat_ClustPair_NJ[min_idx_j][m] - Mean_DistMat_ClustPair_NJ[min_idx_i][min_idx_j]) / 2
			Mean_DistMat_ClustPair_NJ[m][no_of_taxa_clust] = Mean_DistMat_ClustPair_NJ[no_of_taxa_clust][m]

			if (METHOD_USED == M_NJ_ST_XL):
				## comment - sourya
				#XLVal_DistMat_ClustPair_NJ[no_of_taxa_clust][m] = (XLVal_DistMat_ClustPair_NJ[min_idx_i][m] + XLVal_DistMat_ClustPair_NJ[min_idx_j][m] - XLVal_DistMat_ClustPair_NJ[min_idx_i][min_idx_j]) / 2.0
				# add - sourya
				XLVal_DistMat_ClustPair_NJ[no_of_taxa_clust][m] = max(XLVal_DistMat_ClustPair_NJ[min_idx_i][m], XLVal_DistMat_ClustPair_NJ[min_idx_j][m]) - 1
				#XLVal_DistMat_ClustPair_NJ[no_of_taxa_clust][m] = max(XLVal_DistMat_ClustPair_NJ[min_idx_i][m], XLVal_DistMat_ClustPair_NJ[min_idx_j][m], XLVal_DistMat_ClustPair_NJ[min_idx_i][min_idx_j]) - 1
				# end add - sourya
				XLVal_DistMat_ClustPair_NJ[m][no_of_taxa_clust] = XLVal_DistMat_ClustPair_NJ[no_of_taxa_clust][m]

		# now remove the rows and columns corresponding to min_idx_i and min_idx_j
		Mean_DistMat_ClustPair_NJ = numpy.delete(Mean_DistMat_ClustPair_NJ, (min_idx_i), axis=0)	# delete the row
		Mean_DistMat_ClustPair_NJ = numpy.delete(Mean_DistMat_ClustPair_NJ, (min_idx_i), axis=1)	# delete the column
		Mean_DistMat_ClustPair_NJ = numpy.delete(Mean_DistMat_ClustPair_NJ, (min_idx_j - 1), axis=0)	# delete the row
		Mean_DistMat_ClustPair_NJ = numpy.delete(Mean_DistMat_ClustPair_NJ, (min_idx_j - 1), axis=1)	# delete the column

		if (METHOD_USED == M_NJ_ST_XL):
			# update the LCA rank matrix as well
			XLVal_DistMat_ClustPair_NJ = numpy.delete(XLVal_DistMat_ClustPair_NJ, (min_idx_i), axis=0)	# delete the row
			XLVal_DistMat_ClustPair_NJ = numpy.delete(XLVal_DistMat_ClustPair_NJ, (min_idx_i), axis=1)	# delete the column
			XLVal_DistMat_ClustPair_NJ = numpy.delete(XLVal_DistMat_ClustPair_NJ, (min_idx_j - 1), axis=0)	# delete the row
			XLVal_DistMat_ClustPair_NJ = numpy.delete(XLVal_DistMat_ClustPair_NJ, (min_idx_j - 1), axis=1)	# delete the column

		# clear Norm_Mean_DistMat_ClustPair_NJ
		Norm_Mean_DistMat_ClustPair_NJ = numpy.delete(Norm_Mean_DistMat_ClustPair_NJ, (min_idx_i), axis=0)	# delete the row
		Norm_Mean_DistMat_ClustPair_NJ = numpy.delete(Norm_Mean_DistMat_ClustPair_NJ, (min_idx_i), axis=1)	# delete the column
		Norm_Mean_DistMat_ClustPair_NJ.fill(0)
		
		if (METHOD_USED == M_NJ_ST_XL):
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
