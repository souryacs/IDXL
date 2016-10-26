import Header
from Header import *
import UtilFunc
from UtilFunc import *

#--------------------------------------------------------
"""
this function is a shortcut to obtain the normalized expression 
used in the agglomerative clustering proposed in this code
as various methods are experimented, corresponding various forms of 
agglomerative clustering is tried
"""
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
	fp.write('\n\n\n ==>>>> printing contents of ' + str(inp_str) + ' ---- ')
	for i in range(N):
		fp.write('\n ' + str(i) + '--' + str(TaxaList[i]) + '--->>')
		if (i > 0):
			for j in range(i):
				fp.write(' ' + str(inp_data[i][j]))
	fp.close()

#----------------------------------------------------
"""
this function fills the input distance matrix using couplet based average internode count
"""
def Fill_DistMat_InternodeCount(DistMat, ntaxa):
	for i in range(ntaxa - 1):
		for j in range(i+1, ntaxa):
			"""
			target key of a couplet (if exists)
			"""
			l = (i, j)
			
			if l in TaxaPair_Reln_Dict:
				"""
				there exists a valid entry in the couplet dictionary
				"""
				DistMat[i][j] = TaxaPair_Reln_Dict[l]._GetAvgSumLevel()
				DistMat[j][i] = DistMat[i][j]
			else:
				"""
				there exists no valid couplet
				so, set the distance matrix values as (-1)
				"""
				DistMat[i][j] = -1
				DistMat[j][i] = DistMat[i][j]

	return

#----------------------------------------------------
"""
if excess gene count information is used, this 
function fills the distance matrix using either average or filtered average XL measure
depeding on the method and the type of distance matrix
"""
def Fill_DistMat_ExcessGeneCount(DistMat, DIST_MAT_TYPE, ntaxa):
	for i in range(ntaxa - 1):
		for j in range(i+1, ntaxa):
			"""
			target key of a couplet (if exists)
			"""
			l = (i, j)
			
			if l in TaxaPair_Reln_Dict:
				"""
				there exists a valid entry in the couplet dictionary
				"""
				if (DIST_MAT_TYPE == 1):
					DistMat[i][j] = TaxaPair_Reln_Dict[l]._GetAvgXLVal()
				elif (DIST_MAT_TYPE == 2):
					DistMat[i][j] = (TaxaPair_Reln_Dict[l]._GetAvgXLVal() + TaxaPair_Reln_Dict[l]._GetMultiModeXLVal()) / 2.0 
				elif (DIST_MAT_TYPE == 3):
					DistMat[i][j] = (TaxaPair_Reln_Dict[l]._GetAvgXLVal() + \
						TaxaPair_Reln_Dict[l]._MedianXLVal() + TaxaPair_Reln_Dict[l]._GetMultiModeXLVal()) / 3.0
				
				#symmetric property of the distance matrix
				DistMat[j][i] = DistMat[i][j]
			
			else:
				"""
				there exists no valid couplet
				so, set the distance matrix values as (-1)
				"""
				DistMat[i][j] = -1
				DistMat[j][i] = DistMat[i][j]	# symmetric property

	return

#---------------------------------------------
"""
computing the row wise sum for individual taxca clusters
"""
def ComputeSumRowsDistMat(sum_list, nclust, DistMat, Output_Text_File, inpstr):
	for i in range(nclust):
		t1 = 0
		for j in range(nclust):
			if (DistMat[i][j] >= 0):	# add the condition - sourya
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
			NormMat[i][j] = ObtainNormalizedVal(DistMat[i][j], (sum_list[i] - DistMat[i][j]), (sum_list[j] - DistMat[i][j]))
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
# ************************ Minimum finding for the method PNJSTXL / IDXL *********************
#---------------------------------------------
"""
this is the oldest version of the function 
results of this version are stored in the PNJSTXL folders of respective datasets
"""
"""
finds the couplet containing minimum entries of relative distance measures
when the product of internode count and excess gene count based measures are used
"""
def Find_Unique_Min_PNJSTXL_Version1(DM_ID, Norm_DM_ID, DM_XL, Norm_DM_XL, nclust, clust_species_list):
	
	target_val = (Norm_DM_ID[0][1] * Norm_DM_XL[0][1])
	min_idx_i = 0
	min_idx_j = 1
	for i in range(nclust - 1):
		for j in range(i+1, nclust):
			if (i == j):
				continue
			if (DM_ID[i][j] >= 0) and (DM_XL[i][j] >= 0):
				if ((Norm_DM_ID[i][j] * Norm_DM_XL[i][j]) > target_val):
					target_val = (Norm_DM_ID[i][j] * Norm_DM_XL[i][j])
					min_idx_i = i
					min_idx_j = j
				elif (FlEq((Norm_DM_ID[i][j] * Norm_DM_XL[i][j]), target_val) == True):
					"""
					equal value with respect to the earlier minimum
					here we agglomerate the clusters according to the following condition:
					1) if DM_ID[i][j] is strictly lower then use the new cluster pair
					2) if DM_ID[i][j] = DM_ID[min_idx_i][min_idx_j], 
					and DM_XL[i][j] is strictly lower then use the new cluster pair
					3) if both internode count and XL measures are identical for the above mentioned clusters, 
					agglomerate the new cluster pair, if they have higher sum of cardinality
					"""
					if (DM_ID[i][j] < DM_ID[min_idx_i][min_idx_j]):
						min_idx_i = i
						min_idx_j = j
					else:
						if (FlEq(DM_ID[i][j], DM_ID[min_idx_i][min_idx_j]) == True):
							if (DM_XL[i][j] < DM_XL[min_idx_i][min_idx_j]):
								min_idx_i = i
								min_idx_j = j
							else:
								if (FlEq(DM_XL[i][j], DM_XL[min_idx_i][min_idx_j]) == True):
									# higher sum of cardinality of the new cluster pair
									if (len(clust_species_list[i]) + len(clust_species_list[j])) \
										> (len(clust_species_list[min_idx_i]) + len(clust_species_list[min_idx_j])):
										min_idx_i = i
										min_idx_j = j
						
	return min_idx_i, min_idx_j

#---------------------------------------------
"""
this is the second version of the function 
here also, internode count and excess gene count based measures are used
but the positional rank of excess gene count is checked to remove a few error cases
"""
def Find_Unique_Min_PNJSTXL_Version2(DM_ID, Norm_DM_ID, DM_XL, Norm_DM_XL, nclust, clust_species_list):
	"""
	at first, use the DM_XL (containing the absolute values of excess gene leaf count for individual couplets)
	to create a sorted array
	"""
	Sorted_DM_XL = []

	for i in range(nclust - 1):
		for j in range(i+1, nclust):
			if (i == j):
				continue
			if (DM_ID[i][j] >= 0) and (DM_XL[i][j] >= 0):	# check for valid entry
				Sorted_DM_XL.append(DM_XL[i][j])
	
	"""
	find the unique elements of array
	we have found that using set based unique operation is highly efficient
	"""
	Sorted_DM_XL = list(set(Sorted_DM_XL))
	
	# sort the list in ascending order
	Sorted_DM_XL.sort()
	# no of elements in the sorted array
	len_Sorted_DM_XL = len(Sorted_DM_XL)
	
	"""
	this is a flag variable which is set when the first valid distance element is found
	"""
	tgt_found = 0
	
	for i in range(nclust - 1):
		for j in range(i+1, nclust):
			if (i == j):
				continue
			if (DM_ID[i][j] >= 0) and (DM_XL[i][j] >= 0):
				"""
				check the positional rank of the XL based distance matrix entry (absolute)
				"""
				curr_DM_XL = DM_XL[i][j]
				rank_curr_DM_XL = Sorted_DM_XL.index(DM_XL[i][j])
				if (rank_curr_DM_XL <= int(round(len_Sorted_DM_XL / 2))):	#condition - sourya
					if (tgt_found == 0):
						tgt_found = 1
						target_val = (Norm_DM_ID[i][j] * Norm_DM_XL[i][j])
						min_idx_i = i
						min_idx_j = j
					else:
						if ((Norm_DM_ID[i][j] * Norm_DM_XL[i][j]) > target_val):
							target_val = (Norm_DM_ID[i][j] * Norm_DM_XL[i][j])
							min_idx_i = i
							min_idx_j = j
						elif (FlEq((Norm_DM_ID[i][j] * Norm_DM_XL[i][j]), target_val) == True):
							"""
							equal value with respect to the earlier minimum
							here we agglomerate the clusters according to the following condition:
							1) if DM_ID[i][j] is strictly lower then use the new cluster pair
							2) if DM_ID[i][j] = DM_ID[min_idx_i][min_idx_j], 
							and DM_XL[i][j] is strictly lower then use the new cluster pair
							3) if both internode count and XL measures are identical for the above mentioned clusters, 
							agglomerate the new cluster pair, if they have higher sum of cardinality
							"""
							if (DM_ID[i][j] < DM_ID[min_idx_i][min_idx_j]):
								min_idx_i = i
								min_idx_j = j
							else:
								if (FlEq(DM_ID[i][j], DM_ID[min_idx_i][min_idx_j]) == True):
									if (DM_XL[i][j] < DM_XL[min_idx_i][min_idx_j]):
										min_idx_i = i
										min_idx_j = j
									else:
										if (FlEq(DM_XL[i][j], DM_XL[min_idx_i][min_idx_j]) == True):
											# higher sum of cardinality of the new cluster pair
											if (len(clust_species_list[i]) + len(clust_species_list[j])) \
												> (len(clust_species_list[min_idx_i]) + len(clust_species_list[min_idx_j])):
												min_idx_i = i
												min_idx_j = j

	return min_idx_i, min_idx_j

#---------------------------------------------
"""
this is the third version of the PNJSTXL / IDXL method
here, the positional rank of XL based distance matrix is used
also, zero mean unit variance of the relative ID and XL distance matrices are employed
"""
def Find_Unique_Min_PNJSTXL_Version3(DM_ID, Norm_DM_ID, DM_XL, Norm_DM_XL, nclust, TaxaList, outfile):

	if (DEBUG_LEVEL >= 2):
		fp = open(outfile, 'a')
		fp.write('\n\n\n ==>>>> printing contents of Correlation between NDM_ID and NDM_XL --->>> \n\n')

	"""
	at first, use the DM_XL (containing the absolute values of excess gene leaf count for individual couplets)
	to create a sorted array
	"""
	Sorted_DM_XL = []

	"""
	create two arrays containing the relative ID and XL based distance
	"""
	List_Norm_DM_ID = []
	List_Norm_DM_XL = []
	for i in range(nclust - 1):
		for j in range(i+1, nclust):
			if (i == j):
				continue
			if (DM_XL[i][j] >= 0) and (DM_ID[i][j] >= 0):	# check for valid entry
				List_Norm_DM_ID.append(Norm_DM_ID[i][j])
				List_Norm_DM_XL.append(Norm_DM_XL[i][j])
				Sorted_DM_XL.append(DM_XL[i][j])
	
	"""
	find the unique elements of array
	we have found that using set based unique operation is highly efficient
	"""
	Sorted_DM_XL = list(set(Sorted_DM_XL))
	
	# sort the list in ascending order
	Sorted_DM_XL.sort()
	# no of elements in the sorted array
	len_Sorted_DM_XL = len(Sorted_DM_XL)
	
	"""
	compute the mean and standard deviations of these two arrays
	"""
	mean_NDM_ID = Compute_Mean(List_Norm_DM_ID)
	mean_NDM_XL = Compute_Mean(List_Norm_DM_XL)
	stdev_NDM_ID = Pop_StDev(List_Norm_DM_ID)
	stdev_NDM_XL = Pop_StDev(List_Norm_DM_XL)
	
	if (DEBUG_LEVEL >= 2): 
		fp.write('\n mean_NDM_ID: ' + str(mean_NDM_ID) + ' mean_NDM_XL: ' + str(mean_NDM_XL) + \
			'  stdev_NDM_ID: ' + str(stdev_NDM_ID) + '  stdev_NDM_XL: ' + str(stdev_NDM_XL))
	
	"""
	this is a flag variable which is set when the first valid distance element is found
	"""
	tgt_found = 0
	
	for i in range(nclust):
		if (DEBUG_LEVEL >= 2): 
			fp.write('\n ' + str(i) + '--' + str(TaxaList[i]) + '--->>')
		if (i > 0):
			for j in range(i):
				if ((DM_ID[i][j] >= 0) and (DM_XL[i][j] >= 0)):	# corresponding to a valid entry
					if (FlEq(stdev_NDM_ID, 0) == True) or (FlEq(stdev_NDM_XL, 0) == True):
						"""
						here, standard deviation values of one or both of the relative 
						distance matrices are 0, 
						so, only the product of relative internode count 
						and excess gene count measures are employed
						"""
						corr_val = Norm_DM_ID[i][j] * Norm_DM_XL[i][j]
						if (DEBUG_LEVEL >= 2): 
							fp.write(' ' + str(corr_val))	#print the correlation value
						"""
						check the positional rank of the XL based distance matrix entry (absolute)
						"""
						curr_DM_XL = DM_XL[i][j]
						rank_curr_DM_XL = Sorted_DM_XL.index(DM_XL[i][j])
						if (rank_curr_DM_XL <= int(round(len_Sorted_DM_XL / 2))):	#sourya
							if (tgt_found == 0):
								tgt_found = 1
								min_DM_XL = curr_DM_XL
								max_idx_i = i
								max_idx_j = j
							else:
								if (curr_DM_XL < min_DM_XL):
									min_DM_XL = curr_DM_XL
									max_idx_i = i
									max_idx_j = j
								elif (FlEq(curr_DM_XL, min_DM_XL) == True):
									"""
									for equality case, check the ID (absolute) values
									"""
									if (DM_ID[i][j] < DM_ID[max_idx_i][max_idx_j]):
										max_idx_i = i
										max_idx_j = j
					
					else:	#if (stdev_NDM_ID > 0) and (stdev_NDM_XL > 0):	# sourya
						corr_val = ((mean_NDM_ID - Norm_DM_ID[i][j]) / stdev_NDM_ID) + ((mean_NDM_XL - Norm_DM_XL[i][j]) / stdev_NDM_XL)
						if (DEBUG_LEVEL >= 2): 
							fp.write(' ' + str(corr_val))	#print the correlation value
						"""
						as both stdev values are > 0, processing would be done by the correlation values
						find the index (rank) of the current XL val 
						process the element only if the XL value is comparatively low
						"""
						rank_curr_DM_XL = Sorted_DM_XL.index(DM_XL[i][j])
						if (rank_curr_DM_XL <= int(round(len_Sorted_DM_XL / 2))):	#sourya
							if (tgt_found == 0):
								tgt_found = 1
								max_corr_val = corr_val
								max_idx_i = i
								max_idx_j = j
							else:
								if (corr_val > max_corr_val):
									max_corr_val = corr_val
									max_idx_i = i
									max_idx_j = j
								elif (FlEq(corr_val, max_corr_val) == True):
									"""
									for equality case, first check whether the relative internode count is lower
									"""
									if (Norm_DM_ID[i][j] < Norm_DM_ID[max_idx_i][max_idx_j]):
										max_idx_i = i
										max_idx_j = j
									else:
										if (Norm_DM_XL[i][j] < Norm_DM_XL[max_idx_i][max_idx_j]):
											max_idx_i = i
											max_idx_j = j
				else:
					corr_val = -2
					if (DEBUG_LEVEL >= 2): 
						fp.write(' ' + str(corr_val))
	
	if (DEBUG_LEVEL >= 2): 
		if (FlEq(stdev_NDM_ID, 0) == True) or (FlEq(stdev_NDM_XL, 0) == True):
			fp.write('\n Final max correlation index i: ' + str(max_idx_i) + ' j: ' + str(max_idx_j) + \
				'  Min XL: ' + str(min_DM_XL))
		else:	#if (stdev_NDM_ID > 0) and (stdev_NDM_XL > 0):	# sourya
			fp.write('\n Final max correlation index i: ' + str(max_idx_i) + ' j: ' + str(max_idx_j) + \
				'  Max correlation: ' + str(max_corr_val))
			
	if (DEBUG_LEVEL >= 2):
		fp.close()
	
	"""
	reset the lists employed
	"""
	Sorted_DM_XL = []
	List_Norm_DM_ID = []
	List_Norm_DM_XL = []
	
	return max_idx_j, max_idx_i

#---------------------------------------------
"""
this is the fourth version of the PNJSTXL / IDXL method
here, zero mean unit variance of the relative ID and XL distance matrices, and also 
the absolute distance matrices, are employed
"""
def Find_Unique_Min_PNJSTXL_Version4(DM_ID, Norm_DM_ID, DM_XL, Norm_DM_XL, nclust, TaxaList, outfile):

	if (DEBUG_LEVEL >= 2):
		fp = open(outfile, 'a')
		fp.write('\n\n\n ==>>>> printing contents of Correlation between NDM_ID and NDM_XL --->>> \n\n')

	"""
	create two arrays containing the relative ID and XL based distance
	"""
	List_Norm_DM_ID = []
	List_Norm_DM_XL = []
	List_DM_XL = []
	for i in range(nclust - 1):
		for j in range(i+1, nclust):
			if (i == j):
				continue
			if (DM_XL[i][j] >= 0) and (DM_ID[i][j] >= 0):	# check for valid entry
				List_Norm_DM_ID.append(Norm_DM_ID[i][j])
				List_Norm_DM_XL.append(Norm_DM_XL[i][j])
				List_DM_XL.append(DM_XL[i][j])
				
	"""
	compute the mean and standard deviations of these two arrays
	"""
	mean_NDM_ID = Compute_Mean(List_Norm_DM_ID)
	mean_NDM_XL = Compute_Mean(List_Norm_DM_XL)
	mean_DM_XL = Compute_Mean(List_DM_XL)
	stdev_NDM_ID = Pop_StDev(List_Norm_DM_ID)
	stdev_NDM_XL = Pop_StDev(List_Norm_DM_XL)
	stdev_DM_XL = Pop_StDev(List_DM_XL)
	
	if (DEBUG_LEVEL >= 2): 
		fp.write('\n mean_NDM_ID: ' + str(mean_NDM_ID) + ' mean_NDM_XL: ' + str(mean_NDM_XL) + ' mean_DM_XL: ' + str(mean_DM_XL) + \
			'  stdev_NDM_ID: ' + str(stdev_NDM_ID) + '  stdev_NDM_XL: ' + str(stdev_NDM_XL) + ' stdev_DM_XL: ' + str(stdev_DM_XL))
	
	"""
	this is a flag variable which is set when the first valid distance element is found
	"""
	tgt_found = 0
	
	for i in range(nclust):
		if (DEBUG_LEVEL >= 2):
			fp.write('\n ' + str(i) + '--' + str(TaxaList[i]) + '--->>')
		if (i > 0):
			for j in range(i):
				if ((DM_ID[i][j] >= 0) and (DM_XL[i][j] >= 0)):	# condition add - sourya
					if (FlEq(stdev_NDM_ID, 0) == True) or (FlEq(stdev_DM_XL, 0) == True) or (FlEq(stdev_NDM_XL, 0) == True):
						"""
						here, standard deviation values of one or both of the relative 
						distance matrices are 0, 
						so, only the product of relative internode count 
						and excess gene count measures are employed
						"""
						corr_val = Norm_DM_ID[i][j] * Norm_DM_XL[i][j]
						if (DEBUG_LEVEL >= 2):
							fp.write(' ' + str(corr_val))	#print the correlation value
						curr_DM_XL = DM_XL[i][j]
						if (tgt_found == 0):
							tgt_found = 1
							min_DM_XL = curr_DM_XL
							max_idx_i = i
							max_idx_j = j
						else:
							if (curr_DM_XL < min_DM_XL):
								min_DM_XL = curr_DM_XL
								max_idx_i = i
								max_idx_j = j
							elif (FlEq(curr_DM_XL, min_DM_XL) == True):
								"""
								for equality case, check the ID (absolute) values
								"""
								if (DM_ID[i][j] < DM_ID[max_idx_i][max_idx_j]):
									max_idx_i = i
									max_idx_j = j
					else:	#if (stdev_NDM_ID > 0) and (stdev_NDM_XL > 0):	# sourya
						"""
						as both stdev values are > 0, processing would be done by the correlation values 
						computed using both the relative XL and ID distance values
						and the absolute XL values
						"""
						corr_val = ((mean_NDM_ID - Norm_DM_ID[i][j]) / stdev_NDM_ID) \
							+ ((mean_NDM_XL - Norm_DM_XL[i][j]) / stdev_NDM_XL) + ((mean_DM_XL - DM_XL[i][j]) / stdev_DM_XL)
						if (DEBUG_LEVEL >= 2):
							fp.write(' ' + str(corr_val))	#print the correlation value
						if (tgt_found == 0):
							tgt_found = 1
							max_corr_val = corr_val
							max_idx_i = i
							max_idx_j = j
						else:
							if (corr_val > max_corr_val):
								max_corr_val = corr_val
								max_idx_i = i
								max_idx_j = j
							elif (FlEq(corr_val, max_corr_val) == True):
								"""
								for equality case, first check whether the relative internode count is lower
								"""
								if (Norm_DM_ID[i][j] < Norm_DM_ID[max_idx_i][max_idx_j]):
									max_idx_i = i
									max_idx_j = j
								else:
									if (Norm_DM_XL[i][j] < Norm_DM_XL[max_idx_i][max_idx_j]):
										max_idx_i = i
										max_idx_j = j
				else:
					corr_val = -2
					if (DEBUG_LEVEL >= 2):
						fp.write(' ' + str(corr_val))
	
	if (DEBUG_LEVEL >= 2): 
		if (FlEq(stdev_NDM_ID, 0) == True) or (FlEq(stdev_DM_XL, 0) == True) or (FlEq(stdev_NDM_XL, 0) == True):
			fp.write('\n Final max correlation index i: ' + str(max_idx_i) + ' j: ' + str(max_idx_j) + \
				'  Min XL: ' + str(min_DM_XL))
		else:	#if (stdev_NDM_ID > 0) and (stdev_NDM_XL > 0):	# sourya
			fp.write('\n Final max correlation index i: ' + str(max_idx_i) + ' j: ' + str(max_idx_j) + \
				'  Max correlation: ' + str(max_corr_val))
			
	if (DEBUG_LEVEL >= 2):
		fp.close()
	
	"""
	reset the lists employed
	"""
	List_Norm_DM_ID = []
	List_Norm_DM_XL = []
	List_DM_XL = []
	
	return max_idx_j, max_idx_i

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
this function refines the initial species tree (in terms of a star network) to find the species tree
it does using agglomerative clustering (NJ principle)
the distance metric employed for NJ algorithm 
consists of the couplet based internode count and excess gene leaf count
"""
def Form_Species_Tree_NJ_Cluster(Curr_tree, Output_Text_File, METHOD_USED, DIST_MAT_TYPE):

	"""
	initially we have N clusters for N taxa, where individual clusters are isolated taxon
	agglomerating technique introduces a speciation node which contains two taxa as its children
	"""
	no_of_taxa_clust = len(COMPLETE_INPUT_TAXA_LIST)

	"""
	initialize the taxa clusters
	copying the taxa list is done since initial clusters contain single species  
	"""
	clust_species_list = []
	for i in range(len(COMPLETE_INPUT_TAXA_LIST)):
		subl = []
		subl.append(COMPLETE_INPUT_TAXA_LIST[i])
		clust_species_list.append(subl)

	"""
	allocate a 2D square matrix of dimension N X N
	for a pair of taxa clusters Cx and Cy, 
	it contains the couplet based internode count distance values
	"""
	Mean_DistMat_ClustPair_NJ = numpy.zeros((no_of_taxa_clust, no_of_taxa_clust), dtype=numpy.float)

	"""
	allocate one new square matrix which will contain the 
	NJ based relative internode count for individual couplets
	"""
	Norm_Mean_DistMat_ClustPair_NJ = numpy.zeros((no_of_taxa_clust, no_of_taxa_clust), dtype=numpy.float)

	"""
	Similarly allocate N X N matrix to store the XL measure for individual couplets
	"""
	XLVal_DistMat_ClustPair_NJ = numpy.zeros((no_of_taxa_clust, no_of_taxa_clust), dtype=numpy.float)
	
	"""
	allocate one new square matrix which will contain the relative XL distance during NJ based iterations
	"""
	Norm_XLVal_DistMat_ClustPair_NJ = numpy.zeros((no_of_taxa_clust, no_of_taxa_clust), dtype=numpy.float)

	"""
	fill input distance matrix using internode count measure
	"""
	Fill_DistMat_InternodeCount(Mean_DistMat_ClustPair_NJ, no_of_taxa_clust)

	"""
	fill input distance matrix with Xl measure
	"""
	Fill_DistMat_ExcessGeneCount(XLVal_DistMat_ClustPair_NJ, DIST_MAT_TYPE, no_of_taxa_clust)

	#--------------------------------------------------------
	"""
	loop to execute the agglomerative clustering
	"""
	while(no_of_taxa_clust > 2): 
		if (DEBUG_LEVEL >= 2):
			fp = open(Output_Text_File, 'a')
			fp.write('\n\n\n *********************************** \n iteration start --- number of clusters: ' + str(no_of_taxa_clust))
			fp.write('\n clust_species_list : ' + str(clust_species_list))
			fp.close()
			PrintMatrixContent(no_of_taxa_clust, clust_species_list, Mean_DistMat_ClustPair_NJ, \
				'Mean_DistMat_ClustPair_NJ', Output_Text_File)
			PrintMatrixContent(no_of_taxa_clust, clust_species_list, XLVal_DistMat_ClustPair_NJ, \
				'XLVal_DistMat_ClustPair_NJ', Output_Text_File)

		"""
		for individual cluster Cx, it contains XL(Cx, :) - sum of extra lineages considering the cluster pair 
		(Cx, Cy) for all other clusters Cy
		"""
		sum_DistMat_Clust = []
		ComputeSumRowsDistMat(sum_DistMat_Clust, no_of_taxa_clust, \
			Mean_DistMat_ClustPair_NJ, Output_Text_File, 'sum_DistMat_Clust')
		sum_XLVal_Clust = []
		ComputeSumRowsDistMat(sum_XLVal_Clust, no_of_taxa_clust, \
			XLVal_DistMat_ClustPair_NJ, Output_Text_File, 'sum_XLVal_Clust')
		
		"""
		fill the relative distance matrix entries
		"""
		FillNJNormalizeMatrix(Norm_Mean_DistMat_ClustPair_NJ, \
			Mean_DistMat_ClustPair_NJ, sum_DistMat_Clust, no_of_taxa_clust)
		FillNJNormalizeMatrix(Norm_XLVal_DistMat_ClustPair_NJ, \
			XLVal_DistMat_ClustPair_NJ, sum_XLVal_Clust, no_of_taxa_clust)

		if (DEBUG_LEVEL >= 2):
			PrintMatrixContent(no_of_taxa_clust, clust_species_list, Norm_Mean_DistMat_ClustPair_NJ, \
				'Norm_Mean_DistMat_ClustPair_NJ', Output_Text_File)
			if 1:	#(METHOD_USED == NJSTXL) or (METHOD_USED == MedNJSTXL) or (METHOD_USED == ModeNJSTXL):
				PrintMatrixContent(no_of_taxa_clust, clust_species_list, Norm_XLVal_DistMat_ClustPair_NJ, \
					'Norm_XLVal_DistMat_ClustPair_NJ', Output_Text_File)
			
		"""
		find the cluster pairs having minimum distance values
		the distance measure consists of internode count and excess gene leaf count
		"""
		if (METHOD_USED == 1):
			min_idx_i, min_idx_j = Find_Unique_Min_PNJSTXL_Version1(Mean_DistMat_ClustPair_NJ, Norm_Mean_DistMat_ClustPair_NJ, \
				XLVal_DistMat_ClustPair_NJ, Norm_XLVal_DistMat_ClustPair_NJ, no_of_taxa_clust, clust_species_list)
		elif (METHOD_USED == 2):
			min_idx_i, min_idx_j = Find_Unique_Min_PNJSTXL_Version2(Mean_DistMat_ClustPair_NJ, Norm_Mean_DistMat_ClustPair_NJ, \
				XLVal_DistMat_ClustPair_NJ, Norm_XLVal_DistMat_ClustPair_NJ, no_of_taxa_clust, clust_species_list)
		elif (METHOD_USED == 3):
			min_idx_i, min_idx_j = Find_Unique_Min_PNJSTXL_Version3(Mean_DistMat_ClustPair_NJ, Norm_Mean_DistMat_ClustPair_NJ, \
				XLVal_DistMat_ClustPair_NJ, Norm_XLVal_DistMat_ClustPair_NJ, no_of_taxa_clust, clust_species_list, Output_Text_File)
		else:
			min_idx_i, min_idx_j = Find_Unique_Min_PNJSTXL_Version4(Mean_DistMat_ClustPair_NJ, Norm_Mean_DistMat_ClustPair_NJ, \
				XLVal_DistMat_ClustPair_NJ, Norm_XLVal_DistMat_ClustPair_NJ, no_of_taxa_clust, clust_species_list, Output_Text_File)
		#---------------------------------------------------------------      
		"""
		note down the taxa list in these two indices (min_idx_i and min_idx_j) 
		of the clust_species_list
		"""
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
		"""
		adjust the Mean_DistMat_ClustPair_NJ by inserting one new row and column corresponding to the new cluster
		and then deleting the information of earlier two clusters
		"""
		"""
		first append one row
		"""
		Mean_DistMat_ClustPair_NJ = numpy.vstack((Mean_DistMat_ClustPair_NJ, \
			numpy.zeros((1, no_of_taxa_clust), dtype=numpy.float)))
		"""
		then append one column
		"""
		Mean_DistMat_ClustPair_NJ = numpy.hstack((Mean_DistMat_ClustPair_NJ, \
			numpy.zeros((no_of_taxa_clust + 1, 1), dtype=numpy.float)))
		"""
		now reshape the distance matrix
		"""
		Mean_DistMat_ClustPair_NJ = numpy.reshape(Mean_DistMat_ClustPair_NJ, \
			((no_of_taxa_clust + 1), (no_of_taxa_clust + 1)), order='C')
		
		"""
		apply these operations on the excess gene count matrix as well
		"""
		XLVal_DistMat_ClustPair_NJ = numpy.vstack((XLVal_DistMat_ClustPair_NJ, \
			numpy.zeros((1, no_of_taxa_clust), dtype=numpy.float)))
		XLVal_DistMat_ClustPair_NJ = numpy.hstack((XLVal_DistMat_ClustPair_NJ, \
			numpy.zeros((no_of_taxa_clust + 1, 1), dtype=numpy.float)))
		XLVal_DistMat_ClustPair_NJ = numpy.reshape(XLVal_DistMat_ClustPair_NJ, \
			((no_of_taxa_clust + 1), (no_of_taxa_clust + 1)), order='C')
		
		"""
		add taxa_list as a new element of clust_species_list
		"""
		clust_species_list.append(taxa_list)          
		
		"""
		now recompute the entries of this new row and column (which is indexed by no_of_taxa_clust), according to the NJ principle
		compute Mean_DistMat_ClustPair_NJ[no_of_taxa_clust][m] entries 
		where m != min_idx_i and m != min_idx_j
		"""
		for m in range(no_of_taxa_clust):
			if (m == min_idx_i) or (m == min_idx_j):
				continue
			
			"""
			update the distance matrix storing internode count
			maintain the symmetric property of the distance matrix
			"""
			if (Mean_DistMat_ClustPair_NJ[min_idx_i][m] >= 0) and (Mean_DistMat_ClustPair_NJ[min_idx_j][m] < 0):
				Mean_DistMat_ClustPair_NJ[no_of_taxa_clust][m] = \
					Mean_DistMat_ClustPair_NJ[min_idx_i][m] - (Mean_DistMat_ClustPair_NJ[min_idx_i][min_idx_j] / 2.0)
			elif (Mean_DistMat_ClustPair_NJ[min_idx_i][m] < 0) and (Mean_DistMat_ClustPair_NJ[min_idx_j][m] >= 0):
				Mean_DistMat_ClustPair_NJ[no_of_taxa_clust][m] = \
					Mean_DistMat_ClustPair_NJ[min_idx_j][m] - (Mean_DistMat_ClustPair_NJ[min_idx_i][min_idx_j] / 2.0)
			else:
				Mean_DistMat_ClustPair_NJ[no_of_taxa_clust][m] = (Mean_DistMat_ClustPair_NJ[min_idx_i][m] + \
					Mean_DistMat_ClustPair_NJ[min_idx_j][m] - Mean_DistMat_ClustPair_NJ[min_idx_i][min_idx_j]) / 2.0
			
			# maintain symmetric property
			Mean_DistMat_ClustPair_NJ[m][no_of_taxa_clust] = Mean_DistMat_ClustPair_NJ[no_of_taxa_clust][m]

			"""
			update the distance matrix storing excess gene leaf count
			maintain the symmetric property of the distance matrix
			"""
			if (XLVal_DistMat_ClustPair_NJ[min_idx_i][m] >= 0) and (XLVal_DistMat_ClustPair_NJ[min_idx_j][m] < 0):
				XLVal_DistMat_ClustPair_NJ[no_of_taxa_clust][m] = XLVal_DistMat_ClustPair_NJ[min_idx_i][m]
			elif (XLVal_DistMat_ClustPair_NJ[min_idx_i][m] < 0) and (XLVal_DistMat_ClustPair_NJ[min_idx_j][m] >= 0):
				XLVal_DistMat_ClustPair_NJ[no_of_taxa_clust][m] = XLVal_DistMat_ClustPair_NJ[min_idx_j][m]
			else:
				"""
				simple averaging of XL values are used
				"""
				XLVal_DistMat_ClustPair_NJ[no_of_taxa_clust][m] = \
					(XLVal_DistMat_ClustPair_NJ[min_idx_i][m] + XLVal_DistMat_ClustPair_NJ[min_idx_j][m]) / 2.0
			
			# maintain symmetric property
			XLVal_DistMat_ClustPair_NJ[m][no_of_taxa_clust] = XLVal_DistMat_ClustPair_NJ[no_of_taxa_clust][m]

		"""
		now remove the rows and columns corresponding to min_idx_i and min_idx_j
		"""
		Mean_DistMat_ClustPair_NJ = numpy.delete(Mean_DistMat_ClustPair_NJ, (min_idx_i), axis=0)	# delete the row
		Mean_DistMat_ClustPair_NJ = numpy.delete(Mean_DistMat_ClustPair_NJ, (min_idx_i), axis=1)	# delete the column
		Mean_DistMat_ClustPair_NJ = numpy.delete(Mean_DistMat_ClustPair_NJ, (min_idx_j - 1), axis=0)	# delete the row
		Mean_DistMat_ClustPair_NJ = numpy.delete(Mean_DistMat_ClustPair_NJ, (min_idx_j - 1), axis=1)	# delete the column

		"""
		update the excess gene count matrix as well
		"""
		XLVal_DistMat_ClustPair_NJ = numpy.delete(XLVal_DistMat_ClustPair_NJ, (min_idx_i), axis=0)	# delete the row
		XLVal_DistMat_ClustPair_NJ = numpy.delete(XLVal_DistMat_ClustPair_NJ, (min_idx_i), axis=1)	# delete the column
		XLVal_DistMat_ClustPair_NJ = numpy.delete(XLVal_DistMat_ClustPair_NJ, (min_idx_j - 1), axis=0)	# delete the row
		XLVal_DistMat_ClustPair_NJ = numpy.delete(XLVal_DistMat_ClustPair_NJ, (min_idx_j - 1), axis=1)	# delete the column

		"""
		clear Norm_Mean_DistMat_ClustPair_NJ
		"""
		Norm_Mean_DistMat_ClustPair_NJ = numpy.delete(Norm_Mean_DistMat_ClustPair_NJ, (min_idx_i), axis=0)	# delete the row
		Norm_Mean_DistMat_ClustPair_NJ = numpy.delete(Norm_Mean_DistMat_ClustPair_NJ, (min_idx_i), axis=1)	# delete the column
		Norm_Mean_DistMat_ClustPair_NJ.fill(0)
		
		"""
		clear norm excess gene count matrix
		"""
		Norm_XLVal_DistMat_ClustPair_NJ = numpy.delete(Norm_XLVal_DistMat_ClustPair_NJ, (min_idx_i), axis=0)	# delete the row
		Norm_XLVal_DistMat_ClustPair_NJ = numpy.delete(Norm_XLVal_DistMat_ClustPair_NJ, (min_idx_i), axis=1)	# delete the column
		Norm_XLVal_DistMat_ClustPair_NJ.fill(0)    
		
		"""
		remove individual clusters' taxa information from the clust_species_list
		"""
		clust_species_list.pop(min_idx_i)
		clust_species_list.pop(min_idx_j - 1)
		
		# decrement the number of clusters considered
		no_of_taxa_clust = no_of_taxa_clust - 1

	return

