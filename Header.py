#!/usr/bin/env python

"""
program for species tree estimation
computes four criteria - 
1) Level / internode count information
2) Accumulated Rank statistics 
"""

import dendropy
from optparse import OptionParser
#import math
import time
import os
import numpy
import sys
#import matplotlib.pyplot as plt	# add - debug - sourya

""" 
this dictionary defines the taxa pair relations
each entry is indexed by two nodes 
"""
TaxaPair_Reln_Dict = dict()

# these variables store the respective parameters of input gene trees
COMPLETE_INPUT_TAXA_LIST = []

# this is the debug level
# set for printing the necessary information
DEBUG_LEVEL = 0

"""
variables depicting the method employed for species tree construction
first two methods use average statistics of the coalescence rank or the internode count
last two methods employ the mode statistics
"""
"""
accumulated internode count with average statistics
plus average of excess gene count measure 
"""
NJSTXL = 3
"""
accumulated internode count with average statistics
plus minimum (median, average) of excess gene count measure 
"""
MedNJSTXL = 1

"""
product of average internode count and excess lineage information
to infer species tree
"""
ProdNJSTXL = 2

# variables used to denote whether we use traditional NJ method
# or use a variant of it, namely the agglomerative clustering
TRADITIONAL_NJ = 1
AGGLO_CLUST = 2

# add - sourya
MODE_PERCENT = 0.5

MODE_BIN_COUNT = 40

# this list corresponds to various rank merging algorithm
SIMPLE_SUM_RANK = 1
MEAN_RECIPROCAL_RANK = 2

##-----------------------------------------------------
""" 
this class defines the connectivity relationship between a pair of taxa
initially the information are obtained from the input source trees
"""
class Reln_TaxaPair(object):  
	def __init__(self):
		# this is the count of trees for which the couplet is supported
		self.tree_support_count = 0        
		# this list contains the number of levels (tree internodees) between individual couplets
		# computed for all the gene trees
		self.sum_internode_count = 0
		# this is the extra lineage count list for this couplet
		self.XL_val_list = []
		"""
		this is a variable containing the binned average of the XL values
		of very high frequency
		initially the value is set as -1, to signify that the computation is not done
		once the computation (for a couplet) is done, the value is subsequently used and returned
		"""
		self.binned_avg_XL = -1
		
		self.avg_XL = -1
		self.median_XL = -1
		
	# this function adds the count of tree according to the support of 
	# corresponding couplet in the input tree
	def _IncrSupportTreeCount(self):
		self.tree_support_count = self.tree_support_count + 1
	
	# this function returns the number of trees supporting the couplet
	def _GetSupportTreeCount(self):
		return self.tree_support_count        
				
	def _AddXLVal(self, val):
		self.XL_val_list.append(val)

	def _GetAvgXLVal(self):
		if (self.avg_XL == -1):
			self.avg_XL = (sum(self.XL_val_list) * 1.0) / self.tree_support_count
		return self.avg_XL
	
	def _MedianXLVal(self):
		if (self.median_XL == -1):
			self.median_XL = numpy.median(numpy.array(self.XL_val_list))
		return self.median_XL
	
	#------------------------------------------
	#def _GetMultiModeXLVal(self):
		#candidate_score_sum = 0
		#candidate_freq_sum = 0
		#curr_arr = numpy.array(self.XL_val_list)
		## returns the counts of individual elements
		## array size: max_elem + 1
		#counts = numpy.bincount(curr_arr)
		## remove the zero values 
		#values = numpy.nonzero(counts)[0]
		## mode value and corresponding frequency
		#mode_val = numpy.argmax(counts)
		#mode_count = numpy.max(counts)
		## check for the values having frequency at least half of the maximum frequency
		#for v in values:
			#if (counts[v] >= 0.5 * mode_count):
				#candidate_score_sum = candidate_score_sum + (v * counts[v])
				#candidate_freq_sum = candidate_freq_sum + counts[v]
		#return (candidate_score_sum * 1.0) / candidate_freq_sum        

	#------------------------------------------
	#def _GetMultiModeXLVal_New(self):
		#candidate_score_sum = 0
		#candidate_freq_sum = 0
		#curr_arr = numpy.array(self.XL_val_list)
		#uniqw, inverse = numpy.unique(curr_arr, return_inverse=True)
		#counts = numpy.bincount(inverse)
		## remove the zero values 
		#values = numpy.nonzero(counts)[0]
		## mode value and corresponding frequency
		#mode_count = numpy.max(counts)
		
		##print '*****************'
		##print 'curr_arr: ', curr_arr
		##print 'uniqw: ', uniqw
		##print 'inverse: ', inverse
		##print 'counts: ', counts
		##print 'values: ', values
		##print 'mode_count: ', mode_count
		
		## check for the values having frequency at least half of the maximum frequency
		#for v in values:
			#if (counts[v] >= MODE_PERCENT * mode_count):
				##candidate_score_sum = candidate_score_sum + (v * counts[v])	# comment - sourya
				#candidate_score_sum = candidate_score_sum + (uniqw[v] * counts[v])	# add - sourya
				#candidate_freq_sum = candidate_freq_sum + counts[v]
		#return (candidate_score_sum * 1.0) / candidate_freq_sum        

	#------------------------------------------
	def _GetMultiModeXLVal(self, Output_Text_File=None):
		if (self.binned_avg_XL == -1):
			
			Bin_Width = (1.0 / MODE_BIN_COUNT)
			len_list = [0] * MODE_BIN_COUNT
			
			if Output_Text_File is not None:
				fp = open(Output_Text_File, 'a') 
			
			# sort the XL list
			self.XL_val_list.sort()
			
			for j in range(len(self.XL_val_list)):
				curr_xl_val = self.XL_val_list[j]
				bin_idx = int(curr_xl_val / Bin_Width)
				if (bin_idx == MODE_BIN_COUNT):
					bin_idx = bin_idx - 1
				len_list[bin_idx] = len_list[bin_idx] + 1
			
			if Output_Text_File is not None:
				for i in range(MODE_BIN_COUNT):
					fp.write('\n bin idx: ' + str(i) + ' len:  ' + str(len_list[i]))
			
			# this is the maximum length of a particular bin
			# corresponding to max frequency
			max_freq = max(len_list)
			
			if Output_Text_File is not None:
				fp.write('\n Max freq: ' + str(max_freq))
			
			num = 0
			denom = 0
			for i in range(MODE_BIN_COUNT):
				if (len_list[i] >= (MODE_PERCENT * max_freq)):
					list_start_idx = sum(len_list[:i])
					list_end_idx = list_start_idx + len_list[i] - 1
					value_sum = sum(self.XL_val_list[list_start_idx:(list_end_idx+1)])
					num = num + value_sum
					denom = denom + len_list[i]
					if Output_Text_File is not None:
						fp.write('\n Included bin idx: ' + str(i) + ' starting point: ' + str(list_start_idx) \
							+ 'ending point: ' + str(list_end_idx) + ' sum: ' + str(value_sum))
			
			self.binned_avg_XL = (num / denom)
			
			if Output_Text_File is not None:
				fp.write('\n Final binned average XL: ' + str(self.binned_avg_XL))
				fp.close()
			
		return self.binned_avg_XL
	#------------------------------------------

	def _AddLevel(self, val):
		self.sum_internode_count = self.sum_internode_count + val

	def _GetAvgSumLevel(self):
		return (self.sum_internode_count * 1.0) / self.tree_support_count
					
	#def _GetMultiModeSumLevel(self):
		#candidate_score_sum = 0
		#candidate_freq_sum = 0
		#curr_arr = numpy.array(self.sum_internode_count)
		## returns the counts of individual elements
		## array size: max_elem + 1
		#counts = numpy.bincount(curr_arr)
		## remove the zero values 
		#values = numpy.nonzero(counts)[0]
		## mode value and corresponding frequency
		#mode_val = numpy.argmax(counts)
		#mode_count = numpy.max(counts)
		## check for the values having frequency at least half of the maximum frequency
		#for v in values:
			#if (counts[v] >= 0.5 * mode_count):
				#candidate_score_sum = candidate_score_sum + (v * counts[v])
				#candidate_freq_sum = candidate_freq_sum + counts[v]
		#return (candidate_score_sum * 1.0) / candidate_freq_sum
																		
	# this function prints information for the current couplet
	def _PrintTaxaPairRelnInfo(self, key, out_text_file, METHOD_USED):
		fp = open(out_text_file, 'a')    
		fp.write('\n taxa pair key: ' + str(key))
		fp.write('\n supporting number of trees: ' + str(self._GetSupportTreeCount()))
		fp.write('\n *** average sum of internode count : ' + str(self._GetAvgSumLevel()))    
		fp.write('\n *** average XL val : ' + str(self._GetAvgXLVal()))   
		fp.write('\n *** median XL val : ' + str(self._MedianXLVal()))   
		fp.write('\n *** 50 percent mode XL val : ' + str(self._GetMultiModeXLVal()))   
		fp.close()
			
		## sourya - debug
		#if (key[0] == 'HOM' and key[1] == 'TAR') or (key[0] == 'MYO' and key[1] == 'TUR'):
			#print 'printing the file'
			#print 'current directory: ', os.getcwd()
			
			##fig1 = plt.figure()
			##n1, bins1, patches1 = plt.hist(self.sum_internode_count, 37, normed=0)	#, facecolor='green', alpha=0.75)
			##xlabel_str = 'I_G(' + str(key[0]) + ',' + str(key[1]) + ')'
			##plt.xlabel(xlabel_str)
			##plt.ylabel('Frequency')
			##plt.title('Histogram of the internode count across gene trees for the couplet ' + str(key[0]) + ' and ' + str(key[1]))
			##plt.grid(True)
			##plt.tight_layout()
			##fig1.set_size_inches(10, 6)
			##figname = 'internode_count_' + str(key[0]) + '_' + str(key[1]) + '.jpg'
			##print 'figname: ', figname
			##plt.savefig(figname)
			##plt.clf()	# clear the current figure memory
			##plt.close()	# close the figure window
			

			#fig2 = plt.figure()
			#n2, bins2, patches2 = plt.hist(self.XL_val_list, 37, normed=0)	#, facecolor='green', alpha=0.75)
			#xlabel_str = 'X_G(' + str(key[0]) + ',' + str(key[1]) + ')'
			#plt.xlabel(xlabel_str)
			#plt.ylabel('Frequency')
			#plt.title('Histogram of the extra leaf count across gene trees for the couplet ' + str(key[0]) + ' and ' + str(key[1]))
			#plt.grid(True)
			#plt.tight_layout()
			#fig2.set_size_inches(10, 6)
			#figname = 'extra_leaf_' + str(key[0]) + '_' + str(key[1]) + '.jpg'
			#print 'figname: ', figname
			#plt.savefig(figname)      
			#plt.clf()	# clear the current figure memory
			#plt.close()	# close the figure window
			
		## end sourya - debug
			
		
    
      