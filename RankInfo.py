import Header
from Header import *  
import UtilFunc
from UtilFunc import *

#--------------------------------------------------------
"""
this function adds the following two measures for individual couplets:
1) Internode count 
2) Excess gene leaf count
Both measures use the LCA node for the corresponding couplet
"""
def Compute_Internode_ExcessGeneLeaf(taxa_count, xl_val, curr_node_level, node1_level, node2_level, node1_idx, node2_idx):

	"""
	sourya - check between using normalized internode count or simple integer values 
	Note: We use the normalized internode count
	"""
	#internode_count = ((node1_level - curr_node_level) + (node2_level - curr_node_level)) - 1
	internode_count = ((node1_level - curr_node_level) + (node2_level - curr_node_level) - 1) * 1.0 / taxa_count
	
	"""
	create the couplet key in the dictionary provided it does not exist
	the key has a form (idx1, idx2) where idx1 < idx2
	idx is computed with respect to the COMPLETE_INPUT_TAXA_LIST
	"""
	if (node1_idx < node2_idx):
		target_key = (node1_idx, node2_idx)
	else:
		target_key = (node2_idx, node1_idx)
	if target_key not in TaxaPair_Reln_Dict:
		TaxaPair_Reln_Dict.setdefault(target_key, Reln_TaxaPair())
	
	"""
	add the count of supporting input trees for this couplet
	"""
	TaxaPair_Reln_Dict[target_key]._IncrSupportTreeCount()
	"""
	add the internode count measure with respect to the current input tree
	"""
	TaxaPair_Reln_Dict[target_key]._AddLevel(internode_count)
	"""
	add the XL measure with respect to the current input tree
	"""
	TaxaPair_Reln_Dict[target_key]._AddXLVal(xl_val)
	
	return

#--------------------------------------------------------
"""
this function derives couplet relations belonging to one tree
that is provided as an input argument to this function
"""
def DeriveCoupletRelations(Curr_tree):

	Curr_tree_taxa_count = len(Curr_tree.infer_taxa().labels())

	"""
	traverse the internal nodes of the tree in postorder fashion
	for each internal node x, we would find the couplets whose LCA (lowest common ancestor) node is x
	"""
	for curr_node in Curr_tree.postorder_internal_node_iter():
		"""
		level of current internal node
		used to compute the internode count measure
		"""
		curr_node_level = curr_node.level()
		
		"""
		sourya - check between using normalized excess gene count or simple integer values 
		Currently the normalized XL value is used since it is applicable 
		for the gene trees carrying partially overlapping taxa subsets
		"""
		#xl_val = len(curr_node.leaf_nodes()) - 2
		xl_val = ((len(curr_node.leaf_nodes()) - 2) * 1.0 ) / Curr_tree_taxa_count
		
		"""
		list the leaf and internal children of the current node
		"""
		curr_node_child_leaf_nodes = []
		curr_node_child_internal_nodes = []
		for x in curr_node.child_nodes():
			if (x.is_leaf() == True):
				curr_node_child_leaf_nodes.append(x)
			else:
				curr_node_child_internal_nodes.append(x)
		
		"""
		pair of leaf nodes will be related by sibling relations
		"""
		if (len(curr_node_child_leaf_nodes) > 1):
			for i in range(len(curr_node_child_leaf_nodes) - 1):
				node1 = curr_node_child_leaf_nodes[i]
				node1_level = node1.level()
				node1_idx = COMPLETE_INPUT_TAXA_LIST.index(node1.taxon.label)
				for j in range(i+1, len(curr_node_child_leaf_nodes)):
					node2 = curr_node_child_leaf_nodes[j]
					node2_level = node2.level()
					node2_idx = COMPLETE_INPUT_TAXA_LIST.index(node2.taxon.label)
					Compute_Internode_ExcessGeneLeaf(Curr_tree_taxa_count, xl_val, curr_node_level, \
						node1_level, node2_level, node1_idx, node2_idx)
		
		"""
		one leaf node (direct descendant) and another leaf node (under one internal node)
		will be related by ancestor / descendant relations
		"""
		if (len(curr_node_child_leaf_nodes) > 0) and (len(curr_node_child_internal_nodes) > 0):
			for p in curr_node_child_leaf_nodes:
				node1 = p
				node1_level = node1.level()
				node1_idx = COMPLETE_INPUT_TAXA_LIST.index(node1.taxon.label)
				for q in curr_node_child_internal_nodes:
					for r in q.leaf_nodes():
						node2 = r
						node2_level = node2.level()
						node2_idx = COMPLETE_INPUT_TAXA_LIST.index(node2.taxon.label)
						Compute_Internode_ExcessGeneLeaf(Curr_tree_taxa_count, xl_val, curr_node_level, \
							node1_level, node2_level, node1_idx, node2_idx)
		
		"""
		finally a pair of leaf nodes which are descendant of internal nodes are processed
		"""
		if (len(curr_node_child_internal_nodes) > 1):
			for i in range(len(curr_node_child_internal_nodes) - 1):
				for p in curr_node_child_internal_nodes[i].leaf_nodes():
					node1 = p
					node1_level = node1.level()
					node1_idx = COMPLETE_INPUT_TAXA_LIST.index(node1.taxon.label)
					for j in range(i+1, len(curr_node_child_internal_nodes)):
						for q in curr_node_child_internal_nodes[j].leaf_nodes():
							node2 = q
							node2_level = node2.level()
							node2_idx = COMPLETE_INPUT_TAXA_LIST.index(node2.taxon.label)
							Compute_Internode_ExcessGeneLeaf(Curr_tree_taxa_count, xl_val, curr_node_level, \
								node1_level, node2_level, node1_idx, node2_idx)
		
	return
