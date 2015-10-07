import Header
from Header import *  
import UtilFunc
from UtilFunc import *

#--------------------------------------------------------
# this function defines accumulated branch count between two nodes
# via the path through their LCA node
def DefineAccBranchXL(xl_val, curr_node_level, node1, node2):
	"""
	*** important 
	sourya - instead of using the internode count, we are using the number of branches
	which is 1 more than the internode count
	somehow, enforcing greater value increases performance
	"""
	sum_of_branch_count = ((node1.level() - curr_node_level) + (node2.level() - curr_node_level))

	key1 = (node1.taxon.label, node2.taxon.label)
	key2 = (node2.taxon.label, node1.taxon.label)
	if key1 in TaxaPair_Reln_Dict:
		TaxaPair_Reln_Dict[key1]._IncrSupportTreeCount()
		TaxaPair_Reln_Dict[key1]._AddLevel(sum_of_branch_count)
		TaxaPair_Reln_Dict[key1]._AddXLVal(xl_val)
	elif key2 in TaxaPair_Reln_Dict:
		TaxaPair_Reln_Dict[key2]._IncrSupportTreeCount()
		TaxaPair_Reln_Dict[key2]._AddLevel(sum_of_branch_count)
		TaxaPair_Reln_Dict[key2]._AddXLVal(xl_val)
	else:
		TaxaPair_Reln_Dict.setdefault(key1, Reln_TaxaPair())
		TaxaPair_Reln_Dict[key1]._IncrSupportTreeCount()
		TaxaPair_Reln_Dict[key1]._AddLevel(sum_of_branch_count)
		TaxaPair_Reln_Dict[key1]._AddXLVal(xl_val)

	return

#--------------------------------------------------------
# this function defines accumulated branch count between two nodes
# via the path through their LCA node
def DefineAccBranch(curr_node_level, node1, node2):
	"""
	*** important 
	sourya - instead of using the internode count, we are using the number of branches
	which is 1 more than the internode count
	somehow, enforcing greater value increases performance
	"""	
	sum_of_branch_count = ((node1.level() - curr_node_level) + (node2.level() - curr_node_level))

	key1 = (node1.taxon.label, node2.taxon.label)
	key2 = (node2.taxon.label, node1.taxon.label)
	if key1 in TaxaPair_Reln_Dict:
		TaxaPair_Reln_Dict[key1]._IncrSupportTreeCount()
		TaxaPair_Reln_Dict[key1]._AddLevel(sum_of_branch_count)
	elif key2 in TaxaPair_Reln_Dict:
		TaxaPair_Reln_Dict[key2]._IncrSupportTreeCount()
		TaxaPair_Reln_Dict[key2]._AddLevel(sum_of_branch_count)
	else:
		TaxaPair_Reln_Dict.setdefault(key1, Reln_TaxaPair())
		TaxaPair_Reln_Dict[key1]._IncrSupportTreeCount()
		TaxaPair_Reln_Dict[key1]._AddLevel(sum_of_branch_count)

	return

#--------------------------------------------------------
# this function derives coupket relations belonging to one tree
# that is provided as an input argument to this function
def DeriveCoupletRelations(Curr_tree, METHOD_USED):

	# traverse the internal nodes of the tree in postorder fashion
	for curr_node in Curr_tree.postorder_internal_node_iter():
		# compute the rank associated with this node
		curr_node_level = curr_node.level()
		xl_val = len(curr_node.leaf_nodes()) - 2
		
		# list the leaf and internal children of the current node
		curr_node_child_leaf_nodes = []
		curr_node_child_internal_nodes = []
		for x in curr_node.child_nodes():
			if (x.is_leaf() == True):
				curr_node_child_leaf_nodes.append(x)
			else:
				curr_node_child_internal_nodes.append(x)
		
		# pair of leaf nodes will be related by sibling relations
		if (len(curr_node_child_leaf_nodes) > 1):
			for i in range(len(curr_node_child_leaf_nodes) - 1):
				for j in range(i+1, len(curr_node_child_leaf_nodes)):
					if (METHOD_USED == NJ_ST) or (METHOD_USED == M_NJ_ST):
						DefineAccBranch(curr_node_level, curr_node_child_leaf_nodes[i], curr_node_child_leaf_nodes[j])
					elif (METHOD_USED == M_NJ_ST_XL):
						DefineAccBranchXL(xl_val, curr_node_level, curr_node_child_leaf_nodes[i], curr_node_child_leaf_nodes[j])
		
		# one leaf node (direct descendant) and another leaf node (under one internal node)
		# will be related by ancestor / descendant relations
		if (len(curr_node_child_leaf_nodes) > 0) and (len(curr_node_child_internal_nodes) > 0):
			for p in curr_node_child_leaf_nodes:
				for q in curr_node_child_internal_nodes:
					for r in q.leaf_nodes():
						if (METHOD_USED == NJ_ST) or (METHOD_USED == M_NJ_ST):
							DefineAccBranch(curr_node_level, p, r)
						elif (METHOD_USED == M_NJ_ST_XL):
							DefineAccBranchXL(xl_val, curr_node_level, p, r)
				
		# finally a pair of leaf nodes which are descendant of internal nodes will be related by NO_EDGE relation
		if (len(curr_node_child_internal_nodes) > 1):
			for i in range(len(curr_node_child_internal_nodes) - 1):
				for j in range(i+1, len(curr_node_child_internal_nodes)):
					for p in curr_node_child_internal_nodes[i].leaf_nodes():
						for q in curr_node_child_internal_nodes[j].leaf_nodes():
							if (METHOD_USED == NJ_ST) or (METHOD_USED == M_NJ_ST):
								DefineAccBranch(curr_node_level, p, q)
							elif (METHOD_USED == M_NJ_ST_XL):
								DefineAccBranchXL(xl_val, curr_node_level, p, q)
		
	return
