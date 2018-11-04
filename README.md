*********************************
IDXL (Species tree estimation using Internode Distance and eXcess gene Lineage)
*********************************

IDXL is a python based tool for computing species tree from a set of incongruent gene trees 
with Incomplete Lineage Sorting (ILS). Following measures between individual couplets are compurted 
for species tree estimation.

A) Internode distance between individual couplets (proposed in NJst approach (Liu et. al. 2011))

B) A novel couplet based measure, termed as the "Excess gene lineage", computed for all the input gene trees.

These measures are used to form respective distance matrices, which are subsequently used for NJ based species tree construction.

Method Description
-----------------------

Input
-----------

A collection of gene trees (formed by sampling individual genes from a group of taxa) with overlapping taxa set, having topological incongruence due to Incomplete Lineage Sorting (ILS). Gene trees may or may not be weighted; our species tree estimation method does not consider the branch lengths of the tree.

Input gene trees can be either in NEWICK format or in NEXUS format. However, all the gene trees should have identical input formats. They should be placed in a standard tree list file, according to the syntax of NEXUS or NEWICK formats. Such a tree list text file is to be provided as an input of this executable.

Output
----------

A species tree covering all the taxa of the gene trees. Output species tree is generated in the NEWICK format.

Installation
--------------

IDXL is developed in Linux Systems (Ubuntu 14.04), using Python 2.7. It is tested and meant for systems having linux OS (Fedora / Ubuntu). User should download the zipped archieve from GitHub / clone this repository. In addition, they require 
to install the following packages:

1) Python 2.7 (available in Ubuntu, by default) 

Note: The current version does not support python 3. We plan to release a new version supporting python 3.

2) Dendropy 3.12.0 ( available on the link: https://pythonhosted.org/DendroPy/ ) 

Note: there is a new release of Dendropy 4.0 but we have used 3.12.0 for the implementation. Future release would support Dendropy latest version.

3) Numpy ( available on the link: http://www.numpy.org/ )

User should install Numpy using pip, containing the latest version. 


Command line options
----------------------

The file IDXL.py is the main executable file. User should execute this file as:

python IDXL.py [options]

with the following options:

-h, --help show this help message and exit

-I INP_FILENAME, --INPFILE=INP_FILENAME

	name of the input file containing gene trees

-O OUT_FILENAME, --OUTFILE=OUT_FILENAME

	name of the file to contain output species tree

-p INP_FILE_FORMAT, --inpform=INP_FILE_FORMAT

	1 - input file format is NEWICK (default)
	2 - input file format is NEXUS

-m METHOD_TYPE, --method=METHOD_TYPE

	Recommended and default: 1. User should use the default setting.
		
-d DISTMAT_TYPE, --distmat=DISTMAT_TYPE

	Recommended and default: 1. User should use the default setting.		

-r TAXON_NAME, --ROOT=TAXON_NAME

	User can specify a taxon name (string) which will be used to root the output species tree. 
	Useful to benchmark the output species tree with any reference rooted species tree.


Example of a command (followed for the results published in the manuscript)

python IDXL.py -I source_tree_input.txt -O output_file_name -p1 -r "sample_taxon" (to be replaced by a taxon name to be used as a root)

Output files
-------------

If output_file_name is provided, corresponding directory is created (if not exists) and a separate output log file 
'Complete_Desription.txt' is also created within that directory.

If -O option is not specified, a folder named IDXL_M[METHOD_TYPE]_D[DISTMAT_TYPE] is created in the same directory 
containing the input trees. A file named "outtree_Newick.tre" stores the final tree.


Citation
---------

Upon using this package, user should cite the following article:

Sourya Bhattacharyya, Jayanta Mukherjee, IDXL: Species Tree Inference Using Internode Distance and Excess Gene Leaf Count, Journal of Molecular Evolution (Springer), volume 85, issue 1-2, pp. 57-78, 2017, DOI: 10.1007/s00239-017-9807-7


For any queries, please contact
------------------------------

Sourya Bhattacharyya

La Jolla Institute of Allergy and Immunology (LIAI)

La Jolla, CA 92037, USA

sourya.bhatta@gmail.com

Jayanta Mukherjee 

Department of Computer Science and Engineering 

Indian Institute of Technology Kharagpur 

WB 721302, India 

jay@cse.iitkgp.ernet.in




