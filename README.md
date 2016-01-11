*********************************
MNJSTXL & PNJSTXL
*********************************

******************************
Species tree estimation using internode count and excess gene leaves information
******************************

MNJSTXL & PNJSTXL is a python based tool for computing species tree from a set of incongruent gene trees 
with Incomplete Lineage Sorting (ILS). Following measures between individual couplets are compurted 
for species tree estimation.

A) Internode count between individual couplets (proposed in NJst approach (Liu et. al. 2011))

B) Excess gene leaves count between individual couplets, computed for all the input gene trees.

These measures are used to form respective distance matrices, which will be subsequently used for NJ 
based species tree construction.

Description
-----------------------

Input
-----------

A collection of gene trees (formed by sampling individual genes from a group of taxa) with overlapping taxa set, having topological incongruence 
due to Incomplete Lineage Sorting (ILS). Gene trees may or may not be weighted; our species tree estimation method 
does not consider the branch lengths of the tree.

Input gene trees can be either in NEWICK format or in NEXUS format. 
However, all the gene trees should have identical input formats. They should be placed in a 
standard tree list file, according to the syntax of NEXUS or NEWICK formats. Such a tree list 
text file is to be provided as an input of this executable.

Output
--------

A species tree covering all the taxa of the gene trees. Output species tree 
is generated in the NEWICK format.

Methods implemented
------------------------------

We have implemented following 2 kinds of methods for species tree estimation methods.

A) MNJstXL (or MedNJstXL): Internode count and excess gene leaf count between individual couplets are 
computed with respect to input gene trees. Two separate distance matrices are constructed next. One contains 
the couplet based average internode count measure. Another matrix contains the minimum of 
median and average of excess gene leaf count, computed for all the gene trees. These two distance matrices 
are used separately based on the principle of lowest agregated rank, to construct the species tree using NJ.
This method is termed as MedNJstXL or MNJstXL.

B) PNJstXL (or ProdNJstXL): Here, we use the product of internode count and excess gene leaf count to construct a single distance matrix, 
for NJ based species tree construction.

*********************************
Dependencies
*********************************

These methods is developed in Linux Systems (Ubuntu 14.04), using Python 2.7. It is tested and meant for systems 
having linux OS (Fedora / Ubuntu).

User needs to install following before using this package:

1) Python 2.7 (available in Ubuntu, by default) 

Note: We have not tested the code on Python 3. Any user having Python 3 environment need to 
check the correct execution of our code, and optionally needs to upgrade it accordingly.
We plan to support Python 3 environment in some future release.

2) Dendropy 3.12.0 ( available on the link: https://pythonhosted.org/DendroPy/ ) 

Note: there is a new release of Dendropy 4.0 but we have used 3.12.0 for the implementation. We 
did not upgrade the code for Dendropy 4.0 support, so any user having this new version of Dendropy 
might need to check the functionalities of MNJstXL and possibly upgrade / replace / edit few 
dendropy related functions. So, we recommend users to use the earlier version of Dendropy, to avoid any conflict.

Support for Dendropy 4 and corresponding update of code will be done in a future release.

3) Numpy ( available on the link: http://www.numpy.org/ )

User can install Numpy using pip (python software downloader tool) module, which contains the latest 
Numpy module in it. We found that Numpy module in the traditional Apt-Get repository is of lower version.

***************
Command line options
****************

After downloading the current archieve in a zipped format, extracting the archieve reveals a file MNJSTXL.py which 
is the main executable file. Assuming current directory contains this python file, 
following command is to be executed:

./MNJSTXL [options]

Details of the options are mentioned below:

-h, --help show this help message and exit

-I INP_FILENAME, --INPFILE=INP_FILENAME

                name of the input file containing gene trees

-O OUT_FILENAME, --OUTFILE=OUT_FILENAME

                name of the output file to contain target species tree

-p INP_FILE_FORMAT, --inpform=INP_FILE_FORMAT

                1 - input file format is NEWICK (default)
                2 - input file format is NEXUS

-m METHOD_TYPE, --method=METHOD_TYPE

                1 - MNJstXL  (Default Method)  

                2 - PNJstXL 

-r TAXON_NAME, --ROOT=TAXON_NAME

		User can specify a taxon name to root the output species tree with that specified taxon.


Example of a command (followed for the results published in the manuscript)

./MNJSTXL -I source_tree_input.txt -p1 -m2

command descriptions:

1) -I specifies the input filename

2) source_tree_input.txt : contains the input collection of gene trees

3) -p option is for specifying the input tree format input file contains the trees in NEWICK format, 
as specified by the option (-p1) (1 stands for newick)

4) -m option is used to specify the species tree construction method. 
Here 2 is used to denote that PNJstXL is employed. The value can vary between 1 and 2.

In addition, the package contains another option: -O 'output_file_name'

Here, user can specify the output file name containing the derived species tree file.

If no such option is provided, our method performs the following operations:

If m = 1, a directory “MedNJSTXL” is created within the same directory containing the input treelist file. 
Within this new created directory, one file 'outtree_newick.tre' is created, which contains the derived species tree. 
Another text file named 'Complete_Desription.txt' is created, which contains execution and timing information 
for the method. For m = 2,  a directory 'ProdNJSTXL' is created within 
the same directory containing the input treelist file. Above mentioned files within the new directory are 
generated as per the execution.

*********************************
For any queries, please contact
*********************************

Sourya Bhattacharyya 
Department of Computer Science and Engineering
Indian Institute of Technology Kharagpur
<sourya.bhatta@gmail.com>



