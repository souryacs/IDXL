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

Method implemented
---------------------

First we create the following two distance matrices:

I = matrix of average internode distance for individual couplets, as proposed in the method NJst (Liu et. al. 2011)

X = matrix of excess gene lineage for individual couplets (mean value for individual couplets)

      
      







We have implemented following 2 kinds of methods for species tree estimation methods.

A) MNJstXL (or MedNJstXL): 

This approach computes the   separate 

Let 
      
      
      
      

Neighbor Joining (NJ) based species tree construction methods, compute the relative distance of individual 
couplets, with respect to all other couplets, in terms of the given distance matrix.

Let, I(x,y) = relative internode count of the couplet (x,y) (with respect to other couplets (x,z) and (y,z) for 
all other taxon z belonging to the set of input gene trees)
 
 Similarly, let X(x,y) = relative excess gene leaf count of the couplet (x,y) (with respect to other couplets (x,z) and (y,z) for 
all other taxon z belonging to the set of input gene trees)

MNJstXL computes the ranks (position in the sorted list of ascending order) of I(x,y) and X(x,y) and computes the 
aggregate of these ranks. Suppose, r(x,y) denotes the agregated rank for the couplet (x,y).

The couplet (x,y) having minimum r(x,y) gets selected for agglomeration in the current iteration.

B) PNJstXL (or ProdNJstXL): 

Here, we employ the product of I(x,y) and X(x,y) as a criterion for the selection of couplet. 
A couplet (x,y) which is a possible candidate for agglomeration, should have low values of both of these measures.
So, their product would also be very low.

So, the couplet (x,y) having the minimum (I(x,y) * X(x,y)) gets selected for agglomeration in the current iteration.

This method produces better performance.

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

If m = 1, a directory “MNJSTXL” is created within the same directory containing the input treelist file. 
Within this new created directory, one file 'outtree_newick.tre' is created, which contains the derived species tree. 
Another text file named 'Complete_Desription.txt' is created, which contains execution and timing information 
for the method. 

For m = 2,  a directory 'PNJSTXL' is created within 
the same directory containing the input treelist file. Above mentioned files within the new directory are 
generated as per the execution.

*********************************
For any queries, please contact
*********************************

Sourya Bhattacharyya 
Department of Computer Science and Engineering
Indian Institute of Technology Kharagpur
<sourya.bhatta@gmail.com>



