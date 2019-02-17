README for biSBM v1.2


R wrapper by Dr. Luca Marotta, marotta.luca@gmail.com. Thank you!


Daniel Larremore
Harvard School of Public Health
July 29, 2014
http://danlarremore.com/bipartiteSBM
daniel.larremore@gmail.com

biSBM - a method for community detection in bipartite networks, based on the publication: 
	"Efficiently inferring community structure in bipartite networks"
	Daniel B. Larremore, Aaron Clauset, and Abigail Z. Jacobs 
	Physical Review E 90(1), 012805 (2014).
	http://danlarremore.com/pdf/2014_LCJ_EfficientlyInferringCommunityStructureInBipartiteNetworks_PRE.pdf

/***** MATLAB vs R vs C++ *****/

This file explains how to use the R version of the biSBM code. If you prefer using a wrapper in MATLAB or pure C++, see README_MATLAB.txt or README_cplusplus.txt. 

/***** FILES *****/

You'll need these files in your R directory:
	biSBM.h
	biSBM_R.cpp
	biSBM.R
	-----
	example.R
	southernWomen.edgelist
	southernWomen.types
	test.edgelist
	test.types


/***** NOTES *****/

The first three files are the actual implementation: biSBM_R.cpp contains the code for the C++ function called by R; biSBM.h contains all the other C++ functions that are shared by all implementations; finally biSBM.R is the R interface to the C++ workhorse. 

It has the same arguments as the Matlab version and like the Matlab version it accepts as input either the NxN adjacency matrix or the Ex2 edgelist. Put these files in the same folder and sourcing from within R the biSBM.R it should compile and load the function biSBM (you'll need to install the Rcpp package in order to get everything work). 
The R interface still needs some make-up (a few more checks on the input data have to be added), nonetheless I'm sending it to you 'cause I didn't have the chance to test it on any other computer but mine. 

The file example.R is just a very minimalistic guide to the use of the function contained in biSBM.R, and then there are the edgelist and type files for the SouthernWomen network. I had to modify the edgelist because the one in your data contains every edge two times and my implemention currently works properly only with the Ex2 edgelist.

Please don't hesitate to let me know if something is not clear, or the implementation does not run.