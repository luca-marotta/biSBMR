# biSBMR
R wrapper for the bipartite Stochastic block model, based on the publication by Daniel B. Larremore, Aaron Clauset, 
and Abigail Z. Jacobs, "Efficiently inferring community structure in bipartite networks", Physical Review E 90(1), 012805 (2014). 

Matlab wrappers and pure C++ implementations are available on Daniel's personal webpage: http://danlarremore.com/bipartiteSBM/index.html

# Unix R wrapper

### Files 

	biSBM.h
	biSBM_R.cpp
	biSBM.R
	-----
	example.R
	southernWomen.edgelist
	southernWomen.types
	test.edgelist
	test.types


### Notes

The first three files are the actual implementation: biSBM_R.cpp contains the code for the C++ function called by R; biSBM.h contains all the other C++ functions that are shared by all implementations; finally biSBM.R is the R interface to the C++ workhorse. 

It has the same arguments as the Matlab version and like the Matlab version it accepts as input either the NxN adjacency matrix or the Ex2 edgelist. Put these files in the same folder and sourcing from within R the biSBM.R it should compile and load the function biSBM (you'll need to install the Rcpp package in order to get everything work). 

The file example.R is just a very minimalistic guide to the use of the function contained in biSBM.R, and then there are the edgelist and type files for the SouthernWomen network.
