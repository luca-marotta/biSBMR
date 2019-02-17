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
	example_R.R
	southernWomen.edgelist
	southernWomen.types
	test.edgelist
	test.types


### Notes
The first three files are the actual implementation: biSBM_R.cpp contains the code for the C++ function called by R; biSBM.h contains all the other C++ functions that are shared by all implementations; finally biSBM.R is the R interface to the C++ workhorse. 

It has the same arguments as the Matlab version and like the Matlab version it accepts as input either the NxN adjacency matrix or the Ex2 edgelist. Put these files in the same folder and sourcing from within R the biSBM.R it should compile and load the function biSBM (you'll need to install the Rcpp package in order to get everything work). 

The file example.R is just a very minimalistic guide to the use of the function contained in biSBM.R, and then there are the edgelist and type files for the SouthernWomen network.

# Windows R wrapper

### Files 

	biSBMWin.h
	biSBMWin_R.cpp
	biSBMWin.R
	biSBM.bat
	-----
	example_R.R
	southernWomen.edgelist
	southernWomen.types
	test.edgelist
	test.types
### Notes

The following steps to make the wrapper work on Windows were tested on Windows 7 / Windows 10 with R 3.3 and Rtools >= 3.3

1. Ensure you installed the R package Rcpp.
2. Install Rtools from https://cran.r-project.org/bin/windows/Rtools (latest version tested: Rtools34)
   1. Allow the installer to edit your PATH. 
   2. Manually Add R to your PATH (see below for details). 
3. Give the file biSBM.bat read and execute permissions (right-click on it, select properties, and edit permissions in the security tab).
4. Run the example script example_R.R to compile the shared library. Make sure to comment the first source command 
##### About Windows PATH

Usually your R installation on Windows is located in the path C:\Program Files\R\R-vers-num\bin, where vers-num is your R release number. In order to call R from MS-DOS you will have to add this path to your system's PATH environment variable. 
To do this, go to Control Panel --> System --> Change Settings and select the Advanced tab. 
Here click on Environment Variables. In the lower panel that will open (System Variables),click on the variable Path and then the button Edit. Add the path to R in your system and click OK.
