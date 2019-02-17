#### compiling the shared library
library(Rcpp);
system("export PKG_CXXFLAGS=`Rscript -e \"Rcpp:::CxxFlags()\"` && R CMD SHLIB biSBM_R.cpp", intern=T);

# INPUTS:
# data     - NxN adjacency matrix *or* Ex2 edge list
# nodeTypes   - Nx1 or 1xN vector of node types. [0,1] or [1,2].
# ka,kb       - Number of communities for types 0,1 respectively (or 1,2).
# deg.corr - boolean: 0 uncorrected biSBM, 1 degree-corrected biSBM. Defaults to 1.
# iterations - positive integer number of Kernighan-Lin iterations. Defaults to 10.
# 
# OUTPUTS:
# g - Nx1 maximum likelihood partition vector


biSBM <- function(data, nodeType, ka, kb, deg.corr=1, iter=10){
  
  
  ### checking ka, kb, deg.corr and iter are positive
  if(any(c(ka<0, kb<0, deg.corr<0, iter<0))) stop("ka, kb, deg.corr and iter must be positive");
  
  ### deg.corr can be just 0 or 1
  if(deg.corr>1 || deg.corr<0) stop("Use deg.corr=1 for the degree corrected SBM and 0 otherwise");
  
  ### loading the shared library
  dyn.load("biSBM_R.so");
  if(!is.loaded("rFunction")) stop("Something is wrong with the C++ shared library");
  
  ### if the adjacency matrix is supplied, transform it in edgelist
  ###
  if(typeof(data)=="list") data <- vapply(data, as.double, numeric(nrow(data)));
  if(ncol(data)>2){
    nn <- ncol(data); ##total nodes
    data[lower.tri(data)] <- 0;
    nodesA <- which(rowSums(data)==0)[1] - 1;
    el1 <- rep(1:nodesA, rowSums(data[1:nodesA, ]));
    el2 <- unlist(apply(data[1:nodesA,], 1, function(x) which(x!=0)));
    data <- cbind(el1, el2);
  }
  
  if(min(data)==0){
    warning("Network indexing starts from 0", call.=F, immediate.=T);
    network <- network + 1;
    
  } 
  
  ## checking dimensions and type for nodeType
  if(typeof(nodeType)=="list") nodeType <- vapply(nodeType, as.double, numeric(nrow(nodeType)));
  if(is.atomic(nodeType)) dim(nodeType) <- c(length(nodeType), 1);
  if((sum(dim(nodeType)) - abs(diff(dim(nodeType))))!=2) stop("nodeType is not a column or row vector");
  if(!is.double(nodeType)) nodeType <- as.double(nodeType);
  
  if(length(nodeType)!= length(unique(as.numeric(data)))) stop("The number of network's nodes and node types are different. Please, 
          check the connected components in your network");
  
  
  ### calling C++ algorithm and returning the partition
  res <-integer(length(nodeType));
  .C("rFunction", as.numeric(data), nrow(data), nodeType, length(nodeType), as.integer(ka), 
     as.integer(kb), as.integer(deg.corr), as.integer(iter), res=res)$res + 1;
}