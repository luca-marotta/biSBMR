##### example for biSBM.R #####
### loading and compiling biSBM function
setwd("R")
if(Sys.info()["sysname"]=="Windows") {
  source("biSBMWin.R") } else {
    source("biSBM.R")
  }

##### TEST 1 - Southern Women #####
edges <- read.table("southernWomen.edgelist");  ## edgelist
types <- read.table("southernWomen.types");     ## nodeType
### call to the function: by default deg.corr=1 and iter = 10.
result <- biSBM(data = edges, nodeType = types, ka = 2, kb = 3, deg.corr = 1, iter = 3);
##### Should produce output: #####
# Calling biSBM with the following parameters.
# KA:  2
# KB:	3
# NA:	18
# NB:	14
# Type A Communities: 1,2,
# Type B Communities: 3,4,5,
# DegreeCorrect:	1
# Edges:	89
# >1, 0, -367.658
# >2, 0, -367.658
# >3, 0, -367.658
# Final Score: -367.658

### result is an object of class biSBMObj. The summary() generic function prints info
### about the optimization: the number of nodes in each node set, the node ids split by community
### and the score.
summary(result)

### the vector containing the partition can be extracted with the following function:
getPartition(result)

### To get the final likelihood score use the following function:
getScore(result)


##### TEST 2 - Difficult Case (from paper) #####
edges <- read.table("test.edgeList");  ## edgelist
edges = edges[,1:2]; ## get rid of the edge weights (which are all 1).
types <- read.table("test.types");     ## nodeType
g <- biSBM(data = edges, nodeType = types, ka = 2, kb = 3, deg.corr = 1, iter = 3);

summary(g)
##### Should produce output: #####
# Calling biSBM with the following parameters.
# KA:  2
# KB:	3
# NA:	700
# NB:	300
# Type A Communities: 1,2,
# Type B Communities: 3,4,5,
# DegreeCorrect:	1
# Edges:	7696
# >1, 1, -65525.8
# >2, 1, -65525.8
# >3, 2, -65525.8
# Final Score: -65525.8

