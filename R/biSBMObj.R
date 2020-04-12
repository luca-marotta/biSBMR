############# biSBMObj class constructor ################à

new_biSBMObj <- function(partition, score, nodeType){
  structure(partition, class = "biSBMObj", score = score,  nodeType = nodeType)
}

summary.biSBMObj <- function(x){
  cat("**** Bipartite Stochastic bloc model object ****\n\n\n")

  ### summary for input network
  nodeType <- attr(x, "nodeType")
  cat("Size of node sets for input network:\n")
  setSize <- table(nodeType)
  names(setSize) <- c("type 1", "type 2")
  print(setSize)

  ### summary for communities
  cat("Communities type 1:\n")
  seq_type1 <- seq_len(setSize["type 1"])
  communities_type1 <- split(x = seq_type1, f = x[seq_type1])
  print(communities_type1)

  cat("Communities type 2:\n")
  seq_type2 <- (setSize["type 1"] + 1):(length(nodeType))
  communities_type2 <- split(x = seq_type2, f = x[seq_type2])
  print(communities_type2)

  score <- attr(x, "score")
  print(paste0("Likelihood score is: ", score))

  invisible(list(partition = x, score = score, communities = list(type1 = communities_type1, type2 = communities_type2)))
}

getPartition <- function(x){
  UseMethod("getPartition")
}

getPartition.biSBMObj <- function(x){
  attributes(x) <- NULL
  x
}

getScore <- function(x){
  UseMethod("getScore")
}

getScore.biSBMObj <- function(x){
  attr(x, "score")
}
