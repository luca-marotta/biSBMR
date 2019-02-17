// biSBM v1.2
//
// R wrapper by Dr. Luca Marotta, marotta.luca@gmail.com. Thank you!
//
//Daniel Larremore
//Harvard School of Public Health
//July 29, 2014
//http://danlarremore.com/bipartiteSBM
//daniel.larremore@gmail.com
//
//biSBM - a method for community detection in bipartite networks, based on the publication:
//"Efficiently inferring community structure in bipartite networks"
//Daniel B. Larremore, Aaron Clauset, and Abigail Z. Jacobs
//Physical Review E 90(1), 012805 (2014).
//http://danlarremore.com/pdf/2014_LCJ_EfficientlyInferringCommunityStructureInBipartiteNetworks_PRE.pdf
//
// Please do not distribute without contacting the author above at daniel.larremore@gmail.com
// If a bug is located within the code, please contact the author, to correct the official version!
//
// This code is based on code written by Brian Karrer for the stochastic block model, http://arxiv.org/abs/1104.3590
// You can download that code at http://www-personal.umich.edu/~mejn/dcbm/
//

#include <Rcpp.h> 
#include <R.h> 
#include <cstdio>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <math.h>
#include <fstream>
#include <limits>
#include <sys/stat.h>
#include <unistd.h>
#include <ctime>
#include <vector>
#include "biSBM.h"

//using namespace Rcpp;
using namespace std;

/********** biSBM **********/
extern "C" 
{ 
  void rFunction(double* el, int* Nel, double* nodes, int* Nnodes, int *nA, int *nB, int *dc, int *iter, int* g)
  {
    srandom(time(NULL));
    /*
     * 0 edgelist
     * 1 types
     * 2 KA
     * 3 KB
     * 4 DegreeCorrect
     * 5 KLSteps
     */
    int KA, KB;
    double *types;
    //double *edges;
    double *e;
     Types.clear();
     Comms.clear();
    
    /**********  0 edgelist  **********/
     size_t mrows =(long int) *Nel;
     e = el;
    
     GetTheNetworkEdges(e,mrows);
    
    
    /**********  1 types  **********/
    //     types = mxGetPr(prhs[1]);
    //     Nodes = mxGetM(prhs[1]); // GLOBAL - Number of nodes
    types = nodes;
    Nodes = *Nnodes;
    int maxTypes = 0;
    for (unsigned int i=0; i<Nodes; ++i)
    {
      Types.push_back(types[i]-1); // GLOBAL - Vertex types
      maxTypes = std::max(maxTypes,Types[i]+1);
    }
    std::vector<int> tally (maxTypes+1,0);
    for (unsigned int i=0; i<Types.size(); ++i)
    {
      tally[Types[i]]++;
    }
    /**********  2 KA  **********/
     KA = *nA;
     Rcpp::Rcout<<"KA "<<KA<<std::endl;
    /**********  3 KB  **********/
     KB = *nB;
     Rcpp::Rcout<<"KB "<<KB<<std::endl;
    /**********  4 isDegreeCorrect  **********/
    isDegreeCorrect = *dc;
    /**********  5 KLPerNetwork  **********/
    KLPerNetwork = *iter;
    /**********  Output  **********/
    int counter = 0;
    std::vector<int> commlist;
    for (unsigned int q=0; q<KA; ++q)
    {
      commlist.push_back(counter);
      counter++;
    }
    Comms.push_back(commlist);
    commlist.clear();
    for (unsigned int q=0; q<KB; ++q)
    {
      commlist.push_back(counter);
      counter++;
    }
    Comms.push_back(commlist);
    commlist.clear();
    MaxComms = counter;
    
    Rcpp::Rcout<<"\nCalling biSBM with the following parameters.\n";
    Rcpp::Rcout<<"KA:\t"<<KA<<std::endl;
    Rcpp::Rcout<<"KB:\t"<<KB<<std::endl;
    Rcpp::Rcout<<"NA:\t"<<tally[0]<<std::endl;
    Rcpp::Rcout<<"NB:\t"<<tally[1]<<std::endl;
    Rcpp::Rcout<<"Type A Communities: ";
    unsigned int i=0;
    for (unsigned int j=0; j<Comms[i].size(); ++j) {
      Rcpp::Rcout<<Comms[i][j]+1<<",";
    }
    Rcpp::Rcout<<"\n";
    Rcpp::Rcout<<"Type B Communities: ";
    i=1;
    for (unsigned int j=0; j<Comms[i].size(); ++j) {
      Rcpp::Rcout<<Comms[i][j]+1<<",";
    }
    Rcpp::Rcout<<"\n";
    Rcpp::Rcout<<"DegreeCorrect:\t"<<isDegreeCorrect<<std::endl;
    Rcpp::Rcout<<"Edges:\t"<<mrows<<std::endl;
    //mexEvalString("drawnow;");
    
    /**********  Call the biSBM subroutine.  **********/
    biSBM();
    
    /**********  Put the output into g.  **********/
    for (unsigned int i=0; i<Nodes; ++i)
    {
      g[i] = BestState[i];
      
    }
    
   }
  
} 
