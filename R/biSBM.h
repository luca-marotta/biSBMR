// biSBM v1.2
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <math.h>
#include <fstream>
#include <limits>
#include <sys/stat.h>
#include <unistd.h>
#include <ctime>

using namespace std;

// /* const values to control maximum sizes */
const int constComms = 1000;
const long int constNodes = 1000000;
const int MAXEDGES = 10000000;  // this is the maximum number of edges

/* empty global declarations */
long int Nodes;
unsigned int MaxComms;
bool isDegreeCorrect; // 0 means don't correct for the degrees, 1 means do correct for the degrees.
static vector<vector<int> > Comms;
static vector<int> Types;

/* Number of random initializations (default) */
unsigned int KLPerNetwork = 1;

/********** GLOBAL VARIABLES **********/
bool TrueCommsAvailable = 0; // (default) Changes to 1 if user passes in true comms.
bool InitializationOption = 0; // (default) May be changed to 1 by user.

// OVERALL GENERAL VARIABLES NEEDED FOR ALL OF THE CODE
long int *AdjList[constNodes];
long int LastEmpty[constNodes];
long int Degree[constNodes];  // Degree of nodes in the network
long int EdgeList[MAXEDGES][2]; // The first is the maximum number of edges to read in

// FOR KL
unsigned int CurrentState[constNodes];
int BestState[constNodes];
int ChangeSet[constNodes];
int UpdateIndex[constNodes];
int TrueState[constNodes]; // This records the true communities if they exist read in from the file

double TwiceEdges = 0;
double MaxScore = 0;

int BestCommVertices[constComms];
int BestCommStubs[constComms];
int BestEdgeMatrix[constComms][constComms];

int CurrentCommVertices[constComms];
int CurrentCommStubs[constComms];
int CurrentEdgeMatrix[constComms][constComms];

int AdjustmentFrom[constComms];
int AdjustmentDestination[constComms];
int TempNeighborSet[2][constComms];  // the first entry lists the comm and the second entry lists the number of edges to that comm
int NeighborSet[2][constComms]; // this is what we record and use
int SelfEdgeCounter = 0; // this records self-edges to make sure that they are counted correctly
unsigned int ActualDiffComms = 0; // this records the number of different comms in neighborhood

// For reporting best state
int SavedState[constNodes];
int SavedCommVertices[constComms];
int SavedCommStubs[constComms];
int SavedEdgeMatrix[constComms][constComms];
double NMIValue = 0;
double VIValue = 0;
double HighestScore = 0;

/* function declarations */
void freegraph(); // gets rid of the graph pointers at end
void GetTheNetworkEdges(double[], int mwSize);
void RunKL(); // runs Kernighan-Lin once.
void Initialize(); // initializes the data structures for KL
double ComputeInitialScore(); // computes the initial score after initialization
void ComputeNeighborSet(int vertex, int option);  // computes the neighbor set for a vertex, using either best or currentstates
double ComputeProposal(int vertex, int from, int destination); // computes the value of the particular change
void UpdateMatrices(int vertex, int option, int from, int destination); // this updates either the best
//or current matrices based on moving the vertex from from to destination
double LogFunction(double x); // this returns x*log(x) and zero if x=0
void PrintResults(string); // prints out the resulting graph for now
double ComputeVI();
double ComputeNMI();
double Entropy(int entoption);
double JointEntropy();
void biSBM(void);

void biSBM()
{
  HighestScore = -numeric_limits<double>::max( );
  VIValue = 0;
  NMIValue = 0;
  time_t startTime = time(NULL);
  unsigned int i,j,k;
   for(j=0; j < KLPerNetwork; j++)
   {
     RunKL();
     //KL,dt,L:
     cout<<">"<<j+1<<", "<<difftime(time(NULL),startTime)<<", "<<max(MaxScore,HighestScore)<<"\n";
    if(MaxScore >= HighestScore)
    {
      HighestScore = MaxScore;
      if(TrueCommsAvailable == 1)
      {
	VIValue = ComputeVI();
	NMIValue = ComputeNMI();
	cout << "VI Value: " << VIValue << " NMI Value: " << NMIValue << endl;
      }
      for(i=0; i < MaxComms; i++)
      {
	SavedCommVertices[i] = BestCommVertices[i];
	SavedCommStubs[i] = BestCommStubs[i];
	for(k=0; k < MaxComms; k++)
	  SavedEdgeMatrix[i][k] = BestEdgeMatrix[i][k];
      }
      for(i=0; i < Nodes; i++)
      {
	SavedState[i] = BestState[i];
      }
    }
  }
  
  // because PrintResults are written for best values we copy them
  // back over from the saved values which are the best ones.
  for(i=0; i < MaxComms; i++)
  {
    BestCommVertices[i] = SavedCommVertices[i];
    BestCommStubs[i] = SavedCommStubs[i];
    for(k=0; k < MaxComms; k++)
      BestEdgeMatrix[i][k] = SavedEdgeMatrix[i][k];
  }
  for(i=0; i < Nodes; i++)
  {
    BestState[i] = SavedState[i];
  }
  cout << "Final Score: " << ComputeInitialScore() << endl;
  
  freegraph();
}

//void GetTheNetworkEdges(double *e, int mrows)
// NOTE: OVERLOADED. CPP uses GetTheNetworkEdges(string);
void GetTheNetworkEdges(double *e, size_t mrows)
{
  unsigned int counter = 0;
  unsigned int ii;
  for (ii=0; ii<mrows; ++ii)
  {
    EdgeList[ii][0] = (long int)e[ii]-1;
    EdgeList[ii][1] = (long int)e[mrows+ii]-1;
    counter = counter+1;
    Nodes = max(Nodes,(long int)e[ii]);
    Nodes = max(Nodes,(long int)e[mrows+ii]);
  }
  TwiceEdges = 2*counter; // GLOBAL
  
  if(TrueCommsAvailable == 1)
  {
    // DBL
    //        InputFile2.open(trueCommsName.c_str());
    //        if (!InputFile2)
    //        {
      //            cout << "Error in opening file";
      //            cin.get();
      //            return;
      //        }
      //
      //        for(i=0; i < Nodes; i++)
      //            TrueState[i] = -1;
      //
      //        while(std::getline(InputFile2, lineread)) // Read line by line
      //        {
	//            buffer = new char [lineread.size()+1];
	//            strcpy(buffer, lineread.c_str());
	//            //  cout << buffer << endl;
	//            sscanf(buffer, "%ld", &entry1);
	//            //  sscanf(buffer, "n%d,%*[^,],%d", &ignore, &entry1); //%*s
	//            // entry1 = entry1+1;
	//            //sscanf(buffer, "%d %d", &entry1, &entry2);
	//            // TrueState[entry1-1] = entry2;
	//            TrueState[counter2] = entry1;
	//
	//            counter2 = counter2+1;
	//            delete[] buffer;
	//        }
	//        InputFile2.close();
	//
	//        for(i=0; i < Nodes; i++)
	//        {
	  //            if((TrueState[i] == -1) || (TrueState[i] >= MaxComms))
	  //            {
	    //                cout << "STOP A VERTEX WAS NOT LABELED OR WAS LABELED INCORRECTLY." << TrueState[i] << " "  << i << endl;
	    //                cin.get();
	    //            }
	    //        }
  }
  // We start the degree values and LastEmpty all at zero
  for(unsigned int i=0; i<Nodes; ++i)
  {
    Degree[i] = 0;
    LastEmpty[i] = 0;
  }
  // First we count the degrees by scanning through the list once
  for(unsigned int i=0; i<counter; ++i)
  {
    Degree[EdgeList[i][0]]++;
    Degree[EdgeList[i][1]]++;
  }
  // Now we make space in the adjacency lists
  for(unsigned int i=0; i<Nodes; ++i)
  {
    AdjList[i] = new long int [Degree[i]];
  }
  // Now we read the edges into the adjacency lists utilizing lastempty to put them into
  // the proper spots
  for(unsigned int i=0; i<counter; ++i)
  {
    AdjList[EdgeList[i][0]][LastEmpty[EdgeList[i][0]]] = EdgeList[i][1];
    LastEmpty[EdgeList[i][0]]++;
    
    AdjList[EdgeList[i][1]][LastEmpty[EdgeList[i][1]]] = EdgeList[i][0];
    LastEmpty[EdgeList[i][1]]++;
  }
  return;
}

// This function deletes the graph from memory.
void freegraph()
{
  long int i;//-Wunused,j;
  
  for(i=0; i < Nodes; i++)
    delete [] AdjList[i];
  
  return;
}

void RunKL()
{
  unsigned int i,j,k;
  int MaxIndex = 1;
  double CurrentScore;  // records the current log-likelihood
  int MaxVertex;  // this records the index of the largest vertex ratio found so far
  double MaxRatio;  // records the value of the ratio, actually it's the log of the ratio
  int MaxPriority; // records the community that the vertex wants to go to
  long int tempvertex = 0;
  
  double prevMaxScore = -numeric_limits<double>::max( );
  long double tolerance = 0.00000001; // this prevents loops due to numerical errors.
  
  double ProposalRatio;
  double value;
  int Priority;
  
  Initialize();
  
  // This returns the log of the initial score
  MaxScore = ComputeInitialScore();
  int iter = 0;
  int maxIter = 100;
  while(MaxScore >= prevMaxScore + tolerance && iter < maxIter)
  {
    iter++;
    // cout << "MAX SCORE IS: " << MaxScore << endl;
    // we start with everything equal to the best values
    CurrentScore = MaxScore;
    prevMaxScore = MaxScore;
    MaxIndex = -1;
    
    // ChangeSet records which vertices are able to move, in that they haven't
    // already moved during this KL step.  Update index will tell when the vertex
    // was chosen to move.
    for(i=0; i < Nodes; i++)
    {
      CurrentState[i] = BestState[i];
      ChangeSet[i] = i;
      UpdateIndex[i] = -1;
    }
    
    for(i=0; i < MaxComms; i++)
    {
      CurrentCommVertices[i] = BestCommVertices[i];
      CurrentCommStubs[i] = BestCommStubs[i];
      for(j=0; j < MaxComms; j++)
	CurrentEdgeMatrix[i][j] = BestEdgeMatrix[i][j];
    }
    
    // This loop moves each vertex once
    // Note that we DONT reinitialize changeset as this is unnecessary
    // This would make it a factor of 2 slower.
    for(i=0; i < Nodes; i++)
    {
      MaxVertex = 0;
      MaxRatio = -numeric_limits<double>::max( );
      MaxPriority = 0;
      // This loop selects which vertex to move
      for(j=0; j < Nodes-i; j++)
      {
	// get proposal and proposal ratio for ChangeSet[j]
	Priority = 0;
	ProposalRatio = -numeric_limits<double>::max( );
	// we first compute the neighbor set of the vertex, this is fixed
	// and the same for every change,
	// computing this first makes this more efficient
	// zero indicates run with current communities
	ComputeNeighborSet(ChangeSet[j], 0);
	
	// DanLarremore Modification:
	// We actually don't want to try all the comms, because some of them are forbidden.
	// We only get to choose from the entries of Comms[Types[j]].
	
	for (unsigned int q=0; q < Comms[Types[ChangeSet[j]]].size(); ++q)
	{
	  k = Comms[Types[ChangeSet[j]]][q];
	  // we compute the value of vertex ChangeSet[j] going to k
	  // we DONT allow a vertex to remain where it was
	  // This is essential to enforce so that it will go downhill and not be greedy
	  if(k != CurrentState[ChangeSet[j]])
	  {
	    value = ComputeProposal(ChangeSet[j], CurrentState[ChangeSet[j]], k);
	    if(value > ProposalRatio)
	    {
	      Priority = k;
	      ProposalRatio = value;
	    }
	  }
	}
	
	// check whether its higher than what you already have as the max KL move
	if(ProposalRatio > MaxRatio)
	{
	  MaxVertex = j;  // Note this is not the vertex j, but the vertex given by ChangeSet[j]
	  MaxRatio = ProposalRatio;
	  MaxPriority = Priority;
	}
      }
      // now we move it, first recording the current neighbors so that
      // we can update the matrices properly
      ComputeNeighborSet(ChangeSet[MaxVertex], 0);
      // This updates the matrices to represent the vertices new state
      UpdateMatrices(ChangeSet[MaxVertex], 0, CurrentState[ChangeSet[MaxVertex]], MaxPriority);
      CurrentState[ChangeSet[MaxVertex]] = MaxPriority;
      // we are using logs so we add the maxratio to the current score for the new score
      CurrentScore = CurrentScore + MaxRatio;
      UpdateIndex[ChangeSet[MaxVertex]] = i;
      // we switch it with the last element of changeset, removing it from further consideration
      // until we have moved the other vertices
      tempvertex = ChangeSet[MaxVertex];
      ChangeSet[MaxVertex] = ChangeSet[Nodes-i-1];
      ChangeSet[Nodes-i-1] = tempvertex;
      
      // now if the new state is better than the previous best state we record this
      // MaxIndex in combination with UpdateIndex
      // telling us where we are in the movement of vertices
      if(CurrentScore > MaxScore)
      {
	MaxScore = CurrentScore;
	MaxIndex = i;
      }
    }
    
    
    // now we update BestState if a change resulted in a higher maximum
    // by implementing all the changes found above
    
    // There is a potential for speeding this up here.
    if(MaxIndex != -1)
    {
      for(i=0; i < Nodes; i++)
      {
	// we only make the changes to beststate that happened before or equal to maxindex
	// no other vertex is updated
	// fortunately the update order is irrelevant to the final result so
	// we can just do it this way
	// Because we force all moves to be different, these updates are all a switch of community
	if(UpdateIndex[i] <= MaxIndex)
	{
	  // the option 1 does update corresponding to the best states / matrices
	  ComputeNeighborSet(i, 1);
	  UpdateMatrices(i, 1, BestState[i], CurrentState[i]); // 1 does best matrix update
	  BestState[i] = CurrentState[i];
	}
      }
    }
    
  }
  if (iter==maxIter)
  {
    cout << "...maxIterations on this KL run." << endl;
  }
  
  return;
}

// This starts off from a random initial condition
void Initialize()
{
  unsigned int i, j;
  int neighbor;
  int sum;
  
  for(i=0; i < MaxComms; i++)
  {
    BestCommVertices[i] = 0;
    BestCommStubs[i] = 0;
    for(j=0; j < MaxComms; j++)
    {
      BestEdgeMatrix[i][j] = 0;
    }
  }
  
  for(i=0; i < Nodes; i++)
  {
    // BestState[i] = int(numgen.nextDouble(MaxComms));   // REPLACERNG, should return 0 to MaxComms-1 in integer
    // DanLarremore Modification:
    // The initialized communities must be constrained to respect types.
    // cout << i << "," << Types[i] << endl;
    BestState[i] = Comms[Types[i]][0] + (random() % Comms[Types[i]].size() );
    
    if(InitializationOption == 1)
      BestState[i] = TrueState[i];
    BestCommVertices[BestState[i]]++;
    BestCommStubs[BestState[i]] += Degree[i];
  }
  
  // We are going to double count all edges and then divide two
  for(i=0; i < Nodes; i++)
  {
    for(j=0; j < Degree[i]; j++)
    {
      neighbor = AdjList[i][j];
      BestEdgeMatrix[BestState[i]][BestState[neighbor]]++;
      // the following statement prevents us from quadruple counting same comm edges.
      if(BestState[neighbor] != BestState[i])
	BestEdgeMatrix[BestState[neighbor]][BestState[i]]++;
    }
  }
  
  sum = 0;
  // we get rid of the double-counting
  for(i=0; i < MaxComms; i++)
  {
    for(j=0; j < MaxComms; j++)
    {
      BestEdgeMatrix[i][j] = BestEdgeMatrix[i][j]/2;
      if(i != j)
	sum = sum + BestEdgeMatrix[i][j];
      if(i == j)
	sum = sum + 2*BestEdgeMatrix[i][i];
    }
  }
  
  //    cout << "The starting best edge matrix encodes: " << sum << " twice edges." << endl;
  //    for(i=0; i < MaxComms; i++)
  //    {
    //        for(j=0; j < MaxComms; j++)
    //        {
      //            if(i==j)
      //                cout << 2*BestEdgeMatrix[i][j]/TwiceEdges << " ";
      //            if(i!=j)
      //                cout << BestEdgeMatrix[i][j]/TwiceEdges << " ";
      //        }
      //        cout << endl;
      //    }
      //
      return;
}

double ComputeInitialScore()
{
  // For the running of the KL algorithm itself this does not matter as all we use
  // are whether the score increases
  // We will want this when we compare different initializations
  
  // this actually returns 1/2 the unnormalized log-likelihood listed in the paper
  
  unsigned int i,j;
  double sum = 0;
  
  for(i=0; i < MaxComms; i++)
  {
    if(isDegreeCorrect)
    {
      sum = sum - LogFunction(BestCommStubs[i]);
    }
    if(!isDegreeCorrect)
    {
      if(BestCommVertices[i] != 0)
	sum = sum - double(BestCommStubs[i])*log(BestCommVertices[i]);
    }
    for(j=i; j < MaxComms; j++)
    {
      if(j != i)
      {
	sum = sum + LogFunction(BestEdgeMatrix[i][j]);
      }
      if(i==j)
      {
	sum = sum + .5*LogFunction(2*BestEdgeMatrix[i][j]);
      }
    }
  }
  
  return sum;
}

// We compute this using the current comm matrices
// We avoid the potential pitfalls of huge intermediate numbers by adding logs together.  So we treat 0 log 0 as 0.
// We return 0 for degree zero vertices (which really shouldn't be sent into the program
// in the first place.)
// We also return 0 for from = destination cause there is no change then.
// Here we use base e.  It returns the log of the actual value.
// Again this is half of the change in the unnormalized log-likelihood listed in the paper
double ComputeProposal(int vertex, int from, int destination)
{
  unsigned int i;//-Wunused, j, k;
  double ratiovalue = 0;
  int fromcount = 0;
  int destcount = 0;
  
  double help1;
  double help2;
  double help3;
  
  if(from == destination)
    return 0;
  
  // if the degree of the vertex is zero we know nothing about it
    // in this case we don't ever change its community
    // at the end we put all degree zeroes into their own group
    if(isDegreeCorrect)
    {
      if(Degree[vertex] == 0)
	return 0;
    }
    
    // we first add up all the cross-terms (between communities that are not from / destination)
    for(i=0; i < ActualDiffComms; i++)
    {
      // we lost NeighborSet[1][i] edges to NeighborSet[0][i] from the from comm
      // we gain the same amount in the destination comm
      // IFF the comms were not from and destination
      if((NeighborSet[0][i] != from) && (NeighborSet[0][i] != destination))
      {
	// do update NOTE: each community mcc' gets updated once if it had edges switch out
	// which is correct, remembering that mcc' is symmetric and we only count c < c' here
	
	help1 = double(CurrentEdgeMatrix[from][NeighborSet[0][i]]);
	help2 = double(CurrentEdgeMatrix[destination][NeighborSet[0][i]]);
	help3 = double(NeighborSet[1][i]);
	
	ratiovalue = ratiovalue + LogFunction(help1-help3) - LogFunction(help1);
	ratiovalue = ratiovalue + LogFunction(help2+help3) - LogFunction(help2);
      }
      
      if(NeighborSet[0][i] == from)
	fromcount = NeighborSet[1][i];
      
      if(NeighborSet[0][i] == destination)
	destcount = NeighborSet[1][i];
    }
    
    // now we add in the term corresponding to from / dest
    help1 = double(CurrentEdgeMatrix[from][destination]);
    help2 = double(fromcount-destcount);
    ratiovalue = ratiovalue + LogFunction(help1 + help2) - LogFunction(help1);
    
    // now we add in the terms corresponding to from
    if(isDegreeCorrect)
    {
      help1 = double(CurrentCommStubs[from]);
      help2 = double(Degree[vertex]);
      ratiovalue = ratiovalue - LogFunction(help1 - help2) + LogFunction(help1);
    }
    if(!isDegreeCorrect)
    {
      help1 = double(CurrentCommStubs[from]);
      help2 = double(Degree[vertex]);
      if(help1 - help2 != 0)
	ratiovalue = ratiovalue - (help1-help2)*log(double(CurrentCommVertices[from]-1));
      if(help1 != 0)
	ratiovalue = ratiovalue + help1*log(double(CurrentCommVertices[from]));
    }
    
    // now we do from/from
    help1 = double(2*CurrentEdgeMatrix[from][from]);
    help2 = double(2*SelfEdgeCounter + 2*fromcount);
    ratiovalue = ratiovalue + .5*LogFunction(help1 - help2) - .5*LogFunction(help1);
    
    // now we add in the terms corresponding to dest
    if(isDegreeCorrect)
    {
      help1 = double(CurrentCommStubs[destination]);
      help2 = double(Degree[vertex]);
      ratiovalue = ratiovalue - LogFunction(help1 + help2) + LogFunction(help1);
    }
    if(!isDegreeCorrect)
    {
      help1 = double(CurrentCommStubs[destination]);
      help2 = double(Degree[vertex]);
      if(help1 + help2 != 0)
	ratiovalue = ratiovalue - (help1+help2)*log(double(CurrentCommVertices[destination]+1));
      if(help1 != 0)
	ratiovalue = ratiovalue + help1*log(double(CurrentCommVertices[destination]));
    }
    
    // and now dest/dest
    help1 = double(2*CurrentEdgeMatrix[destination][destination]);
    help2 = double(2*SelfEdgeCounter + 2*destcount);
    ratiovalue = ratiovalue + .5*LogFunction(help1 + help2) - .5*LogFunction(help1);
    
    return ratiovalue;
}

void ComputeNeighborSet(int vertex, int option)
{
  unsigned int i;//-Wunused,j;
  int neighbor;
  
  SelfEdgeCounter = 0;
  
  for(i=0; i < MaxComms; i++)
  {
    TempNeighborSet[0][i] = i;
    TempNeighborSet[1][i] = 0;
    NeighborSet[0][i] = i;
    NeighborSet[1][i] = 0;
  }
  
  // NOTE SINCE A SELF-EDGE SHOWS UP TWICE IN THE ADJLIST THIS DOUBLE
  // COUNTS THESE EDGES, WE RECORD THE NUMBER OF TIMES THIS HAPPENS
  // IN A SEPARATE VARIABLE AND THEN DIVIDE BY TWO
  for(i=0; i < Degree[vertex]; i++)
  {
    neighbor = AdjList[vertex][i];
    if(neighbor != vertex)
    {
      if(option == 0)
	TempNeighborSet[1][CurrentState[neighbor]]++;
      
      if(option == 1)
	TempNeighborSet[1][BestState[neighbor]]++;
    }
    if(neighbor == vertex)
      SelfEdgeCounter++;
  }
  
  SelfEdgeCounter = SelfEdgeCounter/2;
  
  ActualDiffComms = 0;
  for(i=0; i < MaxComms; i++)
  {
    if(TempNeighborSet[1][i] != 0)
    {
      NeighborSet[0][ActualDiffComms] = TempNeighborSet[0][i];
      NeighborSet[1][ActualDiffComms] = TempNeighborSet[1][i];
      ActualDiffComms++;
    }
  }
  
  return;
}

void UpdateMatrices(int vertex, int option, int from, int destination)
{
  unsigned int i;//-Wunused,j;
  int fromcount = 0;
  int destcount = 0;
  
  if(option == 0)
  {
    CurrentCommVertices[from]--;
    CurrentCommVertices[destination]++;
    CurrentCommStubs[from] -= Degree[vertex];
    CurrentCommStubs[destination] += Degree[vertex];
    
    for(i=0; i < ActualDiffComms; i++)
    {
      if((NeighborSet[0][i] != from) && (NeighborSet[0][i] != destination))
      {
	// do update NOTE: each community mcc' gets updated once if it had edges switch out
	// which is correct, remembering that mcc' is symmetric and we only count c < c' here
	CurrentEdgeMatrix[from][NeighborSet[0][i]] -= NeighborSet[1][i];
	CurrentEdgeMatrix[NeighborSet[0][i]][from] -= NeighborSet[1][i];
	
	CurrentEdgeMatrix[destination][NeighborSet[0][i]] += NeighborSet[1][i];
	CurrentEdgeMatrix[NeighborSet[0][i]][destination] += NeighborSet[1][i];
      }
      
      if(NeighborSet[0][i] == from)
	fromcount = NeighborSet[1][i];
      
      if(NeighborSet[0][i] == destination)
	destcount = NeighborSet[1][i];
    }
    
    CurrentEdgeMatrix[from][from] -= (SelfEdgeCounter + fromcount);
    CurrentEdgeMatrix[destination][destination] += (SelfEdgeCounter + destcount);
    CurrentEdgeMatrix[from][destination] += (fromcount - destcount);
    CurrentEdgeMatrix[destination][from] += (fromcount - destcount);
  }
  
  if(option == 1)
  {
    BestCommVertices[from]--;
    BestCommVertices[destination]++;
    BestCommStubs[from] -= Degree[vertex];
    BestCommStubs[destination] += Degree[vertex];
    
    for(i=0; i < ActualDiffComms; i++)
    {
      if((NeighborSet[0][i] != from) && (NeighborSet[0][i] != destination))
      {
	// do update NOTE: each community mcc' gets updated once if it had edges switch out
	// which is correct, remembering that mcc' is symmetric and we only count c < c' here
	BestEdgeMatrix[from][NeighborSet[0][i]] -= NeighborSet[1][i];
	BestEdgeMatrix[NeighborSet[0][i]][from] -= NeighborSet[1][i];
	
	BestEdgeMatrix[destination][NeighborSet[0][i]] += NeighborSet[1][i];
	BestEdgeMatrix[NeighborSet[0][i]][destination] += NeighborSet[1][i];
      }
      
      if(NeighborSet[0][i] == from)
	fromcount = NeighborSet[1][i];
      
      if(NeighborSet[0][i] == destination)
	destcount = NeighborSet[1][i];
    }
    
    BestEdgeMatrix[from][from] -= (SelfEdgeCounter + fromcount);
    BestEdgeMatrix[destination][destination] += (SelfEdgeCounter + destcount);
    BestEdgeMatrix[from][destination] += (fromcount - destcount);
    BestEdgeMatrix[destination][from] += (fromcount - destcount);
  }
  
  return;
}

// This function returns zero if x = 0, otherwise it returns x*log(x)
double LogFunction(double x)
{
  if(x < 0)
  {
    cout << "x = " << x << endl;
    //mexErrMsgTxt("SOMETHING WRONG HAS OCCURRED...STOP! x is below zero");
  }
  
  if(x == 0)
  {
    return 0;
  }
  
  return x*log(x);
}

// We do not normalize VI here.
double ComputeVI()
{
  double EntropyA;
  double EntropyB;
  double EntropyAB;
  
  EntropyA = Entropy(0); // 0 calls for best state
  EntropyB = Entropy(1); // 1 calls for true state
  EntropyAB = JointEntropy(); // does joint for best / true
  
  return 2*EntropyAB-EntropyA-EntropyB;
}

double ComputeNMI()
{
  double EntropyA;
  double EntropyB;
  double EntropyAB;
  
  EntropyA = Entropy(0);
  EntropyB = Entropy(1);
  EntropyAB = JointEntropy();
  
  return 2*(EntropyA+EntropyB-EntropyAB)/(EntropyA+EntropyB);
}

double Entropy(int entoption)
{
  double Ent = 0;
  
  unsigned int i, j;//-Wunused, k;
  double *Ni;
  
  Ni = new double [MaxComms];
  
  for(i = 0; i < MaxComms; i++)
  {
    Ni[i] = 0;
  }
  
  for(j=0; j < Nodes; j++)
  {
    if(entoption == 0)
      Ni[BestState[j]]++;
    if(entoption == 1)
      Ni[TrueState[j]]++;
  }
  
  // NOTE WE RETURN THE ENTROPY IN LOG BASE 2
  for(i=0; i < MaxComms; i++)
  {
    if(Ni[i] != 0)
    {
      Ent = Ent - Ni[i]/double(Nodes)*log(Ni[i]/double(Nodes))/log(2);
    }
  }
  
  delete [] Ni;
  
  return Ent;
}

// Calculates the joint entropy
double JointEntropy()
{
  unsigned int i, j, l;
  double JointEnt = 0;
  
  double Nij[constComms][constComms];
  
  // This rapidly calculates Nij in a simple fashion.
  for(i=0; i < MaxComms; i++)
  {
    for(j=0; j < MaxComms; j++)
    {
      Nij[i][j] = 0;
    }
  }
  
  for(l=0; l < Nodes; l++)
  {
    Nij[BestState[l]][TrueState[l]]++;
  }
  
  JointEnt = 0;
  for(i=0; i < MaxComms; i++)
  {
    for(j = 0; j < MaxComms; j++)
    {
      if(Nij[i][j] != 0)
      {
	// divide by log 2 to convert to base 2.
	JointEnt = JointEnt - Nij[i][j]/double(Nodes)*log(Nij[i][j]/double(Nodes))/log(2);
      }
    }
  }
  
  return JointEnt;
}

/***** CPP ONLY *****/
void GetTheVertexTypes(string fileName)
{
    ifstream InputFile;
    string lineread;
    char *buffer;
    int entry = 0;
    
    InputFile.open(fileName.c_str());
    if (!InputFile)
    {
        cout << "Error in opening file";
        printf("\n\n%s\n\n",fileName.c_str());
        exit(0);
        cin.get();
        return;
    }
    
    while(std::getline(InputFile, lineread)) // Read line by line
    {
        buffer = new char [lineread.size()+1];
        strcpy(buffer, lineread.c_str());
        sscanf(buffer, "%d", &entry);
        Types.push_back(entry-1);
        delete[] buffer;
    }
    InputFile.close();
    
    int maxTypes = 0;
    for (unsigned int i=0; i<Types.size(); ++i)
    {
        maxTypes = max(maxTypes,Types[i]);
    }
    vector<int> tally (maxTypes+1,0);
    for (unsigned int i=0; i<Types.size(); ++i)
    {
        tally[Types[i]]++;
    }
    for (unsigned int i=0; i<tally.size(); ++i)
    {
        cout << "TYPE " << i << " : " << tally[i] << endl;
    }
    Nodes = max(Nodes,long(Types.size()));
}
// NOTE: OVERLOADED. R uses GetTheNetworkEdges(double,size_t)
void GetTheNetworkEdges(string fileName)
{
    long int i;//-Wunused,j;
    ifstream InputFile;
    ifstream InputFile2;
    // Modify here for a different file name.
    //string fileName = "matrix.txt";
    // This string is unused, but left declared.
    string fileName2 = "polBlogsSymmLargestFComms.txt";
    string lineread;
    char *buffer;
    long int entry1 = 0;
    long int entry2 = 0;
    int counter = 0;
    int ignore = 0;
    string ignore3;
    
    InputFile.open(fileName.c_str());
    if (!InputFile)
    {
        cout << "Error in opening file";
        printf("\n\n%s\n\n",fileName.c_str());
        exit(0);
        cin.get();
        return;
    }
    
    while(std::getline(InputFile, lineread)) // Read line by line
    {
        buffer = new char [lineread.size()+1];
        strcpy(buffer, lineread.c_str());
        sscanf(buffer, "%ld %ld %d", &entry1, &entry2, &ignore);
        //  sscanf(buffer, "n%d,n%d,%f,%d", &entry1, &entry2, &ignore2, &ignore);
        // cout << entry1 << " " << entry2 << " " << ignore2 << " " << ignore << endl;
        // cin.get();
        EdgeList[counter][0] = entry1-1;
        EdgeList[counter][1] = entry2-1;
        // cout << entry1 << " " << entry2 << endl;
        // cin.get();
        counter = counter+1;
        Nodes = max(Nodes,entry1);
        Nodes = max(Nodes,entry2);
        delete[] buffer;
    }
    InputFile.close();
    
    printf("NODES: %li\n",Nodes);
    
    TwiceEdges = 2*counter;
    
    
    // We start the degree values and LastEmpty all at zero
    for(i=0; i < Nodes; i++)
    {
        Degree[i] = 0;
        LastEmpty[i] = 0;
    }
    // First we count the degrees by scanning through the list once
    for(i=0; i < counter; i++)
    {
        Degree[EdgeList[i][0]]++;
        Degree[EdgeList[i][1]]++;
    }
    // Now we make space in the adjacency lists
    for(i=0; i < Nodes; i++)
    {
        AdjList[i] = new long int [Degree[i]];
    }
    // Now we read the edges into the adjacency lists utilizing lastempty to put them into
    // the proper spots
    for(i=0; i < counter; i++)
    {
        AdjList[EdgeList[i][0]][LastEmpty[EdgeList[i][0]]] = EdgeList[i][1];
        LastEmpty[EdgeList[i][0]]++;
        
        AdjList[EdgeList[i][1]][LastEmpty[EdgeList[i][1]]] = EdgeList[i][0];
        LastEmpty[EdgeList[i][1]]++;
    }
    cout << "Read in twice edges = " << TwiceEdges << endl;
    
    return;
}
void PrintResults(string saveLocation)
{
    chdir(saveLocation.c_str());
    string saveName;
    string saveName2;
    
    int i;//-Wunused, j, k;
    ofstream myfile;
    
    if (isDegreeCorrect)
    {
        saveName = "biDCSBMcomms1.tsv";
        saveName2 = "biDCSBMcomms1.score";
    }
    else
    {
        saveName = "biSBMcomms1.tsv";
        saveName2 = "biSBMcomms1.score";
    }
    int ctr = 1;
    while (access(saveName.c_str(),F_OK)==0)
    {
        ctr++;
        char c[20];
        sprintf(c, "%d", ctr);
        string intStr;
        intStr = string(c);
        if (isDegreeCorrect)
        {
            saveName = std::string("biDCSBMcomms") + c + std::string(".tsv");
            saveName2 = std::string("biDCSBMcomms") + c + std::string(".score");
        }
        else
        {
            saveName = std::string("biSBMcomms") + c + std::string(".tsv");
            saveName2 = std::string("biSBMcomms") + c + std::string(".score");
        }
    }
    
    myfile.open(saveName.c_str());
    for(i=0; i < Nodes; i++)
    {
        myfile << BestState[i] << endl;
    }
    myfile.close();
    
    myfile.open(saveName2.c_str());
    myfile << ComputeInitialScore() << endl;
    myfile.close();
    
    return;
}
