#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <ctime>
#include <stdlib.h>
#include <assert.h>
#include <bitset>
#include <map>
#include <climits>
#include <sstream>
#include <iomanip>
#include <set>
#include <queue>

#include "Md5.h"
#include "Profiler.h"

using namespace std;

typedef unsigned int uint;

//////////////////////// Constants //////////////////////////////////
// Construction of hash functions.
const uint LARGENUMBER = 20000000;
const uint primeMax = 20000000;
const uint primeMin = 10000000;
const uint domain = 15773;

// Parameter for random generation function.
int seed = 123;

// Files for min-wise hash functions and superHashTemplates.
const string hashFunctionsFile = "hashFunctions_100";
const string superHashTemplatesFile = "superHashTemplates_";

// Parameters for superHash construction.
// Assumption: numThresholds must be the same with the size of the settings.
const uint numThresholds = 5;
const double thresholds[ 5 ] = { 0.98, 0.90, 0.80, 0.70 , 0.6 };
const uint numSuperhashForTheta[ 5 ] = { 2, 9, 20, 25, 30 };
const uint numMinhashPerSuperhashForTheta[ 5 ] = { 40, 21, 14, 9, 6 };

// Similarity measure during the Filtering (of candidates) step.
const uint similarityMeasure = 0; // 0: Jaccard, 1: intersection.

// Parameters for number of candidates.
const uint pairsThreshold = 100;
const uint maxCandidates = 300; // before filtering.
const uint numCandidates = 200; // after filtering.

////////////////////// Class declarations /////////////////////////
class Edge {
  public:
   uint _docID;
   double _similarity;
   // Constructors.
   Edge() {
     _similarity = 0.0f;
     _docID = 0;
   }
   Edge(const uint inDocID, const double inSimilarity) {
     _docID = inDocID;
     _similarity = inSimilarity;
   }
};

class edgeComparator {
  public:
    bool operator() (const Edge& a, const Edge& b) const {
      return (a._similarity > b._similarity);
    }  
};

class edgeComparatorLargestFirst {
  public:
    bool operator() (const Edge& a, const Edge& b) const {
      return (a._similarity < b._similarity);
    }  
};

class Node {
  public:
   // Candidate edges before filtering (by computing full score).
   vector<Edge> _candidates;

   // Interface for adding/removing edges in Node.
   // Largest/(Smallest) similarity edge at the top.
   priority_queue<Edge, vector<Edge>, edgeComparatorLargestFirst> _edgeHeap;  

   // Add edge in the priority queue.
   void addEdge(const Edge& edge) {
     _edgeHeap.push(edge); // push edge
   }   

   // Constructor.
   Node() {}
};

class minHash {
  public:
   uint _docid;
   uint _numTerms;
   vector<unsigned short> _minH;
   minHash() {};
   minHash(const uint inDocid, const uint inNumTerms, vector<unsigned short> inMinH) {
     _docid = inDocid;
     _numTerms = inNumTerms; 
     _minH = inMinH;
   }
};

class Graph {
  public: 
    vector<Node> _nodes;
    Graph(const uint sz) { _nodes.resize(sz, Node()); 
    //TOGO
    ones = 0;
    pairs = 0;
    samples = 0;
}
    vector<uint> GreedyNearestNeighborTSP(const uint baseDocID);
    void FilterCandidates(const uint typeOfSimilarity, vector<minHash>& minHashes);

    // helper TOGO
    profiler one, pair, sample;
    uint ones, pairs, samples;
};

// Version with MD5().
class superHash {
  public:
   uint _docID;
   string _superHash;
   minHash* _minH;
   superHash(string inSuperHash, minHash* inMinHashPtr, uint inDocID) {
     _superHash = inSuperHash;
     _minH = inMinHashPtr;
     _docID = inDocID;
   }
};

bool superHashComparator (const superHash i, const superHash j) { 
  return (i._superHash < j._superHash); 
}

////////////////////// Method Signatures /////////////////////////
vector<uint> getPrimesUpTo(const uint max);
void keepPrimesMoreThan(const uint primeMIN,
                        vector<uint>& primes);
void randomlySelectKNumbers(const uint k,
                            vector<uint>&primes);
void printVec(const vector<uint> vec);
void printUshortVec(const vector<unsigned short> vec);
unsigned short hashF(const uint termId,
                     const uint ithHashFunction,
                     const uint domain,
                     const vector<uint>& hashFunctions);
void writeOut(const vector<uint> allPrimes,
              string hashFunctionsFile);
vector<uint> loadHashFunctions(string hashFunctionsFile);
minHash createMinHashes(const uint docID,
                        const vector<uint>& tids,
                        const vector<uint>& hashFunctions);
vector<vector<uint> > createDocs(const uint numDocs,
                                 const uint minDocSz,
                                 const uint maxDocSz,
                                 const uint maxTid);
vector<uint> createDoc(const uint minDocSz,
                       const uint maxDocSz,
                       const uint maxTid);
uint randInRange(const uint min,
                 const uint max);
void printDocs(const vector<vector<uint> > vec);
vector<unsigned short> TokenizeStringToUShort(string str);
unsigned short StrToUShort(const string number_in_string);
vector<uint> randomlySelectKNumbersUpTo(const uint k,
                                        const uint max);
void writeoutSuperHashTemplates(const uint superHashSize,
                                const uint numSuperHashes,
                                const uint numHashFunctions,
                                const string file);
vector<vector<unsigned short> > loadSuperHashTemplates(const string file);
Graph createSparseGraph(vector<minHash>& minHashes,
                        const uint baseDocID);
void printBin(unsigned long i);
vector<minHash> createMinHashesForKDocs(const uint k,
                                        vector<uint>& hashFuncs);
double computeJaccardSimilarity(minHash* a, 
                                minHash* b);
double computeIntersection(minHash* a,
                           minHash*b );
void selectNearestNeighborCandidates(const vector<superHash>& superHashVec,
                                     Graph& sparseGraph, 
                                     const double score,
                                     set<uint>& idxOfNodesLeft);
void findRangeWithSameSuperHash(const vector<superHash>& superHashVec,
                                const uint start, uint& end);
void addEdgesToSimilarDocs(const vector<superHash>& superHashVec,
                           const uint startIdx,
                           const uint endIdx,
                           const double score,
                           Graph& sparseGraph,
                           set<uint>& idxOfNodesLeft);
Graph constructArtificialGraph(const uint numNodes);
void selectKRandom(vector<uint>& randNodes,
                   int numElementsToSample);
string UshortToStr(const unsigned short);
superHash createSuperHash_MD5(minHash& minHashes, 
                              const vector<unsigned short>& superHashTemplate,
                              const uint baseDocID);
void createSuperhashTemplatesGivenSetup(const uint numMinhash,
                                        const string basefile);
uint randInRangeMS(const uint min, const uint max, int* seed);

////////////////////// Methods' implementation /////////////////////////

// For each node of the sparse graph, compute the nodes (or edges) with the
// highest similarity by computing the similarity (given as argument), which 
// is based on the min-hashes (which are also given).
// Assumption: each node of the sparse graph has _candidates that correspond
// to some expected score value based on the setup of the superhashes.
// Assumption: the similarityMeasure is defined in constants!
void Graph::FilterCandidates(const uint typeOfSimilarity, vector<minHash>& minHashes) {

  // TOGO STATS
  unsigned long long uniqueCand = 0;
  unsigned long long cntCand = 0;
  for (uint i = 0; i < _nodes.size(); ++i) {
    cntCand += _nodes[i]._candidates.size();
  }
  double avgCand = ((double) cntCand / (double) _nodes.size());
  cout << "______ Avg Cand number: " << setprecision(20) << avgCand << " for #nodes: " << _nodes.size() << " and total candidates: " << cntCand << endl;
  // ENED OF TODO

  // For each node of the sparse graph, if candidate (edges) vector not empty,
  // calculate the score based on the given similarity.
  for (uint i = 0; i < _nodes.size(); ++i) {
    if (_nodes[i]._candidates.size() == 0) continue;  // corner case for nodes without candidates.
    // Map used for deduplication.
    map<uint, bool> uniqueDest;
    map<uint, bool>::iterator it;
    
    // Vector maintaining edges with final score.
    vector<Edge> tmpEdges; 

    // For each candidate (edge or doc pair) of the given node, compute score.
    for (uint j = 0; j < _nodes[i]._candidates.size(); ++j) {
      uint destDocID = _nodes[i]._candidates[j]._docID;

      // Deduplication.
      it = uniqueDest.find(destDocID);
      if (it != uniqueDest.end()) continue; // skip node, since we have already selected it.
      uniqueDest[ destDocID ] = true;

      // Compute score.
      if (typeOfSimilarity == similarityMeasure) {  // 0: Jaccard, 1: intersection.
        _nodes[i]._candidates[j]._similarity = computeJaccardSimilarity(&minHashes[i], &minHashes[ destDocID ]);
      } else { // Intersection is estimated: Jaccard*(|A| + |B|) / (1 + Jacccard), where |A| #terms.
        uint iSz = minHashes[i]._numTerms;
        uint jSz = minHashes[destDocID]._numTerms;
        double tmpScore = _nodes[i]._candidates[j]._similarity;
        _nodes[i]._candidates[j]._similarity = tmpScore*((double) iSz + (double) jSz) / (1.0f + tmpScore);  
      }
      // Add final candidate to vector.
      tmpEdges.push_back( Edge(destDocID, _nodes[i]._candidates[j]._similarity) );
    }
  
    // Sort final candidates by decreasing score and keep the numCandidates ones.
    sort(tmpEdges.begin(), tmpEdges.end(), edgeComparator());

    // For each candidate, push it into the priority queue of the corresponding node.
    uint candidatesSize = min((uint)tmpEdges.size(), numCandidates); // pick the candidateSize.
    for (uint k = 0; k < candidatesSize; ++k) {
      _nodes[i].addEdge( tmpEdges[k] );
    }

// TOGO
uniqueCand += _nodes[i]._edgeHeap.size();

    // Erase candidates for this node.
    _nodes[i]._candidates.clear();
  }

// TOGO
  double avgUniqCand = ((double) uniqueCand / (double) _nodes.size());
  cout << "___AFTER___ Avg Cand number: " << setprecision(20) << avgUniqCand << " for #nodes: " << _nodes.size() << " and total candidates: " << uniqueCand << endl;
// END OF TOGO
}

// Given the sparseGraph, run TSP to obtain an ordering.
// 1. Start with the largest similarity candidate.
// 2. Add the partial path by the edge picked by step 1.
// 3. Assuming the current node has more candidates, visit them
//    exhaustively by highest similarity order and add the 
//    corresponding path appropriately.
// 4. When the candidates of a node are exhausted, run a heuristic
//    that selects the node that has the highest similarity candidate
//    as the node to proceed.
// 5. In case there are only isolated nodes left (i.e., no candidates),
//    exhaustively pick a node from the remaining nodes.
// The TSP method implements the approach proposed in paper:
// "Inverted file compression through document identifier reassignment"
//  by W. Shieh, T. Chen, J. Shann and C. Chung.
// Note that we provide as argument the baseDocID, so we can create the mapping
// in range [x, y), since so far we have been working in range [0, y - x).
vector<uint> Graph::GreedyNearestNeighborTSP(const uint baseDocID) {

  profiler p1;  // TOGO
  p1.StartPause(); // TOGO
  p1.Continue(); // TOGO
  // Vector to store the docid assignment (oldDocid -> newDocid).
  vector<uint> docidAssignment(_nodes.size(), 0);
  p1.PauseMicroseconds(); // tOGO

  profiler p2;  // TOGO
  p2.StartPause(); // TOGO
  p2.Continue(); // TOGO

// ADDON
  set<uint> remainingNodes;
  set<uint>::iterator it;
// END OF ADDON

  // Select the edge with the highest similarity in the sparseGraph.
  Edge largestSimilarityEdge(0, 0.0f);
  uint source = 0;  // will be used as the source of the current edge.
  for (uint i = 0; i < _nodes.size(); ++i) {
// ADDON
//    remainingNodes.insert(i);
// end of ADDON

    if (!_nodes[i]._edgeHeap.empty() &&
       (_nodes[i]._edgeHeap.top()._similarity > largestSimilarityEdge._similarity)) {
      largestSimilarityEdge._similarity = _nodes[i]._edgeHeap.top()._similarity;
      largestSimilarityEdge._docID = _nodes[i]._edgeHeap.top()._docID;
      source = i;
    }
  }

  p2.PauseMicroseconds(); // tOGO
  profiler p34; //tOGO
  profiler p35;  // TOGO
  p34.StartPause(); // TOGO
  p35.StartPause(); // TOGO
  // Start with a partial path with a vertex/node from the largestSimilarityEdge.
  docidAssignment[0] = source;

  vector<bool> visited(_nodes.size(), false);
  visited[ source ] = true;
  uint numNodesSelected = 1; 
  uint dest = largestSimilarityEdge._docID;

// ADDON
//cout << "size of remaining " << remainingNodes.size() << endl;
//   remainingNodes.erase(source);
//cout << "after size of remaining " << remainingNodes.size() << endl;
// end of addon

  profiler p3;  // TOGO
  p3.StartPause(); // TOGO
  p3.Continue(); // TOGO

  // While there are still unvisited nodes continue.
  while (numNodesSelected < _nodes.size()) {
    // select the edge from the current node with the maximal similarity.
    bool foundNextNode = false;
    uint curNode = 0;
    while (!_nodes[ source ]._edgeHeap.empty()) {
      curNode = _nodes[ source ]._edgeHeap.top()._docID;
      if (visited[ curNode ]) { // if already traversed the node, remove edge.
        _nodes[ source ]._edgeHeap.pop();
      } else {  // next node to visit found and stored in curNode.
        foundNextNode = true;  
        break;
      }
    }

  p34.Continue(); //TOGO

    // if not found, select node with the largest remaining similarity.
    if (!foundNextNode) {
      largestSimilarityEdge._similarity = 0.0f;  // re-use edge.
      largestSimilarityEdge._docID = 0;
      bool foundNextNodeHeuristic = false;
///*
      for (uint i = 0; i < _nodes.size(); ++i) { 
        if (!visited[ i ] &&
            !_nodes[i]._edgeHeap.empty() &&
            (_nodes[i]._edgeHeap.top()._similarity > largestSimilarityEdge._similarity)) {
          largestSimilarityEdge._similarity = _nodes[i]._edgeHeap.top()._similarity;
          largestSimilarityEdge._docID = i;
          curNode = i;
          foundNextNodeHeuristic = true;
        }
      }
//*/

/*
// ADDON
      for (it = remainingNodes.begin(); it != remainingNodes.end(); ++it) { 
        uint cnt = *it;
        if (!visited[ cnt ] &&
            !_nodes[cnt]._edgeHeap.empty() &&
            (_nodes[cnt]._edgeHeap.top()._similarity > largestSimilarityEdge._similarity)) {
          largestSimilarityEdge._similarity = _nodes[cnt]._edgeHeap.top()._similarity;
          largestSimilarityEdge._docID = cnt;
          curNode = cnt;
          foundNextNodeHeuristic = true;
// tOGO
//cout << " -- curNode in next largest sim: " << cnt << endl;
        }
      }
// END OF ADDON
*/
  p34.PauseMicroseconds(); // TOGO
  p35.Continue(); // TOGO
///*
      // In case that there are no nodes with edges left.
      if (!foundNextNodeHeuristic) {
        curNode = 0;
        while (curNode < _nodes.size() && visited[ curNode ]) {
          ++curNode;
        }
      }
//*/
/*
// ADDON
      // In case that there are no nodes with edges left.
      if (!foundNextNodeHeuristic) {
        curNode = 0;
        for (it = remainingNodes.begin(); it != remainingNodes.end(); ++it) { 
          if (!visited[*it]) {
            curNode = *it;
//togo
//cout << " --- found isolated node: " << curNode << endl;
            break;
          }
        }
      }
// ADDON
////*/
  p35.PauseMicroseconds(); // TOGO 
    }
    // The next node has already been selected.
    visited[ curNode ] = true;
    docidAssignment[ numNodesSelected ] = curNode;
    ++numNodesSelected;  
    source = curNode;

// ADDON
//cout << "before removing sourcenode: " << source << " size: " << remainingNodes.size() << endl;
//   remainingNodes.erase(source);
//cout << "removing: " << source << " and newsize: " << remainingNodes.size() << endl;
// end of addon

  }

  p3.PauseMicroseconds(); // TOGO

  profiler p4;  // TOGO
  p4.StartPause(); // TOGO
  p4.Continue(); // TOGO
  // Reconstruct the actual docids in the mapping by adding the baseDocID.
  for (uint i = 0; i < docidAssignment.size(); ++i) {
    docidAssignment[i] += baseDocID;
  }

  p4.PauseMicroseconds(); //TOGO

  cout << setprecision(20) << "P1 init: " << p1.GetTime() << endl;
  cout << setprecision(20) << "P2 find largest edge: " << p2.GetTime() << endl;
  cout << setprecision(20) << "P3 core: " << p3.GetTime() << endl;
  cout << setprecision(20) << "P34 : " << p34.GetTime() << endl;
  cout << setprecision(20) << "P35 no nodes with edges left: " << p35.GetTime() << endl;
  cout << setprecision(20) << "P4 reconstruction: " << p4.GetTime() << endl;

  return docidAssignment;
}

// Given the minHashes, add similar pairs based on the superhashes
// that are constructed based on the setup in constants (corresponding)
// to specific threhold of similarity and return the sparse graph.
Graph createSparseGraph(vector<minHash>& minHashes, const uint baseDocID) {
  assert(minHashes.size() > 0);
 
  // Initialize a graph representation.
  Graph sparseGraph(minHashes.size());

  // TOGO
  sparseGraph.one.StartPause();
  sparseGraph.pair.StartPause();
  sparseGraph.sample.StartPause();
  // ENDOFTOGO

  // Maintain the nodes-documents that do no have enough candidates (at first all).
  // These nodes are going to participate in rounds until enough neighbors
  // have been discovered. In each round there is a threshold of similarity
  // (when constructing superhashes) is reduced.
  set<uint> idxOfNodesLeft;
  set<uint>::iterator it;
  for (uint i = 0; i < minHashes.size(); ++i) {
    idxOfNodesLeft.insert(i);
  }

  // Profilers.
  profiler profSHCreation;
  profiler profSorting;
  profiler profSelectNearest;
  profiler profFiltering;
  // Start Profilers.
  profFiltering.StartPause();
  profSHCreation.StartPause();
  profSorting.StartPause();
  profSelectNearest.StartPause();

  // Foreach threshold, (a) load the corresponding superHash Templates, 
  // (b) create the superHashes given the setup in constants, (c) sort 
  // superHashes and (d) add similar pairs.
  for (unsigned short j = 0; j < numThresholds; ++j) {
    // Load the corresponding superHashTemplates.
    const string input = superHashTemplatesFile + UshortToStr(j);
    vector<vector<unsigned short> > superHashTemplates = loadSuperHashTemplates(input);
 
    // Current threshold.
    const double curTheta = thresholds[ j ];
    cout << "Theta: " << curTheta << endl;

profiler p1; // TOGO
p1.StartPause(); // TOGO
p1.Continue();  //TOGO

    // Foreach super hash template.
    uint numSuperHashes = superHashTemplates.size();
    for (uint i = 0; i < numSuperHashes; ++i) {
      vector<superHash> superHashVec;  // maintain super hashes of this round.

profSHCreation.Continue();
      // Foreach doc that does not have enough neighbors discovered (candidates), 
      // create the super hash according to the i-th super hash template.
      set<uint>::iterator it = idxOfNodesLeft.begin();
      for (it = idxOfNodesLeft.begin(); it != idxOfNodesLeft.end(); ++it) {
        superHashVec.push_back(createSuperHash_MD5(minHashes[*it], superHashTemplates[i], baseDocID));
      }
profSHCreation.PauseMicroseconds();

profSorting.Continue();
      // Sort superHashes by superHash value.
      sort(superHashVec.begin(), superHashVec.end(), superHashComparator);
profSorting.PauseMicroseconds();

      cout << "# SH: " << i << endl; //TOGO

profSelectNearest.Continue();
      // Select nearest neighbor candidates for the sparse graph.
      selectNearestNeighborCandidates(superHashVec, sparseGraph, curTheta, idxOfNodesLeft);
profSelectNearest.PauseMicroseconds();
    }
p1.PauseMicroseconds(); // TOGO
//cout << setprecision(20) << "time: " << p1.GetTime() << endl; // TOGO
  }

profFiltering.Continue();
  // Filter nearest neighbors (candidates) by computing the similariy (given) from the min-hashes.
  sparseGraph.FilterCandidates(similarityMeasure, minHashes);
profFiltering.StartPause();

cout << "Ones: " << sparseGraph.ones << endl;
cout << "Pairs: " << sparseGraph.pairs << endl;
cout << "Samples: " << sparseGraph.samples << endl;

cout << setprecision(20);
cout << "Time for SH creation: " << profSHCreation.GetTime() << endl;
cout << "Time for SH sorting: " << profSorting.GetTime() << endl;
cout << "Time for candidate selection: " << profSelectNearest.GetTime() << endl;
cout << "Time for filtering candidates: " << profFiltering.GetTime() << endl;

// run tSP
profiler profTSP;
profTSP.StartPause();
profTSP.Continue();
vector<uint> docidAssignment = sparseGraph.GreedyNearestNeighborTSP(baseDocID);
profTSP.PauseMicroseconds();
cout << "Time for TSP: " << profTSP.GetTime() << endl;

cout << "---------------------" << endl;
cout << "Time for singles: " << sparseGraph.one.GetTime() << endl;
cout << "Time for pairs: " << sparseGraph.pair.GetTime() << endl;
cout << "Time for sampling: " << sparseGraph.sample.GetTime() << endl;

  return sparseGraph;
}

// Foreach threhold - corresponding superhashes, select the nearest neighbors of similar documents
// by finding blocks of same superhashes (after the superhashes have been sorted) and add these
// neighbors to the sparse graph.
// Assumption: superHash vector of size more than 0.
// Given the superhashes for the nodes-docs that do not have enough neighbors (idxOfNodesLeft),
// the sparse graph (so we can add edges of candidates), the current score based on the superhash
// construction setting and the nodes without enough candidates, add candidates by finding nodes
// with the same superhash values.
void selectNearestNeighborCandidates(const vector<superHash>& superHashVec,
                                     Graph& sparseGraph,
                                     const double score,
                                     set<uint>& idxOfNodesLeft) {
  assert(superHashVec.size() > 0);
  uint start = 0;
  uint end = 0;
  uint last = superHashVec.size();
  vector<superHash>::const_iterator low = superHashVec.begin(), up;

  while (low != superHashVec.end()) {
    start = low - superHashVec.begin();
    up = upper_bound(low, superHashVec.end(), superHashVec[ start ], superHashComparator);
    end = up - 1 - superHashVec.begin();
    low = up;

    // Add edges to document-pairs with similarity == score, only if more than 2 docs have same superhash.
    if (end - start > 0) {
      addEdgesToSimilarDocs(superHashVec, start, end, score, sparseGraph, idxOfNodesLeft);
    }
  }
/*
  // Older version for finding range with the same superhashes.
  uint curStartIdx = 0;
  uint curEndIdx = 0;
  uint lastIdx = superHashVec.size();

  while (curStartIdx != lastIdx) {
    // Find docID ranges with same superHashes.
    findRangeWithSameSuperHash(superHashVec, curStartIdx, curEndIdx);

    // Add edges to document-pairs with similarity == score.
    addEdgesToSimilarDocs(superHashVec, curStartIdx, curEndIdx, score, sparseGraph, idxOfNodesLeft);

    ++curEndIdx;  // next start idx.
    curStartIdx = curEndIdx;  // set curStartIdx.
    assert(curEndIdx == curStartIdx);
  } */
}

// Given the superHashVec, the start and end idx of docs with the same superhash, the sparsegraph
// and the nodes participating in the round (idxOfNodesLeft), (i) do nothing if there is only one docid,
// (ii) add all document pairs as candidates (iii) add a sample of doc pairs as candidates. 
// TOCHECK!! (ii and iii).
void addEdgesToSimilarDocs(const vector<superHash>& superHashVec,
                           const uint startIdx,
                           const uint endIdx,
                           const double score,
                           Graph& sparseGraph,
                           set<uint>& idxOfNodesLeft) {
  uint numDocsWithSameSuperHash = endIdx - startIdx;
  assert(numDocsWithSameSuperHash > 0);

//TOGO
//cout << "range: " << numDocsWithSameSuperHash << " starid: " << startIdx << " - " << endIdx << endl;
 
  if (numDocsWithSameSuperHash == 1) {  // two documents with the same superhash.
++sparseGraph.ones; // TOGO
sparseGraph.one.Continue(); // TOGO
    uint sourceDocid = superHashVec[ startIdx ]._docID;
    // In case the current Node has already enough candidates, we skip the process.
    if (sparseGraph._nodes[sourceDocid]._candidates.size() > maxCandidates) {
      idxOfNodesLeft.erase(sourceDocid);
      return;
    }
    uint destinationDocid = superHashVec[ endIdx ]._docID;
    Edge edgeSource(destinationDocid, score);
    sparseGraph._nodes[sourceDocid]._candidates.push_back(edgeSource);

sparseGraph.one.PauseMicroseconds(); // TOGO
    // Add edge with the opposite direction.
    //Edge edgeDestination(sourceDocid, score);
    //sparseGraph._nodes[destinationDocid]._candidates.push_back(edgeDestination);
  } else if (numDocsWithSameSuperHash < pairsThreshold) {

/* sanity check
uint cnt = 0;
for (uint i = 0; i < superHashVec[startIdx]._minH->_minH.size(); ++i) {
  unsigned short sH1 = superHashVec[startIdx]._minH->_minH[i];
  unsigned short sH2 = superHashVec[startIdx+1]._minH->_minH[i];
  cout << i << " minh[0]: " << sH1 << " minh[1]: " << sH2 << endl;
  if (sH1 == sH2) ++cnt;
}
cout << cnt << " out of :" << superHashVec[startIdx]._minH->_minH.size() << endl;
exit(0);
*/

++sparseGraph.pairs; // TOGO
sparseGraph.pair.Continue(); // TOGO
    for (uint i = startIdx; i < endIdx; ++i) { 
      uint sourceDocid = superHashVec[i]._docID;
      // In case the current Node has already enough candidates, we skip the process.
      if (sparseGraph._nodes[sourceDocid]._candidates.size() > maxCandidates) {
        idxOfNodesLeft.erase(sourceDocid);
        continue;
      }

      for (uint j = startIdx; j < endIdx; ++j) {
        if (i == j) continue; 
        uint destinationDocid = superHashVec[j]._docID;
        Edge edgeSource(destinationDocid, score);
        sparseGraph._nodes[sourceDocid]._candidates.push_back(edgeSource); 

        // Add edge with the opposite direction.
        //Edge edgeDestination(sourceDocid, score);
        //sparseGraph._nodes[destinationDocid]._candidates.push_back(edgeDestination);

        // In case the current Node has already enough candidates, we skip the process.
        // Note that we also need another check inside the second for loop.
        if (sparseGraph._nodes[sourceDocid]._candidates.size() > maxCandidates) {
          idxOfNodesLeft.erase(sourceDocid);
          break;
        }
      }
    }  

sparseGraph.pair.PauseMicroseconds(); // TOGO
  } else {  // compute sample of scores for each node/document.
++sparseGraph.samples; //TOGO
sparseGraph.sample.Continue(); // TOGO
    assert(numDocsWithSameSuperHash >= pairsThreshold);
    for (uint i = startIdx; i < endIdx; ++i) {
      uint sourceDocid = superHashVec[i]._docID;
      // In case the current Node has already enough candidates, we skip the process.
      if (sparseGraph._nodes[sourceDocid]._candidates.size() > maxCandidates) {
        idxOfNodesLeft.erase(sourceDocid);
        continue;
      }
  
      // Select K random nodes.
      vector<uint> randNodes(numDocsWithSameSuperHash, startIdx);
      for (uint cnt = 0; cnt < randNodes.size(); ++cnt) {
        randNodes[ cnt ] += cnt;
        assert(randNodes[cnt] >= startIdx);
        assert(randNodes[cnt] < endIdx);
      }
       
      uint idxOfSameElement = i - startIdx;
      assert(idxOfSameElement < randNodes.size());
      assert(randNodes[ idxOfSameElement ] == i);
      randNodes.erase(randNodes.begin() + idxOfSameElement); // remove the element same with i.

      // Compute how many edges to add to the node.
      int numEdgesToAdd = min((int) numDocsWithSameSuperHash - 1, abs((int) (maxCandidates - sparseGraph._nodes[sourceDocid]._candidates.size()) ) );
      selectKRandom(randNodes, numEdgesToAdd); //abs((int) (maxCandidates - sparseGraph._nodes[sourceDocid]._candidates.size()) ));
      uint randNode = 0; 

      for (uint j = 0; j < randNodes.size(); ++j) {
        randNode = randNodes[ j ];       
        uint destinationDocid = superHashVec[randNode]._docID;
        Edge edgeSource(destinationDocid, score);
        sparseGraph._nodes[sourceDocid]._candidates.push_back(edgeSource);

        // Add edge with the opposite direction.
        //Edge edgeDestination(sourceDocid, score);
        //sparseGraph._nodes[destinationDocid]._candidates.push_back(edgeDestination);

        // In case the current Node has already enough candidates, we skip the process.
        if (sparseGraph._nodes[sourceDocid]._candidates.size() > maxCandidates) {
          idxOfNodesLeft.erase(sourceDocid);
          break;
        }
      }
    }
sparseGraph.sample.PauseMicroseconds(); // TOGO
  }
}

// Given the superHash vector, the current start and end indexes, find the range [start, end],
// such that all docIDs in range have the same superhash.
// Assumption: startIdx equal to endIdx (endIdx is updated here), while startIdx is stable.
void findRangeWithSameSuperHash(const vector<superHash>& superHashVec, const uint startIdx, uint& endIdx) {
  assert(startIdx == endIdx);
  assert(startIdx < superHashVec.size());

  string curSuperHash = superHashVec[startIdx]._superHash;
  for (uint i = startIdx + 1; i < superHashVec.size(); ++i) {
    if (curSuperHash == superHashVec[i]._superHash) {
      ++endIdx;
    } else {
      break;
    }
  }
  // In case, we check the last entry, the endIdx is already set to startIdx.
  assert(startIdx <= endIdx);
}

double computeIntersection(minHash* a, minHash*b ) {
  uint intersectionSz = 0;
  uint sz = a->_minH.size();
  assert(sz == b->_minH.size());
  for (uint i = 0; i < sz; ++i) {
    if (a->_minH[i] == b->_minH[i]) {
      ++intersectionSz;
    }
  }
  return intersectionSz;
}

double computeJaccardSimilarity(minHash* a, minHash* b) {
  uint intersectionSz = 0;
  uint sz = a->_minH.size();
  assert(sz == b->_minH.size());
  assert(sz > 0);
  for (uint i = 0; i < sz; ++i) {
    if (a->_minH[i] == b->_minH[i]) {
      ++intersectionSz;
    }
  }
  uint unionSz = sz;
  double jaccardScore = (double) intersectionSz / (double) sz;
  return jaccardScore;
}

// Given a vector containing the signature min-hashes, a superHashTemplate and a baseDocID,
// create a superHash by concatenating min-hashes accordingly, hashing the concatenated
// result with MD5() and returning the superHash.
// Note that the baseDocID is removed from any docID and must be reconstructed after TSP.
superHash createSuperHash_MD5(minHash& minHashes, 
                              const vector<unsigned short>& superHashTemplate,
                              const uint baseDocID) {
/* // Benchmark for 100k artificial docs for all thresholds, total: 3.88 minutes.
  // String version.
  // Concatenate several min-hash signatures.
  string superHashStr;
  for (uint i = 0; i < superHashTemplate.size(); ++i) {
    superHashStr += UshortToStr( minHashes._minH[i] );
  }
  // Apply MD5() to the concatenation.
  MD5 md5;

  // Create new superhash by reducing the docID in range [0, x).
  assert(minHashes._docid >= baseDocID);
  uint newDocID = minHashes._docid - baseDocID;  
  superHash newSuperHash( md5.digestString(&superHashStr[0]), &minHashes, newDocID);

  return newSuperHash;
*/

///*
  // Bechmark for 100k artificial docs for all thresholds, total: 0.6 minutes.
  assert(minHashes._minH.size() > 0);
  size_t sZ = superHashTemplate.size()*sizeof(minHashes._minH[0]); 
  unsigned char superHashAr[ sZ + 1];  // +1 for the termination '\0'.
  memset((void*) &superHashAr, 0, sizeof(superHashAr));
  superHashAr[ sZ ] = '\0';

  // Byte version.
  for (uint i = 0, j = 0; i < superHashTemplate.size(); ++i, j += 2) {
    memcpy((void*) &superHashAr[j], (void*) &minHashes._minH[i], sizeof(minHashes._minH[i]));
  }

  // Apply MD5() to the concatenation.
  MD5 md5;

  // Create new superhash by reducing the docID in range [0, x).
  assert(minHashes._docid >= baseDocID);
  uint newDocID = minHashes._docid - baseDocID; // Note that we may reduce 1 for the final code if the docIDs start from 1.
  //superHash newSuperHash( md5.digestString((char*) &superHashAr[0]), &minHashes, newDocID); 
  superHash newSuperHash( md5.digestMemory(&superHashAr[0], sizeof(superHashAr)), &minHashes, newDocID); 
  return newSuperHash;
//*/
}

// Given a unsigned short number, return it in string.
string UshortToStr(const unsigned short number) {
  stringstream ss;
  ss << number;
  //cout << number << " -> " << ss.str() << endl;
  return ss.str();
}

// Given the docID, the vector with the termIDs, the min-wise hash functions,
// compute the min-hashes and return the minhash object.
inline minHash createMinHashes(const uint docID,
                        const vector<uint>& tids,
                        const vector<uint>& hashFunctions) {
  vector<unsigned short> minHashes(hashFunctions.size(), USHRT_MAX);
  for (uint i = 0; i < hashFunctions.size(); ++i) {
    for (uint tidIdx = 0; tidIdx < tids.size(); ++tidIdx) {
      // Note, plus 1 in order to ensure that termId always start with 1!
      minHashes[i] = min(hashF(tids[tidIdx] + 1, i, domain, hashFunctions), minHashes[i]);
    }
  }
  minHash tmpMinHash(docID, tids.size(), minHashes);
  return tmpMinHash;
}

// Given the termID, the i-th min-wise hash function, the domain and the min-wise
// hash functions, return the results of the computation.
inline unsigned short hashF(const uint termId,
          const uint ithHashFunction,
          const uint domain,
          const vector<uint>& hashFunctions) {
  return ((termId*hashFunctions[ithHashFunction])%domain);
}


// Given the filename with the superHashTemplates, load the corresponding
// templates in a vector of ushorts and return it.
vector<vector<unsigned short> > loadSuperHashTemplates(const string file) {
  vector<vector<unsigned short> > superHashTemplates;
  ifstream fH(file.c_str());
  if (fH.good()) {
    while (fH.good()) {
      string line;
      getline(fH, line);
      fH >> ws;
      vector<unsigned short> tmpSuperHashTemplate = TokenizeStringToUShort(line);
      superHashTemplates.push_back(tmpSuperHashTemplate);
    }
  }  
  fH.close();
  // Reporting.
  cout << "Total superHashTemplates loaded " << superHashTemplates.size() << " from file " << file << endl;
  return superHashTemplates;
}

// Given the number of nodes, construct an artificial graph and return it.
Graph constructArtificialGraph(const uint numNodes) {
  // Initialize a graph representation.
  Graph sparseGraph(numNodes);

  // Temporary edge.
  Edge tmpEdge(1, 15);
  Edge tmpEdge1(4, 3);

  // Adding edges artificially.
  sparseGraph._nodes[0].addEdge(tmpEdge);
  sparseGraph._nodes[0].addEdge(tmpEdge1);
  cout << "adding edge 0->1 (15)" << endl;
  cout << "adding edge 0->4 (3)" << endl;
  
  Edge tmpEdge2(0, 15);
  Edge tmpEdge3(2, 1);
  Edge tmpEdge3_5(4, 11);
  sparseGraph._nodes[1].addEdge(tmpEdge2);
  sparseGraph._nodes[1].addEdge(tmpEdge3);
  sparseGraph._nodes[1].addEdge(tmpEdge3_5);
  cout << "adding edge 1->0 (15)" << endl;
  cout << "adding edge 1->2 (1)" << endl;
  cout << "adding edge 1->4 (11)" << endl;

  Edge tmpEdge4(4, 7);
  Edge tmpEdge5(1, 1);
  sparseGraph._nodes[2].addEdge(tmpEdge4);
  sparseGraph._nodes[2].addEdge(tmpEdge5);
  cout << "adding edgea 2->4 (7)" << endl;
  cout << "adding edge 2->1 (1)" << endl;

  Edge tmpEdge6(0, 3);
  Edge tmpEdge7(2, 7);
  Edge tmpEdge8(1, 11);
  sparseGraph._nodes[4].addEdge(tmpEdge6);
  sparseGraph._nodes[4].addEdge(tmpEdge7);
  sparseGraph._nodes[4].addEdge(tmpEdge8);
  cout << "adding edgea 4->0 (3)" << endl;
  cout << "adding edge 4->2 (7)" << endl;
  cout << "adding edge 4->1 (11)" << endl;

  return sparseGraph;
}

//TODO(dimopoulos): complete impl.
void printBin(unsigned long num) {
  vector<bool> bin;
  for (size_t i = 0; i < sizeof(unsigned long)*8; ++i) {
    bin.push_back(((num>>i)&1));
    num >>= 1;
  }
  for (int i = bin.size() - 1; i >= 0; --i) {
    cout << bin[i];
  }
  cout << endl;
}

// Given a string, tokenize it and return a vector of uints.
vector<unsigned short> TokenizeStringToUShort(string str) {
  vector<unsigned short> tokens;
  istringstream iss( str );
  do {
    string token;
    iss >> token;
    if (token=="") continue; //cout << "fix this please!"<<endl; // TODO fix
      tokens.push_back( StrToUShort(token) );
    } while (iss);
  return tokens;
}

// Given a string, return the number in ushort.
unsigned short StrToUShort(const string number_in_string) {
  unsigned short number = 0;
  istringstream(number_in_string) >> number;
  return number;
}

// Given the number of minhashes, the base file and assuming that the corresponding
// setup about #superhashes and #ofconcatenatedMinhashesPerShard for each threshold
// are stored in constant tables, write out each superhash template accordingly.
void createSuperhashTemplatesGivenSetup(const uint numMinhash, const string basefile) {
  // Reporting.
  cout << "Create superHash templates for " << numMinhash << " minhashes on basefile: "
       << basefile << endl;

  for (unsigned short i = 0; i < numThresholds; ++i) {
    string output = basefile + UshortToStr(i);
    writeoutSuperHashTemplates( numMinhashPerSuperhashForTheta[ i ], numSuperhashForTheta[ i ], 
                                numMinhash, output);
    // Reporting.
    cout << "Threshold: " << thresholds[i] << " #sh: " << numSuperhashForTheta[i] 
         << " #concatenatedMinhashPerSH: " << numMinhashPerSuperhashForTheta[i] 
         << " created and written in fie: " << output << endl;
  }
}

// Given the number of min-hash to concatenate per superhash (superHashSize),
// the number of superHash templates to create, the number of min-hashes and 
// the file to output the superhash templates, construct super hashes accordingly.
void writeoutSuperHashTemplates(const uint superHashSize, const uint numSuperHashes,
                                const uint numHashFunctions, const string file) {
  // Create index vector of size number of hash functions.
  vector<uint> hashIdx(numHashFunctions, 0);
  for (uint i = 0; i < hashIdx.size(); ++i) {
    hashIdx[i] = i;
  }
  ofstream outFh(file.c_str());
  // Create a superHash template for each superHash and write it out.
  for (uint i = 0; i < numSuperHashes; ++i) {
    vector<uint> tmp = hashIdx;
    randomlySelectKNumbers(superHashSize, tmp);
    for (uint j = 0; j < tmp.size(); ++j) {
      outFh << tmp[j] << " ";
    }
    outFh << endl;
  }
  outFh.close();
  // Reporting.
  cout << numSuperHashes << " superHashTemplates of size " << superHashSize << " on "
       << numHashFunctions << " minHashFunctions were written out in file " << file << endl;
}

// Given the vector with min-wise hash functions and k, the number of documents, create
// k documents and foreach one, construct the min-hashes and return them as a vector.
vector<minHash> createMinHashesForKDocs(const uint k, vector<uint>& hashFuncs) {
  assert(k > 0);
  // Create k docs.
  vector<vector<uint> > docs = createDocs(k, 10, 100, 300);
  //printDocs(docs);

  // TODO(dimopoulos): resize it ? it is costly
  vector<minHash> minHashVec; // maintain all minhashes
  for (uint docID = 0; docID < docs.size(); ++docID) {
    minHashVec.push_back(createMinHashes(docID, docs[docID], hashFuncs));
  }
  return minHashVec;
}

// Given the number of documents, the minimum document size (in terms of termIDs), the maximum
// document size and the maximum termID, create this many docs and return them in a vector
// of vector<uint>.
vector<vector<uint> > createDocs(const uint numDocs,
                                 const uint minDocSz,
                                 const uint maxDocSz,
                                 const uint maxTid) {
  vector<vector<uint> > docs;
  for (uint i = 0; i < numDocs; ++i) {
    vector<uint> tmpDoc = createDoc(minDocSz, maxDocSz, maxTid);
    docs.push_back(tmpDoc);
  }
  // Reporting.
  cout << "The number of docs created is " << numDocs << " with maxTermId: " << maxTid << endl;
  return docs;
}

vector<uint> createDoc(const uint minDocSz,
                       const uint maxDocSz,
                       const uint maxTid) {
  uint curDocSz = randInRange(minDocSz, maxDocSz);
  vector<uint> doc(curDocSz, 0);
  vector<bool> tidBitset(maxTid, 0);
  uint randTid = 0;
  for (uint i = 0; i < curDocSz; ++i) {
    do {
      randTid = randInRange(1, maxTid);
    } while (tidBitset[randTid]);
    tidBitset[randTid] = 1;
    doc[i] = randTid;
  }
  return doc;
}

// Random generator function by Torsten.
uint randInRangeMS(const uint min, const uint max, int* seed) {
  assert(min<max);
  assert(min>=0);
  int lo;
  int hi;
  int test;

  hi=(*seed)/127773;
  lo=(*seed) % 127773;
  test=16807*lo-2836*hi;
  if (test>0) *seed=test;
  else *seed=test+2147483647;
  return (min + (uint)((max-min+1)*((double)(*seed)/(double)2147483647))) ;  
}

// Given the minimum and the maximum unsigned integer, return an integer at random in the given range.
uint randInRange(const uint min, const uint max) {
  if (min == max) return min;
  assert(min<max);
  assert(min>=0);
  return (min + (rand()%(max-min + 1)));
}

// Given the filename with the min-wise hash functions, load it and return a vector of uints.
vector<uint> loadHashFunctions(string hashFunctionsFile) {
  FILE* fp = fopen(hashFunctionsFile.c_str(), "r");
  if (fp == NULL) {
    cout << "Could not open file: " << hashFunctionsFile << ". Exiting." << endl;
  }
  uint hashTmp = 0;
  vector<uint> hashes;
  while (fscanf(fp, "%u\n", &hashTmp) != EOF) {
    hashes.push_back(hashTmp);
  }
  fclose(fp);

  // Reporting.
  cout << "The number of hash functions loaded from file: " << hashFunctionsFile << " is " << hashes.size() << endl;

  return hashes;
}

// Given the vector of min-wise hash functions and a filename, write out the hash functions to the file.
void writeOut(const vector<uint> allPrimes, string hashFunctionsFile) {
  ofstream outfH(hashFunctionsFile.c_str());
  for (size_t i = 0; i < allPrimes.size(); ++i) {
    outfH << allPrimes[i] << endl;
  }
  outfH.close();
  // Reporting.
  cout << "The hashes of size: " << allPrimes.size() << " were writen out to file: "
       << hashFunctionsFile << endl;
}

// bitset size must be the same with max in arg list.
vector<uint> getPrimesUpTo(const uint max) {
  uint pos = 2;  // current position inspected
  bitset<primeMax> isPrime;
  isPrime.set(0);
  isPrime.set(1);
  vector<uint> primes;
  primes.push_back(2);

  while (pos < max) {  // loop until position equals max
    // clear off all multiples of current prime until max
    for (size_t i = pos*2; i < max; i+= pos) {
      isPrime.set(i);
    }

    // find next prime or next unset bit
    uint nextPrime = 1;
    for (size_t i = pos+1; i < max; ++i) {
      if (!isPrime[i]) {
        nextPrime = i;
        break;
      }
    }
 
    // check for termination condition
    if (nextPrime == 1) {
      break;
    } 
    pos = nextPrime;
    primes.push_back(nextPrime);
  } 
  return primes;
}

// Given a vector<uint>, maintain all the numbers more than the given value primeMin.
void keepPrimesMoreThan(const uint primeMIN, vector<uint>& primes) {
  uint borderIdx = 0;
  while (primes.size() > borderIdx && primes[borderIdx] < primeMin) {
    ++borderIdx;
  } // find last prime before range
  --borderIdx;
  primes.erase(primes.begin(), primes.begin() + borderIdx); // remove all other primes less than primeMin
}

// Given a parameter k and a vector of unsigned ints, permute them and return k of them.
void randomlySelectKNumbers(const uint k, vector<uint>&primes) {
  assert(primes.size() >= k);
  random_shuffle(primes.begin(), primes.end());
  primes.resize(k);
}

// Given a vector of uint, print it out.
void printVec(const vector<uint> vec) {
  for (size_t i = 0; i < vec.size(); ++i) {
    cout << "idx: " << i << " val: " << vec[i] << endl;
  }
}

// Given a vector<ushort>, print out the vector.
void printUshortVec(const vector<unsigned short> vec) {
  for (size_t i = 0; i < vec.size(); ++i) {
    cout << "idx: " << i << " val: " << vec[i] << endl;
  }
}

// Given a vector with termIDs, print it out.
void printDocs(const vector<vector<uint> > vec) {
  for (size_t i = 0; i < vec.size(); ++i) {
    cout << "docid: " << i << " with tids: ";
    for (size_t j = 0; j < vec[i].size(); ++j) {
      cout << vec[i][j] << " ";
    }
    cout << endl;
  }
}

// Given a parameter k and the maximum value max, randomly select k number with value less than max.
vector<uint> randomlySelectKNumbersUpTo(const uint k, const uint max) {
  assert(k>0);
  assert(max>0);
  vector<uint> numVec(k, 0);
  for (uint i = 0; i < k; ++i) {
    numVec[i] = rand()%max;
  } 
  return numVec;
}

// Given a vector of numbers and the number of elements to maintain, select k at random and return them.
void selectKRandom(vector<uint>& randNodes, int numElementsToSample) {
  uint randNumber = 0;
  for (int i = 0; i < numElementsToSample; ++i) {
    randNumber = randInRange(i, randNodes.size()); // OLD rand() generator.
    //randNumber = randInRangeMS(i, randNodes.size(), &seed); // New pseudorandom generator.
    uint tmp = randNodes[ randNumber ];
    randNodes[ randNumber ] = randNodes[ i ];
    randNodes[ i ] = tmp;
  }
  randNodes.resize(numElementsToSample);
}
