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

using namespace std;

const uint LARGENUMBER = 20000000;
const uint primeMax = 20000000;
const uint primeMin = 10000000;
const uint numHashFunctions = 100;  // # min-hashes
const uint numSuperHashes = 5; //20; # super hashes
const uint superHashSize = 3; // glue 3 min-hash to get a SH
const uint domain = 15773;
const string hashFunctionsFile = "hashFunctions_100";
const string superHashTemplatesFile = "superHashTemplates";
const uint similarityMeasure = 0; // 0: Jaccard, 1: intersection.
const uint numThresholds = 2; //6;
const double thresholds[ 2 ] = { 0.1, 0.05 }; //0.99, 0.95, 0.90, 0.80, 0.65 };
const uint maxCandidates = 5; //100;

typedef unsigned int uint;

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
   uint _docID;
   uint _degree;

   // Interface for adding/removing edges in Node.
   // Largest/Smallest similarity edge at the top.
   priority_queue<Edge, vector<Edge>, edgeComparatorLargestFirst> _edgeHeap;  

   // Given a parameter (global) maxCandidates per Node, return if we should find
   // more high similarity documents or not.
   bool hasEnoughCandidates() {
     return (_degree >= maxCandidates);
   }
 
   // Given the maxCapacity of edges per node, return how many edges we could potentially add.
   int numRemainingEdges(const uint& maxCapacity) {
     assert(maxCapacity > _degree);
     return (maxCapacity - _degree);
   }

   // Add edge in the priority queue and update the degree.
   void addEdge(const Edge& edge, const uint& maxCapacity) {
     // Note: left-to-right evaluation is guaranteed in &&/||.
     if (_edgeHeap.size() < maxCapacity) {
       _edgeHeap.push(edge); // push edge
       ++_degree;  // increase degree
     }
   }   

   // Constructor.
   Node() { 
     _degree = 0; 
     _docID = 0;
   }
};

class Graph {
  public: 
    vector<Node> _nodes;
    Graph(const uint sz) { _nodes.resize(sz, Node()); }
    vector<uint> GreedyNearestNeighborTSP();
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

class superHash {
  public:
   unsigned long _superHash;
   minHash* _minH;
   superHash(unsigned long inSuperHash, minHash* inMinHashPtr) {
     _superHash = inSuperHash;
     _minH = inMinHashPtr;
   }
   superHash(minHash* inMinHashPtr) {
     _superHash = 0;
     _minH = inMinHashPtr;
   }
};

bool superHashComparator (const superHash i, const superHash j) { 
  return (i._superHash<j._superHash); 
}

////////////////////// Method Signatures /////////////////////////
vector<uint> getPrimesUpTo(const uint max);
void keepPrimesMoreThan(const uint primeMIN, vector<uint>& primes);
void randomlySelectKNumbers(const uint k, vector<uint>&primes);
void printVec(const vector<uint> vec);
void printUshortVec(const vector<unsigned short> vec);
unsigned short hashF(const uint termId, const uint ithHashFunction,
          const uint domain, const vector<uint>& hashFunctions);
void writeOut(const vector<uint> allPrimes, string hashFunctionsFile);
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
uint randInRange(const uint min, const uint max);
void printDocs(const vector<vector<uint> > vec);
vector<unsigned short> TokenizeStringToUShort(string str);
unsigned short StrToUShort(const string number_in_string);
vector<uint> randomlySelectKNumbersUpTo(const uint k, const uint max);
void writeoutSuperHashTemplates(const uint superHashSize, const uint numSuperHashes,
                                const uint numHashFunctions, const string file);
vector<vector<unsigned short> > loadSuperHashTemplates(const string file);
Graph createSparseGraph(vector<minHash>& minHashes, 
                        vector<vector<unsigned short> >& superHashTemplates);
superHash createSuperHash(minHash& minHashes, 
                          const vector<unsigned short>& superHashTemplate);
void printBin(unsigned long i);
vector<minHash> createMinHashesForKDocs(const uint k, vector<uint>& hashFuncs);
double computeJaccardSimilarity(minHash* a, minHash* b);
double computeIntersection(minHash* a, minHash*b );
void selectNearestNeighborCandidates(const vector<superHash> superHashVec, Graph& sparseGraph); 
void findRangeWithSameSuperHash(const vector<superHash> superHashVec, const uint start, uint& end);
void addEdgesToSimilarDocs(const vector<superHash> superHashVec,
                           const uint startIdx,
                           const uint endIdx,
                           const double upperBoundThreshold,
                           const double lowerBoundThreshold,
                           Graph& sparseGraph);
Graph constructArtificialGraph(const uint numNodes);
void selectKRandom(vector<uint>& randNodes, int numElementsToSample);

////////////////////// Driver /////////////////////////
int main() {
  srand(time(0));
/*
  // 1. Create hash functions for min-wise hashing.
  vector<uint> allPrimes = getPrimesUpTo(primeMax);
  keepPrimesMoreThan(primeMin, allPrimes);
  randomlySelectKNumbers(numHashFunctions, allPrimes);
  // 2. Writing out the hash functions in file hashFunctionsFile.
  writeOut(allPrimes, hashFunctionsFile);
  // 3. Loading the hash functions from file hashFunctionsFile.
  vector<uint> hashFuncs = loadHashFunctions(hashFunctionsFile);

  // Print primes.
  printVec(allPrimes);
*/

  vector<uint> hashFuncs = loadHashFunctions(hashFunctionsFile);
  //vector<uint> doc = createDoc(20, 100, 100);
  //uint docid = 12;
  //printVec(doc);

  //vector<minHash> minHashDS;  // maintain all minhashes
  //minHashDS.push_back(createMinHashes(docid, doc, hashFuncs));
  //printUshortVec(minHashes);

  // 4. Create and write out superhash templates.
  //writeoutSuperHashTemplates(superHashSize, numSuperHashes, numHashFunctions, superHashTemplatesFile);
  // 5. Load superhash templates.
  vector<vector<unsigned short> > superHashTemplates = loadSuperHashTemplates(superHashTemplatesFile);

  // 6. Generate signatures for documents. 
  vector<minHash> minHashDS = createMinHashesForKDocs(10000, hashFuncs);
  cout << "###############################" << endl;
  Graph sparseGraph = createSparseGraph(minHashDS, superHashTemplates);

  // Testing artificial graph creation.
  //Graph sparseGraph = constructArtificialGraph(10);
  vector<uint> docidAssignment = sparseGraph.GreedyNearestNeighborTSP();
  printVec(docidAssignment);  

  return 0;
}

////////////////////// Methods' implementation /////////////////////////

vector<uint> Graph::GreedyNearestNeighborTSP() {
  // Vector to store the docid assignment (oldDocid -> newDocid).
  vector<uint> docidAssignment(_nodes.size(), 0);

  // Select the edge with the highest similarity in the sparseGraph.
  Edge largestSimilarityEdge(0, 0.0f);
  uint source = 0;  // will be used as the source of the current edge.
  for (uint i = 0; i < _nodes.size(); ++i) {
    if (!_nodes[i]._edgeHeap.empty() &&
       (_nodes[i]._edgeHeap.top()._similarity > largestSimilarityEdge._similarity)) {
      largestSimilarityEdge._similarity = _nodes[i]._edgeHeap.top()._similarity;
      largestSimilarityEdge._docID = _nodes[i]._edgeHeap.top()._docID;
      source = i;
    }
  }
  // Start with a partial path with a vertex/node from the largestSimilarityEdge.
  docidAssignment[0] = source;

  vector<bool> visited(_nodes.size(), false);
  visited[ source ] = true;
  uint numNodesSelected = 1; 
  uint dest = largestSimilarityEdge._docID;

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

    // if not found, select node with the largest remaining similarity.
    if (!foundNextNode) {
      largestSimilarityEdge._similarity = 0.0f;  // re-use edge.
      largestSimilarityEdge._docID = 0;
      bool foundNextNodeHeuristic = false;
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

      // In case that there are no nodes with edges left.
      if (!foundNextNodeHeuristic) {
        curNode = 0;
        while (curNode < _nodes.size() && visited[ curNode ]) {
          ++curNode;
        }
      } 
    }

    // The next node has already been selected.
    visited[ curNode ] = true;
    docidAssignment[ numNodesSelected ] = curNode;
    ++numNodesSelected;  
    source = curNode;
  }

  return docidAssignment;
}

Graph createSparseGraph(vector<minHash>& minHashes, 
                        vector<vector<unsigned short> >& superHashTemplates) {
  assert(superHashTemplates.size() > 0);
  assert(minHashes.size() > 0);
 
  // Initialize a graph representation.
  Graph sparseGraph(minHashes.size());

  // Foreach super hash template.
  uint numSuperHashes = superHashTemplates.size();
  for (uint i = 0; i < numSuperHashes; ++i) {
    vector<superHash> superHashVec;  // maintain super hashes of this round.

    // Foreach doc create the super hash according to the i-th super hash template.
    for (uint docID = 0; docID < minHashes.size(); ++docID) {
      superHashVec.push_back(createSuperHash(minHashes[docID], superHashTemplates[i]));
    }
    // Sort superHashes by superHash value.
    sort(superHashVec.begin(), superHashVec.end(), superHashComparator);

    cout << "# SH: " << i << endl; //TOGO

    // Add edges in the graph.
    selectNearestNeighborCandidates(superHashVec, sparseGraph);
  }
  return sparseGraph;
}

// Foreach threhold, select the nearest neighbors of similar documents.
// Assumption: superHash vector of size more than 0.
void selectNearestNeighborCandidates(const vector<superHash> superHashVec, Graph& sparseGraph) {
  assert(superHashVec.size() > 0);
  assert(numThresholds > 1);

  uint curStartIdx = 0;
  uint curEndIdx = 0;
  uint lastIdx = superHashVec.size();

  while (curStartIdx != lastIdx) {
    // Find docID ranges with same superHashes.
    findRangeWithSameSuperHash(superHashVec, curStartIdx, curEndIdx);

    // In each iteration we only want to add edges between docs with similarity in range (UP, LB].
    for (uint numThreshold = 0; numThreshold < numThresholds; ++numThreshold) {
      const double upperBoundThreshold = (numThreshold == 0) ? 1.0f : thresholds[numThreshold - 1];
      const double lowerBoundThreshold = thresholds[numThreshold];

 //     cout << " threshold range (" << upperBoundThreshold << ", " << lowerBoundThreshold << "]" << endl; //TOGO

      // Add edges to document-nodes with high similarity.
      addEdgesToSimilarDocs(superHashVec, curStartIdx, curEndIdx, upperBoundThreshold,
                            lowerBoundThreshold, sparseGraph);
    }

    ++curEndIdx;  // next start idx.
    curStartIdx = curEndIdx;  // set curStartIdx.
    assert(curEndIdx == curStartIdx);
  }
}

// Given the superHashVec, the start and end idx of docs with the same superhash and the 
// sparsegraph, (i) do nothing if there is only one docid, (ii) compute all pairwise similarity
// scores if #docsWithSameSuperHash is < pairwiseComputationThreshold, (iii) compute similarity
// on a sample of docs.
// TODO(dimopoulos): add to global variables, the pairwiwseComputationThreshold and maxCapacity.
void addEdgesToSimilarDocs(const vector<superHash> superHashVec,
                           const uint startIdx,
                           const uint endIdx,
                           const double upperBoundThreshold,
                           const double lowerBoundThreshold,
                           Graph& sparseGraph) {
  uint numDocsWithSameSuperHash = endIdx - startIdx;
  if (numDocsWithSameSuperHash == 0) return;

//TOGO
//cout << "range: " << numDocsWithSameSuperHash << " starid: " << startIdx << " - " << endIdx << endl;

  uint pairwiseComputationThreshold = 100; // TODO(globals)
  if (numDocsWithSameSuperHash == 1) {
      uint sourceDocid = superHashVec[ startIdx ]._minH->_docid;
      // In case the current Node has already enough candidates, we skip the process.
      if (sparseGraph._nodes[sourceDocid].hasEnoughCandidates()) {
        cout << "docid: " << sourceDocid << " has enough candidates!! " << endl;
        return;
      }
      uint destinationDocid = superHashVec[ endIdx ]._minH->_docid;
      double score = computeJaccardSimilarity(superHashVec[ startIdx ]._minH, superHashVec[ endIdx ]._minH);
      // Intersection is estimated: Jaccard*(|A| + |B|) / (1 + Jacccard), where |A| #terms.
      if (similarityMeasure == 1) { // intersection.
        uint iSz = superHashVec[ startIdx ]._minH->_numTerms;
        uint jSz = superHashVec[ endIdx ]._minH->_numTerms;
        score = score*((double) iSz + (double) jSz) / (1.0f + score);  
      }
     
      if (score >= lowerBoundThreshold && score < upperBoundThreshold) {
        Edge edgeSource(sourceDocid, score);
        sparseGraph._nodes[sourceDocid].addEdge(edgeSource, maxCandidates); 

        //Edge edgeDestination(destinationDocid, score);
        //sparseGraph._nodes[destinationDocid].addEdge(edgeDestination, maxCandidates);
      }
  } else if (numDocsWithSameSuperHash < pairwiseComputationThreshold) {
    for (uint i = startIdx; i < endIdx; ++i) { // compute all scores for each node.
      uint sourceDocid = superHashVec[i]._minH->_docid;
      // In case the current Node has already enough candidates, we skip the process.
      if (sparseGraph._nodes[sourceDocid].hasEnoughCandidates()) {
        cout << "docid: " << sourceDocid << " has enough candidates!! " << endl;
        continue;
      }

      for (uint j = startIdx; j < endIdx; ++j) {
        if (i == j) continue; 
        uint destinationDocid = superHashVec[j]._minH->_docid;
        double score = computeJaccardSimilarity(superHashVec[i]._minH, superHashVec[j]._minH);
        // Intersection is estimated: Jaccard*(|A| + |B|) / (1 + Jacccard), where |A| #terms.
        if (similarityMeasure == 1) { // intersection.
          uint iSz = superHashVec[i]._minH->_numTerms;
          uint jSz = superHashVec[j]._minH->_numTerms;
          score = score*((double) iSz + (double) jSz) / (1.0f + score);  
        }
        if (score >= lowerBoundThreshold && score < upperBoundThreshold) {
          Edge edgeSource(sourceDocid, score);
          sparseGraph._nodes[sourceDocid].addEdge(edgeSource, maxCandidates); 
          //Edge edgeDestination(destinationDocid, score);
          //sparseGraph._nodes[destinationDocid].addEdge(edgeDestination, maxCandidates);
        }
      }
    }  
  } else {  // compute sample of scores for each node/document.
    for (uint i = startIdx; i < endIdx; ++i) {
      uint sourceDocid = superHashVec[i]._minH->_docid;
      // In case the current Node has already enough candidates, we skip the process.
      if (sparseGraph._nodes[sourceDocid].hasEnoughCandidates()) continue;
  
      // Select K random nodes.
      vector<uint> randNodes((endIdx - startIdx), startIdx);
      for (uint cnt = 0; cnt < randNodes.size(); ++cnt) {
        randNodes[ cnt ] += cnt;
      }
      randNodes.erase(randNodes.begin() + i);  // remove the element same with i.
      selectKRandom(randNodes, sparseGraph._nodes[sourceDocid].numRemainingEdges(maxCandidates));
      uint randNode = 0; 

      for (uint j = 0; j < randNodes.size(); ++j) {
        randNode = randNodes[ j ];       
        uint destinationDocid = superHashVec[randNode]._minH->_docid;
        double score = computeJaccardSimilarity(superHashVec[i]._minH, superHashVec[randNode]._minH);
        if (similarityMeasure == 1) { // intersection.
          uint iSz = superHashVec[i]._minH->_numTerms;
          uint randNodeSz = superHashVec[randNode]._minH->_numTerms;
          score = score*((double) iSz + (double) randNodeSz) / (1.0f + score);  
        }
        if (score >= lowerBoundThreshold && score < upperBoundThreshold) {
          Edge edgeSource(sourceDocid, score);
          sparseGraph._nodes[sourceDocid].addEdge(edgeSource, maxCandidates);
          //Edge edgeDestination(destinationDocid, score);
          //sparseGraph._nodes[destinationDocid].addEdge(edgeDestination, maxCandidates); 
        }
      }
    }
  }
}

// Given the superHash vector, the current start and end indexes, find the range [start, end],
// such that all docIDs in range have the same superhash.
// Assumption: startIdx equal to endIdx (endIdx is updated here), while startIdx is stable.
void findRangeWithSameSuperHash(const vector<superHash> superHashVec, const uint startIdx, uint& endIdx) {
  assert(startIdx == endIdx);
  assert(startIdx < superHashVec.size());

  unsigned long curSuperHash = superHashVec[startIdx]._superHash;
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
  for (uint i = 0; i < sz; ++i) {
    if (a->_minH[i] == b->_minH[i]) {
      ++intersectionSz;
    }
  }
  uint unionSz = sz;
  double jaccardScore = (double) intersectionSz / (double) sz;
  return jaccardScore;
}

// Assumption: superHash of size 8 bytes and each minHash of size 2 bytes.
superHash createSuperHash(minHash& minHashes, 
                          const vector<unsigned short>& superHashTemplate) {
  superHash tmpSH(&minHashes);  // superhash is set to 0 by the constructor.
  assert(superHashSize < 4);  // since SH is 8 bytes, we can only glue 4 min-hash of 8 bytes.
  uint i = 0;
  for (; i < superHashSize - 1; ++i) {
    tmpSH._superHash |= minHashes._minH[i];
    tmpSH._superHash <<= 16;
  }
  assert(numSuperHashes - 1 == i);
  tmpSH._superHash |= minHashes._minH[ numSuperHashes - 1 ];
  return tmpSH;
}

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

inline unsigned short hashF(const uint termId,
          const uint ithHashFunction,
          const uint domain,
          const vector<uint>& hashFunctions) {
  return ((termId*hashFunctions[ithHashFunction])%domain);
}

vector<minHash> createMinHashesForKDocs(const uint k, vector<uint>& hashFuncs) {
  assert(k > 0);
  // Create k docs.
  vector<vector<uint> > docs = createDocs(k, 3, 100, 10000);
  //printDocs(docs);

  // TODO(dimopoulos): resize it ? it is costly
  vector<minHash> minHashVec; // maintain all minhashes
  for (uint docID = 0; docID < docs.size(); ++docID) {
    minHashVec.push_back(createMinHashes(docID, docs[docID], hashFuncs));
  }
  return minHashVec;
}

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

Graph constructArtificialGraph(const uint numNodes) {
  // Initialize a graph representation.
  Graph sparseGraph(numNodes);

  // Temporary edge.
  Edge tmpEdge(1, 15);
  Edge tmpEdge1(4, 3);

  // Adding edges artificially.
  sparseGraph._nodes[0].addEdge(tmpEdge, maxCandidates);
  sparseGraph._nodes[0].addEdge(tmpEdge1, maxCandidates);
  cout << "adding edge 0->1 (15)" << endl;
  cout << "adding edge 0->4 (3)" << endl;
  
  Edge tmpEdge2(0, 15);
  Edge tmpEdge3(2, 1);
  Edge tmpEdge3_5(4, 11);
  sparseGraph._nodes[1].addEdge(tmpEdge2, maxCandidates);
  sparseGraph._nodes[1].addEdge(tmpEdge3, maxCandidates);
  sparseGraph._nodes[1].addEdge(tmpEdge3_5, maxCandidates);
  cout << "adding edge 1->0 (15)" << endl;
  cout << "adding edge 1->2 (1)" << endl;
  cout << "adding edge 1->4 (11)" << endl;

  Edge tmpEdge4(4, 7);
  Edge tmpEdge5(1, 1);
  sparseGraph._nodes[2].addEdge(tmpEdge4, maxCandidates);
  sparseGraph._nodes[2].addEdge(tmpEdge5, maxCandidates);
  cout << "adding edgea 2->4 (7)" << endl;
  cout << "adding edge 2->1 (1)" << endl;

  Edge tmpEdge6(0, 3);
  Edge tmpEdge7(2, 7);
  Edge tmpEdge8(1, 11);
  sparseGraph._nodes[4].addEdge(tmpEdge6, maxCandidates);
  sparseGraph._nodes[4].addEdge(tmpEdge7, maxCandidates);
  sparseGraph._nodes[4].addEdge(tmpEdge8, maxCandidates);
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

unsigned short StrToUShort(const string number_in_string) {
  unsigned short number = 0;
  istringstream(number_in_string) >> number;
  return number;
}

// Construct super hashes.
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

uint randInRange(const uint min, const uint max) {
  assert(min<max);
  assert(min>=0);
  return (min + (rand()%(max-min + 1)));
}


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

void keepPrimesMoreThan(const uint primeMIN, vector<uint>& primes) {
  uint borderIdx = 0;
  while (primes.size() > borderIdx && primes[borderIdx] < primeMin) {
    ++borderIdx;
  } // find last prime before range
  --borderIdx;
  primes.erase(primes.begin(), primes.begin() + borderIdx); // remove all other primes less than primeMin
}

void randomlySelectKNumbers(const uint k, vector<uint>&primes) {
  random_shuffle(primes.begin(), primes.end());
  primes.resize(k);
}

void printVec(const vector<uint> vec) {
  for (size_t i = 0; i < vec.size(); ++i) {
    cout << "idx: " << i << " val: " << vec[i] << endl;
  }
}

void printUshortVec(const vector<unsigned short> vec) {
  for (size_t i = 0; i < vec.size(); ++i) {
    cout << "idx: " << i << " val: " << vec[i] << endl;
  }
}

void printDocs(const vector<vector<uint> > vec) {
  for (size_t i = 0; i < vec.size(); ++i) {
    cout << "docid: " << i << " with tids: ";
    for (size_t j = 0; j < vec[i].size(); ++j) {
      cout << vec[i][j] << " ";
    }
    cout << endl;
  }
}

vector<uint> randomlySelectKNumbersUpTo(const uint k, const uint max) {
  assert(k>0);
  assert(max>0);
  vector<uint> numVec(k, 0);
  for (uint i = 0; i < k; ++i) {
    numVec[i] = rand()%max;
  } 
  return numVec;
}

void selectKRandom(vector<uint>& randNodes, int numElementsToSample) {
  uint randNumber = 0;
  for (int i = 0; i < numElementsToSample; ++i) {
    randNumber = randInRange(i, randNodes.size());
    uint tmp = randNodes[ randNumber ];
    randNodes[ randNumber ] = randNodes[ i ];
    randNodes[ i ] = tmp;
  }
  randNodes.resize(numElementsToSample);
}
