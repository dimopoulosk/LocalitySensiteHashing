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
const uint numHashFunctions = 50;
const uint numSuperHashes = 20;
const uint superHashSize = 3;
const uint numMinSH = 2;
const uint domain = 15773;
const string hashFunctionsFile = "hashFunctions";
const string superHashTemplatesFile = "superHashTemplates";

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

class Node {
  public:
   uint _docID;
   uint _degree;

   // Interface for adding/removing edges in Node.
   // Smallest similarity edge at the top.
   priority_queue<Edge, vector<Edge>, edgeComparator> _edgeHeap;  

   // Add edge in the priority queue and update the degree.
   void addEdge(const Edge& edge, const uint& maxCapacity) {
     // Note: left-to-right evaluation is guaranteed in &&/||.
     if (_edgeHeap.size() < maxCapacity || 
        ((_edgeHeap.size() == maxCapacity) && (edge._similarity > _edgeHeap.top()._similarity))) {
       _edgeHeap.push(edge); // push edge
       ++_degree;  // increase degree
       capacityCheck(maxCapacity); // check priority queue capacity
     }
     assert(_degree == _edgeHeap.size());
   }   

   // Remove smallest similarity edges until size of priority queue is maxCapacity.
   void capacityCheck(const uint maxCapacity) {
     if (_edgeHeap.size() > maxCapacity) {
       while (_edgeHeap.size() > maxCapacity) {
         _edgeHeap.pop();
         --_degree;
       }
     }
     assert(_degree == _edgeHeap.size());
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
};

class minHash {
  public:
   uint _docid;
   vector<unsigned short> _minH;
   minHash() {};
   minHash(const uint inDocid, vector<unsigned short> inMinH) {
     _docid = inDocid; 
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
unsigned short hash(const uint termId, const uint ithHashFunction,
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
void createSparseGraph(vector<minHash>& minHashes, 
                       vector<vector<unsigned short> >& superHashTemplates);
superHash createSuperHash(minHash& minHashes, 
                          const vector<unsigned short>& superHashTemplate);
void printBin(unsigned long i);
vector<minHash> createMinHashesForKDocs(const uint k, vector<uint>& hashFuncs);
double computeJaccardSimilarity(minHash* a, minHash* b);
double computeIntersection(minHash* a, minHash*b );
void addEdges(const vector<superHash> superHashVec, Graph& sparseGraph); 
void findRangeWithSameSuperHash(const vector<superHash> superHashVec, const uint start, uint& end);
void addEdgesToSimilarDocs(const vector<superHash> superHashVec,
                           const uint startIdx,
                           const uint endIdx,
                           Graph& sparseGraph);
////////////////////// Driver /////////////////////////
int main() {
  srand(time(0));
/*
  vector<uint> allPrimes = getPrimesUpTo(primeMax);
  keepPrimesMoreThan(primeMin, allPrimes);
  randomlySelectKNumbers(numHashFunctions, allPrimes);
  writeOut(allPrimes, hashFunctionsFile);
  printVec(allPrimes);
  vector<uint> hashFuncs = loadHashFunctions(hashFunctionsFile);
*/

  vector<uint> hashFuncs = loadHashFunctions(hashFunctionsFile);
  //vector<uint> doc = createDoc(20, 100, 100);
  //uint docid = 12;
  //printVec(doc);

  //vector<minHash> minHashDS;  // maintain all minhashes
  //minHashDS.push_back(createMinHashes(docid, doc, hashFuncs));
  //printUshortVec(minHashes);

  //writeoutSuperHashTemplates(superHashSize, numSuperHashes, numHashFunctions, superHashTemplatesFile);
  vector<vector<unsigned short> > superHashTemplates = loadSuperHashTemplates(superHashTemplatesFile);

  vector<minHash> minHashDS = createMinHashesForKDocs(100000, hashFuncs);
  cout << "###############################" << endl;
  createSparseGraph(minHashDS, superHashTemplates);

  return 0;
}

////////////////////// Methods' implementation /////////////////////////

// TODO(dimopoulos): change method's name.
void createSparseGraph(vector<minHash>& minHashes, 
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

    // Add edges in the graph.
    addEdges(superHashVec, sparseGraph);
  }
  // return sparseGraph;
}

// TODO(dimopoulos): implementation.
// Assumption: superHash vector of size more than 0.
void addEdges(const vector<superHash> superHashVec, Graph& sparseGraph) {
  assert(superHashVec.size() > 0);
  uint curStartIdx = 0;
  uint curEndIdx = 0;
  uint lastIdx = superHashVec.size();

  while (curStartIdx != lastIdx) {
    // Find docID ranges with same superHashes.
    findRangeWithSameSuperHash(superHashVec, curStartIdx, curEndIdx);

    // Add edges to document-nodes with high similarity.
    addEdgesToSimilarDocs(superHashVec, curStartIdx, curEndIdx, sparseGraph);

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
                           Graph& sparseGraph) {
  uint numDocsWithSameSuperHash = endIdx - startIdx;
  if (numDocsWithSameSuperHash == 0) return;
cout << "range: " << numDocsWithSameSuperHash << " starid: " << startIdx << " - " << endIdx << endl;

  uint pairwiseComputationThreshold = 100;
  if (numDocsWithSameSuperHash < pairwiseComputationThreshold) {
    for (uint i = startIdx; i < endIdx - 1; ++i) { // compute all scores.
      for (uint j = i + 1; j < endIdx; ++j) {
        uint sourceDocid = superHashVec[i]._minH->_docid;
        uint destinationDocid =superHashVec[j]._minH->_docid;
        double score = computeJaccardSimilarity(superHashVec[i]._minH, superHashVec[j]._minH);
        Edge edgeSource(sourceDocid, score);
        Edge edgeDestination(destinationDocid, score);
        sparseGraph._nodes[sourceDocid].addEdge(edgeSource, 100);
        sparseGraph._nodes[destinationDocid].addEdge(edgeDestination, 1000);
      }
    }  
  } else {  // compute sample of scores for each node/document.
    for (uint i = startIdx; i < endIdx - 1; ++i) {
      uint sourceDocid = superHashVec[i]._minH->_docid;
      set<uint> selectedNodes;
  
      uint numIterations = 0;
      uint maxIterations = 1000;  // TODO(dimopoulos): move to globals.
      while (numIterations < maxIterations) {  // 
        uint randNode = randInRange(startIdx + 1, endIdx);  // select a Node in range.
        if ((selectedNodes.count(randNode) == 0) && (randNode != i)) {
          selectedNodes.insert(randNode);

          uint sourceDocid = superHashVec[i]._minH->_docid;
          uint destinationDocid =superHashVec[randNode]._minH->_docid;
          double score = computeJaccardSimilarity(superHashVec[i]._minH, superHashVec[randNode]._minH);
          Edge edgeSource(sourceDocid, score);
          Edge edgeDestination(destinationDocid, score);
          sparseGraph._nodes[sourceDocid].addEdge(edgeSource, 100);
          sparseGraph._nodes[destinationDocid].addEdge(edgeDestination, 1000);
        }
        ++numIterations;
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
  assert(superHashSize > numMinSH); // size of each SH must be at least numMinSH 
  for (uint i = 0; i < superHashSize - 1; ++i) {
    tmpSH._superHash |= minHashes._minH[i];
    tmpSH._superHash <<= 16;
  }
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
      minHashes[i] = min(hash(tids[tidIdx] + 1, i, domain, hashFunctions), minHashes[i]);
    }
  }
  minHash tmpMinHash(docID, minHashes);
  return tmpMinHash;
}

inline unsigned short hash(const uint termId,
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
