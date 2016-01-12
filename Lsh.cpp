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

#include "Lsh.h"

using namespace std;

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
  // createSuperhashTemplatesGivenSetup(100, "superHashTemplates_");

  // 5. Generate signatures for documents.
  vector<minHash> minHashDS = createMinHashesForKDocs(300000, hashFuncs);

  // 6. Create sparseGraph by constructing superHashes and adding edges of similar pairs.
  Graph sparseGraph = createSparseGraph(minHashDS, 0);

  // Testing artificial graph creation.
  // Graph sparseGraph = constructArtificialGraph(10);

  // 7. Run TSP on the sparse graph.
  //vector<uint> docidAssignment = sparseGraph.GreedyNearestNeighborTSP();
  //printVec(docidAssignment);  

  return 0;
}
