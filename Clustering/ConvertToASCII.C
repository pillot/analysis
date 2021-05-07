#include <fstream>
#include <iostream>
#include <vector>

#include "DataFormatsMCH/Digit.h"
#include "MCHBase/ClusterBlock.h"

using namespace std;
using namespace o2::mch;

//_________________________________________________________________________________________________
void ReadNextEvent(ifstream& inFile, std::vector<ClusterStruct>& clusters);
void WriteClusters(std::vector<ClusterStruct>& clusters, ofstream& outFile);

//_________________________________________________________________________________________________
void ConvertToASCII(string inFileName, string outFileName)
{
  /// read the binary file of clusters and convert it to ASCII format

  // open files
  ifstream inFile(inFileName, ios::binary);
  if (!inFile.is_open()) {
    cout << "fail opening file" << endl;
    return;
  }
  ofstream outFile(outFileName, ios::out);
  if (!outFile.is_open()) {
    return;
  }

  std::vector<ClusterStruct> clusters{};
  while (inFile.peek() != EOF) {
    ReadNextEvent(inFile, clusters);
    WriteClusters(clusters, outFile);
  }

  inFile.close();
  outFile.close();
}

//_________________________________________________________________________________________________
void ReadNextEvent(ifstream& inFile, std::vector<ClusterStruct>& clusters)
{
  /// read the next event in the input file

  // get the number of clusters
  int nClusters(-1);
  inFile.read(reinterpret_cast<char*>(&nClusters), sizeof(int));
  if (inFile.fail()) {
    throw std::length_error("invalid input");
  }

  // get the number of associated digits
  int nDigits(-1);
  inFile.read(reinterpret_cast<char*>(&nDigits), sizeof(int));
  if (inFile.fail()) {
    throw std::length_error("invalid input");
  }

  if (nClusters < 0 || nDigits < 0) {
    throw std::length_error("invalid input");
  }

  // fill clusters in O2 format, if any
  if (nClusters > 0) {
    clusters.resize(nClusters);
    inFile.read(reinterpret_cast<char*>(&clusters[0]), nClusters * sizeof(ClusterStruct));
    if (inFile.fail()) {
      throw std::length_error("invalid input");
    }
  }

  // skip the digits if any
  if (nDigits > 0) {
    inFile.seekg(nDigits * sizeof(Digit), std::ios::cur);
    if (inFile.fail()) {
      throw std::length_error("invalid input");
    }
  }
}

//_________________________________________________________________________________________________
void WriteClusters(std::vector<ClusterStruct>& clusters, ofstream& outFile)
{
  /// write the clusters in the output file
  outFile << clusters.size() << " clusters:" << endl;
  for (const auto& cluster : clusters) {
    outFile << cluster << endl;
  }
}
