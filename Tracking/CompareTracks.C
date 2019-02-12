#include <fstream>
#include <iostream>

#include "MCHBase/ClusterBlock.h"
#include "MCHBase/TrackBlock.h"

using namespace std;
using namespace o2::mch;

bool AreTrackParamCompatible(const TrackParamStruct& param1, const TrackParamStruct& param2, double precision);
void PrintResiduals(const TrackParamStruct& param1, const TrackParamStruct& param2);

//_________________________________________________________________________________________________
void CompareTracks(string inFileName1, string inFileName2, double precision = 1.e-4, bool printAll = false)
{
  /// Compare the tracks stored in the 2 binary files
  
  // open files
  ifstream inFile1(inFileName1,ios::binary);
  ifstream inFile2(inFileName2,ios::binary);
  if (!inFile1.is_open() || !inFile2.is_open()) {
    cout << "fail opening files" << endl;
    return;
  }

  int event1(0), event2(0);
  int size1(0), size2(0);
  int nTracks1(0), nTracks2(0);
  TrackParamStruct trackParam1{}, trackParam2{};
  int nClusters1(0), nClusters2(0);
  ClusterStruct cluster1{}, cluster2{};
  bool areCompatible(false);
  int nDifferences(0);

  while (!inFile1.eof() || !inFile2.eof()) {

    // find next non-empty event in both files
    size1 = 0;
    while (size1 == 0 && !inFile1.eof()) {
      inFile1.read(reinterpret_cast<char*>(&event1),sizeof(int));
      inFile1.read(reinterpret_cast<char*>(&size1),sizeof(int));
    }
    size2 = 0;
    while (size2 == 0 && !inFile2.eof()) {
      inFile2.read(reinterpret_cast<char*>(&event2),sizeof(int));
      inFile2.read(reinterpret_cast<char*>(&size2),sizeof(int));
    }
    if (size1 == 0 && size2 == 0) {
      continue;
    }
    if (event1 != event2) {
      cout << "different event number" << endl;
      return;
    }
    if (size1 != size2) {
      cout << "different event size" << endl;
      return;
    }

    // get the number of tracks
    inFile1.read(reinterpret_cast<char*>(&nTracks1),sizeof(int));
    inFile2.read(reinterpret_cast<char*>(&nTracks2),sizeof(int));
    if (nTracks1 != nTracks2) {
      cout << "different number of tracks" << endl;
      return;
    }

    for (Int_t iTrack = 0; iTrack < nTracks1; iTrack++) {

      // compare the track parameters
      inFile1.read(reinterpret_cast<char*>(&trackParam1),sizeof(TrackParamStruct));
      inFile2.read(reinterpret_cast<char*>(&trackParam2),sizeof(TrackParamStruct));
      areCompatible = AreTrackParamCompatible(trackParam1, trackParam2, precision);
      if (!areCompatible) ++nDifferences;
      if (printAll || !areCompatible) {
        PrintResiduals(trackParam1, trackParam2);
      }

      // get the number of clusters
      inFile1.read(reinterpret_cast<char*>(&nClusters1),sizeof(int));
      inFile2.read(reinterpret_cast<char*>(&nClusters2),sizeof(int));
      if (nClusters1 != nClusters2) {
        cout << "different number of clusters" << endl;
        return;
      }

      for (Int_t iCl = 0; iCl < nClusters1; iCl++) {

        // compare the clusters
        inFile1.read(reinterpret_cast<char*>(&cluster1),sizeof(ClusterStruct));
        inFile2.read(reinterpret_cast<char*>(&cluster2),sizeof(ClusterStruct));
        if (cluster1.uid != cluster2.uid) {
          cout << "different clusters" << endl;
          return;
        }

      }

    }

  }

  cout << "number of different track parameters = " << nDifferences << endl;

  inFile1.close();
  inFile2.close();

}

//_________________________________________________________________________________________________
bool AreTrackParamCompatible(const TrackParamStruct& param1, const TrackParamStruct& param2, double precision)
{
  /// compare track parameters within precision
  return (abs(param2.x - param1.x) <= precision &&
          abs(param2.y - param1.y) <= precision &&
          abs(param2.z - param1.z) <= precision &&
          abs(param2.px - param1.px) <= precision &&
          abs(param2.py - param1.py) <= precision &&
          abs(param2.pz - param1.pz) <= precision &&
          param2.sign == param1.sign);
}

//_________________________________________________________________________________________________
void PrintResiduals(const TrackParamStruct& param1, const TrackParamStruct& param2)
{
  /// print param2 - param1
  cout << "{dx = " << param2.x - param1.x
  << ", dy = " << param2.y - param1.y
  << ", dz = " << param2.z - param1.z
  << ", dpx = " << (param2.px - param1.px) << " (" << 100.* (param2.px - param1.px) / param1.px << "\%)"
  << ", dpy = " << (param2.py - param1.py) << " (" << 100.* (param2.py - param1.py) / param1.py << "\%)"
  << ", dpz = " << (param2.pz - param1.pz) << " (" << 100.* (param2.pz - param1.pz) / param1.pz << "\%)"
  << ", dsign = " << param2.sign - param1.sign
  << "}" << endl;
}

