#include <fstream>
#include <iostream>
#include <vector>

#include "MCHBase/ClusterBlock.h"
#include "MCHBase/TrackBlock.h"

using namespace std;
using namespace o2::mch;

//_________________________________________________________________________________________________
struct TrackStruct {
  TrackParamStruct paramAtVertex{};
  double dca = 0.;
  double rAbs = 0.;
  TrackParamStruct paramAt1stCluster{};
  double chi2 = 0.;
  std::vector<ClusterStruct> clusters{};
};

//_________________________________________________________________________________________________
int ReadNextEvent(ifstream& inFile, int version, std::vector<TrackStruct>& tracks, int& size);
void ReadTrack(ifstream& inFile, TrackStruct& track, int version);
void WriteTracks(std::vector<TrackStruct>& tracks, ofstream& outFile);

//_________________________________________________________________________________________________
void ConvertToASCII(string trackFileName, int versionFile)
{
  /// read the binary file of tracks and convert it to ASCII format

  // open files
  ifstream inFileTrack(trackFileName, ios::binary);
  if (!inFileTrack.is_open()) {
    cout << "fail opening file" << endl;
    return;
  }
  ofstream outFileTrack("O2Tracks.txt", ios::out);
  if (!outFileTrack.is_open()) {
    return;
  }

  std::vector<TrackStruct> tracks{};
  while (true) {

    // get tracks
    int size(0);
    int event = ReadNextEvent(inFileTrack, versionFile, tracks, size);
    if (event < 0) {
      // reaching end of file
      break;
    }

    // write tracks
    outFileTrack << "------ event " << event << " ------" << endl;
    outFileTrack << "size = " << size << endl;
    WriteTracks(tracks, outFileTrack);
  }

  inFileTrack.close();
  outFileTrack.close();
}

//_________________________________________________________________________________________________
int ReadNextEvent(ifstream& inFile, int version, std::vector<TrackStruct>& tracks, int& size)
{
  /// read the next event in the input file

  int event(-1);
  if (!inFile.read(reinterpret_cast<char*>(&event), sizeof(int))) {
    // reaching end of file
    return -1;
  }

  inFile.read(reinterpret_cast<char*>(&size), sizeof(int));
  if (size == 0) {
    // empty event
    tracks.clear();
    return event;
  }

  int nTracks(0);
  inFile.read(reinterpret_cast<char*>(&nTracks), sizeof(int));
  tracks.resize(nTracks);

  for (Int_t iTrack = 0; iTrack < nTracks; ++iTrack) {
    ReadTrack(inFile, tracks[iTrack], version);
  }

  return event;
}

//_________________________________________________________________________________________________
void ReadTrack(ifstream& inFile, TrackStruct& track, int version)
{
  /// read one track from the input file

  if (version > 2) {
    inFile.read(reinterpret_cast<char*>(&(track.paramAtVertex)), sizeof(TrackParamStruct));
    inFile.read(reinterpret_cast<char*>(&(track.dca)), sizeof(double));
    inFile.read(reinterpret_cast<char*>(&(track.rAbs)), sizeof(double));
  }

  inFile.read(reinterpret_cast<char*>(&(track.paramAt1stCluster)), sizeof(TrackParamStruct));

  if (version > 1) {
    inFile.read(reinterpret_cast<char*>(&(track.chi2)), sizeof(double));
  }

  int nClusters(0);
  inFile.read(reinterpret_cast<char*>(&nClusters), sizeof(int));
  track.clusters.resize(nClusters);
  for (Int_t iCl = 0; iCl < nClusters; ++iCl) {
    inFile.read(reinterpret_cast<char*>(&(track.clusters[iCl])), sizeof(ClusterStruct));
  }
}

//_________________________________________________________________________________________________
void WriteTracks(std::vector<TrackStruct>& tracks, ofstream& outFile)
{
  /// write the tracks in the output file
  outFile << "nTracks = " << tracks.size() << endl;
  for (const auto& track : tracks) {
    outFile << "{x = " << track.paramAtVertex.x << ", y = " << track.paramAtVertex.y << ", z = " << track.paramAtVertex.z
            << ", px = " << track.paramAtVertex.px << ", py = " << track.paramAtVertex.py << ", pz = " << track.paramAtVertex.pz
            << ", sign = " << track.paramAtVertex.sign << "}" << endl;
    outFile << track.dca << endl;
    outFile << track.rAbs << endl;
    outFile << "{x = " << track.paramAt1stCluster.x << ", y = " << track.paramAt1stCluster.y << ", z = " << track.paramAt1stCluster.z
            << ", px = " << track.paramAt1stCluster.px << ", py = " << track.paramAt1stCluster.py << ", pz = " << track.paramAt1stCluster.pz
            << ", sign = " << track.paramAt1stCluster.sign << "}" << endl;
    outFile << track.chi2 << endl;
    outFile << "nClusters = " << track.clusters.size() << endl;
    for (const auto& cluster : track.clusters) {
      outFile << "{x = " << cluster.x << ", y = " << cluster.y << ", z = " << cluster.z << ", ex = " << cluster.ex
              << ", ey = " << cluster.ey << ", uid = " << cluster.uid << "}" << endl;
    }
  }
}
