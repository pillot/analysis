#include <fstream>
#include <iostream>
#include <vector>

#include "DataFormatsMCH/TrackMCH.h"
#include "DataFormatsMCH/ClusterBlock.h"
#include "MCHBase/TrackBlock.h"

#include "/Users/PILLOT/Work/Alice/Macros/Tracking/TrackMCHv1.h"

using namespace std;
using namespace o2::mch;

//_________________________________________________________________________________________________
struct ClusterStructV1 {
  float x;      ///< cluster position along x
  float y;      ///< cluster position along y
  float z;      ///< cluster position along z
  float ex;     ///< cluster resolution along x
  float ey;     ///< cluster resolution along y
  uint32_t uid; ///< cluster unique ID
};

//_________________________________________________________________________________________________
struct TrackAtVtxStruct {
  TrackParamStruct paramAtVertex{};
  double dca = 0.;
  double rAbs = 0.;
  int mchTrackIdx = 0;
};

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
template <class T>
int ReadNextEventV5(ifstream& inFile, std::vector<TrackStruct>& tracks, int& size);
template <typename T>
void FillTrack(TrackStruct& track, const TrackAtVtxStruct* trackAtVtx, const T* mchTrack, std::vector<ClusterStruct>& clusters);
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
    int event(-1);
    if (versionFile < 5) {
      event = ReadNextEvent(inFileTrack, versionFile, tracks, size);
    } else if (versionFile == 5) {
      event = ReadNextEventV5<TrackMCHv1>(inFileTrack, tracks, size);
    } else {
      event = ReadNextEventV5<TrackMCH>(inFileTrack, tracks, size);
    }
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
    if (version < 4) {
      inFile.read(reinterpret_cast<char*>(&(track.clusters[iCl])), sizeof(ClusterStructV1));
    } else {
      inFile.read(reinterpret_cast<char*>(&(track.clusters[iCl])), sizeof(ClusterStruct));
    }
  }
}

//_________________________________________________________________________________________________
template <class T>
int ReadNextEventV5(ifstream& inFile, std::vector<TrackStruct>& tracks, int& size)
{
  /// read the next event in the input file

  static int event = -1;
  ++event;

  tracks.clear();

  int nTracksAtVtx(-1);
  if (!inFile.read(reinterpret_cast<char*>(&nTracksAtVtx), sizeof(int))) {
    size = 0;
    return -1; // reaching end of file
  }
  int nMCHTracks(-1);
  inFile.read(reinterpret_cast<char*>(&nMCHTracks), sizeof(int));
  int nClusters(-1);
  inFile.read(reinterpret_cast<char*>(&nClusters), sizeof(int));
  if (nTracksAtVtx == 0 && nMCHTracks == 0) {
    size = 3 * sizeof(int);
    return event; // event is empty
  }

  // read the tracks at vertex, the MCH tracks and the attached clusters
  std::vector<TrackAtVtxStruct> tracksAtVtx(nTracksAtVtx);
  inFile.read(reinterpret_cast<char*>(tracksAtVtx.data()), nTracksAtVtx * sizeof(TrackAtVtxStruct));
  std::vector<T> mchTracks(nMCHTracks);
  inFile.read(reinterpret_cast<char*>(mchTracks.data()), nMCHTracks * sizeof(T));
  std::vector<ClusterStruct> clusters(nClusters);
  inFile.read(reinterpret_cast<char*>(clusters.data()), nClusters * sizeof(ClusterStruct));

  if (nTracksAtVtx > 0) {

    // fill the internal track structure based on the provided tracks at vertex
    tracks.reserve(nTracksAtVtx);
    for (const auto& trackAtVtx : tracksAtVtx) {
      tracks.emplace_back();
      FillTrack(tracks.back(), &trackAtVtx, (nMCHTracks > 0) ? &(mchTracks[trackAtVtx.mchTrackIdx]) : nullptr, clusters);
    }

  } else {

    // fill the internal track structure based on the MCH tracks extrapolated to the provided vertex
    tracks.reserve(nMCHTracks);
    for (const auto& mchTrack : mchTracks) {
      tracks.emplace_back();
      FillTrack(tracks.back(), nullptr, &mchTrack, clusters);
    }
  }

  size = 3 * sizeof(int) + nTracksAtVtx * sizeof(TrackAtVtxStruct) + nMCHTracks * sizeof(T) + nClusters * sizeof(ClusterStruct);
  return event;
}

//_________________________________________________________________________________________________
template <typename T>
void FillTrack(TrackStruct& track, const TrackAtVtxStruct* trackAtVtx, const T* mchTrack, std::vector<ClusterStruct>& clusters)
{
  /// fill the internal track structure from the provided informations

  if (trackAtVtx) {
    track.paramAtVertex = trackAtVtx->paramAtVertex;
    track.dca = trackAtVtx->dca;
    track.rAbs = trackAtVtx->rAbs;
  }

  if (mchTrack) {
    track.chi2 = mchTrack->getChi2();
    track.paramAt1stCluster.x = mchTrack->getX();
    track.paramAt1stCluster.y = mchTrack->getY();
    track.paramAt1stCluster.z = mchTrack->getZ();
    track.paramAt1stCluster.px = mchTrack->getPx();
    track.paramAt1stCluster.py = mchTrack->getPy();
    track.paramAt1stCluster.pz = mchTrack->getPz();
    track.paramAt1stCluster.sign = mchTrack->getSign();
    track.clusters.insert(track.clusters.end(),
                          std::make_move_iterator(clusters.begin() + mchTrack->getFirstClusterIdx()),
                          std::make_move_iterator(clusters.begin() + mchTrack->getLastClusterIdx() + 1));
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
