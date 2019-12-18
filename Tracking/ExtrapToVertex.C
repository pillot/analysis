#include <fstream>
#include <iostream>
#include <vector>

#include <TMath.h>
#include <TGeoManager.h>
#include <TGeoGlobalMagField.h>

#include "AliCDBManager.h"
#include "AliGRPManager.h"
#include "AliGRPObject.h"
#include "AliMUONCDB.h"
#include "AliMUONTrackExtrap.h"
#include "AliMUONTrackParam.h"

#include "DetectorsBase/GeometryManager.h"
#include "Field/MagneticField.h"

#include "MCHBase/ClusterBlock.h"
#include "MCHBase/TrackBlock.h"

using namespace std;
using namespace o2::mch;

//_________________________________________________________________________________________________
struct VertexStruct {
  double x;
  double y;
  double z;
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
bool LoadOCDB();
bool SetMagField();
int ReadNextEvent(ifstream& inFile, VertexStruct& vertices);
int ReadNextEvent(ifstream& inFile, int version, std::vector<TrackStruct>& tracks);
void ReadTrack(ifstream& inFile, TrackStruct& track, int version);
void ExtrapTrackToVertex(TrackStruct& track, VertexStruct& vertex);
void WriteTracks(std::vector<TrackStruct>& tracks, ofstream& outFile);

//_________________________________________________________________________________________________
void ExtrapToVertex(string trackFileName, int versionFile, string vtxFileName)
{
  /// extrapolate the tracks in trackFileName to the corresponding vertex with AliRoot code
  /// and save them with additional parameters at vertex in a new binary file

  if (versionFile > 2) {
    cout << "tracks are already extrapolated to vertex. Not doing it again" << endl;
    return;
  }

  // prepare track extrapolation
  if (!LoadOCDB()) {
    cout << "fail loading OCDB objects for track extrapolation" << endl;
    return;
  }

  // open files
  ifstream inFileTrack(trackFileName, ios::binary);
  ifstream inFileVtx(vtxFileName, ios::binary);
  if (!inFileTrack.is_open() || !inFileVtx.is_open()) {
    cout << "fail opening files" << endl;
    return;
  }
  ofstream outFileTrack("O2Tracks.vtx", ios::out | ios::binary);
  if (!outFileTrack.is_open()) {
    return;
  }

  std::vector<TrackStruct> tracks{};
  VertexStruct vertex{};
  while (true) {

    // get vertex and tracks
    int event1 = ReadNextEvent(inFileTrack, versionFile, tracks);
    int event2 = ReadNextEvent(inFileVtx, vertex);

    if (event1 < 0 && event2 < 0) {
      // reaching end of both files
      break;
    }
    if (event1 != event2) {
      cout << "inconsistent files" << endl;
      break;
    }

    // extrapolate tracks to vertex
    for (auto& track : tracks) {
      ExtrapTrackToVertex(track, vertex);
    }

    // save tracks with parameters at vertex
    // write the tracks in the second binary file
    outFileTrack.write((char*)&event1, sizeof(int));
    WriteTracks(tracks, outFileTrack);
  }

  inFileTrack.close();
  inFileVtx.close();
  outFileTrack.close();
}

//_________________________________________________________________________________________________
bool LoadOCDB()
{
  /// access OCDB and prepare track extrapolation to vertex

  AliCDBManager* man = AliCDBManager::Instance();
  if (gSystem->AccessPathName("OCDB.root", kFileExists) == 0) {
    man->SetDefaultStorage("local:///dev/null");
    man->SetSnapshotMode("OCDB.root");
  } else {
    man->SetDefaultStorage("local://./OCDB");
  }
  man->SetRun(295584);

  if (!SetMagField()) {
//  if (!AliMUONCDB::LoadField()) {
    return false;
  }
  AliMUONTrackExtrap::SetField();

  o2::base::GeometryManager::loadGeometry("O2geometry.root");
//  o2::base::GeometryManager::loadGeometry("AliRootgeometry.root", "ALICE");
  if (!gGeoManager) {
    return false;
  }

  return true;
}

//_________________________________________________________________________________________________
bool SetMagField()
{
  /// Set the magnetic field using O2 maps and GRP info
  
  AliGRPManager grpMan;
  if (!grpMan.ReadGRPEntry()) {
    Error("SetMagField", "failed to load GRP Data from OCDB");
    return false;
  }
  
  const AliGRPObject *grpData = grpMan.GetGRPData();
  if (!grpData) {
    Error("SetMagField", "GRP Data is not loaded");
    return false;
  }
  
  float l3Current = grpData->GetL3Current((AliGRPObject::Stats)0);
  if (l3Current == AliGRPObject::GetInvalidFloat()) {
    Error("SetMagField", "GRP/GRP/Data entry:  missing value for the L3 current !");
    return false;
  }

  char l3Polarity = grpData->GetL3Polarity();
  if (l3Polarity == AliGRPObject::GetInvalidChar()) {
    Error("SetMagField", "GRP/GRP/Data entry:  missing value for the L3 polarity !");
    return false;
  }
  
  float diCurrent = grpData->GetDipoleCurrent((AliGRPObject::Stats)0);
  if (diCurrent == AliGRPObject::GetInvalidFloat()) {
    Error("SetMagField", "GRP/GRP/Data entry:  missing value for the dipole current !");
    return false;
  }

  char diPolarity = grpData->GetDipolePolarity();
  if (diPolarity == AliGRPObject::GetInvalidChar()) {
    Error("SetMagField", "GRP/GRP/Data entry:  missing value for the dipole polarity !");
    return false;
  }
  
  float beamEnergy = grpData->GetBeamEnergy();
  if (beamEnergy==AliGRPObject::GetInvalidFloat()) {
    Error("SetMagField", "GRP/GRP/Data entry:  missing value for the beam energy !");
    return false;
  }

  TString beamType = grpData->GetBeamType();
  if (beamType==AliGRPObject::GetInvalidString()) {
    Error("SetMagField", "GRP/GRP/Data entry:  missing value for the beam type !");
    return false;
  }
  
  Info("SetMagField", "l3Current = %f, diCurrent = %f", TMath::Abs(l3Current) * (l3Polarity ? -1:1),
                                                        TMath::Abs(diCurrent) * (diPolarity ? -1:1));
  
  auto field = o2::field::MagneticField::createFieldMap(TMath::Abs(l3Current) * (l3Polarity ? -1:1),
                                                        TMath::Abs(diCurrent) * (diPolarity ? -1:1),
                                                        o2::field::MagneticField::kConvLHC, false, beamEnergy, beamType.Data(),
                                                        "$(O2_ROOT)/share/Common/maps/mfchebKGI_sym.root");
  TGeoGlobalMagField::Instance()->SetField(field);
  TGeoGlobalMagField::Instance()->Lock();

  return true;
}

//_________________________________________________________________________________________________
int ReadNextEvent(ifstream& inFile, VertexStruct& vertex)
{
  /// read the next vertex in the input file
  /// return the event number or -1 when reaching end of file

  int event(-1);
  if (!inFile.read(reinterpret_cast<char*>(&event), sizeof(int))) {
    // reaching end of file
    return -1;
  }

  // read the vertex
  inFile.read(reinterpret_cast<char*>(&vertex), sizeof(VertexStruct));

  return event;
}

//_________________________________________________________________________________________________
int ReadNextEvent(ifstream& inFile, int version, std::vector<TrackStruct>& tracks)
{
  /// read the next event in the input file

  int event(-1);
  if (!inFile.read(reinterpret_cast<char*>(&event), sizeof(int))) {
    // reaching end of file
    return -1;
  }
  
  int size(0);
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
void ExtrapTrackToVertex(TrackStruct& track, VertexStruct& vertex)
{
  /// compute the track parameters at vertex, at DCA and at the end of the absorber

  // convert parameters at first cluster in MUON format
  AliMUONTrackParam trackParam;
  trackParam.SetNonBendingCoor(track.paramAt1stCluster.x);
  trackParam.SetBendingCoor(track.paramAt1stCluster.y);
  trackParam.SetZ(track.paramAt1stCluster.z);
  trackParam.SetNonBendingSlope(track.paramAt1stCluster.px / track.paramAt1stCluster.pz);
  trackParam.SetBendingSlope(track.paramAt1stCluster.py / track.paramAt1stCluster.pz);
  trackParam.SetInverseBendingMomentum(track.paramAt1stCluster.sign / TMath::Sqrt(track.paramAt1stCluster.py * track.paramAt1stCluster.py + track.paramAt1stCluster.pz * track.paramAt1stCluster.pz));

  // extrapolate to vertex
  AliMUONTrackParam trackParamAtVertex(trackParam);
  AliMUONTrackExtrap::ExtrapToVertex(&trackParamAtVertex, vertex.x, vertex.y, vertex.z, 0., 0.);
  track.paramAtVertex.x = trackParamAtVertex.GetNonBendingCoor();
  track.paramAtVertex.y = trackParamAtVertex.GetBendingCoor();
  track.paramAtVertex.z = trackParamAtVertex.GetZ();
  track.paramAtVertex.px = trackParamAtVertex.Px();
  track.paramAtVertex.py = trackParamAtVertex.Py();
  track.paramAtVertex.pz = trackParamAtVertex.Pz();
  track.paramAtVertex.sign = trackParamAtVertex.GetCharge();

  // extrapolate to DCA
  AliMUONTrackParam trackParamAtDCA(trackParam);
  AliMUONTrackExtrap::ExtrapToVertexWithoutBranson(&trackParamAtDCA, vertex.z);
  double dcaX = trackParamAtDCA.GetNonBendingCoor() - vertex.x;
  double dcaY = trackParamAtDCA.GetBendingCoor() - vertex.y;
  track.dca = TMath::Sqrt(dcaX*dcaX + dcaY*dcaY);
 
  // extrapolate to the end of the absorber
  AliMUONTrackExtrap::ExtrapToZ(&trackParam, -505.);
  double xAbs = trackParam.GetNonBendingCoor();
  double yAbs = trackParam.GetBendingCoor();
  track.rAbs = TMath::Sqrt(xAbs*xAbs + yAbs*yAbs);
}

//_________________________________________________________________________________________________
void WriteTracks(std::vector<TrackStruct>& tracks, ofstream& outFile)
{
  /// write the tracks in the output file

  int size(0); // total number of bytes requested to store the event, excluding event number and total size
  int nTracks(tracks.size());
  if (nTracks == 0) {
    outFile.write((char*)&size, sizeof(int));
    return;
  } else {
    size = sizeof(int);
  }
  for (const auto& track : tracks) {
    size += 2 * sizeof(TrackParamStruct) + 3 * sizeof(double) + sizeof(int) + track.clusters.size() * sizeof(ClusterStruct);
  }
  outFile.write((char*)&size, sizeof(int));

  outFile.write((char*)&nTracks, sizeof(int));

  for (const auto& track : tracks) {
    outFile.write((char*)&(track.paramAtVertex), sizeof(TrackParamStruct));
    outFile.write((char*)&(track.dca), sizeof(double));
    outFile.write((char*)&(track.rAbs), sizeof(double));
    outFile.write((char*)&(track.paramAt1stCluster), sizeof(TrackParamStruct));
    outFile.write((char*)&(track.chi2), sizeof(double));
    int nClusters = track.clusters.size();
    outFile.write((char*)&(nClusters), sizeof(int));
    for (const auto& cluster : track.clusters) {
      outFile.write((char*)&cluster, sizeof(ClusterStruct));
    }
  }
}
