#include <fstream>
#include <iostream>
#include <vector>

#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TGeoGlobalMagField.h>

#include "AliCDBManager.h"
#include "AliGRPManager.h"
#include "AliGRPObject.h"

#include "AliESDEvent.h"
#include "AliESDMuonTrack.h"

#include "AliMUONCDB.h"
#include "AliMUONConstants.h"
#include "AliMUONTrack.h"
#include "AliMUONTrackParam.h"
#include "AliMUONTrackExtrap.h"
#include "AliMUONVCluster.h"
#include "AliMUONESDInterface.h"

#include "Field/MagneticField.h"

#include "DataFormatsMCH/TrackMCH.h"
#include "DataFormatsMCH/ClusterBlock.h"
#include "MCHBase/TrackBlock.h"

using namespace std;
using namespace o2::mch;

struct TrackAtVtxStruct {
  TrackParamStruct paramAtVertex{};
  double dca = 0.;
  double rAbs = 0.;
  int mchTrackIdx = 0;
};

bool PrepareRefitting(int runNumber);
bool SetMagField();
void ConvertTracks(AliESDEvent& esd, bool refit, std::vector<TrackAtVtxStruct>& o2TracksAtVtx,
                   std::vector<TrackMCH>& o2Tracks, std::vector<ClusterStruct>& o2Clusters);

//_________________________________________________________________________________________________
void ConvertESDTrack(TString esdFileName, TString outFileName = "AliESDTracks.out", bool refit = false, int LastEvent = -1)
{
  /// convert ESD tracks+clusters into O2 structures
  /// saved in a binary file with the following format:
  ///
  /// number of tracks at vertex in event 1 (= 0 if tracks refitted)
  /// number of MCH tracks in event 1
  /// number of associated clusters in event 1
  /// list of TrackAtVtxStruct (unless tracks are refitted)
  /// list of TrackMCH
  /// list of ClusterStruct
  /// number of tracks at vertex in event 2 (= 0)
  /// ...
  /// if the tracks are refitted, the extrapolation to vertex is left to another code
  /// if the tracks are not refitted, the parameters at MID as set to 0 as they cannot be computed

  // open the ESD file
  TFile* esdFile = TFile::Open(esdFileName);
  if (!esdFile || !esdFile->IsOpen()) {
    Error("ConvertESDTrack", "opening ESD file %s failed", esdFileName.Data());
    return;
  }
  AliESDEvent* esd = new AliESDEvent();
  TTree* tree = (TTree*) esdFile->Get("esdTree");
  if (!tree) {
    Error("ConvertESDTrack", "no ESD tree found");
    return;
  }
  esd->ReadFromTree(tree);
  
  // open output file
  ofstream outFile(outFileName.Data(),ios::out | ios::binary);
  if (!outFile.is_open()) {
    return;
  }

  // prepare refitting if needed
  if (refit && (tree->GetEvent(0) <= 0 || !PrepareRefitting(esd->GetRunNumber()))) {
    Error("ConvertESDTrack", "cannot refit tracks");
    return;
  }

  std::vector<TrackAtVtxStruct> o2TracksAtVtx{};
  std::vector<TrackMCH> o2Tracks{};
  std::vector<ClusterStruct> o2Clusters{};

  int nevents = (LastEvent >= 0) ? TMath::Min(LastEvent + 1, (int)tree->GetEntries()) : (int)tree->GetEntries();
  for (int iEvent = 0; iEvent < nevents; iEvent++) {

    // get the ESD event
    if (tree->GetEvent(iEvent) <= 0) {
      Error("ConvertESDTrack", "no ESD object found for event %d", iEvent);
      return;
    }

    // convert the (refitted) tracks in O2 format
    ConvertTracks(*esd, refit, o2TracksAtVtx, o2Tracks, o2Clusters);

    // write the tracks in the binary file
    int nTracksAtVtx = o2TracksAtVtx.size();
    outFile.write(reinterpret_cast<char*>(&nTracksAtVtx), sizeof(int));
    int nTracks = o2Tracks.size();
    outFile.write(reinterpret_cast<char*>(&nTracks), sizeof(int));
    int nClusters = o2Clusters.size();
    outFile.write(reinterpret_cast<char*>(&nClusters), sizeof(int));
    outFile.write(reinterpret_cast<char*>(o2TracksAtVtx.data()), o2TracksAtVtx.size() * sizeof(TrackAtVtxStruct));
    outFile.write(reinterpret_cast<char*>(o2Tracks.data()), o2Tracks.size() * sizeof(TrackMCH));
    outFile.write(reinterpret_cast<char*>(o2Clusters.data()), o2Clusters.size() * sizeof(ClusterStruct));
  }

  esdFile->Close();
  outFile.close();
}

//_________________________________________________________________________________________________
bool PrepareRefitting(int runNumber)
{
  /// prepare the tracker with access to recoParam (from OCDB) and magnetic field (from OCDB/GRP + O2 maps)

  AliCDBManager* man = AliCDBManager::Instance();
  if (gSystem->AccessPathName("OCDB.root", kFileExists) == 0) {
    man->SetDefaultStorage("local:///dev/null");
    man->SetSnapshotMode("OCDB.root");
  } else {
    man->SetDefaultStorage("local://./OCDB");
  }
  man->SetRun(runNumber);

  if (!SetMagField()) {
//  if (!AliMUONCDB::LoadField()) {
    return false;
  }

  AliMUONESDInterface::ResetTracker();

  return true;
}

//_________________________________________________________________________________________________
bool SetMagField()
{
  /// set the magnetic field using O2 maps and GRP info

  AliGRPManager grpMan;
  if (!grpMan.ReadGRPEntry()) {
    Error("SetMagField", "failed to load GRP Data from OCDB");
    return false;
  }

  const AliGRPObject* grpData = grpMan.GetGRPData();
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
  if (beamEnergy == AliGRPObject::GetInvalidFloat()) {
    Error("SetMagField", "GRP/GRP/Data entry:  missing value for the beam energy !");
    return false;
  }

  TString beamType = grpData->GetBeamType();
  if (beamType == AliGRPObject::GetInvalidString()) {
    Error("SetMagField", "GRP/GRP/Data entry:  missing value for the beam type !");
    return false;
  }

  Info("SetMagField", "l3Current = %f, diCurrent = %f", TMath::Abs(l3Current) * (l3Polarity ? -1 : 1),
       TMath::Abs(diCurrent) * (diPolarity ? -1 : 1));

  auto field = o2::field::MagneticField::createFieldMap(TMath::Abs(l3Current) * (l3Polarity ? -1 : 1),
                                                        TMath::Abs(diCurrent) * (diPolarity ? -1 : 1),
                                                        o2::field::MagneticField::kConvLHC, false, beamEnergy, beamType.Data(),
                                                        "$(O2_ROOT)/share/Common/maps/mfchebKGI_sym.root");
  TGeoGlobalMagField::Instance()->SetField(field);
  TGeoGlobalMagField::Instance()->Lock();

  return true;
}

//_________________________________________________________________________________________________
void ConvertTracks(AliESDEvent& esd, bool refit, std::vector<TrackAtVtxStruct>& o2TracksAtVtx,
                   std::vector<TrackMCH>& o2Tracks, std::vector<ClusterStruct>& o2Clusters)
{
  /// write the tracks in O2 format in the output file, refitting them before if requested
  /// if the tracks are refitted, the extrapolation to vertex is left to another code

  static AliMUONTrack muonTrack{};

  o2TracksAtVtx.clear();
  o2Tracks.clear();
  o2Clusters.clear();

  for (int iTrack = 0; iTrack < esd.GetNumberOfMuonTracks(); ++iTrack) {

    AliESDMuonTrack* esdTrack = esd.GetMuonTrack(iTrack);
    if (!esdTrack->ContainTrackerData()) {
      continue;
    }

    if (!refit) {
      double dcaX = esdTrack->GetNonBendingCoorAtDCA() - esdTrack->GetNonBendingCoor();
      double dcaY = esdTrack->GetBendingCoorAtDCA() - esdTrack->GetBendingCoor();
      o2TracksAtVtx.push_back({{esdTrack->GetNonBendingCoor(), esdTrack->GetBendingCoor(), esdTrack->GetZ(),
                                esdTrack->Px(), esdTrack->Py(), esdTrack->Pz(), esdTrack->Charge()},
                               TMath::Sqrt(dcaX * dcaX + dcaY * dcaY),
                               esdTrack->GetRAtAbsorberEnd(),
                               static_cast<int>(o2Tracks.size())});
    }

    AliMUONESDInterface::ESDToMUON(*esdTrack, muonTrack, refit);
    AliMUONTrackParam paramAtMID{};
    if (refit) {
      paramAtMID = *static_cast<AliMUONTrackParam*>(muonTrack.GetTrackParamAtCluster()->Last());
      AliMUONTrackExtrap::ExtrapToZCov(&paramAtMID, AliMUONConstants::MuonFilterZEnd());
      AliMUONTrackExtrap::AddMCSEffect(&paramAtMID, AliMUONConstants::MuonFilterZEnd() - AliMUONConstants::MuonFilterZBeg(), AliMUONConstants::MuonFilterX0());
      AliMUONTrackExtrap::ExtrapToZCov(&paramAtMID, AliMUONConstants::DefaultChamberZ(AliMUONConstants::NTrackingCh()));
    }

    AliMUONTrackParam* param = static_cast<AliMUONTrackParam*>(muonTrack.GetTrackParamAtCluster()->First());
    o2Tracks.emplace_back(param->GetZ(), param->GetParameters(), param->GetCovariances(),
                          muonTrack.GetGlobalChi2(), o2Clusters.size(), muonTrack.GetNClusters(),
                          paramAtMID.GetZ(), paramAtMID.GetParameters(), paramAtMID.GetCovariances());

    for (int iCl = 0; iCl < muonTrack.GetNClusters(); ++iCl) {
      AliMUONVCluster* cluster = static_cast<AliMUONTrackParam*>(muonTrack.GetTrackParamAtCluster()->UncheckedAt(iCl))->GetClusterPtr();
      o2Clusters.push_back({static_cast<float>(cluster->GetX()), static_cast<float>(cluster->GetY()),
                            static_cast<float>(cluster->GetZ()), static_cast<float>(cluster->GetErrX()),
                            static_cast<float>(cluster->GetErrY()), cluster->GetUniqueID(), 0, 0});
    }
  }
}
