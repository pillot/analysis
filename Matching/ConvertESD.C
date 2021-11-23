#include <limits>
#include <vector>
#include <map>

#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TGeoGlobalMagField.h>
#include <TMatrixD.h>

#include "AliCDBManager.h"
#include "AliGRPManager.h"
#include "AliGRPObject.h"
#include "AliGeomManager.h"

#include "AliESDEvent.h"
#include "AliESDMuonTrack.h"

#include "AliMpCDB.h"

#include "AliMUONCDB.h"
#include "AliMUONConstants.h"
#include "AliMUONESDInterface.h"
#include "AliMUONVTrackReconstructor.h"
#include "AliMUONVTrackStore.h"
#include "AliMUONTrack.h"
#include "AliMUONTrackParam.h"
#include "AliMUONTrackExtrap.h"
#include "AliMUONVCluster.h"
#include "AliMUONGeometryTransformer.h"
#include "AliMUONTriggerCircuit.h"
#include "AliMUONVTriggerStore.h"
#include "AliMUONVTriggerTrackStore.h"
#include "AliMUONLocalTrigger.h"
#include "AliMUONTriggerTrack.h"
#include "AliMUONVDigitStore.h"
#include "AliMUONTrackHitPattern.h"

#include "Field/MagneticField.h"

#include "ReconstructionDataFormats/TrackMCHMID.h"
#include "DataFormatsMCH/ROFRecord.h"
#include "DataFormatsMCH/TrackMCH.h"
#include "DataFormatsMID/ROFRecord.h"
#include "DataFormatsMID/Track.h"
#include "DataFormatsMCH/Cluster.h"

using namespace std;
using namespace o2;

AliMUONTriggerCircuit* triggerCircuit = nullptr;
AliMUONTrackHitPattern* trackHitPattern = nullptr;

bool Prepare(int runNumber);
bool SetMagField();
void LoadEvent(AliESDEvent& esd, AliMUONVTrackStore& trackStore, AliMUONVTriggerStore& trigStore);
void MatchTracks(AliMUONVTrackStore& trackStore, std::map<uint32_t, AliMUONTrackParam>& paramsAtMID,
                 const AliMUONVTriggerTrackStore& trigTrackStore, const AliMUONVTriggerStore& trigStore);
void ConvertTracks(const AliMUONVTrackStore& trackStore, std::map<uint32_t, AliMUONTrackParam>& paramsAtMID,
                   std::vector<mch::TrackMCH>& mchTracks, std::vector<mch::Cluster>& mchClusters,
                   std::map<uint32_t, int>& mchTracksIdx);
void ConvertTriggers(const AliMUONVTriggerTrackStore& trigTrackStore, const AliMUONVTriggerStore& trigStore,
                     std::vector<mid::Track>& midTracks, std::map<uint32_t, int>& midTracksIdx);
void ConvertMatched(const AliMUONVTrackStore& trackStore, const AliMUONVTriggerStore& trigStore,
                    const std::map<uint32_t, int>& mchTracksIdx, const std::map<uint32_t, int>& midTracksIdx,
                    const InteractionRecord& midIR, std::vector<dataformats::TrackMCHMID>& muonTracks);

//_________________________________________________________________________________________________
void ConvertESD(TString esdFileName, int nROFsPerTF)
{
  /// convert ESD info into O2 structures saved in a root files

  // open the ESD file
  TFile* esdFile = TFile::Open(esdFileName);
  if (!esdFile || !esdFile->IsOpen()) {
    Error("ConvertESD", "opening ESD file %s failed", esdFileName.Data());
    return;
  }
  AliESDEvent* esd = new AliESDEvent();
  TTree* tree = (TTree*) esdFile->Get("esdTree");
  if (!tree) {
    Error("ConvertESD", "no ESD tree found");
    return;
  }
  esd->ReadFromTree(tree);

  // prepare output files
  TFile* outFileMCH = TFile::Open("ESDMCHTracks.root", "RECREATE");
  TTree* outTreeMCH = new TTree("o2sim", "Tree MCH Standalone Tracks");
  std::vector<mch::TrackMCH> mchTracks{};
  outTreeMCH->Branch("tracks", &mchTracks);
  std::vector<mch::ROFRecord> mchROFs{};
  outTreeMCH->Branch("trackrofs", &mchROFs);
  std::vector<mch::Cluster> mchClusters{};
  outTreeMCH->Branch("trackclusters", &mchClusters);
  TFile* outFileMID = TFile::Open("ESDMIDTracks.root", "RECREATE");
  TTree* outTreeMID = new TTree("midreco", "midreco");
  std::vector<mid::Track> midTracks{};
  outTreeMID->Branch("MIDTrack", &midTracks);
  std::vector<mid::ROFRecord> midROFs{};
  outTreeMID->Branch("MIDTrackROF", &midROFs);
  TFile* outFileMUON = TFile::Open("ESDMUONTracks.root", "RECREATE");
  TTree* outTreeMUON = new TTree("o2sim", "Tree Matched MCH-MID Tracks");
  std::vector<dataformats::TrackMCHMID> muonTracks{};
  outTreeMUON->Branch("tracks", &muonTracks);

  // prepare the reconstruction tools
  if (tree->GetEvent(0) <= 0 || !Prepare(esd->GetRunNumber())) {
    Error("ConvertESD", "failed to setup the reconstruction");
    return;
  }

  AliMUONVTrackStore* trackStore = AliMUONESDInterface::NewTrackStore();
  AliMUONVTriggerStore* trigStore = AliMUONESDInterface::NewTriggerStore();
  AliMUONVTriggerTrackStore* trigTrackStore = AliMUONESDInterface::NewTriggerTrackStore();
  std::map<uint32_t, AliMUONTrackParam> paramsAtMID{};
  std::map<uint32_t, int> mchTracksIdx{};
  std::map<uint32_t, int> midTracksIdx{};

  int nevents = (int)tree->GetEntries();
  int nROFsInCurrentTF = 0;
  for (int iEvent = 0; iEvent < nevents; iEvent++) {

    printf("processing event... %d %%\r", 100 * (iEvent + 1) / nevents);

    // get the ESD event
    if (tree->GetEvent(iEvent) <= 0) {
      Error("ConvertESD", "no ESD object found for event %d", iEvent);
      return;
    }

    // get the tracker and trigger info in MUON format
    LoadEvent(*esd, *trackStore, *trigStore);

    // convert tracker tracks in O2 format
    int firstMCHIdx = mchTracks.size();
    ConvertTracks(*trackStore, paramsAtMID, mchTracks, mchClusters, mchTracksIdx);
    mchROFs.emplace_back(InteractionRecord{0, static_cast<uint32_t>(iEvent)}, firstMCHIdx, trackStore->GetSize());

    // reconstruct the trigger tracks
    trigTrackStore->Clear();
    AliMUONESDInterface::GetTracker()->EventReconstructTrigger(*triggerCircuit, *trigStore, *trigTrackStore);

    // convert them in O2 format
    size_t firstMIDIdx = midTracks.size();
    ConvertTriggers(*trigTrackStore, *trigStore, midTracks, midTracksIdx);
    midROFs.emplace_back(InteractionRecord{0, static_cast<uint32_t>(iEvent)}, mid::EventType::Standard, firstMIDIdx,
                         static_cast<size_t>(trigTrackStore->GetSize()));

    // match tracker and trigger tracks
    MatchTracks(*trackStore, paramsAtMID, *trigTrackStore, *trigStore);

    // convert them in O2 format
    ConvertMatched(*trackStore, *trigStore, mchTracksIdx, midTracksIdx, midROFs.back().interactionRecord, muonTracks);

    // fill the output files if the TF is full and reset the containers
    if (++nROFsInCurrentTF == nROFsPerTF) {
      outTreeMCH->Fill();
      outTreeMID->Fill();
      outTreeMUON->Fill();
      mchROFs.clear();
      mchTracks.clear();
      mchClusters.clear();
      midROFs.clear();
      midTracks.clear();
      muonTracks.clear();
      nROFsInCurrentTF = 0;
    }
  }

  // fill the output files with the last TF if not already done
  if (nROFsInCurrentTF > 0) {
    outTreeMCH->Fill();
    outTreeMID->Fill();
    outTreeMUON->Fill();
  }
  printf("processing event... done \n");

  outFileMCH->cd();
  outTreeMCH->SetEntries();
  outTreeMCH->Write();
  outFileMCH->Close();
  delete outFileMCH;
  outFileMID->cd();
  outTreeMID->SetEntries();
  outTreeMID->Write();
  outFileMID->Close();
  delete outFileMID;
  outFileMUON->cd();
  outTreeMUON->SetEntries();
  outTreeMUON->Write();
  outFileMUON->Close();
  delete outFileMUON;
  esdFile->Close();
}

//_________________________________________________________________________________________________
bool Prepare(int runNumber)
{
  /// prepare the reconstruction tools with necessary objects:
  /// - recoParam (from OCDB)
  /// - magnetic field (from OCDB/GRP + O2 maps)
  /// - geometry (from OCDB)
  /// - mapping (from OCDB)

  // setup the OCDB
  AliCDBManager* man = AliCDBManager::Instance();
  if (gSystem->AccessPathName("OCDB.root", kFileExists) == 0) {
    man->SetDefaultStorage("local:///dev/null");
    man->SetSnapshotMode("OCDB.root");
  } else {
    man->SetDefaultStorage("local://./OCDB");
  }
  man->SetRun(runNumber);

  // load the magnetic field
  if (!SetMagField()) {
//  if (!AliMUONCDB::LoadField()) {
    return false;
  }

  // load the geometry
  AliGeomManager::LoadGeometry();
  if (!AliGeomManager::GetGeometry() || !AliGeomManager::ApplyAlignObjsFromCDB("MUON")) {
    return false;
  }

  // load the mapping
  if (!AliMpCDB::LoadDDLStore()) {
    return false;
  }


  // prepare the trigger track reconstructor
  AliMUONESDInterface::ResetTracker();
  AliMUONGeometryTransformer* transformer = new AliMUONGeometryTransformer();
  transformer->LoadGeometryData();
  triggerCircuit = new AliMUONTriggerCircuit(transformer);

  AliMUONVDigitStore* digitStore = AliMUONESDInterface::NewDigitStore();
  trackHitPattern = new AliMUONTrackHitPattern(AliMUONESDInterface::GetTracker()->GetRecoParam(), *transformer, *digitStore, nullptr);
  delete digitStore; // assume the digitStore is never used when we use AliMUONTrackHitPattern

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
void LoadEvent(AliESDEvent& esd, AliMUONVTrackStore& trackStore, AliMUONVTriggerStore& trigStore)
{
  /// get the tracker and trigger info in MUON format
  trackStore.Clear();
  trigStore.Clear();
  for (int iTrack = 0; iTrack < esd.GetNumberOfMuonTracks(); ++iTrack) {
    AliESDMuonTrack* esdTrack = esd.GetMuonTrack(iTrack);
    if (esdTrack->ContainTrackerData()) {
      AliMUONESDInterface::Add(*esdTrack, trackStore, true);
    }
    if (esdTrack->ContainTriggerData()) {
      AliMUONESDInterface::Add(*esdTrack, trigStore);
    }
  }
}

//_________________________________________________________________________________________________
void MatchTracks(AliMUONVTrackStore& trackStore, std::map<uint32_t, AliMUONTrackParam>& paramsAtMID,
                 const AliMUONVTriggerTrackStore& trigTrackStore, const AliMUONVTriggerStore& trigStore)
{
  /// match the tracker and trigger tracks
  AliMUONTrack* track = nullptr;
  TIter next(trackStore.CreateIterator());
  while ((track = static_cast<AliMUONTrack*>(next()))) {
    trackHitPattern->MatchTriggerTrack(track, paramsAtMID.find(track->GetUniqueID())->second, trigTrackStore, trigStore);
  }
}

//_________________________________________________________________________________________________
void ConvertTracks(const AliMUONVTrackStore& trackStore, std::map<uint32_t, AliMUONTrackParam>& paramsAtMID,
                   std::vector<mch::TrackMCH>& mchTracks, std::vector<mch::Cluster>& mchClusters,
                   std::map<uint32_t, int>& mchTracksIdx)
{
  /// convert MUON tracks in O2 format

  paramsAtMID.clear();
  mchTracksIdx.clear();
  int mchTrackIdx = mchTracks.size() - 1;

  AliMUONTrack* track = nullptr;
  TIter next(trackStore.CreateIterator());
  while ((track = static_cast<AliMUONTrack*>(next()))) {

    AliMUONTrackParam paramAtMID(*static_cast<AliMUONTrackParam*>(track->GetTrackParamAtCluster()->Last()));
    AliMUONTrackExtrap::ExtrapToZCov(&paramAtMID, AliMUONConstants::MuonFilterZEnd());
    AliMUONTrackExtrap::AddMCSEffect(&paramAtMID, AliMUONConstants::MuonFilterZEnd() - AliMUONConstants::MuonFilterZBeg(), AliMUONConstants::MuonFilterX0());
    AliMUONTrackExtrap::ExtrapToZCov(&paramAtMID, AliMUONConstants::DefaultChamberZ(AliMUONConstants::NTrackingCh()));
    paramsAtMID.emplace(track->GetUniqueID(), paramAtMID);

    AliMUONTrackParam* param = static_cast<AliMUONTrackParam*>(track->GetTrackParamAtCluster()->First());
    mchTracks.emplace_back(param->GetZ(), param->GetParameters(), param->GetCovariances(),
                           track->GetGlobalChi2(), mchClusters.size(), track->GetNClusters(),
                           paramAtMID.GetZ(), paramAtMID.GetParameters(), paramAtMID.GetCovariances());

    for (int iCl = 0; iCl < track->GetNClusters(); ++iCl) {
      AliMUONVCluster* cluster = static_cast<AliMUONTrackParam*>(track->GetTrackParamAtCluster()->UncheckedAt(iCl))->GetClusterPtr();
      mchClusters.push_back({static_cast<float>(cluster->GetX()), static_cast<float>(cluster->GetY()),
                             static_cast<float>(cluster->GetZ()), static_cast<float>(cluster->GetErrX()),
                             static_cast<float>(cluster->GetErrY()), cluster->GetUniqueID(), 0, 0});
    }

    mchTracksIdx[track->GetUniqueID()] = ++mchTrackIdx;
  }
}

//_________________________________________________________________________________________________
void ConvertTriggers(const AliMUONVTriggerTrackStore& trigTrackStore, const AliMUONVTriggerStore& trigStore,
                     std::vector<mid::Track>& midTracks, std::map<uint32_t, int>& midTracksIdx)
{
  /// convert MUON trigger tracks in O2 format

  midTracksIdx.clear();
  int midTrackIdx = midTracks.size() - 1;

  AliMUONTriggerTrack* triggerTrack = nullptr;
  TIter next(trigTrackStore.CreateIterator());
  while ((triggerTrack = static_cast<AliMUONTriggerTrack*>(next()))) {

    auto& midTrack = midTracks.emplace_back();

    midTrack.setPosition(triggerTrack->GetX11(), triggerTrack->GetY11(), triggerTrack->GetZ11());

    midTrack.setDirection(triggerTrack->GetSlopeX(), triggerTrack->GetSlopeY(), 1.);

    const TMatrixD& trigCov = triggerTrack->GetCovariances();
    midTrack.setCovarianceParameters(trigCov(0, 0), trigCov(1, 1), std::numeric_limits<float>::max(),
                                     trigCov(2, 2), 0., trigCov(1, 2));

    AliMUONLocalTrigger* locTrg = trigStore.FindLocal(triggerTrack->GetLoTrgNum());
    if (locTrg->LoHpt() > 0) {
      midTrack.setChi2(1.);
    } else if (locTrg->LoLpt() > 0) {
      midTrack.setChi2(2.);
    } else {
      midTrack.setChi2(3.);
    }

    midTracksIdx[triggerTrack->GetUniqueID()] = ++midTrackIdx;
  }
}

//_________________________________________________________________________________________________
void ConvertMatched(const AliMUONVTrackStore& trackStore, const AliMUONVTriggerStore& trigStore,
                    const std::map<uint32_t, int>& mchTracksIdx, const std::map<uint32_t, int>& midTracksIdx,
                    const InteractionRecord& midIR, std::vector<dataformats::TrackMCHMID>& muonTracks)
{
  /// convert MUON matched tracks in O2 format
  AliMUONTrack* track = nullptr;
  TIter next(trackStore.CreateIterator());
  while ((track = static_cast<AliMUONTrack*>(next()))) {
    if (track->GetMatchTrigger() > 0) {
      muonTracks.emplace_back(mchTracksIdx.find(track->GetUniqueID())->second,
                              midTracksIdx.find(trigStore.FindLocal(track->LoCircuit())->GetUniqueID())->second,
                              midIR, track->GetChi2MatchTrigger());
    }
  }
}
