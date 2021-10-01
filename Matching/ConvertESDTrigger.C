#include <limits>
#include <vector>

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
#include "AliMUONESDInterface.h"
#include "AliMUONVTrackReconstructor.h"
#include "AliMUONGeometryTransformer.h"
#include "AliMUONTriggerCircuit.h"
#include "AliMUONVTriggerStore.h"
#include "AliMUONVTriggerTrackStore.h"
#include "AliMUONLocalTrigger.h"
#include "AliMUONTriggerTrack.h"

#include "Field/MagneticField.h"

#include "DataFormatsMID/ROFRecord.h"
#include "DataFormatsMID/Track.h"

using namespace std;
using namespace o2::mid;

AliMUONTriggerCircuit* triggerCircuit = nullptr;

bool Prepare(int runNumber);
bool SetMagField();
void GetTriggers(AliESDEvent& esd, AliMUONVTriggerStore& trigStore);
void ConvertTriggers(const AliMUONVTriggerTrackStore& trigTrackStore, const AliMUONVTriggerStore& trigStore,
                     std::vector<Track>& midTracks);

//_________________________________________________________________________________________________
void ConvertESDTrigger(TString esdFileName, TString outFileName = "MIDTracks.root")
{
  /// convert ESD trigger info into O2 structures saved in a root file

  // open the ESD file
  TFile* esdFile = TFile::Open(esdFileName);
  if (!esdFile || !esdFile->IsOpen()) {
    Error("ConvertESDTrigger", "opening ESD file %s failed", esdFileName.Data());
    return;
  }
  AliESDEvent* esd = new AliESDEvent();
  TTree* tree = (TTree*) esdFile->Get("esdTree");
  if (!tree) {
    Error("ConvertESDTrigger", "no ESD tree found");
    return;
  }
  esd->ReadFromTree(tree);

  std::vector<Track> midTracks{};
  std::vector<ROFRecord> midROFs(1);

  // prepare output file
  TFile* outFile = TFile::Open(outFileName, "RECREATE");
  TTree* outTree = new TTree("midreco", "midreco");
  TBranch* outTracks = outTree->Branch("MIDTrack", &midTracks);
  TBranch* outROFs = outTree->Branch("MIDTrackROF", &midROFs);

  // prepare the reconstruction tools
  if (tree->GetEvent(0) <= 0 || !Prepare(esd->GetRunNumber())) {
    Error("ConvertESDTrigger", "failed to setup the reconstruction");
    return;
  }

  AliMUONVTriggerStore* trigStore = AliMUONESDInterface::NewTriggerStore();
  AliMUONVTriggerTrackStore* trigTrackStore = AliMUONESDInterface::NewTriggerTrackStore();

  int nevents = (int)tree->GetEntries();
  for (int iEvent = 0; iEvent < nevents; iEvent++) {

    // get the ESD event
    if (tree->GetEvent(iEvent) <= 0) {
      Error("ConvertESDTrigger", "no ESD object found for event %d", iEvent);
      return;
    }

    // get the trigger info
    trigStore->Clear();
    GetTriggers(*esd, *trigStore);

    // reconstruct the trigger tracks
    trigTrackStore->Clear();
    AliMUONESDInterface::GetTracker()->EventReconstructTrigger(*triggerCircuit, *trigStore, *trigTrackStore);

    // convert in O2 format
    midTracks.clear();
    ConvertTriggers(*trigTrackStore, *trigStore, midTracks);
    midROFs[0] = ROFRecord({0, static_cast<uint32_t>(iEvent)}, EventType::Standard, 0, trigTrackStore->GetSize());

    // fill the output file
    outTree->Fill();
  }

  outFile->cd();
  outTree->SetEntries();
  outTree->Write();
  outFile->Close();
  delete outFile;
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
void GetTriggers(AliESDEvent& esd, AliMUONVTriggerStore& trigStore)
{
  /// get the trigger info
  for (int iTrack = 0; iTrack < esd.GetNumberOfMuonTracks(); ++iTrack) {
    AliESDMuonTrack* esdTrack = esd.GetMuonTrack(iTrack);
    if (esdTrack->ContainTriggerData()) {
      AliMUONESDInterface::Add(*esdTrack, trigStore);
    }
  }
}

//_________________________________________________________________________________________________
void ConvertTriggers(const AliMUONVTriggerTrackStore& trigTrackStore, const AliMUONVTriggerStore& trigStore,
                     std::vector<Track>& midTracks)
{
  /// convert MUON trigger tracks in O2 format

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
  }
}
