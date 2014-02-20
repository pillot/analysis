#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TFile.h"
#include "TTree.h"
#include <Riostream.h>
#include <TGeoManager.h>
#include <TROOT.h>

// STEER includes
#include "AliCDBManager.h"
#include "AliESDEvent.h"
#include "AliESDMuonTrack.h"

// MUON includes
#include "AliMUONCDB.h"
#include "AliMUONConstants.h"
#include "AliMUONTrackParam.h"
#include "AliMUONTrackExtrap.h"
#include "AliMUONESDInterface.h"
#endif

void CheckExtrapToVertex(Int_t nevents = 1)
{
  
  //Reset ROOT and connect tree file
  gROOT->Reset();
  
  // open the ESD file
  TFile* esdFile = TFile::Open("AliESDs.root");
  if (!esdFile || !esdFile->IsOpen()) {
    Error("CheckExtrapToVertex", "opening ESD file %s failed", "AliESDs.root");
    return;
  }
  AliESDEvent* esd = new AliESDEvent();
  TTree* tree = (TTree*) esdFile->Get("esdTree");
  if (!tree) {
    Error("CheckExtrapToVertex", "no ESD tree found");
    return;
  }
  esd->ReadFromTree(tree);
  
  // get run number
  if (tree->GetEvent(0) <= 0) {
    Error("CheckExtrapToVertex", "no ESD object found for event 0");
    return;
  }
  Int_t runNumber = esd->GetRunNumber();
  
  // Import TGeo geometry (needed by AliMUONTrackExtrap::ExtrapToVertex)
  if (!gGeoManager) {
    TGeoManager::Import("geometry.root");
    if (!gGeoManager) {
      Error("CheckExtrapToVertex", "getting geometry from file %s failed", "generated/galice.root");
      return;
    }
  }
  
  // load necessary data from OCDB
  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  AliCDBManager::Instance()->SetRun(runNumber);
  if (!AliMUONCDB::LoadField()) return;
  
  // set the magnetic field for track extrapolations
  AliMUONTrackExtrap::SetField();
  
  // Loop over events
  nevents = TMath::Min(nevents,(Int_t)tree->GetEntries());
  for (Int_t iEvent = 0; iEvent < nevents; iEvent++) {
    
    // get the event summary data
    if (tree->GetEvent(iEvent) <= 0) {
      Error("CheckExtrapToVertex", "no ESD object found for event %d", iEvent);
      return;
    }
    
    Int_t nTracks = (Int_t)esd->GetNumberOfMuonTracks() ; 
    
    // loop over all reconstructed tracks (also first track of combination)
    for (Int_t iTrack = 0; iTrack <  nTracks;  iTrack++) {
      
      // extrapolate to vertex
      AliMUONTrackParam trackParam;
      AliMUONESDInterface::GetParamAtFirstCluster(*(esd->GetMuonTrack(iTrack)), trackParam);
      trackParam.Print("FULL");
      AliMUONESDInterface::GetParamCov(*(esd->GetMuonTrack(iTrack)), trackParam);
      cout<<"extrap to vertex uncorrected"<<endl;
      AliMUONTrackExtrap::ExtrapToVertexUncorrected(&trackParam, 0.);
      //trackParam.GetCovariances().Print();
      AliMUONTrackExtrap::ExtrapToZCov(&trackParam, AliMUONConstants::AbsZEnd()+415.);
      trackParam.GetCovariances().Print();
      //AliMUONTrackExtrap::ExtrapToZCov(&trackParam, AliMUONConstants::AbsZEnd()-415.);
      //trackParam.GetCovariances().Print();
      
      AliMUONESDInterface::GetParamAtFirstCluster(*(esd->GetMuonTrack(iTrack)), trackParam);
      AliMUONESDInterface::GetParamCov(*(esd->GetMuonTrack(iTrack)), trackParam);
      cout<<"extrap to vertex"<<endl;
      AliMUONTrackExtrap::ExtrapToVertex(&trackParam, 0., 0., 0., 0., 0.);
      //trackParam.GetCovariances().Print();
      AliMUONTrackExtrap::ExtrapToZCov(&trackParam, AliMUONConstants::AbsZEnd()+415.);
      trackParam.GetCovariances().Print();
      //AliMUONTrackExtrap::ExtrapToZCov(&trackParam, AliMUONConstants::AbsZEnd()-415.);
      //trackParam.GetCovariances().Print();
      
    }
    
  }
  
}

