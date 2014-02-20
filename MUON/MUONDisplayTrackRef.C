#if !defined(__CINT__) || defined(__MAKECINT__)
// ROOT includes
#include <Riostream.h>
#include <TROOT.h>
#include <TFile.h>
#include <TSystem.h>

// STEER includes
#include "AliLog.h"
#include "AliMCEventHandler.h"
#include "AliCDBManager.h"
#include "AliGeomManager.h"
#include "AliMCEvent.h"

// MUON includes
#include "AliMUONCDB.h"
#include "AliMUONRecoCheck.h"
#include "AliMUONVTrackStore.h"
#include "AliMUONTrack.h"
#include "AliMUONTrackParam.h"
#include "AliMUONVCluster.h"
#include "AliMUONVDigit.h"
#include "AliMUONQAMappingCheck.h"
#include "AliMUONVTrackerData.h"

#include "mapping/AliMpVSegmentation.h"
#include "mapping/AliMpSegmentation.h"
#include "mapping/AliMpPad.h"

#endif

/// \ingroup macros
/// \file MUONDisplayTrackRef.C
///
/// \author Ph. Pillot, Subatech, March. 2012
///
/// Macro to display trackRef in mchview format
///
/// To display the results: mchview --use MUONDisplayTrackRef.root

//-----------------------------------------------------------------------
void MUONDisplayTrackRef()
{
  
  // disable printout of AliMCEvent
  AliLog::SetClassDebugLevel("AliMCEvent",-1);
  
  // reset ROOT and connect tree file
  gROOT->Reset();
  
  // link to simulated tracks
  AliMCEventHandler mcEventHandler;
  mcEventHandler.SetInputPath("");
  mcEventHandler.InitIO("");
  AliMUONRecoCheck rc(0x0, &mcEventHandler);
  
  // ocdb access
  AliCDBManager::Instance()->SetDefaultStorage("raw://");
  AliCDBManager::Instance()->SetSpecificStorage("MUON/Align/Data","alien://folder=/alice/simulation/2008/v4-15-Release/Full");
//  AliCDBManager::Instance()->SetSpecificStorage("MUON/Align/Data","alien://folder=/alice/simulation/2008/v4-15-Release/Ideal");
  mcEventHandler.BeginEvent(0);
  Int_t run = mcEventHandler.MCEvent()->Header()->GetRun();
  mcEventHandler.FinishEvent();
  AliCDBManager::Instance()->SetRun(run);
  
  // load geometry
  AliGeomManager::LoadGeometry();
  if (!AliGeomManager::GetGeometry()) return;
  AliGeomManager::ApplyAlignObjsFromCDB("MUON");
  
  // load mapping
  if (!AliMUONCDB::LoadMapping()) return;
  
  // helper class to store cluster location
  AliMUONQAMappingCheck qaMappingStore(run);
  
  // loop over Events
  Int_t nEvents = mcEventHandler.GetTree()->GetEntries();
  for (Int_t iEvent = 0; iEvent < nEvents; iEvent++) {
    
    if ((iEvent+1)%1000 == 0) {
      printf("event %d\n", iEvent+1);
      cout<<flush;
    }
    
    if (!mcEventHandler.BeginEvent(iEvent)) continue;
    qaMappingStore.NewEvent();
    
    // get simulated tracks
    rc.ResetStores();
    AliMUONVTrackStore* trackRefStore = rc.TrackRefs(iEvent);
    
    // loop over trackRefs
    TIter next(trackRefStore->CreateIterator());
    AliMUONTrack* trackRef = 0x0;
    while ( ( trackRef = static_cast<AliMUONTrack*>(next()) ) ) {
      
      // loop over Hits
      for (Int_t iCl = 0; iCl < trackRef->GetNClusters(); iCl++) {
	
	AliMUONVCluster* cluster = static_cast<AliMUONTrackParam*>
	  (trackRef->GetTrackParamAtCluster()->UncheckedAt(iCl))->GetClusterPtr();
	
	// add digits and charge to hit so that it is accepted
	Int_t deId = cluster->GetDetElemId();
	const AliMpVSegmentation* seg0 = AliMpSegmentation::Instance()->GetMpSegmentation(deId,AliMp::GetCathodType(0));
	const AliMpVSegmentation* seg1 = AliMpSegmentation::Instance()->GetMpSegmentation(deId,AliMp::GetCathodType(1));
	AliMpPad pad0 = seg0->PadByIndices(seg0->MaxPadIndexX()/2,seg0->MaxPadIndexY()/2);
	AliMpPad pad1 = seg1->PadByIndices(seg1->MaxPadIndexX()/2,seg1->MaxPadIndexY()/2);
	cluster->AddDigitId(AliMUONVDigit::BuildUniqueID(deId,pad0.GetManuId(),pad0.GetManuChannel(),0));
	cluster->AddDigitId(AliMUONVDigit::BuildUniqueID(deId,pad1.GetManuId(),pad1.GetManuChannel(),1));
	cluster->SetCharge(10);
	
	// store cluster position
	qaMappingStore.Store(*cluster);
	
      }
     
    }
    
    // clean memory
    mcEventHandler.FinishEvent();
    
  }
  
  /// make a trackerData
  AliMUONVTrackerData *trackerData = qaMappingStore.CreateData("trackRef");
  
  // save results
  TFile *f = new TFile("MUONDisplayTrackRef.root", "RECREATE");
  trackerData->Write();
  f->Close();
  
  // draw results
  gSystem->Exec("mchview --use MUONDisplayTrackRef.root");
  
}

