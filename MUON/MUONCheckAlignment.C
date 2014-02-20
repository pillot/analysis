#if !defined(__CINT__) || defined(__MAKECINT__)
// ROOT includes
#include <Riostream.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TGeoManager.h>
#include <TSystem.h>
#include <TH2F.h>

// STEER includes
#include "AliCDBManager.h"
#include "AliGeomManager.h"
#include "AliESDEvent.h"
#include "AliESDMuonTrack.h"
#include "AliESDMuonCluster.h"

// MUON includes
#include "AliMUONCDB.h"
#include "AliMUON2DMap.h"
#include "AliMUONCalibParamND.h"
#include "AliMUONGeometryTransformer.h"
#include "AliMUONTrackerData.h"

#include "mapping/AliMpConstants.h"
#include "mapping/AliMpDDLStore.h"
#include "mapping/AliMpDetElement.h"
#include "mapping/AliMpVSegmentation.h"
#include "mapping/AliMpSegmentation.h"
#include "mapping/AliMpPad.h"


#endif

/// \ingroup macros
/// \file MUONCheckAlignment.C
///
/// \author Ph. Pillot, Subatech, March. 2012
///
/// Macro to compare 2 alignment files:
/// - Align[0] is the alignment used for reconstruction
/// - Align[1] and Align[2] are the alignments to be compared


TString align[3] = {
  "alien://folder=/alice/data/2011/OCDB",
//  "alien://folder=/alice/simulation/2008/v4-15-Release/Residual",
  "alien://folder=/alice/data/2011/OCDB",
  "alien://folder=/alice/cern.ch/user/j/jcastill/ReAligni00pbpb11CDB2"
//  "alien://folder=/alice/cern.ch/user/j/jcastill/LHC11hMisAlignCDB"
};


//-----------------------------------------------------------------------
void MUONCheckAlignment()
{
  
  // reset ROOT and connect tree file
  gROOT->Reset();
  
  // link to ESD MUON tracks
  TFile *esdFile = TFile::Open("AliESDs.root"); // open the file
  if (!esdFile || !esdFile->IsOpen()) return;
  TTree *esdTree = (TTree*) esdFile->Get("esdTree"); // get the tree
  if (!esdTree) return;
  AliESDEvent esdEvent;
  esdEvent.ReadFromTree(esdTree); // link ESD event to the tree
  
  // ocdb access
  AliCDBManager* cdbm = AliCDBManager::Instance();
  cdbm->SetDefaultStorage("alien://folder=/alice/data/2011/OCDB");
  esdTree->GetEvent(0);
  cdbm->SetRun(esdEvent.GetRunNumber());
  
  // load mapping
  if (!AliMUONCDB::LoadMapping()) return;
  
  // get geometry transformers
  AliMUONGeometryTransformer geoTransformer[3];
  for (Int_t i = 0; i < 3; i++) {
    cdbm->UnloadFromCache("GRP/Geometry/Data");
    cdbm->UnloadFromCache("MUON/Align/Data");
    AliGeomManager::GetGeometry()->UnlockGeometry();
    AliGeomManager::LoadGeometry();
    if (!AliGeomManager::GetGeometry()) return;
    cdbm->SetSpecificStorage("MUON/Align/Data",align[i].Data());
    AliGeomManager::ApplyAlignObjsFromCDB("MUON");
    geoTransformer[i].LoadGeometryData();
  }
  
  // store for cluster shifts
  AliMUON2DMap shiftStore[3] = {AliMUON2DMap(kTRUE), AliMUON2DMap(kTRUE), AliMUON2DMap(kTRUE)};
  
  TH2F *h = new TH2F("h","h",300,0.,300.,100.,0.,100.);
  h->SetDirectory(0);
  
  // loop over Events
  Int_t nEvents = esdTree->GetEntries();
  for (Int_t iEvent = 0; iEvent < nEvents; iEvent++) {
    
    if ((iEvent+1)%1000 == 0) {
      printf("event %d\n", iEvent+1);
      cout<<flush;
    }
    
    if (esdTree->GetEvent(iEvent) <= 0) continue;
    
    // loop over tracks
    Int_t nTracks = esdEvent.GetNumberOfMuonTracks();
    for (Int_t iTr = 0; iTr < nTracks; iTr++) {
      
      AliESDMuonTrack *esdTrack = esdEvent.GetMuonTrack(iTr);
      if (!esdTrack->ContainTrackerData()) continue;
      
      // loop over clusters
      for (Int_t iCl = 0; iCl < esdTrack->GetNClusters(); iCl++) {
	
	AliESDMuonCluster* cluster = static_cast<AliESDMuonCluster*>(esdTrack->GetClusters().UncheckedAt(iCl));
	Int_t deId = cluster->GetDetElemId();
	
	// position from reconstruction
	Double_t x0 = cluster->GetX();
	Double_t y0 = cluster->GetY();
	Double_t z0 = cluster->GetZ();
	
	// local position
	Double_t xl, yl, zl;
	geoTransformer[0].Global2Local(deId,x0,y0,z0,xl,yl,zl);
	
	// position with first alignment
	Double_t x1, y1, z1;
	geoTransformer[1].Local2Global(deId,xl,yl,zl,x1,y1,z1);
	
	// position with second alignment
	Double_t x2, y2, z2;
	geoTransformer[2].Local2Global(deId,xl,yl,zl,x2,y2,z2);
	
	// cluster shifts
	Double_t shift[3] = {x2-x1, y2-y1, z2-z1};
	
	//printf("DE%d: %f, %f, %f\n",deId, dx, dy, dz);
	if (deId == 904) h->Fill(xl, yl);
	
	// loop over cathods
	for (Int_t icath = 0; icath < 2; icath++) {
	  
	  // cluster location
	  const AliMpVSegmentation* seg = AliMpSegmentation::Instance()->GetMpSegmentation(deId,AliMp::GetCathodType(icath));
	  AliMpPad pad = seg->PadByPosition(xl,yl);
	  if (!pad.IsValid()) continue;
	  Int_t manuId = pad.GetManuId();
	  
	  // store cluster shifts
	  for (Int_t ishift = 0; ishift < 3; ishift++) {
	    
	    AliMUONVCalibParam* p = static_cast<AliMUONVCalibParam*>(shiftStore[ishift].FindObject(deId,manuId));
	    if (!p) {
	      p = new AliMUONCalibParamND(5,AliMpConstants::ManuNofChannels(),deId,manuId,0.);
	      p->SetValueAsDouble(0,3,AliMpDDLStore::Instance()->GetDetElement(deId)->NofChannelsInManu(manuId));
	      p->SetValueAsDouble(0,4,nEvents);
	      shiftStore[ishift].Add(p);
	    }
	    p->SetValueAsDouble(0,0,p->ValueAsDouble(0,0)+shift[ishift]);
	    p->SetValueAsDouble(0,1,p->ValueAsDouble(0,1)+shift[ishift]*shift[ishift]);
	    p->SetValueAsDouble(0,2,p->ValueAsDouble(0,2)+1.);
	    
	  }
	  
	}
	
      }
      
    }
    
  }
  
  // create tracker data
  AliMUONTrackerData data[3] = {
    AliMUONTrackerData("dx","dx",shiftStore[0]),
    AliMUONTrackerData("dy","dy",shiftStore[1]),
    AliMUONTrackerData("dz","dz",shiftStore[2])
  };
  
  // save results
  TFile *f = TFile::Open("MUONCheckAlignment.root", "RECREATE");
  for (Int_t ishift = 0; ishift < 3; ishift++) {
    data[ishift].SetInternalDimensionName(0,"mean");
    data[ishift].SetInternalDimensionName(1,"sigma");
    data[ishift].Write();
  }
  h->Write();
  f->Close();
  
  // draw results
  gSystem->Exec("mchview --use MUONCheckAlignment.root");
  
}

