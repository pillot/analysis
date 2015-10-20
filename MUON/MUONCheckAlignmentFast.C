#if !defined(__CINT__) || defined(__MAKECINT__)
// ROOT includes
#include <Riostream.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TGeoManager.h>
#include <TSystem.h>

// STEER includes
#include "AliCDBManager.h"
#include "AliGeomManager.h"
#include "AliESDEvent.h"
#include "AliESDMuonTrack.h"
#include "AliESDMuonCluster.h"

// MUON includes
#include "AliMUONCDB.h"
#include "AliMUON2DMap.h"
#include "AliMUONConstants.h"
#include "AliMUONCalibParamND.h"
#include "AliMUONGeometryTransformer.h"
#include "AliMUONTrackerData.h"

#include "AliMpConstants.h"
#include "AliMpDEIterator.h"
#include "AliMpDetElement.h"
#include "AliMpVSegmentation.h"
#include "AliMpSegmentation.h"
#include "AliMpVPadIterator.h"
#include "AliMpPad.h"

#endif


/// \ingroup macros
/// \file MUONCheckAlignment.C
///
/// \author Ph. Pillot, Subatech, March. 2012
///
/// Macro to compare 2 alignment files


//Int_t run = 167818;
//Int_t run = 169099; // LHC11h
//Int_t run = 192732; // LHC12h
//Int_t run = 197388; // LHC13f
Int_t run = 221128; // LHC15c

TString alignStorage[2] = {
//  "alien://folder=/alice/cern.ch/user/p/ppillot/OCDB2012",
  "alien://folder=/alice/data/2015/OCDB",
//  "alien://folder=/alice/simulation/2008/v4-15-Release/Ideal",
//  "alien://folder=/alice/cern.ch/user/p/ppillot/OCDB_PbPbSim",
//  "alien://folder=/alice/cern.ch/user/j/jcastill/LHC11hMisAlignCDB3"
//  "alien://folder=/alice/cern.ch/user/j/jcastill/pbpb11wrk/LHC11hMisAlignCDB4",
//  "alien://folder=/alice/data/2011/OCDB"
//  "alien://folder=/alice/cern.ch/user/j/jcastill/ReAligni00pbpb11CDB2",
//  "alien://folder=/alice/cern.ch/user/j/jcastill/LHC11hMisAlignCDB3"
//  "alien://folder=/alice/simulation/2008/v4-15-Release/Residual"
//  "alien://folder=/alice/simulation/2008/v4-15-Release/Full"
  "local://."
};

// disable the display at the channel level to get an output file
// of a reasonable size that fit in a mail
Bool_t disableChannelLevel = kFALSE;

// show the areas covered with the first alignment but not with the second
Bool_t showDeadAreas = kFALSE;


//-----------------------------------------------------------------------
void MUONCheckAlignmentFast()
{
  
  // reset ROOT and connect tree file
  gROOT->Reset();
  
  // ocdb access
  AliCDBManager* cdbm = AliCDBManager::Instance();
  cdbm->SetDefaultStorage("alien://folder=/alice/data/2015/OCDB");
  cdbm->SetRun(run);
  
  // load mapping
  if (!AliMUONCDB::LoadMapping()) return;
  
  // get geometry transformers
  AliMUONGeometryTransformer geoTransformer[2];
  for (Int_t i = 0; i < 2; i++) {
    cdbm->UnloadFromCache("GRP/Geometry/Data");
    cdbm->UnloadFromCache("MUON/Align/Data");
    AliGeomManager::GetGeometry()->UnlockGeometry();
    AliGeomManager::LoadGeometry();
    if (!AliGeomManager::GetGeometry()) return;
    cdbm->SetSpecificStorage("MUON/Align/Data",alignStorage[i].Data());
    AliGeomManager::ApplyAlignObjsFromCDB("MUON");
    geoTransformer[i].LoadGeometryData();
  }
  
  // store for cluster shifts
  AliMUON2DMap shiftStore(kTRUE);
  
  // loop over chamber
  for (Int_t iCh = 0; iCh < AliMUONConstants::NTrackingCh(); iCh++) {
    
    // loop over DEs
    AliMpDEIterator nextDE;
    nextDE.First(iCh);
    while (!nextDE.IsDone()) {
      
      Int_t deId = nextDE.CurrentDE()->GetId();
      
      // loop over cathods
      for (Int_t icath = 0; icath < 2; icath++) {
	
	const AliMpVSegmentation* seg = AliMpSegmentation::Instance()->GetMpSegmentation(deId,AliMp::GetCathodType(icath));
	
	// loop over pads
	AliMpVPadIterator *nextPad = seg->CreateIterator();
	nextPad->First();
	while (!nextPad->IsDone()) {
	  
	  AliMpPad pad = nextPad->CurrentItem();
	  Int_t manuId = pad.GetManuId();
	  Int_t manuChannel = pad.GetManuChannel();
	  
	  // local position
	  Double_t xl = pad.GetPositionX();
	  Double_t yl = pad.GetPositionY();
	  Double_t zl = 0.;
	  
	  // position with first alignment
	  Double_t x1, y1, z1;
	  geoTransformer[0].Local2Global(deId,xl,yl,zl,x1,y1,z1);
	  
	  // position with second alignment
	  Double_t x2, y2, z2;
	  geoTransformer[1].Local2Global(deId,xl,yl,zl,x2,y2,z2);
	  
	  // pad shift
	  Double_t dx = x2 - x1;
	  Double_t dy = y2 - y1;
	  Double_t dz = z2 - z1;
	  
	  // move the pad according to the difference between the 2 geometries
	  if (showDeadAreas) {
	    Double_t x0, y0, z0;
	    geoTransformer[0].Local2Global(deId,xl,yl,zl,x0,y0,z0);
	    
	    Double_t xl2, yl2, zl2;
	    geoTransformer[1].Global2Local(deId,x0,y0,z0,xl2,yl2,zl2);
	    
	    AliMpPad pad2 = seg->PadByPosition(xl2,yl2);
	    if (!pad2.IsValid()) {
	      nextPad->Next();
	      continue;
	    }
	    dx = dy = dz = -1;
	  }
	  
	  // store pad shifts
	  AliMUONVCalibParam* p = static_cast<AliMUONVCalibParam*>(shiftStore.FindObject(deId,manuId));
	  if (!p) {
	    p = new AliMUONCalibParamND(3,AliMpConstants::ManuNofChannels(),deId,manuId,0.);
	    shiftStore.Add(p);
	  }
	  p->SetValueAsDouble(manuChannel,0,dx);
	  p->SetValueAsDouble(manuChannel,1,dy);
	  p->SetValueAsDouble(manuChannel,2,dz);
	  
	  nextPad->Next();
	}
	
	delete nextPad;
	
      }
      
      nextDE.Next();
    }
    
  }
  
  // create tracker data
  AliMUONTrackerData data("shifts","shifts",3,kTRUE);
  data.SetDimensionName(0,"dx"); // max shift in x
  data.SetDimensionName(1,"dy"); // max shift in y
  data.SetDimensionName(2,"dz"); // max shift in z
  if (disableChannelLevel) data.DisableChannelLevel();
  data.Add(shiftStore);
  
  // save results
  TFile *f = TFile::Open("MUONCheckAlignmentFast.root", "RECREATE");
  data.Write();
  f->Close();
  
  // draw results
  gSystem->Exec("mchview --use MUONCheckAlignmentFast.root");
  
}

