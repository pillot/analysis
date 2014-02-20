/*
 *  EfficiencyMapFromStatusMap.C
 *  aliroot_dev
 *
 *  Created by philippe pillot on 24/06/13.
 *  Copyright 2013 SUBATECH. All rights reserved.
 *
 */


#include <Riostream.h>
#include <list>

#include <TFile.h>
#include <TObjArray.h>
#include <TList.h>
#include <TSystem.h>
#include <TParameter.h>

#include "AliCDBManager.h"

#include "AliMUONCDB.h"
#include "AliMUONDigitCalibrator.h"
#include "AliMUONRejectList.h"
#include "AliMUONTrackerData.h"

#include "mapping/AliMpSegmentation.h"
#include "mapping/AliMpDDLStore.h"
#include "mapping/AliMpManuIterator.h"
#include "mapping/AliMpConstants.h"
#include "mapping/AliMpDetElement.h"
#include "mapping/AliMpVSegmentation.h"
#include "mapping/AliMpPad.h"


// Largest pixel size in such a way that every pad in non-bending and
// bending planes together can be divided in an integer number of pixel.
// The pixel size is divided by 2 in st1&2 because the cathod planes are
// shifted by half a pad the one relatively to the other in both directions.
Double_t pixelSize[5][2] = {{0.105, 0.21}, {0.375, 0.25}, {0.357143, 0.5}, {0.357143, 0.5}, {0.357143, 0.5}};


//---------------------------------------------------------------------------
void EfficiencyMapFromStatusMap(Int_t runNumber)
{
  /// Compute an efficiency map according to pad status for the given run (from raw OCDB)
  /// Pad OK     --> efficiency from  1 to 2 according to the status on the other cathod
  /// Pad not OK --> efficiency from -1 to 0 according to the status on the other cathod
  /// Manu efficiency from 0 to 1 according to the fraction of pixel OK on both cathods together
  /*
   aliroot -l
   .include $ALICE_ROOT/MUON
   .L $WORK/Macros/MuonEfficiency/EfficiencyMapFromStatusMap.C+
  */
  
  // load mapping if not already done
  AliCDBManager* man = AliCDBManager::Instance();
  if (!man->IsDefaultStorageSet()) man->SetDefaultStorage("raw://");
  if (man->GetRun() < 0) man->SetRun(runNumber);
  if (!AliMpSegmentation::Instance(kFALSE) || !AliMpDDLStore::Instance(kFALSE)) {
    if (!AliMUONCDB::LoadMapping()) return;
  }
  
  // the digit calibrator will determine the pad status
  AliMUONDigitCalibrator digitCalibrator(runNumber);
  
  AliMUONRejectList effMaps;
  std::list<double> dimx[5][2];
  std::list<double> dimy[5][2];
  std::list<int> npadx[5];
  std::list<int> npady[5];
  
  // loop over Manus
  Int_t iManu = 0;
  AliMpManuIterator manuIt;
  Int_t deId, manuId;
  while (manuIt.Next(deId, manuId)) {
    
    if ((++iManu)%100 == 0) cout << "\rprocessed Manus: " << iManu << flush;
    
    Int_t stationId = (deId/100-1)/2;
    Int_t nPixelsInManu = 0;
    Int_t nPixelsOkInManu = 0;
    
    // loop over channels
    for (Int_t iChannel = 0; iChannel < AliMpConstants::ManuNofChannels(); iChannel++) {
      
      // skip non-connected pads
      AliMpDetElement* de = AliMpDDLStore::Instance()->GetDetElement(deId);
      if (!de || !de->IsConnectedChannel(manuId,iChannel)) continue;
      
      // get the current pad
      const AliMpVSegmentation* seg = AliMpSegmentation::Instance()->GetMpSegmentationByElectronics(deId, manuId);  
      AliMpPad pad = seg->PadByLocation(manuId, iChannel, kFALSE);
      dimx[stationId][seg->PlaneType()].push_back(2.*pad.GetDimensionX());
      dimy[stationId][seg->PlaneType()].push_back(2.*pad.GetDimensionY());
      
      // get the number of pixels in this pad
      Int_t nPixelsx = (Int_t) (2. * pad.GetDimensionX() / pixelSize[stationId][0] + 0.5);
      Int_t nPixelsy = (Int_t) (2. * pad.GetDimensionY() / pixelSize[stationId][1] + 0.5);
      Int_t nPixels = nPixelsx * nPixelsy;
      Int_t nPixelsOk = 0;
      nPixelsInManu += 2 * nPixels;
      
      // set the efficiency of the current pad
      Float_t eff = -1.;
      if (digitCalibrator.IsValidDigit(deId, manuId, iChannel)) {
	nPixelsOkInManu += nPixels;
	eff = 1.;
      }
      
      // get the segmentation on the other cathod
      AliMp::CathodType otherCathod = AliMp::OtherCathodType(de->GetCathodType(seg->PlaneType()));
      const AliMpVSegmentation* seg2 = AliMpSegmentation::Instance()->GetMpSegmentation(deId, otherCathod);
      
      TList *padx = new TList[nPixelsy];
      for (Int_t iy = 0; iy < nPixelsy; iy++) padx[iy].SetOwner(kTRUE);
      TList *pady = new TList[nPixelsx];
      for (Int_t ix = 0; ix < nPixelsx; ix++) pady[ix].SetOwner(kTRUE);
      
      // look at overlapping pads on the other cathod
      for (Int_t ix = 0; ix < nPixelsx; ix++) {
	
	Double_t pixelPosx = pad.GetPositionX() - pad.GetDimensionX() + (ix+0.5) * pixelSize[stationId][0];
	
	for (Int_t iy = 0; iy < nPixelsy; iy++) {
	  
	  Double_t pixelPosy = pad.GetPositionY() - pad.GetDimensionY() + (iy+0.5) * pixelSize[stationId][1];
	  
	  // get the pad below this pixel
	  AliMpPad pad2 = seg2->PadByPosition(pixelPosx, pixelPosy, kFALSE);
	  if (!pad2.IsValid()) continue;
	  Int_t manuId2 = pad2.GetManuId();
	  Int_t channelId2 = pad2.GetManuChannel();
	  TString padId = Form("%u", ((manuId2 << 16) | channelId2));
	  
	  TObject* o = padx[iy].FindObject(padId.Data());
	  if (o) static_cast<TParameter<int>*>(o)->SetVal(static_cast<TParameter<int>*>(o)->GetVal()+1);
	  else padx[iy].AddLast(new TParameter<int>(padId.Data(), 1));
	  
	  o = pady[ix].FindObject(padId.Data());
	  if (o) static_cast<TParameter<int>*>(o)->SetVal(static_cast<TParameter<int>*>(o)->GetVal()+1);
	  else pady[ix].AddLast(new TParameter<int>(padId.Data(), 1));
	  
	  // check if this pad is OK
	  if (digitCalibrator.IsValidDigit(deId, manuId2, channelId2)) nPixelsOk++;
	  
	}
	
      }
      
      TParameter<int> *p = 0x0;
      for (Int_t iy = 0; iy < nPixelsy; iy++) {
	TIter nextPadx(&(padx[iy]));
	while ((p = static_cast<TParameter<int>*>(nextPadx()))) npadx[stationId].push_back(p->GetVal());
      }
      for (Int_t ix = 0; ix < nPixelsx; ix++) {
	TIter nextPady(&pady[ix]);
	while ((p = static_cast<TParameter<int>*>(nextPady()))) npady[stationId].push_back(p->GetVal());
      }
      delete[] padx;
      delete[] pady;
      
      // modify the pad efficiency according to the status of the other cathod
      eff += ((Float_t)nPixelsOk) / ((Float_t)nPixels);
      nPixelsOkInManu += nPixelsOk;
      
      // register the efficiency of the current pad
      effMaps.SetChannelProbability(deId, manuId, iChannel, eff);
      
    }
    
    // register the efficiency of the current manu
    effMaps.SetManuProbability(deId, manuId, ((Float_t)nPixelsOkInManu) / ((Float_t)nPixelsInManu));
    
  }
  
  cout << "\rprocessed Manus: " << iManu << endl;
  
  std::list<double>::iterator it;
  std::list<int>::iterator it2;
  TString plane[2] = {"bending", "non-bending"};
  for (Int_t iSt = 0; iSt < 5; iSt++) {
    std::cout << "\n\nStation " << iSt+1 <<"\n";
    for (Int_t iPlane = 0; iPlane < 2; iPlane++) {
      dimx[iSt][iPlane].sort();
      dimx[iSt][iPlane].unique();
      std::cout << "\n" << plane[iPlane].Data() << " plane: number of pixels in non-bending dimension of pads\n";
      for (it=dimx[iSt][iPlane].begin(); it!=dimx[iSt][iPlane].end(); ++it) std::cout << *it/pixelSize[iSt][0] << '\n';
    }
    npadx[iSt].sort();
    npadx[iSt].unique();
    std::cout << "\nnumber of times a pad is found in non-bending dimension of pads\n";
    for (it2=npadx[iSt].begin(); it2!=npadx[iSt].end(); ++it2) std::cout << *it2 << '\n';
    for (Int_t iPlane = 0; iPlane < 2; iPlane++) {
      dimy[iSt][iPlane].sort();
      dimy[iSt][iPlane].unique();
      std::cout << "\n" << plane[iPlane].Data() << " plane: number of pixels in bending dimension of pads\n";
      for (it=dimy[iSt][iPlane].begin(); it!=dimy[iSt][iPlane].end(); ++it) std::cout << *it/pixelSize[iSt][1] << '\n';
    }
    npady[iSt].sort();
    npady[iSt].unique();
    std::cout << "\nnumber of times a pad is found in bending dimension of pads\n";
    for (it2=npady[iSt].begin(); it2!=npady[iSt].end(); ++it2) std::cout << *it2 << '\n';
  }
  
  // put the efficiency maps in a tracker data for display
  AliMUONTrackerData effMapsDisp("efficiencies", "efficiency maps", effMaps);
  effMapsDisp.SetDimensionName(0,"efficiency");
  
  // save efficiency maps
  TList eff;
  eff.SetName("map");
  eff.SetOwner(kFALSE);
  eff.AddLast(&effMaps);
  TFile *outFile = TFile::Open("EfficiencyMapFromStatusMap.root", "RECREATE");
  if (outFile && outFile->IsOpen()) {
    eff.Write(0x0, TObject::kOverwrite);
    effMapsDisp.Write(0x0, TObject::kOverwrite);
    outFile->Close();
  }
  
  // display efficiency maps
  gSystem->Exec("mchview --use EfficiencyMapFromStatusMap.root");
  
}

