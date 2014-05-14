//
//  DrawOccupancy.C
//  aliroot_dev
//
//  Created by philippe pillot on 30/04/2014.
//  Copyright (c) 2014 Philippe Pillot. All rights reserved.
//

#include <stdio.h>

#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>

#include "AliCDBManager.h"

#include "AliMUONCDB.h"
#include "AliMUONCalibrationData.h"
#include "AliMUONRecoParam.h"
#include "AliMUONDigitCalibrator.h"
#include "AliMUONConstants.h"
#include "AliMUONVDigit.h"
#include "AliMUONVDigitStore.h"
#include "AliMUON2DMap.h"
#include "AliMUONVCalibParam.h"
#include "AliMUONCalibParamND.h"
#include "AliMUONTrackerData.h"

#include "AliMpConstants.h"
#include "AliMpCDB.h"
#include "AliMpDEIterator.h"
#include "AliMpVSegmentation.h"
#include "AliMpSegmentation.h"
#include "AliMpVPadIterator.h"
#include "AliMpPad.h"


void GetNWorkingPads(Int_t run, Int_t nPadsTot[10], AliMUON2DMap &occStore);
void CountFiredPads(AliMUONVDigitStore *digitStore, Int_t nPadsCurrent[10], AliMUON2DMap &occStore, Long64_t nEvents);


//------------------------------------------------------------------
void DrawOccupancy(const char *digitFileName, Int_t run)
{
  /// compute detector occupancy from reconstructed digits
  /*
   .x $ALICE_ROOT/MUON/rootlogon.C
   .x $WORK/Macros/PreClustering/DrawOccupancy.C+
   */
  
  AliMUON2DMap occStore(kTRUE);
  
  // get the number of working pads
  Int_t nPadsTot[10];
  GetNWorkingPads(run, nPadsTot, occStore);
  for (Int_t iCh = 0; iCh < 10; ++iCh)
    printf("ch%02d: nPadsTot = %d\n", iCh+1, nPadsTot[iCh]);
  
  // read digits
  TFile* digitFile = TFile::Open(digitFileName);
  if (!digitFile) return;
  TTree* treeD = static_cast<TTree*>(digitFile->Get("TreeD"));
  if (!treeD) return;
  AliMUONVDigitStore* digitStore = AliMUONVDigitStore::Create(*treeD);
  digitStore->Connect(*treeD);
  
  Int_t nPads[10];
  TH1F *hOcc[10];
  for (Int_t iCh = 0; iCh < 10; ++iCh) {
    hOcc[iCh] = new TH1F(Form("hOccCh%d",iCh+1), Form("occupancy of chamber %d;occupancy;#events",iCh+1), 100, 0., 10.);
    hOcc[iCh]->SetDirectory(0);
    nPads[iCh] = 0;
  }
  
  // loop over events
  Int_t nPadsCurrent[10];
  Long64_t nEvents = treeD->GetEntries();
  for (Long64_t iEv = 0; iEv < nEvents; ++iEv) {
    
    treeD->GetEntry(iEv);
    
    CountFiredPads(digitStore, nPadsCurrent, occStore, nEvents);
    for (Int_t iCh = 0; iCh < 10; ++iCh) {
      nPads[iCh] += nPadsCurrent[iCh];
      hOcc[iCh]->Fill(100.*nPadsCurrent[iCh]/nPadsTot[iCh]);
    }
    
    digitStore->Clear();
    
  }
  
  // print results
  for (Int_t iCh = 0; iCh < 10; ++iCh)
    printf("ch%02d: <nPads> = %07.2f --> occ = %04.2f%%\n", iCh+1, ((Float_t)nPads[iCh])/((Float_t)nEvents),
           100.*((Float_t)nPads[iCh])/((Float_t)nEvents)/((Float_t)nPadsTot[iCh]));
  
  AliMUONTrackerData occupancy("occupancy", "occupancy", 1, kFALSE);
  occupancy.SetDimensionName(0, "occ");
  occupancy.Add(occStore);
  
  // display results
  TCanvas *cOcc = new TCanvas("cOcc", "occupancy per chamber", 1200, 500);
  cOcc->Divide(5,2);
  for (Int_t iCh = 0; iCh < 10; ++iCh) {
    cOcc->cd(iCh+1);
    gPad->SetLogy();
    hOcc[iCh]->Draw();
  }
  
  // save occupancy to a file
  TFile *outFile = TFile::Open("occupancy.root", "UPDATE");
  if (outFile && outFile->IsOpen()) {
    occupancy.Write(0x0, TObject::kOverwrite);
    cOcc->Write(0x0, TObject::kOverwrite);
    outFile->Close();
  }
  
  // display occupancy
  gSystem->Exec("mchview --use occupancy.root");
  
}


//------------------------------------------------------------------
void GetNWorkingPads(Int_t run, Int_t nPadsTot[10], AliMUON2DMap &occStore)
{
  /// get the numbe of working pads
  
  // load OCDB
  AliCDBManager* man = AliCDBManager::Instance();
  man->SetDefaultStorage("raw://");
  man->SetRun(run);
  if (!AliMUONCDB::LoadMapping()) return;
  AliMUONRecoParam *recoParam = AliMUONCDB::LoadRecoParam();
  
  // prepare digit calibrator
  AliMUONCalibrationData *calibrationData = new AliMUONCalibrationData(run);
  AliMUONDigitCalibrator *digitCalibrator = new AliMUONDigitCalibrator(*calibrationData,recoParam);
  
  // loop over DEs
  for (Int_t iCh = 0; iCh < AliMUONConstants::NTrackingCh(); iCh++) {
    
    nPadsTot[iCh] = 0;
    
    AliMpDEIterator deIt;
    deIt.First(iCh);
    while (!deIt.IsDone()) {
      
      Int_t deId = deIt.CurrentDEId();
      
      const AliMpVSegmentation* seg[2] =
      { AliMpSegmentation::Instance()->GetMpSegmentation(deId,AliMp::kCath0),
        AliMpSegmentation::Instance()->GetMpSegmentation(deId,AliMp::kCath1)
      };
      
      // loop over cathods
      for (Int_t iCath = 0; iCath < 2; iCath++) {
        
        // loop over pads
        Int_t iPad = 0;
        AliMpVPadIterator* padIt = seg[iCath]->CreateIterator();
        padIt->First();
        while (!padIt->IsDone()) {
          
          AliMpPad pad = padIt->CurrentItem();
          Int_t manuId = pad.GetManuId();
          Int_t manuChannel = pad.GetManuChannel();
          
          if (digitCalibrator->IsValidDigit(deId, manuId, manuChannel)) {
            
            // make sure all valid pads are taken into account
            AliMUONVCalibParam *c = static_cast<AliMUONVCalibParam*>(occStore.FindObject(deId, manuId));
            if (!c) {
              c = new AliMUONCalibParamND(1, AliMpConstants::ManuNofChannels(), deId, manuId, AliMUONVCalibParam::InvalidFloatValue());
              occStore.Add(c);
            }
            c->SetValueAsDouble(manuChannel, 0, 0.);
            
            // count valid digits
            ++nPadsTot[iCh];
            
          }
          
          padIt->Next();
          
        }
        
      }
      
      deIt.Next();
      
    }
    
  }
  
}


//------------------------------------------------------------------
void CountFiredPads(AliMUONVDigitStore *digitStore, Int_t nPadsCurrent[10], AliMUON2DMap &occStore, Long64_t nEvents)
{
  /// print digits occupancy per chamber
  
  for (Int_t iCh = 0; iCh < 10; ++iCh) nPadsCurrent[iCh] = 0;
  
  // loop over digits
  AliMUONVDigit* digit = 0x0;
  TIter nextDigit(digitStore->CreateIterator());
  while ((digit = static_cast<AliMUONVDigit*>(nextDigit()))) {
    
    if (digit->Charge() <= 0) continue;
    
    Int_t deId = digit->DetElemId();
    Int_t chId = deId/100-1;
    if (chId < 0 || chId > 9) continue;
    Int_t manuId = digit->ManuId();
    Int_t manuChannel = digit->ManuChannel();
    
    // count fired digits
    AliMUONVCalibParam *c = static_cast<AliMUONVCalibParam*>(occStore.FindObject(deId, manuId));
    if (!c) {
      c = new AliMUONCalibParamND(1, AliMpConstants::ManuNofChannels(), deId, manuId, AliMUONVCalibParam::InvalidFloatValue());
      occStore.Add(c);
    }
    Double_t occ = c->ValueAsDouble(manuChannel, 0);
    if (occ >= AliMUONVCalibParam::InvalidFloatValue()) {
      printf("invalid pad fired!\n");
      occ = 0.;
    }
    c->SetValueAsDouble(manuChannel, 0, occ + 1./nEvents);
    
    ++nPadsCurrent[chId];
    
  }
  
}

