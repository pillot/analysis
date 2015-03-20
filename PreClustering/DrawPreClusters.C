//
//  DrawPreClusters.C
//  aliroot_dev
//
//  Created by philippe pillot on 17/02/2014.
//  Copyright (c) 2014 Philippe Pillot. All rights reserved.
//

#include <stdio.h>

#include <TFile.h>
#include <TTree.h>
#include <TExMap.h>
#include <TSystem.h>

#include "AliCodeTimer.h"
#include "AliCDBManager.h"

#include "AliMUONCDB.h"
#include "AliMUONVCluster.h"
#include "AliMUONVClusterStore.h"
#include "AliMUON2DMap.h"
#include "AliMUONVDigit.h"
#include "AliMUONCalibParamNI.h"
#include "AliMUONTrackerData.h"

#include "AliMpConstants.h"


const Int_t nDEs = 157; // 156 + 1 since 1st slote cannot be used (see comment below)
Int_t iDEmax = 0; // index must start from 1 because TExMap::GetValue(...) return 0 if key not found
TExMap deIndices;
Int_t iCluster[nDEs];


//------------------------------------------------------------------
void DrawPreClusters(const char *clusterFileName, Int_t iEvent, const char *outFileName)
{
  /// Draw the preclusters in given event
  /*
   .x $ALICE_ROOT/MUON/rootlogon.C
   .x $WORK/Macros/PreClustering/DrawPreClusters.C+(...)
   */
  
  AliCodeTimerAutoGeneral("",0);
  
  // load mapping
  AliCDBManager* man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  man->SetRun(0);
  if (!AliMUONCDB::LoadMapping()) return;
  
  // read clusters
  TFile* clusterFile = new TFile(clusterFileName);
  if (!clusterFile) return;
  TTree* treeR = static_cast<TTree*>(clusterFile->Get("TreeR"));
  if (!treeR) return;
  AliMUONVClusterStore *clusterStore = AliMUONVClusterStore::Create(*treeR);
  clusterStore->Connect(*treeR);
  
  Long64_t nEvents = treeR->GetEntries();
  if (iEvent < 0 || iEvent >= nEvents) {
    printf("choose an event in the range 0-%lld\n", nEvents-1);
    return;
  }
  
  treeR->GetEntry(iEvent);
  
  AliMUON2DMap digitStore(kTRUE);
  
  // loop over clusters
  AliMUONVCluster* cluster = 0x0;
  TIter nextCluster(clusterStore->CreateIterator());
  while ((cluster = static_cast<AliMUONVCluster*>(nextCluster()))) {
    
    // get the DE index
    Int_t deId = cluster->GetDetElemId();
    Int_t iDE = deIndices.GetValue(deId);
    if (iDE == 0) {
      deIndices.Add(deId, ++iDEmax);
      iCluster[iDEmax] = 1;
      iDE = iDEmax;
    }
    
    // loop over attached digits
    for (Int_t iDigit = 0; iDigit < cluster->GetNDigits(); iDigit++) {
      
      UInt_t digitId = cluster->GetDigitId(iDigit);
      Int_t manuId = AliMUONVDigit::ManuId(digitId);
      Int_t manuChannel = AliMUONVDigit::ManuChannel(digitId);
      
      // register the digit
      AliMUONVCalibParam* c = static_cast<AliMUONVCalibParam*>(digitStore.FindObject(deId, manuId));
      if (!c) {
        c = new AliMUONCalibParamNI(1, AliMpConstants::ManuNofChannels(), deId, manuId);
        digitStore.Add(c);
      }
      c->SetValueAsInt(manuChannel, 0, iCluster[iDE]);
      
    }
    
    ++iCluster[iDE];
    
  }
  
  clusterStore->Clear();
  
  // create the tracker data
  AliMUONTrackerData digitData("preclusters", "preclusters", 1, kTRUE);
  digitData.SetDimensionName(0, "index");
  digitData.Add(digitStore);
  
  // save it to a file
  TFile *outFile = TFile::Open(outFileName, "UPDATE");
  if (outFile && outFile->IsOpen()) {
    digitData.Write(0x0, TObject::kOverwrite);
    outFile->Close();
  }
  
  AliCodeTimer::Instance()->Print();
  
  // draw it
  //gSystem->Exec(Form("mchview --use %s", outFileName));
  
}

