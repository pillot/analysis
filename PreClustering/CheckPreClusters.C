//
//  CheckPreClusters.C
//  aliroot_dev
//
//  Created by philippe pillot on 24/01/2014.
//  Copyright (c) 2014 Philippe Pillot. All rights reserved.
//

#include <stdio.h>
#include <list>

#include <TFile.h>
#include <TTree.h>
#include <TExMap.h>

#include "AliCodeTimer.h"

#include "AliMUONVDigit.h"
#include "AliMUONVDigitStore.h"
#include "AliMUONVCluster.h"
#include "AliMUONVClusterStore.h"

const Int_t nDEs = 157; // 156 + 1 since 1st slote cannot be used (see comment below)
Int_t iDEmax = 0; // index must start from 1 because TExMap::GetValue(...) return 0 if key not found
Int_t deIds[nDEs];
TExMap deIndices;

void FindDigits(AliMUONVClusterStore *clusterStore, AliMUONVDigitStore *digitStore, std::list<UInt_t> *multiUsedDigitId, std::list<UInt_t> *missingDigitId);
void FindExtraDigits(AliMUONVDigitStore *digitStore, std::list<UInt_t> *extraDigitId);
void PrintDigits(std::list<UInt_t> *digitId, const Char_t* type);

//------------------------------------------------------------------
void CheckPreClusters(const char *clusterFileName, const char *digitFileName)
{
  /// Check that the preclusters contains all and only once the digits
  
  // read clusters
  TFile* clusterFile = new TFile(clusterFileName);
  if (!clusterFile) return;
  TTree* treeR = static_cast<TTree*>(clusterFile->Get("TreeR"));
  if (!treeR) return;
  AliMUONVClusterStore *clusterStore = AliMUONVClusterStore::Create(*treeR);
  clusterStore->Connect(*treeR);
  
  // read digits
  TFile* digitFile = TFile::Open(digitFileName);
  if (!digitFile) return;
  TTree* treeD = static_cast<TTree*>(digitFile->Get("TreeD"));
  if (!treeD) return;
  AliMUONVDigitStore* digitStore = AliMUONVDigitStore::Create(*treeD);
  digitStore->Connect(*treeD);
  
  std::list<UInt_t> multiUsedDigitId[nDEs];
  std::list<UInt_t> missingDigitId[nDEs];
  std::list<UInt_t> extraDigitId[nDEs];
  
  // loop over events
  Long64_t nEvents = treeR->GetEntries();
  if (treeD->GetEntries() != nEvents) {
    printf("Warning: not the same number of events in the cluster and the digit trees --> try with the smallest one\n");
    nEvents = TMath::Min(nEvents, treeD->GetEntries());
  }
  for (Long64_t iEv = 0; iEv < nEvents; ++iEv) {
    
    treeR->GetEntry(iEv);
    treeD->GetEntry(iEv);
    
    FindDigits(clusterStore, digitStore, multiUsedDigitId, missingDigitId);
    FindExtraDigits(digitStore, extraDigitId);
    
    Bool_t print = kFALSE;
    for (Int_t iDE = 1; iDE <= iDEmax; iDE++) {
      if (missingDigitId[iDE].size() > 0 || multiUsedDigitId[iDE].size() > 0 || extraDigitId[iDE].size() > 0) {
        print = kTRUE;
        break;
      }
    }
    if (print) {
      printf("Event %lld:\n", iEv);
      PrintDigits(missingDigitId, "Missing digits");
      PrintDigits(multiUsedDigitId, "Multi-used digits");
      PrintDigits(extraDigitId, "Extra digits");
    }
    
    clusterStore->Clear();
    digitStore->Clear();
    
  }
  
  AliCodeTimer::Instance()->Print();
  
}

//------------------------------------------------------------------
void FindDigits(AliMUONVClusterStore *clusterStore, AliMUONVDigitStore *digitStore, std::list<UInt_t> *multiUsedDigitId, std::list<UInt_t> *missingDigitId)
{
  /// find digits used in several preclusters or missing (?)
  
  AliCodeTimerAutoGeneral("",0);
  
  // empty the lists of digits
  for (Int_t iDE = 1; iDE < nDEs; iDE++) {
    multiUsedDigitId[iDE].clear();
    missingDigitId[iDE].clear();
  }
  
  // loop over clusters
  AliMUONVCluster* cluster = 0x0;
  TIter nextCluster(clusterStore->CreateIterator());
  while ((cluster = static_cast<AliMUONVCluster*>(nextCluster()))) {
    
    // get the DE index
    Int_t deId = cluster->GetDetElemId();
    Int_t iDE = deIndices.GetValue(deId);
    if (iDE == 0) {
      deIndices.Add(deId, ++iDEmax);
      deIds[iDEmax] = deId;
      iDE = iDEmax;
    }
    
    // loop over attached digits
    for (Int_t iDigit = 0; iDigit < cluster->GetNDigits(); iDigit++) {
      
      UInt_t digitId = cluster->GetDigitId(iDigit);
      AliMUONVDigit *digit = digitStore->FindObject(deId, AliMUONVDigit::ManuId(digitId), AliMUONVDigit::ManuChannel(digitId), AliMUONVDigit::Cathode(digitId));
      
      if (!digit || digit->Charge() <= 0) missingDigitId[iDE].push_back(digitId);
      else if (digit->IsUsed()) multiUsedDigitId[iDE].push_back(digitId);
      else digit->Used(kTRUE);
      
    }
    
  }
  
}

//------------------------------------------------------------------
void FindExtraDigits(AliMUONVDigitStore *digitStore, std::list<UInt_t> *extraDigitId)
{
  /// print digits never used in pre-clusters
  
  AliCodeTimerAutoGeneral("",0);
  
  // empty the lists of digits
  for (Int_t iDE = 1; iDE < nDEs; iDE++) extraDigitId[iDE].clear();
  
  // loop over digits
  AliMUONVDigit* digit = 0x0;
  TIter nextDigit(digitStore->CreateTrackerIterator());
  while ((digit = static_cast<AliMUONVDigit*>(nextDigit()))) {
    
    if (digit->Charge() <= 0 || digit->IsUsed()) continue;
    
    // get the DE index
    Int_t deId = digit->DetElemId();
    Int_t iDE = deIndices.GetValue(deId);
    if (iDE == 0) {
      deIndices.Add(deId, ++iDEmax);
      deIds[iDEmax] = deId;
      iDE = iDEmax;
    }
    
    extraDigitId[iDE].push_back(digit->GetUniqueID());
    
  }
  
}

//------------------------------------------------------------------
void PrintDigits(std::list<UInt_t> *digitId, const Char_t* type)
{
  /// print the digits in the list of given type
  
  AliCodeTimerAutoGeneral("",0);
  
  // loop over DEs
  for (Int_t iDE = 1; iDE <= iDEmax; iDE++) {
    
    Int_t nDigits = digitId[iDE].size();
    if (nDigits == 0) continue;
    
    printf("\t%s in DE %d = %d digitIDs = (", type, deIds[iDE], nDigits);
    
    for (std::list<UInt_t>::iterator digitIdIt = digitId[iDE].begin(); digitIdIt != digitId[iDE].end(); ++digitIdIt) {
      
      if (digitIdIt == digitId[iDE].begin()) printf("%u", *digitIdIt);
      else printf(", %u", *digitIdIt);
      
    }
    
    printf(")\n");
    
  }
  
}
