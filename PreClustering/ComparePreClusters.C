//
//  ComparePreClusters.C
//  aliroot_dev
//
//  Created by philippe pillot on 22/01/2014.
//  Copyright (c) 2014 Philippe Pillot. All rights reserved.
//

#include <stdio.h>
#include <list>

#include <TFile.h>
#include <TTree.h>
#include <TExMap.h>

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

struct preCluster {
  AliMUONVCluster *cluster; // link to the corresponding cluster
  std::list<UInt_t> digitId; // list of digit Ids
  Bool_t compareMe; // kFALSE if already associate to an identical precluster
};

const Int_t nDEs = 157; // 156 + 1 since 1st slote cannot be used (see comment below)
Int_t iDEmax = 0; // index must start from 1 because TExMap::GetValue(...) return 0 if key not found
Int_t deIds[nDEs];
TExMap deIndices;

Int_t LoadPreClusters(AliMUONVClusterStore *clusterStore, std::list<preCluster> *preclusters);
Bool_t GetDifferences(std::list<preCluster> *preclusters1, std::list<preCluster> *preclusters2,
                      std::list<preCluster> *preclustersDiff1, std::list<preCluster> *preclustersDiff2);
Bool_t AreIdentical(preCluster &cl1, preCluster &cl2);
void PrintDifferences(std::list<preCluster> *preclustersDiff1, std::list<preCluster> *preclustersDiff2,
                      Int_t &integratedNDiff1, Int_t &integratedNDiff2);
void DrawDifferences(std::list<preCluster> *preclustersDiff1, std::list<preCluster> *preclustersDiff2, const char *outFileName);

//------------------------------------------------------------------
void ComparePreClusters(const char *clusterFileName1, const char *clusterFileName2, Long64_t iEvent = -2)
{
  /// Compare the preclusters between the two input files
  /// iEvent <= -2: do not draw
  /// iEvent == -1: draw for all concerned events
  /// iEvent >=  0: draw only for this event (if difference found)
  
  // load mapping
  AliCDBManager* man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALIROOT_OCDB_ROOT/OCDB");
  man->SetRun(0);
  if (!AliMUONCDB::LoadMapping()) return;
  
  // input clusters 1
  TFile* clusterFile1 = new TFile(clusterFileName1);
  if (!clusterFile1) return;
  TTree* treeR1 = static_cast<TTree*>(clusterFile1->Get("TreeR"));
  if (!treeR1) return;
  AliMUONVClusterStore *clusterStore1 = AliMUONVClusterStore::Create(*treeR1);
  clusterStore1->Connect(*treeR1);
  
  // input clusters 2
  TFile* clusterFile2 = new TFile(clusterFileName2);
  if (!clusterFile2) return;
  TTree* treeR2 = static_cast<TTree*>(clusterFile2->Get("TreeR"));
  if (!treeR2) return;
  AliMUONVClusterStore *clusterStore2 = AliMUONVClusterStore::Create(*treeR2);
  clusterStore2->Connect(*treeR2);
  
  // out directory for displays
  if (iEvent >= -1 && gSystem->AccessPathName("displays")) gSystem->Exec(Form("mkdir displays"));
  
  std::list<preCluster> preclusters1[nDEs];
  std::list<preCluster> preclusters2[nDEs];
  std::list<preCluster> preclustersDiff1[nDEs];
  std::list<preCluster> preclustersDiff2[nDEs];
  Int_t integratedNDiff1 = 0;
  Int_t integratedNDiff2 = 0;
  
  // loop over events
  Long64_t nEvents = treeR1->GetEntries();
  if (treeR2->GetEntries() != nEvents) {
    printf("Warning: not the same number of events in the two cluster trees --> try with the smallest one\n");
    nEvents = TMath::Min(nEvents, treeR2->GetEntries());
  }
  for (Long64_t iEv = 0; iEv < nEvents; ++iEv) {
    
    printf("Event %lld:\n", iEv);
    
    treeR1->GetEntry(iEv);
    treeR2->GetEntry(iEv);
    
    Int_t nPreclusters1 = LoadPreClusters(clusterStore1, preclusters1);
    Int_t nPreclusters2 = LoadPreClusters(clusterStore2, preclusters2);
    printf("\tTotal number of preclusters = %d / %d\n", nPreclusters1, nPreclusters2);
    
    Bool_t differenceFound = GetDifferences(preclusters1, preclusters2, preclustersDiff1, preclustersDiff2);
    
    PrintDifferences(preclustersDiff1, preclustersDiff2, integratedNDiff1, integratedNDiff2);
    
    if (differenceFound && (iEvent == -1 || iEv == iEvent))
      DrawDifferences(preclustersDiff1, preclustersDiff2, Form("displays/%lld.root", iEv));
    
    clusterStore1->Clear();
    clusterStore2->Clear();
    
  }
  
  printf("\nIntegrated number of different preclusters = %d / %d\n\n", integratedNDiff1, integratedNDiff2);
  
  AliCodeTimer::Instance()->Print();
  
}

//------------------------------------------------------------------
Int_t LoadPreClusters(AliMUONVClusterStore *clusterStore, std::list<preCluster> *preclusters)
{
  /// fill the mpDE structure with reconstructed digits
  
  AliCodeTimerAutoGeneral("",0);
  
  Int_t nPreclusters = 0;
  
  // empty the lists of preclusters
  for (Int_t iDE = 1; iDE < nDEs; iDE++) preclusters[iDE].clear();
  
  // loop over clusters
  AliMUONVCluster* cluster = 0x0;
  TIter nextCluster(clusterStore->CreateIterator());
  while ((cluster = static_cast<AliMUONVCluster*>(nextCluster()))) {
    
    preCluster cl;
    
    cl.cluster = cluster;
    cl.compareMe = kTRUE;
    
    for (Int_t iDigit = 0; iDigit < cluster->GetNDigits(); iDigit++)
      cl.digitId.push_back(cluster->GetDigitId(iDigit));
    cl.digitId.sort();
    
    Int_t deId = cluster->GetDetElemId();
    Int_t iDE = deIndices.GetValue(deId);
    if (iDE == 0) {
      deIndices.Add(deId, ++iDEmax);
      deIds[iDEmax] = deId;
      iDE = iDEmax;
    }
    
    preclusters[iDE].push_back(cl);
    nPreclusters++;
    
  }
  
  return nPreclusters;
  
}

//------------------------------------------------------------------
Bool_t GetDifferences(std::list<preCluster> *preclusters1, std::list<preCluster> *preclusters2,
                      std::list<preCluster> *preclustersDiff1, std::list<preCluster> *preclustersDiff2)
{
  /// For every DE: find in each list the preclusters which are different from the ones in the other list
  
  AliCodeTimerAutoGeneral("",0);
  
  Bool_t diffenceFound = kFALSE;
  
  // loop over DEs
  for (Int_t iDE = 1; iDE <= iDEmax; iDE++) {
    
    preclustersDiff1[iDE].clear();
    preclustersDiff2[iDE].clear();
    
    // loop over preclusters in the first list
    for (std::list<preCluster>::iterator clusterIt1 = preclusters1[iDE].begin(); clusterIt1 != preclusters1[iDE].end(); ++clusterIt1) {
      
      if (!(*clusterIt1).compareMe) continue;
      (*clusterIt1).compareMe = kFALSE;
      
      Bool_t matchFound = kFALSE;
      
      // loop over preclusters in the second list
      for (std::list<preCluster>::iterator clusterIt2 = preclusters2[iDE].begin(); clusterIt2 != preclusters2[iDE].end(); ++clusterIt2) {
        
        if (!(*clusterIt2).compareMe) continue;
        
        if (AreIdentical(*clusterIt1, *clusterIt2)) {
          
          (*clusterIt2).compareMe = kFALSE;
          
          matchFound = kTRUE;
          break;
          
        }
        
      }
      
      if (!matchFound) {
        preclustersDiff1[iDE].push_back((*clusterIt1));
        diffenceFound = kTRUE;
      }
      
    }
    
    // loop over preclusters in the second list
    for (std::list<preCluster>::iterator clusterIt2 = preclusters2[iDE].begin(); clusterIt2 != preclusters2[iDE].end(); ++clusterIt2) {
      
      if ((*clusterIt2).compareMe) {
        preclustersDiff2[iDE].push_back((*clusterIt2));
        diffenceFound = kTRUE;
      }
      
    }
    
  }
  
  return diffenceFound;
  
}

//------------------------------------------------------------------
Bool_t AreIdentical(preCluster &cl1, preCluster &cl2)
{
  /// compare 2 preclusters. Return kTRUE if they are identical
  
  if (cl1.digitId.size() != cl2.digitId.size()) return kFALSE;
  
  // loop over digitIds in both lists
  std::list<UInt_t>::iterator digitIdIt1 = cl1.digitId.begin();
  std::list<UInt_t>::iterator digitIdIt2 = cl2.digitId.begin();
  while (digitIdIt1 != cl1.digitId.end()) {
    
    if (*digitIdIt1 != *digitIdIt2) return kFALSE;
    
    digitIdIt1++;
    digitIdIt2++;
    
  }
  
  return kTRUE;
  
}

//------------------------------------------------------------------
void PrintDifferences(std::list<preCluster> *preclustersDiff1, std::list<preCluster> *preclustersDiff2,
                      Int_t &integratedNDiff1, Int_t &integratedNDiff2)
{
  /// print the preclusters found to be different in each file
  
  AliCodeTimerAutoGeneral("",0);
  
  Int_t nDiffTot1 = 0;
  Int_t nDiffTot2 = 0;
  
  // loop over DEs
  for (Int_t iDE = 1; iDE <= iDEmax; iDE++) {
    
    Int_t nDiff1 = preclustersDiff1[iDE].size();
    Int_t nDiff2 = preclustersDiff2[iDE].size();
    if (nDiff1 == 0 && nDiff2 == 0) continue;
    
    printf("\tDifferences in DE %d:\n", deIds[iDE]);
    
    // loop over preclusters in the first list
    if (nDiff1 > 0) {
      
      printf("\t\tin first file:\n");
      
      for (std::list<preCluster>::iterator clusterIt = preclustersDiff1[iDE].begin(); clusterIt != preclustersDiff1[iDE].end(); ++clusterIt) {
        
        Int_t nDigits = (*clusterIt).digitId.size();
        printf("\t\tnDigits = %d digitIDs = (", nDigits);
        
        for (std::list<UInt_t>::iterator digitIdIt = (*clusterIt).digitId.begin(); digitIdIt != (*clusterIt).digitId.end(); ++digitIdIt) {
          
          if (digitIdIt == (*clusterIt).digitId.begin()) printf("%u", *digitIdIt);
          else printf(", %u", *digitIdIt);
          
        }
        
        printf(")\n");
        
      }
      
    }
    
    // loop over preclusters in the second list
    if (nDiff2 > 0) {
      
      printf("\t\tin second file:\n");
      
      for (std::list<preCluster>::iterator clusterIt = preclustersDiff2[iDE].begin(); clusterIt != preclustersDiff2[iDE].end(); ++clusterIt) {
        
        Int_t nDigits = (*clusterIt).digitId.size();
        printf("\t\tnDigits = %d digitIDs = (", nDigits);
        
        for (std::list<UInt_t>::iterator digitIdIt = (*clusterIt).digitId.begin(); digitIdIt != (*clusterIt).digitId.end(); ++digitIdIt) {
          
          if (digitIdIt == (*clusterIt).digitId.begin()) printf("%u", *digitIdIt);
          else printf(", %u", *digitIdIt);
          
        }
        
        printf(")\n");
        
      }
      
    }
    
    nDiffTot1 += nDiff1;
    nDiffTot2 += nDiff2;
    
  }
  
  printf("\tTotal number of different preclusters = %d / %d\n", nDiffTot1, nDiffTot2);
  
  integratedNDiff1 += nDiffTot1;
  integratedNDiff2 += nDiffTot2;
  
}

//------------------------------------------------------------------
void DrawDifferences(std::list<preCluster> *preclustersDiff1, std::list<preCluster> *preclustersDiff2, const char *outFileName)
{
  /// draw the preclusters found to be different in each file
  
  AliCodeTimerAutoGeneral("",0);
  
  AliMUON2DMap digitStore1(kTRUE);
  AliMUON2DMap digitStore2(kTRUE);
  
  // loop over DEs
  for (Int_t iDE = 1; iDE <= iDEmax; iDE++) {
    
    Int_t nDiff1 = preclustersDiff1[iDE].size();
    Int_t nDiff2 = preclustersDiff2[iDE].size();
    if (nDiff1 == 0 && nDiff2 == 0) continue;
    
    // loop over preclusters in the first list
    if (nDiff1 > 0) {
      
      Int_t iCluster = 1;
      
      for (std::list<preCluster>::iterator clusterIt = preclustersDiff1[iDE].begin(); clusterIt != preclustersDiff1[iDE].end(); ++clusterIt) {
        
        // loop over attached digits
        for (std::list<UInt_t>::iterator digitIdIt = (*clusterIt).digitId.begin(); digitIdIt != (*clusterIt).digitId.end(); ++digitIdIt) {
          
          Int_t manuId = AliMUONVDigit::ManuId(*digitIdIt);
          Int_t manuChannel = AliMUONVDigit::ManuChannel(*digitIdIt);
          
          // register the digit
          AliMUONVCalibParam* c = static_cast<AliMUONVCalibParam*>(digitStore1.FindObject(deIds[iDE], manuId));
          if (!c) {
            c = new AliMUONCalibParamNI(1, AliMpConstants::ManuNofChannels(), deIds[iDE], manuId);
            digitStore1.Add(c);
          }
          c->SetValueAsInt(manuChannel, 0, iCluster);
          
        }
        
        ++iCluster;
        
      }
      
    }
    
    // loop over preclusters in the second list
    if (nDiff2 > 0) {
      
      Int_t iCluster = 1;
      
      for (std::list<preCluster>::iterator clusterIt = preclustersDiff2[iDE].begin(); clusterIt != preclustersDiff2[iDE].end(); ++clusterIt) {
        
        for (std::list<UInt_t>::iterator digitIdIt = (*clusterIt).digitId.begin(); digitIdIt != (*clusterIt).digitId.end(); ++digitIdIt) {
          
          Int_t manuId = AliMUONVDigit::ManuId(*digitIdIt);
          Int_t manuChannel = AliMUONVDigit::ManuChannel(*digitIdIt);
          
          // register the digit
          AliMUONVCalibParam* c = static_cast<AliMUONVCalibParam*>(digitStore2.FindObject(deIds[iDE], manuId));
          if (!c) {
            c = new AliMUONCalibParamNI(1, AliMpConstants::ManuNofChannels(), deIds[iDE], manuId);
            digitStore2.Add(c);
          }
          c->SetValueAsInt(manuChannel, 0, iCluster);
          
        }
        
        ++iCluster;
        
      }
      
    }
    
  }
  
  // create the tracker data
  AliMUONTrackerData digitData1("preclusters1", "preclusters1", 1, kTRUE);
  digitData1.SetDimensionName(0, "index");
  digitData1.Add(digitStore1);
  AliMUONTrackerData digitData2("preclusters2", "preclusters2", 1, kTRUE);
  digitData2.SetDimensionName(0, "index");
  digitData2.Add(digitStore2);
  
  // save it to a file
  TFile *outFile = TFile::Open(outFileName, "UPDATE");
  if (outFile && outFile->IsOpen()) {
    digitData1.Write(0x0, TObject::kOverwrite);
    digitData2.Write(0x0, TObject::kOverwrite);
    outFile->Close();
  }
  
}
