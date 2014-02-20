//
//  PreClusterizerV1.C
//  aliroot_dev
//
//  Created by philippe pillot on 14/01/2014.
//  Copyright (c) 2014 Philippe Pillot. All rights reserved.
//

#include <stdio.h>
#include <list>

#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>

#include "AliCDBManager.h"
#include "AliCodeTimer.h"

#include "AliMUONCDB.h"
#include "AliMUONVCalibParam.h"
#include "AliMUONCalibrationData.h"
#include "AliMUONConstants.h"
#include "AliMUONVDigit.h"
#include "AliMUONVDigitStore.h"
#include "AliMUONVCluster.h"
#include "AliMUONClusterStoreV2.h"
#include "AliMUON2DMap.h"

#include "AliMpDEIterator.h"
#include "AliMpVSegmentation.h"
#include "AliMpSegmentation.h"
#include "AliMpVPadIterator.h"
#include "AliMpPad.h"

#include "/Users/pillot/Work/Alice/Macros/PreClustering/O2Muon.h"

//#define VERBOSEMODE

struct mpPad {
  Int_t nNeighbours; // number of neighbours
  Int_t neighbours[10]; // indices of neighbours in array stored in mpDE
  Float_t area[2][2]; // 2D area
  AliMUONVDigit *digit; // pointer to the associated digit
  Bool_t useMe; // kFALSE if no digit attached or already visited
};

struct mpDE {
  Int_t id; // unique ID
  Int_t iPlanes[2]; // plane type corresponding to both cathods
  Int_t nPads[2]; // number of pads on each plane
  mpPad *pads[2]; // array of pads on each plane
  TExMap padIndices[2]; // indices of pads from their ID
  std::list<Int_t> iDigits[2]; // indices of reconstructed digits on each plane
};

struct preCluster {
  std::list<mpPad*> pads;
  Float_t area[2][2]; // 2D area containing the precluster
  Bool_t useMe; // kFALSE if precluster already merged to another one
  Bool_t storeMe; // kTRUE if precluster to be saved (merging result)
};


const Int_t nDEs = 156;
mpDE mpDEs[nDEs];
TExMap deIndices;


void CreateMapping(AliMUONVStore* neighbours);
void FindNeighbours(Int_t iDE, Int_t iPlane);
Int_t LoadDigits(AliMUONVDigitStore* digitStore);
void ResetPads();
Int_t PreClusterizeFIFO(std::list<preCluster> preClusters[nDEs][2]);
Int_t PreClusterizeRecursive(std::list<preCluster> preClusters[nDEs][2]);
void AddPad(Int_t iDE, Int_t iPlane, mpPad &pad, preCluster &cl);
Int_t MergePreClusters(std::list<preCluster> preClusters[nDEs][2]);
void MergePreClusters(preCluster &cl, std::list<preCluster> preClusters[2], Int_t iPlane, preCluster *&mergedCl);
void MergePreClusters(preCluster &cl1, preCluster &cl2);
Bool_t AreOverlapping(preCluster &cl1, preCluster &cl2, Float_t precision);
Bool_t AreOverlapping(Float_t area1[2][2], Float_t area2[2][2], Float_t precision);
void StorePreClusters(std::list<preCluster> preClusters[nDEs][2], AliMUONVClusterStore *clusterStore);


//------------------------------------------------------------------
void PreClusterizerV1(const char *digitFileName, const char *clusterFileName, const char *method)
{
  /// Pre-clusterizer at DE level using "recursive" of "FIFO" method
  /// Each pad as an index in the pad array and knows the indices of its neighbours
  /*
   .x $ALICE_ROOT/MUON/rootlogon.C
   .L $WORK/Macros/PreClustering/O2Muon.C+
   .x $WORK/Macros/PreClustering/PreClusterizerV1.C+
  */
  
  if (strcmp(method, "recursive") && strcmp(method, "FIFO") && strcmp(method, "both")) {
    printf("Choose either \"recursive\", \"FIFO\" or \"both\" method.\n");
    printf("If you choose \"both\" only \"recursive\" results will be stored.\n");
    return;
  }
  
  /*
  // create digits
  O2Muon digitizer("raw://");
  //  digitizer.makeDigitFile("/Users/pillot/Work/Alice/Work/Data/2012/LHC12h/raw/190147/local/12000190147002.10.root");
  digitizer.makeDigitFile("/Users/pillot/Work/Alice/Work/Data/2013/LHC13f/raw/196474/13000196474082.15filtered.root");
  
  // standard clusterizing
  O2Muon clusterizer("raw://");
  digitizer.makeClustering("digits.root", "PRECLUSTER", "precluster.log", 196474);
  return;
  */
  // load mapping locally and create the local structure
  AliCDBManager* man = AliCDBManager::Instance();
//  man->SetDefaultStorage("local:///Users/pillot/Work/Alice/Work/Data/2012/LHC12h/raw/190147/saf/CDBMirror/alice/data/2012/OCDB");
//  man->SetRun(190147);
  man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  man->SetRun(0);
  if (!AliMUONCDB::LoadMapping()) return;
  //AliMUONVStore* neighbours = AliMUONCalibrationData::CreateNeighbours(0);
  //CreateMapping(neighbours);
  CreateMapping(0x0);
  //delete neighbours;
  
  // read digits
  TFile* digitFile = TFile::Open(digitFileName);
  if (!digitFile) return;
  TTree* treeD = static_cast<TTree*>(digitFile->Get("TreeD"));
  if (!treeD) return;
  AliMUONVDigitStore* digitStore = AliMUONVDigitStore::Create(*treeD);
  digitStore->Connect(*treeD);
  
  // output clusters
  TFile* clusterFile = new TFile(clusterFileName,"RECREATE");
  TTree* treeR = new TTree("TreeR","Clusters");
  AliMUONVClusterStore *clusterStore = new AliMUONClusterStoreV2();
  clusterStore->Connect(*treeR, kTRUE);
  
  // loop over events
  std::list<preCluster> preClusters[nDEs][2];
  Float_t nDigitsMean = 0.;
  Float_t nPreclustersMean = 0.;
  Float_t nPreClustersMergedMean = 0.;
  Long64_t nEvents = treeD->GetEntries();
  for (Long64_t iEv = 0; iEv < nEvents; ++iEv) {
    
#ifdef VERBOSEMODE
    printf("Event %lld:\n", iEv);
#endif
    
    treeD->GetEntry(iEv);
    nDigitsMean += LoadDigits(digitStore);
    
    if (!strcmp(method, "FIFO") || !strcmp(method, "both")) {
      nPreclustersMean += PreClusterizeFIFO(preClusters);
      nPreClustersMergedMean += MergePreClusters(preClusters);
    }
    
    if (!strcmp(method, "both")) {
      ResetPads();
      nDigitsMean += LoadDigits(digitStore);
    }
    
    if (!strcmp(method, "recursive") || !strcmp(method, "both")) {
      nPreclustersMean += PreClusterizeRecursive(preClusters);
      nPreClustersMergedMean += MergePreClusters(preClusters);
    }
    
    StorePreClusters(preClusters, clusterStore);
    treeR->Fill();
    
    ResetPads();
    digitStore->Clear();
    
  }
  
  if (!strcmp(method, "both")) {
    nDigitsMean /= 2.;
    nPreclustersMean /= 2.;
    nPreClustersMergedMean /= 2.;
  }
  
  printf("\n");
  printf("Average number of digits = %f\n", nDigitsMean/nEvents);
  printf("Average number of preclusters before merging = %f\n", nPreclustersMean/nEvents);
  printf("Average number of preclusters after merging = %f\n", nPreClustersMergedMean/nEvents);
  printf("\n");
  
  // save clusters
  clusterFile->cd();
  treeR->Write();
  clusterFile->Close();
  delete clusterFile;
  
  AliCodeTimer::Instance()->Print();
  
}

//------------------------------------------------------------------
void CreateMapping(AliMUONVStore* neighbours)
{
  /// fill the mpDE and mpPad structures once for all
  
  AliCodeTimerAutoGeneral("",0);
  
  MemInfo_t memBefore, memAfter;
  ProcInfo_t procBefore, procAfter;
  gSystem->GetMemInfo(&memBefore);
  gSystem->GetProcInfo(&procBefore);
  
  // loop over DEs
  Int_t iDE = 0;
  Int_t packValue2 = -1, manuId2 = -1, manuChannel2 = -1;
  UInt_t padId = 0;
  for (Int_t iCh = 0; iCh < AliMUONConstants::NTrackingCh(); iCh++) {
    AliMpDEIterator deIt;
    deIt.First(iCh);
    while (!deIt.IsDone()) {
      
      Int_t deId = deIt.CurrentDEId();
      mpDEs[iDE].id = deId;
      deIndices.Add(deId, iDE);
      
      const AliMpVSegmentation* seg[2] =
      { AliMpSegmentation::Instance()->GetMpSegmentation(deId,AliMp::kCath0),
        AliMpSegmentation::Instance()->GetMpSegmentation(deId,AliMp::kCath1)
      };
      
      // loop over cathods
      for (Int_t iCath = 0; iCath < 2; iCath++) {
        
        Int_t iPlane = seg[iCath]->PlaneType();
        mpDEs[iDE].iPlanes[iCath] = iPlane;
        Int_t nPads = seg[iCath]->NofPads();
        mpDEs[iDE].nPads[iPlane] = nPads;
        mpDEs[iDE].pads[iPlane] = new mpPad[nPads];
        
        // 1st loop over pads to associate an index to their Id
        Int_t iPad = 0;
        AliMpVPadIterator* padIt = seg[iCath]->CreateIterator();
        padIt->First();
        while (!padIt->IsDone()) {
          
          AliMpPad pad = padIt->CurrentItem();
          
          mpDEs[iDE].pads[iPlane][iPad].nNeighbours = 0;
          mpDEs[iDE].pads[iPlane][iPad].area[0][0] = pad.GetPositionX() - pad.GetDimensionX();
          mpDEs[iDE].pads[iPlane][iPad].area[0][1] = pad.GetPositionX() + pad.GetDimensionX();
          mpDEs[iDE].pads[iPlane][iPad].area[1][0] = pad.GetPositionY() - pad.GetDimensionY();
          mpDEs[iDE].pads[iPlane][iPad].area[1][1] = pad.GetPositionY() + pad.GetDimensionY();
          mpDEs[iDE].pads[iPlane][iPad].digit = 0x0;
          mpDEs[iDE].pads[iPlane][iPad].useMe = kFALSE;
          
          padId = AliMUONVDigit::BuildUniqueID(deId, pad.GetManuId(), pad.GetManuChannel(), iCath);
          mpDEs[iDE].padIndices[iPlane].Add(padId, iPad);
          
          iPad++;
          padIt->Next();
          
        }
        
        if (neighbours) {
          
          // 2nd loop over pads to fill the neighbours
          iPad = 0;
          padIt->First();
          while (!padIt->IsDone()) {
            
            AliMpPad pad = padIt->CurrentItem();
            AliMUONVCalibParam* calibParam = static_cast<AliMUONVCalibParam*>(neighbours->FindObject(deId, pad.GetManuId()));
            
            for (Int_t iNeighbour = 1; iNeighbour < calibParam->Dimension(); iNeighbour++) {
              
              packValue2 = calibParam->ValueAsInt(pad.GetManuChannel(),iNeighbour);
              if (packValue2 <= 0) continue;
              
              calibParam->UnpackValue(packValue2, manuId2, manuChannel2);
              padId = AliMUONVDigit::BuildUniqueID(deId, manuId2, manuChannel2, iCath);
              
              mpDEs[iDE].pads[iPlane][iPad].neighbours[mpDEs[iDE].pads[iPlane][iPad].nNeighbours++] = mpDEs[iDE].padIndices[iPlane].GetValue(padId);
              
            }
            
            iPad++;
            padIt->Next();
            
          }
          
        } else FindNeighbours(iDE, iPlane);
        
      }
      
      iDE++;
      deIt.Next();
      
    }
    
  }
  
  gSystem->GetMemInfo(&memAfter);
  gSystem->GetProcInfo(&procAfter);
  printf("Memory comsumption for mapping: %d MB\n", memAfter.fMemUsed - memBefore.fMemUsed);
  printf("Memory comsumption for mapping: %d MB\n", memBefore.fMemFree - memAfter.fMemFree);
  printf("Memory comsumption for mapping: %ld MB\n", (procAfter.fMemResident - procBefore.fMemResident + 500) / 1000);
  printf("Memory comsumption for mapping: %ld MB\n", (procAfter.fMemVirtual - procBefore.fMemVirtual + 500) / 1000);
  
}

//------------------------------------------------------------------
void FindNeighbours(Int_t iDE, Int_t iPlane)
{
  /// fill the mpPad neighbours structures of given DE/plane
  
  AliCodeTimerAutoGeneral("",0);
  
  // loop over pads
  for (Int_t iPad1 = 0; iPad1 < mpDEs[iDE].nPads[iPlane]; ++iPad1) {
    
    // loop over next pads to find the neighbours
    for (Int_t iPad2 = iPad1+1; iPad2 < mpDEs[iDE].nPads[iPlane]; ++iPad2) {
      
      if (AreOverlapping(mpDEs[iDE].pads[iPlane][iPad1].area, mpDEs[iDE].pads[iPlane][iPad2].area, 1.e-4)) {
        
        mpDEs[iDE].pads[iPlane][iPad1].neighbours[mpDEs[iDE].pads[iPlane][iPad1].nNeighbours++] = iPad2;
        mpDEs[iDE].pads[iPlane][iPad2].neighbours[mpDEs[iDE].pads[iPlane][iPad2].nNeighbours++] = iPad1;
        
      }
      
    }
    
  }
  
}

//------------------------------------------------------------------
Int_t LoadDigits(AliMUONVDigitStore* digitStore)
{
  /// fill the mpDE structure with reconstructed digits
  
  AliCodeTimerAutoGeneral("",0);
  
  Int_t nDigits = 0;
  
  // loop over digits
  AliMUONVDigit* digit = 0x0;
  TIter nextDigit(digitStore->CreateIterator());
  while ((digit = static_cast<AliMUONVDigit*>(nextDigit()))) {
    
    if (digit->Charge() <= 0) continue;
    
    Int_t iDE = deIndices.GetValue(digit->DetElemId());
    Int_t iPlane = mpDEs[iDE].iPlanes[digit->Cathode()];
    Int_t iPad = mpDEs[iDE].padIndices[iPlane].GetValue(digit->GetUniqueID());
    
    // attach the digit to this pad
    mpDEs[iDE].pads[iPlane][iPad].digit = digit;
    mpDEs[iDE].pads[iPlane][iPad].useMe = kTRUE;
    mpDEs[iDE].iDigits[iPlane].push_back(iPad);
    ++nDigits;
    
  }
  
#ifdef VERBOSEMODE
  printf("\tTotal number of digits = %d\n", nDigits);
#endif
  
  return nDigits;
  
}

//------------------------------------------------------------------
void ResetPads()
{
  /// reset digit information in concerned pads
  
  AliCodeTimerAutoGeneral("",0);
  
  // loop over DEs
  for (Int_t iDE = 0; iDE < nDEs; iDE++) {
    
    // loop over planes
    for (Int_t iPlane = 0; iPlane < 2; iPlane++) {
      
      // loop over fired digits
      for (std::list<Int_t>::iterator digitIt = mpDEs[iDE].iDigits[iPlane].begin(); digitIt != mpDEs[iDE].iDigits[iPlane].end(); ++digitIt) {
        
        mpDEs[iDE].pads[iPlane][*digitIt].digit = 0x0;
        mpDEs[iDE].pads[iPlane][*digitIt].useMe = kFALSE;
        
      }
      
      // clear list of digits
      mpDEs[iDE].iDigits[iPlane].clear();
      
    }
    
  }
  
}

//------------------------------------------------------------------
Int_t PreClusterizeFIFO(std::list<preCluster> preClusters[nDEs][2])
{
  /// preclusterize both planes of every DE using FIFO algorithm
  
  AliCodeTimerAutoGeneral("",0);
  
  Int_t nPreclustersTot = 0;
  
  // loop over DEs
  for (Int_t iDE = 0; iDE < nDEs; iDE++) {
    
    // loop over planes
    for (Int_t iPlane = 0; iPlane < 2; iPlane++) {
      
      preClusters[iDE][iPlane].clear();
      
      // loop over fired digits
      for (std::list<Int_t>::iterator digitIt = mpDEs[iDE].iDigits[iPlane].begin(); digitIt != mpDEs[iDE].iDigits[iPlane].end(); ++digitIt) {
        
        if (mpDEs[iDE].pads[iPlane][*digitIt].useMe) {
          
          // create the precluster
          preCluster cl;
          cl.pads.push_back(&mpDEs[iDE].pads[iPlane][*digitIt]);
          cl.area[0][0] = mpDEs[iDE].pads[iPlane][*digitIt].area[0][0];
          cl.area[0][1] = mpDEs[iDE].pads[iPlane][*digitIt].area[0][1];
          cl.area[1][0] = mpDEs[iDE].pads[iPlane][*digitIt].area[1][0];
          cl.area[1][1] = mpDEs[iDE].pads[iPlane][*digitIt].area[1][1];
          cl.useMe = kTRUE;
          cl.storeMe = kFALSE;
          mpDEs[iDE].pads[iPlane][*digitIt].useMe = kFALSE;
          
          // loop over all pads of the precluster
          for (std::list<mpPad*>::iterator padIt = cl.pads.begin(); padIt != cl.pads.end(); ++padIt) {
            
            // loop over their neighbours
            for (Int_t iNeighbour = 0; iNeighbour < (*padIt)->nNeighbours; iNeighbour++) {
              
              Int_t iPad = (*padIt)->neighbours[iNeighbour];
              if (mpDEs[iDE].pads[iPlane][iPad].useMe) {
                
                // add the pad to the precluster
                cl.pads.push_back(&mpDEs[iDE].pads[iPlane][iPad]);
                if (mpDEs[iDE].pads[iPlane][iPad].area[0][0] < cl.area[0][0]) cl.area[0][0] = mpDEs[iDE].pads[iPlane][iPad].area[0][0];
                if (mpDEs[iDE].pads[iPlane][iPad].area[0][1] > cl.area[0][1]) cl.area[0][1] = mpDEs[iDE].pads[iPlane][iPad].area[0][1];
                if (mpDEs[iDE].pads[iPlane][iPad].area[1][0] < cl.area[1][0]) cl.area[1][0] = mpDEs[iDE].pads[iPlane][iPad].area[1][0];
                if (mpDEs[iDE].pads[iPlane][iPad].area[1][1] > cl.area[1][1]) cl.area[1][1] = mpDEs[iDE].pads[iPlane][iPad].area[1][1];
                mpDEs[iDE].pads[iPlane][iPad].useMe = kFALSE;
                
              }
              
            }
            
          }
          
          // add the precluster to the list
          preClusters[iDE][iPlane].push_back(cl);
          ++nPreclustersTot;
          
        }
        
      }
      
    }
    
  }
  
#ifdef VERBOSEMODE
  printf("\tTotal number of preclusters before merging = %d\n", nPreclustersTot);
#endif
  
  return nPreclustersTot;
  
}

//------------------------------------------------------------------
Int_t PreClusterizeRecursive(std::list<preCluster> preClusters[nDEs][2])
{
  /// preclusterize both planes of every DE using recursive algorithm
  
  AliCodeTimerAutoGeneral("",0);
  
  Int_t nPreclustersTot = 0;
  
  // loop over DEs
  for (Int_t iDE = 0; iDE < nDEs; iDE++) {
    
    // loop over planes
    for (Int_t iPlane = 0; iPlane < 2; iPlane++) {
      
      preClusters[iDE][iPlane].clear();
      
      // loop over fired digits
      for (std::list<Int_t>::iterator digitIt = mpDEs[iDE].iDigits[iPlane].begin(); digitIt != mpDEs[iDE].iDigits[iPlane].end(); ++digitIt) {
        
        if (mpDEs[iDE].pads[iPlane][*digitIt].useMe) {
          
          // create the precluster
          preCluster cl;
          cl.area[0][0] = 1.e6;
          cl.area[0][1] = -1.e6;
          cl.area[1][0] = 1.e6;
          cl.area[1][1] = -1.e6;
          cl.useMe = kTRUE;
          cl.storeMe = kFALSE;
          
          // add the pad and its fired neighbours recusively
          AddPad(iDE, iPlane, mpDEs[iDE].pads[iPlane][*digitIt], cl);
          
          // add the precluster to the list
          preClusters[iDE][iPlane].push_back(cl);
          ++nPreclustersTot;
          
        }
        
      }
      
    }
    
  }
  
#ifdef VERBOSEMODE
  printf("\tTotal number of preclusters before merging = %d\n", nPreclustersTot);
#endif
  
  return nPreclustersTot;
  
}

//------------------------------------------------------------------
void AddPad(Int_t iDE, Int_t iPlane, mpPad &pad, preCluster &cl)
{
  /// add the given mpPad and its fired neighbours (recursive method)
  
  // add the given pad
  cl.pads.push_back(&pad);
  if (pad.area[0][0] < cl.area[0][0]) cl.area[0][0] = pad.area[0][0];
  if (pad.area[0][1] > cl.area[0][1]) cl.area[0][1] = pad.area[0][1];
  if (pad.area[1][0] < cl.area[1][0]) cl.area[1][0] = pad.area[1][0];
  if (pad.area[1][1] > cl.area[1][1]) cl.area[1][1] = pad.area[1][1];
  pad.useMe = kFALSE;
  
  // loop over its neighbours
  for (Int_t iNeighbour = 0; iNeighbour < pad.nNeighbours; iNeighbour++) {
    
    if (mpDEs[iDE].pads[iPlane][pad.neighbours[iNeighbour]].useMe) {
      
      // add the pad to the precluster
      AddPad(iDE, iPlane, mpDEs[iDE].pads[iPlane][pad.neighbours[iNeighbour]], cl);
      
    }
    
  }
  
}

//------------------------------------------------------------------
Int_t MergePreClusters(std::list<preCluster> preClusters[nDEs][2])
{
  /// merge overlapping preclusters on every DE
  
  AliCodeTimerAutoGeneral("",0);
  
  Int_t nPreClustersMergedTot = 0;
  
  // loop over DEs
  for (Int_t iDE = 0; iDE < nDEs; iDE++) {
    
    // loop over preclusters of one plane
    for (std::list<preCluster>::iterator clusterIt = preClusters[iDE][0].begin(); clusterIt != preClusters[iDE][0].end(); ++clusterIt) {
      
      if (!(*clusterIt).useMe) continue;
      
      (*clusterIt).useMe = kFALSE;
      
      // look for overlapping preclusters in the other plane
      preCluster *mergedCl(0x0);
      MergePreClusters(*clusterIt, preClusters[iDE], 1, mergedCl);
      
      // add the current one
      if (!mergedCl) mergedCl = &(*clusterIt);
      else MergePreClusters(*mergedCl, *clusterIt);
      
      mergedCl->storeMe = kTRUE;
      
      ++nPreClustersMergedTot;
      
    }
    
    // loop over preclusters of the other plane
    for (std::list<preCluster>::iterator clusterIt = preClusters[iDE][1].begin(); clusterIt != preClusters[iDE][1].end(); ++clusterIt) {
      
      if (!(*clusterIt).useMe) continue;
      
      // all remaining preclusters have to be stored
      (*clusterIt).storeMe = kTRUE;
      
      ++nPreClustersMergedTot;
      
    }
    
  }
  
#ifdef VERBOSEMODE
  printf("\tTotal number of preclusters after merging = %d\n", nPreClustersMergedTot);
#endif
  
  return nPreClustersMergedTot;
  
}

//------------------------------------------------------------------
void MergePreClusters(preCluster &cl, std::list<preCluster> preClusters[2], Int_t iPlane, preCluster *&mergedCl)
{
  /// merge preclusters on the given plane overlapping with the given one (recursive method)
  
  // loop over preclusters in the given plane
  for (std::list<preCluster>::iterator clusterIt = preClusters[iPlane].begin(); clusterIt != preClusters[iPlane].end(); ++clusterIt) {
    
    if (!(*clusterIt).useMe) continue;
    
    if (AreOverlapping(cl.area, (*clusterIt).area, -1.e-4) && AreOverlapping(cl, *clusterIt, -1.e-4)) {
      
      (*clusterIt).useMe = kFALSE;
      
      // look for new overlapping preclusters in the other plane
      MergePreClusters(*clusterIt, preClusters, (iPlane+1)%2, mergedCl);
      
      // store overlapping preclusters and merge them
      if (!mergedCl) mergedCl = &(*clusterIt);
      else MergePreClusters(*mergedCl, *clusterIt);
      
    }
    
  }
  
}

//------------------------------------------------------------------
void MergePreClusters(preCluster &cl1, preCluster &cl2)
{
  /// merge precluster2 into precluster1
  
  cl1.pads.splice(cl1.pads.end(), cl2.pads);
  /*
  if (cl2.area[0][0] < cl1.area[0][0]) cl1.area[0][0] = cl2.area[0][0];
  if (cl2.area[0][1] > cl1.area[0][1]) cl1.area[0][1] = cl2.area[0][1];
  if (cl2.area[1][0] < cl1.area[1][0]) cl1.area[1][0] = cl2.area[1][0];
  if (cl2.area[1][1] > cl1.area[1][1]) cl1.area[1][1] = cl2.area[1][1];
  */
}

//___________________________________________________________________________
Bool_t AreOverlapping(preCluster &cl1, preCluster &cl2, Float_t precision)
{
  /// check if the two preclusters overlap
  /// precision in cm: positive = increase pad size / negative = decrease pad size
  
  // loop over all pads of the precluster1
  for (std::list<mpPad*>::iterator padIt1 = cl1.pads.begin(); padIt1 != cl1.pads.end(); ++padIt1) {
    
    // loop over all pads of the precluster2
    for (std::list<mpPad*>::iterator padIt2 = cl2.pads.begin(); padIt2 != cl2.pads.end(); ++padIt2) {
      
      if (AreOverlapping((*padIt1)->area, (*padIt2)->area, precision)) return kTRUE;
      
    }
    
  }
  
  return kFALSE;
  
}

//------------------------------------------------------------------
Bool_t AreOverlapping(Float_t area1[2][2], Float_t area2[2][2], Float_t precision)
{
  /// check if the two areas overlap
  /// precision in cm: positive = increase pad size / negative = decrease pad size
  
  return (!((area1[0][0] - area2[0][1] > precision) || (area2[0][0] - area1[0][1] > precision) ||
            (area1[1][0] - area2[1][1] > precision) || (area2[1][0] - area1[1][1] > precision)));
  
}

//------------------------------------------------------------------
void StorePreClusters(std::list<preCluster> preClusters[nDEs][2], AliMUONVClusterStore *clusterStore)
{
  /// store the preclusters into standard cluster store
  
  AliCodeTimerAutoGeneral("",0);
  
  clusterStore->Clear();
  Int_t clusterIndex = 0;
  
  // loop over DEs
  for (Int_t iDE = 0; iDE < nDEs; iDE++) {
    
    // loop over planes
    for (Int_t iPlane = 0; iPlane < 2; iPlane++) {
      
      // loop over preclusters
      for (std::list<preCluster>::iterator clusterIt = preClusters[iDE][iPlane].begin(); clusterIt != preClusters[iDE][iPlane].end(); ++clusterIt) {
        
        if (!(*clusterIt).storeMe) continue;
        
        AliMUONVCluster *cluster = clusterStore->Add(mpDEs[iDE].id/100-1, mpDEs[iDE].id, clusterIndex++);
        
        // loop over pads
        for (std::list<mpPad*>::iterator padIt = (*clusterIt).pads.begin(); padIt != (*clusterIt).pads.end(); ++padIt) {
          
          cluster->AddDigitId((*padIt)->digit->GetUniqueID());
          
        }
        
      }
      
    }
    
  }
  
}

