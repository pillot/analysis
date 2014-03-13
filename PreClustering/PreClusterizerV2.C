//
//  PreClusterizerV2.C
//  aliroot_dev
//
//  Created by philippe pillot on 11/03/2014.
//  Copyright (c) 2014 Philippe Pillot. All rights reserved.
//

#include <stdio.h>
#include <vector>

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
  mpPad *pads; // array of pads on both planes
  TExMap padIndices[2]; // indices of pads from their ID
  Int_t nFiredPads[2]; // number of fired pads on each plane
  std::vector<mpPad*> firedPads[2]; // indices of fired pads on each plane
  Int_t nOrderedPads[2]; // current number of fired pads in the following arrays
  std::vector<mpPad*> orderedPads[2]; // indices of fired pads ordered after preclustering and merging
};

struct preCluster {
  Int_t firstPad; // index of first associated pad in the orderedPads array
  Int_t lastPad; // index of last associated pad in the orderedPads array
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
Int_t PreClusterizeFIFO(std::vector<preCluster*> preClusters[nDEs][2], Int_t nPreClusters[nDEs][2]);
Int_t PreClusterizeRecursive(std::vector<preCluster*> preClusters[nDEs][2], Int_t nPreClusters[nDEs][2]);
void AddPad(Int_t iDE, mpPad *pad, preCluster &cl);
Int_t MergePreClusters(std::vector<preCluster*> preClusters[nDEs][2], Int_t nPreClusters[nDEs][2]);
void MergePreClusters(preCluster &cl, std::vector<preCluster*> preClusters[nDEs][2], Int_t nPreClusters[nDEs][2],
                      Int_t iDE, Int_t iPlane, preCluster *&mergedCl);
preCluster* UsePreClusters(preCluster *cl, Int_t iDE);
void MergePreClusters(preCluster &cl1, preCluster &cl2, Int_t iDE);
Bool_t AreOverlapping(preCluster &cl1, preCluster &cl2, Int_t iDE, Float_t precision);
Bool_t AreOverlapping(Float_t area1[2][2], Float_t area2[2][2], Float_t precision);
void StorePreClusters(std::vector<preCluster*> preClusters[nDEs][2], Int_t nPreClusters[nDEs][2], AliMUONVClusterStore *clusterStore);


//------------------------------------------------------------------
void PreClusterizerV2(const char *digitFileName, const char *clusterFileName, const char *method)
{
  /// Pre-clusterizer at DE level using "recursive" of "FIFO" method
  /// Each pad as an index in the pad array and knows the indices of its neighbours
  /// Compared to v1, arrays of pads and pre-clusters are reused and eventually extended but never deleted between events
  /*
   .x $ALICE_ROOT/MUON/rootlogon.C
   .L $WORK/Macros/PreClustering/O2Muon.C+
   .x $WORK/Macros/PreClustering/PreClusterizerV2.C+
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
  
  // prepare storage of preclusters
  Int_t nPreClusters[nDEs][2];
  std::vector<preCluster*> preClusters[nDEs][2];
  for (Int_t iDE = 0; iDE < nDEs; iDE++) {
    for (Int_t iPlane = 0; iPlane < 2; iPlane++) {
      preClusters[iDE][iPlane].reserve(100);
    }
  }
  
  // loop over events
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
      nPreclustersMean += PreClusterizeFIFO(preClusters, nPreClusters);
      nPreClustersMergedMean += MergePreClusters(preClusters, nPreClusters);
    }
    
    if (!strcmp(method, "both")) {
      ResetPads();
      nDigitsMean += LoadDigits(digitStore);
    }
    
    if (!strcmp(method, "recursive") || !strcmp(method, "both")) {
      nPreclustersMean += PreClusterizeRecursive(preClusters, nPreClusters);
      nPreClustersMergedMean += MergePreClusters(preClusters, nPreClusters);
    }
    
    StorePreClusters(preClusters, nPreClusters, clusterStore);
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
  
  // clean memory
  for (Int_t iDE = 0; iDE < nDEs; iDE++) {
    for (Int_t iPlane = 0; iPlane < 2; iPlane++) {
      for (Int_t iCluster = 0; iCluster < nPreClusters[iDE][iPlane]; iCluster++) {
        delete preClusters[iDE][iPlane][iCluster];
      }
    }
  }
  
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
      
      mpDEs[iDE].nPads[seg[0]->PlaneType()] = seg[0]->NofPads();
      mpDEs[iDE].nPads[seg[1]->PlaneType()] = seg[1]->NofPads();
      
      mpDEs[iDE].pads = new mpPad[mpDEs[iDE].nPads[0]+mpDEs[iDE].nPads[1]];
      mpDEs[iDE].nOrderedPads[0] = 0;
      mpDEs[iDE].orderedPads[0].reserve((mpDEs[iDE].nPads[0]/10+mpDEs[iDE].nPads[1]/10)); // 10% occupancy
      mpDEs[iDE].nOrderedPads[1] = 0;
      mpDEs[iDE].orderedPads[1].reserve((mpDEs[iDE].nPads[0]/10+mpDEs[iDE].nPads[1]/10)); // 10% occupancy
      
      // loop over cathods
      for (Int_t iCath = 0; iCath < 2; iCath++) {
        
        Int_t iPlane = seg[iCath]->PlaneType();
        mpDEs[iDE].iPlanes[iCath] = iPlane;
        mpDEs[iDE].nFiredPads[iPlane] = 0;
        mpDEs[iDE].firedPads[iPlane].reserve(mpDEs[iDE].nPads[iPlane]/10); // 10% occupancy
        
        // 1st loop over pads to associate an index to their Id
        Int_t iPad = (iPlane == 0) ? 0 : mpDEs[iDE].nPads[0];
        AliMpVPadIterator* padIt = seg[iCath]->CreateIterator();
        padIt->First();
        while (!padIt->IsDone()) {
          
          AliMpPad pad = padIt->CurrentItem();
          
          mpDEs[iDE].pads[iPad].nNeighbours = 0;
          mpDEs[iDE].pads[iPad].area[0][0] = pad.GetPositionX() - pad.GetDimensionX();
          mpDEs[iDE].pads[iPad].area[0][1] = pad.GetPositionX() + pad.GetDimensionX();
          mpDEs[iDE].pads[iPad].area[1][0] = pad.GetPositionY() - pad.GetDimensionY();
          mpDEs[iDE].pads[iPad].area[1][1] = pad.GetPositionY() + pad.GetDimensionY();
          mpDEs[iDE].pads[iPad].digit = 0x0;
          mpDEs[iDE].pads[iPad].useMe = kFALSE;
          
          padId = AliMUONVDigit::BuildUniqueID(deId, pad.GetManuId(), pad.GetManuChannel(), iCath);
          mpDEs[iDE].padIndices[iPlane].Add(padId, iPad);
          
          iPad++;
          padIt->Next();
          
        }
        
        if (neighbours) {
          
          // 2nd loop over pads to fill the neighbours
          iPad = (iPlane == 0) ? 0 : mpDEs[iDE].nPads[0];
          padIt->First();
          while (!padIt->IsDone()) {
            
            AliMpPad pad = padIt->CurrentItem();
            AliMUONVCalibParam* calibParam = static_cast<AliMUONVCalibParam*>(neighbours->FindObject(deId, pad.GetManuId()));
            
            for (Int_t iNeighbour = 1; iNeighbour < calibParam->Dimension(); iNeighbour++) {
              
              packValue2 = calibParam->ValueAsInt(pad.GetManuChannel(),iNeighbour);
              if (packValue2 <= 0) continue;
              
              calibParam->UnpackValue(packValue2, manuId2, manuChannel2);
              padId = AliMUONVDigit::BuildUniqueID(deId, manuId2, manuChannel2, iCath);
              
              mpDEs[iDE].pads[iPad].neighbours[mpDEs[iDE].pads[iPad].nNeighbours++] = mpDEs[iDE].padIndices[iPlane].GetValue(padId);
              
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
  
  Int_t firstPad = (iPlane == 0) ? 0 : mpDEs[iDE].nPads[0];
  Int_t lastPad = firstPad + mpDEs[iDE].nPads[iPlane];
  
  // loop over pads
  for (Int_t iPad1 = firstPad; iPad1 < lastPad; ++iPad1) {
    
    // loop over next pads to find the neighbours
    for (Int_t iPad2 = iPad1+1; iPad2 < lastPad; ++iPad2) {
      
      if (AreOverlapping(mpDEs[iDE].pads[iPad1].area, mpDEs[iDE].pads[iPad2].area, 1.e-4)) {
        
        mpDEs[iDE].pads[iPad1].neighbours[mpDEs[iDE].pads[iPad1].nNeighbours++] = iPad2;
        mpDEs[iDE].pads[iPad2].neighbours[mpDEs[iDE].pads[iPad2].nNeighbours++] = iPad1;
        
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
    mpPad *pad = &mpDEs[iDE].pads[mpDEs[iDE].padIndices[iPlane].GetValue(digit->GetUniqueID())];
    
    // attach the digit to this pad
    pad->digit = digit;
    pad->useMe = kTRUE;
    if (mpDEs[iDE].nFiredPads[iPlane] >= static_cast<Int_t>(mpDEs[iDE].firedPads[iPlane].size()))
      mpDEs[iDE].firedPads[iPlane].push_back(pad);
    else mpDEs[iDE].firedPads[iPlane][mpDEs[iDE].nFiredPads[iPlane]] = pad;
    ++mpDEs[iDE].nFiredPads[iPlane];
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
      
      // loop over fired pads
      for (Int_t iPad = 0; iPad < mpDEs[iDE].nFiredPads[iPlane]; ++iPad) {
        
        mpDEs[iDE].firedPads[iPlane][iPad]->digit = 0x0;
        mpDEs[iDE].firedPads[iPlane][iPad]->useMe = kFALSE;
        
      }
      
      // clear number of fired pads
      mpDEs[iDE].nFiredPads[iPlane] = 0;
      
    }
    
    // clear ordered number of fired pads
    mpDEs[iDE].nOrderedPads[0] = 0;
    mpDEs[iDE].nOrderedPads[1] = 0;
    
  }
  
}

//------------------------------------------------------------------
Int_t PreClusterizeFIFO(std::vector<preCluster*> preClusters[nDEs][2], Int_t nPreClusters[nDEs][2])
{
  /// preclusterize both planes of every DE using FIFO algorithm
  
  AliCodeTimerAutoGeneral("",0);
  
  Int_t nPreclustersTot(0);
  preCluster *cl(0x0);
  mpPad *pad(0x0), *pad2(0x0);
  
  // loop over DEs
  for (Int_t iDE = 0; iDE < nDEs; iDE++) {
    
    // loop over planes
    for (Int_t iPlane = 0; iPlane < 2; iPlane++) {
      
      nPreClusters[iDE][iPlane] = 0;
      
      // loop over fired pads
      for (Int_t iPad = 0; iPad < mpDEs[iDE].nFiredPads[iPlane]; ++iPad) {
        
        pad = mpDEs[iDE].firedPads[iPlane][iPad];
        
        if (pad->useMe) {
          
          // create the precluster if needed
          if (nPreClusters[iDE][iPlane] >= static_cast<Int_t>(preClusters[iDE][iPlane].size()))
            preClusters[iDE][iPlane].push_back(new preCluster);
          
          // get the precluster
          cl = preClusters[iDE][iPlane][nPreClusters[iDE][iPlane]];
          ++nPreClusters[iDE][iPlane];
          ++nPreclustersTot;
          
          // reset its content
          if (mpDEs[iDE].nOrderedPads[0] >= static_cast<Int_t>(mpDEs[iDE].orderedPads[0].size()))
            mpDEs[iDE].orderedPads[0].push_back(pad);
          else mpDEs[iDE].orderedPads[0][mpDEs[iDE].nOrderedPads[0]] = pad;
          cl->firstPad = mpDEs[iDE].nOrderedPads[0];
          cl->lastPad = mpDEs[iDE].nOrderedPads[0];
          ++mpDEs[iDE].nOrderedPads[0];
          cl->area[0][0] = pad->area[0][0];
          cl->area[0][1] = pad->area[0][1];
          cl->area[1][0] = pad->area[1][0];
          cl->area[1][1] = pad->area[1][1];
          cl->useMe = kTRUE;
          cl->storeMe = kFALSE;
          
          pad->useMe = kFALSE;
          
          // loop over all pads of the precluster
          for (Int_t iOrderPad = cl->firstPad; iOrderPad <= cl->lastPad; ++iOrderPad) {
            
            pad = mpDEs[iDE].orderedPads[0][iOrderPad];
            
            // loop over their neighbours
            for (Int_t iNeighbour = 0; iNeighbour < pad->nNeighbours; iNeighbour++) {
              
              pad2 = &mpDEs[iDE].pads[pad->neighbours[iNeighbour]];
              
              if (pad2->useMe) {
                
                // add the pad to the precluster
                if (mpDEs[iDE].nOrderedPads[0] >= static_cast<Int_t>(mpDEs[iDE].orderedPads[0].size()))
                  mpDEs[iDE].orderedPads[0].push_back(pad2);
                else mpDEs[iDE].orderedPads[0][mpDEs[iDE].nOrderedPads[0]] = pad2;
                cl->lastPad = mpDEs[iDE].nOrderedPads[0];
                ++mpDEs[iDE].nOrderedPads[0];
                if (pad2->area[0][0] < cl->area[0][0]) cl->area[0][0] = pad2->area[0][0];
                if (pad2->area[0][1] > cl->area[0][1]) cl->area[0][1] = pad2->area[0][1];
                if (pad2->area[1][0] < cl->area[1][0]) cl->area[1][0] = pad2->area[1][0];
                if (pad2->area[1][1] > cl->area[1][1]) cl->area[1][1] = pad2->area[1][1];
                
                pad2->useMe = kFALSE;
                
              }
              
            }
            
          }
          
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
Int_t PreClusterizeRecursive(std::vector<preCluster*> preClusters[nDEs][2], Int_t nPreClusters[nDEs][2])
{
  /// preclusterize both planes of every DE using recursive algorithm
  
  AliCodeTimerAutoGeneral("",0);
  
  Int_t nPreclustersTot = 0;
  preCluster *cl(0x0);
  
  // loop over DEs
  for (Int_t iDE = 0; iDE < nDEs; iDE++) {
    
    // loop over planes
    for (Int_t iPlane = 0; iPlane < 2; iPlane++) {
      
      nPreClusters[iDE][iPlane] = 0;
      
      // loop over fired pads
      for (Int_t iPad = 0; iPad < mpDEs[iDE].nFiredPads[iPlane]; ++iPad) {
        
        if (mpDEs[iDE].firedPads[iPlane][iPad]->useMe) {
          
          // create the precluster if needed
          if (nPreClusters[iDE][iPlane] >= static_cast<Int_t>(preClusters[iDE][iPlane].size()))
            preClusters[iDE][iPlane].push_back(new preCluster);
          
          // get the precluster
          cl = preClusters[iDE][iPlane][nPreClusters[iDE][iPlane]];
          ++nPreClusters[iDE][iPlane];
          ++nPreclustersTot;
          
          // reset its content
          cl->area[0][0] = 1.e6;
          cl->area[0][1] = -1.e6;
          cl->area[1][0] = 1.e6;
          cl->area[1][1] = -1.e6;
          cl->useMe = kTRUE;
          cl->storeMe = kFALSE;
          
          // add the pad and its fired neighbours recusively
          cl->firstPad = mpDEs[iDE].nOrderedPads[0];
          AddPad(iDE, mpDEs[iDE].firedPads[iPlane][iPad], *cl);
          
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
void AddPad(Int_t iDE, mpPad *pad, preCluster &cl)
{
  /// add the given mpPad and its fired neighbours (recursive method)
  
  // add the given pad
  if (mpDEs[iDE].nOrderedPads[0] >= static_cast<Int_t>(mpDEs[iDE].orderedPads[0].size()))
    mpDEs[iDE].orderedPads[0].push_back(pad);
  else mpDEs[iDE].orderedPads[0][mpDEs[iDE].nOrderedPads[0]] = pad;
  cl.lastPad = mpDEs[iDE].nOrderedPads[0];
  ++mpDEs[iDE].nOrderedPads[0];
  if (pad->area[0][0] < cl.area[0][0]) cl.area[0][0] = pad->area[0][0];
  if (pad->area[0][1] > cl.area[0][1]) cl.area[0][1] = pad->area[0][1];
  if (pad->area[1][0] < cl.area[1][0]) cl.area[1][0] = pad->area[1][0];
  if (pad->area[1][1] > cl.area[1][1]) cl.area[1][1] = pad->area[1][1];
  
  pad->useMe = kFALSE;
  
  // loop over its neighbours
  for (Int_t iNeighbour = 0; iNeighbour < pad->nNeighbours; iNeighbour++) {
    
    if (mpDEs[iDE].pads[pad->neighbours[iNeighbour]].useMe) {
      
      // add the pad to the precluster
      AddPad(iDE, &mpDEs[iDE].pads[pad->neighbours[iNeighbour]], cl);
      
    }
    
  }
  
}

//------------------------------------------------------------------
Int_t MergePreClusters(std::vector<preCluster*> preClusters[nDEs][2], Int_t nPreClusters[nDEs][2])
{
  /// merge overlapping preclusters on every DE
  
  AliCodeTimerAutoGeneral("",0);
  
  Int_t nPreClustersMergedTot = 0;
  
  // loop over DEs
  for (Int_t iDE = 0; iDE < nDEs; iDE++) {
    
    // loop over preclusters of one plane
    for (Int_t iCluster = 0; iCluster < nPreClusters[iDE][0]; ++iCluster) {
      
      if (!preClusters[iDE][0][iCluster]->useMe) continue;
      
      preClusters[iDE][0][iCluster]->useMe = kFALSE;
      
      // look for overlapping preclusters in the other plane
      preCluster *mergedCl(0x0);
      MergePreClusters(*preClusters[iDE][0][iCluster], preClusters, nPreClusters, iDE, 1, mergedCl);
      
      // add the current one
      if (!mergedCl) mergedCl = UsePreClusters(preClusters[iDE][0][iCluster], iDE);
      else MergePreClusters(*mergedCl, *preClusters[iDE][0][iCluster], iDE);
      
      ++nPreClustersMergedTot;
      
    }
    
    // loop over preclusters of the other plane
    for (Int_t iCluster = 0; iCluster < nPreClusters[iDE][1]; ++iCluster) {
      
      if (!preClusters[iDE][1][iCluster]->useMe) continue;
      
      // all remaining preclusters have to be stored
      UsePreClusters(preClusters[iDE][1][iCluster], iDE);
      
      ++nPreClustersMergedTot;
      
    }
    
  }
  
#ifdef VERBOSEMODE
  printf("\tTotal number of preclusters after merging = %d\n", nPreClustersMergedTot);
#endif
  
  return nPreClustersMergedTot;
  
}

//------------------------------------------------------------------
void MergePreClusters(preCluster &cl, std::vector<preCluster*> preClusters[nDEs][2], Int_t nPreClusters[nDEs][2],
                      Int_t iDE, Int_t iPlane, preCluster *&mergedCl)
{
  /// merge preclusters on the given plane overlapping with the given one (recursive method)
  
  // loop over preclusters in the given plane
  for (Int_t iCluster = 0; iCluster < nPreClusters[iDE][iPlane]; ++iCluster) {
    
    if (!preClusters[iDE][iPlane][iCluster]->useMe) continue;
    
    if (AreOverlapping(cl.area, preClusters[iDE][iPlane][iCluster]->area, -1.e-4) &&
        AreOverlapping(cl, *preClusters[iDE][iPlane][iCluster], iDE, -1.e-4)) {
      
      preClusters[iDE][iPlane][iCluster]->useMe = kFALSE;
      
      // look for new overlapping preclusters in the other plane
      MergePreClusters(*preClusters[iDE][iPlane][iCluster], preClusters, nPreClusters, iDE, (iPlane+1)%2, mergedCl);
      
      // store overlapping preclusters and merge them
      if (!mergedCl) mergedCl = UsePreClusters(preClusters[iDE][iPlane][iCluster], iDE);
      else MergePreClusters(*mergedCl, *preClusters[iDE][iPlane][iCluster], iDE);
      
    }
    
  }
  
}

//------------------------------------------------------------------
preCluster* UsePreClusters(preCluster *cl, Int_t iDE)
{
  /// use this precluster as a new merged precluster
  
  Int_t firstPad = mpDEs[iDE].nOrderedPads[1];
  
  // move the fired pads
  for (Int_t iOrderPad = cl->firstPad; iOrderPad <= cl->lastPad; ++iOrderPad) {
    
    if (mpDEs[iDE].nOrderedPads[1] >= static_cast<Int_t>(mpDEs[iDE].orderedPads[1].size()))
      mpDEs[iDE].orderedPads[1].push_back(mpDEs[iDE].orderedPads[0][iOrderPad]);
    else mpDEs[iDE].orderedPads[1][mpDEs[iDE].nOrderedPads[1]] = mpDEs[iDE].orderedPads[0][iOrderPad];
    
    ++mpDEs[iDE].nOrderedPads[1];
    
  }
  
  cl->firstPad = firstPad;
  cl->lastPad = mpDEs[iDE].nOrderedPads[1] - 1;
  
  cl->storeMe = kTRUE;
  
  return cl;
  
}

//------------------------------------------------------------------
void MergePreClusters(preCluster &cl1, preCluster &cl2, Int_t iDE)
{
  /// merge precluster2 into precluster1
  
  // move the fired pads
  for (Int_t iOrderPad = cl2.firstPad; iOrderPad <= cl2.lastPad; ++iOrderPad) {
    
    if (mpDEs[iDE].nOrderedPads[1] >= static_cast<Int_t>(mpDEs[iDE].orderedPads[1].size()))
      mpDEs[iDE].orderedPads[1].push_back(mpDEs[iDE].orderedPads[0][iOrderPad]);
    else mpDEs[iDE].orderedPads[1][mpDEs[iDE].nOrderedPads[1]] = mpDEs[iDE].orderedPads[0][iOrderPad];
    
    ++mpDEs[iDE].nOrderedPads[1];
    
  }
  
  cl1.lastPad = mpDEs[iDE].nOrderedPads[1] - 1;
  
}

//___________________________________________________________________________
Bool_t AreOverlapping(preCluster &cl1, preCluster &cl2, Int_t iDE, Float_t precision)
{
  /// check if the two preclusters overlap
  /// precision in cm: positive = increase pad size / negative = decrease pad size
  
  // loop over all pads of the precluster1
  for (Int_t iOrderPad1 = cl1.firstPad; iOrderPad1 <= cl1.lastPad; ++iOrderPad1) {
    
    // loop over all pads of the precluster2
    for (Int_t iOrderPad2 = cl2.firstPad; iOrderPad2 <= cl2.lastPad; ++iOrderPad2) {
      
      if (AreOverlapping(mpDEs[iDE].orderedPads[0][iOrderPad1]->area,
                         mpDEs[iDE].orderedPads[0][iOrderPad2]->area, precision)) return kTRUE;
      
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
void StorePreClusters(std::vector<preCluster*> preClusters[nDEs][2], Int_t nPreClusters[nDEs][2], AliMUONVClusterStore *clusterStore)
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
      for (Int_t iCluster = 0; iCluster < nPreClusters[iDE][iPlane]; ++iCluster) {
        
        if (!preClusters[iDE][iPlane][iCluster]->storeMe) continue;
        
        AliMUONVCluster *cluster = clusterStore->Add(mpDEs[iDE].id/100-1, mpDEs[iDE].id, clusterIndex++);
        
        // loop over pads
        for (Int_t iOrderPad = preClusters[iDE][iPlane][iCluster]->firstPad; iOrderPad <= preClusters[iDE][iPlane][iCluster]->lastPad; ++iOrderPad) {
          
          cluster->AddDigitId(mpDEs[iDE].orderedPads[1][iOrderPad]->digit->GetUniqueID());
          
        }
        
      }
      
    }
    
  }
  
}

