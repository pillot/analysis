//
//  PreClusterizerV2.cxx
//  aliroot_dev
//
//  Created by philippe pillot on 13/03/2014.
//  Copyright (c) 2014 Philippe Pillot. All rights reserved.
//

#include <stdio.h>

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

#include "AliMpCDB.h"
#include "AliMpDEIterator.h"
#include "AliMpVSegmentation.h"
#include "AliMpSegmentation.h"
#include "AliMpVPadIterator.h"
#include "AliMpPad.h"

#include "PreClusterizerV2.h"

//#define VERBOSEMODE

//------------------------------------------------------------------
PreClusterizerV2::PreClusterizerV2(const char* ocdbPath, Int_t run)
: TObject(),
deIndices()
{
  
  // load mapping locally and create the local structure
  AliCDBManager* man = AliCDBManager::Instance();
  man->SetDefaultStorage(ocdbPath);
  man->SetRun(run);
  
  if (!AliMUONCDB::LoadMapping()) return;
  
  //AliMUONVStore* neighbours = AliMUONCalibrationData::CreateNeighbours(0);
  //CreateMapping(neighbours);
  CreateMapping(0x0);
  //delete neighbours;
  
  AliCodeTimer::Instance()->Print();
  
}

//------------------------------------------------------------------
PreClusterizerV2::~PreClusterizerV2()
{
  
  for (UChar_t iDE = 0; iDE < nDEs; ++iDE) delete[] mpDEs[iDE].pads;
  
  AliMpCDB::UnloadAll();
  
}

//------------------------------------------------------------------
void PreClusterizerV2::PreClusterize(const char *digitFileName, const char *clusterFileName, const char *method)
{
  /// Pre-clusterizer at DE level using "recursive" of "FIFO" method
  /// Each pad as an index in the pad array and knows the indices of its neighbours
  /// Compared to v1, arrays of pads and pre-clusters are reused and eventually extended but never deleted between events
  
  if (strcmp(method, "recursive") && strcmp(method, "FIFO") && strcmp(method, "both")) {
    printf("Choose either \"recursive\", \"FIFO\" or \"both\" method.\n");
    printf("If you choose \"both\" only \"recursive\" results will be stored.\n");
    return;
  }
  
  AliCodeTimer::Instance()->Reset();
  
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
  UShort_t nPreClusters[nDEs][2];
  std::vector<preCluster*> preClusters[nDEs][2];
  for (UChar_t iDE = 0; iDE < nDEs; ++iDE) {
    for (UChar_t iPlane = 0; iPlane < 2; ++iPlane) {
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
  
  // clean memory
  for (UChar_t iDE = 0; iDE < nDEs; ++iDE) {
    for (UChar_t iPlane = 0; iPlane < 2; ++iPlane) {
      for (UShort_t iCluster = 0; iCluster < static_cast<UShort_t>(preClusters[iDE][iPlane].size()); ++iCluster) {
        delete preClusters[iDE][iPlane][iCluster];
      }
    }
  }
  delete clusterStore;
  delete clusterFile;
  delete digitStore;
  
  AliCodeTimer::Instance()->Print();
  
}

//------------------------------------------------------------------
void PreClusterizerV2::CreateMapping(AliMUONVStore* neighbours)
{
  /// fill the mpDE and mpPad structures once for all
  
  AliCodeTimerAutoGeneral("",0);
  
  MemInfo_t memBefore, memAfter;
  ProcInfo_t procBefore, procAfter;
  gSystem->GetMemInfo(&memBefore);
  gSystem->GetProcInfo(&procBefore);
  
  // loop over DEs
  UChar_t iDE = 0;
  Int_t packValue2 = -1, manuId2 = -1, manuChannel2 = -1;
  UInt_t padId = 0;
  for (Int_t iCh = 0; iCh < AliMUONConstants::NTrackingCh(); ++iCh) {
    AliMpDEIterator deIt;
    deIt.First(iCh);
    while (!deIt.IsDone()) {
      
      Int_t deId = deIt.CurrentDEId();
      mpDE &de = mpDEs[iDE];
      de.id = deId;
      deIndices.Add(deId, iDE);
      
      const AliMpVSegmentation* seg[2] =
      { AliMpSegmentation::Instance()->GetMpSegmentation(deId,AliMp::kCath0),
        AliMpSegmentation::Instance()->GetMpSegmentation(deId,AliMp::kCath1)
      };
      
      de.nPads[seg[0]->PlaneType()] = seg[0]->NofPads();
      de.nPads[seg[1]->PlaneType()] = seg[1]->NofPads();
      
      de.pads = new mpPad[de.nPads[0]+de.nPads[1]];
      de.nOrderedPads[0] = 0;
      de.orderedPads[0].reserve((de.nPads[0]/10+de.nPads[1]/10)); // 10% occupancy
      de.nOrderedPads[1] = 0;
      de.orderedPads[1].reserve((de.nPads[0]/10+de.nPads[1]/10)); // 10% occupancy
      
      // loop over cathods
      for (Int_t iCath = 0; iCath < 2; ++iCath) {
        
        UChar_t iPlane = seg[iCath]->PlaneType();
        de.iPlanes[iCath] = iPlane;
        de.nFiredPads[iPlane] = 0;
        de.firedPads[iPlane].reserve(de.nPads[iPlane]/10); // 10% occupancy
        
        // 1st loop over pads to associate an index to their Id
        UShort_t iPad = (iPlane == 0) ? 0 : de.nPads[0];
        AliMpVPadIterator* padIt = seg[iCath]->CreateIterator();
        padIt->First();
        while (!padIt->IsDone()) {
          
          AliMpPad pad = padIt->CurrentItem();
          
          de.pads[iPad].nNeighbours = 0;
          de.pads[iPad].area[0][0] = pad.GetPositionX() - pad.GetDimensionX();
          de.pads[iPad].area[0][1] = pad.GetPositionX() + pad.GetDimensionX();
          de.pads[iPad].area[1][0] = pad.GetPositionY() - pad.GetDimensionY();
          de.pads[iPad].area[1][1] = pad.GetPositionY() + pad.GetDimensionY();
          de.pads[iPad].digit = 0x0;
          de.pads[iPad].useMe = kFALSE;
          
          padId = AliMUONVDigit::BuildUniqueID(deId, pad.GetManuId(), pad.GetManuChannel(), iCath);
          de.padIndices[iPlane].Add(padId, iPad);
          
          ++iPad;
          padIt->Next();
          
        }
        
        if (neighbours) {
          
          // 2nd loop over pads to fill the neighbours
          iPad = (iPlane == 0) ? 0 : de.nPads[0];
          padIt->First();
          while (!padIt->IsDone()) {
            
            AliMpPad pad = padIt->CurrentItem();
            AliMUONVCalibParam* calibParam = static_cast<AliMUONVCalibParam*>(neighbours->FindObject(deId, pad.GetManuId()));
            
            for (Int_t iNeighbour = 1; iNeighbour < calibParam->Dimension(); ++iNeighbour) {
              
              packValue2 = calibParam->ValueAsInt(pad.GetManuChannel(),iNeighbour);
              if (packValue2 <= 0) continue;
              
              calibParam->UnpackValue(packValue2, manuId2, manuChannel2);
              padId = AliMUONVDigit::BuildUniqueID(deId, manuId2, manuChannel2, iCath);
              
              de.pads[iPad].neighbours[de.pads[iPad].nNeighbours] = de.padIndices[iPlane].GetValue(padId);
              ++de.pads[iPad].nNeighbours;
              
            }
            
            ++iPad;
            padIt->Next();
            
          }
          
        } else FindNeighbours(de, iPlane);
        
        delete padIt;
        
      }
      
      ++iDE;
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
void PreClusterizerV2::FindNeighbours(mpDE &de, UChar_t iPlane)
{
  /// fill the mpPad neighbours structures of given DE/plane
  
  AliCodeTimerAutoGeneral("",0);
  
  UShort_t firstPad = (iPlane == 0) ? 0 : de.nPads[0];
  UShort_t lastPad = firstPad + de.nPads[iPlane];
  
  // loop over pads
  for (UShort_t iPad1 = firstPad; iPad1 < lastPad; ++iPad1) {
    
    // loop over next pads to find the neighbours
    for (UShort_t iPad2 = iPad1+1; iPad2 < lastPad; ++iPad2) {
      
      if (AreOverlapping(de.pads[iPad1].area, de.pads[iPad2].area, 1.e-4)) {
        
        de.pads[iPad1].neighbours[de.pads[iPad1].nNeighbours] = iPad2;
        ++de.pads[iPad1].nNeighbours;
        de.pads[iPad2].neighbours[de.pads[iPad2].nNeighbours] = iPad1;
        ++de.pads[iPad2].nNeighbours;
        
      }
      
    }
    
  }
  
}

//------------------------------------------------------------------
Int_t PreClusterizerV2::LoadDigits(AliMUONVDigitStore* digitStore)
{
  /// fill the mpDE structure with reconstructed digits
  
  AliCodeTimerAutoGeneral("",0);
  
  Int_t nDigits(0);
  UChar_t iPlane(0);
  UShort_t iPad(0);
  
  // loop over digits
  AliMUONVDigit* digit = 0x0;
  TIter nextDigit(digitStore->CreateIterator());
  while ((digit = static_cast<AliMUONVDigit*>(nextDigit()))) {
    
    if (digit->Charge() <= 0) continue;
    
    mpDE &de(mpDEs[deIndices.GetValue(digit->DetElemId())]);
    iPlane = de.iPlanes[digit->Cathode()];
    iPad = de.padIndices[iPlane].GetValue(digit->GetUniqueID());
    
    // attach the digit to this pad
    de.pads[iPad].digit = digit;
    de.pads[iPad].useMe = kTRUE;
    if (de.nFiredPads[iPlane] < static_cast<UShort_t>(de.firedPads[iPlane].size()))
      de.firedPads[iPlane][de.nFiredPads[iPlane]] = iPad;
    else de.firedPads[iPlane].push_back(iPad);
    ++de.nFiredPads[iPlane];
    ++nDigits;
    
  }
  
#ifdef VERBOSEMODE
  printf("\tTotal number of digits = %d\n", nDigits);
#endif
  
  return nDigits;
  
}

//------------------------------------------------------------------
void PreClusterizerV2::ResetPads()
{
  /// reset digit information in concerned pads
  
  AliCodeTimerAutoGeneral("",0);
  
  mpPad *pad(0x0);
  
  // loop over DEs
  for (UChar_t iDE = 0; iDE < nDEs; ++iDE) {
    
    mpDE &de(mpDEs[iDE]);
    
    // loop over planes
    for (UChar_t iPlane = 0; iPlane < 2; ++iPlane) {
      
      // loop over fired pads
      for (UShort_t iFiredPad = 0; iFiredPad < de.nFiredPads[iPlane]; ++iFiredPad) {
        
        pad = &de.pads[de.firedPads[iPlane][iFiredPad]];
        pad->digit = 0x0;
        pad->useMe = kFALSE;
        
      }
      
      // clear number of fired pads
      de.nFiredPads[iPlane] = 0;
      
    }
    
    // clear ordered number of fired pads
    de.nOrderedPads[0] = 0;
    de.nOrderedPads[1] = 0;
    
  }
  
}

//------------------------------------------------------------------
Int_t PreClusterizerV2::PreClusterizeFIFO(std::vector<preCluster*> preClusters[nDEs][2], UShort_t nPreClusters[nDEs][2])
{
  /// preclusterize both planes of every DE using FIFO algorithm
  
  AliCodeTimerAutoGeneral("",0);
  
  Int_t nPreclustersTot(0);
  preCluster *cl(0x0);
  UShort_t iPad(0);
  
  // loop over DEs
  for (UChar_t iDE = 0; iDE < nDEs; ++iDE) {
    
    mpDE &de(mpDEs[iDE]);
    
    // loop over planes
    for (UChar_t iPlane = 0; iPlane < 2; ++iPlane) {
      
      nPreClusters[iDE][iPlane] = 0;
      
      // loop over fired pads
      for (UShort_t iFiredPad = 0; iFiredPad < de.nFiredPads[iPlane]; ++iFiredPad) {
        
        iPad = de.firedPads[iPlane][iFiredPad];
        
        if (de.pads[iPad].useMe) {
          
          // create the precluster if needed
          if (nPreClusters[iDE][iPlane] >= static_cast<UShort_t>(preClusters[iDE][iPlane].size()))
            preClusters[iDE][iPlane].push_back(new preCluster);
          
          // get the precluster
          cl = preClusters[iDE][iPlane][nPreClusters[iDE][iPlane]];
          ++nPreClusters[iDE][iPlane];
          ++nPreclustersTot;
          
          // reset its content
          if (de.nOrderedPads[0] < static_cast<UShort_t>(de.orderedPads[0].size()))
            de.orderedPads[0][de.nOrderedPads[0]] = iPad;
          else de.orderedPads[0].push_back(iPad);
          cl->firstPad = de.nOrderedPads[0];
          cl->lastPad = de.nOrderedPads[0];
          ++de.nOrderedPads[0];
          mpPad &pad(de.pads[iPad]);
          cl->area[0][0] = pad.area[0][0];
          cl->area[0][1] = pad.area[0][1];
          cl->area[1][0] = pad.area[1][0];
          cl->area[1][1] = pad.area[1][1];
          cl->useMe = kTRUE;
          cl->storeMe = kFALSE;
          
          pad.useMe = kFALSE;
          
          // loop over all pads of the precluster
          for (UShort_t iOrderPad = cl->firstPad; iOrderPad <= cl->lastPad; ++iOrderPad) {
            
            mpPad &pad1(de.pads[de.orderedPads[0][iOrderPad]]);
            
            // loop over their neighbours
            for (UShort_t iNeighbour = 0; iNeighbour < pad1.nNeighbours; ++iNeighbour) {
              
              iPad = pad1.neighbours[iNeighbour];
              
              if (de.pads[iPad].useMe) {
                
                // add the pad to the precluster
                if (de.nOrderedPads[0] < static_cast<UShort_t>(de.orderedPads[0].size()))
                  de.orderedPads[0][de.nOrderedPads[0]] = iPad;
                else de.orderedPads[0].push_back(iPad);
                cl->lastPad = de.nOrderedPads[0];
                ++de.nOrderedPads[0];
                mpPad &pad2(de.pads[iPad]);
                if (pad2.area[0][0] < cl->area[0][0]) cl->area[0][0] = pad2.area[0][0];
                if (pad2.area[0][1] > cl->area[0][1]) cl->area[0][1] = pad2.area[0][1];
                if (pad2.area[1][0] < cl->area[1][0]) cl->area[1][0] = pad2.area[1][0];
                if (pad2.area[1][1] > cl->area[1][1]) cl->area[1][1] = pad2.area[1][1];
                
                pad2.useMe = kFALSE;
                
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
Int_t PreClusterizerV2::PreClusterizeRecursive(std::vector<preCluster*> preClusters[nDEs][2], UShort_t nPreClusters[nDEs][2])
{
  /// preclusterize both planes of every DE using recursive algorithm
  
  AliCodeTimerAutoGeneral("",0);
  
  Int_t nPreclustersTot = 0;
  preCluster *cl(0x0);
  UShort_t iPad(0);
  
  // loop over DEs
  for (UChar_t iDE = 0; iDE < nDEs; ++iDE) {
    
    mpDE &de(mpDEs[iDE]);
    
    // loop over planes
    for (UChar_t iPlane = 0; iPlane < 2; ++iPlane) {
      
      nPreClusters[iDE][iPlane] = 0;
      
      // loop over fired pads
      for (UShort_t iFiredPad = 0; iFiredPad < de.nFiredPads[iPlane]; ++iFiredPad) {
        
        iPad = de.firedPads[iPlane][iFiredPad];
        
        if (de.pads[iPad].useMe) {
          
          // create the precluster if needed
          if (nPreClusters[iDE][iPlane] >= static_cast<UShort_t>(preClusters[iDE][iPlane].size()))
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
          cl->firstPad = de.nOrderedPads[0];
          AddPad(de, iPad, *cl);
          
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
void PreClusterizerV2::AddPad(mpDE &de, UShort_t iPad, preCluster &cl)
{
  /// add the given mpPad and its fired neighbours (recursive method)
  
  // add the given pad
  mpPad &pad(de.pads[iPad]);
  if (de.nOrderedPads[0] < static_cast<UShort_t>(de.orderedPads[0].size()))
    de.orderedPads[0][de.nOrderedPads[0]] = iPad;
  else de.orderedPads[0].push_back(iPad);
  cl.lastPad = de.nOrderedPads[0];
  ++de.nOrderedPads[0];
  if (pad.area[0][0] < cl.area[0][0]) cl.area[0][0] = pad.area[0][0];
  if (pad.area[0][1] > cl.area[0][1]) cl.area[0][1] = pad.area[0][1];
  if (pad.area[1][0] < cl.area[1][0]) cl.area[1][0] = pad.area[1][0];
  if (pad.area[1][1] > cl.area[1][1]) cl.area[1][1] = pad.area[1][1];
  
  pad.useMe = kFALSE;
  
  // loop over its neighbours
  for (UShort_t iNeighbour = 0; iNeighbour < pad.nNeighbours; ++iNeighbour) {
    
    if (de.pads[pad.neighbours[iNeighbour]].useMe) {
      
      // add the pad to the precluster
      AddPad(de, pad.neighbours[iNeighbour], cl);
      
    }
    
  }
  
}

//------------------------------------------------------------------
Int_t PreClusterizerV2::MergePreClusters(std::vector<preCluster*> preClusters[nDEs][2], UShort_t nPreClusters[nDEs][2])
{
  /// merge overlapping preclusters on every DE
  
  AliCodeTimerAutoGeneral("",0);
  
  Int_t nPreClustersMergedTot = 0;
  preCluster *cl(0x0);
  
  // loop over DEs
  for (UChar_t iDE = 0; iDE < nDEs; ++iDE) {
    
    mpDE &de(mpDEs[iDE]);
    
    // loop over preclusters of one plane
    for (UShort_t iCluster = 0; iCluster < nPreClusters[iDE][0]; ++iCluster) {
      
      if (!preClusters[iDE][0][iCluster]->useMe) continue;
      
      cl = preClusters[iDE][0][iCluster];
      cl->useMe = kFALSE;
      
      // look for overlapping preclusters in the other plane
      preCluster *mergedCl(0x0);
      MergePreClusters(*cl, preClusters[iDE], nPreClusters[iDE], de, 1, mergedCl);
      
      // add the current one
      if (!mergedCl) mergedCl = UsePreClusters(cl, de);
      else MergePreClusters(*mergedCl, *cl, de);
      
      ++nPreClustersMergedTot;
      
    }
    
    // loop over preclusters of the other plane
    for (UShort_t iCluster = 0; iCluster < nPreClusters[iDE][1]; ++iCluster) {
      
      if (!preClusters[iDE][1][iCluster]->useMe) continue;
      
      // all remaining preclusters have to be stored
      UsePreClusters(preClusters[iDE][1][iCluster], de);
      
      ++nPreClustersMergedTot;
      
    }
    
  }
  
#ifdef VERBOSEMODE
  printf("\tTotal number of preclusters after merging = %d\n", nPreClustersMergedTot);
#endif
  
  return nPreClustersMergedTot;
  
}

//------------------------------------------------------------------
void PreClusterizerV2::MergePreClusters(preCluster &cl, std::vector<preCluster*> preClusters[2], UShort_t nPreClusters[2],
                                        mpDE &de, UChar_t iPlane, preCluster *&mergedCl)
{
  /// merge preclusters on the given plane overlapping with the given one (recursive method)
  
  preCluster *cl2(0x0);
  
  // loop over preclusters in the given plane
  for (UShort_t iCluster = 0; iCluster < nPreClusters[iPlane]; ++iCluster) {
    
    if (!preClusters[iPlane][iCluster]->useMe) continue;
    
    cl2 = preClusters[iPlane][iCluster];
    if (AreOverlapping(cl.area, cl2->area, -1.e-4) && AreOverlapping(cl, *cl2, de, -1.e-4)) {
      
      cl2->useMe = kFALSE;
      
      // look for new overlapping preclusters in the other plane
      MergePreClusters(*cl2, preClusters, nPreClusters, de, (iPlane+1)%2, mergedCl);
      
      // store overlapping preclusters and merge them
      if (!mergedCl) mergedCl = UsePreClusters(cl2, de);
      else MergePreClusters(*mergedCl, *cl2, de);
      
    }
    
  }
  
}

//------------------------------------------------------------------
preCluster* PreClusterizerV2::UsePreClusters(preCluster *cl, mpDE &de)
{
  /// use this precluster as a new merged precluster
  
  UShort_t firstPad = de.nOrderedPads[1];
  
  // move the fired pads
  for (UShort_t iOrderPad = cl->firstPad; iOrderPad <= cl->lastPad; ++iOrderPad) {
    
    if (de.nOrderedPads[1] < static_cast<UShort_t>(de.orderedPads[1].size()))
      de.orderedPads[1][de.nOrderedPads[1]] = de.orderedPads[0][iOrderPad];
    else de.orderedPads[1].push_back(de.orderedPads[0][iOrderPad]);
    
    ++de.nOrderedPads[1];
    
  }
  
  cl->firstPad = firstPad;
  cl->lastPad = de.nOrderedPads[1] - 1;
  
  cl->storeMe = kTRUE;
  
  return cl;
  
}

//------------------------------------------------------------------
void PreClusterizerV2::MergePreClusters(preCluster &cl1, preCluster &cl2, mpDE &de)
{
  /// merge precluster2 into precluster1
  
  // move the fired pads
  for (UShort_t iOrderPad = cl2.firstPad; iOrderPad <= cl2.lastPad; ++iOrderPad) {
    
    if (de.nOrderedPads[1] < static_cast<UShort_t>(de.orderedPads[1].size()))
      de.orderedPads[1][de.nOrderedPads[1]] = de.orderedPads[0][iOrderPad];
    else de.orderedPads[1].push_back(de.orderedPads[0][iOrderPad]);
    
    ++de.nOrderedPads[1];
    
  }
  
  cl1.lastPad = de.nOrderedPads[1] - 1;
  
}

//___________________________________________________________________________
Bool_t PreClusterizerV2::AreOverlapping(preCluster &cl1, preCluster &cl2, mpDE &de, Float_t precision)
{
  /// check if the two preclusters overlap
  /// precision in cm: positive = increase pad size / negative = decrease pad size
  
  // loop over all pads of the precluster1
  for (UShort_t iOrderPad1 = cl1.firstPad; iOrderPad1 <= cl1.lastPad; ++iOrderPad1) {
    
    // loop over all pads of the precluster2
    for (UShort_t iOrderPad2 = cl2.firstPad; iOrderPad2 <= cl2.lastPad; ++iOrderPad2) {
      
      if (AreOverlapping(de.pads[de.orderedPads[0][iOrderPad1]].area,
                         de.pads[de.orderedPads[0][iOrderPad2]].area, precision)) return kTRUE;
      
    }
    
  }
  
  return kFALSE;
  
}

//------------------------------------------------------------------
Bool_t PreClusterizerV2::AreOverlapping(Float_t area1[2][2], Float_t area2[2][2], Float_t precision)
{
  /// check if the two areas overlap
  /// precision in cm: positive = increase pad size / negative = decrease pad size
  
  if (area1[0][0] - area2[0][1] > precision) return kFALSE;
  if (area2[0][0] - area1[0][1] > precision) return kFALSE;
  if (area1[1][0] - area2[1][1] > precision) return kFALSE;
  if (area2[1][0] - area1[1][1] > precision) return kFALSE;
  
  return kTRUE;
  
}

//------------------------------------------------------------------
void PreClusterizerV2::StorePreClusters(std::vector<preCluster*> preClusters[nDEs][2], UShort_t nPreClusters[nDEs][2],
                                        AliMUONVClusterStore *clusterStore)
{
  /// store the preclusters into standard cluster store
  
  AliCodeTimerAutoGeneral("",0);
  
  clusterStore->Clear();
  Int_t clusterIndex(0);
  preCluster *cl(0x0);
  AliMUONVCluster *cluster(0x0);
  
  // loop over DEs
  for (UChar_t iDE = 0; iDE < nDEs; ++iDE) {
    
    mpDE &de(mpDEs[iDE]);
    
    // temporary array of IDs
    UInt_t *digitsId = new UInt_t[de.nOrderedPads[1]];
    
    // loop over planes
    for (UChar_t iPlane = 0; iPlane < 2; ++iPlane) {
      
      // loop over preclusters
      for (UShort_t iCluster = 0; iCluster < nPreClusters[iDE][iPlane]; ++iCluster) {
        
        if (!preClusters[iDE][iPlane][iCluster]->storeMe) continue;
        
        cl = preClusters[iDE][iPlane][iCluster];
        cluster = clusterStore->Add(de.id/100-1, de.id, clusterIndex);
        ++clusterIndex;
        
        // loop over pads
        Int_t nDigits(0);
        for (UShort_t iOrderPad = cl->firstPad; iOrderPad <= cl->lastPad; ++iOrderPad) {
          
          digitsId[nDigits] = de.pads[de.orderedPads[1][iOrderPad]].digit->GetUniqueID();
          ++nDigits;
          
        }
        
        cluster->SetDigitsId(nDigits, digitsId);
        
      }
      
    }
    
    // clean memory
    delete[] digitsId;
    
  }
  
}

