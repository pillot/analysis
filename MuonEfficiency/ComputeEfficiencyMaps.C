/*
 *  ComputeEfficiencyMaps.C
 *  aliroot_dev
 *
 *  Created by philippe pillot on 17/06/13.
 *  Copyright 2013 Subatech. All rights reserved.
 *
 */


#include <Riostream.h>

#include <TFile.h>
#include <TString.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TIterator.h>
#include <TList.h>
#include <TSystem.h>

#include "AliCDBManager.h"
#include "AliCounterCollection.h"

#include "AliMUONCDB.h"
#include "AliMUONRejectList.h"
#include "AliMUONTrackerData.h"

#include "mapping/AliMpSegmentation.h"
#include "mapping/AliMpDDLStore.h"
#include "mapping/AliMpDEIterator.h"
#include "mapping/AliMpBusPatch.h"
#include "mapping/AliMpManuIterator.h"
#include "mapping/AliMpConstants.h"


void ResetEfficiencies(Float_t chamberEff[10][2], AliMUONRejectList effMaps[2]);
void ComputeDEEfficiency(AliCounterCollection &nClusters, AliMUONRejectList effMaps[2]);
void ComputeBusPatchEfficiency(AliCounterCollection &nClusters, AliMUONRejectList effMaps[2]);
void ComputeManuEfficiency(AliCounterCollection &nClusters, AliMUONRejectList effMaps[2]);
void ComputeChannelEfficiency(AliCounterCollection &nClusters, AliMUONRejectList effMaps[2]);
void ComputeEfficiency(Float_t nExpected, Float_t nAccepted, Float_t efficiency[2]);


//---------------------------------------------------------------------------
void ComputeEfficiencyMaps()
{
  /// compute the efficiency maps and associated errors per Manu, Bus Patch and Detection Element
  /*
   aliroot -l
   .include $ALICE_ROOT/MUON
   .x $WORK/Macros/MuonEfficiency/ComputeEfficiencyMaps.C+
  */
  
  // load mapping locally if not already done
  AliCDBManager* man = AliCDBManager::Instance();
  if (!man->IsDefaultStorageSet()) man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  if (man->GetRun() < 0) man->SetRun(0);
  if (!AliMpSegmentation::Instance(kFALSE) || !AliMpDDLStore::Instance(kFALSE)) {
    if (!AliMUONCDB::LoadMapping()) return;
  }
  
  // get input counters
  TFile *inFile = new TFile("AnalysisResults.root", "READ");
  if (!inFile || !inFile->IsOpen()) {
    printf("cannot open the file AnalysisResults.root\n");
    return;
  }
  AliCounterCollection *nClusters = static_cast<AliCounterCollection*>(inFile->FindObjectAny("Clusters"));
  if (!nClusters) {
    printf("cannot find the counter collection \"Clusters\"\n");
    return;
  }
  
  // compute the efficiency and associated error at every levels
  Float_t chamberEff[10][2];
  AliMUONRejectList effMaps[2];
  ResetEfficiencies(chamberEff, effMaps);
  ComputeDEEfficiency(*nClusters, effMaps);
  ComputeBusPatchEfficiency(*nClusters, effMaps);
  ComputeManuEfficiency(*nClusters, effMaps);
  ComputeChannelEfficiency(*nClusters, effMaps);
  
  // close input file
  inFile->Close();
  
  // put the efficiency maps in a tracker data for display
  AliMUONTrackerData effMapsDisp("efficiencies", "efficiency maps", effMaps[0]);
  effMapsDisp.SetDimensionName(0,"efficiency");
  AliMUONTrackerData errMapsDisp("uncertainties", "uncertainty maps", effMaps[1]);
  errMapsDisp.SetDimensionName(0,"uncertainty");
  
  // save efficiency maps
  TList eff;
  eff.SetName("maps");
  eff.SetOwner(kFALSE);
  eff.AddLast(&(effMaps[0]));
  eff.AddLast(&(effMaps[1]));
  TFile *outFile = TFile::Open("EfficiencyMaps.root", "RECREATE");
  if (outFile && outFile->IsOpen()) {
    eff.Write(0x0, TObject::kOverwrite);
    effMapsDisp.Write(0x0, TObject::kOverwrite);
    errMapsDisp.Write(0x0, TObject::kOverwrite);
    outFile->Close();
  }
  
  // display efficiency maps
  gSystem->Exec("mchview --use EfficiencyMaps.root");
  
}

//---------------------------------------------------------------------------
void ResetEfficiencies(Float_t chamberEff[10][2], AliMUONRejectList effMaps[2])
{
  /// set the efficiency and associated error to default values at every levels
  
  // at chamber level
  for (Int_t iChamber = 0; iChamber < 10; iChamber++) {
    
    chamberEff[iChamber][0] = -1.;
    chamberEff[iChamber][1] = 2.;
    
    // at Detection Element level
    AliMpDEIterator deIt;
    deIt.First(iChamber);
    while (!deIt.IsDone()) {
      Int_t deId = deIt.CurrentDEId();
      effMaps[0].SetDetectionElementProbability(deId, -1.);
      effMaps[1].SetDetectionElementProbability(deId, 2.);
      deIt.Next();
    }
    
  }
  
  // at Bus Patch level
  TIterator* bpIt = AliMpDDLStore::Instance()->CreateBusPatchIterator();
  AliMpBusPatch* bp = 0x0;
  while ((bp = static_cast<AliMpBusPatch*>(bpIt->Next()))) {
    Int_t bpId = bp->GetId();
    effMaps[0].SetBusPatchProbability(bpId, -1.);
    effMaps[1].SetBusPatchProbability(bpId, 2.);
  }
  delete bpIt;
  
  // at Manu level
  AliMpManuIterator manuIt;
  Int_t deId, manuId;
  while (manuIt.Next(deId, manuId)) {
    
    effMaps[0].SetManuProbability(deId, manuId, -1.);
    effMaps[1].SetManuProbability(deId, manuId, 2.);
    
    // at channel level
    for (Int_t iChannel = 0; iChannel < AliMpConstants::ManuNofChannels(); iChannel++) {
      effMaps[0].SetChannelProbability(deId, manuId, iChannel, -1.);
      effMaps[1].SetChannelProbability(deId, manuId, iChannel, 2.);
    }
    
  }
  
}

//---------------------------------------------------------------------------
void ComputeDEEfficiency(AliCounterCollection &nClusters, AliMUONRejectList effMaps[2])
{
  /// compute the efficiency and associated error per Detection Elements
  
  // get numbers of clusters per Detection Elements
  TH1D *ExpectedDE = nClusters.Get("DE", "Cluster:Expected");
  ExpectedDE->SetName("expected");
  TH1D *AcceptedDE = nClusters.Get("DE", "Cluster:Accepted");
  AcceptedDE->SetName("accepted");
  Int_t nDEs = ExpectedDE->GetNbinsX();
  
  // loop over fired Detection Elements
  for (Int_t iDE = 1; iDE <= nDEs; iDE++) {
    
    // get Detection Element Id
    TString deKey = ExpectedDE->GetXaxis()->GetBinLabel(iDE);
    Int_t deId = deKey.Atoi();
    
    // get numbers of clusters in this Detection Element (! clusters are added to both cathodes)
    Float_t nExpected = 0.5 * static_cast<Float_t>(ExpectedDE->GetBinContent(iDE));
    if (nExpected < 1.) continue;
    Float_t nAccepted = 0.5 * static_cast<Float_t>(AcceptedDE->GetBinContent(iDE));
    
    // compute efficiency
    Float_t efficiency[2];
    ComputeEfficiency(nExpected, nAccepted, efficiency);
    
    // register the results
    effMaps[0].SetDetectionElementProbability(deId, efficiency[0]);
    effMaps[1].SetDetectionElementProbability(deId, efficiency[1]);
    
  }
  
  // clean memory
  delete ExpectedDE;
  delete AcceptedDE;
  
}

//---------------------------------------------------------------------------
void ComputeBusPatchEfficiency(AliCounterCollection &nClusters, AliMUONRejectList effMaps[2])
{
  /// compute the efficiency and associated error per Bus Patch
  
  // get numbers of clusters per Bus Patchs
  TH1D *ExpectedBusPatch = nClusters.Get("BusPatch", "Cluster:Expected");
  ExpectedBusPatch->SetName("expected");
  TH1D *AcceptedBusPatch = nClusters.Get("BusPatch", "Cluster:Accepted");
  AcceptedBusPatch->SetName("accepted");
  Int_t nBusPatchs = ExpectedBusPatch->GetNbinsX();
  
  // loop over fired Bus Patchs
  for (Int_t iBusPatch = 1; iBusPatch <= nBusPatchs; iBusPatch++) {
    
    // get Bus Patch Id
    TString bpKey = ExpectedBusPatch->GetXaxis()->GetBinLabel(iBusPatch);
    Int_t bpId = bpKey.Atoi();
    
    // get numbers of clusters in this Bus Patch
    Float_t nExpected = static_cast<Float_t>(ExpectedBusPatch->GetBinContent(iBusPatch));
    if (nExpected < 1.) continue;
    Float_t nAccepted = static_cast<Float_t>(AcceptedBusPatch->GetBinContent(iBusPatch));
    
    // compute efficiency
    Float_t efficiency[2];
    ComputeEfficiency(nExpected, nAccepted, efficiency);
    
    // register the results
    effMaps[0].SetBusPatchProbability(bpId, efficiency[0]);
    effMaps[1].SetBusPatchProbability(bpId, efficiency[1]);
    
  }
  
  // clean memory
  delete ExpectedBusPatch;
  delete AcceptedBusPatch;
  
}

//---------------------------------------------------------------------------
void ComputeManuEfficiency(AliCounterCollection &nClusters, AliMUONRejectList effMaps[2])
{
  /// compute the efficiency and associated error per Manu
  
  // get numbers of clusters for each DE/manu combination
  TH2D *ExpectedManuVsDE = nClusters.Get("Manu", "DE", "Cluster:Expected");
  ExpectedManuVsDE->SetName("expected");
  TH2D *AcceptedManuVsDE = nClusters.Get("Manu", "DE", "Cluster:Accepted");
  AcceptedManuVsDE->SetName("accepted");
  Int_t nDEs = ExpectedManuVsDE->GetNbinsX();
  Int_t nManus = ExpectedManuVsDE->GetNbinsY();
  
  // loop over fired Detection Elements
  for (Int_t iDE = 1; iDE <= nDEs; iDE++) {
    
    // get Detection Element Id
    TString deKey = ExpectedManuVsDE->GetXaxis()->GetBinLabel(iDE);
    Int_t deId = deKey.Atoi();
    
    // loop over fired Manus
    for (Int_t iManu = 1; iManu <= nManus; iManu++) {
      
      // get Manu Id
      TString manuKey = ExpectedManuVsDE->GetYaxis()->GetBinLabel(iManu);
      Int_t manuId = manuKey.Atoi();
      
      // get the numbers of clusters in this Manu
      Float_t nExpected = static_cast<Float_t>(ExpectedManuVsDE->GetBinContent(iDE, iManu));
      if (nExpected < 1.) continue;
      Float_t nAccepted = static_cast<Float_t>(AcceptedManuVsDE->GetBinContent(iDE, iManu));
      
      // compute efficiency
      Float_t efficiency[2];
      ComputeEfficiency(nExpected, nAccepted, efficiency);
      
      // register the results
      effMaps[0].SetManuProbability(deId, manuId, efficiency[0]);
      effMaps[1].SetManuProbability(deId, manuId, efficiency[1]);
      
    }
    
  }
  
  // clean memory
  delete ExpectedManuVsDE;
  delete AcceptedManuVsDE;
  
}

//---------------------------------------------------------------------------
void ComputeChannelEfficiency(AliCounterCollection &nClusters, AliMUONRejectList effMaps[2])
{
  /// compute the efficiency and associated error per channel
  
  // list of fired DEs
  TObjArray *deKeys = nClusters.GetKeyWords("DE").Tokenize(",");
  Int_t nDEs = deKeys->GetEntriesFast();
  
  // loop over fired DEs
  for (Int_t iDE = 0; iDE < nDEs; iDE++) {
    
    // get DE Id
    TString deKey = static_cast<TObjString*>(deKeys->UncheckedAt(iDE))->GetString();
    Int_t deId = deKey.Atoi();
    
    // get numbers of clusters in this DE for each manu/channel combination
    TH2D *ExpectedChannelVsManu = nClusters.Get("channel", "Manu", Form("Cluster:Expected/DE:%s", deKey.Data()));
    ExpectedChannelVsManu->SetName("expected");
    TH2D *AcceptedChannelVsManu = nClusters.Get("channel", "Manu", Form("Cluster:Accepted/DE:%s", deKey.Data()));
    AcceptedChannelVsManu->SetName("accepted");
    Int_t nManus = ExpectedChannelVsManu->GetNbinsX();
    Int_t nChannels = ExpectedChannelVsManu->GetNbinsY();
    
    // loop over fired Manus
    for (Int_t iManu = 1; iManu <= nManus; iManu++) {
      
      // get Manu Id
      TString manuKey = ExpectedChannelVsManu->GetXaxis()->GetBinLabel(iManu);
      Int_t manuId = manuKey.Atoi();
      
      // loop over fired channels
      for (Int_t iChannel = 1; iChannel <= nChannels; iChannel++) {
	
	// get channel Id
	TString channelKey = ExpectedChannelVsManu->GetYaxis()->GetBinLabel(iChannel);
	Int_t channelId = channelKey.Atoi();
	
	// get the numbers of clusters in this channel
	Float_t nExpected = static_cast<Float_t>(ExpectedChannelVsManu->GetBinContent(iManu, iChannel));
	if (nExpected < 1.) continue;
	Float_t nAccepted = static_cast<Float_t>(AcceptedChannelVsManu->GetBinContent(iManu, iChannel));
	
	// compute efficiency
	Float_t efficiency[2];
	ComputeEfficiency(nExpected, nAccepted, efficiency);
	
	// register the results
	effMaps[0].SetChannelProbability(deId, manuId, channelId, efficiency[0]);
	effMaps[1].SetChannelProbability(deId, manuId, channelId, efficiency[1]);
	
      }
      
    }
    
    // clean memory
    delete ExpectedChannelVsManu;
    delete AcceptedChannelVsManu;
    
  }
  
  // clean memory
  delete deKeys;
  
}

//---------------------------------------------------------------------------
void ComputeEfficiency(Float_t nExpected, Float_t nAccepted, Float_t efficiency[2])
{
  /// compute the efficiency and associated error
  
  if (nExpected > 0.) {
    
    efficiency[0] = nAccepted / nExpected;
    efficiency[1] = TMath::Max(1./nExpected, TMath::Sqrt(efficiency[0]*TMath::Abs(1.-efficiency[0])/nExpected));
  
  } else {
    
    efficiency[0] = -1.;
    efficiency[1] = 2.;
    
  }
  
}

