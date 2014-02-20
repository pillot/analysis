/*
 *  checkEfficiency.C
 *  aliroot
 *
 *  Created by Philippe Pillot on 09/12/08.
 *  Copyright 2008 Subatech. All rights reserved.
 *
 */


#if !defined(__CINT__) || defined(__MAKECINT__)
// ROOT includes
#include "TTree.h"
#include "TBranch.h"
#include "TObjArray.h"
#include "TLorentzVector.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TParticle.h"
#include "TTree.h"
#include <Riostream.h>
#include <TROOT.h>

// STEER includes
#include "AliLog.h"
#include "AliRunLoader.h"
#include "AliLoader.h"

// MUON includes
#include "AliMUONConstants.h"
#include "AliMUONTrack.h"
#include "AliMUONVTrackStore.h"
#include "AliMUONTrackParam.h"
#include "AliMUONRecoCheck.h"
#include "AliMUONVCluster.h"
#include "AliMUONVClusterStore.h"
#endif

/// \ingroup macros
/// \file checkEfficiency.C
/// \brief Macro checkEfficiency.C
///
/// \author Ph. Pillot, Subatech, December 2008

Double_t sigmaCut = 5;

void checkEfficiency(Bool_t useLabel = kFALSE)
{
  //Reset ROOT and connect tree file
  gROOT->Reset();
  
  // File for histograms and histogram booking
  TFile *histoFile = new TFile("checkEfficiency.root", "RECREATE");
  TH1F* hESDnTrackRefsTotPerCh = new TH1F("hESDnTrackRefsTotPerCh", "number of trackRefs per chamber;chamber ID", 10, 0., 10.);
  hESDnTrackRefsTotPerCh->SetFillColor(kRed);
  TH1F* hESDnMatchedTrackRefsTotPerCh = new TH1F("hESDnMatchedTrackRefsTotPerCh", "number of trackRefs matched with cluster per chamber;chamber ID", 10, 0., 10.);
  hESDnMatchedTrackRefsTotPerCh->SetFillColor(kRed);
  TH1F* hESDFracMatchedTrackRefsTotPerCh = new TH1F("hESDFracMatchedTrackRefsTotPerCh", "fraction of trackRefs matched with cluster per chamber;chamber ID", 10, 0., 10.);
  hESDFracMatchedTrackRefsTotPerCh->SetFillColor(kRed);
  
  TH1F* hESDnClustersTotPerCh = new TH1F("hESDnClustersTotPerCh", "number of clusters per chamber;chamber ID", 10, 0., 10.);
  hESDnClustersTotPerCh->SetFillColor(kRed);
  TH1F* hESDnMatchedClustersTotPerCh = new TH1F("hESDnMatchedClustersTotPerCh", "number of clusters matched with trackRef per chamber;chamber ID", 10, 0., 10.);
  hESDnMatchedClustersTotPerCh->SetFillColor(kRed);
  TH1F* hESDFracMatchedClustersTotPerCh = new TH1F("hESDFracMatchedClustersTotPerCh", "fraction of clusters matched with trackRef per chamber;chamber ID", 10, 0., 10.);
  hESDFracMatchedClustersTotPerCh->SetFillColor(kRed);
  
  TH1F* hESDnClustersTotOvernTrackRefTotPerCh = new TH1F("hESDnClustersTotOvernTrackRefTotPerCh", "number of clusters over number of trackRef per chamber;chamber ID", 10, 0., 10.);
  hESDnClustersTotOvernTrackRefTotPerCh->SetFillColor(kRed);
  
  TH2F* hESDunmatchedTrackRefsMap[10];
  TH2F* hESDunmatchedClustersMap[10];
  for (Int_t i = 0; i < 10; i++) {
    Float_t rMax = AliMUONConstants::Rmax(i/2);
    hESDunmatchedClustersMap[i] = new TH2F(Form("hESDunmatchedClustersMap%d",i+1), Form("unmatched cluster position distribution in chamber %d",i+1),2000, -rMax, rMax, 2000, -rMax, rMax);
    hESDunmatchedTrackRefsMap[i] = new TH2F(Form("hESDunmatchedTrackRefsMap%d",i+1), Form("unmatched trackRef position distribution in chamber %d",i+1),2000, -rMax, rMax, 2000, -rMax, rMax);
  }
  
  AliRunLoader * rl = AliRunLoader::Open("galice.root","MUONLoader");
  AliLoader* MUONLoader = rl->GetDetectorLoader("MUON");
  MUONLoader->LoadRecPoints("READ");   
  
//  AliMUONRecoCheck rc("AliESDs.root", "./generated/");
  AliMUONRecoCheck rc("AliESDs.root", "./");
  
  Int_t nGoodLabelC = 0, nGoodLabelT = 0, nBadLabelC = 0, nBadLabelT = 0, nNoLabelC = 0, nNoLabelT = 0;
  Int_t nEventUsed = 0;
  
  // Loop over events
  Int_t nevents = rc.NumberOfEvents();
  for (Int_t iEvent = 0; iEvent < nevents; iEvent++) {
    
    if (!(rl->GetEvent(iEvent) == 0)) continue;
    TTree* treeR = MUONLoader->TreeR();
    AliMUONVClusterStore* clusterStore = AliMUONVClusterStore::Create(*treeR);
    if ( clusterStore != 0x0 ) {
      clusterStore->Clear();
      clusterStore->Connect(*treeR);
      treeR->GetEvent(0);
    }
    
    AliMUONVTrackStore* TrackRefStore = rc.TrackRefs(iEvent);
    
    // remove shower events
    //if (TrackRefStore->GetSize() > 2) {
    //  delete clusterStore;
    //  continue;
    //}
    
    nEventUsed++;
    
    // ------------------- TrackRef -------------------
    TIter nextTrackRef(TrackRefStore->CreateIterator());
    AliMUONTrack* simulatedTrack;
    while ( ( simulatedTrack = static_cast<AliMUONTrack*>(nextTrackRef()) ) ) {
      
      Int_t label = simulatedTrack->GetUniqueID();
      
      AliMUONTrackParam* trackParam = static_cast<AliMUONTrackParam*>(simulatedTrack->GetTrackParamAtCluster()->First());
      while (trackParam) {
	
	AliMUONVCluster *trackRef = trackParam->GetClusterPtr();
	Int_t chId = trackRef->GetChamberId();
	Int_t deId = trackRef->GetDetElemId();
	Double_t bestChi2 =  2*sigmaCut*sigmaCut + 1.;
	Int_t bestLabel = -1;
	
	hESDnTrackRefsTotPerCh->Fill(chId);
	
	TIter nextCluster(clusterStore->CreateChamberIterator(chId,chId));
	AliMUONVCluster *cluster;
	Bool_t matched = kFALSE;
	while ( ( cluster = static_cast<AliMUONVCluster*>(nextCluster()) ) ) {
	  
	  if (cluster->GetDetElemId() != deId) continue;
	  
	  if (useLabel) { // with label
	    
	    if (cluster->GetMCLabel() == label) matched = kTRUE;
	    
	  } else { // without label
	    
	    Double_t deltaX = cluster->GetX() - trackRef->GetX();
	    Double_t deltaY = cluster->GetY() - trackRef->GetY();
	    Double_t chi2 = deltaX * deltaX / cluster->GetErrX2() + deltaY * deltaY / cluster->GetErrY2();
	    if (chi2 <= 2*sigmaCut*sigmaCut) {
	      matched = kTRUE;
	      if (chi2<bestChi2) {
		bestChi2 = chi2;
		bestLabel = cluster->GetMCLabel();
	      }
	    }
	    
	  }
	  
	  if (useLabel && matched) break;
	  
	}
	
	if (matched) {
	  hESDnMatchedTrackRefsTotPerCh->Fill(chId);
	  if (bestLabel == label) nGoodLabelT++;
	  else if (bestLabel < 0) nNoLabelT++;
	  else nBadLabelT++;
	} else hESDunmatchedTrackRefsMap[chId]->Fill(trackRef->GetX(), trackRef->GetY());
	
	trackParam = static_cast<AliMUONTrackParam*>(simulatedTrack->GetTrackParamAtCluster()->After(trackParam));
      }
      
    }
    
    // ------------------- reconstructed clusters -------------------
    TIter nextCluster(clusterStore->CreateIterator());
    AliMUONVCluster *cluster;
    while ( ( cluster = static_cast<AliMUONVCluster*>(nextCluster()) ) ) {
      
      Int_t chId = cluster->GetChamberId();
      Int_t deId = cluster->GetDetElemId();
      Int_t label = cluster->GetMCLabel();
      Double_t bestChi2 =  2*sigmaCut*sigmaCut + 1.;
      Int_t bestLabel = -1;
      
      hESDnClustersTotPerCh->Fill(chId);
      
      TIter nextTrackRef(TrackRefStore->CreateIterator());
      AliMUONTrack* simulatedTrack;
      Bool_t matched = kFALSE;
      while ( ( simulatedTrack = static_cast<AliMUONTrack*>(nextTrackRef()) ) ) {
	
	if (useLabel) { // with label
	  
	  if (Int_t(simulatedTrack->GetUniqueID()) == label) matched = kTRUE;
	  
	} else { // without label
	  
	  AliMUONTrackParam* trackParam = static_cast<AliMUONTrackParam*>(simulatedTrack->GetTrackParamAtCluster()->First());
	  while (trackParam) {
	    
	    AliMUONVCluster *trackRef = trackParam->GetClusterPtr();
	    if (trackRef->GetDetElemId() == deId) {
	      Double_t deltaX = cluster->GetX() - trackRef->GetX();
	      Double_t deltaY = cluster->GetY() - trackRef->GetY();
	      Double_t chi2 = deltaX * deltaX / cluster->GetErrX2() + deltaY * deltaY / cluster->GetErrY2();
	      if (chi2 <= 2*sigmaCut*sigmaCut) {
		matched = kTRUE;
		if (chi2<bestChi2) {
		  bestChi2 = chi2;
		  bestLabel = (Int_t) simulatedTrack->GetUniqueID();
		}
	      }
	    }
	    
	    trackParam = static_cast<AliMUONTrackParam*>(simulatedTrack->GetTrackParamAtCluster()->After(trackParam));
	  }
	  
	}
	
	if (useLabel && matched) break;
	
      }
      
      if (matched) {
	hESDnMatchedClustersTotPerCh->Fill(chId);
	if (bestLabel == label) nGoodLabelC++;
	else if (label < 0) nNoLabelC++;
	else nBadLabelC++;
      } else hESDunmatchedClustersMap[chId]->Fill(cluster->GetX(), cluster->GetY());
      
    }
    
    delete clusterStore;
  } // for (Int_t iEvent = FirstEvent;
  
  MUONLoader->UnloadRecPoints();
  
  hESDFracMatchedTrackRefsTotPerCh->Add(hESDnMatchedTrackRefsTotPerCh);
  hESDFracMatchedTrackRefsTotPerCh->Divide(hESDnTrackRefsTotPerCh);
  
  hESDFracMatchedClustersTotPerCh->Add(hESDnMatchedClustersTotPerCh);
  hESDFracMatchedClustersTotPerCh->Divide(hESDnClustersTotPerCh);
  
  hESDnClustersTotOvernTrackRefTotPerCh->Add(hESDnClustersTotPerCh);
  hESDnClustersTotOvernTrackRefTotPerCh->Divide(hESDnTrackRefsTotPerCh);
  
  histoFile->Write();
  histoFile->Close();
  
  cout<<"Warning: the comparison between the 2 methods require to identify the closest cluster or TrackRef"<<endl;
  cout<<"         --> it is meaningless in case of unknown misalignment"<<endl;
  cout<<endl;
  cout<<"Number of events used: "<<nEventUsed<<endl;
  cout<<"Matching TrackRefs with clusters:"<<endl;
  cout<<"cluster label: "<<nGoodLabelT<<" good, "<<nBadLabelT<<" bad, "<<nNoLabelT<<" not labeled"<<endl;
  cout<<"Matching clusters with TrackRefs:"<<endl;
  cout<<"cluster label: "<<nGoodLabelC<<" good, "<<nBadLabelC<<" bad, "<<nNoLabelC<<" not labeled"<<endl;
  
}

