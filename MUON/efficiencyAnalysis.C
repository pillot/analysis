/*
 *  efficiencyAnalysis.C
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
/// \file efficiencyAnalysis.C
/// \brief Macro efficiencyAnalysis.C
///
/// \author Ph. Pillot, Subatech, December 2008

Double_t sigmaCut = 5;

Bool_t IsReconstructible(Int_t chMatched[20], Int_t nCluster);

void efficiencyAnalysis()
{
  //Reset ROOT and connect tree file
  gROOT->Reset();
  
  // File for histograms and histogram booking
  TFile *histoFile = new TFile("efficiencyAnalysis.root", "RECREATE");
  TH1F* hESDnTrackRefsTotPerCh = new TH1F("hESDnTrackRefsTotPerCh", "number of trackRefs per chamber;chamber ID", 10, 0., 10.);
  hESDnTrackRefsTotPerCh->SetFillColor(kRed);
  TH1F* hESDnMatchedTrackRefsTotPerCh = new TH1F("hESDnMatchedTrackRefsTotPerCh", "number of trackRefs matched with cluster per chamber;chamber ID", 10, 0., 10.);
  hESDnMatchedTrackRefsTotPerCh->SetFillColor(kRed);
  TH1F* hESDFracMatchedTrackRefsTotPerCh = new TH1F("hESDFracMatchedTrackRefsTotPerCh", "fraction of trackRefs matched with cluster per chamber;chamber ID", 10, 0., 10.);
  hESDFracMatchedTrackRefsTotPerCh->SetFillColor(kRed);
  
  TH1F* hESDnTrackRefsPerCh = new TH1F("hESDnTrackRefsPerCh", "number of trackRef per chamber per reconstructible track;chamber ID", 10, 0., 10.);
  hESDnTrackRefsPerCh->SetFillColor(kRed);
  TH1F* hESDnMatchedTrackRefsPerCh = new TH1F("hESDnMatchedTrackRefsPerCh", "number of trackRefs matched with cluster per chamber per reconstructible track;chamber ID", 10, 0., 10.);
  hESDnMatchedTrackRefsPerCh->SetFillColor(kRed);
  TH1F* hESDFracMatchedTrackRefsPerCh = new TH1F("hESDFracMatchedTrackRefsPerCh", "fraction of trackRefs matched with cluster per chamber per reconstructible track;chamber ID", 10, 0., 10.);
  hESDFracMatchedTrackRefsPerCh->SetFillColor(kRed);
  TH1F* hESDnMatchedTrackRefsRecoPerCh = new TH1F("hESDnMatchedTrackRefsRecoPerCh", "number of trackRefs matched with cluster per chamber per still reconstructible track;chamber ID", 10, 0., 10.);
  hESDnMatchedTrackRefsRecoPerCh->SetFillColor(kRed);
  TH1F* hESDFracMatchedTrackRefsRecoPerCh = new TH1F("hESDFracMatchedTrackRefsRecoPerCh", "fraction of trackRefs matched with cluster per chamber per still reconstructible track;chamber ID", 10, 0., 10.);
  hESDFracMatchedTrackRefsRecoPerCh->SetFillColor(kRed);
  TH1F* hESDnExpectClustersPerCh = new TH1F("hESDnExpectClustersPerCh", "number of expected clusters per chamber per reconstructible track;chamber ID", 10, 0., 10.);
  hESDnExpectClustersPerCh->SetFillColor(kRed);
  
  TH1F* hESDnClustersPerCh = new TH1F("hESDnClustersPerCh", "number of clusters per chamber per reconstructed track;chamber ID", 10, 0., 10.);
  hESDnClustersPerCh->SetFillColor(kRed);
  TH1F* hESDnMatchedClustersPerCh = new TH1F("hESDnMatchedClustersPerCh", "number of clusters matched with trackRef per chamber per reconstructed track;chamber ID", 10, 0., 10.);
  hESDnMatchedClustersPerCh->SetFillColor(kRed);
  TH1F* hESDFracMatchedClustersPerCh = new TH1F("hESDFracMatchedClustersPerCh", "fraction of clusters matched with trackRef per chamber per reconstructed track;chamber ID", 10, 0., 10.);
  hESDFracMatchedClustersPerCh->SetFillColor(kRed);
  
  TH1F* hESDnClustersTotPerCh = new TH1F("hESDnClustersTotPerCh", "number of clusters per chamber;chamber ID", 10, 0., 10.);
  hESDnClustersTotPerCh->SetFillColor(kRed);
  TH1F* hESDnMatchedClustersTotPerCh = new TH1F("hESDnMatchedClustersTotPerCh", "number of clusters matched with trackRef per chamber;chamber ID", 10, 0., 10.);
  hESDnMatchedClustersTotPerCh->SetFillColor(kRed);
  TH1F* hESDFracMatchedClustersTotPerCh = new TH1F("hESDFracMatchedClustersTotPerCh", "fraction of clusters matched with trackRef per chamber;chamber ID", 10, 0., 10.);
  hESDFracMatchedClustersTotPerCh->SetFillColor(kRed);
  TH1F* hESDnClustersTotOvernTrackRefTotPerCh = new TH1F("hESDnClustersTotOvernTrackRefTotPerCh", "number of clusters over number of trackRef per chamber;chamber ID", 10, 0., 10.);
  hESDnClustersTotOvernTrackRefTotPerCh->SetFillColor(kRed);
  
  TH1F* hESDResidualXInCh[10];
  TH1F* hESDResidualYInCh[10];
  TH2F* hESDunmatchedTrackRefsMap[10];
  TH2F* hESDunmatchedClustersMap[10];
  for (Int_t i = 0; i < 10; i++) {
    hESDResidualXInCh[i] = new TH1F(Form("hESDResidualXInCh%d",i+1), Form("cluster-track residual-X distribution in chamber %d",i+1), 1000, -5., 5.);
    hESDResidualYInCh[i] = new TH1F(Form("hESDResidualYInCh%d",i+1), Form("cluster-track residual-Y distribution in chamber %d",i+1), 1000, -1., 1.);
    Float_t rMax = AliMUONConstants::Rmax(i/2);
    hESDunmatchedClustersMap[i] = new TH2F(Form("hESDunmatchedClustersMap%d",i+1), Form("unmatched cluster position distribution in chamber %d",i+1),1000, -rMax, rMax, 1000, -rMax, rMax);
    hESDunmatchedTrackRefsMap[i] = new TH2F(Form("hESDunmatchedTrackRefsMap%d",i+1), Form("unmatched trackRef position distribution in chamber %d",i+1),1000, -rMax, rMax, 1000, -rMax, rMax);
  }
  
  Int_t nSimuTracksTot = 0, nRecTracksTot = 0;
  
  AliRunLoader * rl = AliRunLoader::Open("galice.root","MUONLoader");
  AliLoader* MUONLoader = rl->GetDetectorLoader("MUON");
  MUONLoader->LoadRecPoints("READ");   
  
  AliMUONRecoCheck rc("AliESDs.root", "./generated/");
  
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
    
    // ------------------- TrackRef -------------------
    AliMUONVTrackStore* TrackRefStore = rc.TrackRefs(iEvent);
    if (TrackRefStore->GetSize() > 2) {
      delete clusterStore;
      continue;
    }
    TIter nextTrackRef(TrackRefStore->CreateIterator());
    AliMUONTrack* simulatedTrack;
    while ( ( simulatedTrack = static_cast<AliMUONTrack*>(nextTrackRef()) ) ) {
      
      AliMUONTrackParam* trackParam = static_cast<AliMUONTrackParam*>(simulatedTrack->GetTrackParamAtCluster()->First());
      while (trackParam) {
	
	AliMUONVCluster *trackRef = trackParam->GetClusterPtr();
	Int_t chId = trackRef->GetChamberId();
	Int_t deId = trackRef->GetDetElemId();
	
	hESDnTrackRefsTotPerCh->Fill(chId);
	
	TIter nextCluster(clusterStore->CreateChamberIterator(chId,chId));
	AliMUONVCluster *cluster;
	Bool_t matched = kFALSE;
	while ( ( cluster = static_cast<AliMUONVCluster*>(nextCluster()) ) ) {
	  
	  if (cluster->GetDetElemId() != deId) continue;
	  
	  Double_t deltaX = cluster->GetX() - trackRef->GetX();
	  Double_t deltaY = cluster->GetY() - trackRef->GetY();
	  Double_t chi2 = deltaX * deltaX / cluster->GetErrX2() + deltaY * deltaY / cluster->GetErrY2();
	  if (chi2 > 2*sigmaCut*sigmaCut) continue;
	  
	  hESDnMatchedTrackRefsTotPerCh->Fill(chId);
	  hESDResidualXInCh[chId]->Fill(deltaX);
	  hESDResidualYInCh[chId]->Fill(deltaY);
	  
	  matched = kTRUE;
	  break;
	}
	
	if (!matched) hESDunmatchedTrackRefsMap[chId]->Fill(trackRef->GetX(), trackRef->GetY());
	
	trackParam = static_cast<AliMUONTrackParam*>(simulatedTrack->GetTrackParamAtCluster()->After(trackParam));
      }
    }
    
    // ------------------- reconstructible tracks -------------------
    AliMUONVTrackStore* reconstructibleTrackStore = rc.ReconstructibleTracks(iEvent);
    TIter nextReconstructibleTrack(reconstructibleTrackStore->CreateIterator());
    AliMUONTrack* reconstructibleTrack;
    while ( ( reconstructibleTrack = static_cast<AliMUONTrack*>(nextReconstructibleTrack()) ) ) {
      
      Int_t chMatched[20];
      Int_t nCluster = 0;
      nSimuTracksTot++;
      
      AliMUONTrackParam* trackParam = static_cast<AliMUONTrackParam*>(reconstructibleTrack->GetTrackParamAtCluster()->First());
      while (trackParam) {
	
	AliMUONVCluster *trackRef = trackParam->GetClusterPtr();
	Int_t chId = trackRef->GetChamberId();
	Int_t deId = trackRef->GetDetElemId();
	chMatched[nCluster] = -1;
	
	hESDnTrackRefsPerCh->Fill(chId);
	
	TIter nextCluster(clusterStore->CreateChamberIterator(chId,chId));
	AliMUONVCluster *cluster;
	while ( ( cluster = static_cast<AliMUONVCluster*>(nextCluster()) ) ) {
	  
	  if (cluster->GetDetElemId() != deId) continue;
	  
	  Double_t deltaX = cluster->GetX() - trackRef->GetX();
	  Double_t deltaY = cluster->GetY() - trackRef->GetY();
	  Double_t chi2 = deltaX * deltaX / cluster->GetErrX2() + deltaY * deltaY / cluster->GetErrY2();
	  if (chi2 > 2*sigmaCut*sigmaCut) continue;
	  
	  chMatched[nCluster] = chId;
	  hESDnMatchedTrackRefsPerCh->Fill(chId);
	  
	  break;
	}
	
	nCluster++;
	trackParam = static_cast<AliMUONTrackParam*>(reconstructibleTrack->GetTrackParamAtCluster()->After(trackParam));
      }
      
      if (IsReconstructible(chMatched, nCluster))
	for (Int_t i = 0; i < nCluster; i++)
	  if (chMatched[i] > -1) hESDnMatchedTrackRefsRecoPerCh->Fill(chMatched[i]);
      
    }
    
    // ------------------- reconstructed tracks -------------------
    AliMUONVTrackStore* reconstructedTrackStore = rc.ReconstructedTracks(iEvent, kFALSE);
    TIter nextReconstructedTrack(reconstructedTrackStore->CreateIterator());
    AliMUONTrack* reconstructedTrack;
    while ( ( reconstructedTrack = static_cast<AliMUONTrack*>(nextReconstructedTrack()) ) ) {
      /*
      // match with simulated tracks
      Int_t nMatchClusters = 0;
      if (!rc.FindCompatibleTrack(*reconstructedTrack, *TrackRefStore, nMatchClusters, kFALSE, sigmaCut)) continue;
      */
      nRecTracksTot++;
      AliMUONTrackParam* trackParam = static_cast<AliMUONTrackParam*>(reconstructedTrack->GetTrackParamAtCluster()->First());
      while (trackParam) {
	
	AliMUONVCluster *cluster = trackParam->GetClusterPtr();
	Int_t chId = cluster->GetChamberId();
	Int_t deId = cluster->GetDetElemId();
	
	hESDnClustersPerCh->Fill(chId);
	
	TIter nextTrackRef(TrackRefStore->CreateIterator());
	AliMUONTrack* simulatedTrack;
	Bool_t matched = kFALSE;
	while ( ( simulatedTrack = static_cast<AliMUONTrack*>(nextTrackRef()) ) ) {
	  
	  AliMUONTrackParam* trackParam = static_cast<AliMUONTrackParam*>(simulatedTrack->GetTrackParamAtCluster()->First());
	  while (trackParam) {
	    
	    AliMUONVCluster *trackRef = trackParam->GetClusterPtr();
	    if (trackRef->GetDetElemId() == deId) {
	      Double_t deltaX = cluster->GetX() - trackRef->GetX();
	      Double_t deltaY = cluster->GetY() - trackRef->GetY();
	      Double_t chi2 = deltaX * deltaX / cluster->GetErrX2() + deltaY * deltaY / cluster->GetErrY2();
	      if (chi2 <= 2*sigmaCut*sigmaCut) {
		hESDnMatchedClustersPerCh->Fill(chId);
		matched = kTRUE;
		break;
	      }
	    }
	    
	    trackParam = static_cast<AliMUONTrackParam*>(simulatedTrack->GetTrackParamAtCluster()->After(trackParam));
	  }
	  
	  if (matched) break;
	}
	
	trackParam = static_cast<AliMUONTrackParam*>(reconstructedTrack->GetTrackParamAtCluster()->After(trackParam));
      }
      
    }
    
    // ------------------- reconstructed clusters -------------------
    TIter nextCluster(clusterStore->CreateIterator());
    AliMUONVCluster *cluster;
    while ( ( cluster = static_cast<AliMUONVCluster*>(nextCluster()) ) ) {
      
      Int_t chId = cluster->GetChamberId();
      Int_t deId = cluster->GetDetElemId();
      
      hESDnClustersTotPerCh->Fill(chId);
      
      TIter nextTrackRef(TrackRefStore->CreateIterator());
      AliMUONTrack* simulatedTrack;
      Bool_t matched = kFALSE;
      while ( ( simulatedTrack = static_cast<AliMUONTrack*>(nextTrackRef()) ) ) {
	
	AliMUONTrackParam* trackParam = static_cast<AliMUONTrackParam*>(simulatedTrack->GetTrackParamAtCluster()->First());
	while (trackParam) {
	  
	  AliMUONVCluster *trackRef = trackParam->GetClusterPtr();
	  if (trackRef->GetDetElemId() == deId) {
	    Double_t deltaX = cluster->GetX() - trackRef->GetX();
	    Double_t deltaY = cluster->GetY() - trackRef->GetY();
	    Double_t chi2 = deltaX * deltaX / cluster->GetErrX2() + deltaY * deltaY / cluster->GetErrY2();
	    if (chi2 <= 2*sigmaCut*sigmaCut) {
	      hESDnMatchedClustersTotPerCh->Fill(chId);
	      matched = kTRUE;
	      break;
	    }
	  }
	  
	  trackParam = static_cast<AliMUONTrackParam*>(simulatedTrack->GetTrackParamAtCluster()->After(trackParam));
	}
	
	if (matched) break;
      }
      
      if (!matched) hESDunmatchedClustersMap[chId]->Fill(cluster->GetX(), cluster->GetY());
      
    }
    
    delete clusterStore;
  } // for (Int_t iEvent = FirstEvent;
  
  MUONLoader->UnloadRecPoints();
  
  hESDFracMatchedTrackRefsTotPerCh->Add(hESDnMatchedTrackRefsTotPerCh);
  hESDFracMatchedTrackRefsTotPerCh->Divide(hESDnTrackRefsTotPerCh);
  
  hESDnTrackRefsPerCh->Scale(1./nSimuTracksTot);
  hESDnExpectClustersPerCh->Add(hESDnTrackRefsPerCh);
  hESDnExpectClustersPerCh->Multiply(hESDFracMatchedTrackRefsTotPerCh);
  hESDnMatchedTrackRefsPerCh->Scale(1./nSimuTracksTot);
  hESDFracMatchedTrackRefsPerCh->Add(hESDnMatchedTrackRefsPerCh);
  hESDFracMatchedTrackRefsPerCh->Divide(hESDnTrackRefsPerCh);
  hESDnMatchedTrackRefsRecoPerCh->Scale(1./nSimuTracksTot);
  hESDFracMatchedTrackRefsRecoPerCh->Add(hESDnMatchedTrackRefsRecoPerCh);
  hESDFracMatchedTrackRefsRecoPerCh->Divide(hESDnTrackRefsPerCh);
  
  hESDnClustersPerCh->Scale(1./nRecTracksTot);
  hESDnMatchedClustersPerCh->Scale(1./nRecTracksTot);
  hESDFracMatchedClustersPerCh->Add(hESDnMatchedClustersPerCh);
  hESDFracMatchedClustersPerCh->Divide(hESDnClustersPerCh);
  hESDFracMatchedClustersTotPerCh->Add(hESDnMatchedClustersTotPerCh);
  hESDFracMatchedClustersTotPerCh->Divide(hESDnClustersTotPerCh);
  hESDnClustersTotOvernTrackRefTotPerCh->Add(hESDnClustersTotPerCh);
  hESDnClustersTotOvernTrackRefTotPerCh->Divide(hESDnTrackRefsTotPerCh);
  
  histoFile->Write();
  histoFile->Close();
  
}

Bool_t IsReconstructible(Int_t chMatched[20], Int_t nCluster)
{
  // Check track reconstructibility
  
  Int_t n1 = 0, n2 = 0, n3 = 0, n45 = 0;
  for (Int_t i = 0; i < nCluster; i++) {
    switch (chMatched[i]) {
      case 0:
      case 1:
	n1++;
	break;
      case 2:
      case 3:
	n2++;
	break;
      case 4:
      case 5:
	n3++;
	break;
      case 6:
      case 7:
      case 8:
      case 9:
	n45++;
	break;
    }
  }
  
  if (n1 > 0 && n2 > 0 && n3 > 0 && n45 > 2) return kTRUE;
  
  return kFALSE;
  
}

