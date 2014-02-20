/*
 *  CheckAcceptance.C
 *  aliroot
 *
 *  Created by Philippe Pillot on 10/08/09.
 *  Copyright 2009 SUBATECH. All rights reserved.
 *
 */

#include <Riostream.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TGeoGlobalMagField.h>

#include "AliMagF.h"

#include "AliMUONConstants.h"
#include "AliMUONTrackParam.h"
#include "AliMUONTrackExtrap.h"

void CheckAcceptance(Int_t nTracksToBeUsed)
{
  /// return the acceptance of the spectrometer according to the requested stations.
  /// (tracks are generated at the last station)
  
  // parameters
  gRandom->SetSeed(0);
  const Bool_t requestStation[5] = {kTRUE, kTRUE, kTRUE, kTRUE, kTRUE};
  const Double_t maxNonBendingSlope = 0.8;
  const Double_t maxBendingSlope = 1.5;
  const Double_t minBendingMomentum = 0.1;
  const Double_t maxBendingMomentum = 2;
  const Bool_t fieldON = kTRUE;
  
  // set  mag field for track extrapolations
  // waiting for mag field in CDB 
  if (!TGeoGlobalMagField::Instance()->GetField()) {
    printf("Loading field map...\n");
    AliMagF* field;
    if (fieldON) field = new AliMagF("Maps","Maps", -1., -1., AliMagF::k5kG);
    else field = new AliMagF("Maps","Maps", 0., 0., AliMagF::k5kG);
    TGeoGlobalMagField::Instance()->SetField(field);
  }
  AliMUONTrackExtrap::SetField();
  
  // initialize histograms
  TFile* f = new TFile("acceptance.root","RECREATE");
  TH1F* hBendingSlope0 = new TH1F("hBendingSlope0", "bending slope distribution at first chamber", 100, -maxBendingSlope, maxBendingSlope);
  TH1F* hBendingSlope10 = new TH1F("hBendingSlope10", "bending slope distribution at last chamber", 100, -maxBendingSlope, maxBendingSlope);
  TH1F* hBendingMomentum = new TH1F("hBendingMomentum", "bending momentum distribution", 100, minBendingMomentum, maxBendingMomentum);
  TH2F* hBendingMomentumVsBendingSlope0 = new TH2F("hBendingMomentumVsBendingSlope0", "bending momentum vs. bending slope distribution at first chamber",
						   100, -maxBendingSlope, maxBendingSlope, 100, minBendingMomentum, maxBendingMomentum);
  TH1F* hNonBendingImpact = new TH1F("hNonBendingImpact", "non bending impact parameter distribution", 100, -1000, 1000);
  TH1F* hBendingImpact = new TH1F("hBendingImpact", "bending impact parameter distribution", 100, -1000, 1000);
  TH2F* hBendingMomentumVsBendingImpact = new TH2F("hBendingMomentumVsBendingImpact", "bending momentum vs. bending impact parameter distribution",
						   100, -1000, 1000, 100, minBendingMomentum, maxBendingMomentum);
  TH1F* hBendingMomentumTrackingFailed = new TH1F("hBendingMomentumTrackingFailed", "bending momentum distribution when the tracking failed", 100, minBendingMomentum, maxBendingMomentum);

  
  // Get starting station
  Int_t st0 = 0;
  if (requestStation[4]) st0 = 4;
  else if (requestStation[3]) st0 = 3;
  else return;
  
  // produce tracks to fill acceptance histograms
  Int_t ntrials = 0;
  Int_t nAcceptedTracks = 0;
  AliMUONTrackParam param;
  while (nAcceptedTracks < nTracksToBeUsed) {
    
    if ((ntrials)%1000 == 0) cout<<"\rProcessing track "<<ntrials<<"   (accepted = "<<nAcceptedTracks<<")"<<flush;
    
    ntrials++;
    
    // Set starting parameters
    param.SetNonBendingCoor(gRandom->Uniform(-AliMUONConstants::Rmax(st0), AliMUONConstants::Rmax(st0)));
    param.SetBendingCoor(gRandom->Uniform(-AliMUONConstants::Rmax(st0), AliMUONConstants::Rmax(st0)));
    param.SetZ(AliMUONConstants::DefaultChamberZ(2*st0));
    param.SetNonBendingSlope(gRandom->Uniform(-maxNonBendingSlope, maxNonBendingSlope));
    Double_t bendingSlope = gRandom->Uniform(-maxBendingSlope, maxBendingSlope);
    param.SetBendingSlope(bendingSlope);
    param.SetInverseBendingMomentum(1./TMath::Sign(gRandom->Uniform(minBendingMomentum, maxBendingMomentum), gRandom->Uniform(-1,1)));
    
    // Check if the track is within acceptance of the requested stations
    Bool_t keepTrack = kTRUE;
    for (Int_t st = st0-1; st >= 0; st--) {
      
      if (!requestStation[st]) continue;
      
      if (!AliMUONTrackExtrap::ExtrapToZ(&param, AliMUONConstants::DefaultChamberZ(2*st)))
	hBendingMomentumTrackingFailed->Fill(1./TMath::Abs(param.GetInverseBendingMomentum()));
      
      if (TMath::Abs(param.GetNonBendingCoor()) > AliMUONConstants::Rmax(st) ||
	  TMath::Abs(param.GetBendingCoor()) > AliMUONConstants::Rmax(st)) {
	keepTrack = kFALSE;
	break;
      }
      
    }
    
    if (!keepTrack) continue;
    nAcceptedTracks++;
    
    // fill histograms with accepted tracks at last chamber
    hBendingSlope10->Fill(bendingSlope);
    
    // fill histograms with accepted tracks at first chamber
    hBendingSlope0->Fill(param.GetBendingSlope());
    Double_t bendingMomentum = 1./TMath::Abs(param.GetInverseBendingMomentum());
    hBendingMomentum->Fill(bendingMomentum);
    hBendingMomentumVsBendingSlope0->Fill(param.GetBendingSlope(), bendingMomentum);
    
    // fill histograms with accepted tracks extrapolated to vertex
    AliMUONTrackExtrap::ExtrapToZ(&param, 0);
    hNonBendingImpact->Fill(param.GetNonBendingCoor());
    hBendingImpact->Fill(param.GetBendingCoor());
    hBendingMomentumVsBendingImpact->Fill(param.GetBendingCoor(), bendingMomentum);
    
  }
  
  // save results
  f->cd();
  hBendingSlope0->Write();
  hBendingSlope10->Write();
  hBendingMomentum->Write();
  hBendingMomentumVsBendingSlope0->Write();
  hNonBendingImpact->Write();
  hBendingImpact->Write();
  hBendingMomentumVsBendingImpact->Write();
  hBendingMomentumTrackingFailed->Write();
  f->Close();
  
  cout<<"\rnumber of trials = "<<ntrials<<"   (accepted tracks = "<<nAcceptedTracks<<")"<<endl;
  
}

//###########################################################################################
void CheckAcceptance2(Int_t nTracksToBeUsed)
{
  /// return the acceptance of the spectrometer according to the requested stations.
  /// (tracks are generated at the vertex)
  
  // parameters
  gRandom->SetSeed(0);
  const Bool_t requestStation[5] = {kTRUE, kTRUE, kTRUE, kTRUE, kTRUE};
  const Double_t maxNonBendingSlope = 0.4;
  const Double_t maxBendingSlope = 0.5;
  const Double_t minBendingMomentum = 0.1;
  const Double_t maxBendingMomentum = 2;
  const Double_t maxImpactParameter = 200;
  const Bool_t fieldON = kTRUE;
  
  // set  mag field for track extrapolations
  // waiting for mag field in CDB 
  if (!TGeoGlobalMagField::Instance()->GetField()) {
    printf("Loading field map...\n");
    AliMagF* field;
    if (fieldON) field = new AliMagF("Maps","Maps",1.,1.,AliMagF::k5kG);
    else field = new AliMagF("Maps","Maps",0.,0.,AliMagF::k5kG);
    TGeoGlobalMagField::Instance()->SetField(field);
  }
  AliMUONTrackExtrap::SetField();
  
  // initialize histograms
  TFile* f = new TFile("acceptance.root","RECREATE");
  TH1F* hBendingSlope0 = new TH1F("hBendingSlope0", "bending slope distribution at vertex", 100, -maxBendingSlope, maxBendingSlope);
  TH1F* hBendingSlope10 = new TH1F("hBendingSlope10", "bending slope distribution at last chamber", 100, -maxBendingSlope, maxBendingSlope);
  TH1F* hBendingMomentum = new TH1F("hBendingMomentum", "bending momentum distribution", 100, minBendingMomentum, maxBendingMomentum);
  TH2F* hBendingMomentumVsBendingSlope0 = new TH2F("hBendingMomentumVsBendingSlope0", "bending momentum vs. bending slope distribution at vertex",
						   100, -maxBendingSlope, maxBendingSlope, 100, minBendingMomentum, maxBendingMomentum);
  TH2F* hBendingMomentumVsBendingSlope10 = new TH2F("hBendingMomentumVsBendingSlope10", "bending momentum vs. bending slope distribution at last chamber",
						   100, -1.5, 1.5, 100, minBendingMomentum, maxBendingMomentum);
  TH1F* hBendingImpact = new TH1F("hBendingImpact", "bending impact parameter distribution", 100, -1000, 1000);
  TH2F* hBendingMomentumVsBendingImpact = new TH2F("hBendingMomentumVsBendingImpact", "bending momentum vs. bending impact parameter distribution",
						   100, -1000, 1000, 100, minBendingMomentum, maxBendingMomentum);
  TH1F* hBendingMomentumTrackingFailed = new TH1F("hBendingMomentumTrackingFailed", "bending momentum distribution when the tracking failed", 100, minBendingMomentum, maxBendingMomentum);
  
  
  // produce tracks to fill acceptance histograms
  Int_t ntrials = 0;
  Int_t nAcceptedTracks = 0;
  AliMUONTrackParam param;
  while (nAcceptedTracks < nTracksToBeUsed) {
    
    if ((ntrials)%1000 == 0) cout<<"\rProcessing track "<<ntrials<<"   (accepted = "<<nAcceptedTracks<<")"<<flush;
    
    ntrials++;
    
    // Set starting parameters
    if (maxImpactParameter > 0) {
      param.SetNonBendingCoor(gRandom->Uniform(-maxImpactParameter, maxImpactParameter));
      param.SetBendingCoor(gRandom->Uniform(-maxImpactParameter, maxImpactParameter));
    } else {
      param.SetNonBendingCoor(0.);
      param.SetBendingCoor(0.);
    }
    Double_t bendingImpact = param.GetBendingCoor();
    param.SetZ(0.);
    param.SetNonBendingSlope(gRandom->Uniform(-maxNonBendingSlope, maxNonBendingSlope));
    Double_t bendingSlope = gRandom->Uniform(-maxBendingSlope, maxBendingSlope);
    param.SetBendingSlope(bendingSlope);
    Double_t bendingMomentum = gRandom->Uniform(minBendingMomentum, maxBendingMomentum);
    param.SetInverseBendingMomentum(1./TMath::Sign(bendingMomentum, gRandom->Uniform(-1.,1.)));
    
    // Check if the track is within acceptance of the requested stations
    Bool_t keepTrack = kTRUE;
    for (Int_t st = 0; st < 5; st++) {
      
      if (!requestStation[st]) continue;
      
      if (!AliMUONTrackExtrap::ExtrapToZ(&param, AliMUONConstants::DefaultChamberZ(2*st)))
	hBendingMomentumTrackingFailed->Fill(1./TMath::Abs(param.GetInverseBendingMomentum()));
      
      if (TMath::Abs(param.GetNonBendingCoor()) > AliMUONConstants::Rmax(st) ||
	  TMath::Abs(param.GetBendingCoor()) > AliMUONConstants::Rmax(st)) {
	keepTrack = kFALSE;
	break;
      }
      
    }
    
    if (!keepTrack) continue;
    nAcceptedTracks++;
    
    // fill histograms with accepted tracks at last chamber
    hBendingSlope10->Fill(param.GetBendingSlope());
    hBendingMomentumVsBendingSlope10->Fill(param.GetBendingSlope(), 1./TMath::Abs(param.GetInverseBendingMomentum()));
    
    // fill histograms with accepted tracks at vertex
    hBendingSlope0->Fill(bendingSlope);
    hBendingMomentum->Fill(bendingMomentum);
    hBendingMomentumVsBendingSlope0->Fill(bendingSlope, bendingMomentum);
    hBendingImpact->Fill(bendingImpact);
    hBendingMomentumVsBendingImpact->Fill(bendingImpact, bendingMomentum);
    
  }
  
  // save results
  f->cd();
  hBendingSlope0->Write();
  hBendingSlope10->Write();
  hBendingMomentum->Write();
  hBendingMomentumVsBendingSlope0->Write();
  hBendingMomentumVsBendingSlope10->Write();
  hBendingImpact->Write();
  hBendingMomentumVsBendingImpact->Write();
  hBendingMomentumTrackingFailed->Write();
  f->Close();
  
  cout<<"\rnumber of trials = "<<ntrials<<"   (accepted tracks = "<<nAcceptedTracks<<")"<<endl;
  
}
