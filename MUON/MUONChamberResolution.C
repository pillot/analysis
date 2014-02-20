/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

/// \ingroup macros
/// \file MUONChamberResolution.C
/// \brief Macro to compute the resolution of clusters attached to tracks
///
/// \author Philippe Pillot, SUBATECH

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TStopwatch.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TMultiGraph.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <Riostream.h>
#include <TString.h>

// STEER includes
#include "AliESDEvent.h"
#include "AliESDMuonTrack.h"
#include "AliCDBManager.h"

// MUON includes
#include "AliMUONCDB.h"
#include "AliMUONRecoParam.h"
#include "AliMUONESDInterface.h"
#include "AliMUONVTrackReconstructor.h"
#include "AliMUONTrack.h"
#include "AliMUONTrackParam.h"
#include "AliMUONTrackExtrap.h"
#include "AliMUONVCluster.h"
#include "AliMUONConstants.h"
#endif

void FillResiduals(Int_t extrapMode, Int_t nevents, TTree* esdTree, AliESDEvent* esd, Double_t clusterResNB[10], Double_t clusterResB[10]);
void ComputeResolution(Bool_t correctForSystematics, Double_t* clusterResNB = 0x0, Double_t* clusterResB = 0x0, Double_t* clusterResNBErr = 0x0, Double_t* clusterResBErr = 0x0);
TTree* GetESDTree(TFile *esdFile);
void SetClusterResolution(AliMUONTrack& track, Double_t clusterResNB[10], Double_t clusterResB[10]);
void Zoom(TH1* h, Double_t fractionCut = 0.01);
void ZoomLeft(TH1* h, Double_t fractionCut = 0.02);
void ZoomRight(TH1* h, Double_t fractionCut = 0.02);
void GetMean(TH1* h, Double_t& mean, Double_t& meanErr, TGraphErrors* g = 0x0, Int_t i = 0, Double_t x = 0, Bool_t zoom = kTRUE);
void GetRMS(TH1* h, Double_t& rms, Double_t& rmsErr, TGraphErrors* g = 0x0, Int_t i = 0, Double_t x = 0, Bool_t zoom = kTRUE);
void FillSigmaClusterVsP(TH2* hIn, TH2* hOut, TGraphErrors* g, Bool_t zoom = kTRUE);

//-----------------------------------------------------------------------
void MUONChamberResolution(Int_t nSteps = 5, Int_t extrapMode = 1, Bool_t correctForSystematics = kTRUE, Int_t nevents = -1,
			   const char* esdFileNameIn = "AliESDs.root", const char* ocdbPath = "local://$ALICE_ROOT/OCDB")
{
  /// Compute the cluster resolution by studying cluster-track residual, deconvoluting from track resolution
  /// if extrapMode == 0: extrapolate from the closest cluster
  /// if extrapMode == 1: extrapolate from the previous cluster except between stations 2-3-4
  
  nSteps = TMath::Max(nSteps,1);
  if (extrapMode != 0 && extrapMode != 1) {
    Error("MUONChamberResolution", "wrong extrapolation mode");
    return;
  }
  
  // open the ESD file and tree
  TFile* esdFile = TFile::Open(esdFileNameIn);
  TTree* esdTree = GetESDTree(esdFile);
  
  // connect ESD event to the ESD tree
  AliESDEvent* esd = new AliESDEvent();
  esd->ReadFromTree(esdTree);
  
  // get run number
  if (esdTree->GetEvent(0) <= 0) {
    Error("MUONChamberResolution", "no ESD object found for event 0");
    return;
  }
  Int_t runNumber = esd->GetRunNumber();
  
  // load necessary data from OCDB
  AliCDBManager::Instance()->SetDefaultStorage(ocdbPath);
  AliCDBManager::Instance()->SetRun(runNumber);
  if (!AliMUONCDB::LoadField()) return;
  AliMUONRecoParam* recoParam = AliMUONCDB::LoadRecoParam();
  if (!recoParam) return;
  
  // reset tracker for computing track parameters at cluster
  AliMUONESDInterface::ResetTracker(recoParam);
  
  // set default chamber resolution
  //Double_t clusterResNB[10] = {0.059, 0.059, 0.071, 0.070, 0.070, 0.071, 0.070, 0.070, 0.069, 0.069};
  //Double_t clusterResB[10] = {0.0056, 0.0047, 0.0073, 0.0074, 0.0104, 0.0101, 0.0079, 0.0072, 0.0073, 0.0076};
  //Double_t clusterResNB[10] = {0.061, 0.058, 0.072, 0.074, 0.070, 0.071, 0.071, 0.071, 0.070, 0.072};
  //Double_t clusterResB[10]  = {0.0069, 0.0062, 0.0085, 0.0081, 0.0100, 0.0095, 0.0118, 0.0091, 0.0096, 0.0106};
  //chamberResolution_clostestAndPrev-St23.root:
  //Double_t clusterResNB[10] = {0.061, 0.058, 0.071, 0.074, 0.070, 0.071, 0.071, 0.071, 0.070, 0.072};
  //Double_t clusterResB[10]  = {0.0069, 0.0062, 0.0090, 0.0078, 0.0102, 0.0094, 0.0106, 0.0094, 0.0096, 0.0105};
  //chamberResolution_clostestAndPrev.root:
  //Double_t clusterResNB[10] = {0.061, 0.058, 0.071, 0.074, 0.070, 0.071, 0.071, 0.071, 0.070, 0.072};
  //Double_t clusterResB[10]  = {0.0069, 0.0062, 0.0090, 0.0078, 0.0100, 0.0096, 0.0111, 0.0093, 0.0096, 0.0105};
  //chamberResolution_Prev-St23.root:
  //Double_t clusterResNB[10] = {0.061, 0.058, 0.072, 0.074, 0.070, 0.071, 0.071, 0.071, 0.070, 0.072};
  //Double_t clusterResB[10]  = {0.0069, 0.0062, 0.0083, 0.0082, 0.0102, 0.0094, 0.0106, 0.0094, 0.0096, 0.0105};
  //chamberResolution_mode1_pAbove5GeV.root:
  //Double_t clusterResNB[10] = {0.061, 0.058, 0.071, 0.074, 0.070, 0.071, 0.071, 0.071, 0.071, 0.072};
  //Double_t clusterResB[10]  = {0.0063, 0.0059, 0.0081, 0.0081, 0.0101, 0.0093, 0.0100, 0.0091, 0.0094, 0.0099};
  // systematics
  //Double_t clusterResNB[10] = {0.061, 0.058, 0.218, 0.074, 0.070, 0.124, 0.071, 0.071, 0.071, 0.311};
  //Double_t clusterResB[10]  = {0.0063, 0.0059, 0.1, 0.0081, 0.0101, 0.05, 0.0100, 0.1, 0.0094, 0.0099};
  //Double_t clusterResNB[10] = {0.061, 0.058, 0.512, 0.074, 0.070, 0.126, 0.071, 0.071, 0.070, 0.311};
  //Double_t clusterResB[10]  = {0.0069, 0.0062, 0.1, 0.0082, 0.0102, 0.2, 0.0106, 0.1, 0.0096, 0.0105};
  Double_t clusterResNB[10];
  Double_t clusterResB[10];
  Double_t clusterResNBErr[10];
  Double_t clusterResBErr[10];
  for (Int_t i = 0; i < AliMUONConstants::NTrackingCh(); i++) {
    clusterResNB[i] = recoParam->GetDefaultNonBendingReso(i);
    clusterResB[i] = recoParam->GetDefaultBendingReso(i);
    clusterResNBErr[i] = 0.;
    clusterResBErr[i] = 0.;
  }
  
  // output graphs
  TMultiGraph* mgClusterResXVsStep = new TMultiGraph("mgClusterResXVsStep","cluster X-resolution versus step;step;#sigma_{X} (cm)");
  TMultiGraph* mgClusterResYVsStep = new TMultiGraph("mgClusterResYVsStep","cluster Y-resolution versus step;step;#sigma_{Y} (cm)");
  TGraphErrors* gClusterResXVsStep[10];
  TGraphErrors* gClusterResYVsStep[10];
  for (Int_t i = 0; i < AliMUONConstants::NTrackingCh(); i++) {
    gClusterResXVsStep[i] = new TGraphErrors(nSteps+1);
    gClusterResXVsStep[i]->SetName(Form("gResX_ch%d",i+1));
    gClusterResXVsStep[i]->SetMarkerStyle(kFullDotMedium);
    gClusterResXVsStep[i]->SetMarkerColor(i+1+i/9);
    gClusterResXVsStep[i]->SetPoint(0, 0, clusterResNB[i]);
    gClusterResXVsStep[i]->SetPointError(0, 0., clusterResNBErr[i]);
    mgClusterResXVsStep->Add(gClusterResXVsStep[i],"lp");
    
    gClusterResYVsStep[i] = new TGraphErrors(nSteps+1);
    gClusterResYVsStep[i]->SetName(Form("gResY_ch%d",i+1));
    gClusterResYVsStep[i]->SetMarkerStyle(kFullDotMedium);
    gClusterResYVsStep[i]->SetMarkerColor(i+1+i/9);
    gClusterResYVsStep[i]->SetPoint(0, 0, clusterResB[i]);
    gClusterResYVsStep[i]->SetPointError(0, 0., clusterResBErr[i]);
    mgClusterResYVsStep->Add(gClusterResYVsStep[i],"lp");
  }
  
  // timer start...
  TStopwatch timer;
  
  // loop over step
  for (Int_t iStep = 0; iStep < nSteps; iStep++) {
    cout<<"step "<<iStep+1<<"/"<<nSteps<<endl;
    
    // refit the tracks and fill residual distributions
    FillResiduals(extrapMode, nevents, esdTree, esd, clusterResNB, clusterResB);
    
    // compute new clusters resolution from residuals
    ComputeResolution(correctForSystematics, clusterResNB, clusterResB, clusterResNBErr, clusterResBErr);
    
    // fill graphs
    for (Int_t i = 0; i < AliMUONConstants::NTrackingCh(); i++) {
      gClusterResXVsStep[i]->SetPoint(iStep+1, iStep+1, clusterResNB[i]);
      gClusterResXVsStep[i]->SetPointError(iStep+1, 0., clusterResNBErr[i]);
      gClusterResYVsStep[i]->SetPoint(iStep+1, iStep+1, clusterResB[i]);
      gClusterResYVsStep[i]->SetPointError(iStep+1, 0., clusterResBErr[i]);
    }
    
  }
  
  // ...timer stop
  timer.Stop();
  cout<<endl<<"processing time: R:"<<timer.RealTime()<<" C:"<<timer.CpuTime()<<endl<<endl;
  
  // display
  TCanvas* convergence = new TCanvas("convergence","convergence");
  convergence->Divide(1,2);
  convergence->cd(1);
  mgClusterResXVsStep->Draw("ap");
  convergence->cd(2);
  mgClusterResYVsStep->Draw("ap");
  
  // save results
  TFile* outFile = TFile::Open("chamberResolution.root","UPDATE");
  outFile->cd();
  mgClusterResXVsStep->Write();
  mgClusterResYVsStep->Write();
  convergence->Write();
  outFile->Close();
  
  // print results
  printf("\nchamber resolution:\n");
  printf(" - non-bending:");
  for (Int_t i = 0; i < AliMUONConstants::NTrackingCh(); i++) printf((i==0)?" %5.3f":", %5.3f",clusterResNB[i]);
  printf("\n -     bending:");
  for (Int_t i = 0; i < AliMUONConstants::NTrackingCh(); i++) printf((i==0)?" %6.4f":", %6.4f",clusterResB[i]);
  printf("\n\n");
  
  // free memory
  esdFile->Close();
  delete esd;
  
}

//-----------------------------------------------------------------------
void FillResiduals(Int_t extrapMode, Int_t nevents, TTree* esdTree, AliESDEvent* esd, Double_t clusterResNB[10], Double_t clusterResB[10])
{
  /// Fill the cluster-track residuals
  
  // get tracker to refit
  AliMUONVTrackReconstructor* tracker = AliMUONESDInterface::GetTracker();
  
  // find the highest chamber resolution and set histogram bins
  const AliMUONRecoParam* recoParam = tracker->GetRecoParam();
  Double_t maxSigma[2] = {-1., -1.};
  for (Int_t i = 0; i < 10; i++) {
    if (recoParam->GetDefaultNonBendingReso(i) > maxSigma[0]) maxSigma[0] = recoParam->GetDefaultNonBendingReso(i);
    if (recoParam->GetDefaultBendingReso(i) > maxSigma[1]) maxSigma[1] = recoParam->GetDefaultBendingReso(i);
  }
  const Int_t nBins = 2000;
  const Int_t nSigma = 10;
  Double_t maxRes[2] = {nSigma*maxSigma[0], nSigma*maxSigma[1]};
  
  // define bining in momentum
  const Int_t pNBins = 20;
  const Double_t pEdges[2] = {0., 300.};
  
  // output
  TFile* outFile = TFile::Open("chamberResolution.root","RECREATE");
  outFile->mkdir("residuals","residuals");
  outFile->cd("residuals");
  TH1F* hResidualXInCh_ClusterIn[AliMUONConstants::NTrackingCh()];
  TH1F* hResidualYInCh_ClusterIn[AliMUONConstants::NTrackingCh()];
  TH1F* hResidualXInCh_ClusterOut[AliMUONConstants::NTrackingCh()];
  TH1F* hResidualYInCh_ClusterOut[AliMUONConstants::NTrackingCh()];
  for (Int_t i = 0; i < AliMUONConstants::NTrackingCh(); i++) {
    hResidualXInCh_ClusterIn[i] = new TH1F(Form("hResidualXInCh%d_ClusterIn",i+1), Form("cluster-track residual-X distribution in chamber %d (cluster attached to the track);#Delta_{X} (cm)",i+1), nBins, -maxRes[0], maxRes[0]);
    hResidualYInCh_ClusterIn[i] = new TH1F(Form("hResidualYInCh%d_ClusterIn",i+1), Form("cluster-track residual-Y distribution in chamber %d (cluster attached to the track);#Delta_{Y} (cm)",i+1), nBins, -maxRes[1], maxRes[1]);
    hResidualXInCh_ClusterOut[i] = new TH1F(Form("hResidualXInCh%d_ClusterOut",i+1), Form("cluster-track residual-X distribution in chamber %d (cluster not attached to the track);#Delta_{X} (cm)",i+1), nBins, -2.*maxRes[0], 2.*maxRes[0]);
    hResidualYInCh_ClusterOut[i] = new TH1F(Form("hResidualYInCh%d_ClusterOut",i+1), Form("cluster-track residual-Y distribution in chamber %d (cluster not attached to the track);#Delta_{Y} (cm)",i+1), nBins, -2.*maxRes[1], 2.*maxRes[1]);
  }
  
  TH2F *hTrackResXPerCh = new TH2F("hTrackResXPerCh","track #sigma_{X} per Ch;chamber ID;#sigma_{X} (cm)", 10, 0.5, 10.5, nBins, 0., maxRes[0]);
  TH2F *hTrackResYPerCh = new TH2F("hTrackResYPerCh","track #sigma_{Y} per Ch;chamber ID;#sigma_{Y} (cm)", 10, 0.5, 10.5, nBins, 0., maxRes[1]);
  
  TH2F *hMCSXPerCh = new TH2F("hMCSXPerCh","MCS X-dispersion of extrapolated track per Ch;chamber ID;#sigma_{X} (cm)", 10, 0.5, 10.5, nBins, 0., 0.2);
  TH2F *hMCSYPerCh = new TH2F("hMCSYPerCh","MCS Y-dispersion of extrapolated track per Ch;chamber ID;#sigma_{Y} (cm)", 10, 0.5, 10.5, nBins, 0., 0.2);
  
  TH2F *hClusterRes2XPerCh = new TH2F("hClusterRes2XPerCh","cluster #sigma_{X}^{2} per Ch;chamber ID;#sigma_{X}^{2} (cm^{2})", 10, 0.5, 10.5, nSigma*nBins, -0.1*maxRes[0]*maxRes[0], maxRes[0]*maxRes[0]);
  TH2F *hClusterRes2YPerCh = new TH2F("hClusterRes2YPerCh","cluster #sigma_{Y}^{2} per Ch;chamber ID;#sigma_{Y}^{2} (cm^{2})", 10, 0.5, 10.5, nSigma*nBins, -0.1*maxRes[1]*maxRes[1], maxRes[1]*maxRes[1]);
  
  outFile->mkdir("residualsVsP","residualsVsP");
  outFile->cd("residualsVsP");
  TH2F* hResidualXInChVsP_ClusterIn[AliMUONConstants::NTrackingCh()];
  TH2F* hResidualYInChVsP_ClusterIn[AliMUONConstants::NTrackingCh()];
  TH2F* hResidualXInChVsP_ClusterOut[AliMUONConstants::NTrackingCh()];
  TH2F* hResidualYInChVsP_ClusterOut[AliMUONConstants::NTrackingCh()];
  for (Int_t i = 0; i < AliMUONConstants::NTrackingCh(); i++) {
    hResidualXInChVsP_ClusterIn[i] = new TH2F(Form("hResidualXInCh%dVsP_ClusterIn",i+1), Form("cluster-track residual-X distribution in chamber %d versus momentum (cluster attached to the track);p (GeV/c^{2});#Delta_{X} (cm)",i+1), pNBins, pEdges[0], pEdges[1], nBins, -maxRes[0], maxRes[0]);
    hResidualYInChVsP_ClusterIn[i] = new TH2F(Form("hResidualYInCh%dVsP_ClusterIn",i+1), Form("cluster-track residual-Y distribution in chamber %d versus momentum (cluster attached to the track);p (GeV/c^{2});#Delta_{Y} (cm)",i+1), pNBins, pEdges[0], pEdges[1], nBins, -maxRes[1], maxRes[1]);
    hResidualXInChVsP_ClusterOut[i] = new TH2F(Form("hResidualXInCh%dVsP_ClusterOut",i+1), Form("cluster-track residual-X distribution in chamber %d versus momentum (cluster not attached to the track);p (GeV/c^{2});#Delta_{X} (cm)",i+1), pNBins, pEdges[0], pEdges[1], nBins, -2.*maxRes[0], 2.*maxRes[0]);
    hResidualYInChVsP_ClusterOut[i] = new TH2F(Form("hResidualYInCh%dVsP_ClusterOut",i+1), Form("cluster-track residual-Y distribution in chamber %d versus momentum (cluster not attached to the track);p (GeV/c^{2});#Delta_{Y} (cm)",i+1), pNBins, pEdges[0], pEdges[1], nBins, -2.*maxRes[1], 2.*maxRes[1]);
  }
  
  // Loop over ESD events
  if (nevents > 0) nevents = TMath::Min(nevents,(Int_t)esdTree->GetEntries());
  else nevents = (Int_t)esdTree->GetEntries();
  for (Int_t iEvent = 0; iEvent < nevents; iEvent++) {
    if ((iEvent+1)%10 == 0) cout<<"\rEvent processing... "<<iEvent+1<<flush;
    
    // get the ESD of current event
    esdTree->GetEvent(iEvent);
    if (!esd) {
      Error("FillResiduals", "no ESD object found for event %d", iEvent);
      return;
    }
    
    // loop over the ESD tracks
    Int_t nTracks = (Int_t)esd->GetNumberOfMuonTracks() ;
    for (Int_t iTrack = 0; iTrack <  nTracks;  iTrack++) {
      
      // get the ESD track
      AliESDMuonTrack* esdTrack = esd->GetMuonTrack(iTrack);
      
      // skip ghost tracks
      if (!esdTrack->ContainTrackerData()) continue;
      
      // skip low momentum tracks
      if (esdTrack->PUncorrected() < 5.) continue;
      
      // get the corresponding MUON track
      AliMUONTrack track;
      AliMUONESDInterface::ESDToMUON(*esdTrack, track, kFALSE);
      
      // change the cluster resolution
      SetClusterResolution(track, clusterResNB, clusterResB);
      
      // refit the track
      if (!tracker->RefitTrack(track, kFALSE)) break;
      
      // save track unchanged
      AliMUONTrack referenceTrack(track);
      Double_t p = static_cast<AliMUONTrackParam*>(track.GetTrackParamAtCluster()->First())->P();
      
      // loop over clusters
      Int_t nClusters = track.GetNClusters();
      for (Int_t iCluster=0; iCluster<nClusters; iCluster++) {
	
	// Get current, previous and next trackParams
	AliMUONTrackParam* trackParam = static_cast<AliMUONTrackParam*>(track.GetTrackParamAtCluster()->UncheckedAt(iCluster));
	AliMUONTrackParam* previousTrackParam = static_cast<AliMUONTrackParam*>(track.GetTrackParamAtCluster()->Before(trackParam));
	AliMUONTrackParam* nextTrackParam = static_cast<AliMUONTrackParam*>(track.GetTrackParamAtCluster()->After(trackParam));
	
	// save current trackParam and remove it from the track
	AliMUONTrackParam currentTrackParam(*trackParam);
	track.RemoveTrackParamAtCluster(trackParam);
	
	// get cluster info
	AliMUONVCluster* cluster = currentTrackParam.GetClusterPtr();
	Int_t chId = cluster->GetChamberId();
	
	// make sure the track has another cluster in the same station and can still be refitted
	Bool_t refit = track.IsValid( 1 << (chId/2) );
	if (refit) {
	  
	  // refit the track and proceed if everything goes fine
	  if (tracker->RefitTrack(track, kFALSE)) {
	    
	    // compute residuals
	    AliMUONTrackParam* referenceTrackParam = static_cast<AliMUONTrackParam*>(referenceTrack.GetTrackParamAtCluster()->UncheckedAt(iCluster));
	    Double_t deltaX = cluster->GetX() - referenceTrackParam->GetNonBendingCoor();
	    Double_t deltaY = cluster->GetY() - referenceTrackParam->GetBendingCoor();
	    
	    // fill histograms of residuals with cluster still attached to the track
	    hResidualXInCh_ClusterIn[chId]->Fill(deltaX);
	    hResidualYInCh_ClusterIn[chId]->Fill(deltaY);
	    hResidualXInChVsP_ClusterIn[chId]->Fill(p, deltaX);
	    hResidualYInChVsP_ClusterIn[chId]->Fill(p, deltaY);
	    
	    // find the track parameters closest to the current cluster position
	    Double_t dZWithPrevious = (previousTrackParam) ? TMath::Abs(previousTrackParam->GetClusterPtr()->GetZ() - cluster->GetZ()) : FLT_MAX;
	    Int_t previousChId = (previousTrackParam) ? previousTrackParam->GetClusterPtr()->GetChamberId() : -1;
	    Double_t dZWithNext = (nextTrackParam) ? TMath::Abs(nextTrackParam->GetClusterPtr()->GetZ() - cluster->GetZ()) : FLT_MAX;
	    AliMUONTrackParam* startingTrackParam = 0x0;
	    if ((extrapMode == 0 && dZWithPrevious < dZWithNext) ||
		(extrapMode == 1 && previousTrackParam && !(chId/2 == 2 && previousChId/2 == 1) &&
		 !(chId/2 == 3 && previousChId/2 == 2))) startingTrackParam = previousTrackParam;
	    else startingTrackParam = nextTrackParam;
	    
	    // reset current parameters
	    currentTrackParam.SetParameters(startingTrackParam->GetParameters());
	    currentTrackParam.SetZ(startingTrackParam->GetZ());
	    currentTrackParam.SetCovariances(startingTrackParam->GetCovariances());
	    currentTrackParam.ResetPropagator();
	    
	    // extrapolate to the current cluster position and fill histograms of residuals if everything goes fine
	    if (AliMUONTrackExtrap::ExtrapToZCov(&currentTrackParam, currentTrackParam.GetClusterPtr()->GetZ(), kTRUE)) {
	      
	      // compute MCS dispersion if starting from next chamber
	      TMatrixD mcsCov(5,5);
	      if (startingTrackParam == nextTrackParam && chId == 0) {
		AliMUONTrackParam trackParamForMCS;
		trackParamForMCS.SetParameters(nextTrackParam->GetParameters());
		AliMUONTrackExtrap::AddMCSEffect(&trackParamForMCS,AliMUONConstants::ChamberThicknessInX0(nextTrackParam->GetClusterPtr()->GetChamberId()),-1.);
		const TMatrixD &propagator = currentTrackParam.GetPropagator();
		TMatrixD tmp(trackParamForMCS.GetCovariances(),TMatrixD::kMultTranspose,propagator);
		mcsCov.Mult(propagator,tmp);
	      } else mcsCov.Zero();
	      
	      // compute residuals
	      Double_t trackResX2 = currentTrackParam.GetCovariances()(0,0) + mcsCov(0,0);
	      Double_t trackResY2 = currentTrackParam.GetCovariances()(2,2) + mcsCov(2,2);
	      deltaX = cluster->GetX() - currentTrackParam.GetNonBendingCoor();
	      deltaY = cluster->GetY() - currentTrackParam.GetBendingCoor();
	      
	      // fill histograms with cluster not attached to the track
	      hResidualXInCh_ClusterOut[chId]->Fill(deltaX);
	      hResidualYInCh_ClusterOut[chId]->Fill(deltaY);
	      hResidualXInChVsP_ClusterOut[chId]->Fill(p, deltaX);
	      hResidualYInChVsP_ClusterOut[chId]->Fill(p, deltaY);
	      hTrackResXPerCh->Fill(chId+1,TMath::Sqrt(trackResX2));
	      hTrackResYPerCh->Fill(chId+1,TMath::Sqrt(trackResY2));
	      hMCSXPerCh->Fill(chId+1,TMath::Sqrt(mcsCov(0,0)));
	      hMCSYPerCh->Fill(chId+1,TMath::Sqrt(mcsCov(2,2)));
	      hClusterRes2XPerCh->Fill(chId+1, deltaX*deltaX - trackResX2);
	      hClusterRes2YPerCh->Fill(chId+1, deltaY*deltaY - trackResY2);
	    }
	    
	  }
	  
	}
	
	// restore the track
	track.AddTrackParamAtCluster(currentTrackParam, *(currentTrackParam.GetClusterPtr()), kTRUE);
	
      }
      
    }
    
  }
  cout<<"\rEvent processing... "<<nevents<<" done"<<endl;
  
  // fill output
  outFile->Write();
  outFile->Close();
  
}

//-----------------------------------------------------------------------
void ComputeResolution(Bool_t correctForSystematics, Double_t* clusterResNB, Double_t* clusterResB, Double_t* clusterResNBErr, Double_t* clusterResBErr)
{
  /// Compute cluster resolution
  
  // output
  TFile* outFile = TFile::Open("chamberResolution.root","UPDATE");
  if (!outFile || !outFile->IsOpen()) {
    cout<<"you must first fill the residuals"<<endl;
    return;
  }
  
  TH1F* hResidualXInCh_ClusterIn[AliMUONConstants::NTrackingCh()];
  TH1F* hResidualYInCh_ClusterIn[AliMUONConstants::NTrackingCh()];
  TH1F* hResidualXInCh_ClusterOut[AliMUONConstants::NTrackingCh()];
  TH1F* hResidualYInCh_ClusterOut[AliMUONConstants::NTrackingCh()];
  for (Int_t i = 0; i < AliMUONConstants::NTrackingCh(); i++) {
    hResidualXInCh_ClusterIn[i] = static_cast<TH1F*>(outFile->FindObjectAny(Form("hResidualXInCh%d_ClusterIn",i+1)));
    hResidualYInCh_ClusterIn[i] = static_cast<TH1F*>(outFile->FindObjectAny(Form("hResidualYInCh%d_ClusterIn",i+1)));
    hResidualXInCh_ClusterOut[i] = static_cast<TH1F*>(outFile->FindObjectAny(Form("hResidualXInCh%d_ClusterOut",i+1)));
    hResidualYInCh_ClusterOut[i] = static_cast<TH1F*>(outFile->FindObjectAny(Form("hResidualYInCh%d_ClusterOut",i+1)));
  }
  
  TH2F* hResidualXInChVsP_ClusterIn[AliMUONConstants::NTrackingCh()];
  TH2F* hResidualYInChVsP_ClusterIn[AliMUONConstants::NTrackingCh()];
  TH2F* hResidualXInChVsP_ClusterOut[AliMUONConstants::NTrackingCh()];
  TH2F* hResidualYInChVsP_ClusterOut[AliMUONConstants::NTrackingCh()];
  for (Int_t i = 0; i < AliMUONConstants::NTrackingCh(); i++) {
    hResidualXInChVsP_ClusterIn[i] = static_cast<TH2F*>(outFile->FindObjectAny(Form("hResidualXInCh%dVsP_ClusterIn",i+1)));
    hResidualYInChVsP_ClusterIn[i] = static_cast<TH2F*>(outFile->FindObjectAny(Form("hResidualYInCh%dVsP_ClusterIn",i+1)));
    hResidualXInChVsP_ClusterOut[i] = static_cast<TH2F*>(outFile->FindObjectAny(Form("hResidualXInCh%dVsP_ClusterOut",i+1)));
    hResidualYInChVsP_ClusterOut[i] = static_cast<TH2F*>(outFile->FindObjectAny(Form("hResidualYInCh%dVsP_ClusterOut",i+1)));
  }
  
  TGraphErrors* gResidualXPerChMean_ClusterIn = static_cast<TGraphErrors*>(outFile->FindObjectAny("gResidualXPerChMean_ClusterIn"));
  if (!gResidualXPerChMean_ClusterIn) {
    gResidualXPerChMean_ClusterIn = new TGraphErrors(AliMUONConstants::NTrackingCh());
    gResidualXPerChMean_ClusterIn->SetName("gResidualXPerChMean_ClusterIn");
    gResidualXPerChMean_ClusterIn->SetTitle("cluster-track residual-X per Ch: mean (cluster in);chamber ID;<#Delta_{X}> (cm)");
    gResidualXPerChMean_ClusterIn->SetMarkerStyle(kFullDotLarge);
  }
  TGraphErrors* gResidualYPerChMean_ClusterIn = static_cast<TGraphErrors*>(outFile->FindObjectAny("gResidualYPerChMean_ClusterIn"));
  if (!gResidualYPerChMean_ClusterIn) {
    gResidualYPerChMean_ClusterIn = new TGraphErrors(AliMUONConstants::NTrackingCh());
    gResidualYPerChMean_ClusterIn->SetName("gResidualYPerChMean_ClusterIn");
    gResidualYPerChMean_ClusterIn->SetTitle("cluster-track residual-Y per Ch: mean (cluster in);chamber ID;<#Delta_{Y}> (cm)");
    gResidualYPerChMean_ClusterIn->SetMarkerStyle(kFullDotLarge);
  }
  TGraphErrors* gResidualXPerChMean_ClusterOut = static_cast<TGraphErrors*>(outFile->FindObjectAny("gResidualXPerChMean_ClusterOut"));
  if (!gResidualXPerChMean_ClusterOut) {
    gResidualXPerChMean_ClusterOut = new TGraphErrors(AliMUONConstants::NTrackingCh());
    gResidualXPerChMean_ClusterOut->SetName("gResidualXPerChMean_ClusterOut");
    gResidualXPerChMean_ClusterOut->SetTitle("cluster-track residual-X per Ch: mean (cluster out);chamber ID;<#Delta_{X}> (cm)");
    gResidualXPerChMean_ClusterOut->SetMarkerStyle(kFullDotLarge);
  }
  TGraphErrors* gResidualYPerChMean_ClusterOut = static_cast<TGraphErrors*>(outFile->FindObjectAny("gResidualYPerChMean_ClusterOut"));
  if (!gResidualYPerChMean_ClusterOut) {
    gResidualYPerChMean_ClusterOut = new TGraphErrors(AliMUONConstants::NTrackingCh());
    gResidualYPerChMean_ClusterOut->SetName("gResidualYPerChMean_ClusterOut");
    gResidualYPerChMean_ClusterOut->SetTitle("cluster-track residual-Y per Ch: mean (cluster out);chamber ID;<#Delta_{Y}> (cm)");
    gResidualYPerChMean_ClusterOut->SetMarkerStyle(kFullDotLarge);
  }
  TGraphErrors* gResidualXPerChSigma_ClusterIn = static_cast<TGraphErrors*>(outFile->FindObjectAny("gResidualXPerChSigma_ClusterIn"));
  if (!gResidualXPerChSigma_ClusterIn) {
    gResidualXPerChSigma_ClusterIn = new TGraphErrors(AliMUONConstants::NTrackingCh());
    gResidualXPerChSigma_ClusterIn->SetName("gResidualXPerChSigma_ClusterIn");
    gResidualXPerChSigma_ClusterIn->SetTitle("cluster-track residual-X per Ch: sigma (cluster in);chamber ID;#sigma_{X} (cm)");
    gResidualXPerChSigma_ClusterIn->SetMarkerStyle(kFullDotLarge);
  }
  TGraphErrors* gResidualYPerChSigma_ClusterIn = static_cast<TGraphErrors*>(outFile->FindObjectAny("gResidualYPerChSigma_ClusterIn"));
  if (!gResidualYPerChSigma_ClusterIn) {
    gResidualYPerChSigma_ClusterIn = new TGraphErrors(AliMUONConstants::NTrackingCh());
    gResidualYPerChSigma_ClusterIn->SetName("gResidualYPerChSigma_ClusterIn");
    gResidualYPerChSigma_ClusterIn->SetTitle("cluster-track residual-Y per Ch: sigma (cluster in);chamber ID;#sigma_{Y} (cm)");
    gResidualYPerChSigma_ClusterIn->SetMarkerStyle(kFullDotLarge);
  }
  TGraphErrors* gResidualXPerChSigma_ClusterOut = static_cast<TGraphErrors*>(outFile->FindObjectAny("gResidualXPerChSigma_ClusterOut"));
  if (!gResidualXPerChSigma_ClusterOut) {
    gResidualXPerChSigma_ClusterOut = new TGraphErrors(AliMUONConstants::NTrackingCh());
    gResidualXPerChSigma_ClusterOut->SetName("gResidualXPerChSigma_ClusterOut");
    gResidualXPerChSigma_ClusterOut->SetTitle("cluster-track residual-X per Ch: sigma (cluster out);chamber ID;#sigma_{X} (cm)");
    gResidualXPerChSigma_ClusterOut->SetMarkerStyle(kFullDotLarge);
  }
  TGraphErrors* gResidualYPerChSigma_ClusterOut = static_cast<TGraphErrors*>(outFile->FindObjectAny("gResidualYPerChSigma_ClusterOut"));
  if (!gResidualYPerChSigma_ClusterOut) {
    gResidualYPerChSigma_ClusterOut = new TGraphErrors(AliMUONConstants::NTrackingCh());
    gResidualYPerChSigma_ClusterOut->SetName("gResidualYPerChSigma_ClusterOut");
    gResidualYPerChSigma_ClusterOut->SetTitle("cluster-track residual-Y per Ch: sigma (cluster out);chamber ID;#sigma_{Y} (cm)");
    gResidualYPerChSigma_ClusterOut->SetMarkerStyle(kFullDotLarge);
  }
  TGraphErrors* gResidualXPerChDispersion_ClusterOut = static_cast<TGraphErrors*>(outFile->FindObjectAny("gResidualXPerChDispersion_ClusterOut"));
  if (!gResidualXPerChDispersion_ClusterOut) {
    gResidualXPerChDispersion_ClusterOut = new TGraphErrors(AliMUONConstants::NTrackingCh());
    gResidualXPerChDispersion_ClusterOut->SetName("gResidualXPerChDispersion_ClusterOut");
    gResidualXPerChDispersion_ClusterOut->SetTitle("cluster-track residual-X per Ch: dispersion (cluster out);chamber ID;#sigma_{X} (cm)");
    gResidualXPerChDispersion_ClusterOut->SetMarkerStyle(kFullDotLarge);
  }
  TGraphErrors* gResidualYPerChDispersion_ClusterOut = static_cast<TGraphErrors*>(outFile->FindObjectAny("gResidualYPerChDispersion_ClusterOut"));
  if (!gResidualYPerChDispersion_ClusterOut) {
    gResidualYPerChDispersion_ClusterOut = new TGraphErrors(AliMUONConstants::NTrackingCh());
    gResidualYPerChDispersion_ClusterOut->SetName("gResidualYPerChDispersion_ClusterOut");
    gResidualYPerChDispersion_ClusterOut->SetTitle("cluster-track residual-Y per Ch: dispersion (cluster out);chamber ID;#sigma_{Y} (cm)");
    gResidualYPerChDispersion_ClusterOut->SetMarkerStyle(kFullDotLarge);
  }
  TGraphErrors* gCombinedResidualXPerChSigma = static_cast<TGraphErrors*>(outFile->FindObjectAny("gCombinedResidualXPerChSigma"));
  if (!gCombinedResidualXPerChSigma) {
    gCombinedResidualXPerChSigma = new TGraphErrors(AliMUONConstants::NTrackingCh());
    gCombinedResidualXPerChSigma->SetName("gCombinedResidualXPerChSigma");
    gCombinedResidualXPerChSigma->SetTitle("combined cluster-track residual-X per Ch: sigma;chamber ID;#sigma_{X} (cm)");
    gCombinedResidualXPerChSigma->SetMarkerStyle(kFullDotLarge);
  }
  TGraphErrors* gCombinedResidualYPerChSigma = static_cast<TGraphErrors*>(outFile->FindObjectAny("gCombinedResidualYPerChSigma"));
  if (!gCombinedResidualYPerChSigma) {
    gCombinedResidualYPerChSigma = new TGraphErrors(AliMUONConstants::NTrackingCh());
    gCombinedResidualYPerChSigma->SetName("gCombinedResidualYPerChSigma");
    gCombinedResidualYPerChSigma->SetTitle("combined cluster-track residual-Y per Ch: sigma;chamber ID;#sigma_{Y} (cm)");
    gCombinedResidualYPerChSigma->SetMarkerStyle(kFullDotLarge);
  }
  
  TMultiGraph* mgCombinedResidualXSigmaVsP = new TMultiGraph("mgCombinedResidualXSigmaVsP","cluster X-resolution versus momentum;p (GeV/c^{2});#sigma_{X} (cm)");
  TMultiGraph* mgCombinedResidualYSigmaVsP = new TMultiGraph("mgCombinedResidualYSigmaVsP","cluster Y-resolution versus momentum;p (GeV/c^{2});#sigma_{Y} (cm)");
  TGraphErrors* gCombinedResidualXSigmaVsP[10];
  TGraphErrors* gCombinedResidualYSigmaVsP[10];
  for (Int_t i = 0; i < AliMUONConstants::NTrackingCh(); i++) {
    gCombinedResidualXSigmaVsP[i] = new TGraphErrors(hResidualXInChVsP_ClusterIn[i]->GetNbinsX());
    gCombinedResidualXSigmaVsP[i]->SetName(Form("gResXVsP_ch%d",i+1));
    gCombinedResidualXSigmaVsP[i]->SetMarkerStyle(kFullDotMedium);
    gCombinedResidualXSigmaVsP[i]->SetMarkerColor(i+1+i/9);
    mgCombinedResidualXSigmaVsP->Add(gCombinedResidualXSigmaVsP[i],"p");
    
    gCombinedResidualYSigmaVsP[i] = new TGraphErrors(hResidualYInChVsP_ClusterIn[i]->GetNbinsX());
    gCombinedResidualYSigmaVsP[i]->SetName(Form("gResYVsP_ch%d",i+1));
    gCombinedResidualYSigmaVsP[i]->SetMarkerStyle(kFullDotMedium);
    gCombinedResidualYSigmaVsP[i]->SetMarkerColor(i+1+i/9);
    mgCombinedResidualYSigmaVsP->Add(gCombinedResidualYSigmaVsP[i],"p");
  }
  
  TH2F *hTrackResXPerCh = static_cast<TH2F*>(outFile->FindObjectAny("hTrackResXPerCh"));
  TH2F *hTrackResYPerCh = static_cast<TH2F*>(outFile->FindObjectAny("hTrackResYPerCh"));
  
  TH2F *hMCSXPerCh = static_cast<TH2F*>(outFile->FindObjectAny("hMCSXPerCh"));
  TH2F *hMCSYPerCh = static_cast<TH2F*>(outFile->FindObjectAny("hMCSYPerCh"));
  
  TGraphErrors* gTrackResXPerCh = static_cast<TGraphErrors*>(outFile->FindObjectAny("gTrackResXPerCh"));
  if (!gTrackResXPerCh) {
    gTrackResXPerCh = new TGraphErrors(AliMUONConstants::NTrackingCh());
    gTrackResXPerCh->SetName("gTrackResXPerCh");
    gTrackResXPerCh->SetTitle("track <#sigma_{X}> per Ch;chamber ID;<#sigma_{X}> (cm)");
    gTrackResXPerCh->SetMarkerStyle(kFullDotLarge);
  }
  TGraphErrors* gTrackResYPerCh = static_cast<TGraphErrors*>(outFile->FindObjectAny("gTrackResYPerCh"));
  if (!gTrackResYPerCh) {
    gTrackResYPerCh = new TGraphErrors(AliMUONConstants::NTrackingCh());
    gTrackResYPerCh->SetName("gTrackResYPerCh");
    gTrackResYPerCh->SetTitle("track <#sigma_{Y}> per Ch;chamber ID;<#sigma_{Y}> (cm)");
    gTrackResYPerCh->SetMarkerStyle(kFullDotLarge);
  }
  TGraphErrors* gMCSXPerCh = static_cast<TGraphErrors*>(outFile->FindObjectAny("gMCSXPerCh"));
  if (!gMCSXPerCh) {
    gMCSXPerCh = new TGraphErrors(AliMUONConstants::NTrackingCh());
    gMCSXPerCh->SetName("gMCSXPerCh");
    gMCSXPerCh->SetTitle("MCS X-dispersion of extrapolated track per Ch;chamber ID;<#sigma_{X}> (cm)");
    gMCSXPerCh->SetMarkerStyle(kFullDotLarge);
  }
  TGraphErrors* gMCSYPerCh = static_cast<TGraphErrors*>(outFile->FindObjectAny("gMCSYPerCh"));
  if (!gMCSYPerCh) {
    gMCSYPerCh = new TGraphErrors(AliMUONConstants::NTrackingCh());
    gMCSYPerCh->SetName("gMCSYPerCh");
    gMCSYPerCh->SetTitle("MCS Y-dispersion of extrapolated track per Ch;chamber ID;<#sigma_{Y}> (cm)");
    gMCSYPerCh->SetMarkerStyle(kFullDotLarge);
  }
  TGraphErrors* gClusterResXPerCh = static_cast<TGraphErrors*>(outFile->FindObjectAny("gClusterResXPerCh"));
  if (!gClusterResXPerCh) {
    gClusterResXPerCh = new TGraphErrors(AliMUONConstants::NTrackingCh());
    gClusterResXPerCh->SetName("gClusterResXPerCh");
    gClusterResXPerCh->SetTitle("cluster #sigma_{X} per Ch;chamber ID;#sigma_{X} (cm)");
    gClusterResXPerCh->SetMarkerStyle(kFullDotLarge);
  }
  TGraphErrors* gClusterResYPerCh = static_cast<TGraphErrors*>(outFile->FindObjectAny("gClusterResYPerCh"));
  if (!gClusterResYPerCh) {
    gClusterResYPerCh = new TGraphErrors(AliMUONConstants::NTrackingCh());
    gClusterResYPerCh->SetName("gClusterResYPerCh");
    gClusterResYPerCh->SetTitle("cluster #sigma_{Y} per Ch;chamber ID;#sigma_{Y} (cm)");
    gClusterResYPerCh->SetMarkerStyle(kFullDotLarge);
  }
  
  TH2F *hClusterRes2XPerCh = static_cast<TH2F*>(outFile->FindObjectAny("hClusterRes2XPerCh"));
  TH2F *hClusterRes2YPerCh = static_cast<TH2F*>(outFile->FindObjectAny("hClusterRes2YPerCh"));
  
  TGraphErrors* gCalcClusterResXPerCh = static_cast<TGraphErrors*>(outFile->FindObjectAny("gCalcClusterResXPerCh"));
  if (!gCalcClusterResXPerCh) {
    gCalcClusterResXPerCh = new TGraphErrors(AliMUONConstants::NTrackingCh());
    gCalcClusterResXPerCh->SetName("gCalcClusterResXPerCh");
    gCalcClusterResXPerCh->SetTitle("calculated cluster #sigma_{X} per Ch;chamber ID;#sigma_{X} (cm)");
    gCalcClusterResXPerCh->SetMarkerStyle(kFullDotLarge);
  }
  TGraphErrors* gCalcClusterResYPerCh = static_cast<TGraphErrors*>(outFile->FindObjectAny("gCalcClusterResYPerCh"));
  if (!gCalcClusterResYPerCh) {
    gCalcClusterResYPerCh = new TGraphErrors(AliMUONConstants::NTrackingCh());
    gCalcClusterResYPerCh->SetName("gCalcClusterResYPerCh");
    gCalcClusterResYPerCh->SetTitle("calculated cluster #sigma_{Y} per Ch;chamber ID;#sigma_{Y} (cm)");
    gCalcClusterResYPerCh->SetMarkerStyle(kFullDotLarge);
  }
  
  // compute residual mean and dispersion
  Double_t meanIn, meanInErr, meanOut, meanOutErr, sigmaIn, sigmaInErr, sigmaOut, sigmaOutErr;
  Double_t sigmaTrack, sigmaTrackErr, sigmaMCS, sigmaMCSErr, clusterRes, clusterResErr, sigmaCluster, sigmaClusterErr;
  for (Int_t i = 0; i < AliMUONConstants::NTrackingCh(); i++) {
    
    // non-bending direction
    GetMean(hResidualXInCh_ClusterIn[i], meanIn, meanInErr, gResidualXPerChMean_ClusterIn, i, i+1);
    GetRMS(hResidualXInCh_ClusterIn[i], sigmaIn, sigmaInErr, gResidualXPerChSigma_ClusterIn, i, i+1);
    
    GetMean(hResidualXInCh_ClusterOut[i], meanOut, meanOutErr, gResidualXPerChMean_ClusterOut, i, i+1);
    GetRMS(hResidualXInCh_ClusterOut[i], sigmaOut, sigmaOutErr, gResidualXPerChSigma_ClusterOut, i, i+1);
    
    if (correctForSystematics) {
      sigmaIn = TMath::Sqrt(sigmaIn*sigmaIn + meanIn*meanIn);
      sigmaInErr = TMath::Sqrt(sigmaIn*sigmaIn*sigmaInErr*sigmaInErr + meanIn*meanIn*meanInErr*meanInErr) / sigmaIn;
      sigmaOut = TMath::Sqrt(sigmaOut*sigmaOut + meanOut*meanOut);
      sigmaOutErr = TMath::Sqrt(sigmaOut*sigmaOut*sigmaOutErr*sigmaOutErr + meanOut*meanOut*meanOutErr*meanOutErr) / sigmaOut;
    }
    gResidualXPerChDispersion_ClusterOut->SetPoint(i, i+1, sigmaOut);
    gResidualXPerChDispersion_ClusterOut->SetPointError(i, 0., sigmaOutErr);
    
    TH1D *tmp = hTrackResXPerCh->ProjectionY("tmp",i+1,i+1,"e");
    GetMean(tmp, sigmaTrack, sigmaTrackErr, gTrackResXPerCh, i, i+1, kFALSE);
    delete tmp;
    
    tmp = hMCSXPerCh->ProjectionY("tmp",i+1,i+1,"e");
    GetMean(tmp, sigmaMCS, sigmaMCSErr, gMCSXPerCh, i, i+1, kFALSE);
    delete tmp;
    
    clusterRes = TMath::Sqrt(sigmaIn*sigmaOut);
    clusterResErr = 0.5 * TMath::Sqrt(sigmaInErr*sigmaInErr*sigmaOut*sigmaOut + sigmaIn*sigmaIn*sigmaOutErr*sigmaOutErr) / clusterRes;
    gCombinedResidualXPerChSigma->SetPoint(i, i+1, clusterRes);
    gCombinedResidualXPerChSigma->SetPointError(i, 0., clusterResErr);
    if (clusterResNB) clusterResNB[i] = clusterRes;
    if (clusterResNBErr) clusterResNBErr[i] = clusterResErr;
    
    sigmaCluster = sigmaOut*sigmaOut - sigmaTrack*sigmaTrack;
    if (sigmaCluster > 0.) {
      sigmaCluster = TMath::Sqrt(sigmaCluster);
      sigmaClusterErr = TMath::Sqrt(sigmaOut*sigmaOut*sigmaOutErr*sigmaOutErr + sigmaTrack*sigmaTrack*sigmaTrackErr*sigmaTrackErr) / sigmaCluster;
    } else {
      sigmaCluster = 0.;
      sigmaClusterErr = 0.;
    }
    gClusterResXPerCh->SetPoint(i, i+1, sigmaCluster);
    gClusterResXPerCh->SetPointError(i, 0., sigmaClusterErr);
    
    tmp = hClusterRes2XPerCh->ProjectionY("tmp",i+1,i+1,"e");
    ZoomRight(tmp);
    clusterRes = tmp->GetMean();
    if (clusterRes > 0.) {
      gCalcClusterResXPerCh->SetPoint(i, i+1, TMath::Sqrt(clusterRes));
      gCalcClusterResXPerCh->SetPointError(i, 0., 0.5 * tmp->GetMeanError() / TMath::Sqrt(clusterRes));
    } else {
      gCalcClusterResXPerCh->SetPoint(i, i+1, 0.);
      gCalcClusterResXPerCh->SetPointError(i, 0., 0.);
    }
    delete tmp;
    
    FillSigmaClusterVsP(hResidualXInChVsP_ClusterIn[i], hResidualXInChVsP_ClusterOut[i], gCombinedResidualXSigmaVsP[i]);
    
    // bending direction
    GetMean(hResidualYInCh_ClusterIn[i], meanIn, meanInErr, gResidualYPerChMean_ClusterIn, i, i+1);
    GetRMS(hResidualYInCh_ClusterIn[i], sigmaIn, sigmaInErr, gResidualYPerChSigma_ClusterIn, i, i+1);
    
    GetMean(hResidualYInCh_ClusterOut[i], meanOut, meanOutErr, gResidualYPerChMean_ClusterOut, i, i+1);
    GetRMS(hResidualYInCh_ClusterOut[i], sigmaOut, sigmaOutErr, gResidualYPerChSigma_ClusterOut, i, i+1);
    
    if (correctForSystematics) {
      sigmaIn = TMath::Sqrt(sigmaIn*sigmaIn + meanIn*meanIn);
      sigmaInErr = TMath::Sqrt(sigmaIn*sigmaIn*sigmaInErr*sigmaInErr + meanIn*meanIn*meanInErr*meanInErr) / sigmaIn;
      sigmaOut = TMath::Sqrt(sigmaOut*sigmaOut + meanOut*meanOut);
      sigmaOutErr = TMath::Sqrt(sigmaOut*sigmaOut*sigmaOutErr*sigmaOutErr + meanOut*meanOut*meanOutErr*meanOutErr) / sigmaOut;
    }
    gResidualYPerChDispersion_ClusterOut->SetPoint(i, i+1, sigmaOut);
    gResidualYPerChDispersion_ClusterOut->SetPointError(i, 0., sigmaOutErr);
    
    tmp = hTrackResYPerCh->ProjectionY("tmp",i+1,i+1,"e");
    GetMean(tmp, sigmaTrack, sigmaTrackErr, gTrackResYPerCh, i, i+1, kFALSE);
    delete tmp;
    
    tmp = hMCSYPerCh->ProjectionY("tmp",i+1,i+1,"e");
    GetMean(tmp, sigmaMCS, sigmaMCSErr, gMCSYPerCh, i, i+1, kFALSE);
    delete tmp;
    
    clusterRes = TMath::Sqrt(sigmaIn*sigmaOut);
    clusterResErr = 0.5 * TMath::Sqrt(sigmaInErr*sigmaInErr*sigmaOut*sigmaOut + sigmaIn*sigmaIn*sigmaOutErr*sigmaOutErr) / clusterRes;
    gCombinedResidualYPerChSigma->SetPoint(i, i+1, clusterRes);
    gCombinedResidualYPerChSigma->SetPointError(i, 0., clusterResErr);
    if (clusterResB) clusterResB[i] = clusterRes;
    if (clusterResBErr) clusterResBErr[i] = clusterResErr;
    
    sigmaCluster = sigmaOut*sigmaOut - sigmaTrack*sigmaTrack;
    if (sigmaCluster > 0.) {
      sigmaCluster = TMath::Sqrt(sigmaCluster);
      sigmaClusterErr = TMath::Sqrt(sigmaOut*sigmaOut*sigmaOutErr*sigmaOutErr + sigmaTrack*sigmaTrack*sigmaTrackErr*sigmaTrackErr) / sigmaCluster;
    } else {
      sigmaCluster = 0.;
      sigmaClusterErr = 0.;
    }
    gClusterResYPerCh->SetPoint(i, i+1, sigmaCluster);
    gClusterResYPerCh->SetPointError(i, 0., sigmaClusterErr);
    
    tmp = hClusterRes2YPerCh->ProjectionY("tmp",i+1,i+1,"e");
    ZoomRight(tmp);
    clusterRes = tmp->GetMean();
    if (clusterRes > 0.) {
      gCalcClusterResYPerCh->SetPoint(i, i+1, TMath::Sqrt(clusterRes));
      gCalcClusterResYPerCh->SetPointError(i, 0., 0.5 * tmp->GetMeanError() / TMath::Sqrt(clusterRes));
    } else {
      gCalcClusterResYPerCh->SetPoint(i, i+1, 0.);
      gCalcClusterResYPerCh->SetPointError(i, 0., 0.);
    }
    delete tmp;
    
    FillSigmaClusterVsP(hResidualYInChVsP_ClusterIn[i], hResidualYInChVsP_ClusterOut[i], gCombinedResidualYSigmaVsP[i]);
    
  }
  
  // display
  TCanvas* summary1 = new TCanvas("summary1","summary1");
  summary1->Divide(2,2);
  summary1->cd(1);
  gResidualXPerChMean_ClusterOut->Draw("ap");
  gResidualXPerChMean_ClusterOut->SetMarkerColor(2);
  gResidualXPerChMean_ClusterOut->SetLineColor(2);
  gResidualXPerChMean_ClusterIn->Draw("p");
  gResidualXPerChMean_ClusterIn->SetMarkerColor(4);
  gResidualXPerChMean_ClusterIn->SetLineColor(4);
  summary1->cd(2);
  gResidualYPerChMean_ClusterOut->Draw("ap");
  gResidualYPerChMean_ClusterOut->SetMarkerColor(2);
  gResidualYPerChMean_ClusterOut->SetLineColor(2);
  gResidualYPerChMean_ClusterIn->Draw("p");
  gResidualYPerChMean_ClusterIn->SetMarkerColor(4);
  gResidualYPerChMean_ClusterIn->SetLineColor(4);
  summary1->cd(3);
  gResidualXPerChSigma_ClusterOut->Draw("ap");
  gResidualXPerChSigma_ClusterOut->SetMinimum(0.);
  gResidualXPerChSigma_ClusterOut->SetMarkerColor(2);
  gResidualXPerChSigma_ClusterOut->SetLineColor(2);
  gResidualXPerChSigma_ClusterIn->Draw("p");
  gResidualXPerChSigma_ClusterIn->SetMarkerColor(4);
  gResidualXPerChSigma_ClusterIn->SetLineColor(4);
  gMCSXPerCh->Draw("p");
  gMCSXPerCh->SetMarkerColor(5);
  gMCSXPerCh->SetLineColor(5);
  gCombinedResidualXPerChSigma->Draw("p");
  gCombinedResidualXPerChSigma->SetMarkerColor(3);
  gCombinedResidualXPerChSigma->SetLineColor(3);
  summary1->cd(4);
  gResidualYPerChSigma_ClusterOut->Draw("ap");
  gResidualYPerChSigma_ClusterOut->SetMinimum(0.);
  gResidualYPerChSigma_ClusterOut->SetMarkerColor(2);
  gResidualYPerChSigma_ClusterOut->SetLineColor(2);
  gResidualYPerChSigma_ClusterIn->Draw("p");
  gResidualYPerChSigma_ClusterIn->SetMarkerColor(4);
  gResidualYPerChSigma_ClusterIn->SetLineColor(4);
  gMCSYPerCh->Draw("p");
  gMCSYPerCh->SetMarkerColor(5);
  gMCSYPerCh->SetLineColor(5);
  gCombinedResidualYPerChSigma->Draw("p");
  gCombinedResidualYPerChSigma->SetMarkerColor(3);
  gCombinedResidualYPerChSigma->SetLineColor(3);
  summary1->Update();
  
  TCanvas* summary2 = new TCanvas("summary2","summary2");
  summary2->Divide(2,2);
  summary2->cd(1);
  gResidualXPerChDispersion_ClusterOut->Draw("ap");
  gResidualXPerChDispersion_ClusterOut->SetMinimum(0.);
  gResidualXPerChDispersion_ClusterOut->SetMarkerColor(2);
  gResidualXPerChDispersion_ClusterOut->SetLineColor(2);
  gMCSXPerCh->Draw("p");
  gTrackResXPerCh->Draw("p");
  gTrackResXPerCh->SetMarkerColor(4);
  gTrackResXPerCh->SetLineColor(4);
  gClusterResXPerCh->Draw("p");
  summary2->cd(2);
  gResidualYPerChDispersion_ClusterOut->Draw("ap");
  gResidualYPerChDispersion_ClusterOut->SetMinimum(0.);
  gResidualYPerChDispersion_ClusterOut->SetMarkerColor(2);
  gResidualYPerChDispersion_ClusterOut->SetLineColor(2);
  gMCSYPerCh->Draw("p");
  gTrackResYPerCh->Draw("p");
  gTrackResYPerCh->SetMarkerColor(4);
  gTrackResYPerCh->SetLineColor(4);
  gClusterResYPerCh->Draw("p");
  summary2->cd(3);
  gCombinedResidualXPerChSigma->Draw("ap");
  gCombinedResidualXPerChSigma->SetMinimum(0.);
  gClusterResXPerCh->Draw("p");
  gCalcClusterResXPerCh->Draw("p");
  gCalcClusterResXPerCh->SetMarkerColor(6);
  gCalcClusterResXPerCh->SetLineColor(6);
  summary2->cd(4);
  gCombinedResidualYPerChSigma->Draw("ap");
  gCombinedResidualYPerChSigma->SetMinimum(0.);
  gClusterResYPerCh->Draw("p");
  gCalcClusterResYPerCh->Draw("p");
  gCalcClusterResYPerCh->SetMarkerColor(6);
  gCalcClusterResYPerCh->SetLineColor(6);
  summary2->Update();
  
  TCanvas* resVsP = new TCanvas("resVsP","resVsP");
  resVsP->Divide(1,2);
  resVsP->cd(1);
  mgCombinedResidualXSigmaVsP->Draw("ap");
  resVsP->cd(2);
  mgCombinedResidualYSigmaVsP->Draw("ap");
  resVsP->Update();
  
  // fill output
  outFile->cd();
  gResidualXPerChMean_ClusterIn->Write(0,TObject::kOverwrite);
  gResidualXPerChMean_ClusterOut->Write(0,TObject::kOverwrite);
  gResidualXPerChSigma_ClusterIn->Write(0,TObject::kOverwrite);
  gResidualXPerChSigma_ClusterOut->Write(0,TObject::kOverwrite);
  gResidualXPerChDispersion_ClusterOut->Write(0,TObject::kOverwrite);
  gCombinedResidualXPerChSigma->Write(0,TObject::kOverwrite);
  gTrackResXPerCh->Write(0,TObject::kOverwrite);
  gMCSXPerCh->Write(0,TObject::kOverwrite);
  gClusterResXPerCh->Write(0,TObject::kOverwrite);
  gCalcClusterResXPerCh->Write(0,TObject::kOverwrite);
  gResidualYPerChMean_ClusterIn->Write(0,TObject::kOverwrite);
  gResidualYPerChMean_ClusterOut->Write(0,TObject::kOverwrite);
  gResidualYPerChSigma_ClusterIn->Write(0,TObject::kOverwrite);
  gResidualYPerChSigma_ClusterOut->Write(0,TObject::kOverwrite);
  gResidualYPerChDispersion_ClusterOut->Write(0,TObject::kOverwrite);
  gCombinedResidualYPerChSigma->Write(0,TObject::kOverwrite);
  gTrackResYPerCh->Write(0,TObject::kOverwrite);
  gMCSYPerCh->Write(0,TObject::kOverwrite);
  gClusterResYPerCh->Write(0,TObject::kOverwrite);
  gCalcClusterResYPerCh->Write(0,TObject::kOverwrite);
  mgCombinedResidualXSigmaVsP->Write(0,TObject::kOverwrite);
  mgCombinedResidualYSigmaVsP->Write(0,TObject::kOverwrite);
  summary1->Write(0,TObject::kOverwrite);
  summary2->Write(0,TObject::kOverwrite);
  resVsP->Write(0,TObject::kOverwrite);
  outFile->Close();
  
}

//-----------------------------------------------------------------------
TTree* GetESDTree(TFile *esdFile)
{
  /// Check that the file is properly open
  /// Return pointer to the ESD Tree
  
  if (!esdFile || !esdFile->IsOpen()) {
    Error("GetESDTree", "opening ESD file failed");
    exit(-1);
  }
  
  TTree* tree = (TTree*) esdFile->Get("esdTree");
  if (!tree) {
    Error("GetESDTree", "no ESD tree found");
    exit(-1);
  }
  
  return tree;
  
}

//-----------------------------------------------------------------------
void SetClusterResolution(AliMUONTrack& track, Double_t clusterResNB[10], Double_t clusterResB[10])
{
  /// Reset the clusters resolution from the ones given in argument
  Int_t nClusters = track.GetNClusters();
  for (Int_t iCluster=0; iCluster<nClusters; iCluster++) {
    AliMUONVCluster* cl = static_cast<AliMUONTrackParam*>(track.GetTrackParamAtCluster()->UncheckedAt(iCluster))->GetClusterPtr();
    cl->SetErrXY(clusterResNB[cl->GetChamberId()],clusterResB[cl->GetChamberId()]);
  }
}

//-----------------------------------------------------------------------
void Zoom(TH1* h, Double_t fractionCut)
{
  /// Reduce the range of the histogram by removing a given fration of the statistic at each edge
  ZoomLeft(h, fractionCut);
  ZoomRight(h, fractionCut);
}

//-----------------------------------------------------------------------
void ZoomLeft(TH1* h, Double_t fractionCut)
{
  /// Reduce the range of the histogram by removing a given fration of the statistic on the left side
  Int_t maxEventsCut = fractionCut * h->GetEntries();
  Int_t nBins = h->GetNbinsX();
  
  // set low edge  
  Int_t minBin;
  Int_t eventsCut = 0;
  for (minBin = 1; minBin <= nBins; minBin++) {
    eventsCut += h->GetBinContent(minBin);
    if (eventsCut > maxEventsCut) break;
  }
  
  // set new axis range
  h->GetXaxis()->SetRange(minBin, h->GetXaxis()->GetLast());
}

//-----------------------------------------------------------------------
void ZoomRight(TH1* h, Double_t fractionCut)
{
  /// Reduce the range of the histogram by removing a given fration of the statistic on the right side
  Int_t maxEventsCut = fractionCut * h->GetEntries();
  Int_t nBins = h->GetNbinsX();
  
  // set high edge
  Int_t maxBin;
  Int_t eventsCut = 0;
  for (maxBin = nBins; maxBin >= 1; maxBin--) {
    eventsCut += h->GetBinContent(maxBin);
    if (eventsCut > maxEventsCut) break;
  }
  
  // set new axis range
  h->GetXaxis()->SetRange(h->GetXaxis()->GetFirst(), maxBin);
}

//-----------------------------------------------------------------------
void GetMean(TH1* h, Double_t& mean, Double_t& meanErr, TGraphErrors* g, Int_t i, Double_t x, Bool_t zoom)
{
  /// Fill graph with the mean value of the histogram and the corresponding error (zooming if required)
  Int_t firstBin = h->GetXaxis()->GetFirst();
  Int_t lastBin = h->GetXaxis()->GetLast();
  //if (zoom) h->GetXaxis()->SetRangeUser(-3.*h->GetRMS(), 3.*h->GetRMS());
  if (zoom) Zoom(h);
  mean = h->GetMean();
  meanErr = h->GetMeanError();
  if (g) {
    g->SetPoint(i, x, mean);
    g->SetPointError(i, 0., meanErr);
  }
  if (zoom) h->GetXaxis()->SetRange(firstBin,lastBin);
}

//-----------------------------------------------------------------------
void GetRMS(TH1* h, Double_t& rms, Double_t& rmsErr, TGraphErrors* g, Int_t i, Double_t x, Bool_t zoom)
{
  /// Return the RMS of the histogram and the corresponding error (zooming if required) and fill graph if !=0x0
  Int_t firstBin = h->GetXaxis()->GetFirst();
  Int_t lastBin = h->GetXaxis()->GetLast();
  //if (zoom) h->GetXaxis()->SetRangeUser(-3.*h->GetRMS(), 3.*h->GetRMS());
  if (zoom) Zoom(h);
  rms = h->GetRMS();
  rmsErr = h->GetRMSError();
  if (g) {
    g->SetPoint(i, x, rms);
    g->SetPointError(i, 0., rmsErr);
  }
  if (zoom) h->GetXaxis()->SetRange(firstBin,lastBin);
}

//-----------------------------------------------------------------------
void FillSigmaClusterVsP(TH2* hIn, TH2* hOut, TGraphErrors* g, Bool_t zoom)
{
  /// Fill graph with cluster resolution from combined residuals with cluster in/out (zooming if required)
  Double_t sigmaIn, sigmaInErr, sigmaOut, sigmaOutErr, clusterRes, clusterResErr;
  for (Int_t j = 1; j <= hIn->GetNbinsX(); j++) {
    TH1D* tmp = hIn->ProjectionY("tmp",j,j,"e");
    GetRMS(tmp, sigmaIn, sigmaInErr, 0x0, 0, 0., zoom);
    delete tmp;
    tmp = hOut->ProjectionY("tmp",j,j,"e");
    GetRMS(tmp, sigmaOut, sigmaOutErr, 0x0, 0, 0., zoom);
    delete tmp;
    Double_t p = 0.5 * (hIn->GetBinLowEdge(j) + hIn->GetBinLowEdge(j+1));
    Double_t pErr = p - hIn->GetBinLowEdge(j);
    clusterRes = TMath::Sqrt(sigmaIn*sigmaOut);
    clusterResErr = 0.5 * TMath::Sqrt(sigmaInErr*sigmaInErr*sigmaOut*sigmaOut + sigmaIn*sigmaIn*sigmaOutErr*sigmaOutErr) / clusterRes;
    g->SetPoint(j, p, clusterRes);
    g->SetPointError(j, pErr, clusterResErr);
  }
}

