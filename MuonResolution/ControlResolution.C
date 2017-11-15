/*
 *  ControlResolution.C
 *
 *  Created by Philippe Pillot on 14/11/17.
 *  Copyright 2017 SUBATECH
 *  This software is made available under the terms of the GNU GPL 3.0
 *
 */

#if !defined(__CINT__) || defined(__MAKECINT__)

#include <TROOT.h>
#include <TMath.h>
#include <TFile.h>
#include <TObjArray.h>
#include <TString.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TGraphAsymmErrors.h>
#include <TF1.h>
#include <TCanvas.h>
#include "AliMuonTrackSmearing.h"

#endif

void DrawExpectedRes(TCanvas *cRes, Bool_t atFirstCluster);
void DrawResolutionFromData(TCanvas *cRes, TCanvas *cResX, TCanvas *cResY, TString fileResMeas, Bool_t atFirstCluster);
void DrawResolutionFromMC(TCanvas *cRes, TCanvas *cResX, TCanvas *cResY, TString fileResSim,
                          Bool_t atFirstCluster, Bool_t refitResVsP, Bool_t drawResPVsP, Int_t rebinP);
void DrawResolutionFromSmearing(TCanvas *cRes, TString fileResQA, Bool_t atFirstCluster, Bool_t drawResPVsP, Int_t rebinP);
Double_t GetSigma(TH1 *h, Double_t *sigmaErr = 0x0);
Double_t GetSigmaCrystalBall(TH1 *h, Double_t *sigmaErr = 0x0);
Double_t GetSigmaGaus(TH1 *h, Double_t *sigmaErr = 0x0);
Double_t GetSigmaBreitWigner(TH1 *h, Double_t *sigmaErr = 0x0);
void FitGausResVsMom(TH2 *h, TGraphAsymmErrors *gSigma, Int_t rebinP);
Double_t langaufun(Double_t *x, Double_t *par);
Double_t GausInvXFun(Double_t *x, Double_t *par);
void FitPResVsP(TH2 *h, TGraphAsymmErrors *gSigma, Int_t rebinP, Double_t pSwitch, Bool_t drawResPVsP);

// smearing class with settings used in MC
AliMuonTrackSmearing muonTrackSmearing(AliMuonTrackSmearing::kCrystalBall);

//-----------------------------------------------------------------------
void ControlResolution(TString fileResQA = "AnalysisResults.root", TString fileResMeas = "", TString fileResSim = "",
                       Bool_t atFirstCluster = kFALSE, Bool_t refitResVsP = kTRUE, Bool_t drawResPVsP = kFALSE,
                       Int_t rebinP = 2)
{
  /// Extract the resolutions vs p from the results of the smearing QA task
  /// and compare them to the expected shape from the smearing settings
  /// and to the results from the resolution and performance tasks if provided.
  /// If refitResVsP = kTRUE: recompute the resolution vs p from the results of
  /// the performance task with the same function as for the smeared results.
  
  setbuf(stdout, NULL);
  /*
  muonTrackSmearing.SetSigmaxChSt1(0.000519);
  muonTrackSmearing.SetSigmayChSt1(0.000064);
  muonTrackSmearing.SetSigmayCh(0.000105);
  muonTrackSmearing.SetCrystalBallParams(3.671484,4.590817,2.453060,1.715218,2.145035,1.782320);
  */
  // momentum and slope resolution versus p
  TCanvas *cRes = new TCanvas("cRes","cRes",10,10,1200,300);
  cRes->Divide(3,1);
  
  // cluster-track and cluster-trackRef residuals
  TCanvas *cResX = 0x0, *cResY = 0x0;
  if (!fileResMeas.IsNull() || !fileResSim.IsNull()) {
    cResX = new TCanvas("cResX","cResX",10,10,1200,600);
    cResX->Divide(6,3);
    cResY = new TCanvas("cResY","cResY",10,10,1200,600);
    cResY->Divide(6,3);
  }
  
  DrawResolutionFromData(cRes, cResX, cResY, fileResMeas, atFirstCluster);
  DrawResolutionFromMC(cRes, cResX, cResY, fileResSim, atFirstCluster, refitResVsP, drawResPVsP, rebinP);
  DrawExpectedRes(cRes, atFirstCluster);
  DrawResolutionFromSmearing(cRes, fileResQA, atFirstCluster, drawResPVsP, rebinP);

}

//-----------------------------------------------------------------------
void DrawExpectedRes(TCanvas *cRes, Bool_t atFirstCluster)
{
  /// Draw the expected resolutions vs p from the smearing settings
  
  // momentum relative resolution versus p
  Double_t where = atFirstCluster ? 1. : -1.;
  TF1 *fPResVsP02 = new TF1("fPResVsP02", &muonTrackSmearing, &AliMuonTrackSmearing::PResVsP, 0.,1000.,2, "AliMuonTrackSmearing", "PResVsP");
  fPResVsP02->SetParameters(1.,where);
  fPResVsP02->SetNpx(10000);
  fPResVsP02->SetLineColor(15);
  TF1 *fPResVsP23 = new TF1("fPResVsP23", &muonTrackSmearing, &AliMuonTrackSmearing::PResVsP, 0.,1000.,2, "AliMuonTrackSmearing", "PResVsP");
  fPResVsP23->SetParameters(2.5,where);
  fPResVsP23->SetNpx(10000);
  fPResVsP23->SetLineColor(2);
  TF1 *fPResVsP310 = new TF1("fPResVsP310", &muonTrackSmearing, &AliMuonTrackSmearing::PResVsP, 0.,1000.,2, "AliMuonTrackSmearing", "PResVsP");
  fPResVsP310->SetParameters(6.,where);
  fPResVsP310->SetNpx(10000);
  fPResVsP310->SetLineColor(6);
  
  // slope resolution versus p
  TF1 *fSlopeXResVsP02 = new TF1("fSlopeXResVsP02", &muonTrackSmearing, &AliMuonTrackSmearing::SlopeResVsP, 0.,1000.,3, "AliMuonTrackSmearing", "SlopeResVsP");
  fSlopeXResVsP02->SetParameters(1.5,-1.,where);
  fSlopeXResVsP02->SetNpx(10000);
  fSlopeXResVsP02->SetLineColor(15);
  TF1 *fSlopeXResVsP23 = new TF1("fSlopeXResVsP23", &muonTrackSmearing, &AliMuonTrackSmearing::SlopeResVsP, 0.,1000.,3, "AliMuonTrackSmearing", "SlopeResVsP");
  fSlopeXResVsP23->SetParameters(2.5,-1.,where);
  fSlopeXResVsP23->SetNpx(10000);
  fSlopeXResVsP23->SetLineColor(2);
  TF1 *fSlopeXResVsP310 = new TF1("fSlopeXResVsP310", &muonTrackSmearing, &AliMuonTrackSmearing::SlopeResVsP, 0.,1000.,3, "AliMuonTrackSmearing", "SlopeResVsP");
  fSlopeXResVsP310->SetParameters(6.,-1.,where);
  fSlopeXResVsP310->SetNpx(10000);
  fSlopeXResVsP310->SetLineColor(6);
  TF1 *fSlopeYResVsP02 = new TF1("fSlopeYResVsP02", &muonTrackSmearing, &AliMuonTrackSmearing::SlopeResVsP, 0.,1000.,3, "AliMuonTrackSmearing", "SlopeResVsP");
  fSlopeYResVsP02->SetParameters(1.5,1.,where);
  fSlopeYResVsP02->SetNpx(10000);
  fSlopeYResVsP02->SetLineColor(15);
  TF1 *fSlopeYResVsP23 = new TF1("fSlopeYResVsP23", &muonTrackSmearing, &AliMuonTrackSmearing::SlopeResVsP, 0.,1000.,3, "AliMuonTrackSmearing", "SlopeResVsP");
  fSlopeYResVsP23->SetParameters(2.5,1.,where);
  fSlopeYResVsP23->SetNpx(10000);
  fSlopeYResVsP23->SetLineColor(2);
  TF1 *fSlopeYResVsP310 = new TF1("fSlopeYResVsP310", &muonTrackSmearing, &AliMuonTrackSmearing::SlopeResVsP, 0.,1000.,3, "AliMuonTrackSmearing", "SlopeResVsP");
  fSlopeYResVsP310->SetParameters(6.,1.,where);
  fSlopeYResVsP310->SetNpx(10000);
  fSlopeYResVsP310->SetLineColor(6);
  
  // draw momentum and slope resolution versus p
  gROOT->SetSelectedPad(cRes->cd(1));
  gPad->SetLogz();
  if (gPad->GetListOfPrimitives()->GetEntries() > 0) {
    fPResVsP02->Draw("same");
  } else {
    fPResVsP02->GetYaxis()->SetRangeUser(0., 10.);
    fPResVsP02->Draw();
  }
  fPResVsP23->Draw("same");
  fPResVsP310->Draw("same");
  gROOT->SetSelectedPad(cRes->cd(2));
  gPad->SetLogz();
  if (gPad->GetListOfPrimitives()->GetEntries() > 0)
    fSlopeXResVsP02->Draw("same");
  else {
    fSlopeXResVsP02->GetYaxis()->SetRangeUser(0., 0.02);
    fSlopeXResVsP02->Draw();
  }
  fSlopeXResVsP23->Draw("same");
  fSlopeXResVsP310->Draw("same");
  gROOT->SetSelectedPad(cRes->cd(3));
  gPad->SetLogz();
  if (gPad->GetListOfPrimitives()->GetEntries() > 0)
    fSlopeYResVsP02->Draw("same");
  else {
    fSlopeYResVsP02->GetYaxis()->SetRangeUser(0., 0.02);
    fSlopeYResVsP02->Draw();
  }
  fSlopeYResVsP23->Draw("same");
  fSlopeYResVsP310->Draw("same");
  
}

//-----------------------------------------------------------------------
void DrawResolutionFromData(TCanvas *cRes, TCanvas *cResX, TCanvas *cResY, TString fileResMeas, Bool_t atFirstCluster)
{
  /// Draw the resolutions vs p and cluster resolution from the task AliAnalysisTaskMuonResolution
  
  if (fileResMeas.IsNull()) return;
  
  TCanvas *dummy = new TCanvas("dummy","dummy",10,10,68,1);
  
  TH2F *hPRes = 0x0, *hSlopeXRes = 0x0, *hSlopeYRes = 0x0;
  TH1D *hXResSt_clIn[5], *hXRes_clIn = 0x0, *hXResSt_clOut[5], *hXRes_clOut = 0x0;
  TH1D *hYResSt_clIn[5], *hYRes_clIn = 0x0, *hYResSt_clOut[5], *hYRes_clOut = 0x0;
  for (Int_t iSt = 0; iSt < 5; ++iSt) {
    hXResSt_clIn[iSt] = 0x0;
    hXResSt_clOut[iSt] = 0x0;
    hYResSt_clIn[iSt] = 0x0;
    hYResSt_clOut[iSt] = 0x0;
  }
  
  TFile *fResMeas = TFile::Open(fileResMeas.Data(),"READ");
  if (!fResMeas || !fResMeas->IsOpen()) return;
  
  // get the momentum and slope measured resolution to control their parameterization
  TObjArray* trackResList = static_cast<TObjArray*>(fResMeas->FindObjectAny("TrackRes"));
  if (!trackResList) return;
  if (atFirstCluster) {
    hPRes = static_cast<TH2F*>(trackResList->FindObject("hUncorrPRes"));
    hSlopeXRes = static_cast<TH2F*>(trackResList->FindObject("hUncorrSlopeXRes"));
    hSlopeYRes = static_cast<TH2F*>(trackResList->FindObject("hUncorrSlopeYRes"));
  } else {
    hPRes = static_cast<TH2F*>(trackResList->FindObject("hPRes"));
    hSlopeXRes = static_cast<TH2F*>(trackResList->FindObject("hSlopeXRes"));
    hSlopeYRes = static_cast<TH2F*>(trackResList->FindObject("hSlopeYRes"));
  }
  if (!hPRes || !hSlopeXRes || !hSlopeYRes) return;
  hPRes->SetDirectory(0);
  hSlopeXRes->SetDirectory(0);
  hSlopeYRes->SetDirectory(0);
  
  // get the cluster-track residuals to control their parameterization
  TObjArray* residualsList = static_cast<TObjArray*>(fResMeas->FindObjectAny("Residuals"));
  if (!residualsList) return;
  TH2F *hResidualXPerCh_ClusterIn = static_cast<TH2F*>(residualsList->FindObject("hResidualXPerCh_ClusterIn"));
  TH2F *hResidualXPerCh_ClusterOut = static_cast<TH2F*>(residualsList->FindObject("hResidualXPerCh_ClusterOut"));
  TH2F *hResidualYPerCh_ClusterIn = static_cast<TH2F*>(residualsList->FindObject("hResidualYPerCh_ClusterIn"));
  TH2F *hResidualYPerCh_ClusterOut = static_cast<TH2F*>(residualsList->FindObject("hResidualYPerCh_ClusterOut"));
  if (!hResidualXPerCh_ClusterIn || !hResidualXPerCh_ClusterOut || !hResidualYPerCh_ClusterIn || !hResidualYPerCh_ClusterOut) return;
  printf("\nmeasured resolution (m):\n");
  dummy->cd();
  for (Int_t iSt = 0; iSt < 5; ++iSt) {
    hXResSt_clIn[iSt] = hResidualXPerCh_ClusterIn->ProjectionY(Form("hXResSt%d_clIn",iSt+1),2*iSt+1,2*iSt+2,"e");
    hXResSt_clIn[iSt]->SetDirectory(0);
    hXResSt_clOut[iSt] = hResidualXPerCh_ClusterOut->ProjectionY(Form("hXResSt%d_clOut",iSt+1),2*iSt+1,2*iSt+2,"e");
    hXResSt_clOut[iSt]->SetDirectory(0);
    printf("- sigma x st%d = %.6f\n", iSt+1, 0.01*TMath::Sqrt(GetSigma(hXResSt_clIn[iSt])*GetSigma(hXResSt_clOut[iSt])));
    hYResSt_clIn[iSt] = hResidualYPerCh_ClusterIn->ProjectionY(Form("hYResSt%d_clIn",iSt+1),2*iSt+1,2*iSt+2,"e");
    hYResSt_clIn[iSt]->SetDirectory(0);
    hYResSt_clOut[iSt] = hResidualYPerCh_ClusterOut->ProjectionY(Form("hYResSt%d_clOut",iSt+1),2*iSt+1,2*iSt+2,"e");
    hYResSt_clOut[iSt]->SetDirectory(0);
    printf("- sigma y st%d = %.6f\n", iSt+1, 0.01*TMath::Sqrt(GetSigma(hYResSt_clIn[iSt])*GetSigma(hYResSt_clOut[iSt])));
  }
  hXRes_clIn = hResidualXPerCh_ClusterIn->ProjectionY("hXRes_clIn",1,10,"e");
  hXRes_clIn->SetDirectory(0);
  hXRes_clOut = hResidualXPerCh_ClusterOut->ProjectionY("hXRes_clOut",1,10,"e");
  hXRes_clOut->SetDirectory(0);
  printf("- sigma x = %.6f\n", 0.01*TMath::Sqrt(GetSigma(hXRes_clIn)*GetSigma(hXRes_clOut)));
  hYRes_clIn = hResidualYPerCh_ClusterIn->ProjectionY("hYRes_clIn",1,10,"e");
  hYRes_clIn->SetDirectory(0);
  hYRes_clOut = hResidualYPerCh_ClusterOut->ProjectionY("hYRes_clOut",1,10,"e");
  hYRes_clOut->SetDirectory(0);
  printf("- sigma y = %.6f\n\n", 0.01*TMath::Sqrt(GetSigma(hYRes_clIn)*GetSigma(hYRes_clOut)));
  
  fResMeas->Close();
  
  // draw momentum and slope resolution versus p
  gROOT->SetSelectedPad(cRes->cd(1));
  gPad->SetLogz();
  hPRes->Draw((gPad->GetListOfPrimitives()->GetEntries() > 0) ? "samecolz" : "colz");
  gROOT->SetSelectedPad(cRes->cd(2));
  gPad->SetLogz();
  hSlopeXRes->Draw((gPad->GetListOfPrimitives()->GetEntries() > 0) ? "samecolz" : "colz");
  gROOT->SetSelectedPad(cRes->cd(3));
  gPad->SetLogz();
  hSlopeYRes->Draw((gPad->GetListOfPrimitives()->GetEntries() > 0) ? "samecolz" : "colz");
  
  // draw cluster-track x-residuals
  for (Int_t iSt = 0; iSt < 5; ++iSt) {
    gROOT->SetSelectedPad(cResX->cd(iSt+1));
    gPad->SetLogy();
    hXResSt_clIn[iSt]->Draw();
    gROOT->SetSelectedPad(cResX->cd(iSt+7));
    gPad->SetLogy();
    hXResSt_clOut[iSt]->Draw();
  }
  gROOT->SetSelectedPad(cResX->cd(6));
  gPad->SetLogy();
  hXRes_clIn->Draw();
  gROOT->SetSelectedPad(cResX->cd(12));
  gPad->SetLogy();
  hXRes_clOut->Draw();
  
  // draw cluster-track y-residuals
  for (Int_t iSt = 0; iSt < 5; ++iSt) {
    gROOT->SetSelectedPad(cResY->cd(iSt+1));
    gPad->SetLogy();
    hYResSt_clIn[iSt]->Draw();
    gROOT->SetSelectedPad(cResY->cd(iSt+7));
    gPad->SetLogy();
    hYResSt_clOut[iSt]->Draw();
  }
  gROOT->SetSelectedPad(cResY->cd(6));
  gPad->SetLogy();
  hYRes_clIn->Draw();
  gROOT->SetSelectedPad(cResY->cd(12));
  gPad->SetLogy();
  hYRes_clOut->Draw();
  
}

//-----------------------------------------------------------------------
void DrawResolutionFromMC(TCanvas *cRes, TCanvas *cResX, TCanvas *cResY, TString fileResSim,
                          Bool_t atFirstCluster, Bool_t refitResVsP, Bool_t drawResPVsP, Int_t rebinP)
{
  /// Draw the resolutions vs p and cluster resolution from the task AliAnalysisTaskMuonPerformance

  if (fileResSim.IsNull()) return;
  
  TCanvas *dummy = new TCanvas("dummy","dummy",10,10,68,1);
  
  // read ouput of the task AliAnalysisTaskMuonPerformance
  TGraphAsymmErrors *gSigmaResPVsP = 0x0, *gSigmaResSlopeXVsP = 0x0, *gSigmaResSlopeYVsP = 0x0;
  TGraphAsymmErrors *gSigmaResPVsP2 = 0x0, *gSigmaResSlopeXVsP2 = 0x0, *gSigmaResSlopeYVsP2 = 0x0;
  TGraphAsymmErrors *gSigmaResPAtVtxVsPIn02deg = 0x0, *gSigmaResPAtVtxVsPIn23deg = 0x0, *gSigmaResPAtVtxVsPIn310deg = 0x0;
  TH1D *hXResSt[5], *hXRes = 0x0;
  TH1D *hYResSt[5], *hYRes = 0x0;
  for (Int_t iSt = 0; iSt < 5; ++iSt) {
    hXResSt[iSt] = 0x0;
    hYResSt[iSt] = 0x0;
  }
  
  TFile *fResSim = TFile::Open(fileResSim.Data(),"READ");
  if (!fResSim || !fResSim->IsOpen()) return;
  
  TObjArray* trackerResolutionList = static_cast<TObjArray*>(fResSim->FindObjectAny("TrackerResolution"));
  if (!trackerResolutionList) return;
  
  if (!refitResVsP) {
    
    // get the momentum and slope resolution from simulation
    if (atFirstCluster) {
      TObjArray* momentumAtFirstClList = static_cast<TObjArray*>(fResSim->FindObjectAny("MomentumAtFirstCl"));
      TObjArray* slopeAtFirstClList = static_cast<TObjArray*>(fResSim->FindObjectAny("SlopeAtFirstCl"));
      if (!momentumAtFirstClList || !slopeAtFirstClList) return;
      gSigmaResPVsP = static_cast<TGraphAsymmErrors*>(momentumAtFirstClList->FindObject("gSigmaResPAt1stClVsP"));
      gSigmaResSlopeXVsP = static_cast<TGraphAsymmErrors*>(slopeAtFirstClList->FindObject("gSigmaResSlopeXAt1stClVsP"));
      gSigmaResSlopeYVsP = static_cast<TGraphAsymmErrors*>(slopeAtFirstClList->FindObject("gSigmaResSlopeYAt1stClVsP"));
    } else {
      TObjArray* momentumAtVtxList = static_cast<TObjArray*>(fResSim->FindObjectAny("MomentumAtVtx"));
      TObjArray* slopeAtVtxList = static_cast<TObjArray*>(fResSim->FindObjectAny("SlopeAtVtx"));
      if (!momentumAtVtxList || !slopeAtVtxList) return;
      gSigmaResPVsP = static_cast<TGraphAsymmErrors*>(momentumAtVtxList->FindObject("gSigmaResPAtVtxVsP"));
      gSigmaResSlopeXVsP = static_cast<TGraphAsymmErrors*>(slopeAtVtxList->FindObject("gSigmaResSlopeXAtVtxVsP"));
      gSigmaResSlopeYVsP = static_cast<TGraphAsymmErrors*>(slopeAtVtxList->FindObject("gSigmaResSlopeYAtVtxVsP"));
    }
    if (!gSigmaResPVsP || !gSigmaResSlopeXVsP || !gSigmaResSlopeYVsP) return;
    
  } else {
    
    if (atFirstCluster) {
      
      // recompute momentum resolution at 1st cluster versus p
      TH2F *hResPAt1stClVsP = static_cast<TH2F*>(trackerResolutionList->FindObject("hResPAt1stClVsP"));
      if (!hResPAt1stClVsP) return;
      dummy->cd();
      gSigmaResPVsP2 = new TGraphAsymmErrors(hResPAt1stClVsP->GetNbinsX()/rebinP);
      gSigmaResPVsP2->SetName("gSigmaResPAt1stClVsP2");
      gSigmaResPVsP2->SetTitle("#sigma_{p}/p at 1st cluster versus p;p (GeV/c);#sigma_{p}/p (%)");
      gSigmaResPVsP2->SetLineColor(6);
      FitPResVsP(hResPAt1stClVsP, gSigmaResPVsP2, rebinP, -1., drawResPVsP);
      
      // recompute slope resolution at 1st cluster versus p
      TH2F *hResSlopeXAt1stClVsP = static_cast<TH2F*>(trackerResolutionList->FindObject("hResSlopeXAt1stClVsP"));
      TH2F *hResSlopeYAt1stClVsP = static_cast<TH2F*>(trackerResolutionList->FindObject("hResSlopeYAt1stClVsP"));
      if (!hResSlopeXAt1stClVsP || !hResSlopeYAt1stClVsP) return;
      dummy->cd();
      gSigmaResSlopeXVsP2 = new TGraphAsymmErrors(hResSlopeXAt1stClVsP->GetNbinsX()/rebinP);
      gSigmaResSlopeXVsP2->SetName("gSigmaResSlopeXAt1stClVsP2");
      gSigmaResSlopeXVsP2->SetTitle("#sigma_{slope_{X}} at 1st cluster versus p;p (GeV/c);#sigma_{slope_{X}}");
      gSigmaResSlopeXVsP2->SetLineColor(2);
      FitGausResVsMom(hResSlopeXAt1stClVsP, gSigmaResSlopeXVsP2, rebinP);
      gSigmaResSlopeYVsP2 = new TGraphAsymmErrors(hResSlopeYAt1stClVsP->GetNbinsX()/rebinP);
      gSigmaResSlopeYVsP2->SetName("gSigmaResSlopeYAt1stClVsP2");
      gSigmaResSlopeYVsP2->SetTitle("#sigma_{slope_{Y}} at 1st cluster versus p;p (GeV/c);#sigma_{slope_{Y}}");
      gSigmaResSlopeYVsP2->SetLineColor(2);
      FitGausResVsMom(hResSlopeYAt1stClVsP, gSigmaResSlopeYVsP2, rebinP);
      
    } else {
      
      // recompute momentum resolution at vertex versus p
      TH2F *hResPAtVtxVsPIn02deg = static_cast<TH2F*>(trackerResolutionList->FindObject("hResPAtVtxVsPIn02degMC"));
      TH2F *hResPAtVtxVsPIn23deg = static_cast<TH2F*>(trackerResolutionList->FindObject("hResPAtVtxVsPIn23deg"));
      TH2F *hResPAtVtxVsPIn310deg = static_cast<TH2F*>(trackerResolutionList->FindObject("hResPAtVtxVsPIn310deg"));
      if (!hResPAtVtxVsPIn02deg || !hResPAtVtxVsPIn23deg || !hResPAtVtxVsPIn310deg) return;
      dummy->cd();
      gSigmaResPAtVtxVsPIn02deg = new TGraphAsymmErrors(hResPAtVtxVsPIn02deg->GetNbinsX()/rebinP);
      gSigmaResPAtVtxVsPIn02deg->SetName("gSigmaResPAtVtxVsPIn02deg");
      gSigmaResPAtVtxVsPIn02deg->SetTitle("#sigma_{p}/p at vertex versus p in [0,2[ deg MC;p (GeV/c);#sigma_{p}/p (%)");
      gSigmaResPAtVtxVsPIn02deg->SetLineColor(15);
      FitPResVsP(hResPAtVtxVsPIn02deg, gSigmaResPAtVtxVsPIn02deg, rebinP, 300., drawResPVsP);
      gSigmaResPAtVtxVsPIn23deg = new TGraphAsymmErrors(hResPAtVtxVsPIn23deg->GetNbinsX()/rebinP);
      gSigmaResPAtVtxVsPIn23deg->SetName("gSigmaResPAtVtxVsPIn23deg");
      gSigmaResPAtVtxVsPIn23deg->SetTitle("#sigma_{p}/p at vertex versus p in [2,3[ deg;p (GeV/c);#sigma_{p}/p (%)");
      gSigmaResPAtVtxVsPIn23deg->SetLineColor(2);
      FitPResVsP(hResPAtVtxVsPIn23deg, gSigmaResPAtVtxVsPIn23deg, rebinP, 220., drawResPVsP);
      gSigmaResPAtVtxVsPIn310deg = new TGraphAsymmErrors(hResPAtVtxVsPIn310deg->GetNbinsX()/rebinP);
      gSigmaResPAtVtxVsPIn310deg->SetName("gSigmaResPAtVtxVsPIn310deg");
      gSigmaResPAtVtxVsPIn310deg->SetTitle("#sigma_{p}/p at vertex versus p in [3,10[ deg;p (GeV/c);#sigma_{p}/p (%)");
      gSigmaResPAtVtxVsPIn310deg->SetLineColor(6);
      FitPResVsP(hResPAtVtxVsPIn310deg, gSigmaResPAtVtxVsPIn310deg, rebinP, 160., drawResPVsP);
      
      // recompute slope resolution at vertex versus p
      TH2F *hResSlopeXAtVtxVsP = static_cast<TH2F*>(trackerResolutionList->FindObject("hResSlopeXAtVtxVsP"));
      TH2F *hResSlopeYAtVtxVsP = static_cast<TH2F*>(trackerResolutionList->FindObject("hResSlopeYAtVtxVsP"));
      if (!hResSlopeXAtVtxVsP || !hResSlopeYAtVtxVsP) return;
      dummy->cd();
      gSigmaResSlopeXVsP2 = new TGraphAsymmErrors(hResSlopeXAtVtxVsP->GetNbinsX()/rebinP);
      gSigmaResSlopeXVsP2->SetName("gSigmaResSlopeXAtVtxVsP2");
      gSigmaResSlopeXVsP2->SetTitle("#sigma_{slope_{X}} at vertex versus p;p (GeV/c);#sigma_{slope_{X}}");
      gSigmaResSlopeXVsP2->SetLineColor(2);
      FitGausResVsMom(hResSlopeXAtVtxVsP, gSigmaResSlopeXVsP2, rebinP);
      gSigmaResSlopeYVsP2 = new TGraphAsymmErrors(hResSlopeYAtVtxVsP->GetNbinsX()/rebinP);
      gSigmaResSlopeYVsP2->SetName("gSigmaResSlopeYAtVtxVsP2");
      gSigmaResSlopeYVsP2->SetTitle("#sigma_{slope_{Y}} at vertex versus p;p (GeV/c);#sigma_{slope_{Y}}");
      gSigmaResSlopeYVsP2->SetLineColor(2);
      FitGausResVsMom(hResSlopeYAtVtxVsP, gSigmaResSlopeYVsP2, rebinP);
      
    }
    
  }
  
  // get the cluster-trackRef residuals to control their parameterization
  TH2F *hResClXVsCh = static_cast<TH2F*>(trackerResolutionList->FindObject("hResClXVsCh"));
  TH2F *hResClYVsCh = static_cast<TH2F*>(trackerResolutionList->FindObject("hResClYVsCh"));
  if (!hResClXVsCh || !hResClYVsCh) return;
  printf("true resolution (m):\n");
  dummy->cd();
  for (Int_t iSt = 0; iSt < 5; ++iSt) {
    hXResSt[iSt] = hResClXVsCh->ProjectionY(Form("hXResSt%d",iSt+1),2*iSt+1,2*iSt+2,"e");
    hXResSt[iSt]->SetDirectory(0);
    printf("- sigma x st%d = %.6f\n", iSt+1, 0.01*GetSigma(hXResSt[iSt]));
    hYResSt[iSt] = hResClYVsCh->ProjectionY(Form("hYResSt%d",iSt+1),2*iSt+1,2*iSt+2,"e");
    hYResSt[iSt]->SetDirectory(0);
    printf("- sigma y st%d = %.6f\n", iSt+1, 0.01*GetSigma(hYResSt[iSt]));
  }
  hXRes = hResClXVsCh->ProjectionY("hXRes",1,10,"e");
  hXRes->SetDirectory(0);
  printf("- sigma x = %.6f\n", 0.01*GetSigma(hXRes));
  hYRes = hResClYVsCh->ProjectionY("hYRes",1,10,"e");
  hYRes->SetDirectory(0);
  printf("- sigma y = %.6f\n\n", 0.01*GetSigma(hYRes));
  
  fResSim->Close();
  
  // draw momentum and slope resolution versus p
  gROOT->SetSelectedPad(cRes->cd(1));
  gPad->SetLogz();
  if (gSigmaResPVsP) gSigmaResPVsP->Draw((gPad->GetListOfPrimitives()->GetEntries() > 0) ? "p" : "ap");
  if (gSigmaResPVsP2) gSigmaResPVsP2->Draw((gPad->GetListOfPrimitives()->GetEntries() > 0) ? "p" : "ap");
  if (gSigmaResPAtVtxVsPIn02deg) gSigmaResPAtVtxVsPIn02deg->Draw((gPad->GetListOfPrimitives()->GetEntries() > 0) ? "p" : "ap");
  if (gSigmaResPAtVtxVsPIn23deg) gSigmaResPAtVtxVsPIn23deg->Draw("p");
  if (gSigmaResPAtVtxVsPIn310deg) gSigmaResPAtVtxVsPIn310deg->Draw("p");
  gROOT->SetSelectedPad(cRes->cd(2));
  gPad->SetLogz();
  if (gSigmaResSlopeXVsP) gSigmaResSlopeXVsP->Draw((gPad->GetListOfPrimitives()->GetEntries() > 0) ? "p" : "ap");
  if (gSigmaResSlopeXVsP2) gSigmaResSlopeXVsP2->Draw((gPad->GetListOfPrimitives()->GetEntries() > 0) ? "p" : "ap");
  gROOT->SetSelectedPad(cRes->cd(3));
  gPad->SetLogz();
  if (gSigmaResSlopeYVsP) gSigmaResSlopeYVsP->Draw((gPad->GetListOfPrimitives()->GetEntries() > 0) ? "p" : "ap");
  if (gSigmaResSlopeYVsP2) gSigmaResSlopeYVsP2->Draw((gPad->GetListOfPrimitives()->GetEntries() > 0) ? "p" : "ap");
  
  // draw cluster-trackRef x-residuals
  for (Int_t iSt = 0; iSt < 5; ++iSt) {
    gROOT->SetSelectedPad(cResX->cd(iSt+13));
    gPad->SetLogy();
    hXResSt[iSt]->Draw();
  }
  gROOT->SetSelectedPad(cResX->cd(18));
  gPad->SetLogy();
  hXRes->Draw();
  
  // draw cluster-trackRef y-residuals
  for (Int_t iSt = 0; iSt < 5; ++iSt) {
    gROOT->SetSelectedPad(cResY->cd(iSt+13));
    gPad->SetLogy();
    hYResSt[iSt]->Draw();
  }
  gROOT->SetSelectedPad(cResY->cd(18));
  gPad->SetLogy();
  hYRes->Draw();
  
}

//-----------------------------------------------------------------------
void DrawResolutionFromSmearing(TCanvas *cRes, TString fileResQA, Bool_t atFirstCluster, Bool_t drawResPVsP, Int_t rebinP)
{
  /// Draw the resolutions vs p from the task AliTaskMuonTrackSmearingQA
  
  if (atFirstCluster) return; // no QA plots at first cluster yet
  
  if (fileResQA.IsNull()) return;
  
  TCanvas *dummy = new TCanvas("dummy","dummy",10,10,68,1);
  
  TGraphAsymmErrors /**gSigmaResPVsP = 0x0, */*gSigmaResSlopeXVsP = 0x0, *gSigmaResSlopeYVsP = 0x0;
  TGraphAsymmErrors *gSigmaResPAtVtxVsPIn02deg = 0x0, *gSigmaResPAtVtxVsPIn23deg = 0x0, *gSigmaResPAtVtxVsPIn310deg = 0x0;
  
  TFile *fResQA = TFile::Open(fileResQA.Data(),"READ");
  if (!fResQA || !fResQA->IsOpen()) return;
  
  TObjArray* resHistoList = static_cast<TObjArray*>(fResQA->FindObjectAny("ResHistos"));
  if (!resHistoList) return;
  
  if (atFirstCluster) {
    /*
    // compute momentum resolution versus p
     TH2F *hPPP = static_cast<TH2F*>(resHistoList->FindObject("hPPP"));
     if (!hPPP) return;
    dummy->cd();
    gSigmaResPVsP = new TGraphAsymmErrors(hPPP->GetNbinsX()/rebinP);
    gSigmaResPVsP->SetName("gSigmaResPAt1stClVsP3");
    gSigmaResPVsP->SetTitle("#sigma_{p}/p at 1st cluster versus p;p (GeV/c);#sigma_{p}/p (%)");
    gSigmaResPVsP->SetLineColor(9);
    FitPResVsP(hPPP, gSigmaResPVsP, rebinP, -1., drawResPVsP);
    
    // compute slope resolution versus p
     TH2F *hXXX = static_cast<TH2F*>(resHistoList->FindObject("hXXX"));
     TH2F *hYYY = static_cast<TH2F*>(resHistoList->FindObject("hYYY"));
     if (!hXXX || !hYYY) return;
    dummy->cd();
    gSigmaResSlopeXVsP = new TGraphAsymmErrors(hXXX->GetNbinsX()/rebinP);
    gSigmaResSlopeXVsP->SetName("gSigmaResSlopeXAt1stClVsP3");
    gSigmaResSlopeXVsP->SetTitle("#sigma_{slope_{X}} at 1st cluster versus p;p (GeV/c);#sigma_{slope_{X}}");
    gSigmaResSlopeXVsP->SetLineColor(4);
    FitGausResVsMom(hXXX, gSigmaResSlopeXVsP, rebinP);
    gSigmaResSlopeYVsP = new TGraphAsymmErrors(hYYY->GetNbinsX()/rebinP);
    gSigmaResSlopeYVsP->SetName("gSigmaResSlopeYAt1stClVsP3");
    gSigmaResSlopeYVsP->SetTitle("#sigma_{slope_{Y}} at 1st cluster versus p;p (GeV/c);#sigma_{slope_{Y}}");
    gSigmaResSlopeYVsP->SetLineColor(4);
    FitGausResVsMom(hYYY, gSigmaResSlopeYVsP, rebinP);
    */
  } else {
    
    // compute momentum resolution versus p
    TH2F *hResPAtVtxVsPIn02degMC = static_cast<TH2F*>(resHistoList->FindObject("hResPAtVtxVsPIn02deg"));
    TH2F *hResPAtVtxVsPIn23deg = static_cast<TH2F*>(resHistoList->FindObject("hResPAtVtxVsPIn23deg"));
    TH2F *hResPAtVtxVsPIn310deg = static_cast<TH2F*>(resHistoList->FindObject("hResPAtVtxVsPIn310deg"));
    if (!hResPAtVtxVsPIn02degMC || !hResPAtVtxVsPIn23deg || !hResPAtVtxVsPIn310deg) return;
    dummy->cd();
    gSigmaResPAtVtxVsPIn02deg = new TGraphAsymmErrors(hResPAtVtxVsPIn02degMC->GetNbinsX()/rebinP);
    gSigmaResPAtVtxVsPIn02deg->SetName("gSigmaResPAtVtxVsPIn02deg2");
    gSigmaResPAtVtxVsPIn02deg->SetTitle("#sigma_{p}/p at vertex versus p in [0,2[ deg MC;p (GeV/c);#sigma_{p}/p (%)");
    gSigmaResPAtVtxVsPIn02deg->SetLineColor(11);
    FitPResVsP(hResPAtVtxVsPIn02degMC, gSigmaResPAtVtxVsPIn02deg, rebinP, 300., drawResPVsP);
    gSigmaResPAtVtxVsPIn23deg = new TGraphAsymmErrors(hResPAtVtxVsPIn23deg->GetNbinsX()/rebinP);
    gSigmaResPAtVtxVsPIn23deg->SetName("gSigmaResPAtVtxVsPIn23deg2");
    gSigmaResPAtVtxVsPIn23deg->SetTitle("#sigma_{p}/p at vertex versus p in [2,3[ deg;p (GeV/c);#sigma_{p}/p (%)");
    gSigmaResPAtVtxVsPIn23deg->SetLineColor(4);
    FitPResVsP(hResPAtVtxVsPIn23deg, gSigmaResPAtVtxVsPIn23deg, rebinP, 220., drawResPVsP);
    gSigmaResPAtVtxVsPIn310deg = new TGraphAsymmErrors(hResPAtVtxVsPIn310deg->GetNbinsX()/rebinP);
    gSigmaResPAtVtxVsPIn310deg->SetName("gSigmaResPAtVtxVsPIn310deg2");
    gSigmaResPAtVtxVsPIn310deg->SetTitle("#sigma_{p}/p at vertex versus p in [3,10[ deg;p (GeV/c);#sigma_{p}/p (%)");
    gSigmaResPAtVtxVsPIn310deg->SetLineColor(9);
    FitPResVsP(hResPAtVtxVsPIn310deg, gSigmaResPAtVtxVsPIn310deg, rebinP, 160., drawResPVsP);
    
    // compute slope resolution versus p
    TH2F *hResSlopeXAtVtxVsP = static_cast<TH2F*>(resHistoList->FindObject("hResSlopeXAtVtxVsP"));
    TH2F *hResSlopeYAtVtxVsP = static_cast<TH2F*>(resHistoList->FindObject("hResSlopeYAtVtxVsP"));
    if (!hResSlopeXAtVtxVsP || !hResSlopeYAtVtxVsP) return;
    dummy->cd();
    gSigmaResSlopeXVsP = new TGraphAsymmErrors(hResSlopeXAtVtxVsP->GetNbinsX()/rebinP);
    gSigmaResSlopeXVsP->SetName("gSigmaResSlopeXAtVtxVsP3");
    gSigmaResSlopeXVsP->SetTitle("#sigma_{slope_{X}} at vertex versus p;p (GeV/c);#sigma_{slope_{X}}");
    gSigmaResSlopeXVsP->SetLineColor(4);
    FitGausResVsMom(hResSlopeXAtVtxVsP, gSigmaResSlopeXVsP, rebinP);
    gSigmaResSlopeYVsP = new TGraphAsymmErrors(hResSlopeYAtVtxVsP->GetNbinsX()/rebinP);
    gSigmaResSlopeYVsP->SetName("gSigmaResSlopeYAtVtxVsP3");
    gSigmaResSlopeYVsP->SetTitle("#sigma_{slope_{Y}} at vertex versus p;p (GeV/c);#sigma_{slope_{Y}}");
    gSigmaResSlopeYVsP->SetLineColor(4);
    FitGausResVsMom(hResSlopeYAtVtxVsP, gSigmaResSlopeYVsP, rebinP);
    
  }
  
  fResQA->Close();
  
  // draw momentum and slope resolution versus p
  gROOT->SetSelectedPad(cRes->cd(1));
  gPad->SetLogz();
  //if (gSigmaResPVsP) gSigmaResPVsP->Draw((gPad->GetListOfPrimitives()->GetEntries() > 0) ? "p" : "ap");
  if (gSigmaResPAtVtxVsPIn02deg) gSigmaResPAtVtxVsPIn02deg->Draw((gPad->GetListOfPrimitives()->GetEntries() > 0) ? "p" : "ap");
  if (gSigmaResPAtVtxVsPIn23deg) gSigmaResPAtVtxVsPIn23deg->Draw("p");
  if (gSigmaResPAtVtxVsPIn310deg) gSigmaResPAtVtxVsPIn310deg->Draw("p");
  gROOT->SetSelectedPad(cRes->cd(2));
  gPad->SetLogz();
  if (gSigmaResSlopeXVsP) gSigmaResSlopeXVsP->Draw((gPad->GetListOfPrimitives()->GetEntries() > 0) ? "p" : "ap");
  gROOT->SetSelectedPad(cRes->cd(3));
  gPad->SetLogz();
  if (gSigmaResSlopeYVsP) gSigmaResSlopeYVsP->Draw((gPad->GetListOfPrimitives()->GetEntries() > 0) ? "p" : "ap");

}

//-----------------------------------------------------------------------
Double_t GetSigma(TH1 *h, Double_t *sigmaErr)
{
  /// get the dispersion of the histo
  switch (muonTrackSmearing.GetChosenFunc()) {
    case AliMuonTrackSmearing::kCrystalBall:
      return GetSigmaCrystalBall(h, sigmaErr);
      break;
    case AliMuonTrackSmearing::kBreitWigner:
      return GetSigmaBreitWigner(h, sigmaErr);
      break;
    case AliMuonTrackSmearing::kGaus:
      return GetSigmaGaus(h, sigmaErr);
      break;
    default:
      printf("Unknown option chosen!\n");
      return -1.;
      break;
  }
}

//-----------------------------------------------------------------------
Double_t GetSigmaCrystalBall(TH1 *h, Double_t *sigmaErr)
{
  /// get the dispersion of the histo
  
  static TF1 *fCrystalBall = new TF1("CrystalBall", &muonTrackSmearing, &AliMuonTrackSmearing::CrystalBallSymmetric, -2.5, 2.5, 5, "AliMuonTrackSmearing", "CrystalBallSymmetric");

  if (h->GetEntries() < 10.) {
    if (sigmaErr) *sigmaErr = 0.;
    return 0.;
  }
  
  Double_t sigmaTrk = muonTrackSmearing.GetSigmaTrk();
  Double_t sigmaTrkCut = muonTrackSmearing.GetSigmaTrkCut();
  
  // first fit
  Double_t xMin = -sigmaTrkCut*sigmaTrk*100.;
  Double_t xMax = sigmaTrkCut*sigmaTrk*100.;
  fCrystalBall->SetRange(xMin, xMax);
  fCrystalBall->SetParameters(h->GetEntries(), 0., 0.1, 2., 1.5);
  fCrystalBall->SetParLimits(1, xMin, xMax);
  fCrystalBall->SetParLimits(2, 0., 1.);
  fCrystalBall->FixParameter(3, 1.e6);
  h->Fit(fCrystalBall, "WWRNQ");
  
  // rebin histo
  Int_t rebin = static_cast<Int_t>(TMath::Min(0.1*h->GetNbinsX(),TMath::Max(0.3*fCrystalBall->GetParameter(2)/h->GetBinWidth(1),1.)));
  while (h->GetNbinsX()%rebin!=0) rebin--;
  h->Rebin(rebin);
  
  // second fit
  fCrystalBall->ReleaseParameter(3);
  fCrystalBall->SetParameter(3, 2.);
  fCrystalBall->SetParameter(4, 1.5);
  h->Fit(fCrystalBall,"RQ");
  
  if (!TString(h->GetName()).Contains("In")) printf("alpha = %f; n = %f\n", fCrystalBall->GetParameter(3), fCrystalBall->GetParameter(4));
  if (sigmaErr) *sigmaErr = fCrystalBall->GetParError(2);
  return fCrystalBall->GetParameter(2);
  
}

//-----------------------------------------------------------------------
Double_t GetSigmaGaus(TH1 *h, Double_t *sigmaErr)
{
  /// get the dispersion of the histo
  
  static TF1 *fGaus = new TF1("fGaus","gaus");
  
  if (h->GetEntries() < 10.) {
    if (sigmaErr) *sigmaErr = 0.;
    return 0.;
  }
  
  // first fit
  Double_t xMin = h->GetXaxis()->GetXmin();
  Double_t xMax = h->GetXaxis()->GetXmax();
  fGaus->SetRange(xMin, xMax);
  fGaus->SetParameters(h->GetEntries(), 0., 0.1);
  fGaus->SetParLimits(1, xMin, xMax);
  h->Fit("fGaus", "WWNQ");
  
  // rebin histo
  Int_t rebin = static_cast<Int_t>(TMath::Min(0.1*h->GetNbinsX(),TMath::Max(0.3*fGaus->GetParameter(2)/h->GetBinWidth(1),1.)));
  while (h->GetNbinsX()%rebin!=0) rebin--;
  h->Rebin(rebin);
  
  // second fit
  xMin = TMath::Max(fGaus->GetParameter(1)-10.*fGaus->GetParameter(2), h->GetXaxis()->GetXmin());
  xMax = TMath::Min(fGaus->GetParameter(1)+10.*fGaus->GetParameter(2), h->GetXaxis()->GetXmax());
  fGaus->SetRange(xMin, xMax);
  fGaus->SetParLimits(1, xMin, xMax);
  h->Fit("fGaus","RQ");
  
  if (sigmaErr) *sigmaErr = fGaus->GetParError(2);
  return fGaus->GetParameter(2);
  
}

//-----------------------------------------------------------------------
Double_t GetSigmaBreitWigner(TH1 *h, Double_t *sigmaErr)
{
  /// get the dispersion of the histo
  
  static TF1 *fBreitWigner = new TF1("fBreitWigner","[0]*TMath::BreitWigner(x,[1],[2])",-10.,10.);
  
  if (h->GetEntries() < 10.) {
    if (sigmaErr) *sigmaErr = 0.;
    return 0.;
  }
  
  // first fit
  Double_t xMin = h->GetXaxis()->GetXmin()/2.;
  Double_t xMax = h->GetXaxis()->GetXmax()/2.;
  fBreitWigner->SetRange(xMin, xMax);
  fBreitWigner->SetParameters(h->GetEntries(), 0., 0.01);
  fBreitWigner->SetParLimits(1, xMin, xMax);
  h->Fit("fBreitWigner", "WWNQ");
  
  // rebin histo
  Int_t rebin = static_cast<Int_t>(TMath::Min(0.1*h->GetNbinsX(),TMath::Max(0.1*fBreitWigner->GetParameter(2)/h->GetBinWidth(1),1.)));
  while (h->GetNbinsX()%rebin!=0) rebin--;
  h->Rebin(rebin);
  
  // second fit
  xMin = TMath::Max(fBreitWigner->GetParameter(1)-10.*fBreitWigner->GetParameter(2), h->GetXaxis()->GetXmin()/2.);
  xMax = TMath::Min(fBreitWigner->GetParameter(1)+10.*fBreitWigner->GetParameter(2), h->GetXaxis()->GetXmax()/2.);
  fBreitWigner->SetRange(xMin, xMax);
  fBreitWigner->SetParLimits(1, xMin, xMax);
  h->Fit("fBreitWigner","RQI");
  
  if (sigmaErr) *sigmaErr = fBreitWigner->GetParError(2);
  return fBreitWigner->GetParameter(2);
  
}

//-----------------------------------------------------------------------
void FitGausResVsMom(TH2 *h, TGraphAsymmErrors *gSigma, Int_t rebinP)
{
  /// generic function to fit residuals versus momentum with a gaussian
  
  for (Int_t i = rebinP; i <= h->GetNbinsX(); i+=rebinP) {
    
    if (i+rebinP <= h->GetNbinsX()) printf("compute slope resolution versus p... %d/%d\r",i/rebinP,h->GetNbinsX()/rebinP);
    else printf("compute slope resolution versus p... %d/%d\n",i/rebinP,h->GetNbinsX()/rebinP);
    
    TH1D *tmp = h->ProjectionY(Form("%s_%d",h->GetName(),i/rebinP),i-rebinP+1,i,"e");
    tmp->SetDirectory(0);
    if (tmp->GetEntries() < 10.) {
      delete tmp;
      continue;
    }
    
    //TCanvas *c = new TCanvas(Form("c%s_%d",h->GetName(),i/rebinP),Form("c%s_%d",h->GetName(),i/rebinP));
    //tmp->Draw();
    
    // get fit results and fill graph
    Double_t sigmaErr = 0.;
    Double_t sigma = GetSigmaGaus(tmp, &sigmaErr);
    h->GetXaxis()->SetRange(i-rebinP+1,i);
    Double_t p = (tmp->GetEntries() > 0) ? h->GetMean() : 0.5 * (h->GetXaxis()->GetBinLowEdge(i-rebinP+1) + h->GetXaxis()->GetBinLowEdge(i+1));
    h->GetXaxis()->SetRange();
    Double_t pErr[2] = {p-h->GetXaxis()->GetBinLowEdge(i-rebinP+1), h->GetXaxis()->GetBinLowEdge(i+1)-p};
    gSigma->SetPoint(i/rebinP-1, p, sigma);
    gSigma->SetPointError(i/rebinP-1, pErr[0], pErr[1], sigmaErr, sigmaErr);
    
    // clean memory
    //c->Update();
    delete tmp;
    
  }
  
}

//-----------------------------------------------------------------------
Double_t langaufun(Double_t *x, Double_t *par)
{
  /// Landau convoluted with a Gaussian
  //
  //Fit parameters:
  //par[0]=Width (scale) parameter of Landau density
  //par[1]=Most Probable (MP, location) parameter of Landau density
  //par[2]=Total area (integral -inf to inf, normalization constant)
  //par[3]=Width (sigma) of convoluted Gaussian function
  //
  //In the Landau distribution (represented by the CERNLIB approximation),
  //the maximum is located at x=-0.22278298 with the location parameter=0.
  //This shift is corrected within this function, so that the actual
  //maximum is identical to the MP parameter.
  
  // Numeric constants
  Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  Double_t mpshift  = -0.22278298;       // Landau maximum location
  
  // MP shift correction
  Double_t mpc = par[1] - mpshift * par[0];
  //return par[2]*TMath::Landau(-x[0],mpc,par[0])/par[0];
  //return par[2]*TMath::Gaus(x[0],par[1],par[3])*invsq2pi/par[3];
  
  // Range of convolution integral
  Double_t xlow = x[0] - 5. * par[3];
  Double_t xupp = x[0] + 5. * par[3];
  Double_t step = TMath::Min(0.2*par[0],0.1*par[3]);
  
  // Convolution integral of Landau and Gaussian by sum
  Double_t xx = xlow + 0.5 * step;
  Double_t sum = 0.0;
  while (xx < xupp) {
    sum += TMath::Landau(-xx,mpc,par[0]) * TMath::Gaus(x[0],xx,par[3]);
    xx += step;
  }
 
  return (par[2] * step * sum * invsq2pi / par[3] / par[0]);
}

//-----------------------------------------------------------------------
Double_t GausInvXFun(Double_t *x, Double_t *par)
{
  /// norm * k/p/p * gaus(k/p, k/p0, k/p0/p0*sigmap)
  /// with k = factor to convert momentum into angular deviation
  
  Double_t p0 = par[1];
  Double_t p = x[0] + p0 + par[3];
  Double_t s0 = muonTrackSmearing.PToThetaDev(p0);
  Double_t s = muonTrackSmearing.PToThetaDev(p);
  Double_t dsdp = s/p;
  
  return par[0] * s/p * TMath::Gaus(s, s0, s0/p0*par[2], kTRUE);
}

//-----------------------------------------------------------------------
void FitPResVsP(TH2 *h, TGraphAsymmErrors *gSigma, Int_t rebinP, Double_t pSwitch, Bool_t drawResPVsP)
{
  /// generic function to fit momentum residuals versus momentum
  /// if p < pSwitch: use a landau convoluted with a gaussian
  /// if p > pSwitch: use a gaus(k/p) distribution
  
  static TF1 *fGaus2 = new TF1("fGaus2","gausn");
  static TF1 *fLandauGaus = 0x0;
  if (!fLandauGaus) {
    fLandauGaus = new TF1("fLandauGaus",langaufun,h->GetYaxis()->GetBinLowEdge(1),h->GetYaxis()->GetBinLowEdge(h->GetNbinsY()+1),4);
    fLandauGaus->SetNpx(10000);
  }
  static TF1 *fGausInvP = new TF1("fGausInvP",GausInvXFun, h->GetYaxis()->GetBinLowEdge(1),h->GetYaxis()->GetBinLowEdge(h->GetNbinsY()+1),4);
  
  for (Int_t i = rebinP; i <= h->GetNbinsX(); i+=rebinP) {
    
    if (i+rebinP <= h->GetNbinsX()) printf("compute momentum resolution versus p... %d/%d\r",i/rebinP,h->GetNbinsX()/rebinP);
    else printf("compute momentum resolution versus p... %d/%d\n",i/rebinP,h->GetNbinsX()/rebinP);
    
    TH1D *tmp = h->ProjectionY(Form("%s_%d",h->GetName(),i/rebinP),i-rebinP+1,i,"e");
    tmp->SetDirectory(0);
    
    // draw histogram if required
    TCanvas *c = 0x0;
    if (drawResPVsP) {
      
      TString mc = "";
      c = static_cast<TCanvas*>(gROOT->GetListOfCanvases()->FindObject(Form("c%s_%d",h->GetName(),i/rebinP)));
      if (!c) {
        c = static_cast<TCanvas*>(gROOT->GetListOfCanvases()->FindObject(Form("c%sMC_%d",h->GetName(),i/rebinP)));
        mc = "MC";
      }
      if (c) {
        TH1D *tmp0 = static_cast<TH1D*>(c->GetListOfPrimitives()->FindObject(Form("%s%s_%d",h->GetName(),mc.Data(),i/rebinP)));
        tmp->Scale(tmp0->GetSumOfWeights()*tmp0->GetBinWidth(1)/tmp->GetSumOfWeights()/tmp->GetBinWidth(1));
        c->cd();
        tmp->SetLineColor(4);
        fLandauGaus->SetLineColor(4);
        fGausInvP->SetLineColor(4);
        tmp->Draw("sames");
      } else {
        c = new TCanvas(Form("c%s_%d",h->GetName(),i/rebinP),Form("c%s_%d",h->GetName(),i/rebinP));
        tmp->SetLineColor(2);
        fLandauGaus->SetLineColor(2);
        fGausInvP->SetLineColor(2);
        tmp->Draw();
      }
      
    }
    
    // get mean momentum in this bin
    h->GetXaxis()->SetRange(i-rebinP+1,i);
    Double_t p = (tmp->GetEntries() > 0) ? h->GetMean() : 0.5 * (h->GetXaxis()->GetBinLowEdge(i-rebinP+1) + h->GetXaxis()->GetBinLowEdge(i+1));
    h->GetXaxis()->SetRange();
    
    // extract the resolution if enough entries
    Double_t sigmaP = 0., sigmaPErr = 0.;
    if (tmp->GetEntries() > 50.) {
      
      // rebin histo
      Double_t sigma = tmp->GetRMS();
      Int_t rebin = static_cast<Int_t>(TMath::Min(0.1*tmp->GetNbinsX(),TMath::Max(0.1*sigma/tmp->GetBinWidth(1),1.)));
      while (tmp->GetNbinsX()%rebin!=0) rebin--;
      tmp->Rebin(rebin);
      tmp->Scale(1./rebin);
      
      // first fit
      fGaus2->SetParameters(tmp->GetSumOfWeights()*tmp->GetBinWidth(1), 0., tmp->GetRMS());
      Double_t xMin = TMath::Max(-TMath::Max(sigma,2.*tmp->GetBinWidth(1)), tmp->GetXaxis()->GetXmin());
      Double_t xMax = TMath::Min(TMath::Max(sigma,2.*tmp->GetBinWidth(1)), tmp->GetXaxis()->GetXmax());
      if (xMin < tmp->GetXaxis()->GetXmax() && xMax > tmp->GetXaxis()->GetXmin()) fGaus2->SetRange(xMin, xMax);
      tmp->Fit("fGaus2", "RQ");
      
      // rebin histo
      sigma = fGaus2->GetParameter(2);
      rebin = static_cast<Int_t>(TMath::Min(0.1*tmp->GetNbinsX(),TMath::Max(0.5*sigma/tmp->GetBinWidth(1),1.)));
      while (tmp->GetNbinsX()%rebin!=0) rebin--;
      tmp->Rebin(rebin);
      tmp->Scale(1./rebin);
      
      // switch between fitting functions
      if (p < pSwitch) {
        
        // second fit
        Double_t mean = fGaus2->GetParameter(1);
        fLandauGaus->SetParameters(0.25*sigma*TMath::Sqrt(8.*log(2.)), mean, fGaus2->GetParameter(0), 0.5*sigma);
        fLandauGaus->SetParLimits(0, 0.0025*sigma*TMath::Sqrt(8.*log(2.)), 1000.);
        fLandauGaus->SetParLimits(3, 0., 2.*sigma);
        xMin = TMath::Max(mean-50.*sigma, tmp->GetXaxis()->GetXmin());
        xMax = TMath::Min(mean+10.*sigma, tmp->GetXaxis()->GetXmax());
        if (xMin < tmp->GetXaxis()->GetXmax() && xMax > tmp->GetXaxis()->GetXmin()) fLandauGaus->SetRange(xMin, xMax);
        tmp->Fit("fLandauGaus","RQ");
        
        // third fit
        xMin = TMath::Max(mean-4.*sigma, tmp->GetXaxis()->GetXmin());
        xMax = TMath::Min(mean+3.*sigma, tmp->GetXaxis()->GetXmax());
        if (xMin < tmp->GetXaxis()->GetXmax() && xMax > tmp->GetXaxis()->GetXmin()) fLandauGaus->SetRange(xMin, xMax);
        tmp->Fit("fLandauGaus","RQ");
        
        // get fit results
        Double_t fwhm = 4.*fLandauGaus->GetParameter(0);
        sigma = fLandauGaus->GetParameter(3);
        sigmaP = TMath::Sqrt(sigma*sigma + fwhm*fwhm/(8.*log(2.)));
        Double_t fwhmErr = fLandauGaus->GetParError(0);
        Double_t sigmaErr = fLandauGaus->GetParError(3);
        sigmaPErr = TMath::Sqrt(sigma*sigma*sigmaErr*sigmaErr + fwhm*fwhm*fwhmErr*fwhmErr/(64.*log(2.)*log(2.))) / sigmaP;
        
      } else {
        
        // second fit
        Double_t mean = fGaus2->GetParameter(1);
        fGausInvP->SetParameters(fGaus2->GetParameter(0), p, sigma, 0.);
        Double_t xMin = TMath::Max(mean-TMath::Max(2.*sigma,2.*tmp->GetBinWidth(1)), tmp->GetXaxis()->GetXmin());
        Double_t xMax = TMath::Min(mean+TMath::Max(3.*sigma,2.*tmp->GetBinWidth(1)), tmp->GetXaxis()->GetXmax());
        if (xMin < tmp->GetXaxis()->GetXmax() && xMax > tmp->GetXaxis()->GetXmin()) fGausInvP->SetRange(xMin, xMax);
        tmp->Fit("fGausInvP","RIQ");
        
        // get fit results
        sigmaP = fGausInvP->GetParameter(2);
        sigmaPErr = fGausInvP->GetParError(2);
        
      }
      
    }
    
    // fill graph
    Double_t pErr[2] = {p-h->GetXaxis()->GetBinLowEdge(i-rebinP+1), h->GetXaxis()->GetBinLowEdge(i+1)-p};
    gSigma->SetPoint(i/rebinP-1, p, 100.*sigmaP/p);
    gSigma->SetPointError(i/rebinP-1, pErr[0], pErr[1], 100.*sigmaPErr/p, 100.*sigmaPErr/p);
    
    // clean memory
    if (drawResPVsP) c->Update();
    else delete tmp;
    
  }
  
}

