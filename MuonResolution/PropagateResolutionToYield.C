//
//  PropagateResolutionToYield.C
//  aliroot_dev
//
//  Created by philippe pillot on 24/07/2014.
//  Copyright (c) 2014 Philippe Pillot. All rights reserved.
//


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
#include <TRandom3.h>
#include <TCanvas.h>
#include <TDatabasePDG.h>

#endif

enum {kCrystalBall, kBreitWigner, kGaus};
Int_t chosenFunc = kCrystalBall;

Bool_t tuneData = kTRUE;
Double_t sigmaxChSt1 = -1.; // used to compute uncertainty on the slope at vtx
Double_t sigmayChSt1 = -1.; // used to compute uncertainty on the slope at vtx
Double_t sigmayCh = -1.; // used to compute uncertainty on the momentum at vtx
Double_t tailxChSt1[2] = {0., 0.}; // CB tails for station 1 resolution dispersion along x
Double_t tailyChSt1[2] = {0., 0.}; // CB tails for station 1 resolution dispersion along y
Double_t tailyCh[2] = {0., 0.}; // CB tails for spectrometer resolution dispersion along y

void InitSigmas()
{
  /// average (non-)bending resolution at the chamber level (m)
  /// to get the resolution at the station level it has to be divided by sqrt(2)
  
  if (tuneData) { // --- DATA ---
    
    switch ( chosenFunc ) {
      case kCrystalBall:
        // Crystal Ball fit
        sigmaxChSt1 = 0.000401;
        tailxChSt1[0] = 2.017293;
        tailxChSt1[1] = 1.890778;
        sigmayChSt1 = 0.000110;
        tailyChSt1[0] = 1.514588;
        tailyChSt1[1] = 1.938707;
        sigmayCh = 0.000153;
        tailyCh[0] = 1.280116;
        tailyCh[1] = 2.239019;
        break;
      case kBreitWigner:
        // Breit Wigner fit (correspond to gamma values = FWHM/2.)
        sigmaxChSt1 = 0.000550;
        sigmayChSt1 = 0.000190;
        sigmayCh = 0.000280;
        break;
      case kGaus:
        // gaussian fit
        sigmaxChSt1 = 0.000419;
        sigmayChSt1 = 0.000130;
        sigmayCh = 0.000207;
        break;
      default:
        break;
    }
    
  } else { // --- SIM ---
    
    switch ( chosenFunc ) {
      case kCrystalBall:
        // Crystal Ball fit
        sigmaxChSt1 = 0.000519;
        tailxChSt1[0] = 3.671484;
        tailxChSt1[1] = 4.590817;
        sigmayChSt1 = 0.000064;
        tailyChSt1[0] = 2.453060;
        tailyChSt1[1] = 1.715218;
        sigmayCh = 0.000105;
        tailyCh[0] = 2.145035;
        tailyCh[1] = 1.782320;
        break;
      case kBreitWigner:
        // Breit Wigner fit (correspond to gamma values = FWHM/2.)
        sigmaxChSt1 = 0.000665;
        sigmayChSt1 = 0.000166;
        sigmayCh = 0.000224;
        break;
      case kGaus:
        // gaussian fit
        sigmaxChSt1 = 0.000506;
        sigmayChSt1 = 0.000065;
        sigmayCh = 0.000112;
        break;
      default:
        break;
    }
    
  }
  
}

Double_t sigmaTrk = 0.002; // chamber resolution used during tracking
Double_t sigmaTrkCut = 4.; // sigma cut used during tracking

// tune the parameterization of MCS and energy loss to fit the momentum and angular resolution given by:
// - kTRUE:  the Kalman filter
// - kFALSE: the performance task
Bool_t tuneKalman = kFALSE;

// draw resolution at first cluster
Bool_t atFirstCluster = kFALSE;

/*
 absorber parameters from geometry for a muon of 10 GeV (using AliMUONTrackExtrap)
 0 < theta < 2 (yAbsEnd=16cm):
  x/x0 = 828.913425
  zB = 382.142805
  eLoss = 9.82487194539394793e+00
  sigmaELoss = 0.355647
  sigmaELoss = 0.621332 (linear sum)
 2 < theta < 3 (yAbsEnd=20cm):
  x/x0 = 134.193184
  zB = 471.874591
  eLoss = 3.00994537855345357e+00
  sigmaELoss = 0.087540
  sigmaELoss = 0.183265 (linear sum)
 3 < theta < 10 (yAbsEnd=50cm):
  x/x0 = 55.905367
  zB = 446.183741
  eLoss = 2.52037057051776614e+00
  sigmaELoss = 0.067273
  sigmaELoss = 0.150743 (linear sum)
*/

// branson plan position
Double_t zB02 = 3.82;
Double_t zB23 = 4.72;
Double_t zB310 = 4.46;

// minimum momentum to cross the whole spectrometer
Double_t pAccCut = 0.5;

// pT/eta/phi ranges for generation
//Double_t pTRange[2] = {0.8,50.};
//Double_t pTRange[2] = {1.,100.};
Double_t pTRange[2] = {1.,15.};
Double_t npTBinPerGeV = 2.;
Double_t etaRange[2] = {-4.2,-2.3};
Double_t nEtaBinPerUnit = 50.;
Double_t phiRange[2] = {0.,360.};
Double_t nPhiBinPerUnit = 0.25;

// number of events to generate
Int_t nEvents = 100000;

// generate uniform in pT/eta and weight by the distribution instead of generating according to the distributions
Bool_t generateUniform = kTRUE;

// generate in pT/y instead of pT/eta
Bool_t generateY = kTRUE;

// draw Delta_p histograms versus p
Bool_t drawResPVsP = kTRUE;

Double_t DnDpT( const Double_t *x, const Double_t *par );
Double_t DnDeta( const Double_t *x, const Double_t *par );
Double_t PResVsP( const Double_t *x, const Double_t *par );
Double_t SlopeResVsP( const Double_t *x, const Double_t *par );
Double_t PToThetaDev(Double_t p);
Double_t ThetaDevToP(Double_t thetaDev);
Double_t FWHMELoss2(Double_t p, Double_t theta);
Double_t SigmaThetaDevFromMCS2(Double_t p);
Double_t SigmaThetaDevFromRes2();
Double_t SigmaSlopeFromMCSInAbs2(Double_t p, Double_t theta);
Double_t SigmaSlopeFromMCSInCh2(Double_t p, Bool_t at1stCl, Double_t zB);
//Double_t SigmaPosFromRes2(Bool_t bendingDir, Bool_t at1stCl, Double_t zB);
Double_t SigmaSlopeFromRes2(Bool_t bendingDir, Bool_t at1stCl, Double_t zB);
Double_t ELoss(Double_t p, Double_t theta);
Double_t ELossFluctuation2(Double_t p, Double_t rhoZoverA);
Double_t MCS2(Double_t p, Double_t dZ, Double_t x0);
Double_t GetSigma(TH1 *h, Double_t *sigmaErr = 0x0);
Double_t GetSigmaCrystalBall(TH1 *h, Double_t *sigmaErr = 0x0);
Double_t GetSigmaGaus(TH1 *h, Double_t *sigmaErr = 0x0);
Double_t GetSigmaBreitWigner(TH1 *h, Double_t *sigmaErr = 0x0);
void FitGausResVsMom(TH2 *h, TGraphAsymmErrors *gSigma, Int_t rebinP);
Double_t langaufun(Double_t *x, Double_t *par);
Double_t GausInvXFun(Double_t *x, Double_t *par);
void FitPResVsP(TH2 *h, TGraphAsymmErrors *gSigma, Int_t rebinP, Double_t pSwitch);
Double_t CrystalBallSymmetric(Double_t *x,Double_t *par);
Double_t GenRndGaus(Double_t mean, Double_t sigma);
Double_t GenRndBreitWigner(Double_t mean, Double_t sigma, Double_t max);
Double_t GenRndCrystalBall(Double_t mean, Double_t sigma, Double_t tail[2], Double_t max);


//-----------------------------------------------------------------------
void PropagateResolutionToYield(TString fileResMeas = "results.root",
                                TString fileResSim = "AnalysisResults.root",
                                Bool_t checkResVsP = kFALSE)
{
  /// propagate the momentum resolution to an uncertainty on the muon yield vs pT
  
  InitSigmas();
  if (sigmaxChSt1 < 0. || sigmayChSt1 < 0. || sigmayCh < 0.) return;
  
  gRandom->SetSeed(0);
  
  setbuf(stdout, NULL);
  TCanvas *dummy = new TCanvas("dummy","dummy",10,10,68,1);
  
  // momentum relative resolution versus p
  TF1 *fPResVsP02 = new TF1("fPResVsP02",PResVsP,0.,1000.,1);
  fPResVsP02->SetParameter(0,1);
  fPResVsP02->SetNpx(10000);
  fPResVsP02->SetLineColor(15);
  TF1 *fPResVsP23 = new TF1("fPResVsP23",PResVsP,0.,1000.,1);
  fPResVsP23->SetParameter(0,2.5);
  fPResVsP23->SetNpx(10000);
  fPResVsP23->SetLineColor(2);
  TF1 *fPResVsP310 = new TF1("fPResVsP310",PResVsP,0.,1000.,1);
  fPResVsP310->SetParameter(0,6);
  fPResVsP310->SetNpx(10000);
  fPResVsP310->SetLineColor(6);
  // slope resolution versus p
  TF1 *fSlopeXResVsP02 = new TF1("fSlopeXResVsP02",SlopeResVsP,0.,1000.,2);
  fSlopeXResVsP02->SetParameters(1.5,-1.);
  fSlopeXResVsP02->SetNpx(10000);
  fSlopeXResVsP02->SetLineColor(15);
  TF1 *fSlopeXResVsP23 = new TF1("fSlopeXResVsP23",SlopeResVsP,0.,1000.,2);
  fSlopeXResVsP23->SetParameters(2.5,-1.);
  fSlopeXResVsP23->SetNpx(10000);
  fSlopeXResVsP23->SetLineColor(2);
  TF1 *fSlopeXResVsP310 = new TF1("fSlopeXResVsP310",SlopeResVsP,0.,1000.,2);
  fSlopeXResVsP310->SetParameters(6,-1.);
  fSlopeXResVsP310->SetNpx(10000);
  fSlopeXResVsP310->SetLineColor(6);
  TF1 *fSlopeYResVsP02 = new TF1("fSlopeYResVsP02",SlopeResVsP,0.,1000.,2);
  fSlopeYResVsP02->SetParameters(1.5,1.);
  fSlopeYResVsP02->SetNpx(10000);
  fSlopeYResVsP02->SetLineColor(15);
  TF1 *fSlopeYResVsP23 = new TF1("fSlopeYResVsP23",SlopeResVsP,0.,1000.,2);
  fSlopeYResVsP23->SetParameters(2.5,1.);
  fSlopeYResVsP23->SetNpx(10000);
  fSlopeYResVsP23->SetLineColor(2);
  TF1 *fSlopeYResVsP310 = new TF1("fSlopeYResVsP310",SlopeResVsP,0.,1000.,2);
  fSlopeYResVsP310->SetParameters(6,1.);
  fSlopeYResVsP310->SetNpx(10000);
  fSlopeYResVsP310->SetLineColor(6);
  
  // read ouput of the task AliAnalysisTaskMuonResolution
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
  if (fResMeas) {
    
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
  }
  
  // read ouput of the task AliAnalysisTaskMuonPerformance
  TGraphAsymmErrors *gSigmaResPAtVtxVsP = 0x0, *gSigmaResSlopeXAtVtxVsP = 0x0, *gSigmaResSlopeYAtVtxVsP = 0x0;
  TGraphAsymmErrors *gSigmaResPAtVtxVsPIn02deg = 0x0, *gSigmaResPAtVtxVsPIn23deg = 0x0, *gSigmaResPAtVtxVsPIn310deg = 0x0;
  TGraphAsymmErrors *gSigmaResSlopeXAtVtxVsP2 = 0x0, *gSigmaResSlopeYAtVtxVsP2 = 0x0;
  TH1D *hXResSt[5], *hXRes = 0x0;
  TH1D *hYResSt[5], *hYRes = 0x0;
  for (Int_t iSt = 0; iSt < 5; ++iSt) {
    hXResSt[iSt] = 0x0;
    hYResSt[iSt] = 0x0;
  }
  TFile *fResSim = TFile::Open(fileResSim.Data(),"READ");
  if (fResSim) {
    
    // get the momentum and slope resolution from simulation to control their parameterization
    if (atFirstCluster) {
      TObjArray* momentumAtVtxList = static_cast<TObjArray*>(fResSim->FindObjectAny("MomentumAtFirstCl"));
      TObjArray* slopeAtVtxList = static_cast<TObjArray*>(fResSim->FindObjectAny("SlopeAtFirstCl"));
      if (!momentumAtVtxList || !slopeAtVtxList) return;
      gSigmaResPAtVtxVsP = static_cast<TGraphAsymmErrors*>(momentumAtVtxList->FindObject("gSigmaResPAt1stClVsP"));
      gSigmaResSlopeXAtVtxVsP = static_cast<TGraphAsymmErrors*>(slopeAtVtxList->FindObject("gSigmaResSlopeXAt1stClVsP"));
      gSigmaResSlopeYAtVtxVsP = static_cast<TGraphAsymmErrors*>(slopeAtVtxList->FindObject("gSigmaResSlopeYAt1stClVsP"));
    } else {
      TObjArray* momentumAtVtxList = static_cast<TObjArray*>(fResSim->FindObjectAny("MomentumAtVtx"));
      TObjArray* slopeAtVtxList = static_cast<TObjArray*>(fResSim->FindObjectAny("SlopeAtVtx"));
      if (!momentumAtVtxList || !slopeAtVtxList) return;
      gSigmaResPAtVtxVsP = static_cast<TGraphAsymmErrors*>(momentumAtVtxList->FindObject("gSigmaResPAtVtxVsP"));
      gSigmaResSlopeXAtVtxVsP = static_cast<TGraphAsymmErrors*>(slopeAtVtxList->FindObject("gSigmaResSlopeXAtVtxVsP"));
      gSigmaResSlopeYAtVtxVsP = static_cast<TGraphAsymmErrors*>(slopeAtVtxList->FindObject("gSigmaResSlopeYAtVtxVsP"));
    }
    if (!gSigmaResPAtVtxVsP || !gSigmaResSlopeXAtVtxVsP || !gSigmaResSlopeYAtVtxVsP) return;
    
    // get the cluster-trackRef residuals to control their parameterization
    TObjArray* trackerResolutionList = static_cast<TObjArray*>(fResSim->FindObjectAny("TrackerResolution"));
    if (!trackerResolutionList) return;
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
    
    if (checkResVsP) {
      
      if (atFirstCluster) {
        
        // recompute momentum resolution at 1st cluster versus p
        TH2F *hResPAtVtxVsPIn310deg = static_cast<TH2F*>(trackerResolutionList->FindObject("hResPAt1stClVsP"));
        if (!hResPAtVtxVsPIn310deg) return;
        dummy->cd();
        Int_t rebinP = 2;
        gSigmaResPAtVtxVsPIn310deg = new TGraphAsymmErrors(hResPAtVtxVsPIn310deg->GetNbinsX()/rebinP);
        gSigmaResPAtVtxVsPIn310deg->SetName("gSigmaResPAt1stClVsP");
        gSigmaResPAtVtxVsPIn310deg->SetTitle("#sigma_{p}/p at 1st cluster versus p;p (GeV/c);#sigma_{p}/p (%)");
        gSigmaResPAtVtxVsPIn310deg->SetLineColor(6);
        FitPResVsP(hResPAtVtxVsPIn310deg, gSigmaResPAtVtxVsPIn310deg, rebinP, -1.);
        
        // recompute slope resolution at 1st cluster versus p
        TH2F *hResSlopeXAtVtxVsP = static_cast<TH2F*>(trackerResolutionList->FindObject("hResSlopeXAt1stClVsP"));
        TH2F *hResSlopeYAtVtxVsP = static_cast<TH2F*>(trackerResolutionList->FindObject("hResSlopeYAt1stClVsP"));
        if (!hResSlopeXAtVtxVsP || !hResSlopeYAtVtxVsP) return;
        dummy->cd();
        gSigmaResSlopeXAtVtxVsP2 = new TGraphAsymmErrors(hResSlopeXAtVtxVsP->GetNbinsX()/rebinP);
        gSigmaResSlopeXAtVtxVsP2->SetName("gSigmaResSlopeXAt1stClVsP2");
        gSigmaResSlopeXAtVtxVsP2->SetTitle("#sigma_{slope_{X}} at 1st cluster versus p;p (GeV/c);#sigma_{slope_{X}}");
        gSigmaResSlopeXAtVtxVsP2->SetLineColor(2);
        FitGausResVsMom(hResSlopeXAtVtxVsP, gSigmaResSlopeXAtVtxVsP2, rebinP);
        gSigmaResSlopeYAtVtxVsP2 = new TGraphAsymmErrors(hResSlopeYAtVtxVsP->GetNbinsX()/rebinP);
        gSigmaResSlopeYAtVtxVsP2->SetName("gSigmaResSlopeYAt1stClVsP2");
        gSigmaResSlopeYAtVtxVsP2->SetTitle("#sigma_{slope_{Y}} at 1st cluster versus p;p (GeV/c);#sigma_{slope_{Y}}");
        gSigmaResSlopeYAtVtxVsP2->SetLineColor(2);
        FitGausResVsMom(hResSlopeYAtVtxVsP, gSigmaResSlopeYAtVtxVsP2, rebinP);
        
      } else {
        
        // recompute momentum resolution at vertex versus p
        TH2F *hResPAtVtxVsPIn02deg = static_cast<TH2F*>(trackerResolutionList->FindObject("hResPAtVtxVsPIn02degMC"));
        TH2F *hResPAtVtxVsPIn23deg = static_cast<TH2F*>(trackerResolutionList->FindObject("hResPAtVtxVsPIn23deg"));
        TH2F *hResPAtVtxVsPIn310deg = static_cast<TH2F*>(trackerResolutionList->FindObject("hResPAtVtxVsPIn310deg"));
        if (!hResPAtVtxVsPIn23deg || !hResPAtVtxVsPIn310deg) return;
        dummy->cd();
        Int_t rebinP = 2;
        gSigmaResPAtVtxVsPIn02deg = new TGraphAsymmErrors(hResPAtVtxVsPIn02deg->GetNbinsX()/rebinP);
        gSigmaResPAtVtxVsPIn02deg->SetName("gSigmaResPAtVtxVsPIn02deg");
        gSigmaResPAtVtxVsPIn02deg->SetTitle("#sigma_{p}/p at vertex versus p in [0,2[ deg;p (GeV/c);#sigma_{p}/p (%)");
        gSigmaResPAtVtxVsPIn02deg->SetLineColor(15);
        FitPResVsP(hResPAtVtxVsPIn02deg, gSigmaResPAtVtxVsPIn02deg, rebinP, 300.);
        gSigmaResPAtVtxVsPIn23deg = new TGraphAsymmErrors(hResPAtVtxVsPIn23deg->GetNbinsX()/rebinP);
        gSigmaResPAtVtxVsPIn23deg->SetName("gSigmaResPAtVtxVsPIn23deg");
        gSigmaResPAtVtxVsPIn23deg->SetTitle("#sigma_{p}/p at vertex versus p in [2,3[ deg;p (GeV/c);#sigma_{p}/p (%)");
        gSigmaResPAtVtxVsPIn23deg->SetLineColor(2);
        FitPResVsP(hResPAtVtxVsPIn23deg, gSigmaResPAtVtxVsPIn23deg, rebinP, 220.);
        gSigmaResPAtVtxVsPIn310deg = new TGraphAsymmErrors(hResPAtVtxVsPIn310deg->GetNbinsX()/rebinP);
        gSigmaResPAtVtxVsPIn310deg->SetName("gSigmaResPAtVtxVsPIn310deg");
        gSigmaResPAtVtxVsPIn310deg->SetTitle("#sigma_{p}/p at vertex versus p in [3,10[ deg;p (GeV/c);#sigma_{p}/p (%)");
        gSigmaResPAtVtxVsPIn310deg->SetLineColor(6);
        FitPResVsP(hResPAtVtxVsPIn310deg, gSigmaResPAtVtxVsPIn310deg, rebinP, 160.);
        
        // recompute slope resolution at vertex versus p
        TH2F *hResSlopeXAtVtxVsP = static_cast<TH2F*>(trackerResolutionList->FindObject("hResSlopeXAtVtxVsP"));
        TH2F *hResSlopeYAtVtxVsP = static_cast<TH2F*>(trackerResolutionList->FindObject("hResSlopeYAtVtxVsP"));
        if (!hResSlopeXAtVtxVsP || !hResSlopeYAtVtxVsP) return;
        dummy->cd();
        gSigmaResSlopeXAtVtxVsP2 = new TGraphAsymmErrors(hResSlopeXAtVtxVsP->GetNbinsX()/rebinP);
        gSigmaResSlopeXAtVtxVsP2->SetName("gSigmaResSlopeXAtVtxVsP2");
        gSigmaResSlopeXAtVtxVsP2->SetTitle("#sigma_{slope_{X}} at vertex versus p;p (GeV/c);#sigma_{slope_{X}}");
        gSigmaResSlopeXAtVtxVsP2->SetLineColor(2);
        FitGausResVsMom(hResSlopeXAtVtxVsP, gSigmaResSlopeXAtVtxVsP2, rebinP);
        gSigmaResSlopeYAtVtxVsP2 = new TGraphAsymmErrors(hResSlopeYAtVtxVsP->GetNbinsX()/rebinP);
        gSigmaResSlopeYAtVtxVsP2->SetName("gSigmaResSlopeYAtVtxVsP2");
        gSigmaResSlopeYAtVtxVsP2->SetTitle("#sigma_{slope_{Y}} at vertex versus p;p (GeV/c);#sigma_{slope_{Y}}");
        gSigmaResSlopeYAtVtxVsP2->SetLineColor(2);
        FitGausResVsMom(hResSlopeYAtVtxVsP, gSigmaResSlopeYAtVtxVsP2, rebinP);
        
      }
      
    }
    
    fResSim->Close();
  }
  
  // draw parameterization of momentum and slope versus p
  TCanvas *cRes = new TCanvas("cRes","cRes",10,10,1200,300);
  cRes->Divide(3,1);
  gROOT->SetSelectedPad(cRes->cd(1));
  gPad->SetLogz();
  if (hPRes) hPRes->Draw("colz");
  if (atFirstCluster && gSigmaResPAtVtxVsP) {
    gSigmaResPAtVtxVsP->Draw(hPRes?"p":"ap");
    if (gSigmaResPAtVtxVsPIn310deg) gSigmaResPAtVtxVsPIn310deg->Draw("p");
  } else if (!atFirstCluster && gSigmaResPAtVtxVsPIn02deg && gSigmaResPAtVtxVsPIn23deg && gSigmaResPAtVtxVsPIn310deg) {
    gSigmaResPAtVtxVsPIn02deg->Draw(hPRes?"p":"ap");
    gSigmaResPAtVtxVsPIn23deg->Draw("p");
    gSigmaResPAtVtxVsPIn310deg->Draw("p");
  }
  if (hPRes || (atFirstCluster && gSigmaResPAtVtxVsP) || (!atFirstCluster && gSigmaResPAtVtxVsPIn02deg && gSigmaResPAtVtxVsPIn23deg && gSigmaResPAtVtxVsPIn310deg))
    fPResVsP02->Draw("same");
  else {
    fPResVsP02->GetYaxis()->SetRangeUser(0., 10.);
    fPResVsP02->Draw();
  }
  fPResVsP23->Draw("same");
  fPResVsP310->Draw("same");
  gROOT->SetSelectedPad(cRes->cd(2));
  gPad->SetLogz();
  if (hSlopeXRes) hSlopeXRes->Draw("colz");
  if (atFirstCluster && gSigmaResSlopeXAtVtxVsP) {
    gSigmaResSlopeXAtVtxVsP->Draw(hPRes?"p":"ap");
    if (gSigmaResSlopeXAtVtxVsP2) gSigmaResSlopeXAtVtxVsP2->Draw("p");
  } else if (!atFirstCluster && gSigmaResSlopeXAtVtxVsP2) gSigmaResSlopeXAtVtxVsP2->Draw(hPRes?"p":"ap");
  if (hSlopeXRes || (atFirstCluster && gSigmaResSlopeXAtVtxVsP) || (!atFirstCluster && gSigmaResSlopeXAtVtxVsP2))
    fSlopeXResVsP02->Draw("same");
  else {
    fSlopeXResVsP02->GetYaxis()->SetRangeUser(0., 0.02);
    fSlopeXResVsP02->Draw();
  }
  fSlopeXResVsP23->Draw("same");
  fSlopeXResVsP310->Draw("same");
  gROOT->SetSelectedPad(cRes->cd(3));
  gPad->SetLogz();
  if (hSlopeYRes) hSlopeYRes->Draw("colz");
  if (atFirstCluster && gSigmaResSlopeYAtVtxVsP) {
    gSigmaResSlopeYAtVtxVsP->Draw(hPRes?"p":"ap");
    if (gSigmaResSlopeYAtVtxVsP2) gSigmaResSlopeYAtVtxVsP2->Draw("p");
  } else if (!atFirstCluster && gSigmaResSlopeYAtVtxVsP2) gSigmaResSlopeYAtVtxVsP2->Draw(hPRes?"p":"ap");
  if (hSlopeYRes || (atFirstCluster && gSigmaResSlopeYAtVtxVsP) || (!atFirstCluster && gSigmaResSlopeYAtVtxVsP2))
    fSlopeYResVsP02->Draw("same");
  else {
    fSlopeYResVsP02->GetYaxis()->SetRangeUser(0., 0.02);
    fSlopeYResVsP02->Draw();
  }
  fSlopeYResVsP23->Draw("same");
  fSlopeYResVsP310->Draw("same");
  
  // draw cluster-track and cluster-trackRef x-residuals
  if (hXResSt_clIn[0] || hXRes_clIn || hXResSt_clOut[0] || hXRes_clOut || hXResSt[0] || hXRes) {
    TCanvas *cResX = new TCanvas("cResX","cResX",10,10,1200,600);
    cResX->Divide(6,3);
    for (Int_t iSt = 0; iSt < 5; ++iSt) {
      gROOT->SetSelectedPad(cResX->cd(iSt+1));
      gPad->SetLogy();
      if (hXResSt_clIn[iSt]) hXResSt_clIn[iSt]->Draw();
      gROOT->SetSelectedPad(cResX->cd(iSt+7));
      gPad->SetLogy();
      if (hXResSt_clOut[iSt]) hXResSt_clOut[iSt]->Draw();
      gROOT->SetSelectedPad(cResX->cd(iSt+13));
      gPad->SetLogy();
      if (hXResSt[iSt]) hXResSt[iSt]->Draw();
    }
    gROOT->SetSelectedPad(cResX->cd(6));
    gPad->SetLogy();
    if (hXRes_clIn) hXRes_clIn->Draw();
    gROOT->SetSelectedPad(cResX->cd(12));
    gPad->SetLogy();
    if (hXRes_clOut) hXRes_clOut->Draw();
    gROOT->SetSelectedPad(cResX->cd(18));
    gPad->SetLogy();
    if (hXRes) hXRes->Draw();
  }
  
  // draw cluster-track and cluster-trackRef y-residuals
  if (hYResSt_clIn[0] || hYRes_clIn || hYResSt_clOut[0] || hYRes_clOut || hYResSt[0] || hYRes) {
    TCanvas *cResY = new TCanvas("cResY","cResY",10,10,1200,600);
    cResY->Divide(6,3);
    for (Int_t iSt = 0; iSt < 5; ++iSt) {
      gROOT->SetSelectedPad(cResY->cd(iSt+1));
      gPad->SetLogy();
      if (hYResSt_clIn[iSt]) hYResSt_clIn[iSt]->Draw();
      gROOT->SetSelectedPad(cResY->cd(iSt+7));
      gPad->SetLogy();
      if (hYResSt_clOut[iSt]) hYResSt_clOut[iSt]->Draw();
      gROOT->SetSelectedPad(cResY->cd(iSt+13));
      gPad->SetLogy();
      if (hYResSt[iSt]) hYResSt[iSt]->Draw();
    }
    gROOT->SetSelectedPad(cResY->cd(6));
    gPad->SetLogy();
    if (hYRes_clIn) hYRes_clIn->Draw();
    gROOT->SetSelectedPad(cResY->cd(12));
    gPad->SetLogy();
    if (hYRes_clOut) hYRes_clOut->Draw();
    gROOT->SetSelectedPad(cResY->cd(18));
    gPad->SetLogy();
    if (hYRes) hYRes->Draw();
  }
  
  // generation functions
  TF1 *fpTGen = new TF1("fpTGen",DnDpT,pTRange[0],pTRange[1],1);
  fpTGen->SetParameter(0,1.);
  fpTGen->SetNpx(10000);
  Double_t pTInt = fpTGen->Integral(pTRange[0],pTRange[1]);
  Double_t sumwpT = pTInt/(pTRange[1]-pTRange[0]);
  TF1 *fEtaGen = new TF1("fEtaGen",DnDeta,etaRange[0],etaRange[1],1);
  fEtaGen->SetParameter(0,1.);
  fEtaGen->SetNpx(10000);
  Double_t etaInt = fEtaGen->Integral(etaRange[0],etaRange[1]);
  Double_t sumwEta = etaInt/(etaRange[1]-etaRange[0]);
  
  // generated and reconstructed histograms
  Int_t npTBins = (Int_t) (npTBinPerGeV*pTRange[1]);
  Double_t pTBinSize = pTRange[1]/npTBins;
  Int_t nEtaBins = (Int_t) (nEtaBinPerUnit*(etaRange[1]-etaRange[0]));
  Double_t etaBinSize = (etaRange[1]-etaRange[0])/nEtaBins;
  Int_t nPhiBins = (Int_t) (nPhiBinPerUnit*(phiRange[1]-phiRange[0]));
  TH1F *hpTGen = new TH1F("hpTGen","hpTGen",2*npTBins,0.,2.*pTRange[1]);
  hpTGen->Sumw2();
  TH1F *hpGen = new TH1F("hpGen","hpGen",5*npTBins,0.,50.*pTRange[1]);
  hpGen->Sumw2();
  TH1F *hetaGen = new TH1F("hetaGen","hetaGen",nEtaBins+0.4*nEtaBinPerUnit,etaRange[0]-0.2,etaRange[1]+0.2);
  hetaGen->Sumw2();
  TH1F *hyGen = new TH1F("hyGen","hyGen",nEtaBins+0.4*nEtaBinPerUnit,etaRange[0]-0.2,etaRange[1]+0.2);
  hyGen->Sumw2();
  TH1F *hphiGen = new TH1F("hphiGen","hphiGen",nPhiBins,phiRange[0],phiRange[1]);
  hphiGen->Sumw2();
  TH1F *hpTRec = new TH1F("hpTRec","hpTRec",2*npTBins,0.,2.*pTRange[1]);
  hpTRec->Sumw2();
  TH1F *hpRec = new TH1F("hpRec","hpRec",5*npTBins,0.,50.*pTRange[1]);
  hpRec->Sumw2();
  TH1F *hetaRec = new TH1F("hetaRec","hetaRec",nEtaBins+0.4*nEtaBinPerUnit,etaRange[0]-0.2,etaRange[1]+0.2);
  hetaRec->Sumw2();
  TH1F *hyRec = new TH1F("hyRec","hyRec",nEtaBins+0.4*nEtaBinPerUnit,etaRange[0]-0.2,etaRange[1]+0.2);
  hyRec->Sumw2();
  TH1F *hphiRec = new TH1F("hphiRec","hphiRec",nPhiBins,phiRange[0],phiRange[1]);
  hphiRec->Sumw2();
  
  // resolution histograms
  Double_t deltaPAtVtxEdges[2];
  if (atFirstCluster) {
    deltaPAtVtxEdges[0] = -5. - 0.0005*(20.*pTRange[1])*(20.*pTRange[1]);
    deltaPAtVtxEdges[1] = 5. + 0.0005*(20.*pTRange[1])*(20.*pTRange[1]);
  } else {
    deltaPAtVtxEdges[0] = -20. - 0.0005*(20.*pTRange[1])*(20.*pTRange[1]);
    deltaPAtVtxEdges[1] = 5. + 0.0005*(20.*pTRange[1])*(20.*pTRange[1]);
  }
  Int_t deltaPAtVtxNBins = (Int_t)(10.*(deltaPAtVtxEdges[1]-deltaPAtVtxEdges[0]));
  if (deltaPAtVtxNBins > 1000) deltaPAtVtxNBins = (deltaPAtVtxNBins+200)/200*200;
  else if (deltaPAtVtxNBins > 100) deltaPAtVtxNBins = (deltaPAtVtxNBins+20)/20*20;
  TH2F *hResPAtVtxVsPIn02degMC = new TH2F("hResPAtVtxVsPIn02deg","#Delta_{p} at vertex versus p in [0,2[ deg MC;p (GeV/c);#Delta_{p} (GeV/c)",2*npTBins,0.,20.*pTRange[1],deltaPAtVtxNBins,deltaPAtVtxEdges[0],deltaPAtVtxEdges[1]);
  hResPAtVtxVsPIn02degMC->SetDirectory(0);
  hResPAtVtxVsPIn02degMC->Sumw2();
  TH2F *hResPAtVtxVsPIn23deg = new TH2F("hResPAtVtxVsPIn23deg","#Delta_{p} at vertex versus p in [2,3[ deg;p (GeV/c);#Delta_{p} (GeV/c)",2*npTBins,0.,20.*pTRange[1],deltaPAtVtxNBins,deltaPAtVtxEdges[0],deltaPAtVtxEdges[1]);
  hResPAtVtxVsPIn23deg->SetDirectory(0);
  hResPAtVtxVsPIn23deg->Sumw2();
  TH2F *hResPAtVtxVsPIn310deg = (atFirstCluster) ? new TH2F("hResPAt1stClVsP","#Delta_{p} at 1st cluster versus p;p (GeV/c);#Delta_{p} (GeV/c)",2*npTBins,0.,20.*pTRange[1],deltaPAtVtxNBins,deltaPAtVtxEdges[0],deltaPAtVtxEdges[1]) : new TH2F("hResPAtVtxVsPIn310deg","#Delta_{p} at vertex versus p in [3,10[ deg;p (GeV/c);#Delta_{p} (GeV/c)",2*npTBins,0.,20.*pTRange[1],deltaPAtVtxNBins,deltaPAtVtxEdges[0],deltaPAtVtxEdges[1]);
  hResPAtVtxVsPIn310deg->SetDirectory(0);
  hResPAtVtxVsPIn310deg->Sumw2();
  TH2F *hResPtAtVtxVsPt = new TH2F("hResPtAtVtxVsPt","#Delta_{p_{t}} at vertex versus p_{t};p_{t} (GeV/c);#Delta_{p_{t}} (GeV/c)",npTBins,0.,pTRange[1],deltaPAtVtxNBins,deltaPAtVtxEdges[0]/10.,deltaPAtVtxEdges[1]/10.);
  hResPtAtVtxVsPt->SetDirectory(0);
  hResPtAtVtxVsPt->Sumw2();
  Double_t deltaSlopeAtVtxEdges[2];
  if (atFirstCluster) {
    deltaSlopeAtVtxEdges[0] = -0.01;
    deltaSlopeAtVtxEdges[1] = 0.01;
  } else {
    deltaSlopeAtVtxEdges[0] = -0.05;
    deltaSlopeAtVtxEdges[1] = 0.05;
  }
  Int_t deltaSlopeAtVtxNBins = 1000;
  TH2F *hResSlopeXAtVtxVsP = (atFirstCluster) ? new TH2F("hResSlopeXAt1stClVsP","#Delta_{slope_{X}} at 1st cluster versus p;p (GeV/c);#Delta_{slope_{X}}",2*npTBins,0.,20.*pTRange[1], deltaSlopeAtVtxNBins, deltaSlopeAtVtxEdges[0], deltaSlopeAtVtxEdges[1]) : new TH2F("hResSlopeXAtVtxVsP","#Delta_{slope_{X}} at vertex versus p;p (GeV/c);#Delta_{slope_{X}}",2*npTBins,0.,20.*pTRange[1], deltaSlopeAtVtxNBins, deltaSlopeAtVtxEdges[0], deltaSlopeAtVtxEdges[1]);
  hResSlopeXAtVtxVsP->SetDirectory(0);
  hResSlopeXAtVtxVsP->Sumw2();
  TH2F *hResSlopeYAtVtxVsP = (atFirstCluster) ? new TH2F("hResSlopeYAt1stVsP","#Delta_{slope_{Y}} at 1st cluster versus p;p (GeV/c);#Delta_{slope_{Y}}",2*npTBins,0.,20.*pTRange[1], deltaSlopeAtVtxNBins, deltaSlopeAtVtxEdges[0], deltaSlopeAtVtxEdges[1]) : new TH2F("hResSlopeYAtVtxVsP","#Delta_{slope_{Y}} at vertex versus p;p (GeV/c);#Delta_{slope_{Y}}",2*npTBins,0.,20.*pTRange[1], deltaSlopeAtVtxNBins, deltaSlopeAtVtxEdges[0], deltaSlopeAtVtxEdges[1]);
  hResSlopeYAtVtxVsP->SetDirectory(0);
  hResSlopeYAtVtxVsP->Sumw2();
  Double_t deltaEtaAtVtxEdges[2] = {-0.5, 0.5};
  Int_t deltaEtaAtVtxNBins = 1000;
  TH2F *hResEtaAtVtxVsP = new TH2F("hResEtaAtVtxVsP","#Delta_{eta} at vertex versus p;p (GeV/c);#Delta_{eta}",2*npTBins,0.,20.*pTRange[1], deltaEtaAtVtxNBins, deltaEtaAtVtxEdges[0], deltaEtaAtVtxEdges[1]);
  hResEtaAtVtxVsP->SetDirectory(0);
  hResEtaAtVtxVsP->Sumw2();
  TH2F *hResPhiAtVtxVsP = new TH2F("hResPhiAtVtxVsP","#Delta_{phi} at vertex versus p;p (GeV/c);#Delta_{phi}",2*npTBins,0.,20.*pTRange[1], deltaEtaAtVtxNBins, deltaEtaAtVtxEdges[0], deltaEtaAtVtxEdges[1]);
  hResPhiAtVtxVsP->SetDirectory(0);
  hResPhiAtVtxVsP->Sumw2();
  
  // loop over events
  Int_t iEvent = 0;
  Double_t muMass = TDatabasePDG::Instance()->GetParticle("mu-")->Mass(); // GeV
  while (++iEvent <= nEvents) {
    
    printf("processing event... %.0f%%%s", 100.*iEvent/nEvents, (iEvent < nEvents) ? "\r" : "\n");
    
    // generate a muon pT/eta/phi
    Double_t pT, eta, phi, w;
    Double_t par[1] = {1.};
    if (generateUniform) {
      pT = gRandom->Uniform(pTRange[0],pTRange[1]);
      eta = gRandom->Uniform(etaRange[0],etaRange[1]);
      w = DnDpT(&pT,par)/sumwpT*DnDeta(&eta,par)/sumwEta;
    } else {
      pT = fpTGen->GetRandom();
      eta = fEtaGen->GetRandom();
      w = 1.;
    }
    phi = gRandom->Uniform(phiRange[0],phiRange[1])*TMath::DegToRad();
    
    // compute eta (y) if generation is done vs y (eta)
    Double_t y;
    Double_t mT = TMath::Sqrt(muMass*muMass + pT*pT);
    if (generateY) {
      y = eta;
      Double_t pZ = mT*TMath::SinH(y);
      eta = TMath::ASinH(pZ/pT);
    } else {
      Double_t pZ = pT*TMath::SinH(eta);
      y = TMath::ASinH(pZ/mT);
    }
    
    // compute total momentum
    Double_t p = pT*TMath::CosH(eta);
    
    // fill generated histograms
    hpGen->Fill(p,w);
    hpTGen->Fill(pT,w);
    hetaGen->Fill(eta,w);
    hyGen->Fill(y,w);
    hphiGen->Fill(phi*TMath::RadToDeg(),w);
    
    // apply acceptance cut
    if (p < pAccCut) continue;
    
    // get energy loss and fluctuation at different angle at this momentum
    Double_t eLoss02 = ELoss(p,1.5); // = 110cm eLoss_Common + 305 cm eLoss_tungsten
    Double_t eLoss23 = ELoss(p,2.5); // = 378cm eLoss_Common + 37 cm eLoss_tungsten
    Double_t eLoss310 = ELoss(p,6.); // = 378cm eLoss_Common + 37 cm eLoss_steel
    Double_t eLossC = (eLoss02 - 305./37.*eLoss23) / (110 - 378.*305./37);
    Double_t eLossW = (eLoss02 - 110.*eLossC) / 305.;
    Double_t eLossS = (eLoss310 - 378.*eLossC) / 37.;
    Double_t fwhmELoss02 = TMath::Sqrt(FWHMELoss2(p, 1.5));
    Double_t fwhmELoss23 = TMath::Sqrt(FWHMELoss2(p, 2.5));
    Double_t fwhmELoss310 = TMath::Sqrt(FWHMELoss2(p, 6.));
    Double_t fwhmELossC = (fwhmELoss02 - 305./37.*fwhmELoss23) / (110 - 378.*305./37);
    Double_t fwhmELossW = (fwhmELoss02 - 110.*fwhmELossC) / 305.;
    Double_t fwhmELossS = (fwhmELoss310 - 378.*fwhmELossC) / 37.;
    
    // get mean energy loss and Branson plane for this eta value
    Double_t theta = 2.*TMath::ATan(TMath::Exp(eta))*TMath::RadToDeg();
    Double_t eLoss, fwhmELoss, zB;
    if (theta < 2.) {
      eLoss = eLoss02;
      fwhmELoss = fwhmELoss02;
      zB = zB02;
    } else if (theta < 3.) {
      eLoss = eLoss23;
      fwhmELoss = fwhmELoss23;
      zB = zB23;
    } else {
      eLoss = eLoss310;
      fwhmELoss = fwhmELoss310;
      zB = zB310;
    }
    
    // MCS in the absorber
    Double_t pHalfLoss = p-0.5*TMath::Min(p-1.e-6,eLoss);
    Double_t slopeX = TMath::Cos(phi)/TMath::SinH(eta);
    Double_t slopeY = TMath::Sin(phi)/TMath::SinH(eta);
    Double_t sigmaMCSAbs = TMath::Sqrt(SigmaSlopeFromMCSInAbs2((theta < 2.) ? p : pHalfLoss, theta));
    Double_t slopeXMCS = GenRndGaus(slopeX,sigmaMCSAbs);
    Double_t slopeYMCS = GenRndGaus(slopeY,sigmaMCSAbs);
    
    // track position at the end of the absorber
    Double_t slopeXAbs, slopeYAbs;
    if (theta < 2) { // 305 cm of tungsten from 200 cm to 505 cm
      slopeXAbs = slopeXMCS + 2./TMath::Sqrt(3.) * (slopeXMCS - slopeX) * zB23 / (zB23 - 2.) * 3.05 / 5.05;
      slopeYAbs = slopeYMCS + 2./TMath::Sqrt(3.) * (slopeYMCS - slopeY) * zB23 / (zB23 - 2.) * 3.05 / 5.05;
    } else {
      slopeXAbs = slopeXMCS + 2./TMath::Sqrt(3.) * (slopeXMCS - slopeX) * zB / (zB - 0.9) * 4.15 / 5.05;
      slopeYAbs = slopeYMCS + 2./TMath::Sqrt(3.) * (slopeYMCS - slopeY) * zB / (zB - 0.9) * 4.15 / 5.05;
    }
    Double_t thetaAbs = TMath::ATan(TMath::Sqrt(slopeXAbs*slopeXAbs + slopeYAbs*slopeYAbs))*TMath::RadToDeg();
    
    // energy loss according to the "real" path in the absorber (should be done in 3D...)
    if ((theta > 2. && thetaAbs < 2.) || (theta > 3. && thetaAbs < 3.)) {
      
      Double_t lAt2out = (thetaAbs < 2.) ? 415.*TMath::Power((2.-theta)/(thetaAbs-theta),2./3.) : 415.;
      Double_t lAt2 = lAt2out*(505.-415.)/(505.-lAt2out);
      Double_t lAt3out = (theta > 3.) ? 415.*TMath::Power((3.-theta)/(thetaAbs-theta),2./3.) : 0.;
      Double_t lAt3 = lAt3out*(505.-415.)/(505.-lAt3out);
      
      Double_t lInC = TMath::Max(TMath::Min(lAt2,378.),110.);
      Double_t lInS = TMath::Max(lAt3-378.,0.);
      
      eLoss = (415.-lInC-lInS)*eLossW + lInS*eLossS + lInC*eLossC;
      fwhmELoss = (415.-lInC-lInS)*fwhmELossW + lInS*fwhmELossS + lInC*fwhmELossC;
      
    } else if ((theta < 2. && thetaAbs > 2.) || (theta < 3. && thetaAbs > 3.)) {
      
      Double_t lInMaterial = (theta < 2.) ? 305. : 415.;
      Double_t lAt2out = (theta < 2.) ? lInMaterial*TMath::Power((2.-theta)/(thetaAbs-theta),2./3.) : 0.;
      Double_t lAt2 = (415.-lInMaterial) + lAt2out*(505.-lInMaterial)/(505.-lAt2out);
      Double_t lAt3out = (thetaAbs > 3.) ? lInMaterial*TMath::Power((3.-theta)/(thetaAbs-theta),2./3.) : lInMaterial;
      Double_t lAt3 = (415.-lInMaterial) + lAt3out*(505.-lInMaterial)/(505.-lAt3out);
      
      Double_t lInC = 110. + TMath::Max(378.-TMath::Max(lAt2,110.),0.);
      Double_t lInS = 415.-TMath::Max(lAt3,378.);
      
      eLoss = (415.-lInC-lInS)*eLossW + lInS*eLossS + lInC*eLossC;
      fwhmELoss = (415.-lInC-lInS)*fwhmELossW + lInS*fwhmELossS + lInC*fwhmELossC;
      
    }
    
    // energy loss in the absorber
    Double_t dp = -1.;
    while (dp < 0) dp = gRandom->Landau(eLoss+0.22278298*0.25*fwhmELoss,0.25*fwhmELoss);
    Double_t pAbsEnd = p - dp;
    if (pAbsEnd < pAccCut) continue;
    
    // Branson plane according to thetaAbs
    Double_t zBRec;
    if (thetaAbs < 2.) zBRec = zB02;
    else if (thetaAbs < 3.) zBRec = zB23;
    else zBRec = zB310;
    
    // compute reconstructed slopes
    Double_t sigmaMCSCh2 = SigmaSlopeFromMCSInCh2(pAbsEnd, kFALSE, zBRec);
    Double_t sigmaResX2 = SigmaSlopeFromRes2(kFALSE, kFALSE, zBRec);
    Double_t sigmaResY2 = SigmaSlopeFromRes2(kTRUE, kFALSE, zBRec);
    Double_t slopeXRec = 0., slopeYRec = 0.;
    if (chosenFunc == kGaus) {
      slopeXRec = GenRndGaus(slopeXMCS,TMath::Sqrt(sigmaMCSCh2+sigmaResX2));
      slopeYRec = GenRndGaus(slopeYMCS,TMath::Sqrt(sigmaMCSCh2+sigmaResY2));
    } else {
      Double_t slopeXMCSCh = GenRndGaus(slopeXMCS,TMath::Sqrt(sigmaMCSCh2));
      Double_t slopeYMCSCh = GenRndGaus(slopeYMCS,TMath::Sqrt(sigmaMCSCh2));
      Double_t sigmaResX = TMath::Sqrt(sigmaResX2);
      Double_t sigmaResY = TMath::Sqrt(sigmaResY2);
      if (chosenFunc == kBreitWigner) {
        slopeXRec = GenRndBreitWigner(slopeXMCSCh, sigmaResX, sigmaTrkCut * sigmaResX / sigmaxChSt1 * sigmaTrk);
        slopeYRec = GenRndBreitWigner(slopeYMCSCh, sigmaResY, sigmaTrkCut * sigmaResY / sigmayChSt1 * sigmaTrk);
      } else if (chosenFunc == kCrystalBall) {
        slopeXRec = GenRndCrystalBall(slopeXMCSCh, sigmaResX, tailxChSt1, sigmaTrkCut * sigmaResX / sigmaxChSt1 * sigmaTrk);
        slopeYRec = GenRndCrystalBall(slopeYMCSCh, sigmaResY, tailyChSt1, sigmaTrkCut * sigmaResY / sigmayChSt1 * sigmaTrk);
      }
    }
    /*
    // track slope dispersion at first cluster
    Double_t sigmaMCSCh2 = SigmaSlopeFromMCSInCh2(pAbsEnd, kFALSE, zBRec);
    Double_t sigmaSlopeX2 = SigmaSlopeFromRes2(kFALSE, kFALSE, zBRec);
    Double_t sigmaSlopeY2 = SigmaSlopeFromRes2(kTRUE, kFALSE, zBRec);
    Double_t dSlopeX = 0., dSlopeY = 0.;
    if (chosenFunc == kGaus) {
      dSlopeX = GenRndGaus(0.,TMath::Sqrt(sigmaMCSCh2+sigmaSlopeX2));
      dSlopeY = GenRndGaus(0.,TMath::Sqrt(sigmaMCSCh2+sigmaSlopeY2));
    } else {
      Double_t sigmaMCSCh = TMath::Sqrt(sigmaMCSCh2);
      Double_t dSlopeXMCS = GenRndGaus(0.,sigmaMCSCh);
      Double_t dSlopeYMCS = GenRndGaus(0.,sigmaMCSCh);
      Double_t sigmaSlopeX = TMath::Sqrt(sigmaSlopeX2);
      Double_t sigmaSlopeY = TMath::Sqrt(sigmaSlopeX2);
      if (chosenFunc == kBreitWigner) {
        dSlopeX = GenRndBreitWigner(dSlopeXMCS, sigmaSlopeX, sigmaTrkCut * sigmaSlopeX / sigmaxChSt1 * sigmaTrk);
        dSlopeY = GenRndBreitWigner(dSlopeYMCS, sigmaSlopeY, sigmaTrkCut * sigmaSlopeY / sigmayChSt1 * sigmaTrk);
      } else if (chosenFunc == kCrystalBall) {
        dSlopeX = GenRndCrystalBall(dSlopeXMCS, sigmaSlopeX, tailxChSt1, sigmaTrkCut * sigmaSlopeX / sigmaxChSt1 * sigmaTrk);
        dSlopeY = GenRndCrystalBall(dSlopeYMCS, sigmaSlopeY, tailyChSt1, sigmaTrkCut * sigmaSlopeY / sigmayChSt1 * sigmaTrk);
      }
    }
    
    // track position dispersion at first cluster
    Double_t sigmaResX = TMath::Sqrt(SigmaPosFromRes2(kFALSE, kFALSE, zBRec));
    Double_t sigmaResY = TMath::Sqrt(SigmaPosFromRes2(kTRUE, kFALSE, zBRec));
    Double_t dPosX = 0., dPosY = 0.;
    if (chosenFunc == kGaus) {
      dPosX = GenRndGaus(0.,sigmaResX);
      dPosY = GenRndGaus(0.,sigmaResY);
    } else if (chosenFunc == kBreitWigner) {
      dPosX = GenRndBreitWigner(0., sigmaResX, sigmaTrkCut * sigmaResX / sigmaxChSt1 * sigmaTrk);
      dPosY = GenRndBreitWigner(0., sigmaResY, sigmaTrkCut * sigmaResY / sigmayChSt1 * sigmaTrk);
    } else if (chosenFunc == kCrystalBall) {
      dPosX = GenRndCrystalBall(0., sigmaResX, tailxChSt1, sigmaTrkCut * sigmaResX / sigmaxChSt1 * sigmaTrk);
      dPosY = GenRndCrystalBall(0., sigmaResY, tailyChSt1, sigmaTrkCut * sigmaResY / sigmayChSt1 * sigmaTrk);
    }
    
    // compute reconstructed slopes at vertex
    Double_t slopeXRec = slopeXMCS + (dSlopeX*(5.36-zBRec) + dPosX)/zBRec;
    Double_t slopeYRec = slopeYMCS + (dSlopeY*(5.36-zBRec) + dPosY)/zBRec;
    */
    // compute reconstructed momentum at first cluster
    Double_t pAbsEndRec = 0.;
    if (chosenFunc == kGaus) {
      Double_t sigmaThetaDev = TMath::Sqrt(SigmaThetaDevFromMCS2(pAbsEnd) + SigmaThetaDevFromRes2());
      pAbsEndRec = ThetaDevToP(GenRndGaus(PToThetaDev(pAbsEnd),sigmaThetaDev));
    } else {
      Double_t thetaDevMCS = GenRndGaus(PToThetaDev(pAbsEnd),TMath::Sqrt(SigmaThetaDevFromMCS2(pAbsEnd)));
      Double_t sigmaThetaRes = TMath::Sqrt(SigmaThetaDevFromRes2());
      if (chosenFunc == kBreitWigner) {
        pAbsEndRec = ThetaDevToP(GenRndBreitWigner(thetaDevMCS, sigmaThetaRes, sigmaTrkCut * sigmaThetaRes / sigmayCh * sigmaTrk));
      } else if (chosenFunc == kCrystalBall) {
        pAbsEndRec = ThetaDevToP(GenRndCrystalBall(thetaDevMCS, sigmaThetaRes, tailyCh, sigmaTrkCut * sigmaThetaRes / sigmayCh * sigmaTrk));
      }
    }
    if (pAbsEndRec < 0.) pAbsEndRec = 1.e6;
    
    // compute reconstructed kinematics
    Double_t phiRec = TMath::ATan2(slopeYRec,slopeXRec) + TMath::Pi();
    Double_t etaRec = -TMath::ASinH(1./TMath::Sqrt(slopeXRec*slopeXRec + slopeYRec*slopeYRec));
    Double_t pRec = pAbsEndRec + ELoss(pAbsEndRec, thetaAbs);
    Double_t pTRec = pRec/TMath::CosH(etaRec);
    Double_t pZRec = pTRec*TMath::SinH(etaRec);
    Double_t mTRec = TMath::Sqrt(muMass*muMass + pTRec*pTRec);
    Double_t yRec = TMath::ASinH(pZRec/mTRec);
    
    // fill resolution histograms
    if (atFirstCluster) {
      hResPAtVtxVsPIn310deg->Fill(pAbsEnd,pAbsEndRec-pAbsEnd,w);
//      hResSlopeXAtVtxVsP->Fill(pAbsEnd,dSlopeX,w);
//      hResSlopeYAtVtxVsP->Fill(pAbsEnd,dSlopeY,w);
    } else {
      if (theta < 2.) hResPAtVtxVsPIn02degMC->Fill(p,pRec-p,w);
      if (thetaAbs > 2. && thetaAbs < 3.) hResPAtVtxVsPIn23deg->Fill(p,pRec-p,w);
      else if (thetaAbs >= 3. && thetaAbs < 10.) hResPAtVtxVsPIn310deg->Fill(p,pRec-p,w);
      hResSlopeXAtVtxVsP->Fill(p,slopeXRec-slopeX,w);
      hResSlopeYAtVtxVsP->Fill(p,slopeYRec-slopeY,w);
    }
    hResPtAtVtxVsPt->Fill(pT,pTRec-pT,w);
    hResEtaAtVtxVsP->Fill(p,etaRec-eta,w);
    Double_t dPhi = phiRec-phi;
    if (dPhi < -TMath::Pi()) dPhi += 2.*TMath::Pi();
    else if (dPhi > TMath::Pi()) dPhi -= 2.*TMath::Pi();
    hResPhiAtVtxVsP->Fill(p,dPhi,w);
    
    // fill reconstructed histograms
    if (thetaAbs < 2. || thetaAbs > 10. || etaRec < -4. || etaRec > -2.5) continue;
    hpRec->Fill(pRec,w);
    hpTRec->Fill(pTRec,w);
    hetaRec->Fill(etaRec,w);
    hyRec->Fill(yRec,w);
    hphiRec->Fill(phiRec*TMath::RadToDeg(),w);
    
  }
  
  // draw generated histograms
  TCanvas *cGen = new TCanvas("cGen","cGen",10,10,600,600);
  cGen->Divide(2,2);
  gROOT->SetSelectedPad(cGen->cd(1));
  gPad->SetLogy();
  hpTGen->Draw();
  fpTGen->SetParameter(0,nEvents*pTBinSize/pTInt);
  fpTGen->Draw("same");
  gROOT->SetSelectedPad(cGen->cd(2));
  gPad->SetLogy();
  hetaGen->Draw();
  hyGen->SetLineColor(9);
  hyGen->Draw("same");
  fEtaGen->SetParameter(0,nEvents*etaBinSize/etaInt);
  fEtaGen->Draw("same");
  gROOT->SetSelectedPad(cGen->cd(3));
  gPad->SetLogy();
  hpGen->Draw();
  gROOT->SetSelectedPad(cGen->cd(4));
  gPad->SetLogy();
  hphiGen->Draw();
  
  // draw reconstructed histograms
  TCanvas *cRec = new TCanvas("cRec","cRec",10,10,600,600);
  cRec->Divide(2,2);
  gROOT->SetSelectedPad(cRec->cd(1));
  gPad->SetLogy();
  hpTRec->Draw();
  gROOT->SetSelectedPad(cRec->cd(2));
  gPad->SetLogy();
  hetaRec->Draw();
  hyRec->SetLineColor(9);
  hyRec->Draw("same");
  gROOT->SetSelectedPad(cRec->cd(3));
  gPad->SetLogy();
  hpRec->Draw();
  gROOT->SetSelectedPad(cRec->cd(4));
  gPad->SetLogy();
  hphiRec->Draw();
  
  // draw reconstructed over generated ratios
  TCanvas *cRat = new TCanvas("cRat","cRat",10,10,600,600);
  cRat->Divide(2,2);
  gROOT->SetSelectedPad(cRat->cd(1));
  TH1 *hpTRat = static_cast<TH1*>(hpTRec->Clone());
  hpTRat->SetNameTitle("hpTRat","hpTRat");
  hpTRat->Divide(hpTGen);
  hpTRat->Draw();
  gROOT->SetSelectedPad(cRat->cd(2));
  TH1 *hetaRat = static_cast<TH1*>(hetaRec->Clone());
  hetaRat->SetNameTitle("hetaRat","hetaRat");
  hetaRat->Divide(hetaGen);
  hetaRat->Draw();
  TH1 *hyRat = static_cast<TH1*>(hyRec->Clone());
  hyRat->SetNameTitle("hyRat","hyRat");
  hyRat->Divide(hyGen);
  hyRat->SetLineColor(9);
  hyRat->Draw("same");
  gROOT->SetSelectedPad(cRat->cd(3));
  TH1 *hpRat = static_cast<TH1*>(hpRec->Clone());
  hpRat->SetNameTitle("hpRat","hpRat");
  hpRat->Divide(hpGen);
  hpRat->Draw();
  gROOT->SetSelectedPad(cRat->cd(4));
  TH1 *hphiRat = static_cast<TH1*>(hphiRec->Clone());
  hphiRat->SetNameTitle("hphiRat","hphiRat");
  hphiRat->Divide(hphiGen);
  hphiRat->Draw();
  
  // draw resolution vs p
  TCanvas *cResVsP = new TCanvas("cResVsP","cResVsP",10,10,1200,600);
  cResVsP->Divide(4,2);
  gROOT->SetSelectedPad(cResVsP->cd(1));
  gPad->SetLogz();
  hResPAtVtxVsPIn02degMC->Draw("colz");
  gROOT->SetSelectedPad(cResVsP->cd(2));
  gPad->SetLogz();
  hResPAtVtxVsPIn23deg->Draw("colz");
  gROOT->SetSelectedPad(cResVsP->cd(3));
  gPad->SetLogz();
  hResPAtVtxVsPIn310deg->Draw("colz");
  gROOT->SetSelectedPad(cResVsP->cd(4));
  gPad->SetLogz();
  hResPtAtVtxVsPt->Draw("colz");
  gROOT->SetSelectedPad(cResVsP->cd(5));
  gPad->SetLogz();
  hResSlopeXAtVtxVsP->Draw("colz");
  gROOT->SetSelectedPad(cResVsP->cd(6));
  gPad->SetLogz();
  hResSlopeYAtVtxVsP->Draw("colz");
  gROOT->SetSelectedPad(cResVsP->cd(7));
  gPad->SetLogz();
  hResEtaAtVtxVsP->Draw("colz");
  gROOT->SetSelectedPad(cResVsP->cd(8));
  gPad->SetLogz();
  hResPhiAtVtxVsP->Draw("colz");
  
  if (checkResVsP) {
    
    TGraphAsymmErrors *gSigmaResPAtVtxVsPIn02degMC2 = 0x0, *gSigmaResPAtVtxVsPIn23deg2 = 0x0, *gSigmaResPAtVtxVsPIn310deg2 = 0x0, *gSigmaResSlopeXAtVtxVsP3 = 0x0, *gSigmaResSlopeYAtVtxVsP3 = 0x0;
    
    if (atFirstCluster) {
      
      // compute momentum resolution versus p
      Int_t rebinP = 2;
      dummy->cd();
      gSigmaResPAtVtxVsPIn310deg2 = new TGraphAsymmErrors(hResPAtVtxVsPIn310deg->GetNbinsX()/rebinP);
      gSigmaResPAtVtxVsPIn310deg2->SetName("gSigmaResPAt1stClVsP2");
      gSigmaResPAtVtxVsPIn310deg2->SetTitle("#sigma_{p}/p at 1st cluster versus p;p (GeV/c);#sigma_{p}/p (%)");
      gSigmaResPAtVtxVsPIn310deg2->SetLineColor(9);
      FitPResVsP(hResPAtVtxVsPIn310deg, gSigmaResPAtVtxVsPIn310deg2, rebinP, -1.);
      
      // compute slope resolution versus p
      dummy->cd();
      gSigmaResSlopeXAtVtxVsP3 = new TGraphAsymmErrors(hResSlopeXAtVtxVsP->GetNbinsX()/rebinP);
      gSigmaResSlopeXAtVtxVsP3->SetName("gSigmaResSlopeXAt1stClVsP3");
      gSigmaResSlopeXAtVtxVsP3->SetTitle("#sigma_{slope_{X}} at 1st cluster versus p;p (GeV/c);#sigma_{slope_{X}}");
      gSigmaResSlopeXAtVtxVsP3->SetLineColor(4);
      FitGausResVsMom(hResSlopeXAtVtxVsP, gSigmaResSlopeXAtVtxVsP3, rebinP);
      gSigmaResSlopeYAtVtxVsP3 = new TGraphAsymmErrors(hResSlopeYAtVtxVsP->GetNbinsX()/rebinP);
      gSigmaResSlopeYAtVtxVsP3->SetName("gSigmaResSlopeYAt1stClVsP3");
      gSigmaResSlopeYAtVtxVsP3->SetTitle("#sigma_{slope_{Y}} at 1st cluster versus p;p (GeV/c);#sigma_{slope_{Y}}");
      gSigmaResSlopeYAtVtxVsP3->SetLineColor(4);
      FitGausResVsMom(hResSlopeYAtVtxVsP, gSigmaResSlopeYAtVtxVsP3, rebinP);
      
    } else {
      
      // compute momentum resolution versus p
      Int_t rebinP = 2;
      dummy->cd();
      gSigmaResPAtVtxVsPIn02degMC2 = new TGraphAsymmErrors(hResPAtVtxVsPIn02degMC->GetNbinsX()/rebinP);
      gSigmaResPAtVtxVsPIn02degMC2->SetName("gSigmaResPAtVtxVsPIn02degMC2");
      gSigmaResPAtVtxVsPIn02degMC2->SetTitle("#sigma_{p}/p at vertex versus p in [0,2[ deg MC;p (GeV/c);#sigma_{p}/p (%)");
      gSigmaResPAtVtxVsPIn02degMC2->SetLineColor(11);
      FitPResVsP(hResPAtVtxVsPIn02degMC, gSigmaResPAtVtxVsPIn02degMC2, rebinP, 300.);
      gSigmaResPAtVtxVsPIn23deg2 = new TGraphAsymmErrors(hResPAtVtxVsPIn23deg->GetNbinsX()/rebinP);
      gSigmaResPAtVtxVsPIn23deg2->SetName("gSigmaResPAtVtxVsPIn23deg2");
      gSigmaResPAtVtxVsPIn23deg2->SetTitle("#sigma_{p}/p at vertex versus p in [2,3[ deg;p (GeV/c);#sigma_{p}/p (%)");
      gSigmaResPAtVtxVsPIn23deg2->SetLineColor(4);
      FitPResVsP(hResPAtVtxVsPIn23deg, gSigmaResPAtVtxVsPIn23deg2, rebinP, 220.);
      gSigmaResPAtVtxVsPIn310deg2 = new TGraphAsymmErrors(hResPAtVtxVsPIn310deg->GetNbinsX()/rebinP);
      gSigmaResPAtVtxVsPIn310deg2->SetName("gSigmaResPAtVtxVsPIn310deg2");
      gSigmaResPAtVtxVsPIn310deg2->SetTitle("#sigma_{p}/p at vertex versus p in [3,10[ deg;p (GeV/c);#sigma_{p}/p (%)");
      gSigmaResPAtVtxVsPIn310deg2->SetLineColor(9);
      FitPResVsP(hResPAtVtxVsPIn310deg, gSigmaResPAtVtxVsPIn310deg2, rebinP, 160.);
      
      // compute slope resolution versus p
      dummy->cd();
      gSigmaResSlopeXAtVtxVsP3 = new TGraphAsymmErrors(hResSlopeXAtVtxVsP->GetNbinsX()/rebinP);
      gSigmaResSlopeXAtVtxVsP3->SetName("gSigmaResSlopeXAtVtxVsP3");
      gSigmaResSlopeXAtVtxVsP3->SetTitle("#sigma_{slope_{X}} at vertex versus p;p (GeV/c);#sigma_{slope_{X}}");
      gSigmaResSlopeXAtVtxVsP3->SetLineColor(4);
      FitGausResVsMom(hResSlopeXAtVtxVsP, gSigmaResSlopeXAtVtxVsP3, rebinP);
      gSigmaResSlopeYAtVtxVsP3 = new TGraphAsymmErrors(hResSlopeYAtVtxVsP->GetNbinsX()/rebinP);
      gSigmaResSlopeYAtVtxVsP3->SetName("gSigmaResSlopeYAtVtxVsP3");
      gSigmaResSlopeYAtVtxVsP3->SetTitle("#sigma_{slope_{Y}} at vertex versus p;p (GeV/c);#sigma_{slope_{Y}}");
      gSigmaResSlopeYAtVtxVsP3->SetLineColor(4);
      FitGausResVsMom(hResSlopeYAtVtxVsP, gSigmaResSlopeYAtVtxVsP3, rebinP);
      
    }
    
    // draw resolution vs p
    gROOT->SetSelectedPad(cRes->cd(1));
    if (gSigmaResPAtVtxVsPIn02degMC2) gSigmaResPAtVtxVsPIn02degMC2->Draw("p");
    if (gSigmaResPAtVtxVsPIn23deg2) gSigmaResPAtVtxVsPIn23deg2->Draw("p");
    gSigmaResPAtVtxVsPIn310deg2->Draw("p");
    gROOT->SetSelectedPad(cRes->cd(2));
    if (gSigmaResSlopeXAtVtxVsP3) gSigmaResSlopeXAtVtxVsP3->Draw("p");
    gROOT->SetSelectedPad(cRes->cd(3));
    if (gSigmaResSlopeYAtVtxVsP3) gSigmaResSlopeYAtVtxVsP3->Draw("p");
    
  }
  
  // save results
  TFile *fResults = TFile::Open("ResToYield.root","RECREATE");
  hpGen->Write();
  hpTGen->Write();
  hetaGen->Write();
  hyGen->Write();
  hphiGen->Write();
  hpRec->Write();
  hpTRec->Write();
  hetaRec->Write();
  hyRec->Write();
  hphiRec->Write();
  fResults->Close();
  
}
/*
//-----------------------------------------------------------------------
Double_t DnDpT( const Double_t *x, const Double_t *par )
{
  // muon pT
  Double_t pT = *x;
  Double_t p[6] = {371.665, 0.845642, 0.56192, 9.34859, 0.000474519, -0.851091}; // LHC13de tune1
//  Double_t p[6] = {1., 0., 0., 0., 0., 0.}; // flat
  return par[0] * p[0] * (1. / TMath::Power(p[1] + TMath::Power(pT,p[2]), p[3]) + p[4] * TMath::Exp(p[5]*pT));
}

//-----------------------------------------------------------------------
Double_t DnDeta( const Double_t *x, const Double_t *par )
{
  // muon eta
  Double_t eta = *x;
  Double_t p[8] = {0.777922, 1., 0., -0.0184202, 0., -0.00107081, 0., 2.}; // LHC13de tune1
//  Double_t p[8] = {1., 1., 0., 0., 0., 0., 0., 1.}; // flat
  Double_t arg = eta/p[7];
  return par[0] * p[0] * (p[1] * (1. + p[2]*eta + p[3]*eta*eta + p[4]*eta*eta*eta + p[5]*eta*eta*eta*eta) + p[6]*TMath::Exp(-0.5*arg*arg));
}
*/
//-----------------------------------------------------------------------
Double_t DnDpT( const Double_t *px, const Double_t *par )
{
  // muon pT FONLL
  
  Double_t x=*px;
  if( x < 1. ) return 0.;
  Float_t p0,p1,p2,p3, p4, p5, p6, p7;
  p0 =  4.13770e+02; //
  p1 =  9.99639e-02; //
  p2 = -1.75012e-02; //
  p3 = -1.79602e-01; //
  p4 =  3.84611e+00; //
  p5 =  2.69159e+05; //
  p6 =  1.07241e+04; //
  p7 =  3.34421e+01; //
  
  // norm
  Float_t norm = p0;
  //slope
  Float_t slope = p1 * ( 1. - TMath::Exp( p2*x)) + p3;
  //double
  Float_t den = 1./TMath::Power(x,p4);
  // pol
  Float_t pol = p5 + p6*x+p7*x*x;
  
  return par[0] * norm * TMath::Exp( slope*x ) * den * pol;
}

//-----------------------------------------------------------------------
Double_t DnDeta( const Double_t *py, const Double_t *par )
{
  // muon y FONLL
  
  Float_t yswitch = -1.;
  Double_t x = *py;
  x *= yswitch;
  //pol4 only valid in y= -4;-2.5
  Float_t p0,p1,p2,p3, p4, p5, p6, p7, p8;
  p0 =  1.28729e+06; //
  p1 = -1.54430e+05; //
  p2 = -4.31467e+04; //
  p3 = -4.74633e+03; //
  p4 =  2.24207e+02; //
  p5 =  2.54699e+02; //
  p6 =  4.89765e+01; //
  p7 =  8.83499e+00; //
  p8 = -3.15046e+00; //
  return par[0] *
  (p0 +
   p1* x +
   p2* x*x +
   p3* x*x*x +
   p4* x*x*x*x +
   p5* x*x*x*x*x +
   p6* x*x*x*x*x*x +
   p7* x*x*x*x*x*x*x +
   p8* x*x*x*x*x*x*x*x);
}

//-----------------------------------------------------------------------
Double_t PResVsP( const Double_t *x, const Double_t *par )
{
  // momentum resolution versus p
  Double_t p = *x;
  Double_t theta = par[0];
  Double_t dp = ELoss(p, theta);
  Double_t pCorr = TMath::Max(atFirstCluster ? p : p - dp, pAccCut);
  Double_t thetaDev = PToThetaDev(pCorr);
  Double_t sigmaThetaDev2 = SigmaThetaDevFromRes2();
  if (chosenFunc == kBreitWigner) sigmaThetaDev2 /= 2.*log(2.); // FWHM = 2 * gamma = 2 * srqt(2*ln(2)) * sigma
  sigmaThetaDev2 += SigmaThetaDevFromMCS2(pCorr);
  if (atFirstCluster) return 100.*TMath::Sqrt(sigmaThetaDev2)/thetaDev;
  else {
    Double_t sigmaELoss2 = FWHMELoss2(p, theta) / (8.*log(2.)); // gaussian: fwmh = 2 * srqt(2*ln(2)) * sigma
    return 100.*TMath::Sqrt(sigmaELoss2 + sigmaThetaDev2/thetaDev/thetaDev*pCorr*pCorr)/p;
  }
}

//-----------------------------------------------------------------------
Double_t SlopeResVsP( const Double_t *x, const Double_t *par )
{
  // slope resolution versus p
  Double_t p = *x;
  Double_t theta = par[0];
  Double_t dp = ELoss(p, theta);
  Double_t zB;
  if (theta < 2.) zB = zB02;
  else if (theta < 3.) zB = zB23;
  else zB = zB310;
  Double_t pCorr = TMath::Max(atFirstCluster ? p : p - dp, pAccCut);
  Double_t sigmaMCSCh2 = SigmaSlopeFromMCSInCh2(pCorr, atFirstCluster, zB);
  Double_t sigmaRes2 = SigmaSlopeFromRes2((par[1] > 0), atFirstCluster, zB);
  if (chosenFunc == kBreitWigner) sigmaRes2 /= 2.*log(2.); // FWHM = 2 * gamma = 2 * srqt(2*ln(2)) * sigma
  if (atFirstCluster) return TMath::Sqrt(sigmaMCSCh2+sigmaRes2);
  else {
    Double_t pHalfCorr = p-0.5*TMath::Min(p-1.e-6,dp);
    Double_t sigmaMCSAbs2 = SigmaSlopeFromMCSInAbs2(pHalfCorr,theta);
    return TMath::Sqrt(sigmaMCSAbs2+sigmaMCSCh2+sigmaRes2);
  }
}

//-----------------------------------------------------------------------
Double_t PToThetaDev(Double_t p)
{
  /// deviation angle for a track of a given momentum in the spectrometer
  return 3./p*0.3; // BL = 3Tm; GeV/c = 10^9/3.10^8
}

//-----------------------------------------------------------------------
Double_t ThetaDevToP(Double_t thetaDev)
{
  /// deviation angle for a track of a given momentum in the spectrometer
  return 3./thetaDev*0.3; // BL = 3Tm; GeV/c = 10^9/3.10^8
}

//-----------------------------------------------------------------------
Double_t FWHMELoss2(Double_t p, Double_t theta)
{
  // energy loss dispersion of muon in the absorber
  Double_t rhoZoverA; // effective rho * Z / A of the absorber (Carbone = 2.21*0.5; Be = 0.44*1.85)
  if (tuneKalman) {
    Double_t corr = 2.; // FWHM are summed quadratically instead of linearly in AliMUONTrackExtrap --> bug to be corrected
    if (theta < 2.) rhoZoverA = 3.3*corr;
    else if (theta < 3.) rhoZoverA = 0.81*corr;
    else rhoZoverA = 0.62*corr;
  } else {
    if (theta < 2.) rhoZoverA = 11.;
    else if (theta < 3.) rhoZoverA = 2.7;
    else rhoZoverA = 1.7;
  }
  return ELossFluctuation2(p, rhoZoverA);
}

//-----------------------------------------------------------------------
Double_t SigmaThetaDevFromMCS2(Double_t p)
{
  /// resolution of the deviation angle for a track of a given momentum in the spectrometer
  return tuneKalman ? 4.*MCS2(p, 0.035, 1.) : 2*MCS2(p, 0.075, 1.) + 4.*MCS2(p, 0.035, 1.);
}

//-----------------------------------------------------------------------
Double_t SigmaThetaDevFromRes2()
{
  /// resolution of the deviation angle for a track of a given momentum in the spectrometer
  static Double_t magicFactor = 1./2.2; // to get thetaDev resolution^2 from averaged chamber resolution^2
  return magicFactor*sigmayCh*sigmayCh;
}

//-----------------------------------------------------------------------
Double_t SigmaSlopeFromMCSInAbs2(Double_t p, Double_t theta)
{
  /// slope dispersion due to the MCS in the front absorber
  Double_t x0;
//  if (theta < 2.) x0 = (tuneKalman) ? 1. : 0.5;
//  else x0 = (tuneKalman) ? 35.4 : 19.3; // X0 Carbone = 42.7/2.21; Be = 65.2/1.85
  if (theta < 2.) x0 = 1.;
  else x0 = 35.4; // X0 Carbone = 42.7/2.21; Be = 65.2/1.85
  return MCS2(p, 415., x0) / 4.; // sigma reduced by ~2 by the Branson correction
}

//-----------------------------------------------------------------------
Double_t SigmaSlopeFromMCSInCh2(Double_t p, Bool_t at1stCl, Double_t zB)
{
  /// slope dispersion due to the MCS in chambers
  Double_t sigmaMCS2 = (tuneKalman) ? 2.*MCS2(p, 0.065, 1.) : 2.*MCS2(p, 0.065, 1.) + 2.*MCS2(p, 0.075, 1.);
  if (at1stCl) return sigmaMCS2;
  else return sigmaMCS2*(5.36-zB)*(5.36-zB)/zB/zB; // propagate to the Branson plane then to the vertex
}
/*
//-----------------------------------------------------------------------
Double_t SigmaSlopeFromMCSInCh2(Double_t p, Bool_t at1stCl, Double_t zB)
{
  /// slope dispersion due to the MCS in chambers
  return (tuneKalman) ? 2.*MCS2(p, 0.065, 1.) : 2.*MCS2(p, 0.065, 1.) + 2.*MCS2(p, 0.075, 1.);
}
*/
//-----------------------------------------------------------------------
Double_t SigmaSlopeFromRes2(Bool_t bendingDir, Bool_t at1stCl, Double_t zB)
{
  /// slope dispersion due to chamber resolution
  static Double_t magicFactorX = 1./23.; // to get angular resolution^2 from station 1 x-resolution^2
  static Double_t magicFactorY = 1./3.; // to get angular resolution^2 from station 1 y-resolution^2
  Double_t sigmaSt12, sigmaSlope2;
  if (bendingDir) {
    sigmaSt12 = sigmayChSt1*sigmayChSt1/2.;
    sigmaSlope2 = magicFactorY*sigmaSt12;
  } else {
    sigmaSt12 = sigmaxChSt1*sigmaxChSt1/2.;
    sigmaSlope2 = magicFactorX*sigmaSt12;
  }
  if (at1stCl) return sigmaSlope2;
  else return (sigmaSlope2*(5.36-zB)*(5.36-zB) + sigmaSt12)/zB/zB; // propagate to the Branson plane then to the vertex
}
/*
//-----------------------------------------------------------------------
Double_t SigmaPosFromRes2(Bool_t bendingDir, Bool_t at1stCl, Double_t zB)
{
  /// slope dispersion due to chamber resolution
  if (bendingDir) return sigmayChSt1*sigmayChSt1/2.;
  else return sigmaxChSt1*sigmaxChSt1/2.;
}

//-----------------------------------------------------------------------
Double_t SigmaSlopeFromRes2(Bool_t bendingDir, Bool_t at1stCl, Double_t zB)
{
  /// slope dispersion due to chamber resolution
  static Double_t magicFactorX = 1./23.; // to get angular resolution^2 from station 1 x-resolution^2
  static Double_t magicFactorY = 1./3.; // to get angular resolution^2 from station 1 y-resolution^2
  if (bendingDir) return magicFactorY*sigmayChSt1*sigmayChSt1/2.;
  else return magicFactorX*sigmaxChSt1*sigmaxChSt1/2.;
}
*/
//-----------------------------------------------------------------------
Double_t ELoss(Double_t p, Double_t theta)
{
  // Returns the total momentum energy loss of muon in the absorber
  if (theta < 2.) return 9.3 + 0.05*p;
  else if (theta < 3.) return 2.9 + 0.0134*p;
  else return 2.4 + 0.007*p;
}

//-----------------------------------------------------------------------
Double_t ELossFluctuation2(Double_t p, Double_t rhoZoverA)
{
  // Returns the total momentum energy loss fluctuation of muon in the absorber
  Double_t muMass = TDatabasePDG::Instance()->GetParticle("mu-")->Mass(); // GeV
  Double_t k = 0.307075e-3; // GeV.g^-1.cm^2
  Double_t pathLength = 415.;
  Double_t p2=p*p;
  Double_t beta2=p2/(p2 + muMass*muMass);
  Double_t fwhm = 2. * k * rhoZoverA * pathLength / beta2; // FWHM of the energy loss Landau distribution
  return fwhm*fwhm;
}

//-----------------------------------------------------------------------
Double_t MCS2(Double_t p, Double_t dZ, Double_t x0)
{
  /// Return the angular dispersion square due to multiple Coulomb scattering
  /// through a material of thickness "dZ" and of radiation length "x0"
  Double_t muMass = TDatabasePDG::Instance()->GetParticle("mu-")->Mass(); // GeV
  Double_t p2=p*p;
  Double_t betap=p2/TMath::Sqrt(p2 + muMass*muMass);
  Double_t theta02 = 0.0136 / betap * (1 + 0.038 * TMath::Log(dZ/x0));
  return theta02 * theta02 * dZ / x0;
}

//-----------------------------------------------------------------------
Double_t GetSigma(TH1 *h, Double_t *sigmaErr)
{
  /// get the dispersion of the histo
  switch (chosenFunc) {
    case kCrystalBall:
      return GetSigmaCrystalBall(h, sigmaErr);
      break;
    case kBreitWigner:
      return GetSigmaBreitWigner(h, sigmaErr);
      break;
    case kGaus:
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
  
  static TF1 *fCrystalBall = new TF1("CrystalBall",CrystalBallSymmetric,-2.5,2.5,5);
  
  if (h->GetEntries() < 10.) {
    if (sigmaErr) *sigmaErr = 0.;
    return 0.;
  }
  
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
  Double_t s0 = PToThetaDev(p0);
  Double_t s = PToThetaDev(p);
  Double_t dsdp = s/p;
  
  return par[0] * s/p * TMath::Gaus(s, s0, s0/p0*par[2], kTRUE);
}

//-----------------------------------------------------------------------
void FitPResVsP(TH2 *h, TGraphAsymmErrors *gSigma, Int_t rebinP, Double_t pSwitch)
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

//-----------------------------------------------------------------------
Double_t CrystalBallSymmetric(Double_t *x,Double_t *par)
{
  ///par[0] = Normalization
  ///par[1] = mean
  ///par[2] = sigma
  ///par[3] = alpha = alpha'
  ///par[4] = n = n'
  
  Double_t t = fabs((x[0]-par[1])/par[2]);
  
  Double_t absAlpha = fabs(par[3]);
  Double_t a = TMath::Power(par[4]/absAlpha,par[4])*exp(-0.5*absAlpha*absAlpha);
  Double_t b = par[4]/absAlpha - absAlpha;
  
  if (t < absAlpha) return par[0]*(exp(-0.5*t*t)); // gaussian core
  else return par[0]*(a/TMath::Power(b + t, par[4])); //left and right tails
  
}

//-----------------------------------------------------------------------
Double_t GenRndGaus(Double_t mean, Double_t sigma)
{
  /// generate a random number following a gaussian distribution
  return gRandom->Gaus(mean,sigma);
}

//-----------------------------------------------------------------------
Double_t GenRndBreitWigner(Double_t mean, Double_t sigma, Double_t max)
{
  /// generate a random number following a Breit-Wigner distribution in the range mean ± max
  Double_t val = 0.;
  do val = gRandom->BreitWigner(mean,sigma);
  while (TMath::Abs(val-mean) > max);
  return val;
}

//-----------------------------------------------------------------------
Double_t GenRndCrystalBall(Double_t mean, Double_t sigma, Double_t tail[2], Double_t max)
{
  /// generate a random number following a Crystal Ball distribution in the range mean ± max
  static TF1 *fCrystalBall = 0x0;
  if (!fCrystalBall) {
    fCrystalBall = new TF1("CrystalBall2",CrystalBallSymmetric,-2.,2.,5);
    fCrystalBall->SetNpx(1000);
  }
  fCrystalBall->SetParameters(1.,mean,sigma,tail[0],tail[1]);
  fCrystalBall->SetRange(mean-max,mean+max);
  return fCrystalBall->GetRandom();
}

