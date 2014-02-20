/*
 *  JPsiInPbPb2010_v6.C
 *  aliroot_dev
 *
 *  Created by philippe pillot on 20/09/11.
 *  Copyright 2011 SUBATECH. All rights reserved.
 *
 */

#if !defined(__CINT__) || defined(__MAKECINT__)

#include <Riostream.h>
#include "TCanvas.h"
#include "TH1.h"
#include "TMath.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TString.h"
#include "TASImage.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TLine.h"
#include "TROOT.h"

#endif

// y=-1: loop over the 3 choices below
// y=0: -4<y<-2.5
// y=1: -4<y<-3.25
// y=2: -3.25<y<-2

// pt=-1: loop over the 3 choices below
// pt=0: pt>0
// pt=1: 0<pt<3GeV/c
// pt=2: pt>3GeV/c

// branching ratio
Double_t BR=0.059;
Double_t eBR=0.01; // 1%

// luminosity in p-p
Double_t eppLumi=0.043; // 3.1% + 3% du R-factor = 4.3% 
//Double_t eppLumi=0.0;  // if already included in the pp result

Bool_t fineDisplay = kFALSE;

Bool_t PRL = kTRUE;
//Bool_t PRL = kFALSE;


Float_t gCanvasWidth = 1024;
Float_t gCanvasHeigh = 768;

Int_t gMarkerStylePhenForw = 22; // 22 triangle plein pointe haut
Float_t gMarkerSizePhenForw = 1.9;
Int_t gMarkerStylePhenMid = 23; // 23 triangle plein pointe bas
Float_t gMarkerSizePhenMid = 1.8;
Int_t gMarkerStyleAliceForw = 21; // 21 carre plein,  25 carre vide
Float_t gMarkerSizeAliceForw = 1.9;

void SetGraphNameTitleVsCent(Int_t nCentBins, TString *&gName, TString *&gTitle);
void SetGraphLabelsVsCent(Int_t nCentBins, Double_t *&x, Double_t *&ex, Double_t *&bin, TString *&label,
			  Double_t *&NPart, Double_t *&eNPart, Double_t *&dNchdEta, Double_t *&edNchdEta,
			  Double_t *&dNchdEtaNorm, Double_t *&edNchdEtaNorm);
void SetNormVar(Int_t pt, Int_t y, Double_t &sigmaJPsipp, Double_t &esigmaJPsipp);
void SetNormVarVsCent(Int_t pt, Int_t y, Int_t nCentBins, Double_t *&nMB, Double_t *&TAA, Double_t *&eTAA, Double_t *&AccEff, Double_t *&eAccEff);
void SetSystVarVsCent(Int_t nCentBins, Double_t *&eGen, Double_t *&eRec, Double_t *&eTrk, Double_t *&eTrg);
void SetNJPsiVsCent(Double_t **&nJPsi, Double_t **&enJPsi, Int_t &nTests, Int_t &nCentBins);
void SetNJPsiVsCenty1(Double_t **&nJPsi, Double_t **&enJPsi, Int_t &nTests, Int_t &nCentBins);
void SetNJPsiVsCenty2(Double_t **&nJPsi, Double_t **&enJPsi, Int_t &nTests, Int_t &nCentBins);
void SetNJPsiVsCentpt1(Double_t **&nJPsi, Double_t **&enJPsi, Int_t &nTests, Int_t &nCentBins);
void SetNJPsiVsCentpt2(Double_t **&nJPsi, Double_t **&enJPsi, Int_t &nTests, Int_t &nCentBins);
void SetNJPsiVsCentpt2y1(Double_t **&nJPsi, Double_t **&enJPsi, Int_t &nTests, Int_t &nCentBins);
void SetNJPsiVsCentpt2y2(Double_t **&nJPsi, Double_t **&enJPsi, Int_t &nTests, Int_t &nCentBins);
void ComputeNJpsiMean(Int_t pt, Int_t y, Double_t *&nJPsiMean, Double_t *&enJPsiMean, Double_t *&sysnJPsiMean, Int_t &nCentBins, Bool_t weightedMeanRMS);
void ComputeRAA(Int_t pt, Int_t y, Int_t nCentBins, Double_t *nJPsi, Double_t *enJPsi, Double_t *sysnJPsi, Double_t *&renJPsi,
		Double_t *&RAA, Double_t *&eRAA, Double_t *&eRAA2, Double_t &eRAACorr, Bool_t print = kFALSE);
void BuildRAAGraphs(Int_t nCentBins, Double_t *renJPsi, Double_t *RAA, Double_t *eRAA, Double_t *eRAA2, Double_t eRAACorr,
		    TGraphErrors *&gRAA, TGraphErrors *&gRAASys, TGraphErrors *&gRAASysAll, Int_t xAxisType, Int_t color = 2, Int_t sysPos = 0);
void BuildForwardPHENIXRAAGraphs(TGraphErrors *&gRAA_PHENIX, TGraphAsymmErrors *&gRAA_PHENIXSys,
				 TGraphErrors *&gRAA_PHENIXSysAll, Int_t xAxisType);
void BuildMidPHENIXRAAGraphs(TGraphErrors *&gRAA_PHENIX2, TGraphErrors *&gRAA_PHENIXSys2,
			     TGraphErrors *&gRAA_PHENIXSysAll2, Int_t xAxisType);
void DisplayRAA(Int_t pt, Int_t y, Int_t nCentBins, TGraphErrors *gRAA, TGraphErrors *gRAASys, TGraphErrors *gRAASysAll, Int_t xAxisType);
void DisplayRAA(TString var, Int_t nCentBins, TGraphErrors *gRAA[3], TGraphErrors *gRAASys[3], TGraphErrors *gRAASysAll[3], Int_t xAxisType);
void DisplayRAAWithPHENIX(Int_t pt, Int_t y, TGraphErrors *gRAA, TGraphErrors *gRAASys, TGraphErrors *gRAASysAll, Int_t xAxisType);
void DisplayRAAWithPHENIX2(Int_t pt, Int_t y, TGraphErrors *gRAA, TGraphErrors *gRAASys, TGraphErrors *gRAASysAll, Int_t xAxisType);
void ALICEseal(TString type, Double_t xPad, Double_t yPad);


//------------------------------------------------------------------------
void JPsiInPbPb2010_v6(Bool_t weightedMeanRMS = kTRUE, Int_t xAxisType = 1, Int_t y = -1, Int_t pt = -1)
{
  
  if (y < -1 || y > 2 || pt < -1 || pt > 2 || (y > 0 && pt == 1)) {
    ::Error("JPsiInPbPb2010_v6", "\ninvalid choice of pt or rapidity bin\n\n");
    exit(0);
  }
  
  Int_t font=42;
  gStyle->SetCanvasColor(10);
  gStyle->SetFrameFillColor(10);
  gStyle->SetStatColor(10);
  gStyle->SetFillColor(10);
  gStyle->SetOptTitle(0);
  gStyle->SetTitleFillColor(10);
  gStyle->SetTextFont(font);
  gStyle->SetLabelFont(font,"x");
  gStyle->SetTitleFont(font,"x");
  gStyle->SetLabelFont(font,"y");
  gStyle->SetTitleFont(font,"y");
  gROOT->SetStyle("Plain");
  gROOT->ForceStyle();
  Int_t color[3] = {2, 4, kGreen+2};
  
  if (y == -1) {
    
    Int_t nCentBins;
    TGraphErrors *gRAA[3], *gRAASys[3], *gRAASysAll[3];
    
    for (Int_t i=0; i<3; i++) {
      
      Double_t *nJPsiMean, *enJPsiMean, *sysnJPsiMean;
      ComputeNJpsiMean(0, i, nJPsiMean, enJPsiMean, sysnJPsiMean, nCentBins, weightedMeanRMS);
      
      Double_t eRAACorr;
      Double_t *renJPsiMean, *RAA, *eRAA, *eRAA2;
      ComputeRAA(0, i, nCentBins, nJPsiMean, enJPsiMean, sysnJPsiMean, renJPsiMean, RAA, eRAA, eRAA2, eRAACorr, kTRUE);
      
      BuildRAAGraphs(nCentBins, renJPsiMean, RAA, eRAA, eRAA2, eRAACorr, gRAA[i], gRAASys[i], gRAASysAll[i], xAxisType, color[i], i);
      
    }
    
    DisplayRAA("y", nCentBins, gRAA, gRAASys, gRAASysAll, xAxisType);
    
  }
  
  if (pt == -1) {
    
    Int_t nCentBins;
    TGraphErrors *gRAA[3], *gRAASys[3], *gRAASysAll[3];
    
    for (Int_t i=0; i<3; i++) {
      
      Double_t *nJPsiMean, *enJPsiMean, *sysnJPsiMean;
      ComputeNJpsiMean(i, 0, nJPsiMean, enJPsiMean, sysnJPsiMean, nCentBins, weightedMeanRMS);
      
      Double_t eRAACorr;
      Double_t *renJPsiMean, *RAA, *eRAA, *eRAA2;
      ComputeRAA(i, 0, nCentBins, nJPsiMean, enJPsiMean, sysnJPsiMean, renJPsiMean, RAA, eRAA, eRAA2, eRAACorr, kTRUE);
      
      BuildRAAGraphs(nCentBins, renJPsiMean, RAA, eRAA, eRAA2, eRAACorr, gRAA[i], gRAASys[i], gRAASysAll[i], xAxisType, color[i], i);
      
    }
    
    DisplayRAA("pt", nCentBins, gRAA, gRAASys, gRAASysAll, xAxisType);
    
  }
  
  if (y == -1) {
    
    Int_t nCentBins;
    
    for (Int_t i=1; i<3; i++) {
      
      Double_t *nJPsiMean, *enJPsiMean, *sysnJPsiMean;
      ComputeNJpsiMean(2, i, nJPsiMean, enJPsiMean, sysnJPsiMean, nCentBins, weightedMeanRMS);
      
      Double_t eRAACorr;
      Double_t *renJPsiMean, *RAA, *eRAA, *eRAA2;
      ComputeRAA(2, i, nCentBins, nJPsiMean, enJPsiMean, sysnJPsiMean, renJPsiMean, RAA, eRAA, eRAA2, eRAACorr, kTRUE);
      
    }
    
  }
  
  if (y<0) y = 0;
  if (pt<0) pt = 0;
  
  Int_t nCentBins;
  Double_t *nJPsiMean, *enJPsiMean, *sysnJPsiMean;
  ComputeNJpsiMean(pt, y, nJPsiMean, enJPsiMean, sysnJPsiMean, nCentBins, weightedMeanRMS);
  
  Double_t eRAACorr;
  Double_t *renJPsiMean, *RAA, *eRAA, *eRAA2;
  ComputeRAA(pt, y, nCentBins, nJPsiMean, enJPsiMean, sysnJPsiMean, renJPsiMean, RAA, eRAA, eRAA2, eRAACorr, kTRUE);
  
  TGraphErrors *gRAA, *gRAASys, *gRAASysAll;
  BuildRAAGraphs(nCentBins, renJPsiMean, RAA, eRAA, eRAA2, eRAACorr, gRAA, gRAASys, gRAASysAll, xAxisType);
  
  DisplayRAA(pt, y, nCentBins, gRAA, gRAASys, gRAASysAll, xAxisType);
  
  DisplayRAAWithPHENIX(pt, y, gRAA, gRAASys, gRAASysAll, xAxisType);
  
  DisplayRAAWithPHENIX2(pt, y, gRAA, gRAASys, gRAASysAll, xAxisType);
  
}


//------------------------------------------------------------------------
void SetGraphNameTitleVsCent(Int_t nCentBins, TString *&gName, TString *&gTitle)
{
  // define graph name and title according to the number of centrality bins
  
  if (nCentBins == 1) {
    
    gName = new TString[1];
    gName[0] = "080";
    gTitle = new TString[1];
    gTitle[0] = "0-80%";
    
  } else if (nCentBins == 4) {
    
    gName = new TString[4];
    gName[0] = "080"; gName[1] = "020"; gName[2] = "2040"; gName[3] = "4080";
    gTitle = new TString[4];
    gTitle[0] = "0-80%"; gTitle[1] = "0-20%"; gTitle[2] = "20-40%"; gTitle[3] = "40-80%";
    
  } else if (nCentBins == 6) {
    
    gName = new TString[6];
    gName[0] = "080"; gName[1] = "010"; gName[2] = "1020"; gName[3] = "2030"; gName[4] = "3050"; gName[5] = "5080";
    gTitle = new TString[6];
    gTitle[0] = "0-80%"; gTitle[1] = "0-10%"; gTitle[2] = "10-20%"; gTitle[3] = "20-30%"; gTitle[4] = "30-50%"; gTitle[5] = "50-80%";
    
  } else {
    
    ::Error("SetGraphNameTitleVsCent", "invalid number of centrality bins");
    exit(0);
    
  }
  
}


//------------------------------------------------------------------------
void SetGraphLabelsVsCent(Int_t nCentBins, Double_t *&x, Double_t *&ex, Double_t *&bin, TString *&label,
			  Double_t *&NPart, Double_t *&eNPart, Double_t *&dNchdEta, Double_t *&edNchdEta,
			  Double_t *&dNchdEtaNorm, Double_t *&edNchdEtaNorm)
{
  // define graph x-axis labels according to the number of centrality bins
  
  if (nCentBins == 1) {
    
    // plot versus centrality bin
    x = new Double_t[1];
    x[0] = 100.;
    ex = new Double_t[1];
    ex[0] = 1.;
    bin = new Double_t[2];
    bin[0] = 0.; bin[1] = 100.;
    label = new TString[2];
    label[0] = "0-80%"; label[1] = "0-80%";
    
    // plot versus NPart
    NPart = new Double_t[1]; // weighted using Alberica's macro
    NPart[0] = 264.0;
    eNPart = new Double_t[1];
    eNPart[0] = 4.5; 
    
    // plot versus dNch/deta
    dNchdEta = new Double_t[1]; 
    dNchdEta[0] =  1000;
    // put a 0% error for now
    edNchdEta = new Double_t[1];
    Double_t errPC = 0.;
    edNchdEta[0] = dNchdEta[0]*errPC;
    
    // plot versus dNch/deta /Npart/2
    dNchdEtaNorm =new Double_t[1]; 
    dNchdEtaNorm[0] =  7.42;
    // put a 0% error for now
    edNchdEtaNorm = new Double_t[1];
    errPC = 0.;
    edNchdEtaNorm[0] = dNchdEtaNorm[0]*errPC;
    
  } else if (nCentBins == 4) {
    
    // plot versus centrality bin
    x = new Double_t[4];
    x[0] = 100.; x[1] = 10.; x[2] = 30.; x[3] = 60.;
    ex = new Double_t[4];
    ex[0] = 1.; ex[1] = 10.; ex[2] = 10.; ex[3] = 20.;
    bin = new Double_t[5];
    bin[0] = 0.; bin[1] = 20.; bin[2] = 60.; bin[3] = 80.; bin[4] = 100.;
    label = new TString[5];
    label[0] = "0-80%"; label[1] = "0-20%"; label[2] = "20-40%"; label[3] = "40-80%"; label[4] = "80-100%";
    
    // plot versus NPart
    NPart = new Double_t[4]; // weighted using Alberica's macro
    NPart[0] = 264.0; NPart[1] = 323.6; NPart[2] = 168.2; NPart[3] = 69.6;
    eNPart = new Double_t[4];
    eNPart[0] = 4.5; eNPart[1] = 3.6; eNPart[2] = 3.2; eNPart[3] = 2.8;
    
    // plot versus dNch/deta
    dNchdEta =new Double_t[4]; 
    dNchdEta[0] = 1000; dNchdEta[1] = 1250; dNchdEta[2] = 580; dNchdEta[3] = 200; 
    // put a 0% error for now
    edNchdEta = new Double_t[4];
    Double_t errPC = 0.;
    for (Int_t i = 1; i < 4; i++) edNchdEta[i] = dNchdEta[i]*errPC;
    
    // plot versus dNch/deta /Npart/2
    dNchdEtaNorm =new Double_t[4]; 
    dNchdEtaNorm[0] =  7.42; dNchdEtaNorm[1] = 7.78; dNchdEtaNorm[2] = 6.86; dNchdEtaNorm[3] = 5.9;     
    // put a 0% error for now
    edNchdEtaNorm = new Double_t[4];
    errPC = 0.;
    for (Int_t i = 1; i < 4; i++) edNchdEtaNorm[i] = dNchdEtaNorm[i]*errPC;
    
  } else if (nCentBins == 6) {
    
    // plot versus centrality bin
    x = new Double_t[6];
    x[0] = 100.; x[1] = 5.; x[2] = 15.; x[3] = 25.; x[4] = 40.; x[5] = 65.;
    ex = new Double_t[6];
    ex[0] = 1.; ex[1] = 5.; ex[2] = 5.; ex[3] = 5.; ex[4] = 10.; ex[5] = 15.;
    bin = new Double_t[7];
    bin[0] = 0.; bin[1] = 20.; bin[2] = 50.; bin[3] = 70.; bin[4] = 80.; bin[5] = 90.; bin[6] = 100.;
    label = new TString[7];
    label[0] = "0-80%"; label[1] = "0-10%"; label[2] = "10-20%"; label[3] = "20-30%"; label[4] = "30-50%"; label[5] = "50-80%"; label[6] = "80-100%";
    
    // plot versus NPart
    NPart = new Double_t[6]; // weighted using Alberica's macro
    NPart[0] = 264.0; NPart[1] = 360.7; NPart[2] = 264.0; NPart[3] = 189.4; NPart[4] = 117.3; NPart[5] = 46.5;
    eNPart = new Double_t[6];
    eNPart[0] = 4.5; eNPart[1] = 3.6; eNPart[2] = 3.7; eNPart[3] = 3.5; eNPart[4] = 3.4; eNPart[5] = 1.7;
    
    // plot versus dNch/deta
    dNchdEta =new Double_t[6]; 
    dNchdEta[0] = 1000; dNchdEta[1] = 1463; dNchdEta[2] = 979; dNchdEta[3] = 658; dNchdEta[4] = 369; dNchdEta[5] = 110; 
    //dNchdEta[0] = 1000; dNchdEta[1] = 1326.; dNchdEta[2] = 966.; dNchdEta[3] = 649.; dNchdEta[4] = 127.; dNchdEta[5] = 51.;
    // put a 0% error for now
    edNchdEta = new Double_t[6];
    Double_t errPC = 0.;
    edNchdEta[0] = 50.;  edNchdEta[1] = 60.; edNchdEta[2] = 37.; edNchdEta[3] = 23.; edNchdEta[4] = 13.; edNchdEta[5] = 5.; 
    
    
    // plot versus dNch/deta /Npart/2
    dNchdEtaNorm =new Double_t[6];
    dNchdEtaNorm[0] = 7.42; dNchdEtaNorm[1] = 8.1; dNchdEtaNorm[2] = 7.6; dNchdEtaNorm[3] = 7.0; dNchdEtaNorm[4] = 6.5; dNchdEtaNorm[5] = 5.5;
    // put a 0% error for now
    edNchdEtaNorm = new Double_t[6];
    errPC = 0.;
    for (Int_t i = 1; i < 6; i++) edNchdEtaNorm[i] = dNchdEtaNorm[i]*errPC;
    
  } else {
    
    ::Error("SetGraphLabelsVsCent", "invalid number of centrality bins");
    exit(0);
    
  }
  
}


//------------------------------------------------------------------------
void SetNormVar(Int_t pt, Int_t y, Double_t &sigmaJPsipp, Double_t &esigmaJPsipp)
{
  // set the normalization variables according to rapidity or pt bin
  
  Double_t ppStatErr = 0.;
  Double_t ppSysErr = 0.;
  
  //(syst. uncertainty is evaluated from all the contributions, including Luminosity. Polarization is not included.):
  //luminosity is 3.1%.
  
  if (y == 0 && pt == 0) { // all pt, y
    
    // old --- 3.46 +- 0.13(stat) +- 0.32(syst) +- 0.28(syst. lumi) mub
    // measured: sigma_J/psi sigma(all pT, 2.5<y<4):   3.24 +- 0.13(stat) +- 0.32(syst) mub
    // Cross Section(1M) = 3.343 +/- 0.132 (Stat) +/- 0.261 (SystTot not incl.Lumi) +/- 0.104(Lumi) +/- 0.281 (TotSystErr) (mub)
    sigmaJPsipp = 3.343;
    ppStatErr = 0.132;
    ppSysErr = TMath::Sqrt(0.281*0.281-eppLumi*sigmaJPsipp*eppLumi*sigmaJPsipp);
    
  } else if (y == 1 && pt == 0) { // all pt, -4<y<-3.25 
    
    // old --- 1.50+-0.08(stat)+-0.13(syst)+-0.12(syst. lumi) mub
    // sigma(all pT, 3.25<y<4): 1.38 +- 0.07(stat) +- 0.13(syst) mub
    // Cross Section(1M) = 1.423 +/- 0.071 (Stat) +/- 0.111 (SystTot not incl.Lumi) +/- 0.044(Lumi) +/- 0.120 (TotSystErr)(mub)
    sigmaJPsipp = 1.423;
    ppStatErr = 0.071;
    ppSysErr = TMath::Sqrt(0.12*0.12-eppLumi*sigmaJPsipp*eppLumi*sigmaJPsipp);
    
  } else if (y == 2 && pt == 0) { // all pt, -3.25<y<-2
    
    // old --- 1.89+-0.11(stat)+-0.17(syst)+-0.15(syst. lumi) mub
    // sigma(all pT, 2.5<y<3.25): 1.86 +- 0.10(stat) +- 0.17(syst) mub
    // Cross Section(1M) = 1.922 +/- 0.108 (Stat) +/- 0.150 (SystTot not incl.Lumi) +/- 0.060(Lumi) +/- 0.162 (TotSystErr)(mub)
    sigmaJPsipp = 1.922;
    ppStatErr = 0.108;
    ppSysErr = TMath::Sqrt(0.162*0.162-eppLumi*sigmaJPsipp*eppLumi*sigmaJPsipp);
    
  } else if (y == 0 && pt == 1) { // pt<3GeV/c
    
    // sigma (pT<3GeV, 2.5<y<4) = 2.63+-0.11(stat)+-0.24(syst)+-0.21(syst. lumi) mub
    sigmaJPsipp = 2.63;
    ppStatErr = 0.11;
    ppSysErr = TMath::Sqrt(0.24*0.24-eppLumi*sigmaJPsipp*eppLumi*sigmaJPsipp);
    
  } else if (y == 0 && pt == 2) { // pt>3GeV/c
    
    // sigma (pT>3 GeV, 2.5<y<4) = 0.87+-0.06(stat)+-0.09(syst)+-0.07(syst. lumi) mub
    sigmaJPsipp = 0.87;
    ppStatErr = 0.06;
    ppSysErr = TMath::Sqrt(0.09*0.09-eppLumi*sigmaJPsipp*eppLumi*sigmaJPsipp);
    
  } else if (y == 1 && pt == 2) { // pt>3GeV/c and -4<y<-3.25
    
    // old --- 0.337+-0.030(Stat)+-0.031(Syst)+-0.027(Lumi) mub
    // sigma(pT>3, 3.25<y<4): 0.33 +- 0.03(stat) +- 0.03(syst) mub 
    // Cross Section(1M) = 0.344 +/- 0.029 (Stat) +/- 0.027 (SystTot not incl.Lumi) +/- 0.011(Lumi) +/- 0.029 (TotSystErr) (mub)
    sigmaJPsipp = 0.344;
    ppStatErr = 0.029;
    ppSysErr =  TMath::Sqrt(0.029*0.029-eppLumi*sigmaJPsipp*eppLumi*sigmaJPsipp);
    
  } else if (y == 2 && pt == 2) { // pt>3GeV/c and -3.25<y<-2.5
    
    // old --- 0.485+-0.042(Stat)+-0.044(Syst)+-0.039(Lumi) mub
    // sigma(pT>3, 2.5<y<3.25): 0.48 +- 0.04(stat) +- 0.05(syst) mub
    // Cross Section(1M) = 0.496 +/- 0.042 (Stat) +/- 0.039 (SystTot not incl.Lumi) +/- 0.015(Lumi) +/- 0.042 (TotSystErr) (mub)
    sigmaJPsipp = 0.496;
    ppStatErr = 0.042;
    ppSysErr = TMath::Sqrt(0.042*0.042-eppLumi*sigmaJPsipp*eppLumi*sigmaJPsipp);
    
  }
  
  // relative error**2 on p-p cross-section (without uncertainty from luminosity)
  Double_t esigmaJPsipp2 = (ppStatErr*ppStatErr + ppSysErr*ppSysErr) / (sigmaJPsipp*sigmaJPsipp);
  
  // remove relative error on the branching ratio
  esigmaJPsipp2 -= eBR*eBR;
  
  // relative error on p-p cross-section
  esigmaJPsipp = TMath::Sqrt(esigmaJPsipp2);
  
}


//------------------------------------------------------------------------
void SetNormVarVsCent(Int_t pt, Int_t y, Int_t nCentBins, Double_t *&nMB, Double_t *&TAA, Double_t *&eTAA, Double_t *&AccEff, Double_t *&eAccEff)
{
  // set the normalization variables according to the number of centrality bins
  
  if (nCentBins == 1) {
    
    // Christophe pass2
    nMB = new Double_t[1];
    nMB[0] = 17650876.;
    
    TAA = new Double_t[1]; // mb^-1
    TAA[0] = 7.04;
    eTAA = new Double_t[1];
    eTAA[0] = 0.27;
    
    AccEff = new Double_t[1];
    eAccEff = new Double_t[1];
    if (y == 0 && pt == 0) { // all pt, y      
      //AccEff = 0.1877;
      //eAccEff = 0.0005/AccEff;
      AccEff[0] = 0.176436;
      eAccEff[0] = 0.000746;
    } else if (y == 1 && pt == 0) { // all pt, -4<y<-3.25 
      //AccEff = 0.186;
      //eAccEff = 0.0007/AccEff;
      AccEff[0] = 0.170522; 
      eAccEff[0] = 0.001121;
    } else if (y == 2 && pt == 0) { // all pt, -3.25<y<-2.5
      //AccEff = 0.1888;
      //eAccEff = 0.0006/AccEff;
      AccEff[0] = 0.180863;
      eAccEff[0] = 0.000999;
    } else if (y == 0 && pt == 1) { // pt<3GeV/c
      //AccEff = 0.1751;
      //eAccEff = 0.0005/AccEff;
      AccEff[0] = 0.164124;
      eAccEff[0] =  0.0008471;
    } else if (y == 0 && pt == 2) { // pt>3GeV/c
      //AccEff = 0.2204;
      //eAccEff = 0.0009/AccEff;
      AccEff[0] =  0.208509;
      eAccEff[0] =  0.001541;
    } else if (y == 1 && pt == 2) { // pt>3GeV/c and -4<y<-3.25
      //AccEff = 0.2131;
      //eAccEff = 0.0014/AccEff;
      AccEff[0] = 0.197948;
      eAccEff[0] =   0.002297;
    } else if (y == 2 && pt == 2) { // pt>3GeV/c and -3.25<y<-2.5
      //AccEff = 0.2258;
      //eAccEff = 0.0012/AccEff;
      AccEff[0] =  0.216383;
      eAccEff[0] = 0.002075;
    }
    
  } else if (nCentBins == 4) {
    
    // Christophe pass2
    nMB = new Double_t[4];
    nMB[0] = 17650876.; nMB[1] = 4400805.; nMB[2] = 4415524.5; nMB[3] = 8834546.5;
    
    TAA = new Double_t[4]; // mb^-1
    TAA[0] = 7.04; TAA[1] = 18.93; TAA[2] = 6.86; TAA[3] = 1.2;
    eTAA = new Double_t[4];
    eTAA[0] = 0.27; eTAA[1] = 0.74; eTAA[2] = 0.28; eTAA[3] = 0.07;
    
    AccEff = new Double_t[4];
    eAccEff = new Double_t[4];
    if (y == 0 && pt == 0) { // all pt, y   
      AccEff[0] = 0.176436;      AccEff[1] = 1;      AccEff[2] = 1;      AccEff[3] = 1;
      eAccEff[0] = 0.000746;      eAccEff[1] = 1;      eAccEff[2] = 1;      eAccEff[3] = 1;
    } else if (y == 1 && pt == 0) { // all pt, -4<y<-3.25 
      AccEff[0] = 0.170522;      AccEff[1] = 1;      AccEff[2] = 1;      AccEff[3] = 1;
      eAccEff[0] = 0.001121;      eAccEff[1] = 1;      eAccEff[2] = 1;      eAccEff[3] = 1;
    } else if (y == 2 && pt == 0) { // all pt, -3.25<y<-2.5
      AccEff[0] = 0.180863;      AccEff[1] = 1;      AccEff[2] = 1;      AccEff[3] = 1;
      eAccEff[0] = 0.000999;      eAccEff[1] = 1;      eAccEff[2] = 1;      eAccEff[3] = 1;
    } else if (y == 0 && pt == 1) { // pt<3GeV/c
      AccEff[0] = 0.164124;      AccEff[1] = 1;      AccEff[2] = 1;      AccEff[3] = 1;
      eAccEff[0] = 0.0008471;      eAccEff[1] = 1;      eAccEff[2] = 1;      eAccEff[3] = 1;
    } else if (y == 0 && pt == 2) { // pt>3GeV/c
      AccEff[0] = 0.208509; AccEff[1] = 0.206747;  AccEff[2] = 0.211615;  AccEff[3] = 0.216682;
      eAccEff[0] = 0.001541;  eAccEff[1] = 0.001155;  eAccEff[2] = 0.001166;      eAccEff[3] = 0.000833;
    } else if (y == 1 && pt == 2) { // pt>3GeV/c and -4<y<-3.25
      AccEff[0] = 0.197948;      AccEff[1] = 1;      AccEff[2] = 1;      AccEff[3] = 1;
      eAccEff[0] = 0.002297;      eAccEff[1] = 1;      eAccEff[2] = 1;      eAccEff[3] = 1;
    } else if (y == 2 && pt == 2) { // pt>3GeV/c and -3.25<y<-2.5
      AccEff[0] = 0.216383;      AccEff[1] = 1;      AccEff[2] = 1;      AccEff[3] = 1;
      eAccEff[0] = 0.002075;      eAccEff[1] = 1;      eAccEff[2] = 1;      eAccEff[3] = 1;
    }
    
  } else if (nCentBins == 6) {
    
    // Christophe pass2
    nMB = new Double_t[6];
    nMB[0] = 17650876.; nMB[1] = 2204316.; nMB[2] = 2196489.; nMB[3] = 2204114.; nMB[4] = 4422821.; nMB[5] = 6623136.;
    
    TAA = new Double_t[6]; // mb^-1
    TAA[0] = 7.04; TAA[1] = 23.48; TAA[2] = 14.43; TAA[3] = 8.74; TAA[4] = 3.87; TAA[5] = 0.72;
    eTAA = new Double_t[6];
    eTAA[0] = 0.27; eTAA[1] = 0.97; eTAA[2] = 0.57; eTAA[3] = 0.37; eTAA[4] = 0.18; eTAA[5] = 0.05;
    
    AccEff = new Double_t[6];
    eAccEff = new Double_t[6];
    if (y == 0 && pt == 0) { // all pt, y    pt>0GeV/c    -4<y<-2.5 
      AccEff[0] = 0.176436; AccEff[1] = 0.171510;   AccEff[2] = 0.177077;    AccEff[3] = 0.179810;  AccEff[4] = 0.182709;   AccEff[5] = 0.185980;  
      eAccEff[0] = 0.000746;  eAccEff[1] = 0.000783;   eAccEff[2] = 0.000797;  eAccEff[3] = 0.000802;  eAccEff[4] =0.000570;  eAccEff[5] = 0.000470;
    } else if (y == 1 && pt == 0) { // all pt, -4<y<-3.25 
      AccEff[0] =0.170522 ;  AccEff[1] = 0.164712;  AccEff[2] = 0.171817;  AccEff[3] = 0.173653;  AccEff[4] = 0.178297;  AccEff[5] = 0.181119;
      eAccEff[0] = 0.001121;  eAccEff[1] = 0.001172; eAccEff[2] = 0.001202; eAccEff[3] = 0.001202;  eAccEff[4] = 0.000861;  eAccEff[5] = 0.000708;
    } else if (y == 2 && pt == 0) { // all pt, -3.25<y<-2.5
      AccEff[0] = 0.180863;   AccEff[1] = 0.176598;      AccEff[2] = 0.180986;    AccEff[3] = 0.184464;  AccEff[4] = 0.186008;  AccEff[5] = 0.189620;
      eAccEff[0] = 0.000999;  eAccEff[1] = 0.001050;  eAccEff[2] = 0.001063;  eAccEff[3] = 0.001077;  eAccEff[4] = 0.000761;  eAccEff[5] = 0.000627;
    } else if (y == 0 && pt == 1) { // pt<3GeV/c  -4<y<-2.5 
      AccEff[0] = 0.164124;  AccEff[1] = 0.159009;  AccEff[2] = 0.164658;  AccEff[3] = 0.167625;  AccEff[4] = 0.170913;  AccEff[5] = 0.173978;
      eAccEff[0] = 0.000847;  eAccEff[1] = 0.000886;  eAccEff[2] = 0.000904;  eAccEff[3] = 0.000911;  eAccEff[4] = 0.000649;  eAccEff[5] = 0.000535;
    } else if (y == 0 && pt == 2) { // pt>3GeV/c  -4<y<-2.5 
      AccEff[0] =  0.208509;  AccEff[1] = 0.204061;  AccEff[2] = 0.209445;  AccEff[3] = 0.211624;  AccEff[4] = 0.213393;  AccEff[5] = 0.217178;
      eAccEff[0] = 0.001541;  eAccEff[1] = 0.001620;  eAccEff[2] = 0.001645;  eAccEff[3] = 0.001654;  eAccEff[4] = 0.001169;  eAccEff[5] = 0.000963;
    } else if (y == 1 && pt == 2) { // pt>3GeV/c and -4<y<-3.25
      AccEff[0] = 0.197948;  AccEff[1] = 0.192132;  AccEff[2] = 0.200183;  AccEff[3] = 0.198991;  AccEff[4] = 0.206409;  AccEff[5] = 0.209419;
      eAccEff[0] = 0.002297;  eAccEff[1] = 0.002403;  eAccEff[2] = 0.002470;  eAccEff[3] = 0.002452;  eAccEff[4] = 0.001758;  eAccEff[5] = 0.001449;
    } else if (y == 2 && pt == 2) { // pt>3GeV/c and -3.25<y<-2.5
      AccEff[0] = 0.216383;  AccEff[1] = 0.212990;  AccEff[2] = 0.216265;  AccEff[3] = 0.221077;  AccEff[4] = 0.218613;  AccEff[5] = 0.222938;
      eAccEff[0] = 0.002075;  eAccEff[1] = 0.002189; eAccEff[2] = 0.002203;  eAccEff[3] = 0.002235;  eAccEff[4] = 0.001565;  eAccEff[5] = 0.001288;
    }
    
  } else {
    
    ::Error("SetNormVarVsCent", "invalid number of centrality bins");
    exit(0);
    
  }
  
}


//------------------------------------------------------------------------
void SetSystVarVsCent(Int_t nCentBins, Double_t *&eGen, Double_t *&eRec, Double_t *&eTrk, Double_t *&eTrg)
{
  // set the sytematic errors (in percent) according to the number of centrality bins
  
  if (nCentBins == 1) { // 0-80, all
    
    eGen = new Double_t[2];
    eGen[0] = 0.; eGen[1] = 0.05;
    eRec = new Double_t[2];
    eRec[0] = 0.; eRec[1] = 0.02;
    eTrk = new Double_t[2];
    eTrk[0] = 0.005; eTrk[1] = 0.05;
    eTrg = new Double_t[2];
    eTrg[0] = 0.015; eTrg[1] = 0.04;
    
  } else if (nCentBins == 4) { // 0-80, 0-20, 20-40, 40-80, all
    
    eGen = new Double_t[5];
    eGen[0] = 0.; eGen[1] = 0.; eGen[2] = 0.; eGen[3] = 0.; eGen[4] = 0.05;
    eRec = new Double_t[5];
    eRec[0] = 0.; eRec[1] = 0.; eRec[2] = 0.; eRec[3] = 0.; eRec[4] = 0.02;
    eTrk = new Double_t[5];
    eTrk[0] = 0.005; eTrk[1] = 0.01; eTrk[2] = 0.005; eTrk[3] = 0., eTrk[4] = 0.05;
    eTrg = new Double_t[5];
    eTrg[0] = 0.015; eTrg[1] = 0.02; eTrg[2] = 0.01; eTrg[3] = 0.; eTrg[4] = 0.04;
    
  } else if (nCentBins == 6) { // 0-80, 0-10, 10-20, 20-30, 30-50, 50-80, all
    
    eGen = new Double_t[7];
    eGen[0] = 0.; eGen[1] = 0.; eGen[2] = 0.; eGen[3] = 0.; eGen[4] = 0.; eGen[5] = 0.; eGen[6] = 0.05;
    eRec = new Double_t[7];
    eRec[0] = 0.; eRec[1] = 0.; eRec[2] = 0.; eRec[3] = 0.; eRec[4] = 0.; eRec[5] = 0.; eRec[6] = 0.02;
    eTrk = new Double_t[7];
    eTrk[0] = 0.005; eTrk[1] = 0.01; eTrk[2] = 0.005; eTrk[3] = 0.; eTrk[4] = 0.; eTrk[5] = 0.; eTrk[6] = 0.05;
    eTrg = new Double_t[7];
    eTrg[0] = 0.015; eTrg[1] = 0.02; eTrg[2] = 0.015; eTrg[3] = 0.01; eTrg[4] = 0.005; eTrg[5] = 0.; eTrg[6] = 0.04;
    
  } else {
    
    ::Error("SetSystVarVsCent", "invalid number of centrality bins");
    exit(0);
    
  }
  
}


//------------------------------------------------------------------------
void SetNJPsiVsCent(Double_t **&nJPsi, Double_t **&enJPsi, Int_t &nTests, Int_t &nCentBins)
{
  // return the number of JPsi and errors for each test and each centrality bins
  // as well as the number of tests and of centrality bins
  
  const Int_t nT = 24;
  const Int_t nC = 6;
  
  Double_t n[nT][nC] = {
    
    // ------ CB1 ------
    
    // CB tail fixed to simulations (pure J/Œ®) in all cases (Œ±=0.98 n=5.2) M and œÉ fixed to 0-80% fit
    {2186.5, 870.356, 538.051, 353.358, 298.849, 117.837},
    // CB tail fixed to simulations (pure J/Œ®) in all cases (Œ±=0.98 n=5.2) M and œÉ fixed to 0-80% fit: œÉ + 1 sigma
    {2322.79, 923.378, 563.304, 387.22, 316.586, 121.781},
    // CB tail fixed to simulations (pure J/Œ®) in all cases (Œ±=0.98 n=5.2) M and œÉ fixed to 0-80% fit: œÉ - 1 sigma
    {2037.66, 812.948, 509.08, 318.194, 278.79, 112.881},
    /*
     // CB tail fixed to simulations (pure J/Œ®) in all cases (Œ±=0.98 n=5.2) M free and œÉ fixed (74 instead of 76) to 0-80% fit
     {2148.17, 858.215, 534.377, 353.892, 293.84, 119.819},
     // CB tail fixed to simulations (pure J/Œ®) in all cases (Œ±=0.98 n=5.2) M and œÉ fixed to 0-80% fit: M + 1 sigma
     {2130.62, 846.493, 523.07, 349.381, 290.515, 113.487},
     // CB tail fixed to simulations (pure J/Œ®) in all cases (Œ±=0.98 n=5.2) M and œÉ fixed to 0-80% fit: M - 1 sigma
     {2158.51, 861.681, 535.576, 338.947, 295.812, 118.912},
     */
    // CB tail fixed to pp LHC10g in all cases (Œ±=1.15 n=1.59) M and œÉ fix to 0-80%
    {2423.41, 963.731, 595.516, 391.634, 333.495, 130.448},
    // CB tail fixed to pp LHC10g in all cases (Œ±=1.15 n=1.59) M and œÉ fix to 0-80% + 1 sigma
    {2556.59, 1015.76, 618.779, 426.221, 350.372, 134.072},
    // CB tail fixed to pp LHC10g in all cases (Œ±=1.15 n=1.59) M and œÉ fix to 0-80% - 1 sigma
    {2274.16, 906.158, 567.631, 355.242, 313.489, 125.633},
    
    // CB tail fixed to pp ‚Äúold‚Äù MC in all cases (Œ±=1.15 n=3.6) M and œÉ fix to 0-80%
    {2190.66, 871.788, 536.225, 355.211, 300.749, 118.163},
    // CB tail fixed to pp ‚Äúold‚Äù MC in all cases (Œ±=1.15 n=3.6) M and œÉ fix to 0-80% + 1 sigma
    {2318.81, 922.092, 559.358, 388.048, 316.572, 121.677},
    // CB tail fixed to pp ‚Äúold‚Äù MC in all cases (Œ±=1.15 n=3.6) M and œÉ fix to 0-80% - 1 sigma
    {2048.48, 816.592, 509.062, 320.778, 282.295, 113.539},
    
    // ------ CB2 ------
    /*
     // CB2 tails fixed to pure J/Œ® simulation (pass2, Philippe) M and œÉ fixed to 0-80% fit
     {2234.45, 890.814, 548.486, 361.476, 305.59, 119.888},
     // CB2 tails fixed to pure J/Œ® simulation (pass2, Philippe) M and œÉ fixed to 0-80% fit + 1 sigma
     {2378.96, 946.687, 575.999, 396.942, 324.773, 123.878},
     // CB2 tails fixed to pure J/Œ® simulation (pass2, Philippe) M and œÉ fixed to 0-80% fit - 1 sigma
     {2079.44, 830.903, 517.584, 325.431, 284.551, 115.001},
     */
    // CB2 tails fixed to embedding (pass1, Javier) in all cases M and œÉ fixed to 0-80% fit
    {2256.64, 917.674, 563.277, 365.872, 309.336, 119.935},
    // CB2 tails fixed to embedding (pass1, Javier) in all cases M and œÉ fixed to 0-80% fit + 1 sigma
    {2414.61, 980.297, 593.982, 405.8, 329.895, 124.208},
    // CB2 tails fixed to embedding (pass1, Javier) in all cases M and œÉ fixed to 0-80% fit - 1 sigma
    {2084.15, 849.629, 528.584, 325.278, 286.167, 114.463},
    
    // ------ mixing + CB1 ------
    
    // Background from mixing, CB (Œ±=0.98  n=5.2) + 1 expo. + range [2;4] M and œÉ fixed to 0-80% fit
    {2390.45, 1004.35, 577.393, 434.207, 305.023, 114.492},
    // Background from mixing, CB (Œ±=0.98  n=5.2) + 1 expo. + range [2;4] M and œÉ fixed to 0-80% fit + 1 sigma
    {2515.28, 1060.17, 600.255, 461.875, 319.111, 115.488},
    // Background from mixing, CB (Œ±=0.98  n=5.2) + 1 expo. + range [2;4] M and œÉ fixed to 0-80% fit - 1 sigma
    {2251.76, 944.053, 550.814, 403.699, 288.8, 112.326},
    
    // Background from mixing, CB (Œ±=1.15  n=1.59) + 1 expo. + range [2;4.4] M and œÉ fixed to 0-80% fit
    {2617.16, 1098.85, 631.025, 469.823, 335.205, 122.868},
    // Background from mixing, CB (Œ±=1.15  n=1.59) + 1 expo. + range [2;4.4] M and œÉ fixed to 0-80% fit + 1 sigma
    {2724.08, 1141.5, 649.254, 493.803, 346.905, 123.353},
    // Background from mixing, CB (Œ±=1.15  n=1.59) + 1 expo. + range [2;4.4] M and œÉ fixed to 0-80% fit - 1 sigma
    {2479.27, 1040.09, 605.974, 439.355, 319.317, 121.16},
    
    // Background from mixing, CB (Œ±=0.98  n=5.2) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit
    {2317.06, 981.918, 552.516, 396.647, 276.668, 113.474},
    // Background from mixing, CB (Œ±=0.98  n=5.2) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit + 1 sigma
    {2451.1, 1041.49, 575.784, 429.452, 288.718, 115.2},
    // Background from mixing, CB (Œ±=0.98  n=5.2) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit - 1 sigma
    {2170.11, 917.927, 525.474, 362.086, 262.256, 110.66},
    
    // Background from mixing, CB (Œ±=1.15  n=1.59) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit
    {2567.51, 1089.19, 612.267, 437.791, 307.2, 122.519},
    // Background from mixing, CB (Œ±=1.15  n=1.59) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit + 1 sigma
    {2697.85, 1147.66, 633.521, 468.561, 318.118, 123.672},
    // Background from mixing, CB (Œ±=1.15  n=1.59) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit - 1 sigma
    {2420.72, 1024.95, 586.377, 402.142, 293.268, 120.195},
    
  };
  
  Double_t en[nT][nC] = {
    
    // ------ CB1 ------
    
    // CB tail fixed to simulations (pure J/Œ®) in all cases (Œ±=0.98 n=5.2) M and œÉ fixed to 0-80% fit
    {134.44, 99.7141, 67.4565, 46.4516, 35.6592, 15.1885},
    // CB tail fixed to simulations (pure J/Œ®) in all cases (Œ±=0.98 n=5.2) M and œÉ fixed to 0-80% fit: œÉ + 1 sigma
    {143.043, 106.206, 71.7208, 52.6941, 37.8042, 15.8515},
    // CB tail fixed to simulations (pure J/Œ®) in all cases (Œ±=0.98 n=5.2) M and œÉ fixed to 0-80% fit: œÉ - 1 sigma
    {125.742, 92.9699, 62.4126, 46.1791, 33.4652, 14.4957},
    /*
     // CB tail fixed to simulations (pure J/Œ®) in all cases (Œ±=0.98 n=5.2) M free and œÉ fixed (74 instead of 76) to 0-80% fit
     {132.712, 98.6414, 66.888, 48.9598, 35.3893, 15.361},
     // CB tail fixed to simulations (pure J/Œ®) in all cases (Œ±=0.98 n=5.2) M and œÉ fixed to 0-80% fit: M + 1 sigma
     {130.162, 97.3748, 64.8484, 46.6657, 34.783, 14.673},
     // CB tail fixed to simulations (pure J/Œ®) in all cases (Œ±=0.98 n=5.2) M and œÉ fixed to 0-80% fit: M - 1 sigma
     {133.27, 98.7987, 66.8767, 46.7818, 35.4391, 15.1923},
     */
    // CB tail fixed to pp LHC10g in all cases (Œ±=1.15 n=1.59) M and œÉ fix to 0-80%
    {150.18, 102.969, 74.9456, 53.0865, 39.6457, 17.1552},
    // CB tail fixed to pp LHC10g in all cases (Œ±=1.15 n=1.59) M and œÉ fix to 0-80% + 1 sigma
    {156.832, 117.185, 79.1499, 56.7817, 41.7559, 17.534},
    // CB tail fixed to pp LHC10g in all cases (Œ±=1.15 n=1.59) M and œÉ fix to 0-80% - 1 sigma
    {141.516, 97.9642, 70.617, 50.2681, 37.4772, 16.1694},
    
    // CB tail fixed to pp ‚Äúold‚Äù MC in all cases (Œ±=1.15 n=3.6) M and œÉ fix to 0-80%
    {134.833, 100.037, 67.6587, 48.1231, 35.7797, 15.2072},
    // CB tail fixed to pp ‚Äúold‚Äù MC in all cases (Œ±=1.15 n=3.6) M and œÉ fix to 0-80% + 1 sigma
    {142.993, 105.125, 71.703, 50.3564, 37.7807, 15.82},
    // CB tail fixed to pp ‚Äúold‚Äù MC in all cases (Œ±=1.15 n=3.6) M and œÉ fix to 0-80% - 1 sigma
    {126.486, 93.7361, 63.5367, 45.6999, 33.6804, 14.5534},
    
    // ------ CB2 ------
    /*
     // CB2 tails fixed to pure J/Œ® simulation (pass2, Philippe) M and œÉ fixed to 0-80% fit
     {129.946, 101.876, 68.281, 49.4867, 36.4153, 15.542},
     // CB2 tails fixed to pure J/Œ® simulation (pass2, Philippe) M and œÉ fixed to 0-80% fit + 1 sigma
     {146.379, 107.281, 73.4044, 51.4688, 38.6727, 16.2362},
     // CB2 tails fixed to pure J/Œ® simulation (pass2, Philippe) M and œÉ fixed to 0-80% fit - 1 sigma
     {128.197, 94.9917, 64.4061, 47.6612, 34.1405, 14.8109},
     */
    // CB2 tails fixed to embedding (pass1, Javier) in all cases M and œÉ fixed to 0-80% fit
    {138.767, 100.812, 71.1305, 50.3989, 36.7773, 15.5368},
    // CB2 tails fixed to embedding (pass1, Javier) in all cases M and œÉ fixed to 0-80% fit + 1 sigma
    {147.639, 106.537, 76.282, 55.3495, 39.2381, 16.2746},
    // CB2 tails fixed to embedding (pass1, Javier) in all cases M and œÉ fixed to 0-80% fit - 1 sigma
    {128.627, 97.2074, 65.9661, 44.6864, 34.2315, 14.2289},
    
    // ------ mixing + CB1 ------
    
    // Background from mixing, CB (Œ±=0.98  n=5.2) + 1 expo. + range [2;4] M and œÉ fixed to 0-80% fit
    {132.406, 99.3915, 66.5593, 41.7434, 36.7759, 14.1805},
    // Background from mixing, CB (Œ±=0.98  n=5.2) + 1 expo. + range [2;4] M and œÉ fixed to 0-80% fit + 1 sigma
    {139.876, 104.762, 70.3044, 43.7173, 38.2638, 14.5188},
    // Background from mixing, CB (Œ±=0.98  n=5.2) + 1 expo. + range [2;4] M and œÉ fixed to 0-80% fit - 1 sigma
    {124.85, 93.0023, 62.7949, 39.6663, 34.3464, 13.8469},
    
    // Background from mixing, CB (Œ±=1.15  n=1.59) + 1 expo. + range [2;4.4] M and œÉ fixed to 0-80% fit
    {145.803, 102.715, 72.6466, 45.4445, 41.092, 15.2794},
    // Background from mixing, CB (Œ±=1.15  n=1.59) + 1 expo. + range [2;4.4] M and œÉ fixed to 0-80% fit + 1 sigma
    {152.66, 106.749, 75.7527, 47.1222, 43.8861, 15.539},
    // Background from mixing, CB (Œ±=1.15  n=1.59) + 1 expo. + range [2;4.4] M and œÉ fixed to 0-80% fit - 1 sigma
    {137.994, 103.03, 69.1256, 43.3965, 38.2456, 14.9228},
    
    // Background from mixing, CB (Œ±=0.98  n=5.2) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit
    {135.404, 99.9435, 67.9387, 46.361, 36.05, 14.4447},
    // Background from mixing, CB (Œ±=0.98  n=5.2) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit + 1 sigma
    {143.574, 106.309, 71.9753, 49.4461, 38.0861, 14.8958},
    // Background from mixing, CB (Œ±=0.98  n=5.2) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit - 1 sigma
    {124.967, 93.8037, 63.6633, 43.4203, 33.9627, 13.9508},
    
    // Background from mixing, CB (Œ±=1.15  n=1.59) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit
    {150.289, 111.255, 75.3876, 51.8039, 40.0182, 15.7053},
    // Background from mixing, CB (Œ±=1.15  n=1.59) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit + 1 sigma
    {158.097, 117.811, 79.3478, 48.5908, 41.9958, 16.0971},
    // Background from mixing, CB (Œ±=1.15  n=1.59) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit - 1 sigma
    {141.255, 104.512, 70.9752, 48.564, 37.9425, 15.2593}
    
  };
  
  // return the values
  nTests = nT;
  nCentBins = nC;
  nJPsi = new Double_t*[nT];
  enJPsi = new Double_t*[nT];
  for (Int_t i=0; i<nT; i++) {
    nJPsi[i] = new Double_t[nC];
    enJPsi[i] = new Double_t[nC];
    for (Int_t j=0; j<nC; j++) {
      nJPsi[i][j] = n[i][j];
      enJPsi[i][j] = en[i][j];
    }
  }
  
}


//------------------------------------------------------------------------
void SetNJPsiVsCenty1(Double_t **&nJPsi, Double_t **&enJPsi, Int_t &nTests, Int_t &nCentBins)
{
  // return the number of JPsi in -4<y<-3.25 and errors for each test and each centrality bins
  // as well as the number of tests and of centrality bins
  
  const Int_t nT = 24;
  const Int_t nC = 6;
  
  Double_t n[nT][nC] = {
    
    // ------ CB1 ------
    
    // CB tail fixed to simulations (pure J/Œ®) in all cases (Œ±=0.98 n=5.2) M and œÉ fixed to 0-80% fit
    {851.313, 347.622, 241.577, 128.61, 130.668, 46.2461},
    // CB tail fixed to simulations (pure J/Œ®) in all cases (Œ±=0.98 n=5.2) M and œÉ fixed to 0-80% fit: œÉ + 1 sigma
    {927.498, 374.328, 264.225, 144.545, 141.788, 50.0376},
    // CB tail fixed to simulations (pure J/Œ®) in all cases (Œ±=0.98 n=5.2) M and œÉ fixed to 0-80% fit: œÉ - 1 sigma
    {764.184, 316.959, 216.192, 110.752, 117.005, 41.5492},
    
    // CB tail fixed to pp LHC10g in all cases (Œ±=1.15 n=1.59) M and œÉ fix to 0-80%
    {931.55, 380.479, 264.67, 139.946, 143.558, 50.7102},
    // CB tail fixed to pp LHC10g in all cases (Œ±=1.15 n=1.59) M and œÉ fix to 0-80% + 1 sigma
    {1011.76, 408.094, 288.549, 157.196, 155.439, 54.7571},
    // CB tail fixed to pp LHC10g in all cases (Œ±=1.15 n=1.59) M and œÉ fix to 0-80% - 1 sigma
    {835.504, 346.966, 236.817, 119.901, 128.328, 45.3675},
    
    // CB tail fixed to pp ‚Äúold‚Äù MC in all cases (Œ±=1.15 n=3.6) M and œÉ fix to 0-80%
    {851.087, 345.748, 242.246, 128.852, 131.455, 46.5601},
    // CB tail fixed to pp ‚Äúold‚Äù MC in all cases (Œ±=1.15 n=3.6) M and œÉ fix to 0-80% + 1 sigma
    {922.776, 370.809, 263.561, 143.996, 141.843, 50.0143},
    // CB tail fixed to pp ‚Äúold‚Äù MC in all cases (Œ±=1.15 n=3.6) M and œÉ fix to 0-80% - 1 sigma
    {767.759, 316.33, 218.076, 111.629, 118.514, 42.0832},
    
    // ------ CB2 ------
    
    // CB2 tails fixed to embedding (pass1, Javier) in all cases M and œÉ fixed to 0-80% fit
    {872.744, 353.963, 270.366, 134.703, 133.648, 46.8638},
    // CB2 tails fixed to embedding (pass1, Javier) in all cases M and œÉ fixed to 0-80% fit + 1 sigma
    {951.228, 380.019, 296.861, 151.628, 145.308, 49.9457},
    // CB2 tails fixed to embedding (pass1, Javier) in all cases M and œÉ fixed to 0-80% fit - 1 sigma
    {783.65, 322.995, 241.347, 116.493, 119.639, 42.546},
    
    // ------ mixing + CB1 ------
    
    // Background from mixing, CB (Œ±=0.98  n=5.2) + 1 expo. + range [2;4] M and œÉ fixed to 0-80% fit
    {926.365, 392.909, 258.001, 151.332, 112.933, 40.0402},
    // Background from mixing, CB (Œ±=0.98  n=5.2) + 1 expo. + range [2;4] M and œÉ fixed to 0-80% fit + 1 sigma
    {987.889, 413.545, 276.861, 165.083, 116.22, 48.2809},
    // Background from mixing, CB (Œ±=0.98  n=5.2) + 1 expo. + range [2;4] M and œÉ fixed to 0-80% fit - 1 sigma
    {844.358, 341.319, 236.532, 134.932, 106.901, 37.6027},
    
    // Background from mixing, CB (Œ±=1.15  n=1.59) + 1 expo. + range [2;4.4] M and œÉ fixed to 0-80% fit
    {1004.35, 423.823, 281.136, 162.331, 121.547, 48.7732},
    // Background from mixing, CB (Œ±=1.15  n=1.59) + 1 expo. + range [2;4.4] M and œÉ fixed to 0-80% fit + 1 sigma
    {1055.74, 441.03, 297.709, 174.537, 124.176, 50.8707},
    // Background from mixing, CB (Œ±=1.15  n=1.59) + 1 expo. + range [2;4.4] M and œÉ fixed to 0-80% fit - 1 sigma
    {938.828, 401.972, 262.045, 147.719, 116.533, 45.9685},
    
    // Background from mixing, CB (Œ±=0.98  n=5.2) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit
    {926.659, 368.414, 244.107, 140.613, 116.83, 39.3194},
    // Background from mixing, CB (Œ±=0.98  n=5.2) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit + 1 sigma
    {999.169, 388.283, 262.366, 154.069, 123.187, 40.6816},
    // Background from mixing, CB (Œ±=0.98  n=5.2) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit - 1 sigma
    {843.147, 343.63, 223.317, 124.57, 107.952, 37.0964},
    
    // Background from mixing, CB (Œ±=1.15  n=1.59) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit
    {1021.68, 399.274, 269.091, 151.798, 128.592, 40.8967},
    // Background from mixing, CB (Œ±=1.15  n=1.59) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit + 1 sigma
    {1067.64, 414.214, 283.934, 162.58, 133.843, 41.8852},
    // Background from mixing, CB (Œ±=1.15  n=1.59) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit - 1 sigma
    {953.588, 380.736, 252.272, 139.132, 121.454, 39.3473}
    
  };
  
  Double_t en[nT][nC] = {
    
    // ------ CB1 ------
    
    // CB tail fixed to simulations (pure J/Œ®) in all cases (Œ±=0.98 n=5.2) M and œÉ fixed to 0-80% fit
    {84.5089, 61.7555, 42.8026, 24.6292, 19.6246, 9.9129},
    // CB tail fixed to simulations (pure J/Œ®) in all cases (Œ±=0.98 n=5.2) M and œÉ fixed to 0-80% fit: œÉ + 1 sigma
    {83.5203, 60.9336, 43.1178, 26.7767, 20.9081, 10.7243},
    // CB tail fixed to simulations (pure J/Œ®) in all cases (Œ±=0.98 n=5.2) M and œÉ fixed to 0-80% fit: œÉ - 1 sigma
    {78.3402, 51.0368, 38.1606, 22.3876, 17.731, 9.17286},
    
    // CB tail fixed to pp LHC10g in all cases (Œ±=1.15 n=1.59) M and œÉ fix to 0-80%
    {91.7787, 67.0182, 46.323, 27.0049, 21.5739, 10.8857},
    // CB tail fixed to pp LHC10g in all cases (Œ±=1.15 n=1.59) M and œÉ fix to 0-80% + 1 sigma
    {93.147, 72.8853, 48.0653, 29.2693, 23.2403, 11.6885},
    // CB tail fixed to pp LHC10g in all cases (Œ±=1.15 n=1.59) M and œÉ fix to 0-80% - 1 sigma
    {84.6779, 61.4215, 39.5494, 24.5536, 19.7213, 10.0468},
    
    // CB tail fixed to pp ‚Äúold‚Äù MC in all cases (Œ±=1.15 n=3.6) M and œÉ fix to 0-80%
    {84.9281, 62.2076, 42.9312, 24.6431, 19.1524, 9.88199},
    // CB tail fixed to pp ‚Äúold‚Äù MC in all cases (Œ±=1.15 n=3.6) M and œÉ fix to 0-80% + 1 sigma
    {83.7909, 57.2514, 43.8155, 26.674, 21.1447, 10.6982},
    // CB tail fixed to pp ‚Äúold‚Äù MC in all cases (Œ±=1.15 n=3.6) M and œÉ fix to 0-80% - 1 sigma
    {78.9805, 57.5877, 40.6451, 22.4922, 18.0364, 9.23347},
    
    // ------ CB2 ------
    
    // CB2 tails fixed to embedding (pass1, Javier) in all cases M and œÉ fixed to 0-80% fit
    {88.3829, 65.5713, 44.6215, 25.4999, 19.8701, 9.66005},
    // CB2 tails fixed to embedding (pass1, Javier) in all cases M and œÉ fixed to 0-80% fit + 1 sigma
    {95.216, 71.4089, 47.0777, 27.7461, 21.4406, 10.3331},
    // CB2 tails fixed to embedding (pass1, Javier) in all cases M and œÉ fixed to 0-80% fit - 1 sigma
    {82.777, 60.5635, 42.7591, 23.1687, 18.1802, 9.3845},
    
    // ------ mixing + CB1 ------
    
    // Background from mixing, CB (Œ±=0.98  n=5.2) + 1 expo. + range [2;4] M and œÉ fixed to 0-80% fit
    {79.5522, 49.1628, 38.1189, 23.2394, 18.2457, 9.67064},
    // Background from mixing, CB (Œ±=0.98  n=5.2) + 1 expo. + range [2;4] M and œÉ fixed to 0-80% fit + 1 sigma
    {76.4307, 52.1056, 41.2529, 24.7964, 19.187, 9.94317},
    // Background from mixing, CB (Œ±=0.98  n=5.2) + 1 expo. + range [2;4] M and œÉ fixed to 0-80% fit - 1 sigma
    {70.2505, 48.4915, 34.9666, 21.5677, 17.1367, 8.9897},
    
    // Background from mixing, CB (Œ±=1.15  n=1.59) + 1 expo. + range [2;4.4] M and œÉ fixed to 0-80% fit
    {77.5898, 53.0944, 41.4721, 25.1044, 19.6977, 10.1954},
    // Background from mixing, CB (Œ±=1.15  n=1.59) + 1 expo. + range [2;4.4] M and œÉ fixed to 0-80% fit + 1 sigma
    {81.8426, 55.6943, 44.2593, 26.4888, 20.4872, 10.6849},
    // Background from mixing, CB (Œ±=1.15  n=1.59) + 1 expo. + range [2;4.4] M and œÉ fixed to 0-80% fit - 1 sigma
    {79.6728, 50.3125, 38.5846, 23.6248, 18.7583, 9.62442},
    
    // Background from mixing, CB (Œ±=0.98  n=5.2) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit
    {77.9401, 51.5338, 39.5595, 24.143, 20.6023, 9.54801},
    // Background from mixing, CB (Œ±=0.98  n=5.2) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit + 1 sigma
    {108.096, 54.9658, 42.7954, 25.8447, 22.2185, 10.0978},
    // Background from mixing, CB (Œ±=0.98  n=5.2) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit - 1 sigma
    {71.2215, 47.8578, 36.2199, 25.8097, 18.9063, 8.90175},
    
    // Background from mixing, CB (Œ±=1.15  n=1.59) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit
    {129.364, 55.986, 43.6404, 26.2299, 22.7736, 10.3478},
    // Background from mixing, CB (Œ±=1.15  n=1.59) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit + 1 sigma
    {80.3408, 58.6775, 46.283, 27.573, 23.8571, 10.7716},
    // Background from mixing, CB (Œ±=1.15  n=1.59) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit - 1 sigma
    {80.5299, 53.114, 40.9029, 24.7801, 21.3977, 9.8575}
    
  };
  
  // return the values
  nTests = nT;
  nCentBins = nC;
  nJPsi = new Double_t*[nT];
  enJPsi = new Double_t*[nT];
  for (Int_t i=0; i<nT; i++) {
    nJPsi[i] = new Double_t[nC];
    enJPsi[i] = new Double_t[nC];
    for (Int_t j=0; j<nC; j++) {
      nJPsi[i][j] = n[i][j];
      enJPsi[i][j] = en[i][j];
    }
  }
  
}


//------------------------------------------------------------------------
void SetNJPsiVsCenty2(Double_t **&nJPsi, Double_t **&enJPsi, Int_t &nTests, Int_t &nCentBins)
{
  // return the number of JPsi in -3.25<y<-2 and errors for each test and each centrality bins
  // as well as the number of tests and of centrality bins
  
  const Int_t nT = 24;
  const Int_t nC = 6;
  
  Double_t n[nT][nC] = {
    
    // ------ CB1 ------
    
    // CB tail fixed to simulations (pure J/Œ®) in all cases (Œ±=0.98 n=5.2) M and œÉ fixed to 0-80% fit
    {1258.48, 496.111, 294.859, 219.912, 170.073, 77.6667},
    // CB tail fixed to simulations (pure J/Œ®) in all cases (Œ±=0.98 n=5.2) M and œÉ fixed to 0-80% fit: œÉ + 1 sigma
    {1353.69, 534.215, 308.953, 244.134, 182.173, 81.9705},
    // CB tail fixed to simulations (pure J/Œ®) in all cases (Œ±=0.98 n=5.2) M and œÉ fixed to 0-80% fit: œÉ - 1 sigma
    {1153.64, 454.328, 277.823, 194.953, 156.434, 71.9037},
    
    // CB tail fixed to pp LHC10g in all cases (Œ±=1.15 n=1.59) M and œÉ fix to 0-80%
    {1404.6, 552.772, 328.234, 245.091, 191.441, 86.6497},
    // CB tail fixed to pp LHC10g in all cases (Œ±=1.15 n=1.59) M and œÉ fix to 0-80% + 1 sigma
    {1509.21, 594.938, 342.247, 272.938, 204.325, 91.2391},
    // CB tail fixed to pp LHC10g in all cases (Œ±=1.15 n=1.59) M and œÉ fix to 0-80% - 1 sigma
    {1285.43, 505.42, 309.763, 216.152, 175.705, 80.0336},
    
    // CB tail fixed to pp ‚Äúold‚Äù MC in all cases (Œ±=1.15 n=3.6) M and œÉ fix to 0-80%
    {1252.55, 490.995, 290.368, 222.527, 172.808, 77.2557},
    // CB tail fixed to pp ‚Äúold‚Äù MC in all cases (Œ±=1.15 n=3.6) M and œÉ fix to 0-80% + 1 sigma
    {1349.02, 530.432, 302.821, 248.187, 184.921, 81.1428},
    // CB tail fixed to pp ‚Äúold‚Äù MC in all cases (Œ±=1.15 n=3.6) M and œÉ fix to 0-80% - 1 sigma
    {1143.87, 447.419, 274.201, 195.825, 158.238, 71.3205},
    
    // ------ CB2 ------
    
    // CB2 tails fixed to embedding (pass1, Javier) in all cases M and œÉ fixed to 0-80% fit
    {1303.15, 521.919, 308.063, 224.449, 176.463, 79.4745},
    // CB2 tails fixed to embedding (pass1, Javier) in all cases M and œÉ fixed to 0-80% fit + 1 sigma
    {1410.88, 565.406, 324.29, 251.931, 190.108, 84.2883},
    // CB2 tails fixed to embedding (pass1, Javier) in all cases M and œÉ fixed to 0-80% fit - 1 sigma
    {1184.03, 474.003, 288.206, 196.444, 160.982, 73.0032},
    
    // ------ mixing + CB1 ------
    
    // Background from mixing, CB (Œ±=0.98  n=5.2) + 1 expo. + range [2;4] M and œÉ fixed to 0-80% fit
    {1537.08, 653.002, 347.278, 286.188, 189.417, 77.3785},
    // Background from mixing, CB (Œ±=0.98  n=5.2) + 1 expo. + range [2;4] M and œÉ fixed to 0-80% fit + 1 sigma
    {1636.52, 699.02, 362.154, 308.065, 199.865, 77.9798},
    // Background from mixing, CB (Œ±=0.98  n=5.2) + 1 expo. + range [2;4] M and œÉ fixed to 0-80% fit - 1 sigma
    {1423.82, 602.273, 328.434, 261.604, 177.121, 74.7773},
    
    // Background from mixing, CB (Œ±=1.15  n=1.59) + 1 expo. + range [2;4.4] M and œÉ fixed to 0-80% fit
    {1670.63, 708.782, 377.211, 306.891, 207.639, 83.2952},
    // Background from mixing, CB (Œ±=1.15  n=1.59) + 1 expo. + range [2;4.4] M and œÉ fixed to 0-80% fit + 1 sigma
    {1767.2, 753.68, 390.318, 328.542, 217.002, 85.9694},
    // Background from mixing, CB (Œ±=1.15  n=1.59) + 1 expo. + range [2;4.4] M and œÉ fixed to 0-80% fit - 1 sigma
    {1557.89, 658.157, 359.697, 282.316, 196.006, 80.6893},
    
    // Background from mixing, CB (Œ±=0.98  n=5.2) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit
    {1434.17, 603.773, 328.633, 258.801, 170.183, 76.2157},
    // Background from mixing, CB (Œ±=0.98  n=5.2) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit + 1 sigma
    {1541.58, 652.564, 344.54, 286.376, 180.121, 77.5662},
    // Background from mixing, CB (Œ±=0.98  n=5.2) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit - 1 sigma
    {1314.23, 550.81, 308.709, 229.947, 158.275, 72.7966},
    
    // Background from mixing, CB (Œ±=1.15  n=1.59) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit
    {1581.83, 665.837, 363.264, 283.689, 189.273, 82.3346},
    // Background from mixing, CB (Œ±=1.15  n=1.59) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit + 1 sigma
    {1697.52, 719, 378.993, 314.976, 198.845, 83.3879},
    // Background from mixing, CB (Œ±=1.15  n=1.59) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit - 1 sigma
    {1446.61, 606.124, 341.565, 250.363, 176.148, 78.5929}
    
  };
  
  Double_t en[nT][nC] = {
    
    // ------ CB1 ------
    
    // CB tail fixed to simulations (pure J/Œ®) in all cases (Œ±=0.98 n=5.2) M and œÉ fixed to 0-80% fit
    {99.0675, 78.324, 52.6912, 36.2296, 28.051, 12.0125},
    // CB tail fixed to simulations (pure J/Œ®) in all cases (Œ±=0.98 n=5.2) M and œÉ fixed to 0-80% fit: œÉ + 1 sigma
    {113.817, 84.5837, 56.7495, 42.8258, 30.1112, 12.5386},
    // CB tail fixed to simulations (pure J/Œ®) in all cases (Œ±=0.98 n=5.2) M and œÉ fixed to 0-80% fit: œÉ - 1 sigma
    {97.0984, 71.0488, 48.5993, 39.712, 25.9429, 11.2484},
    
    // CB tail fixed to pp LHC10g in all cases (Œ±=1.15 n=1.59) M and œÉ fix to 0-80%
    {112.347, 85.9821, 59.2827, 44.2922, 31.3817, 13.4013},
    // CB tail fixed to pp LHC10g in all cases (Œ±=1.15 n=1.59) M and œÉ fix to 0-80% + 1 sigma
    {127.418, 94.832, 63.5985, 46.6767, 33.6571, 14.1857},
    // CB tail fixed to pp LHC10g in all cases (Œ±=1.15 n=1.59) M and œÉ fix to 0-80% - 1 sigma
    {108.586, 80.6752, 54.4401, 38.0442, 28.9829, 12.2751},
    
    {105.396, 78.1783, 52.5944, 42.6021, 28.0472, 11.7979},
    // CB tail fixed to pp ‚Äúold‚Äù MC in all cases (Œ±=1.15 n=3.6) M and œÉ fix to 0-80%
    {114.129, 84.7699, 56.8648, 43.1828, 30.222, 12.7287},
    // CB tail fixed to pp ‚Äúold‚Äù MC in all cases (Œ±=1.15 n=3.6) M and œÉ fix to 0-80% + 1 sigma
    {96.4429, 71.4576, 48.2414, 41.8117, 25.7891, 11.1673},
    // CB tail fixed to pp ‚Äúold‚Äù MC in all cases (Œ±=1.15 n=3.6) M and œÉ fix to 0-80% - 1 sigma
    
    // ------ CB2 ------
    
    // CB2 tails fixed to embedding (pass1, Javier) in all cases M and œÉ fixed to 0-80% fit
    {109.383, 82.3468, 55.3138, 37.3597, 29.0779, 12.3144},
    // CB2 tails fixed to embedding (pass1, Javier) in all cases M and œÉ fixed to 0-80% fit + 1 sigma
    {118.808, 89.5739, 60.0435, 44.0272, 31.4207, 13.1572},
    // CB2 tails fixed to embedding (pass1, Javier) in all cases M and œÉ fixed to 0-80% fit - 1 sigma
    {99.7437, 74.9715, 50.5393, 34.6954, 26.6614, 11.4893},
    
    // ------ mixing + CB1 ------
    
    // Background from mixing, CB (Œ±=0.98  n=5.2) + 1 expo. + range [2;4] M and œÉ fixed to 0-80% fit
    {111.601, 82.9199, 56.1825, 35.0825, 30.1519, 12.4015},
    // Background from mixing, CB (Œ±=0.98  n=5.2) + 1 expo. + range [2;4] M and œÉ fixed to 0-80% fit + 1 sigma
    {119.256, 88.4596, 60.3849, 37.1722, 32.3893, 12.9342},
    // Background from mixing, CB (Œ±=0.98  n=5.2) + 1 expo. + range [2;4] M and œÉ fixed to 0-80% fit - 1 sigma
    {103.647, 76.8446, 51.992, 32.9402, 27.9632, 11.7794},
    
    // Background from mixing, CB (Œ±=1.15  n=1.59) + 1 expo. + range [2;4.4] M and œÉ fixed to 0-80% fit
    {121.555, 90.2505, 60.8451, 37.8623, 32.8333, 13.3782},
    // Background from mixing, CB (Œ±=1.15  n=1.59) + 1 expo. + range [2;4.4] M and œÉ fixed to 0-80% fit + 1 sigma
    {129.172, 96.0708, 64.8843, 40.0378, 35.0714, 14.2941},
    // Background from mixing, CB (Œ±=1.15  n=1.59) + 1 expo. + range [2;4.4] M and œÉ fixed to 0-80% fit - 1 sigma
    {113.586, 84.181, 56.7582, 35.7372, 30.6454, 12.7858},
    
    // Background from mixing, CB (Œ±=0.98  n=5.2) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit
    {111.651, 82.8247, 55.5843, 38.1371, 29.639, 11.9804},
    // Background from mixing, CB (Œ±=0.98  n=5.2) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit + 1 sigma
    {120.408, 89.4553, 59.8861, 41.2773, 31.7938, 12.524},
    // Background from mixing, CB (Œ±=0.98  n=5.2) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit - 1 sigma
    {102.677, 75.9927, 51.1032, 35.0399, 27.4146, 11.3451},
    
    // Background from mixing, CB (Œ±=1.15  n=1.59) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit
    {123.404, 91.5624, 61.4248, 42.2417, 32.6355, 13.0099},
    // Background from mixing, CB (Œ±=1.15  n=1.59) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit + 1 sigma
    {132.921, 98.7914, 66.9546, 46.2925, 35.0799, 13.5465},
    // Background from mixing, CB (Œ±=1.15  n=1.59) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit - 1 sigma
    {113.359, 83.9441, 56.5089, 38.6951, 30.3053, 12.3377}
    
  };
  
  // return the values
  nTests = nT;
  nCentBins = nC;
  nJPsi = new Double_t*[nT];
  enJPsi = new Double_t*[nT];
  for (Int_t i=0; i<nT; i++) {
    nJPsi[i] = new Double_t[nC];
    enJPsi[i] = new Double_t[nC];
    for (Int_t j=0; j<nC; j++) {
      nJPsi[i][j] = n[i][j];
      enJPsi[i][j] = en[i][j];
    }
  }
  
}


//------------------------------------------------------------------------
void SetNJPsiVsCentpt1(Double_t **&nJPsi, Double_t **&enJPsi, Int_t &nTests, Int_t &nCentBins)
{
  // return the number of JPsi with pT<3GeV/c and errors for each test and each centrality bins
  // as well as the number of tests and of centrality bins
  
  const Int_t nT = 24;
  const Int_t nC = 6;
  
  Double_t n[nT][nC] = {
    
    // ------ CB1 ------
    
    // CB tail fixed to simulations (pure J/Œ®) in all cases (Œ±=0.98 n=5.2) M and œÉ fixed to 0-80% fit
    {1617.33, 638.692, 412.956, 286.446, 204.783, 80.2432},
    // CB tail fixed to simulations (pure J/Œ®) in all cases (Œ±=0.98 n=5.2) M and œÉ fixed to 0-80% fit: œÉ + 1 sigma
    {1719.76, 679.383, 431.855, 313.008, 217.208, 83.0669},
    // CB tail fixed to simulations (pure J/Œ®) in all cases (Œ±=0.98 n=5.2) M and œÉ fixed to 0-80% fit: œÉ - 1 sigma
    {1504.54, 594.62, 390.949, 258.444, 190.035, 76.626},
    
    // CB tail fixed to pp LHC10g in all cases (Œ±=1.15 n=1.59) M and œÉ fix to 0-80%
    {1768.42, 695.614, 453.286, 311.778, 225.899, 88.0493},
    // CB tail fixed to pp LHC10g in all cases (Œ±=1.15 n=1.59) M and œÉ fix to 0-80% + 1 sigma
    {1883.72, 741.295, 473.578, 342.59, 240.134, 91.1166},
    // CB tail fixed to pp LHC10g in all cases (Œ±=1.15 n=1.59) M and œÉ fix to 0-80% - 1 sigma
    {1637.71, 645.134, 428.373, 278.965, 208.125, 83.8404},
    
    // CB tail fixed to pp ‚Äúold‚Äù MC in all cases (Œ±=1.15 n=3.6) M and œÉ fix to 0-80%
    {1605.8, 632.717, 408.574, 285.267, 205.023, 79.8484},
    // CB tail fixed to pp ‚Äúold‚Äù MC in all cases (Œ±=1.15 n=3.6) M and œÉ fix to 0-80% + 1 sigma
    {1715.23, 676.842, 428.189, 314.206, 217.68, 82.7439},
    // CB tail fixed to pp ‚Äúold‚Äù MC in all cases (Œ±=1.15 n=3.6) M and œÉ fix to 0-80% - 1 sigma
    {1482.75, 584.244, 384.905, 254.335, 189.236, 75.932},
    
    // ------ CB2 ------
    
    // CB2 tails fixed to embedding (pass1, Javier) in all cases M and œÉ fixed to 0-80% fit
    {1639.59, 639.851, 426.99, 298.479, 208.957, 80.9967},
    // CB2 tails fixed to embedding (pass1, Javier) in all cases M and œÉ fixed to 0-80% fit + 1 sigma
    {1756.89, 686.018, 449.234, 329.691, 223.6, 84.128},
    // CB2 tails fixed to embedding (pass1, Javier) in all cases M and œÉ fixed to 0-80% fit - 1 sigma
    {1509.62, 589.88, 400.982, 266.035, 191.534, 76.8668},
    
    // ------ mixing + CB1 ------
    
    // Background from mixing, CB (Œ±=0.98  n=5.2) + 1 expo. + range [2;4] M and œÉ fixed to 0-80% fit
    {1843.49, 870.53, 464.779, 343.105, 219.528, 78.6781},
    // Background from mixing, CB (Œ±=0.98  n=5.2) + 1 expo. + range [2;4] M and œÉ fixed to 0-80% fit + 1 sigma
    {1952.99, 920.264, 485.626, 350.719, 230.557, 79.768},
    // Background from mixing, CB (Œ±=0.98  n=5.2) + 1 expo. + range [2;4] M and œÉ fixed to 0-80% fit - 1 sigma
    {1720.15, 814.715, 440.19, 316.411, 205.793, 76.7342},
    
    // Background from mixing, CB (Œ±=1.15  n=1.59) + 1 expo. + range [2;4.4] M and œÉ fixed to 0-80% fit
    {2008.43, 940.306, 507.335, 368.258, 240.982, 84.7256},
    // Background from mixing, CB (Œ±=1.15  n=1.59) + 1 expo. + range [2;4.4] M and œÉ fixed to 0-80% fit + 1 sigma
    {2116.9, 987.867, 526.715, 391.673, 251.401, 88.1777},
    // Background from mixing, CB (Œ±=1.15  n=1.59) + 1 expo. + range [2;4.4] M and œÉ fixed to 0-80% fit - 1 sigma
    {1883.09, 885.481, 483.314, 341.481, 226.986, 82.9759},
    
    // Background from mixing, CB (Œ±=0.98  n=5.2) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit
    {1843.49, 870.53, 464.779, 343.105, 219.528, 78.6781},
    // Background from mixing, CB (Œ±=0.98  n=5.2) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit + 1 sigma
    {1952.99, 920.264, 485.626, 350.719, 230.557, 79.768},
    // Background from mixing, CB (Œ±=0.98  n=5.2) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit - 1 sigma
    {1720.15, 814.715, 440.19, 316.411, 205.793, 76.7342},
    
    // Background from mixing, CB (Œ±=1.15  n=1.59) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit
    {1950.55, 826.36, 480.842, 347.269, 213.431, 84.0398},
    // Background from mixing, CB (Œ±=1.15  n=1.59) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit + 1 sigma
    {2069.17, 880.687, 500.275, 371.679, 222.31, 84.9235},
    // Background from mixing, CB (Œ±=1.15  n=1.59) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit - 1 sigma
    {1815.94, 767.048, 456.776, 316.383, 200.958, 82.152}
    
  };
  
  Double_t en[nT][nC] = {
    
    // ------ CB1 ------
    
    // CB tail fixed to simulations (pure J/Œ®) in all cases (Œ±=0.98 n=5.2) M and œÉ fixed to 0-80% fit
    {116.244, 87.292, 59.8398, 38.6591, 30.4151, 12.6413},
    // CB tail fixed to simulations (pure J/Œ®) in all cases (Œ±=0.98 n=5.2) M and œÉ fixed to 0-80% fit: œÉ + 1 sigma
    {124.354, 93.0231, 63.6734, 41.0511, 32.3534, 13.2018},
    // CB tail fixed to simulations (pure J/Œ®) in all cases (Œ±=0.98 n=5.2) M and œÉ fixed to 0-80% fit: œÉ - 1 sigma
    {108.881, 83.2699, 56.0753, 38.5379, 28.4204, 12.0446},
    
    // CB tail fixed to pp LHC10g in all cases (Œ±=1.15 n=1.59) M and œÉ fix to 0-80%
    {120.505, 95.8974, 65.8385, 44.7432, 33.3637, 13.8961},
    // CB tail fixed to pp LHC10g in all cases (Œ±=1.15 n=1.59) M and œÉ fix to 0-80% + 1 sigma
    {140.528, 103.072, 70.1389, 44.9658, 35.9914, 14.5089},
    // CB tail fixed to pp LHC10g in all cases (Œ±=1.15 n=1.59) M and œÉ fix to 0-80% - 1 sigma
    {119.108, 92.4968, 61.8263, 39.6102, 31.0586, 13.2178},
    
    // CB tail fixed to pp ‚Äúold‚Äù MC in all cases (Œ±=1.15 n=3.6) M and œÉ fix to 0-80%
    {114.234, 86.4449, 59.6023, 38.4807, 30.2915, 12.5794},
    // CB tail fixed to pp ‚Äúold‚Äù MC in all cases (Œ±=1.15 n=3.6) M and œÉ fix to 0-80% + 1 sigma
    {124.361, 92.9996, 63.7255, 41.0809, 32.3709, 13.1726},
    // CB tail fixed to pp ‚Äúold‚Äù MC in all cases (Œ±=1.15 n=3.6) M and œÉ fix to 0-80% - 1 sigma
    {107.574, 82.739, 53.6922, 35.8878, 27.8855, 11.9136},
    
    // ------ CB2 ------
    
    // CB2 tails fixed to embedding (pass1, Javier) in all cases M and œÉ fixed to 0-80% fit
    {115.237, 87.5927, 61.8683, 40.0884, 30.8036, 12.7971},
    // CB2 tails fixed to embedding (pass1, Javier) in all cases M and œÉ fixed to 0-80% fit + 1 sigma
    {127.167, 94.5516, 66.4088, 42.9504, 33.1984, 13.4337},
    // CB2 tails fixed to embedding (pass1, Javier) in all cases M and œÉ fixed to 0-80% fit - 1 sigma
    {108.575, 80.3193, 57.3539, 39.7059, 28.4492, 11.6139},
    
    // ------ mixing + CB1 ------
    
    // Background from mixing, CB (Œ±=0.98  n=5.2) + 1 expo. + range [2;4] M and œÉ fixed to 0-80% fit
    {119.184, 80.027, 60.3255, 37.6009, 31.98, 12.5351},
    // Background from mixing, CB (Œ±=0.98  n=5.2) + 1 expo. + range [2;4] M and œÉ fixed to 0-80% fit + 1 sigma
    {127.051, 84.1004, 64.3381, 41.5284, 34.1253, 12.9775},
    // Background from mixing, CB (Œ±=0.98  n=5.2) + 1 expo. + range [2;4] M and œÉ fixed to 0-80% fit - 1 sigma
    {111.303, 75.6655, 56.4416, 35.5257, 29.7591, 12.0444},
    
    // Background from mixing, CB (Œ±=1.15  n=1.59) + 1 expo. + range [2;4.4] M and œÉ fixed to 0-80% fit
    {130.544, 86.4225, 65.8304, 40.6187, 35.1767, 13.5675},
    // Background from mixing, CB (Œ±=1.15  n=1.59) + 1 expo. + range [2;4.4] M and œÉ fixed to 0-80% fit + 1 sigma
    {138.814, 90.2767, 69.684, 42.5355, 37.4253, 14.0703},
    // Background from mixing, CB (Œ±=1.15  n=1.59) + 1 expo. + range [2;4.4] M and œÉ fixed to 0-80% fit - 1 sigma
    {122.324, 82.1738, 61.8083, 38.5713, 32.8016, 13.1059},
    
    // Background from mixing, CB (Œ±=0.98  n=5.2) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit
    {119.184, 80.027, 60.3255, 37.6009, 31.98, 12.5351},
    // Background from mixing, CB (Œ±=0.98  n=5.2) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit + 1 sigma
    {127.051, 84.1004, 64.3381, 41.5284, 34.1253, 12.9775},
    // Background from mixing, CB (Œ±=0.98  n=5.2) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit - 1 sigma
    {111.303, 75.6655, 56.4416, 35.5257, 29.7591, 12.0444},
    
    // Background from mixing, CB (Œ±=1.15  n=1.59) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit
    {133.823, 99.3407, 67.2928, 41.489, 35.0212, 13.314},
    // Background from mixing, CB (Œ±=1.15  n=1.59) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit + 1 sigma
    {142.088, 105.381, 71.4805, 43.5984, 37.1461, 13.699},
    // Background from mixing, CB (Œ±=1.15  n=1.59) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit - 1 sigma
    {124.994, 92.5088, 62.9521, 43.1157, 32.7367, 12.8589}
    
  };
  
  // return the values
  nTests = nT;
  nCentBins = nC;
  nJPsi = new Double_t*[nT];
  enJPsi = new Double_t*[nT];
  for (Int_t i=0; i<nT; i++) {
    nJPsi[i] = new Double_t[nC];
    enJPsi[i] = new Double_t[nC];
    for (Int_t j=0; j<nC; j++) {
      nJPsi[i][j] = n[i][j];
      enJPsi[i][j] = en[i][j];
    }
  }
  
}


//------------------------------------------------------------------------
void SetNJPsiVsCentpt2(Double_t **&nJPsi, Double_t **&enJPsi, Int_t &nTests, Int_t &nCentBins)
{
  // return the number of JPsi with pT>3GeV/c and errors for each test and each centrality bins
  // as well as the number of tests and of centrality bins
  
  const Int_t nT = 24;
  const Int_t nC = 4;
  
  Double_t n[nT][nC] = {
    
    // ------ CB1 ------
    
    // CB tail fixed to simulations (pure J/Œ®) in all cases (Œ±=0.98 n=5.2) M and œÉ fixed to 0-80% fit
    {514.196, 317.6, 141.249, 59.534},
    // CB tail fixed to simulations (pure J/Œ®) in all cases (Œ±=0.98 n=5.2) M and œÉ fixed to 0-80% fit: œÉ + 1 sigma
    {564.279, 342.963, 160.992, 64.6686},
    // CB tail fixed to simulations (pure J/Œ®) in all cases (Œ±=0.98 n=5.2) M and œÉ fixed to 0-80% fit: œÉ - 1 sigma
    {456.934, 287.6, 119.574, 54.0183},
    
    // CB tail fixed to pp LHC10g in all cases (Œ±=1.15 n=1.59) M and œÉ fix to 0-80%
    {570.799, 353.847, 155.75, 65.539},
    // CB tail fixed to pp LHC10g in all cases (Œ±=1.15 n=1.59) M and œÉ fix to 0-80% + 1 sigma
    {620.573, 378.491, 175.786, 70.8642},
    // CB tail fixed to pp LHC10g in all cases (Œ±=1.15 n=1.59) M and œÉ fix to 0-80% - 1 sigma
    {514.017, 325.026, 133.644, 59.7585},
    
    // CB tail fixed to pp ‚Äúold‚Äù MC in all cases (Œ±=1.15 n=3.6) M and œÉ fix to 0-80%
    {506.696, 312.532, 139.13, 59.2589},
    // CB tail fixed to pp ‚Äúold‚Äù MC in all cases (Œ±=1.15 n=3.6) M and œÉ fix to 0-80% + 1 sigma
    {555.131, 336.727, 158.522, 64.3398},
    // CB tail fixed to pp ‚Äúold‚Äù MC in all cases (Œ±=1.15 n=3.6) M and œÉ fix to 0-80% - 1 sigma
    {452.028, 284.627, 117.917, 53.7017},
    
    // ------ CB2 ------
    
    // CB2 tails fixed to embedding (pass1, Javier) in all cases M and œÉ fixed to 0-80% fit
    {517.476, 329.749, 138.575, 60.764},
    // CB2 tails fixed to embedding (pass1, Javier) in all cases M and œÉ fixed to 0-80% fit + 1 sigma
    {565.948, 354.847, 157.719, 65.5767},
    // CB2 tails fixed to embedding (pass1, Javier) in all cases M and œÉ fixed to 0-80% fit - 1 sigma
    {462.582, 299.761, 118.001, 55.4394},
    
    // ------ mixing + CB1 ------
    
    // Background from mixing, CB (Œ±=0.98  n=5.2) + 1 expo. + range [2;4] M and œÉ fixed to 0-80% fit
    {560.348, 342.714, 159.693, 47.9475},
    // Background from mixing, CB (Œ±=0.98  n=5.2) + 1 expo. + range [2;4] M and œÉ fixed to 0-80% fit + 1 sigma
    {593.961, 360.279, 175.801, 48.838},
    // Background from mixing, CB (Œ±=0.98  n=5.2) + 1 expo. + range [2;4] M and œÉ fixed to 0-80% fit - 1 sigma
    {504.732, 318.327, 141.246, 46.1684},
    
    // Background from mixing, CB (Œ±=1.15  n=1.59) + 1 expo. + range [2;4.4] M and œÉ fixed to 0-80% fit
    {594.647, 365.941, 167.914, 50.1909},
    // Background from mixing, CB (Œ±=1.15  n=1.59) + 1 expo. + range [2;4.4] M and œÉ fixed to 0-80% fit + 1 sigma
    {636.613, 395.716, 184.621, 51.2535},
    // Background from mixing, CB (Œ±=1.15  n=1.59) + 1 expo. + range [2;4.4] M and œÉ fixed to 0-80% fit - 1 sigma
    {552.409, 344.534, 148.699, 48.2471},
    
    // Background from mixing, CB (Œ±=0.98  n=5.2) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit
    {547.054, 340.555, 152.431, 52.0147},
    // Background from mixing, CB (Œ±=0.98  n=5.2) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit + 1 sigma
    {591.127, 365.816, 169.048, 54.7776},
    // Background from mixing, CB (Œ±=0.98  n=5.2) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit - 1 sigma
    {490.97, 309.979, 132.447, 48.5049},
    
    // Background from mixing, CB (Œ±=1.15  n=1.59) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit
    {593.737, 373.982, 162.962, 56.5464},
    // Background from mixing, CB (Œ±=1.15  n=1.59) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit + 1 sigma
    {635.059, 396.231, 179.194, 58.9787},
    // Background from mixing, CB (Œ±=1.15  n=1.59) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit - 1 sigma
    {542.818, 344.383, 143.37, 52.8693}
    
  };
  
  Double_t en[nT][nC] = {
    
    // ------ CB1 ------
    
    // CB tail fixed to simulations (pure J/Œ®) in all cases (Œ±=0.98 n=5.2) M and œÉ fixed to 0-80% fit
    {53.6924, 50.7627, 22.2675, 10.8222},
    // CB tail fixed to simulations (pure J/Œ®) in all cases (Œ±=0.98 n=5.2) M and œÉ fixed to 0-80% fit: œÉ + 1 sigma
    {58.2002, 51.4722, 26.226, 11.5237},
    // CB tail fixed to simulations (pure J/Œ®) in all cases (Œ±=0.98 n=5.2) M and œÉ fixed to 0-80% fit: œÉ - 1 sigma
    {52.0659, 43.1339, 20.1795, 10.0688},
    
    // CB tail fixed to pp LHC10g in all cases (Œ±=1.15 n=1.59) M and œÉ fix to 0-80%
    {59.1012, 52.2545, 24.4976, 12.7822},
    // CB tail fixed to pp LHC10g in all cases (Œ±=1.15 n=1.59) M and œÉ fix to 0-80% + 1 sigma
    {63.5336, 56.1919, 26.3775, 12.5669},
    // CB tail fixed to pp LHC10g in all cases (Œ±=1.15 n=1.59) M and œÉ fix to 0-80% - 1 sigma
    {54.3755, 48.0677, 23.8271, 11.1606},
    
    // CB tail fixed to pp ‚Äúold‚Äù MC in all cases (Œ±=1.15 n=3.6) M and œÉ fix to 0-80%
    {52.9998, 50.1906, 23.3295, 10.7193},
    // CB tail fixed to pp ‚Äúold‚Äù MC in all cases (Œ±=1.15 n=3.6) M and œÉ fix to 0-80% + 1 sigma
    {57.387, 50.7406, 23.8553, 11.4059},
    // CB tail fixed to pp ‚Äúold‚Äù MC in all cases (Œ±=1.15 n=3.6) M and œÉ fix to 0-80% - 1 sigma
    {48.3548, 42.7348, 20.7078, 9.99563},
    
    // ------ CB2 ------
    
    // CB2 tails fixed to embedding (pass1, Javier) in all cases M and œÉ fixed to 0-80% fit
    {53.841, 48.9327, 21.9912, 10.9334},
    // CB2 tails fixed to embedding (pass1, Javier) in all cases M and œÉ fixed to 0-80% fit + 1 sigma
    {58.1757, 52.8906, 23.8051, 11.6236},
    // CB2 tails fixed to embedding (pass1, Javier) in all cases M and œÉ fixed to 0-80% fit - 1 sigma
    {49.238, 44.6889, 20.0628, 10.2016},
    
    // ------ mixing + CB1 ------
    
    // Background from mixing, CB (Œ±=0.98  n=5.2) + 1 expo. + range [2;4] M and œÉ fixed to 0-80% fit
    {51.8373, 46.59, 22.8902, 10.5106},
    // Background from mixing, CB (Œ±=0.98  n=5.2) + 1 expo. + range [2;4] M and œÉ fixed to 0-80% fit + 1 sigma
    {55.144, 48.6379, 24.9313, 11.0098},
    // Background from mixing, CB (Œ±=0.98  n=5.2) + 1 expo. + range [2;4] M and œÉ fixed to 0-80% fit - 1 sigma
    {51.6881, 46.1581, 20.7118, 9.9474},
    
    // Background from mixing, CB (Œ±=1.15  n=1.59) + 1 expo. + range [2;4.4] M and œÉ fixed to 0-80% fit
    {55.031, 48.4778, 24.5431, 11.1799},
    // Background from mixing, CB (Œ±=1.15  n=1.59) + 1 expo. + range [2;4.4] M and œÉ fixed to 0-80% fit + 1 sigma
    {64.007, 57.808, 26.5387, 11.6582},
    // Background from mixing, CB (Œ±=1.15  n=1.59) + 1 expo. + range [2;4.4] M and œÉ fixed to 0-80% fit - 1 sigma
    {51.4835, 45.3433, 22.3669, 10.6455},
    
    // Background from mixing, CB (Œ±=0.98  n=5.2) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit
    {56.4573, 49.7608, 22.6984, 11.138},
    // Background from mixing, CB (Œ±=0.98  n=5.2) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit + 1 sigma
    {59.0077, 54.4117, 24.6202, 12.082},
    // Background from mixing, CB (Œ±=0.98  n=5.2) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit - 1 sigma
    {50.7573, 44.9122, 20.5817, 10.3572},
    
    // Background from mixing, CB (Œ±=1.15  n=1.59) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit
    {59.1045, 52.233, 24.6151, 11.6916},
    // Background from mixing, CB (Œ±=1.15  n=1.59) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit + 1 sigma
    {63.4905, 55.8096, 26.4819, 12.3281},
    // Background from mixing, CB (Œ±=1.15  n=1.59) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit - 1 sigma
    {54.3796, 49.8625, 22.5764, 11.5072}
    
  };
  
  // return the values
  nTests = nT;
  nCentBins = nC;
  nJPsi = new Double_t*[nT];
  enJPsi = new Double_t*[nT];
  for (Int_t i=0; i<nT; i++) {
    nJPsi[i] = new Double_t[nC];
    enJPsi[i] = new Double_t[nC];
    for (Int_t j=0; j<nC; j++) {
      nJPsi[i][j] = n[i][j];
      enJPsi[i][j] = en[i][j];
    }
  }
  
}


//------------------------------------------------------------------------
void SetNJPsiVsCentpt2y1(Double_t **&nJPsi, Double_t **&enJPsi, Int_t &nTests, Int_t &nCentBins)
{
  // return the number of JPsi with pT>3GeV/c and errors for each test and each centrality bins
  // as well as the number of tests and of centrality bins
  // return the number of JPsi in -4<y<-3.25 and errors for each test and each centrality bins
  // as well as the number of tests and of centrality bins
  
  const Int_t nT = 24;
  const Int_t nC = 1;
  
  Double_t n[nT][nC] = {
    
    // ------ CB1 ------
    
    // CB tail fixed to simulations (pure J/Œ®) in all cases (Œ±=0.98 n=5.2) M and œÉ fixed to 0-80% fit
    {120},
    // CB tail fixed to simulations (pure J/Œ®) in all cases (Œ±=0.98 n=5.2) M and œÉ fixed to 0-80% fit: œÉ + 1 sigma
    {139},
    // CB tail fixed to simulations (pure J/Œ®) in all cases (Œ±=0.98 n=5.2) M and œÉ fixed to 0-80% fit: œÉ - 1 sigma
    {94},
    
    // CB tail fixed to pp LHC10g in all cases (Œ±=1.15 n=1.59) M and œÉ fix to 0-80%
    {132},
    // CB tail fixed to pp LHC10g in all cases (Œ±=1.15 n=1.59) M and œÉ fix to 0-80% + 1 sigma
    {150},
    // CB tail fixed to pp LHC10g in all cases (Œ±=1.15 n=1.59) M and œÉ fix to 0-80% - 1 sigma
    {106},
    
    // CB tail fixed to pp ‚Äúold‚Äù MC in all cases (Œ±=1.15 n=3.6) M and œÉ fix to 0-80%
    {119},
    // CB tail fixed to pp ‚Äúold‚Äù MC in all cases (Œ±=1.15 n=3.6) M and œÉ fix to 0-80% + 1 sigma
    {137},
    // CB tail fixed to pp ‚Äúold‚Äù MC in all cases (Œ±=1.15 n=3.6) M and œÉ fix to 0-80% - 1 sigma
    {94},
    
    // ------ CB2 ------
    
    // CB2 tails fixed to embedding (pass1, Javier) in all cases M and œÉ fixed to 0-80% fit
    {115},
    // CB2 tails fixed to embedding (pass1, Javier) in all cases M and œÉ fixed to 0-80% fit + 1 sigma
    {135},
    // CB2 tails fixed to embedding (pass1, Javier) in all cases M and œÉ fixed to 0-80% fit - 1 sigma
    {88},
    
    // ------ mixing + CB1 ------
    
    // Background from mixing, CB (Œ±=0.98  n=5.2) + 1 expo. + range [2;4] M and œÉ fixed to 0-80% fit
    {164},
    // Background from mixing, CB (Œ±=0.98  n=5.2) + 1 expo. + range [2;4] M and œÉ fixed to 0-80% fit + 1 sigma
    {191},
    // Background from mixing, CB (Œ±=0.98  n=5.2) + 1 expo. + range [2;4] M and œÉ fixed to 0-80% fit - 1 sigma
    {138},
    
    // Background from mixing, CB (Œ±=1.15  n=1.59) + 1 expo. + range [2;4.4] M and œÉ fixed to 0-80% fit
    {175},
    // Background from mixing, CB (Œ±=1.15  n=1.59) + 1 expo. + range [2;4.4] M and œÉ fixed to 0-80% fit + 1 sigma
    {200},
    // Background from mixing, CB (Œ±=1.15  n=1.59) + 1 expo. + range [2;4.4] M and œÉ fixed to 0-80% fit - 1 sigma
    {150},
    
    // Background from mixing, CB (Œ±=0.98  n=5.2) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit
    {164},
    // Background from mixing, CB (Œ±=0.98  n=5.2) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit + 1 sigma
    {185},
    // Background from mixing, CB (Œ±=0.98  n=5.2) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit - 1 sigma
    {138},
    
    // Background from mixing, CB (Œ±=1.15  n=1.59) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit
    {177},
    // Background from mixing, CB (Œ±=1.15  n=1.59) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit + 1 sigma
    {194},
    // Background from mixing, CB (Œ±=1.15  n=1.59) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit - 1 sigma
    {152}
    
  };
  
  Double_t en[nT][nC] = {
    
    // ------ CB1 ------
    
    // CB tail fixed to simulations (pure J/Œ®) in all cases (Œ±=0.98 n=5.2) M and œÉ fixed to 0-80% fit
    {28},
    // CB tail fixed to simulations (pure J/Œ®) in all cases (Œ±=0.98 n=5.2) M and œÉ fixed to 0-80% fit: œÉ + 1 sigma
    {32},
    // CB tail fixed to simulations (pure J/Œ®) in all cases (Œ±=0.98 n=5.2) M and œÉ fixed to 0-80% fit: œÉ - 1 sigma
    {24},
    
    // CB tail fixed to pp LHC10g in all cases (Œ±=1.15 n=1.59) M and œÉ fix to 0-80%
    {31},
    // CB tail fixed to pp LHC10g in all cases (Œ±=1.15 n=1.59) M and œÉ fix to 0-80% + 1 sigma
    {35},
    // CB tail fixed to pp LHC10g in all cases (Œ±=1.15 n=1.59) M and œÉ fix to 0-80% - 1 sigma
    {27},
    
    // CB tail fixed to pp ‚Äúold‚Äù MC in all cases (Œ±=1.15 n=3.6) M and œÉ fix to 0-80%
    {28},
    // CB tail fixed to pp ‚Äúold‚Äù MC in all cases (Œ±=1.15 n=3.6) M and œÉ fix to 0-80% + 1 sigma
    {31},
    // CB tail fixed to pp ‚Äúold‚Äù MC in all cases (Œ±=1.15 n=3.6) M and œÉ fix to 0-80% - 1 sigma
    {24},
    
    // ------ CB2 ------
    
    // CB2 tails fixed to embedding (pass1, Javier) in all cases M and œÉ fixed to 0-80% fit
    {27},
    // CB2 tails fixed to embedding (pass1, Javier) in all cases M and œÉ fixed to 0-80% fit + 1 sigma
    {35},
    // CB2 tails fixed to embedding (pass1, Javier) in all cases M and œÉ fixed to 0-80% fit - 1 sigma
    {22},
    
    // ------ mixing + CB1 ------
    
    // Background from mixing, CB (Œ±=0.98  n=5.2) + 1 expo. + range [2;4] M and œÉ fixed to 0-80% fit
    {32},
    // Background from mixing, CB (Œ±=0.98  n=5.2) + 1 expo. + range [2;4] M and œÉ fixed to 0-80% fit + 1 sigma
    {34},
    // Background from mixing, CB (Œ±=0.98  n=5.2) + 1 expo. + range [2;4] M and œÉ fixed to 0-80% fit - 1 sigma
    {27},
    
    // Background from mixing, CB (Œ±=1.15  n=1.59) + 1 expo. + range [2;4.4] M and œÉ fixed to 0-80% fit
    {33},
    // Background from mixing, CB (Œ±=1.15  n=1.59) + 1 expo. + range [2;4.4] M and œÉ fixed to 0-80% fit + 1 sigma
    {34},
    // Background from mixing, CB (Œ±=1.15  n=1.59) + 1 expo. + range [2;4.4] M and œÉ fixed to 0-80% fit - 1 sigma
    {30},
    
    // Background from mixing, CB (Œ±=0.98  n=5.2) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit
    {32},
    // Background from mixing, CB (Œ±=0.98  n=5.2) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit + 1 sigma
    {35},
    // Background from mixing, CB (Œ±=0.98  n=5.2) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit - 1 sigma
    {27},
    
    // Background from mixing, CB (Œ±=1.15  n=1.59) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit
    {33},
    // Background from mixing, CB (Œ±=1.15  n=1.59) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit + 1 sigma
    {36},
    // Background from mixing, CB (Œ±=1.15  n=1.59) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit - 1 sigma
    {30}
    
  };
  
  // return the values
  nTests = nT;
  nCentBins = nC;
  nJPsi = new Double_t*[nT];
  enJPsi = new Double_t*[nT];
  for (Int_t i=0; i<nT; i++) {
    nJPsi[i] = new Double_t[nC];
    enJPsi[i] = new Double_t[nC];
    for (Int_t j=0; j<nC; j++) {
      nJPsi[i][j] = n[i][j];
      enJPsi[i][j] = en[i][j];
    }
  }
  
}

//------------------------------------------------------------------------
void SetNJPsiVsCentpt2y2(Double_t **&nJPsi, Double_t **&enJPsi, Int_t &nTests, Int_t &nCentBins)
{
  // return the number of JPsi with pT>3GeV/c and errors for each test and each centrality bins
  // as well as the number of tests and of centrality bins
  // return the number of JPsi in -3.25<y<-2.5 and errors for each test and each centrality bins
  // as well as the number of tests and of centrality bins
  
  const Int_t nT = 24;
  const Int_t nC = 1;
  
  Double_t n[nT][nC] = {
    
    // ------ CB1 ------
    
    // CB tail fixed to simulations (pure J/Œ®) in all cases (Œ±=0.98 n=5.2) M and œÉ fixed to 0-80% fit
    {369},
    // CB tail fixed to simulations (pure J/Œ®) in all cases (Œ±=0.98 n=5.2) M and œÉ fixed to 0-80% fit: œÉ + 1 sigma
    {409},
    // CB tail fixed to simulations (pure J/Œ®) in all cases (Œ±=0.98 n=5.2) M and œÉ fixed to 0-80% fit: œÉ - 1 sigma
    {320},
    
    // CB tail fixed to pp LHC10g in all cases (Œ±=1.15 n=1.59) M and œÉ fix to 0-80%
    {412},
    // CB tail fixed to pp LHC10g in all cases (Œ±=1.15 n=1.59) M and œÉ fix to 0-80% + 1 sigma
    {452},
    // CB tail fixed to pp LHC10g in all cases (Œ±=1.15 n=1.59) M and œÉ fix to 0-80% - 1 sigma
    {363},
    
    // CB tail fixed to pp ‚Äúold‚Äù MC in all cases (Œ±=1.15 n=3.6) M and œÉ fix to 0-80%
    {367},
    // CB tail fixed to pp ‚Äúold‚Äù MC in all cases (Œ±=1.15 n=3.6) M and œÉ fix to 0-80% + 1 sigma
    {406},
    // CB tail fixed to pp ‚Äúold‚Äù MC in all cases (Œ±=1.15 n=3.6) M and œÉ fix to 0-80% - 1 sigma
    {321},
    
    // ------ CB2 ------
    
    // CB2 tails fixed to embedding (pass1, Javier) in all cases M and œÉ fixed to 0-80% fit
    {375},
    // CB2 tails fixed to embedding (pass1, Javier) in all cases M and œÉ fixed to 0-80% fit + 1 sigma
    {417},
    // CB2 tails fixed to embedding (pass1, Javier) in all cases M and œÉ fixed to 0-80% fit - 1 sigma
    {326},
    
    // ------ mixing + CB1 ------
    
    // Background from mixing, CB (Œ±=0.98  n=5.2) + 1 expo. + range [2;4] M and œÉ fixed to 0-80% fit
    {385},
    // Background from mixing, CB (Œ±=0.98  n=5.2) + 1 expo. + range [2;4] M and œÉ fixed to 0-80% fit + 1 sigma
    {420},
    // Background from mixing, CB (Œ±=0.98  n=5.2) + 1 expo. + range [2;4] M and œÉ fixed to 0-80% fit - 1 sigma
    {344},
    
    // Background from mixing, CB (Œ±=1.15  n=1.59) + 1 expo. + range [2;4.4] M and œÉ fixed to 0-80% fit
    {425},
    // Background from mixing, CB (Œ±=1.15  n=1.59) + 1 expo. + range [2;4.4] M and œÉ fixed to 0-80% fit + 1 sigma
    {458},
    // Background from mixing, CB (Œ±=1.15  n=1.59) + 1 expo. + range [2;4.4] M and œÉ fixed to 0-80% fit - 1 sigma
    {383},
    
    // Background from mixing, CB (Œ±=0.98  n=5.2) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit
    {384},
    // Background from mixing, CB (Œ±=0.98  n=5.2) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit + 1 sigma
    {415},
    // Background from mixing, CB (Œ±=0.98  n=5.2) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit - 1 sigma
    {343},
    
    // Background from mixing, CB (Œ±=1.15  n=1.59) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit
    {418},
    // Background from mixing, CB (Œ±=1.15  n=1.59) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit + 1 sigma
    {447},
    // Background from mixing, CB (Œ±=1.15  n=1.59) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit - 1 sigma
    {380}
    
  };
  
  Double_t en[nT][nC] = {
    
    // ------ CB1 ------
    
    // CB tail fixed to simulations (pure J/Œ®) in all cases (Œ±=0.98 n=5.2) M and œÉ fixed to 0-80% fit
    {44},
    // CB tail fixed to simulations (pure J/Œ®) in all cases (Œ±=0.98 n=5.2) M and œÉ fixed to 0-80% fit: œÉ + 1 sigma
    {49},
    // CB tail fixed to simulations (pure J/Œ®) in all cases (Œ±=0.98 n=5.2) M and œÉ fixed to 0-80% fit: œÉ - 1 sigma
    {40},
    
    // CB tail fixed to pp LHC10g in all cases (Œ±=1.15 n=1.59) M and œÉ fix to 0-80%
    {49},
    // CB tail fixed to pp LHC10g in all cases (Œ±=1.15 n=1.59) M and œÉ fix to 0-80% + 1 sigma
    {53},
    // CB tail fixed to pp LHC10g in all cases (Œ±=1.15 n=1.59) M and œÉ fix to 0-80% - 1 sigma
    {44},
    
    // CB tail fixed to pp ‚Äúold‚Äù MC in all cases (Œ±=1.15 n=3.6) M and œÉ fix to 0-80%
    {44},
    // CB tail fixed to pp ‚Äúold‚Äù MC in all cases (Œ±=1.15 n=3.6) M and œÉ fix to 0-80% + 1 sigma
    {48},
    // CB tail fixed to pp ‚Äúold‚Äù MC in all cases (Œ±=1.15 n=3.6) M and œÉ fix to 0-80% - 1 sigma
    {40},
    
    // ------ CB2 ------
    
    // CB2 tails fixed to embedding (pass1, Javier) in all cases M and œÉ fixed to 0-80% fit
    {45},
    // CB2 tails fixed to embedding (pass1, Javier) in all cases M and œÉ fixed to 0-80% fit + 1 sigma
    {49},
    // CB2 tails fixed to embedding (pass1, Javier) in all cases M and œÉ fixed to 0-80% fit - 1 sigma
    {40},
    
    // ------ mixing + CB1 ------
    
    // Background from mixing, CB (Œ±=0.98  n=5.2) + 1 expo. + range [2;4] M and œÉ fixed to 0-80% fit
    {45},
    // Background from mixing, CB (Œ±=0.98  n=5.2) + 1 expo. + range [2;4] M and œÉ fixed to 0-80% fit + 1 sigma
    {50},
    // Background from mixing, CB (Œ±=0.98  n=5.2) + 1 expo. + range [2;4] M and œÉ fixed to 0-80% fit - 1 sigma
    {38},
    
    // Background from mixing, CB (Œ±=1.15  n=1.59) + 1 expo. + range [2;4.4] M and œÉ fixed to 0-80% fit
    {50},
    // Background from mixing, CB (Œ±=1.15  n=1.59) + 1 expo. + range [2;4.4] M and œÉ fixed to 0-80% fit + 1 sigma
    {54},
    // Background from mixing, CB (Œ±=1.15  n=1.59) + 1 expo. + range [2;4.4] M and œÉ fixed to 0-80% fit - 1 sigma
    {44},
    
    // Background from mixing, CB (Œ±=0.98  n=5.2) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit
    {44},
    // Background from mixing, CB (Œ±=0.98  n=5.2) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit + 1 sigma
    {48},
    // Background from mixing, CB (Œ±=0.98  n=5.2) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit - 1 sigma
    {41},
    
    // Background from mixing, CB (Œ±=1.15  n=1.59) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit
    {48},
    // Background from mixing, CB (Œ±=1.15  n=1.59) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit + 1 sigma
    {52},
    // Background from mixing, CB (Œ±=1.15  n=1.59) + pol1 + range [2;4.4] M and œÉ fixed to 0-80% fit - 1 sigma
    {44}
    
  };
  
  // return the values
  nTests = nT;
  nCentBins = nC;
  nJPsi = new Double_t*[nT];
  enJPsi = new Double_t*[nT];
  for (Int_t i=0; i<nT; i++) {
    nJPsi[i] = new Double_t[nC];
    enJPsi[i] = new Double_t[nC];
    for (Int_t j=0; j<nC; j++) {
      nJPsi[i][j] = n[i][j];
      enJPsi[i][j] = en[i][j];
    }
  }
  
}


//------------------------------------------------------------------------
void ComputeNJpsiMean(Int_t pt, Int_t y, Double_t *&nJPsiMean, Double_t *&enJPsiMean, Double_t *&sysnJPsiMean, Int_t &nCentBins, Bool_t weightedMeanRMS)
{
  // return the average number of JPsi over all the tests and the corresponding
  // statistical and systematic errors for the given rapidity or pt bin
  
  // set results for all tests
  Int_t nTests;
  Double_t **nJPsi;
  Double_t **enJPsi;
  if (y == 0 && pt == 0) {
    printf("\nresults integrated over pt and y:\n\n");
    SetNJPsiVsCent(nJPsi, enJPsi, nTests, nCentBins);
  } else if (y == 1 && pt == 0) {
    printf("\nresults integrated over pt in -4<y<-3.25:\n\n");
    SetNJPsiVsCenty1(nJPsi, enJPsi, nTests, nCentBins);
  } else if (y == 2 && pt == 0) {
    printf("\nresults integrated over pt in -3.25<y<-2.5:\n\n");
    SetNJPsiVsCenty2(nJPsi, enJPsi, nTests, nCentBins);
  } else if (y == 0 && pt == 1) {
    printf("\nresults integrated over y with 0<pt<3GeV/c:\n\n");
    SetNJPsiVsCentpt1(nJPsi, enJPsi, nTests, nCentBins);
  } else if (y == 0 && pt == 2) {
    printf("\nresults integrated over y with pt>3GeV/c:\n\n");
    SetNJPsiVsCentpt2(nJPsi, enJPsi, nTests, nCentBins);
  } else if (y == 1 && pt == 2) {
    printf("\nresults with -4<y<-3.25 and  with pt>3GeV/c:\n\n");
    SetNJPsiVsCentpt2y1(nJPsi, enJPsi, nTests, nCentBins);
  } else if (y == 2 && pt == 2) {
    printf("\nresults with -3.25<y<-2.5 and  with pt>3GeV/c:\n\n");
    SetNJPsiVsCentpt2y2(nJPsi, enJPsi, nTests, nCentBins);
  }
  
  // define histograms
  TGraphErrors **gNJpsi = new TGraphErrors*[nCentBins];
  TH1F **hNJpsi = new TH1F*[nCentBins];
  TH1F **heNJpsi = new TH1F*[nCentBins];
  TString *gName, *gTitle;
  SetGraphNameTitleVsCent(nCentBins, gName, gTitle);
  for(Int_t i=0; i<nCentBins; i++) {
    gNJpsi[i] = new TGraphErrors(nTests);
    gNJpsi[i]->SetName(Form("g%s_pt%d_y%d",gName[i].Data(),pt,y));
    gNJpsi[i]->SetTitle(Form("g%s_pt%d_y%d",gTitle[i].Data(),pt,y));
    hNJpsi[i] = new TH1F(Form("h%s_pt%d_y%d",gName[i].Data(),pt,y), Form("h%s_pt%d_y%d",gTitle[i].Data(),pt,y), 300000, 0., 3000.);
    heNJpsi[i] = new TH1F(Form("he%s_pt%d_y%d",gName[i].Data(),pt,y), Form("he%s_pt%d_y%d",gTitle[i].Data(),pt,y), 50000, 0., 500.);
  }
  
  // loop over tests
  Double_t *enJPsiMean2 = new Double_t[nCentBins];
  for(Int_t j=0; j<nCentBins; j++) enJPsiMean2[j] = 0.;
  Double_t *noSys = new Double_t[nCentBins];
  Double_t *renJPsi, *RAA, *eRAA, *eRAA2;
  Double_t eRAACorr;
  TGraphErrors *g, *gSys, *gSysAll;
  TCanvas *c = new TCanvas(Form("compareRAA_pt%d_y%d",pt,y), Form("compareRAA -- pt%d -- y%d",pt,y));
  gPad->SetFrameBorderMode(0);
  for(Int_t i=0; i<nTests; i++) {
    
    // display current RAA
    ComputeRAA(pt, y, nCentBins, nJPsi[i], enJPsi[i], noSys, renJPsi, RAA, eRAA, eRAA2, eRAACorr);
    if (nCentBins > 1) {
      BuildRAAGraphs(nCentBins, renJPsi, RAA, eRAA, eRAA2, eRAACorr, g, gSys, gSysAll, kTRUE);
      g->SetMarkerStyle(gMarkerStyleAliceForw);
      g->SetMarkerSize(gMarkerSizeAliceForw);
      g->SetMarkerColor(i+1);
      g->SetLineColor(i+1);
      c->cd();
      if (i==0) {
	g->GetYaxis()->SetRangeUser(0., 1.2);
	g->Draw("ap");
      }
      else g->Draw("p");
    }
    
    // loop over centrality classes
    for(Int_t j=0; j<nCentBins; j++) {
      
      Double_t wstat = 1./nJPsi[i][j];
      Double_t w = weightedMeanRMS ? 1./enJPsi[i][j]/enJPsi[i][j]/wstat : 1.;
//      Double_t w = 1./(enJPsi[i][j]*enJPsi[i][j]-nJPsi[i][j]);
      enJPsiMean2[j] += w*w*enJPsi[i][j]*enJPsi[i][j];
      
      // fill histogram
      gNJpsi[j]->SetPoint(i, i, nJPsi[i][j]);
      gNJpsi[j]->SetPointError(i, 0, enJPsi[i][j]);
      hNJpsi[j]->Fill(nJPsi[i][j], w);
      heNJpsi[j]->Fill(enJPsi[i][j], w);
      
    }
    
  }
  
  // compute NJpsi means and errors and display histograms
  nJPsiMean = new Double_t[nCentBins];
  enJPsiMean = new Double_t[nCentBins];
  sysnJPsiMean = new Double_t[nCentBins];
  TCanvas *sys1;
  if (fineDisplay) sys1 = new TCanvas(Form("sys1_pt%d_y%d",pt,y), Form("systematics of NJpsi -- pt%d -- y%d",pt,y), 1200, 180*nCentBins);
  else sys1 = new TCanvas(Form("sys1_pt%d_y%d",pt,y), Form("systematics of NJpsi -- pt%d -- y%d",pt,y), 400, 180*nCentBins);
  gPad->SetFrameBorderMode(0);
  if (fineDisplay) sys1->Divide(3,nCentBins);
  else sys1->Divide(1,nCentBins);
  for(Int_t i=0; i<nCentBins; i++) {
    
    // compute NJpsi means and errors
    nJPsiMean[i] = hNJpsi[i]->GetMean();
    sysnJPsiMean[i] = hNJpsi[i]->GetRMS();
    //enJPsiMean[i] = heNJpsi[i]->GetMean();
    enJPsiMean[i] = TMath::Sqrt(enJPsiMean2[i]*nTests)/hNJpsi[i]->GetSumOfWeights();
    
    // find largest difference from mean value
    Double_t maxDiff = -1.;
    for(Int_t j=0; j<nTests; j++) {
      if (TMath::Abs(nJPsi[j][i]-nJPsiMean[i]) > maxDiff) maxDiff = TMath::Abs(nJPsi[j][i]-nJPsiMean[i]);
    }
    
    // print results
    printf("NJpsi %s = %4.0f ± %3.0f (%4.1f%%) ± %3.0f (%4.1f%%) --- maxDiff = %3.0f (%3.1f RMS)\n",gName[i].Data(),
	   nJPsiMean[i], enJPsiMean[i], 100.*enJPsiMean[i]/nJPsiMean[i],
	   sysnJPsiMean[i], 100.*sysnJPsiMean[i]/nJPsiMean[i],
	   maxDiff, maxDiff / sysnJPsiMean[i]);
    
    // dispay histograms
    if (fineDisplay) sys1->cd(3*i+1);
    else sys1->cd(i+1);
    gPad->SetFrameBorderMode(0);
    gPad->SetTopMargin(0.01);
    gPad->SetBottomMargin(0.01);
    gNJpsi[i]->GetXaxis()->Set(nTests,-0.5,nTests-0.5);
    gNJpsi[i]->GetXaxis()->SetLabelColor(10);
    gNJpsi[i]->GetYaxis()->SetTitle("N_{J/#psi}");
    gNJpsi[i]->GetYaxis()->SetTitleSize(0.1);
    gNJpsi[i]->GetYaxis()->SetTitleOffset(0.4);
    gNJpsi[i]->GetYaxis()->SetLabelSize(0.08);
    gNJpsi[i]->SetMarkerStyle(20);
    gNJpsi[i]->Draw("ap");
    TPaveText *t = new TPaveText(0.2, 0.8, 0.5, 0.95,"NDC");
    t->AddText(0.,0.,gTitle[i].Data());
    t->SetFillStyle(0);
    t->SetBorderSize(0);
    t->SetTextFont(42);
    t->Draw();
    TLine *l0 = new TLine(-0.5, nJPsiMean[i], nTests-0.5, nJPsiMean[i]);
    l0->Draw();
    TLine *l1 = new TLine(-0.5, nJPsiMean[i]-sysnJPsiMean[i], nTests-0.5, nJPsiMean[i]-sysnJPsiMean[i]);
    l1->SetLineStyle(3);
    l1->Draw();
    TLine *l2 = new TLine(-0.5, nJPsiMean[i]+sysnJPsiMean[i], nTests-0.5, nJPsiMean[i]+sysnJPsiMean[i]);
    l2->SetLineStyle(3);
    l2->Draw();
    if (fineDisplay) {
      sys1->cd(3*i+2);
      gPad->SetFrameBorderMode(0);
      hNJpsi[i]->Draw();
      sys1->cd(3*i+3);
      gPad->SetFrameBorderMode(0);
      heNJpsi[i]->Draw();
    }
  }
  printf("\n");
  
}


//------------------------------------------------------------------------
void ComputeRAA(Int_t pt, Int_t y, Int_t nCentBins, Double_t *nJPsi, Double_t *enJPsi, Double_t *sysnJPsi, Double_t *&renJPsi,
		Double_t *&RAA, Double_t *&eRAA, Double_t *&eRAA2, Double_t &eRAACorr, Bool_t print)
{
  // compute the RAA versus centrality in the given pt or rapidity bin
  
  // get graph name in nCentBins centrality bins for printout
  TString *gName, *gTitle;
  SetGraphNameTitleVsCent(nCentBins, gName, gTitle);
  
  // get centrality independent normalization variables
  Double_t sigmaJPsipp, esigmaJPsipp;
  SetNormVar(pt, y, sigmaJPsipp, esigmaJPsipp);
  
  // get normalization variables in nCentBins centrality bins
  Double_t *nMB, *TAA, *eTAA, *AccEff, *eAccEff; 
  SetNormVarVsCent(pt, y, nCentBins, nMB, TAA, eTAA,AccEff, eAccEff);
  
  // get the systematic errors in nCentBins centrality bins
  Double_t *eGen, *eRec, *eTrk, *eTrg;
  SetSystVarVsCent(nCentBins, eGen, eRec, eTrk, eTrg);
  
  // Compute relative statistic & systematic errors
  renJPsi = new Double_t[nCentBins];
  Double_t *rsysnJPsi = new Double_t[nCentBins];
  Double_t *renMB = new Double_t[nCentBins];
  Double_t *reTAA = new Double_t[nCentBins];
  Double_t *reAccEff = new Double_t[nCentBins];
  for(Int_t i=0; i<nCentBins; i++) {    
    renJPsi[i]=enJPsi[i]/nJPsi[i];
    rsysnJPsi[i]=sysnJPsi[i]/nJPsi[i];
    renMB[i]=1./TMath::Sqrt(nMB[i]);
    reTAA[i]=eTAA[i]/TAA[i];
    reAccEff[i]=eAccEff[i]/AccEff[i];
  }
  
  // compute RAA and related errors for each centrality bin
  RAA = new Double_t[nCentBins];
  eRAA = new Double_t[nCentBins];
  eRAA2 = new Double_t[nCentBins];
  for(Int_t i=0; i<nCentBins; i++) { 
    
    RAA[i] = nJPsi[i]/BR/AccEff[i]/nMB[i]/sigmaJPsipp*1000./TAA[i]; // 1000. factor mub to mb 
    
    //printf(" nJPsi[i] %f BR %f AccEff[i] %f nMB[i] %f sigmaJPsipp %f TAA[i] %f \n",nJPsi[i],BR,AccEff[i],nMB[i],sigmaJPsipp,TAA[i]); 
    
    eRAA[i] = RAA[i] * TMath::Sqrt(rsysnJPsi[i]*rsysnJPsi[i] +
				   eTrk[i]*eTrk[i] +
				   eTrg[i]*eTrg[i] +
				   renMB[i]*renMB[i]);
    
    eRAA2[i] = RAA[i] * reTAA[i];
    
    Double_t eRAATot = TMath::Sqrt(eRAA[i]*eRAA[i] + eRAA2[i]*eRAA2[i]);
    
    if (print) printf("RAA %s = %5.3f ± %5.3f ± %5.3f\n", gName[i].Data(), RAA[i], renJPsi[i]*RAA[i], eRAATot);
    
  }
  if (print) printf("\n");
  
  //printf(" eGen %f eRec %f eTrk %f eTrg %f esigmaJPsipp %f eppLumi %f \n",eGen[nCentBins],eRec[nCentBins],eTrk[nCentBins],eTrg[nCentBins],esigmaJPsipp,eppLumi); 
  
  // correlated error**2 for RAA (without error on luminosity)
  Double_t eRAACorr2 = eGen[nCentBins]*eGen[nCentBins] +
		       eRec[nCentBins]*eRec[nCentBins] +
                       eTrk[nCentBins]*eTrk[nCentBins] +
		       eTrg[nCentBins]*eTrg[nCentBins] +
		       //eBR*eBR +
		       esigmaJPsipp*esigmaJPsipp;
  
  if (print) printf("correlated RAA error (wo lumi) = %4.1f%%\n", 100.*TMath::Sqrt(eRAACorr2));
  
  // add relative error on the luminosity
  eRAACorr2 += eppLumi*eppLumi;
  
  // correlated error for RAA
  eRAACorr = TMath::Sqrt(eRAACorr2);
  
  if (print) printf("correlated RAA error = %4.1f%%\n\n", 100.*eRAACorr);
  
}


//------------------------------------------------------------------------
void BuildRAAGraphs(Int_t nCentBins, Double_t *renJPsi, Double_t *RAA, Double_t *eRAA, Double_t *eRAA2, Double_t eRAACorr,
		    TGraphErrors *&gRAA, TGraphErrors *&gRAASys, TGraphErrors *&gRAASysAll, Int_t xAxisType, Int_t color, Int_t sysPos)
{
  // build the TGraphs to display our RAA results
  
  // Get graph x-axis labels according to the number of centrality bins
  Double_t *x, *ex, *bin, *NPart, *eNPart, *dNchdEta, *edNchdEta, *dNchdEtaNorm, *edNchdEtaNorm;;
  TString *label;
  SetGraphLabelsVsCent(nCentBins, x, ex, bin, label, NPart, eNPart,dNchdEta, edNchdEta, dNchdEtaNorm, edNchdEtaNorm);
  
  // fill RAA graphs
  gRAA = new TGraphErrors(nCentBins-1);
  gRAASys = new TGraphErrors(nCentBins-1);
  gRAASysAll = new TGraphErrors(1);
  if (xAxisType==0) {
    gRAASysAll->SetPoint(0,99.-3.*sysPos,1.);
    gRAASysAll->SetPointError(0,1.,eRAACorr);
  } else if (xAxisType==1){
    gRAASysAll->SetPoint(0,395.-12.*sysPos,1.);
    gRAASysAll->SetPointError(0,5.,eRAACorr);
  } else if (xAxisType==2){
    gRAASysAll->SetPoint(0,1505,1.);
    gRAASysAll->SetPointError(0,15,eRAACorr);
  } else if (xAxisType==3){
    gRAASysAll->SetPoint(0,9.4,1.);
    gRAASysAll->SetPointError(0,.1,eRAACorr);
  }
  for(Int_t i=nCentBins-1; i>=1; i--) {
    if (xAxisType==0) {
      gRAA->SetPoint(nCentBins-1-i,100.-x[i],RAA[i]);
      gRAA->SetPointError(nCentBins-1-i,ex[i],renJPsi[i]*RAA[i]);
      gRAASys->SetPoint(nCentBins-1-i,100.-x[i],RAA[i]);
      gRAASys->SetPointError(nCentBins-1-i,1.,TMath::Sqrt(eRAA[i]*eRAA[i] + eRAA2[i]*eRAA2[i]));
    } else if (xAxisType==1){
      gRAA->SetPoint(nCentBins-1-i,NPart[i],RAA[i]);
      gRAA->SetPointError(nCentBins-1-i,eNPart[i],renJPsi[i]*RAA[i]);
      gRAASys->SetPoint(nCentBins-1-i,NPart[i],RAA[i]);
      gRAASys->SetPointError(nCentBins-1-i,5.,TMath::Sqrt(eRAA[i]*eRAA[i] + eRAA2[i]*eRAA2[i]));
    } else if (xAxisType==2){
      gRAA->SetPoint(nCentBins-1-i,dNchdEta[i],RAA[i]);
      gRAA->SetPointError(nCentBins-1-i,edNchdEta[i],renJPsi[i]*RAA[i]);
      gRAASys->SetPoint(nCentBins-1-i,dNchdEta[i],RAA[i]);
      gRAASys->SetPointError(nCentBins-1-i,10,TMath::Sqrt(eRAA[i]*eRAA[i] + eRAA2[i]*eRAA2[i]));
    }else if (xAxisType==3){
      gRAA->SetPoint(nCentBins-1-i,dNchdEtaNorm[i],RAA[i]);
      gRAA->SetPointError(nCentBins-1-i,edNchdEtaNorm[i],renJPsi[i]*RAA[i]);
      gRAASys->SetPoint(nCentBins-1-i,dNchdEtaNorm[i],RAA[i]);
      gRAASys->SetPointError(nCentBins-1-i,0.4,TMath::Sqrt(eRAA[i]*eRAA[i] + eRAA2[i]*eRAA2[i]));
    }
  }
  
  if (xAxisType==0) {
    gRAA->GetXaxis()->Set(nCentBins, bin);
    for(Int_t i=nCentBins; i>=1; i--)
      gRAA->GetXaxis()->SetBinLabel(nCentBins+1-i, label[i].Data());
  } else if (xAxisType==1){
    gRAA->GetXaxis()->Set(40, 0., 400.);
  }else if (xAxisType==2){
    gRAA->GetXaxis()->Set(40, 0., 1550.);
  }else if (xAxisType==3){
    gRAA->GetXaxis()->Set(40, 0., 10.);
  }
  
  // define axis
  if (xAxisType==0) {
    gRAA->GetXaxis()->SetTitle("centrality");
    gRAA->GetXaxis()->SetLabelSize(0.07);
    gRAA->GetXaxis()->SetTitleOffset(1.1);
    gRAA->GetXaxis()->SetTitleSize(0.05);
    gRAA->GetXaxis()->SetTickLength(0.);
  } else if (xAxisType==1){
    gRAA->GetXaxis()->SetTitle("#LT N_{part} #GT");
    //gRAA->GetXaxis()->SetTitle("#LT N_{part} #GT weighted by N_{coll}");
    gRAA->GetXaxis()->SetTitleSize(0.05); 
    gRAA->GetXaxis()->SetTitleOffset(1.0);
    gRAA->GetXaxis()->SetLabelSize(0.05);
  }
  else if (xAxisType==2){
    gRAA->GetXaxis()->SetTitle("dN_{ch}/d#eta |_{ #eta = 0}");
    gRAA->GetXaxis()->SetTitleSize(0.05);
    gRAA->GetXaxis()->SetTitleOffset(1.1);
    gRAA->GetXaxis()->SetLabelSize(0.05);
  }
  else if (xAxisType==3){
    gRAA->GetXaxis()->SetTitle("(dN_{ch}/d#eta)/(#LT N_{part} #GT /2)");
    gRAA->GetXaxis()->SetTitleSize(0.05);
    gRAA->GetXaxis()->SetTitleOffset(1.1);
    gRAA->GetXaxis()->SetLabelSize(0.05);
  }
  gRAA->GetYaxis()->SetTitle("R_{AA}");
  gRAA->GetYaxis()->SetTitleSize(0.05);
  gRAA->GetYaxis()->SetTitleOffset(0.9);
  gRAA->GetYaxis()->SetLabelSize(0.05);
  gRAA->GetYaxis()->SetRangeUser(0.,1.4);
  
  // define display
  gRAA->SetMarkerStyle(gMarkerStyleAliceForw);
  gRAA->SetMarkerSize(gMarkerSizeAliceForw);
  gRAA->SetMarkerColor(color);
  gRAA->SetLineWidth(2);
  gRAA->SetLineColor(color);
  gRAASys->SetLineColor(color);
  gRAASys->SetFillStyle(0);
  gRAASysAll->SetLineWidth(1);
  gRAASysAll->SetFillStyle(3001);
  gRAASysAll->SetFillColor(color);
  
}


//------------------------------------------------------------------------
void BuildForwardPHENIXRAAGraphs(TGraphErrors *&gRAA_PHENIX, TGraphAsymmErrors *&gRAA_PHENIXSys,
				 TGraphErrors *&gRAA_PHENIXSysAll, Int_t xAxisType)
{
  // build the TGraphs to display PHENIX RAA results at forward rapidity
  
  // fill RAA graphs
  Double_t cent_PHENIX[17] = {2.5, 7.5, 12.5, 17.5, 22.5, 27.5, 32.5, 37.5, 42.5, 47.5, 52.5, 57.5, 62.5, 67.5, 72.5, 77.5, 86.};
  Double_t ecent_PHENIX[17] = {2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 6.};
  Double_t NPart_PHENIX[17] = {350.8, 301.7, 255.7, 216.4, 182.4, 152.7, 126.8, 104.2, 84.6, 67.7, 53.2, 41., 30.8, 22.6, 16.1, 11.2, 5.6};
  Double_t eNPart_PHENIX[17] = {3.1, 4.7, 5.4, 5.6, 5.7, 5.9, 5.9, 5.8, 5.6, 5.4, 5., 4.5, 3.9, 3.4, 2.8, 2.2, 0.8};
  Double_t dNchdEta_PHENIX[17] = {687,560,457,372,302,246,197,156,124,95.3,70.9,52.2,37.5,25.6,-100.0,-100.0, -100.0};
  Double_t edNchdEta_PHENIX[17] = {37,28,22,18,16,14,12,11,9.6,8.6,7.6,6.5,5.4,4.5,1,1,1};
  Double_t dNchdEtaNorm_PHENIX[17] = {3.89, 3.73, 3.59, 3.45, 3.34, 3.25, 3.15, 3.05, 2.96, 2.86, 2.7, 2.6, 2.48, 2.33, 1.0, 1.0, 1.0};
  Double_t edNchdEtaNorm_PHENIX[17] = {0.23,0.22,0.21,0.21,0.21,0.22,0.24,0.26,0.28,0.32,0.36,0.41,0.46,0.55,1,1,1};
  Double_t RAA_PHENIX[17] = {0.17, 0.16, 0.20, 0.25, 0.25, 0.35, 0.35, 0.41, 0.52, 0.49, 0.54, 0.80, 0.68, 0.72, 0.91, 1.03, 1.20};
  Double_t eRAA_PHENIX[17] = {0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.05, 0.06, 0.06, 0.08, 0.11, 0.10};
  Double_t syspRAA_PHENIX[17] = {0.04, 0.02, 0.03, 0.04, 0.03, 0.04, 0.04, 0.06, 0.07, 0.07, 0.09, 0.14, 0.13, 0.15, 0.21, 0.26, 0.23};
  Double_t sysmRAA_PHENIX[17] = {0.02, 0.02, 0.03, 0.04, 0.03, 0.04, 0.04, 0.06, 0.07, 0.07, 0.09, 0.14, 0.13, 0.15, 0.21, 0.26, 0.23};
  gRAA_PHENIX = new TGraphErrors(17);
  gRAA_PHENIXSys = new TGraphAsymmErrors(17);
  for (Int_t i=16; i>=0; i--) {
    if (xAxisType==0) {
      gRAA_PHENIX->SetPoint(16-i, 100.-cent_PHENIX[i], RAA_PHENIX[i]);
      gRAA_PHENIX->SetPointError(16-i, ecent_PHENIX[i], eRAA_PHENIX[i]);
      gRAA_PHENIXSys->SetPoint(16-i, 100.-cent_PHENIX[i], RAA_PHENIX[i]);
      gRAA_PHENIXSys->SetPointError(16-i, 1., 1., sysmRAA_PHENIX[i], syspRAA_PHENIX[i]);
    } else if (xAxisType==1){
      gRAA_PHENIX->SetPoint(16-i, NPart_PHENIX[i], RAA_PHENIX[i]);
      gRAA_PHENIX->SetPointError(16-i, eNPart_PHENIX[i], eRAA_PHENIX[i]);
      gRAA_PHENIXSys->SetPoint(16-i, NPart_PHENIX[i], RAA_PHENIX[i]);
      gRAA_PHENIXSys->SetPointError(16-i, 5., 5., sysmRAA_PHENIX[i], syspRAA_PHENIX[i]);
    }else if (xAxisType==2){
      gRAA_PHENIX->SetPoint(16-i, dNchdEta_PHENIX[i], RAA_PHENIX[i]);
      gRAA_PHENIX->SetPointError(16-i, edNchdEta_PHENIX[i], eRAA_PHENIX[i]);
      gRAA_PHENIXSys->SetPoint(16-i, dNchdEta_PHENIX[i], RAA_PHENIX[i]);
      gRAA_PHENIXSys->SetPointError(16-i, 10., 10., sysmRAA_PHENIX[i], syspRAA_PHENIX[i]);
    }else if (xAxisType==3){
      gRAA_PHENIX->SetPoint(16-i, dNchdEtaNorm_PHENIX[i], RAA_PHENIX[i]);
      gRAA_PHENIX->SetPointError(16-i, edNchdEtaNorm_PHENIX[i], eRAA_PHENIX[i]);
      gRAA_PHENIXSys->SetPoint(16-i, dNchdEtaNorm_PHENIX[i], RAA_PHENIX[i]);
      gRAA_PHENIXSys->SetPointError(16-i, .2, .2, sysmRAA_PHENIX[i], syspRAA_PHENIX[i]);
    }
  }
  gRAA_PHENIXSysAll = new TGraphErrors(1);
  if (xAxisType==0) {
    gRAA_PHENIXSysAll->SetPoint(0,96.,1.);
    gRAA_PHENIXSysAll->SetPointError(0,1.,0.092);
  } else if (xAxisType==1){
    gRAA_PHENIXSysAll->SetPoint(0,383.,1.);
    gRAA_PHENIXSysAll->SetPointError(0,5.,0.092);
  }else if (xAxisType==2){
    gRAA_PHENIXSysAll->SetPoint(0,1535,1.);
    gRAA_PHENIXSysAll->SetPointError(0,15,0.092);
  }else if (xAxisType==3){
    gRAA_PHENIXSysAll->SetPoint(0,9.,1.);
    gRAA_PHENIXSysAll->SetPointError(0,0.1,0.092);
  }
  
  // define axis
  if (xAxisType==0) {
    gRAA_PHENIX->GetXaxis()->SetTitle("1 - centrality");
    gRAA_PHENIX->GetXaxis()->Set(100., 0., 100.);
  } else if(xAxisType==1){
    //gRAA->GetXaxis()->SetTitle("<N_{part}>");
    gRAA_PHENIX->GetXaxis()->SetTitle("#LT N_{part} #GT");
    gRAA_PHENIX->GetXaxis()->Set(40, 0., 400.);
  } else if(xAxisType==2){
    gRAA_PHENIX->GetXaxis()->SetTitle("dN_{ch}/d#eta|_{ #eta = 0}");
    gRAA_PHENIX->GetXaxis()->Set(40, 0., 1550.);
  }else if(xAxisType==3){
    gRAA_PHENIX->GetXaxis()->SetTitle("(dN_{ch}/d#eta)/(<N_{part}>/2)");
    gRAA_PHENIX->GetXaxis()->Set(40, 0., 10.);
  }
  
  gRAA_PHENIX->GetXaxis()->SetTitleSize(0.06);
  gRAA_PHENIX->GetXaxis()->SetTitleOffset(0.95);
  gRAA_PHENIX->GetXaxis()->SetLabelSize(0.05);
  
  gRAA_PHENIX->GetYaxis()->SetTitle("R_{AA}");
  gRAA_PHENIX->GetYaxis()->SetTitleSize(0.07);
  gRAA_PHENIX->GetYaxis()->SetTitleOffset(0.9);
  gRAA_PHENIX->GetYaxis()->SetLabelSize(0.05);
  gRAA_PHENIX->GetYaxis()->SetRangeUser(0.,1.4);
  
  // define display
  gRAA_PHENIX->SetMarkerStyle(gMarkerStylePhenForw);
  gRAA_PHENIX->SetMarkerSize(gMarkerSizePhenForw);
  gRAA_PHENIX->SetMarkerColor(4);
  gRAA_PHENIX->SetLineWidth(2);
  gRAA_PHENIX->SetLineColor(4);
  gRAA_PHENIXSys->SetLineColor(4);
  gRAA_PHENIXSys->SetFillStyle(0);
  gRAA_PHENIXSysAll->SetLineWidth(1);
  gRAA_PHENIXSysAll->SetFillStyle(3001);
  gRAA_PHENIXSysAll->SetFillColor(4);
  
}


//------------------------------------------------------------------------
void BuildMidPHENIXRAAGraphs(TGraphErrors *&gRAA_PHENIX2, TGraphErrors *&gRAA_PHENIXSys2,
			     TGraphErrors *&gRAA_PHENIXSysAll2, Int_t xAxisType)
{
  // build the TGraphs to display PHENIX RAA results at forward rapidity
  
  // fill RAA graphs
  Double_t cent_PHENIX2[8] = {2.5, 7.5, 12.5, 17.5, 25., 35., 50., 76.5};
  Double_t ecent_PHENIX2[8] = {2.5, 2.5, 2.5, 2.5, 5., 5., 10., 16.5};
  Double_t NPart_PHENIX2[8] = {351.4, 299.0, 253.9, 215.3, 166.6, 114.2, 58.4, 14.5};
  Double_t eNPart_PHENIX2[8] = {2.9, 3.8, 4.3, 5.3, 5.4, 4.4, 4., 2.5};
  Double_t dNchdEta_PHENIX2[8] = {687,560,457,372,273,181,80,1};
  Double_t edNchdEta_PHENIX2[8] = {37,28,22,18,16,12,10,1};
  Double_t dNchdEtaNorm_PHENIX2[8] = {3.89, 3.73, 3.59, 3.45, 3.30, 3.1, 2.78, 1.0};
  Double_t edNchdEtaNorm_PHENIX2[8] = {0.23,0.22,0.21,0.21,0.3,0.3,0.5,1.0};
  Double_t RAA_PHENIX2[8] = {0.26, 0.34, 0.36, 0.45, 0.58, 0.58, 0.65, 0.74};
  Double_t eRAA_PHENIX2[8] = {0.05, 0.06, 0.06, 0.07, 0.07, 0.08, 0.07, 0.12};
  Double_t sysRAA_PHENIX2[8] = {0.04, 0.05, 0.05, 0.07, 0.08, 0.08, 0.1, 0.21};
  gRAA_PHENIX2 = new TGraphErrors(8);
  gRAA_PHENIXSys2 = new TGraphErrors(8);
  for (Int_t i=7; i>=0; i--) {
    if (xAxisType==0) {
      gRAA_PHENIX2->SetPoint(7-i, 100.-cent_PHENIX2[i], RAA_PHENIX2[i]);
      gRAA_PHENIX2->SetPointError(7-i, ecent_PHENIX2[i], eRAA_PHENIX2[i]);
      gRAA_PHENIXSys2->SetPoint(7-i, 100.-cent_PHENIX2[i], RAA_PHENIX2[i]);
      gRAA_PHENIXSys2->SetPointError(7-i, 1., sysRAA_PHENIX2[i]);
    } else if (xAxisType==1){
      gRAA_PHENIX2->SetPoint(7-i, NPart_PHENIX2[i], RAA_PHENIX2[i]);
      gRAA_PHENIX2->SetPointError(7-i, eNPart_PHENIX2[i], eRAA_PHENIX2[i]);
      gRAA_PHENIXSys2->SetPoint(7-i, NPart_PHENIX2[i], RAA_PHENIX2[i]);
      gRAA_PHENIXSys2->SetPointError(7-i, 5., sysRAA_PHENIX2[i]);
    } else if (xAxisType==2){
      gRAA_PHENIX2->SetPoint(7-i, dNchdEta_PHENIX2[i], RAA_PHENIX2[i]);
      gRAA_PHENIX2->SetPointError(7-i, edNchdEta_PHENIX2[i], eRAA_PHENIX2[i]);
      gRAA_PHENIXSys2->SetPoint(7-i, dNchdEta_PHENIX2[i], RAA_PHENIX2[i]);
      gRAA_PHENIXSys2->SetPointError(7-i, 10, sysRAA_PHENIX2[i]);
    }else if (xAxisType==3){
      gRAA_PHENIX2->SetPoint(7-i, dNchdEtaNorm_PHENIX2[i], RAA_PHENIX2[i]);
      gRAA_PHENIX2->SetPointError(7-i, edNchdEtaNorm_PHENIX2[i], eRAA_PHENIX2[i]);
      gRAA_PHENIXSys2->SetPoint(7-i, dNchdEtaNorm_PHENIX2[i], RAA_PHENIX2[i]);
      gRAA_PHENIXSys2->SetPointError(7-i, .2, sysRAA_PHENIX2[i]);
    }
  }
  gRAA_PHENIXSysAll2 = new TGraphErrors(1);
  if (xAxisType==0) {
    gRAA_PHENIXSysAll2->SetPoint(0,93.,1.);
    gRAA_PHENIXSysAll2->SetPointError(0,1.,0.12);
  } else if (xAxisType==1){
    gRAA_PHENIXSysAll2->SetPoint(0,371.,1.);
    gRAA_PHENIXSysAll2->SetPointError(0,5.,0.12);
  } else if (xAxisType==2){
    gRAA_PHENIXSysAll2->SetPoint(0,1475,1.);
    gRAA_PHENIXSysAll2->SetPointError(0,15,0.12);
  }else if (xAxisType==3){
    gRAA_PHENIXSysAll2->SetPoint(0,9.2,1.);
    gRAA_PHENIXSysAll2->SetPointError(0,.1,0.12);
  }
  
  // define axis
  if (xAxisType==0) {
    gRAA_PHENIX2->GetXaxis()->SetTitle("1 - centrality");
    gRAA_PHENIX2->GetXaxis()->Set(100., 0., 100.);
  } else if (xAxisType==1) {
    //gRAA->GetXaxis()->SetTitle("<N_{part}>");
    gRAA_PHENIX2->GetXaxis()->SetTitle("#LT N_{part}* #GT");
    gRAA_PHENIX2->GetXaxis()->Set(40, 0., 400.);
  } else if(xAxisType==2){
    gRAA_PHENIX2->GetXaxis()->SetTitle("dN_{ch}/d#eta|_{#eta = 0}");
    gRAA_PHENIX2->GetXaxis()->Set(40, 0., 1550.);
  }else if(xAxisType==3){
    gRAA_PHENIX2->GetXaxis()->SetTitle("(dN_{ch}/d#eta)/(<N_{part}>/2)");
    gRAA_PHENIX2->GetXaxis()->Set(40, 0., 10.);
  }
  gRAA_PHENIX2->GetXaxis()->SetTitleSize(0.05);
  gRAA_PHENIX2->GetXaxis()->SetTitleOffset(0.95);
  gRAA_PHENIX2->GetXaxis()->SetLabelSize(0.05);
  
  gRAA_PHENIX2->GetYaxis()->SetTitle("R_{AA}");
  gRAA_PHENIX2->GetYaxis()->SetTitleSize(0.05);
  gRAA_PHENIX2->GetYaxis()->SetTitleOffset(0.9);
  gRAA_PHENIX2->GetYaxis()->SetLabelSize(0.05);
  gRAA_PHENIX2->GetYaxis()->SetRangeUser(0., 1.4);
  
  // define display
  gRAA_PHENIX2->SetMarkerStyle(gMarkerStylePhenMid);
  gRAA_PHENIX2->SetMarkerSize(gMarkerSizePhenMid);
  gRAA_PHENIX2->SetMarkerColor(kGreen+2);
  gRAA_PHENIX2->SetLineWidth(2);
  gRAA_PHENIX2->SetLineColor(kGreen+2);
  gRAA_PHENIXSys2->SetLineColor(kGreen+2);
  gRAA_PHENIXSys2->SetFillStyle(0);
  gRAA_PHENIXSysAll2->SetLineWidth(1);
  gRAA_PHENIXSysAll2->SetFillStyle(3001);
  gRAA_PHENIXSysAll2->SetFillColor(kGreen+2);
  
}


//------------------------------------------------------------------------
void DisplayRAA(Int_t pt, Int_t y, Int_t nCentBins, TGraphErrors *gRAA, TGraphErrors *gRAASys, TGraphErrors *gRAASysAll, Int_t xAxisType)
{
  // display our RAA results
  
  new TCanvas(Form("RAA_pt%d_y%d",pt,y),Form("RAA -- pt%d -- y%d",pt,y),gCanvasWidth,gCanvasHeigh);
  gPad->SetTopMargin(0.05);
  gPad->SetBottomMargin(0.12);
  gPad->SetRightMargin(0.05);
  gRAA->Draw("ap");
  gRAASys->Draw("e2");
  gRAASysAll->Draw("e2");
  TLine *l = 0x0;  
  if (xAxisType==0) l = new TLine(0., 1., 100., 1.);
  else if (xAxisType==1) l =  new TLine(0., 1., 400., 1.);
  else if (xAxisType==2) l =  new TLine(0., 1., 1550., 1.);
  else if (xAxisType==3) l =  new TLine(0., 1., 10., 1.);
  l->SetLineStyle(3);
  l->Draw();
  if (xAxisType==0) {
    if (nCentBins == 6) {
      TLine *l1 = new TLine(20, 0, 20, 0.03);
      l1->Draw();
      TLine *l2 = new TLine(50, 0, 50, 0.03);
      l2->Draw();
      TLine *l3 = new TLine(70, 0, 70, 0.03);
      l3->Draw();
      TLine *l4 = new TLine(80, 0, 80, 0.03);
      l4->Draw();
      TLine *l5 = new TLine(90, 0, 90, 0.03);
      l5->Draw();
    }
  }
  TLegend* tlg = new TLegend(0.13, 0.82, 0.9, 0.90,"","NDC");
  if (y == 0 && pt == 0) tlg->AddEntry(gRAA,"ALICE (Pb-Pb #sqrt{s_{NN}} = 2.76 TeV), 2.5<y<4, p_{T}>0 ","p");
  else if (y == 1 && pt == 0) tlg->AddEntry(gRAA,"inclusive J/#psi in Pb-Pb  #sqrt{s_{NN}} = 2.76 TeV, 3.25<y<4, p_{T}>0","p");
  else if (y == 2 && pt == 0) tlg->AddEntry(gRAA,"inclusive J/#psi in Pb-Pb  #sqrt{s_{NN}} = 2.76 TeV, 2.5<y<3.25, p_{T}>0","p");
  else if (y == 0 && pt == 1) tlg->AddEntry(gRAA,"inclusive J/#psi in Pb-Pb  #sqrt{s_{NN}} = 2.76 TeV, 2.5<y<4, 0<p_{T}<3 GeV/c","p");
  else if (y == 0 && pt == 2) tlg->AddEntry(gRAA,"inclusive J/#psi in Pb-Pb  #sqrt{s_{NN}} = 2.76 TeV, 2.5<y<4, p_{T}>3 GeV/c","p");
  tlg->SetFillStyle(0);
  tlg->SetBorderSize(0);
  tlg->SetTextFont(42);
  tlg->SetMargin(0.05);
  tlg->Draw();
  //ALICEseal("ALICE preliminary", 0.2, 0.2);
  //ALICEseal("", 0.2, 0.2);
  
}


//------------------------------------------------------------------------
void DisplayRAA(TString var, Int_t nCentBins, TGraphErrors *gRAA[3], TGraphErrors *gRAASys[3], TGraphErrors *gRAASysAll[3], Int_t xAxisType)
{
  // display our RAA results
  
  if (var == "y") new TCanvas("RAA_pt0_yAll","RAA -- pt0 -- yAll",gCanvasWidth,gCanvasHeigh);
  else new TCanvas("RAA_ptAll_y0","RAA -- ptAll -- y0",gCanvasWidth,gCanvasHeigh);
  gPad->SetTopMargin(0.05);
  gPad->SetBottomMargin(0.12);
  gPad->SetRightMargin(0.05);
  for (Int_t i=0; i<3; i++) {
    if (i==0) gRAA[i]->Draw("ap");
    else gRAA[i]->Draw("p");
    gRAASys[i]->Draw("e2");
    gRAASysAll[i]->Draw("e2");
  }
  TLine *l = 0x0; 
  if  (xAxisType==0)  l = new TLine(0., 1., 100., 1.);
  else if (xAxisType==1) l = new TLine(0., 1., 400., 1.);
  else if (xAxisType==2) l= new TLine(0., 1., 1550., 1.);
  else if (xAxisType==3) l= new TLine(0., 1., 10., 1.);
  l->SetLineStyle(3);
  l->Draw();
  if (xAxisType==0) {
    if (nCentBins == 6) {
      TLine *l1 = new TLine(20, 0, 20, 0.03);
      l1->Draw();
      TLine *l2 = new TLine(50, 0, 50, 0.03);
      l2->Draw();
      TLine *l3 = new TLine(70, 0, 70, 0.03);
      l3->Draw();
      TLine *l4 = new TLine(80, 0, 80, 0.03);
      l4->Draw();
      TLine *l5 = new TLine(90, 0, 90, 0.03);
      l5->Draw();
    }
  }
  TLegend *lg = new TLegend(0.17, 0.81, 0.92, 0.93,"","NDC");
  lg->AddEntry(gRAA[0],"ALICE (Pb-Pb #sqrt{s_{NN}} = 2.76 TeV), 2.5<y<4, p_{T}>0 ","p");
  if (var == "y") {
    lg->AddEntry(gRAA[1],"ALICE (Pb-Pb #sqrt{s_{NN}} = 2.76 TeV), 3.25<y<4, p_{T}>0 ","p");
    lg->AddEntry(gRAA[2],"ALICE (Pb-Pb #sqrt{s_{NN}} = 2.76 TeV), 2.5<y<3.25, p_{T}>0 ","p");
  } else {
    lg->AddEntry(gRAA[1],"ALICE (Pb-Pb #sqrt{s_{NN}} = 2.76 TeV), 2.5<y<4, 0<p_{T}<3 GeV/c ","p");
    lg->AddEntry(gRAA[2],"ALICE (Pb-Pb #sqrt{s_{NN}} = 2.76 TeV), 2.5<y<4, p_{T}>3 GeV/c ","p");
  }
  lg->SetFillStyle(0);
  lg->SetBorderSize(0);
  lg->SetTextFont(42);
  lg->SetMargin(0.05);
  lg->Draw();
  //ALICEseal("ALICE preliminary", 0.2, 0.2);
  ALICEseal("", 0.2, 0.2);
  
}


//------------------------------------------------------------------------
void DisplayRAAWithPHENIX(Int_t pt, Int_t y, TGraphErrors *gRAA, TGraphErrors *gRAASys, TGraphErrors *gRAASysAll, Int_t xAxisType)
{
  // display our RAA results together with PHENIX results at forward rapidity
  
  // Get PHENIX graphs at forward rapidity
  TGraphErrors *gRAA_PHENIX, *gRAA_PHENIXSysAll;
  TGraphAsymmErrors *gRAA_PHENIXSys;
  BuildForwardPHENIXRAAGraphs(gRAA_PHENIX, gRAA_PHENIXSys, gRAA_PHENIXSysAll, xAxisType);
  
  // plot
  new TCanvas(Form("RAAvsPHENIX_pt%d_y%d",pt,y),Form("RAAvsPHENIX -- pt%d -- y%d",pt,y),gCanvasWidth,gCanvasHeigh);
  gPad->SetTopMargin(0.05);
  gPad->SetBottomMargin(0.12);
  gPad->SetRightMargin(0.05);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.03);
  gPad->SetTopMargin(0.03);
  gPad->SetBottomMargin(0.13);
  
  gRAA_PHENIX->Draw("ap");
  gRAA_PHENIXSys->Draw("e2");
  gRAA_PHENIXSysAll->Draw("e2");
  gRAA->Draw("p");
  gRAASys->Draw("e2");
  gRAASysAll->Draw("e2");
  TLine *l = 0x0; 
  if (xAxisType==0) l = new TLine(0., 1., 100., 1.) ;
  else if (xAxisType==1) l = new TLine(0., 1., 400., 1.);
  else if (xAxisType==2) l = new TLine(0., 1., 1550., 1.);
  else if (xAxisType==3) l = new TLine(0., 1., 10., 1.);
  l->SetLineStyle(3);
  l->Draw();
  
  TLegend *lg = new TLegend(0.19, 0.77, 0.89, 0.95,"","NDC");
  if (y == 0 && pt == 0) lg->AddEntry(gRAA,"ALICE (Pb-Pb #sqrt{s_{NN}} = 2.76 TeV), 2.5<y<4, p_{T}>0 (preliminary)","p");
  else if (y == 1 && pt == 0) lg->AddEntry(gRAA,"ALICE (Pb-Pb #sqrt{s_{NN}} = 2.76 TeV), 3.25<y<4, p_{T}>0 ","p");
  else if (y == 2 && pt == 0) lg->AddEntry(gRAA,"ALICE (Pb-Pb #sqrt{s_{NN}} = 2.76 TeV), 2.5<y<3.25, p_{T}>0 ","p");
  else if (y == 0 && pt == 1) lg->AddEntry(gRAA,"ALICE (Pb-Pb #sqrt{s_{NN}} = 2.76 TeV), 2.5<y<4, 0<p_{T}<3 GeV/c ","p");
  else if (y == 0 && pt == 2) lg->AddEntry(gRAA,"ALICE (Pb-Pb #sqrt{s_{NN}} = 2.76 TeV), 2.5<y<4, p_{T}>3 GeV/c ","p");
  if (PRL) lg->AddEntry(gRAA_PHENIX,"PHENIX (Au-Au #sqrt{s_{NN}} = 0.2 TeV), 1.2<|y|<2.2, p_{T}>0","p");
  else  lg->AddEntry(gRAA_PHENIX,"PHENIX (Au-Au #sqrt{s_{NN}} = 0.2 TeV), 1.2<|y|<2.2, p_{T}>0 (arXiv:1103.6269)","p");
  lg->SetFillStyle(0);
  lg->SetBorderSize(0);
  lg->SetTextFont(42);
  lg->SetMargin(0.05);
  lg->Draw();
  if (xAxisType==1 && !PRL) {
    TPaveText *t3 = new TPaveText(0.2, 0.15, 0.6, 0.21,"NDC");
    t3->AddText("(*) ALICE #LT N_{part} #GT is weighted by N_{coll}");
    t3->SetFillStyle(0);
    t3->SetBorderSize(0);
    t3->SetTextFont(42);
    t3->Draw();
  }
  
}


//------------------------------------------------------------------------
void DisplayRAAWithPHENIX2(Int_t pt, Int_t y, TGraphErrors *gRAA, TGraphErrors *gRAASys, TGraphErrors *gRAASysAll, Int_t xAxisType)
{
  // display our RAA results together with PHENIX results at forward and mid-rapidity
  
  // Get PHENIX graphs at forward rapidity
  TGraphErrors *gRAA_PHENIX, *gRAA_PHENIXSysAll;
  TGraphAsymmErrors *gRAA_PHENIXSys;
  BuildForwardPHENIXRAAGraphs(gRAA_PHENIX, gRAA_PHENIXSys, gRAA_PHENIXSysAll, xAxisType);
  
  // Get PHENIX graphs at mid-rapidity
  TGraphErrors *gRAA_PHENIX2, *gRAA_PHENIXSys2, *gRAA_PHENIXSysAll2;
  BuildMidPHENIXRAAGraphs(gRAA_PHENIX2, gRAA_PHENIXSys2, gRAA_PHENIXSysAll2, xAxisType);
  
  // plot RAA
  new TCanvas(Form("RAAvsPHENIX2_pt%d_y%d",pt,y),Form("RAAvsPHENIX2 -- pt%d -- y%d",pt,y),gCanvasWidth,gCanvasHeigh);
  gPad->SetTopMargin(0.05);
  gPad->SetBottomMargin(0.12);
  gPad->SetRightMargin(0.05);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.03);
  gPad->SetTopMargin(0.03);
  gPad->SetBottomMargin(0.13);
  
  gRAA_PHENIX->Draw("ap");
  gRAA_PHENIXSys->Draw("e2");
  gRAA_PHENIXSysAll->Draw("e2");
  gRAA_PHENIX2->Draw("p");
  gRAA_PHENIXSys2->Draw("e2");
  gRAA_PHENIXSysAll2->Draw("e2");
  gRAA->Draw("p");
  gRAASys->Draw("e2");
  gRAASysAll->Draw("e2");
  TLine *l = 0x0;
  if (xAxisType==0) l = new TLine(0., 1., 100., 1.); 
  else if (xAxisType==1)l = new TLine(0., 1., 400., 1.);
  else if (xAxisType==2)l = new TLine(0., 1., 1550., 1.);
  else if (xAxisType==3)l = new TLine(0., 1., 10., 1.);
  l->SetLineStyle(3);
  l->Draw();
  TLegend *lg = new TLegend(0.17, 0.78, 0.95, 0.93,"","NDC");
  if (y == 0 && pt == 0) lg->AddEntry(gRAA,"ALICE (Pb-Pb #sqrt{s_{NN}} = 2.76 TeV), 2.5<y<4, p_{T}>0 ","p");
  else if (y == 1 && pt == 0) lg->AddEntry(gRAA,"ALICE (Pb-Pb #sqrt{s_{NN}} = 2.76 TeV), 3.25<y<4, p_{T}>0 ","p");
  else if (y == 2 && pt == 0) lg->AddEntry(gRAA,"ALICE (Pb-Pb #sqrt{s_{NN}} = 2.76 TeV), 2.5<y<3.25, p_{T}>0 ","p");
  else if (y == 0 && pt == 1) lg->AddEntry(gRAA,"ALICE (Pb-Pb #sqrt{s_{NN}} = 2.76 TeV), 2.5<y<4, 0<p_{T}<3 GeV/c ","p");
  else if (y == 0 && pt == 2) lg->AddEntry(gRAA,"ALICE (Pb-Pb #sqrt{s_{NN}} = 2.76 TeV), 2.5<y<4, p_{T}>3 GeV/c ","p");
  //lg->AddEntry(gRAA_PHENIX,"PHENIX (Au-Au #sqrt{s_{NN}} = 0.2 TeV), 1.2<|y|<2.2, p_{T}>0 (arXiv:1103.6269)","p");
  //lg->AddEntry(gRAA_PHENIX2,"PHENIX (Au-Au #sqrt{s_{NN}} = 0.2 TeV), |y|<0.35, p_{T}>0 (nucl-ex/0611020)","p");
  lg->AddEntry(gRAA_PHENIX,"PHENIX (Au-Au #sqrt{s_{NN}} = 0.2 TeV), 1.2<|y|<2.2, p_{T}>0","p");
  lg->AddEntry(gRAA_PHENIX2,"PHENIX (Au-Au #sqrt{s_{NN}} = 0.2 TeV), |y|<0.35, p_{T}>0","p");
  lg->SetFillStyle(0);
  lg->SetBorderSize(0);
  lg->SetTextFont(42);
  lg->SetMargin(0.05);
  lg->Draw();
  if (xAxisType==1) {
    TPaveText *t3 = new TPaveText(0.2, 0.15, 0.6, 0.25,"NDC");
    t3->AddText("(*) ALICE #LT N_{part} #GT is weighted by N_{coll}");
    t3->SetFillStyle(0);
    t3->SetBorderSize(0);
    t3->SetTextFont(42);
    t3->Draw();
  }
  
}


//------------------------------------------------------------------------
void ALICEseal(TString type, Double_t xPad, Double_t yPad)
{
  TVirtualPad* currPad = gPad;
  TPad *myPadLogo = new TPad("myPadLogo", "Pad for ALICE Logo",xPad,yPad,xPad+0.17,yPad+0.17);
  myPadLogo->SetFillColor(0);
  myPadLogo->SetBorderMode(0);
  myPadLogo->SetBorderSize(2);
  myPadLogo->SetFrameBorderMode(0);
  myPadLogo->SetLeftMargin(0.0);
  myPadLogo->SetTopMargin(0.0);
  myPadLogo->SetBottomMargin(0.0);
  myPadLogo->SetRightMargin(0.0);
  myPadLogo->Draw();
  myPadLogo->cd();
  TASImage *myAliceLogo = new TASImage("/Users/pillot/Pictures/alice_logo.png");
  myAliceLogo->Draw();
  currPad->cd();
  Double_t x1 = xPad - 0.07, y1 = yPad - 0.06;
  Double_t x2 = x1 + 0.25, y2 = y1 + 0.08;
  TPaveText* t1=new TPaveText(x1,y1,x2,y2,"NDC");
  t1->SetFillStyle(0);
  t1->SetBorderSize(0);
  t1->AddText(0.,0.,Form("%s", type.Data()));
  t1->SetTextColor(kRed);
  t1->SetTextFont(42);
  t1->Draw();
  TPaveText* t2=new TPaveText(x1+0.06,y1-0.06,x2-0.06,y2-0.06,"NDC");
  t2->SetFillStyle(0);
  t2->SetBorderSize(0);
  t2->SetTextColor(kRed);
  t2->SetTextFont(52);
  //TDatime dt;
  //TString today = Form("%02i/%02i/%4i", dt.GetDay(), dt.GetMonth(), dt.GetYear());
  //t2->AddText(0.,0.,today.Data());
  //t2->Draw();
}

