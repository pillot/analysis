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

// y=-1: loop over the 4 choices below
// y=0: -4<y<-2.5
// y=1: -3.25<y<-2.5
// y=2: -4<y<-3.25
// y=3: y bins cent 0-90%

// pt=-1: loop over the 8 choices below
// pt=0: 0<pt<8 GeV/c
// pt=1: 0<pt<2 GeV/c
// pt=2: 2<pt<5 GeV/c
// pt=3: 5<pt<8 GeV/c
// pt=4: pt bins cent 0-90%
// pt=5: pt bins cent 0-20%
// pt=6: pt bins cent 20-40%
// pt=7: pt bins cent 40-90%


// branching ratio
Double_t BR=0.059;
Double_t eBR=0.01; // 1%

// luminosity in p-p
Double_t eppLumi=0.0355; // 1.9% + 3% du R-factor = 3.55%
Double_t eppEff=0.045; // 2% trigger eff + 4% reco eff = 4.5%

// factor normalization nMUL -> nMB
//Double_t nbMUL[10]={17409527., 9492252., 4459693., 2040116., 872857., 347547., 128538., 45868., 16545., 6111.}; // nb MUL AOD101
//Double_t fnorm090 = 27.28;
Double_t errfnorm=0.021;  // error syst of the fnorm factor for 2011 data
//nbMB[i] = nbMUL[i] * ( fnorm090 / (9. * nbMUL[i]/nbMUL[0]) );


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


void SetGraphNameTitleVsCent(Int_t nBins, TString *&gName, TString *&gTitle);
void SetGraphLabelsVsCent(Int_t nBins, Double_t *&x, Double_t *&ex, Double_t *&bin, TString *&label,
						  Double_t *&NPart, Double_t *&eNPart, Double_t *&dNchdEta, Double_t *&edNchdEta,
						  Double_t *&dNchdEtaNorm, Double_t *&edNchdEtaNorm);
void SetGraphLabelsVsPtY(Int_t nBins, Double_t *&x, Double_t *&ex);

void SetNormVar(Int_t pt, Int_t y, Double_t &sigmaJPsipp, Double_t &esigmaJPsippStat, Double_t &esigmaJPsippSyst);
void SetNormVarPtY(Int_t pt, Int_t y, Double_t *&sigmaJPsipp, Double_t *&esigmaJPsippStat, Double_t *&esigmaJPsippSyst);
void SetNormVarVsCent(Int_t pt, Int_t y, Int_t nBins, Double_t *&nMB, Double_t *&TAA, Double_t *&eTAA, Double_t *&AccEff, Double_t *&eAccEff);
void SetNormVarVsPtY(Int_t pt, Int_t y, Double_t &nMB, Double_t &TAA, Double_t &eTAA, Double_t *&AccEff, Double_t *&eAccEff);
void SetSystVarVsCent(Int_t nBins, Double_t *&eGen, Double_t *&eRec, Double_t *&eTrk, Double_t *&eTrg);
void SetSystVarVsPtY(Int_t pt, Int_t y, Double_t *&eGen, Double_t *&eRec, Double_t *&eTrk, Double_t *&eTrg);

void ExtractNJPsiVsCent(Int_t pt, Int_t y, Int_t nBins, Int_t nT, TString SigFct, Double_t **&n, Double_t **&en);
void ExtractNJPsiVsPtY(Int_t pt, Int_t y, Int_t nBins, Int_t nT, TString SigFct, TString Cent, Double_t **&n, Double_t **&en);

void SetNJPsiVsCent(Double_t **&nJPsi, Double_t **&enJPsi, Int_t &nTests, Int_t &nBins);
void SetNJPsiVsCentpt1(Double_t **&nJPsi, Double_t **&enJPsi, Int_t &nTests, Int_t &nBins);
void SetNJPsiVsCentpt2(Double_t **&nJPsi, Double_t **&enJPsi, Int_t &nTests, Int_t &nBins);
void SetNJPsiVsCentpt3(Double_t **&nJPsi, Double_t **&enJPsi, Int_t &nTests, Int_t &nBins);
void SetNJPsiVsCenty1(Double_t **&nJPsi, Double_t **&enJPsi, Int_t &nTests, Int_t &nBins);
void SetNJPsiVsCenty2(Double_t **&nJPsi, Double_t **&enJPsi, Int_t &nTests, Int_t &nBins);

void SetNJPsiVsPt4(Double_t **&nJPsi, Double_t **&enJPsi, Int_t &nTests, Int_t &nBins);
void SetNJPsiVsPt5(Double_t **&nJPsi, Double_t **&enJPsi, Int_t &nTests, Int_t &nBins);
void SetNJPsiVsPt6(Double_t **&nJPsi, Double_t **&enJPsi, Int_t &nTests, Int_t &nBins);
void SetNJPsiVsPt7(Double_t **&nJPsi, Double_t **&enJPsi, Int_t &nTests, Int_t &nBins);
void SetNJPsiVsY3(Double_t **&nJPsi, Double_t **&enJPsi, Int_t &nTests, Int_t &nBins);

void ComputeNJpsiMean(Int_t pt, Int_t y, Double_t *&nJPsiMean, Double_t *&enJPsiMean, Double_t *&sysnJPsiMean, Int_t &nBins, Bool_t weightedMeanRMS);
void ComputeRAA(Int_t pt, Int_t y, Int_t nBins, Double_t *nJPsi, Double_t *enJPsi, Double_t *sysnJPsi, Double_t *&renJPsi,
		Double_t *&RAA, Double_t *&eRAAstat, Double_t *&eRAA, Double_t *&eRAA2, Double_t &eRAACorr, Bool_t print = kFALSE);

void BuildRAAGraphs(Int_t nBins, Double_t *renJPsi, Double_t *RAA, Double_t *eRAAstat, Double_t *eRAA, Double_t *eRAA2, Double_t eRAACorr,
		    TGraphErrors *&gRAA, TGraphErrors *&gRAASys, TGraphErrors *&gRAASysAll, Int_t xAxisType, Int_t color = 2, Int_t sysPos = 0);
void BuildForwardPHENIXRAAGraphs(TGraphErrors *&gRAA_PHENIX, TGraphAsymmErrors *&gRAA_PHENIXSys,
				 TGraphErrors *&gRAA_PHENIXSysAll, Int_t xAxisType);
void BuildMidPHENIXRAAGraphs(TGraphErrors *&gRAA_PHENIX2, TGraphErrors *&gRAA_PHENIXSys2,
			     TGraphErrors *&gRAA_PHENIXSysAll2, Int_t xAxisType);

void DisplayRAA(Int_t pt, Int_t y, Int_t nBins, TGraphErrors *gRAA, TGraphErrors *gRAASys, TGraphErrors *gRAASysAll, Int_t xAxisType);
void DisplayRAA(TString var, Int_t nBins, TGraphErrors *gRAA[8], TGraphErrors *gRAASys[8], TGraphErrors *gRAASysAll[8], Int_t xAxisType);
void DisplayRAAratio(TString var, Int_t nBins, TGraphErrors *gRAA[8], TGraphErrors *gRAASys[8], TGraphErrors *gRAASysAll[8]);

void DisplayRAAWithPHENIX(Int_t pt, Int_t y, TGraphErrors *gRAA, TGraphErrors *gRAASys, TGraphErrors *gRAASysAll, Int_t xAxisType);
void DisplayRAAWithPHENIX2(Int_t pt, Int_t y, TGraphErrors *gRAA, TGraphErrors *gRAASys, TGraphErrors *gRAASysAll, Int_t xAxisType);
void ALICEseal(TString type, Double_t xPad, Double_t yPad);




//------------------------------------------------------------------------
void Raa_2011_12test(Int_t y = -1, Int_t pt = -1, Bool_t weightedMeanRMS = kTRUE, Int_t xAxisType = 1)
{
  
  if (y < -1 || y > 3 || pt < -1 || pt > 8 || (y > 0 && pt > 0) ) {
    ::Error("Raa_2011", "\ninvalid choice of pt or rapidity bin\n\n");
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
  Int_t color[8] = {kRed, kBlue, kGreen+2, kOrange, kViolet, kMagenta-4, kBlue+2, kSpring+10};
  
  
  
  if (y == -1) {
    
    Int_t nBins;
    TGraphErrors *gRAA[4], *gRAASys[4], *gRAASysAll[4];
    
    for (Int_t i=0; i<4; i++) {
      
      Double_t *nJPsiMean, *enJPsiMean, *sysnJPsiMean;
      ComputeNJpsiMean(0, i, nJPsiMean, enJPsiMean, sysnJPsiMean, nBins, weightedMeanRMS);
      
      Double_t eRAACorr;
      Double_t *renJPsiMean, *RAA, *eRAAstat, *eRAA, *eRAA2;
      ComputeRAA(0, i, nBins, nJPsiMean, enJPsiMean, sysnJPsiMean, renJPsiMean, RAA, eRAAstat, eRAA, eRAA2, eRAACorr, kTRUE);
      
      BuildRAAGraphs(nBins, renJPsiMean, RAA, eRAAstat, eRAA, eRAA2, eRAACorr, gRAA[i], gRAASys[i], gRAASysAll[i], xAxisType, color[i], 0); // i
	  
	  TGraphErrors *gRAAi=gRAA[i], *gRAASysi=gRAASys[i], *gRAASysAlli=gRAASysAll[i];
	  DisplayRAA(0, i, nBins, gRAAi, gRAASysi, gRAASysAlli, xAxisType);
	  
	  if (i==0) DisplayRAAWithPHENIX(0, 0, gRAAi, gRAASysi, gRAASysAlli, xAxisType);
    }
    
    DisplayRAA("y", nBins, gRAA, gRAASys, gRAASysAll, xAxisType);
	DisplayRAAratio("y", 7, gRAA, gRAASys, gRAASysAll);
	
  }
  
  if (pt == -1) {
    
    Int_t nBins;
    TGraphErrors *gRAA[8], *gRAASys[8], *gRAASysAll[8];
    
    for (Int_t i=1; i<9; i++) {
      
      Double_t *nJPsiMean, *enJPsiMean, *sysnJPsiMean;
      ComputeNJpsiMean(i, 0, nJPsiMean, enJPsiMean, sysnJPsiMean, nBins, weightedMeanRMS);
      
      Double_t eRAACorr;
      Double_t *renJPsiMean, *RAA, *eRAAstat, *eRAA, *eRAA2;
      ComputeRAA(i, 0, nBins, nJPsiMean, enJPsiMean, sysnJPsiMean, renJPsiMean, RAA, eRAAstat, eRAA, eRAA2, eRAACorr, kTRUE);
      
      BuildRAAGraphs(nBins, renJPsiMean, RAA, eRAAstat, eRAA, eRAA2, eRAACorr, gRAA[i-1], gRAASys[i-1], gRAASysAll[i-1], xAxisType, color[i-1], 0); // i
      
	  TGraphErrors *gRAAi=gRAA[i-1], *gRAASysi=gRAASys[i-1], *gRAASysAlli=gRAASysAll[i-1];
	  DisplayRAA(i, 0, nBins, gRAAi, gRAASysi, gRAASysAlli, xAxisType);
    }
    
    DisplayRAA("pt", nBins, gRAA, gRAASys, gRAASysAll, xAxisType);
	DisplayRAA("cent", nBins, gRAA, gRAASys, gRAASysAll, xAxisType);
	DisplayRAAratio("pt", 7, gRAA, gRAASys, gRAASysAll);
	
  }
  
}



//------------------------------------------------------------------------
void SetGraphNameTitleVsCent(Int_t nBins, TString *&gName, TString *&gTitle)
{
  // define graph name and title according to the number of centrality bins
  
  if (nBins == 8) {
    
    gName = new TString[8];
    gName[0] = "090"; gName[1] = "010"; gName[2] = "1020"; gName[3] = "2030"; gName[4] = "3040"; gName[5] = "4050"; gName[6] = "5060"; gName[7] = "6090";
    gTitle = new TString[8];
    gTitle[0] = "0-90%"; gTitle[1] = "0-10%"; gTitle[2] = "10-20%"; gTitle[3] = "20-30%"; gTitle[4] = "30-40%"; gTitle[5] = "40-50%"; gTitle[6] = "50-60%"; gTitle[7] = "60-90%";
    
  } else if (nBins == 10) {
    
    gName = new TString[10];
    gName[0] = "090"; gName[1] = "010"; gName[2] = "1020"; gName[3] = "2030"; gName[4] = "3040"; gName[5] = "4050"; gName[6] = "5060"; gName[7] = "6070"; gName[8] = "7080"; gName[9] = "8090";
    gTitle = new TString[10];
    gTitle[0] = "0-90%"; gTitle[1] = "0-10%"; gTitle[2] = "10-20%"; gTitle[3] = "20-30%"; gTitle[4] = "30-40%"; gTitle[5] = "40-50%"; gTitle[6] = "50-60%"; gTitle[7] = "60-70%"; gTitle[8] = "70-80%"; gTitle[9] = "80-90%";
    
  } else {
    
    ::Error("SetGraphNameTitleVsCent", "invalid number of centrality bins");
    exit(0);
    
  }
  
}

//------------------------------------------------------------------------
void SetGraphNameTitleVsPtY(Int_t nBins, TString *&gName, TString *&gTitle)
{
  // define graph name and title according to the number of centrality bins
  
  if (nBins == 6) {
    
    gName = new TString[6];
    gName[0] = "y1"; gName[1] = "y2"; gName[2] = "y3"; gName[3] = "y4"; gName[4] = "y5"; gName[5] = "y6";
    gTitle = new TString[6];
    gTitle[0] = "-2.75<y<-2.5"; gTitle[1] = "-3<y<2.75"; gTitle[2] = "-3.25<y<-3"; gTitle[3] = "-3.5<y<-3.25"; gTitle[4] = "-3.75<y<-3.5"; gTitle[5] = "-4<y<-3.75";
    
  } else if (nBins == 7) {
    
    gName = new TString[7];
    gName[0] = "pt1"; gName[1] = "pt2"; gName[2] = "pt3"; gName[3] = "pt4"; gName[4] = "pt5"; gName[5] = "pt6"; gName[6] = "pt7";
    gTitle = new TString[7];
    gTitle[0] = "0<pt<1"; gTitle[1] = "1<pt<2"; gTitle[2] = "2<pt<3"; gTitle[3] = "3<pt<4"; gTitle[4] = "4<pt<5"; gTitle[5] = "5<pt<6"; gTitle[6] = "6<pt<8";
    
  } else if (nBins == 3) {
    
    gName = new TString[3];
    gName[0] = "pt1"; gName[1] = "pt2"; gName[2] = "pt3";
    gTitle = new TString[3];
    gTitle[0] = "0<pt<2"; gTitle[1] = "2<pt<4"; gTitle[2] = "4<pt<8";
    
  } else {
    
    ::Error("SetGraphNameTitleVsPtY", "invalid number of centrality bins");
    exit(0);
    
  }
  
}

//------------------------------------------------------------------------
void SetGraphLabelsVsCent(Int_t nBins, Double_t *&x, Double_t *&ex, Double_t *&bin, TString *&label,
						  Double_t *&NPart, Double_t *&eNPart, Double_t *&dNchdEta, Double_t *&edNchdEta,
						  Double_t *&dNchdEtaNorm, Double_t *&edNchdEtaNorm)
{
  // define graph x-axis labels according to the number of centrality bins
  
  if (nBins == 8) {
    
    // plot versus centrality bin
    x = new Double_t[8];
    x[0] = 100.; x[1] = 5.; x[2] = 15.; x[3] = 25.; x[4] = 35.; x[5] = 45.; x[6] = 55.; x[7] = 75.;
    ex = new Double_t[8];
    ex[0] = 1.; ex[1] = 5.; ex[2] = 5.; ex[3] = 5.; ex[4] = 5.; ex[5] = 5.; ex[6] = 5.; ex[7] = 15.;
    bin = new Double_t[8];
    bin[0] = 0.; bin[1] = 10.; bin[2] = 20.; bin[3] = 30.; bin[4] = 40.; bin[5] = 50.; bin[6] = 60.; bin[7] = 90.;
    label = new TString[9];
    label[0] = "0-90%"; label[1] = "0-10%"; label[2] = "10-20%"; label[3] = "20-30%"; label[4] = "30-40%"; label[5] = "40-50%"; label[6] = "50-60%"; label[7] = "60-90%"; label[8] = "90-100%";	
    // plot versus NPart
    NPart = new Double_t[8]; // weighted using Alberica's macro
    NPart[0] = 124.3; NPart[1] = 356.5; NPart[2] = 260.5; NPart[3] = 186.4; NPart[4] = 128.9; NPart[5] = 85.0; NPart[6] = 52.8; NPart[7] = 17.55;
    eNPart = new Double_t[8];
    eNPart[0] = 2.3; eNPart[1] = 3.6; eNPart[2] = 4.4; eNPart[3] = 3.9; eNPart[4] = 3.3; eNPart[5] = 2.6; eNPart[6] = 2.0; eNPart[7] = 0.72;
	
  } else if (nBins == 10) {
    
    // plot versus centrality bin
    x = new Double_t[10];
    x[0] = 100.; x[1] = 5.; x[2] = 15.; x[3] = 25.; x[4] = 35.; x[5] = 45.; x[6] = 55.; x[7] = 65.; x[8] = 75.; x[9] = 85.;
    ex = new Double_t[10];
    ex[0] = 1.; ex[1] = 5.; ex[2] = 5.; ex[3] = 5.; ex[4] = 5.; ex[5] = 5.; ex[6] = 5.; ex[7] = 5.; ex[8] = 5.; ex[9] = 5.;
    bin = new Double_t[10];
    bin[0] = 0.; bin[1] = 10.; bin[2] = 20.; bin[3] = 30.; bin[4] = 40.; bin[5] = 50.; bin[6] = 60.; bin[7] = 70.; bin[8] = 80.; bin[9] = 90.;
    label = new TString[11];
    label[0] = "0-90%"; label[1] = "0-10%"; label[2] = "10-20%"; label[3] = "20-30%"; label[4] = "30-40%"; label[5] = "40-50%"; label[6] = "50-60%"; label[7] = "60-70%"; label[8] = "70-80%"; label[9] = "80-90%"; label[10] = "90-100%";
	
    // plot versus NPart
    NPart = new Double_t[10]; // weighted using Alberica's macro
    NPart[0] = 124.3; NPart[1] = 356.5; NPart[2] = 260.5; NPart[3] = 186.4; NPart[4] = 128.9; NPart[5] = 85.0; NPart[6] = 52.8; NPart[7] = 30.0; NPart[8] = 15.8; NPart[9] = 7.52;
    eNPart = new Double_t[10];
    eNPart[0] = 2.3; eNPart[1] = 3.6; eNPart[2] = 4.4; eNPart[3] = 3.9; eNPart[4] = 3.3; eNPart[5] = 2.6; eNPart[6] = 2.0; eNPart[7] = 1.3; eNPart[8] = 0.6; eNPart[9] = 0.4;
    
  } else {
    
    ::Error("SetGraphLabelsVsCent", "invalid number of centrality bins");
    exit(0);
    
  }
  
}

//------------------------------------------------------------------------
void SetGraphLabelsVsPtY(Int_t nBins, Double_t *&x, Double_t *&ex)
{
  // define graph x-axis labels according to the number of centrality bins
  
  if (nBins == 6) {
    
    // plot versus y bin
    x = new Double_t[6];
    x[0] = 2.625; x[1] = 2.875; x[2] = 3.125; x[3] = 3.375; x[4] = 3.625; x[5] = 3.875;
    ex = new Double_t[6];
    ex[0] = 0.125; ex[1] = 0.125; ex[2] = 0.125; ex[3] = 0.125; ex[4] = 0.125; ex[5] = 0.125;
  } 
  else if (nBins == 7) {
    
    // plot versus pt bin
    x = new Double_t[7];
    x[0] = 0.5; x[1] = 1.5; x[2] = 2.5; x[3] = 3.5; x[4] = 4.5; x[5] = 5.5; x[6] = 7.;
    ex = new Double_t[7];
    ex[0] = 0.5; ex[1] = 0.5; ex[2] = 0.5; ex[3] = 0.5; ex[4] = 0.5; ex[5] = 0.5; ex[6] = 1.;
  }
  else if (nBins == 3) {
    
    // plot versus pt bin
    x = new Double_t[3];
    x[0] = 1.; x[1] = 3.; x[2] = 6.;
    ex = new Double_t[3];
    ex[0] = 1.; ex[1] = 1.; ex[2] = 2.;
  } 
  else {
    
    ::Error("SetGraphLabelsVsPtY", "invalid number of centrality bins");
    exit(0);
    
  }
  
}





//------------------------------------------------------------------------
void SetNormVar(Int_t pt, Int_t y, Double_t &sigmaJPsipp, Double_t &esigmaJPsippStat, Double_t &esigmaJPsippSyst)
{
  // set the normalization variables according to rapidity or pt bin
  
  Double_t ppStatErr = 0.;
  Double_t ppSysErr = 0.;
  
  //(syst. uncertainty is evaluated from all the contributions, including Luminosity. Polarization is not included.):
  //luminosity is 3%.
  
  if (y == 0 && pt == 0) { // 0<pt<8 ,  -2.5<y<-4

    sigmaJPsipp = 3.343;
    ppStatErr = 0.132;
    ppSysErr = 0.241;  // 0.281
    
  } else if (y == 1 && pt == 0) { // 0<pt<8 , -3.<y<-2.5
    
    sigmaJPsipp = 1.355;
    ppStatErr = 0.135;
    ppSysErr = 0.075;
	// -3.25<y<-2.5
	//sigmaJPsipp = 1.922;
    //ppStatErr = 0.108;
    //ppSysErr = 0.162;
	
  } else if (y == 2 && pt == 0) { // 0<pt<8 , -4<y<-3.5 

	sigmaJPsipp = 0.92;
    ppStatErr = 0.0875;
    ppSysErr = 0.05;
	// -4<y<-3.25
    //sigmaJPsipp = 1.423;
    //ppStatErr = 0.071;
    //ppSysErr = 0.12;
    
  } else if (y == 0 && pt == 1) { // 0<pt<2 ,  -2.5<y<-4
    
    sigmaJPsipp = 1.6275;
    ppStatErr = 0.080;
    ppSysErr = 0.068;
    
  } else if (y == 0 && pt == 2) { // 2<pt<5 ,  -2.5<y<-4
    
    sigmaJPsipp = 1.559;
    ppStatErr = 0.074;
    ppSysErr = 0.058;
    
  } else if (y == 0 && pt == 3) { // 5<pt<8 ,  -2.5<y<-4
    
    sigmaJPsipp = 0.167;
    ppStatErr = 0.020;
    ppSysErr =  0.007;
    
  } 
  
  // relative error**2 on p-p cross-section (without uncertainty from luminosity)
  Double_t esigmaJPsippSyst2 = (ppSysErr*ppSysErr) / (sigmaJPsipp*sigmaJPsipp);
  Double_t esigmaJPsippStat2 = (ppStatErr*ppStatErr) / (sigmaJPsipp*sigmaJPsipp);

  // remove relative error on the branching ratio
  esigmaJPsippSyst2 -= eBR*eBR;
  
  // relative error on p-p cross-section
  esigmaJPsippSyst = TMath::Sqrt(esigmaJPsippSyst2);
  esigmaJPsippStat = TMath::Sqrt(esigmaJPsippStat2);

  
}

//------------------------------------------------------------------------
void SetNormVarPtY(Int_t pt, Int_t y, Double_t *&sigmaJPsipp, Double_t *&esigmaJPsippStat, Double_t *&esigmaJPsippSyst)
{
  // set the normalization variables according to rapidity or pt bin
  
  //(syst. uncertainty is evaluated from all the contributions, including Luminosity. Polarization is not included.):
  //luminosity is 3%.
  
  if (y == 0 && pt > 3 && pt <8) { // pt bins ,  -2.5<y<-4
    
	sigmaJPsipp = new Double_t[7];
	Double_t *ppStatErr = new Double_t[7];
	Double_t *ppSysErr = new Double_t[7];
	
    sigmaJPsipp[0] = 1.5 * 0.380 ; sigmaJPsipp[1] = 1.5 * 0.705 ; sigmaJPsipp[2] = 1.5 * 0.583 ; sigmaJPsipp[3] = 1.5 * 0.321 ; sigmaJPsipp[4] = 1.5 * 0.135 ; sigmaJPsipp[5] = 1.5 * 0.073 ; sigmaJPsipp[6] = 1.5 * 2 * 0.019 ; 
    ppStatErr[0] = 1.5 * 0.033 ; ppStatErr[1] = 1.5 * 0.042 ; ppStatErr[2] = 1.5 * 0.038 ; ppStatErr[3] = 1.5 * 0.027 ; ppStatErr[4] = 1.5 * 0.017 ; ppStatErr[5] = 1.5 * 0.011 ; ppStatErr[6] = 1.5 * 2 * 0.004 ; 
    ppSysErr[0] = 1.5 * 0.021 ; ppSysErr[1] = 1.5 * 0.040 ; ppSysErr[2] = 1.5 * 0.033 ; ppSysErr[3] = 1.5 * 0.018 ; ppSysErr[4] = 1.5 * 0.008 ; ppSysErr[5] = 1.5 * 0.004 ; ppSysErr[6] = 1.5 * 2 * 0.001 ; 
	
	Double_t *esigmaJPsippSyst2 = new Double_t[7];
	Double_t *esigmaJPsippStat2 = new Double_t[7];
	esigmaJPsippSyst = new Double_t[7];
	esigmaJPsippStat = new Double_t[7];
	
	for (Int_t i=0; i<7; i++)
	{
	  // relative error**2 on p-p cross-section (without uncertainty from luminosity)
	  esigmaJPsippSyst2[i] = (ppSysErr[i]*ppSysErr[i]) / (sigmaJPsipp[i]*sigmaJPsipp[i]);
	  esigmaJPsippStat2[i] = (ppStatErr[i]*ppStatErr[i]) / (sigmaJPsipp[i]*sigmaJPsipp[i]);

	  // remove relative error on the branching ratio
	  esigmaJPsippSyst2[i] -= eBR*eBR;
	  
	  // relative error on p-p cross-section
	  esigmaJPsippSyst[i] = TMath::Sqrt(esigmaJPsippSyst2[i]);
	  esigmaJPsippStat[i] = TMath::Sqrt(esigmaJPsippStat2[i]);

	}
	delete esigmaJPsippSyst2;
	delete esigmaJPsippStat2;
	delete ppStatErr;
	delete ppSysErr;
    
  } else if (y == 0 && pt == 8) { // pt bins cent 60-90% ,  -2.5<y<-4
  
	sigmaJPsipp = new Double_t[3];
	Double_t *ppStatErr = new Double_t[3];
	Double_t *ppSysErr = new Double_t[3];
	
    sigmaJPsipp[0] = 1.5 * (0.380 + 0.705); sigmaJPsipp[1] = 1.5 * (0.583 + 0.321) ; sigmaJPsipp[2] = 1.5 * (0.135 + 0.073 + (2*0.019)) ;
    ppStatErr[0] = 0.080 ; ppStatErr[1] = 0.070 ; ppStatErr[2] = 0.032 ;
    ppSysErr[0] = 0.068 ; ppSysErr[1] = 0.056 ; ppSysErr[2] = 0.014 ;
	
	Double_t *esigmaJPsippSyst2 = new Double_t[3];
	Double_t *esigmaJPsippStat2 = new Double_t[3];
	esigmaJPsippSyst = new Double_t[3];
	esigmaJPsippStat = new Double_t[3];
	
	for (Int_t i=0; i<3; i++)
	{
	  // relative error**2 on p-p cross-section (without uncertainty from luminosity)
	  esigmaJPsippSyst2[i] = (ppSysErr[i]*ppSysErr[i]) / (sigmaJPsipp[i]*sigmaJPsipp[i]);
	  esigmaJPsippStat2[i] = (ppStatErr[i]*ppStatErr[i]) / (sigmaJPsipp[i]*sigmaJPsipp[i]);
	  
	  // remove relative error on the branching ratio
	  esigmaJPsippSyst2[i] -= eBR*eBR;
	  
	  // relative error on p-p cross-section
	  esigmaJPsippSyst[i] = TMath::Sqrt(esigmaJPsippSyst2[i]);
	  esigmaJPsippStat[i] = TMath::Sqrt(esigmaJPsippStat2[i]);
	  
	}
	delete esigmaJPsippSyst2;
	delete esigmaJPsippStat2;
	delete ppStatErr;
	delete ppSysErr;
	
  } else if (y == 3 && pt == 0) { // 0<pt<8 ,  y bins
    
	sigmaJPsipp = new Double_t[6];
	Double_t *ppStatErr = new Double_t[6];
	Double_t *ppSysErr = new Double_t[6];
	
    sigmaJPsipp[0] = 0.25 * 3.05 ; sigmaJPsipp[1] = 0.25 * 2.37 ; sigmaJPsipp[2] = 0.25 * 2.26 ; sigmaJPsipp[3] = 0.25 * 2.01 ; sigmaJPsipp[4] = 0.25 * 2.00 ; sigmaJPsipp[5] = 0.25 * 1.68 ;
    ppStatErr[0] = 0.25 * 0.35 ; ppStatErr[1] = 0.25 * 0.19 ; ppStatErr[2] = 0.25 * 0.15 ; ppStatErr[3] = 0.25 * 0.14 ; ppStatErr[4] = 0.25 * 0.16 ; ppStatErr[5] = 0.25 * 0.19 ; 
    ppSysErr[0] = 0.25 * 0.17 ; ppSysErr[1] = 0.25 * 0.13 ; ppSysErr[2] = 0.25 * 0.13 ; ppSysErr[3] = 0.25 * 0.11 ; ppSysErr[4] = 0.25 * 0.11 ; ppSysErr[5] = 0.25 * 0.09 ; 
    
	Double_t *esigmaJPsippSyst2 = new Double_t[6];
	Double_t *esigmaJPsippStat2 = new Double_t[6];
	esigmaJPsippSyst = new Double_t[6];
	esigmaJPsippStat = new Double_t[6];
	
	for (Int_t i=0; i<6; i++)
	{
	  // relative error**2 on p-p cross-section (without uncertainty from luminosity)
	  esigmaJPsippSyst2[i] = (ppSysErr[i]*ppSysErr[i]) / (sigmaJPsipp[i]*sigmaJPsipp[i]);
	  esigmaJPsippStat2[i] = (ppStatErr[i]*ppStatErr[i]) / (sigmaJPsipp[i]*sigmaJPsipp[i]);
	  
	  // remove relative error on the branching ratio
	  esigmaJPsippSyst2[i] -= eBR*eBR;
	  
	  // relative error on p-p cross-section
	  esigmaJPsippSyst[i] = TMath::Sqrt(esigmaJPsippSyst2[i]);
	  esigmaJPsippStat[i] = TMath::Sqrt(esigmaJPsippStat2[i]);
	}
	delete esigmaJPsippSyst2;
	delete esigmaJPsippStat2;
	delete ppStatErr;
	delete ppSysErr;
  }

}


//------------------------------------------------------------------------
void SetNormVarVsCent(Int_t pt, Int_t y, Int_t nBins, Double_t *&nMB, Double_t *&TAA, Double_t *&eTAA, Double_t *&AccEff, Double_t *&eAccEff)
{
  // set the normalization variables according to the number of centrality bins
  
  if (nBins == 8) {
	 
	nMB = new Double_t[8];    
	//nMB[0] = 480292339.; nMB[1] = 53340210.; nMB[2] = 53366740.; nMB[3] = 53364137.; nMB[4] = 53366597.; nMB[5] = 53364308.; nMB[6] = 53368749.; nMB[7] = 160102906.;   // AOD083
	nMB[0] = 475266322.; nMB[1] = 52807369.; nMB[2] = 52807369.; nMB[3] = 52807369.; nMB[4] = 52807369.; nMB[5] = 52807369.; nMB[6] = 52807369.; nMB[7] = 158422107;  // AOD101

    TAA = new Double_t[8]; // mb^-1
    TAA[0] = 6.26377; TAA[1] = 23.4797; TAA[2] = 14.4318; TAA[3] = 8.73769; TAA[4] = 5.02755; TAA[5] = 2.68327; TAA[6] = 1.32884; TAA[7] = 0.312433;
    eTAA = new Double_t[8];
    eTAA[0] = 0.237682; eTAA[1] = 0.972587; eTAA[2] = 0.573289; eTAA[3] = 0.370219; eTAA[4] = 0.22099; eTAA[5] = 0.137073; eTAA[6] = 0.0929536; eTAA[7] = 0.0271638;    
	
	
	// Acc * Eff
    AccEff = new Double_t[8];
    eAccEff = new Double_t[8];
	
    if (y == 1 && pt == 0) { //  0<pt<8     -3<y<-2.5 
      AccEff[0] = 0.085568 ; AccEff[1] = 0.086017 ; AccEff[2] = 0.086043 ; AccEff[3] = 0.088314 ; AccEff[4] = 0.087452 ; AccEff[5] = 0.089372 ; AccEff[6] = 0.089116 ; AccEff[7] = 0.089040 ;  
      eAccEff[0] = 0.000922;  eAccEff[1] = 0.000937;   eAccEff[2] = 0.000939;  eAccEff[3] = 0.000962;  eAccEff[4] =0.000955;  eAccEff[5] = 0.000964; eAccEff[6] = 0.000949; eAccEff[7] = 0.000554;
	  // -3.25<y<-2.5 
	  //AccEff[0] = 0.122992 ; AccEff[1] = 0.120086 ; AccEff[2] = 0.123262 ; AccEff[3] = 0.125173 ; AccEff[4] = 0.125345 ; AccEff[5] = 0.128141 ; AccEff[6] = 0.126402 ; AccEff[7] = 0.127504 ;  
      //eAccEff[0] = 0.000909;  eAccEff[1] = 0.000906;   eAccEff[2] = 0.000922;  eAccEff[3] = 0.000924;  eAccEff[4] =0.000927;  eAccEff[5] = 0.000938; eAccEff[6] = 0.000931; eAccEff[7] = 0.000538;
	}	
	else if (y == 2 && pt == 0) { //  0<pt<8    -4<y<-3.5 
      AccEff[0] = 0.120188 ;  AccEff[1] = 0.116536;  AccEff[2] = 0.120453;  AccEff[3] = 0.123172;  AccEff[4] = 0.123574;  AccEff[5] = 0.122878; AccEff[6] = 0.128529 ; AccEff[7] = 0.128296 ;  
      eAccEff[0] = 0.001316;  eAccEff[1] = 0.001317; eAccEff[2] = 0.001341; eAccEff[3] = 0.001326;  eAccEff[4] = 0.001350;  eAccEff[5] = 0.001347; eAccEff[6] = 0.001384; eAccEff[7] = 0.000787;
	  // -4<y<-3.25
	  //AccEff[0] =0.155556 ;  AccEff[1] = 0.150021;  AccEff[2] = 0.156111;  AccEff[3] = 0.158418;  AccEff[4] = 0.161031;  AccEff[5] = 0.160710; AccEff[6] = 0.165914 ; AccEff[7] = 0.165924 ;  
      //eAccEff[0] = 0.001159;  eAccEff[1] = 0.001159; eAccEff[2] = 0.001179; eAccEff[3] = 0.001169;  eAccEff[4] = 0.001191;  eAccEff[5] = 0.001190; eAccEff[6] = 0.001216; eAccEff[7] = 0.000695;
	}
	else if (y == 0 && pt == 1) { //  0<pt<2    -4<y<-2.5 
      AccEff[0] = 0.128539 ;  AccEff[1] = 0.124928;  AccEff[2] = 0.128148;  AccEff[3] = 0.132239;  AccEff[4] = 0.132922;  AccEff[5] = 0.132807; AccEff[6] = 0.134235 ; AccEff[7] = 0.135017 ;
      eAccEff[0] = 0.001017;  eAccEff[1] = 0.001090; eAccEff[2] = 0.001027; eAccEff[3] = 0.001030;  eAccEff[4] = 0.001046;  eAccEff[5] = 0.001038; eAccEff[6] = 0.001050; eAccEff[7] = 0.000604;
	} 
	else if (y == 0 && pt == 2) { //  2<pt<5    -4<y<-2.5 
      AccEff[0] = 0.133051;   AccEff[1] = 0.128706;      AccEff[2] = 0.134742;    AccEff[3] = 0.134412;  AccEff[4] = 0.136324;  AccEff[5] = 0.137792; AccEff[6] = 0.139714 ; AccEff[7] = 0.140498;
      eAccEff[0] = 0.001039;  eAccEff[1] = 0.001033;  eAccEff[2] = 0.001060;  eAccEff[3] = 0.001050;  eAccEff[4] = 0.001065;  eAccEff[5] = 0.001065; eAccEff[6] = 0.001076; eAccEff[7] = 0.000621;
	} 
	else if (y == 0 && pt == 3) { //  5<pt<8    -4<y<-2.5 
      AccEff[0] = 0.236292;  AccEff[1] = 0.231518;  AccEff[2] = 0.232606;  AccEff[3] = 0.237468;  AccEff[4] = 0.237403;  AccEff[5] = 0.252162; AccEff[6] = 0.247304 ; AccEff[7] = 0.244922;
      eAccEff[0] = 0.002110;  eAccEff[1] = 0.003680;  eAccEff[2] = 0.003762;  eAccEff[3] = 0.003759;  eAccEff[4] = 0.003691;  eAccEff[5] = 0.003968; eAccEff[6] = 0.003911; eAccEff[7] = 0.002198;
	}
	
	
  } else if (nBins == 10) {
    
    nMB = new Double_t[10];
    //nMB[0] = 480292339.; nMB[1] = 53340210.; nMB[2] = 53366740.; nMB[3] = 53364137.; nMB[4] = 53366597.; nMB[5] = 53364308.; nMB[6] = 53368749.; nMB[7] = 53365378.; nMB[8] = 53364434.; nMB[9] = 53373094.;  // AOD083
   	nMB[0] = 475266322.; nMB[1] = 52807369.; nMB[2] = 52807369.; nMB[3] = 52807369.; nMB[4] = 52807369.; nMB[5] = 52807369.; nMB[6] = 52807369.; nMB[7] = 52807369.; nMB[8] = 52807369.; nMB[9] = 52807369.;   // AOD101
	
    TAA = new Double_t[10]; // mb^-1
    TAA[0] = 6.26377; TAA[1] = 23.4797; TAA[2] = 14.4318; TAA[3] = 8.73769; TAA[4] = 5.02755; TAA[5] = 2.68327; TAA[6] = 1.32884; TAA[7] = 0.590485; TAA[8] = 0.239622; TAA[9] = 0.0973951;
    eTAA = new Double_t[10];
    eTAA[0] = 0.237682; eTAA[1] = 0.972587; eTAA[2] = 0.573289; eTAA[3] = 0.370219; eTAA[4] = 0.22099; eTAA[5] = 0.137073; eTAA[6] = 0.0929536; eTAA[7] = 0.0449641; eTAA[8] = 0.0259041; eTAA[9] = 0.00988253;
    
	
    AccEff = new Double_t[10];
    eAccEff = new Double_t[10];
    
	if (y == 0 && pt == 0) { //  0<pt<8    -4<y<-2.5 
      AccEff[0] = 0.136726 ; AccEff[1] = 0.132838 ; AccEff[2] = 0.137253 ; AccEff[3] = 0.139356 ; AccEff[4] = 0.140610 ; AccEff[5] = 0.142014 ; AccEff[6] = 0.143283 ; AccEff[7] = 0.143082 ; AccEff[8] = 0.144013 ; AccEff[9] = 0.144583 ;  
      eAccEff[0] = 0.000727;  eAccEff[1] = 0.000718;   eAccEff[2] = 0.000730;  eAccEff[3] = 0.000728;  eAccEff[4] =0.000737;  eAccEff[5] = 0.000740; eAccEff[6] = 0.000746; eAccEff[7] = 0.000740; eAccEff[8] = 0.000740; eAccEff[9] = 0.000746;
    }
	/*
	else if (y == 0 && pt == 1) { //  0<pt<2    -4<y<-3.25 
      AccEff[0] = 0.128539 ;  AccEff[1] = 0.124928;  AccEff[2] = 0.128148;  AccEff[3] = 0.132239;  AccEff[4] = 0.132922;  AccEff[5] = 0.132807; AccEff[6] = 0.134235 ; AccEff[7] = 0.135454 ; AccEff[8] = 0.135259 ; AccEff[9] = 0.134339 ;
      eAccEff[0] = 0.001017;  eAccEff[1] = 0.001090; eAccEff[2] = 0.001027; eAccEff[3] = 0.001030;  eAccEff[4] = 0.001046;  eAccEff[5] = 0.001038; eAccEff[6] = 0.001050; eAccEff[7] = 0.001124; eAccEff[8] = 0.001118; eAccEff[9] = 0.001130;
	} 
	else if (y == 0 && pt == 2) { //  2<pt<5    -4<y<-3.25 
      AccEff[0] = 0.133051;   AccEff[1] = 0.128706;      AccEff[2] = 0.134742;    AccEff[3] = 0.134412;  AccEff[4] = 0.136324;  AccEff[5] = 0.137792; AccEff[6] = 0.139714 ; AccEff[7] = 0.138584 ; AccEff[8] = 0.140751 ; AccEff[9] = 0.142144 ;
      eAccEff[0] = 0.001039;  eAccEff[1] = 0.001033;  eAccEff[2] = 0.001060;  eAccEff[3] = 0.001050;  eAccEff[4] = 0.001065;  eAccEff[5] = 0.001065; eAccEff[6] = 0.001076; eAccEff[7] = 0.001158; eAccEff[8] = 0.001156; eAccEff[9] = 0.001166;
	} 
	else if (y == 0 && pt == 3) { //  5<pt<8    -4<y<-3.25 
      AccEff[0] = 0.236292;  AccEff[1] = 0.231518;  AccEff[2] = 0.232606;  AccEff[3] = 0.237468;  AccEff[4] = 0.237403;  AccEff[5] = 0.252162; AccEff[6] = 0.247304 ; AccEff[7] = 0.241632 ; AccEff[8] = 0.243625 ; AccEff[9] = 0.249534 ;
      eAccEff[0] = 0.002110;  eAccEff[1] = 0.003680;  eAccEff[2] = 0.003762;  eAccEff[3] = 0.003759;  eAccEff[4] = 0.003691;  eAccEff[5] = 0.003968; eAccEff[6] = 0.003911; eAccEff[7] = 0.004277; eAccEff[8] = 0.004410; eAccEff[9] = 0.004472;
	} 
    */
	
  } else {
    
    ::Error("SetNormVarVsCent", "invalid number of centrality bins");
    exit(0);
    
  }
  
}

//------------------------------------------------------------------------
void SetNormVarVsPtY(Int_t pt, Int_t y, Double_t &nMB, Double_t &TAA, Double_t &eTAA, Double_t *&AccEff, Double_t *&eAccEff)
{
  // set the normalization variables according to the number of centrality bins
  
  if (y == 0 && pt == 4) { // pt bins  -4<y<-2.5   Cent 0-90%
    
    //nMB = 480292339.; // AOD083
	nMB = 475266322.; // AOD101
	
    TAA = 6.26377; // mb^-1
    eTAA = 0.237682;  
	
    AccEff = new Double_t[7];
    eAccEff = new Double_t[7];
	
	AccEff[0] = 0.143316 ; AccEff[1] = 0.121247 ; AccEff[2] = 0.119741 ; AccEff[3] = 0.136591 ; AccEff[4] = 0.172319 ; AccEff[5] = 0.216304 ; AccEff[6] = 0.260633 ;  
	eAccEff[0] = 0.000930;   eAccEff[1] = 0.000621;  eAccEff[2] = 0.000688;  eAccEff[3] =0.000955;  eAccEff[4] = 0.001506; eAccEff[5] = 0.002386; eAccEff[6] = 0.002992;
  }
  else if (y == 0 && pt == 5) { // pt bins  -4<y<-2.5   Cent 0-20%
    
    //nMB = 106706950.; // AOD083
	nMB = 105614738.; // AOD101
	
    TAA = 18.9262; // mb^-1
    eTAA = 0.741423;  
	
    AccEff = new Double_t[7];
    eAccEff = new Double_t[7];
	
	AccEff[0] = 0.140862 ; AccEff[1] = 0.118835 ; AccEff[2] = 0.118610 ; AccEff[3] = 0.134156 ; AccEff[4] = 0.168886 ; AccEff[5] = 0.213416 ; AccEff[6] = 0.257198 ;  
	eAccEff[0] = 0.001336;   eAccEff[1] = 0.000894;  eAccEff[2] = 0.000991;  eAccEff[3] = 0.001375;  eAccEff[4] = 0.002158; eAccEff[5] = 0.003418; eAccEff[6] = 0.004330;
  }  
  else if (y == 0 && pt == 6) { // pt bins  -4<y<-2.5   Cent 20-40%
    
    //nMB = 106730734.; // AOD083
	nMB = 105614738.; // AOD101
	
    TAA = 6.85556; // mb^-1
    eTAA = 0.283436;  
	
    AccEff = new Double_t[7];
    eAccEff = new Double_t[7];
	
	AccEff[0] = 0.147138 ; AccEff[1] = 0.125221 ; AccEff[2] = 0.120628 ; AccEff[3] = 0.139761 ; AccEff[4] = 0.177449 ; AccEff[5] = 0.218353 ; AccEff[6] = 0.263094 ;  
	eAccEff[0] = 0.001369;   eAccEff[1] = 0.000900;  eAccEff[2] = 0.000992;  eAccEff[3] = 0.001403;  eAccEff[4] = 0.002229; eAccEff[5] = 0.003549; eAccEff[6] = 0.004248;
  }
  else if (y == 0 && pt == 7) { // pt bins  -4<y<-2.5   Cent 40-90%
    
    //nMB = 266835963.; // AOD083
	nMB = 264036845.; // AOD101
	
    TAA = 1.033; // mb^-1
    eTAA = 0.05017;      
	
    AccEff = new Double_t[7];
    eAccEff = new Double_t[7];
	
	AccEff[0] = 0.148891 ; AccEff[1] = 0.126241 ; AccEff[2] = 0.124492 ; AccEff[3] = 0.143573 ; AccEff[4] = 0.180765 ; AccEff[5] = 0.228988 ; AccEff[6] = 0.275603 ;  
	eAccEff[0] = 0.001095;   eAccEff[1] = 0.000718;  eAccEff[2] = 0.000802;  eAccEff[3] = 0.001129;  eAccEff[4] = 0.001792; eAccEff[5] = 0.002938; eAccEff[6] = 0.003590;
  }
  else if (y == 0 && pt == 8) { // pt bins  -4<y<-2.5   Cent 60-90%
    
	nMB = 158422107.; // AOD101
	
    TAA = 0.312433; // mb^-1
    eTAA = 0.0271638; 
	
    AccEff = new Double_t[3];
    eAccEff = new Double_t[3];
	
	AccEff[0] = 0.135245 ; AccEff[1] = 0.132686 ; AccEff[2] = 0.208055 ;
	eAccEff[0] = 0.000694;   eAccEff[1] = 0.000758;  eAccEff[2] = 0.001587;
  }
  else if (y == 3 && pt == 0) { // y bins  0<pt<8   Cent 0-90%
    
    //nMB = 480292339.; // AOD083
	nMB = 475266322.; // AOD101
	
    TAA = 6.26377; // mb^-1
    eTAA = 0.237682;  
	
    AccEff = new Double_t[6];
    eAccEff = new Double_t[6];
	
	AccEff[0] = 0.038823 ; AccEff[1] = 0.138947 ; AccEff[2] = 0.205179 ; AccEff[3] = 0.215015 ; AccEff[4] = 0.169379 ; AccEff[5] = 0.064762 ; 
	eAccEff[0] = 0.000320;   eAccEff[1] = 0.000596;  eAccEff[2] = 0.000725;  eAccEff[3] = 0.000774;  eAccEff[4] = 0.000747; eAccEff[5] = 0.000528;
  }
  else {
    
    ::Error("SetNormVarVsPtY", "invalid number of pt or y");
    exit(0);
    
  }
  
  
}


//------------------------------------------------------------------------
void SetSystVarVsCent(Int_t nBins, Double_t *&eGen, Double_t *&eRec, Double_t *&eTrk, Double_t *&eTrg)
{
  // set the sytematic errors (in percent) according to the number of centrality bins
  
  if (nBins == 8) { // 0-90, 0-10, 10-20, 20-30, 30-40, 40-50, 50-60, 60-90, all
    
    eGen = new Double_t[9];
    eGen[0] = 0.; eGen[1] = 0.; eGen[2] = 0.; eGen[3] = 0.; eGen[4] = 0.; eGen[5] = 0.; eGen[6] = 0.; eGen[7] = 0.; eGen[8] = 0.05;
    eRec = new Double_t[9];
    eRec[0] = 0.; eRec[1] = 0.; eRec[2] = 0.; eRec[3] = 0.; eRec[4] = 0.; eRec[5] = 0.; eRec[6] = 0.; eRec[7] = 0.; eRec[8] = 0.02;
    eTrk = new Double_t[9];
    eTrk[0] = 0.005; eTrk[1] = 0.01; eTrk[2] = 0.005; eTrk[3] = 0.; eTrk[4] = 0.; eTrk[5] = 0.; eTrk[6] = 0.; eTrk[7] = 0.; eTrk[8] = 0.06;
    eTrg = new Double_t[9];
    eTrg[0] = 0.015; eTrg[1] = 0.02; eTrg[2] = 0.015; eTrg[3] = 0.01; eTrg[4] = 0.005; eTrg[5] = 0.; eTrg[6] = 0.; eTrg[7] = 0.; eTrg[8] = 0.064;
    
  } else if (nBins == 10) { // 0-90, 0-10, 10-20, 20-30, 30-40, 40-50, 50-60, 60-70, 70-80, 80-90, all
    
    eGen = new Double_t[11];
    eGen[0] = 0.; eGen[1] = 0.; eGen[2] = 0.; eGen[3] = 0.; eGen[4] = 0.; eGen[5] = 0.; eGen[6] = 0.; eGen[7] = 0.; eGen[8] = 0.; eGen[9] = 0.; eGen[10] = 0.05;
    eRec = new Double_t[11];
    eRec[0] = 0.; eRec[1] = 0.; eRec[2] = 0.; eRec[3] = 0.; eRec[4] = 0.; eRec[5] = 0.; eRec[6] = 0.; eRec[7] = 0.; eRec[8] = 0.; eRec[9] = 0.; eRec[10] = 0.02;
    eTrk = new Double_t[11];
    eTrk[0] = 0.005; eTrk[1] = 0.01; eTrk[2] = 0.005; eTrk[3] = 0.; eTrk[4] = 0.; eTrk[5] = 0.; eTrk[6] = 0.; eTrk[7] = 0.; eTrk[8] = 0.; eTrk[9] = 0.; eTrk[10] = 0.06;
    eTrg = new Double_t[11];
    eTrg[0] = 0.015; eTrg[1] = 0.02; eTrg[2] = 0.015; eTrg[3] = 0.01; eTrg[4] = 0.005; eTrg[5] = 0.; eTrg[6] = 0.; eTrg[7] = 0.; eTrg[8] = 0.; eTrg[9] = 0.; eTrg[10] = 0.064;
    
  } else {
    
    ::Error("SetSystVarVsCent", "invalid number of centrality bins");
    exit(0);
    
  }
  
}

//------------------------------------------------------------------------
void SetSystVarVsPtY(Int_t pt, Int_t y, Double_t *&eGen, Double_t *&eRec, Double_t *&eTrk, Double_t *&eTrg)
{
  // set the sytematic errors (in percent) according to the number of centrality bins
  
  if ( (y == 0 && pt == 4) || (y == 3 && pt == 0) ) { // 0-90, all
    
    eGen = new Double_t[2];
    eGen[0] = 0.; eGen[1] = 0.05;
    eRec = new Double_t[2];
    eRec[0] = 0.; eRec[1] = 0.02;
    eTrk = new Double_t[2];
    eTrk[0] = 0.005; eTrk[1] = 0.06;
    eTrg = new Double_t[2];
    eTrg[0] = 0.015; eTrg[1] = 0.064;
  } 
  else if (y == 0 && pt == 5) { // 0-20, all
    
    eGen = new Double_t[2];
    eGen[0] = 0.; eGen[1] = 0.05;
    eRec = new Double_t[2];
    eRec[0] = 0.; eRec[1] = 0.02;
    eTrk = new Double_t[2];
    eTrk[0] = 0.01; eTrk[1] = 0.06;
    eTrg = new Double_t[2];
    eTrg[0] = 0.02; eTrg[1] = 0.064;
  }
  else if (y == 0 && pt == 6) { // 20-40, all
    
    eGen = new Double_t[2];
    eGen[0] = 0.; eGen[1] = 0.05;
    eRec = new Double_t[2];
    eRec[0] = 0.; eRec[1] = 0.02;
    eTrk = new Double_t[2];
    eTrk[0] = 0.; eTrk[1] = 0.06;
    eTrg = new Double_t[2];
    eTrg[0] = 0.01; eTrg[1] = 0.064;
  }
  else if (y == 0 && pt > 6) { // 40-90 or 60-90, all
    
    eGen = new Double_t[2];
    eGen[0] = 0.; eGen[1] = 0.05;
    eRec = new Double_t[2];
    eRec[0] = 0.; eRec[1] = 0.02;
    eTrk = new Double_t[2];
    eTrk[0] = 0.; eTrk[1] = 0.06;
    eTrg = new Double_t[2];
    eTrg[0] = 0.; eTrg[1] = 0.064;
  }
  else {
    
    ::Error("SetSystVarVsPtY", "invalid number of pt or y");
    exit(0);
    
  }
  
}



//------------------------------------------------------------------------
void ExtractNJPsiVsCent(Int_t pt, Int_t y, Int_t nBins, Int_t nT, TString SigFct, Double_t **&n, Double_t **&en)
{
  
  TString centBinName[10] = {"090", "010", "1020", "2030", "3040", "4050", "5060", "6090", "7080", "8090"};
  
  
  TString Study;
  if (y == 0 && pt == 0) {Study="pt_0"; centBinName[7]="6070";}
  else if (y == 1 && pt == 0) Study="y_1";
  else if (y == 2 && pt == 0) Study="y_3";
  else if (y == 0 && pt == 1) Study="pt_1";
  else if (y == 0 && pt == 2) Study="pt_2";
  else if (y == 0 && pt == 3) Study="pt_3";
  
  TString dir = "../../../SignalExtraction/SL";
  Int_t nTests=0;
  
  n = new Double_t*[nT];
  en = new Double_t*[nT];
  for (Int_t i=0; i<nT; i++) {
    n[i] = new Double_t[nBins];
    en[i] = new Double_t[nBins];
  }
  
  // Different tails
  TList* Tails = new TList();
  Tails->Add(new TObjString("embedding2011"));
  Tails->Add(new TObjString("simuJpsi2011"));
  //########################################
  // loop over cases of tails
  //########################################
  TIter nextTails(Tails);
  TObjString* nameTails;
  while (( nameTails = static_cast<TObjString*>(nextTails()) ))
  {
	  
	// Different Msigma
	TList* Msigma = new TList();
	Msigma->Add(new TObjString("0-90"));
	Msigma->Add(new TObjString("psigma"));
	Msigma->Add(new TObjString("msigma"));
	//########################################
	// loop over cases of Mass and sigma
	//########################################
	TIter nextMsigma(Msigma);
	TObjString* nameMsigma;
	while (( nameMsigma = static_cast<TObjString*>(nextMsigma()) ))
	{
	  
	  //########################################
	  // loop over centrality bin
	  //########################################
	  for (Int_t i=0; i<nBins; i++)
	  {
		
		// init
		Double_t nbJpsiRaw=0., nbJpsiMix=0., errnbJpsiRaw=0., errnbJpsiMix=0.;
		
		//###############  Raw
		// open file		
		ifstream inFileRaw( Form( "%s/FitResults_%s_cent_%s.txt", dir.Data(), Study.Data(), centBinName[i].Data() ) );
		TString currTestRaw;
		if (inFileRaw.is_open())
		{
		  
		  while (! inFileRaw.eof() )
		  {
			currTestRaw.ReadLine(inFileRaw,kTRUE);
			if(currTestRaw.IsNull()) continue;
			
			TObjArray* arr = currTestRaw.Tokenize(" ");
			
			// Select Test
			if (arr->At(0)->GetName() != SigFct) continue;
			else if (arr->At(1)->GetName() != nameTails->String()) continue;
			else if (arr->At(3)->GetName() != nameMsigma->String()) continue;
			
			// Get values
			TString strnbJpsi = arr->At(5)->GetName();
			TString strerrnbJpsi = arr->At(7)->GetName();
			nbJpsiRaw = strnbJpsi.Atof();
			errnbJpsiRaw = strerrnbJpsi.Atof();
			
		  }  // loop in file
		  
		  inFileRaw.close();
		  
		} // file is open
		else printf("Input file unknow !!!!!! \n");
		
		
		//###############  Mix
		// open file
		ifstream inFileMix( Form( "%s/FitResults_Mixing_%s_cent_%s.txt", dir.Data(), Study.Data(), centBinName[i].Data() ) );
		TString currTestMix;
		if (inFileMix.is_open())
		{
		  
		  while (! inFileMix.eof() )
		  {
			currTestMix.ReadLine(inFileMix,kTRUE);
			if(currTestMix.IsNull()) continue;
			
			TObjArray* arr = currTestMix.Tokenize(" ");
			
			// Select Test
			if (arr->At(0)->GetName() != SigFct) continue;
			else if (arr->At(1)->GetName() != nameTails->String()) continue;
			else if (arr->At(3)->GetName() != nameMsigma->String()) continue;
			
			// Get values
			TString strnbJpsi = arr->At(5)->GetName();
			TString strerrnbJpsi = arr->At(7)->GetName();
			nbJpsiMix = strnbJpsi.Atof();
			errnbJpsiMix = strerrnbJpsi.Atof();
			
		  }  // loop in file
		  
		  inFileMix.close();
		}
		  
		n[nTests][i] = nbJpsiRaw;
		en[nTests][i] = errnbJpsiRaw;
		n[nTests+6][i] = nbJpsiMix;
		en[nTests+6][i] = errnbJpsiMix;
		
	  } // loop on centrality
	  
	  nTests++;
	  
	} // end loop on Mas and sigma
	
  } // end loop on Tails
  
}

//------------------------------------------------------------------------
void ExtractNJPsiVsPtY(Int_t pt, Int_t y, Int_t nBins, Int_t nT, TString SigFct, TString Cent, Double_t **&n, Double_t **&en)
{
    
  TString Study="pt";
  Int_t shift=4;
  if (y == 3 && pt == 0) Study="y";
  //else if (y == 0 && pt == 4) Study="Pt_Bins_4y25_cent090";
  //else if (y == 0 && pt == 5) Study="Pt_Bins_4y25_cent020";
  //else if (y == 0 && pt == 6) Study="Pt_Bins_4y25_cent2040";
  //else if (y == 0 && pt == 7) Study="Pt_Bins_4y25_cent4090";
  //else if (y == 0 && pt == 8) Study="Pt_Bins_4y25_cent6090";
  
  TString dir = "../../../SignalExtraction/SL";
  Int_t nTests=0;
  
  n = new Double_t*[nT];
  en = new Double_t*[nT];
  for (Int_t i=0; i<nT; i++) {
    n[i] = new Double_t[nBins];
    en[i] = new Double_t[nBins];
  }
  
  // Different tails
  TList* Tails = new TList();
  Tails->Add(new TObjString("embedding2011"));
  Tails->Add(new TObjString("simuJpsi2011"));
  //########################################
  // loop over cases of tails
  //########################################
  TIter nextTails(Tails);
  TObjString* nameTails;
  while (( nameTails = static_cast<TObjString*>(nextTails()) ))
  {
	
	// Different Msigma
	TList* Msigma = new TList();
	Msigma->Add(new TObjString("0-90"));
	Msigma->Add(new TObjString("psigma"));
	Msigma->Add(new TObjString("msigma"));
	//########################################
	// loop over cases of Mass and sigma
	//########################################
	TIter nextMsigma(Msigma);
	TObjString* nameMsigma;
	while (( nameMsigma = static_cast<TObjString*>(nextMsigma()) ))
	{
	  
	  //########################################
	  // loop over centrality bin
	  //########################################
	  for (Int_t i=0; i<nBins; i++)
	  {
		
		// init
		Double_t nbJpsiRaw=0., nbJpsiMix=0., errnbJpsiRaw=0., errnbJpsiMix=0.;
		
		//###############  Raw
		// open file		
		ifstream inFileRaw( Form( "%s/FitResults_%s_%i_cent_%s.txt", dir.Data(), Study.Data(), i+shift, Cent.Data()) );
		TString currTestRaw;
		if (inFileRaw.is_open())
		{
		  
		  while (! inFileRaw.eof() )
		  {
			currTestRaw.ReadLine(inFileRaw,kTRUE);
			if(currTestRaw.IsNull()) continue;
			
			TObjArray* arr = currTestRaw.Tokenize(" ");
			
			// Select Test
			if (arr->At(0)->GetName() != SigFct) continue;
			else if (arr->At(1)->GetName() != nameTails->String()) continue;
			else if (arr->At(3)->GetName() != nameMsigma->String()) continue;
			
			// Get values
			TString strnbJpsi = arr->At(5)->GetName();
			TString strerrnbJpsi = arr->At(7)->GetName();
			nbJpsiRaw = strnbJpsi.Atof();
			errnbJpsiRaw = strerrnbJpsi.Atof();
			
		  }  // loop in file
		  
		  inFileRaw.close();
		  
		} // file is open
		else printf("Input file unknow !!!!!! \n");
		
		
		//###############  Mix
		// open file
		ifstream inFileMix( Form( "%s/FitResults_Mixing_%s_%i_cent_%s.txt", dir.Data(), Study.Data(), i+shift, Cent.Data()) );
		TString currTestMix;
		if (inFileMix.is_open())
		{
		  
		  while (! inFileMix.eof() )
		  {
			currTestMix.ReadLine(inFileMix,kTRUE);
			if(currTestMix.IsNull()) continue;
			
			TObjArray* arr = currTestMix.Tokenize(" ");
			
			// Select Test
			if (arr->At(0)->GetName() != SigFct) continue;
			else if (arr->At(1)->GetName() != nameTails->String()) continue;
			else if (arr->At(3)->GetName() != nameMsigma->String()) continue;
			
			// Get values
			TString strnbJpsi = arr->At(5)->GetName();
			TString strerrnbJpsi = arr->At(7)->GetName();
			nbJpsiMix = strnbJpsi.Atof();
			errnbJpsiMix = strerrnbJpsi.Atof();
			
		  }  // loop in file
		  
		  inFileMix.close();
		}
		
		n[nTests][i] = nbJpsiRaw;
		en[nTests][i] = errnbJpsiRaw;
		n[nTests+6][i] = nbJpsiMix;
		en[nTests+6][i] = errnbJpsiMix;
		
	  } // loop on centrality
	  
	  nTests++;
	  
	} // end loop on Mas and sigma
	
  } // end loop on Tails
  
}


//------------------------------------------------------------------------
void SetNJPsiVsCent(Double_t **&nJPsi, Double_t **&enJPsi, Int_t &nTests, Int_t &nBins)
{
  // return the number of JPsi and errors for each test and each centrality bins
  // as well as the number of tests and of centrality bins
  
  const Int_t nT = 12;
  const Int_t nC = 10;
  
  ExtractNJPsiVsCent(0, 0, nC, nT, "1CB2", nJPsi, enJPsi);
  
  // return the values
  nTests = nT;
  nBins = nC;
  
}

//------------------------------------------------------------------------
void SetNJPsiVsCentpt1(Double_t **&nJPsi, Double_t **&enJPsi, Int_t &nTests, Int_t &nBins)
{
  // return the number of JPsi with pT<3GeV/c and errors for each test and each centrality bins
  // as well as the number of tests and of centrality bins
  
  const Int_t nT = 12;
  const Int_t nC = 8;
  
  ExtractNJPsiVsCent(1, 0, nC, nT, "1CB2", nJPsi, enJPsi);
  
  // return the values
  nTests = nT;
  nBins = nC;  
}


//------------------------------------------------------------------------
void SetNJPsiVsCentpt2(Double_t **&nJPsi, Double_t **&enJPsi, Int_t &nTests, Int_t &nBins)
{
  // return the number of JPsi with pT>3GeV/c and errors for each test and each centrality bins
  // as well as the number of tests and of centrality bins
  
  const Int_t nT = 12;
  const Int_t nC = 8;
  
  ExtractNJPsiVsCent(2, 0, nC, nT, "1CB2", nJPsi, enJPsi);
  
  // return the values
  nTests = nT;
  nBins = nC;  
}


//------------------------------------------------------------------------
void SetNJPsiVsCentpt3(Double_t **&nJPsi, Double_t **&enJPsi, Int_t &nTests, Int_t &nBins)
{
  // return the number of JPsi with pT>3GeV/c and errors for each test and each centrality bins
  // as well as the number of tests and of centrality bins
  
  const Int_t nT = 12;
  const Int_t nC = 8;
  
  ExtractNJPsiVsCent(3, 0, nC, nT, "1CB2", nJPsi, enJPsi);
  
  // return the values
  nTests = nT;
  nBins = nC;  
}


//------------------------------------------------------------------------
void SetNJPsiVsCenty1(Double_t **&nJPsi, Double_t **&enJPsi, Int_t &nTests, Int_t &nBins)
{
  // return the number of JPsi in -4<y<-3.25 and errors for each test and each centrality bins
  // as well as the number of tests and of centrality bins
  
  const Int_t nT = 12;
  const Int_t nC = 8;
 
  ExtractNJPsiVsCent(0, 1, nC, nT, "1CB2", nJPsi, enJPsi);
  
  // return the values
  nTests = nT;
  nBins = nC;  
}


//------------------------------------------------------------------------
void SetNJPsiVsCenty2(Double_t **&nJPsi, Double_t **&enJPsi, Int_t &nTests, Int_t &nBins)
{
  // return the number of JPsi in -3.25<y<-2 and errors for each test and each centrality bins
  // as well as the number of tests and of centrality bins
  
  const Int_t nT = 12;
  const Int_t nC = 8;
    
  ExtractNJPsiVsCent(0, 2, nC, nT, "1CB2", nJPsi, enJPsi);
  
  // return the values
  nTests = nT;
  nBins = nC;  
}


//------------------------------------------------------------------------
void SetNJPsiVsPt4(Double_t **&nJPsi, Double_t **&enJPsi, Int_t &nTests, Int_t &nBins)
{
  // return the number of JPsi with pT<3GeV/c and errors for each test and each centrality bins
  // as well as the number of tests and of centrality bins
  
  const Int_t nT = 12;
  const Int_t nC = 7;
  
  ExtractNJPsiVsPtY(4, 0, nC, nT, "1CB2", "090", nJPsi, enJPsi);
  
  // return the values
  nTests = nT;
  nBins = nC;  
}


//------------------------------------------------------------------------
void SetNJPsiVsPt5(Double_t **&nJPsi, Double_t **&enJPsi, Int_t &nTests, Int_t &nBins)
{
  // return the number of JPsi with pT>3GeV/c and errors for each test and each centrality bins
  // as well as the number of tests and of centrality bins
  
  const Int_t nT = 12;
  const Int_t nC = 7;
  
  ExtractNJPsiVsPtY(5, 0, nC, nT, "1CB2", "020", nJPsi, enJPsi);
  
  // return the values
  nTests = nT;
  nBins = nC;  
}


//------------------------------------------------------------------------
void SetNJPsiVsPt6(Double_t **&nJPsi, Double_t **&enJPsi, Int_t &nTests, Int_t &nBins)
{
  // return the number of JPsi with pT>3GeV/c and errors for each test and each centrality bins
  // as well as the number of tests and of centrality bins
  
  const Int_t nT = 12;
  const Int_t nC = 7;
  
  ExtractNJPsiVsPtY(6, 0, nC, nT, "1CB2", "2040", nJPsi, enJPsi);
  
  // return the values
  nTests = nT;
  nBins = nC;  
}


//------------------------------------------------------------------------
void SetNJPsiVsPt7(Double_t **&nJPsi, Double_t **&enJPsi, Int_t &nTests, Int_t &nBins)
{
  // return the number of JPsi in -4<y<-3.25 and errors for each test and each centrality bins
  // as well as the number of tests and of centrality bins
  
  const Int_t nT = 12;
  const Int_t nC = 7;
  
  ExtractNJPsiVsPtY(7, 0, nC, nT, "1CB2", "4090", nJPsi, enJPsi);
  
  // return the values
  nTests = nT;
  nBins = nC;  
}

//------------------------------------------------------------------------
void SetNJPsiVsPt8(Double_t **&nJPsi, Double_t **&enJPsi, Int_t &nTests, Int_t &nBins)
{
  // return the number of JPsi in -4<y<-3.25 and errors for each test and each centrality bins
  // as well as the number of tests and of centrality bins
  
  const Int_t nT = 12;
  const Int_t nC = 3;
  
  ExtractNJPsiVsPtY(8, 0, nC, nT, "1CB2", "6090", nJPsi, enJPsi);
  
  // return the values
  nTests = nT;
  nBins = nC;  
}


//------------------------------------------------------------------------
void SetNJPsiVsY3(Double_t **&nJPsi, Double_t **&enJPsi, Int_t &nTests, Int_t &nBins)
{
  // return the number of JPsi in -3.25<y<-2 and errors for each test and each centrality bins
  // as well as the number of tests and of centrality bins
  
  const Int_t nT = 12;
  const Int_t nC = 6;
   
  ExtractNJPsiVsPtY(0, 3, nC, nT, "1CB2", "090", nJPsi, enJPsi);
  
  // return the values
  nTests = nT;
  nBins = nC;
}





//------------------------------------------------------------------------
void ComputeNJpsiMean(Int_t pt, Int_t y, Double_t *&nJPsiMean, Double_t *&enJPsiMean, Double_t *&sysnJPsiMean, Int_t &nBins, Bool_t weightedMeanRMS)
{
  // return the average number of JPsi over all the tests and the corresponding
  // statistical and systematic errors for the given rapidity or pt bin
  
  // set results for all tests
  Int_t nTests;
  Double_t **nJPsi;
  Double_t **enJPsi;
  if (y == 0 && pt == 0) {
    printf("\nresults integrated over pt and y:\n\n");
    SetNJPsiVsCent(nJPsi, enJPsi, nTests, nBins);
  } 
  else if (y == 1 && pt == 0) {
    printf("\nresults integrated over pt in -3.25<y<-2.5:\n\n");
    SetNJPsiVsCenty1(nJPsi, enJPsi, nTests, nBins);
  } 
  else if (y == 2 && pt == 0) {
    printf("\nresults integrated over pt in -4<y<-3.25:\n\n");
    SetNJPsiVsCenty2(nJPsi, enJPsi, nTests, nBins);
  }
  else if (y == 3 && pt == 0) {
    printf("\nresults integrated over pt and  Y bins in cent 090%%:\n\n");
    SetNJPsiVsY3(nJPsi, enJPsi, nTests, nBins);
  }
  else if (y == 0 && pt == 1) {
    printf("\nresults integrated over y with 0<pt<2GeV/c:\n\n");
    SetNJPsiVsCentpt1(nJPsi, enJPsi, nTests, nBins);
  } 
  else if (y == 0 && pt == 2) {
    printf("\nresults integrated over y with 2<pt<5GeV/c:\n\n");
    SetNJPsiVsCentpt2(nJPsi, enJPsi, nTests, nBins);
  } 
  else if (y == 0 && pt == 3) {
    printf("\nresults integrated over y with 5<pt<8GeV/c:\n\n");
    SetNJPsiVsCentpt3(nJPsi, enJPsi, nTests, nBins);
  } 
  else if (y == 0 && pt == 4) {
    printf("\nresults integrated over y  and  Pt bins in cent 090%%:\n\n");
    SetNJPsiVsPt4(nJPsi, enJPsi, nTests, nBins);
  }
  else if (y == 0 && pt == 5) {
    printf("\nresults integrated over y  and  Pt bins in cent 020%%:\n\n");
    SetNJPsiVsPt5(nJPsi, enJPsi, nTests, nBins);
  }
  else if (y == 0 && pt == 6) {
    printf("\nresults integrated over y  and  Pt bins in cent 2040%%:\n\n");
    SetNJPsiVsPt6(nJPsi, enJPsi, nTests, nBins);
  }
  else if (y == 0 && pt == 7) {
    printf("\nresults integrated over y  and  Pt bins in cent 4090%%:\n\n");
    SetNJPsiVsPt7(nJPsi, enJPsi, nTests, nBins);
  }
  else if (y == 0 && pt == 8) {
    printf("\nresults integrated over y  and  Pt bins in cent 6090%%:\n\n");
    SetNJPsiVsPt8(nJPsi, enJPsi, nTests, nBins);
  }
  
  
  // define histograms
  TGraphErrors **gNJpsi = new TGraphErrors*[nBins];
  TH1F **hNJpsi = new TH1F*[nBins];
  TH1F **heNJpsi = new TH1F*[nBins];
  TString *gName, *gTitle;
  
  if (y<3 && pt<4) SetGraphNameTitleVsCent(nBins, gName, gTitle);
  else SetGraphNameTitleVsPtY(nBins, gName, gTitle);
  
  for(Int_t i=0; i<nBins; i++) {
    gNJpsi[i] = new TGraphErrors(nTests);
    gNJpsi[i]->SetName(Form("g%s_pt%d_y%d",gName[i].Data(),pt,y));
    gNJpsi[i]->SetTitle(Form("g%s_pt%d_y%d",gTitle[i].Data(),pt,y));
    hNJpsi[i] = new TH1F(Form("h%s_pt%d_y%d",gName[i].Data(),pt,y), Form("h%s_pt%d_y%d",gTitle[i].Data(),pt,y), 4200000, 0., 42000.);
    heNJpsi[i] = new TH1F(Form("he%s_pt%d_y%d",gName[i].Data(),pt,y), Form("he%s_pt%d_y%d",gTitle[i].Data(),pt,y), 55000, 0., 550.);
  }
  
  // loop over tests
  Double_t *enJPsiMean2 = new Double_t[nBins];
  for(Int_t j=0; j<nBins; j++) enJPsiMean2[j] = 0.;
  Double_t *noSys = new Double_t[nBins];
  Double_t *renJPsi, *RAA, *eRAAstat, *eRAA, *eRAA2;
  Double_t eRAACorr;
  TGraphErrors *g, *gSys, *gSysAll;
  TCanvas *c = new TCanvas(Form("compareRAA_pt%d_y%d",pt,y), Form("compareRAA -- pt%d -- y%d",pt,y));
  gPad->SetFrameBorderMode(0);
  for(Int_t i=0; i<nTests; i++) {
    
    // display current RAA
    ComputeRAA(pt, y, nBins, nJPsi[i], enJPsi[i], noSys, renJPsi, RAA, eRAAstat, eRAA, eRAA2, eRAACorr);
    if (nBins > 1) {
      BuildRAAGraphs(nBins, renJPsi, RAA, eRAAstat, eRAA, eRAA2, eRAACorr, g, gSys, gSysAll, kTRUE);
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
    for(Int_t j=0; j<nBins; j++) {
      
      Double_t wstat = 1./nJPsi[i][j];
      Double_t w = weightedMeanRMS ? 1./enJPsi[i][j]/enJPsi[i][j]/wstat : 1.;
      enJPsiMean2[j] += w*w*enJPsi[i][j]*enJPsi[i][j];
      
      // fill histogram
      gNJpsi[j]->SetPoint(i, i, nJPsi[i][j]);
      gNJpsi[j]->SetPointError(i, 0, enJPsi[i][j]);
      hNJpsi[j]->Fill(nJPsi[i][j], w);
      heNJpsi[j]->Fill(enJPsi[i][j], w);
      
    }
    
  }
  
  
  // compute NJpsi means and errors and display histograms
  nJPsiMean = new Double_t[nBins];
  enJPsiMean = new Double_t[nBins];
  sysnJPsiMean = new Double_t[nBins];
  TCanvas *sys1;
  if (fineDisplay) sys1 = new TCanvas(Form("sys1_pt%d_y%d",pt,y), Form("systematics of NJpsi -- pt%d -- y%d",pt,y), 1200, 180*nBins);
  else sys1 = new TCanvas(Form("sys1_pt%d_y%d",pt,y), Form("systematics of NJpsi -- pt%d -- y%d",pt,y), 400, 180*nBins);
  gPad->SetFrameBorderMode(0);
  if (fineDisplay) sys1->Divide(3,nBins);
  else sys1->Divide(1,nBins);
  for(Int_t i=0; i<nBins; i++) {
    
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
    
    // display histograms
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
  
  sys1->SaveAs(Form("Outputs/Syst_pt%i_y%i.pdf",pt,y));

  
}


//------------------------------------------------------------------------
void ComputeRAA(Int_t pt, Int_t y, Int_t nBins, Double_t *nJPsi, Double_t *enJPsi, Double_t *sysnJPsi, Double_t *&renJPsi,
				Double_t *&RAA, Double_t *&eRAAstat, Double_t *&eRAA, Double_t *&eRAA2, Double_t &eRAACorr, Bool_t print)
{
  // compute the RAA versus centrality in the given pt or rapidity bin
  
  // get graph name in nBins centrality bins for printout
  TString *gName, *gTitle;
  
  if (y<3 && pt<4) SetGraphNameTitleVsCent(nBins, gName, gTitle);
  else SetGraphNameTitleVsPtY(nBins, gName, gTitle);
  
  // get centrality independent normalization variables
  Double_t sigmaJPsipp, esigmaJPsippStat, esigmaJPsippSyst;
  Double_t *sigmaJPsippPtY, *esigmaJPsippStatPtY, *esigmaJPsippSystPtY;
  
  if (y<3 && pt<4) SetNormVar(pt, y, sigmaJPsipp, esigmaJPsippStat, esigmaJPsippSyst);
  else SetNormVarPtY(pt, y, sigmaJPsippPtY, esigmaJPsippStatPtY, esigmaJPsippSystPtY);
  
  // get normalization variables in nBins centrality bins
  Double_t *nMB, *TAA, *eTAA, *AccEff, *eAccEff;
  Double_t nMBPtY, TAAPtY, eTAAPtY;
  
  if (y<3 && pt<4) SetNormVarVsCent(pt, y, nBins, nMB, TAA, eTAA, AccEff, eAccEff);
  else SetNormVarVsPtY(pt, y, nMBPtY, TAAPtY, eTAAPtY, AccEff, eAccEff);
  
  // get the systematic errors in nBins centrality bins
  Double_t *eGen, *eRec, *eTrk, *eTrg;
  
  if (y<3 && pt<4) SetSystVarVsCent(nBins, eGen, eRec, eTrk, eTrg);
  else SetSystVarVsPtY(pt, y, eGen, eRec, eTrk, eTrg);
  
  
  
  if (y<3 && pt<4)   //  cases   centrality  bins
  {
  
	// Compute relative statistic & systematic errors
	renJPsi = new Double_t[nBins];
	Double_t *rsysnJPsi = new Double_t[nBins];
	Double_t *renMB = new Double_t[nBins];
	Double_t *reTAA = new Double_t[nBins];
	Double_t *reAccEff = new Double_t[nBins];
	for(Int_t i=0; i<nBins; i++) {    
	  renJPsi[i]=enJPsi[i]/nJPsi[i];
	  rsysnJPsi[i]=sysnJPsi[i]/nJPsi[i];
	  renMB[i]=1./TMath::Sqrt(nMB[i]);
	  reTAA[i]=eTAA[i]/TAA[i];
	  reAccEff[i]=eAccEff[i]/AccEff[i];
	}
	
	// compute RAA and related errors for each centrality bin
	RAA = new Double_t[nBins];
	eRAAstat = new Double_t[nBins];
	eRAA = new Double_t[nBins];
	eRAA2 = new Double_t[nBins];
	for(Int_t i=0; i<nBins; i++) { 
	  
	  RAA[i] = nJPsi[i]/BR/AccEff[i]/nMB[i]/sigmaJPsipp*1000./TAA[i]; // 1000. factor mub to mb 
	  
	  eRAAstat[i] = RAA[i] * renJPsi[i];
	  
	  eRAA[i] = RAA[i] * TMath::Sqrt(rsysnJPsi[i]*rsysnJPsi[i] +
									 eTrk[i]*eTrk[i] +
									 eTrg[i]*eTrg[i] +
									 renMB[i]*renMB[i]);
	  
	  eRAA2[i] = RAA[i] * reTAA[i];
	  
	  Double_t eRAATot = TMath::Sqrt(eRAA[i]*eRAA[i] + eRAA2[i]*eRAA2[i]);
	  
	  if (print) printf("RAA %s = %5.3f ± %5.3f ± %5.3f\n", gName[i].Data(), RAA[i], eRAAstat[i], eRAATot);
	  
	}
	if (print) printf("\n");
	
	//printf(" eGen %f eRec %f eTrk %f eTrg %f esigmaJPsippSyst %f eppLumi %f \n",eGen[nBins],eRec[nBins],eTrk[nBins],eTrg[nBins],esigmaJPsippSyst,eppLumi); 
	
	// correlated error**2 for RAA (without error on luminosity)
	Double_t eRAACorr2 = eGen[nBins]*eGen[nBins] +
						 eRec[nBins]*eRec[nBins] +
						 eTrk[nBins]*eTrk[nBins] +
						 eTrg[nBins]*eTrg[nBins] +
						 //eBR*eBR +
						 esigmaJPsippSyst*esigmaJPsippSyst +
						 esigmaJPsippStat*esigmaJPsippStat +
						 // eff error from pp
						 eppEff*eppEff +
						 //error fnorm
						 errfnorm*errfnorm;
	
	if (y==0 && pt==0) eRAACorr2 = eGen[nBins]*eGen[nBins] +
								   eRec[nBins]*eRec[nBins] +
								   eTrk[nBins]*eTrk[nBins] +
								   eTrg[nBins]*eTrg[nBins] +
								   //eBR*eBR +
								   esigmaJPsippSyst*esigmaJPsippSyst +
								   esigmaJPsippStat*esigmaJPsippStat +
								   // error of eff from pp
								   //eppEff*eppEff +
								   //error fnorm
								   errfnorm*errfnorm;
	
	if (print) printf("correlated RAA error (wo lumi) = %4.1f%%\n", 100.*TMath::Sqrt(eRAACorr2));
	
	// add relative error on the luminosity
	eRAACorr2 += eppLumi*eppLumi;
	
	// correlated error for RAA
	eRAACorr = TMath::Sqrt(eRAACorr2);
	
	if (print) printf("correlated RAA error = %4.1f%%\n\n", 100.*eRAACorr);
	
	
  }
  else    //  cases  Pt  and  Y  bins
  {
	
	// Compute relative statistic & systematic errors
	renJPsi = new Double_t[nBins];
	Double_t *rsysnJPsi = new Double_t[nBins];
	Double_t *reAccEff = new Double_t[nBins];
	for(Int_t i=0; i<nBins; i++) {    
	  renJPsi[i]=enJPsi[i]/nJPsi[i];
	  rsysnJPsi[i]=sysnJPsi[i]/nJPsi[i];
	  reAccEff[i]=eAccEff[i]/AccEff[i];
	}
	Double_t renMB = 1./TMath::Sqrt(nMBPtY);
	Double_t reTAA = eTAAPtY/TAAPtY;
	
	
	// compute RAA and related errors for each centrality bin
	RAA = new Double_t[nBins];
	eRAAstat = new Double_t[nBins];
	eRAA = new Double_t[nBins];
	eRAA2 = new Double_t[nBins];
	for(Int_t i=0; i<nBins; i++) { 
	  
	  RAA[i] = nJPsi[i]/BR/AccEff[i]/nMBPtY/sigmaJPsippPtY[i]*1000./TAAPtY; // 1000. factor mub to mb 
	  
	  eRAAstat[i] = RAA[i] * TMath::Sqrt(renJPsi[i]*renJPsi[i] + esigmaJPsippStatPtY[i]*esigmaJPsippStatPtY[i]);
	  
	  eRAA[i] = RAA[i] * TMath::Sqrt(rsysnJPsi[i]*rsysnJPsi[i] +
									 eTrk[0]*eTrk[0] +
									 eTrg[0]*eTrg[0] +
									 renMB*renMB);
	  
	  eRAA2[i] = RAA[i] * TMath::Sqrt(eGen[1]*eGen[1] +
									  eRec[1]*eRec[1] +
									  eTrk[1]*eTrk[1] +
									  eTrg[1]*eTrg[1] +
									  esigmaJPsippSystPtY[i]*esigmaJPsippSystPtY[i]);       
									  
	  
	  Double_t eRAATot = TMath::Sqrt(eRAA[i]*eRAA[i] + eRAA2[i]*eRAA2[i]);
	  
	  if (print) printf("RAA %s = %5.3f ± %5.3f ± %5.3f\n", gName[i].Data(), RAA[i], eRAAstat[i], eRAATot);
	  
	}
	if (print) printf("\n");
	
	//printf(" eGen %f eRec %f eTrk %f eTrg %f esigmaJPsippSyst %f eppLumi %f \n",eGen[nBins],eRec[nBins],eTrk[nBins],eTrg[nBins],esigmaJPsippSyst,eppLumi); 
	
	// correlated error**2 for RAA (without error on luminosity)
	Double_t eRAACorr2 = reTAA*reTAA +
						 // error of eff from pp
						 eppEff*eppEff + 
						 //eBR*eBR +
						 //error fnorm
						 errfnorm*errfnorm;
	
	if (print) printf("correlated RAA error (wo lumi) = %4.1f%%\n", 100.*TMath::Sqrt(eRAACorr2));
	
	// add relative error on the luminosity
	eRAACorr2 += eppLumi*eppLumi;
	
	// correlated error for RAA
	eRAACorr = TMath::Sqrt(eRAACorr2);
	
	if (print) printf("correlated RAA error = %4.1f%%\n\n", 100.*eRAACorr);
	
  }
  
  
  
}


//------------------------------------------------------------------------
void BuildRAAGraphs(Int_t nBins, Double_t *renJPsi, Double_t *RAA, Double_t *eRAAstat, Double_t *eRAA, Double_t *eRAA2, Double_t eRAACorr,
					TGraphErrors *&gRAA, TGraphErrors *&gRAASys, TGraphErrors *&gRAASysAll, Int_t xAxisType, Int_t color, Int_t sysPos)
{
  // build the TGraphs to display our RAA results
  
  // Get graph x-axis labels according to the number of centrality bins
  Double_t *x, *ex, *bin, *NPart, *eNPart, *dNchdEta, *edNchdEta, *dNchdEtaNorm, *edNchdEtaNorm;;
  TString *label;
  
  if (nBins>7) {
	
	SetGraphLabelsVsCent(nBins, x, ex, bin, label, NPart, eNPart,dNchdEta, edNchdEta, dNchdEtaNorm, edNchdEtaNorm);
	
	// fill RAA graphs
	gRAA = new TGraphErrors(nBins-1);
	gRAASys = new TGraphErrors(nBins-1);
	gRAASysAll = new TGraphErrors(1);
	if (xAxisType==0) {
	  gRAASysAll->SetPoint(0,99.-3.*sysPos,1.);
	  gRAASysAll->SetPointError(0,1.,eRAACorr);
	} 
	else if (xAxisType==1){
	  gRAASysAll->SetPoint(0,395.-12.*sysPos,1.);
	  gRAASysAll->SetPointError(0,5.,eRAACorr);
	} 
	else if (xAxisType==2){
	  gRAASysAll->SetPoint(0,1505,1.);
	  gRAASysAll->SetPointError(0,15,eRAACorr);
	} else if (xAxisType==3){
	  gRAASysAll->SetPoint(0,9.4,1.);
	  gRAASysAll->SetPointError(0,.1,eRAACorr);
	}
	
	for(Int_t i=nBins-1; i>=1; i--) {
	  if (xAxisType==0) {
		gRAA->SetPoint(nBins-1-i,100.-x[i],RAA[i]);
		gRAA->SetPointError(nBins-1-i,ex[i],eRAAstat[i]);
		gRAASys->SetPoint(nBins-1-i,100.-x[i],RAA[i]);
		gRAASys->SetPointError(nBins-1-i,1.,TMath::Sqrt(eRAA[i]*eRAA[i] + eRAA2[i]*eRAA2[i]));
	  } 
	  else if (xAxisType==1){
		gRAA->SetPoint(nBins-1-i,NPart[i],RAA[i]);
		gRAA->SetPointError(nBins-1-i,eNPart[i],eRAAstat[i]);
		gRAASys->SetPoint(nBins-1-i,NPart[i],RAA[i]);
		gRAASys->SetPointError(nBins-1-i,5.,TMath::Sqrt(eRAA[i]*eRAA[i] + eRAA2[i]*eRAA2[i]));
	  } 
	  else if (xAxisType==2){
		gRAA->SetPoint(nBins-1-i,dNchdEta[i],RAA[i]);
		gRAA->SetPointError(nBins-1-i,edNchdEta[i],eRAAstat[i]);
		gRAASys->SetPoint(nBins-1-i,dNchdEta[i],RAA[i]);
		gRAASys->SetPointError(nBins-1-i,10,TMath::Sqrt(eRAA[i]*eRAA[i] + eRAA2[i]*eRAA2[i]));
	  }else if (xAxisType==3){
		gRAA->SetPoint(nBins-1-i,dNchdEtaNorm[i],RAA[i]);
		gRAA->SetPointError(nBins-1-i,edNchdEtaNorm[i],eRAAstat[i]);
		gRAASys->SetPoint(nBins-1-i,dNchdEtaNorm[i],RAA[i]);
		gRAASys->SetPointError(nBins-1-i,0.4,TMath::Sqrt(eRAA[i]*eRAA[i] + eRAA2[i]*eRAA2[i]));
	  }
	}
	
	if (xAxisType==0) {
	  gRAA->GetXaxis()->Set(nBins, bin);
	  for(Int_t i=nBins; i>=1; i--)
		gRAA->GetXaxis()->SetBinLabel(nBins+1-i, label[i].Data());
	} 
	else if (xAxisType==1){
	  gRAA->GetXaxis()->Set(40, 0., 400.);
	}
	else if (xAxisType==2){
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
	} 
	
	else if (xAxisType==1){
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
	gRAA->SetMarkerSize(1.);
	gRAA->SetMarkerColor(color);
	gRAA->SetLineWidth(2);
	gRAA->SetLineColor(color);
	gRAASys->SetLineColor(color);
	gRAASys->SetFillStyle(0);
	gRAASysAll->SetLineWidth(1);
	gRAASysAll->SetFillStyle(3001);
	gRAASysAll->SetFillColor(color);
	
	gRAA->SetTitle("");
  }
  
  
  else {
	
	SetGraphLabelsVsPtY(nBins, x, ex);
	
	// fill RAA graphs
	gRAA = new TGraphErrors(nBins);
	gRAASys = new TGraphErrors(nBins);
	gRAASysAll = new TGraphErrors(1);
	
	if (nBins==6) {
	  gRAASysAll->SetPoint(0,3.98125-3.*sysPos,1);
	  gRAASysAll->SetPointError(0,0.01875,eRAACorr);
	}
	if (nBins==7) {
	  gRAASysAll->SetPoint(0,7.9-3.*sysPos,1);
	  gRAASysAll->SetPointError(0,0.1,eRAACorr);
	}
		
	for(Int_t i=0; i<nBins; i++) {
	  
	  gRAA->SetPoint(i, x[i], RAA[i]);
	  gRAA->SetPointError(i, ex[i], eRAAstat[i]);
	  gRAASys->SetPoint(i, x[i], RAA[i]);
	  gRAASys->SetPointError(i, ex[i], TMath::Sqrt(eRAA[i]*eRAA[i] + eRAA2[i]*eRAA2[i]));
	}

	// define axis
	if (nBins==6) 
	{
	  gRAA->GetXaxis()->Set(40, 2.5, 4.); // 6
	  gRAA->GetXaxis()->SetTitle("y");
	}
	else if (nBins==7) 
	{
	  gRAA->GetXaxis()->Set(40, 0., 8.);  // 8
	  gRAA->GetXaxis()->SetTitle("p_{T} GeV/c");
	}
	
	gRAA->GetXaxis()->SetLabelSize(0.05); // 0.07
	gRAA->GetXaxis()->SetTitleOffset(1.1);
	gRAA->GetXaxis()->SetTitleSize(0.05);
	//gRAA->GetXaxis()->SetTickLength(0.);
	
	gRAA->GetYaxis()->SetTitle("R_{AA}");
	gRAA->GetYaxis()->SetTitleSize(0.05);
	gRAA->GetYaxis()->SetTitleOffset(0.9);
	gRAA->GetYaxis()->SetLabelSize(0.05);
	gRAA->GetYaxis()->SetRangeUser(0.,1.4);

	
	// define display
	gRAA->SetMarkerStyle(gMarkerStyleAliceForw);
	gRAA->SetMarkerSize(1.);
	gRAA->SetMarkerColor(color);
	gRAA->SetLineWidth(2);
	gRAA->SetLineColor(color);
	gRAASys->SetLineColor(color);
	gRAASys->SetFillStyle(0);
	gRAASysAll->SetLineWidth(1);
	gRAASysAll->SetFillStyle(3001);
	gRAASysAll->SetFillColor(color);
	
	gRAA->SetTitle("");
  }
  
}


//------------------------------------------------------------------------
void BuildForwardPHENIXRAAGraphs(TGraphErrors *&gRAA_PHENIX, TGraphAsymmErrors *&gRAA_PHENIXSys,
								 TGraphErrors *&gRAA_PHENIXSysAll, Int_t xAxisType)
{
  // build the TGraphs to display PHENIX RAA results at forward rapidity
  
  // fill RAA graphs
  Double_t cent_PHENIX[17] = {2.5, 7.5, 12.5, 17.5, 22.5, 27.5, 32.5, 37.5, 42.5, 47.5, 52.5, 57.5, 62.5, 67.5, 72.5, 77.5, 86.};
  Double_t ecent_PHENIX[17] = {2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 6.};
  
  
  
  //Double_t NPart_PHENIX[17] = {350.8, 301.7, 255.7, 216.4, 182.4, 152.7, 126.8, 104.2, 84.6, 67.7, 53.2, 41., 30.8, 22.6, 16.1, 11.2, 5.6};
  //Double_t eNPart_PHENIX[17] = {3.1, 4.7, 5.4, 5.6, 5.7, 5.9, 5.9, 5.8, 5.6, 5.4, 5., 4.5, 3.9, 3.4, 2.8, 2.2, 0.8};
  
  Double_t NPart_PHENIX[5] = {356.5,260.5,186.4,117.3,46.5};
  Double_t eNPart_PHENIX[5] = {3.6,4.4,3.9,3.4,1.7};
  
    
  
  
  
  Double_t dNchdEta_PHENIX[17] = {687,560,457,372,302,246,197,156,124,95.3,70.9,52.2,37.5,25.6,-100.0,-100.0, -100.0};
  Double_t edNchdEta_PHENIX[17] = {37,28,22,18,16,14,12,11,9.6,8.6,7.6,6.5,5.4,4.5,1,1,1};
  Double_t dNchdEtaNorm_PHENIX[17] = {3.89, 3.73, 3.59, 3.45, 3.34, 3.25, 3.15, 3.05, 2.96, 2.86, 2.7, 2.6, 2.48, 2.33, 1.0, 1.0, 1.0};
  Double_t edNchdEtaNorm_PHENIX[17] = {0.23,0.22,0.21,0.21,0.21,0.22,0.24,0.26,0.28,0.32,0.36,0.41,0.46,0.55,1,1,1};
  
  
  
  
  //Double_t RAA_PHENIX[17] = {0.17, 0.16, 0.20, 0.25, 0.25, 0.35, 0.35, 0.41, 0.52, 0.49, 0.54, 0.80, 0.68, 0.72, 0.91, 1.03, 1.20};
  //Double_t eRAA_PHENIX[17] = {0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.05, 0.06, 0.06, 0.08, 0.11, 0.10};
  //Double_t syspRAA_PHENIX[17] = {0.04, 0.02, 0.03, 0.04, 0.03, 0.04, 0.04, 0.06, 0.07, 0.07, 0.09, 0.14, 0.13, 0.15, 0.21, 0.26, 0.23};
  //Double_t sysmRAA_PHENIX[17] = {0.02, 0.02, 0.03, 0.04, 0.03, 0.04, 0.04, 0.06, 0.07, 0.07, 0.09, 0.14, 0.13, 0.15, 0.21, 0.26, 0.23};
  
  Double_t RAA_PHENIX[5] = {0.5559455, 0.5171857, 0.5940191, 0.4964814, 0.6814482};
  Double_t eRAA_PHENIX[5] = {0.05856555, 0.06327899, 0.06983689, 0.06074629, 0.08696302};
  Double_t syspRAA_PHENIX[5] = {0.05979968, 0.04133377, 0.07610803, 0.04259053, 0.05724437};
  Double_t sysmRAA_PHENIX[5] = {0.05979968, 0.04133377, 0.07610803, 0.04259053, 0.05724437};
  
  
  gRAA_PHENIX = new TGraphErrors(5);
  gRAA_PHENIXSys = new TGraphAsymmErrors(5);
  for (Int_t i=4; i>=0; i--) {
    if (xAxisType==0) {
      gRAA_PHENIX->SetPoint(16-i, 100.-cent_PHENIX[i], RAA_PHENIX[i]);
      gRAA_PHENIX->SetPointError(16-i, ecent_PHENIX[i], eRAA_PHENIX[i]);
      gRAA_PHENIXSys->SetPoint(16-i, 100.-cent_PHENIX[i], RAA_PHENIX[i]);
      gRAA_PHENIXSys->SetPointError(16-i, 1., 1., sysmRAA_PHENIX[i], syspRAA_PHENIX[i]);
    } else if (xAxisType==1){
      gRAA_PHENIX->SetPoint(4-i, NPart_PHENIX[i], RAA_PHENIX[i]);
      gRAA_PHENIX->SetPointError(4-i, eNPart_PHENIX[i], eRAA_PHENIX[i]);
      gRAA_PHENIXSys->SetPoint(4-i, NPart_PHENIX[i], RAA_PHENIX[i]);
      gRAA_PHENIXSys->SetPointError(4-i, 5., 5., sysmRAA_PHENIX[i], syspRAA_PHENIX[i]);
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
    gRAA_PHENIXSysAll->SetPointError(0,5.,0.1245976);
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
void DisplayRAA(Int_t pt, Int_t y, Int_t nBins, TGraphErrors *gRAA, TGraphErrors *gRAASys, TGraphErrors *gRAASysAll, Int_t xAxisType)
{
  // display our RAA results
  
  TCanvas* c = new TCanvas(Form("RAA_pt%d_y%d",pt,y),Form("RAA -- pt%d -- y%d",pt,y),gCanvasWidth,gCanvasHeigh);
  gPad->SetTopMargin(0.05);
  gPad->SetBottomMargin(0.12);
  gPad->SetRightMargin(0.05);
  gRAA->Draw("ap");
  gRAASys->Draw("e2");
  gRAASysAll->Draw("e2");
  
  TLine *l = 0x0;  
  if (xAxisType==0) l = new TLine(0., 1., 100., 1.);
  else if (xAxisType==1){
	if (y<3 && pt<4) l =  new TLine(0., 1., 400., 1.);
	else if (y==3) l =  new TLine(2.5, 1., 4., 1.);
	else l =  new TLine(0., 1., 8., 1.);
  }
  else if (xAxisType==2) l =  new TLine(0., 1., 1550., 1.);
  else if (xAxisType==3) l =  new TLine(0., 1., 10., 1.);
  l->SetLineStyle(3);
  l->Draw();
  if (xAxisType==0) {
    if (nBins == 6) {
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
	
  
  TPaveText* t = new TPaveText(0.18, 0.83, 0.87, 0.93,"NDC");
  t->AddText(0.,0.,"Inclusive J/#psi in Pb-Pb  #sqrt{s_{NN}} = 2.76 TeV");
  t->SetFillStyle(0);
  t->SetBorderSize(0);
  t->SetTextFont(42);
  t->Draw();
  
  
  TLegend* tlg = new TLegend(0.25, 0.77, 0.9, 0.83,"","NDC");
  if (y == 0 && pt == 0) tlg->AddEntry(gRAA,"ALICE  2.5<y<4, 0<p_{T}<8 GeV/c","p");
  else if (y == 1 && pt == 0) tlg->AddEntry(gRAA,"ALICE  2.5<y<3.25, 0<p_{T}<8 GeV/c","p");
  else if (y == 2 && pt == 0) tlg->AddEntry(gRAA,"ALICE  3.25<y<4, 0<p_{T}<8 GeV/c  ","p");
  else if (y == 3 && pt == 0) tlg->AddEntry(gRAA,"ALICE  0<p_{T}<8 GeV/c, cent 0-90%","p");
  else if (y == 0 && pt == 1) tlg->AddEntry(gRAA,"ALICE  2.5<y<4, 0<p_{T}<2 GeV/c   ","p");
  else if (y == 0 && pt == 2) tlg->AddEntry(gRAA,"ALICE  2.5<y<4, 2<p_{T}<5 GeV/c   ","p");
  else if (y == 0 && pt == 3) tlg->AddEntry(gRAA,"ALICE  2.5<y<4, 5<p_{T}<8 GeV/c   ","p");
  else if (y == 0 && pt == 4) tlg->AddEntry(gRAA,"ALICE  2.5<y<4, cent 0-90%  ","p");
  else if (y == 0 && pt == 5) tlg->AddEntry(gRAA,"ALICE  2.5<y<4, cent 0-20%  ","p");
  else if (y == 0 && pt == 6) tlg->AddEntry(gRAA,"ALICE  2.5<y<4, cent 20-40% ","p");
  else if (y == 0 && pt == 7) tlg->AddEntry(gRAA,"ALICE  2.5<y<4, cent 40-90% ","p");
  else if (y == 0 && pt == 8) tlg->AddEntry(gRAA,"ALICE  2.5<y<4, cent 60-90% ","p");

  tlg->SetFillStyle(0);
  tlg->SetBorderSize(0);
  tlg->SetTextFont(42);
  tlg->SetMargin(0.05);
  tlg->Draw();
  //ALICEseal("ALICE preliminary", 0.2, 0.2);
  //ALICEseal("", 0.2, 0.2);
  
  c->SaveAs(Form("Outputs/Raa_pt%i_y%i.pdf",pt,y));
  //c->SaveAs(Form("Outputs/Raa_pt%i_y%i.C",pt,y));
  //c->SaveAs(Form("Outputs/Raa_pt%i_y%i.root",pt,y));
  
  
}


//------------------------------------------------------------------------
void DisplayRAA(TString var, Int_t nBins, TGraphErrors *gRAA[8], TGraphErrors *gRAASys[8], TGraphErrors *gRAASysAll[8], Int_t xAxisType)
{
  // display our RAA results
  
  Int_t nbegin=0, nend = 3;
  TCanvas *c;
  if (var == "y") c = new TCanvas("RAA_pt0_yAll","RAA -- pt0 -- yAll",gCanvasWidth,gCanvasHeigh);
  else if (var == "pt") {nbegin=0; nend=3; c = new TCanvas("RAA_ptAll_y0","RAA -- ptAll -- y0",gCanvasWidth,gCanvasHeigh);}
  else {nbegin=4; nend=8; c = new TCanvas("RAA_ptBins_centAll_y0","RAA -- ptBins -- centAll -- y0",gCanvasWidth,gCanvasHeigh);}
  
  gPad->SetTopMargin(0.05);
  gPad->SetBottomMargin(0.12);
  gPad->SetRightMargin(0.05);
  
  for (Int_t i=nbegin; i<nend; i++) {
    if (i==nbegin) gRAA[i]->Draw("ap");
    else gRAA[i]->Draw("p");
    gRAASys[i]->Draw("e2");
    
	Double_t yRAASysAll = gRAASysAll[i]->GetErrorY(0);
	if (var == "y" || var == "pt") {
	  gRAASysAll[i]->SetPoint(0,395.-12*(i-nbegin),1.);
	  gRAASysAll[i]->SetPointError(0,5.,yRAASysAll);
	}
	else {
	  gRAASysAll[i]->SetPoint(0,7.9-0.24*(i-nbegin),1.);
	  gRAASysAll[i]->SetPointError(0,0.1,yRAASysAll);
	}
	gRAASysAll[i]->Draw("e2");
  }
  
  TLine *l = 0x0; 
  if  (xAxisType==0)  l = new TLine(0., 1., 100., 1.);
  else if (xAxisType==1) 
  {
	if (var=="y" || var=="pt") l = new TLine(0., 1., 400., 1.);
	else l = new TLine(0., 1., 8., 1.);
  }
  else if (xAxisType==2) l= new TLine(0., 1., 1550., 1.);
  else if (xAxisType==3) l= new TLine(0., 1., 10., 1.);
  l->SetLineStyle(3);
  l->Draw();
  
  if (xAxisType==0) {
    if (nBins == 6) {
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
  if (var == "y") {
	lg->AddEntry(gRAA[0],"ALICE (Pb-Pb #sqrt{s_{NN}} = 2.76 TeV), 0<p_{T}<8, 2.5<y<4 ","p");
    lg->AddEntry(gRAA[1],"ALICE (Pb-Pb #sqrt{s_{NN}} = 2.76 TeV), 0<p_{T}<8, 2.5<y<3.25 ","p");
    lg->AddEntry(gRAA[2],"ALICE (Pb-Pb #sqrt{s_{NN}} = 2.76 TeV), 0<p_{T}<8, 3.25<y<4 ","p");
  } 
  else if (var == "pt") {
    lg->AddEntry(gRAA[0],"ALICE (Pb-Pb #sqrt{s_{NN}} = 2.76 TeV), 2.5<y<4, 0<p_{T}<2 GeV/c ","p");
    lg->AddEntry(gRAA[1],"ALICE (Pb-Pb #sqrt{s_{NN}} = 2.76 TeV), 2.5<y<4, 2<p_{T}<5 GeV/c ","p");
	lg->AddEntry(gRAA[2],"ALICE (Pb-Pb #sqrt{s_{NN}} = 2.76 TeV), 2.5<y<4, 5<p_{T}<8 GeV/c ","p");
  }
  else {
	//lg->AddEntry(gRAA[3],"ALICE (Pb-Pb #sqrt{s_{NN}} = 2.76 TeV), 2.5<y<4, cent 0-90% ","p");
    lg->AddEntry(gRAA[4],"ALICE (Pb-Pb #sqrt{s_{NN}} = 2.76 TeV), 2.5<y<4, cent 0-20% ","p");
    lg->AddEntry(gRAA[5],"ALICE (Pb-Pb #sqrt{s_{NN}} = 2.76 TeV), 2.5<y<4, cent 20-40% ","p");
	lg->AddEntry(gRAA[6],"ALICE (Pb-Pb #sqrt{s_{NN}} = 2.76 TeV), 2.5<y<4, cent 40-90% ","p");
	lg->AddEntry(gRAA[7],"ALICE (Pb-Pb #sqrt{s_{NN}} = 2.76 TeV), 2.5<y<4, cent 60-90% ","p");
  }
  lg->SetFillStyle(0);
  lg->SetBorderSize(0);
  lg->SetTextFont(42);
  lg->SetMargin(0.05);
  lg->Draw();
  //ALICEseal("ALICE preliminary", 0.2, 0.2);
  //ALICEseal("", 0.2, 0.2);
  
  c->SaveAs(Form("Outputs/Raa_All_%s.pdf", var.Data()));
  
}

//------------------------------------------------------------------------
void DisplayRAAratio(TString var, Int_t nBins, TGraphErrors *gRAA[8], TGraphErrors *gRAASys[8], TGraphErrors *gRAASysAll[8])
{
  
  
  if (var == "y") 
  {
	Double_t *x, *ex, *bin, *NPart, *eNPart, *dNchdEta, *edNchdEta, *dNchdEtaNorm, *edNchdEtaNorm;;
	TString *label;
    SetGraphLabelsVsCent(nBins+1, x, ex, bin, label, NPart, eNPart,dNchdEta, edNchdEta, dNchdEtaNorm, edNchdEtaNorm);
	  
	
	Double_t **yRAA = new Double_t*[2];
	Double_t **eyRAA = new Double_t*[2];
	Double_t **yRAAsys = new Double_t*[2];
	Double_t *yRAAsysAll = new Double_t[2];
	
	Double_t **nMBevt = new Double_t*[2];
	Double_t **TAAcent = new Double_t*[2];
	Double_t **eTAAcent = new Double_t*[2];
	Double_t sigmaJPsipp[2], esigmaJPsippStat[2], esigmaJPsippSyst[2];
	
	TGraphErrors *gRAAratio = new TGraphErrors(nBins);
	TGraphErrors *gRAASysratio = new TGraphErrors(nBins);
	TGraphErrors *gRAASysAllratio = new TGraphErrors(1);
	
	
	// Get values of Raa
	for (Int_t i=1; i<3; i++) 
	{
	  yRAA[i-1] = new Double_t[nBins];
	  yRAA[i-1] = gRAA[i]->GetY();
	  eyRAA[i-1] = new Double_t[nBins];
	  eyRAA[i-1] = gRAA[i]->GetEY();
	  yRAAsys[i-1] = new Double_t[nBins];
	  yRAAsys[i-1] = gRAASys[i]->GetEY();
	  yRAAsysAll[i-1] = gRAASysAll[i]->GetErrorY(0);
	  
	  // get normalization variables in nBins centrality bins
	  Double_t *nMB, *TAA, *eTAA, *AccEff, *eAccEff;
	  SetNormVarVsCent(0, i, nBins+1, nMB, TAA, eTAA, AccEff, eAccEff);
	  
	  nMBevt[i-1] = new Double_t[nBins+1];
	  nMBevt[i-1] = nMB;
	  TAAcent[i-1] = new Double_t[nBins+1];
	  TAAcent[i-1] = TAA;
	  eTAAcent[i-1] = new Double_t[nBins+1];
	  eTAAcent[i-1] = eTAA;
	  
	  // get centrality independent normalization variables
	  SetNormVar(0, i, sigmaJPsipp[i-1], esigmaJPsippStat[i-1], esigmaJPsippSyst[i-1]);
	  
	}

	// Computation and Build graph
	for (Int_t j=0; j<nBins; j++)
	{
	  Double_t yRAAratio = yRAA[1][j] / yRAA[0][j];
	  Double_t reRAA1 = eyRAA[0][j] / yRAA[0][j];
	  Double_t reRAA2 = eyRAA[1][j] / yRAA[1][j];
	  Double_t eraaRatio = TMath::Sqrt( reRAA1*reRAA1 + reRAA2*reRAA2 );
	  
	  Double_t renMB1= 1. / TMath::Sqrt(nMBevt[0][nBins-j]);
	  Double_t reTAA1= eTAAcent[0][nBins-j] / TAAcent[0][nBins-j]; 
	  Double_t yRAAsys1 = yRAA[0][j] * TMath::Sqrt( (yRAAsys[0][j]*yRAAsys[0][j]) / (yRAA[0][j]*yRAA[0][j]) - renMB1*renMB1 - reTAA1*reTAA1 );
	  Double_t renMB2= 1. / TMath::Sqrt(nMBevt[1][nBins-j]);
	  Double_t reTAA2= eTAAcent[1][nBins-j] / TAAcent[1][nBins-j]; 
	  Double_t yRAAsys2 = yRAA[1][j] * TMath::Sqrt( (yRAAsys[1][j]*yRAAsys[1][j]) / (yRAA[1][j]*yRAA[1][j]) - renMB2*renMB2 - reTAA2*reTAA2 );
	  Double_t yRAAsysRatio = TMath::Sqrt( yRAAsys1*yRAAsys1 + yRAAsys2*yRAAsys2 ) * yRAAratio;
	  	  
	  gRAAratio->SetPoint(j, NPart[nBins-j], yRAAratio);
	  gRAAratio->SetPointError(j, eNPart[nBins-j], eraaRatio);
	  gRAASysratio->SetPoint(j, NPart[nBins-j], yRAAratio);
	  gRAASysratio->SetPointError(j, 5., yRAAsysRatio);
	}
		
	Double_t yRAAsysAll1 = TMath::Sqrt( yRAAsysAll[0]*yRAAsysAll[0] - eppEff*eppEff - errfnorm*errfnorm - eppLumi*eppLumi + esigmaJPsippSyst[1]*esigmaJPsippSyst[1] );
	Double_t yRAAsysAllRatio = yRAAsysAll1;
	gRAASysAllratio->SetPoint(0,395.,1.);
	gRAASysAllratio->SetPointError(0,5.,yRAAsysAllRatio);
	
	
	// display our RAA results
	Int_t color = kRed;
	gRAAratio->GetXaxis()->Set(40, 0., 400.);
	gRAAratio->SetTitle("");
	
	// define axis
	gRAAratio->GetXaxis()->SetTitle("#LT N_{part} #GT");
	gRAAratio->GetXaxis()->SetTitleSize(0.05); 
	gRAAratio->GetXaxis()->SetTitleOffset(1.0);
	gRAAratio->GetXaxis()->SetLabelSize(0.05);
	
	gRAAratio->GetYaxis()->SetTitle("Ratio R_{AA}");
	gRAAratio->GetYaxis()->SetTitleSize(0.05);
	gRAAratio->GetYaxis()->SetTitleOffset(0.9);
	gRAAratio->GetYaxis()->SetLabelSize(0.05);
	gRAAratio->GetYaxis()->SetRangeUser(0.5,1.25);
	
	// define display
	gRAAratio->SetMarkerStyle(gMarkerStyleAliceForw);
	gRAAratio->SetMarkerSize(1.);
	gRAAratio->SetMarkerColor(color);
	gRAAratio->SetLineWidth(2);
	gRAAratio->SetLineColor(color);
	gRAASysratio->SetLineColor(color);
	gRAASysratio->SetFillStyle(0);
	gRAASysAllratio->SetLineWidth(1);
	gRAASysAllratio->SetFillStyle(3001);
	gRAASysAllratio->SetFillColor(color);
	

	TCanvas *c1 = new TCanvas("RAA_pt0_yAll_Ratio","RAA -- pt0 -- yAll -- Ratio",gCanvasWidth,gCanvasHeigh);
	gPad->SetTopMargin(0.05);
	gPad->SetBottomMargin(0.12);
	gPad->SetRightMargin(0.05);
	
	gRAAratio->Draw("ap");
    gRAASysratio->Draw("e2");
    gRAASysAllratio->Draw("e2");
	
	TLine *l =  new TLine(0., 1., 400., 1.);
	l->SetLineStyle(3);
	l->Draw();
	
	TLegend *lg = new TLegend(0.17, 0.81, 0.92, 0.93,"","NDC");
	lg->AddEntry(gRAAratio," 3.25<y<4 / 2.5<y<3.25 ","p");
	lg->SetFillStyle(0);
	lg->SetBorderSize(0);
	lg->SetTextFont(42);
	lg->SetMargin(0.05);
	lg->Draw();
	
	c1->SaveAs("Outputs/Raa_pt0_yAll_Ratio.pdf");
	
  }
  
  
  
  else if (var == "pt") 
  {
	Double_t *x, *ex, *bin, *NPart, *eNPart, *dNchdEta, *edNchdEta, *dNchdEtaNorm, *edNchdEtaNorm;;
	TString *label;
    SetGraphLabelsVsCent(nBins+1, x, ex, bin, label, NPart, eNPart,dNchdEta, edNchdEta, dNchdEtaNorm, edNchdEtaNorm);
	
	
	Double_t **yRAA = new Double_t*[3];
	Double_t **eyRAA = new Double_t*[3];
	Double_t **yRAAsys = new Double_t*[3];
	Double_t *yRAAsysAll = new Double_t[3];
	
	Double_t **nMBevt = new Double_t*[3];
	Double_t **TAAcent = new Double_t*[3];
	Double_t **eTAAcent = new Double_t*[3];
	Double_t sigmaJPsipp[3], esigmaJPsippStat[3], esigmaJPsippSyst[3];
	
	
	// Get values of Raa
	for (Int_t i=0; i<3; i++) 
	{
	  yRAA[i] = new Double_t[nBins];
	  yRAA[i] = gRAA[i]->GetY();
	  eyRAA[i] = new Double_t[nBins];
	  eyRAA[i] = gRAA[i]->GetEY();
	  yRAAsys[i] = new Double_t[nBins];
	  yRAAsys[i] = gRAASys[i]->GetEY();
	  yRAAsysAll[i] = gRAASysAll[i]->GetErrorY(0);
	  
	  // get normalization variables in nBins centrality bins
	  Double_t *nMB, *TAA, *eTAA, *AccEff, *eAccEff;
	  SetNormVarVsCent(i, 0, nBins+1, nMB, TAA, eTAA, AccEff, eAccEff);
	  
	  nMBevt[i] = new Double_t[nBins+1];
	  nMBevt[i] = nMB;
	  TAAcent[i] = new Double_t[nBins+1];
	  TAAcent[i] = TAA;
	  eTAAcent[i] = new Double_t[nBins+1];
	  eTAAcent[i] = eTAA;
	  
	  // get centrality independent normalization variables
	  SetNormVar(i, 0, sigmaJPsipp[i], esigmaJPsippStat[i], esigmaJPsippSyst[i]);
	  
	}
	
	
	TGraphErrors *gRAAratio = new TGraphErrors(nBins);
	TGraphErrors *gRAASysratio = new TGraphErrors(nBins);
	TGraphErrors *gRAASysAllratio = new TGraphErrors(1);
	
	TGraphErrors *gRAAratio1 = new TGraphErrors(nBins);
	TGraphErrors *gRAASysratio1 = new TGraphErrors(nBins);
	TGraphErrors *gRAASysAllratio1 = new TGraphErrors(1);

	TGraphErrors *gRAAratio2 = new TGraphErrors(nBins);
	TGraphErrors *gRAASysratio2 = new TGraphErrors(nBins);
	TGraphErrors *gRAASysAllratio2 = new TGraphErrors(1);

	
	// Computation and Build graph
	for (Int_t j=0; j<nBins; j++)
	{
	  Double_t yRAAratio = yRAA[1][j] / yRAA[0][j];
	  Double_t reRAA1 = eyRAA[0][j] / yRAA[0][j];
	  Double_t reRAA2 = eyRAA[1][j] / yRAA[1][j];
	  Double_t eraaRatio = TMath::Sqrt( reRAA1*reRAA1 + reRAA2*reRAA2 );
	  
	  Double_t renMB1= 1. / TMath::Sqrt(nMBevt[0][nBins-j]);
	  Double_t reTAA1= eTAAcent[0][nBins-j] / TAAcent[0][nBins-j]; 
	  Double_t yRAAsys1 = yRAA[0][j] * TMath::Sqrt( (yRAAsys[0][j]*yRAAsys[0][j]) / (yRAA[0][j]*yRAA[0][j]) - renMB1*renMB1 - reTAA1*reTAA1 );
	  Double_t renMB2= 1. / TMath::Sqrt(nMBevt[1][nBins-j]);
	  Double_t reTAA2= eTAAcent[1][nBins-j] / TAAcent[1][nBins-j]; 
	  Double_t yRAAsys2 = yRAA[1][j] * TMath::Sqrt( (yRAAsys[1][j]*yRAAsys[1][j]) / (yRAA[1][j]*yRAA[1][j]) - renMB2*renMB2 - reTAA2*reTAA2 );
	  Double_t yRAAsysRatio = TMath::Sqrt( yRAAsys1*yRAAsys1 + yRAAsys2*yRAAsys2 ) * yRAAratio;
	  
	  gRAAratio->SetPoint(j, NPart[nBins-j], yRAAratio);
	  gRAAratio->SetPointError(j, eNPart[nBins-j], eraaRatio);
	  gRAASysratio->SetPoint(j, NPart[nBins-j], yRAAratio);
	  gRAASysratio->SetPointError(j, 5., yRAAsysRatio);
	}
	
	for (Int_t j=0; j<nBins; j++)
	{
	  Double_t yRAAratio = yRAA[2][j] / yRAA[1][j];
	  Double_t reRAA1 = eyRAA[1][j] / yRAA[1][j];
	  Double_t reRAA2 = eyRAA[2][j] / yRAA[2][j];
	  Double_t eraaRatio = TMath::Sqrt( reRAA1*reRAA1 + reRAA2*reRAA2 );
	  
	  Double_t renMB1= 1. / TMath::Sqrt(nMBevt[1][nBins-j]);
	  Double_t reTAA1= eTAAcent[1][nBins-j] / TAAcent[1][nBins-j]; 
	  Double_t yRAAsys1 = yRAA[1][j] * TMath::Sqrt( (yRAAsys[1][j]*yRAAsys[1][j]) / (yRAA[1][j]*yRAA[1][j]) - renMB1*renMB1 - reTAA1*reTAA1 );
	  Double_t renMB2= 1. / TMath::Sqrt(nMBevt[2][nBins-j]);
	  Double_t reTAA2= eTAAcent[2][nBins-j] / TAAcent[2][nBins-j]; 
	  Double_t yRAAsys2 = yRAA[2][j] * TMath::Sqrt( (yRAAsys[2][j]*yRAAsys[2][j]) / (yRAA[2][j]*yRAA[2][j]) - renMB2*renMB2 - reTAA2*reTAA2 );
	  Double_t yRAAsysRatio = TMath::Sqrt( yRAAsys1*yRAAsys1 + yRAAsys2*yRAAsys2 ) * yRAAratio;
	  
	  gRAAratio1->SetPoint(j, NPart[nBins-j], yRAAratio);
	  gRAAratio1->SetPointError(j, eNPart[nBins-j], eraaRatio);
	  gRAASysratio1->SetPoint(j, NPart[nBins-j], yRAAratio);
	  gRAASysratio1->SetPointError(j, 5., yRAAsysRatio);
	}
	
	for (Int_t j=0; j<nBins; j++)
	{
	  Double_t yRAAratio = yRAA[2][j] / yRAA[0][j];
	  Double_t reRAA1 = eyRAA[0][j] / yRAA[0][j];
	  Double_t reRAA2 = eyRAA[2][j] / yRAA[2][j];
	  Double_t eraaRatio = TMath::Sqrt( reRAA1*reRAA1 + reRAA2*reRAA2 );
	  
	  Double_t renMB1= 1. / TMath::Sqrt(nMBevt[0][nBins-j]);
	  Double_t reTAA1= eTAAcent[0][nBins-j] / TAAcent[0][nBins-j]; 
	  Double_t yRAAsys1 = yRAA[0][j] * TMath::Sqrt( (yRAAsys[0][j]*yRAAsys[0][j]) / (yRAA[0][j]*yRAA[0][j]) - renMB1*renMB1 - reTAA1*reTAA1 );
	  Double_t renMB2= 1. / TMath::Sqrt(nMBevt[2][nBins-j]);
	  Double_t reTAA2= eTAAcent[2][nBins-j] / TAAcent[2][nBins-j]; 
	  Double_t yRAAsys2 = yRAA[2][j] * TMath::Sqrt( (yRAAsys[2][j]*yRAAsys[2][j]) / (yRAA[2][j]*yRAA[2][j]) - renMB2*renMB2 - reTAA2*reTAA2 );
	  Double_t yRAAsysRatio = TMath::Sqrt( yRAAsys1*yRAAsys1 + yRAAsys2*yRAAsys2 ) * yRAAratio;
	  
	  gRAAratio2->SetPoint(j, NPart[nBins-j], yRAAratio);
	  gRAAratio2->SetPointError(j, eNPart[nBins-j], eraaRatio);
	  gRAASysratio2->SetPoint(j, NPart[nBins-j], yRAAratio);
	  gRAASysratio2->SetPointError(j, 5., yRAAsysRatio);
	}
	
	
	Double_t yRAAsysAllRatio = TMath::Sqrt( yRAAsysAll[0]*yRAAsysAll[0] - eppEff*eppEff - errfnorm*errfnorm - eppLumi*eppLumi + esigmaJPsippSyst[1]*esigmaJPsippSyst[1] );
	gRAASysAllratio->SetPoint(0,395.,1.);
	gRAASysAllratio->SetPointError(0,5.,yRAAsysAllRatio);
	Double_t yRAAsysAllRatio1 = TMath::Sqrt( yRAAsysAll[1]*yRAAsysAll[1] - eppEff*eppEff - errfnorm*errfnorm - eppLumi*eppLumi + esigmaJPsippSyst[2]*esigmaJPsippSyst[2] );
	gRAASysAllratio1->SetPoint(0,371.,1.);
	gRAASysAllratio1->SetPointError(0,5.,yRAAsysAllRatio1);
	Double_t yRAAsysAllRatio2 = TMath::Sqrt( yRAAsysAll[0]*yRAAsysAll[0] - eppEff*eppEff - errfnorm*errfnorm - eppLumi*eppLumi + esigmaJPsippSyst[2]*esigmaJPsippSyst[2] );
	gRAASysAllratio2->SetPoint(0,383.,1.);
	gRAASysAllratio2->SetPointError(0,5.,yRAAsysAllRatio2);
	
	
	// display our RAA results
	gRAAratio->GetXaxis()->Set(40, 0., 400.);
	gRAAratio->SetTitle("");
	
	// define axis
	gRAAratio->GetXaxis()->SetTitle("#LT N_{part} #GT");
	gRAAratio->GetXaxis()->SetTitleSize(0.05); 
	gRAAratio->GetXaxis()->SetTitleOffset(1.0);
	gRAAratio->GetXaxis()->SetLabelSize(0.05);
	
	gRAAratio->GetYaxis()->SetTitle("Ratio R_{AA}");
	gRAAratio->GetYaxis()->SetTitleSize(0.05);
	gRAAratio->GetYaxis()->SetTitleOffset(0.9);
	gRAAratio->GetYaxis()->SetLabelSize(0.05);
	gRAAratio->GetYaxis()->SetRangeUser(0.3,1.4);
	
	// define display
	gRAAratio->SetMarkerStyle(gMarkerStyleAliceForw);
	gRAAratio->SetMarkerSize(1.);
	gRAAratio->SetMarkerColor(kRed);
	gRAAratio->SetLineWidth(2);
	gRAAratio->SetLineColor(kRed);
	gRAASysratio->SetLineColor(kRed);
	gRAASysratio->SetFillStyle(0);
	gRAASysAllratio->SetLineWidth(1);
	gRAASysAllratio->SetFillStyle(3001);
	gRAASysAllratio->SetFillColor(kRed);
	
	gRAAratio1->SetMarkerStyle(gMarkerStyleAliceForw);
	gRAAratio1->SetMarkerSize(1.);
	gRAAratio1->SetMarkerColor(kBlue);
	gRAAratio1->SetLineWidth(2);
	gRAAratio1->SetLineColor(kBlue);
	gRAASysratio1->SetLineColor(kBlue);
	gRAASysratio1->SetFillStyle(0);
	gRAASysAllratio1->SetLineWidth(1);
	gRAASysAllratio1->SetFillStyle(3001);
	gRAASysAllratio1->SetFillColor(kBlue);
	
	gRAAratio2->SetMarkerStyle(gMarkerStyleAliceForw);
	gRAAratio2->SetMarkerSize(1.);
	gRAAratio2->SetMarkerColor(kGreen+2);
	gRAAratio2->SetLineWidth(2);
	gRAAratio2->SetLineColor(kGreen+2);
	gRAASysratio2->SetLineColor(kGreen+2);
	gRAASysratio2->SetFillStyle(0);
	gRAASysAllratio2->SetLineWidth(1);
	gRAASysAllratio2->SetFillStyle(3001);
	gRAASysAllratio2->SetFillColor(kGreen+2);
	
	
	TCanvas *c2 = new TCanvas("RAA_ptAll_y0_Ratio","RAA -- ptAll -- y0 -- Ratio",gCanvasWidth,gCanvasHeigh);
	gPad->SetTopMargin(0.05);
	gPad->SetBottomMargin(0.12);
	gPad->SetRightMargin(0.05);
	
	gRAAratio->Draw("ap");
    gRAASysratio->Draw("e2");
    gRAASysAllratio->Draw("e2");
	gRAAratio1->Draw("p");
    gRAASysratio1->Draw("e2");
    gRAASysAllratio1->Draw("e2");
	gRAAratio2->Draw("p");
    gRAASysratio2->Draw("e2");
    gRAASysAllratio2->Draw("e2");
	
	TLine *l =  new TLine(0., 1., 400., 1.);
	l->SetLineStyle(3);
	l->Draw();
	
	TLegend *lg = new TLegend(0.55, 0.73, 0.9, 0.93,"","NDC");
	lg->AddEntry(gRAAratio," 2<pt<5 / 0<pt<2 ","p");
	lg->AddEntry(gRAAratio2," 5<pt<8 / 0<pt<2 ","p");
	lg->AddEntry(gRAAratio1," 5<pt<8 / 2<pt<5 ","p");
	lg->SetFillStyle(0);
	lg->SetBorderSize(0);
	lg->SetTextFont(42);
	lg->SetMargin(0.05);
	lg->Draw();
	
	c2->SaveAs("Outputs/Raa_ptAll_y0_Ratio.pdf");
	
  }
  
  /*
  else if (var == "ptBins") 
  {
	Double_t *x, *ex, *bin, *NPart, *eNPart, *dNchdEta, *edNchdEta, *dNchdEtaNorm, *edNchdEtaNorm;;
	TString *label;
    SetGraphLabelsVsCent(nBins+1, x, ex, bin, label, NPart, eNPart,dNchdEta, edNchdEta, dNchdEtaNorm, edNchdEtaNorm);
	
	
	Double_t **yRAA = new Double_t*[3];
	Double_t **eyRAA = new Double_t*[3];
	Double_t **yRAAsys = new Double_t*[3];
	Double_t *yRAAsysAll = new Double_t[3];
	
	Double_t **nMBevt = new Double_t*[3];
	Double_t **TAAcent = new Double_t*[3];
	Double_t **eTAAcent = new Double_t*[3];
	Double_t sigmaJPsipp[3], esigmaJPsippSyst[3];
	
	
	// Get values of Raa
	for (Int_t i=5; i<8; i++) 
	{
	  yRAA[i-5] = new Double_t[nBins];
	  yRAA[i-5] = gRAA[i]->GetY();
	  eyRAA[i-5] = new Double_t[nBins];
	  eyRAA[i-5] = gRAA[i]->GetEY();
	  yRAAsys[i-5] = new Double_t[nBins];
	  yRAAsys[i-5] = gRAASys[i]->GetEY();
	  yRAAsysAll[i-5] = gRAASysAll[i]->GetErrorY(0);
	  
	  // get normalization variables in nBins centrality bins
	  Double_t *nMB, *TAA, *eTAA, *AccEff, *eAccEff;
	  SetNormVarVsCent(i, 0, nBins+1, nMB, TAA, eTAA, AccEff, eAccEff);
	  
	  nMBevt[i-5] = new Double_t[nBins+1];
	  nMBevt[i-5] = nMB;
	  TAAcent[i-5] = new Double_t[nBins+1];
	  TAAcent[i-5] = TAA;
	  eTAAcent[i-5] = new Double_t[nBins+1];
	  eTAAcent[i-5] = eTAA;
	  
	  // get centrality independent normalization variables
	  SetNormVar(i, 0, sigmaJPsipp[i-1], esigmaJPsippSyst[i-1]);
	  
	}
	
	
	TGraphErrors *gRAAratio = new TGraphErrors(nBins);
	TGraphErrors *gRAASysratio = new TGraphErrors(nBins);
	TGraphErrors *gRAASysAllratio = new TGraphErrors(1);
	
	TGraphErrors *gRAAratio1 = new TGraphErrors(nBins);
	TGraphErrors *gRAASysratio1 = new TGraphErrors(nBins);
	TGraphErrors *gRAASysAllratio1 = new TGraphErrors(1);
	
	TGraphErrors *gRAAratio2 = new TGraphErrors(nBins);
	TGraphErrors *gRAASysratio2 = new TGraphErrors(nBins);
	TGraphErrors *gRAASysAllratio2 = new TGraphErrors(1);
	
	
	// Computation and Build graph
	for (Int_t j=0; j<nBins; j++)
	{
	  Double_t yRAAratio = yRAA[1][j] / yRAA[0][j];
	  Double_t reRAA1 = eyRAA[0][j] / yRAA[0][j];
	  Double_t reRAA2 = eyRAA[1][j] / yRAA[1][j];
	  Double_t eraaRatio = TMath::Sqrt( reRAA1*reRAA1 + reRAA2*reRAA2 );
	  
	  Double_t renMB1= 1. / TMath::Sqrt(nMBevt[0][nBins-j]);
	  Double_t reTAA1= eTAAcent[0][nBins-j] / TAAcent[0][nBins-j]; 
	  Double_t yRAAsys1 = yRAA[0][j] * TMath::Sqrt( (yRAAsys[0][j]*yRAAsys[0][j]) / (yRAA[0][j]*yRAA[0][j]) - renMB1*renMB1 - reTAA1*reTAA1 );
	  Double_t renMB2= 1. / TMath::Sqrt(nMBevt[1][nBins-j]);
	  Double_t reTAA2= eTAAcent[1][nBins-j] / TAAcent[1][nBins-j]; 
	  Double_t yRAAsys2 = yRAA[1][j] * TMath::Sqrt( (yRAAsys[1][j]*yRAAsys[1][j]) / (yRAA[1][j]*yRAA[1][j]) - renMB2*renMB2 - reTAA2*reTAA2 );
	  Double_t yRAAsysRatio = TMath::Sqrt( yRAAsys1*yRAAsys1 + yRAAsys2*yRAAsys2 ) * yRAAratio;
	  
	  gRAAratio->SetPoint(j, NPart[nBins-j], yRAAratio);
	  gRAAratio->SetPointError(j, eNPart[nBins-j], eraaRatio);
	  gRAASysratio->SetPoint(j, NPart[nBins-j], yRAAratio);
	  gRAASysratio->SetPointError(j, 5., yRAAsysRatio);
	}
	
	for (Int_t j=0; j<nBins; j++)
	{
	  Double_t yRAAratio = yRAA[2][j] / yRAA[1][j];
	  Double_t reRAA1 = eyRAA[1][j] / yRAA[1][j];
	  Double_t reRAA2 = eyRAA[2][j] / yRAA[2][j];
	  Double_t eraaRatio = TMath::Sqrt( reRAA1*reRAA1 + reRAA2*reRAA2 );
	  
	  Double_t renMB1= 1. / TMath::Sqrt(nMBevt[1][nBins-j]);
	  Double_t reTAA1= eTAAcent[1][nBins-j] / TAAcent[1][nBins-j]; 
	  Double_t yRAAsys1 = yRAA[1][j] * TMath::Sqrt( (yRAAsys[1][j]*yRAAsys[1][j]) / (yRAA[1][j]*yRAA[1][j]) - renMB1*renMB1 - reTAA1*reTAA1 );
	  Double_t renMB2= 1. / TMath::Sqrt(nMBevt[2][nBins-j]);
	  Double_t reTAA2= eTAAcent[2][nBins-j] / TAAcent[2][nBins-j]; 
	  Double_t yRAAsys2 = yRAA[2][j] * TMath::Sqrt( (yRAAsys[2][j]*yRAAsys[2][j]) / (yRAA[2][j]*yRAA[2][j]) - renMB2*renMB2 - reTAA2*reTAA2 );
	  Double_t yRAAsysRatio = TMath::Sqrt( yRAAsys1*yRAAsys1 + yRAAsys2*yRAAsys2 ) * yRAAratio;
	  
	  gRAAratio1->SetPoint(j, NPart[nBins-j], yRAAratio);
	  gRAAratio1->SetPointError(j, eNPart[nBins-j], eraaRatio);
	  gRAASysratio1->SetPoint(j, NPart[nBins-j], yRAAratio);
	  gRAASysratio1->SetPointError(j, 5., yRAAsysRatio);
	}
	
	for (Int_t j=0; j<nBins; j++)
	{
	  Double_t yRAAratio = yRAA[2][j] / yRAA[0][j];
	  Double_t reRAA1 = eyRAA[0][j] / yRAA[0][j];
	  Double_t reRAA2 = eyRAA[2][j] / yRAA[2][j];
	  Double_t eraaRatio = TMath::Sqrt( reRAA1*reRAA1 + reRAA2*reRAA2 );
	  
	  Double_t renMB1= 1. / TMath::Sqrt(nMBevt[0][nBins-j]);
	  Double_t reTAA1= eTAAcent[0][nBins-j] / TAAcent[0][nBins-j]; 
	  Double_t yRAAsys1 = yRAA[0][j] * TMath::Sqrt( (yRAAsys[0][j]*yRAAsys[0][j]) / (yRAA[0][j]*yRAA[0][j]) - renMB1*renMB1 - reTAA1*reTAA1 );
	  Double_t renMB2= 1. / TMath::Sqrt(nMBevt[2][nBins-j]);
	  Double_t reTAA2= eTAAcent[2][nBins-j] / TAAcent[2][nBins-j]; 
	  Double_t yRAAsys2 = yRAA[2][j] * TMath::Sqrt( (yRAAsys[2][j]*yRAAsys[2][j]) / (yRAA[2][j]*yRAA[2][j]) - renMB2*renMB2 - reTAA2*reTAA2 );
	  Double_t yRAAsysRatio = TMath::Sqrt( yRAAsys1*yRAAsys1 + yRAAsys2*yRAAsys2 ) * yRAAratio;
	  
	  gRAAratio2->SetPoint(j, NPart[nBins-j], yRAAratio);
	  gRAAratio2->SetPointError(j, eNPart[nBins-j], eraaRatio);
	  gRAASysratio2->SetPoint(j, NPart[nBins-j], yRAAratio);
	  gRAASysratio2->SetPointError(j, 5., yRAAsysRatio);
	}
	
	
	Double_t yRAAsysAllRatio = TMath::Sqrt( yRAAsysAll[0]*yRAAsysAll[0] - eppEff*eppEff - errfnorm*errfnorm - eppLumi*eppLumi + esigmaJPsippSyst[1]*esigmaJPsippSyst[1] );
	gRAASysAllratio->SetPoint(0,395.,1.);
	gRAASysAllratio->SetPointError(0,5.,yRAAsysAllRatio);
	Double_t yRAAsysAllRatio1 = TMath::Sqrt( yRAAsysAll[1]*yRAAsysAll[1] - eppEff*eppEff - errfnorm*errfnorm - eppLumi*eppLumi + esigmaJPsippSyst[2]*esigmaJPsippSyst[2] );
	gRAASysAllratio1->SetPoint(0,383.,1.);
	gRAASysAllratio1->SetPointError(0,5.,yRAAsysAllRatio1);
	Double_t yRAAsysAllRatio2 = TMath::Sqrt( yRAAsysAll[0]*yRAAsysAll[0] - eppEff*eppEff - errfnorm*errfnorm - eppLumi*eppLumi + esigmaJPsippSyst[2]*esigmaJPsippSyst[2] );
	gRAASysAllratio2->SetPoint(0,371.,1.);
	gRAASysAllratio2->SetPointError(0,5.,yRAAsysAllRatio2);
	
	
	// display our RAA results
	gRAAratio->GetXaxis()->Set(40, 0., 400.);
	gRAAratio->SetTitle("");
	
	// define axis
	gRAAratio->GetXaxis()->SetTitle("#LT N_{part} #GT");
	gRAAratio->GetXaxis()->SetTitleSize(0.05); 
	gRAAratio->GetXaxis()->SetTitleOffset(1.0);
	gRAAratio->GetXaxis()->SetLabelSize(0.05);
	
	gRAAratio->GetYaxis()->SetTitle("Ratio R_{AA}");
	gRAAratio->GetYaxis()->SetTitleSize(0.05);
	gRAAratio->GetYaxis()->SetTitleOffset(0.9);
	gRAAratio->GetYaxis()->SetLabelSize(0.05);
	gRAAratio->GetYaxis()->SetRangeUser(0.5,1.25);
	
	// define display
	gRAAratio->SetMarkerStyle(gMarkerStyleAliceForw);
	gRAAratio->SetMarkerSize(1.);
	gRAAratio->SetMarkerColor(kRed);
	gRAAratio->SetLineWidth(2);
	gRAAratio->SetLineColor(kRed);
	gRAASysratio->SetLineColor(kRed);
	gRAASysratio->SetFillStyle(0);
	gRAASysAllratio->SetLineWidth(1);
	gRAASysAllratio->SetFillStyle(3001);
	gRAASysAllratio->SetFillColor(kRed);
	
	gRAAratio1->SetMarkerStyle(gMarkerStyleAliceForw);
	gRAAratio1->SetMarkerSize(1.);
	gRAAratio1->SetMarkerColor(kBlue);
	gRAAratio1->SetLineWidth(2);
	gRAAratio1->SetLineColor(kBlue);
	gRAASysratio1->SetLineColor(kBlue);
	gRAASysratio1->SetFillStyle(0);
	gRAASysAllratio1->SetLineWidth(1);
	gRAASysAllratio1->SetFillStyle(3001);
	gRAASysAllratio1->SetFillColor(kBlue);
	
	gRAAratio2->SetMarkerStyle(gMarkerStyleAliceForw);
	gRAAratio2->SetMarkerSize(1.);
	gRAAratio2->SetMarkerColor(kGreen+2);
	gRAAratio2->SetLineWidth(2);
	gRAAratio2->SetLineColor(kGreen+2);
	gRAASysratio2->SetLineColor(kGreen+2);
	gRAASysratio2->SetFillStyle(0);
	gRAASysAllratio2->SetLineWidth(1);
	gRAASysAllratio2->SetFillStyle(3001);
	gRAASysAllratio2->SetFillColor(kGreen+2);
	
	
	TCanvas *c2 = new TCanvas("RAA_centAll_y0","RAA -- centAll -- y0",gCanvasWidth,gCanvasHeigh);
	gPad->SetTopMargin(0.05);
	gPad->SetBottomMargin(0.12);
	gPad->SetRightMargin(0.05);
	
	gRAAratio->Draw("ap");
    gRAASysratio->Draw("e2");
    gRAASysAllratio->Draw("e2");
	gRAAratio1->Draw("p");
    gRAASysratio1->Draw("e2");
    gRAASysAllratio1->Draw("e2");
	gRAAratio2->Draw("p");
    gRAASysratio2->Draw("e2");
    gRAASysAllratio2->Draw("e2");
	
	TLine *l =  new TLine(0., 1., 400., 1.);
	l->SetLineStyle(3);
	l->Draw();
	
	TLegend *lg = new TLegend(0.17, 0.75, 0.92, 0.97,"","NDC");
	lg->AddEntry(gRAAratio," 2<pt<5 / 0<pt<2 ","p");
	lg->AddEntry(gRAAratio," 5<pt<8 / 2<pt<5 ","p");
	lg->AddEntry(gRAAratio," 5<pt<8 / 0<pt<2 ","p");
	lg->SetFillStyle(0);
	lg->SetBorderSize(0);
	lg->SetTextFont(42);
	lg->SetMargin(0.05);
	lg->Draw();
	
	c2->SaveAs("Outputs/Raa_centAll_y0_Ratio.pdf");

  }
   */
  
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
  if (y == 0 && pt == 0) lg->AddEntry(gRAA,"ALICE 2011 (Pb-Pb #sqrt{s_{NN}} = 2.76 TeV), 2.5<y<4, p_{T}>0 ","p");
  else if (y == 1 && pt == 0) lg->AddEntry(gRAA,"ALICE (Pb-Pb #sqrt{s_{NN}} = 2.76 TeV), 3.25<y<4, p_{T}>0 ","p");
  else if (y == 2 && pt == 0) lg->AddEntry(gRAA,"ALICE (Pb-Pb #sqrt{s_{NN}} = 2.76 TeV), 2.5<y<3.25, p_{T}>0 ","p");
  else if (y == 0 && pt == 1) lg->AddEntry(gRAA,"ALICE (Pb-Pb #sqrt{s_{NN}} = 2.76 TeV), 2.5<y<4, 0<p_{T}<3 GeV/c ","p");
  else if (y == 0 && pt == 2) lg->AddEntry(gRAA,"ALICE (Pb-Pb #sqrt{s_{NN}} = 2.76 TeV), 2.5<y<4, p_{T}>3 GeV/c ","p");
  if (PRL) lg->AddEntry(gRAA_PHENIX,"ALICE 2010 (Pb-Pb #sqrt{s_{NN}} = 2.76 TeV), 2.5<y<4, p_{T}>0 ","p");
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
  
  gRAA_PHENIX->SetTitle("");
  gRAA_PHENIXSys->SetTitle("");
  gRAA_PHENIXSysAll->SetTitle("");
  gRAA->SetTitle("");
  gRAASys->SetTitle("");
  gRAASysAll->SetTitle("");
  
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
  TASImage *myAliceLogo = new TASImage("alice_logo.png");
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
  
  TDatime dt;
  TString today = Form("%02i/%02i/%4i", dt.GetDay(), dt.GetMonth(), dt.GetYear());
  t2->AddText(0.,0.,today.Data());
  t2->Draw();
}

