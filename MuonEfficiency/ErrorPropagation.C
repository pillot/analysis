/*
 *  ErrorPropagation.C
 *  aliroot_dev
 *
 *  Created by philippe pillot on 18/11/13.
 *  Copyright 2013 Subatech. All rights reserved.
 *
 */

#if !defined(__CINT__) || defined(__MAKECINT__)

#include <TMath.h>
#include <TH1D.h>
#include <TRandom2.h>
#include <TCanvas.h>
#include <TEfficiency.h>

#endif

const Int_t knMode = 6;
const Int_t kRefMode = 1; // can go from 0 to 2+knMode
TString sigmaModeMeas[2] = {"variance:      ", "level/2 tails: "};
TString sigmaModeCalc[knMode] = {"variance:      ", "bayesian:      ", "bayesian short:", "Wilson:        ", "ClopperPearson:", "FeldmanCousins:"};

void PrintError(Double_t eff, Double_t sigmaMeas[2][3], Double_t sigmaCalc[knMode][3]);

//----------------------------------------------------------------------------
void ErrorPropagation()
{
  
  Bool_t bayesian = kTRUE; // bayesian estimator = (k+a)/(n+a+b). Otherwise use frequentist estimator = k/n
  Double_t level = 0.682689492137; // 1-sigma confidence level
  Double_t a = 1.; // prior for k
  Double_t b = 1.; // prior for n-k
  Double_t quantiles[2];
  quantiles[0] = 0.5*(1.-level);
  quantiles[1] = 0.5*(1.+level);
  Int_t n = 10000000;
  
  // measurements
//  Double_t tot[10] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
//  Double_t pass[10] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
  Double_t tot[10] = {1., 1., 1., 1., 1., 1., 1., 1., 1., 1.};
  Double_t pass[10] = {1., 1., 0., 0., 1., 0., 1., 0., 0., 1.};
//  Double_t tot[10] = {10., 10., 10., 10., 10., 10., 10., 10., 10., 10.};
//  Double_t pass[10] = {10., 10., 0., 0., 10., 0., 10., 0., 0., 10.};
//  Double_t tot[10] = {100., 100., 100., 100., 100., 100., 100., 100., 100., 100.};
//  Double_t pass[10] = {100., 100., 0., 0., 100., 0., 100., 0., 0., 100.};
//  Double_t tot[10] = {1000., 1000., 1000., 1000., 1000., 1000., 1000., 1000., 1000., 1000.};
//  Double_t pass[10] = {800., 900., 900., 700., 600., 600., 800., 700., 800., 900.};
  
  // ranges for simulation
  Double_t range[10][2];
  for (Int_t iCh = 0; iCh < 10; iCh++) {
    range[iCh][0] = TMath::Max(0., TEfficiency::Bayesian((Int_t)tot[iCh], (Int_t)pass[iCh], 0.99, a, b, kFALSE));
    range[iCh][1] = TMath::Min(1., TEfficiency::Bayesian((Int_t)tot[iCh], (Int_t)pass[iCh], 0.99, a, b, kTRUE));
    printf("simulation range for ch%d = [%f, %f]\n", iCh+1, range[iCh][0], range[iCh][1]);
  }
  
  // efficiency per chamber
  Double_t effCh[10], varCh[10][2], nCh[10][2];
  for (Int_t iCh = 0; iCh < 10; iCh++) {
    if (bayesian) effCh[iCh] = (pass[iCh]+a)/(tot[iCh]+a+b);
    else effCh[iCh] = (tot[iCh] > 0.) ? pass[iCh]/tot[iCh] : 1.;
    for (Int_t i = 0; i < 2; i++) {
      varCh[iCh][i] = 0.;
      nCh[iCh][i] = 0.;
    }
  }
  
  TH1D *hEffCh[10];
  for (Int_t iCh = 0; iCh < 10; iCh++) {
    hEffCh[iCh] = new TH1D(Form("hEffCh%d",iCh+1), Form("chamber %d efficiency",iCh+1), 10000, 0., 1.);
    hEffCh[iCh]->Sumw2();
  }
  
  // efficiency per station
  Double_t effSt[6], varSt[6][2], nSt[6][2];
  for (Int_t iSt = 0; iSt < 6; iSt++) {
    if (iSt < 5) effSt[iSt] = 1.-(1.-effCh[2*iSt])*(1.-effCh[2*iSt+1]);
    else effSt[5] = effCh[6]*effCh[7]*effCh[8] + effCh[6]*effCh[7]*effCh[9] + effCh[6]*effCh[8]*effCh[9] + effCh[7]*effCh[8]*effCh[9] - 3.*effCh[6]*effCh[7]*effCh[8]*effCh[9];
    for (Int_t i = 0; i < 2; i++) {
      varSt[iSt][i] = 0.;
      nSt[iSt][i] = 0.;
    }
    
  }
  
  TH1D *hEffSt[6];
  for (Int_t iSt = 0; iSt < 6; iSt++) {
    if (iSt < 5) hEffSt[iSt] = new TH1D(Form("hEffSt%d",iSt+1), Form("station %d efficiency",iSt+1), 10000, 0., 1.);
    else hEffSt[iSt] = new TH1D("hEffSt45", "station 4-5 efficiency", 10000, 0., 1.);
    hEffSt[iSt]->Sumw2();
  }
  
  // spectrometer efficiency
  Double_t effSpectro = effSt[0] * effSt[1] * effSt[2] * effSt[3] * effSt[4];
  Double_t varSpectro[2] = {0., 0.};
  Double_t nSpectro[2] = {0., 0.};
  
  TH1D *hEffSpectro = new TH1D("hEffSpectro", "spectrometer efficiency", 10000, 0., 1.);
  hEffSpectro->Sumw2();
  
  Double_t effSpectro45 = effSt[0] * effSt[1] * effSt[2] * effSt[5];
  Double_t varSpectro45[2] = {0., 0.};
  Double_t nSpectro45[2] = {0., 0.};
  
  TH1D *hEffSpectro45 = new TH1D("hEffSpectro45", "spectrometer efficiency (3/4ch in st45)", 10000, 0., 1.);
  hEffSpectro45->Sumw2();
  
  gRandom->SetSeed(0);
  
  for (Int_t iEv = 0; iEv < n; iEv++) {
    
    // efficiency distribution per chamber
    Double_t eCh[10], wCh[10];
    for (Int_t iCh = 0; iCh < 10; iCh++) {
      
      eCh[iCh] = gRandom->Uniform(range[iCh][0], range[iCh][1]);
      wCh[iCh] = TMath::BetaDist(eCh[iCh], pass[iCh]+a, tot[iCh]-pass[iCh]+b);
      
      hEffCh[iCh]->Fill(eCh[iCh], wCh[iCh]);
      
      Int_t i = (eCh[iCh] < effCh[iCh]) ? 0 : 1;
      varCh[iCh][i] += wCh[iCh]*(eCh[iCh]-effCh[iCh])*(eCh[iCh]-effCh[iCh]);
      nCh[iCh][i] += wCh[iCh];
      
    }
    
    // efficiency distribution per station
    Double_t eSt[10], wSt[10];
    for (Int_t iSt = 0; iSt < 6; iSt++) {
      
      if (iSt < 5) {
	
	eSt[iSt] = 1.-(1.-eCh[2*iSt])*(1.-eCh[2*iSt+1]);
	wSt[iSt] = wCh[2*iSt]*wCh[2*iSt+1];
	
      } else {
	
	eSt[5] = eCh[6]*eCh[7]*eCh[8] + eCh[6]*eCh[7]*eCh[9] + eCh[6]*eCh[8]*eCh[9] + eCh[7]*eCh[8]*eCh[9] - 3.*eCh[6]*eCh[7]*eCh[8]*eCh[9];
	wSt[5] = wCh[6]*wCh[7]*wCh[8]*wCh[9];
	
      }
      
      hEffSt[iSt]->Fill(eSt[iSt], wSt[iSt]);
      
      Int_t i = (eSt[iSt] < effSt[iSt]) ? 0 : 1;
      varSt[iSt][i] += wSt[iSt]*(eSt[iSt]-effSt[iSt])*(eSt[iSt]-effSt[iSt]);
      nSt[iSt][i] += wSt[iSt];
      
    }
    
    // spectrometer efficiency
    Double_t e = eSt[0]*eSt[1]*eSt[2]*eSt[3]*eSt[4];
    Double_t w = wSt[0]*wSt[1]*wSt[2]*wSt[3]*wSt[4];
    
    hEffSpectro->Fill(e, w);
    
    Int_t i = (e < effSpectro) ? 0 : 1;
    varSpectro[i] += w*(e-effSpectro)*(e-effSpectro);
    nSpectro[i] += w;
    
    // spectrometer efficiency 3/4ch in st45
    e = eSt[0]*eSt[1]*eSt[2]*eSt[5];
    w = wSt[0]*wSt[1]*wSt[2]*wSt[5];
    
    hEffSpectro45->Fill(e, w);
    
    i = (e < effSpectro45) ? 0 : 1;
    varSpectro45[i] += w*(e-effSpectro45)*(e-effSpectro45);
    nSpectro45[i] += w;
    
  }
  
  // efficiency per chamber
  Double_t sigmaChMeas[10][2][3];
  Double_t sigmaChCalc[10][knMode][3];
  for (Int_t iCh = 0; iCh < 10; iCh++) {
    
    sigmaChMeas[iCh][0][0] = (nCh[iCh][0] > 0) ? sqrt(varCh[iCh][0] / nCh[iCh][0]) : 0.;
    sigmaChMeas[iCh][0][1] = (nCh[iCh][1] > 0) ? sqrt(varCh[iCh][1] / nCh[iCh][1]) : 0.;
    sigmaChMeas[iCh][0][2] = (nCh[iCh][0]+nCh[iCh][1] > 0) ? sqrt((varCh[iCh][0]+varCh[iCh][1]) / (nCh[iCh][0]+nCh[iCh][1])) : 0.;
    
    hEffCh[iCh]->GetQuantiles(2, sigmaChMeas[iCh][1], quantiles);
    sigmaChMeas[iCh][1][0] = TMath::Max(0., effCh[iCh] - sigmaChMeas[iCh][1][0]);
    sigmaChMeas[iCh][1][1] = TMath::Max(0., sigmaChMeas[iCh][1][1] - effCh[iCh]);
    
    sigmaChCalc[iCh][0][0] = sigmaChMeas[iCh][0][0];
    sigmaChCalc[iCh][0][1] = sigmaChMeas[iCh][0][1];
    sigmaChCalc[iCh][0][2] = sigmaChMeas[iCh][0][2];
    
    sigmaChCalc[iCh][1][0] = TMath::Max(0., effCh[iCh] - TEfficiency::Bayesian((Int_t)tot[iCh], (Int_t)pass[iCh], level, a, b, kFALSE));
    sigmaChCalc[iCh][1][1] = TMath::Max(0., TEfficiency::Bayesian((Int_t)tot[iCh], (Int_t)pass[iCh], level, a, b, kTRUE) - effCh[iCh]);
    
    sigmaChCalc[iCh][2][0] = effCh[iCh] - TEfficiency::Bayesian((Int_t)tot[iCh], (Int_t)pass[iCh], level, a, b, kFALSE, kTRUE);
    sigmaChCalc[iCh][2][1] = TEfficiency::Bayesian((Int_t)tot[iCh], (Int_t)pass[iCh], level, a, b, kTRUE, kTRUE) - effCh[iCh];
    
    sigmaChCalc[iCh][3][0] = effCh[iCh] - TEfficiency::Wilson((Int_t)tot[iCh], (Int_t)pass[iCh], level, kFALSE);
    sigmaChCalc[iCh][3][1] = TEfficiency::Wilson((Int_t)tot[iCh], (Int_t)pass[iCh], level, kTRUE) - effCh[iCh];
    
    sigmaChCalc[iCh][4][0] = effCh[iCh] - TEfficiency::ClopperPearson((Int_t)tot[iCh], (Int_t)pass[iCh], level, kFALSE);
    sigmaChCalc[iCh][4][1] = TEfficiency::ClopperPearson((Int_t)tot[iCh], (Int_t)pass[iCh], level, kTRUE) - effCh[iCh];
    
    sigmaChCalc[iCh][5][0] = effCh[iCh] - TEfficiency::FeldmanCousins((Int_t)tot[iCh], (Int_t)pass[iCh], level, kFALSE);
    sigmaChCalc[iCh][5][1] = TEfficiency::FeldmanCousins((Int_t)tot[iCh], (Int_t)pass[iCh], level, kTRUE) - effCh[iCh];
    
  }
  
  // print results per chamber
  for (Int_t iCh = 0; iCh < 10; iCh++) {
    printf("\nefficiency of chamber %d:\n", iCh+1);
    PrintError(effCh[iCh], sigmaChMeas[iCh], sigmaChCalc[iCh]);
  }
  
  // efficiency per station
  Double_t sigmaStMeas[6][2][3];
  Double_t sigmaStCalc[6][knMode][3];
  for (Int_t iSt = 0; iSt < 6; iSt++) {
    
    sigmaStMeas[iSt][0][0] = (nSt[iSt][0] > 0) ? sqrt(varSt[iSt][0] / nSt[iSt][0]) : 0.;
    sigmaStMeas[iSt][0][1] = (nSt[iSt][1] > 0) ? sqrt(varSt[iSt][1] / nSt[iSt][1]) : 0.;
    sigmaStMeas[iSt][0][2] = (nSt[iSt][0]+nSt[iSt][1] > 0) ? sqrt((varSt[iSt][0]+varSt[iSt][1]) / (nSt[iSt][0]+nSt[iSt][1])) : 0.;
    
    hEffSt[iSt]->GetQuantiles(2, sigmaStMeas[iSt][1], quantiles);
    sigmaStMeas[iSt][1][0] = TMath::Max(0., effSt[iSt] - sigmaStMeas[iSt][1][0]);
    sigmaStMeas[iSt][1][1] = TMath::Max(0., sigmaStMeas[iSt][1][1] - effSt[iSt]);
    
    Double_t d1, d2, d3, d4, d12, d13, d14, d23, d24, d34, d123, d124, d134, d234, d1234;
    if (iSt < 5) {
      
      d1 = (1. - effCh[2*iSt+1]); d1 *= d1;
      d2 = (1. - effCh[2*iSt]); d2 *= d2;
      
    } else {
      
      d1 = effCh[7]*effCh[8] + effCh[7]*effCh[9] + effCh[8]*effCh[9] - 3.*effCh[7]*effCh[8]*effCh[9]; d1 *= d1;
      d2 = effCh[6]*effCh[8] + effCh[6]*effCh[9] + effCh[8]*effCh[9] - 3.*effCh[6]*effCh[8]*effCh[9]; d2 *= d2;
      d3 = effCh[6]*effCh[7] + effCh[6]*effCh[9] + effCh[7]*effCh[9] - 3.*effCh[6]*effCh[7]*effCh[9]; d3 *= d3;
      d4 = effCh[6]*effCh[7] + effCh[6]*effCh[8] + effCh[7]*effCh[8] - 3.*effCh[6]*effCh[7]*effCh[8]; d4 *= d4;
      d12 = effCh[8] + effCh[9] - 3.*effCh[8]*effCh[9]; d12 *= d12;
      d13 = effCh[7] + effCh[9] - 3.*effCh[7]*effCh[9]; d13 *= d13;
      d14 = effCh[7] + effCh[8] - 3.*effCh[7]*effCh[8]; d14 *= d14;
      d23 = effCh[6] + effCh[9] - 3.*effCh[6]*effCh[9]; d23 *= d23;
      d24 = effCh[6] + effCh[8] - 3.*effCh[6]*effCh[8]; d24 *= d24;
      d34 = effCh[6] + effCh[7] - 3.*effCh[6]*effCh[7]; d34 *= d34;
      d123 = 1. - 3.*effCh[9]; d123 *= d123;
      d124 = 1. - 3.*effCh[8]; d124 *= d124;
      d134 = 1. - 3.*effCh[7]; d134 *= d134;
      d234 = 1. - 3.*effCh[6]; d234 *= d234;
      d1234 = - 3.; d1234 *= d1234;
      
    }
    
    for (Int_t iMode = 0; iMode < knMode; iMode++) {
      
      for (Int_t i = 0; (iMode == 0 && i < 3) || i < 2; i++) {
	
	if (iSt < 5) {
	  
	  Double_t s1 = sigmaChCalc[2*iSt][iMode][i] * sigmaChCalc[2*iSt][iMode][i];
	  Double_t s2 = sigmaChCalc[2*iSt+1][iMode][i] * sigmaChCalc[2*iSt+1][iMode][i];
	  sigmaStCalc[iSt][iMode][i] = TMath::Sqrt(d1*s1 + d2*s2 + s1*s2);
	  
	} else {
	  
	  Double_t s1 = sigmaChCalc[6][iMode][i] * sigmaChCalc[6][iMode][i];
	  Double_t s2 = sigmaChCalc[7][iMode][i] * sigmaChCalc[7][iMode][i];
	  Double_t s3 = sigmaChCalc[8][iMode][i] * sigmaChCalc[8][iMode][i];
	  Double_t s4 = sigmaChCalc[9][iMode][i] * sigmaChCalc[9][iMode][i];
	  sigmaStCalc[iSt][iMode][i] = TMath::Sqrt(d1*s1 + d2*s2 + d3*s3 + d4*s4
                                                   + d12*s1*s2 + d13*s1*s3 + d14*s1*s4 + d23*s2*s3 + d24*s2*s4 + d34*s3*s4
                                                   + d123*s1*s2*s3 + d124*s1*s2*s4 + d134*s1*s3*s4 + d234*s2*s3*s4
                                                   + d1234*s1*s2*s3*s4);
	  
	}
	
      }
      
    }
    
  }
  
  // print results per station
  for (Int_t iSt = 0; iSt < 6; iSt++) {
    if (iSt < 5) printf("\nefficiency of station %d:\n", iSt+1);
    else  printf("\nefficiency of station 4-5:\n");
    PrintError(effSt[iSt], sigmaStMeas[iSt], sigmaStCalc[iSt]);
  }
  
  // spectrometer efficiency
  Double_t sigmaSpectroMeas[2][3];
  Double_t sigmaSpectroCalc[knMode][3];
  Double_t sigmaSpectro45Meas[2][3];
  Double_t sigmaSpectro45Calc[knMode][3];
  
  sigmaSpectroMeas[0][0] = (nSpectro[0] > 0) ? sqrt(varSpectro[0] / nSpectro[0]) : 0.;
  sigmaSpectroMeas[0][1] = (nSpectro[1] > 0) ? sqrt(varSpectro[1] / nSpectro[1]) : 0.;
  sigmaSpectroMeas[0][2] = sqrt((varSpectro[0]+varSpectro[1]) / (nSpectro[0]+nSpectro[1]));
  
  sigmaSpectro45Meas[0][0] = (nSpectro45[0] > 0) ? sqrt(varSpectro45[0] / nSpectro45[0]) : 0.;
  sigmaSpectro45Meas[0][1] = (nSpectro45[1] > 0) ? sqrt(varSpectro45[1] / nSpectro45[1]) : 0.;
  sigmaSpectro45Meas[0][2] = sqrt((varSpectro45[0]+varSpectro45[1]) / (nSpectro45[0]+nSpectro45[1]));
  
  hEffSpectro->GetQuantiles(2, sigmaSpectroMeas[1], quantiles);
  sigmaSpectroMeas[1][0] = TMath::Max(0., effSpectro - sigmaSpectroMeas[1][0]);
  sigmaSpectroMeas[1][1] = TMath::Max(0., sigmaSpectroMeas[1][1] - effSpectro);
  
  hEffSpectro45->GetQuantiles(2, sigmaSpectro45Meas[1], quantiles);
  sigmaSpectro45Meas[1][0] = TMath::Max(0., effSpectro45 - sigmaSpectro45Meas[1][0]);
  sigmaSpectro45Meas[1][1] = TMath::Max(0., sigmaSpectro45Meas[1][1] - effSpectro45);
  
  Double_t de[6][2];
  for (Int_t iSt = 0; iSt < 6; iSt++) de[iSt][0] = effSt[iSt]*effSt[iSt];
  
  for (Int_t iMode = 0; iMode < knMode; iMode++) {
    
    for (Int_t i = 0; (iMode == 0 && i < 3) || i < 2; i++) {
      
      for (Int_t iSt = 0; iSt < 6; iSt++) de[iSt][1] = sigmaStCalc[iSt][iMode][i]*sigmaStCalc[iSt][iMode][i];
      
      sigmaSpectroCalc[iMode][i] = 0.;
      
      for (Int_t j = 1; j < 32; j++) {
	Double_t sigmaAdd = 1.;
	for (Int_t iSt = 0; iSt < 5; iSt++) sigmaAdd *= de[iSt][TESTBIT(j,iSt)];
	sigmaSpectroCalc[iMode][i] += sigmaAdd;
      }
      
      sigmaSpectroCalc[iMode][i] = TMath::Sqrt(sigmaSpectroCalc[iMode][i]);
      
      sigmaSpectro45Calc[iMode][i] = 0.;
      
      for (Int_t j = 1; j < 16; j++) {
	Double_t sigmaAdd = de[5][TESTBIT(j,3)];
	for (Int_t iSt = 0; iSt < 3; iSt++) sigmaAdd *= de[iSt][TESTBIT(j,iSt)];
	sigmaSpectro45Calc[iMode][i] += sigmaAdd;
      }
      
      sigmaSpectro45Calc[iMode][i] = TMath::Sqrt(sigmaSpectro45Calc[iMode][i]);
      
    }
    
  }
  
  // print spectrometer results
  printf("\nspectrometer efficiency:\n");
  PrintError(effSpectro, sigmaSpectroMeas, sigmaSpectroCalc);
  
  printf("\nspectrometer efficiency (3/4ch in st45):\n");
  PrintError(effSpectro45, sigmaSpectro45Meas, sigmaSpectro45Calc);
  
  // plot efficiency per chamber
  TCanvas *cEffCh = new TCanvas("cEffCh", "efficiency per chamber", 1000, 400);
  cEffCh->Divide(5,2);
  for (Int_t iCh = 0; iCh < 10; iCh++) {
    cEffCh->cd(iCh+1);
    hEffCh[iCh]->Rebin(10);
    hEffCh[iCh]->Draw();
  }
  
  // plot efficiency per station
  TCanvas *cEffSt = new TCanvas("cEffSt", "efficiency per station", 600, 400);
  cEffSt->Divide(3,2);
  for (Int_t iSt = 0; iSt < 6; iSt++) {
    cEffSt->cd(iSt+1);
    hEffSt[iSt]->Rebin(10);
    hEffSt[iSt]->Draw();
  }
  
  // plot spectrometer efficiency
  TCanvas *cEffSpectro = new TCanvas("cEffSpectro", "spectrometer efficiency", 400, 200);
  cEffSpectro->Divide(2,1);
  cEffSpectro->cd(1);
  hEffSpectro->Rebin(10);
  hEffSpectro->Draw();
  cEffSpectro->cd(2);
  hEffSpectro45->Rebin(10);
  hEffSpectro45->Draw();
  
}

//----------------------------------------------------------------------------
void PrintError(Double_t eff, Double_t sigmaMeas[2][3], Double_t sigmaCalc[knMode][3])
{
  /// print results per chamber
  
  if (kRefMode < 0 || kRefMode > knMode+2) {
    printf("bad reference\n");
    return;
  }
  
  for (Int_t iMode = 0; iMode < 2; iMode++) {
    
    Double_t diff[2];
    if (kRefMode <= 2) {
      diff[0] = (sigmaMeas[kRefMode][0] > 0.) ? 100.*(sigmaMeas[iMode][0]-sigmaMeas[kRefMode][0])/sigmaMeas[kRefMode][0] : 0.;
      diff[1] = (sigmaMeas[kRefMode][1] > 0.) ? 100.*(sigmaMeas[iMode][1]-sigmaMeas[kRefMode][1])/sigmaMeas[kRefMode][1] : 0.;
    } else {
      diff[0] = (sigmaCalc[kRefMode-2][0] > 0.) ? 100.*(sigmaMeas[iMode][0]-sigmaCalc[kRefMode-2][0])/sigmaCalc[kRefMode-2][0] : 0.;
      diff[1] = (sigmaCalc[kRefMode-2][1] > 0.) ? 100.*(sigmaMeas[iMode][1]-sigmaCalc[kRefMode-2][1])/sigmaCalc[kRefMode-2][1] : 0.;
    }
    
    if (iMode == kRefMode) printf("%s %f \e[34m+ %f (  ref  ) - %f (  ref  )\e[0m", sigmaModeMeas[iMode].Data(), eff, sigmaMeas[iMode][1], sigmaMeas[iMode][0]);
    else if ((iMode == 0 && kRefMode == 2) || (iMode == 1 && kRefMode == 3)) printf("%s %f \e[34m+ %f (%6.2f%%) - %f (%6.2f%%)\e[0m", sigmaModeMeas[iMode].Data(), eff, sigmaMeas[iMode][1], diff[1], sigmaMeas[iMode][0], diff[0]);
    else printf("%s %f + %f (%6.2f%%) - %f (%6.2f%%)", sigmaModeMeas[iMode].Data(), eff, sigmaMeas[iMode][1], diff[1], sigmaMeas[iMode][0], diff[0]);
    if (iMode == 0) printf(" | \e[31m± %f (  ref  )\e[0m\n", sigmaMeas[iMode][2]);
    else printf("\n");
    
  }
  
  printf("calculated errors:\n");
  for (Int_t iMode = 0; iMode < knMode; iMode++) {
    
    Double_t diff[3];
    if (kRefMode <= 2) {
      diff[0] = (sigmaMeas[kRefMode][0] > 0.) ? 100.*(sigmaCalc[iMode][0]-sigmaMeas[kRefMode][0])/sigmaMeas[kRefMode][0] : 0.;
      diff[1] = (sigmaMeas[kRefMode][1] > 0.) ? 100.*(sigmaCalc[iMode][1]-sigmaMeas[kRefMode][1])/sigmaMeas[kRefMode][1] : 0.;
    } else {
      diff[0] = (sigmaCalc[kRefMode-2][0] > 0.) ? 100.*(sigmaCalc[iMode][0]-sigmaCalc[kRefMode-2][0])/sigmaCalc[kRefMode-2][0] : 0.;
      diff[1] = (sigmaCalc[kRefMode-2][1] > 0.) ? 100.*(sigmaCalc[iMode][1]-sigmaCalc[kRefMode-2][1])/sigmaCalc[kRefMode-2][1] : 0.;
    }
    if (iMode == 0) diff[2] = (sigmaMeas[0][2] > 0.) ? 100.*(sigmaCalc[iMode][2]-sigmaMeas[0][2])/sigmaMeas[0][2] : 0.;
    
    if (iMode == kRefMode-2) printf("%s %f \e[34m+ %f (  ref  ) - %f (  ref  )\e[0m", sigmaModeCalc[iMode].Data(), eff, sigmaCalc[iMode][1], sigmaCalc[iMode][0]);
    else if ((iMode == 0 && kRefMode == 0) || (iMode == 1 && kRefMode == 1)) printf("%s %f \e[34m+ %f (%6.2f%%) - %f (%6.2f%%)\e[0m", sigmaModeCalc[iMode].Data(), eff, sigmaCalc[iMode][1], diff[1], sigmaCalc[iMode][0], diff[0]);
    else printf("%s %f + %f (%6.2f%%) - %f (%6.2f%%)", sigmaModeCalc[iMode].Data(), eff, sigmaCalc[iMode][1], diff[1], sigmaCalc[iMode][0], diff[0]);
    if (iMode == 0) printf(" | \e[31m± %f (%6.2f%%)\e[0m\n", sigmaCalc[iMode][2], diff[2]);
    else printf("\n");
    
  }
  
}

