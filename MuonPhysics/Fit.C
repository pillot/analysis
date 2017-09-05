/*
 *  Fit.C
 *  aliroot_dev
 *
 *  Created by philippe pillot on 15/10/12.
 *  Copyright 2012 Subatech. All rights reserved.
 *
 */

#if !defined(__CINT__) || defined(__MAKECINT__)

#include <Riostream.h>
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TH1F.h"
#include "TF1.h"
#include "TMath.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TInterpreter.h"
#include "TSystem.h"

#endif

//------------------------------------------------------------------------------
Double_t BackgroundVWG(Double_t *x, Double_t *par)
{
  // gaussian variable width
  Double_t sigma = par[2]+par[3]*((x[0]-par[1])/par[1]);
  return par[0]*TMath::Exp(-(x[0]-par[1])*(x[0]-par[1])/(2.*sigma*sigma));
  
}

//------------------------------------------------------------------------------
Double_t CrystalBallExtended(Double_t *x,Double_t *par)
{
  //par[0] = Normalization
  //par[1] = mean
  //par[2] = sigma
  //par[3] = alpha
  //par[4] = n
  //par[5] = alpha'
  //par[6] = n'  
  
  Double_t t = (x[0]-par[1])/par[2];
  if (par[2] < 0) t = -t;
  
  Double_t absAlpha = fabs((Double_t)par[3]);
  Double_t absAlpha2 = fabs((Double_t)par[5]);
  
  if (t >= -absAlpha && t < absAlpha2) // gaussian core
  {
    return par[0]*(exp(-0.5*t*t));
  }
  
  if (t < -absAlpha) //left tail
  {
    Double_t a =  TMath::Power(par[4]/absAlpha,par[4])*exp(-0.5*absAlpha*absAlpha);
    Double_t b = par[4]/absAlpha - absAlpha;
    return par[0]*(a/TMath::Power(b - t, par[4]));
  }
  
  if (t >= absAlpha2) //right tail
  {
    
    Double_t c =  TMath::Power(par[6]/absAlpha2,par[6])*exp(-0.5*absAlpha2*absAlpha2);
    Double_t d = par[6]/absAlpha2 - absAlpha2;
    return par[0]*(c/TMath::Power(d + t, par[6]));
  }
  
  return 0. ; 
} 

//------------------------------------------------------------------------------
Double_t Gaus(Double_t *x, Double_t *par)
{
  // gaussian
  return par[0]/TMath::Sqrt(2.*TMath::Pi())/par[2]*TMath::Exp(-(x[0]-par[1])*(x[0]-par[1])/(2.*par[2]*par[2]));
  
}

//------------------------------------------------------------------------------
Double_t Exp(Double_t *x, Double_t *par)
{
  // exponential
  return par[0]*TMath::Exp(par[1]*x[0]);
  
}

//------------------------------------------------------------------------------
Double_t Pow(Double_t *x, Double_t *par)
{
  // power law
  return par[0]*TMath::Power(x[0],par[1]);
  
}

//------------------------------------------------------------------------------
Double_t fitFunctionVWG(Double_t *x, Double_t *par)
{
  if (x[0] > 2.9 && x[0] < 3.3) TF1::RejectPoint();
  return BackgroundVWG(x, par);
}

//------------------------------------------------------------------------------
Double_t fitFunctionCB2VWG(Double_t *x, Double_t *par)
{
  return BackgroundVWG(x, par) + CrystalBallExtended(x, &par[4]);
}

//------------------------------------------------------------------------------
Double_t fitFunction2CB2VWG(Double_t *x, Double_t *par)
{
  Double_t par2[7] = {
    par[11],
    3.68609+(par[5]-3.096916)/3.096916*3.68609,
    par[6]/3.096916*3.68609,
    par[7],
    par[8],
    par[9],
    par[10]
  };
  return BackgroundVWG(x, par) + CrystalBallExtended(x, &par[4]) + CrystalBallExtended(x, par2);
}

//------------------------------------------------------------------------------
Double_t fitFunction3GausExp(Double_t *x, Double_t *par)
{
  Double_t par2[3] = {par[4]*par[7], 10.02326+(par[5]-9.4603)/9.4603*10.02326, par[6]/9.4603*10.02326};
  Double_t par3[3] = {par[4]*par[7]*par[8], 10.3552+(par[5]-9.4603)/9.4603*10.3552, par[6]/9.4603*10.3552};
  return Exp(x, par) + Exp(x, &par[2]) + Gaus(x, &par[4]) + Gaus(x, par2) + Gaus(x, par3);
}

//------------------------------------------------------------------------------
Double_t fitFunction3GausPow(Double_t *x, Double_t *par)
{
  Double_t par2[3] = {par[4]*par[7], 10.02326+(par[5]-9.4603)/9.4603*10.02326, par[6]/9.4603*10.02326};
  Double_t par3[3] = {par[4]*par[7]*par[8], 10.3552+(par[5]-9.4603)/9.4603*10.3552, par[6]/9.4603*10.3552};
  return Pow(x, par) + Pow(x, &par[2]) + Gaus(x, &par[4]) + Gaus(x, par2) + Gaus(x, par3);
}

//------------------------------------------------------------------------------
void Fit(TString fileName = "AnalysisResults.root", TString containerName = "Histograms", TString collision = "pp")
{
  /// Fit the JPsi and Upsilon mass distribution
  
  collision.ToUpper();
  
  gROOT->SetStyle("Plain");
  gROOT->ForceStyle();
  gStyle->SetOptFit(1);
  
  // prepare environment
  if (!gInterpreter->IsLoaded("$WORK/Macros/Facilities/runTaskFacilities.C")) {
    gROOT->LoadMacro("$WORK/Macros/Facilities/runTaskFacilities.C");
    gROOT->ProcessLineFast("LoadAlirootLocally(\"\",\"include\",\"\",\"\");");
  }
  
  // open file
  TFile* file = TFile::Open(fileName.Data(),"READ");
  if (!file || !file->IsOpen()) return;
  
  // get histo
  TObjArray* histos = static_cast<TObjArray*>(file->FindObjectAny(containerName.Data()));
  if (!histos) return;
  TH1F* h = static_cast<TH1F*>(histos->FindObject("hMass"));
  if (!h) return;
  TH1F* hJPsi = (TH1F*)h->Clone("JPsi");
  hJPsi->SetDirectory(0);
  TH1F* hUpsilon = (TH1F*)h->Clone("hUpsilon");
  hUpsilon->SetDirectory(0);
  Bool_t isW = kFALSE;
  TH1F *hPtMuPlus = 0x0, *hPtMuMinus = 0x0, *hPt = 0x0, *hPtPlusOverMinus = 0x0;
  h = static_cast<TH1F*>(histos->FindObject("hPtMuPlus"));
  if (h) {
    isW = kTRUE;
    hPtMuPlus = (TH1F*)h->Clone("hptMuPlus");
    hPtMuPlus->Rebin(25);
    hPtMuPlus->Sumw2();
    hPtMuPlus->SetDirectory(0);
    h = static_cast<TH1F*>(histos->FindObject("hPtMuMinus"));
    if (!h) return;
    hPtMuMinus = (TH1F*)h->Clone("hptMuMinus");
    hPtMuMinus->Rebin(25);
    hPtMuMinus->Sumw2();
    hPtMuMinus->SetDirectory(0);
    hPt = (TH1F*)hPtMuPlus->Clone("hpt");
    hPt->Add(hPtMuMinus);
    hPt->SetTitle("#mu^{#pm} transverse momentum distribution");
    hPt->SetDirectory(0);
    hPtPlusOverMinus = (TH1F*)hPtMuPlus->Clone("hptMuPlusOverMuMinus");
    hPtPlusOverMinus->Divide(hPtMuMinus);
    hPtPlusOverMinus->SetTitle("#mu^{+} / #mu^{-} transverse momentum distribution ratio");
    hPtPlusOverMinus->SetDirectory(0);
  }
  
  // close file
  file->Close();
  delete file;
  
  // plot histo
  TCanvas* c = new TCanvas();
  c->SetWindowSize(900,400);
  c->Divide(2,1);
  
  c->cd(1);
  gPad->SetLogy();
  hJPsi->Rebin(5);
  hJPsi->GetXaxis()->SetRangeUser(2.,5.);
  hJPsi->Draw("e0");
  
  // BKG fit function
  /*
  TF1 *fitVWG = new TF1("fitVWG",fitFunctionVWG,2.4,4.5,4);
  fitVWG->SetParameter(0, 100.);
  fitVWG->SetParameter(1, 1.80243);
  fitVWG->SetParameter(2, 0.586137);
  fitVWG->SetParameter(3, 0.268486);
  hJPsi->Fit(fitVWG,"R");
  fitVWG->Draw("same");
  */
  /*
  // JPsi fit function
  TF1 *fitFctJPsi = new TF1("fitFctJPsi", CrystalBallExtended, 0., 4.9, 7);
  fitFctJPsi->SetLineColor(2);
  fitFctJPsi->SetParNames("kPsi","mPsi","sPsi","alPsi","nlPsi","auPsi","nuPsi");
  fitFctJPsi->SetParameter(0, 5000.); 
  fitFctJPsi->SetParameter(1, 3.); 
  fitFctJPsi->SetParameter(2, 0.08);
  fitFctJPsi->SetParameter(3, 1.);   
  fitFctJPsi->SetParLimits(3, 0., 5.);
  fitFctJPsi->SetParameter(4, 5.);   
  fitFctJPsi->SetParLimits(4, 0., 10.);
  fitFctJPsi->SetParameter(5, 2.);   
  fitFctJPsi->SetParLimits(5, 0., 5.);
  fitFctJPsi->SetParameter(6, 3.);
  fitFctJPsi->SetParLimits(6, 0., 10.);
  hJPsi->Fit(fitFctJPsi,"RIM","e0");
  fitFctJPsi->Draw("same");
  */
  
  TF1 *fitFctCB2VWG = 0x0;
  if (collision == "PBPB") {
    fitFctCB2VWG = new TF1("fitFctCB2VWG", fitFunctionCB2VWG, 2.2, 5., 11);
    fitFctCB2VWG->SetLineColor(2);
    fitFctCB2VWG->SetParNames("kVWG","mVWG","sVWG1","sVWG2","kPsi","mPsi","sPsi","alPsi","nlPsi","auPsi","nuPsi");
    fitFctCB2VWG->SetParameter(0, 10000.);
    fitFctCB2VWG->SetParameter(1, 1.9);
    fitFctCB2VWG->SetParameter(2, 0.5);
    fitFctCB2VWG->SetParLimits(2, 0., 100.);
    fitFctCB2VWG->SetParameter(3, 0.3);
    fitFctCB2VWG->SetParLimits(3, 0., 100.);
    fitFctCB2VWG->SetParameter(4, 100.); 
    fitFctCB2VWG->SetParameter(5, 3.13); 
    fitFctCB2VWG->SetParLimits(5, 3.08, 3.2);  
    fitFctCB2VWG->SetParameter(6, 0.08);
    fitFctCB2VWG->SetParLimits(6, 0.05, 0.15);
    fitFctCB2VWG->FixParameter(7, 1.12);
    fitFctCB2VWG->FixParameter(8, 3.96);
    fitFctCB2VWG->FixParameter(9, 2.5);
    fitFctCB2VWG->FixParameter(10, 2.48);
    hJPsi->Fit(fitFctCB2VWG,"IER","e0");
    fitFctCB2VWG->Draw("same");
  } else {
    fitFctCB2VWG = new TF1("fitFctCB2VWG", fitFunction2CB2VWG, 2., 5., 12);
    fitFctCB2VWG->SetLineColor(2);
    fitFctCB2VWG->SetParNames("kVWG","mVWG","sVWG1","sVWG2","kPsi","mPsi","sPsi","alPsi","nlPsi","auPsi","nuPsi");
    fitFctCB2VWG->SetParName(11, "kPsi'");
    fitFctCB2VWG->SetParameter(0, 10000.);
    fitFctCB2VWG->SetParameter(1, 1.9);
    //fitFctCB2VWG->SetParameter(1, 0.9);
    fitFctCB2VWG->SetParameter(2, 0.5);
    //fitFctCB2VWG->SetParameter(2, 0.3);
    fitFctCB2VWG->SetParLimits(2, 0., 100.);
    fitFctCB2VWG->SetParameter(3, 0.3);
    fitFctCB2VWG->SetParLimits(3, 0., 100.);
    fitFctCB2VWG->SetParameter(4, 100.); 
    fitFctCB2VWG->SetParameter(5, 3.13); 
    fitFctCB2VWG->SetParLimits(5, 3.08, 3.2);  
    fitFctCB2VWG->SetParameter(6, 0.08);
    fitFctCB2VWG->SetParLimits(6, 0.05, 0.15);
//    fitFctCB2VWG->FixParameter(7, 0.93);
//    fitFctCB2VWG->FixParameter(8, 5.59);
//    fitFctCB2VWG->FixParameter(9, 2.32);
//    fitFctCB2VWG->FixParameter(10, 3.39);
    fitFctCB2VWG->FixParameter(7, 0.984);
    fitFctCB2VWG->FixParameter(8, 5.839);
    fitFctCB2VWG->FixParameter(9, 1.972);
    fitFctCB2VWG->FixParameter(10, 3.444);
    fitFctCB2VWG->SetParameter(11, 10.);
    hJPsi->Fit(fitFctCB2VWG,"IER","e0");
    fitFctCB2VWG->Draw("same");
  }
  
  TF1 *fitFctCB2 = new TF1("fitFctCB2", CrystalBallExtended, 2.2, 5., 7);
  fitFctCB2->SetLineColor(4);
  fitFctCB2->SetParNames("kPsi","mPsi","sPsi","alPsi","nlPsi","auPsi","nuPsi");
  fitFctCB2->FixParameter(0, fitFctCB2VWG->GetParameter(4)); 
  fitFctCB2->FixParameter(1, fitFctCB2VWG->GetParameter(5)); 
  fitFctCB2->FixParameter(2, fitFctCB2VWG->GetParameter(6));
  fitFctCB2->FixParameter(3, fitFctCB2VWG->GetParameter(7));   
  fitFctCB2->FixParameter(4, fitFctCB2VWG->GetParameter(8));   
  fitFctCB2->FixParameter(5, fitFctCB2VWG->GetParameter(9));   
  fitFctCB2->FixParameter(6, fitFctCB2VWG->GetParameter(10));
  fitFctCB2->Draw("same");
  printf("\n --> nJPsi = %f\n\n", fitFctCB2->Integral(0., 100.)/hJPsi->GetBinWidth(1));
  
  if (collision != "PBPB") {
    TF1 *fitFctCB22 = new TF1("fitFctCB22", CrystalBallExtended, 2.2, 5., 7);
    fitFctCB22->SetLineColor(4);
    fitFctCB22->SetParNames("kPsi'","mPsi'","sPsi'","alPsi'","nlPsi'","auPsi'","nuPsi'");
    fitFctCB22->FixParameter(0, fitFctCB2VWG->GetParameter(11)); 
    fitFctCB22->FixParameter(1, 3.68609+(fitFctCB2VWG->GetParameter(5)-3.096916)/3.096916*3.68609); 
    fitFctCB22->FixParameter(2, fitFctCB2VWG->GetParameter(6)/3.096916*3.68609);
    fitFctCB22->FixParameter(3, fitFctCB2VWG->GetParameter(7));   
    fitFctCB22->FixParameter(4, fitFctCB2VWG->GetParameter(8));   
    fitFctCB22->FixParameter(5, fitFctCB2VWG->GetParameter(9));   
    fitFctCB22->FixParameter(6, fitFctCB2VWG->GetParameter(10));
    fitFctCB22->Draw("same");
    printf(" --> nJPsi' = %f\n\n", fitFctCB22->Integral(0., 100.)/hJPsi->GetBinWidth(1));
  }
  
  c->cd(2);
  gPad->SetLogy();
  hUpsilon->Rebin(10);
  hUpsilon->GetXaxis()->SetRangeUser(7.,13.);
  hUpsilon->Draw("e0");
  
  // Upsilon fit function
//  TF1 *fitFctGausExp = new TF1("fitFctGausExp", fitFunction3GausPow, 6.5, 13.5, 9);
  TF1 *fitFctGausExp = new TF1("fitFctGausExp", fitFunction3GausExp, 6., 16., 9);
  fitFctGausExp->SetLineColor(2);
  fitFctGausExp->SetParNames("kExp1","eExp1","kExp2","eExp2","kUpsilon","mUpsilon","sUpsilon","kUpsilon'","fUpsilon\"");
  fitFctGausExp->SetParameter(0, 100.);
  fitFctGausExp->SetParLimits(0, 0., 1.e9);
  fitFctGausExp->SetParameter(1, -1.1);
  fitFctGausExp->SetParLimits(1, -100., 0.);
  fitFctGausExp->SetParameter(2, 10.);
  fitFctGausExp->SetParLimits(2, 0., 1.e9);
  fitFctGausExp->SetParameter(3, -1.1);
  fitFctGausExp->SetParLimits(3, -10., 0.);
  fitFctGausExp->SetParameter(4, 10.); 
  fitFctGausExp->SetParameter(5, 9.46); 
  fitFctGausExp->SetParLimits(5, 9.3,9.6);  
  fitFctGausExp->SetParameter(6, 0.15);
  fitFctGausExp->SetParLimits(6, 0.05, 0.30);
  fitFctGausExp->SetParameter(7, 0.5); 
  fitFctGausExp->SetParLimits(7, 0.001, 1.);
  if (collision == "PBPB") {
    fitFctGausExp->FixParameter(8, 0.5); 
  } else {
    fitFctGausExp->SetParameter(8, 0.5);
    fitFctGausExp->SetParLimits(8, 0.001, 1.);
  }
  hUpsilon->Fit(fitFctGausExp,"IER","e0");
  fitFctGausExp->Draw("same");
  
  TF1 *fitFctGaus = new TF1("fitFctGaus", Gaus, 5., 15., 3);
  fitFctGaus->SetLineColor(4);
  fitFctGaus->SetParNames("kUpsilon","mUpsilon","sUpsilon");
  fitFctGaus->FixParameter(0, fitFctGausExp->GetParameter(4)); 
  fitFctGaus->FixParameter(1, fitFctGausExp->GetParameter(5)); 
  fitFctGaus->FixParameter(2, fitFctGausExp->GetParameter(6));
  fitFctGaus->Draw("same");
  printf("\n --> nUpsilon = %f\n\n", fitFctGaus->Integral(0., 20.)/hUpsilon->GetBinWidth(1));
  
  TF1 *fitFctGaus2 = new TF1("fitFctGaus2", Gaus, 5., 15., 3);
  fitFctGaus2->SetLineColor(4);
  fitFctGaus2->SetParNames("kUpsilon'","mUpsilon'","sUpsilon'");
  fitFctGaus2->FixParameter(0, fitFctGausExp->GetParameter(4)*fitFctGausExp->GetParameter(7)); 
  fitFctGaus2->FixParameter(1, 10.02326+(fitFctGausExp->GetParameter(5)-9.4603)/9.4603*10.02326); 
  fitFctGaus2->FixParameter(2, fitFctGausExp->GetParameter(6)/9.4603*10.02326);
  fitFctGaus2->Draw("same");
  printf(" --> nUpsilon' = %f\n\n", fitFctGaus2->Integral(0., 20.)/hUpsilon->GetBinWidth(1));
  
  TF1 *fitFctGaus3 = new TF1("fitFctGaus3", Gaus, 5., 15., 3);
  fitFctGaus3->SetLineColor(4);
  fitFctGaus3->SetParNames("kUpsilon\"","mUpsilon\"","sUpsilon\"");
  fitFctGaus3->FixParameter(0, fitFctGausExp->GetParameter(4)*fitFctGausExp->GetParameter(7)*fitFctGausExp->GetParameter(8)); 
  fitFctGaus3->FixParameter(1, 10.3552+(fitFctGausExp->GetParameter(5)-9.4603)/9.4603*10.3552); 
  fitFctGaus3->FixParameter(2, fitFctGausExp->GetParameter(6)/9.4603*10.3552);
  fitFctGaus3->Draw("same");
  printf(" --> nUpsilon\" = %f\n\n", fitFctGaus3->Integral(0., 20.)/hUpsilon->GetBinWidth(1));
  
  // W
  if (isW) {
    TCanvas* c2 = new TCanvas();
    c2->SetWindowSize(800,800);
    c2->Divide(2,2);
    c2->cd(1);
    gPad->SetLogy();
    hPtMuPlus->Draw();
    c2->cd(2);
    gPad->SetLogy();
    hPtMuMinus->Draw();
    c2->cd(3);
    gPad->SetLogy();
    hPt->Draw();
    c2->cd(4);
    hPtPlusOverMinus->GetYaxis()->SetRangeUser(0., 1.2);
    hPtPlusOverMinus->Draw();
  }
  
}

