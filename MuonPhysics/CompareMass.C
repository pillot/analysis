//
//  CompareMass.C
//
//  Created by philippe pillot on 22/06/2015.
//  Copyright (c) 2015 Philippe Pillot. All rights reserved.
//

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
#include "TPaveText.h"

#endif

//Double_t fractionOfAfter[5] = {0.,0.15,0.25,0.5,0.66};
//Double_t fractionOfAfter[5] = {0.,0.42,0.49,0.5,0.66};
//Double_t fractionOfAfter[5] = {0.,0.39,0.63,0.79,1.};
//Double_t fractionOfAfter[5] = {0.,0.536,0.73,0.834,1.};
Double_t fractionOfAfter[5] = {0.,0.434,0.673,0.80,1.};


Double_t CrystalBallExtended(Double_t *x,Double_t *par);

//------------------------------------------------------------------------------
void CompareMass(TString fileName_before = "AnalysisResults.root", TString fileName_after = "AnalysisResults.root")
{
  /// compare mass plots, fit them and mix them with given fractions
  
  gROOT->SetStyle("Plain");
  gROOT->ForceStyle();
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  
  Bool_t sameFile = (fileName_before==fileName_after);
  //  TString extBefore = sameFile ? "_before" : "_after";
  //  TString extBefore = "_after";
  TString extBefore = "";
  //  TString extAfter = "_after";
  TString extAfter = "";
  
  // open files
  TFile* file_before = TFile::Open(fileName_before.Data(),"READ");
  if (!file_before || !file_before->IsOpen()) return;
  TFile* file_after = file_before;
  if (!sameFile) {
    file_after = TFile::Open(fileName_after.Data(),"READ");
    if (!file_after || !file_after->IsOpen()) return;
  }
  
  // histo before
  TObjArray* histo_before = static_cast<TObjArray*>(file_before->FindObjectAny(Form("Histograms%s",extBefore.Data())));
  if (!histo_before) return;
  
  // histo after
  TObjArray* histo_after = static_cast<TObjArray*>(file_after->FindObjectAny(Form("Histograms%s",extAfter.Data())));
  if (!histo_after) return;
  
  // draw differences
  TString sMass = "hMass";
  TCanvas* cHist = new TCanvas("cHist", "cHist");
  gPad->SetLogy();
  TH1F* h_before = static_cast<TH1F*>(histo_before->FindObject("hMass"));
  if (!h_before) return;
  TH1F* h_after = static_cast<TH1F*>(histo_after->FindObject("hMass"));
  if (!h_after) return;
  h_before->Draw("e0");
  h_after->Draw("e0sames");
  h_after->SetLineColor(2);
  
  // fit
  Bool_t psi = (h_before->GetMean() < 5.);
  TF1 *fitFct = new TF1("fitFct", CrystalBallExtended, psi ? 0. : 5., psi ? 5. : 15., 7);
  if (psi) {
    fitFct->SetParNames("kPsi","mPsi","sPsi","alPsi","nlPsi","auPsi","nuPsi");
    fitFct->SetParameter(1, 3.);
    fitFct->SetParameter(2, 0.08);
  }
  else {
    fitFct->SetParNames("kUps","mUps","sUps","alUps","nlUps","auUps","nuUps");
    fitFct->SetParameter(1, 10.);
    fitFct->SetParameter(2, 0.120);
  }
  fitFct->SetParameter(0, 5000.);
  fitFct->SetParameter(3, 1.);
  fitFct->SetParLimits(3, 0., 10.);
  fitFct->SetParameter(4, 5.);
  fitFct->SetParLimits(4, 0., 10.);
  fitFct->SetParameter(5, 2.);
  fitFct->SetParLimits(5, 0., 10.);
  fitFct->SetParameter(6, 3.);
  fitFct->SetParLimits(6, 0., 10.);
  fitFct->SetLineColor(4);
  h_before->Fit(fitFct,"RIML","e0sames");
  fitFct->DrawClone("same");
  Double_t sigma_before = fitFct->GetParameter(2)*1000.;
  fitFct->SetNpx(1000);
  fitFct->SetLineColor(2);
  h_after->Fit(fitFct,"RIML","e0sames");
  fitFct->DrawClone("same");
  Double_t sigma_after = fitFct->GetParameter(2)*1000.;
  
  // results
  TPaveText *txt = new TPaveText(0.15, 0.75, 0.45, 0.85, "NBNDC");
  txt->AddText(Form("sigma = %3.0f MeV/c^{2}",sigma_before))->SetTextColor(4);
  txt->AddText(Form("sigma = %3.0f MeV/c^{2}",sigma_after))->SetTextColor(2);
  txt->Draw("same");
  
  // draw and fit mixture
  Int_t color[5] = {1, 4, 2, 6, 8};
  TCanvas* cMix = new TCanvas("cMix", "cMix");
  gPad->SetLogy();
  TPaveText *txtMix = new TPaveText(0.15, 0.75, 0.45, 0.85, "NBNDC");
  fitFct->SetParameter(0, 10000);
  fitFct->SetLineWidth(2);
  for (Int_t i = 0; i < 5; ++i) {
    TH1F* h_Mixb = static_cast<TH1F*>(h_before->Clone("hMassMix"));
    h_Mixb->Sumw2();
    h_Mixb->Scale(1.-fractionOfAfter[i], "width");
    TH1F* h_Mixa = static_cast<TH1F*>(h_after->Clone("hMassMixb"));
    h_Mixa->Sumw2();
    h_Mixa->Scale(fractionOfAfter[i], "width");
    h_Mixb->Add(h_Mixa);
    h_Mixb->SetLineColor(color[i]);
    h_Mixb->Draw(i == 0 ? "e0" : "e0sames");
    fitFct->SetLineColor(color[i]);
    //if (i > 0) for (Int_t j = 3; j < 7; ++j) fitFct->FixParameter(j, fitFct->GetParameter(j));
    h_Mixb->Fit(fitFct,"RIML","e0sames");
    fitFct->DrawClone("same");
    Double_t sigma_mix = fitFct->GetParameter(2)*1000.;
    txtMix->AddText(Form("f = %4.2f: sigma = %3.0f MeV/c^{2}", fractionOfAfter[i], sigma_mix))->SetTextColor(color[i]);
    Double_t nEv = h_Mixb->Integral("width");
    Double_t nEvFit = fitFct->Integral(0., 20.);
    printf("nEv = %d; fit = %d; diff = %f%%\n", (Int_t)nEv, (Int_t)nEvFit, 100.*(nEvFit-nEv)/nEv);
  }
  txtMix->Draw("same");
  
  
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
  if (par[3] < 0) t = -t;
  
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

