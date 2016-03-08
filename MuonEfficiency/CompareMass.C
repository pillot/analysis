//
//  CompareMass.C
//  aliphysics-dev
//
//  Created by philippe pillot on 07/03/2016.
//  Copyright © 2016 Philippe Pillot. All rights reserved.
//

#if !defined(__CINT__) || defined(__MAKECINT__)

#include <Riostream.h>
#include "TFile.h"
#include "TObjArray.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TH1F.h"
#include "TF1.h"
#include "TMath.h"
#include "TStyle.h"
#include "TPaveStats.h"
#include "TPaveText.h"

#endif

Double_t CrystalBallExtended(Double_t *x,Double_t *par);
void SetFctParam(TF1 *fitFct);

//------------------------------------------------------------------------------
void CompareMass(TString fileName1, TString fileName2)
{
  /// compare Acc*Eff results between file1 and file2
  
  gStyle->SetOptStat(1110);
  gStyle->SetFillColor(0);
  
  TString fileName[2] = {fileName1.Data(), fileName2.Data()};
  Int_t nFiles = 2;
  
  // get results
  TFile* outFile[2];
  TObjArray *MassVspT[2];
  TObjArray *MassVsy[2];
  for (Int_t i = 0; i < nFiles; i++) {
    outFile[i] = TFile::Open(fileName[i].Data(),"READ");
    if (!outFile[i] || !outFile[i]->IsOpen()) return;
    MassVspT[i] = static_cast<TObjArray*>(outFile[i]->FindObjectAny("MassVspT"));
    MassVsy[i] = static_cast<TObjArray*>(outFile[i]->FindObjectAny("MassVsy"));
    if (!MassVspT[i] || !MassVsy[i]) return;
  }
  
  // fit function
  TF1 *fitFct = new TF1("fitFct", CrystalBallExtended, 0., 6., 7);
  
  // draw histos vs pT
  Int_t color[2] = {1, 2};
  Int_t npT = MassVspT[0]->GetEntries();
  if (MassVspT[1]->GetEntries() != npT) return;
  Int_t nPadsx = TMath::CeilNint(TMath::Sqrt(npT));
  Int_t nPadsy = TMath::FloorNint(TMath::Sqrt(npT));
  TCanvas* cpT = new TCanvas("massVtpT", "mass vs pT", TMath::Max(300*nPadsx, 1200), TMath::Max(300*nPadsy, 1200));
  cpT->Divide(nPadsx,nPadsy);
  for (Int_t ipT = 0; ipT < npT; ++ipT) {
    cpT->cd(ipT+1);
    gPad->SetLogy();
    
    ((TH1F*)MassVspT[0]->UncheckedAt(ipT))->SetLineColor(color[0]);
    fitFct->SetLineColor(color[0]);
    SetFctParam(fitFct);
    ((TH1F*)MassVspT[0]->UncheckedAt(ipT))->Fit(fitFct,"RIML+","e0");
    gPad->Update();
    TPaveStats *stats = (TPaveStats*)gPad->GetPrimitive("stats");
    stats->SetName("stats1");
    stats->SetX1NDC(0.6);
    stats->SetX2NDC(0.8);
    stats->SetY1NDC(0.8);
    stats->SetY2NDC(0.9);
    stats->SetTextColor(color[0]);
    TPaveText *txt = new TPaveText(0.8, 0.8, 0.9, 0.87, "NBNDC");
    txt->AddText(Form("m %5.3f",fitFct->GetParameter(1)))->SetTextColor(color[0]);
    txt->AddText(Form("s %5.4f",fitFct->GetParameter(2)))->SetTextColor(color[0]);
    txt->Draw("same");
    
    ((TH1F*)MassVspT[1]->UncheckedAt(ipT))->SetLineColor(color[1]);
    fitFct->SetLineColor(color[1]);
    SetFctParam(fitFct);
    ((TH1F*)MassVspT[1]->UncheckedAt(ipT))->Fit(fitFct,"RIML+","e0sames");
    gPad->Update();
    stats = (TPaveStats*)gPad->GetPrimitive("stats");
    stats->SetName("stats2");
    stats->SetX1NDC(0.6);
    stats->SetX2NDC(0.8);
    stats->SetY1NDC(0.7);
    stats->SetY2NDC(0.8);
    stats->SetTextColor(color[1]);
    txt = new TPaveText(0.8, 0.7, 0.9, 0.77, "NBNDC");
    txt->AddText(Form("m %5.3f",fitFct->GetParameter(1)))->SetTextColor(color[1]);
    txt->AddText(Form("s %5.4f",fitFct->GetParameter(2)))->SetTextColor(color[1]);
    txt->Draw("same");
    
    gPad->Modified();
  }
  
  // draw histos vs y
  Int_t ny = MassVsy[0]->GetEntries();
  if (MassVsy[1]->GetEntries() != ny) return;
  nPadsx = TMath::CeilNint(TMath::Sqrt(ny));
  nPadsy = TMath::FloorNint(TMath::Sqrt(ny));
  TCanvas* cy = new TCanvas("massVty", "mass vs y", TMath::Max(300*nPadsx, 1200), TMath::Max(300*nPadsy, 1200));
  cy->Divide(nPadsx,nPadsy);
  for (Int_t iy = 0; iy < ny; ++iy) {
    cy->cd(iy+1);
    gPad->SetLogy();
    
    ((TH1F*)MassVsy[0]->UncheckedAt(iy))->SetLineColor(color[0]);
    fitFct->SetLineColor(color[0]);
    SetFctParam(fitFct);
    ((TH1F*)MassVsy[0]->UncheckedAt(iy))->Fit(fitFct,"RIML+","e0");
    gPad->Update();
    TPaveStats *stats = (TPaveStats*)gPad->GetPrimitive("stats");
    stats->SetName("stats1");
    stats->SetX1NDC(0.6);
    stats->SetX2NDC(0.8);
    stats->SetY1NDC(0.8);
    stats->SetY2NDC(0.9);
    stats->SetTextColor(color[0]);
    TPaveText *txt = new TPaveText(0.8, 0.8, 0.9, 0.9, "NBNDC");
    txt->AddText(Form("m %5.3f",fitFct->GetParameter(1)))->SetTextColor(color[0]);
    txt->AddText(Form("s %5.4f",fitFct->GetParameter(2)))->SetTextColor(color[0]);
    txt->Draw("same");
    
    ((TH1F*)MassVsy[1]->UncheckedAt(iy))->SetLineColor(color[1]);
    fitFct->SetLineColor(color[1]);
    SetFctParam(fitFct);
    ((TH1F*)MassVsy[1]->UncheckedAt(iy))->Fit(fitFct,"RIML+","e0sames");
    gPad->Update();
    stats = (TPaveStats*)gPad->GetPrimitive("stats");
    stats->SetName("stats2");
    stats->SetX1NDC(0.6);
    stats->SetX2NDC(0.8);
    stats->SetY1NDC(0.7);
    stats->SetY2NDC(0.8);
    stats->SetTextColor(color[1]);
    txt = new TPaveText(0.8, 0.7, 0.9, 0.8, "NBNDC");
    txt->AddText(Form("m %5.3f",fitFct->GetParameter(1)))->SetTextColor(color[1]);
    txt->AddText(Form("s %5.4f",fitFct->GetParameter(2)))->SetTextColor(color[1]);
    txt->Draw("same");
    
    gPad->Modified();
  }
  
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

//------------------------------------------------------------------------------
void SetFctParam(TF1 *fitFct)
{
  fitFct->SetParameter(0, 5000.);
  fitFct->SetParameter(1, 3.);
  fitFct->SetParameter(2, 0.06);
  fitFct->SetParameter(3, 1.);
  fitFct->SetParLimits(3, 0., 10.);
  fitFct->SetParameter(4, 5.);
  fitFct->SetParLimits(4, 0., 10.);
  fitFct->SetParameter(5, 2.);
  fitFct->SetParLimits(5, 0., 10.);
  fitFct->SetParameter(6, 3.);
  fitFct->SetParLimits(6, 0., 10.);
}

