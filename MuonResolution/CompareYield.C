//
//  CompareYield.C
//  aliroot_dev
//
//  Created by philippe pillot on 19/08/2014.
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


//-----------------------------------------------------------------------
void CompareYield(TString file0 = "ResToYield.root", TString file1 = "ResToYield.root", TString file2 = "")
{
  /// compare reconstructed yields
  
  const Int_t nFiles = file2.IsNull() ? 2 : 3;
  TString fileName[nFiles];
  fileName[0] = file0;
  fileName[1] = file1;
  if (nFiles > 2) fileName[2] = file2;
  
  TH1 *hpGen[nFiles], *hpTGen[nFiles], *hetaGen[nFiles], *hyGen[nFiles], *hphiGen[nFiles];
  TH1 *hpRec[nFiles], *hpTRec[nFiles], *hetaRec[nFiles], *hyRec[nFiles], *hphiRec[nFiles];
  for (Int_t iFile = 0; iFile < nFiles; ++iFile) {
    TFile *fResults = TFile::Open(fileName[iFile].Data(),"READ");
    if (!fResults) return;
    hpGen[iFile] = static_cast<TH1*>(fResults->FindObjectAny("hpGen"));
    hpTGen[iFile] = static_cast<TH1*>(fResults->FindObjectAny("hpTGen"));
    hetaGen[iFile] = static_cast<TH1*>(fResults->FindObjectAny("hetaGen"));
    hyGen[iFile] = static_cast<TH1*>(fResults->FindObjectAny("hyGen"));
    hphiGen[iFile] = static_cast<TH1*>(fResults->FindObjectAny("hphiGen"));
    hpRec[iFile] = static_cast<TH1*>(fResults->FindObjectAny("hpRec"));
    hpTRec[iFile] = static_cast<TH1*>(fResults->FindObjectAny("hpTRec"));
    hetaRec[iFile] = static_cast<TH1*>(fResults->FindObjectAny("hetaRec"));
    hyRec[iFile] = static_cast<TH1*>(fResults->FindObjectAny("hyRec"));
    hphiRec[iFile] = static_cast<TH1*>(fResults->FindObjectAny("hphiRec"));
    if (!hpGen[iFile] || !hpTGen[iFile] || !hetaGen[iFile] || !hyGen[iFile] || !hphiGen[iFile] ||
        !hpRec[iFile] || !hpTRec[iFile] || !hetaRec[iFile] || !hyRec[iFile] || !hphiRec[iFile]) return;
    hpGen[iFile]->SetDirectory(0);
    hpTGen[iFile]->SetDirectory(0);
    hetaGen[iFile]->SetDirectory(0);
    hyGen[iFile]->SetDirectory(0);
    hphiGen[iFile]->SetDirectory(0);
    hpRec[iFile]->SetDirectory(0);
    hpTRec[iFile]->SetDirectory(0);
    hetaRec[iFile]->SetDirectory(0);
    hyRec[iFile]->SetDirectory(0);
    hphiRec[iFile]->SetDirectory(0);
    fResults->Close();
  }
  
  // drawing option
  TString drawOpt[nFiles];
  Int_t color[nFiles];
  drawOpt[0] = "";
  color[0] = 2;
  for (Int_t iFile = 1; iFile < nFiles; ++iFile) {
    drawOpt[iFile] = "same";
    color[iFile] = color[iFile-1] + 2;
  }
  
  // draw generated histograms
  TCanvas *cGen = new TCanvas("cGen","cGen",10,10,600,600);
  cGen->Divide(2,2);
  gROOT->SetSelectedPad(cGen->cd(1));
  gPad->SetLogy();
  for (Int_t iFile = 0; iFile < nFiles; ++iFile) {
    hpTGen[iFile]->SetLineColor(color[iFile]);
    hpTGen[iFile]->Draw(drawOpt[iFile].Data());
  }
  gROOT->SetSelectedPad(cGen->cd(2));
  gPad->SetLogy();
  for (Int_t iFile = 0; iFile < nFiles; ++iFile) {
    hetaGen[iFile]->SetLineColor(color[iFile]);
    hetaGen[iFile]->Draw(drawOpt[iFile].Data());
    hyGen[iFile]->SetLineColor(7+color[iFile]);
    hyGen[iFile]->Draw("same");
  }
  gROOT->SetSelectedPad(cGen->cd(3));
  gPad->SetLogy();
  for (Int_t iFile = 0; iFile < nFiles; ++iFile) {
    hpGen[iFile]->SetLineColor(color[iFile]);
    hpGen[iFile]->Draw(drawOpt[iFile].Data());
  }
  gROOT->SetSelectedPad(cGen->cd(4));
  gPad->SetLogy();
  for (Int_t iFile = 0; iFile < nFiles; ++iFile) {
    hphiGen[iFile]->SetLineColor(color[iFile]);
    hphiGen[iFile]->Draw(drawOpt[iFile].Data());
  }
  
  // draw reconstructed histograms
  TCanvas *cRec = new TCanvas("cRec","cRec",10,10,600,600);
  cRec->Divide(2,2);
  gROOT->SetSelectedPad(cRec->cd(1));
  gPad->SetLogy();
  for (Int_t iFile = 0; iFile < nFiles; ++iFile) {
    hpTRec[iFile]->SetLineColor(color[iFile]);
    hpTRec[iFile]->Draw(drawOpt[iFile].Data());
  }
  gROOT->SetSelectedPad(cRec->cd(2));
  gPad->SetLogy();
  for (Int_t iFile = 0; iFile < nFiles; ++iFile) {
    hetaRec[iFile]->SetLineColor(color[iFile]);
    hetaRec[iFile]->Draw(drawOpt[iFile].Data());
    hyRec[iFile]->SetLineColor(7+color[iFile]);
    hyRec[iFile]->Draw("same");
  }
  gROOT->SetSelectedPad(cRec->cd(3));
  gPad->SetLogy();
  for (Int_t iFile = 0; iFile < nFiles; ++iFile) {
    hpRec[iFile]->SetLineColor(color[iFile]);
    hpRec[iFile]->Draw(drawOpt[iFile].Data());
  }
  gROOT->SetSelectedPad(cRec->cd(4));
  gPad->SetLogy();
  for (Int_t iFile = 0; iFile < nFiles; ++iFile) {
    hphiRec[iFile]->SetLineColor(color[iFile]);
    hphiRec[iFile]->Draw(drawOpt[iFile].Data());
  }
  
  // draw ratios of reconstructed histograms
  TCanvas *cRat = new TCanvas("cRat","cRat",10,10,600,600);
  cRat->Divide(2,2);
  gROOT->SetSelectedPad(cRat->cd(1));
  for (Int_t iFile = 1; iFile < nFiles; ++iFile) {
    TH1 *hpTRat = static_cast<TH1*>(hpTRec[iFile]->Clone());
    hpTRat->SetNameTitle(Form("hpTRat%d",iFile),Form("hpTRat%d",iFile));
    hpTRat->Divide(hpTRec[0]);
    hpTRat->SetLineColor(color[iFile]);
    hpTRat->Draw(drawOpt[iFile-1].Data());
  }
  gROOT->SetSelectedPad(cRat->cd(2));
  for (Int_t iFile = 1; iFile < nFiles; ++iFile) {
    TH1 *hetaRat = static_cast<TH1*>(hetaRec[iFile]->Clone());
    hetaRat->SetNameTitle(Form("hetaRat%d",iFile),Form("hetaRat%d",iFile));
    hetaRat->Divide(hetaRec[0]);
    hetaRat->SetLineColor(color[iFile]);
    hetaRat->Draw(drawOpt[iFile-1].Data());
    TH1 *hyRat = static_cast<TH1*>(hyRec[iFile]->Clone());
    hyRat->SetNameTitle(Form("hyRat%d",iFile),Form("hyRat%d",iFile));
    hyRat->Divide(hyRec[0]);
    hyRat->SetLineColor(7+color[iFile]);
    hyRat->Draw("same");
  }
  gROOT->SetSelectedPad(cRat->cd(3));
  for (Int_t iFile = 1; iFile < nFiles; ++iFile) {
    TH1 *hpRat = static_cast<TH1*>(hpRec[iFile]->Clone());
    hpRat->SetNameTitle(Form("hpRat%d",iFile),Form("hpRat%d",iFile));
    hpRat->Divide(hpRec[0]);
    hpRat->SetLineColor(color[iFile]);
    hpRat->Draw(drawOpt[iFile-1].Data());
  }
  gROOT->SetSelectedPad(cRat->cd(4));
  for (Int_t iFile = 1; iFile < nFiles; ++iFile) {
    TH1 *hphiRat = static_cast<TH1*>(hphiRec[iFile]->Clone());
    hphiRat->SetNameTitle(Form("hphiRat%d",iFile),Form("hphiRat%d",iFile));
    hphiRat->Divide(hphiRec[0]);
    hphiRat->SetLineColor(color[iFile]);
    hphiRat->Draw(drawOpt[iFile-1].Data());
  }
  
}

