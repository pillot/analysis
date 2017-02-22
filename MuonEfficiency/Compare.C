/*
 *  Compare.C
 *  aliroot_dev
 *
 *  Created by philippe pillot on 13/09/11.
 *  Copyright 2011 SUBATECH. All rights reserved.
 *
 */

#if !defined(__CINT__) || defined(__MAKECINT__)

#include <Riostream.h>
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include <TFile.h>
#include <TH1F.h>
#include <TMath.h>
#include <TStyle.h>
#include <TString.h>

#endif

void FillRatio(TGraphAsymmErrors *g, TGraphAsymmErrors *ref, TGraphAsymmErrors *rat);

//--------------------------------------------------------------------
void Compare(TString fileName0, TString fileNames, Bool_t vsCent = kFALSE)
{
  /// compare Acc*Eff results between fileO and other(s)
  /// in case of file list, histos in files preceded by '*' will be drawn in blue
  
  gStyle->SetOptStat(11);
  gStyle->SetFillColor(0);
  
  Int_t color[1000];
  Double_t markerSize[1000];
  TString fileName[1000];
  color[0] = 1;
  markerSize[0] = 0.2;
  fileName[0] = fileName0;
  Int_t nFiles = 1;
  if (fileNames.EndsWith(".root")) {
    color[nFiles] = 2;
    markerSize[nFiles] = 0.2;
    fileName[nFiles++] = fileNames;
  }
  else {
    ifstream inFile(fileNames.Data());
    if (!inFile.is_open()) {
      printf("cannot open file %s\n", fileNames.Data());
      return;
    }
    while (!inFile.eof()) {
      TString fName;
      fName.ReadLine(inFile,kTRUE);
      if (fName.IsNull() || fName.BeginsWith("#") || !fName.EndsWith(".root")) continue;
      if (fName.BeginsWith("*")) {
        color[nFiles] = 4;
        markerSize[nFiles] = 0.4;
        fName.Remove(TString::kLeading,'*');
      } else {
        color[nFiles] = 2;
        markerSize[nFiles] = 0.2;
      }
      fileName[nFiles++] = fName;
    }
  }
  
  const Int_t nData = 7;
  TString dataGen[nData] = {"hPtGen", "hYGen", "hPtGenMu", "hYGenMu", "hGenPtSummary", "hGenYSummary", "hGenSummary"};
  TString dataRec[nData] = {"hPtRec", "hYRec", "hPtRecMu", "hYRecMu", "hRecPtSummary", "hRecYSummary", "hRecSummary"};
  TString dataAcc[nData] = {"hPtAcc", "hYAcc", "hPtAccMu", "hYAccMu", "hAccPtSummary", "hAccYSummary", "hAccSummary"};
  TString cName[nData] = {"cPt", "cY", "cPtMu", "cYMu", "cPtSummary", "cYSummary", "cSummary"};
  TString cTitle[nData] = {"versus pT (J/Psi)", "versus y (J/Psi)", "versus pT (single-mu)",
    "versus y (single-mu)", "summary versus pt", "summary versus y", "summary versus pt/y"};
  
  // get results
  TFile* outFile[1000];
  TH1F* hDataGen[1000][nData];
  TH1F* hDataRec[1000][nData];
  TH1F* hDataAcc[1000][nData];
  Int_t nDataCentVsPt = 0, nDataCentVsY = 0;
  TGraphAsymmErrors *gDataAccCent[1000];
  TGraphAsymmErrors *gDataAccCentVsPt[1000][100];
  TGraphAsymmErrors *gDataAccCentVsY[1000][100];
  for (Int_t i = 0; i < nFiles; i++) {
    outFile[i] = TFile::Open(fileName[i].Data(),"READ");
    if (!outFile[i] || !outFile[i]->IsOpen()) return;
    for (Int_t j = 0; j < nData; j++) {
      hDataGen[i][j] = static_cast<TH1F*>(outFile[i]->FindObjectAny(dataGen[j].Data()));
      if (j < 4) hDataGen[i][j]->Sumw2();
      hDataRec[i][j] = static_cast<TH1F*>(outFile[i]->FindObjectAny(dataRec[j].Data()));
      if (j < 4) hDataRec[i][j]->Sumw2();
      hDataAcc[i][j] = static_cast<TH1F*>(outFile[i]->FindObjectAny(dataAcc[j].Data()));
    }
    if (vsCent) {
      gDataAccCent[i] = static_cast<TGraphAsymmErrors*>(outFile[i]->FindObjectAny("accEffVsCent_pt0_y0"));
      Int_t j = -1;
      do {
        ++j;
        gDataAccCentVsPt[i][j] = static_cast<TGraphAsymmErrors*>(outFile[i]->FindObjectAny(Form("accEffVsCent_pt%d_y0", j+1)));
      } while (gDataAccCentVsPt[i][j]);
      if (i == 0) nDataCentVsPt = j;
      else if (j != nDataCentVsPt) {
        printf("inconsistent pT bins for acceff vs centrality\n");
        return;
      }
      j = -1;
      do {
        ++j;
        gDataAccCentVsY[i][j] = static_cast<TGraphAsymmErrors*>(outFile[i]->FindObjectAny(Form("accEffVsCent_pt0_y%d", j+1)));
      } while (gDataAccCentVsY[i][j]);
      if (i == 0) nDataCentVsY = j;
      else if (j != nDataCentVsY) {
        printf("inconsistent y bins for acceff vs centrality\n");
        return;
      }
    }
  }
  
  // compute ratio
  TH1F* hGenRat[1000][nData];
  TH1F* hRecRat[1000][nData];
  TH1F* hAccRat[1000][nData];
  TGraphAsymmErrors *gDataAccRatCent[1000];
  TGraphAsymmErrors *gDataAccRatCentVsPt[1000][100];
  TGraphAsymmErrors *gDataAccRatCentVsY[1000][100];
  for (Int_t i = 1; i < nFiles; i++) {
    for (Int_t j = 0; j < nData; j++) {
      if (!hDataGen[i][j]) continue;
      hGenRat[i-1][j] = (TH1F*) hDataGen[i][j]->Clone(Form("%sRat",dataGen[j].Data()));
      hGenRat[i-1][j]->SetTitle("ratio");
      hGenRat[i-1][j]->Divide(hDataGen[0][j]);
      hRecRat[i-1][j] = (TH1F*) hDataRec[i][j]->Clone(Form("%sRat",dataRec[j].Data()));
      hRecRat[i-1][j]->SetTitle("ratio");
      hRecRat[i-1][j]->Divide(hDataRec[0][j]);
      hAccRat[i-1][j] = (TH1F*) hDataAcc[i][j]->Clone(Form("%sRat",dataAcc[j].Data()));
      hAccRat[i-1][j]->SetTitle("ratio");
      hAccRat[i-1][j]->Divide(hDataAcc[0][j]);
    }
    if (vsCent) {
      if (gDataAccCent[i] && gDataAccCent[0]) {
        gDataAccRatCent[i-1] = new TGraphAsymmErrors(gDataAccCent[0]->GetN());
        gDataAccRatCent[i-1]->SetNameTitle("accEffRatVsCent_pt0_y0", gDataAccCent[0]->GetTitle());
        FillRatio(gDataAccCent[i], gDataAccCent[0], gDataAccRatCent[i-1]);
      } else gDataAccRatCent[i-1] = 0x0;
      for (Int_t j = 0; j < nDataCentVsPt; j++) {
        if (gDataAccCentVsPt[i][j] && gDataAccCentVsPt[0][j]) {
          gDataAccRatCentVsPt[i-1][j] = new TGraphAsymmErrors(gDataAccCentVsPt[0][j]->GetN());
          gDataAccRatCentVsPt[i-1][j]->SetNameTitle(Form("accEffRatVsCent_pt%d_y0",j+1), gDataAccCentVsPt[0][j]->GetTitle());
          FillRatio(gDataAccCentVsPt[i][j], gDataAccCentVsPt[0][j], gDataAccRatCentVsPt[i-1][j]);
        } else gDataAccRatCentVsPt[i-1][j] = 0x0;
      }
      for (Int_t j = 0; j < nDataCentVsY; j++) {
        if (gDataAccCentVsY[i][j] && gDataAccCentVsY[0][j]) {
          gDataAccRatCentVsY[i-1][j] = new TGraphAsymmErrors(gDataAccCentVsY[0][j]->GetN());
          gDataAccRatCentVsY[i-1][j]->SetNameTitle(Form("accEffRatVsCent_pt0_y%d",j+1), gDataAccCentVsY[0][j]->GetTitle());
          FillRatio(gDataAccCentVsY[i][j], gDataAccCentVsY[0][j], gDataAccRatCentVsY[i-1][j]);
        } else gDataAccRatCentVsY[i-1][j] = 0x0;
      }
    }
  }
  
  /*
  // results from Javier (embedding)
  Double_t pTBinLowEdge[3] = {0., 3., 10.};
  Double_t yBinLowEdge[3] = {-4., -3.25, -2.5};
  Double_t acc[3][3] = {{18.8, 17.8, 19.6}, {17.0, 16.2, 17.8}, {21.7, 20.3, 22.9}};
  Double_t erracc[3][3] = {{0.1, 0.2, 0.2}, {0.2, 0.2, 0.2}, {0.2, 0.3, 0.3}};
  TH1F* hAccSummaryJ = new TH1F("hAccSummaryJ","Acc * Eff versus pt/y bins", 9, -0.5, 8.5);
  for (Int_t ipt = 0; ipt <3; ipt++) {
    for (Int_t iy = 0; iy <3; iy++) {
      TString label = "";
      label += (ipt == 0) ? Form("%g<pt<%g",pTBinLowEdge[0],pTBinLowEdge[2]) : Form("%g<pt<%g",pTBinLowEdge[ipt-1],pTBinLowEdge[ipt]);
      label += (iy == 0) ? Form(" / %g<y<%g",yBinLowEdge[0],yBinLowEdge[2]) : Form(" / %g<y<%g",yBinLowEdge[iy-1],yBinLowEdge[iy]);
      hAccSummaryJ->SetBinContent(ipt*3+iy+1, acc[ipt][iy]/100.);
      hAccSummaryJ->SetBinError(ipt*3+iy+1, erracc[ipt][iy]/100.);
      hAccSummaryJ->GetXaxis()->SetBinLabel(ipt*3+iy+1, label.Data());
    }
  }
  TH1F* hAccSummaryJRat = hAccSummaryJ->Clone("hAccSummaryJRat");
  hAccSummaryJRat->Divide(hAccSummary1);
  */
  
  // draw histos
  Int_t markerStyle = 21;
  Double_t scale[1000];
  for (Int_t j = 0; j < nData; j++) {
    TCanvas* c = new TCanvas(cName[j].Data(), cTitle[j].Data(), 900, 600);
    c->Divide(3,2);
    c->cd(1);
    gPad->SetLogy();
    if (!hDataGen[0][j]) continue;
    hDataGen[0][j]->SetLineColor(color[0]);
    hDataGen[0][j]->SetMarkerStyle(markerStyle);
    hDataGen[0][j]->SetMarkerSize(markerSize[0]);
    hDataGen[0][j]->SetMarkerColor(color[0]);
    hDataGen[0][j]->Draw("e0");
    for (Int_t i = 1; i < nFiles; i++) {
      if (j < 4) scale[i-1] = ((Double_t)hDataGen[0][j]->GetEntries())/((Double_t)hDataGen[i][j]->GetEntries());
      else scale[i-1] = ((Double_t)hDataGen[0][j]->GetBinContent(1))/((Double_t)hDataGen[i][j]->GetBinContent(1));
      hDataGen[i][j]->Scale(scale[i-1]);
      hDataGen[i][j]->SetLineColor(color[i]);
      hDataGen[i][j]->SetMarkerStyle(markerStyle);
      hDataGen[i][j]->SetMarkerSize(markerSize[i]);
      hDataGen[i][j]->SetMarkerColor(color[i]);
      hDataGen[i][j]->Draw("e0same");
    }
    c->cd(2);
    gPad->SetLogy();
    hDataRec[0][j]->SetLineColor(color[0]);
    hDataRec[0][j]->SetMarkerStyle(markerStyle);
    hDataRec[0][j]->SetMarkerSize(markerSize[0]);
    hDataRec[0][j]->SetMarkerColor(color[0]);
    hDataRec[0][j]->Draw("e0");
    for (Int_t i = 1; i < nFiles; i++) {
      hDataRec[i][j]->Scale(scale[i-1]);
      hDataRec[i][j]->SetLineColor(color[i]);
      hDataRec[i][j]->SetMarkerStyle(markerStyle);
      hDataRec[i][j]->SetMarkerSize(markerSize[i]);
      hDataRec[i][j]->SetMarkerColor(color[i]);
      hDataRec[i][j]->Draw("e0same");
    }
    c->cd(3);
    hDataAcc[0][j]->SetLineColor(color[0]);
    hDataAcc[0][j]->SetMarkerStyle(markerStyle);
    hDataAcc[0][j]->SetMarkerSize(markerSize[0]);
    hDataAcc[0][j]->SetMarkerColor(color[0]);
    hDataAcc[0][j]->Draw("e0");
    for (Int_t i = 1; i < nFiles; i++) {
      hDataAcc[i][j]->SetLineColor(color[i]);
      hDataAcc[i][j]->SetMarkerStyle(markerStyle);
      hDataAcc[i][j]->SetMarkerSize(markerSize[i]);
      hDataAcc[i][j]->SetMarkerColor(color[i]);
      hDataAcc[i][j]->Draw("e0same");
    }
    c->cd(4);
    for (Int_t i = 0; i < nFiles-1; i++) {
      hGenRat[i][j]->SetLineColor(color[i+1]);
      hGenRat[i][j]->SetMarkerStyle(markerStyle);
      hGenRat[i][j]->SetMarkerSize(markerSize[i+1]);
      hGenRat[i][j]->SetMarkerColor(color[i+1]);
      hGenRat[i][j]->Scale(scale[i]);
      if (i == 0) hGenRat[i][j]->Draw((j<4)?"e0":"histpx0");
      else hGenRat[i][j]->Draw((j<4)?"e0same":"histpx0same");
    }
    c->cd(5);
    for (Int_t i = 0; i < nFiles-1; i++) {
      hRecRat[i][j]->SetLineColor(color[i+1]);
      hRecRat[i][j]->SetMarkerStyle(markerStyle);
      hRecRat[i][j]->SetMarkerSize(markerSize[i+1]);
      hRecRat[i][j]->SetMarkerColor(color[i+1]);
      hRecRat[i][j]->Scale(scale[i]);
      if (i == 0) hRecRat[i][j]->Draw((j<4)?"e0":"histpx0");
      else hRecRat[i][j]->Draw((j<4)?"e0same":"histpx0same");
    }
    c->cd(6);
    for (Int_t i = 0; i < nFiles-1; i++) {
      hAccRat[i][j]->SetLineColor(color[i+1]);
      hAccRat[i][j]->SetMarkerStyle(markerStyle);
      hAccRat[i][j]->SetMarkerSize(markerSize[i+1]);
      hAccRat[i][j]->SetMarkerColor(color[i+1]);
      if (i == 0) hAccRat[i][j]->Draw((j<4)?"e0":"histpx0");
      else hAccRat[i][j]->Draw((j<4)?"e0same":"histpx0same");
    }
  }
  
  // draw histos vs pT differently for the note
  TH1 *h;
  TCanvas* cpT = new TCanvas("CompareAccEffVspT", "CompareAccEffVspT", 900, 600);
  cpT->Divide(2,2);
  cpT->cd(1);
  gPad->SetLogy();
  if (hDataGen[0][0]) {
    h = static_cast<TH1*>(hDataGen[0][0]->Clone());
    h->SetStats(0);
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetXaxis()->SetLabelSize(0.05);
    h->GetYaxis()->SetTitleSize(0.05);
    h->GetYaxis()->SetTitleOffset(0.8);
    h->GetYaxis()->SetLabelSize(0.05);
    h->Draw("e0");
    for (Int_t i = 1; i < nFiles; i++) hDataGen[i][0]->Draw("e0same");
  }
  cpT->cd(2);
  if (hDataAcc[0][4]) {
    h = static_cast<TH1*>(hDataAcc[0][4]->Clone());
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetXaxis()->SetTitleOffset(0.);
    h->GetXaxis()->SetLabelSize(0.05);
    h->GetYaxis()->SetTitleSize(0.05);
    h->GetYaxis()->SetLabelSize(0.05);
    h->Draw("e0");
    for (Int_t i = 1; i < nFiles; i++) hDataAcc[i][4]->Draw("e0same");
  }
  cpT->cd(3);
  if (hGenRat[0][0]) {
    h = static_cast<TH1*>(hGenRat[0][0]->Clone());
    h->SetStats(0);
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetXaxis()->SetLabelSize(0.05);
    h->GetYaxis()->SetTitleSize(0.05);
    h->GetYaxis()->SetTitleOffset(0.8);
    h->GetYaxis()->SetLabelSize(0.05);
    h->Draw("e0");
    for (Int_t i = 1; i < nFiles-1; i++) hGenRat[i][0]->Draw("e0same");
  }
  cpT->cd(4);
  gPad->SetGridy();
  if (hAccRat[0][4]) {
    h = static_cast<TH1*>(hAccRat[0][4]->Clone());
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetXaxis()->SetTitleOffset(0.);
    h->GetXaxis()->SetLabelSize(0.05);
    h->GetYaxis()->SetTitleSize(0.05);
    h->GetYaxis()->SetLabelSize(0.05);
    h->Draw("histpx0");
    for (Int_t i = 1; i < nFiles-1; i++) hAccRat[i][4]->Draw("histpx0same");
  }
  
  // draw histos vs y differently for the note
  TCanvas* cy = new TCanvas("CompareAccEffVsy", "CompareAccEffVsy", 900, 600);
  cy->Divide(2,2);
  cy->cd(1);
  if (hDataGen[0][1]) {
    h = static_cast<TH1*>(hDataGen[0][1]->Clone());
    h->SetStats(0);
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetXaxis()->SetLabelSize(0.05);
    h->GetYaxis()->SetTitleSize(0.05);
    h->GetYaxis()->SetTitleOffset(0.8);
    h->GetYaxis()->SetLabelSize(0.05);
    h->Draw("e0");
    for (Int_t i = 1; i < nFiles; i++) hDataGen[i][1]->Draw("e0same");
  }
  cy->cd(2);
  if (hDataAcc[0][5]) {
    h = static_cast<TH1*>(hDataAcc[0][5]->Clone());
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetXaxis()->SetTitleOffset(0.);
    h->GetXaxis()->SetLabelSize(0.05);
    h->GetYaxis()->SetTitleSize(0.05);
    h->GetYaxis()->SetLabelSize(0.05);
    h->Draw("e0");
    for (Int_t i = 1; i < nFiles; i++) hDataAcc[i][5]->Draw("e0same");
  }
  cy->cd(3);
  if (hGenRat[0][1]) {
    h = static_cast<TH1*>(hGenRat[0][1]->Clone());
    h->SetStats(0);
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetXaxis()->SetLabelSize(0.05);
    h->GetYaxis()->SetTitleSize(0.05);
    h->GetYaxis()->SetTitleOffset(0.8);
    h->GetYaxis()->SetLabelSize(0.05);
    h->Draw("e0");
    for (Int_t i = 1; i < nFiles-1; i++) hGenRat[i][1]->Draw("e0same");
  }
  cy->cd(4);
  gPad->SetGridy();
  if (hAccRat[0][5]) {
    h = static_cast<TH1*>(hAccRat[0][5]->Clone());
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetXaxis()->SetTitleOffset(0.);
    h->GetXaxis()->SetLabelSize(0.05);
    h->GetYaxis()->SetTitleSize(0.05);
    h->GetYaxis()->SetLabelSize(0.05);
    h->Draw("histpx0");
    for (Int_t i = 1; i < nFiles-1; i++) hAccRat[i][5]->Draw("histpx0same");
  }
  
  // draw histos vs pT/y differently for the note
  TCanvas* cpTy = new TCanvas("CompareAccEffVspTy", "CompareAccEffVspTy", 900, 300);
  cpTy->Divide(2,1);
  cpTy->cd(1);
  if (hDataAcc[0][6]) {
    h = static_cast<TH1*>(hDataAcc[0][6]->Clone());
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetXaxis()->SetTitleOffset(0.);
    h->GetXaxis()->SetLabelSize(0.05);
    h->GetYaxis()->SetTitleSize(0.05);
    h->GetYaxis()->SetLabelSize(0.05);
    h->Draw("e0");
    for (Int_t i = 1; i < nFiles; i++) hDataAcc[i][6]->Draw("e0same");
  }
  cpTy->cd(2);
  gPad->SetGridy();
  if (hAccRat[0][6]) {
    h = static_cast<TH1*>(hAccRat[0][6]->Clone());
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetXaxis()->SetTitleOffset(0.);
    h->GetXaxis()->SetLabelSize(0.05);
    h->GetYaxis()->SetTitleSize(0.05);
    h->GetYaxis()->SetLabelSize(0.05);
    h->Draw("histpx0");
    for (Int_t i = 1; i < nFiles-1; i++) hAccRat[i][6]->Draw("histpx0same");
  }
  
  // draw histos vs cent
  if (vsCent) {
    TGraphAsymmErrors *g;
    TCanvas* ccent = new TCanvas("CompareAccEffVsCent", "CompareAccEffVsCent", 900, 300);
    ccent->Divide(2,1);
    ccent->cd(1);
    for (Int_t i = 0; i < nFiles; i++) if (gDataAccCent[i]) {
      gDataAccCent[i]->SetLineColor(color[i]);
      gDataAccCent[i]->SetMarkerStyle(markerStyle);
      gDataAccCent[i]->SetMarkerSize(markerSize[i]);
      gDataAccCent[i]->SetMarkerColor(color[i]);
      if (i == 0) {
        gDataAccCent[0]->GetXaxis()->SetTitleSize(0.05);
        gDataAccCent[0]->GetXaxis()->SetTitleFont(42);
        gDataAccCent[0]->GetXaxis()->SetTitleOffset(0.8);
        gDataAccCent[0]->GetXaxis()->SetLabelSize(0.05);
        gDataAccCent[0]->GetXaxis()->SetLabelFont(42);
        gDataAccCent[0]->GetYaxis()->SetTitleSize(0.05);
        gDataAccCent[0]->GetYaxis()->SetLabelSize(0.05);
        gDataAccCent[0]->GetYaxis()->SetLabelFont(42);
        gDataAccCent[0]->Draw("ap");
      } else gDataAccCent[i]->Draw("p");
    }
    ccent->cd(2);
    gPad->SetGridy();
    for (Int_t i = 0; i < nFiles-1; i++) if (gDataAccRatCent[i]) {
      gDataAccRatCent[i]->SetLineColor(color[i+1]);
      gDataAccRatCent[i]->SetMarkerStyle(markerStyle);
      gDataAccRatCent[i]->SetMarkerSize(markerSize[i+1]);
      gDataAccRatCent[i]->SetMarkerColor(color[i+1]);
      if (i == 0) {
        gDataAccRatCent[0]->GetXaxis()->SetTitle("centrality");
        gDataAccRatCent[0]->GetXaxis()->SetTitleSize(0.05);
        gDataAccRatCent[0]->GetXaxis()->SetTitleOffset(0.8);
        gDataAccRatCent[0]->GetXaxis()->SetLabelSize(0.05);
        gDataAccRatCent[0]->GetYaxis()->SetTitle("ratio");
        gDataAccRatCent[0]->GetYaxis()->SetTitleSize(0.05);
        gDataAccRatCent[0]->GetYaxis()->SetLabelSize(0.05);
        gDataAccRatCent[0]->Draw("ap");
      } else gDataAccRatCent[i]->Draw("p");
    }
  }
  
  // draw histos vs cent vs pT
  if (vsCent && nDataCentVsPt > 0) {
    TGraphAsymmErrors *g;
    TCanvas* ccentpt = new TCanvas("CompareAccEffVsCentVspT", "CompareAccEffVsCentVspT", 1200, 600);
    ccentpt->Divide(nDataCentVsPt,2);
    for (Int_t j = 0; j < nDataCentVsPt; ++j) {
      ccentpt->cd(j+1);
      for (Int_t i = 0; i < nFiles; i++) if (gDataAccCentVsPt[i][j]) {
        gDataAccCentVsPt[i][j]->SetLineColor(color[i]);
        gDataAccCentVsPt[i][j]->SetMarkerStyle(markerStyle);
        gDataAccCentVsPt[i][j]->SetMarkerSize(markerSize[i]);
        gDataAccCentVsPt[i][j]->SetMarkerColor(color[i]);
        if (i == 0) {
          gDataAccCentVsPt[0][j]->GetXaxis()->SetTitleSize(0.05);
          gDataAccCentVsPt[0][j]->GetXaxis()->SetTitleFont(42);
          gDataAccCentVsPt[0][j]->GetXaxis()->SetTitleOffset(0.8);
          gDataAccCentVsPt[0][j]->GetXaxis()->SetLabelSize(0.05);
          gDataAccCentVsPt[0][j]->GetXaxis()->SetLabelFont(42);
          gDataAccCentVsPt[0][j]->GetYaxis()->SetTitleSize(0.05);
          gDataAccCentVsPt[0][j]->GetYaxis()->SetLabelSize(0.05);
          gDataAccCentVsPt[0][j]->GetYaxis()->SetLabelFont(42);
          gDataAccCentVsPt[0][j]->Draw("ap");
        } else gDataAccCentVsPt[i][j]->Draw("p");
      }
      ccentpt->cd(j+nDataCentVsPt+1);
      gPad->SetGridy();
      for (Int_t i = 0; i < nFiles-1; i++) if (gDataAccRatCentVsPt[i][j]) {
        gDataAccRatCentVsPt[i][j]->SetLineColor(color[i+1]);
        gDataAccRatCentVsPt[i][j]->SetMarkerStyle(markerStyle);
        gDataAccRatCentVsPt[i][j]->SetMarkerSize(markerSize[i+1]);
        gDataAccRatCentVsPt[i][j]->SetMarkerColor(color[i+1]);
        if (i == 0) {
          gDataAccRatCentVsPt[0][j]->GetXaxis()->SetTitle("centrality");
          gDataAccRatCentVsPt[0][j]->GetXaxis()->SetTitleSize(0.05);
          gDataAccRatCentVsPt[0][j]->GetXaxis()->SetTitleOffset(0.8);
          gDataAccRatCentVsPt[0][j]->GetXaxis()->SetLabelSize(0.05);
          gDataAccRatCentVsPt[0][j]->GetYaxis()->SetTitle("ratio");
          gDataAccRatCentVsPt[0][j]->GetYaxis()->SetTitleSize(0.05);
          gDataAccRatCentVsPt[0][j]->GetYaxis()->SetLabelSize(0.05);
          gDataAccRatCentVsPt[0][j]->Draw("ap");
        } else gDataAccRatCentVsPt[i][j]->Draw("p");
      }
    }
  }
  
  // draw histos vs cent vs y
  if (vsCent && nDataCentVsY > 0) {
    TGraphAsymmErrors *g;
    TCanvas* ccenty = new TCanvas("CompareAccEffVsCentVsY", "CompareAccEffVsCentVsY", 1200, 600);
    ccenty->Divide(nDataCentVsY,2);
    for (Int_t j = 0; j < nDataCentVsY; ++j) {
      ccenty->cd(j+1);
      for (Int_t i = 0; i < nFiles; i++) if (gDataAccCentVsY[i][j]) {
        gDataAccCentVsY[i][j]->SetLineColor(color[i]);
        gDataAccCentVsY[i][j]->SetMarkerStyle(markerStyle);
        gDataAccCentVsY[i][j]->SetMarkerSize(markerSize[i]);
        gDataAccCentVsY[i][j]->SetMarkerColor(color[i]);
        if (i == 0) {
          gDataAccCentVsY[0][j]->GetXaxis()->SetTitleSize(0.05);
          gDataAccCentVsY[0][j]->GetXaxis()->SetTitleFont(42);
          gDataAccCentVsY[0][j]->GetXaxis()->SetTitleOffset(0.8);
          gDataAccCentVsY[0][j]->GetXaxis()->SetLabelSize(0.05);
          gDataAccCentVsY[0][j]->GetXaxis()->SetLabelFont(42);
          gDataAccCentVsY[0][j]->GetYaxis()->SetTitleSize(0.05);
          gDataAccCentVsY[0][j]->GetYaxis()->SetLabelSize(0.05);
          gDataAccCentVsY[0][j]->GetYaxis()->SetLabelFont(42);
          gDataAccCentVsY[0][j]->Draw("ap");
        } else gDataAccCentVsY[i][j]->Draw("p");
      }
      ccenty->cd(j+nDataCentVsY+1);
      gPad->SetGridy();
      for (Int_t i = 0; i < nFiles-1; i++) if (gDataAccRatCentVsY[i][j]) {
        gDataAccRatCentVsY[i][j]->SetLineColor(color[i+1]);
        gDataAccRatCentVsY[i][j]->SetMarkerStyle(markerStyle);
        gDataAccRatCentVsY[i][j]->SetMarkerSize(markerSize[i+1]);
        gDataAccRatCentVsY[i][j]->SetMarkerColor(color[i+1]);
        if (i == 0) {
          gDataAccRatCentVsY[0][j]->GetXaxis()->SetTitle("centrality");
          gDataAccRatCentVsY[0][j]->GetXaxis()->SetTitleSize(0.05);
          gDataAccRatCentVsY[0][j]->GetXaxis()->SetTitleOffset(0.8);
          gDataAccRatCentVsY[0][j]->GetXaxis()->SetLabelSize(0.05);
          gDataAccRatCentVsY[0][j]->GetYaxis()->SetTitle("ratio");
          gDataAccRatCentVsY[0][j]->GetYaxis()->SetTitleSize(0.05);
          gDataAccRatCentVsY[0][j]->GetYaxis()->SetLabelSize(0.05);
          gDataAccRatCentVsY[0][j]->Draw("ap");
        } else gDataAccRatCentVsY[i][j]->Draw("p");
      }
    }
  }
  
}

//--------------------------------------------------------------------
void FillRatio(TGraphAsymmErrors *g, TGraphAsymmErrors *ref, TGraphAsymmErrors *rat)
{
  /// make the ratio g/ref
  
  Int_t n = (g) ? g->GetN() : -1;
  
  if (!g || !ref || ref->GetN() != n || !rat || rat->GetN() != n) {
    printf("inconsistent graphs\n");
    return;
  }
  
  Double_t x, acc[2];
  for (Int_t i = 0; i < n; ++i) {
    g->GetPoint(i,x,acc[1]);
    ref->GetPoint(i,x,acc[0]);
    rat->SetPoint(i,x,(acc[0] > 0.) ? acc[1]/acc[0] : 0.);
    rat->SetPointError(i,0.,0.,g->GetErrorYlow(i),g->GetErrorYhigh(i)); // take one but should partly cancel
  }
  
  rat->GetXaxis()->Set(n, -0.5, n-0.5);
  rat->GetXaxis()->SetNdivisions(1,kFALSE);
  for (Int_t i = 0; i < n; ++i) {
    rat->GetXaxis()->SetBinLabel(i+1, g->GetXaxis()->GetBinLabel(i+1));
  }

}


