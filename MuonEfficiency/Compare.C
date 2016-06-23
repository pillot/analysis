/*
 *  Compare.C
 *  aliroot_dev
 *
 *  Created by philippe pillot on 13/09/11.
 *  Copyright 2011 SUBATECH. All rights reserved.
 *
 */

void Compare(TString fileName0, TString fileNames)
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
  }
  
  // compute ratio
  TH1F* hGenRat[1000][nData];
  TH1F* hRecRat[1000][nData];
  TH1F* hAccRat[1000][nData];
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
  TCanvas* cpT = new TCanvas("CompareAccEffVspT", "CompareAccEffVspT", 900, 600);
  cpT->Divide(2,2);
  cpT->cd(1);
  gPad->SetLogy();
  if (!hDataGen[0][0]) continue;
  TH1 *h = static_cast<TH1*>(hDataGen[0][0]->Clone());
  h->SetStats(0);
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetXaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetTitleOffset(0.8);
  h->GetYaxis()->SetLabelSize(0.05);
  h->Draw("e0");
  for (Int_t i = 1; i < nFiles; i++) hDataGen[i][0]->Draw("e0same");
  cpT->cd(2);
  if (!hDataAcc[0][4]) continue;
  h = static_cast<TH1*>(hDataAcc[0][4]->Clone());
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetXaxis()->SetTitleOffset(0.);
  h->GetXaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetLabelSize(0.05);
  h->Draw("e0");
  for (Int_t i = 1; i < nFiles; i++) hDataAcc[i][4]->Draw("e0same");
  cpT->cd(3);
  if (!hGenRat[0][0]) continue;
  h = static_cast<TH1*>(hGenRat[0][0]->Clone());
  h->SetStats(0);
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetXaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetTitleOffset(0.8);
  h->GetYaxis()->SetLabelSize(0.05);
  h->Draw("e0");
  for (Int_t i = 1; i < nFiles-1; i++) hGenRat[i][0]->Draw("e0same");
  cpT->cd(4);
  gPad->SetGridy();
  if (!hAccRat[0][4]) continue;
  h = static_cast<TH1*>(hAccRat[0][4]->Clone());
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetXaxis()->SetTitleOffset(0.);
  h->GetXaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetLabelSize(0.05);
  h->Draw("histpx0");
  for (Int_t i = 1; i < nFiles-1; i++) hAccRat[i][4]->Draw("histpx0same");
  
  // draw histos vs y differently for the note
  TCanvas* cy = new TCanvas("CompareAccEffVsy", "CompareAccEffVsy", 900, 600);
  cy->Divide(2,2);
  cy->cd(1);
  if (!hDataGen[0][1]) continue;
  TH1 *h = static_cast<TH1*>(hDataGen[0][1]->Clone());
  h->SetStats(0);
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetXaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetTitleOffset(0.8);
  h->GetYaxis()->SetLabelSize(0.05);
  h->Draw("e0");
  for (Int_t i = 1; i < nFiles; i++) hDataGen[i][1]->Draw("e0same");
  cy->cd(2);
  if (!hDataAcc[0][5]) continue;
  h = static_cast<TH1*>(hDataAcc[0][5]->Clone());
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetXaxis()->SetTitleOffset(0.);
  h->GetXaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetLabelSize(0.05);
  h->Draw("e0");
  for (Int_t i = 1; i < nFiles; i++) hDataAcc[i][5]->Draw("e0same");
  cy->cd(3);
  if (!hGenRat[0][1]) continue;
  h = static_cast<TH1*>(hGenRat[0][1]->Clone());
  h->SetStats(0);
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetXaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetTitleOffset(0.8);
  h->GetYaxis()->SetLabelSize(0.05);
  h->Draw("e0");
  for (Int_t i = 1; i < nFiles-1; i++) hGenRat[i][1]->Draw("e0same");
  cy->cd(4);
  gPad->SetGridy();
  if (!hAccRat[0][5]) continue;
  h = static_cast<TH1*>(hAccRat[0][5]->Clone());
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetXaxis()->SetTitleOffset(0.);
  h->GetXaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetLabelSize(0.05);
  h->Draw("histpx0");
  for (Int_t i = 1; i < nFiles-1; i++) hAccRat[i][5]->Draw("histpx0same");
  
}

