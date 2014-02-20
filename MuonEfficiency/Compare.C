/*
 *  Compare.C
 *  aliroot_dev
 *
 *  Created by philippe pillot on 13/09/11.
 *  Copyright 2011 SUBATECH. All rights reserved.
 *
 */

void Compare(TString fileName1, TString fileName2, TString fileName3 = "", TString fileName4 = "", TString fileName5 = "")
{
  /// compare Acc*Eff results between file1 and file2
  
  gStyle->SetOptStat(11);
  gStyle->SetFillColor(0);
  
  TString fileName[5] = {fileName1.Data(), fileName2.Data(), fileName3.Data(), fileName4.Data(), fileName5.Data()};
  Int_t nFiles = 2;
  if (!fileName3.IsNull()) {
    nFiles++;
    if (!fileName4.IsNull()) {
      nFiles++;
      if (!fileName5.IsNull()) {
	nFiles++;
      }
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
  TFile* outFile[5];
  TH1F* hDataGen[5][nData];
  TH1F* hDataRec[5][nData];
  TH1F* hDataAcc[5][nData];
  for (Int_t i = 0; i < nFiles; i++) {
    outFile[i] = TFile::Open(fileName[i].Data(),"READ");
    if (!outFile[i] || !outFile[i]->IsOpen()) return;
    for (Int_t j = 0; j < nData; j++) {
      hDataGen[i][j] = static_cast<TH1F*>(outFile[i]->FindObjectAny(dataGen[j].Data()));
      hDataRec[i][j] = static_cast<TH1F*>(outFile[i]->FindObjectAny(dataRec[j].Data()));
      hDataAcc[i][j] = static_cast<TH1F*>(outFile[i]->FindObjectAny(dataAcc[j].Data()));
    }
  }
  
  // compute ratio
  TH1F* hGenRat[4][nData];
  TH1F* hRecRat[4][nData];
  TH1F* hAccRat[4][nData];
  for (Int_t i = 1; i < nFiles; i++) {
    for (Int_t j = 0; j < nData; j++) {
      if (!hDataGen[i][j]) continue;
      hGenRat[i-1][j] = (TH1F*) hDataGen[i][j]->Clone(Form("%sRat",dataGen[j].Data()));
      hGenRat[i-1][j]->SetTitle("ratio");
      hGenRat[i-1][j]->Sumw2();
      hGenRat[i-1][j]->Divide(hDataGen[0][j]);
      hRecRat[i-1][j] = (TH1F*) hDataRec[i][j]->Clone(Form("%sRat",dataRec[j].Data()));
      hRecRat[i-1][j]->SetTitle("ratio");
      hRecRat[i-1][j]->Sumw2();
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
  Double_t scale[4] = {1., 1., 1., 1.};
  Int_t color[4] = {2, 4, 3, 6};
  for (Int_t j = 0; j < nData; j++) {
    TCanvas* c = new TCanvas(cName[j].Data(), cTitle[j].Data(), 900, 600);
    c->Divide(3,2);
    c->cd(1);
    gPad->SetLogy();
    if (!hDataGen[0][j]) continue;
    hDataGen[0][j]->SetLineColor(1);
    hDataGen[0][j]->Draw();
    for (Int_t i = 1; i < nFiles; i++) {
      hDataGen[i][j]->Sumw2();
      if (j < 4) scale[i-1] = ((Double_t)hDataGen[0][j]->GetEntries())/((Double_t)hDataGen[i][j]->GetEntries());
      else scale[i-1] = ((Double_t)hDataGen[0][j]->GetBinContent(1))/((Double_t)hDataGen[i][j]->GetBinContent(1));
      hDataGen[i][j]->Scale(scale[i-1]);
      hDataGen[i][j]->SetLineColor(color[i-1]);
      hDataGen[i][j]->Draw("histsame");
    }
    c->cd(2);
    gPad->SetLogy();
    hDataRec[0][j]->SetLineColor(1);
    hDataRec[0][j]->Draw();
    for (Int_t i = 1; i < nFiles; i++) {
      hDataRec[i][j]->Sumw2();
      hDataRec[i][j]->Scale(scale[i-1]);
      hDataRec[i][j]->SetLineColor(color[i-1]);
      hDataRec[i][j]->Draw("histsame");
    }
    c->cd(3);
    hDataAcc[0][j]->SetLineColor(1);
    hDataAcc[0][j]->SetMarkerColor(1);
    hDataAcc[0][j]->Draw("e0");
    for (Int_t i = 1; i < nFiles; i++) {
      hDataAcc[i][j]->SetLineColor(color[i-1]);
      hDataAcc[i][j]->SetMarkerColor(color[i-1]);
      hDataAcc[i][j]->Draw("e0same");
    }
    c->cd(4);
    for (Int_t i = 0; i < nFiles-1; i++) {
      hGenRat[i][j]->SetLineColor(color[i]);
      hGenRat[i][j]->SetMarkerStyle(21);
      hGenRat[i][j]->SetMarkerSize(0.6);
      hGenRat[i][j]->SetMarkerColor(color[i]);
      hGenRat[i][j]->Scale(scale[i]);
      if (i == 0) hGenRat[i][j]->Draw("e0");
      else hGenRat[i][j]->Draw("e0same");
    }
    c->cd(5);
    for (Int_t i = 0; i < nFiles-1; i++) {
      hRecRat[i][j]->SetLineColor(color[i]);
      hRecRat[i][j]->SetMarkerStyle(21);
      hRecRat[i][j]->SetMarkerSize(0.6);
      hRecRat[i][j]->SetMarkerColor(color[i]);
      hRecRat[i][j]->Scale(scale[i]);
      if (i == 0) hRecRat[i][j]->Draw("e0");
      else hRecRat[i][j]->Draw("e0same");
    }
    c->cd(6);
    for (Int_t i = 0; i < nFiles-1; i++) {
      hAccRat[i][j]->SetLineColor(color[i]);
      hAccRat[i][j]->SetMarkerStyle(21);
      hAccRat[i][j]->SetMarkerSize(0.6);
      hAccRat[i][j]->SetMarkerColor(color[i]);
      if (i == 0) hAccRat[i][j]->Draw("e0");
      else hAccRat[i][j]->Draw("e0same");
    }
  }
  
}

