/*
 *  BasicCompare.C
 *  aliroot_dev
 *
 *  Created by philippe pillot on 09/09/11.
 *  Copyright 2011 SUBATECH. All rights reserved.
 *
 */

void BasicCompare(TString fileName1 = "AnalysisResults.root", TString fileName2 = "AnalysisResults.root")
{
  /// compare results before/after refit if fileName_before==fileName_after
  /// or after refit in both files if they are different
  
  // prepare environment
  gROOT->LoadMacro("$ALICE/Macros/Facilities/runTaskFacilities.C");
  LoadAlirootLocally("PWG3base", "", "");
  
  // open files
  TFile* file1 = TFile::Open(fileName1.Data(),"READ");
  if (!file1 || !file1->IsOpen()) return;
  TFile* file2 = TFile::Open(fileName2.Data(),"READ");
  if (!file2 || !file2->IsOpen()) return;
  
  // histo 1
  TObjArray* histo1 = static_cast<TObjArray*>(file1->FindObjectAny("Histograms"));
  if (!histo1) return;
  
  // histo 2
  TObjArray* histo2 = static_cast<TObjArray*>(file2->FindObjectAny("Histograms"));
  if (!histo2) return;
  
  // draw differences
  TString sHist[11] = {"hChi2", "hNClustersPerTrack", "hNChamberHitPerTrack", "hDCA", "hPt", "hRapidity", "hMass", "hPUncorrected", "hRAbs", "hDCAX", "hDCAY"};
  TCanvas* cHist = new TCanvas("cHist", "cHist", 1000, 800);
  cHist->Divide(4,3);
  TCanvas* cDiff = new TCanvas("cDiff", "cDiff", 1000, 800);
  cDiff->Divide(4,3);
  TCanvas* cRatio = new TCanvas("cRatio", "cRatio", 1000, 800);
  cRatio->Divide(4,3);
  for (Int_t i=0; i<11; i++) {
    cHist->cd(i+1);
    gPad->SetLogy();
    TH1F* h1 = static_cast<TH1F*>(histo1->FindObject(sHist[i].Data()));
    if (!h1) continue;
    TH1F* h2 = static_cast<TH1F*>(histo2->FindObject(sHist[i].Data()));
    h1->Scale(1./h1->GetEntries());
    h1->Draw();
    h2->Scale(1./h2->GetEntries());
    h2->Draw("sames");
    h2->SetLineColor(2);
    cDiff->cd(i+1);
    TH1F* h_diff = h1->Clone();
    h_diff->Add(h2, -1.);
    h_diff->Draw();
    h_diff->SetLineColor(4);
    cRatio->cd(i+1);
    TH1F* h_ratio = h1->Clone();
    h_ratio->Divide(h2);
    h_ratio->Draw();
    h_ratio->SetLineColor(4);
  }
  
}

