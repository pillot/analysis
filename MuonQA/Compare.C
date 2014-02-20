/*
 *  Compare.C
 *  aliroot_dev
 *
 *  Created by philippe pillot on 23/02/11.
 *  Copyright 2011 SUBATECH. All rights reserved.
 *
 */

void Compare(TString fileName_before = "AnalysisResults.root", TString fileName_after = "AnalysisResults.root")
{
  /// compare results before/after refit if fileName_before==fileName_after
  /// or after refit in both files if they are different
  
  gStyle->SetFillColor(0);
  
  // prepare environment
  gROOT->LoadMacro("/Users/pillot/Work/Alice/Work/Macros/Facilities/runTaskFacilities.C");
  LoadAlirootLocally("","","");
  
  
  // open files
  TFile* file_before = TFile::Open(fileName_before.Data(),"READ");
  if (!file_before || !file_before->IsOpen()) return;
  TFile* file_after = TFile::Open(fileName_after.Data(),"READ");
  if (!file_after || !file_after->IsOpen()) return;
  
  // histo before
  TObjArray* histo_before = static_cast<TObjArray*>(file_before->FindObjectAny("general1"));
  if (!histo_before) return;
  TObjArray* histo2_before = static_cast<TObjArray*>(file_before->FindObjectAny("general2"));
  if (!histo2_before) return;
  
  // histo after
  TObjArray* histo_after = static_cast<TObjArray*>(file_after->FindObjectAny("general1"));
  if (!histo_after) return;
  TObjArray* histo2_after = static_cast<TObjArray*>(file_after->FindObjectAny("general2"));
  if (!histo2_after) return;
  
  // counters before
  AliCounterCollection* counters_before = static_cast<AliCounterCollection*>(file_before->FindObjectAny("eventCounters"));
  if (!counters_before) return;
  counters_before->Sort();
  AliCounterCollection* tracks_before = static_cast<AliCounterCollection*>(file_before->FindObjectAny("trackCounters"));
  if (!tracks_before) return;
  tracks_before->Sort();
  
  // counters after
  AliCounterCollection* counters_after = static_cast<AliCounterCollection*>(file_after->FindObjectAny("eventCounters"));
  if (!counters_after) return;
  counters_after->Sort();
  AliCounterCollection* tracks_after = static_cast<AliCounterCollection*>(file_after->FindObjectAny("trackCounters"));
  if (!tracks_after) return;
  tracks_after->Sort();
  
  // number of events for normalization
  Double_t nEvent_before = counters_before->GetSum("event:any/selected:any");
  //TH1F* hChi2_before = static_cast<TH1F*>(histo_before->FindObject("hChi2"));
  //Double_t nEvent_before = hChi2_before->GetEntries();
  cout<<"- nEvent in file 1 = "<<nEvent_before<<endl;
  Double_t nEvent_after = counters_after->GetSum("event:any/selected:any");
  //TH1F* hChi2_after = static_cast<TH1F*>(histo_after->FindObject("hChi2"));
  //Double_t nEvent_after = hChi2_after->GetEntries();
  cout<<"- nEvent in file 2 = "<<nEvent_after<<endl;
  Double_t norm = nEvent_before/nEvent_after;
  cout<<"--> normalization factor = "<<norm<<endl;
  
  // draw plots
  TString sHist[6] = {"hChi2", "hNClustersPerTrack", "hNChamberHitPerTrack", "hDCA", "hPt", "hRapidity"};
  TCanvas* cHist = new TCanvas("cHist", "cHist", 1000, 800);
  cHist->Divide(3,2);
  for (Int_t i=0; i<6; i++) {
    cHist->cd(i+1);
    gPad->SetLogy();
    TH1* h_before = static_cast<TH1*>(histo_before->FindObject(sHist[i].Data()));
    TH1* h_after = static_cast<TH1*>(histo_after->FindObject(sHist[i].Data())));
    h_after->Scale(norm);
    h_before->Draw();
    h_before->SetLineColor(4);
    h_after->Draw("sames");
    h_after->SetLineColor(2);
  }
  
  TString sHist2[5] = {"hClusterChargePerChMean", "hClusterSizePerChMean", "hNClustersPerCh", "hClusterChargePerChSigma", "hClusterSizePerChSigma"};
  TCanvas* cHist2 = new TCanvas("cHist2", "cHist2", 1000, 800);
  cHist2->Divide(3,2);
  for (Int_t i=0; i<5; i++) {
    cHist2->cd(i+1);
    TH1* h_before = static_cast<TH1*>(histo2_before->FindObject(sHist2[i].Data()));
    TH1* h_after = static_cast<TH1*>(histo2_after->FindObject(sHist2[i].Data())));
    h_before->Draw();
    h_before->SetLineColor(4);
    h_after->Draw("sames");
    h_after->SetLineColor(2);
    h_after->SetMarkerColor(2);
  }
  
  TString sHist3[6] = {"track", "trigger", "v0mult", "track", "trigger", "v0mult"};
  TString sSel3[6] = {"", "track:trackeronly,matched", "track:trackeronly,matched", "selected:yes", "track:trackeronly,matched/selected:yes", "track:trackeronly,matched/selected:yes"};
  TCanvas* cHist3 = new TCanvas("cHist3", "cHist3", 1200, 600);
  cHist3->Divide(3,2);
  for (Int_t i=0; i<6; i++) {
    cHist3->cd(i+1);
    gPad->SetLogy();
    TH1* h_before = tracks_before->Get(sHist3[i].Data(),sSel3[i].Data());
    TH1* h_after = tracks_after->Get(sHist3[i].Data(),sSel3[i].Data());
    h_after->Scale(norm);
    h_before->Draw();
    h_before->SetLineColor(4);
    h_after->Draw("sames");
    h_after->SetLineColor(2);
    h_after->SetMarkerColor(2);
  }
  
  TString sHist4[6] = {"trigger", "trigger", "v0mult", "trigger", "trigger", "v0mult"};
  TString sSel4[6] = {"", "event:muon", "", "selected:yes", "event:muon/selected:yes", "selected:yes"};
  TCanvas* cHist4 = new TCanvas("cHist4", "cHist4", 1200, 600);
  cHist4->Divide(3,2);
  for (Int_t i=0; i<6; i++) {
    cHist4->cd(i+1);
    gPad->SetLogy();
    TH1* h_before = counters_before->Get(sHist4[i].Data(),sSel4[i].Data());
    TH1* h_after = counters_after->Get(sHist4[i].Data(),sSel4[i].Data());
    h_after->Scale(norm);
    h_before->Draw();
    h_before->SetLineColor(4);
    h_after->Draw("sames");
    h_after->SetLineColor(2);
    h_after->SetMarkerColor(2);
  }
  
}

