/*
 *  CompareTrack.C
 *  aliroot_dev
 *
 *  Created by philippe pillot on 18/04/13.
 *  Copyright 2013 SUBATECH. All rights reserved.
 *
 */

void CompareTrack(TString fileName_before = "AnalysisResults.root", TString fileName_after = "AnalysisResults.root")
{
  /// compare selected tracks versus run
  
  // prepare environment
  gROOT->LoadMacro("$ALICE/Macros/Facilities/runTaskFacilities.C");
  LoadAlirootLocally("", "", "");
  
  // open files
  TFile* file_before = TFile::Open(fileName_before.Data(),"READ");
  if (!file_before || !file_before->IsOpen()) return;
  TFile* file_after = TFile::Open(fileName_after.Data(),"READ");
  if (!file_after || !file_after->IsOpen()) return;
  
  // counters before
  AliCounterCollection* counters_before = static_cast<AliCounterCollection*>(file_before->FindObjectAny("trackCounters"));
  if (!counters_before) return;
  
  // counters after
  AliCounterCollection* counters_after = static_cast<AliCounterCollection*>(file_after->FindObjectAny("trackCounters"));
  if (!counters_after) return;
  
  // histos before
  counters_before->Sort("run",1);
  TH1D *b = counters_before->Draw("run","track:matched/trig:low/rabs:yes/eta:yes/pdca:yes/pt:1gev");
  
  // histos after
  counters_after->Sort("run",1);
  TH1D *a = counters_after->Draw("run","track:matched/trig:low/rabs:yes/eta:yes/pdca:yes/pt:1gev");
  
  // ratio
  TH1D *r = a->Clone("ratio");
  r->Sumw2();
  r->Divide(b);
  
  // draw
  TCanvas* cRatio = new TCanvas("cRatio", "cRatio", 1000, 400);
  cRatio->Divide(1,2);
  cRatio->cd(1);
  b->SetLineColor(2);
  b->Draw("h");
  a->Draw("hsame");
  cRatio->cd(2);
  r->Draw("e0");
  
}

