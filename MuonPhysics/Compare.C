/*
 *  Compare.C
 *  aliroot
 *
 *  Created by Philippe Pillot on 03/11/10.
 *  Copyright 2010 SUBATECH. All rights reserved.
 *
 */

void Compare(TString fileName_before = "AnalysisResults.root", TString fileName_after = "AnalysisResults.root")
{
  /// compare results before/after refit if fileName_before==fileName_after
  /// or after refit in both files if they are different
  
  gROOT->SetStyle("Plain");
  gROOT->ForceStyle();
  gStyle->SetOptFit(1);
  
  // prepare environment
  TString extraLibs = "RAWDatabase:CDB:STEER:MUONcore:MUONmapping:MUONcalib:MUONgeometry:MUONtrigger:MUONraw:MUONbase:MUONrec";
  LoadAlirootLocally(extraLibs);
  
  Bool_t sameFile = (fileName_before==fileName_after);
//  TString extBefore = sameFile ? "_before" : "_after";
  TString extBefore = "_phys";
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
  
  // counters before
  AliCounterCollection* counters_before = static_cast<AliCounterCollection*>(file_before->FindObjectAny(Form("trackCounters%s",extBefore.Data())));
  if (!counters_before) return;
  
  // counters after
  AliCounterCollection* counters_after = static_cast<AliCounterCollection*>(file_after->FindObjectAny(Form("trackCounters%s",extAfter.Data())));
  if (!counters_after) return;
  
  // number of events before
  AliCounterCollection* event_before = static_cast<AliCounterCollection*>(file_before->FindObjectAny(Form("eventCounters%s",extBefore.Data())));
  if (!event_before) return;
  Double_t nEventsBefore = event_before->GetSum("event:all");
  printf("\n# events before = %d\n",(Int_t)nEventsBefore);
  
  // number of events after
  AliCounterCollection* event_after = static_cast<AliCounterCollection*>(file_after->FindObjectAny(Form("eventCounters%s",extAfter.Data())));
  if (!event_after) return;
  Double_t nEventsAfter = event_after->GetSum("event:all");
  printf("# events after = %d\n",(Int_t)nEventsAfter);
  
  // draw differences
  const Int_t nHists = 16;
  TString sHist[nHists] = {"hChi2", "hNClustersPerTrack", "hNChamberHitPerTrack", "hNTracks", "hDCA", "hPt", "hRapidity", "hNTrigAll", "hMass", "hPUncorrected", "hRAbs", "hNMatchTracks", "hPtMuPlus", "hPtMuMinus", "hChi2Trig", "hMultSelect"};
  TCanvas* cHist = new TCanvas("cHist", "cHist", 1200, 1100);
  cHist->Divide(4,4);
  TLegend *lHist = new TLegend(0.5,0.65,0.9,0.8);
  lHist->SetFillStyle(0);
  lHist->SetBorderSize(0);
  TCanvas* cDiff = new TCanvas("cDiff", "cDiff", 1200, 1100);
  cDiff->Divide(4,4);
  TLegend *lDiff = new TLegend(0.5,0.65,0.92,0.75);
  lDiff->SetFillStyle(0);
  lDiff->SetBorderSize(0);
  TCanvas* cRatio = new TCanvas("cRatio", "cRatio", 1200, 1100);
  cRatio->Divide(4,4);
  TLegend *lRatio = new TLegend(0.5,0.65,0.9,0.75);
  lRatio->SetFillStyle(0);
  lRatio->SetBorderSize(0);
  for (Int_t i=0; i<nHists; i++) {
    cHist->cd(i+1);
    gPad->SetLogy();
    TH1F* h_before = static_cast<TH1F*>(histo_before->FindObject(sHist[i].Data()));
    if (!h_before) continue;
    //if (i == 15) h_before->Rebin(10);
    //h_before->Scale(1./nEventsBefore);
    TH1F* h_after = static_cast<TH1F*>(histo_after->FindObject(sHist[i].Data())));
    if (!h_after) continue;
    //if (i == 15) h_after->Rebin(10);
    //h_after->Scale(1./nEventsAfter);
    h_before->Draw();
    h_after->Draw("sames");
    h_after->SetLineColor(2);
    if (i == 0) {
      lHist->AddEntry(h_before,"before","l");
      lHist->AddEntry(h_after,"after","l");
      lHist->Draw("same");
    }
    cDiff->cd(i+1);
    TH1F* h_diff = h_before->Clone();
    h_diff->Add(h_after, -1.);
    h_diff->Draw();
    h_diff->SetLineColor(4);
    if (i == 0) {
      lDiff->AddEntry(h_diff,"before - after","l");
      lDiff->Draw("same");
    }
    cRatio->cd(i+1);
    TH1F* h_ratio = h_after->Clone();
    h_ratio->Divide(h_before);
    h_ratio->Draw();
    h_ratio->SetLineColor(4);
    if (i == 0) {
      lRatio->AddEntry(h_ratio,"after / before","l");
      lRatio->Draw("same");
    }
  }
  
  // print statistic
  Double_t b,a;
  printf("\n");
  printf("track       before     after     diff\n");
  b = counters_before->GetSum(); a = counters_after->GetSum();
  printf("total:   %9d %9d %7.2f %%\n", (Int_t)b, (Int_t)a, 100.*(a-b)/b);
  b = counters_before->GetSum("track:matched"); a = counters_after->GetSum("track:matched");
  printf("match:   %9d %9d %7.2f %%\n", (Int_t)b, (Int_t)a, 100.*(a-b)/b);
  b = counters_before->GetSum("track:matched/rabs:yes/eta:yes"); a = counters_after->GetSum("track:matched/rabs:yes/eta:yes");
  printf("+ acc:   %9d %9d %7.2f %%\n", (Int_t)b, (Int_t)a, 100.*(a-b)/b);
  b = counters_before->GetSum("track:matched/rabs:yes/eta:yes/pdca:yes"); a = counters_after->GetSum("track:matched/rabs:yes/eta:yes/pdca:yes");
  printf("+ DCA:   %9d %9d %7.2f %%\n", (Int_t)b, (Int_t)a, 100.*(a-b)/b);
  b = counters_before->GetSum("track:matched/rabs:yes/eta:yes/pdca:yes/pt:0.5GeV"); a = counters_after->GetSum("track:matched/rabs:yes/eta:yes/pdca:yes/pt:0.5GeV");
  printf(">0.5GeV: %9d %9d %7.2f %%\n", (Int_t)b, (Int_t)a, 100.*(a-b)/b);
  b = counters_before->GetSum("track:matched/rabs:yes/eta:yes/pdca:yes/pt:1GeV"); a = counters_after->GetSum("track:matched/rabs:yes/eta:yes/pdca:yes/pt:1GeV");
  printf(">1GeV:   %9d %9d %7.2f %%\n", (Int_t)b, (Int_t)a, 100.*(a-b)/b);
  b = counters_before->GetSum("track:matched/rabs:yes/eta:yes/pdca:yes/pt:2GeV"); a = counters_after->GetSum("track:matched/rabs:yes/eta:yes/pdca:yes/pt:2GeV");
  printf(">2GeV:   %9d %9d %7.2f %%\n", (Int_t)b, (Int_t)a, 100.*(a-b)/b);
  b = counters_before->GetSum("track:matched/rabs:yes/eta:yes/pdca:yes/trig:low"); a = counters_after->GetSum("track:matched/rabs:yes/eta:yes/pdca:yes/trig:low");
  printf("lowPt:   %9d %9d %7.2f %%\n", (Int_t)b, (Int_t)a, 100.*(a-b)/b);
  b = counters_before->GetSum("track:matched/rabs:yes/eta:yes/pdca:yes/trig:high"); a = counters_after->GetSum("track:matched/rabs:yes/eta:yes/pdca:yes/trig:high");
  printf("highPt:  %9d %9d %7.2f %%\n", (Int_t)b, (Int_t)a, 100.*(a-b)/b);
  printf("\n");
  /*
  // print statistic
  Double_t b,a;
  printf("\n");
  printf("track       before     after     diff\n");
  b = counters_before->GetSum(); a = counters_after->GetSum();
  printf("total:   %9d %9d %7.2f %%\n", (Int_t)b, (Int_t)a, 100.*(a-b)/b);
  b = counters_before->GetSum("track:matched"); a = counters_after->GetSum("track:matched");
  printf("match:   %9d %9d %7.2f %%\n", (Int_t)b, (Int_t)a, 100.*(a-b)/b);
  b = counters_before->GetSum("track:matched/rabs:yes/eta:yes"); a = counters_after->GetSum("track:matched/rabs:yes/eta:yes");
  printf("+ acc:   %9d %9d %7.2f %%\n", (Int_t)b, (Int_t)a, 100.*(a-b)/b);
  b = counters_before->GetSum("track:matched/rabs:yes/eta:yes/pdca:yes"); a = counters_after->GetSum("track:matched/rabs:yes/eta:yes/pdca:yes");
  printf("+ DCA:   %9d %9d %7.2f %%\n", (Int_t)b, (Int_t)a, 100.*(a-b)/b);
  b = counters_before->GetSum("track:matched/rabs:yes/eta:yes/pdca:yes/chi2:yes"); a = counters_after->GetSum("track:matched/rabs:yes/eta:yes/pdca:yes/chi2:yes");
  printf("+ chi2:   %9d %9d %7.2f %%\n", (Int_t)b, (Int_t)a, 100.*(a-b)/b);
  b = counters_before->GetSum("track:matched/rabs:yes/eta:yes/pdca:yes/chi2:yes/pt:0.5GeV"); a = counters_after->GetSum("track:matched/rabs:yes/eta:yes/pdca:yes/chi2:yes/pt:0.5GeV");
  printf(">0.5GeV: %9d %9d %7.2f %%\n", (Int_t)b, (Int_t)a, 100.*(a-b)/b);
  b = counters_before->GetSum("track:matched/rabs:yes/eta:yes/pdca:yes/chi2:yes/pt:1GeV"); a = counters_after->GetSum("track:matched/rabs:yes/eta:yes/pdca:yes/chi2:yes/pt:1GeV");
  printf(">1GeV:   %9d %9d %7.2f %%\n", (Int_t)b, (Int_t)a, 100.*(a-b)/b);
  b = counters_before->GetSum("track:matched/rabs:yes/eta:yes/pdca:yes/chi2:yes/pt:2GeV"); a = counters_after->GetSum("track:matched/rabs:yes/eta:yes/pdca:yes/chi2:yes/pt:2GeV");
  printf(">2GeV:   %9d %9d %7.2f %%\n", (Int_t)b, (Int_t)a, 100.*(a-b)/b);
  b = counters_before->GetSum("track:matched/rabs:yes/eta:yes/pdca:yes/chi2:yes/trig:low"); a = counters_after->GetSum("track:matched/rabs:yes/eta:yes/pdca:yes/chi2:yes/trig:low");
  printf("lowPt:   %9d %9d %7.2f %%\n", (Int_t)b, (Int_t)a, 100.*(a-b)/b);
  b = counters_before->GetSum("track:matched/rabs:yes/eta:yes/pdca:yes/chi2:yes/trig:high"); a = counters_after->GetSum("track:matched/rabs:yes/eta:yes/pdca:yes/chi2:yes/trig:high");
  printf("highPt:  %9d %9d %7.2f %%\n", (Int_t)b, (Int_t)a, 100.*(a-b)/b);
  printf("\n");
  */
}

//______________________________________________________________________________
void LoadAlirootLocally(TString& extraLibs)
{
  /// Load libraries locally
  
  // Load common libraries
  gSystem->Load("libVMC");
  gSystem->Load("libTree.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libMinuit.so");
  gSystem->Load("libXMLParser.so");
  gSystem->Load("libGui.so");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  
  // Load additional libraries
  gSystem->Load("libProofPlayer");
  TObjArray* libs = extraLibs.Tokenize(":");
  for (Int_t i = 0; i < libs->GetEntriesFast(); i++)
    gSystem->Load(Form("lib%s",static_cast<TObjString*>(libs->UncheckedAt(i))->GetName()));
  delete libs;
  
}

