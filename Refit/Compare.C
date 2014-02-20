/*
 *  Compare.C
 *  aliroot
 *
 *  Created by Philippe Pillot on 03/11/10.
 *  Copyright 2010 SUBATECH. All rights reserved.
 *
 */

TString trigType = "trig:any";

//------------------------------------------------------------------------------
void ComputeEfficiency(Double_t gen, Double_t rec, Double_t &eff, Double_t &err)
{
  /// compute the efficiency and the corresponding error
  if (gen > 0.) {
    eff = rec/gen;
    err = 100.*TMath::Max(1./gen, TMath::Sqrt(eff*TMath::Abs(1.-eff)/gen));
    eff *= 100.;
  } else {
    eff = 100.;
    err = 100.;
  }
}

//------------------------------------------------------------------------------
void FillEfficiency(TGraphErrors &g, Int_t i, Double_t x, Double_t gen, Double_t rec)
{
  /// compute the efficiency and the corresponding error and fill the graph
  Double_t eff, err;
  ComputeEfficiency(gen, rec, eff, err);
  g.SetPoint(i, x, eff);
  g.SetPointError(i, 0., err);
}

//------------------------------------------------------------------------------
void Compare(TString fileName_before = "AnalysisResults.root", TString fileName_after = "AnalysisResults.root",
	     Bool_t file1AsRef = kFALSE, Bool_t showBadTracks = kFALSE)
{
  /// compare results before/after refit if fileName_before==fileName_after
  /// or after refit in both files if they are different
  
  gROOT->SetStyle("Plain");
  gROOT->ForceStyle();
  gStyle->SetOptFit(1);
  
  // prepare environment
  gROOT->LoadMacro("$WORK/Macros/Facilities/runTaskFacilities.C");
  TString extraLibs = "RAWDatabase:CDB:STEER:MUONcore:MUONmapping:MUONcalib:MUONgeometry:MUONtrigger:MUONraw:MUONbase:MUONrec";
  LoadAlirootLocally(extraLibs,"","");
  
  Bool_t sameFile = (fileName_before==fileName_after);
  
  // open files
  TFile* file_before = TFile::Open(fileName_before.Data(),"READ");
  if (!file_before || !file_before->IsOpen()) return;
  TFile* file_after = file_before;
  if (!sameFile) {
    file_after = TFile::Open(fileName_after.Data(),"READ");
    if (!file_after || !file_after->IsOpen()) return;
  }
  
  TString afterSuf = showBadTracks ? "bad" : "after";
  
  // histo before
  TObjArray* histo_before;
  if (sameFile && !file1AsRef) histo_before = static_cast<TObjArray*>(file_before->FindObjectAny("Histograms_before"));
  else histo_before = static_cast<TObjArray*>(file_before->FindObjectAny(Form("Histograms_%s",afterSuf.Data())));
  if (!histo_before) return;
  
  // histo after
  TObjArray* histo_after = static_cast<TObjArray*>(file_after->FindObjectAny(Form("Histograms_%s",afterSuf.Data())));
  if (!histo_after) return;
  
  // counters before
  AliCounterCollection* counters_before;
  if (sameFile && !file1AsRef) counters_before = static_cast<AliCounterCollection*>(file_before->FindObjectAny("trackCounters_before"));
  else counters_before = static_cast<AliCounterCollection*>(file_before->FindObjectAny(Form("trackCounters_%s",afterSuf.Data())));
  if (!counters_before) return;
  
  // counters after
  AliCounterCollection* counters_after = static_cast<AliCounterCollection*>(file_after->FindObjectAny(Form("trackCounters_%s",afterSuf.Data())));
  if (!counters_after) return;
  
  // draw differences
  TString sHist[12] = {"hChi2", "hNClustersPerTrack", "hNChamberHitPerTrack", "hNTracks", "hDCA", "hPt", "hRapidity", "hNTrigAll", "hMass", "hPUncorrected", "hRAbs", "hNMatchTracks"};
  TCanvas* cHist = new TCanvas("cHist", "cHist", 1200, 800);
  cHist->Divide(4,3);
  TLegend *lHist = new TLegend(0.5,0.67,0.9,0.8);
  lHist->SetFillStyle(0);
  lHist->SetBorderSize(0);
  TCanvas* cDiff = new TCanvas("cDiff", "cDiff", 1200, 800);
  cDiff->Divide(4,3);
  TLegend *lDiff = new TLegend(0.5,0.7,0.9,0.8);
  lDiff->SetFillStyle(0);
  lDiff->SetBorderSize(0);
  TCanvas* cRatio = new TCanvas("cRatio", "cRatio", 1200, 800);
  cRatio->Divide(4,3);
  TLegend *lRatio = new TLegend(0.5,0.7,0.9,0.8);
  lRatio->SetFillStyle(0);
  lRatio->SetBorderSize(0);
  for (Int_t i=0; i<12; i++) {
    cHist->cd(i+1);
    gPad->SetLogy();
    TH1F* h_before = static_cast<TH1F*>(histo_before->FindObject(sHist[i].Data()));
    if (!h_before) continue;
    TH1F* h_after = static_cast<TH1F*>(histo_after->FindObject(sHist[i].Data())));
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
  
  // draw DCA
  TString sDCA[2] = {"hDCAX", "hDCAY"};
  TCanvas* cDCA = new TCanvas("cDCA", "cDCA", 800, 400);
  cDCA->Divide(2,1);
  for (Int_t i=0; i<2; i++) {
    cDCA->cd(i+1);
    gPad->SetLogy();
    TH1F* h_before = static_cast<TH1F*>(histo_before->FindObject(sDCA[i].Data()));
    if (!h_before) continue;
    TH1F* h_after = static_cast<TH1F*>(histo_after->FindObject(sDCA[i].Data())));
    h_before->Draw();
    h_before->GetXaxis()->SetRangeUser(-10.,10.);
    h_before->Fit("gaus","","same",-5.,5.);
    h_after->Draw("sames");
    h_after->SetLineColor(2);
    h_after->Fit("gaus","","same",-5.,5.);
    lHist->DrawClone("same");
  }
  
  // print statistic
  Double_t b,a,eff,err;
  printf("\n");
  printf("track       before     after        diff\n");
  b = counters_before->GetSum(); a = counters_after->GetSum();
  ComputeEfficiency(b, a, eff, err);
  printf("total:   %9d %9d %7.2f ± %4.2f %%\n", (Int_t)b, (Int_t)a, eff-100., err);
  b = counters_before->GetSum("track:matched"); a = counters_after->GetSum("track:matched");
  ComputeEfficiency(b, a, eff, err);
  printf("match:   %9d %9d %7.2f ± %4.2f %%\n", (Int_t)b, (Int_t)a, eff-100., err);
  b = counters_before->GetSum("track:matched/rabs:yes/eta:yes"); a = counters_after->GetSum("track:matched/rabs:yes/eta:yes");
  ComputeEfficiency(b, a, eff, err);
  printf("+ acc:   %9d %9d %7.2f ± %4.2f %%\n", (Int_t)b, (Int_t)a, eff-100., err);
  b = counters_before->GetSum("track:matched/rabs:yes/eta:yes/pdca:yes"); a = counters_after->GetSum("track:matched/rabs:yes/eta:yes/pdca:yes");
  ComputeEfficiency(b, a, eff, err);
  printf("+ DCA:   %9d %9d %7.2f ± %4.2f %%\n", (Int_t)b, (Int_t)a, eff-100., err);
  b = counters_before->GetSum("track:matched/rabs:yes/eta:yes/pdca:yes/pt:any-0.5GeV"); a = counters_after->GetSum("track:matched/rabs:yes/eta:yes/pdca:yes/pt:any-0.5GeV");
  ComputeEfficiency(b, a, eff, err);
  printf("<0.5GeV: %9d %9d %7.2f ± %4.2f %%\n", (Int_t)b, (Int_t)a, eff-100., err);
  b = counters_before->GetSum("track:matched/rabs:yes/eta:yes/pdca:yes/pt:0.5GeV"); a = counters_after->GetSum("track:matched/rabs:yes/eta:yes/pdca:yes/pt:0.5GeV");
  ComputeEfficiency(b, a, eff, err);
  printf(">0.5GeV: %9d %9d %7.2f ± %4.2f %%\n", (Int_t)b, (Int_t)a, eff-100., err);
  b = counters_before->GetSum("track:matched/rabs:yes/eta:yes/pdca:yes/pt:1GeV"); a = counters_after->GetSum("track:matched/rabs:yes/eta:yes/pdca:yes/pt:1GeV");
  ComputeEfficiency(b, a, eff, err);
  printf(">1GeV:   %9d %9d %7.2f ± %4.2f %%\n", (Int_t)b, (Int_t)a, eff-100., err);
  b = counters_before->GetSum("track:matched/rabs:yes/eta:yes/pdca:yes/pt:2GeV"); a = counters_after->GetSum("track:matched/rabs:yes/eta:yes/pdca:yes/pt:2GeV");
  ComputeEfficiency(b, a, eff, err);
  printf(">2GeV:   %9d %9d %7.2f ± %4.2f %%\n", (Int_t)b, (Int_t)a, eff-100., err);
  b = counters_before->GetSum("track:matched/rabs:yes/eta:yes/pdca:yes/pt:4GeV"); a = counters_after->GetSum("track:matched/rabs:yes/eta:yes/pdca:yes/pt:4GeV");
  ComputeEfficiency(b, a, eff, err);
  printf(">4GeV:   %9d %9d %7.2f ± %4.2f %%\n", (Int_t)b, (Int_t)a, eff-100., err);
  b = counters_before->GetSum("track:matched/rabs:yes/eta:yes/pdca:yes/trig:low"); a = counters_after->GetSum("track:matched/rabs:yes/eta:yes/pdca:yes/trig:low");
  ComputeEfficiency(b, a, eff, err);
  printf("lowPt:   %9d %9d %7.2f ± %4.2f %%\n", (Int_t)b, (Int_t)a, eff-100., err);
  b = counters_before->GetSum("track:matched/rabs:yes/eta:yes/pdca:yes/trig:low/pt:1GeV"); a = counters_after->GetSum("track:matched/rabs:yes/eta:yes/pdca:yes/trig:low/pt:1GeV");
  ComputeEfficiency(b, a, eff, err);
  printf(">1GeV:   %9d %9d %7.2f ± %4.2f %%\n", (Int_t)b, (Int_t)a, eff-100., err);
  b = counters_before->GetSum("track:matched/rabs:yes/eta:yes/pdca:yes/trig:high"); a = counters_after->GetSum("track:matched/rabs:yes/eta:yes/pdca:yes/trig:high");
  ComputeEfficiency(b, a, eff, err);
  printf("highPt:  %9d %9d %7.2f ± %4.2f %%\n", (Int_t)b, (Int_t)a, eff-100., err);
  b = counters_before->GetSum("track:matched/rabs:yes/eta:yes/pdca:yes/trig:high/pt:4GeV"); a = counters_after->GetSum("track:matched/rabs:yes/eta:yes/pdca:yes/trig:high/pt:4GeV");
  ComputeEfficiency(b, a, eff, err);
  printf(">4GeV:   %9d %9d %7.2f ± %4.2f %%\n", (Int_t)b, (Int_t)a, eff-100., err);
  printf("\n");
  
  // plot efficiency trends
  gStyle->SetFillColor(0);
  /*
  new TCanvas;
  TGraph* gOldAlignSelect = new TGraph(6);
  gOldAlignSelect->SetPoint(0,4.,99.7);
  gOldAlignSelect->SetPoint(1,2.,98.7);
  gOldAlignSelect->SetPoint(2,1.,97.2);
  gOldAlignSelect->SetPoint(3,0.65,96.);
  gOldAlignSelect->SetPoint(4,0.5,94.8);
  gOldAlignSelect->SetPoint(5,0.3,92.3);
  gOldAlignSelect->Draw("apc*");
  gOldAlignSelect->SetTitle("efficiency versus resolution: sigma cut = 4, alignment pass1/pass2");
  TGraph* gNewAlignSelect = new TGraph(6);
  gNewAlignSelect->SetPoint(0,4.,99.7);
  gNewAlignSelect->SetPoint(1,2.,98.7);
  gNewAlignSelect->SetPoint(2,1.,97.3);
  gNewAlignSelect->SetPoint(3,0.65,96.05);
  gNewAlignSelect->SetPoint(4,0.5,95.);
  gNewAlignSelect->SetPoint(5,0.3,92.7);
  gNewAlignSelect->Draw("samepc*");
  gNewAlignSelect->SetMarkerColor(2);
  gNewAlignSelect->SetLineColor(2);
  */
//  const Int_t nTestRes = 7;
//  Double_t chRes[nTestRes] = {0.25, 0.5, 0.5, 0.75, 1., 2., 4.};
//  TString sChRes[nTestRes] = {"0.75-0.25mm", "0.75-0.5mm", "1-0.5mm", "1-0.75mm", "1-1mm", "2-2mm", "4-4mm"};
  const Int_t nTestRes = 8;
  Double_t chRes[nTestRes] = {0.25, 0.5, 0.75, 1., 1.5, 2., 3., 4.};
  TString sChRes[nTestRes] = {"0.75-0.25mm", "0.75-0.5mm", "1-0.75mm", "1-1mm", "1.5-1.5mm", "2-2mm", "3-3mm", "4-4mm"};
//  const Int_t nTestRes = 6;
//  Double_t chRes[nTestRes] = {0.25, 0.5, 0.75, 1., 2., 4.};
//  TString sChRes[nTestRes] = {"0.75-0.25mm", "0.75-0.5mm", "1-0.75mm", "1-1mm", "2-2mm", "4-4mm"};
  Int_t nPadXRes = (nTestRes+1)/2;
//  const Int_t nTestSigTrk = 3;
//  Double_t sigmaCut[nTestSigTrk] = {2., 3., 4.};
//  TString sSigmaCut[nTestSigTrk] = {"2sigma", "3sigma", "4sigma"};
  const Int_t nTestSigTrk = 5;
  Double_t sigmaCut[nTestSigTrk] = {2., 3., 4., 5., 6.};
  TString sSigmaCut[nTestSigTrk] = {"2sigma", "3sigma", "4sigma", "5sigma", "6sigma"};
  Int_t nPadXSigTrk = (nTestSigTrk+1)/2;
  Double_t ranges[2];
  if (showBadTracks) {ranges[0] = 0.; ranges[1] = 25.;}
  else {ranges[0] = 75.; ranges[1] = 100.5;}
  
  // open all files and get the counters
  AliCounterCollection* countersRef[nTestRes][nTestSigTrk];
  AliCounterCollection* counters[nTestRes][nTestSigTrk];
  TH1F* hPUncorrectedRef[nTestRes][nTestSigTrk];
  TH1F* hPUncorrected[nTestRes][nTestSigTrk];
  for (Int_t i=0; i<nTestRes; i++) {
    for (Int_t j=0; j<nTestSigTrk; j++) {
      countersRef[i][j] = NULL;
      counters[i][j] = NULL;
      hPUncorrectedRef[i][j] = NULL;
      hPUncorrected[i][j] = NULL;
      TFile* file = TFile::Open(Form("%s_%s_4sigma/AnalysisResults.root",sChRes[i].Data(),sSigmaCut[j].Data()),"READ");
      if (!file || !file->IsOpen()) continue;
      if (file1AsRef) {
	countersRef[i][j] = counters_before;
	hPUncorrectedRef[i][j] = static_cast<TH1F*>(histo_before->FindObject("hPUncorrected"));
      } else {
	countersRef[i][j] = static_cast<AliCounterCollection*>(file->FindObjectAny("trackCounters_before"));
	hPUncorrectedRef[i][j] = static_cast<TH1F*>(static_cast<TObjArray*>(file->FindObjectAny("Histograms_before"))->FindObject("hPUncorrected"));
      }
      counters[i][j] = static_cast<AliCounterCollection*>(file->FindObjectAny(Form("trackCounters_%s",afterSuf.Data())));
      hPUncorrected[i][j] = static_cast<TH1F*>(static_cast<TObjArray*>(file->FindObjectAny(Form("Histograms_%s",afterSuf.Data())))->FindObject("hPUncorrected"));
    }
  }
  
  // Draw efficiency versus chamber resolution for different track selection and sigma cut
  TCanvas* cEffVsSigmaCut = new TCanvas("cEffVsSigmaCut", "cEffVsSigmaCut", 300.*nPadXRes, 600);
  cEffVsSigmaCut->Divide(nPadXRes,2);
  TLegend *lEff = new TLegend(0.6,0.4,1.,0.6);
  lEff->SetFillStyle(0);
  lEff->SetBorderSize(0);
  for (Int_t i=0; i<nTestRes; i++) {
    cEffVsSigmaCut->cd(i+1);
    TGraphErrors* g1 = new TGraphErrors(nTestSigTrk);
    TGraphErrors* g2 = new TGraphErrors(nTestSigTrk);
    TGraphErrors* g3 = new TGraphErrors(nTestSigTrk);
    TGraphErrors* g4 = new TGraphErrors(nTestSigTrk);
    for (Int_t j=0; j<nTestSigTrk; j++) {
      if (!counters[i][j]) continue;
      Double_t nAllRef = countersRef[i][j]->GetSum();
      Double_t nTrigRef = countersRef[i][j]->GetSum(Form("track:matched/%s",trigType.Data()));
      Double_t nTrigAccRef = countersRef[i][j]->GetSum(Form("track:matched/%s/rabs:yes/eta:yes",trigType.Data()));
      Double_t nTrigAccDCARef = countersRef[i][j]->GetSum(Form("track:matched/%s/rabs:yes/eta:yes/pdca:yes",trigType.Data()));
      Double_t nAll = counters[i][j]->GetSum();
      Double_t nTrig = counters[i][j]->GetSum(Form("track:matched/%s",trigType.Data()));
      Double_t nTrigAcc = counters[i][j]->GetSum(Form("track:matched/%s/rabs:yes/eta:yes",trigType.Data()));
      Double_t nTrigAccDCA = counters[i][j]->GetSum(Form("track:matched/%s/rabs:yes/eta:yes/pdca:yes",trigType.Data()));
      FillEfficiency(*g1, j, sigmaCut[j], nAllRef, nAll);
      FillEfficiency(*g2, j, sigmaCut[j], nTrigRef, nTrig);
      FillEfficiency(*g3, j, sigmaCut[j], nTrigAccRef, nTrigAcc);
      FillEfficiency(*g4, j, sigmaCut[j], nTrigAccDCARef, nTrigAccDCA);
    }
    g1->Draw("apl*");
    g1->SetTitle(Form("efficiency vs sigma cut: chamber resolution = %4.2f;#sigma cut;eff.", chRes[i]));
    g1->GetYaxis()->SetRangeUser(ranges[0], ranges[1]);
    g1->GetYaxis()->SetLabelSize(0.05);
    g1->GetXaxis()->SetLabelSize(0.05);
    g1->GetXaxis()->SetTitleSize(0.05);
    g1->GetXaxis()->SetTitleOffset(0.9);
    g1->GetYaxis()->SetLabelSize(0.05);
    g1->GetYaxis()->SetTitleSize(0.05);
    g2->Draw("samepl*");
    g2->SetMarkerColor(2);
    g2->SetLineColor(2);
    g3->Draw("samepl*");
    g3->SetMarkerColor(4);
    g3->SetLineColor(4);
    g4->Draw("samepl*");
    g4->SetMarkerColor(8);
    g4->SetLineColor(8);
  }
  cEffVsSigmaCut->cd(2*nPadXRes-1);
  lEff->AddEntry(g1,"all","l");
  lEff->AddEntry(g2,"+ match","l");
  lEff->AddEntry(g3,"+ acc","l");
  lEff->AddEntry(g4,"+ pDCA","l");
  lEff->Draw("same");
  
  // Draw efficiency versus sigma cut for different track selection and chamber resolution
  TCanvas* cEffVsChRes = new TCanvas("cEffVsChRes", "cEffVsChRes", 300.*nPadXSigTrk, 600);
  cEffVsChRes->Divide(nPadXSigTrk,2);
  for (Int_t j=0; j<nTestSigTrk; j++) {
    cEffVsChRes->cd(j+1);
    TGraphErrors* g1 = new TGraphErrors(nTestRes);
    TGraphErrors* g2 = new TGraphErrors(nTestRes);
    TGraphErrors* g3 = new TGraphErrors(nTestRes);
    TGraphErrors* g4 = new TGraphErrors(nTestRes);
    for (Int_t i=0; i<nTestRes; i++) {
      if (!counters[i][j]) continue;
      Double_t nAllRef = countersRef[i][j]->GetSum();
      Double_t nTrigRef = countersRef[i][j]->GetSum(Form("track:matched/%s",trigType.Data()));
      Double_t nTrigAccRef = countersRef[i][j]->GetSum(Form("track:matched/%s/rabs:yes/eta:yes",trigType.Data()));
      Double_t nTrigAccDCARef = countersRef[i][j]->GetSum(Form("track:matched/%s/rabs:yes/eta:yes/pdca:yes",trigType.Data()));
      Double_t nAll = counters[i][j]->GetSum();
      Double_t nTrig = counters[i][j]->GetSum(Form("track:matched/%s",trigType.Data()));
      Double_t nTrigAcc = counters[i][j]->GetSum(Form("track:matched/%s/rabs:yes/eta:yes",trigType.Data()));
      Double_t nTrigAccDCA = counters[i][j]->GetSum(Form("track:matched/%s/rabs:yes/eta:yes/pdca:yes",trigType.Data()));
      FillEfficiency(*g1, i, chRes[i], nAllRef, nAll);
      FillEfficiency(*g2, i, chRes[i], nTrigRef, nTrig);
      FillEfficiency(*g3, i, chRes[i], nTrigAccRef, nTrigAcc);
      FillEfficiency(*g4, i, chRes[i], nTrigAccDCARef, nTrigAccDCA);
    }
    g1->Draw("apl*");
    g1->SetTitle(Form("efficiency vs chamber resolution: sigma cut = %d;#sigma_{ch};eff.", (Int_t)sigmaCut[j]));
    g1->GetYaxis()->SetRangeUser(ranges[0], ranges[1]);
    g1->GetXaxis()->SetLabelSize(0.05);
    g1->GetXaxis()->SetTitleSize(0.05);
    g1->GetXaxis()->SetTitleOffset(0.9);
    g1->GetYaxis()->SetLabelSize(0.05);
    g1->GetYaxis()->SetTitleSize(0.05);
    g2->Draw("samepl*");
    g2->SetMarkerColor(2);
    g2->SetLineColor(2);
    g3->Draw("samepl*");
    g3->SetMarkerColor(4);
    g3->SetLineColor(4);
    g4->Draw("samepl*");
    g4->SetMarkerColor(8);
    g4->SetLineColor(8);
  }
  cEffVsChRes->cd(2*nPadXSigTrk-1);
  lEff->DrawClone("same");
  
  // Draw efficiency versus chamber resolution for different pt range and sigma cut
  TCanvas* cEffptVsSigmaCut = new TCanvas("cEffptVsSigmaCut", "cEffptVsSigmaCut", 300.*nPadXRes, 600);
  cEffptVsSigmaCut->Divide(nPadXRes,2);
  TLegend *lEffpt = new TLegend(0.45,0.4,0.9,0.6);
  lEffpt->SetFillStyle(0);
  lEffpt->SetBorderSize(0);
  for (Int_t i=0; i<nTestRes; i++) {
    cEffptVsSigmaCut->cd(i+1);
    TGraphErrors* g1 = new TGraphErrors(nTestSigTrk);
    TGraphErrors* g2 = new TGraphErrors(nTestSigTrk);
    TGraphErrors* g3 = new TGraphErrors(nTestSigTrk);
    TGraphErrors* g4 = new TGraphErrors(nTestSigTrk);
    TGraphErrors* g5 = new TGraphErrors(nTestSigTrk);
    for (Int_t j=0; j<nTestSigTrk; j++) {
      if (!counters[i][j]) continue;
      Double_t nPt0Ref = countersRef[i][j]->GetSum(Form("track:matched/%s/rabs:yes/eta:yes/pdca:yes",trigType.Data()));
      Double_t nPt05Ref = countersRef[i][j]->GetSum(Form("track:matched/%s/rabs:yes/eta:yes/pdca:yes/pt:0.5GeV",trigType.Data()));
      Double_t nPt1Ref = countersRef[i][j]->GetSum(Form("track:matched/%s/rabs:yes/eta:yes/pdca:yes/pt:1GeV",trigType.Data()));
      Double_t nPt2Ref = countersRef[i][j]->GetSum(Form("track:matched/%s/rabs:yes/eta:yes/pdca:yes/pt:2GeV",trigType.Data()));
      Double_t nPt4Ref = countersRef[i][j]->GetSum(Form("track:matched/%s/rabs:yes/eta:yes/pdca:yes/pt:4GeV",trigType.Data()));
      Double_t nPt0 = counters[i][j]->GetSum(Form("track:matched/%s/rabs:yes/eta:yes/pdca:yes",trigType.Data()));
      Double_t nPt05 = counters[i][j]->GetSum(Form("track:matched/%s/rabs:yes/eta:yes/pdca:yes/pt:0.5GeV",trigType.Data()));
      Double_t nPt1 = counters[i][j]->GetSum(Form("track:matched/%s/rabs:yes/eta:yes/pdca:yes/pt:1GeV",trigType.Data()));
      Double_t nPt2 = counters[i][j]->GetSum(Form("track:matched/%s/rabs:yes/eta:yes/pdca:yes/pt:2GeV",trigType.Data()));
      Double_t nPt4 = counters[i][j]->GetSum(Form("track:matched/%s/rabs:yes/eta:yes/pdca:yes/pt:4GeV",trigType.Data()));
      FillEfficiency(*g1, j, sigmaCut[j], nPt0Ref-nPt05Ref, nPt0-nPt05);
      FillEfficiency(*g2, j, sigmaCut[j], nPt05Ref-nPt1Ref, nPt05-nPt1);
      FillEfficiency(*g3, j, sigmaCut[j], nPt1Ref-nPt2Ref, nPt1-nPt2);
      FillEfficiency(*g4, j, sigmaCut[j], nPt2Ref-nPt4Ref, nPt2-nPt4);
      FillEfficiency(*g5, j, sigmaCut[j], nPt4Ref, nPt4);
    }
    g1->Draw("apl*");
    g1->SetTitle(Form("efficiency vs sigma cut: chamber resolution = %4.2f;#sigma cut;eff.", chRes[i]));
    g1->GetYaxis()->SetRangeUser(ranges[0], ranges[1]);
    g1->GetXaxis()->SetLabelSize(0.05);
    g1->GetXaxis()->SetTitleSize(0.05);
    g1->GetXaxis()->SetTitleOffset(0.9);
    g1->GetYaxis()->SetLabelSize(0.05);
    g1->GetYaxis()->SetTitleSize(0.05);
    g2->Draw("samepl*");
    g2->SetMarkerColor(2);
    g2->SetLineColor(2);
    g3->Draw("samepl*");
    g3->SetMarkerColor(6);
    g3->SetLineColor(6);
    g4->Draw("samepl*");
    g4->SetMarkerColor(4);
    g4->SetLineColor(4);
    g5->Draw("samepl*");
    g5->SetMarkerColor(8);
    g5->SetLineColor(8);
  }
  cEffptVsSigmaCut->cd(2*nPadXRes-1);
  lEffpt->AddEntry(g1,"p_{T} < 0.5 GeV/c","l");
  lEffpt->AddEntry(g2,"0.5 < p_{T} < 1 GeV/c","l");
  lEffpt->AddEntry(g3,"1 < p_{T} < 2 GeV/c","l");
  lEffpt->AddEntry(g4,"2 < p_{T} < 4 GeV/c","l");
  lEffpt->AddEntry(g5,"p_{T} > 4 GeV/c","l");
  lEffpt->Draw("same");
  
  // Draw efficiency versus sigma cut for different pt range and chamber resolution
  TCanvas* cEffptVsChRes = new TCanvas("cEffptVsChRes", "cEffptVsChRes", 300.*nPadXSigTrk, 600);
  cEffptVsChRes->Divide(nPadXSigTrk,2);
  for (Int_t j=0; j<nTestSigTrk; j++) {
    cEffptVsChRes->cd(j+1);
    TGraphErrors* g1 = new TGraphErrors(nTestRes);
    TGraphErrors* g2 = new TGraphErrors(nTestRes);
    TGraphErrors* g3 = new TGraphErrors(nTestRes);
    TGraphErrors* g4 = new TGraphErrors(nTestRes);
    TGraphErrors* g5 = new TGraphErrors(nTestRes);
    for (Int_t i=0; i<nTestRes; i++) {
      if (!counters[i][j]) continue;
      Double_t nPt0Ref = countersRef[i][j]->GetSum(Form("track:matched/%s/rabs:yes/eta:yes/pdca:yes",trigType.Data()));
      Double_t nPt05Ref = countersRef[i][j]->GetSum(Form("track:matched/%s/rabs:yes/eta:yes/pdca:yes/pt:0.5GeV",trigType.Data()));
      Double_t nPt1Ref = countersRef[i][j]->GetSum(Form("track:matched/%s/rabs:yes/eta:yes/pdca:yes/pt:1GeV",trigType.Data()));
      Double_t nPt2Ref = countersRef[i][j]->GetSum(Form("track:matched/%s/rabs:yes/eta:yes/pdca:yes/pt:2GeV",trigType.Data()));
      Double_t nPt4Ref = countersRef[i][j]->GetSum(Form("track:matched/%s/rabs:yes/eta:yes/pdca:yes/pt:4GeV",trigType.Data()));
      Double_t nPt0 = counters[i][j]->GetSum(Form("track:matched/%s/rabs:yes/eta:yes/pdca:yes",trigType.Data()));
      Double_t nPt05 = counters[i][j]->GetSum(Form("track:matched/%s/rabs:yes/eta:yes/pdca:yes/pt:0.5GeV",trigType.Data()));
      Double_t nPt1 = counters[i][j]->GetSum(Form("track:matched/%s/rabs:yes/eta:yes/pdca:yes/pt:1GeV",trigType.Data()));
      Double_t nPt2 = counters[i][j]->GetSum(Form("track:matched/%s/rabs:yes/eta:yes/pdca:yes/pt:2GeV",trigType.Data()));
      Double_t nPt4 = counters[i][j]->GetSum(Form("track:matched/%s/rabs:yes/eta:yes/pdca:yes/pt:4GeV",trigType.Data()));
      FillEfficiency(*g1, i, chRes[i], nPt0Ref-nPt05Ref, nPt0-nPt05);
      FillEfficiency(*g2, i, chRes[i], nPt05Ref-nPt1Ref, nPt05-nPt1);
      FillEfficiency(*g3, i, chRes[i], nPt1Ref-nPt2Ref, nPt1-nPt2);
      FillEfficiency(*g4, i, chRes[i], nPt2Ref-nPt4Ref, nPt2-nPt4);
      FillEfficiency(*g5, i, chRes[i], nPt4Ref, nPt4);
    }
    g1->Draw("apl*");
    g1->SetTitle(Form("efficiency vs chamber resolution: sigma cut = %d;#sigma_{ch};eff.", (Int_t)sigmaCut[j]));
    g1->GetYaxis()->SetRangeUser(ranges[0], ranges[1]);
    g1->GetXaxis()->SetLabelSize(0.05);
    g1->GetXaxis()->SetTitleSize(0.05);
    g1->GetXaxis()->SetTitleOffset(0.9);
    g1->GetYaxis()->SetLabelSize(0.05);
    g1->GetYaxis()->SetTitleSize(0.05);
    g2->Draw("samepl*");
    g2->SetMarkerColor(2);
    g2->SetLineColor(2);
    g3->Draw("samepl*");
    g3->SetMarkerColor(6);
    g3->SetLineColor(6);
    g4->Draw("samepl*");
    g4->SetMarkerColor(4);
    g4->SetLineColor(4);
    g5->Draw("samepl*");
    g5->SetMarkerColor(8);
    g5->SetLineColor(8);
  }
  cEffptVsChRes->cd(2*nPadXSigTrk-1);
  lEffpt->DrawClone("same");
  
  // Draw efficiency versus chamber resolution for different uncorrected-p range and sigma cut
  TCanvas* cEffpVsSigmaCut = new TCanvas("cEffpVsSigmaCut", "cEffpVsSigmaCut", 300.*nPadXRes, 600);
  cEffpVsSigmaCut->Divide(nPadXRes,2);
  TLegend *lEffp = new TLegend(0.45,0.4,0.9,0.6);
  lEffp->SetFillStyle(0);
  lEffp->SetBorderSize(0);
  for (Int_t i=0; i<nTestRes; i++) {
    cEffpVsSigmaCut->cd(i+1);
    TGraphErrors* g1 = new TGraphErrors(nTestSigTrk);
    TGraphErrors* g2 = new TGraphErrors(nTestSigTrk);
    TGraphErrors* g3 = new TGraphErrors(nTestSigTrk);
    TGraphErrors* g4 = new TGraphErrors(nTestSigTrk);
    TGraphErrors* g5 = new TGraphErrors(nTestSigTrk);
    for (Int_t j=0; j<nTestSigTrk; j++) {
      if (!hPUncorrected[i][j]) continue;
      Double_t nP0Ref = hPUncorrectedRef[i][j]->Integral(1,5);
      Double_t nP5Ref = hPUncorrectedRef[i][j]->Integral(6,10);
      Double_t nP10Ref = hPUncorrectedRef[i][j]->Integral(11,20);
      Double_t nP20Ref = hPUncorrectedRef[i][j]->Integral(21,40);
      Double_t nP40Ref = hPUncorrectedRef[i][j]->Integral(41,300);
      Double_t nP0 = hPUncorrected[i][j]->Integral(1,5);
      Double_t nP5 = hPUncorrected[i][j]->Integral(6,10);
      Double_t nP10 = hPUncorrected[i][j]->Integral(11,20);
      Double_t nP20 = hPUncorrected[i][j]->Integral(21,40);
      Double_t nP40 = hPUncorrected[i][j]->Integral(41,300);
      FillEfficiency(*g1, j, sigmaCut[j], nP0Ref, nP0);
      FillEfficiency(*g2, j, sigmaCut[j], nP5Ref, nP5);
      FillEfficiency(*g3, j, sigmaCut[j], nP10Ref, nP10);
      FillEfficiency(*g4, j, sigmaCut[j], nP20Ref, nP20);
      FillEfficiency(*g5, j, sigmaCut[j], nP40Ref, nP40);
    }
    g1->Draw("apl*");
    g1->SetTitle(Form("efficiency vs sigma cut: chamber resolution = %4.2f;#sigma cut;eff.", chRes[i]));
    g1->GetYaxis()->SetRangeUser(ranges[0], ranges[1]);
    g1->GetXaxis()->SetLabelSize(0.05);
    g1->GetXaxis()->SetTitleSize(0.05);
    g1->GetXaxis()->SetTitleOffset(0.9);
    g1->GetYaxis()->SetLabelSize(0.05);
    g1->GetYaxis()->SetTitleSize(0.05);
    g2->Draw("samepl*");
    g2->SetMarkerColor(2);
    g2->SetLineColor(2);
    g3->Draw("samepl*");
    g3->SetMarkerColor(6);
    g3->SetLineColor(6);
    g4->Draw("samepl*");
    g4->SetMarkerColor(4);
    g4->SetLineColor(4);
    g5->Draw("samepl*");
    g5->SetMarkerColor(8);
    g5->SetLineColor(8);
  }
  cEffpVsSigmaCut->cd(2*nPadXRes-1);
  lEffp->AddEntry(g1,"p < 5 GeV/c","l");
  lEffp->AddEntry(g2,"5 < p < 10 GeV/c","l");
  lEffp->AddEntry(g3,"10 < p < 20 GeV/c","l");
  lEffp->AddEntry(g4,"20 < p < 40 GeV/c","l");
  lEffp->AddEntry(g5,"p > 40 GeV/c","l");
  lEffp->Draw("same");
  
  // Draw efficiency versus sigma cut for different uncorrected-p range and chamber resolution
  TCanvas* cEffpVsChRes = new TCanvas("cEffpVsChRes", "cEffpVsChRes", 300.*nPadXSigTrk, 600);
  cEffpVsChRes->Divide(nPadXSigTrk,2);
  for (Int_t j=0; j<nTestSigTrk; j++) {
    cEffpVsChRes->cd(j+1);
    TGraphErrors* g1 = new TGraphErrors(nTestRes);
    TGraphErrors* g2 = new TGraphErrors(nTestRes);
    TGraphErrors* g3 = new TGraphErrors(nTestRes);
    TGraphErrors* g4 = new TGraphErrors(nTestRes);
    TGraphErrors* g5 = new TGraphErrors(nTestRes);
    for (Int_t i=0; i<nTestRes; i++) {
      if (!hPUncorrected[i][j]) continue;
      Double_t nP0Ref = hPUncorrectedRef[i][j]->Integral(1,5);
      Double_t nP5Ref = hPUncorrectedRef[i][j]->Integral(6,10);
      Double_t nP10Ref = hPUncorrectedRef[i][j]->Integral(11,20);
      Double_t nP20Ref = hPUncorrectedRef[i][j]->Integral(21,40);
      Double_t nP40Ref = hPUncorrectedRef[i][j]->Integral(41,300);
      Double_t nP0 = hPUncorrected[i][j]->Integral(1,5);
      Double_t nP5 = hPUncorrected[i][j]->Integral(6,10);
      Double_t nP10 = hPUncorrected[i][j]->Integral(11,20);
      Double_t nP20 = hPUncorrected[i][j]->Integral(21,40);
      Double_t nP40 = hPUncorrected[i][j]->Integral(41,300);
      FillEfficiency(*g1, i, chRes[i], nP0Ref, nP0);
      FillEfficiency(*g2, i, chRes[i], nP5Ref, nP5);
      FillEfficiency(*g3, i, chRes[i], nP10Ref, nP10);
      FillEfficiency(*g4, i, chRes[i], nP20Ref, nP20);
      FillEfficiency(*g5, i, chRes[i], nP40Ref, nP40);
    }
    g1->Draw("apl*");
    g1->SetTitle(Form("efficiency vs chamber resolution: sigma cut = %d;#sigma_{ch};eff.", (Int_t)sigmaCut[j]));
    g1->GetYaxis()->SetRangeUser(ranges[0], ranges[1]);
    g1->GetXaxis()->SetLabelSize(0.05);
    g1->GetXaxis()->SetTitleSize(0.05);
    g1->GetXaxis()->SetTitleOffset(0.9);
    g1->GetYaxis()->SetLabelSize(0.05);
    g1->GetYaxis()->SetTitleSize(0.05);
    g2->Draw("samepl*");
    g2->SetMarkerColor(2);
    g2->SetLineColor(2);
    g3->Draw("samepl*");
    g3->SetMarkerColor(6);
    g3->SetLineColor(6);
    g4->Draw("samepl*");
    g4->SetMarkerColor(4);
    g4->SetLineColor(4);
    g5->Draw("samepl*");
    g5->SetMarkerColor(8);
    g5->SetLineColor(8);
  }
  cEffpVsChRes->cd(2*nPadXSigTrk-1);
  lEffp->DrawClone("same");
  
  const Int_t nTestSigTrg = 5;
  Double_t sigmaCutTrig[nTestSigTrg] = {2., 3., 4., 5., 6.};
  TString sSigmaCutTrig[nTestSigTrg] = {"2sigma", "3sigma", "4sigma", "5sigma", "6sigma"};
//  const Int_t nTestSigTrg = 3;
//  Double_t sigmaCutTrig[nTestSigTrg] = {2., 3., 4.};
//  TString sSigmaCutTrig[nTestSigTrg] = {"2sigma", "3sigma", "4sigma"};
  Double_t rangesTrig[2];
  if (showBadTracks) {rangesTrig[0] = 0.; rangesTrig[1] = 25.;}
  else {rangesTrig[0] = 75.; rangesTrig[1] = 100.5;}
  
  // open all files and get the counters for trigger
  AliCounterCollection* countersRefTrig[nTestSigTrg];
  AliCounterCollection* countersTrig[nTestSigTrg];
  TH1F* hPUncorrectedRefTrig[nTestSigTrg];
  TH1F* hPUncorrectedTrig[nTestSigTrg];
  for (Int_t i=0; i<nTestSigTrg; i++) {
    countersRefTrig[i] = NULL;
    countersTrig[i] = NULL;
    hPUncorrectedRefTrig[i] = NULL;
    hPUncorrectedTrig[i] = NULL;
    TFile* file = TFile::Open(Form("2-2mm_4sigma_%s/AnalysisResults.root",sSigmaCutTrig[i].Data()),"READ");
//    TFile* file = TFile::Open(Form("1-1mm_4sigma_%s/AnalysisResults.root",sSigmaCutTrig[i].Data()),"READ");
//    TFile* file = TFile::Open(Form("2-2mm_4sigma_%s_Less5TrgTracks/AnalysisResults.root",sSigmaCutTrig[i].Data()),"READ");
//    TFile* file = TFile::Open(Form("2-2mm_4sigma_%s_More5TrgTracks/AnalysisResults.root",sSigmaCutTrig[i].Data()),"READ");
//    TFile* file = TFile::Open(Form("2-2mm_4sigma_%s_woAdjacentTrgTracks/AnalysisResults.root",sSigmaCutTrig[i].Data()),"READ");
//    TFile* file = TFile::Open(Form("2-2mm_4sigma_%s_wo2AdjacentTrgTracks/AnalysisResults.root",sSigmaCutTrig[i].Data()),"READ");
    if (!file || !file->IsOpen()) continue;
    if (file1AsRef) {
      countersRefTrig[i] = counters_before;
      hPUncorrectedRefTrig[i] = static_cast<TH1F*>(histo_before->FindObject("hPUncorrected"));
    } else {
      countersRefTrig[i] = static_cast<AliCounterCollection*>(file->FindObjectAny("trackCounters_before"));
      hPUncorrectedRefTrig[i] = static_cast<TH1F*>(static_cast<TObjArray*>(file->FindObjectAny("Histograms_before"))->FindObject("hPUncorrected"));
    }
    countersTrig[i] = static_cast<AliCounterCollection*>(file->FindObjectAny(Form("trackCounters_%s",afterSuf.Data())));
    hPUncorrectedTrig[i] = static_cast<TH1F*>(static_cast<TObjArray*>(file->FindObjectAny(Form("Histograms_%s",afterSuf.Data())))->FindObject("hPUncorrected"));
  }
  // Draw efficiency versus chamber resolution for different track selection and sigma cut
  TCanvas* cEffVsSigmaCutTrig = new TCanvas("cEffVsSigmaCutTrig", "cEffVsSigmaCutTrig", 600, 600);
  cEffVsSigmaCutTrig->Divide(2,2);
  
  cEffVsSigmaCutTrig->cd(1);
  TGraphErrors* g2 = new TGraphErrors(nTestSigTrg);
  TGraphErrors* g3 = new TGraphErrors(nTestSigTrg);
  TGraphErrors* g4 = new TGraphErrors(nTestSigTrg);
  TLegend *lEffTrig = new TLegend(0.6,0.42,1.,0.58);
  lEffTrig->SetFillStyle(0);
  lEffTrig->SetBorderSize(0);
  for (Int_t i=0; i<nTestSigTrg; i++) {
    if (!countersTrig[i]) continue;
    Double_t nTrigRefTrig = countersRefTrig[i]->GetSum(Form("track:matched/%s",trigType.Data()));
    Double_t nTrigAccRefTrig = countersRefTrig[i]->GetSum(Form("track:matched/%s/rabs:yes/eta:yes",trigType.Data()));
    Double_t nTrigAccDCARefTrig = countersRefTrig[i]->GetSum(Form("track:matched/%s/rabs:yes/eta:yes/pdca:yes",trigType.Data()));
    Double_t nTrig = countersTrig[i]->GetSum(Form("track:matched/%s",trigType.Data()));
    Double_t nTrigAcc = countersTrig[i]->GetSum(Form("track:matched/%s/rabs:yes/eta:yes",trigType.Data()));
    Double_t nTrigAccDCA = countersTrig[i]->GetSum(Form("track:matched/%s/rabs:yes/eta:yes/pdca:yes",trigType.Data()));
    FillEfficiency(*g2, i, sigmaCutTrig[i], nTrigRefTrig, nTrig);
    FillEfficiency(*g3, i, sigmaCutTrig[i], nTrigAccRefTrig, nTrigAcc);
    FillEfficiency(*g4, i, sigmaCutTrig[i], nTrigAccDCARefTrig, nTrigAccDCA);
  }
  g2->Draw("apl*");
  g2->SetTitle("efficiency vs trigger sigma cut for different track selection;trig #sigma cut;eff.");
  g2->GetYaxis()->SetRangeUser(rangesTrig[0], rangesTrig[1]);
  g2->GetXaxis()->SetLabelSize(0.05);
  g2->GetXaxis()->SetTitleSize(0.05);
  g2->GetXaxis()->SetTitleOffset(0.9);
  g2->GetYaxis()->SetLabelSize(0.05);
  g2->GetYaxis()->SetTitleSize(0.05);
  g2->SetMarkerColor(2);
  g2->SetLineColor(2);
  g3->Draw("samepl*");
  g3->SetMarkerColor(4);
  g3->SetLineColor(4);
  g4->Draw("samepl*");
  g4->SetMarkerColor(8);
  g4->SetLineColor(8);
  lEffTrig->AddEntry(g2,"match","l");
  lEffTrig->AddEntry(g3,"+ acc","l");
  lEffTrig->AddEntry(g4,"+ pDCA","l");
  lEffTrig->Draw("same");
  
  cEffVsSigmaCutTrig->cd(2);
  TGraphErrors* g1 = new TGraphErrors(nTestSigTrg);
  TGraphErrors* g2 = new TGraphErrors(nTestSigTrg);
  TGraphErrors* g3 = new TGraphErrors(nTestSigTrg);
  TGraphErrors* g4 = new TGraphErrors(nTestSigTrg);
  TGraphErrors* g5 = new TGraphErrors(nTestSigTrg);
  for (Int_t i=0; i<nTestSigTrg; i++) {
    if (!countersTrig[i]) continue;
    Double_t nPt0RefTrig = countersRefTrig[i]->GetSum(Form("track:matched/%s/rabs:yes/eta:yes/pdca:yes",trigType.Data()));
    Double_t nPt05RefTrig = countersRefTrig[i]->GetSum(Form("track:matched/%s/rabs:yes/eta:yes/pdca:yes/pt:0.5GeV",trigType.Data()));
    Double_t nPt1RefTrig = countersRefTrig[i]->GetSum(Form("track:matched/%s/rabs:yes/eta:yes/pdca:yes/pt:1GeV",trigType.Data()));
    Double_t nPt2RefTrig = countersRefTrig[i]->GetSum(Form("track:matched/%s/rabs:yes/eta:yes/pdca:yes/pt:2GeV",trigType.Data()));
    Double_t nPt4RefTrig = countersRefTrig[i]->GetSum(Form("track:matched/%s/rabs:yes/eta:yes/pdca:yes/pt:4GeV",trigType.Data()));
    Double_t nPt0 = countersTrig[i]->GetSum(Form("track:matched/%s/rabs:yes/eta:yes/pdca:yes",trigType.Data()));
    Double_t nPt05 = countersTrig[i]->GetSum(Form("track:matched/%s/rabs:yes/eta:yes/pdca:yes/pt:0.5GeV",trigType.Data()));
    Double_t nPt1 = countersTrig[i]->GetSum(Form("track:matched/%s/rabs:yes/eta:yes/pdca:yes/pt:1GeV",trigType.Data()));
    Double_t nPt2 = countersTrig[i]->GetSum(Form("track:matched/%s/rabs:yes/eta:yes/pdca:yes/pt:2GeV",trigType.Data()));
    Double_t nPt4 = countersTrig[i]->GetSum(Form("track:matched/%s/rabs:yes/eta:yes/pdca:yes/pt:4GeV",trigType.Data()));
    FillEfficiency(*g1, i, sigmaCutTrig[i], nPt0RefTrig-nPt05RefTrig, nPt0-nPt05);
    FillEfficiency(*g2, i, sigmaCutTrig[i], nPt05RefTrig-nPt1RefTrig, nPt05-nPt1);
    FillEfficiency(*g3, i, sigmaCutTrig[i], nPt1RefTrig-nPt2RefTrig, nPt1-nPt2);
    FillEfficiency(*g4, i, sigmaCutTrig[i], nPt2RefTrig-nPt4RefTrig, nPt2-nPt4);
    FillEfficiency(*g5, i, sigmaCutTrig[i], nPt4RefTrig, nPt4);
  }
  g1->Draw("apl*");
  g1->SetTitle("efficiency vs trigger sigma cut for different pt;trig #sigma cut;eff.");
  g1->GetYaxis()->SetRangeUser(rangesTrig[0], rangesTrig[1]);
  g1->GetXaxis()->SetLabelSize(0.05);
  g1->GetXaxis()->SetTitleSize(0.05);
  g1->GetXaxis()->SetTitleOffset(0.9);
  g1->GetYaxis()->SetLabelSize(0.05);
  g1->GetYaxis()->SetTitleSize(0.05);
  g2->Draw("samepl*");
  g2->SetMarkerColor(2);
  g2->SetLineColor(2);
  g3->Draw("samepl*");
  g3->SetMarkerColor(6);
  g3->SetLineColor(6);
  g4->Draw("samepl*");
  g4->SetMarkerColor(4);
  g4->SetLineColor(4);
  g5->Draw("samepl*");
  g5->SetMarkerColor(8);
  g5->SetLineColor(8);
  lEffpt->DrawClone("same");
  
  cEffVsSigmaCutTrig->cd(3);
  TGraphErrors* g1 = new TGraphErrors(nTestSigTrg);
  TGraphErrors* g2 = new TGraphErrors(nTestSigTrg);
  TGraphErrors* g3 = new TGraphErrors(nTestSigTrg);
  TGraphErrors* g4 = new TGraphErrors(nTestSigTrg);
  TGraphErrors* g5 = new TGraphErrors(nTestSigTrg);
  for (Int_t i=0; i<nTestSigTrg; i++) {
    if (!hPUncorrectedTrig[i]) continue;
    Double_t nP0RefTrig = hPUncorrectedRefTrig[i]->Integral(1,5);
    Double_t nP5RefTrig = hPUncorrectedRefTrig[i]->Integral(6,10);
    Double_t nP10RefTrig = hPUncorrectedRefTrig[i]->Integral(11,20);
    Double_t nP20RefTrig = hPUncorrectedRefTrig[i]->Integral(21,40);
    Double_t nP40RefTrig = hPUncorrectedRefTrig[i]->Integral(41,300);
    Double_t nP0 = hPUncorrectedTrig[i]->Integral(1,5);
    Double_t nP5 = hPUncorrectedTrig[i]->Integral(6,10);
    Double_t nP10 = hPUncorrectedTrig[i]->Integral(11,20);
    Double_t nP20 = hPUncorrectedTrig[i]->Integral(21,40);
    Double_t nP40 = hPUncorrectedTrig[i]->Integral(41,300);
    FillEfficiency(*g1, i, sigmaCutTrig[i], nP0RefTrig, nP0);
    FillEfficiency(*g2, i, sigmaCutTrig[i], nP5RefTrig, nP5);
    FillEfficiency(*g3, i, sigmaCutTrig[i], nP10RefTrig, nP10);
    FillEfficiency(*g4, i, sigmaCutTrig[i], nP20RefTrig, nP20);
    FillEfficiency(*g5, i, sigmaCutTrig[i], nP40RefTrig, nP40);
  }
  g1->Draw("apl*");
  g1->SetTitle("efficiency vs trigger sigma cut for different p;trig #sigma cut;eff.");
  g1->GetYaxis()->SetRangeUser(rangesTrig[0], rangesTrig[1]);
  g1->GetXaxis()->SetLabelSize(0.05);
  g1->GetXaxis()->SetTitleSize(0.05);
  g1->GetXaxis()->SetTitleOffset(0.9);
  g1->GetYaxis()->SetLabelSize(0.05);
  g1->GetYaxis()->SetTitleSize(0.05);
  g2->Draw("samepl*");
  g2->SetMarkerColor(2);
  g2->SetLineColor(2);
  g3->Draw("samepl*");
  g3->SetMarkerColor(6);
  g3->SetLineColor(6);
  g4->Draw("samepl*");
  g4->SetMarkerColor(4);
  g4->SetLineColor(4);
  g5->Draw("samepl*");
  g5->SetMarkerColor(8);
  g5->SetLineColor(8);
  lEffp->DrawClone("same");
  
  cEffVsSigmaCutTrig->cd(4);
  TGraphErrors* g2 = new TGraphErrors(nTestSigTrg);
  TGraphErrors* g3 = new TGraphErrors(nTestSigTrg);
  TGraphErrors* g4 = new TGraphErrors(nTestSigTrg);
  TLegend *lEffTrigTrig = new TLegend(0.6,0.42,1.,0.58);
  lEffTrigTrig->SetFillStyle(0);
  lEffTrigTrig->SetBorderSize(0);
  for (Int_t i=0; i<nTestSigTrg; i++) {
    if (!countersTrig[i]) continue;
    Double_t nTrigRefTrig = countersRefTrig[i]->GetSum("track:matched/rabs:yes/eta:yes/pdca:yes");
    Double_t nTrigLowRefTrig = countersRefTrig[i]->GetSum("track:matched/rabs:yes/eta:yes/pdca:yes/trig:low");
    Double_t nTrigHighRefTrig = countersRefTrig[i]->GetSum("track:matched/rabs:yes/eta:yes/pdca:yes/trig:high");
    Double_t nTrig = countersTrig[i]->GetSum("track:matched/rabs:yes/eta:yes/pdca:yes");
    Double_t nTrigLow = countersTrig[i]->GetSum("track:matched/rabs:yes/eta:yes/pdca:yes/trig:low");
    Double_t nTrigHigh = countersTrig[i]->GetSum("track:matched/rabs:yes/eta:yes/pdca:yes/trig:high");
    FillEfficiency(*g2, i, sigmaCutTrig[i], nTrigRefTrig-nTrigLowRefTrig, nTrig-nTrigLow);
    FillEfficiency(*g3, i, sigmaCutTrig[i], nTrigLowRefTrig-nTrigHighRefTrig, nTrigLow-nTrigHigh);
    FillEfficiency(*g4, i, sigmaCutTrig[i], nTrigHighRefTrig, nTrigHigh);
  }
  g2->Draw("apl*");
  g2->SetTitle("efficiency vs trigger sigma cut for different trigger level;trig #sigma cut;eff.");
  g2->GetYaxis()->SetRangeUser(rangesTrig[0], rangesTrig[1]);
  g2->GetXaxis()->SetLabelSize(0.05);
  g2->GetXaxis()->SetTitleSize(0.05);
  g2->GetXaxis()->SetTitleOffset(0.9);
  g2->GetYaxis()->SetLabelSize(0.05);
  g2->GetYaxis()->SetTitleSize(0.05);
  g2->SetMarkerColor(2);
  g2->SetLineColor(2);
  g3->Draw("samepl*");
  g3->SetMarkerColor(4);
  g3->SetLineColor(4);
  g4->Draw("samepl*");
  g4->SetMarkerColor(8);
  g4->SetLineColor(8);
  lEffTrigTrig->AddEntry(g2,"All p_{T}","l");
  lEffTrigTrig->AddEntry(g3,"Low p_{T}","l");
  lEffTrigTrig->AddEntry(g4,"High p_{T}","l");
  lEffTrigTrig->Draw("same");
  
}

