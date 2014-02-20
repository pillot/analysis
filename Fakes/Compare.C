/*
 *  Compare.C
 *  aliroot
 *
 *  Created by Philippe Pillot on 08/11/10.
 *  Copyright 2010 SUBATECH. All rights reserved.
 *
 */

void Compare(TString fileName_before = "AnalysisResults.root", TString fileName_after = "AnalysisResults.root")
{
  /// compare results before/after refit if fileName_before==fileName_after
  /// or after refit in both files if they are different
  
  // prepare environment
  TString extraLibs = "RAWDatabase:CDB:STEER:MUONcore:MUONmapping:MUONcalib:MUONgeometry:MUONtrigger:MUONraw:MUONbase:MUONrec:PWG3base";
  LoadAlirootLocally(extraLibs);
  
  Bool_t sameFile = (fileName_before==fileName_after);
  
  // open files
  TFile* file_before = TFile::Open(fileName_before.Data(),"READ");
  if (!file_before || !file_before->IsOpen()) return;
  TFile* file_after = file_before;
  if (!sameFile) {
    file_after = TFile::Open(fileName_after.Data(),"READ");
    if (!file_after || !file_after->IsOpen()) return;
  }
  
  // histo before
  TObjArray* histo_before;
  if (sameFile) histo_before = static_cast<TObjArray*>(file_before->FindObjectAny("Histograms_before"));
  else histo_before = static_cast<TObjArray*>(file_before->FindObjectAny("Histograms_after"));
  if (!histo_before) return;
  
  // histo after
  TObjArray* histo_after = static_cast<TObjArray*>(file_after->FindObjectAny("Histograms_after"));
  if (!histo_after) return;
  
  // counters before
  AliCounterCollection* counters_before;
  if (sameFile) counters_before = static_cast<AliCounterCollection*>(file_before->FindObjectAny("trackCounters_before"));
  else counters_before = static_cast<AliCounterCollection*>(file_before->FindObjectAny("trackCounters_after"));
  if (!counters_before) return;
  
  // counters after
  AliCounterCollection* counters_after = static_cast<AliCounterCollection*>(file_after->FindObjectAny("trackCounters_after"));
  if (!counters_after) return;
  
  // histo fakes before
  TObjArray* histoFakes_before = static_cast<TObjArray*>(file_before->FindObjectAny("histos"));
  if (!histoFakes_before) return;
  
  // histo fakes after
  TObjArray* histoFakes_after = static_cast<TObjArray*>(file_after->FindObjectAny("histos"));
  if (!histoFakes_after) return;
  
  // counters fakes before
  AliCounterCollection* countersFakes_before = static_cast<AliCounterCollection*>(file_before->FindObjectAny("trackCounters"));
  if (!countersFakes_before) return;
  
  // counters fakes after
  AliCounterCollection* countersFakes_after = static_cast<AliCounterCollection*>(file_after->FindObjectAny("trackCounters"));
  if (!countersFakes_after) return;
  
  // draw differences in selected tracks
  TString sHist[8] = {"hChi2", "hNClustersPerTrack", "hNChamberHitPerTrack", "hDCA", "hPt", "hRapidity", "hMass", "hPUncorrected"};
  TCanvas* cHist = new TCanvas("cHist", "cHist", 1000, 900);
  cHist->Divide(3,3);
  TCanvas* cDiff = new TCanvas("cDiff", "cDiff", 1000, 900);
  cDiff->Divide(3,3);
  TCanvas* cRatio = new TCanvas("cRatio", "cRatio", 1000, 900);
  cRatio->Divide(3,3);
  for (Int_t i=0; i<8; i++) {
    cHist->cd(i+1);
    gPad->SetLogy();
    TH1F* h_before = static_cast<TH1F*>(histo_before->FindObject(sHist[i].Data()));
    if (!h_before) continue;
    TH1F* h_after = static_cast<TH1F*>(histo_after->FindObject(sHist[i].Data())));
    h_before->Draw();
    h_after->Draw("sames");
    h_after->SetLineColor(2);
    cDiff->cd(i+1);
    TH1F* h_diff = h_before->Clone();
    h_diff->Add(h_after, -1.);
    h_diff->Draw();
    h_diff->SetLineColor(4);
    cRatio->cd(i+1);
    TH1F* h_ratio = h_after->Clone();
    h_ratio->Divide(h_before);
    h_ratio->Draw();
    h_ratio->SetLineColor(4);
  }
  
  // draw differences in fake tracks
  TString sFakes[8] = {"hNumberOfClustersF", "hChi2PerDofF", "hPF", "hEtaF", "hNumberOfChamberHitF", "hDCAF", "hPtF", "hPhiF"};
  TCanvas* cFakes = new TCanvas("cFakes", "cFakes", 1200, 600);
  cFakes->Divide(4,2);
  for (Int_t i=0; i<8; i++) {
    cFakes->cd(i+1);
    gPad->SetLogy();
    TH1F* h_before = static_cast<TH1F*>(histoFakes_before->FindObject(sFakes[i].Data()));
    if (!h_before) continue;
    TH1F* h_after = static_cast<TH1F*>(histoFakes_after->FindObject(sFakes[i].Data())));
    h_before->Draw();
    h_before->SetLineColor(1);
    h_before->SetFillStyle(0);
    h_after->Draw("sames");
    h_after->SetLineColor(2);
    h_after->SetFillStyle(0);
  }
  
  // print statistic of all tracks
  Double_t b,a;
  printf("\n");
  printf("statistic of all tracks:\n");
  printf("track       before     after     diff\n");
  b = counters_before->GetSum(); a = counters_after->GetSum();
  printf("total:   %9d %9d %7.2f %%\n", (Int_t)b, (Int_t)a, 100.*(a-b)/b);
  b = counters_before->GetSum("track:matched"); a = counters_after->GetSum("track:matched");
  printf("match:   %9d %9d %7.2f %%\n", (Int_t)b, (Int_t)a, 100.*(a-b)/b);
  b = counters_before->GetSum("track:matched/rabs:yes/eta:yes"); a = counters_after->GetSum("track:matched/rabs:yes/eta:yes");
  printf("selec:   %9d %9d %7.2f %%\n", (Int_t)b, (Int_t)a, 100.*(a-b)/b);
  b = counters_before->GetSum("track:matched/rabs:yes/eta:yes/pt:0.5GeV"); a = counters_after->GetSum("track:matched/rabs:yes/eta:yes/pt:0.5GeV");
  printf(">0.5GeV: %9d %9d %7.2f %%\n", (Int_t)b, (Int_t)a, 100.*(a-b)/b);
  b = counters_before->GetSum("track:matched/rabs:yes/eta:yes/pt:1GeV"); a = counters_after->GetSum("track:matched/rabs:yes/eta:yes/pt:1GeV");
  printf(">1GeV:   %9d %9d %7.2f %%\n", (Int_t)b, (Int_t)a, 100.*(a-b)/b);
  b = counters_before->GetSum("track:matched/rabs:yes/eta:yes/pt:2GeV"); a = counters_after->GetSum("track:matched/rabs:yes/eta:yes/pt:2GeV");
  printf(">2GeV:   %9d %9d %7.2f %%\n", (Int_t)b, (Int_t)a, 100.*(a-b)/b);
  printf("\n");
  
  // print statistic of fakes
  printf("\n");
  printf("statistic of good tracks:\n");
  printf("track       before     after     diff\n");
  b = countersFakes_before->GetSum("track:matched"); a = countersFakes_after->GetSum("track:matched");
  printf("total:   %9d %9d %7.2f %%\n", (Int_t)b, (Int_t)a, 100.*(a-b)/b);
  b = countersFakes_before->GetSum("track:matched/trig:yes"); a = countersFakes_after->GetSum("track:matched/trig:yes");
  printf("match:   %9d %9d %7.2f %%\n", (Int_t)b, (Int_t)a, 100.*(a-b)/b);
  b = countersFakes_before->GetSum("track:matched/trig:yes/acc:in"); a = countersFakes_after->GetSum("track:matched/trig:yes/acc:in");
  printf("selec:   %9d %9d %7.2f %%\n", (Int_t)b, (Int_t)a, 100.*(a-b)/b);
  printf("\n");
  
  // print statistic of fakes
  printf("\n");
  printf("statistic of fake tracks:\n");
  printf("track       before     after     diff\n");
  b = countersFakes_before->GetSum("track:fake"); a = countersFakes_after->GetSum("track:fake");
  printf("total:   %9d %9d %7.2f %%\n", (Int_t)b, (Int_t)a, 100.*(a-b)/b);
  b = countersFakes_before->GetSum("track:fake/trig:yes"); a = countersFakes_after->GetSum("track:fake/trig:yes");
  printf("match:   %9d %9d %7.2f %%\n", (Int_t)b, (Int_t)a, 100.*(a-b)/b);
  b = countersFakes_before->GetSum("track:fake/trig:yes/acc:in"); a = countersFakes_after->GetSum("track:fake/trig:yes/acc:in");
  printf("selec:   %9d %9d %7.2f %%\n", (Int_t)b, (Int_t)a, 100.*(a-b)/b);
  printf("\n");
  
  // Draw efficiency trends
  gStyle->SetFillColor(0);
  
  Double_t chRes[6] = {0.3, 0.5, 0.65, 1., 2., 4.};
  TString sChRes[6] = {"0.5-0.3mm", "0.5mm", "0.65mm", "1mm", "2mm", "4mm"};
  Double_t sigmaCut[5] = {2., 3., 4., 5., 6.};
  TString sSigmaCut[5] = {"_2sigma", "_3sigma", "_4sigma", "_5sigma", "_6sigma"};
  Double_t ranges[2] = {0., 100.5};
  
  // open all files and get the counters
  AliCounterCollection* counters[6][5];
  for (Int_t i=0; i<6; i++) {
    for (Int_t j=0; j<5; j++) {
      counters[i][j] = NULL;
      TFile* file = TFile::Open(Form("%s%s/AnalysisResults.root",sChRes[i].Data(),sSigmaCut[j].Data()),"READ");
      if (!file || !file->IsOpen()) continue;
      counters[i][j] = static_cast<AliCounterCollection*>(file->FindObjectAny("trackCounters"));
    }
  }
  
  Int_t iRef = 5, jRef = 4;
  Double_t nAllMRef, nTrigMRef, nTrigAccMRef;
  Double_t nAllFRef, nTrigFRef, nTrigAccFRef;
  if (counters[iRef][jRef]) {
    
    nAllMRef = counters[iRef][jRef]->GetSum("track:matched");
    nTrigMRef = counters[iRef][jRef]->GetSum("track:matched/trig:yes");
    nTrigAccMRef = counters[iRef][jRef]->GetSum("track:matched/trig:yes/acc:in");
    nAllFRef = counters[iRef][jRef]->GetSum("track:fake");
    nTrigFRef = counters[iRef][jRef]->GetSum("track:fake/trig:yes");
    nTrigAccFRef = counters[iRef][jRef]->GetSum("track:fake/trig:yes/acc:in");
    
    // Draw efficiency versus chamber resolution for different track selection and sigma cut
    TCanvas* cEffVsSigmaCut = new TCanvas("cEffVsSigmaCut", "cEffVsSigmaCut", 1000, 600);
    cEffVsSigmaCut->Divide(3,2);
    for (Int_t i=0; i<6; i++) {
      cEffVsSigmaCut->cd(i+1);
      TGraph* g1 = new TGraph(5);
      TGraph* g2 = new TGraph(5);
      TGraph* g3 = new TGraph(5);
      TGraph* g4 = new TGraph(5);
      TGraph* g5 = new TGraph(5);
      TGraph* g6 = new TGraph(5);
      for (Int_t j=0; j<5; j++) {
	if (!counters[i][j]) continue;
	Double_t nAllM = counters[i][j]->GetSum("track:matched");
	Double_t nTrigM = counters[i][j]->GetSum("track:matched/trig:yes");
	Double_t nTrigAccM = counters[i][j]->GetSum("track:matched/trig:yes/acc:in");
	Double_t nAllF = counters[i][j]->GetSum("track:fake");
	Double_t nTrigF = counters[i][j]->GetSum("track:fake/trig:yes");
	Double_t nTrigAccF = counters[i][j]->GetSum("track:fake/trig:yes/acc:in");
	g1->SetPoint(j, sigmaCut[j], 100. + 100.*(nAllM-nAllMRef)/nAllMRef);
	g2->SetPoint(j, sigmaCut[j], 100. + 100.*(nTrigM-nTrigMRef)/nTrigMRef);
	g3->SetPoint(j, sigmaCut[j], 100. + 100.*(nTrigAccM-nTrigAccMRef)/nTrigAccMRef);
	g4->SetPoint(j, sigmaCut[j], - 100.*(nAllF-nAllFRef)/nAllFRef);
	g5->SetPoint(j, sigmaCut[j], - 100.*(nTrigF-nTrigFRef)/nTrigFRef);
	g6->SetPoint(j, sigmaCut[j], - 100.*(nTrigAccF-nTrigAccFRef)/nTrigAccFRef);
      }
      g1->Draw("apc*");
      g1->SetTitle(Form("efficiency vs sigma cut: chamber resolution = %4.2f;#sigma cut;eff.", chRes[i]));
      g1->GetYaxis()->SetRangeUser(ranges[0], ranges[1]);
      g1->GetYaxis()->SetLabelSize(0.05);
      g1->GetXaxis()->SetLabelSize(0.05);
      g1->GetXaxis()->SetTitleSize(0.05);
      g1->GetXaxis()->SetTitleOffset(0.9);
      g1->GetYaxis()->SetLabelSize(0.05);
      g1->GetYaxis()->SetTitleSize(0.05);
//      g2->Draw("samepc*");
      g2->SetMarkerColor(2);
      g2->SetLineColor(2);
//      g3->Draw("samepc*");
      g3->SetMarkerColor(4);
      g3->SetLineColor(4);
      g4->Draw("samepc*");
      g4->SetLineStyle(2);
//      g5->Draw("samepc*");
      g5->SetMarkerColor(2);
      g5->SetLineColor(2);
      g5->SetLineStyle(2);
//      g6->Draw("samepc*");
      g6->SetMarkerColor(4);
      g6->SetLineColor(4);
      g6->SetLineStyle(2);
    }
    
    // Draw efficiency versus sigma cut for different track selection and chamber resolution
    TCanvas* cEffVsChRes = new TCanvas("cEffVsChRes", "cEffVsChRes", 1000, 600);
    cEffVsChRes->Divide(3,2);
    for (Int_t j=0; j<5; j++) {
      cEffVsChRes->cd(j+1);
      TGraph* g1 = new TGraph(6);
      TGraph* g2 = new TGraph(6);
      TGraph* g3 = new TGraph(6);
      TGraph* g4 = new TGraph(6);
      TGraph* g5 = new TGraph(6);
      TGraph* g6 = new TGraph(6);
      for (Int_t i=0; i<6; i++) {
	if (!counters[i][j]) continue;
	Double_t nAllM = counters[i][j]->GetSum("track:matched");
	Double_t nTrigM = counters[i][j]->GetSum("track:matched/trig:yes");
	Double_t nTrigAccM = counters[i][j]->GetSum("track:matched/trig:yes/acc:in");
	Double_t nAllF = counters[i][j]->GetSum("track:fake");
	Double_t nTrigF = counters[i][j]->GetSum("track:fake/trig:yes");
	Double_t nTrigAccF = counters[i][j]->GetSum("track:fake/trig:yes/acc:in");
	g1->SetPoint(i, chRes[i], 100. + 100.*(nAllM-nAllMRef)/nAllMRef);
	g2->SetPoint(i, chRes[i], 100. + 100.*(nTrigM-nTrigMRef)/nTrigMRef);
	g3->SetPoint(i, chRes[i], 100. + 100.*(nTrigAccM-nTrigAccMRef)/nTrigAccMRef);
	g4->SetPoint(i, chRes[i], - 100.*(nAllF-nAllFRef)/nAllFRef);
	g5->SetPoint(i, chRes[i], - 100.*(nTrigF-nTrigFRef)/nTrigFRef);
	g6->SetPoint(i, chRes[i], - 100.*(nTrigAccF-nTrigAccFRef)/nTrigAccFRef);
      }
      g1->Draw("apc*");
      g1->SetTitle(Form("efficiency vs chamber resolution: sigma cut = %d;#sigma_{ch};eff.", (Int_t)sigmaCut[j]));
      g1->GetYaxis()->SetRangeUser(ranges[0], ranges[1]);
      g1->GetXaxis()->SetLabelSize(0.05);
      g1->GetXaxis()->SetTitleSize(0.05);
      g1->GetXaxis()->SetTitleOffset(0.9);
      g1->GetYaxis()->SetLabelSize(0.05);
      g1->GetYaxis()->SetTitleSize(0.05);
//      g2->Draw("samepc*");
      g2->SetMarkerColor(2);
      g2->SetLineColor(2);
//      g3->Draw("samepc*");
      g3->SetMarkerColor(4);
      g3->SetLineColor(4);
      g4->Draw("samepc*");
      g4->SetLineStyle(2);
//      g5->Draw("samepc*");
      g5->SetMarkerColor(2);
      g5->SetLineColor(2);
      g5->SetLineStyle(2);
//      g6->Draw("samepc*");
      g6->SetMarkerColor(4);
      g6->SetLineColor(4);
      g6->SetLineStyle(2);
    }
    
  }
  
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

