/*
 *  plotMuonFakes.C
 *  aliroot
 *
 *  Created by philippe pillot on 10/02/11.
 *  Copyright 2011 SUBATECH. All rights reserved.
 *
 */


void plotMuonFakes(TString fileName = "AnalysisResults.root", Int_t run = -1)
{
  /// display statistic of fakes
  
  TString srun = (run < 0) ? "run:any" : Form("run:%d", run);
  
  // prepare environment
  gROOT->LoadMacro("/Users/pillot/Work/Alice/Work/Macros/Facilities/runTaskFacilities.C");
  LoadAlirootLocally("", "", "");
  
  // open file
  TFile* file = TFile::Open(fileName.Data(),"READ");
  if (!file || !file->IsOpen()) return;
  
  // get counters
  AliCounterCollection* trackCounters = static_cast<AliCounterCollection*>(file->FindObjectAny("trackCounters_CombiOCDB"));
  if (!trackCounters) return;
  AliCounterCollection* pairCounters = static_cast<AliCounterCollection*>(file->FindObjectAny("pairCounters_CombiOCDB"));
  if (!pairCounters) return;
  
  // get statistics of singles versus run
  TH1D* hTrackablesVsRun = trackCounters->Get("run",Form("track:reconstructible/%s",srun.Data()));
  TH1D* hMatchedVsRun = trackCounters->Get("run",Form("track:matched/%s",srun.Data()));
  TGraphAsymmErrors* gEffVsRun = new TGraphAsymmErrors(hMatchedVsRun,hTrackablesVsRun);
  TH1D* hTracksVsRun = trackCounters->Get("run",Form("track:reconstructed/%s",srun.Data()));
  TH1D* hFakesVsRun = trackCounters->Get("run",Form("track:fake/%s",srun.Data()));
  TGraphAsymmErrors* gRatioFakesVsRun = new TGraphAsymmErrors(hFakesVsRun,hTracksVsRun);
  TH1D* hDecaysVsRun = trackCounters->Get("run",Form("track:decay/%s",srun.Data()));
  TGraphAsymmErrors* gRatioDecaysVsRun = new TGraphAsymmErrors(hDecaysVsRun,hTracksVsRun);
  
  // get statistics of singles versus centrality
  TH1D* hTrackables = trackCounters->Get("cent",Form("track:reconstructible/%s",srun.Data()));
  TH1D* hMatched = trackCounters->Get("cent",Form("track:matched/%s",srun.Data()));
  TGraphAsymmErrors* gEff = new TGraphAsymmErrors(hMatched,hTrackables);
  TH1D* hMatchTrackable = hMatched->Clone();
  hMatchTrackable->Add(trackCounters->Get("cent",Form("track:matchedyet/%s",srun.Data())),-1.);
  TGraphAsymmErrors* gEffTrackable = new TGraphAsymmErrors(hMatchTrackable,hTrackables);
  TH1D* hTracks = trackCounters->Get("cent",Form("track:reconstructed/%s",srun.Data()));
  TH1D* hFakes = trackCounters->Get("cent",Form("track:fake/%s",srun.Data()));
  TGraphAsymmErrors* gRatioFakes = new TGraphAsymmErrors(hFakes,hTracks);
  TH1D* hDecays = trackCounters->Get("cent",Form("track:decay/%s",srun.Data()));
  TGraphAsymmErrors* gRatioDecays = new TGraphAsymmErrors(hDecays,hTracks);
  TH1D* hTracks_trig = trackCounters->Get("cent",Form("track:reconstructed/trig:yes/%s",srun.Data()));
  TH1D* hFakes_trig = trackCounters->Get("cent",Form("track:fake/trig:yes/%s",srun.Data()));
  TGraphAsymmErrors* gRatioFakes_trig = new TGraphAsymmErrors(hFakes_trig,hTracks_trig);
  TH1D* hDecays_trig = trackCounters->Get("cent",Form("track:decay/trig:yes/%s",srun.Data()));
  TGraphAsymmErrors* gRatioDecays_trig = new TGraphAsymmErrors(hDecays_trig,hTracks_trig);
  TH1D* hTracks_trig_acc = trackCounters->Get("cent",Form("track:reconstructed/trig:yes/acc:in/%s",srun.Data()));
  TH1D* hFakes_trig_acc = trackCounters->Get("cent",Form("track:fake/trig:yes/acc:in/%s",srun.Data()));
  TGraphAsymmErrors* gRatioFakes_trig_acc = new TGraphAsymmErrors(hFakes_trig_acc,hTracks_trig_acc);
  TH1D* hDecays_trig_acc = trackCounters->Get("cent",Form("track:decay/trig:yes/acc:in/%s",srun.Data()));
  TGraphAsymmErrors* gRatioDecays_trig_acc = new TGraphAsymmErrors(hDecays_trig_acc,hTracks_trig_acc);
  
  // get statistics of pairs versus centrality
  TH1D* hPairs = pairCounters->Get("cent",Form("pair:reconstructed/%s",srun.Data()));
  TH1D* h12Fakes = pairCounters->Get("cent",Form("pair:1fake,2fakes/%s",srun.Data()));
  TGraphAsymmErrors* gRatioP = new TGraphAsymmErrors(h12Fakes,hPairs);
  TH1D* hPairs_trig = pairCounters->Get("cent",Form("pair:reconstructed/trig:2/%s",srun.Data()));
  TH1D* h12Fakes_trig = pairCounters->Get("cent",Form("pair:1fake,2fakes/trig:2/%s",srun.Data()));
  TGraphAsymmErrors* gRatioP_trig = new TGraphAsymmErrors(h12Fakes_trig,hPairs_trig);
  TH1D* hPairs_trig_acc = pairCounters->Get("cent",Form("pair:reconstructed/trig:2/acc:in/%s",srun.Data()));
  TH1D* h12Fakes_trig_acc = pairCounters->Get("cent",Form("pair:1fake,2fakes/trig:2/acc:in/%s",srun.Data()));
  TGraphAsymmErrors* gRatioP_trig_acc = new TGraphAsymmErrors(h12Fakes_trig_acc,hPairs_trig_acc);
  
  // draw statistics of singles versus run
  TCanvas* cTracksVsRun = new TCanvas("cTracksVsRun", "cTracksVsRun", 900, 600);
  cTracksVsRun->SetFillColor(0);
  cTracksVsRun->Divide(3,2);
  cTracksVsRun->cd(1);
  gPad->SetLogy();
  hTracksVsRun->Draw();
  hFakesVsRun->SetLineColor(2);
  hFakesVsRun->SetFillColor(2);
  hFakesVsRun->SetFillStyle(3017);
  hFakesVsRun->Draw("sames");
  cTracksVsRun->cd(2);
  gPad->SetLogy();
  hTracksVsRun->Draw();
  hDecaysVsRun->SetLineColor(2);
  hDecaysVsRun->SetFillColor(2);
  hDecaysVsRun->SetFillStyle(3017);
  hDecaysVsRun->Draw("sames");
  cTracksVsRun->cd(3);
  gPad->SetLogy();
  hTrackablesVsRun->Draw();
  hMatchedVsRun->SetLineColor(2);
  hMatchedVsRun->Draw("sames");
  cTracksVsRun->cd(4);
  gRatioFakesVsRun->SetHistogram((TH1F*)hTracksVsRun->Clone());
  gRatioFakesVsRun->GetYaxis()->SetTitle("fakes / total");
  gRatioFakesVsRun->Draw("ap");
  cTracksVsRun->cd(5);
  gRatioDecaysVsRun->SetHistogram((TH1F*)hTracksVsRun->Clone());
  gRatioDecaysVsRun->GetYaxis()->SetTitle("decays / total");
  gRatioDecaysVsRun->Draw("ap");
  cTracksVsRun->cd(6);
  gEffVsRun->SetHistogram((TH1F*)hTrackablesVsRun->Clone());
  gEffVsRun->GetYaxis()->SetTitle("efficiency");
  gEffVsRun->Draw("ap");
  
  // draw statistics of fakes versus centrality
  TCanvas* cTracks = new TCanvas("cFakes", "cFakes", 900, 600);
  cTracks->SetFillColor(0);
  cTracks->Divide(3,2);
  cTracks->cd(1);
  gPad->SetLogy();
  hTracks->GetXaxis()->SetTitle("centrality class");
  hTracks->Draw();
  hFakes->SetLineColor(2);
  hFakes->SetFillColor(2);
  hFakes->SetFillStyle(3017);
  hFakes->Draw("sames");
  cTracks->cd(2);
  gPad->SetLogy();
  hTracks_trig->GetXaxis()->SetTitle("centrality class");
  hTracks_trig->Draw();
  hFakes_trig->SetLineColor(2);
  hFakes_trig->SetFillColor(2);
  hFakes_trig->SetFillStyle(3017);
  hFakes_trig->Draw("sames");
  cTracks->cd(3);
  gPad->SetLogy();
  hTracks_trig_acc->GetXaxis()->SetTitle("centrality class");
  hTracks_trig_acc->Draw();
  hFakes_trig_acc->SetLineColor(2);
  hFakes_trig_acc->SetFillColor(2);
  hFakes_trig_acc->SetFillStyle(3017);
  hFakes_trig_acc->Draw("sames");
  cTracks->cd(4);
  gRatioFakes->SetHistogram((TH1F*)hTracks->Clone());
  gRatioFakes->GetYaxis()->SetTitle("fakes / total");
  gRatioFakes->Draw("ap");
  cTracks->cd(5);
  gRatioFakes_trig->SetHistogram((TH1F*)hTracks_trig->Clone());
  gRatioFakes_trig->GetYaxis()->SetTitle("fakes (trig) / total");
  gRatioFakes_trig->Draw("ap");
  cTracks->cd(6);
  gRatioFakes_trig_acc->SetHistogram((TH1F*)hTracks_trig_acc->Clone());
  gRatioFakes_trig_acc->GetYaxis()->SetTitle("fakes (trig+acc) / total");
  gRatioFakes_trig_acc->Draw("ap");
  
  // draw statistics of decays versus centrality
  TCanvas* cTracks = new TCanvas("cDecays", "cDecays", 900, 600);
  cTracks->SetFillColor(0);
  cTracks->Divide(3,2);
  cTracks->cd(1);
  gPad->SetLogy();
  hTracks->GetXaxis()->SetTitle("centrality class");
  hTracks->Draw();
  hDecays->SetLineColor(2);
  hDecays->SetFillColor(2);
  hDecays->SetFillStyle(3017);
  hDecays->Draw("sames");
  cTracks->cd(2);
  gPad->SetLogy();
  hTracks_trig->GetXaxis()->SetTitle("centrality class");
  hTracks_trig->Draw();
  hDecays_trig->SetLineColor(2);
  hDecays_trig->SetFillColor(2);
  hDecays_trig->SetFillStyle(3017);
  hDecays_trig->Draw("sames");
  cTracks->cd(3);
  gPad->SetLogy();
  hTracks_trig_acc->GetXaxis()->SetTitle("centrality class");
  hTracks_trig_acc->Draw();
  hDecays_trig_acc->SetLineColor(2);
  hDecays_trig_acc->SetFillColor(2);
  hDecays_trig_acc->SetFillStyle(3017);
  hDecays_trig_acc->Draw("sames");
  cTracks->cd(4);
  gRatioDecays->SetHistogram((TH1F*)hTracks->Clone());
  gRatioDecays->GetYaxis()->SetTitle("decays / total");
  gRatioDecays->Draw("ap");
  cTracks->cd(5);
  gRatioDecays_trig->SetHistogram((TH1F*)hTracks_trig->Clone());
  gRatioDecays_trig->GetYaxis()->SetTitle("decays (trig) / total");
  gRatioDecays_trig->Draw("ap");
  cTracks->cd(6);
  gRatioDecays_trig_acc->SetHistogram((TH1F*)hTracks_trig_acc->Clone());
  gRatioDecays_trig_acc->GetYaxis()->SetTitle("decays (trig+acc) / total");
  gRatioDecays_trig_acc->Draw("ap");
  
  // draw statistics of pairs versus centrality
  TCanvas* cPairs = new TCanvas("cPairs", "cPairs", 900, 600);
  cPairs->SetFillColor(0);
  cPairs->Divide(3,2);
  cPairs->cd(1);
  gPad->SetLogy();
  hPairs->GetXaxis()->SetTitle("centrality class");
  hPairs->Draw();
  h12Fakes->SetLineColor(2);
  h12Fakes->SetFillColor(2);
  h12Fakes->SetFillStyle(3017);
  h12Fakes->Draw("sames");
  cPairs->cd(2);
  gPad->SetLogy();
  hPairs_trig->GetXaxis()->SetTitle("centrality class");
  hPairs_trig->Draw();
  h12Fakes_trig->SetLineColor(2);
  h12Fakes_trig->SetFillColor(2);
  h12Fakes_trig->SetFillStyle(3017);
  h12Fakes_trig->Draw("sames");
  cPairs->cd(3);
  gPad->SetLogy();
  hPairs_trig_acc->GetXaxis()->SetTitle("centrality class");
  hPairs_trig_acc->Draw();
  h12Fakes_trig_acc->SetLineColor(2);
  h12Fakes_trig_acc->SetFillColor(2);
  h12Fakes_trig_acc->SetFillStyle(3017);
  h12Fakes_trig_acc->Draw("sames");
  cPairs->cd(4);
  gRatioP->SetHistogram((TH1F*)hPairs->Clone());
  gRatioP->GetYaxis()->SetTitle("fakes / total");
  gRatioP->Draw("ap");
  cPairs->cd(5);
  gRatioP_trig->SetHistogram((TH1F*)hPairs_trig->Clone());
  gRatioP_trig->GetYaxis()->SetTitle("fakes (trig) / total");
  gRatioP_trig->Draw("ap");
  cPairs->cd(6);
  gRatioP_trig_acc->SetHistogram((TH1F*)hPairs_trig_acc->Clone());
  gRatioP_trig_acc->GetYaxis()->SetTitle("fakes (trig+acc) / total");
  gRatioP_trig_acc->Draw("ap");
  
  // draw efficiency versus centrality
  TCanvas* cEff = new TCanvas("cEff", "cEff", 600, 300);
  cEff->SetFillColor(0);
  cEff->Divide(2,1);
  cEff->cd(1);
  gPad->SetLogy();
  hTrackables->GetXaxis()->SetTitle("centrality class");
  hTrackables->Draw();
  hMatched->SetLineColor(2);
  hMatched->Draw("sames");
  //hMatchTrackable->SetLineColor(6);
  //hMatchTrackable->Draw("sames");
  cEff->cd(2);
  gEff->SetHistogram((TH1F*)hTrackables->Clone());
  gEff->GetYaxis()->SetTitle("efficiency");
  gEff->Draw("ap");
  //gEffTrackable->SetHistogram((TH1F*)hTrackables->Clone());
  //gEffTrackable->SetLineColor(4);
  //gEffTrackable->Draw("p");
  
}