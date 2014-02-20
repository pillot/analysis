/*
 *  plotMuonPerformance.C
 *  aliroot_dev
 *
 *  Created by philippe pillot on 05/04/11.
 *  Copyright 2011 Subatech. All rights reserved.
 *
 */

#include <TH1.h>
#include <TH2.h>
#include <TAxis.h>
#include <TString.h>
#include <Riostream.h>
#include <TFile.h>
#include <TList.h>
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TMath.h>
#include <TArrayD.h>
#include <TStyle.h>

void plotMuonEfficiencyVsRun(TString runList = "");
void integratedEfficiency();

//---------------------------------------------------------------------------
void plotMuonPerformance(TString runList = "")
{
  /// plot muon performances (like efficiency) integrated and versus run
  /// !!! to be compiled because of CINT !!!
  
  plotMuonEfficiencyVsRun(runList);
  integratedEfficiency();
  
}

//---------------------------------------------------------------------------
void plotMuonEfficiencyVsRun(TString runList)
{
  /// plot the tracking efficiency versus run
  /// !!! to be compiled because of CINT !!!
  
  gStyle->SetFillColor(0);
  
  // output graphs
  TGraphAsymmErrors *trackingEffVsRun = new TGraphAsymmErrors;
  
  TArrayD chambersEff(11);
  TArrayD chambersEffErr[2];
  chambersEffErr[0].Set(11);
  chambersEffErr[1].Set(11);
  
  ifstream inFile(runList.Data());
  if (!inFile.is_open()) {
    printf("cannot open file %s\n",runList.Data());
    return;
  }
  
  // loop over runs
  Int_t irun = 0;
  TString currRun;
  TList runs;
  runs.SetOwner();
  while (!inFile.eof()) {
    
    // get current run number
    currRun.ReadLine(inFile,kTRUE);
    if(currRun.IsNull()) continue;
    
    // get input hists
    TH1F *effSummary = 0x0;
    TFile *file = new TFile(Form("runs/%d/AnalysisResults.root",currRun.Atoi()), "read");
//    TFile *file = new TFile(Form("runs/%d/AnalysisResults_wcut.root",currRun.Atoi()), "read");
    TList *effList = 0x0;
    if (!file || !file->IsOpen()) {
      printf("cannot open file runs/%d/AnalysisResults.root\n",currRun.Atoi());
    } else {
      effList = static_cast<TList*>(file->FindObjectAny("Efficiency"));
      if (!effList) {
	printf("cannot find Efficiency object in runs/%d/AnalysisResults.root\n",currRun.Atoi());
      } else {
	effSummary = static_cast<TH1F*>(effList->FindObject("effSummary"));
	if (!effSummary) {
	  printf("cannot find Efficiency object in runs/%d/AnalysisResults.root\n",currRun.Atoi());
	}
      }
    }
    
    // fill graph
    if (effSummary) {
      printf("%s %f ± %f\n",currRun.Data(), effSummary->GetBinContent(2), effSummary->GetBinError(2));
      trackingEffVsRun->SetPoint(irun,irun,effSummary->GetBinContent(2));
      trackingEffVsRun->SetPointError(irun++,0.,0.,effSummary->GetBinError(2),effSummary->GetBinError(2));
    } else {
      trackingEffVsRun->SetPoint(irun,irun,-1);
      trackingEffVsRun->SetPointError(irun++,0.,0.,0.,0.);
    }
    runs.AddLast(new TObjString(currRun));
    
    file->Close();
  }
  inFile.close();
  
  // set bin labels
  trackingEffVsRun->GetXaxis()->Set(irun, -0.5, irun-0.5);
  TIter nextRun(&runs);
  TObjString *srun = 0x0;
  irun = 1;
  while ((srun = static_cast<TObjString*>(nextRun())))
    trackingEffVsRun->GetXaxis()->SetBinLabel(irun++, srun->GetName());
  
  // display
  new TCanvas("cRealTrackingEffVsRun", "Real tracking efficiency versus run",1000,400);
  trackingEffVsRun->SetName("realTrackingEffVsRun");
  trackingEffVsRun->SetTitle("Real tracking efficiency versus run");
  trackingEffVsRun->SetLineStyle(1);
  trackingEffVsRun->SetLineColor(1); 
  trackingEffVsRun->SetMarkerStyle(20);
  trackingEffVsRun->SetMarkerSize(0.7);
  trackingEffVsRun->SetMarkerColor(4);
  trackingEffVsRun->GetXaxis()->SetTitle("Run #");
  trackingEffVsRun->GetXaxis()->SetLabelFont(22);
  trackingEffVsRun->GetXaxis()->SetLabelSize(0.05);
  trackingEffVsRun->GetXaxis()->SetTitleFont(22);
  trackingEffVsRun->GetXaxis()->SetTitleSize(0.05);
  trackingEffVsRun->GetYaxis()->SetTitle("Efficiency");
  trackingEffVsRun->GetYaxis()->SetLabelFont(22);
  trackingEffVsRun->GetYaxis()->SetLabelSize(0.06);
  trackingEffVsRun->GetYaxis()->SetTitleFont(22);
  trackingEffVsRun->GetYaxis()->SetTitleSize(0.05);
  trackingEffVsRun->GetYaxis()->SetTitleOffset(0.7);
  trackingEffVsRun->GetYaxis()->SetRangeUser(0.5, 1.05);
  trackingEffVsRun->Draw("ap");
  
  // save output
  TFile* file = new TFile("performance.root","update");
  trackingEffVsRun->Write(0x0, TObject::kOverwrite);
  file->Close();
  
}

//---------------------------------------------------------------------------
void integratedEfficiency()
{
  /// compute the integrated tracking efficiency in 0-80%
  /// !!! to be compiled because of CINT !!!
  
  // get input hists
  TFile *file = new TFile("AnalysisResults.root", "read");
  if (!file || !file->IsOpen()) {
    printf("cannot open file ./AnalysisResults.root\n");
    return;
  }
  TList *effList = static_cast<TList*>(file->FindObjectAny("Efficiency"));
  TH1F *effSummary = static_cast<TH1F*>(effList->FindObject("effSummary"));
  
  // print results
  cout << "Real tracking efficiency : " << effSummary->GetBinContent(2) << " ± " << effSummary->GetBinError(2) << endl;
  
  // close input file
  file->Close();
  
}

//---------------------------------------------------------------------------
void ratio(TString measEff, TString realEff)
{
  /// compute the integrated tracking efficiency in 0-80%
  /// !!! to be compiled because of CINT !!!
  
  gStyle->SetFillColor(0);
  
  Bool_t measVsReal = measEff.Contains("efficiency");
  
  // get measured efficiencies
  TFile *fEff = TFile::Open(measEff.Data());
  if (!fEff || !fEff->IsOpen()) {
    printf("cannot open file %s\n", measEff.Data());
    return;
  }
  TGraphAsymmErrors *gEff = measVsReal ?
    static_cast<TGraphAsymmErrors*>(fEff->FindObjectAny("trackingEffVsRun")) :
    static_cast<TGraphAsymmErrors*>(fEff->FindObjectAny("realTrackingEffVsRun"));
  if (!gEff) {
    printf("measured efficiencies per run not found\n");
    return;
  }
  
  // get real efficiencies
  TFile *fPer = TFile::Open(realEff.Data());
  if (!fPer || !fPer->IsOpen()) {
    printf("cannot open file %s\n", realEff.Data());
    return;
  }
  TGraphAsymmErrors *gPer = static_cast<TGraphAsymmErrors*>(fPer->FindObjectAny("realTrackingEffVsRun"));
  if (!gPer) {
    printf("real efficiencies per run not found\n");
    return;
  }
  
  // compute ratio measured/real
  Int_t nBins = gEff->GetN();
  TGraphErrors *ratio = new TGraphErrors(nBins);
  Double_t x,effM,effR,effMErr,effRErr,rat,ratErr;
  for (Int_t i = 0; i < nBins; i++) {
    gEff->GetPoint(i,x,effM);
    gPer->GetPoint(i,x,effR);
    effMErr = gEff->GetErrorY(i);
    effRErr = gPer->GetErrorY(i);
    if (effM > 0. && effR > 0.) {
      rat = effM/effR;
      ratErr = rat*TMath::Sqrt(effMErr*effMErr/effM*effM + effRErr*effRErr/effR*effR);
    } else {
      rat = 0.;
      ratErr = 0.;
    }
    ratio->SetPoint(i,x,rat);
    ratio->SetPointError(i,0.,ratErr);
  }
  
  
  // set bin labels
  ratio->GetXaxis()->Set(nBins, -0.5, nBins-0.5);
  for (Int_t i = 1; i <= nBins; i++)
    ratio->GetXaxis()->SetBinLabel(i, gEff->GetXaxis()->GetBinLabel(i));
  
  // display ratio
  TCanvas *c = measVsReal ?
    new TCanvas("cMeasOverReal", "measured over real efficiency versus run",1200,600) :
    new TCanvas("cNewOverOld", "new over old efficiency versus run",1200,600);
  c->Divide(1,2);
  c->cd(1);
  gEff->SetMarkerColor(2);
  gEff->SetTitle("tracking efficiency versus run");
  gEff->DrawClone("ap");
  gPer->SetMarkerColor(4);
  gPer->DrawClone("psame");
  c->cd(2);
  if (measVsReal) {
    ratio->SetName("MeasOverReal");
    ratio->SetTitle("measured over real efficiency versus run");
    ratio->GetYaxis()->SetTitle("measured / real efficiency");
  } else {
    ratio->SetName("NewOverOld");
    ratio->SetTitle("new over old efficiency versus run");
    ratio->GetYaxis()->SetTitle("new / old efficiency");
  }
  ratio->SetLineStyle(1);
  ratio->SetLineColor(1); 
  ratio->SetMarkerStyle(20);
  ratio->SetMarkerSize(0.7);
  ratio->SetMarkerColor(1);
  ratio->GetXaxis()->SetTitle("Run #");
  ratio->GetXaxis()->SetLabelFont(22);
  ratio->GetXaxis()->SetLabelSize(0.05);
  ratio->GetXaxis()->SetTitleFont(22);
  ratio->GetXaxis()->SetTitleSize(0.05);
  ratio->GetYaxis()->SetLabelFont(22);
  ratio->GetYaxis()->SetLabelSize(0.06);
  ratio->GetYaxis()->SetTitleFont(22);
  ratio->GetYaxis()->SetTitleSize(0.05);
  ratio->GetYaxis()->SetTitleOffset(0.7);
  ratio->SetMinimum(1.);
  ratio->SetMaximum(1.4);
  ratio->Draw("ap");
  
  // close input files
  fEff->Close();
  fPer->Close();
  
}

