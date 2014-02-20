/*
 *  ComputeTriggerEff.C
 *  aliroot_dev
 *
 *  Created by philippe pillot on 24/11/11.
 *  Copyright 2011 SUBATECH. All rights reserved.
 *
 */

#if !defined(__CINT__) || defined(__MAKECINT__)

#include <TROOT.h>
#include <TStyle.h>
#include <TMath.h>
#include <TString.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH2D.h>
#include <TH1D.h>

#include "AliCounterCollection.h"

#endif

Double_t intrinsicTrigEff = 1.;
//Double_t intrinsicTrigEff = 0.76;
//Double_t intrinsicTrigEff = 0.44;

void controlPlots(AliCounterCollection *eventCounters, AliCounterCollection *trigCounters,
		  TString data, TH2D *&hOccVsCentPerEvent, TH1D *&hOccPerNTrig);
void efficiencyPlots(TH2D *hOccVsCentPerEvent, TH1D *hOccPerNTrig);

void ComputeTriggerEff(TString fileNameMB = "AnalysisResults.root", TString fileNameJPsi = "AnalysisResults.root")
{
  /// compute the trigger efficiency from the local board occupancies
  /// to be compiled --> you must load the environment before:
  /*
  gROOT->LoadMacro("$WORK/Macros/Facilities/runTaskFacilities.C");
  LoadAlirootLocally("PWG3base", "", "");
  */
  
  gStyle->SetOptStat(0);
  
  // open file containing MB data
  TFile* fileMB = TFile::Open(fileNameMB.Data(),"READ");
  if (!fileMB || !fileMB->IsOpen()) {
    Error("ComputeTriggerEff", "cannot open MB file");
    return;
  }
  
  // event counters in MB data
  AliCounterCollection *eventCountersMB = static_cast<AliCounterCollection*>(fileMB->FindObjectAny("eventCounters"));
  if (!eventCountersMB) {
    Error("ComputeTriggerEff", "cannot find eventCounters object for MB");
    return;
  }
  printf("\n# MB events = %d\n",(Int_t)eventCountersMB->GetSum("event:all"));
  
  // trigger counters in MB data
  AliCounterCollection *trigCountersMB = static_cast<AliCounterCollection*>(fileMB->FindObjectAny("trigCounters"));
  if (!trigCountersMB) {
    Error("ComputeTriggerEff", "cannot find trigCounters object for MB");
    return;
  }
  trigCountersMB->Sort("board", kTRUE);
  trigCountersMB->Sort("ntrig", kTRUE);
  printf("# trigger track in MB events = %d\n",(Int_t)trigCountersMB->GetSum());
  
  // display control plots
  TH2D *hOccVsCentPerEvent = 0x0;
  TH1D *hDummy1 = 0x0;
  controlPlots(eventCountersMB, trigCountersMB, "MB", hOccVsCentPerEvent, hDummy1);
  
  // open file containing JPsi data
  TFile* fileJPsi = TFile::Open(fileNameJPsi.Data(),"READ");
  if (!fileJPsi || !fileJPsi->IsOpen()) {
    Error("ComputeTriggerEff", "cannot open JPsi file");
    return;
  }
  
  // event counters in JPsi data
  AliCounterCollection *eventCountersJPsi = static_cast<AliCounterCollection*>(fileJPsi->FindObjectAny("eventCounters"));
  if (!eventCountersJPsi) {
    Error("ComputeTriggerEff", "cannot find eventCounters object for JPsi");
    return;
  }
  printf("\n# JPsi events = %d\n",(Int_t)eventCountersJPsi->GetSum("event:all"));
  
  // trigger counters in JPsi data
  AliCounterCollection *trigCountersJPsi = static_cast<AliCounterCollection*>(fileJPsi->FindObjectAny("trigCounters"));
  if (!trigCountersJPsi) {
    Error("ComputeTriggerEff", "cannot find trigCounters object for JPsi");
    return;
  }
  trigCountersJPsi->Sort("board", kTRUE);
  trigCountersJPsi->Sort("ntrig", kTRUE);
  printf("# trigger track in JPsi events = %d\n",(Int_t)trigCountersJPsi->GetSum());
  
  // display control plots
  TH2D *hDummy2 = 0x0;
  TH1D *hOccPerNTrig = 0x0;
  controlPlots(eventCountersJPsi, trigCountersJPsi, "JPsi", hDummy2, hOccPerNTrig);
  
  // compute and display trigger efficiency
  efficiencyPlots(hOccVsCentPerEvent, hOccPerNTrig);
  
}

//--------------------------------------------------------------------------
void controlPlots(AliCounterCollection *eventCounters, AliCounterCollection *trigCounters,
		  TString data, TH2D *&hOccVsCentPerEvent, TH1D *&hOccPerNTrig)
{
  /// draw some control plots
  
  // get number of events per centrality class
  TH1D *hNEventVsCent = eventCounters->Get("cent","event:all");
  hNEventVsCent->SetNameTitle(Form("hN%sVsCent",data.Data()),Form("hN%sVsCent",data.Data()));
  TArrayI nEventVsCent(hNEventVsCent->GetXaxis()->GetNbins());
  for (Int_t icent = 1; icent <= hNEventVsCent->GetXaxis()->GetNbins(); icent++)
    nEventVsCent[icent-1] = (Int_t) TMath::Max(hNEventVsCent->GetBinContent(icent), 1.);
  new TCanvas;
  hNEventVsCent->Draw();
  
  // plot the local board occupancy vs centrality per events
  hOccVsCentPerEvent = trigCounters->Get("cent","board","");
  hOccVsCentPerEvent->SetNameTitle(Form("hoccVsCentPer%s",data.Data()),Form("hoccVsCentPer%s",data.Data()));
  // normalize to the number of MB events per centrality bin
  for (Int_t icent = 1; icent <= hOccVsCentPerEvent->GetYaxis()->GetNbins(); icent++) {
    for (Int_t iboard = 1; iboard <= hOccVsCentPerEvent->GetXaxis()->GetNbins(); iboard++) {
      hOccVsCentPerEvent->SetBinContent(iboard, icent, hOccVsCentPerEvent->GetBinContent(iboard, icent) / nEventVsCent[icent-1]);
    }
  }
  new TCanvas;
  hOccVsCentPerEvent->Draw("colz");
  
  // get the number of trigger tracks per centrality class and per value of number of trigger tracks per event
  TH2D *hNTrigVsCent = trigCounters->Get("ntrig","cent","");
  hNTrigVsCent->SetNameTitle(Form("hNTrigVsCent_%s",data.Data()),Form("hNTrigVsCent_%s",data.Data()));
  TArrayI nTrigVsCent(hNTrigVsCent->GetXaxis()->GetNbins());
  TArrayI nTrigVsNTrig(hNTrigVsCent->GetYaxis()->GetNbins());
  for (Int_t intrig = 1; intrig <= hNTrigVsCent->GetYaxis()->GetNbins(); intrig++) {
    nTrigVsNTrig[intrig-1] = (Int_t) TMath::Max(hNTrigVsCent->GetBinContent(1, intrig), 1.);
    TString label = hNTrigVsCent->GetYaxis()->GetBinLabel(intrig);
    Int_t nTrig = label.Atoi();
    for (Int_t icent = 1; icent <= hNTrigVsCent->GetXaxis()->GetNbins(); icent++) {
      nTrigVsCent[icent-1] += (Int_t) hNTrigVsCent->GetBinContent(icent, intrig);
      hNTrigVsCent->SetBinContent(icent, intrig, hNTrigVsCent->GetBinContent(icent, intrig) / nTrig);
    }
  }
  for (Int_t icent = 1; icent <= hNTrigVsCent->GetXaxis()->GetNbins(); icent++)
    nTrigVsCent[icent-1] = TMath::Max(nTrigVsCent[icent-1], 1);
  new TCanvas;
  hNTrigVsCent->Draw("colz");
  
  // plot the local board occupancy vs centrality normalized to the number of trigger tracks
  TH2D *hOccVsCentPerNTrig = trigCounters->Get("cent","board","");
  hOccVsCentPerNTrig->SetNameTitle(Form("hOccVsCentPerNTrig_%s",data.Data()),Form("hOccVsCentPerNTrig_%s",data.Data()));
  // normalize to the number of MB events per centrality bin
  for (Int_t icent = 1; icent <= hOccVsCentPerNTrig->GetYaxis()->GetNbins(); icent++) {
    for (Int_t iboard = 1; iboard <= hOccVsCentPerNTrig->GetXaxis()->GetNbins(); iboard++) {
      hOccVsCentPerNTrig->SetBinContent(iboard, icent, hOccVsCentPerNTrig->GetBinContent(iboard, icent) / nTrigVsCent[icent-1]);
    }
  }
  new TCanvas;
  hOccVsCentPerNTrig->Draw("colz");
  
  // plot the local board occupancy vs number of trigger tracks normalized to the number of trigger tracks
  TH2D *hOccVsNTrig = trigCounters->Get("ntrig","board","");
  hOccVsNTrig->SetNameTitle(Form("hOccVsNTrig_%s",data.Data()),Form("hOccVsNTrig_%s",data.Data()));
  // normalize to the number of trigger tracks
  for (Int_t intrig = 1; intrig <= hOccVsNTrig->GetYaxis()->GetNbins(); intrig++) {
    for (Int_t iboard = 1; iboard <= hOccVsNTrig->GetXaxis()->GetNbins(); iboard++) {
      hOccVsNTrig->SetBinContent(iboard, intrig, hOccVsNTrig->GetBinContent(iboard, intrig) / nTrigVsNTrig[intrig-1]);
    }
  }
  new TCanvas;
  hOccVsNTrig->Draw("colz");
  
  // plot the local board occupancy per trigger tracks
  Double_t nTrigTot = trigCounters->GetSum();
  hOccPerNTrig = trigCounters->Get("board","");
  hOccPerNTrig->SetNameTitle(Form("hOccPerNTrig_%s",data.Data()),Form("hOccPerNTrig_%s",data.Data()));
  if (nTrigTot > 0.) hOccPerNTrig->Scale(1./nTrigTot);
  new TCanvas;
  hOccPerNTrig->Draw();
  
}

//--------------------------------------------------------------------------
void efficiencyPlots(TH2D *hOccVsCentPerEvent, TH1D *hOccPerNTrig)
{
  /// compute and display trigger efficiency
  
  // plot the local board occupancy per trigger tracks corrected for the trigger efficiency versus centrality
  TH2D *hCoorOcc = (TH2D*) hOccVsCentPerEvent->Clone("hCoorOcc");
  hCoorOcc->SetNameTitle("hCoorOcc","hCoorOcc");
  for (Int_t icent = 1; icent <= hCoorOcc->GetYaxis()->GetNbins(); icent++) {
    for (Int_t iboard = 1; iboard <= hOccPerNTrig->GetXaxis()->GetNbins(); iboard++) {
      hCoorOcc->SetBinContent(iboard, icent, hOccPerNTrig->GetBinContent(iboard) * (1.-hOccVsCentPerEvent->GetBinContent(iboard, icent)));
    }
  }
  new TCanvas;
  hCoorOcc->Draw("colz");
  
  // plot the trigger efficiency versus centrality
  TH1D *hEffVsCent = hCoorOcc->ProjectionY("hEffVsCent", 1, hCoorOcc->GetXaxis()->GetNbins());
  hEffVsCent->SetNameTitle("hEffVsCent","hEffVsCent");
  hEffVsCent->Scale(intrinsicTrigEff);
  new TCanvas;
  hEffVsCent->Draw();
  
}

