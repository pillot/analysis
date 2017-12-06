/*  ScanHV.C
*
*  Created by Philippe Pillot on 01/12/17.
*  Copyright 2017 SUBATECH
*  This software is made available under the terms of the GNU GPL 3.0
*
*/

#include <Riostream.h>
#include <TROOT.h>
#include <TString.h>
#include <TCanvas.h>
#include <TMultiGraph.h>
#include <TGraph.h>
#include <TRegexp.h>
#include <TAxis.h>
#include <TTimeStamp.h>
#include <TLine.h>
#include <TList.h>

#include "AliMUONTrackerHV.h"

Double_t yRange[2] = {-1., 1700.};

Int_t nRuns = 0;
UInt_t (*runBoundaries)[3] = 0x0;

void GetRunBoundaries(TGraph *gRunBoudaries, TString &runList);
Bool_t Find1400VIssues(TGraph *g);
void Print1400VIssues(TString dcsName, UInt_t t0, UInt_t t1);
TString FindRuns(UInt_t t0, UInt_t t1);
void DrawRunBoudaries(TCanvas *c2);

//----------------------------------------------------------------------------
void ScanHV(TString runList, TString ocdbPath = "raw://")
{
  /// Scan the HV of every sectors to check for "1400V issues"
  
  AliMUONTrackerHV hv(runList.Data(), ocdbPath.Data());
  
  hv.Plot("",kFALSE,kTRUE);
  
  TCanvas *c = static_cast<TCanvas*>(gROOT->GetListOfCanvases()->Last());
  if (!c) return;
  TMultiGraph *mg = dynamic_cast<TMultiGraph*>(c->GetListOfPrimitives()->At(1));
  if (!mg) return;
  
  TCanvas *c2 = new TCanvas("HVPerCh", "HV status per chamber", 1200, 1000);
  c2->Divide(2,5);
  
  TGraph *gRunBoudaries = static_cast<TGraph*>(mg->GetListOfGraphs()->FindObject("runBoundaries"));
  GetRunBoundaries(gRunBoudaries, runList);
  
  TIter nextg(mg->GetListOfGraphs());
  TGraph *g = 0x0;
  TRegexp reChName("Chamber..");
  Bool_t issues = kFALSE;
  printf("\n------ list of issues ------\n");
  while ((g = static_cast<TGraph*>(nextg()))) {
    TString gName = g->GetName();
    TString chName = gName(reChName);
    Int_t iCh = chName.Remove(0,7).Atoi();
    if (iCh > 0) {
      c2->cd(iCh);
      if (gPad->GetListOfPrimitives()->GetEntries() == 0) {
        g->SetTitle(TString::Format("chamber %d", iCh));
        g->Draw("alp");
        g->GetXaxis()->SetTimeDisplay(1);
        g->GetXaxis()->SetTimeFormat("%d/%m %H:%M");
        g->GetXaxis()->SetTimeOffset(0,"gmt");
        g->GetXaxis()->SetNdivisions(505);
        g->SetMinimum(yRange[0]);
        g->SetMaximum(yRange[1]);
      } else g->Draw("lp");
      if (Find1400VIssues(g)) {
        issues = kTRUE;
        printf("----------------------------\n");
      }
    }
  }
  if (!issues) printf("----------------------------\n");
  printf("\n");

  DrawRunBoudaries(c2);
  
}

//----------------------------------------------------------------------------
void GetRunBoundaries(TGraph *gRunBoudaries, TString &runList)
{
  /// return the time boundaries of each run
  
  nRuns = 0;
  if (runBoundaries) delete[] runBoundaries;
  runBoundaries = 0x0;
  
  if (!gRunBoudaries) return;
  
  Double_t *t = gRunBoudaries->GetX();
  Int_t n = gRunBoudaries->GetN();
  if (n < 2 || n%2 != 0) return;
  
  ifstream inFile(runList.Data());
  if (!inFile.is_open()) return;
  
  nRuns = n/2;
  runBoundaries = new UInt_t[nRuns][3];
  
  TString run;
  Bool_t ok = kTRUE;
  for (Int_t i = 0; i < nRuns; ++i) {
    
    if (inFile.eof()) {
      printf("!!! inconsistent run list !!!\n");
      ok = kFALSE;
      break;
    }
    
    run.ReadLine(inFile,kTRUE);
    if (run.IsNull()) {
      printf("!!! invalid run !!!\n");
      ok = kFALSE;
      break;
    }
    
    runBoundaries[i][0] = run.Atoi();
    runBoundaries[i][1] = t[2*i];
    runBoundaries[i][2] = t[2*i+1];
    
  }
  
  run.ReadLine(inFile,kTRUE);
  if (!run.IsNull() || !inFile.eof()) {
    printf("!!! inconsistent run list !!!\n");
    ok = kFALSE;
  }
  
  inFile.close();
  
  if (!ok) {
    nRuns = 0;
    delete[] runBoundaries;
    runBoundaries = 0x0;
  }

}

//----------------------------------------------------------------------------
Bool_t Find1400VIssues(TGraph *g)
{
  /// look for HV struggling at ~1400V.
  /// return kTRUE if problem found
  
  static const UInt_t minStrugglingTime = 60; // s
  
  Bool_t issues = kFALSE;
  
  Double_t *t = g->GetX();
  Double_t *hv = g->GetY();
  Int_t n = g->GetN();
  
  UInt_t t0 = 0;
  UInt_t dt = 0;
  for (Int_t i = 0; i < n; ++i) {
    if (hv[i] > 1250. && hv[i] < 1550.) {
      if (t0 < 1) t0 = t[i];
      else dt = t[i] - t0;
    } else {
      if (dt > minStrugglingTime) {
        issues = kTRUE;
        Print1400VIssues(g->GetName(), t0, t0+dt);
      }
      t0 = dt = 0;
    }
  }
  if (dt > minStrugglingTime) {
    issues = kTRUE;
    Print1400VIssues(g->GetName(), t0, t0+dt);
  }
  
  return issues;
  
}

//----------------------------------------------------------------------------
void Print1400VIssues(TString dcsName, UInt_t t0, UInt_t t1)
{
  /// print "1400V issues"
  TTimeStamp st0(t0);
  TTimeStamp st1(t1);
  TString runs = FindRuns(t0, t1);
  printf("Problem found for %s between %s and %s --> run(s) %s\n", dcsName.Data(), st0.AsString("s"), st1.AsString("s"), runs.Data());
}

//----------------------------------------------------------------------------
TString FindRuns(UInt_t t0, UInt_t t1)
{
  /// Find the list of affected runs in this time range
  
  TString affectedRuns = "";
  
  if (!runBoundaries) return affectedRuns;
  //Int_t n = (Int_t)(sizeof(runBoundaries)/sizeof(UInt_t[3]));

  for (Int_t i = 0; i < nRuns; ++i) {
    
    if (runBoundaries[i][2] <= t0 || runBoundaries[i][1] >= t1) continue;
    
    affectedRuns += runBoundaries[i][0];
    affectedRuns += ",";
    
  }
  
  affectedRuns.Remove(TString::kTrailing, ',');
  return affectedRuns;
  
}

//----------------------------------------------------------------------------
void DrawRunBoudaries(TCanvas *c2)
{
  /// Draw the run time boundaries
  
  if (!runBoundaries) return;

  for (Int_t i = 0; i < nRuns; ++i) {
    
    TLine startRunLine(runBoundaries[i][1],yRange[0],runBoundaries[i][1],yRange[1]);
    startRunLine.SetUniqueID(runBoundaries[i][0]);
    startRunLine.SetLineColor(4);
    startRunLine.SetLineWidth(1);
    
    TLine endRunLine(runBoundaries[i][2],yRange[0],runBoundaries[i][2],yRange[1]);
    endRunLine.SetUniqueID(runBoundaries[i][0]);
    endRunLine.SetLineColor(2);
    endRunLine.SetLineWidth(1);
    
    for (Int_t j = 1; j <= 10; ++j) {
      c2->cd(j);
      startRunLine.Clone()->Draw();
      endRunLine.Clone()->Draw();
    }
    
  }
  
}
