/*  ScanHV.C
*
*  Created by Philippe Pillot on 01/12/17.
*  Copyright 2017 SUBATECH
*  This software is made available under the terms of the GNU GPL 3.0
*
*/

#include <Riostream.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TString.h>
#include <TCanvas.h>
#include <TMultiGraph.h>
#include <TGraph.h>
#include <TRegexp.h>
#include <TAxis.h>
#include <TTimeStamp.h>
#include <TLine.h>
#include <TList.h>

#include "AliMpDCSNamer.h"
#include "AliMUONTrackerHV.h"
#include "AliCDBManager.h"
#include "AliMUONCDB.h"
#include "AliMUONRecoParam.h"

Double_t yRange[2] = {-1., 1700.};

Int_t nRuns = 0;
UInt_t (*runBoundaries)[3] = 0x0;
Double_t hvMin = 1210.;
Double_t hvLimits[10];

void GetRunBoundaries(TGraph *gRunBoudaries, TString &runList);
void GetHVLimits(TString &runList, TString &ocdbPath);
Bool_t Find1400VIssues(TGraph *g, Int_t iCh);
void Print1400VIssues(TString dcsName, UInt_t t0, UInt_t t1, Bool_t first);
TString FindRuns(TString dcsName, UInt_t t0, UInt_t t1);
Bool_t isKnown(Int_t run, TString dcsName);
void DrawRunBoudaries(TCanvas *c2);

//----------------------------------------------------------------------------
void ScanHV(TString runList, TString ocdbPath = "raw://", Bool_t allIssues = kFALSE)
{
  /// Scan the HV of every sectors to check for "1400V issues".
  /// if allIssues = kTRUE, look for all cases where HV < RecoParam limit more than 15s
  
  AliMUONTrackerHV hv(runList.Data(), ocdbPath.Data());
  
  hv.Plot("");
  
  TCanvas *c = static_cast<TCanvas*>(gROOT->GetListOfCanvases()->Last());
  if (!c) return;
  TMultiGraph *mg = dynamic_cast<TMultiGraph*>(c->GetListOfPrimitives()->At(1));
  if (!mg) return;
  
  TCanvas *c2 = new TCanvas("HVPerCh", "HV status per chamber", 1200, 1000);
  c2->Divide(2,5);
  
  TGraph *gRunBoudaries = static_cast<TGraph*>(mg->GetListOfGraphs()->FindObject("runBoundaries"));
  GetRunBoundaries(gRunBoudaries, runList);
  
  GetHVLimits(runList, ocdbPath);
  if (allIssues) hvMin = -1.;
  
  if (!gSystem->AccessPathName("CheckHVLogs")) {
    printf("\n\e[0;31m!!! The search for know issues is done from log files in CheckHVLogs. !!!\e[0m\n");
    printf("\e[0;31m!!! If you want to recreate the logs (because of e.g. an OCDB update) !!!\e[0m\n");
    printf("\e[0;31m!!! you have to erase them first.                                     !!!\e[0m\n");
  }
  
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
      if (Find1400VIssues(g, iCh-1)) {
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
void GetHVLimits(TString &runList, TString &ocdbPath)
{
  /// get the HV limits per chamber from the RecoParam
  
  for (Int_t i = 0; i < 10; ++i) hvLimits[i] = -1.;
  
  ifstream inFile(runList.Data());
  if (!inFile.is_open()) return;
  TString run;
  do run.ReadLine(inFile,kTRUE);
  while (run.IsNull() && !inFile.eof());
  if (!run.IsDec()) return;
  inFile.close();
  
  AliCDBManager* cdbm = AliCDBManager::Instance();
  cdbm->SetDefaultStorage(ocdbPath.Data());
  cdbm->SetRun(run.Atoi());
  AliMUONRecoParam* recoParam = AliMUONCDB::LoadRecoParam();
  if (!recoParam) return;
  
  for (Int_t i = 0; i < 10; ++i) {
    hvLimits[i] = recoParam->HVLimit(i);
    printf("hvLimits[%d] = %f\n",i,hvLimits[i]);
  }

}

//----------------------------------------------------------------------------
Bool_t Find1400VIssues(TGraph *g, Int_t iCh)
{
  /// look for HV struggling at ~1400V.
  /// return kTRUE if problem found
  
  static const UInt_t minStrugglingTime = 15; // s
  
  Bool_t issues = kFALSE;
  
  Double_t *t = g->GetX();
  Double_t *hv = g->GetY();
  Int_t n = g->GetN();
  
  UInt_t t0 = 0;
  UInt_t dt = 0;
  for (Int_t i = 0; i < n; ++i) {
    if (hv[i] > hvMin && hv[i] < hvLimits[iCh]) {
      if (t0 < 1) t0 = t[i];
      else dt = t[i] - t0;
    } else {
      if (t0 > 0) dt = t[i] - t0;
      if (dt > minStrugglingTime) {
        Print1400VIssues(g->GetName(), t0, t0+dt, !issues);
        issues = kTRUE;
      }
      t0 = dt = 0;
    }
  }
  if (dt > minStrugglingTime) {
    Print1400VIssues(g->GetName(), t0, t0+dt, !issues);
    issues = kTRUE;
  }
  
  return issues;
  
}

//----------------------------------------------------------------------------
void Print1400VIssues(TString dcsName, UInt_t t0, UInt_t t1, Bool_t first)
{
  /// print "1400V issues"
  if (first) printf("Problem found for %s:\n", dcsName.Data());
  TTimeStamp st0(t0);
  TTimeStamp st1(t1);
  TString runs = FindRuns(dcsName, t0, t1);
  printf("- between %s and %s --> run(s) %s\n", st0.AsString("s"), st1.AsString("s"), runs.Data());
}

//----------------------------------------------------------------------------
TString FindRuns(TString dcsName, UInt_t t0, UInt_t t1)
{
  /// Find the list of affected runs in this time range
  
  TString affectedRuns = "";
  
  if (!runBoundaries) return affectedRuns;
  //Int_t n = (Int_t)(sizeof(runBoundaries)/sizeof(UInt_t[3]));

  for (Int_t i = 0; i < nRuns; ++i) {
    
    if (runBoundaries[i][2] <= t0 || runBoundaries[i][1] >= t1) continue;
    
    if (isKnown(runBoundaries[i][0], dcsName)) affectedRuns += runBoundaries[i][0];
    else affectedRuns += TString::Format("\e[0;31;103m%d\e[0m",runBoundaries[i][0]);
    affectedRuns += ",";
    
  }
  
  affectedRuns.Remove(TString::kTrailing, ',');
  return affectedRuns;
  
}

//----------------------------------------------------------------------------
Bool_t isKnown(Int_t run, TString dcsName)
{
  /// return kTRUE if the problematic sector is already identified
  
  static AliMpDCSNamer dcsHelper;
  
  if (gSystem->AccessPathName("CheckHVLogs")) gSystem->Exec("mkdir -p CheckHVLogs");
  
  TString log(TString::Format("CheckHVLogs/%d.log",run));
  if (gSystem->AccessPathName(log.Data()))
    gROOT->ProcessLineSync(TString::Format("AliMUONCDB::CheckHV(%d); > %s",run,log.Data()));
  
  TString dcsAlias = dcsHelper.DCSAliasFromName(dcsName);
  return (!gSystem->Exec(TString::Format("grep -c %s %s > /dev/null",dcsAlias.Data(),log.Data())));

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
