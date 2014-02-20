/*
 *  plotTriggerRatio.C
 *  aliroot_dev
 *
 *  Created by philippe pillot on 16/10/12.
 *  Copyright 2012 Subatech. All rights reserved.
 *
 */


#if !defined(__CINT__) || defined(__MAKECINT__)

#include <Riostream.h>
#include "TFile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TH1F.h"
#include "TROOT.h"
#include "TInterpreter.h"
#include "TGraphAsymmErrors.h"
#include "AliCFContainer.h"
#include "AliCFGridSparse.h"
#include "AliCFEffGrid.h"

#endif

void plotTriggerRatio(TString inputFile = "AnalysisResults.root")
{
  /// plot the ratios LpT/ApT and HpT/LpT
  /*
  gROOT->LoadMacro("$WORK/Macros/Facilities/runTaskFacilities.C");
  LoadAlirootLocally("CORRFW","","")
  .x $WORK/Macros/MuonPerformance/plotTriggerRatio.C+
  */
  
  gROOT->SetStyle("Plain");
  gROOT->ForceStyle();
  
  // prepare environment
  if (!gInterpreter->IsLoaded("$WORK/Macros/Facilities/runTaskFacilities.C")) {
    gROOT->LoadMacro("$WORK/Macros/Facilities/runTaskFacilities.C");
    TString extraLibs = "CORRFW";
    gROOT->ProcessLineFast(Form("LoadAlirootLocally(\"%s\",\"\",\"\");",extraLibs.Data()));
  }
  
  TFile *file = new TFile(inputFile.Data(), "read");
  if (!file || !file->IsOpen()) {
    printf("cannot open file ./AnalysisResults.root\n");
    return;
  }
  AliCFContainer* cfContainer = static_cast<AliCFContainer*>(file->FindObjectAny("EffContainer"));
  if (!cfContainer) return;
  AliCFEffGrid* efficiency = new AliCFEffGrid("eff","",*cfContainer);
  if (!efficiency) return;
  efficiency->CalculateEfficiency(0, 1);
  
  efficiency->GetNum()->SetRangeUser(5, 1., 1.); // tracker
  efficiency->GetNum()->SetRangeUser(8, 1., 3.); // match MC
//  efficiency->GetNum()->SetRangeUser(9, 1., 4.); // MC trig level
  efficiency->GetNum()->SetRangeUser(1, -3.999, -2.501); // eta
  efficiency->GetNum()->SetRangeUser(3, 1., 2.); // theta
//  efficiency->GetNum()->SetRangeUser(10, 0., 10.); // centrality
  TH1 *hAll = efficiency->GetNum()->Project(0); // project over pT
  
  efficiency->GetNum()->SetRangeUser(6, 2., 4.); // trigger
  TH1 *hApT = efficiency->GetNum()->Project(0); // project over pT
  
  efficiency->GetNum()->SetRangeUser(6, 3., 4.); // trigger
  TH1 *hLpT = efficiency->GetNum()->Project(0); // project over pT
  
  efficiency->GetNum()->SetRangeUser(6, 4., 4.); // trigger
  TH1 *hHpT = efficiency->GetNum()->Project(0); // project over pT
  
  TGraphAsymmErrors *hLpTOverAll = new TGraphAsymmErrors(hLpT, hAll);
  hLpTOverAll->SetNameTitle("hLpTOverAll","LpT/All");
  TGraphAsymmErrors *hLpTOverApT = new TGraphAsymmErrors(hLpT, hApT);
  hLpTOverApT->SetNameTitle("hLpTOverApT","LpT/ApT");
  TGraphAsymmErrors *hHpTOverApT = new TGraphAsymmErrors(hHpT, hApT);
  hHpTOverApT->SetNameTitle("hHpTOverApT","HpT/ApT");
  TGraphAsymmErrors *hHpTOverLpT = new TGraphAsymmErrors(hHpT, hLpT);
  hHpTOverLpT->SetNameTitle("hHpTOverLpT","HpT/LpT");
  
  TCanvas *c1 = new TCanvas("cTrigger","cTrigger",10,10,1200,300);
  c1->Divide(4,1);
  c1->cd(1);
  gPad->SetLogy();
  hApT->DrawClone();
  c1->cd(2);
  gPad->SetLogy();
  hLpT->DrawClone();
  c1->cd(3);
  gPad->SetLogy();
  hHpT->DrawClone();
  c1->cd(4);
  gPad->SetLogy();
  hAll->DrawClone();
  
  TCanvas *c2 = new TCanvas("cTriggerRatio","cTriggerRatio",10,10,1200,300);
  c2->Divide(4,1);
  c2->cd(1);
  hLpTOverApT->DrawClone("ap");
  c2->cd(2);
  hHpTOverApT->DrawClone("ap");
  c2->cd(3);
  hHpTOverLpT->DrawClone("ap");
  c2->cd(4);
  hLpTOverAll->DrawClone("ap");
  /*
  // clean memory
  delete hLpTOverApT;
  delete hHpTOverApT;
  delete hHpTOverLpT;
  delete efficiency;
  file->Close();
  delete file;
  */
}

