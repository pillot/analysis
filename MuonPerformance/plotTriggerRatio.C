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

AliCFEffGrid* GetCFEff(TString inputFile);
void GetTriggerRatio(AliCFEffGrid* efficiency, TString var, TGraphAsymmErrors *hIpTOverAll[3],
                     TGraphAsymmErrors *hIpTOverJpT[3]);
void DrawTriggerRatio(TGraphAsymmErrors *hIpTOverAll[3], TGraphAsymmErrors *hIpTOverJpT[3]);
void CompareTriggerRatio(TGraphAsymmErrors *hIpTOverX1[3], TGraphAsymmErrors *hIpTOverX2[3]);

//-------------------------------------------------------------------------------------
void plotTriggerRatio(TString inputFile = "AnalysisResults.root", TString inputFile2 = "")
{
  /// plot the ratios LpT/ApT, HpT/LpT, ...
  /*
   .include $ALICE_ROOT/include
   .include $ALICE_PHYSICS/include
   .x $WORK/Macros/MuonPerformance/plotTriggerRatio.C+
  */
  
  gROOT->SetStyle("Plain");
  gROOT->ForceStyle();
  
  // prepare environment
  if (!gInterpreter->IsLoaded("$WORK/Macros/Facilities/runTaskFacilities.C")) {
    gROOT->LoadMacro("$WORK/Macros/Facilities/runTaskFacilities.C");
    gROOT->ProcessLineFast("LoadAlirootLocally(\"\",\"include\",\"\",\"\");");
  }
  
  AliCFEffGrid *efficiency = GetCFEff(inputFile);
  if (!efficiency) return;
  
  if (inputFile2.IsNull()) {
    
    TGraphAsymmErrors *hIpTOverAll[3];
    TGraphAsymmErrors *hIpTOverJpT[3];
    
    GetTriggerRatio(efficiency, "pT", hIpTOverAll, hIpTOverJpT);
    DrawTriggerRatio(hIpTOverAll, hIpTOverJpT);
    
    GetTriggerRatio(efficiency, "eta", hIpTOverAll, hIpTOverJpT);
    DrawTriggerRatio(hIpTOverAll, hIpTOverJpT);
    
    GetTriggerRatio(efficiency, "phi", hIpTOverAll, hIpTOverJpT);
    DrawTriggerRatio(hIpTOverAll, hIpTOverJpT);
    
  } else {
    
    AliCFEffGrid *efficiency2 = GetCFEff(inputFile2);
    if (!efficiency2) return;
    
    //efficiency->GetNum()->SetRangeUser(4, 1., 1.); // charge
    //efficiency2->GetNum()->SetRangeUser(4, 1., 1.); // charge
    
    TGraphAsymmErrors *hIpTOverAll[3];
    TGraphAsymmErrors *hIpTOverJpT[3];
    TGraphAsymmErrors *hIpTOverAll2[3];
    TGraphAsymmErrors *hIpTOverJpT2[3];
    
    GetTriggerRatio(efficiency, "pT", hIpTOverAll, hIpTOverJpT);
    GetTriggerRatio(efficiency2, "pT", hIpTOverAll2, hIpTOverJpT2);
    CompareTriggerRatio(hIpTOverAll, hIpTOverAll2);
    CompareTriggerRatio(hIpTOverJpT, hIpTOverJpT2);
    
    GetTriggerRatio(efficiency, "eta", hIpTOverAll, hIpTOverJpT);
    GetTriggerRatio(efficiency2, "eta", hIpTOverAll2, hIpTOverJpT2);
    CompareTriggerRatio(hIpTOverAll, hIpTOverAll2);
    CompareTriggerRatio(hIpTOverJpT, hIpTOverJpT2);
    
    GetTriggerRatio(efficiency, "phi", hIpTOverAll, hIpTOverJpT);
    GetTriggerRatio(efficiency2, "phi", hIpTOverAll2, hIpTOverJpT2);
    CompareTriggerRatio(hIpTOverAll, hIpTOverAll2);
    CompareTriggerRatio(hIpTOverJpT, hIpTOverJpT2);
    
  }
  
}

//-------------------------------------------------------------------------------------
AliCFEffGrid* GetCFEff(TString inputFile)
{
  /// get the efficiency object and define some ranges for projection
  
  TFile *file = new TFile(inputFile.Data(), "read");
  if (!file || !file->IsOpen()) {
    printf("cannot open file %s\n", inputFile.Data());
    return 0x0;
  }
  
  AliCFContainer* cfContainer = static_cast<AliCFContainer*>(file->FindObjectAny("EffContainer"));
  if (!cfContainer) return 0x0;
  AliCFEffGrid* efficiency = new AliCFEffGrid("eff","",*cfContainer);
  if (!efficiency) return 0x0;
  efficiency->CalculateEfficiency(0, 1);
  
  // track selection
  //efficiency->GetNum()->SetRangeUser(0, 2., 15.); // pT
  efficiency->GetNum()->SetRangeUser(5, 1., 1.); // tracker
  efficiency->GetNum()->SetRangeUser(8, 1., 3.); // match MC
  //  efficiency->GetNum()->SetRangeUser(9, 1., 4.); // MC trig level
  efficiency->GetNum()->SetRangeUser(1, -3.999, -2.501); // eta
  efficiency->GetNum()->SetRangeUser(3, 1., 2.); // theta_abs
  //  efficiency->GetNum()->SetRangeUser(10, 0., 10.); // centrality
  
  file->Close();
  
  return efficiency;
  
}

//-------------------------------------------------------------------------------------
void GetTriggerRatio(AliCFEffGrid* efficiency, TString var,
                     TGraphAsymmErrors *hIpTOverAll[3],
                     TGraphAsymmErrors *hIpTOverJpT[3])
{
  /// get all the ratios
  
  Int_t iVar;
  if (var == "pT") iVar = 0;
  else if (var == "eta") iVar = 1;
  else if (var == "phi") iVar = 2;
  else return;
  
  efficiency->GetNum()->GetAxis(6)->SetRange();
  TH1 *hAll = efficiency->GetNum()->Project(iVar);
  
  efficiency->GetNum()->SetRangeUser(6, 2., 4.);
  TH1 *hApT = efficiency->GetNum()->Project(iVar);
  
  efficiency->GetNum()->SetRangeUser(6, 3., 4.);
  TH1 *hLpT = efficiency->GetNum()->Project(iVar);
  
  efficiency->GetNum()->SetRangeUser(6, 4., 4.);
  TH1 *hHpT = efficiency->GetNum()->Project(iVar);
  
  hIpTOverAll[0] = new TGraphAsymmErrors(hApT, hAll);
  hIpTOverAll[0]->SetNameTitle(Form("hApTOverAllvs%s",var.Data()),Form("ApT/All vs %s",var.Data()));
  
  hIpTOverAll[1] = new TGraphAsymmErrors(hLpT, hAll);
  hIpTOverAll[1]->SetNameTitle(Form("hLpTOverAllvs%s",var.Data()),Form("LpT/All vs %s",var.Data()));
  
  hIpTOverAll[2] = new TGraphAsymmErrors(hHpT, hAll);
  hIpTOverAll[2]->SetNameTitle(Form("hHpTOverAllvs%s",var.Data()),Form("HpT/All vs %s",var.Data()));
  
  hIpTOverJpT[0] = new TGraphAsymmErrors(hLpT, hApT);
  hIpTOverJpT[0]->SetNameTitle(Form("hLpTOverApTvs%s",var.Data()),Form("LpT/ApT vs %s",var.Data()));
  
  hIpTOverJpT[1] = new TGraphAsymmErrors(hHpT, hApT);
  hIpTOverJpT[1]->SetNameTitle(Form("hHpTOverApTvs%s",var.Data()),Form("HpT/ApT vs %s",var.Data()));
  
  hIpTOverJpT[2] = new TGraphAsymmErrors(hHpT, hLpT);
  hIpTOverJpT[2]->SetNameTitle(Form("hHpTOverLpTvs%s",var.Data()),Form("HpT/LpT vs %s",var.Data()));
  
}

//-------------------------------------------------------------------------------------
void DrawTriggerRatio(TGraphAsymmErrors *hIpTOverAll[3], TGraphAsymmErrors *hIpTOverJpT[3])
{
  /// get all the ratios
  
  TString var;
  TString hName = hIpTOverAll[0]->GetName();
  if (hName.Contains("vspT")) var = "pT";
  else if (hName.Contains("vseta")) var = "eta";
  else if (hName.Contains("vsphi")) var = "phi";
  else return;
  
  TCanvas *c = new TCanvas(Form("cTriggerRatiovs%s",var.Data()),Form("Trigger ratio vs %s",var.Data()),10,10,900,600);
  c->Divide(3,2);
  for (Int_t i = 0; i < 3; ++i) {
    c->cd(i+1);
    hIpTOverAll[i]->Draw("ap");
    c->cd(i+4);
    hIpTOverJpT[i]->Draw("ap");
  }
  
}

//-------------------------------------------------------------------------------------
void CompareTriggerRatio(TGraphAsymmErrors *hIpTOverX1[3], TGraphAsymmErrors *hIpTOverX2[3])
{
  /// get all the ratios
  
  TString var;
  TString hName = hIpTOverX1[0]->GetName();
  if (hName.Contains("vspT")) var = "pT";
  else if (hName.Contains("vseta")) var = "eta";
  else if (hName.Contains("vsphi")) var = "phi";
  else return;
  
  TString den = "XpT";
  if (hName.Contains("OverAll")) den = "All";
  
  TGraphAsymmErrors *hRatio[3];
  for (Int_t i = 0; i < 3; ++i) {
    
    Int_t n = hIpTOverX1[i]->GetN();
    hRatio[i] = new TGraphAsymmErrors(n);
    hRatio[i]->SetNameTitle(Form("%sRatio",hIpTOverX1[i]->GetName()),Form("ratio of %s",hIpTOverX1[i]->GetTitle()));
    
    const Double_t *x = hIpTOverX1[i]->GetX();
    const Double_t *y1 = hIpTOverX1[i]->GetY();
    const Double_t *eyl1 = hIpTOverX1[i]->GetEYlow();
    const Double_t *eyh1 = hIpTOverX1[i]->GetEYhigh();
    const Double_t *y2 = hIpTOverX2[i]->GetY();
    const Double_t *eyl2 = hIpTOverX2[i]->GetEYlow();
    const Double_t *eyh2 = hIpTOverX2[i]->GetEYhigh();
    
    for (Int_t j = 0; j < n; ++j) {
      
      Double_t y = (y1[j] > 0.) ? y2[j]/y1[j] : 0.;
      hRatio[i]->SetPoint(j, x[j], y);
      
      Double_t eyl = 0., eyh = 0.;
      if (y1[j] > 0. && y2[j] > 0.) {
        eyl = y * TMath::Sqrt(eyl2[j]*eyl2[j]/y2[j]/y2[j] + eyh1[j]*eyh1[j]/y1[j]/y1[j]);
        eyh = y * TMath::Sqrt(eyh2[j]*eyh2[j]/y2[j]/y2[j] + eyl1[j]*eyl1[j]/y1[j]/y1[j]);
      }
      hRatio[i]->SetPointError(j, 0., 0., eyl, eyh);
      
    }
    
  }
  
  TCanvas *c = new TCanvas(Form("cTriggerRatioOver%svs%s",den.Data(),var.Data()),Form("Trigger ratio over %s vs %s",den.Data(),var.Data()),10,10,900,600);
  c->Divide(3,2);
  for (Int_t i = 0; i < 3; ++i) {
    c->cd(i+1);
    hIpTOverX1[i]->Draw("ap");
    hIpTOverX2[i]->SetLineColor(2);
    hIpTOverX2[i]->Draw("p");
    c->cd(i+4);
    hRatio[i]->Draw("ap");
  }
  
}

