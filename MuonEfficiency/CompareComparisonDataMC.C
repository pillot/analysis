//
//  CompareComparisonDataMC.C
//  aliphysics-dev
//
//  Created by philippe pillot on 10/05/2016.
//  Copyright (c) 2016 Philippe Pillot. All rights reserved.
//

#include <TAxis.h>
#include <TString.h>
#include <TObjString.h>
#include <Riostream.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include <TMath.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TPRegexp.h>

void CompareComparisonDataMC(TString input)
{
  /// input is a text file with the format: file label
  
  ifstream inFile(input.Data());
  if (!inFile.is_open()) {
    printf("cannot open file %s\n", input.Data());
    return;
  }
  
  TCanvas *cRatioCent = new TCanvas("cRatioCent","cRatioCent",1200,400);
  cRatioCent->SetGridy();
  TCanvas *cRatioRun = new TCanvas("cRatioRun","cRatioRun",1200,400);
  cRatioRun->SetGridy();
  TCanvas *cRatioY = new TCanvas("cRatioY","cRatioY",1200,400);
  cRatioY->SetGridy();
  TCanvas *cRatioPt = new TCanvas("cRatioPt","cRatioPt",1200,400);
  cRatioPt->SetGridy();
  TCanvas *cRatioPhi = new TCanvas("cRatioPhi","cRatioPhi",1200,400);
  cRatioPhi->SetGridy();
  
  TLegend *legend = new TLegend(0.8, 0.8, 0.95, 0.95);
  legend->SetTextSize(0.06);
  legend->SetMargin(0.1);
  
  Int_t color = 0;
  
  TString line;
  TPRegexp re("( |\t).*");
  while (! inFile.eof() ) {
    
    line.ReadLine(inFile,kTRUE);
    if (line.IsNull() || line.BeginsWith("#")) continue;
    
    TString label = line(re);
    TString fileName = line.ReplaceAll(label.Data(),"");
    label.Remove(TString::kLeading,' ');
    label.Remove(TString::kLeading,'\t');
    if (fileName.IsNull() || label.IsNull()) {
      printf("incorrect list of files and/or legends\n");
      return;
    }
    
    TFile *file = new TFile(fileName.Data(), "read");
    if (!file || !file->IsOpen()) {
      printf("cannot open file %s \n",fileName.Data());
      return;
    }
    
    TObjArray *globalRatios = static_cast<TObjArray*>(file->FindObjectAny("GlobalEffRatios"));
    if (!globalRatios) {
      printf("GlobalEffRatios not found\n");
      return;
    }
    
    TGraphAsymmErrors *ratioCent = static_cast<TGraphAsymmErrors*>(globalRatios->FindObject("RatioEffVsCent"));
    if (!ratioCent) {
      printf("RatioEffVsCent not found\n");
      return;
    }
    
    TGraphAsymmErrors *ratioRun = static_cast<TGraphAsymmErrors*>(globalRatios->FindObject("RatioEffVsRun"));
    if (!ratioRun) {
      printf("RatioEffVsRun not found\n");
      return;
    }
    
    TGraphAsymmErrors *ratioY = static_cast<TGraphAsymmErrors*>(globalRatios->FindObject("RatioEffVsY"));
    if (!ratioY) {
      printf("RatioEffVsY not found\n");
      return;
    }
    
    TGraphAsymmErrors *ratioPt = static_cast<TGraphAsymmErrors*>(globalRatios->FindObject("RatioEffVsPt"));
    if (!ratioPt) {
      printf("RatioEffVsPt not found\n");
      return;
    }
    
    TGraphAsymmErrors *ratioPhi = static_cast<TGraphAsymmErrors*>(globalRatios->FindObject("RatioEffVsPhi"));
    if (!ratioPhi) {
      printf("RatioEffVsPhi not found\n");
      return;
    }
    
    ++color;
    if (color == 5 || color == 10) ++color;
    
    cRatioCent->cd();
    ratioCent->SetMarkerColor(color);
    ratioCent->SetLineColor(color);
    if (color == 1) {
      ratioCent->SetMinimum(0.95);
      ratioCent->SetMaximum(1.05);
      ratioCent->GetXaxis()->SetLabelSize(0.06);
      ratioCent->GetYaxis()->SetLabelSize(0.06);
      ratioCent->Draw("ap");
    } else ratioCent->Draw("p");
    
    cRatioRun->cd();
    ratioRun->SetMarkerColor(color);
    ratioRun->SetLineColor(color);
    if (color == 1) {
      ratioRun->SetMinimum(0.9);
      ratioRun->SetMaximum(1.1);
      ratioRun->GetXaxis()->SetLabelSize(0.06);
      ratioRun->GetYaxis()->SetLabelSize(0.06);
      ratioRun->Draw("ap");
    } else ratioRun->Draw("p");
    
    cRatioY->cd();
    ratioY->SetMarkerColor(color);
    ratioY->SetLineColor(color);
    if (color == 1) {
      ratioY->SetMinimum(0.95);
      ratioY->SetMaximum(1.05);
      ratioY->GetXaxis()->SetLabelSize(0.06);
      ratioY->GetYaxis()->SetLabelSize(0.06);
      ratioY->Draw("ap");
    } else ratioY->Draw("p");
    
    cRatioPt->cd();
    ratioPt->SetMarkerColor(color);
    ratioPt->SetLineColor(color);
    if (color == 1) {
      ratioPt->SetMinimum(0.95);
      ratioPt->SetMaximum(1.05);
      ratioPt->GetXaxis()->SetLabelSize(0.06);
      ratioPt->GetYaxis()->SetLabelSize(0.06);
      ratioPt->Draw("ap");
    } else ratioPt->Draw("p");
    
    cRatioPhi->cd();
    ratioPhi->SetMarkerColor(color);
    ratioPhi->SetLineColor(color);
    if (color == 1) {
      ratioPhi->SetMinimum(0.92);
      ratioPhi->SetMaximum(1.08);
      ratioPhi->GetXaxis()->SetLabelSize(0.06);
      ratioPhi->GetYaxis()->SetLabelSize(0.06);
      ratioPhi->Draw("ap");
    } else ratioPhi->Draw("p");
    
    legend->AddEntry(ratioPhi, label.Data(), "ep");
    
  }
  
  inFile.close();
  
  cRatioPhi->cd();
  legend->Draw("same");
  
}
