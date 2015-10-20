//
//  CompareMuonDistributions.C
//  aliphysics-dev
//
//  Created by philippe pillot on 14/08/2015.
//  Copyright (c) 2015 Philippe Pillot. All rights reserved.
//

#include <TH1.h>
#include <TString.h>
#include <Riostream.h>
#include <TFile.h>
#include <TList.h>
#include <TCanvas.h>
#include <THashList.h>
#include <TParameter.h>

const Int_t nHist = 3;
TString sRes[nHist] = {"fHistPt", "fHistY", "fHistPhi"};
THashList *runWeights = 0x0;

void LoadRunWeights(TString fileName);
void AddHisto(TString sfile[2], TH1 *hRes[nHist][3], Double_t weight);

//______________________________________________________________________________
void CompareMuonDistributions(TString dir1, TString dir2, TString fileNameWeights = "")
{
  /// compare reconstructed muon distributions used for efficiency measurements
  
  // load run weights
  if (!fileNameWeights.IsNull()) {
    LoadRunWeights(fileNameWeights);
    if (!runWeights) return;
  }
  
  TH1 *hRes[nHist][3];
  for (Int_t i = 0; i < nHist; i++) for (Int_t j = 0; j < 3; j++) hRes[i][j] = 0x0;
  
  // get results
  if (runWeights) {
    TIter next(runWeights);
    TParameter<Double_t> *weight = 0x0;
    while ((weight = static_cast<TParameter<Double_t>*>(next()))) {
      TString sfile[2];
      sfile[0] = Form("%s/runs/%s/AnalysisResults.root", dir1.Data(), weight->GetName());
      sfile[1] = Form("%s/runs/%s/AnalysisResults.root", dir2.Data(), weight->GetName());
      AddHisto(sfile, hRes, weight->GetVal());
    }
  } else {
    TString sfile[2];
    sfile[0] = Form("%s/AnalysisResults.root", dir1.Data());
    sfile[1] = Form("%s/AnalysisResults.root", dir2.Data());
    AddHisto(sfile, hRes, 1.);
  }
  
  // compute ratios
  for (Int_t i = 0; i < nHist && hRes[i][0] && hRes[i][1]; i++) {
    hRes[i][2] = static_cast<TH1*>(hRes[i][1]->Clone(Form("%sOver1",hRes[i][1]->GetName())));
    hRes[i][2]->Divide(hRes[i][0]);
  }
  
  // display results at the reconstruction level
  TCanvas *cRec = new TCanvas("cRec", "cRec", 1200, 800);
  cRec->Divide(nHist,2);
  for (Int_t i = 0; i < nHist; i++) {
    cRec->cd(i+1);
    if (i == 0) gPad->SetLogy();
    for (Int_t j = 0; j < 2 && hRes[i][j]; j++) {
      hRes[i][j]->SetLineColor(2*j+2);
      hRes[i][j]->Draw((j == 0) ? "" : "sames");
    }
    cRec->cd(i+nHist+1);
    if (hRes[i][2]) {
      hRes[i][2]->Draw();
      hRes[i][2]->SetMinimum(0.5);
      hRes[i][2]->SetMaximum(1.5);
    }
  }
  
}

//______________________________________________________________________________
void AddHisto(TString sfile[2], TH1 *hRes[nHist][3], Double_t weight)
{
  /// get or add histograms with given weight
  
  // get results
  for (Int_t j = 0; j < 2; j++) {
    TFile *file = TFile::Open(sfile[j].Data(),"READ");
    if (!file || !file->IsOpen()) {
      printf("cannot open file\n");
      return;
    }
    if (file && file->IsOpen()) {
      TList *list = static_cast<TList*>(file->FindObjectAny("ExtraHistos"));
      if (!list) {
        printf("cannot find histograms\n");
        return;
      }
      for (Int_t i = 0; i < nHist; i++) {
        if (!hRes[i][j]) {
          hRes[i][j] = static_cast<TH1*>(list->FindObject(sRes[i].Data())->Clone(Form("%s%d",sRes[i].Data(),j+1)));
          if (hRes[i][j]) {
            hRes[i][j]->SetDirectory(0);
            hRes[i][j]->Sumw2();
            Double_t nEntries = static_cast<Double_t>(hRes[i][j]->GetEntries());
            if (nEntries > 0.) hRes[i][j]->Scale(weight/nEntries);
          }
        } else {
          TH1* h = static_cast<TH1*>(list->FindObject(sRes[i].Data()));
          Double_t nEntries = static_cast<Double_t>(h->GetEntries());
          hRes[i][j]->Add(h, weight/nEntries);
        }
      }
    }
    file->Close();
  }
  
}

//______________________________________________________________________________
void LoadRunWeights(TString fileName)
{
  /// Set the number of interested events per run
  /// (used to weight the acc*eff correction integrated
  /// over run for any pt/y/centrality bins)
  
  if (runWeights) return;
  
  ifstream inFile(fileName.Data());
  if (!inFile.is_open()) {
    printf("cannot open file %s\n", fileName.Data());
    return;
  }
  
  runWeights = new THashList(1000);
  runWeights->SetOwner();
  
  TString line;
  while (! inFile.eof() ) {
    
    line.ReadLine(inFile,kTRUE);
    if(line.IsNull()) continue;
    
    TObjArray *param = line.Tokenize(" ");
    if (param->GetEntries() != 2) {
      printf("bad input line %s", line.Data());
      continue;
    }
    
    Int_t run = ((TObjString*)param->UncheckedAt(0))->String().Atoi();
    if (run < 0) {
      printf("invalid run number: %d", run);
      continue;
    }
    
    Float_t weight = ((TObjString*)param->UncheckedAt(1))->String().Atof();
    if (weight <= 0.) {
      printf("invalid weight: %g", weight);
      continue;
    }
    
    runWeights->Add(new TParameter<Double_t>(Form("%d",run), weight));
    
    delete param;
  }
  
  inFile.close();
  
}


