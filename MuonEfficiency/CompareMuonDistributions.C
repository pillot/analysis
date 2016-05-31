//
//  CompareMuonDistributions.C
//  aliphysics-dev
//
//  Created by philippe pillot on 14/08/2015.
//  Copyright (c) 2015 Philippe Pillot. All rights reserved.
//

#include <TH1.h>
#include <THnSparse.h>
#include <TAxis.h>
#include <TString.h>
#include <Riostream.h>
#include <TFile.h>
#include <TList.h>
#include <TCanvas.h>
#include <THashList.h>
#include <TParameter.h>
#include <TLegend.h>
#include "AliCounterCollection.h"

// kinematics range for histogram projection (pT, y, phi, charge)
const Float_t kineRange[4][2] = {{-999., 999.}, {-999., 999.}, {-999., 999.}, {-999., 999.}};

// centrality range for each input file
const Float_t centRange[2][2] = {{-999., 999.}, {-999., 999.}};

// in case the ouput containers of the efficiency task have an extension in their name (like in the train)
TString extension[2] = {"",""};

// normalize to the total number of events instead of to the number of entries in the histograms
Bool_t eventNorm = kFALSE;

const Int_t nHist = 3;
TString sRes[nHist] = {"fHistPt", "fHistY", "fHistPhi"};
THashList *runWeights = 0x0;

void LoadRunWeights(TString fileName);
void AddHisto(TString sfile[2], TH1 *hRes[nHist][3], Double_t weight);
void AddHistoProj(TString sfile[2], TH1 *hProj[4][3], Double_t weight);
void SetKineRange(THnSparse& hKine);
void SetCentRange(THnSparse& hKine, const Float_t centRange[2]);

//______________________________________________________________________________
void CompareMuonDistributions(TString dir1, TString dir2, TString fileNameWeights = "")
{
  /// compare reconstructed muon distributions used for efficiency measurements
  /*
   .include $ALICE_ROOT/include
   .x $WORK/Macros/MuonEfficiency/CompareMuonDistributions.C+(...)
  */
  
  // load run weights
  if (!fileNameWeights.IsNull()) {
    LoadRunWeights(fileNameWeights);
    if (!runWeights) return;
  }
  
  TH1 *hRes[nHist][3];
  for (Int_t i = 0; i < nHist; i++) for (Int_t j = 0; j < 3; j++) hRes[i][j] = 0x0;
  TH1 *hProj[4][3];
  for (Int_t i = 0; i < 4; i++) for (Int_t j = 0; j < 3; j++) hProj[i][j] = 0x0;
  
  // get results
  if (runWeights) {
    TIter next(runWeights);
    TParameter<Double_t> *weight = 0x0;
    while ((weight = static_cast<TParameter<Double_t>*>(next()))) {
      TString sfile[2];
      sfile[0] = Form("%s/runs/%s/AnalysisResults.root", dir1.Data(), weight->GetName());
//      sfile[0] = Form("%s/runs/%s/QAresults.root", dir1.Data(), weight->GetName());
      sfile[1] = Form("%s/runs/%s/AnalysisResults.root", dir2.Data(), weight->GetName());
      AddHisto(sfile, hRes, weight->GetVal());
      AddHistoProj(sfile, hProj, weight->GetVal());
    }
  } else {
    TString sfile[2];
    sfile[0] = Form("%s/AnalysisResults.root", dir1.Data());
    sfile[1] = Form("%s/AnalysisResults.root", dir2.Data());
//    sfile[1] = Form("%s/QAresults.root", dir2.Data());
    AddHisto(sfile, hRes, 1.);
    AddHistoProj(sfile, hProj, 1.);
  }
  
  // compute ratios
  for (Int_t i = 0; i < nHist && hRes[i][0] && hRes[i][1]; i++) {
    hRes[i][2] = static_cast<TH1*>(hRes[i][1]->Clone(Form("%sOver1",hRes[i][1]->GetName())));
    hRes[i][2]->Divide(hRes[i][0]);
  }
  for (Int_t i = 0; i < 4 && hProj[i][0] && hProj[i][1]; i++) {
    hProj[i][2] = static_cast<TH1*>(hProj[i][1]->Clone(Form("%sOver1",hProj[i][1]->GetName())));
    hProj[i][2]->Divide(hProj[i][0]);
  }
  
  // display results at the reconstruction level
  TLegend *lRes = new TLegend(0.5,0.55,0.85,0.75);
  TCanvas *cRec = new TCanvas("cRec", "cRec", 1200, 800);
  cRec->Divide(nHist,2);
  for (Int_t i = 0; i < nHist; i++) {
    cRec->cd(i+1);
    if (i == 0) gPad->SetLogy();
    for (Int_t j = 0; j < 2 && hRes[i][j]; j++) {
      if (i == 0) lRes->AddEntry(hRes[i][j],Form("dir%d",j+1),"l");
      hRes[i][j]->SetLineColor(2*j+2);
      hRes[i][j]->Draw((j == 0) ? "" : "sames");
    }
    if (i == 0) lRes->Draw("same");
    cRec->cd(i+nHist+1);
    if (hRes[i][2]) {
      hRes[i][2]->Draw();
      hRes[i][2]->SetMinimum(0.5);
      hRes[i][2]->SetMaximum(1.5);
    }
  }
  TCanvas *cRecProj = new TCanvas("cRecProj", "cRecProj", 1600, 800);
  cRecProj->Divide(4,2);
  for (Int_t i = 0; i < 4; i++) {
    cRecProj->cd(i+1);
    if (i == 0) gPad->SetLogy();
    for (Int_t j = 0; j < 2 && hProj[i][j]; j++) {
      hProj[i][j]->SetLineColor(2*j+2);
      hProj[i][j]->Draw((j == 0) ? "" : "sames");
    }
    cRecProj->cd(i+5);
    if (hProj[i][2]) {
      hProj[i][2]->Draw();
      hProj[i][2]->SetMinimum(0.5);
      hProj[i][2]->SetMaximum(1.5);
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
      TList *list = static_cast<TList*>(file->FindObjectAny(Form("ExtraHistos%s",extension[j].Data())));
      if (!list) {
        printf("cannot find histograms\n");
        return;
      }
      AliCounterCollection *eventCounters = static_cast<AliCounterCollection*>(file->FindObjectAny(Form("EventsCounters%s",extension[j].Data())));
      if (!eventCounters) {
        printf("cannot find event counters\n");
        return;
      }
      Double_t norm = eventCounters->GetSum();
//      if (j == 1) norm *= 5./2.;
      for (Int_t i = 0; i < nHist; i++) {
        if (!hRes[i][j]) {
          hRes[i][j] = static_cast<TH1*>(list->FindObject(sRes[i].Data())->Clone(Form("%s%d",sRes[i].Data(),j+1)));
          if (hRes[i][j]) {
            hRes[i][j]->SetDirectory(0);
            hRes[i][j]->Sumw2();
            if (!eventNorm) norm = static_cast<Double_t>(hRes[i][j]->GetEntries());
            if (norm > 0.) hRes[i][j]->Scale(weight/norm);
          }
        } else {
          TH1* h = static_cast<TH1*>(list->FindObject(sRes[i].Data()));
          if (!eventNorm) norm = static_cast<Double_t>(h->GetEntries());
          if (norm > 0.) hRes[i][j]->Add(h, weight/norm);
        }
      }
    }
    file->Close();
  }
  
}
//______________________________________________________________________________
void AddHistoProj(TString sfile[2], TH1 *hProj[4][3], Double_t weight)
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
      TList *list = static_cast<TList*>(file->FindObjectAny(Form("ExtraHistos%s",extension[j].Data())));
      if (!list) {
        printf("cannot find histograms\n");
        return;
      }
      THnSparse *hKine = static_cast<THnSparse*>(list->FindObject("hKine"));
      if (!hKine) return;
/*      if (j == 0) {
        kineRange[3][0] = 1.;
        kineRange[3][1] = 1.;
      } else {
        kineRange[3][0] = -1.;
        kineRange[3][1] = -1.;
      }
*/      SetKineRange(*hKine);
      SetCentRange(*hKine, centRange[j]);
      AliCounterCollection *eventCounters = static_cast<AliCounterCollection*>(file->FindObjectAny(Form("EventsCounters%s",extension[j].Data())));
      if (!eventCounters) {
        printf("cannot find event counters\n");
        return;
      }
      Double_t norm = eventCounters->GetSum();
//      if (j == 1) norm *= 5./2.;
      for (Int_t i = 0; i < 4; i++) {
        if (!hProj[i][j]) {
          hProj[i][j] = hKine->Projection(i+1,"eo");
          if (hProj[i][j]) {
            hProj[i][j]->SetDirectory(0);
            if (!eventNorm) norm = static_cast<Double_t>(hProj[i][j]->GetEntries());
            if (norm > 0.) hProj[i][j]->Scale(weight/norm);
          }
        } else {
          TH1* h = hKine->Projection(i+1,"eo");
          if (!eventNorm) norm = static_cast<Double_t>(h->GetEntries());
          if (norm > 0.) hProj[i][j]->Add(h, weight/norm);
          delete h;
        }
      }
    }
    file->Close();
  }
  
}

//______________________________________________________________________________
void SetKineRange(THnSparse& hKine)
{
  /// Sets the kinematics range for histogram projection
  /// If the range exceeds the minimum (maximum) it includes the underflow (overflow)
  /// in the integration (same behaviour as if the range is not set)
  
  for (Int_t i = 0; i < 4; i++) {
    TAxis *a = hKine.GetAxis(i+1);
    Int_t lowBin = a->FindBin(kineRange[i][0]);
    Int_t upBin = a->FindBin(kineRange[i][1]);
    a->SetRange(lowBin, upBin);
  }
  
}

//______________________________________________________________________________
void SetCentRange(THnSparse& hKine, const Float_t centRange[2])
{
  /// Sets the centrality range for histogram projection
  /// If the range exceeds the minimum (maximum) it includes the underflow (overflow)
  /// in the integration (same behaviour as if the range is not set)
  
  TAxis *a = hKine.GetAxis(0);
  Int_t lowBin = a->FindBin(centRange[0]);
  Int_t upBin = a->FindBin(centRange[1]);
  a->SetRange(lowBin, upBin);
  
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


