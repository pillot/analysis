//
//  CompareDiff.C
//  aliroot_dev
//
//  Created by philippe pillot on 13/07/2014.
//  Copyright (c) 2014 Philippe Pillot. All rights reserved.
//

void CompareDiff()
{
  /// compare diff histos between files
  
  gROOT->SetStyle("Plain");
  gROOT->ForceStyle();
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  /*
  const Int_t nFiles = 3;
  TString fileName[nFiles] = {
    "/Users/pillot/Work/Alice/Data/2011/LHC11h/pass2_muon/JPsi/Diff_MULorMLL_MULorMLL_TrgSign_2.8-3.3.root",
    "/Users/pillot/Work/Alice/Data/2011/LHC11h/pass2_muon/JPsipDCA/Diff_MULorMLL_MULorMLL_TrgSign_2.8-3.3.root",
    "/Users/pillot/Work/Alice/Sim/LHC11h/EmbeddingPass2/JPsi/Diff_MB_MB_TrgSign.root"
  };
  TString legend[nFiles] = {
    "w/o pDCA",
    "w/   pDCA",
    "JPsi simu"
  };
  */
  /*
  const Int_t nFiles = 3;
  TString fileName[nFiles] = {
    "/Users/pillot/Work/Alice/Data/2011/LHC11h/pass2_muon/JPsipDCA/Diff_MULorMLL_MULorMLL_TrgSign_2.8-3.3.root",
    "/Users/pillot/Work/Alice/Data/2011/LHC11h/pass2_muon/JPsipDCA/DiffJPsi_MULorMLL_MULorMLL_TrgSign.root",
    "/Users/pillot/Work/Alice/Sim/LHC11h/EmbeddingPass2/JPsi/Diff_MB_MB_TrgSign.root"
  };
  TString legend[nFiles] = {
    "2.8-3.3",
    "JPsi data",
    "JPsi simu"
  };
  */
  /*
  const Int_t nFiles = 2;
  TString fileName[4] = {
    "/Users/pillot/Work/Alice/Data/2011/LHC11h/pass2_muon/JPsipDCA/Diff_MULorMLL_MUL_2.8-3.3.root",
    "/Users/pillot/Work/Alice/Data/2011/LHC11h/pass2_muon/JPsipDCA/Diff_MULorMLL_MUL_TrgSign_2.8-3.3.root",
    "/Users/pillot/Work/Alice/Data/2011/LHC11h/pass2_muon/JPsipDCA/DiffJPsi_MULorMLL_MULorMLL_TrgSign.root",
    "/Users/pillot/Work/Alice/Sim/LHC11h/EmbeddingPass2/JPsi/Diff_MB_MB_TrgSign.root"
  };
  TString legend[4] = {
    "MUL",
    "MUL+Trg",
    "JPsi",
    "simu"
  };
  */
  /*
  const Int_t nFiles = 2;
  TString fileName[nFiles] = {
    "/Users/pillot/Work/Alice/Data/2011/LHC11h/pass2_muon/JPsipDCA/Diff_MULorMLL_MULorMLL_TrgSign_9-10.root",
    "/Users/pillot/Work/Alice/Data/2011/LHC11h/pass2_muon/JPsipDCA/Diff_MULorMLL_MULorMLL_TrgSign_2.8-3.3.root"
  };
  TString legend[nFiles] = {
    "9-10",
    "2.8-3.3"
  };
  */
  /*
  const Int_t nFiles = 2;
  TString fileName[nFiles] = {
    "/Users/pillot/Work/Alice/Data/2011/LHC11h/pass2_muon/JPsipDCA/Diff_MULorMLL_MULorMLL_TrgSign_2.8-3.3.root",
    "/Users/pillot/Work/Alice/Data/2011/LHC11h/pass2_muon/JPsipDCA/Diff_MULorMLL_MULorMLL_TrgFake_2.8-3.3.root"
  };
  TString legend[nFiles] = {
    "MUL+Trg",
    "Fake Trg"
  };
  */
  /*
  const Int_t nFiles = 3;
  TString fileName[nFiles] = {
    "/Users/pillot/Work/Alice/Data/2011/LHC11h/pass2_muon/JPsipDCA/Diff_MULorMLL_MULorMLL_TrgSign_2.8-3.3.root",
    "/Users/pillot/Work/Alice/Data/2011/LHC11h/pass2_muon/JPsipDCA/Diff_MULorMLL_MULorMLL_TrgSign_2.8-3.3_PP.root",
    "/Users/pillot/Work/Alice/Data/2011/LHC11h/pass2_muon/JPsipDCA/Diff_MULorMLL_MULorMLL_TrgSign_2.8-3.3_MM.root"
  };
  TString legend[nFiles] = {
    "PM",
    "PP",
    "MM"
  };
  */
  /*
  const Int_t nFiles = 2;
  TString fileName[nFiles] = {
    "/Users/pillot/Work/Alice/Data/2011/LHC11h/pass2_muon/JPsi/Diff_MULorMLL_MULorMLL_TrgFake_2.8-3.3.root",
    "/Users/pillot/Work/Alice/Data/2011/LHC11h/pass2_muon/JPsipDCA/Diff_MULorMLL_MULorMLL_TrgFake_2.8-3.3.root"
  };
  TString legend[nFiles] = {
    "w/o pDCA",
    "w/   pDCA"
  };
  */
  const Int_t nFiles = 4;
  TString fileName[nFiles] = {
    "/Users/pillot/Work/Alice/Sim/LHC11h/EmbeddingPass2/JPsi/Diff_MB_MB_TrgSign.root",
    "/Users/pillot/Work/Alice/Sim/LHC11h/EmbeddingPass2/JPsi/Diff_MB_relabeled_MB_TrgSign_relabeled.root",
    "/Users/pillot/Work/Alice/Sim/LHC11h/EmbeddingPass2/JPsi/Diff_MB_relabeled_10sigma_MB_TrgSign_relabeled_10sigma.root",
    "/Users/pillot/Work/Alice/Sim/LHC11h/EmbeddingPass2/JPsi/Diff_MB_relabeled_25sigma_MB_TrgSign_relabeled_25sigma.root"
  };
  TString legend[nFiles] = {
    "label reco",
    "relabeled",
    "relabeled10",
    "relabeled25"
  };
  
  Int_t lineColor[5] = {1, 2, 4, 8, 6};
  
  TString histoName[16] = {
    "hDiffVsCent",
    "hDiffVsCentpT1",
    "hDiffVsCentpT2",
    "hDiffVsCentpT3",
    "hDiffVspTCent0",
    "hDiffVspTCent1",
    "hDiffVspTCent2",
    "hDiffVspTCent3",
    "hDiffVspTCent4",
    "hDiffVsCenty1",
    "hDiffVsCenty2",
    "hDiffVsCenty3",
    "hDiffVsyCent0",
    "hDiffVsyCent1",
    "hDiffVsyCent2",
    "hDiffVsyCent4"
  };
  
  TLegend *l = new TLegend(0.7,0.7,0.9,0.9);
  l->SetFillStyle(0);
  l->SetBorderSize(0);
  
  TH1F *h[nFiles][16];
  for (Int_t i = 0; i < nFiles; ++i) {
    TFile *file = TFile::Open(fileName[i].Data(),"READ");
    for (Int_t j = 0; j < 16; ++j) {
      h[i][j] = static_cast<TH1F*>(file->FindObjectAny(histoName[j].Data()));
    }
    l->AddEntry(h[i][0],legend[i].Data(),"l");
  }
  
  for (Int_t j = 0; j < 16; ++j) {
    TCanvas* c = new TCanvas(Form("c%s",histoName[j].Data()), h[0][j]->GetTitle(), 800, 400);
    for (Int_t i = 0; i < nFiles; ++i) {
      h[i][j]->SetLineColor(lineColor[i]);
      h[i][j]->GetXaxis()->SetLabelSize(0.08);
      h[i][j]->GetYaxis()->SetLabelSize(0.06);
      h[i][j]->GetYaxis()->SetRangeUser(-3.5,0.);
      h[i][j]->Draw((i==0)?"":"same");
    }
    l->DrawClone("same");
  }
  
}

