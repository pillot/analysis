#include "TFile.h"
#include "THnSparse.h"
#include "TH1D.h"
#include "TH2.h"
#include "TProfile.h"
#include "TCanvas.h"

void CompareMeanTracklets(TString fileNameData1, TString fileNameData2, Bool_t corr1 = kFALSE, Bool_t corr2 = kFALSE)
{
  
  TString fileName[2] = {fileNameData1.Data(), fileNameData2.Data()};
  Bool_t corr[2] = {corr1, corr2};
  THnSparse *fhNtrk[2];
  for (Int_t i = 0; i < 2; ++i) {
    TFile *file = new TFile(fileName[i].Data(), "read");
    if (!file || !file->IsOpen()) {
      printf("cannot open file %s \n",fileName[i].Data());
      return;
    }
    fhNtrk[i] = static_cast<THnSparse*>(file->FindObjectAny(corr[i] ? "hNtrkCorr" : "hNtrk"));
    if (!fhNtrk[i]) return;
  }
  
//  fhNtrk[0]->GetAxis(2)->SetRangeUser(244484, 244484);
//  fhNtrk[1]->GetAxis(2)->SetRangeUser(244542, 244542);
//  fhNtrk[0]->GetAxis(2)->SetRangeUser(244340, 244359);
//  fhNtrk[1]->GetAxis(2)->SetRangeUser(244617, 244628);
//  fhNtrk[0]->GetAxis(2)->SetRangeUser(244411, 244418);
//  fhNtrk[1]->GetAxis(2)->SetRangeUser(244411, 244418);
//  fhNtrk[0]->GetAxis(2)->SetRangeUser(244340, 244421);
//  fhNtrk[1]->GetAxis(2)->SetRangeUser(244340, 244421);
//  fhNtrk[0]->GetAxis(2)->SetRangeUser(244453, 244628);
//  fhNtrk[1]->GetAxis(2)->SetRangeUser(244453, 244628);
//  fhNtrk[0]->GetAxis(2)->SetRangeUser(266208, 266208);
//  fhNtrk[1]->GetAxis(2)->SetRangeUser(266208, 266208);
//  fhNtrk[0]->GetAxis(2)->SetRangeUser(266234, 266234);
//  fhNtrk[1]->GetAxis(2)->SetRangeUser(266234, 266234);
  
  TH1D *hNtrk[2];
  for (Int_t i = 0; i < 2; ++i) {
    hNtrk[i] = fhNtrk[i]->Projection(0,"e");
    hNtrk[i]->Scale(1./hNtrk[i]->GetEntries());
    hNtrk[i]->SetTitle("N tracklets distribution");
  }
  TH1D *hrNtrk = new TH1D(*(hNtrk[0]));
  hrNtrk->Divide(hNtrk[1]);
  hrNtrk->SetTitle("ratio");
  
  TH1D *hZvtx[2];
  for (Int_t i = 0; i < 2; ++i) {
    hZvtx[i] = fhNtrk[i]->Projection(1,"e");
    hZvtx[i]->Scale(1./hZvtx[i]->GetEntries());
    hZvtx[i]->SetTitle("Z vertex distribution");
  }
  TH1D *hrZvtx = new TH1D(*(hZvtx[0]));
  hrZvtx->Divide(hZvtx[1]);
  hrZvtx->SetTitle("ratio");
  
  TProfile *pMeanTrkVsZvtx[2];
  for (Int_t i = 0; i < 2; ++i) {
    TH2D *hNtrkVsZvtx = fhNtrk[i]->Projection(0,1,"e");
    pMeanTrkVsZvtx[i] = hNtrkVsZvtx->ProfileX(Form("pMeanTrkVsZvtx_%d",i+1));
    pMeanTrkVsZvtx[i]->SetTitle("<Ntracklets> vs Zvtx");
  }
  TH1 *prMeanTrkVsZvtx = pMeanTrkVsZvtx[0]->ProjectionX("prMeanTrkVsZvtx");
  prMeanTrkVsZvtx->Divide(pMeanTrkVsZvtx[1]);
  prMeanTrkVsZvtx->SetTitle("ratio");
  
  // draw histogram
  TCanvas *c = new TCanvas("cC", "cC", 1200, 600);
  c->Divide(3,2);
  c->cd(1);
  gPad->SetLogy();
  hNtrk[0]->Draw();
  hNtrk[1]->SetLineColor(2);
  hNtrk[1]->Draw("sames");
  c->cd(2);
  gPad->SetLogy();
  hZvtx[0]->Draw();
  hZvtx[1]->SetLineColor(2);
  hZvtx[1]->Draw("sames");
  c->cd(3);
  pMeanTrkVsZvtx[0]->Draw();
  pMeanTrkVsZvtx[1]->SetLineColor(2);
  pMeanTrkVsZvtx[1]->Draw("sames");
  c->cd(4);
  hrNtrk->Draw();
  c->cd(5);
  hrZvtx->Draw();
  c->cd(6);
  prMeanTrkVsZvtx->Draw();
  
}

