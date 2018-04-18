#include "TROOT.h"
#include "TFile.h"
#include "THnSparse.h"
#include "TH1D.h"
#include "TH2.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TMath.h"

void DrawMeanTracklets(TString fileNameData = "AnalysisResults.root", Bool_t corr = kFALSE, Bool_t same = kFALSE)
{
  
  TFile *file = new TFile(fileNameData.Data(), "read");
  if (!file || !file->IsOpen()) {
    printf("cannot open file %s \n",fileNameData.Data());
    return;
  }
  THnSparse *fhNtrk = static_cast<THnSparse*>(file->FindObjectAny(corr ? "hNtrkCorr" : "hNtrk"));
  if (!fhNtrk) return;
  //fhNtrk->GetAxis(2)->SetRangeUser(244351, 244351);
  
  TH2D *hNtrkVsZvtx = fhNtrk->Projection(0,1,"e");
  hNtrkVsZvtx->SetDirectory(0);
  TString pName = corr ? "pMeanTrkVsZvtxCorr" : "pMeanTrkVsZvtx";
  while (gROOT->FindObject(pName.Data())) pName += "0";
  TProfile *pMeanTrkVsZvtx = hNtrkVsZvtx->ProfileX(pName.Data());
  pMeanTrkVsZvtx->SetDirectory(0);
  pName = corr ? "pMeanZvtxVsTrkCorr" : "pMeanZvtxVsTrk";
  while (gROOT->FindObject(pName.Data())) pName += "0";
  TProfile *pMeanZvtxVsTrk = hNtrkVsZvtx->ProfileY(pName.Data());
  pMeanZvtxVsTrk->SetDirectory(0);
  //fhNtrk->GetAxis(0)->SetRangeUser(0, 10);
  TH2D *hNtrkVsRun = fhNtrk->Projection(0,2,"e");
  hNtrkVsRun->SetDirectory(0);
  TProfile *pMeanTrkVsRun = hNtrkVsRun->ProfileX("pMeanTrkVsRun");
  pMeanTrkVsRun->SetDirectory(0);
  //fhNtrk->GetAxis(0)->SetRange();
  TH2D *hZvtxVsRun = fhNtrk->Projection(1,2,"e");
  hZvtxVsRun->SetDirectory(0);
  TProfile *pMeanZvtxVsRun = hZvtxVsRun->ProfileX("pMeanZvtxVsRun",1,-1,"s");
  pMeanZvtxVsRun->SetDirectory(0);
  
  file->Close();
  delete file;
  
  // re-plot results vs run to reassemble them
  TH1D *hMeanTrkVsRun = new TH1D("hMeanTrkVsRun","hMeanTrkVsRun",1000,-0.5,999.5);
  hMeanTrkVsRun->Sumw2();
  TH1D *hMeanZvtxVsRun = new TH1D("hMeanZvtxVsRun","hMeanZvtxVsRun",1000,-0.5,999.5);
  hMeanZvtxVsRun->Sumw2();
  TH1D *hRMSZvtxVsRun = new TH1D("hRMSZvtxVsRun","hRMSZvtxVsRun",1000,-0.5,999.5);
  hRMSZvtxVsRun->Sumw2();
  Int_t iBin = 0;
  for (Int_t i = 0; i <= pMeanTrkVsRun->GetNbinsX(); ++i) {
    Double_t mean = pMeanTrkVsRun->GetBinContent(i);
    if (mean < 1.e-6) continue;
    ++iBin;
    hMeanTrkVsRun->SetBinContent(iBin,mean);
    hMeanTrkVsRun->SetBinError(iBin,pMeanTrkVsRun->GetBinError(i));
    hMeanTrkVsRun->GetXaxis()->SetBinLabel(iBin,Form("%d",(Int_t)(pMeanTrkVsRun->GetXaxis()->GetBinCenter(i))));
    hMeanZvtxVsRun->SetBinContent(iBin,pMeanZvtxVsRun->GetBinContent(i));
    hMeanZvtxVsRun->SetBinError(iBin,pMeanZvtxVsRun->GetBinError(i)/TMath::Sqrt(pMeanZvtxVsRun->GetBinEntries(i)));
    hMeanZvtxVsRun->GetXaxis()->SetBinLabel(iBin,Form("%d",(Int_t)(pMeanZvtxVsRun->GetXaxis()->GetBinCenter(i))));
    hRMSZvtxVsRun->SetBinContent(iBin,pMeanZvtxVsRun->GetBinError(i));
    hRMSZvtxVsRun->SetBinError(iBin,0.);
    hRMSZvtxVsRun->GetXaxis()->SetBinLabel(iBin,Form("%d",(Int_t)(pMeanZvtxVsRun->GetXaxis()->GetBinCenter(i))));
  }
  hMeanTrkVsRun->GetXaxis()->SetRange(1,iBin);
  hMeanZvtxVsRun->GetXaxis()->SetRange(1,iBin);
  hRMSZvtxVsRun->GetXaxis()->SetRange(1,iBin);

  // Ntrk distribution in bins of zVtx
  const Int_t nBinZvtx = 10; // 10 bins of 2 cm
  Double_t dZvtx = (hNtrkVsZvtx->GetXaxis()->GetXmax()-hNtrkVsZvtx->GetXaxis()->GetXmin())/nBinZvtx;
  TH1D *hNtrk[nBinZvtx+1];
  TH1D *hNtrkRatio[nBinZvtx];
  hNtrk[nBinZvtx] = hNtrkVsZvtx->ProjectionY("hNtrk", 0, -1, "e");
  hNtrk[nBinZvtx]->SetTitle(corr?"corrected Ntrk;Ntrk^{corr}":"Ntrk;Ntrk");
  hNtrk[nBinZvtx]->Scale(1./hNtrk[nBinZvtx]->GetEntries());
  for (Int_t iz = 0; iz < nBinZvtx; ++iz) {
    Double_t zMin = hNtrkVsZvtx->GetXaxis()->GetXmin() + dZvtx*iz;
    Double_t zMax = hNtrkVsZvtx->GetXaxis()->GetXmin() + dZvtx*(iz+1);
    hNtrk[iz] = hNtrkVsZvtx->ProjectionY(Form("hNtrk%d",iz+1), hNtrkVsZvtx->GetXaxis()->FindBin(zMin+0.1), hNtrkVsZvtx->GetXaxis()->FindBin(zMax-0.1), "e");
    hNtrk[iz]->SetTitle(Form(corr?"corrected Ntrk in %g < zVtx < %g;Ntrk^{corr}":"Ntrk in %g < zVtx < %g;Ntrk",zMin,zMax));
    hNtrk[iz]->Scale(1./hNtrk[iz]->GetEntries());
    hNtrkRatio[iz] = (TH1D*)hNtrk[iz]->Clone(Form("hNtrkRatio%d",iz+1));
    hNtrkRatio[iz]->SetTitle("ratio to zVtx-integrated distribution");
    hNtrkRatio[iz]->Divide(hNtrk[nBinZvtx]);
  }

  // draw histogram
  TString cName = corr ? "cDCorr" : "cD";
  
  TCanvas *c = 0x0;
  if (!same || !(c = static_cast<TCanvas*>(gROOT->GetListOfCanvases()->FindObject(cName.Data())))) {
    c = new TCanvas(cName.Data(), cName.Data(), 1200, 600);
    c->Divide(4,2);
    same = kFALSE;
  }
  Int_t color = same ? 2 : 4;
  c->cd(1);
  gPad->SetLogy();
  TH1 *h = fhNtrk->Projection(0,"e");
  h->Scale(1./h->GetEntries());
  h->SetLineColor(color);
  h->Draw(same ? "sames" : 0x0);
  c->cd(2);
  gPad->SetLogy();
  h = fhNtrk->Projection(1,"e");
  h->Scale(1./h->GetEntries());
  h->SetLineColor(color);
  h->Draw(same ? "sames" : 0x0);
  c->cd(3);
  gPad->SetLogz();
  hNtrkVsZvtx->Draw("colz");
  c->cd(4);
  pMeanTrkVsZvtx->SetLineColor(color);
  pMeanTrkVsZvtx->Draw(same ? "sames" : 0x0);
  c->cd(5);
  hMeanTrkVsRun->SetLineColor(color);
  hMeanTrkVsRun->Draw(same ? "sames" : 0x0);
  c->cd(6);
  hMeanZvtxVsRun->SetLineColor(color);
  hMeanZvtxVsRun->Draw(same ? "sames" : 0x0);
  c->cd(7);
  hRMSZvtxVsRun->SetLineColor(color);
  hRMSZvtxVsRun->Draw(same ? "sames" : 0x0);
  c->cd(8);
  pMeanZvtxVsTrk->SetLineColor(color);
  pMeanZvtxVsTrk->Draw(same ? "sames" : 0x0);
  
  TCanvas *c2 = new TCanvas("cNtrk", "cNtrk", 10, 10, 800, 800);
  c2->Divide(1,2);
  color = 2;
  for (Int_t iz = 0; iz < nBinZvtx; ++iz) {
    c2->cd(1);
    gPad->SetLogy();
    hNtrk[iz]->SetLineColor(color);
    hNtrk[iz]->Draw((iz==0)?"":"sames");
    c2->cd(2);
    hNtrkRatio[iz]->SetLineColor(color);
    hNtrkRatio[iz]->Draw((iz==0)?"":"same");
    ++color;
    if (color == 5 || color == 10) ++color;
  }

}

