#include "TFile.h"
#include "THnSparse.h"
#include "TH1D.h"
#include "TH2.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TMath.h"

void DrawMeanTrackletsVsMeanNCh(TString fileNameData = "AnalysisResults.root", Bool_t corr = kFALSE)
{
  
  TFile *file = new TFile(fileNameData.Data(), "read");
  if (!file || !file->IsOpen()) {
    printf("cannot open file %s \n",fileNameData.Data());
    return;
  }
  THnSparse *fhNtrk = static_cast<THnSparse*>(file->FindObjectAny(corr ? "hNtrkCorr" : "hNtrk"));
  if (!fhNtrk) return;
  
//  fhNtrk->GetAxis(0)->SetRangeUser(30, 40);
  
  TH2D *hNchVsNtrk = fhNtrk->Projection(3,0,"e");
  hNchVsNtrk->SetDirectory(0);
  TString pName = corr ? "pMeanNchVsNtrkCorr" : "pMeanNchVsNtrk";
  TProfile *pMeanNchVsNtrk = hNchVsNtrk->ProfileX(pName.Data());
  pMeanNchVsNtrk->SetDirectory(0);
  
  TH2D *hNchVsZvtx = fhNtrk->Projection(3,1,"e");
  hNchVsZvtx->SetDirectory(0);
  TProfile *pMeanNchVsZvtx = hNchVsZvtx->ProfileX("pMeanNchVsZvtx");
  pMeanNchVsZvtx->SetDirectory(0);
  
  file->Close();
  delete file;
  
  // draw histograms
  TString cName = corr ? "cDCCorr" : "cDC";
  TCanvas *c = new TCanvas(cName.Data(), cName.Data());
  gPad->SetLogz();
  hNchVsNtrk->Draw("colz");
  pMeanNchVsNtrk->Draw("sames");
  
  TCanvas *c2 = new TCanvas("cNchVsZvtx", "cNchVsZvtx");
  c2->Divide(2,1);
  c2->cd(1);
  gPad->SetLogz();
  hNchVsZvtx->Draw("colz");
  c2->cd(2);
  pMeanNchVsZvtx->Draw("");
  
}

