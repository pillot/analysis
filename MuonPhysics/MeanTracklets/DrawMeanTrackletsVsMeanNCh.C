#include "TFile.h"
#include "THnSparse.h"
#include "TH1D.h"
#include "TH2.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TF1.h"
#include "TFitResult.h"

Int_t trkBin[] = {1, 8, 13, 19, 30, 49, 60, 100}; // binning for <Ntrk> equalized to the maximum
//Int_t trkBin[] = {1, 7, 12, 17, 28, 45, 55, 92}; // binning for <Ntrk> equalized to the average
//Int_t trkBin[] = {1, 6, 9, 13, 21, 34, 42, 70}; // binning for <Ntrk> equalized to the minimum

Double_t MeanNchVsNtrk(Double_t *x, Double_t *par);
Double_t GetMeanMultWithAlpha(TH1D *hTrk, Int_t minTrk, Int_t maxTrk, TF1 *fitAlpha, Double_t &meanMultErr);
Double_t GetMeanMultWithAdHocFunc(TH1D *hTrk, Int_t minTrk, Int_t maxTrk, TF1 *fMeanNchVsNtrk,
                                  TFitResultPtr &rfMeanNchVsNtrk, Double_t &meanMultErrInt);

//---------------------------------------------------------------------------------------------
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
  
  // Nch vs Ntrk
  TH2D *hNchVsNtrk = fhNtrk->Projection(3,0,"e");
  hNchVsNtrk->SetTitle(corr?"Nch versus corrected Ntrk;Ntrk^{corr};Nch":"Nch versus Ntrk;Ntrk;Nch");
  hNchVsNtrk->SetDirectory(0);
  TString pName = corr ? "pMeanNchVsNtrkCorr" : "pMeanNchVsNtrk";
  TProfile *pMeanNchVsNtrk = hNchVsNtrk->ProfileX(pName.Data());
  pMeanNchVsNtrk->SetDirectory(0);
  
  // <Nch> vs Zvtx per tracklet bin
  Int_t nBins = (Int_t)(sizeof(trkBin) / sizeof(Int_t)) - 1;
  TProfile **pMeanNchVsZvtx = new TProfile*[nBins+1];
  for (Int_t i = 0; i < nBins; ++i) {
    fhNtrk->GetAxis(0)->SetRangeUser(trkBin[i], trkBin[i+1]-1);
    TH2D *hNchVsZvtxInBin = fhNtrk->Projection(3,1,"e");
    pMeanNchVsZvtx[i] = hNchVsZvtxInBin->ProfileX(Form("pMeanNchVsZvtx%d",i+1));
    pMeanNchVsZvtx[i]->SetTitle(Form("<Nch>(Zvtx)/<Nch> in %sNtrk bin %d;Zvtx;<Nch>(Zvtx)/<Nch>",corr?"corrected ":"",i+1));
    pMeanNchVsZvtx[i]->SetDirectory(0);
    delete hNchVsZvtxInBin;
  }
  fhNtrk->GetAxis(0)->SetRange();
  TH2D *hNchVsZvtx = fhNtrk->Projection(3,1,"e");
  hNchVsZvtx->SetTitle("Nch versus Zvtx;Zvtx;Nch");
  hNchVsZvtx->SetDirectory(0);
  pMeanNchVsZvtx[nBins] = hNchVsZvtx->ProfileX("pMeanNchVsZvtx");
  pMeanNchVsZvtx[nBins]->SetTitle("<Nch> versus Zvtx;Zvtx;<Nch>");
  pMeanNchVsZvtx[nBins]->SetDirectory(0);

  file->Close();
  delete file;

  // fit Nch vs Ntrk with y = ax
  TF1 *fitAlpha = new TF1("fitAlpha","[0]*x",1.,100.);
  hNchVsNtrk->Fit(fitAlpha,"NR");

  // fit <Nch> vs Ntrk profile with an ad-hoc polynomial function
  TF1 *fMeanNchVsNtrk = new TF1("fMeanNchVsNtrk",MeanNchVsNtrk,1.,100.,5);
  fMeanNchVsNtrk->SetParameters(1., 1., 1., 1., 10.);
  TFitResultPtr rfMeanNchVsNtrk = pMeanNchVsNtrk->Fit(fMeanNchVsNtrk,"NRS");

  // tracklet distribution
  TH1D *hTrk = hNchVsNtrk->ProjectionX("hTrk", 0, -1, "e");

  // multiplicity distribution and mean multiplicity per tracklet bin
  TH1D **hMult = new TH1D*[nBins+1];
  hMult[nBins] = hNchVsNtrk->ProjectionY("hMult", hTrk->FindBin(1), -1, "e");
  hMult[nBins]->SetTitle(Form("Nch distribution in %sNtrk bins;Nch",corr?"corrected ":""));
  TH1D *hMeanMultVsBin = new TH1D("hMeanMultVsBin", Form("<Nch> versus %sNtrk bin;%sNtrk bin;<Nch>",corr?"corrected ":"",corr?"corrected ":""), nBins, 0.5, nBins+0.5);
  TH1D *hMeanMultVsBinWithAlpha = new TH1D("hMeanMultVsBinWithAlpha", "hMeanMultVsBinWithAlpha", nBins, 0.5, nBins+0.5);
  TH1D *hMeanMultVsBinWithAdHocFunc = new TH1D("hMeanMultVsBinWithAdHocFunc", "hMeanMultVsBinWithAdHocFunc", nBins, 0.5, nBins+0.5);
  for (Int_t i = 0; i < nBins; ++i) {

    hMult[i] = hNchVsNtrk->ProjectionY(Form("hMult%d",i+1), hTrk->FindBin(trkBin[i]), hTrk->FindBin(trkBin[i+1]-1), "e");

    // mean multiplicity in each bin from MC truth
    hMeanMultVsBin->SetBinContent(i+1, hMult[i]->GetMean());
    hMeanMultVsBin->SetBinError(i+1, hMult[i]->GetMeanError());
    printf("meanMult = %f\n",hMult[i]->GetMean());

    // mean multiplicity in each bin computed with constant alpha factor
    Double_t meanMultErr = 0.;
    Double_t meanMult = GetMeanMultWithAlpha(hTrk, trkBin[i], trkBin[i+1]-1, fitAlpha, meanMultErr);
    hMeanMultVsBinWithAlpha->SetBinContent(i+1, meanMult);
    hMeanMultVsBinWithAlpha->SetBinError(i+1, meanMultErr);

    // mean multiplicity in each bin computed with the ad-hoc fit of the <Nch> vs Ntrk profile
    meanMult = GetMeanMultWithAdHocFunc(hTrk, trkBin[i], trkBin[i+1]-1, fMeanNchVsNtrk, rfMeanNchVsNtrk, meanMultErr);
    printf("meanMultCalc = %f\n",meanMult);
    hMeanMultVsBinWithAdHocFunc->SetBinContent(i+1, meanMult);
    hMeanMultVsBinWithAdHocFunc->SetBinError(i+1, meanMultErr);

  }

  // ratio of mean multiplicity
  TH1D *hMeanMultVsBinWithAlphaOverMC = (TH1D*)hMeanMultVsBinWithAlpha->Clone("hMeanMultVsBinWithAlphaOverMC");
  hMeanMultVsBinWithAlphaOverMC->Divide(hMeanMultVsBin);
  hMeanMultVsBinWithAlphaOverMC->SetTitle(Form("calculated <Nch> over MC truth versus %sNtrk bin;%sNtrk bin;<Nch>_{calc}/<Nch>_{MC}",corr?"corrected ":"",corr?"corrected ":""));
  TH1D *hMeanMultVsBinWithAdHocFuncOverMC = (TH1D*)hMeanMultVsBinWithAdHocFunc->Clone("hMeanMultVsBinWithAdHocFuncOverMC");
  hMeanMultVsBinWithAdHocFuncOverMC->Divide(hMeanMultVsBin);

  // draw histograms
  TString cName = corr ? "cDCCorr" : "cDC";
  TCanvas *c = new TCanvas(cName.Data(), cName.Data());
  gPad->SetLogz();
  hNchVsNtrk->Draw("colz");
  pMeanNchVsNtrk->SetLineWidth(2);
  pMeanNchVsNtrk->Draw("sames");
  fitAlpha->SetLineWidth(1);
  fitAlpha->SetLineColor(4);
  fitAlpha->Draw("same");
  fMeanNchVsNtrk->SetLineWidth(1);
  fMeanNchVsNtrk->SetLineColor(4);
  fMeanNchVsNtrk->Draw("same");

  TCanvas *c2 = new TCanvas("cNchVsZvtx", "cNchVsZvtx");
  c2->Divide(2,1);
  c2->cd(1);
  gPad->SetLogz();
  hNchVsZvtx->Draw("colz");
  c2->cd(2);
  pMeanNchVsZvtx[nBins]->Draw("");
  
  TCanvas *c22 = new TCanvas("cNchVsZvtxInNtrkBin", "cNchVsZvtxInNtrkBin", 0, 0, 1500, 600);
  c22->Divide((nBins+1)/2,2);
  for (Int_t i = 0; i < nBins; ++i) {
    c22->cd(i+1);
    pMeanNchVsZvtx[i]->Scale(1./hMult[i]->GetMean());
    pMeanNchVsZvtx[i]->GetYaxis()->SetLabelSize(0.05);
    pMeanNchVsZvtx[i]->GetYaxis()->SetTitleOffset(1.5);
    pMeanNchVsZvtx[i]->Draw("");
    gPad->SetGridy();
  }

  TCanvas *c3 = new TCanvas("cMult", "cMult");
  gPad->SetLogy();
  hMult[nBins]->SetLineColor(1);
  hMult[nBins]->Draw();
  for (Int_t i = 0; i < nBins; ++i) {
    hMult[i]->SetLineColor(i+2);
    hMult[i]->Draw("sames");
  }

  TCanvas *c4 = new TCanvas("cMeanMultVsBin", "cMeanMultVsBin");
  c4->Divide(1,2);
  c4->cd(1);
  hMeanMultVsBin->SetLineColor(1);
  hMeanMultVsBin->Draw();
  hMeanMultVsBinWithAlpha->SetLineColor(2);
  hMeanMultVsBinWithAlpha->Draw("esame");
  hMeanMultVsBinWithAdHocFunc->SetLineColor(4);
  hMeanMultVsBinWithAdHocFunc->Draw("esame");
  c4->cd(2);
  hMeanMultVsBinWithAlphaOverMC->SetLineColor(2);
  hMeanMultVsBinWithAlphaOverMC->GetYaxis()->SetLabelSize(0.05);
  hMeanMultVsBinWithAlphaOverMC->Draw();
  hMeanMultVsBinWithAdHocFuncOverMC->SetLineColor(4);
  hMeanMultVsBinWithAdHocFuncOverMC->Draw("same");

}

//---------------------------------------------------------------------------------------------
Double_t MeanNchVsNtrk(Double_t *x, Double_t *par)
{
  // fit function of <Nch> vs Ntrk

  Double_t nTrk = *x;
  Double_t a = par[0];
  Double_t b = par[1];
  Double_t c = par[2];
  Double_t c2 = par[3];
  Double_t nTrk0 = par[4];
  Double_t a2 = a*c/c2*TMath::Power(nTrk0,c-c2);
  Double_t b2 = (a-a*c/c2)*TMath::Power(nTrk0,c) + b;

  if (nTrk < nTrk0) {
    // y = ax^c + b when x < x0
    return a*TMath::Power(nTrk,c) + b;
  } else {
    // y' = a'x^c' + b' when x >= x0
    // with the following continuity conditions at x = x0:
    // 1) y(x0) = y'(x0)
    // 2) dy/dx(x0) = dy'/dx(x0)
    return a2*TMath::Power(nTrk,c2) + b2;
  }
}

//---------------------------------------------------------------------------------------------
Double_t GetMeanMultWithAlpha(TH1D *hTrk, Int_t minTrk, Int_t maxTrk, TF1 *fitAlpha, Double_t &meanMultErr)
{
  // mean multiplicity in the given Ntrk range computed with constant alpha factor

  hTrk->SetAxisRange(minTrk, maxTrk);
  Double_t meanTrk = hTrk->GetMean();
  Double_t meanTrkErr = hTrk->GetMeanError();
  hTrk->GetXaxis()->SetRange();

  Double_t alpha = fitAlpha->GetParameter(0);
  Double_t alphaErr = fitAlpha->GetParError(0);

  meanMultErr = TMath::Sqrt(meanTrk*meanTrk*alphaErr*alphaErr + alpha*alpha*meanTrkErr*meanTrkErr);
  return alpha * meanTrk;
}

//---------------------------------------------------------------------------------------------
Double_t GetMeanMultWithAdHocFunc(TH1D *hTrk, Int_t minTrk, Int_t maxTrk, TF1 *fMeanNchVsNtrk,
                                  TFitResultPtr &rfMeanNchVsNtrk, Double_t &meanMultErrInt)
{
  // mean multiplicity in the given Ntrk range computed with the ad-hoc fit of the <Nch> vs Ntrk profile

  /*
   hTrk->SetAxisRange(minTrk, maxTrk);
   Double_t meanTrk = hTrk->GetMean();
   Double_t meanTrkErr = hTrk->GetMeanError();
   hTrk->GetXaxis()->SetRange();

   Double_t meanMult = fMeanNchVsNtrk->Eval(meanTrk);
   printf("meanMult = %f\n", meanMult);
   meanMultErr = meanTrkErr * meanMult / meanTrk;
   return meanMult;
   */

  Double_t SumNEvMeanMult(0.);
  Double_t SumNEvErr2(0.);
  Double_t SumNEvErr2MeanMult(0.);
  Double_t SumNEvErr2MeanMult2(0.);
  Double_t SumNEv2MeanMultErr2(0.);
  Double_t SumNEv(0.);
  for (Int_t j = minTrk; j <= maxTrk; ++j) {
    Double_t nTrk = j;
    Double_t nEv = hTrk->GetBinContent(hTrk->FindBin(nTrk));
    Double_t nEvErr2 = hTrk->GetBinError(hTrk->FindBin(nTrk)) * hTrk->GetBinError(hTrk->FindBin(nTrk));
    Double_t meanMult = fMeanNchVsNtrk->Eval(nTrk);
    Double_t meanMultErr(0);
    rfMeanNchVsNtrk->GetConfidenceIntervals(1, 1, 1, &nTrk, &meanMultErr, 0.683, kFALSE);
    SumNEvMeanMult += nEv * meanMult;
    SumNEvErr2 += nEvErr2;
    SumNEvErr2MeanMult += nEvErr2 * meanMult;
    SumNEvErr2MeanMult2 += nEvErr2 * meanMult * meanMult;
    SumNEv2MeanMultErr2 += nEv * nEv * meanMultErr * meanMultErr;
    SumNEv += nEv;
  }
  Double_t meanMultInt = (SumNEv > 0.) ? SumNEvMeanMult/SumNEv : 0.;
  meanMultErrInt = (SumNEv > 0.) ? TMath::Sqrt(SumNEvErr2MeanMult2 + meanMultInt*meanMultInt*SumNEvErr2 - 2.*meanMultInt*SumNEvErr2MeanMult + SumNEv2MeanMultErr2)/SumNEv : 0.;
  return meanMultInt;
}

