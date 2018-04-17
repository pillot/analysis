#include "TFile.h"
#include "THnSparse.h"
#include "TH1D.h"
#include "TH2.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TF1.h"
#include "TFitResult.h"

// define the binning in term of (corrected) number of tracklet
Int_t trkBin[] = {1, 8, 13, 19, 30, 49, 60, 100}; // binning for <Ntrk> equalized to the maximum
//Int_t trkBin[] = {1, 7, 12, 17, 28, 45, 55, 92}; // binning for <Ntrk> equalized to the average
//Int_t trkBin[] = {1, 6, 9, 13, 21, 34, 42, 70}; // binning for <Ntrk> equalized to the minimum

Double_t MeanNchVsNtrk(Double_t *x, Double_t *par);
Double_t RelativeTrkOverRelativeMeanNchWithAlpha(Double_t *x, Double_t *par);
Double_t RelativeTrkOverRelativeMeanNchWithAdHocFunc(Double_t *x, Double_t *par);
Double_t GetMeanMultWithAlpha(TH1D *hTrk, Int_t minTrk, Int_t maxTrk, TF1 *fitAlpha, Double_t &meanMultErr);
Double_t GetMeanMultWithAdHocFunc(TH1D *hTrk, Int_t minTrk, Int_t maxTrk, TF1 *fMeanNchVsNtrk,
                                  TFitResultPtr &rfMeanNchVsNtrk, Double_t &meanMultErrInt);
void DrawMeanTrackletsVsMeanNChFlat(TH2D *hNchVsNtrk, Bool_t corr, TF1 *fitAlpha, TF1 *fMeanNchVsNtrk, TFitResultPtr &rfMeanNchVsNtrk);

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
  
  // Nch vs Ntrk in bins of Zvtx and integrated
  const Int_t nBinZvtx = 10; // 10 bins of 2 cm
  Double_t dZvtx = (fhNtrk->GetAxis(1)->GetXmax()-fhNtrk->GetAxis(1)->GetXmin())/nBinZvtx;
  TH2D *hNchVsNtrk[nBinZvtx+1];
  TProfile *pMeanNchVsNtrk[nBinZvtx+1];
  for (Int_t iz = 0; iz < nBinZvtx; ++iz) {
    Double_t zMin = fhNtrk->GetAxis(1)->GetXmin() + dZvtx*iz;
    Double_t zMax = fhNtrk->GetAxis(1)->GetXmin() + dZvtx*(iz+1);
    fhNtrk->GetAxis(1)->SetRangeUser(zMin, zMax);
    hNchVsNtrk[iz] = fhNtrk->Projection(3,0,"e");
    hNchVsNtrk[iz]->SetTitle(Form(corr?"Nch versus corrected Ntrk in %g < zVtx < %g;Ntrk^{corr};Nch":"Nch versus Ntrk;Ntrk;Nch",zMin,zMax));
    hNchVsNtrk[iz]->SetDirectory(0);
    pMeanNchVsNtrk[iz] = hNchVsNtrk[iz]->ProfileX(Form(corr ? "pMeanNchVsNtrkCorr%d" : "pMeanNchVsNtrk%d",iz+1));
    pMeanNchVsNtrk[iz]->SetDirectory(0);
  }
  fhNtrk->GetAxis(1)->SetRange();
  hNchVsNtrk[nBinZvtx] = fhNtrk->Projection(3,0,"e");
  hNchVsNtrk[nBinZvtx]->SetTitle(corr?"Nch versus corrected Ntrk;Ntrk^{corr};Nch":"Nch versus Ntrk;Ntrk;Nch");
  hNchVsNtrk[nBinZvtx]->SetDirectory(0);
  pMeanNchVsNtrk[nBinZvtx] = hNchVsNtrk[nBinZvtx]->ProfileX(corr ? "pMeanNchVsNtrkCorr" : "pMeanNchVsNtrk");
  pMeanNchVsNtrk[nBinZvtx]->SetDirectory(0);
  TProfile *pMeanNtrkVsNch = hNchVsNtrk[nBinZvtx]->ProfileY(corr ? "pMeanNtrCorrVsNch" : "pMeanNtrkVsNch");
  pMeanNtrkVsNch->SetDirectory(0);

  // <Nch> vs Zvtx per tracklet bin and integrated
  Int_t nBins = (Int_t)(sizeof(trkBin) / sizeof(Int_t)) - 1;
  TProfile **pMeanNchVsZvtx = new TProfile*[nBins+1];
  TH1D **hRMSNchVsZvtx = new TH1D*[nBins+1];
  for (Int_t i = 0; i < nBins; ++i) {
    fhNtrk->GetAxis(0)->SetRangeUser(trkBin[i], trkBin[i+1]-1);
    TH2D *hNchVsZvtxInBin = fhNtrk->Projection(3,1,"e");
    pMeanNchVsZvtx[i] = hNchVsZvtxInBin->ProfileX(Form("pMeanNchVsZvtx%d",i+1));
    pMeanNchVsZvtx[i]->SetTitle(Form("<Nch>(Zvtx)/<Nch> in %sNtrk bin %d;Zvtx;<Nch>(Zvtx)/<Nch>",corr?"corrected ":"",i+1));
    pMeanNchVsZvtx[i]->SetDirectory(0);
    hRMSNchVsZvtx[i] = new TH1D(Form("hRMSNchVsZvtx%d",i+1), Form("Nch RMS(Zvtx)/RMS in %sNtrk bin %d;Zvtx;Nch RMS(Zvtx)/RMS",corr?"corrected ":"",i+1), pMeanNchVsZvtx[i]->GetNbinsX(), pMeanNchVsZvtx[i]->GetXaxis()->GetXmin(), pMeanNchVsZvtx[i]->GetXaxis()->GetXmax());
    for (Int_t iz = 1; iz <= pMeanNchVsZvtx[i]->GetNbinsX(); ++iz) {
      hRMSNchVsZvtx[i]->SetBinContent(iz, pMeanNchVsZvtx[i]->GetBinError(iz)*TMath::Sqrt(pMeanNchVsZvtx[i]->GetBinEntries(iz)));
      hRMSNchVsZvtx[i]->SetBinError(iz, 0.);
    }
    hRMSNchVsZvtx[i]->SetDirectory(0);
    delete hNchVsZvtxInBin;
  }
  fhNtrk->GetAxis(0)->SetRange();
  TH2D *hNchVsZvtx = fhNtrk->Projection(3,1,"e");
  hNchVsZvtx->SetTitle("Nch versus Zvtx;Zvtx;Nch");
  hNchVsZvtx->SetDirectory(0);
  pMeanNchVsZvtx[nBins] = hNchVsZvtx->ProfileX("pMeanNchVsZvtx");
  pMeanNchVsZvtx[nBins]->SetTitle("<Nch> versus Zvtx;Zvtx;<Nch>");
  pMeanNchVsZvtx[nBins]->SetDirectory(0);
  hRMSNchVsZvtx[nBins] = new TH1D("hRMSNchVsZvtx", "Nch dispersion vs Zvtx;Zvtx;Nch dispersion", pMeanNchVsZvtx[nBins]->GetNbinsX(), pMeanNchVsZvtx[nBins]->GetXaxis()->GetXmin(), pMeanNchVsZvtx[nBins]->GetXaxis()->GetXmax());
  for (Int_t iz = 1; iz <= pMeanNchVsZvtx[nBins]->GetNbinsX(); ++iz) {
    hRMSNchVsZvtx[nBins]->SetBinContent(iz, pMeanNchVsZvtx[nBins]->GetBinError(iz)*TMath::Sqrt(pMeanNchVsZvtx[nBins]->GetBinEntries(iz)));
    hRMSNchVsZvtx[nBins]->SetBinError(iz, 0.);
  }
  hRMSNchVsZvtx[nBins]->SetDirectory(0);

  file->Close();
  delete file;

  // fit Nch vs Ntrk with y = ax
  TF1 *fitAlpha[nBinZvtx+1];
  for (Int_t iz = 0; iz < nBinZvtx+1; ++iz) {
    fitAlpha[iz] = new TF1(Form("fitAlpha%d",iz+1),"[0]*x",1.,100.);
    hNchVsNtrk[iz]->Fit(fitAlpha[iz],"NR");
  }

  // fit <Nch> vs Ntrk profile with an ad-hoc polynomial function
  TF1 *fMeanNchVsNtrk[nBinZvtx+1];
  TFitResultPtr rfMeanNchVsNtrk[nBinZvtx+1];
  for (Int_t iz = 0; iz < nBinZvtx+1; ++iz) {
    fMeanNchVsNtrk[iz] = new TF1(Form("fMeanNchVsNtrk%d",iz+1),MeanNchVsNtrk,1.,100.,5);
    fMeanNchVsNtrk[iz]->SetParameters(1., 1., 1., 1., 10.);
    fMeanNchVsNtrk[iz]->SetParLimits(4, 2., 50.);
    rfMeanNchVsNtrk[iz] = pMeanNchVsNtrk[iz]->Fit(fMeanNchVsNtrk[iz],"NRS");
  }

  // tracklet distribution and integrated <Ntrk> for Ntrk >= 1 versus Zvtx
  TH1D *hTrk[nBinZvtx+1];
  Double_t meanNtrkInt[nBinZvtx+1];
  for (Int_t iz = 0; iz < nBinZvtx+1; ++iz) {
    hTrk[iz] = hNchVsNtrk[iz]->ProjectionX(Form("hTrk%d",iz+1), 0, -1, "e");
    hTrk[iz]->SetAxisRange(1, hTrk[iz]->GetXaxis()->GetXmax());
    meanNtrkInt[iz] = hTrk[iz]->GetMean();
    printf(corr?"integrated <NtrkCorr> = %f\n":"<Ntrk> = %f\n",meanNtrkInt[iz]);
    hTrk[iz]->GetXaxis()->SetRange();
  }

  // integrated <Nch> for Ntrk >= 1 versus Zvtx
  Double_t meanNchInt[nBinZvtx+1];
  Double_t meanNchIntWithAlpha[nBinZvtx+1];
  Double_t meanNchIntWithAdHocFunc[nBinZvtx+1];
  for (Int_t iz = 0; iz < nBinZvtx+1; ++iz) {
    TH1D *hTmp = hNchVsNtrk[iz]->ProjectionY("hTmp", hTrk[iz]->FindBin(1), -1, "e");
    meanNchInt[iz] = hTmp->GetMean();
    printf(corr?"integrated <Nch> = %f\n":"<Nch> = %f\n",meanNchInt[iz]);
    delete hTmp;
    Double_t err;
    meanNchIntWithAlpha[iz] = GetMeanMultWithAlpha(hTrk[iz], 1., hTrk[iz]->GetBinCenter(hTrk[iz]->GetNbinsX()), fitAlpha[iz], err);
    meanNchIntWithAdHocFunc[iz] = GetMeanMultWithAdHocFunc(hTrk[iz], 1., hTrk[iz]->GetBinCenter(hTrk[iz]->GetNbinsX()), fMeanNchVsNtrk[iz], rfMeanNchVsNtrk[iz], err);
  }

  // relative Ntrk over relative <Nch> ((Ntrk^i/<Ntrk>) / (<Nch>^i/<Nch>)) versus Zvtx
  TH1D *hRelativeTrkOverRelativeMeanNch[nBinZvtx+1];
  TF1 *fRelativeTrkOverRelativeMeanNchWithAlpha[nBinZvtx+1];
  TF1 *fRelativeTrkOverRelativeMeanNchWithAdHocFunc[nBinZvtx+1];
  for (Int_t iz = 0; iz < nBinZvtx+1; ++iz) {
    hRelativeTrkOverRelativeMeanNch[iz] = new TH1D(Form("hRelativeTrk%sOverRelativeMeanNch%d",corr?"Corr":"",iz+1), Form("relative Ntrk%s over relative <Nch>;Ntrk^{%s};(Ntrk%s^{i}/<Ntrk%s>)/(<Nch>^{i}/<Nch>)",corr?"Corr":"",corr?"corr":"",corr?"Corr":"",corr?"Corr":""), hTrk[iz]->GetNbinsX(), hTrk[iz]->GetXaxis()->GetXmin(), hTrk[iz]->GetXaxis()->GetXmax());
    for (Int_t i = 1; i <= hTrk[iz]->GetNbinsX(); ++i) {
      if (pMeanNchVsNtrk[iz]->GetBinContent(i) > 0. && pMeanNchVsNtrk[iz]->GetBinError(i) > 0.) {
        Double_t relativeTrkOverRelativeMeanNch = hTrk[iz]->GetBinCenter(i)/meanNtrkInt[iz]*meanNchInt[iz]/pMeanNchVsNtrk[iz]->GetBinContent(i);
        hRelativeTrkOverRelativeMeanNch[iz]->SetBinContent(i,relativeTrkOverRelativeMeanNch);
        hRelativeTrkOverRelativeMeanNch[iz]->SetBinError(i,pMeanNchVsNtrk[iz]->GetBinError(i)/pMeanNchVsNtrk[iz]->GetBinContent(i)*relativeTrkOverRelativeMeanNch);
      }
      fRelativeTrkOverRelativeMeanNchWithAlpha[iz] = new TF1(Form("fRelativeTrkOverRelativeMeanNchWithAlpha%d",iz+1), RelativeTrkOverRelativeMeanNchWithAlpha,1.,100.,3);
      fRelativeTrkOverRelativeMeanNchWithAlpha[iz]->SetParameters(fitAlpha[iz]->GetParameter(0), meanNtrkInt[iz], meanNchIntWithAlpha[iz]);
      fRelativeTrkOverRelativeMeanNchWithAdHocFunc[iz] = new TF1(Form("fRelativeTrkOverRelativeMeanNchWithAdHocFunc%d",iz+1), RelativeTrkOverRelativeMeanNchWithAdHocFunc,1.,100.,7);
      fRelativeTrkOverRelativeMeanNchWithAdHocFunc[iz]->SetParameters(fMeanNchVsNtrk[iz]->GetParameter(0), fMeanNchVsNtrk[iz]->GetParameter(1), fMeanNchVsNtrk[iz]->GetParameter(2), fMeanNchVsNtrk[iz]->GetParameter(3), fMeanNchVsNtrk[iz]->GetParameter(4), meanNtrkInt[iz], meanNchIntWithAdHocFunc[iz]);
    }
  }

  // integrate <Ntrk> for Nch >= 1
  TH1D *hTmp = hNchVsNtrk[nBinZvtx]->ProjectionX("hTmp", hNchVsNtrk[nBinZvtx]->GetYaxis()->FindBin(1), -1, "e");
  Double_t meanNtrkIntNchPos = hTmp->GetMean();
  delete hTmp;

  // Nch distribution and integrated <Nch> for Nch >= 1
  TH1D *hMultInt = hNchVsNtrk[nBinZvtx]->ProjectionY("hMultInt", 0, -1, "e");
  hMultInt->SetAxisRange(1, hMultInt->GetXaxis()->GetXmax());
  Double_t meanNchIntNchPos = hMultInt->GetMean();
  hMultInt->GetXaxis()->SetRange();

  // relative <Ntrk> over relative Nch ((<Ntrk>^i/<Ntrk>) / (Nch^i/<Nch>))
  TH1D *hRelativeMeanTrkOverRelativeNch = new TH1D(Form("hRelativeMeanTrk%sOverRelativeNch",corr?"Corr":""), Form("relative <Ntrk%s> over relative Nch;Nch;(<Ntrk%s>^{i}/<Ntrk%s>)/(Nch^{i}/<Nch>)",corr?"Corr":"",corr?"Corr":"",corr?"Corr":""), hMultInt->GetNbinsX(), hMultInt->GetXaxis()->GetXmin(), hMultInt->GetXaxis()->GetXmax());
  for (Int_t i = 2; i <= hMultInt->GetNbinsX(); ++i) {
    if (pMeanNtrkVsNch->GetBinContent(i) > 0. && pMeanNtrkVsNch->GetBinError(i) > 0.) {
      Double_t relativeMeanTrkOverRelativeNch = pMeanNtrkVsNch->GetBinContent(i)/meanNtrkIntNchPos*meanNchIntNchPos/hMultInt->GetBinCenter(i);
      hRelativeMeanTrkOverRelativeNch->SetBinContent(i,relativeMeanTrkOverRelativeNch);
      hRelativeMeanTrkOverRelativeNch->SetBinError(i,pMeanNtrkVsNch->GetBinError(i)/pMeanNtrkVsNch->GetBinContent(i)*relativeMeanTrkOverRelativeNch);
    }
  }

  // multiplicity distribution and mean multiplicity per tracklet bin
  TH1D **hMult = new TH1D*[nBins+1];
  hMult[nBins] = hNchVsNtrk[nBinZvtx]->ProjectionY("hMult", hTrk[nBinZvtx]->FindBin(1), -1, "e");
  hMult[nBins]->SetTitle(Form("Nch distribution in %sNtrk bins;Nch",corr?"corrected ":""));
  TH1D *hMeanMultVsBin[nBinZvtx+1];
  TH1D *hMeanMultVsBinWithAlpha[nBinZvtx+1][2];
  TH1D *hMeanMultVsBinWithAdHocFunc[nBinZvtx+1][2];
  for (Int_t iz = 0; iz < nBinZvtx+1; ++iz) {
    hMeanMultVsBin[iz] = new TH1D(Form("hMeanMultVsBin%d",iz+1), Form("<Nch> versus %sNtrk bin;%sNtrk bin;<Nch>",corr?"corrected ":"",corr?"corrected ":""), nBins, 0.5, nBins+0.5);
    hMeanMultVsBinWithAlpha[iz][0] = new TH1D(Form("hMeanMultVsBinWithAlphaFromZvtxBin%d",iz+1), Form("hMeanMultVsBinWithAlphaFromZvtxBin%d",iz+1), nBins, 0.5, nBins+0.5);
    hMeanMultVsBinWithAdHocFunc[iz][0] = new TH1D(Form("hMeanMultVsBinWithAdHocFuncFromZvtxBin%d",iz+1), Form("hMeanMultVsBinWithAdHocFuncFromZvtxBin%d",iz+1), nBins, 0.5, nBins+0.5);
    hMeanMultVsBinWithAlpha[iz][1] = new TH1D(Form("hMeanMultVsBinWithAlpha%d",iz+1), Form("hMeanMultVsBinWithAlpha%d",iz+1), nBins, 0.5, nBins+0.5);
    hMeanMultVsBinWithAdHocFunc[iz][1] = new TH1D(Form("hMeanMultVsBinWithAdHocFunc%d",iz+1), Form("hMeanMultVsBinWithAdHocFunc%d",iz+1), nBins, 0.5, nBins+0.5);
  }
  for (Int_t i = 0; i < nBins; ++i) {

    hMult[i] = hNchVsNtrk[nBinZvtx]->ProjectionY(Form("hMult%d",i+1), hTrk[nBinZvtx]->FindBin(trkBin[i]), hTrk[nBinZvtx]->FindBin(trkBin[i+1]-1), "e");

    for (Int_t iz = 0; iz < nBinZvtx+1; ++iz) {

      // mean Ntrk in zVtx bin "iz" and Ntrk bin "i"
      hTrk[iz]->SetAxisRange(trkBin[i], trkBin[i+1]-1);
      printf(corr?"<NtrkCorr> = %f":"<Ntrk> = %f",hTrk[iz]->GetMean());
      hTrk[iz]->GetXaxis()->SetRange();

      // mean multiplicity in zVtx bin "iz" and Ntrk bin "i" from MC truth
      hTmp = hNchVsNtrk[iz]->ProjectionY("hTmp", hTrk[iz]->FindBin(trkBin[i]), hTrk[iz]->FindBin(trkBin[i+1]-1), "e");
      hMeanMultVsBin[iz]->SetBinContent(i+1, hTmp->GetMean());
      hMeanMultVsBin[iz]->SetBinError(i+1, hTmp->GetMeanError());
      printf(" <Nch>(Calc) = %f",hTmp->GetMean());
      delete hTmp;

      // mean multiplicity integrated over zVtx in Ntrk bin "i" computed with constant alpha factor from Zvtx bin "iz"
      Double_t meanMultErr = 0.;
      Double_t meanMult = GetMeanMultWithAlpha(hTrk[nBinZvtx], trkBin[i], trkBin[i+1]-1, fitAlpha[iz], meanMultErr);
      printf(" (%f)",meanMult);
      hMeanMultVsBinWithAlpha[iz][0]->SetBinContent(i+1, meanMult);
      hMeanMultVsBinWithAlpha[iz][0]->SetBinError(i+1, meanMultErr);

      // mean multiplicity in zVtx bin "iz" and Ntrk bin "i" computed with constant alpha factor integrated over zVtx
      meanMultErr = 0.;
      meanMult = GetMeanMultWithAlpha(hTrk[iz], trkBin[i], trkBin[i+1]-1, fitAlpha[nBinZvtx], meanMultErr);
      printf(" (%f)",meanMult);
      hMeanMultVsBinWithAlpha[iz][1]->SetBinContent(i+1, meanMult);
      hMeanMultVsBinWithAlpha[iz][1]->SetBinError(i+1, meanMultErr);

      // mean multiplicity integrated over zVtx in Ntrk bin "i" computed with the ad-hoc fit of the <Nch> vs Ntrk profile from Zvtx bin "iz"
      meanMult = GetMeanMultWithAdHocFunc(hTrk[nBinZvtx], trkBin[i], trkBin[i+1]-1, fMeanNchVsNtrk[iz], rfMeanNchVsNtrk[iz], meanMultErr);
      printf(" (%f)",meanMult);
      hMeanMultVsBinWithAdHocFunc[iz][0]->SetBinContent(i+1, meanMult);
      hMeanMultVsBinWithAdHocFunc[iz][0]->SetBinError(i+1, meanMultErr);

      // mean multiplicity in zVtx bin "iz" and Ntrk bin "i" computed with the ad-hoc fit of the <Nch> vs Ntrk profile integrated over zVtx
      meanMult = GetMeanMultWithAdHocFunc(hTrk[iz], trkBin[i], trkBin[i+1]-1, fMeanNchVsNtrk[nBinZvtx], rfMeanNchVsNtrk[nBinZvtx], meanMultErr);
      printf(" (%f)\n",meanMult);
      hMeanMultVsBinWithAdHocFunc[iz][1]->SetBinContent(i+1, meanMult);
      hMeanMultVsBinWithAdHocFunc[iz][1]->SetBinError(i+1, meanMultErr);

    }

  }

  // ratio of mean multiplicity
  TH1D *hMeanMultVsBinWithAlphaOverMC[nBinZvtx+1][2];
  TH1D *hMeanMultVsBinWithAdHocFuncOverMC[nBinZvtx+1][2];
  for (Int_t iz = 0; iz < nBinZvtx+1; ++iz) {

    hMeanMultVsBinWithAlphaOverMC[iz][0] = (TH1D*)hMeanMultVsBinWithAlpha[iz][0]->Clone(Form("hMeanMultVsBinWithAlphaFromZvtxBin%dOverMC",iz+1));
    hMeanMultVsBinWithAlphaOverMC[iz][0]->Divide(hMeanMultVsBin[nBinZvtx]);
    hMeanMultVsBinWithAlphaOverMC[iz][0]->SetTitle(Form("calculated <Nch> over MC truth versus %sNtrk bin;%sNtrk bin;<Nch>_{calc}/<Nch>_{MC}",corr?"corrected ":"",corr?"corrected ":""));
    hMeanMultVsBinWithAdHocFuncOverMC[iz][0] = (TH1D*)hMeanMultVsBinWithAdHocFunc[iz][0]->Clone(Form("hMeanMultVsBinWithAdHocFuncFromZvtxBin%dOverMC",iz+1));
    hMeanMultVsBinWithAdHocFuncOverMC[iz][0]->Divide(hMeanMultVsBin[nBinZvtx]);

    hMeanMultVsBinWithAlphaOverMC[iz][1] = (TH1D*)hMeanMultVsBinWithAlpha[iz][1]->Clone(Form("hMeanMultVsBinWithAlphaOverMC%d",iz+1));
    hMeanMultVsBinWithAlphaOverMC[iz][1]->Divide(hMeanMultVsBin[iz]);
    hMeanMultVsBinWithAlphaOverMC[iz][1]->SetTitle(Form("calculated <Nch> over MC truth versus %sNtrk bin;%sNtrk bin;<Nch>_{calc}/<Nch>_{MC}",corr?"corrected ":"",corr?"corrected ":""));
    hMeanMultVsBinWithAdHocFuncOverMC[iz][1] = (TH1D*)hMeanMultVsBinWithAdHocFunc[iz][1]->Clone(Form("hMeanMultVsBinWithAdHocFuncOverMC%d",iz+1));
    hMeanMultVsBinWithAdHocFuncOverMC[iz][1]->Divide(hMeanMultVsBin[iz]);

  }

  // draw histograms
  TString cName = corr ? "cDCCorr" : "cDC";
  TCanvas *c = new TCanvas(cName.Data(), cName.Data());
  gPad->SetLogz();
  hNchVsNtrk[nBinZvtx]->Draw("colz");
  pMeanNchVsNtrk[nBinZvtx]->SetLineWidth(2);
  pMeanNchVsNtrk[nBinZvtx]->Draw("sames");
  fitAlpha[nBinZvtx]->SetLineWidth(1);
  fitAlpha[nBinZvtx]->SetLineColor(4);
  fitAlpha[nBinZvtx]->Draw("same");
  fMeanNchVsNtrk[nBinZvtx]->SetLineWidth(1);
  fMeanNchVsNtrk[nBinZvtx]->SetLineColor(4);
  fMeanNchVsNtrk[nBinZvtx]->Draw("same");

  TString cName1 = corr ? "cDCCorrInZvtxBin" : "cDCInZvtxBin";
  TCanvas *c1 = new TCanvas(cName1.Data(), cName1.Data(), 0, 0, 1500, 600);
  c1->Divide((nBinZvtx+1)/2,2);
  for (Int_t iz = 0; iz < nBinZvtx; ++iz) {
    c1->cd(iz+1);
    gPad->SetLogz();
    hNchVsNtrk[iz]->Draw("colz");
    pMeanNchVsNtrk[iz]->SetLineWidth(2);
    pMeanNchVsNtrk[iz]->Draw("sames");
    fitAlpha[iz]->SetLineWidth(1);
    fitAlpha[iz]->SetLineColor(4);
    fitAlpha[iz]->Draw("same");
    fMeanNchVsNtrk[iz]->SetLineWidth(1);
    fMeanNchVsNtrk[iz]->SetLineColor(4);
    fMeanNchVsNtrk[iz]->Draw("same");
  }

  TCanvas *c2 = new TCanvas("cNchVsZvtx", "cNchVsZvtx");
  c2->Divide(3,1);
  c2->cd(1);
  gPad->SetLogz();
  hNchVsZvtx->Draw("colz");
  c2->cd(2);
  pMeanNchVsZvtx[nBins]->Draw("");
  c2->cd(3);
  hRMSNchVsZvtx[nBins]->SetMarkerStyle(20);
  hRMSNchVsZvtx[nBins]->Draw("p");

  TCanvas *c22 = new TCanvas("cMeanNchVsZvtxInNtrkBin", "cMeanNchVsZvtxInNtrkBin", 0, 0, 1500, 600);
  c22->Divide((nBins+1)/2,2);
  for (Int_t i = 0; i < nBins; ++i) {
    c22->cd(i+1);
    pMeanNchVsZvtx[i]->Scale(1./hMult[i]->GetMean());
    pMeanNchVsZvtx[i]->GetYaxis()->SetLabelSize(0.05);
    pMeanNchVsZvtx[i]->GetYaxis()->SetTitleOffset(1.5);
    pMeanNchVsZvtx[i]->Draw("");
    gPad->SetGridy();
  }

  TCanvas *c23 = new TCanvas("cRMSNchVsZvtxInNtrkBin", "cRMSNchVsZvtxInNtrkBin", 0, 0, 1500, 600);
  c23->Divide((nBins+1)/2,2);
  for (Int_t i = 0; i < nBins; ++i) {
    c23->cd(i+1);
    hRMSNchVsZvtx[i]->Scale(1./hMult[i]->GetRMS());
    hRMSNchVsZvtx[i]->GetYaxis()->SetLabelSize(0.05);
    hRMSNchVsZvtx[i]->GetYaxis()->SetTitleOffset(1.5);
    hRMSNchVsZvtx[i]->SetMarkerStyle(20);
    hRMSNchVsZvtx[i]->Draw("p");
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
  hMeanMultVsBin[nBinZvtx]->SetLineColor(1);
  hMeanMultVsBin[nBinZvtx]->Draw();
  hMeanMultVsBinWithAlpha[nBinZvtx][0]->SetLineColor(2);
  hMeanMultVsBinWithAlpha[nBinZvtx][0]->Draw("esame");
  hMeanMultVsBinWithAdHocFunc[nBinZvtx][0]->SetLineColor(4);
  hMeanMultVsBinWithAdHocFunc[nBinZvtx][0]->Draw("esame");
  c4->cd(2);
  hMeanMultVsBinWithAlphaOverMC[nBinZvtx][0]->SetLineColor(2);
  hMeanMultVsBinWithAlphaOverMC[nBinZvtx][0]->GetYaxis()->SetLabelSize(0.05);
  hMeanMultVsBinWithAlphaOverMC[nBinZvtx][0]->Draw();
  hMeanMultVsBinWithAdHocFuncOverMC[nBinZvtx][0]->SetLineColor(4);
  hMeanMultVsBinWithAdHocFuncOverMC[nBinZvtx][0]->Draw("same");

  TCanvas *c42 = new TCanvas("cMeanMultVsBinFromZvtxBin", "cMeanMultVsBinFromZvtxBin");
  c42->Divide(1,2);
  c42->cd(1);
  hMeanMultVsBin[nBinZvtx]->SetLineColor(1);
  hMeanMultVsBin[nBinZvtx]->Draw();
  for (Int_t iz = 0; iz < nBinZvtx; ++iz) {
    c42->cd(1);
    hMeanMultVsBinWithAlpha[iz][0]->SetLineColor(2);
    hMeanMultVsBinWithAlpha[iz][0]->Draw("esame");
    hMeanMultVsBinWithAdHocFunc[iz][0]->SetLineColor(4);
    hMeanMultVsBinWithAdHocFunc[iz][0]->Draw("esame");
    c42->cd(2);
    hMeanMultVsBinWithAlphaOverMC[iz][0]->SetLineColor(2);
    hMeanMultVsBinWithAlphaOverMC[iz][0]->GetYaxis()->SetLabelSize(0.05);
    hMeanMultVsBinWithAlphaOverMC[iz][0]->Draw((iz == 0) ? "" : "same");
    hMeanMultVsBinWithAdHocFuncOverMC[iz][0]->SetLineColor(4);
    hMeanMultVsBinWithAdHocFuncOverMC[iz][0]->Draw("same");
  }

  TCanvas *c43 = new TCanvas("cMeanMultVsBinInZvtxBin", "cMeanMultVsBinInZvtxBin");
  c43->Divide(1,2);
  for (Int_t iz = 0; iz < nBinZvtx; ++iz) {
    c43->cd(1);
    hMeanMultVsBin[iz]->SetLineColor(1);
    hMeanMultVsBin[iz]->Draw((iz == 0) ? "" : "same");
    hMeanMultVsBinWithAlpha[iz][1]->SetLineColor(2);
    hMeanMultVsBinWithAlpha[iz][1]->Draw("esame");
    hMeanMultVsBinWithAdHocFunc[iz][1]->SetLineColor(4);
    hMeanMultVsBinWithAdHocFunc[iz][1]->Draw("esame");
    c43->cd(2);
    hMeanMultVsBinWithAlphaOverMC[iz][1]->SetLineColor(2);
    hMeanMultVsBinWithAlphaOverMC[iz][1]->GetYaxis()->SetLabelSize(0.05);
    hMeanMultVsBinWithAlphaOverMC[iz][1]->Draw((iz == 0) ? "" : "same");
    hMeanMultVsBinWithAdHocFuncOverMC[iz][1]->SetLineColor(4);
    hMeanMultVsBinWithAdHocFuncOverMC[iz][1]->Draw("same");
  }

  TCanvas *c7 = new TCanvas("cRelativeTrkOverRelativeMeanNch", "cRelativeTrkOverRelativeMeanNch");
  Int_t color = 2;
  for (Int_t iz = 0; iz < nBinZvtx; ++iz) {
    hRelativeTrkOverRelativeMeanNch[iz]->SetLineColor(color);
    hRelativeTrkOverRelativeMeanNch[iz]->Draw((iz == 0) ? "" : "same");
    fRelativeTrkOverRelativeMeanNchWithAlpha[iz]->SetLineWidth(1);
    fRelativeTrkOverRelativeMeanNchWithAlpha[iz]->SetLineColor(color);
    fRelativeTrkOverRelativeMeanNchWithAlpha[iz]->Draw("same");
    fRelativeTrkOverRelativeMeanNchWithAdHocFunc[iz]->SetLineWidth(1);
    fRelativeTrkOverRelativeMeanNchWithAdHocFunc[iz]->SetLineColor(color);
    fRelativeTrkOverRelativeMeanNchWithAdHocFunc[iz]->Draw("same");
    ++color;
    if (color == 5 || color == 10) ++color;
  }
  hRelativeTrkOverRelativeMeanNch[nBinZvtx]->SetLineColor(1);
  hRelativeTrkOverRelativeMeanNch[nBinZvtx]->SetLineWidth(2);
  hRelativeTrkOverRelativeMeanNch[nBinZvtx]->Draw("same");
  fRelativeTrkOverRelativeMeanNchWithAlpha[nBinZvtx]->SetLineWidth(2);
  fRelativeTrkOverRelativeMeanNchWithAlpha[nBinZvtx]->SetLineColor(1);
  fRelativeTrkOverRelativeMeanNchWithAlpha[nBinZvtx]->Draw("same");
  fRelativeTrkOverRelativeMeanNchWithAdHocFunc[nBinZvtx]->SetLineWidth(2);
  fRelativeTrkOverRelativeMeanNchWithAdHocFunc[nBinZvtx]->SetLineColor(1);
  fRelativeTrkOverRelativeMeanNchWithAdHocFunc[nBinZvtx]->Draw("same");

  TCanvas *c9 = new TCanvas("cRelativeMeanTrkOverRelativeNch", "cRelativeMeanTrkOverRelativeNch");
  hRelativeMeanTrkOverRelativeNch->SetLineColor(1);
  hRelativeMeanTrkOverRelativeNch->SetLineWidth(2);
  hRelativeMeanTrkOverRelativeNch->Draw();

  DrawMeanTrackletsVsMeanNChFlat(hNchVsNtrk[nBinZvtx], corr, fitAlpha[nBinZvtx], fMeanNchVsNtrk[nBinZvtx], rfMeanNchVsNtrk[nBinZvtx]);

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
Double_t RelativeTrkOverRelativeMeanNchWithAlpha(Double_t *x, Double_t *par)
{
  // function to compute (Ntrk^i/<Ntrk>) / (<Nch>^i/<Nch>)
  // with MeanNchVsNtrk fitted with y = alpha * x

  Double_t nTrk = *x;
  Double_t meanNch = par[0] * nTrk;
  Double_t meanNtrkInt = par[1];
  Double_t meanNchInt = par[2];

  return nTrk/meanNtrkInt*meanNchInt/meanNch;
}

//---------------------------------------------------------------------------------------------
Double_t RelativeTrkOverRelativeMeanNchWithAdHocFunc(Double_t *x, Double_t *par)
{
  // function to compute (Ntrk^i/<Ntrk>) / (<Nch>^i/<Nch>)
  // with MeanNchVsNtrk fitted with ad-hoc function

  Double_t nTrk = *x;
  Double_t meanNch = MeanNchVsNtrk(x,par);
  Double_t meanNtrkInt = par[5];
  Double_t meanNchInt = par[6];

  return nTrk/meanNtrkInt*meanNchInt/meanNch;
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
  meanMultErrInt = meanTrkErr * meanMult / meanTrk;
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

//---------------------------------------------------------------------------------------------
void DrawMeanTrackletsVsMeanNChFlat(TH2D *hNchVsNtrk, Bool_t corr, TF1 *fitAlpha, TF1 *fMeanNchVsNtrk, TFitResultPtr &rfMeanNchVsNtrk)
{
  /// Study the Nch vs Ntrk correlation in case of a flat Nch distribution

  // reweight the Nch vs Ntrk 2D-correlation to correspond to a flat Nch distribution
  TH1D *hMultTot = hNchVsNtrk->ProjectionY("hMultTot", 0, -1, "e");
  TH2D *hNchVsNtrkFlat = (TH2D*)hNchVsNtrk->Clone("hNchVsNtrkFlat");
  for (Int_t j = 1; j <= hNchVsNtrkFlat->GetNbinsY(); ++j) {
    Double_t nEv = hMultTot->GetBinContent(j);
    if (nEv > 0.) for (Int_t i = 1; i <= hNchVsNtrkFlat->GetNbinsX(); ++i) {
      hNchVsNtrkFlat->SetBinContent(i,j,hNchVsNtrkFlat->GetBinContent(i,j)/nEv);
      hNchVsNtrkFlat->SetBinError(i,j,hNchVsNtrkFlat->GetBinError(i,j)/nEv);
    }
  }
  delete hMultTot;

  // <nCh> vs Ntrk profile from reweighted correlation
  TProfile *pMeanNchVsNtrkFlat = hNchVsNtrkFlat->ProfileX(corr ? "pMeanNchVsNtrkFlatCorr" : "pMeanNchVsNtrkFlat");

  // fit reweighted Nch vs Ntrk with y = ax
  TF1 *fitAlphaFlat = new TF1("fitAlphaFlat","[0]*x",1.,100.);
  pMeanNchVsNtrkFlat->Fit(fitAlphaFlat,"NR");

  // fit reweighted <Nch> vs Ntrk profile with an ad-hoc polynomial function
  TF1 *fMeanNchVsNtrkFlat = new TF1("fMeanNchVsNtrkFlat",MeanNchVsNtrk,1.,100.,5);
  fMeanNchVsNtrkFlat->SetParameters(1., 1., 1., 1., 10.);
  fMeanNchVsNtrkFlat->SetParLimits(4, 2., 40.);
  TFitResultPtr rfMeanNchVsNtrkFlat = pMeanNchVsNtrkFlat->Fit(fMeanNchVsNtrkFlat,"NRS");

  // reweigthed tracklet distribution and integrated <Ntrk>
  Double_t maxNtrkFlat = 30.; // beyond this point we don't have enought stat to get reliable results
  TH1D *hTrkFlat = hNchVsNtrkFlat->ProjectionX("hTrkFlat", 0, -1, "e");
  hTrkFlat->SetAxisRange(1, maxNtrkFlat);
  Double_t meanNtrkIntFlat = hTrkFlat->GetMean();
  printf(corr?"integrated <NtrkCorr> = %f\n":"<Ntrk> = %f\n",meanNtrkIntFlat);
  hTrkFlat->GetXaxis()->SetRange();

  // integrated <Nch> from reweigthed distribution
  TH1D *hTmp = hNchVsNtrkFlat->ProjectionY("hTmp", hTrkFlat->FindBin(1), hTrkFlat->FindBin(maxNtrkFlat), "e");
  Double_t meanNchIntFlat = hTmp->GetMean();
  printf(corr?"integrated <Nch> = %f\n":"<Nch> = %f\n",meanNchIntFlat);
  delete hTmp;
  Double_t err;
  Double_t meanNchIntFlatWithAlpha = GetMeanMultWithAlpha(hTrkFlat, 1., maxNtrkFlat, fitAlphaFlat, err);
  Double_t meanNchIntFlatWithAdHocFunc = GetMeanMultWithAdHocFunc(hTrkFlat, 1., maxNtrkFlat, fMeanNchVsNtrkFlat, rfMeanNchVsNtrkFlat, err);

  // relative Ntrk over relative <Nch> ((Ntrk^i/<Ntrk>) / (<Nch>^i/<Nch>)) from reweigthed correlated
  TH1D *hRelativeTrkOverRelativeMeanNchFlat = new TH1D(Form("hRelativeTrk%sOverRelativeMeanNchFlat",corr?"Corr":""), Form("relative Ntrk%s over relative <Nch>;Ntrk^{%s};(Ntrk%s^{i}/<Ntrk%s>)/(<Nch>^{i}/<Nch>)",corr?"corr":"",corr?"Corr":"",corr?"Corr":"",corr?"Corr":""), hTrkFlat->GetNbinsX(), hTrkFlat->GetXaxis()->GetXmin(), hTrkFlat->GetXaxis()->GetXmax());
  for (Int_t i = 1; i <= hTrkFlat->GetNbinsX(); ++i) {
    if (pMeanNchVsNtrkFlat->GetBinContent(i) > 0. && pMeanNchVsNtrkFlat->GetBinError(i) > 0.) {
      Double_t relativeTrkOverRelativeMeanNch = hTrkFlat->GetBinCenter(i)/meanNtrkIntFlat*meanNchIntFlat/pMeanNchVsNtrkFlat->GetBinContent(i);
      hRelativeTrkOverRelativeMeanNchFlat->SetBinContent(i,relativeTrkOverRelativeMeanNch);
      hRelativeTrkOverRelativeMeanNchFlat->SetBinError(i,pMeanNchVsNtrkFlat->GetBinError(i)/pMeanNchVsNtrkFlat->GetBinContent(i)*relativeTrkOverRelativeMeanNch);
    }
  }
  TF1 *fRelativeTrkOverRelativeMeanNchFlatWithAlpha = new TF1("fRelativeTrkOverRelativeMeanNchFlatWithAlpha", RelativeTrkOverRelativeMeanNchWithAlpha,1.,100.,3);
  fRelativeTrkOverRelativeMeanNchFlatWithAlpha->SetParameters(fitAlphaFlat->GetParameter(0), meanNtrkIntFlat, meanNchIntFlatWithAlpha);
  TF1 *fRelativeTrkOverRelativeMeanNchFlatWithAdHocFunc = new TF1("fRelativeTrkOverRelativeMeanNchFlatWithAdHocFunc", RelativeTrkOverRelativeMeanNchWithAdHocFunc,1.,100.,7);
  fRelativeTrkOverRelativeMeanNchFlatWithAdHocFunc->SetParameters(fMeanNchVsNtrkFlat->GetParameter(0), fMeanNchVsNtrkFlat->GetParameter(1), fMeanNchVsNtrkFlat->GetParameter(2), fMeanNchVsNtrkFlat->GetParameter(3), fMeanNchVsNtrkFlat->GetParameter(4), meanNtrkIntFlat, meanNchIntFlatWithAdHocFunc);

  // reweighted multiplicity distribution and mean multiplicity per tracklet bin
  Int_t nBins = (Int_t)(sizeof(trkBin) / sizeof(Int_t)) - 1;
  TH1D **hMultFlat = new TH1D*[nBins+1];
  hMultFlat[nBins] = hNchVsNtrkFlat->ProjectionY("hMultFlat", hTrkFlat->FindBin(0), -1, "e");
  hMultFlat[nBins]->SetTitle(Form("Nch distribution in %sNtrk bins;Nch",corr?"corrected ":""));
  TH1D *hMeanMultFlatVsBin = new TH1D("hMeanMultFlatVsBin", Form("<Nch> versus %sNtrk bin;%sNtrk bin;<Nch>",corr?"corrected ":"",corr?"corrected ":""), nBins, 0.5, nBins+0.5);;
  TH1D *hMeanMultFlatVsBinWithAlpha[2];
  TH1D *hMeanMultFlatVsBinWithAdHocFunc[2];
  hMeanMultFlatVsBinWithAlpha[0] = new TH1D("hMeanMultFlatVsBinWithAlphaFlat", "hMeanMultFlatVsBinWithAlphaFlat", nBins, 0.5, nBins+0.5);
  hMeanMultFlatVsBinWithAdHocFunc[0] = new TH1D("hMeanMultFlatVsBinWithAdHocFuncFlat", "hMeanMultFlatVsBinWithAdHocFuncFlat", nBins, 0.5, nBins+0.5);
  hMeanMultFlatVsBinWithAlpha[1] = new TH1D("hMeanMultFlatVsBinWithAlpha", "hMeanMultFlatVsBinWithAlpha", nBins, 0.5, nBins+0.5);
  hMeanMultFlatVsBinWithAdHocFunc[1] = new TH1D("hMeanMultFlatVsBinWithAdHocFunc", "hMeanMultFlatVsBinWithAdHocFunc", nBins, 0.5, nBins+0.5);
  for (Int_t i = 0; i < nBins; ++i) {

    hMultFlat[i] = hNchVsNtrkFlat->ProjectionY(Form("hMultFlat%d",i+1), hTrkFlat->FindBin(trkBin[i]), hTrkFlat->FindBin(trkBin[i+1]-1), "e");

    // mean Ntrk in zVtx bin "iz" and Ntrk bin "i"
    hTrkFlat->SetAxisRange(trkBin[i], trkBin[i+1]-1);
    printf(corr?"<NtrkCorr> = %f":"<Ntrk> = %f",hTrkFlat->GetMean());
    hTrkFlat->GetXaxis()->SetRange();

    // mean multiplicity in zVtx bin "iz" and Ntrk bin "i" from MC truth
    hMeanMultFlatVsBin->SetBinContent(i+1, hMultFlat[i]->GetMean());
    hMeanMultFlatVsBin->SetBinError(i+1, hMultFlat[i]->GetMeanError());
    printf(" <Nch>(Calc) = %f",hMultFlat[i]->GetMean());

    // mean multiplicity in Ntrk bin "i" computed with constant alpha factor from the weighted correlation
    Double_t meanMultErr = 0.;
    Double_t meanMult = GetMeanMultWithAlpha(hTrkFlat, trkBin[i], trkBin[i+1]-1, fitAlphaFlat, meanMultErr);
    printf(" (%f)",meanMult);
    hMeanMultFlatVsBinWithAlpha[0]->SetBinContent(i+1, meanMult);
    hMeanMultFlatVsBinWithAlpha[0]->SetBinError(i+1, meanMultErr);

    // mean multiplicity in Ntrk bin "i" computed with constant alpha factor from original correlation
    meanMultErr = 0.;
    meanMult = GetMeanMultWithAlpha(hTrkFlat, trkBin[i], trkBin[i+1]-1, fitAlpha, meanMultErr);
    printf(" (%f)",meanMult);
    hMeanMultFlatVsBinWithAlpha[1]->SetBinContent(i+1, meanMult);
    hMeanMultFlatVsBinWithAlpha[1]->SetBinError(i+1, meanMultErr);

    // mean multiplicity in Ntrk bin "i" computed with the ad-hoc fit of the weighted <Nch> vs Ntrk profile
    meanMult = GetMeanMultWithAdHocFunc(hTrkFlat, trkBin[i], trkBin[i+1]-1, fMeanNchVsNtrkFlat, rfMeanNchVsNtrkFlat, meanMultErr);
    printf(" (%f)",meanMult);
    hMeanMultFlatVsBinWithAdHocFunc[0]->SetBinContent(i+1, meanMult);
    hMeanMultFlatVsBinWithAdHocFunc[0]->SetBinError(i+1, meanMultErr);

    // mean multiplicity in Ntrk bin "i" computed with the ad-hoc fit of the original <Nch> vs Ntrk profile
    meanMult = GetMeanMultWithAdHocFunc(hTrkFlat, trkBin[i], trkBin[i+1]-1, fMeanNchVsNtrk, rfMeanNchVsNtrk, meanMultErr);
    printf(" (%f)\n",meanMult);
    hMeanMultFlatVsBinWithAdHocFunc[1]->SetBinContent(i+1, meanMult);
    hMeanMultFlatVsBinWithAdHocFunc[1]->SetBinError(i+1, meanMultErr);

  }

  // ratio of mean multiplicity from reweighted distribution
  TH1D *hMeanMultFlatVsBinWithAlphaOverMC[2];
  TH1D *hMeanMultFlatVsBinWithAdHocFuncOverMC[2];
  hMeanMultFlatVsBinWithAlphaOverMC[0] = (TH1D*)hMeanMultFlatVsBinWithAlpha[0]->Clone("hMeanMultFlatVsBinWithAlphaFlatOverMC");
  hMeanMultFlatVsBinWithAlphaOverMC[0]->Divide(hMeanMultFlatVsBin);
  hMeanMultFlatVsBinWithAlphaOverMC[0]->SetTitle(Form("calculated <Nch> over MC truth versus %sNtrk bin;%sNtrk bin;<Nch>_{calc}/<Nch>_{MC}",corr?"corrected ":"",corr?"corrected ":""));
  hMeanMultFlatVsBinWithAdHocFuncOverMC[0] = (TH1D*)hMeanMultFlatVsBinWithAdHocFunc[0]->Clone("hMeanMultFlatVsBinWithAdHocFuncFlatOverMC");
  hMeanMultFlatVsBinWithAdHocFuncOverMC[0]->Divide(hMeanMultFlatVsBin);
  hMeanMultFlatVsBinWithAlphaOverMC[1] = (TH1D*)hMeanMultFlatVsBinWithAlpha[1]->Clone("hMeanMultFlatVsBinWithAlphaOverMC");
  hMeanMultFlatVsBinWithAlphaOverMC[1]->Divide(hMeanMultFlatVsBin);
  hMeanMultFlatVsBinWithAlphaOverMC[1]->SetTitle(Form("calculated <Nch> over MC truth versus %sNtrk bin;%sNtrk bin;<Nch>_{calc}/<Nch>_{MC}",corr?"corrected ":"",corr?"corrected ":""));
  hMeanMultFlatVsBinWithAdHocFuncOverMC[1] = (TH1D*)hMeanMultFlatVsBinWithAdHocFunc[1]->Clone("hMeanMultFlatVsBinWithAdHocFuncOverMC");
  hMeanMultFlatVsBinWithAdHocFuncOverMC[1]->Divide(hMeanMultFlatVsBin);

  // draw histograms
  TCanvas *c5 = new TCanvas("cMultFlat", "cMultFlat");
  gPad->SetLogy();
  hMultFlat[nBins]->SetLineColor(1);
  hMultFlat[nBins]->Draw();
  for (Int_t i = 0; i < nBins; ++i) {
    hMultFlat[i]->SetLineColor(i+2);
    hMultFlat[i]->Draw("sames");
  }

  TCanvas *c51 = new TCanvas("cNchVsNtrkFlat", "cNchVsNtrkFlat");
  gPad->SetLogz();
  hNchVsNtrkFlat->Draw("colz");
  pMeanNchVsNtrkFlat->SetLineWidth(2);
  pMeanNchVsNtrkFlat->Draw("sames");
  fitAlphaFlat->SetLineWidth(1);
  fitAlphaFlat->SetLineColor(4);
  fitAlphaFlat->Draw("same");
  fMeanNchVsNtrkFlat->SetLineWidth(1);
  fMeanNchVsNtrkFlat->SetLineColor(4);
  fMeanNchVsNtrkFlat->Draw("same");

  TCanvas *c6 = new TCanvas("cMeanMultFlatVsBin", "cMeanMultFlatVsBin");
  c6->Divide(1,2);
  c6->cd(1);
  hMeanMultFlatVsBin->SetLineColor(1);
  hMeanMultFlatVsBin->Draw();
  hMeanMultFlatVsBinWithAlpha[0]->SetLineColor(2);
  hMeanMultFlatVsBinWithAlpha[0]->Draw("esame");
  hMeanMultFlatVsBinWithAdHocFunc[0]->SetLineColor(4);
  hMeanMultFlatVsBinWithAdHocFunc[0]->Draw("esame");
  c6->cd(2);
  hMeanMultFlatVsBinWithAlphaOverMC[0]->SetLineColor(2);
  hMeanMultFlatVsBinWithAlphaOverMC[0]->GetYaxis()->SetLabelSize(0.05);
  hMeanMultFlatVsBinWithAlphaOverMC[0]->Draw();
  hMeanMultFlatVsBinWithAdHocFuncOverMC[0]->SetLineColor(4);
  hMeanMultFlatVsBinWithAdHocFuncOverMC[0]->Draw("same");

  TCanvas *c62 = new TCanvas("cMeanMultFlatVsBinFromOriginal", "cMeanMultFlatVsBinFromOriginal");
  c62->Divide(1,2);
  c62->cd(1);
  hMeanMultFlatVsBin->SetLineColor(1);
  hMeanMultFlatVsBin->Draw();
  hMeanMultFlatVsBinWithAlpha[1]->SetLineColor(2);
  hMeanMultFlatVsBinWithAlpha[1]->Draw("esame");
  hMeanMultFlatVsBinWithAdHocFunc[1]->SetLineColor(4);
  hMeanMultFlatVsBinWithAdHocFunc[1]->Draw("esame");
  c62->cd(2);
  hMeanMultFlatVsBinWithAlphaOverMC[1]->SetLineColor(2);
  hMeanMultFlatVsBinWithAlphaOverMC[1]->GetYaxis()->SetLabelSize(0.05);
  hMeanMultFlatVsBinWithAlphaOverMC[1]->Draw();
  hMeanMultFlatVsBinWithAdHocFuncOverMC[1]->SetLineColor(4);
  hMeanMultFlatVsBinWithAdHocFuncOverMC[1]->Draw("same");

  TCanvas *c8 = new TCanvas("cRelativeTrkOverRelativeMeanNchFlat", "cRelativeTrkOverRelativeMeanNchFlat");
  hRelativeTrkOverRelativeMeanNchFlat->SetLineColor(1);
  hRelativeTrkOverRelativeMeanNchFlat->SetLineWidth(2);
  hRelativeTrkOverRelativeMeanNchFlat->Draw();
  fRelativeTrkOverRelativeMeanNchFlatWithAlpha->SetLineWidth(2);
  fRelativeTrkOverRelativeMeanNchFlatWithAlpha->SetLineColor(1);
  fRelativeTrkOverRelativeMeanNchFlatWithAlpha->Draw("same");
  fRelativeTrkOverRelativeMeanNchFlatWithAdHocFunc->SetLineWidth(2);
  fRelativeTrkOverRelativeMeanNchFlatWithAdHocFunc->SetLineColor(1);
  fRelativeTrkOverRelativeMeanNchFlatWithAdHocFunc->Draw("same");

}

