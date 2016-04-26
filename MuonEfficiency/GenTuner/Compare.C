/*
 *  Compare.C
 *  aliroot_dev
 *
 *  Created by philippe pillot on 22/03/13.
 *  Copyright 2013 SUBATECH. All rights reserved.
 *
 */


TF1 *fRes[2][3] = {{0x0,0x0,0x0},{0x0,0x0,0x0}};

//______________________________________________________________________________
Double_t PtRat(const Double_t *x, const Double_t */*p*/)
{
  /// generated pT fit function ratio
  return (fRes[0][1] && fRes[0][0]) ? fRes[0][1]->Eval(*x) / fRes[0][0]->Eval(*x) : 0.;
}

//______________________________________________________________________________
Double_t YRat(const Double_t *x, const Double_t */*p*/)
{
  /// generated y fit function ratio
  return (fRes[1][1] && fRes[1][0]) ? fRes[1][1]->Eval(*x) / fRes[1][0]->Eval(*x) : 0.;
}

//______________________________________________________________________________
Double_t GetLowEdge(TH1 &h)
{
  /// adjust the lower edge of the fit range according to the content of the histogram
  Int_t binAbove0 = h.FindFirstBinAbove(0.);
  if (h.GetBinContent(binAbove0) < 0.1*h.GetBinContent(binAbove0+1)) binAbove0++;
  return h.GetBinLowEdge(binAbove0);
}

//______________________________________________________________________________
Double_t GetUpEdge(TH1 &h)
{
  /// adjust the upper edge of the fit range according to the content of the histogram
  Int_t binAbove0 = h.FindLastBinAbove(0.);
  if (h.GetBinContent(binAbove0) < 0.1*h.GetBinContent(binAbove0-1)) binAbove0--;
  return h.GetBinLowEdge(binAbove0+1);
}

//______________________________________________________________________________
void Compare(TString sfile1, TString sfile2)
{
  /// compare generated histograms and functions
  
  TString sfile[2];
  sfile[0] = sfile1.Data();
  sfile[1] = sfile2.Data();
  TString sRes[8] = {"hPtGen", "hYGen", "hPhiGen", "hSignGen", "hPtRec", "hYRec", "hPhiRec", "hSignRec"};
  TH1 *hRes[8][3];
  for (Int_t i = 0; i < 8; i++) for (Int_t j = 0; j < 3; j++) hRes[i][j] = 0x0;
//  TString sfunc[2][2] = {"fPtFuncMC", "fPtFuncMC", "fYFuncMC", "fYFuncMC"};
  TString sfunc[2][2] = {"fPtFuncNew", "fPtFuncNew", "fYFuncNew", "fYFuncNew"};
  void *func[2] = {PtRat, YRat};
  TParameter<Double_t> *newMuPlusFrac[2] = {0x0, 0x0};
  Double_t nMu[2][2] = {{0.,0.},{0.,0.}};
  
  // get results
  for (Int_t j = 0; j < 2; j++) {
    TFile *file = TFile::Open(sfile[j].Data(),"READ");
    if (!file || !file->IsOpen()) {
      ::Error("cannot open file");
      return;
    }
    if (file && file->IsOpen()) {
      TObjArray *list = static_cast<TObjArray*>(file->FindObjectAny("Histograms"));
      if (!list) {
	::Error("cannot find histograms");
	return;
      }
      for (Int_t i = 0; i < 8; i++) {
	hRes[i][j] = static_cast<TH1*>(list->FindObject(sRes[i].Data())->Clone(Form("%s%d",sRes[i].Data(),j+1)));
	if (hRes[i][j]) hRes[i][j]->SetDirectory(0);
//        if (hRes[i][j]) hRes[i][j]->Rebin(24);
      }
      for (Int_t i = 0; i < 2; i++) {
	TF1 *f = static_cast<TF1*>(file->FindObjectAny(sfunc[i][j].Data()));
	if (f) {
          fRes[i][j] = new TF1(*f);
          fRes[i][j]->SetName(Form("%s%d",f->GetName(),j+1));
	}
      }
      newMuPlusFrac[j] = static_cast<TParameter<Double_t>*>(file->FindObjectAny("newMuPlusFrac"));
      if (hRes[1][j]) nMu[j][0] = hRes[1][j]->Integral();
      if (hRes[5][j]) nMu[j][1] = hRes[5][j]->Integral();
    }
    file->Close();
  }
  
  // normalize histograms
  for (Int_t i = 0; i < 8; i++) {
    if (hRes[i][0] && hRes[i][1]) {
      for (Int_t j = 0; j < 2; j++) {
        Double_t integral;
        if (i < 2 && fRes[i][j]) integral = hRes[i][j]->Integral(hRes[i][j]->FindBin(fRes[i][j]->GetXmin()), hRes[i][j]->FindBin(fRes[i][j]->GetXmax()), "width");
        else if ((i == 0 && !fRes[i][0] && !fRes[i][1]) || i == 4) {
          Double_t min = TMath::Max(GetLowEdge(*hRes[i][0]), GetLowEdge(*hRes[i][1]));
          Double_t max = TMath::Min(GetUpEdge(*hRes[i][0]), GetUpEdge(*hRes[i][1]));
          integral = hRes[i][j]->Integral(hRes[i][j]->FindBin(min), hRes[i][j]->FindBin(max), "width");
        } else integral = hRes[i][j]->Integral("width");
        Double_t norm = (integral != 0.) ? 1./integral : 1.;
        hRes[i][j]->Scale(norm);
      }
    }
  }
  
  // compute ratios
  for (Int_t i = 0; i < 8; i++) {
    if (hRes[i][0] && hRes[i][1]) {
      hRes[i][2] = static_cast<TH1*>(hRes[i][1]->Clone(Form("%sOver1",hRes[i][1]->GetName())));
      hRes[i][2]->Divide(hRes[i][0]);
    }
  }
  for (Int_t i = 0; i < 2 && fRes[i][0] && fRes[i][1]; i++) {
    fRes[i][2] = new TF1(Form("%sOver1",fRes[i][1]->GetName()), func[i], fRes[i][1]->GetXmin(), fRes[i][1]->GetXmax(), 0);
  }
  
  // print results
  TString var[2] = {"pT", "y"};
  for (Int_t i = 0; i < 2; i++) {
    for (Int_t j = 0; j < 2; j++) {
      if (fRes[i][j]) {
        Double_t *param = fRes[i][j]->GetParameters();
        printf("\n%s parameters for single muon generator in file %d:\n", var[i].Data(), j+1);
        printf("Double_t p[%d] = {", fRes[i][j]->GetNpar());
        for (Int_t k = 0; k < fRes[i][j]->GetNpar()-1; k++) printf("%g, ", fRes[i][j]->GetParameter(k));
        printf("%g};\n", fRes[i][j]->GetParameter(fRes[i][j]->GetNpar()-1));
      }
    }
  }
  for (Int_t j = 0; j < 2; j++) {
    if (newMuPlusFrac[j]) {
      printf("\nfraction of mu+ for single muon generator in file %d:\n", j+1);
      printf("Double_t newMuPlusFrac = %f\n", newMuPlusFrac[j]->GetVal());
    }
  }
  for (Int_t j = 0; j < 2; j++) {
    printf("\nnumber of muons in file %d (full gen/reco range):\n", j+1);
    printf("generated: %g ; reconstructed: %g", nMu[j][0], nMu[j][1]);
    if (nMu[j][0] > 0.) printf(" --> ratio = %g", nMu[j][1]/nMu[j][0]);
    printf("\n");
  }
  
  // display results at the generation level
  TCanvas *cGen = new TCanvas("cGen", "cGen", 1200, 600);
  cGen->Divide(4,2);
  for (Int_t i = 0; i < 4; i++) {
    cGen->cd(i+1);
    if (i == 0) gPad->SetLogy();
    for (Int_t j = 0; j < 2 && hRes[i][j]; j++) {
      hRes[i][j]->SetLineColor(2*j+2);
      hRes[i][j]->Draw((j == 0) ? "" : "sames");
      if (i < 2 && fRes[i][j]) {
	fRes[i][j]->SetLineColor(2*j+2);
	fRes[i][j]->Draw("sames");
      }
    }
    cGen->cd(i+5);
    if (hRes[i][2]) {
      hRes[i][2]->Draw();
      hRes[i][2]->SetMinimum(0.8);
      hRes[i][2]->SetMaximum(1.2);
    }
    if (i < 2 && fRes[i][2]) fRes[i][2]->Draw("sames");
  }
  
  // display results at the reconstruction level
  TLegend *lRes = new TLegend(0.5,0.55,0.85,0.75);
  TCanvas *cRec = new TCanvas("cRec", "cRec", 1200, 600);
  cRec->Divide(4,2);
  cRec->cd(1);
  for (Int_t i = 4; i < 8; i++) {
    cRec->cd(i-3);
    if (i == 4) gPad->SetLogy();
    for (Int_t j = 0; j < 2 && hRes[i][j]; j++) {
      if (i == 4) lRes->AddEntry(hRes[i][j],Form("file%d",j+1),"l");
      hRes[i][j]->SetLineColor(2*j+2);
      hRes[i][j]->Draw((j == 0) ? "" : "sames");
    }
    if (i == 4) lRes->Draw("same");
    cRec->cd(i+1);
    if (hRes[i][2]) {
      hRes[i][2]->Draw();
      hRes[i][2]->SetMinimum(0.8);
      hRes[i][2]->SetMaximum(1.2);
    }
  }
  
}

