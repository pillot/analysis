/*
 *  Compare.C
 *  aliroot_dev
 *
 *  Created by philippe pillot on 22/03/13.
 *  Copyright 2013 SUBATECH. All rights reserved.
 *
 */


//______________________________________________________________________________
Double_t Pt(const Double_t *x, const Double_t *p)
{
  /// generated pT fit function
  Double_t pT = *x;
  return p[0] * (1. / TMath::Power(p[1] + TMath::Power(pT,p[2]), p[3]) + p[4] * TMath::Exp(p[5]*pT));
}

//______________________________________________________________________________
Double_t Y(const Double_t *x, const Double_t *p)
{
  /// generated y fit function
  Double_t y = *x;
  Double_t arg = y/p[7];
  return p[0] * (p[1] * (1. + p[2]*y + p[3]*y*y + p[4]*y*y*y + p[5]*y*y*y*y) + p[6]*TMath::Exp(-0.5*arg*arg));
}

//______________________________________________________________________________
Double_t PtRat(const Double_t *x, const Double_t *p)
{
  /// generated pT fit function ratio
  const Double_t *p1 = &(p[0]);
  const Double_t *p2 = &(p[6]);
  return Pt(x,p1) / Pt(x,p2);
}

//______________________________________________________________________________
Double_t YRat(const Double_t *x, const Double_t *p)
{
  /// generated y fit function ratio
  const Double_t *p1 = &(p[0]);
  const Double_t *p2 = &(p[8]);
  return Y(x,p1) / Y(x,p2);
}

//______________________________________________________________________________
void Compare(TString sfile1, TString sfile2)
{
  /// compare generated histograms and functions
  
  TString sfile[2];
  sfile[0] = sfile1.Data();
  sfile[1] = sfile2.Data();
  TString sRes[6] = {"hPtGen", "hYGen", "hPhiGen", "hPtRec", "hYRec", "hPhiRec"};
  TH1 *hRes[6][3];
  for (Int_t i = 0; i < 6; i++) for (Int_t j = 0; j < 3; j++) hRes[i][j] = 0x0;
  TString sfunc[2][2] = {"fPtFuncMC", "fPtFuncMC", "fYFuncMC", "fYFuncMC"};
  void *func[4] = {Pt, Y, PtRat, YRat};
  Int_t nPar[4] = {6, 8, 12, 16};
  TF1 *fRes[2][3];
  for (Int_t i = 0; i < 2; i++) for (Int_t j = 0; j < 3; j++) fRes[i][j] = 0x0;
  
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
      for (Int_t i = 0; i < 6; i++) {
	hRes[i][j] = static_cast<TH1*>(list->FindObject(sRes[i].Data())->Clone(Form("%s%d",sRes[i].Data(),j+1)));
	if (hRes[i][j]) hRes[i][j]->SetDirectory(0);
      }
      for (Int_t i = 0; i < 2; i++) {
	TF1 *f = static_cast<TF1*>(file->FindObjectAny(sfunc[i][j].Data()));
	if (f) {
	  fRes[i][j] = new TF1(Form("%s%d",sfunc[i][j].Data(),j+1), func[i], f->GetXmin(), f->GetXmax(), nPar[i]);
	  fRes[i][j]->SetParameters(f->GetParameters());
	}
      }
    }
    file->Close();
  }
  
  // normalize histograms
  for (Int_t i = 0; i < 6 && hRes[i][0] && hRes[i][1]; i++) {
    for (Int_t j = 0; j < 2; j++) {
      Double_t integral;
      if (i == 0) integral = hRes[i][j]->Integral(hRes[i][j]->FindBin(0.8), hRes[i][j]->FindBin(30.), "width");
      else if (i == 3) integral = hRes[i][j]->Integral(hRes[i][j]->FindBin(4.), hRes[i][j]->FindBin(30.), "width");
      else integral = hRes[i][j]->Integral("width");
      Double_t norm = (integral != 0.) ? 1./integral : 1.;
      hRes[i][j]->Scale(norm);
    }
  }
  
  // compute ratios
  for (Int_t i = 0; i < 6 && hRes[i][0] && hRes[i][1]; i++) {
    hRes[i][2] = static_cast<TH1*>(hRes[i][1]->Clone(Form("%sOver1",hRes[i][1]->GetName())));
    hRes[i][2]->Divide(hRes[i][0]);
  }
  for (Int_t i = 0; i < 2 && fRes[i][0] && fRes[i][1]; i++) {
    fRes[i][2] = new TF1(Form("%sOver1",fRes[i][1]->GetName()), func[i+2], fRes[i][1]->GetXmin(), fRes[i][1]->GetXmax(), nPar[i+2]);
    Double_t *p = new Double_t[nPar[i+2]];
    fRes[i][1]->GetParameters(p);
    fRes[i][0]->GetParameters(&(p[nPar[i]]));
    fRes[i][2]->SetParameters(p);
    delete[] p;
  }
  
  // print results
  TString var[2] = {"pT", "y"};
  for (Int_t i = 0; i < 2; i++) {
    for (Int_t j = 0; j < 2 && fRes[i][j]; j++) {
      Double_t *param = fRes[i][j]->GetParameters();
      printf("\n%s parameters for single muon generator in file %d:\n", var[i].Data(), j+1);
      printf("Double_t p[%d] = {", nPar[i]);
      for (Int_t k = 0; k < nPar[i]-1; k++) printf("%g, ", param[k]);
      printf("%g};\n", param[nPar[i]-1]);
    }
  }
  
  // display results at the generation level
  TCanvas *cGen = new TCanvas("cGen", "cGen", 1200, 800);
  cGen->Divide(3,2);
  for (Int_t i = 0; i < 3; i++) {
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
    cGen->cd(i+4);
    if (hRes[i][2]) {
      hRes[i][2]->Draw();
      hRes[i][2]->SetMinimum(0.8);
      hRes[i][2]->SetMaximum(1.2);
    }
    if (i < 2 && fRes[i][2]) fRes[i][2]->Draw("sames");
  }
  
  // display results at the reconstruction level
  TCanvas *cRec = new TCanvas("cRec", "cRec", 1200, 800);
  cRec->Divide(3,2);
  cRec->cd(1);
  for (Int_t i = 3; i < 6; i++) {
    cRec->cd(i-2);
    if (i == 3) gPad->SetLogy();
    for (Int_t j = 0; j < 2 && hRes[i][j]; j++) {
      hRes[i][j]->SetLineColor(2*j+2);
      hRes[i][j]->Draw((j == 0) ? "" : "sames");
    }
    cRec->cd(i+1);
    if (hRes[i][2]) {
      hRes[i][2]->Draw();
      hRes[i][2]->SetMinimum(0.8);
      hRes[i][2]->SetMaximum(1.2);
    }
  }
  
}

