//
//  CompareNoCent.C
//  aliroot-dev
//
//  Created by philippe pillot on 13/01/2015.
//  Copyright (c) 2015 Philippe Pillot. All rights reserved.
//

Int_t rebin = 2;
Int_t size0, size1, size2;

Double_t massRange[2] = {0.91, 1.19};
//Double_t massRange[2] = {2.81, 3.29};
//Double_t massRange[2] = {9.01, 9.99};
//Double_t massRange[2] = {0., 13.99};
//Double_t massRange[2] = {2.01, 3.99};

Bool_t printRatio = kFALSE;

Bool_t drawHist = kTRUE;
Bool_t drawDiff = kTRUE;
Bool_t drawRatio = kTRUE;

Double_t Getn(TH1F *h);
void ComputeDiffErr(Double_t n[2], Double_t differr[2]);
void PrintStat(Double_t n[2], Double_t differr[2], TString label);
void DrawStat(Double_t differr[2], TH1F *d, Int_t bin, TString label);

void CompareNoCent(TString containerName1 = "cOut_MULorMLL", TString containerName2 = "cOut_MULorMLL_TrgSign3")
{
  /// compare histos between the 2 containers
  
  gROOT->SetStyle("Plain");
  gROOT->ForceStyle();
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  
  gROOT->LoadMacro("$WORK/Macros/Facilities/runTaskFacilities.C");
  LoadAlirootLocally("CORRFW", "", "AliAnalysisTaskJPsi");
  
  TString label1 = containerName1;
  label1.ReplaceAll("cOut_","");
  size1 = TMath::Max(9,label1.Length());
  TString label2 = containerName2;
  label2.ReplaceAll("cOut_","");
  size2 = TMath::Max(9,label2.Length());
  
  TFile* outfile = TFile::Open(Form("Diff_%s_%s_%4.2f-%4.2f.root",label1.Data(),label2.Data(),massRange[0],massRange[1]),"RECREATE");
  
  TFile* file = TFile::Open("Output.root","READ");
  TList* container1 = static_cast<TList*>(file->FindObjectAny(containerName1.Data()));
  TList* container2 = static_cast<TList*>(file->FindObjectAny(containerName2.Data()));
  if (!container1 || !container2) return;
  
  // pT bins
  const Int_t nPtBins = 15;
  Float_t dPtLowEdge[nPtBins] = {0., 0., 1., 2., 3., 4., 5., 6., 7., 8., 10., 12., 14., 16., 18.};
  Float_t dPtUpEdge[nPtBins] = {20., 1., 2., 3., 4., 5., 6., 7., 8., 10., 12., 14., 16., 18., 20.};
  
  // y bins
  const Int_t nYBins = 10;
  Float_t dYLowEdge[nYBins] = { 2.5, 2.5, 3., 3.5, 2.5, 2.75, 3., 3.25, 3.5, 3.75};
  Float_t dYUpEdge[nYBins] = { 4., 3., 3.5, 4., 2.75, 3., 3.25, 3.5, 3.75, 4.};
  
  // legends
  TLegend *lHist = new TLegend(0.3,0.65,0.9,0.8);
  lHist->SetFillStyle(0);
  lHist->SetBorderSize(0);
  TLegend *lDiff = new TLegend(0.1,0.72,0.9,0.8);
  lDiff->SetFillStyle(0);
  lDiff->SetBorderSize(0);
  TLegend *lRatio = new TLegend(0.0,0.22,0.8,0.3);
  lRatio->SetFillStyle(0);
  lRatio->SetBorderSize(0);
  
  // ###### integrated over pT/y ######
  TH1F *hInt1 = new TH1F("hDimuPM_1","hDimuPM",560,0.,14.);
  TH1F *hInt2 = new TH1F("hDimuPM_2","hDimuPM",560,0.,14.);
  TCanvas *cHist = 0x0, *cDiff = 0x0, *cRatio = 0x0;
  TString label0 = "centrality";
  size0 = TMath::Max(9,label0.Length());
  printf("\nstat all integrated\n");
  printf(Form("%%%ds %%%ds %%%ds        diff\n",size0,size1,size2),label0.Data(),label1.Data(),label2.Data());
  TH1F* h1 = static_cast<TH1F*>(container1->FindObject(Form("hDimuPM_y_0_any")));
  if (!h1) continue;
  TH1F* h2 = static_cast<TH1F*>(container2->FindObject(Form("hDimuPM_y_0_any")));
  if (!h2) continue;
  hInt1->Add(h1);
  hInt2->Add(h2);
  Double_t ntot[2], differrtot[2];
  ntot[0] = Getn(hInt1);
  ntot[1] = Getn(hInt2);
  ComputeDiffErr(ntot, differrtot);
  hInt1->Rebin(rebin);
  hInt2->Rebin(rebin);
  if (1) {
    cHist = new TCanvas("cIntHist", "histo", 400, 400);
    gPad->SetLogy();
    hInt1->Draw();
    hInt1->SetLineColor(4);
    hInt2->Draw("sames");
    hInt2->SetLineColor(2);
    lHist->AddEntry(hInt1,label1.Data(),"l");
    lHist->AddEntry(hInt2,label2.Data(),"l");
    lHist->Draw("same");
  }
  if (1) {
    cDiff = new TCanvas("cIntDiff", "diff", 400, 400);
    gPad->SetLogy();
    TH1F* h_diff = static_cast<TH1F*>(hInt1->Clone());
    h_diff->Add(hInt2, -1.);
    h_diff->Draw();
    h_diff->SetLineColor(4);
    lDiff->AddEntry(h_diff,Form("%s - %s",label1.Data(),label2.Data()),"");
    lDiff->Draw("same");
  }
  if (1) {
    cRatio = new TCanvas("cIntRatio", "ratio", 400, 400);
    TH1F* h_ratio = static_cast<TH1F*>(hInt2->Clone());
    h_ratio->Divide(hInt1);
    h_ratio->Draw();
    h_ratio->SetLineColor(4);
    lRatio->AddEntry(h_ratio,Form("%s / %s",label2.Data(),label1.Data()),"");
    lRatio->Draw("same");
  }
  PrintStat(ntot, differrtot, "any");
  
  // ###### plots versus pT ######
  TH1F *hDiffVspT = new TH1F("hDiffVspT","diff vs pT",nPtBins-1,0.,nPtBins-1.);
  if (drawHist) {
    cHist = new TCanvas("cpTHist", "histo vs pT", 1500, 600);
    cHist->Divide(5,3);
  }
  if (drawDiff) {
    cDiff = new TCanvas("cpTDiff", "diff vs pT", 1500, 600);
    cDiff->Divide(5,3);
  }
  if (drawRatio) {
    cRatio = new TCanvas("cpTRatio", "ratio vs pT", 1500, 600);
    cRatio->Divide(5,3);
  }
  label0 = "pT";
  size0 = TMath::Max(9,label0.Length());
  printf("\nstat versus pT\n");
  printf(Form("%%%ds %%%ds %%%ds        diff\n",size0,size1,size2),label0.Data(),label1.Data(),label2.Data());
  for (Int_t ipT = 1; ipT < nPtBins; ++ipT) {
    TH1F* h1 = static_cast<TH1F*>(container1->FindObject(Form("hDimuPM_pt_%d_any", ipT)));
    if (!h1) continue;
    h1 = static_cast<TH1F*>(h1->Clone());
    TH1F* h2 = static_cast<TH1F*>(container2->FindObject(Form("hDimuPM_pt_%d_any", ipT)));
    if (!h2) continue;
    h2 = static_cast<TH1F*>(h2->Clone());
    ntot[0] = Getn(h1);
    ntot[1] = Getn(h2);
    ComputeDiffErr(ntot, differrtot);
    h1->Rebin(rebin);
    h2->Rebin(rebin);
    if (drawHist) {
      gROOT->SetSelectedPad(cHist->cd(ipT));
      gPad->SetLogy();
      h1->Draw();
      h2->Draw("sames");
      h2->SetLineColor(2);
      if (ipT == 1) lHist->DrawClone("same");
    }
    if (drawDiff) {
      gROOT->SetSelectedPad(cDiff->cd(ipT));
      gPad->SetLogy();
      TH1F* h_diff = static_cast<TH1F*>(h1->Clone());
      h_diff->Add(h2, -1.);
      h_diff->Draw();
      h_diff->SetLineColor(4);
      if (ipT == 1) lDiff->DrawClone("same");
    }
    if (drawRatio) {
      gROOT->SetSelectedPad(cRatio->cd(ipT));
      TH1F* h_ratio = static_cast<TH1F*>(h2->Clone());
      h_ratio->Divide(h1);
      h_ratio->Draw();
      h_ratio->SetLineColor(4);
      if (ipT == 1) lRatio->DrawClone("same");
    }
    PrintStat(ntot, differrtot, Form("%3.1f-%3.1f",dPtLowEdge[ipT],dPtUpEdge[ipT]));
    DrawStat(differrtot, hDiffVspT, ipT, Form("%3.1f-%3.1f",dPtLowEdge[ipT],dPtUpEdge[ipT]));
  }
  TCanvas* cDiffVspT = new TCanvas("cDiffVspT", "diff vs pT", 800, 400);
  hDiffVspT->Draw();
  outfile->cd();
  hDiffVspT->Write();
  
  // ###### plots versus y ######
  TH1F *hDiffVsy = new TH1F("hDiffVsy","diff vs y",6,0.,6.);
  if (drawHist) {
    cHist = new TCanvas("cyHist", "histo vs y", 900, 600);
    cHist->Divide(3,2);
  }
  if (drawDiff) {
    cDiff = new TCanvas("cyDiff", "diff vs y", 900, 600);
    cDiff->Divide(3,2);
  }
  if (drawRatio) {
    cRatio = new TCanvas("cyRatio", "ratio vs y", 900, 600);
    cRatio->Divide(3,2);
  }
  label0 = "y";
  size0 = TMath::Max(9,label0.Length());
  printf("\nstat versus y\n");
  printf(Form("%%%ds %%%ds %%%ds        diff\n",size0,size1,size2),label0.Data(),label1.Data(),label2.Data());
  for (Int_t iy = 4; iy < 10; ++iy) {
    TH1F* h1 = static_cast<TH1F*>(container1->FindObject(Form("hDimuPM_y_%d_any", iy)));
    if (!h1) continue;
    h1 = static_cast<TH1F*>(h1->Clone());
    TH1F* h2 = static_cast<TH1F*>(container2->FindObject(Form("hDimuPM_y_%d_any", iy)));
    if (!h2) continue;
    h2 = static_cast<TH1F*>(h2->Clone());
    ntot[0] = Getn(h1);
    ntot[1] = Getn(h2);
    ComputeDiffErr(ntot, differrtot);
    h1->Rebin(rebin);
    h2->Rebin(rebin);
    if (drawHist) {
      gROOT->SetSelectedPad(cHist->cd(iy-3));
      gPad->SetLogy();
      h1->Draw();
      h2->Draw("sames");
      h2->SetLineColor(2);
      if (iy == 4) lHist->DrawClone("same");
    }
    if (drawDiff) {
      gROOT->SetSelectedPad(cDiff->cd(iy-3));
      gPad->SetLogy();
      TH1F* h_diff = static_cast<TH1F*>(h1->Clone());
      h_diff->Add(h2, -1.);
      h_diff->Draw();
      h_diff->SetLineColor(4);
      if (iy == 4) lDiff->DrawClone("same");
    }
    if (drawRatio) {
      gROOT->SetSelectedPad(cRatio->cd(iy-3));
      TH1F* h_ratio = static_cast<TH1F*>(h2->Clone());
      h_ratio->Divide(h1);
      h_ratio->Draw();
      h_ratio->SetLineColor(4);
      if (iy == 4) lRatio->DrawClone("same");
    }
    PrintStat(ntot, differrtot, Form("%3.1f-%3.1f",dYLowEdge[iy],dYUpEdge[iy]));
    DrawStat(differrtot, hDiffVsy, iy-3, Form("%3.1f-%3.1f",dYLowEdge[iy],dYUpEdge[iy]));
  }
  TCanvas* cDiffVsy = new TCanvas("cDiffVsy", "diff vs y", 800, 400);
  hDiffVsy->Draw();
  outfile->cd();
  hDiffVsy->Write();
  
}

Double_t Getn(TH1F *h)
{
  /// get the number of dimuons in 2.8-3.3 GeV/c^2
  
  return h->Integral(h->FindBin(massRange[0]), h->FindBin(massRange[1]));
  
}

void ComputeDiffErr(Double_t n[2], Double_t differr[2])
{
  /// compute the relative difference and associated error between the 2 numbers of dimuons
  
  differr[0] = 0.;
  differr[1] = 1.;
  
  if (n[0] > 0.) {
    
    if (printRatio) differr[0] = -n[1]/n[0];
    else differr[0] = (n[1]-n[0])/n[0];
    
    if (n[1] > 0.)
      differr[1] = TMath::Max(1./TMath::Max(n[0],n[1]),n[1]/n[0]*TMath::Sqrt(TMath::Abs(1./n[1]-1./n[0])));
    else
      differr[1] = 1./n[0];
    
  }
  
}

void PrintStat(Double_t n[2], Double_t differr[2], TString label)
{
  /// print number of dimuons in 2.8-3.3 GeV/c^2
  
  printf(Form("%%%ds %%%dd %%%dd   %%4.2f ± %%4.2f %%%%\n",
              size0,size1,size2), label.Data(), (Int_t)n[0], (Int_t)n[1], 100.*differr[0], 100.*differr[1]);
  
}

void DrawStat(Double_t differr[2], TH1F *d, Int_t bin, TString label)
{
  /// print number of dimuons in 2.8-3.3 GeV/c^2
  
  d->SetBinContent(bin, 100.*differr[0]);
  d->SetBinError(bin, 100.*differr[1]);
  d->GetXaxis()->SetBinLabel(bin, label.Data());
  
}
