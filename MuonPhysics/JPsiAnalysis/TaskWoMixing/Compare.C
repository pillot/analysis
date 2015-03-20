//
//  Compare.C
//  aliroot_dev
//
//  Created by philippe pillot on 07/07/2014.
//  Copyright (c) 2014 Philippe Pillot. All rights reserved.
//

Int_t rebin = 2;
Int_t size0, size1, size2;

Double_t massRange[2] = {2.81, 3.29};
//Double_t massRange[2] = {9.01, 9.99};
//Double_t massRange[2] = {0.01, 13.99};
//Double_t massRange[2] = {2.01, 3.99};

Bool_t printRatio = kFALSE;

Bool_t weight = kTRUE;
Double_t weightVsCent[9] = {
  16280.4,
  10391.7,
  6571.66,
  3640.80,
  2103.33,
  1253.36,
  661.988,
  272.671,
  139.588
};

Bool_t drawHist = kTRUE;
Bool_t drawDiff = kTRUE;
Bool_t drawRatio = kTRUE;

Double_t Getn(TH1F *h);
void ComputeDiffErr(Double_t n[2], Double_t differr[2]);
void ComputeDiffErrW(Double_t n[9][2], Double_t differr[9][2], Int_t centBins[10], Double_t differrw[2]);
void PrintStat(Double_t n[2], Double_t differr[2], TString label);
void PrintStat(Double_t n[2], Double_t differr[2], Double_t differrw[2], TString label);
void DrawStat(Double_t differr[2], TH1F *d, Int_t bin, TString label);

void Compare(TString containerName1 = "cOut_MB", TString containerName2 = "cOut_MB_TrgSign3")
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
  
  // centrality bins
  const Int_t nCent=9;
  TString centBinName[nCent] = {"010", "1020", "2030", "3040", "4050", "5060", "6070", "7080", "8090"};
  
/*  // pT bins
  const Int_t nPtBins = 19;
  Float_t dPtLowEdge[nPtBins] = {0., 0., 2., 5., 0., 1., 2., 3., 4., 5., 6., 8., 2., 0., 3., 0.3, 0.3, 0.3, 0.3};
  Float_t dPtUpEdge[nPtBins] = {8., 2., 5., 8., 1., 2., 3., 4., 5., 6., 8., 20., 8., 3., 8., 1., 8., 2., 3.};
  Int_t ipTmin = 4;
  Int_t ipTmax = 10;
*/  const Int_t nPtBins = 21;
  Float_t dPtLowEdge[nPtBins] = {0., 0., 2., 5., 0., 1., 2., 3., 4., 5., 6., 8., 10., 14., 2., 0., 3., 0.3, 0.3, 0.3, 0.3};
  Float_t dPtUpEdge[nPtBins] = {8., 2., 5., 8., 1., 2., 3., 4., 5., 6., 8., 10., 14., 20., 8., 3., 8., 1., 8., 2., 3.};
  Int_t ipTmin = 4;
  Int_t ipTmax = 13;
  
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
  
  // ###### plots versus centrality ######
  Double_t n[9][2];
  Double_t differr[9][2];
  TH1F *hDiffVsCent = new TH1F("hDiffVsCent","diff vs cent",10,0.,10.);
  TH1F *hInt1 = new TH1F("hDimuPM_cent0-90_1","hDimuPM_cent0-90",560,0.,14.);
  TH1F *hInt2 = new TH1F("hDimuPM_cent0-90_2","hDimuPM_cent0-90",560,0.,14.);
  TCanvas *cHist = 0x0, *cDiff = 0x0, *cRatio = 0x0;
  if (drawHist) {
    cHist = new TCanvas("cCentHist", "histo vs cent", 900, 900);
    cHist->Divide(3,3);
  }
  if (drawDiff) {
    cDiff = new TCanvas("cCentDiff", "diff vs cent", 900, 900);
    cDiff->Divide(3,3);
  }
  if (drawRatio) {
    cRatio = new TCanvas("cCentRatio", "ratio vs cent", 900, 900);
    cRatio->Divide(3,3);
  }
  TString label0 = "centrality";
  size0 = TMath::Max(9,label0.Length());
  printf("\nstat versus centrality\n");
  if (weight) printf(Form("%%%ds %%%ds %%%ds        diff            w-diff\n",size0,size1,size2),label0.Data(),label1.Data(),label2.Data());
  else printf(Form("%%%ds %%%ds %%%ds        diff\n",size0,size1,size2),label0.Data(),label1.Data(),label2.Data());
  for (Int_t i=0; i<nCent; i++) {
    TH1F* h1 = static_cast<TH1F*>(container1->FindObject(Form("hDimuPM_y_0_%s", centBinName[i].Data())));
    if (!h1) continue;
    TH1F* h2 = static_cast<TH1F*>(container2->FindObject(Form("hDimuPM_y_0_%s", centBinName[i].Data())));
    if (!h2) continue;
    n[i][0] = Getn(h1);
    n[i][1] = Getn(h2);
    ComputeDiffErr(n[i], differr[i]);
    hInt1->Add(h1);
    hInt2->Add(h2);
    h1->Rebin(rebin);
    h2->Rebin(rebin);
    if (drawHist) {
      cHist->cd(i+1);
      gPad->SetLogy();
      h1->Draw();
      h2->Draw("sames");
      h2->SetLineColor(2);
      if (i == 0) {
        lHist->AddEntry(h1,label1.Data(),"l");
        lHist->AddEntry(h2,label2.Data(),"l");
        lHist->Draw("same");
      }
    }
    if (drawDiff) {
      cDiff->cd(i+1);
      gPad->SetLogy();
      TH1F* h_diff = static_cast<TH1F*>(h1->Clone());
      h_diff->Add(h2, -1.);
      h_diff->Draw();
      h_diff->SetLineColor(4);
      if (i == 0) {
        lDiff->AddEntry(h_diff,Form("%s - %s",label1.Data(),label2.Data()),"");
        lDiff->Draw("same");
      }
    }
    if (drawRatio) {
      cRatio->cd(i+1);
      TH1F* h_ratio = static_cast<TH1F*>(h2->Clone());
      h_ratio->Divide(h1);
      h_ratio->Draw();
      h_ratio->SetLineColor(4);
      if (i == 0) {
        lRatio->AddEntry(h_ratio,Form("%s / %s",label2.Data(),label1.Data()),"");
        lRatio->Draw("same");
      }
    }
    PrintStat(n[i], differr[i], centBinName[i].Data());
    DrawStat(differr[i], hDiffVsCent, i+1, centBinName[i].Data());
  }
  Double_t ntot[2], differrtot[2], differrtotw[2];
  Int_t centBinstot[10] = {9,  0,  1,  2,  3,  4,  5,  6,  7,  8};
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
    lHist->DrawClone("same");
  }
  if (1) {
    cDiff = new TCanvas("cIntDiff", "diff", 400, 400);
    gPad->SetLogy();
    TH1F* h_diff = static_cast<TH1F*>(hInt1->Clone());
    h_diff->Add(hInt2, -1.);
    h_diff->Draw();
    h_diff->SetLineColor(4);
    lDiff->DrawClone("same");
  }
  if (1) {
    cRatio = new TCanvas("cIntRatio", "ratio", 400, 400);
    TH1F* h_ratio = static_cast<TH1F*>(hInt2->Clone());
    h_ratio->Divide(hInt1);
    h_ratio->Draw();
    h_ratio->SetLineColor(4);
    lRatio->DrawClone("same");
  }
  if (weight) {
    ComputeDiffErrW(n, differr, centBinstot, differrtotw);
    PrintStat(ntot, differrtot, differrtotw, "090");
    DrawStat(differrtotw, hDiffVsCent, 10, "090");
  } else {
    PrintStat(ntot, differrtot, "090");
    DrawStat(differrtot, hDiffVsCent, 10, "090");
  }
  TCanvas* cDiffVsCent = new TCanvas("cDiffVsCent", "diff vs cent", 800, 400);
  hDiffVsCent->Draw();
  outfile->cd();
  hDiffVsCent->Write();
  
  // ###### plots versus centrality in pT bins ######
  for (Int_t ipT = 1; ipT < 4; ++ipT) {
    TH1F *hDiffVsCentpT = new TH1F(Form("hDiffVsCentpT%d",ipT),Form("diff vs cent in %3.1f < pT < %3.1f",dPtLowEdge[ipT],dPtUpEdge[ipT]),7,0.,7.);
    if (drawHist) {
      cHist = new TCanvas(Form("cCentpT%dHist",ipT), Form("histo vs cent in %3.1f < pT < %3.1f",dPtLowEdge[ipT],dPtUpEdge[ipT]), 1200, 600);
      cHist->Divide(4,2);
    }
    if (drawDiff) {
      cDiff = new TCanvas(Form("cCentpT%dDiff",ipT), Form("diff vs cent in %3.1f < pT < %3.1f",dPtLowEdge[ipT],dPtUpEdge[ipT]), 1200, 600);
      cDiff->Divide(4,2);
    }
    if (drawRatio) {
      cRatio = new TCanvas(Form("cCentpT%dRatio",ipT), Form("ratio vs cent in %3.1f < pT < %3.1f",dPtLowEdge[ipT],dPtUpEdge[ipT]), 1200, 600);
      cRatio->Divide(4,2);
    }
    label0 = "centrality";
    size0 = TMath::Max(9,label0.Length());
    printf("\nstat versus centrality in %3.1f < pT < %3.1f\n",dPtLowEdge[ipT],dPtUpEdge[ipT]);
    if (weight) printf(Form("%%%ds %%%ds %%%ds        diff            w-diff\n",size0,size1,size2),label0.Data(),label1.Data(),label2.Data());
    else printf(Form("%%%ds %%%ds %%%ds        diff\n",size0,size1,size2),label0.Data(),label1.Data(),label2.Data());
    for (Int_t i=0; i<7; i++) {
      TH1F* h1 = static_cast<TH1F*>(container1->FindObject(Form("hDimuPM_pt_%d_%s", ipT, centBinName[i].Data())));
      if (!h1) continue;
      TH1F* h2 = static_cast<TH1F*>(container2->FindObject(Form("hDimuPM_pt_%d_%s", ipT, centBinName[i].Data())));
      if (!h2) continue;
      n[i][0] = Getn(h1);
      n[i][1] = Getn(h2);
      ComputeDiffErr(n[i], differr[i]);
      if (i==6) {
        h1 = static_cast<TH1F*>(h1->Clone());
        h1->SetNameTitle(Form("hDimuPM_pt_%d_6090",ipT), Form("hDimuPM_pt_%d_6090",ipT));
        h2 = static_cast<TH1F*>(h2->Clone());
        h2->SetNameTitle(Form("hDimuPM_pt_%d_6090",ipT), Form("hDimuPM_pt_%d_6090",ipT));
        for (Int_t j=7; j<9; j++) {
          TH1F* h1tmp = static_cast<TH1F*>(container1->FindObject(Form("hDimuPM_pt_%d_%s", ipT, centBinName[j].Data())));
          if (!h1tmp) continue;
          TH1F* h2tmp = static_cast<TH1F*>(container2->FindObject(Form("hDimuPM_pt_%d_%s", ipT, centBinName[j].Data())));
          if (!h2tmp) continue;
          n[j][0] = Getn(h1tmp);
          n[j][1] = Getn(h2tmp);
          ComputeDiffErr(n[j], differr[j]);
          h1->Add(h1tmp);
          h2->Add(h2tmp);
        }
        ntot[0] = Getn(h1);
        ntot[1] = Getn(h2);
        ComputeDiffErr(ntot, differrtot);
      }
      h1->Rebin(rebin);
      h2->Rebin(rebin);
      if (drawHist) {
        gROOT->SetSelectedPad(cHist->cd(i+1));
        gPad->SetLogy();
        h1->Draw();
        h2->Draw("sames");
        h2->SetLineColor(2);
        if (i == 0) lHist->DrawClone("same");
      }
      if (drawDiff) {
        gROOT->SetSelectedPad(cDiff->cd(i+1));
        gPad->SetLogy();
        TH1F* h_diff = static_cast<TH1F*>(h1->Clone());
        h_diff->Add(h2, -1.);
        h_diff->Draw();
        h_diff->SetLineColor(4);
        if (i == 0) lDiff->DrawClone("same");
      }
      if (drawRatio) {
        gROOT->SetSelectedPad(cRatio->cd(i+1));
        TH1F* h_ratio = static_cast<TH1F*>(h2->Clone());
        h_ratio->Divide(h1);
        h_ratio->Draw();
        h_ratio->SetLineColor(4);
        if (i == 0) lRatio->DrawClone("same");
      }
      if (i<6) {
        PrintStat(n[i], differr[i], centBinName[i].Data());
        DrawStat(differr[i], hDiffVsCentpT, i+1, centBinName[i].Data());
      } else {
        if (weight) {
          centBinstot[0] = 3; centBinstot[1] = 6; centBinstot[2] = 7; centBinstot[3] = 8;
          ComputeDiffErrW(n, differr, centBinstot, differrtotw);
          PrintStat(ntot, differrtot, differrtotw, "6090");
          DrawStat(differrtotw, hDiffVsCentpT, i+1, "6090");
        } else {
          PrintStat(ntot, differrtot, "6090");
          DrawStat(differrtot, hDiffVsCentpT, i+1, "6090");
        }
      }
    }
    TCanvas* cDiffVsCentpT = new TCanvas(Form("cDiffVsCentpT%d",ipT), Form("diff vs cent in %3.1f < pT < %3.1f",dPtLowEdge[ipT],dPtUpEdge[ipT]), 800, 400);
    hDiffVsCentpT->Draw();
    outfile->cd();
    hDiffVsCentpT->Write();
  }
  
  // ###### plots versus pT in centrality bins ######
  Int_t centBins[5][10] = {
    {2,  0,  1, -1, -1, -1, -1, -1, -1, -1},
    {2,  2,  3, -1, -1, -1, -1, -1, -1, -1},
    {5,  4,  5,  6,  7,  8, -1, -1, -1, -1},
    {4,  0,  1,  2,  3, -1, -1, -1, -1, -1},
    {9,  0,  1,  2,  3,  4,  5,  6,  7,  8}};
  TString centBinLabel[5] = {"020", "2040", "4090", "040", "090"};
  for (Int_t iBin = 0; iBin < 5; ++iBin) {
    TH1F *hDiffVspTCent = new TH1F(Form("hDiffVspTCent%d",iBin),Form("diff vs pT in cent %s%%",centBinLabel[iBin].Data()),ipTmax-ipTmin+1,0.,ipTmax-ipTmin+1);
    if (drawHist) {
      cHist = new TCanvas(Form("cpTCent%dHist",iBin), Form("histo vs pT in cent %s%%",centBinLabel[iBin].Data()), 1200, 900);
      cHist->Divide(4,3);
    }
    if (drawDiff) {
      cDiff = new TCanvas(Form("cpTCent%dDiff",iBin), Form("diff vs pT in cent %s%%",centBinLabel[iBin].Data()), 1200, 900);
      cDiff->Divide(4,3);
    }
    if (drawRatio) {
      cRatio = new TCanvas(Form("cpTCent%dRatio",iBin), Form("ratio vs pT in cent %s%%",centBinLabel[iBin].Data()), 1200, 900);
      cRatio->Divide(4,3);
    }
    label0 = "pT";
    size0 = TMath::Max(9,label0.Length());
    printf("\nstat versus pT in cent %s%%\n",centBinLabel[iBin].Data());
    if (weight) printf(Form("%%%ds %%%ds %%%ds        diff            w-diff\n",size0,size1,size2),label0.Data(),label1.Data(),label2.Data());
    else printf(Form("%%%ds %%%ds %%%ds        diff\n",size0,size1,size2),label0.Data(),label1.Data(),label2.Data());
    for (Int_t ipT = ipTmin; ipT <= ipTmax; ++ipT) {
      TH1F* h1 = static_cast<TH1F*>(container1->FindObject(Form("hDimuPM_pt_%d_%s", ipT, centBinName[centBins[iBin][1]].Data())));
      if (!h1) continue;
      h1 = static_cast<TH1F*>(h1->Clone());
      h1->SetNameTitle(Form("hDimuPM_pt_%d_%s",ipT,centBinLabel[iBin].Data()), Form("hDimuPM_pt_%d_%s",ipT,centBinLabel[iBin].Data()));
      TH1F* h2 = static_cast<TH1F*>(container2->FindObject(Form("hDimuPM_pt_%d_%s", ipT, centBinName[centBins[iBin][1]].Data())));
      if (!h2) continue;
      h2 = static_cast<TH1F*>(h2->Clone());
      h2->SetNameTitle(Form("hDimuPM_pt_%d_%s",ipT,centBinLabel[iBin].Data()), Form("hDimuPM_pt_%d_%s",ipT,centBinLabel[iBin].Data()));
      n[centBins[iBin][1]][0] = Getn(h1);
      n[centBins[iBin][1]][1] = Getn(h2);
      ComputeDiffErr(n[centBins[iBin][1]], differr[centBins[iBin][1]]);
      for (Int_t iCent = 2; iCent <= centBins[iBin][0]; ++iCent) {
        TH1F* h1tmp = static_cast<TH1F*>(container1->FindObject(Form("hDimuPM_pt_%d_%s", ipT, centBinName[centBins[iBin][iCent]].Data()))));
        if (!h1tmp) continue;
        TH1F* h2tmp = static_cast<TH1F*>(container2->FindObject(Form("hDimuPM_pt_%d_%s", ipT, centBinName[centBins[iBin][iCent]].Data()))));
        if (!h2tmp) continue;
        n[centBins[iBin][iCent]][0] = Getn(h1tmp);
        n[centBins[iBin][iCent]][1] = Getn(h2tmp);
        ComputeDiffErr(n[centBins[iBin][iCent]], differr[centBins[iBin][iCent]]);
        h1->Add(h1tmp);
        h2->Add(h2tmp);
      }
      ntot[0] = Getn(h1);
      ntot[1] = Getn(h2);
      ComputeDiffErr(ntot, differrtot);
      h1->Rebin(rebin);
      h2->Rebin(rebin);
      if (drawHist) {
        gROOT->SetSelectedPad(cHist->cd(ipT-ipTmin+1));
        gPad->SetLogy();
        h1->Draw();
        h2->Draw("sames");
        h2->SetLineColor(2);
        if (ipT == ipTmin) lHist->DrawClone("same");
      }
      if (drawDiff) {
        gROOT->SetSelectedPad(cDiff->cd(ipT-ipTmin+1));
        gPad->SetLogy();
        TH1F* h_diff = static_cast<TH1F*>(h1->Clone());
        h_diff->Add(h2, -1.);
        h_diff->Draw();
        h_diff->SetLineColor(4);
        if (ipT == ipTmin) lDiff->DrawClone("same");
      }
      if (drawRatio) {
        gROOT->SetSelectedPad(cRatio->cd(ipT-ipTmin+1));
        TH1F* h_ratio = static_cast<TH1F*>(h2->Clone());
        h_ratio->Divide(h1);
        h_ratio->Draw();
        h_ratio->SetLineColor(4);
        if (ipT == ipTmin) lRatio->DrawClone("same");
      }
      if (weight) {
        ComputeDiffErrW(n, differr, centBins[iBin], differrtotw);
        PrintStat(ntot, differrtot, differrtotw, Form("%3.1f-%3.1f",dPtLowEdge[ipT],dPtUpEdge[ipT]));
        DrawStat(differrtotw, hDiffVspTCent, ipT-3, Form("%3.1f-%3.1f",dPtLowEdge[ipT],dPtUpEdge[ipT]));
      } else {
        PrintStat(ntot, differrtot, Form("%3.1f-%3.1f",dPtLowEdge[ipT],dPtUpEdge[ipT]));
        DrawStat(differrtot, hDiffVspTCent, ipT-3, Form("%3.1f-%3.1f",dPtLowEdge[ipT],dPtUpEdge[ipT]));
      }
    }
    TCanvas* cDiffVspTCent = new TCanvas(Form("cDiffVspTCent%d",iBin), Form("diff vs pT in cent %s%%",centBinLabel[iBin].Data()), 800, 400);
    hDiffVspTCent->Draw();
    outfile->cd();
    hDiffVspTCent->Write();
  }
  
  // ###### plots versus centrality in y bins ######
  for (Int_t iy = 1; iy < 4; ++iy) {
    TH1F *hDiffVsCenty = new TH1F(Form("hDiffVsCenty%d",iy),Form("diff vs cent in %4.2f < y < %4.2f",dYLowEdge[iy],dYUpEdge[iy]),7,0.,7.);
    if (drawHist) {
      cHist = new TCanvas(Form("cCenty%dHist",iy), Form("histo vs cent in %4.2f < y < %4.2f",dYLowEdge[iy],dYUpEdge[iy]), 1200, 600);
      cHist->Divide(4,2);
    }
    if (drawDiff) {
      cDiff = new TCanvas(Form("cCenty%dDiff",iy), Form("diff vs cent in %4.2f < y < %4.2f",dYLowEdge[iy],dYUpEdge[iy]), 1200, 600);
      cDiff->Divide(4,2);
    }
    if (drawRatio) {
      cRatio = new TCanvas(Form("cCenty%dRatio",iy), Form("ratio vs cent in %4.2f < y < %4.2f",dYLowEdge[iy],dYUpEdge[iy]), 1200, 600);
      cRatio->Divide(4,2);
    }
    label0 = "centrality";
    size0 = TMath::Max(9,label0.Length());
    printf("\nstat versus centrality in %4.2f < y < %4.2f\n",dYLowEdge[iy],dYUpEdge[iy]);
    if (weight) printf(Form("%%%ds %%%ds %%%ds        diff            w-diff\n",size0,size1,size2),label0.Data(),label1.Data(),label2.Data());
    else printf(Form("%%%ds %%%ds %%%ds        diff\n",size0,size1,size2),label0.Data(),label1.Data(),label2.Data());
    for (Int_t i=0; i<7; i++) {
      TH1F* h1 = static_cast<TH1F*>(container1->FindObject(Form("hDimuPM_y_%d_%s", iy, centBinName[i].Data())));
      if (!h1) continue;
      TH1F* h2 = static_cast<TH1F*>(container2->FindObject(Form("hDimuPM_y_%d_%s", iy, centBinName[i].Data())));
      if (!h2) continue;
      n[i][0] = Getn(h1);
      n[i][1] = Getn(h2);
      ComputeDiffErr(n[i], differr[i]);
      if (i==6) {
        h1 = static_cast<TH1F*>(h1->Clone());
        h1->SetNameTitle(Form("hDimuPM_y_%d_6090",iy), Form("hDimuPM_y_%d_6090",iy));
        h2 = static_cast<TH1F*>(h2->Clone());
        h2->SetNameTitle(Form("hDimuPM_y_%d_6090",iy), Form("hDimuPM_y_%d_6090",iy));
        for (Int_t j=7; j<9; j++) {
          TH1F* h1tmp = static_cast<TH1F*>(container1->FindObject(Form("hDimuPM_y_%d_%s", iy, centBinName[j].Data())));
          if (!h1tmp) continue;
          TH1F* h2tmp = static_cast<TH1F*>(container2->FindObject(Form("hDimuPM_y_%d_%s", iy, centBinName[j].Data())));
          if (!h2tmp) continue;
          n[j][0] = Getn(h1tmp);
          n[j][1] = Getn(h2tmp);
          ComputeDiffErr(n[j], differr[j]);
          h1->Add(h1tmp);
          h2->Add(h2tmp);
        }
        ntot[0] = Getn(h1);
        ntot[1] = Getn(h2);
        ComputeDiffErr(ntot, differrtot);
      }
      h1->Rebin(rebin);
      h2->Rebin(rebin);
      if (drawHist) {
        gROOT->SetSelectedPad(cHist->cd(i+1));
        gPad->SetLogy();
        h1->Draw();
        h2->Draw("sames");
        h2->SetLineColor(2);
        if (i == 0) lHist->DrawClone("same");
      }
      if (drawDiff) {
        gROOT->SetSelectedPad(cDiff->cd(i+1));
        gPad->SetLogy();
        TH1F* h_diff = static_cast<TH1F*>(h1->Clone());
        h_diff->Add(h2, -1.);
        h_diff->Draw();
        h_diff->SetLineColor(4);
        if (i == 0) lDiff->DrawClone("same");
      }
      if (drawRatio) {
        gROOT->SetSelectedPad(cRatio->cd(i+1));
        TH1F* h_ratio = static_cast<TH1F*>(h2->Clone());
        h_ratio->Divide(h1);
        h_ratio->Draw();
        h_ratio->SetLineColor(4);
        if (i == 0) lRatio->DrawClone("same");
      }
      if (i<6) {
        PrintStat(n[i], differr[i], centBinName[i].Data());
        DrawStat(differr[i], hDiffVsCenty, i+1, centBinName[i].Data());
      } else {
        if (weight) {
          centBinstot[0] = 3; centBinstot[1] = 6; centBinstot[2] = 7; centBinstot[3] = 8;
          ComputeDiffErrW(n, differr, centBinstot, differrtotw);
          PrintStat(ntot, differrtot, differrtotw, "6090");
          DrawStat(differrtotw, hDiffVsCenty, i+1, "6090");
        } else {
          PrintStat(ntot, differrtot, "6090");
          DrawStat(differrtot, hDiffVsCenty, i+1, "6090");
        }
      }
    }
    TCanvas* cDiffVsCenty = new TCanvas(Form("cDiffVsCenty%d",iy), Form("diff vs cent in %4.2f < y < %4.2f",dYLowEdge[iy],dYUpEdge[iy]), 800, 400);
    hDiffVsCenty->Draw();
    outfile->cd();
    hDiffVsCenty->Write();
  }
  
  // ###### plots versus y in centrality bins ######
  for (Int_t iBin = 0; iBin < 5; ++iBin) {
    TH1F *hDiffVsyCent = new TH1F(Form("hDiffVsyCent%d",iBin),Form("diff vs y in cent %s%%",centBinLabel[iBin].Data()),6,0.,6.);
    if (drawHist) {
      cHist = new TCanvas(Form("cyCent%dHist",iBin), Form("histo vs y in cent %s%%",centBinLabel[iBin].Data()), 900, 600);
      cHist->Divide(3,2);
    }
    if (drawDiff) {
      cDiff = new TCanvas(Form("cyCent%dDiff",iBin), Form("diff vs y in cent %s%%",centBinLabel[iBin].Data()), 900, 600);
      cDiff->Divide(3,2);
    }
    if (drawRatio) {
      cRatio = new TCanvas(Form("cyCent%dRatio",iBin), Form("ratio vs y in cent %s%%",centBinLabel[iBin].Data()), 900, 600);
      cRatio->Divide(3,2);
    }
    label0 = "y";
    size0 = TMath::Max(9,label0.Length());
    printf("\nstat versus y in cent %s%%\n",centBinLabel[iBin].Data());
    if (weight) printf(Form("%%%ds %%%ds %%%ds        diff            w-diff\n",size0,size1,size2),label0.Data(),label1.Data(),label2.Data());
    else printf(Form("%%%ds %%%ds %%%ds        diff\n",size0,size1,size2),label0.Data(),label1.Data(),label2.Data());
    for (Int_t iy = 4; iy < 10; ++iy) {
      TH1F* h1 = static_cast<TH1F*>(container1->FindObject(Form("hDimuPM_y_%d_%s", iy, centBinName[centBins[iBin][1]].Data())));
      if (!h1) continue;
      h1 = static_cast<TH1F*>(h1->Clone());
      h1->SetNameTitle(Form("hDimuPM_y_%d_%s",iy,centBinLabel[iBin].Data()), Form("hDimuPM_y_%d_%s",iy,centBinLabel[iBin].Data()));
      TH1F* h2 = static_cast<TH1F*>(container2->FindObject(Form("hDimuPM_y_%d_%s", iy, centBinName[centBins[iBin][1]].Data())));
      if (!h2) continue;
      h2 = static_cast<TH1F*>(h2->Clone());
      h2->SetNameTitle(Form("hDimuPM_y_%d_%s",iy,centBinLabel[iBin].Data()), Form("hDimuPM_y_%d_%s",iy,centBinLabel[iBin].Data()));
      n[centBins[iBin][1]][0] = Getn(h1);
      n[centBins[iBin][1]][1] = Getn(h2);
      ComputeDiffErr(n[centBins[iBin][1]], differr[centBins[iBin][1]]);
      for (Int_t iCent = 2; iCent <= centBins[iBin][0]; ++iCent) {
        TH1F* h1tmp = static_cast<TH1F*>(container1->FindObject(Form("hDimuPM_y_%d_%s", iy, centBinName[centBins[iBin][iCent]].Data()))));
        if (!h1tmp) continue;
        TH1F* h2tmp = static_cast<TH1F*>(container2->FindObject(Form("hDimuPM_y_%d_%s", iy, centBinName[centBins[iBin][iCent]].Data()))));
        if (!h2tmp) continue;
        n[centBins[iBin][iCent]][0] = Getn(h1tmp);
        n[centBins[iBin][iCent]][1] = Getn(h2tmp);
        ComputeDiffErr(n[centBins[iBin][iCent]], differr[centBins[iBin][iCent]]);
        h1->Add(h1tmp);
        h2->Add(h2tmp);
      }
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
      if (weight) {
        ComputeDiffErrW(n, differr, centBins[iBin], differrtotw);
        PrintStat(ntot, differrtot, differrtotw, Form("%4.2f-%4.2f",dYLowEdge[iy],dYUpEdge[iy]));
        DrawStat(differrtotw, hDiffVsyCent, iy-3, Form("%4.2f-%4.2f",dYLowEdge[iy],dYUpEdge[iy]));
      } else {
        PrintStat(ntot, differrtot, Form("%4.2f-%4.2f",dYLowEdge[iy],dYUpEdge[iy]));
        DrawStat(differrtot, hDiffVsyCent, iy-3, Form("%4.2f-%4.2f",dYLowEdge[iy],dYUpEdge[iy]));
      }
    }
    TCanvas* cDiffVsyCent = new TCanvas(Form("cDiffVsyCent%d",iBin), Form("diff vs y in cent %s%%",centBinLabel[iBin].Data()), 800, 400);
    hDiffVsyCent->Draw();
    outfile->cd();
    hDiffVsyCent->Write();
  }
  
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

void ComputeDiffErrW(Double_t n[9][2], Double_t differr[9][2], Int_t centBins[10], Double_t differrw[2])
{
  /// compute the weigthed relative difference and associated error between the 2 numbers of dimuons
  
  differrw[0] = 0.;
  differrw[1] = 0.;
  Double_t sumw = 0.;
  
  for (Int_t iCent = 1; iCent <= centBins[0]; ++iCent) {
    sumw += weightVsCent[centBins[iCent]];
    differrw[0] += weightVsCent[centBins[iCent]]*differr[centBins[iCent]][0];
    differrw[1] += weightVsCent[centBins[iCent]]*weightVsCent[centBins[iCent]]*differr[centBins[iCent]][1]*differr[centBins[iCent]][1];
  }
  
  differrw[0] /= sumw;
  differrw[1] = TMath::Sqrt(differrw[1])/sumw;
  
}

void PrintStat(Double_t n[2], Double_t differr[2], TString label)
{
  /// print number of dimuons in 2.8-3.3 GeV/c^2
  
  printf(Form("%%%ds %%%dd %%%dd   %%4.2f ± %%4.2f %%%%\n",
              size0,size1,size2), label.Data(), (Int_t)n[0], (Int_t)n[1], 100.*differr[0], 100.*differr[1]);
  
}

void PrintStat(Double_t n[2], Double_t differr[2], Double_t differrw[2], TString label)
{
  /// print number of dimuons in 2.8-3.3 GeV/c^2
  
  printf(Form("%%%ds %%%dd %%%dd   %%4.2f ± %%4.2f %%%%   %%4.2f ± %%4.2f %%%%\n",
              size0,size1,size2), label.Data(), (Int_t)n[0], (Int_t)n[1], 100.*differr[0], 100.*differr[1], 100.*differrw[0], 100.*differrw[1]);
  
}

void DrawStat(Double_t differr[2], TH1F *d, Int_t bin, TString label)
{
  /// print number of dimuons in 2.8-3.3 GeV/c^2
  
  d->SetBinContent(bin, 100.*differr[0]);
  d->SetBinError(bin, 100.*differr[1]);
  d->GetXaxis()->SetBinLabel(bin, label.Data());
  
}
