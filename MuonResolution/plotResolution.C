/*
 *  plotResolution.C
 *  aliroot_dev
 *
 *  Created by philippe pillot on 25/04/11.
 *  Copyright 2011 SUBATECH. All rights reserved.
 *
 */

//-------------------------------------------------------------------------
void plotResolution(const char* fileName)
{
  /// plot resolution versus centrality
  
  // resolution per chamber
  for (Int_t i = 1; i <= 10; i++) plotResolutionPerChamber(fileName, i);
  
  // integrated resolution
  plotResolutionPerChamber(fileName, -1);
  
}


//-------------------------------------------------------------------------
void plotResolutionPerChamber(const char* fileName, Int_t ch)
{
  /// plot resolution versus centrality
  
  gStyle->SetFillColor(0);
  gStyle->SetOptStat(0);
  
  // open first file
  TFile* outFile = TFile::Open(fileName,"READ");
  if (!outFile || !outFile->IsOpen()) return;
  
  // get results
  TObjArray* list = static_cast<TObjArray*>(outFile->FindObjectAny("ResidualsVsCent"));
  if (!list) return;
  TH2* hResidualInChVsCent[4] = {0};
  if (ch >= 1 && ch <= 10) {
    hResidualInChVsCent[0] = static_cast<TH2*>(list->FindObject(Form("hResidualXInCh%dVsCent_ClusterIn",ch)));
    hResidualInChVsCent[1] = static_cast<TH2*>(list->FindObject(Form("hResidualYInCh%dVsCent_ClusterIn",ch)));
    hResidualInChVsCent[2] = static_cast<TH2*>(list->FindObject(Form("hResidualXInCh%dVsCent_ClusterOut",ch)));
    hResidualInChVsCent[3] = static_cast<TH2*>(list->FindObject(Form("hResidualYInCh%dVsCent_ClusterOut",ch)));
  } else {
    hResidualInChVsCent[0] = static_cast<TH2*>(list->FindObject("hResidualXInCh1VsCent_ClusterIn"))->Clone();
    hResidualInChVsCent[1] = static_cast<TH2*>(list->FindObject("hResidualYInCh1VsCent_ClusterIn"))->Clone();
    hResidualInChVsCent[2] = static_cast<TH2*>(list->FindObject("hResidualXInCh1VsCent_ClusterOut"))->Clone();
    hResidualInChVsCent[3] = static_cast<TH2*>(list->FindObject("hResidualYInCh1VsCent_ClusterOut"))->Clone();
    for (Int_t i = 2; i < 10; i++) {
      hResidualInChVsCent[0]->Add(static_cast<TH2*>(list->FindObject(Form("hResidualXInCh%dVsCent_ClusterIn",i))));
      hResidualInChVsCent[1]->Add(static_cast<TH2*>(list->FindObject(Form("hResidualYInCh%dVsCent_ClusterIn",i))));
      hResidualInChVsCent[2]->Add(static_cast<TH2*>(list->FindObject(Form("hResidualXInCh%dVsCent_ClusterOut",i))));
      hResidualInChVsCent[3]->Add(static_cast<TH2*>(list->FindObject(Form("hResidualYInCh%dVsCent_ClusterOut",i))));
    }
  }
  
  // display
  TH1D* p;
  TString pName1[4] = {"x4080in", "y4080in", "x4080out", "y4080out"};
  TString pName2[4] = {"x5in", "y5in", "x5out", "y5out"};
  Double_t norm;
  Int_t rebin = 10;
  Double_t range[4] = {0.5, 0.5, 1., 1.};
  TF1* f = new TF1("f", "gaus", -2., 2.);
  Double_t sigma1, sigma2;
  TCanvas* c1;
  if (ch >= 1 && ch <= 10) {
    c1 = new TCanvas(Form("ch%d",ch),Form("ch%d",ch),1200,900);
    printf("\n--- chamber %d:\n", ch);
  }
  else {
    c1 = new TCanvas("integrated","integrated",1200,900);
    printf("\n--- integrated:\n");
  }
  c1->Divide(2,2);
  printf("sigma in 40-80%%\tsigma 0-5%%\tdifference (%%)\n");
  for (Int_t i = 0; i < 4; i++) {
    c1->cd(i+1);
    gPad->SetLogy();
    p = hResidualInChVsCent[i]->ProjectionY(pName1[i].Data(),10,17,"e"); // 40-80%
    p->SetDirectory(0);
    p->Rebin(rebin);
    norm = p->GetEntries();
    p->SetLineColor(4);
    p->GetXaxis()->SetRangeUser(-range[i], range[i]);
    p->Draw();
    p->Fit("f", "NQ");
    sigma1 = f->GetParameter(2);
    f->SetLineColor(4);
    f->SetLineWidth(1);
    f->DrawCopy("same");
    p = hResidualInChVsCent[i]->ProjectionY(pName2[i].Data(),2,2,"e"); // 5%
    p->SetDirectory(0);
    p->Rebin(rebin);
    norm /= p->GetEntries();
    p->SetLineColor(2);
    p->Scale(norm);
    p->Draw("sames");
    p->Fit("f", "NQ");
    sigma2 = f->GetParameter(2);
    f->SetLineColor(2);
    f->SetLineWidth(1);
    f->DrawCopy("same");
    printf("%f\t%f\t%f\n", sigma1, sigma2, 100.*(sigma2-sigma1)/sigma1);
  }
  
  // close file and clean memory
  f->Delete();
  outFile->Close();
  
}

