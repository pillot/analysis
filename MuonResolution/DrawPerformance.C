/*
 *  DrawPerformance.C
 *  aliroot_dev
 *
 *  Created by philippe pillot on 12/05/11.
 *  Copyright 2011 SUBATECH. All rights reserved.
 *
 */

void ALICEseal(TString type, Double_t xPad, Double_t yPad);

//------------------------------------------------------------------------
void DrawPerformance(TString file = "results.root")
{
  gStyle->SetCanvasColor(10);
  gStyle->SetFrameFillColor(10);
  gStyle->SetStatColor(10);
  gStyle->SetFillColor(10);
  gStyle->SetTitleFillColor(10);
  gROOT->SetStyle("Plain");
  gROOT->ForceStyle();
  
  TFile* outFile1 = TFile::Open(file.Data(),"READ");
  if (!outFile1 || !outFile1->IsOpen()) return;
  TObjArray* list1 = static_cast<TObjArray*>(outFile1->FindObjectAny("ChamberRes"));
  if (!list1) return;
  TGraphErrors* gCombinedResidualXSigmaVsCent = static_cast<TGraphErrors*>(list1->FindObject("gCombinedResidualXSigmaVsCent"));
  TGraphErrors* gCombinedResidualYSigmaVsCent = static_cast<TGraphErrors*>(list1->FindObject("gCombinedResidualYSigmaVsCent"));
  
  new TCanvas("resolutionVsCent","resolutionVsCent");
  gPad->SetTopMargin(0.05);
  gPad->SetRightMargin(0.05);
  gCombinedResidualXSigmaVsCent->GetXaxis()->SetTitle("centrality percentile");
  gCombinedResidualXSigmaVsCent->GetXaxis()->SetRangeUser(0.,90.);
  gCombinedResidualXSigmaVsCent->GetXaxis()->SetTitleSize(0.045);
  gCombinedResidualXSigmaVsCent->GetXaxis()->SetTitleOffset(0.99);
  gCombinedResidualXSigmaVsCent->GetXaxis()->SetLabelSize(0.045);
  gCombinedResidualXSigmaVsCent->GetYaxis()->SetTitle("cluster resolution (cm)");
  gCombinedResidualXSigmaVsCent->GetYaxis()->SetTitleSize(0.045);
  gCombinedResidualXSigmaVsCent->GetYaxis()->SetTitleOffset(1.1);
  gCombinedResidualXSigmaVsCent->GetYaxis()->SetLabelSize(0.045);
  gCombinedResidualXSigmaVsCent->GetYaxis()->SetRangeUser(0., 0.12);
  gCombinedResidualXSigmaVsCent->SetMarkerStyle(20);
  gCombinedResidualXSigmaVsCent->SetMarkerColor(2);
  gCombinedResidualXSigmaVsCent->SetLineWidth(2);
  gCombinedResidualXSigmaVsCent->SetLineColor(2);
  gCombinedResidualXSigmaVsCent->Draw("ap");
  gCombinedResidualYSigmaVsCent->SetMarkerStyle(21);
  gCombinedResidualYSigmaVsCent->SetMarkerColor(4);
  gCombinedResidualYSigmaVsCent->SetLineWidth(2);
  gCombinedResidualYSigmaVsCent->SetLineColor(4);
  gCombinedResidualYSigmaVsCent->Draw("p");
  TLegend *lg = new TLegend(0.15, 0.7, 0.8, 0.9,"muon spectrometer resolution:","NDC");
  lg->AddEntry(gCombinedResidualXSigmaVsCent,"non-bending direction","p");
  lg->AddEntry(gCombinedResidualYSigmaVsCent,"bending direction","p");
  lg->SetFillStyle(0);
  lg->SetBorderSize(0);
  lg->SetTextFont(52);
  lg->Draw();
  ALICEseal("ALICE Performance", 0.75, 0.75);
  
}

//------------------------------------------------------------------------
void ALICEseal(TString type, Double_t xPad, Double_t yPad)
{
  TVirtualPad* currPad = gPad;
  TPad *myPadLogo = new TPad("myPadLogo", "Pad for ALICE Logo",xPad,yPad,xPad+0.17,yPad+0.17);
  myPadLogo->SetFillColor(0);
  myPadLogo->SetBorderMode(0);
  myPadLogo->SetBorderSize(2);
  myPadLogo->SetFrameBorderMode(0);
  myPadLogo->SetLeftMargin(0.0);
  myPadLogo->SetTopMargin(0.0);
  myPadLogo->SetBottomMargin(0.0);
  myPadLogo->SetRightMargin(0.0);
  myPadLogo->Draw();
  myPadLogo->cd();
  TASImage *myAliceLogo = new TASImage("/Users/pillot/Pictures/alice_logo.png");
  myAliceLogo->Draw();
  currPad->cd();
  Double_t x1 = xPad - 0.07, y1 = yPad - 0.06;
  Double_t x2 = x1 + 0.25, y2 = y1 + 0.08;
  TPaveText* t1=new TPaveText(x1,y1,x2,y2,"NDC");
  t1->SetFillStyle(0);
  t1->SetBorderSize(0);
  t1->AddText(0.,0.,Form("%s", type.Data()));
  t1->SetTextColor(kRed);
  t1->SetTextFont(42);
  t1->Draw();
  TPaveText* t2=new TPaveText(x1+0.02,y1-0.06,x2-0.02,y2-0.06,"NDC");
  t2->SetFillStyle(0);
  t2->SetBorderSize(0);
  t2->SetTextColor(kRed);
  t2->SetTextFont(52);
  TDatime dt;
  TString today = Form("%02i/%02i/%4i", dt.GetDay(), dt.GetMonth(), dt.GetYear());
  t2->AddText(0.,0.,today.Data());
  t2->AddText(0.,0.,"Pb-Pb #sqrt{s_{NN}} = 2.76 TeV");
  t2->Draw();
}
