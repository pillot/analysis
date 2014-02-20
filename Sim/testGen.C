/*
 *  testGen.C
 *  aliroot_dev
 *
 *  Created by philippe pillot on 21/04/11.
 *  Copyright 2011 Subatech. All rights reserved.
 *
 */

#include "TSystem.h"
#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"

#include "EVGEN/AliGenLib.h"
#include "EVGEN/AliGenMUONlib.h"


typedef Double_t (*GenFunc) (const Double_t*,  const Double_t*);

GenFunc fptPbPb;
Double_t intPtPbPb;
GenFunc fptpp;
Double_t intPtpp;
GenFunc fyPbPb;
Double_t intYPbPb;
GenFunc fypp;
Double_t intYpp;

Double_t fptratio(const Double_t *x,  const Double_t *dummy) {
  Double_t pbpbval = fptPbPb(x,dummy);
  return (pbpbval > 0.) ? intPtPbPb * fptpp(x,dummy) / pbpbval / intPtpp : 0.;
}

Double_t fyratio(const Double_t *x,  const Double_t *dummy) {
  Double_t pbpbval = fyPbPb(x,dummy);
  return (pbpbval > 0.) ? intYPbPb * fypp(x,dummy) / pbpbval / intYpp : 0.;
}

void testGen()
{
  
  // to use it:
  /*
  gSystem->SetIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/include\ -I$ALICE_ROOT -I$ALICE/geant3/TGeant3");
   .x testGen.C+
  */
  
  AliGenLib* pLibrary = new AliGenMUONlib();
  
  fptPbPb = pLibrary->GetPt(AliGenMUONlib::kJpsi,"CDF PbPb 3.94");
  TF1 *ptPbPb = new TF1("ptPbPb", fptPbPb, 0., 10., 0);
  intPtPbPb  = ptPbPb->Integral(0., 10., (Double_t*) 0x0, 1.e-6);
  
  fyPbPb = pLibrary->GetY(AliGenMUONlib::kJpsi,"CDF PbPb 3.94");
  TF1 *yPbPb = new TF1("yPbPb", fyPbPb, -4.2, -2.3, 0);
  intYPbPb  = yPbPb->Integral(-4.2, -2.3, (Double_t*) 0x0, 1.e-6);
  
  fptpp = pLibrary->GetPt(AliGenMUONlib::kJpsi,"PbPb 2.76c11");
  TF1 *ptpp = new TF1("ptPbPb", fptpp, 0., 10., 0);
  intPtpp  = ptpp->Integral(0., 10., (Double_t*) 0x0, 1.e-6);
  
  fypp = pLibrary->GetY(AliGenMUONlib::kJpsi,"PbPb 2.76c11");
  TF1 *ypp = new TF1("yPbPb", fypp, -4.2, -2.3, 0);
  intYpp  = ypp->Integral(-4.2, -2.3, (Double_t*) 0x0, 1.e-6);
  
  TF1 *ptratio = new TF1("ptratio", fptratio, 0., 10., 0);
  
  TF1 *yratio = new TF1("yratio", fyratio, -4.2, -2.3, 0);
  
  TCanvas *c0 = new TCanvas("c0","Canvas 0",400,10,600,700);
  c0->Divide(2,2);
  c0->cd(1);
  gPad->SetLogy();
  ptPbPb->SetLineColor(kBlue);
  ptPbPb->Draw();
  ptPbPb->GetHistogram()->SetXTitle("p_{T} (GeV)");     
  ptPbPb->GetHistogram()->SetYTitle("dN/dpt (a.u.)");     
  ptpp->SetLineColor(kRed);
  ptpp->Draw("same");
  c0->cd(2);
  gPad->SetLogy();
  yPbPb->SetLineColor(kBlue);
  yPbPb->Draw();
  yPbPb->GetHistogram()->SetXTitle("y");   
  yPbPb->GetHistogram()->SetYTitle("dN/dy (a.u.)");     
  ypp->SetLineColor(kRed);
  ypp->Draw("same");
  c0->cd(3);
  ptratio->GetHistogram()->SetXTitle("p_{T} (GeV)");     
  ptratio->GetHistogram()->SetYTitle("ratio (a.u.)");     
  ptratio->SetLineColor(1);
  ptratio->Draw();
  c0->cd(4);
  yratio->GetHistogram()->SetXTitle("p_{T} (GeV)");     
  yratio->GetHistogram()->SetYTitle("ratio (a.u.)");     
  yratio->SetLineColor(1);
  yratio->Draw();
  
}