/*
 *  JPsiAODAnalysis.C
 *  aliroot_dev
 *
 *  Created by philippe pillot on 12/10/11.
 *  Copyright 2011 SUBATECH. All rights reserved.
 *
 */

#if !defined(__CINT__) || defined(__MAKECINT__)

#include <Riostream.h>
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH3F.h"
#include "TF1.h"
#include "TLegend.h"
#include "TMath.h"
#include "TStyle.h"
#include "TGraphAsymmErrors.h"

#include "AliAODEvent.h"
#include "AliAODMCParticle.h"

#endif

// pt/y bining
const Int_t nPtBins = 2;
Double_t pTBinLowEdge[nPtBins+1] = {0., 3., 8.};
const Int_t nYBins = 2;
Double_t yBinLowEdge[nYBins+1] = {-4., -3.25, -2.5};


Double_t CrystalBall(Double_t *x, Double_t *par);
Double_t CrystalBallExtended(Double_t *x, Double_t *par);


//---------------------------------------------------------------------------
void JPsiAODAnalysis(TString runList, Bool_t useCB2 = kTRUE, Bool_t fixTails = kFALSE)
{
  /// Fit the mass ditribution in JPsi MC simulation
  /// Plot the evolution of the parameters versus pt/y
  
  gROOT->SetStyle("Plain");
  gROOT->ForceStyle();
  gStyle->SetOptStat(10);
  gStyle->SetOptFit(2222);
  
  // mass distribution versus pT/y bins
  TH3F *hMass = new TH3F("hMass","J/#psi mass distribution versus pt/y bins;Mass (GeV/c^{2});p_T bin;y bin",
			 60, 0., 6., nPtBins+1, -0.5, nPtBins+0.5, nYBins+1, -0.5, nYBins+0.5);
  // set pt bin labels
  for (Int_t ipt = 0; ipt <= nPtBins; ipt++) {
    TString label = (ipt == 0) ? Form("%g<pt<%g",pTBinLowEdge[0],pTBinLowEdge[nPtBins]) : Form("%g<pt<%g",pTBinLowEdge[ipt-1],pTBinLowEdge[ipt]);
    hMass->GetYaxis()->SetBinLabel(ipt+1, label.Data());
  }
  // set y bin labels
  for (Int_t iy = 0; iy <= nYBins; iy++) {
    TString label = (iy == 0) ? Form("%g<y<%g",yBinLowEdge[0],yBinLowEdge[nYBins]) : Form("%g<y<%g",yBinLowEdge[iy-1],yBinLowEdge[iy]);
    hMass->GetZaxis()->SetBinLabel(iy+1, label.Data());
  }
  
  // fitting functions
  TF1 *fCB = (useCB2) ? new TF1("fCB", CrystalBallExtended, 0, 6, 7) : new TF1("fCB", CrystalBall, 0, 6, 5);
  
  // open run list
  ifstream inFile(runList.Data());
  if (!inFile.is_open()) {
    printf("cannot open file %s\n",runList.Data());
    return;
  }
  
  // loop over runs
  Int_t irun = -1;
  TString currRun;
  TList runs;
  runs.SetOwner();
  while (!inFile.eof()) {
    
    // get current run number
    currRun.ReadLine(inFile,kTRUE);
    if(currRun.IsNull()) continue;
    runs.AddLast(new TObjString(currRun));
    irun++;
    printf("run %s:\n", currRun.Data());
    
    // get input hists
    TFile *file = new TFile(Form("runs/%d/AliAOD.Muons.root",currRun.Atoi()), "read");
    if (!file || !file->IsOpen()) {
      printf("cannot open file runs/%d/AliAOD.Muons.root\n",currRun.Atoi());
      delete file;
      continue;
    }
    TTree *tAOD = (TTree*) file->Get("aodTree");
    
    // loop over AOD events
    AliAODEvent *aodEv = new AliAODEvent();
    aodEv->ReadFromTree(tAOD);
    Int_t nEvents = 0;
    while(tAOD->GetEvent(nEvents++)) {
      
      if(nEvents%10000 == 0) printf("Event %d\n",nEvents);
      
      // loop over reco dimuons
      Int_t ndimu = aodEv->GetNDimuons();
      for(Int_t q=0; q<ndimu; q++){
	
	AliAODDimuon *dimu = aodEv->GetDimuon(q);
	
	if(dimu->Charge()==0 && dimu->Y()>yBinLowEdge[0] && dimu->Y()<yBinLowEdge[nYBins] &&
	   dimu->GetMu(0)->Eta()>-4. && dimu->GetMu(0)->Eta()<-2.5 &&
	   dimu->GetMu(1)->Eta()>-4. && dimu->GetMu(1)->Eta()<-2.5 &&
	   dimu->GetMu(0)->GetRAtAbsorberEnd()>17.6 && dimu->GetMu(0)->GetRAtAbsorberEnd()<89.5 &&
	   dimu->GetMu(1)->GetRAtAbsorberEnd()>17.6 && dimu->GetMu(1)->GetRAtAbsorberEnd()<89.5 &&
	   dimu->GetMu(0)->GetMatchTrigger()>1 && dimu->GetMu(1)->GetMatchTrigger()>1
	   && dimu->GetMu(0)->GetLabel() >= 0 && dimu->GetMu(1)->GetLabel() >= 0
	   ){
	  
	  // pt/y integrated
	  hMass->Fill(dimu->M(), 0., 0.);
	  
	  // y integrated / pt bins
	  for (Int_t ipt = 0; ipt < nPtBins; ipt++) {
	    
	    if (dimu->Pt() > pTBinLowEdge[ipt] && dimu->Pt() < pTBinLowEdge[ipt+1]) {
	      
	      hMass->Fill(dimu->M(), ipt+1., 0.);
	      
	    }
	    
	  }
	  
	  // y bins
	  for (Int_t iy = 0; iy < nYBins; iy++) {
	    
	    if (dimu->Y() > yBinLowEdge[iy] && dimu->Y() < yBinLowEdge[iy+1]) {
	      
	      // pt integrated
	      hMass->Fill(dimu->M(), 0., iy+1.);
	      
	      // pt bins
	      for (Int_t ipt = 0; ipt < nPtBins; ipt++) {
		
		if (dimu->Pt() > pTBinLowEdge[ipt] && dimu->Pt() < pTBinLowEdge[ipt+1]) {
		  
		  hMass->Fill(dimu->M(), ipt+1., iy+1.);
		  
		}
		
	      }
	      
	    }
	    
	  }
	  
	}
	
      }
      
    }
    
    file->Close();
  }
  inFile.close();
  
  // display mass distribution for each pt/y bins
  TCanvas* cMass = new TCanvas("cMass", "Mass versus pT/y bins", 900, 900);
  cMass->Divide(nPtBins+1, nYBins+1);
  for (Int_t ipt = 0; ipt <= nPtBins; ipt++) {
    for (Int_t iy = 0; iy <= nYBins; iy++) {
      cMass->cd(ipt*(nYBins+1) + iy + 1);
      gPad->SetLogy();
      TH1D* p = hMass->ProjectionX(Form("hgen_pt%d_y%d",ipt,iy), ipt+1, ipt+1, iy+1, iy+1, "e");
      p->SetTitle(Form("%s / %s", hMass->GetYaxis()->GetBinLabel(ipt+1), hMass->GetZaxis()->GetBinLabel(iy+1)));
      Double_t nEntries = p->GetEntries();
      if (useCB2) fCB->SetParameters(nEntries, 3., 0.08, 1., 5., 1., 5.);
      else fCB->SetParameters(p->GetEntries(), 3., 0.08, 1., 5.);
      if (fixTails) {
	if (useCB2) {
	  fCB->FixParameter(3, 1.05);
	  fCB->FixParameter(4, 4.2);
	  fCB->FixParameter(5, 1.84);
	  fCB->FixParameter(6, 3.2);
	} else {
	  fCB->FixParameter(3, 0.98);
	  fCB->FixParameter(4, 5.2);
	}
      }
      p->Fit(fCB, "RIM");
      Double_t nJPsi = fCB->Integral(0., 6.);
      Double_t binWidth = p->GetBinWidth(1);
      printf("%26s: chi2/ndf = %4.2f, sigma = %5.3f ± %5.3f, NJPsi = %d ± %d, Entries = %d (%+3.1f%%)\n",
	     p->GetTitle(), fCB->GetChisquare()/fCB->GetNDF(), fCB->GetParameter(2), fCB->GetParError(2),
	     (Int_t)(nJPsi/binWidth), (Int_t)(fCB->IntegralError(0., 6.)/binWidth),
	     (Int_t)nEntries, 100.*(nJPsi/binWidth-nEntries)/nEntries);
    }
  }
  
}


//---------------------------------------------------------------------------
void ShowCBs()
{
  /// display the different CB parametrization used in the analysis
  
  gROOT->SetStyle("Plain");
  gROOT->ForceStyle();
  
  TLegend *l = new TLegend(0.12,0.65,0.52,0.89);
/*  TF1 *fCB1a = new TF1("fCB1a", CrystalBall, 2, 4, 5);
  fCB1a->FixParameter(0, 1.);
  fCB1a->FixParameter(1, 3.1);
  fCB1a->FixParameter(2, 0.08);
  fCB1a->FixParameter(3, 0.98);
  fCB1a->FixParameter(4, 5.2);
  l->AddEntry(fCB1a,"#alpha=0.98, n=5.2","l");
  TF1 *fCB1b = new TF1("fCB1b", CrystalBall, 2, 4, 5);
  fCB1b->FixParameter(0, 1.);
  fCB1b->FixParameter(1, 3.1);
  fCB1b->FixParameter(2, 0.08);
  fCB1b->FixParameter(3, 1.15);
  fCB1b->FixParameter(4, 1.59);
  l->AddEntry(fCB1b,"#alpha=1.15, n=1.59","l");
  TF1 *fCB1c = new TF1("fCB1c", CrystalBall, 2, 4, 5);
  fCB1c->FixParameter(0, 1.);
  fCB1c->FixParameter(1, 3.1);
  fCB1c->FixParameter(2, 0.08);
  fCB1c->FixParameter(3, 1.15);
  fCB1c->FixParameter(4, 3.6);
  l->AddEntry(fCB1c,"#alpha=1.15, n=3.6","l");
*/  TF1 *fCB2a = new TF1("fCB2a", CrystalBallExtended, 2, 4, 7);
  fCB2a->FixParameter(0, 1.);
  fCB2a->FixParameter(1, 3.1);
  fCB2a->FixParameter(2, 0.07);
  fCB2a->FixParameter(3, 1.01);
  fCB2a->FixParameter(4, 4.86);
  fCB2a->FixParameter(5, 2.18);
  fCB2a->FixParameter(6, 2.68);
//  l->AddEntry(fCB2a,Form("#alpha=%4.2f|%4.2f, n=%4.2f|%4.2f (embedding)",fCB2a->GetParameter(3),
//			 fCB2a->GetParameter(5),fCB2a->GetParameter(4),fCB2a->GetParameter(6)),"l");
  TF1 *fCB2a2 = new TF1("fCB2a", CrystalBallExtended, 2, 4, 7);
  fCB2a2->FixParameter(0, 1.);
  fCB2a2->FixParameter(1, 3.1);
  fCB2a2->FixParameter(2, 0.07);
  fCB2a2->FixParameter(3, 1.02015e+00);
  fCB2a2->FixParameter(4, 4.95025e+00);
  fCB2a2->FixParameter(5, 2.13755e+00);
  fCB2a2->FixParameter(6, 2.81845e+00);
  l->AddEntry(fCB2a2,Form("#alpha=%4.2f|%4.2f, n=%4.2f|%4.2f (embedding)",fCB2a2->GetParameter(3),
			 fCB2a2->GetParameter(5),fCB2a2->GetParameter(4),fCB2a2->GetParameter(6)),"l");
  TF1 *fCB2b = new TF1("fCB2b", CrystalBallExtended, 2, 4, 7);
  fCB2b->FixParameter(0, 1.);
  fCB2b->FixParameter(1, 3.1);
  fCB2b->FixParameter(2, 0.07);
  fCB2b->FixParameter(3, 0.993713);
  fCB2b->FixParameter(4, 6.06265);
  fCB2b->FixParameter(5, 2.20168);
  fCB2b->FixParameter(6, 2.9975);
//  l->AddEntry(fCB2b,Form("#alpha=%4.2f|%4.2f, n=%4.2f|%4.2f (pure geo raw)",fCB2b->GetParameter(3),
//			 fCB2b->GetParameter(5),fCB2b->GetParameter(4),fCB2b->GetParameter(6)),"l");
  TF1 *fCB2b2 = new TF1("fCB2b2", CrystalBallExtended, 2, 4, 7);
  fCB2b2->FixParameter(0, 1.);
  fCB2b2->FixParameter(1, 3.1);
  fCB2b2->FixParameter(2, 0.07);
  fCB2b2->FixParameter(3, 1.03604e+00);
  fCB2b2->FixParameter(4, 5.08201e+00);
  fCB2b2->FixParameter(5, 2.21576e+00);
  fCB2b2->FixParameter(6, 2.71227e+00);
  l->AddEntry(fCB2b2,Form("#alpha=%4.2f|%4.2f, n=%4.2f|%4.2f (pure geo raw)",fCB2b2->GetParameter(3),
			 fCB2b2->GetParameter(5),fCB2b2->GetParameter(4),fCB2b2->GetParameter(6)),"l");
  TF1 *fCB2c = new TF1("fCB2c", CrystalBallExtended, 2, 4, 7);
  fCB2c->FixParameter(0, 1.);
  fCB2c->FixParameter(1, 3.1);
  fCB2c->FixParameter(2, 0.07);
  fCB2c->FixParameter(3, 1.01199);
  fCB2c->FixParameter(4, 6.29919);
  fCB2c->FixParameter(5, 2.15545);
  fCB2c->FixParameter(6, 2.52695);
//  l->AddEntry(fCB2c,Form("#alpha=%4.2f|%4.2f, n=%4.2f|%4.2f (pure geo ideal)",fCB2c->GetParameter(3),
//			 fCB2c->GetParameter(5),fCB2c->GetParameter(4),fCB2c->GetParameter(6)),"l");
  TF1 *fCB2c2 = new TF1("fCB2c2", CrystalBallExtended, 2, 4, 7);
  fCB2c2->FixParameter(0, 1.);
  fCB2c2->FixParameter(1, 3.1);
  fCB2c2->FixParameter(2, 0.07);
  fCB2c2->FixParameter(3, 1.05620);
  fCB2c2->FixParameter(4, 5.15412);
  fCB2c2->FixParameter(5, 2.08599);
  fCB2c2->FixParameter(6, 2.76698);
//  l->AddEntry(fCB2c2,Form("#alpha=%4.2f|%4.2f, n=%4.2f|%4.2f (pure geo ideal)",fCB2c2->GetParameter(3),
//			 fCB2c2->GetParameter(5),fCB2c2->GetParameter(4),fCB2c2->GetParameter(6)),"l");
  TF1 *fCB2c3 = new TF1("fCB2c3", CrystalBallExtended, 2, 4, 7);
  fCB2c3->FixParameter(0, 1.);
  fCB2c3->FixParameter(1, 3.1);
  fCB2c3->FixParameter(2, 0.07);
  fCB2c3->FixParameter(3, 1.06160);
  fCB2c3->FixParameter(4, 5.07988);
  fCB2c3->FixParameter(5, 2.17418);
  fCB2c3->FixParameter(6, 2.89152);
  l->AddEntry(fCB2c3,Form("#alpha=%4.2f|%4.2f, n=%4.2f|%4.2f (pure geo ideal)",fCB2c3->GetParameter(3),
			 fCB2c3->GetParameter(5),fCB2c3->GetParameter(4),fCB2c3->GetParameter(6)),"l");
  TF1 *fCB2d = new TF1("fCB2d", CrystalBallExtended, 2, 4, 7);
  fCB2d->FixParameter(0, 1.);
  fCB2d->FixParameter(1, 3.1);
  fCB2d->FixParameter(2, 0.07);
  fCB2d->FixParameter(3, 1.12222);
  fCB2d->FixParameter(4, 2.44219);
  fCB2d->FixParameter(5, 1.91267);
  fCB2d->FixParameter(6, 10.);
  l->AddEntry(fCB2d,Form("#alpha=%4.2f|%4.2f, n=%4.2f|%4.2f (pp [1.8-5])",fCB2d->GetParameter(3),
			 fCB2d->GetParameter(5),fCB2d->GetParameter(4),fCB2d->GetParameter(6)),"l");
  TF1 *fCB2e = new TF1("fCB2e", CrystalBallExtended, 2, 4, 7);
  fCB2e->FixParameter(0, 1.);
  fCB2e->FixParameter(1, 3.1);
  fCB2e->FixParameter(2, 0.07);
  fCB2e->FixParameter(3, 0.998627);
  fCB2e->FixParameter(4, 6.06996);
  fCB2e->FixParameter(5, 1.91119);
  fCB2e->FixParameter(6, 10.);
  l->AddEntry(fCB2e,Form("#alpha=%4.2f|%4.2f, n=%4.2f|%4.2f (pp [2.0-5])",fCB2e->GetParameter(3),
			 fCB2e->GetParameter(5),fCB2e->GetParameter(4),fCB2e->GetParameter(6)),"l");
  TF1 *fCB2f = new TF1("fCB2f", CrystalBallExtended, 2, 4, 7);
  fCB2f->FixParameter(0, 1.);
  fCB2f->FixParameter(1, 3.1);
  fCB2f->FixParameter(2, 0.07);
  fCB2f->FixParameter(3, 0.97417);
  fCB2f->FixParameter(4, 8.55141);
  fCB2f->FixParameter(5, 1.92676);
  fCB2f->FixParameter(6, 10.);
  l->AddEntry(fCB2f,Form("#alpha=%4.2f|%4.2f, n=%4.2f|%4.2f (pp [2.2-5])",fCB2f->GetParameter(3),
			 fCB2f->GetParameter(5),fCB2f->GetParameter(4),fCB2f->GetParameter(6)),"l");
  TF1 *fCB2g = new TF1("fCB2g", CrystalBallExtended, 2, 4, 7);
  fCB2g->FixParameter(0, 1.);
  fCB2g->FixParameter(1, 3.1);
  fCB2g->FixParameter(2, 0.07);
  fCB2g->FixParameter(3, 8.75435e-01);
  fCB2g->FixParameter(4, 9.57296e+00);
  fCB2g->FixParameter(5, 2.32829e+00);
  fCB2g->FixParameter(6, 2.84961e+00);
  l->AddEntry(fCB2g,Form("#alpha=%4.2f|%4.2f, n=%4.2f|%4.2f (UPC)",fCB2g->GetParameter(3),
			 fCB2g->GetParameter(5),fCB2g->GetParameter(4),fCB2g->GetParameter(6)),"l");
  
  new TCanvas;
  gPad->SetLogy();
/*  fCB1a->SetLineColor(2);
  fCB1a->SetLineStyle(1);
  fCB1a->GetYaxis()->SetRangeUser(0.0003, 3.);
  fCB1a->Draw();
  fCB1b->SetLineColor(8);
  fCB1b->SetLineStyle(2);
  fCB1b->Draw("same");
  fCB1c->SetLineColor(4);
  fCB1c->SetLineStyle(3);
  fCB1c->Draw("same");
*/
/*  fCB2a->SetLineColor(4);
  fCB2a->SetLineStyle(1);
  fCB2a->GetYaxis()->SetRangeUser(0.0003, 3.);
  fCB2a->Draw();
*/  fCB2a2->SetLineColor(4);
  fCB2a2->GetYaxis()->SetRangeUser(0.0003, 3.);
  fCB2a2->SetLineStyle(1);
  fCB2a2->Draw();
/*  fCB2b->SetLineColor(2);
  fCB2b->SetLineStyle(2);
  fCB2b->Draw("same");
*/  fCB2b2->SetLineColor(2);
  fCB2b2->SetLineStyle(2);
  fCB2b2->Draw("same");
/*  fCB2c->SetLineColor(6);
  fCB2c->SetLineStyle(3);
  fCB2c->Draw("same");
  fCB2c2->SetLineColor(2);
  fCB2c2->SetLineStyle(3);
  fCB2c2->Draw("same");
*/  fCB2c3->SetLineColor(6);
  fCB2c3->SetLineStyle(3);
  fCB2c3->Draw("same");
  fCB2d->SetLineColor(15);
  fCB2d->SetLineStyle(4);
  fCB2d->Draw("same");
  fCB2e->SetLineColor(8);
  fCB2e->SetLineStyle(5);
  fCB2e->Draw("same");
  fCB2f->SetLineColor(9);
  fCB2f->SetLineStyle(6);
  fCB2f->Draw("same");
  fCB2g->SetLineColor(1);
  fCB2g->SetLineStyle(7);
  fCB2g->Draw("same");
  l->SetMargin(0.1);
  l->SetFillStyle(0);
  l->SetBorderSize(0);
  l->Draw("same");
}


//---------------------------------------------------------------------------
Double_t CrystalBall(Double_t *x, Double_t *par)
{
  //par[0] = Normalization
  //par[1] = mean
  //par[2] = sigma
  //par[3] = alpha
  //par[4] = n 
  
  Double_t t = (x[0]-par[1])/par[2];
  if (par[3] < 0) t = -t;
  
  Double_t absAlpha = fabs((Double_t)par[3]);
  
  if (t >= -absAlpha) // gaussian core
  {
    return par[0]*(exp(-0.5*t*t));
  }
  
  if (t < -absAlpha) //left tail
  {
    Double_t a =  TMath::Power(par[4]/absAlpha,par[4])*exp(-0.5*absAlpha*absAlpha);
    Double_t b = par[4]/absAlpha - absAlpha;
    return par[0]*(a/TMath::Power(b - t, par[4]));
  }
  
  return 0. ; 
} 


//---------------------------------------------------------------------------
Double_t CrystalBallExtended(Double_t *x, Double_t *par)
{
  //par[0] = Normalization
  //par[1] = mean
  //par[2] = sigma
  //par[3] = alpha
  //par[4] = n
  //par[5] = alpha'
  //par[6] = n'  
  
  Double_t t = (x[0]-par[1])/par[2];
  if (par[3] < 0) t = -t;
  
  Double_t absAlpha = fabs((Double_t)par[3]);
  Double_t absAlpha2 = fabs((Double_t)par[5]);
  
  if (t >= -absAlpha && t < absAlpha2) // gaussian core
  {
    return par[0]*(exp(-0.5*t*t));
  }
  
  if (t < -absAlpha) //left tail
  {
    Double_t a =  TMath::Power(par[4]/absAlpha,par[4])*exp(-0.5*absAlpha*absAlpha);
    Double_t b = par[4]/absAlpha - absAlpha;
    return par[0]*(a/TMath::Power(b - t, par[4]));
  }
  
  if (t >= absAlpha2) //right tail
  {
    
    Double_t c =  TMath::Power(par[6]/absAlpha2,par[6])*exp(-0.5*absAlpha2*absAlpha2);
    Double_t d = par[6]/absAlpha2 - absAlpha2;
    return par[0]*(c/TMath::Power(d + t, par[6]));
  }
  
  return 0. ; 
} 
