#if !defined(__CINT__) || defined(__MAKECINT__)

#include <Riostream.h>
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "THnSparse.h"
#include "TMath.h"
#include "TStyle.h"
#include "TGraphAsymmErrors.h"

#include "AliAODEvent.h"
#include "AliAODMCParticle.h"

#endif

// number of dimensions for efficiency histograms
const Int_t nDim = 4;

// pt/y bining
const Int_t nPtBins = 2;
Double_t pTBinLowEdge[nPtBins+1] = {0., 3., 1000000.};
const Int_t nYBins = 2;
Double_t yBinLowEdge[nYBins+1] = {-4., -3.25, -2.5};

// centrality bining
const Int_t nCentBins = 12;
Double_t centRange[2] = {-10., 110.};

// summary histograms
TH1F* hGenSummary = 0x0;
TH1F* hRecSummary = 0x0;
TH1F* hAccSummary = 0x0;

// histograms vs. pt (J/Psi)
TH1F* hPtGen = 0x0;
TH1F* hPtRec = 0x0;

// histograms vs. y (J/Psi)
TH1F* hYGen = 0x0;
TH1F* hYRec = 0x0;

// histograms vs. pt (single mu)
TH1F* hPtGenMu = 0x0;
TH1F* hPtRecMu = 0x0;

// histograms vs. y (single mu)
TH1F* hYGenMu = 0x0;
TH1F* hYRecMu = 0x0;

Int_t FillHisto(TTree *tAOD, Int_t irun, THnSparse* hgen, THnSparse* hrec, Bool_t print);
void  ComputeAccEff(TList& runs, Int_t ptBin, Int_t yBin, TH1D* hgen, TH1D* hrec);

//---------------------------------------------------------------------------
void ValuesAccEff(TString runList){
  
  //----STYLE----
  gStyle->SetPalette(1);
  gStyle->SetDrawOption(0);
  gStyle->SetFillColor(0);
  gStyle->SetOptStat(11);
  
  Int_t nBins[nDim] = {1003, nPtBins+1, nYBins+1, nCentBins};
  Double_t xMin[nDim] = {-0.5, -0.5, -0.5, centRange[0]};
  Double_t xMax[nDim] = {1002.5, nPtBins+0.5, nYBins+0.5, centRange[1]};
  THnSparse *hgen = new THnSparseT<TArrayF>("hgen", "hgen", nDim, nBins, xMin, xMax);
  THnSparse *hrec = new THnSparseT<TArrayF>("hrec", "hrec", nDim, nBins, xMin, xMax);
  
  hGenSummary = new TH1F("hGenSummary","generated J/#psi versus pt/y bins;;a.u.", (nPtBins+1)*(nYBins+1), -0.5, (nPtBins+1)*(nYBins+1)-0.5);
  hRecSummary = new TH1F("hRecSummary","reconstructed J/#psi versus pt/y bins;;a.u.", (nPtBins+1)*(nYBins+1), -0.5, (nPtBins+1)*(nYBins+1)-0.5);
  hAccSummary = new TH1F("hAccSummary","Acc * Eff versus pt/y bins", (nPtBins+1)*(nYBins+1), -0.5, (nPtBins+1)*(nYBins+1)-0.5);
  
  hPtGen = new TH1F("hPtGen","generated J/#psi versus pt;p_{T} (GeV/c);N_{J/#psi} / 0.5 GeV/c", 30, 0., 30.);
  hPtRec = new TH1F("hPtRec","reconstructed J/#psi versus pt;p_{T} (GeV/c);N_{J/#psi} / 0.5 GeV/c", 30, 0., 30.);
  
  hYGen = new TH1F("hYGen","generated J/#psi versus y;y;N_{J/#psi}", 15, -4., -2.5);
  hYRec = new TH1F("hYRec","reconstructed J/#psi versus y;y;N_{J/#psi}", 15, -4., -2.5);
  
  hPtGenMu = new TH1F("hPtGenMu","generated single-#mu versus pt;p_{T} (GeV/c);N_{single-#mu} / 0.5 GeV/c", 30, 0., 30.);
  hPtRecMu = new TH1F("hPtRecMu","reconstructed single-#mu versus pt;p_{T} (GeV/c);N_{single-#mu} / 0.5 GeV/c", 30, 0., 30.);
  
  hYGenMu = new TH1F("hYGenMu","generated single-#mu versus y;y;N_{single-#mu}", 15, -4., -2.5);
  hYRecMu = new TH1F("hYRecMu","reconstructed single-#mu versus y;y;N_{single-#mu}", 15, -4., -2.5);
  
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
  Int_t nEvents = 0;
  while (!inFile.eof()) {
    
    // get current run number
    currRun.ReadLine(inFile,kTRUE);
    if(currRun.IsNull()) continue;
    runs.AddLast(new TObjString(currRun));
    irun++;
    printf("run %s:\n", currRun.Data());
    
    // get input hists
    TFile *file = new TFile(Form("runs/%d/AliAOD.root",currRun.Atoi()), "read");
    if (!file || !file->IsOpen()) {
      printf("cannot open file runs/%d/AliAOD.root\n",currRun.Atoi());
      delete file;
      continue;
    }
    TTree *tAOD = (TTree*) file->Get("aodTree");
    
    // fill histo
    nEvents += FillHisto(tAOD, irun, hgen, hrec, kTRUE);
    
    file->Close();
  }
  inFile.close();
  
  // compute integrated acc*eff correction
  for (Int_t ipt = 0; ipt <= nPtBins; ipt++) {
    for (Int_t iy = 0; iy <= nYBins; iy++) {
      hgen->GetAxis(1)->SetRange(ipt+1,ipt+1);
      hgen->GetAxis(2)->SetRange(iy+1,iy+1);
      //hgen->GetAxis(3)->SetRange(hgen->GetAxis(3)->FindBin(1.),hgen->GetAxis(3)->FindBin(79.));
      //hgen->GetAxis(3)->SetRange(hgen->GetAxis(3)->FindBin(51.),hgen->GetAxis(3)->FindBin(79.));
      TH1D* hg = hgen->Projection(0,"e");
      hg->SetName(Form("hgen_pt%d_y%d",ipt,iy));
      hrec->GetAxis(1)->SetRange(ipt+1,ipt+1);
      hrec->GetAxis(2)->SetRange(iy+1,iy+1);
      //hrec->GetAxis(3)->SetRange(hrec->GetAxis(3)->FindBin(1.),hrec->GetAxis(3)->FindBin(79.));
      //hrec->GetAxis(3)->SetRange(hrec->GetAxis(3)->FindBin(51.),hrec->GetAxis(3)->FindBin(79.));
      TH1D* hr = hrec->Projection(0,"e");
      hr->SetName(Form("hrec_pt%d_y%d",ipt,iy));
      ComputeAccEff(runs, ipt, iy, hg, hr);
    }
  }
  
  // compute acc*eff correction versus centrality
  TGraphAsymmErrors *gCentAcc[nPtBins+1][nYBins+1];
  for (Int_t ipt = 0; ipt <= nPtBins; ipt++) {
    for (Int_t iy = 0; iy <= nYBins; iy++) {
      hgen->GetAxis(0)->SetRange(1003,1003);
      hgen->GetAxis(1)->SetRange(ipt+1,ipt+1);
      hgen->GetAxis(2)->SetRange(iy+1,iy+1);
      TH1D* hg = hgen->Projection(3,"e");
      hrec->GetAxis(0)->SetRange(1003,1003);
      hrec->GetAxis(1)->SetRange(ipt+1,ipt+1);
      hrec->GetAxis(2)->SetRange(iy+1,iy+1);
      TH1D* hr = hrec->Projection(3,"e");
      gCentAcc[ipt][iy] = new TGraphAsymmErrors(hr, hg);
      TString title = "J/#psi Acc * Eff versus centrality";
      title += (ipt == 0) ? Form(" (%g<pt<%g",pTBinLowEdge[0],pTBinLowEdge[nPtBins]) : Form(" (%g<pt<%g",pTBinLowEdge[ipt-1],pTBinLowEdge[ipt]);
      title += (iy == 0) ? Form(" / %g<y<%g)",yBinLowEdge[0],yBinLowEdge[nYBins]) : Form(" / %g<y<%g)",yBinLowEdge[iy-1],yBinLowEdge[iy]);
      gCentAcc[ipt][iy]->SetNameTitle(Form("gCentAcc_pt%d_y%d",ipt,iy), title.Data());
      delete hg;
      delete hr;
    }
  }
  
  // compute acc*eff vs. pt (J/Psi)
  TH1F* hPtAcc = new TH1F("hPtAcc","Acc * Eff versus pt (J/#psi)", hPtGen->GetNbinsX(), hPtGen->GetXaxis()->GetXmin(), hPtGen->GetXaxis()->GetXmax());
  hPtAcc->Sumw2();
  hPtAcc->Divide(hPtRec, hPtGen, 1., 1., "B");
  
  // compute acc*eff vs. y (J/Psi)
  TH1F* hYAcc = new TH1F("hYAcc","Acc * Eff versus y (J/#psi)", hYGen->GetNbinsX(), hYGen->GetXaxis()->GetXmin(), hYGen->GetXaxis()->GetXmax());
  hYAcc->Sumw2();
  hYAcc->Divide(hYRec, hYGen, 1., 1., "B");
  
  // compute acc*eff vs. pt (single mu)
  TH1F* hPtAccMu = new TH1F("hPtAccMu","Acc * Eff versus pt (single-#mu)", hPtGenMu->GetNbinsX(), hPtGenMu->GetXaxis()->GetXmin(), hPtGenMu->GetXaxis()->GetXmax());
  hPtAccMu->Sumw2();
  hPtAccMu->Divide(hPtRecMu, hPtGenMu, 1., 1., "B");
  
  // compute acc*eff vs. y (single mu)
  TH1F* hYAccMu = new TH1F("hYAccMu","Acc * Eff versus y (single-#mu)", hYGenMu->GetNbinsX(), hYGenMu->GetXaxis()->GetXmin(), hYGenMu->GetXaxis()->GetXmax());
  hYAccMu->Sumw2();
  hYAccMu->Divide(hYRecMu, hYGenMu, 1., 1., "B");
  
  // draw histos vs pt (J/Psi)
  TCanvas* cPt = new TCanvas("cPt", "versus pT (J/Psi)", 900, 300);
  cPt->Divide(3,1);
  cPt->cd(1);
  gPad->SetLogy();
  hPtGen->Draw();
  cPt->cd(2);
  gPad->SetLogy();
  hPtRec->Draw();
  cPt->cd(3);
  hPtAcc->SetMarkerStyle(21);
  hPtAcc->SetMarkerSize(0.6);
  hPtAcc->SetStats(kFALSE);
  hPtAcc->Draw("e0");
  
  // draw histos vs y (J/Psi)
  TCanvas* cY = new TCanvas("cY", "versus y (J/Psi)", 900, 300);
  cY->Divide(3,1);
  cY->cd(1);
  hYGen->Draw();
  cY->cd(2);
  hYRec->Draw();
  cY->cd(3);
  hYAcc->SetMarkerStyle(21);
  hYAcc->SetMarkerSize(0.6);
  hYAcc->SetStats(kFALSE);
  hYAcc->Draw("e0");
  
  // draw histos vs pt (single mu)
  TCanvas* cPtMu = new TCanvas("cPtMu", "versus pT (single-mu)", 900, 300);
  cPtMu->Divide(3,1);
  cPtMu->cd(1);
  gPad->SetLogy();
  hPtGenMu->Draw();
  cPtMu->cd(2);
  gPad->SetLogy();
  hPtRecMu->Draw();
  cPtMu->cd(3);
  hPtAccMu->SetMarkerStyle(21);
  hPtAccMu->SetMarkerSize(0.6);
  hPtAccMu->SetStats(kFALSE);
  hPtAccMu->Draw("e0");
  
  // draw histos vs y (single mu)
  TCanvas* cYMu = new TCanvas("cYMu", "versus y (single-mu)", 900, 300);
  cYMu->Divide(3,1);
  cYMu->cd(1);
  hYGenMu->Draw();
  cYMu->cd(2);
  hYRecMu->Draw();
  cYMu->cd(3);
  hYAccMu->SetMarkerStyle(21);
  hYAccMu->SetMarkerSize(0.6);
  hYAccMu->SetStats(kFALSE);
  hYAccMu->Draw("e0");
  
  // draw summary histos
  TCanvas* cSummary = new TCanvas("cSummary", "summary", 900, 300);
  cSummary->Divide(3,1);
  cSummary->cd(1);
  hGenSummary->Sumw2();
  hGenSummary->Scale(1./nEvents);
  hGenSummary->SetStats(kFALSE);
  hGenSummary->Draw("e0");
  cSummary->cd(2);
  hRecSummary->Sumw2();
  hRecSummary->Scale(1./nEvents);
  hRecSummary->SetStats(kFALSE);
  hRecSummary->Draw("e0");
  cSummary->cd(3);
  hAccSummary->SetMarkerStyle(21);
  hAccSummary->SetMarkerSize(0.6);
  hAccSummary->SetStats(kFALSE);
  hAccSummary->Draw("e0");
  
  // draw histos vs centrality
  TCanvas *cCentAcc[nPtBins+1][nYBins+1];
  for (Int_t ipt = 0; ipt <= nPtBins; ipt++) {
    for (Int_t iy = 0; iy <= nYBins; iy++) {
      cCentAcc[ipt][iy] = new TCanvas(Form("cCentAcc_pt%d_y%d",ipt,iy), gCentAcc[ipt][iy]->GetTitle(), 900, 300);
      gCentAcc[ipt][iy]->GetXaxis()->Set(nCentBins, centRange[0], centRange[1]);
      gCentAcc[ipt][iy]->Draw("ap");
    }
  }
  
  // save histos
  TFile* file = new TFile("acceff.root","update");
  hGenSummary->Write(0x0, TObject::kOverwrite);
  hRecSummary->Write(0x0, TObject::kOverwrite);
  hAccSummary->Write(0x0, TObject::kOverwrite);
  cSummary->Write(0x0, TObject::kOverwrite);
  hPtGen->Write(0x0, TObject::kOverwrite);
  hPtRec->Write(0x0, TObject::kOverwrite);
  hPtAcc->Write(0x0, TObject::kOverwrite);
  cPt->Write(0x0, TObject::kOverwrite);
  hYGen->Write(0x0, TObject::kOverwrite);
  hYRec->Write(0x0, TObject::kOverwrite);
  hYAcc->Write(0x0, TObject::kOverwrite);
  cY->Write(0x0, TObject::kOverwrite);
  hPtGenMu->Write(0x0, TObject::kOverwrite);
  hPtRecMu->Write(0x0, TObject::kOverwrite);
  hPtAccMu->Write(0x0, TObject::kOverwrite);
  cPtMu->Write(0x0, TObject::kOverwrite);
  hYGenMu->Write(0x0, TObject::kOverwrite);
  hYRecMu->Write(0x0, TObject::kOverwrite);
  hYAccMu->Write(0x0, TObject::kOverwrite);
  cYMu->Write(0x0, TObject::kOverwrite);
  for (Int_t ipt = 0; ipt <= nPtBins; ipt++) {
    for (Int_t iy = 0; iy <= nYBins; iy++) {
      gCentAcc[ipt][iy]->Write(0x0, TObject::kOverwrite);
      cCentAcc[ipt][iy]->Write(0x0, TObject::kOverwrite);
    }
  }
  file->Close();
}

//---------------------------------------------------------------------------
Int_t FillHisto(TTree *tAOD, Int_t irun, THnSparse* hgen, THnSparse* hrec, Bool_t print) {
  
  Int_t intgen[nPtBins+1][nYBins+1];
  for (Int_t i=0; i<=nPtBins; i++) for (Int_t j=0; j<=nYBins; j++) intgen[i][j] = 0;
  Int_t intrec[3][nPtBins+1][nYBins+1];
  for (Int_t i=0; i<3; i++) for (Int_t j=0; j<=nPtBins; j++) for (Int_t k=0; k<=nYBins; k++) intrec[i][j][k] = 0;
  Double_t bin[nDim];
  
  printf("Entries: %d\n", (Int_t)tAOD->GetEntries());
  
  // get AOD event
  AliAODEvent *aodEv = new AliAODEvent();
  aodEv->ReadFromTree(tAOD);
  TClonesArray *mcarray = dynamic_cast<TClonesArray*>(aodEv->FindListObject(AliAODMCParticle::StdBranchName()));
  Int_t nEvents = 0;
  while(tAOD->GetEvent(nEvents++)) {
    
    if(nEvents%10000 == 0) printf("Event %d\n",nEvents);
    /*
    // remove events with JPsi outside -4.2<y<-2.3
    Bool_t skip = kFALSE;
    for(int ii=0;ii<mcarray->GetEntries();ii++) {
      AliAODMCParticle *mctrack = (AliAODMCParticle*) mcarray->UncheckedAt(ii); 
      if (mctrack->GetPdgCode() == 443 && (mctrack->Y()<-4.2 || mctrack->Y()>-2.3)) {
	skip = kTRUE;
	break;
      }
    }
    if (skip) continue;
    */
    // get the centrality percentile
    Float_t centrality = aodEv->GetCentrality()->GetCentralityPercentileUnchecked("V0M");
    bin[3] = centrality;
    
    // select on centrality
    //if (centrality < 0. || centrality > 80.) continue;
    
    // loop over MC tracks
    for(int ii=0;ii<mcarray->GetEntries();ii++){
      
      AliAODMCParticle *mctrack = (AliAODMCParticle*) mcarray->At(ii); 
      
      // look for muons
      if (TMath::Abs(mctrack->GetPdgCode()) == 13 &&
	  mctrack->Eta()>yBinLowEdge[0] && mctrack->Eta()<yBinLowEdge[nYBins]) {
	
	hPtGenMu->Fill(mctrack->Pt());
	hYGenMu->Fill(mctrack->Y());
	
      }
      
      // look for generated particle
      if(mctrack->IsPrimary() && !mctrack->IsPhysicalPrimary() &&
	 mctrack->Y()>yBinLowEdge[0] && mctrack->Y()<yBinLowEdge[nYBins]) {
	
	// pt/y integrated
	intgen[0][0]++;
	bin[0] = irun; bin[1] = 0.; bin[2] = 0.;
	hgen->Fill(bin);
	for (Int_t itrg = 0; itrg < 3; itrg++) {
	  bin[0] = 1000+itrg;
	  hgen->Fill(bin);
	}
	hPtGen->Fill(mctrack->Pt());
	hYGen->Fill(mctrack->Y());
	
	// y integrated / pt bins
	for (Int_t ipt = 0; ipt < nPtBins; ipt++) {
	  
	  if (mctrack->Pt() >= pTBinLowEdge[ipt] && mctrack->Pt() < pTBinLowEdge[ipt+1]) {
	    
	    intgen[ipt+1][0]++;
	    bin[0] = irun; bin[1] = ipt+1; bin[2] = 0.;
	    hgen->Fill(bin);
	    for (Int_t itrg = 0; itrg < 3; itrg++) {
	      bin[0] = 1000+itrg;
	      hgen->Fill(bin);
	    }
	    
	  }
	  
	}
	
	// y bins
	for (Int_t iy = 0; iy < nYBins; iy++) {
	  
	  if (mctrack->Y() > yBinLowEdge[iy] && mctrack->Y() < yBinLowEdge[iy+1]) {
	    
	    // pt integrated
	    intgen[0][iy+1]++;
	    bin[0] = irun; bin[1] = 0.; bin[2] = iy+1;
	    hgen->Fill(bin);
	    for (Int_t itrg = 0; itrg < 3; itrg++) {
	      bin[0] = 1000+itrg;
	      hgen->Fill(bin);
	    }
	    
	    // pt bins
	    for (Int_t ipt = 0; ipt < nPtBins; ipt++) {
	      
	      if (mctrack->Pt() >= pTBinLowEdge[ipt] && mctrack->Pt() < pTBinLowEdge[ipt+1]) {
		
		intgen[ipt+1][iy+1]++;
		bin[0] = irun; bin[1] = ipt+1; bin[2] = iy+1;
		hgen->Fill(bin);
		for (Int_t itrg = 0; itrg < 3; itrg++) {
		  bin[0] = 1000+itrg;
		  hgen->Fill(bin);
		}
		
	      }
	      
	    }
	    
	  }
	  
	}
	
      }
      
    }
    
    // loop over reco tracks
    Int_t ntracks = aodEv->GetNTracks();
    for(Int_t q=0; q<ntracks; q++){
      
      AliAODTrack *mu = aodEv->GetTrack(q);
      
      if (mu->IsMuonTrack() && mu->GetMatchTrigger()>0 && mu->GetLabel() >= 0 &&
	  mu->Eta()>yBinLowEdge[0] && mu->Eta()<yBinLowEdge[nYBins] &&
	  mu->GetRAtAbsorberEnd()>17.6 && mu->GetRAtAbsorberEnd()<89.5
	  ) {
	
	hPtRecMu->Fill(mu->Pt());
	hYRecMu->Fill(mu->Y());
	
      }
      
    }
    
    // loop over reco dimuons
    Int_t ndimu = aodEv->GetNDimuons();
    for(Int_t q=0; q<ndimu; q++){
      
      AliAODDimuon *dimu = aodEv->GetDimuon(q);
      
      Double_t thetaTrackAbsEnd0 = TMath::ATan(dimu->GetMu(0)->GetRAtAbsorberEnd()/505.) * TMath::RadToDeg();
      Double_t thetaTrackAbsEnd1 = TMath::ATan(dimu->GetMu(1)->GetRAtAbsorberEnd()/505.) * TMath::RadToDeg();
      if(dimu->Charge()==0 && dimu->Y()>yBinLowEdge[0] && dimu->Y()<yBinLowEdge[nYBins] &&
	 dimu->GetMu(0)->Eta()>-4. && dimu->GetMu(0)->Eta()<-2.5 &&
	 dimu->GetMu(1)->Eta()>-4. && dimu->GetMu(1)->Eta()<-2.5 &&
	 //dimu->GetMu(0)->Y()>-4. && dimu->GetMu(0)->Y()<-2.5 &&
	 //dimu->GetMu(1)->Y()>-4. && dimu->GetMu(1)->Y()<-2.5 &&
	 //dimu->GetMu(0)->GetRAtAbsorberEnd()>17.622 && dimu->GetMu(0)->GetRAtAbsorberEnd()<89. &&
	 //dimu->GetMu(1)->GetRAtAbsorberEnd()>17.622 && dimu->GetMu(1)->GetRAtAbsorberEnd()<89. 
	 dimu->GetMu(0)->GetRAtAbsorberEnd()>17.6 && dimu->GetMu(0)->GetRAtAbsorberEnd()<89.5 &&
	 dimu->GetMu(1)->GetRAtAbsorberEnd()>17.6 && dimu->GetMu(1)->GetRAtAbsorberEnd()<89.5 
	 //thetaTrackAbsEnd0>2. && thetaTrackAbsEnd0<10. &&
	 //thetaTrackAbsEnd1>2. && thetaTrackAbsEnd1<10.
	 && dimu->GetMu(0)->GetLabel() >= 0 && dimu->GetMu(1)->GetLabel() >= 0
	 //&& dimu->GetMu(0)->Pt()>0.5 && dimu->GetMu(1)->Pt()>0.5
	 ){
	
	// loop over trigger cases
	Bool_t trigOk[3] = {kTRUE, kFALSE, kFALSE};
	if(dimu->GetMu(0)->GetMatchTrigger()>0 || dimu->GetMu(1)->GetMatchTrigger()>0) trigOk[1] = kTRUE;
	if(dimu->GetMu(0)->GetMatchTrigger()>0 && dimu->GetMu(1)->GetMatchTrigger()>0) trigOk[2] = kTRUE;
	for (Int_t itrg = 0; itrg < 3; itrg++) {
	  
	  if (!trigOk[itrg]) continue;
	  
	  // pt/y integrated
	  intrec[itrg][0][0]++;
	  bin[0] = 1000+itrg; bin[1] = 0.; bin[2] = 0.;
	  hrec->Fill(bin);
	  if (itrg == 2) {
	    bin[0] = irun;
	    hrec->Fill(bin);
	    hPtRec->Fill(dimu->Pt());
	    hYRec->Fill(dimu->Y());
	  }
	  
	  // y integrated / pt bins
	  for (Int_t ipt = 0; ipt < nPtBins; ipt++) {
	    
	    if (dimu->Pt() >= pTBinLowEdge[ipt] && dimu->Pt() < pTBinLowEdge[ipt+1]) {
	      
	      intrec[itrg][ipt+1][0]++;
	      bin[0] = 1000+itrg; bin[1] = ipt+1; bin[2] = 0.;
	      hrec->Fill(bin);
	      if (itrg == 2) {
		bin[0] = irun;
		hrec->Fill(bin);
	      }
	      
	    }
	    
	  }
	  
	  // y bins
	  for (Int_t iy = 0; iy < nYBins; iy++) {
	    
	    if (dimu->Y() > yBinLowEdge[iy] && dimu->Y() < yBinLowEdge[iy+1]) {
	      
	      // pt integrated
	      intrec[itrg][0][iy+1]++;
	      bin[0] = 1000+itrg; bin[1] = 0.; bin[2] = iy+1;
	      hrec->Fill(bin);
	      if (itrg == 2) {
		bin[0] = irun;
		hrec->Fill(bin);
	      }
	      
	      // pt bins
	      for (Int_t ipt = 0; ipt < nPtBins; ipt++) {
		
		if (dimu->Pt() >= pTBinLowEdge[ipt] && dimu->Pt() < pTBinLowEdge[ipt+1]) {
		  
		  intrec[itrg][ipt+1][iy+1]++;
		  bin[0] = 1000+itrg; bin[1] = ipt+1; bin[2] = iy+1;
		  hrec->Fill(bin);
		  if (itrg == 2) {
		    bin[0] = irun;
		    hrec->Fill(bin);
		  }
		  
		}
		
	      }
	      
	    }
	    
	  }
	  
	}
	
      }
      
    }
    
  }
  
  //   //----NUMBERS---- for -2.5 -4
  Double_t eff0 = (Double_t)intrec[0][0][0]/(Double_t)intgen[0][0];
  Double_t eff1 = (Double_t)intrec[1][0][0]/(Double_t)intgen[0][0];
  Double_t eff2 = (Double_t)intrec[2][0][0]/(Double_t)intgen[0][0];
  Double_t erreff0 = TMath::Sqrt((Double_t)intrec[0][0][0])/(Double_t)intgen[0][0];
  Double_t erreff1 = TMath::Sqrt((Double_t)intrec[1][0][0])/(Double_t)intgen[0][0];
  Double_t erreff2 = TMath::Sqrt((Double_t)intrec[2][0][0])/(Double_t)intgen[0][0];
  
  //----WRITING RESULTS----
  if (print) {
    printf("---- VALUES\n");
    printf("---- generated: %d\n",intgen[0][0]);
    printf("---- nomatching: %d\n",intrec[0][0][0]);
    printf("---- 1matching: %d\n",intrec[1][0][0]);
    printf("---- 2matching: %d\n",intrec[2][0][0]);
    printf("\n---- EFFICIENCIES\n");
    printf("---- nomatching required: ( %f +- %f )\n",eff0,erreff0);
    printf("---- 1matching  required: ( %f +- %f )\n",eff1,erreff1);
    printf("---- 2matching  required: ( %f +- %f )\n\n",eff2,erreff2);
  }
  
  return nEvents;
}

//---------------------------------------------------------------------------
void ComputeAccEff(TList& runs, Int_t ptBin, Int_t yBin, TH1D* hgen, TH1D* hrec)
{
  // compute integrated acc*eff correction
  
  TGraphAsymmErrors *acceff = new TGraphAsymmErrors(hrec, hgen);
  
  // set bin labels
  acceff->GetXaxis()->Set(1003, -0.5, 1002.5);
  TIter nextRun(&runs);
  TObjString *srun = 0x0;
  Int_t irun = 1;
  while ((srun = static_cast<TObjString*>(nextRun())))
    acceff->GetXaxis()->SetBinLabel(irun++, srun->GetName());
  acceff->GetXaxis()->SetBinLabel(1001, "All0Match");
  acceff->GetXaxis()->SetBinLabel(1002, "All1Match");
  acceff->GetXaxis()->SetBinLabel(1003, "All2Match");
  acceff->GetXaxis()->SetRange(0, irun-1);
  
  // display
  TString title = "acceptance * efficiency versus run";
  title += (ptBin == 0) ? Form(" (%g<pt<%g",pTBinLowEdge[0],pTBinLowEdge[nPtBins]) : Form(" (%g<pt<%g",pTBinLowEdge[ptBin-1],pTBinLowEdge[ptBin]);
  title += (yBin == 0) ? Form(" / %g<y<%g)",yBinLowEdge[0],yBinLowEdge[nYBins]) : Form(" / %g<y<%g)",yBinLowEdge[yBin-1],yBinLowEdge[yBin]);
  new TCanvas(Form("cAccEffVsRun_pt%d_y%d",ptBin,yBin), title.Data(), 1000, 300);
  acceff->SetNameTitle(Form("accEffVsRun_pt%d_y%d",ptBin,yBin), title.Data());
  acceff->SetLineStyle(1);
  acceff->SetLineColor(1); 
  acceff->SetMarkerStyle(20);
  acceff->SetMarkerSize(0.7);
  acceff->SetMarkerColor(2);
  acceff->GetXaxis()->SetTitle("Run #");
  acceff->GetXaxis()->SetLabelFont(22);
  acceff->GetXaxis()->SetTitleFont(22);
  acceff->GetYaxis()->SetTitle("Acc*Eff");
  acceff->GetYaxis()->SetLabelFont(22);
  acceff->GetYaxis()->SetLabelFont(22);
  acceff->Draw("ap");
  
  // save output
  TFile* file = new TFile("acceff.root","update");
  hgen->Write(0x0, TObject::kOverwrite);
  hrec->Write(0x0, TObject::kOverwrite);
  acceff->Write(0x0, TObject::kOverwrite);
  file->Close();
  
  // print integrated value
  title = "---- Integrated acc*eff";
  title += (ptBin == 0) ? Form(" (%g<pt<%g",pTBinLowEdge[0],pTBinLowEdge[nPtBins]) : Form(" (%g<pt<%g",pTBinLowEdge[ptBin-1],pTBinLowEdge[ptBin]);
  title += (yBin == 0) ? Form(" / %g<y<%g):",yBinLowEdge[0],yBinLowEdge[nYBins]) : Form(" / %g<y<%g):",yBinLowEdge[yBin-1],yBinLowEdge[yBin]);
  printf("\n%s\n",title.Data());
  Int_t lastPoint = acceff->GetN()-1;
  Double_t x,eff;
  acceff->GetPoint(lastPoint-2,x,eff);
  printf("---- no matching required: %f + %f - %f\n", eff, acceff->GetErrorYhigh(lastPoint-2), acceff->GetErrorYlow(lastPoint-2));
  acceff->GetPoint(lastPoint-1,x,eff);
  printf("---- 1 matching  required: %f + %f - %f\n", eff, acceff->GetErrorYhigh(lastPoint-1), acceff->GetErrorYlow(lastPoint-1));
  acceff->GetPoint(lastPoint,x,eff);
  printf("---- 2 matching  required: %f + %f - %f\n", eff, acceff->GetErrorYhigh(lastPoint), acceff->GetErrorYlow(lastPoint));
  
  // fill summary histos
  TString label = "";
  label += (ptBin == 0) ? Form("%g<pt<%g",pTBinLowEdge[0],pTBinLowEdge[nPtBins]) : Form("%g<pt<%g",pTBinLowEdge[ptBin-1],pTBinLowEdge[ptBin]);
  label += (yBin == 0) ? Form(" / %g<y<%g",yBinLowEdge[0],yBinLowEdge[nYBins]) : Form(" / %g<y<%g",yBinLowEdge[yBin-1],yBinLowEdge[yBin]);
  hGenSummary->Fill(ptBin*(nYBins+1)+yBin, hgen->GetBinContent(1003));
  hGenSummary->GetXaxis()->SetBinLabel(ptBin*(nYBins+1)+yBin+1, label.Data());
  hRecSummary->Fill(ptBin*(nYBins+1)+yBin, hrec->GetBinContent(1003));
  hRecSummary->GetXaxis()->SetBinLabel(ptBin*(nYBins+1)+yBin+1, label.Data());
  hAccSummary->SetBinContent(ptBin*(nYBins+1)+yBin+1, eff);
  hAccSummary->SetBinError(ptBin*(nYBins+1)+yBin+1, acceff->GetErrorYhigh(lastPoint));
  hAccSummary->GetXaxis()->SetBinLabel(ptBin*(nYBins+1)+yBin+1, label.Data());
}

