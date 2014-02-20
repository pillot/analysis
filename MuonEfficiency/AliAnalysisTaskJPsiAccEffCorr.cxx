/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include <Riostream.h>
#include <TRandom3.h>

// ROOT includes
#include <TMath.h>
#include <TH1F.h>
#include <TList.h>
#include <THashList.h>
#include <TObjArray.h>
#include <TClonesArray.h>
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>
#include <TFile.h>

// STEER includes
#include "AliAODEvent.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliAODTrack.h"
#include "AliAODDimuon.h"

// ANALYSIS includes
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisDataSlot.h"
#include "AliAnalysisDataContainer.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskJPsiAccEffCorr.h"

// PWG3 includes
#include "AliCounterCollection.h"


ClassImp(AliAnalysisTaskJPsiAccEffCorr)

//________________________________________________________________________
AliAnalysisTaskJPsiAccEffCorr::AliAnalysisTaskJPsiAccEffCorr() :
AliAnalysisTaskSE(), 
fList(0x0),
fEventCounters(0x0),
fJPsiCounters(0x0),
fCentBinLowEdge(0),
fPtBinLowEdge(0),
fYBinLowEdge(0),
fTrigLevel(1),
fMuLowPtCut(-1.)
{
  // Dummy constructor
}

//________________________________________________________________________
AliAnalysisTaskJPsiAccEffCorr::AliAnalysisTaskJPsiAccEffCorr(const char *name) :
AliAnalysisTaskSE(name), 
fList(0x0),
fEventCounters(0x0),
fJPsiCounters(0x0),
fCentBinLowEdge(0),
fPtBinLowEdge(0),
fYBinLowEdge(0),
fTrigLevel(1),
fMuLowPtCut(-1.)
{
  /// Constructor
  
  // set default binning
  Float_t centBinLowEdge[21] = {0., 5., 10., 15., 20., 25., 30., 35., 40., 45., 50., 55., 60., 65., 70., 75., 80., 85., 90., 95., 100.};
  //Float_t centBinLowEdge[11] = {0., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100.};
  SetCentBins((Int_t)(sizeof(centBinLowEdge)/sizeof(Float_t))-1, centBinLowEdge);
  Float_t pTBinLowEdge[3] = {0., 3., 1000000.};
  SetPtBins((Int_t)(sizeof(pTBinLowEdge)/sizeof(Float_t))-1, pTBinLowEdge);
  Float_t yBinLowEdge[3] = {-4., -3.25, -2.5};
  SetYBins((Int_t)(sizeof(yBinLowEdge)/sizeof(Float_t))-1, yBinLowEdge);
  
  // Output slot #1 writes into a TObjArray container
  DefineOutput(1,TObjArray::Class());
  // Output slot #2 writes event counters
  DefineOutput(2,AliCounterCollection::Class());
  // Output slot #3 writes JPsi counters
  DefineOutput(3,AliCounterCollection::Class());
}

//________________________________________________________________________
AliAnalysisTaskJPsiAccEffCorr::~AliAnalysisTaskJPsiAccEffCorr()
{
  /// Destructor
  if (!AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
    delete fList;
    delete fEventCounters;
    delete fJPsiCounters;
  }
}

//___________________________________________________________________________
void AliAnalysisTaskJPsiAccEffCorr::UserCreateOutputObjects()
{
  /// Create histograms and counters
  
  //gRandom->SetSeed(0);
  
  // initialize histos
  fList = new TObjArray(2000);
  fList->SetOwner();
  
  TH1F* hPtGen = new TH1F("hPtGen","generated J/#psi versus pt;p_{T} (GeV/c);N_{J/#psi} / 0.5 GeV/c", 300, 0., 30.);
  fList->AddAtAndExpand(hPtGen, kPtGen);
  TH1F* hPtRec = new TH1F("hPtRec","reconstructed J/#psi versus pt;p_{T} (GeV/c);N_{J/#psi} / 0.5 GeV/c", 300, 0., 30.);
  fList->AddAtAndExpand(hPtRec, kPtRec);
  
  TH1F* hYGen = new TH1F("hYGen","generated J/#psi versus y;y;N_{J/#psi}", 15, -4., -2.5);
  fList->AddAtAndExpand(hYGen, kYGen);
  TH1F* hYRec = new TH1F("hYRec","reconstructed J/#psi versus y;y;N_{J/#psi}", 15, -4., -2.5);
  fList->AddAtAndExpand(hYRec, kYRec);
  
  TH1F* hPtGenMu = new TH1F("hPtGenMu","generated single-#mu versus pt;p_{T} (GeV/c);N_{single-#mu} / 0.5 GeV/c", 300, 0., 30.);
  fList->AddAtAndExpand(hPtGenMu, kPtGenMu);
  TH1F* hPtRecMu = new TH1F("hPtRecMu","reconstructed single-#mu versus pt;p_{T} (GeV/c);N_{single-#mu} / 0.5 GeV/c", 300, 0., 30.);
  fList->AddAtAndExpand(hPtRecMu, kPtRecMu);
  
  TH1F* hYGenMu = new TH1F("hYGenMu","generated single-#mu versus y;y;N_{single-#mu}", 15, -4., -2.5);
  fList->AddAtAndExpand(hYGenMu, kYGenMu);
  TH1F* hYRecMu = new TH1F("hYRecMu","reconstructed single-#mu versus y;y;N_{single-#mu}", 15, -4., -2.5);
  fList->AddAtAndExpand(hYRecMu, kYRecMu);
  
  TH1F* hDzVtx = new TH1F("hDzVtx","vertex resolution;#DeltaZ (cm)", 200, -1., 1.);
  fList->AddAtAndExpand(hDzVtx, kDzVtx);
  
  TH1F* hDzVtxJPsi = new TH1F("hDzVtxJPsi","vertex resolution for event with J/#psi;#DeltaZ (cm)", 200, -1., 1.);
  fList->AddAtAndExpand(hDzVtxJPsi, kDzVtx2);
  
  TH1F* hMass = new TH1F("hMass","invariant mass distribution;Mass (GeV/c^{2};N_{J/#psi} / 0.025 GeV/c)", 200., 0., 5.);
  fList->AddAtAndExpand(hMass, kMass);
  
  // initialize counters
  
  // centrality binning
  TString centbins = "any";
  for (Int_t i=0; i<fCentBinLowEdge.GetSize()-1; i++) centbins += Form("/%g-%g",fCentBinLowEdge[i],fCentBinLowEdge[i+1]);
  
  // pt binning
  TString ptbins = Form("%g-%g",fPtBinLowEdge[0],fPtBinLowEdge[fPtBinLowEdge.GetSize()-1]);
  for (Int_t i=0; i<fPtBinLowEdge.GetSize()-1; i++) ptbins += Form("/%g-%g",fPtBinLowEdge[i],fPtBinLowEdge[i+1]);
  
  // y binning
  TString ybins = Form("%g-%g",fYBinLowEdge[0],fYBinLowEdge[fYBinLowEdge.GetSize()-1]);
  for (Int_t i=0; i<fYBinLowEdge.GetSize()-1; i++) ybins += Form("/%g-%g",fYBinLowEdge[i],fYBinLowEdge[i+1]);
  
  fEventCounters = new AliCounterCollection(GetOutputSlot(2)->GetContainer()->GetName());
  fEventCounters->AddRubric("event", "any");
  fEventCounters->AddRubric("cent", centbins.Data());
  fEventCounters->AddRubric("run", 100000);
  fEventCounters->Init();
  
  fJPsiCounters = new AliCounterCollection(GetOutputSlot(3)->GetContainer()->GetName());
  fJPsiCounters->AddRubric("type", "gen/rec");
  fJPsiCounters->AddRubric("match", "any/1/2");
  fJPsiCounters->AddRubric("ptbin", ptbins.Data());
  fJPsiCounters->AddRubric("ybin", ybins.Data());
  fJPsiCounters->AddRubric("cent", centbins.Data());
  fJPsiCounters->AddRubric("run", 100000);
  fJPsiCounters->Init();
  
  // Post ata at least once per task to ensure data synchronisation (required for merging)
  PostData(1, fList);
  PostData(2, fEventCounters);
  PostData(3, fJPsiCounters);
}

//________________________________________________________________________
void AliAnalysisTaskJPsiAccEffCorr::UserExec(Option_t *)
{
  /// Called for each event
  
  // get AOD event
  AliAODEvent* aod = dynamic_cast<AliAODEvent*>(InputEvent());
  if ( !aod ) return;
  
  // get the centrality percentile
  Float_t centrality = aod->GetCentrality()->GetCentralityPercentileUnchecked("V0M");
  TString centKey = "";
  for (Int_t icent = 0; icent < fCentBinLowEdge.GetSize()-1; icent++)
    if (centrality > fCentBinLowEdge[icent] && centrality <= fCentBinLowEdge[icent+1])
      centKey = Form("%g-%g",fCentBinLowEdge[icent],fCentBinLowEdge[icent+1]);
  
  //if (centrality > 90) return;
  
  TString matchKey[3] = {"any", "1", "2"};
//  Double_t rAbsMin = 17.622;
//  Double_t rAbsMax = 89.;
  Double_t rAbsMin = 17.6;
  Double_t rAbsMax = 89.5;
  
  TClonesArray *mcarray = dynamic_cast<TClonesArray*>(aod->FindListObject(AliAODMCParticle::StdBranchName()));
  /*
  // remove events with JPsi outside -4.2<y<-2.3
  for (Int_t ii=0; ii<mcarray->GetEntries(); ii++) {
    AliAODMCParticle *mctrack = (AliAODMCParticle*) mcarray->UncheckedAt(ii); 
    if (mctrack->GetPdgCode() == 443 && (mctrack->Y()<-4.1 || mctrack->Y()>-2.4)) return;
  }
  */
  // fill event counter
  fEventCounters->Count(Form("event:any/cent:any/run:%d",fCurrentRunNumber));
  if (!centKey.IsNull()) fEventCounters->Count(Form("event:any/cent:%s/run:%d",centKey.Data(),fCurrentRunNumber));
  
  // vertex resolution
  AliAODMCHeader* aodMCHeader = static_cast<AliAODMCHeader*>(aod->FindListObject(AliAODMCHeader::StdBranchName()));
  ((TH1F*)fList->UncheckedAt(kDzVtx))->Fill(aod->GetPrimaryVertex()->GetZ()-aodMCHeader->GetVtxZ());
  
  // loop over MC tracks
  for (Int_t ii=0; ii<mcarray->GetEntries(); ii++) {
    
    AliAODMCParticle *mctrack = (AliAODMCParticle*) mcarray->UncheckedAt(ii); 
    
    // look for muons
    if (TMath::Abs(mctrack->GetPdgCode()) == 13 &&
	mctrack->Eta()>fYBinLowEdge[0] && mctrack->Eta()<fYBinLowEdge[fYBinLowEdge.GetSize()-1]) {
      
      ((TH1F*)fList->UncheckedAt(kPtGenMu))->Fill(mctrack->Pt());
      ((TH1F*)fList->UncheckedAt(kYGenMu))->Fill(mctrack->Y());
      
    }
    
    // look for generated particles
    if(mctrack->IsPrimary() && !mctrack->IsPhysicalPrimary() &&
       mctrack->Y()>fYBinLowEdge[0] && mctrack->Y()<fYBinLowEdge[fYBinLowEdge.GetSize()-1]) {
      
      ((TH1F*)fList->UncheckedAt(kPtGen))->Fill(mctrack->Pt());
      ((TH1F*)fList->UncheckedAt(kYGen))->Fill(mctrack->Y());
      
      // pt bin
      TString ptKey = "";
      for (Int_t ipt = 0; ipt < fPtBinLowEdge.GetSize()-1; ipt++)
	if (mctrack->Pt() >= fPtBinLowEdge[ipt] && mctrack->Pt() < fPtBinLowEdge[ipt+1])
	  ptKey = Form("%g-%g",fPtBinLowEdge[ipt],fPtBinLowEdge[ipt+1]);
      if (ptKey.IsNull()) continue;
      
      // y bin
      TString yKey = "";
      for (Int_t iy = 0; iy < fYBinLowEdge.GetSize()-1; iy++)
	if (mctrack->Y() >= fYBinLowEdge[iy] && mctrack->Y() < fYBinLowEdge[iy+1])
	  yKey = Form("%g-%g",fYBinLowEdge[iy],fYBinLowEdge[iy+1]);
      if (yKey.IsNull()) continue;
      
      for (Int_t itrg = 0; itrg < 3; itrg++) {
	fJPsiCounters->Count(Form("type:gen/match:%s/ptbin:%s/ybin:%s/cent:any/run:%d",matchKey[itrg].Data(),ptKey.Data(),yKey.Data(),fCurrentRunNumber));
	if (!centKey.IsNull()) fJPsiCounters->Count(Form("type:gen/match:%s/ptbin:%s/ybin:%s/cent:%s/run:%d",matchKey[itrg].Data(),ptKey.Data(),yKey.Data(),centKey.Data(),fCurrentRunNumber));
      }
      
    }
    
  }
  
  // loop over reco tracks
  Int_t ntracks = aod->GetNTracks();
  for (Int_t q=0; q<ntracks; q++){
    
    AliAODTrack *mu = aod->GetTrack(q);
    
    //Double_t ptCut = gRandom->Gaus(fMuLowPtCut,0.3);
    
    if (mu->IsMuonTrack() && mu->GetMatchTrigger()>=fTrigLevel && mu->GetLabel() >= 0 &&
	mu->Eta()>fYBinLowEdge[0] && mu->Eta()<fYBinLowEdge[fYBinLowEdge.GetSize()-1] &&
	//mu->Pt()>ptCut &&
	mu->GetRAtAbsorberEnd()>rAbsMin && mu->GetRAtAbsorberEnd()<rAbsMax
	) {
      
      ((TH1F*)fList->UncheckedAt(kPtRecMu))->Fill(mu->Pt());
      ((TH1F*)fList->UncheckedAt(kYRecMu))->Fill(mu->Y());
      
    }
    
  }
  
  // loop over reco dimuons
  Bool_t jpsiFound = kFALSE;
  Int_t ndimu = aod->GetNDimuons();
  for(Int_t q=0; q<ndimu; q++) {
    
    AliAODDimuon *dimu = aod->GetDimuon(q);
    
    //Double_t ptCut1 = gRandom->Gaus(fMuLowPtCut,0.3);
    //Double_t ptCut2 = gRandom->Gaus(fMuLowPtCut,0.3);
    
    //Double_t thetaTrackAbsEnd0 = TMath::ATan(dimu->GetMu(0)->GetRAtAbsorberEnd()/505.) * TMath::RadToDeg();
    //Double_t thetaTrackAbsEnd1 = TMath::ATan(dimu->GetMu(1)->GetRAtAbsorberEnd()/505.) * TMath::RadToDeg();
    if(dimu->Charge()==0 && dimu->Y()>fYBinLowEdge[0] && dimu->Y()<fYBinLowEdge[fYBinLowEdge.GetSize()-1] &&
       dimu->GetMu(0)->Eta()>-4. && dimu->GetMu(0)->Eta()<-2.5 &&
       dimu->GetMu(1)->Eta()>-4. && dimu->GetMu(1)->Eta()<-2.5 &&
       //dimu->GetMu(0)->Pt()>ptCut1 && dimu->GetMu(1)->Pt()>ptCut2 &&
       dimu->GetMu(0)->Pt()>fMuLowPtCut && dimu->GetMu(1)->Pt()>fMuLowPtCut &&
       dimu->GetMu(0)->GetRAtAbsorberEnd()>rAbsMin && dimu->GetMu(0)->GetRAtAbsorberEnd()<rAbsMax &&
       dimu->GetMu(1)->GetRAtAbsorberEnd()>rAbsMin && dimu->GetMu(1)->GetRAtAbsorberEnd()<rAbsMax 
       //thetaTrackAbsEnd0>2. && thetaTrackAbsEnd0<10. &&
       //thetaTrackAbsEnd1>2. && thetaTrackAbsEnd1<10.
       && dimu->GetMu(0)->GetLabel() >= 0 && dimu->GetMu(1)->GetLabel() >= 0
       ){
      
      //if (dimu->GetMu(0)->GetLabel() < 0 || dimu->GetMu(1)->GetLabel() < 0) printf("no label\n");
      /*
      // remove events with muons generated outside 168.5<Theta<178.5
      Bool_t accOk = kTRUE;
      for (Int_t imu = 0; imu < 2; imu++) {
	Int_t label = dimu->GetMu(imu)->GetLabel();
	AliAODMCParticle *mctrack = (label >= 0) ? (AliAODMCParticle*) mcarray->UncheckedAt(label) : 0x0;
	if (mctrack && mctrack->GetPdgCode() == 13 && mctrack->GetMother() >= 0 &&
	    ((AliAODMCParticle*) mcarray->UncheckedAt(mctrack->GetMother()))->GetPdgCode() == 443 &&
	    (mctrack->Theta()*TMath::RadToDeg() < 168.5 || mctrack->Theta()*TMath::RadToDeg() > 178.5)) {
	  printf("theta = %f\n",mctrack->Theta()*TMath::RadToDeg());
	  accOk = kFALSE;
	}
      }
      if (!accOk) {
	printf("--> reject event\n");
	continue;
      }
      */
      Bool_t trigOk[3] = {kTRUE, kFALSE, kFALSE};
      if(dimu->GetMu(0)->GetMatchTrigger()>=fTrigLevel || dimu->GetMu(1)->GetMatchTrigger()>=fTrigLevel) trigOk[1] = kTRUE;
      if(dimu->GetMu(0)->GetMatchTrigger()>=fTrigLevel && dimu->GetMu(1)->GetMatchTrigger()>=fTrigLevel) trigOk[2] = kTRUE;
      
      if (trigOk[2]) {
	((TH1F*)fList->UncheckedAt(kPtRec))->Fill(dimu->Pt());
	((TH1F*)fList->UncheckedAt(kYRec))->Fill(dimu->Y());
	((TH1F*)fList->UncheckedAt(kMass))->Fill(dimu->M());
	//((TH1F*)fList->UncheckedAt(kPtRecMu))->Fill(dimu->GetMu(0)->Pt());
	//((TH1F*)fList->UncheckedAt(kYRecMu))->Fill(dimu->GetMu(0)->Y());
	//((TH1F*)fList->UncheckedAt(kPtRecMu))->Fill(dimu->GetMu(1)->Pt());
	//((TH1F*)fList->UncheckedAt(kYRecMu))->Fill(dimu->GetMu(1)->Y());
	jpsiFound = kTRUE;
      }
      
      // pt bin
      TString ptKey;
      for (Int_t ipt = 0; ipt < fPtBinLowEdge.GetSize()-1; ipt++)
	if (dimu->Pt() >= fPtBinLowEdge[ipt] && dimu->Pt() < fPtBinLowEdge[ipt+1])
	  ptKey = Form("%g-%g",fPtBinLowEdge[ipt],fPtBinLowEdge[ipt+1]);
      if (ptKey.IsNull()) continue;
      
      // y bin
      TString yKey;
      for (Int_t iy = 0; iy < fYBinLowEdge.GetSize()-1; iy++)
	if (dimu->Y() >= fYBinLowEdge[iy] && dimu->Y() < fYBinLowEdge[iy+1])
	  yKey = Form("%g-%g",fYBinLowEdge[iy],fYBinLowEdge[iy+1]);
      if (yKey.IsNull()) continue;
      
      for (Int_t itrg = 0; itrg < 3; itrg++) {
	if (!trigOk[itrg]) continue;
	fJPsiCounters->Count(Form("type:rec/match:%s/ptbin:%s/ybin:%s/cent:any/run:%d",matchKey[itrg].Data(),ptKey.Data(),yKey.Data(),fCurrentRunNumber));
	if (!centKey.IsNull()) fJPsiCounters->Count(Form("type:rec/match:%s/ptbin:%s/ybin:%s/cent:%s/run:%d",matchKey[itrg].Data(),ptKey.Data(),yKey.Data(),centKey.Data(),fCurrentRunNumber));
      }
      
    }
    
  }
  
  // vertex resolution for event with a valid reconstructed JPsi
  if (jpsiFound) ((TH1F*)fList->UncheckedAt(kDzVtx2))->Fill(aod->GetPrimaryVertex()->GetZ()-aodMCHeader->GetVtxZ());
  
  // Post final data. It will be written to a file with option "RECREATE"
  PostData(1, fList);
  PostData(2, fEventCounters);
  PostData(3, fJPsiCounters);
  
}

//________________________________________________________________________
void AliAnalysisTaskJPsiAccEffCorr::Terminate(Option_t *)
{
  /// draw final results
  
  fList = static_cast<TObjArray*> (GetOutputData(1));
  fEventCounters = static_cast<AliCounterCollection*> (GetOutputData(2));
  fEventCounters->Sort("run", kTRUE);
  Int_t nEv = (Int_t) fEventCounters->GetSum();
  if (nEv == 0) {
    AliError("event counter is empty");
    return;
  }
  fJPsiCounters = static_cast<AliCounterCollection*> (GetOutputData(3));
  fJPsiCounters->Sort("run", kTRUE);
  if (((Int_t)fJPsiCounters->GetSum()) == 0) {
    AliError("J/Psi counter is empty");
    return;
  }
  
  // number of events
  printf("\nTotal number of events = %d\n", nEv);
  TH1D *hevVsRun = fEventCounters->Get("run","");
  TH1D *hevVsCent = fEventCounters->Get("cent","");
  TFile* file = new TFile("acceff_new.root","recreate");
  hevVsRun->Write(0x0, TObject::kOverwrite);
  hevVsCent->Write(0x0, TObject::kOverwrite);
  file->Close();
  delete hevVsRun;
  delete hevVsCent;
  
  // compute acc*eff vs. pt (J/Psi)
  TH1F* hPtGen = (TH1F*)fList->UncheckedAt(kPtGen);
  TH1F* hPtRec = (TH1F*)fList->UncheckedAt(kPtRec);
  TH1F* hPtAcc = ComputeAccEff(*hPtGen, *hPtRec, "hPtAcc","Acc * Eff versus pt (J/#psi)");
  
  // compute acc*eff vs. y (J/Psi)
  TH1F* hYGen = (TH1F*)fList->UncheckedAt(kYGen);
  TH1F* hYRec = (TH1F*)fList->UncheckedAt(kYRec);
  TH1F* hYAcc = ComputeAccEff(*hYGen, *hYRec, "hYAcc","Acc * Eff versus y (J/#psi)");
  
  // compute acc*eff vs. pt (single mu)
  TH1F* hPtGenMu = (TH1F*)fList->UncheckedAt(kPtGenMu);
  TH1F* hPtRecMu = (TH1F*)fList->UncheckedAt(kPtRecMu);
  TH1F* hPtAccMu = ComputeAccEff(*hPtGenMu, *hPtRecMu, "hPtAccMu","Acc * Eff versus pt (single-#mu)");
  
  // compute acc*eff vs. y (single mu)
  TH1F* hYGenMu = (TH1F*)fList->UncheckedAt(kYGenMu);
  TH1F* hYRecMu = (TH1F*)fList->UncheckedAt(kYRecMu);
  TH1F* hYAccMu = ComputeAccEff(*hYGenMu, *hYRecMu, "hYAccMu","Acc * Eff versus y (single-#mu)");
  
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
  
  // save histos
  file = new TFile("acceff_new.root","update");
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
  file->Close();
  
  // number of track matching required
  Int_t nMatch = 2;
  
  // Ncoll per centrality bin (width=10) to weight the acc*eff calculation
  Float_t nColl10[8] = {1502.7, 923.26, 558.68, 321.20, 171.67, 85.13, 38.51, 15.78};
  Float_t nColl10binLowEdge[10] = {0., 10., 20., 30., 40., 50., 60., 70., 80., 90.};
  
  // Ncoll per centrality bin (width=5) to weight the acc*eff calculation
  Float_t nColl5[16] = {1686.87, 1319.89, 1031.9, 807.90, 627.99, 483.95, 369.13, 274.03, 199.30, 143.45, 100.54, 68.82, 46.09, 29.70, 18.80, 11.95};
  Float_t nColl5binLowEdge[17] = {0., 5., 10., 15., 20., 25., 30., 35., 40., 45., 50., 55., 60., 65., 70., 75., 80.};
  
  // number of JPsi per centrality/pt/y bin to weight the acc*eff calculation
  Float_t nJPsi[3][3][5] = {{{973., 573., 406., 306., 119.}, {377., 257., 142., 126., 45.}, {586., 327., 260., 182., 79.}},
			    {{785., 450., 323., 215., 81.}, {973., 573., 406., 306., 119.}, {973., 573., 406., 306., 119.}},
			    {{0., 0., 340., 151., 56.}, {973., 573., 406., 306., 119.}, {973., 573., 406., 306., 119.}}};
  Float_t jPsibinLowEdge[3][3][6] = {{{0., 10., 20., 30., 50., 80.}, {0., 10., 20., 30., 50., 80.}, {0., 10., 20., 30., 50., 80.}},
				    {{0., 10., 20., 30., 50., 80.}, {0., 10., 20., 30., 50., 80.}, {0., 10., 20., 30., 50., 80.}},
				    {{-1., -1., 0., 20., 40., 80.}, {0., 10., 20., 30., 50., 80.}, {0., 10., 20., 30., 50., 80.}}};
  
  // acceptance*efficiency for different pt bins
  TH1F* hGenPtSummary = new TH1F("hGenPtSummary","generated J/#psi versus pt bins;;a.u.", fPtBinLowEdge.GetSize(), -0.5, fPtBinLowEdge.GetSize()-0.5);
  TH1F* hRecPtSummary = new TH1F("hRecPtSummary","reconstructed J/#psi versus pt bins;;a.u.", fPtBinLowEdge.GetSize(), -0.5, fPtBinLowEdge.GetSize()-0.5);
  TH1F* hAccPtSummary = new TH1F("hAccPtSummary","Acc * Eff versus pt bins", fPtBinLowEdge.GetSize(), -0.5, fPtBinLowEdge.GetSize()-0.5);
  
  // loop over pt bins+1 (0 = integrated over pt)
  for (Int_t ipt = 0; ipt < fPtBinLowEdge.GetSize(); ipt++) {
    
    // acceptance*efficiency versus run
    DrawAccEffVsRun(ipt, 0);
    
    // acceptance*efficiency versus centrality
    DrawAccEffVsCent(ipt, 0);
    
    // print integrated value
    Double_t ae, aee;
    TString title = "---- Integrated acc*eff";
    title += (ipt == 0) ? Form(" (%g<pt<%g",fPtBinLowEdge[0],fPtBinLowEdge[fPtBinLowEdge.GetSize()-1]) : Form(" (%g<pt<%g",fPtBinLowEdge[ipt-1],fPtBinLowEdge[ipt]);
    title += Form(" / %g<y<%g):",fYBinLowEdge[0],fYBinLowEdge[fYBinLowEdge.GetSize()-1]);
    printf("\n%s\n",title.Data());
    IntegratedAccEff(ipt, 0, nMatch, 0., -1., ae, aee, kTRUE, hGenPtSummary, hRecPtSummary, hAccPtSummary, ipt);
    
    // print integrated value in 0-80%
    IntegratedAccEff(ipt, 0, nMatch, 0., 80., ae, aee);
    printf("   -  %d matching required - 0-80%%: %f ± %f\n", nMatch, ae, aee);
    
    // print integrated value in 0-80% weighted by nColl (bin=10)
    IntegratedAccEff(ipt, 0, nMatch, 8, nColl10, nColl10binLowEdge, ae, aee);
    printf("   -  %d matching required - 0-80%% - weighted per nColl (bin=10): %f ± %f\n", nMatch, ae, aee);
    
    // print integrated value in 0-80% weighted by nColl (bin=5)
    IntegratedAccEff(ipt, 0, nMatch, 16, nColl5, nColl5binLowEdge, ae, aee);
    printf("   -  %d matching required - 0-80%% - weighted per nColl (bin=5): %f ± %f\n", nMatch, ae, aee);
    
    // print integrated value in 0-80% weighted by nJPsi integrated over pt/y
    IntegratedAccEff(ipt, 0, nMatch, 5, nJPsi[0][0], jPsibinLowEdge[0][0], ae, aee);
    printf("   -  %d matching required - 0-80%% - weighted per nJPsi integrated over pt/y: %f ± %f\n", nMatch, ae, aee);
    
  }
  
  // draw summary histos
  TCanvas* cSummary = new TCanvas("cPtSummary", "pTSummary", 900, 300);
  cSummary->Divide(3,1);
  cSummary->cd(1);
  hGenPtSummary->SetStats(kFALSE);
  hGenPtSummary->Draw("e0");
  cSummary->cd(2);
  hRecPtSummary->SetStats(kFALSE);
  hRecPtSummary->Draw("e0");
  cSummary->cd(3);
  hAccPtSummary->SetMarkerStyle(21);
  hAccPtSummary->SetMarkerSize(0.6);
  hAccPtSummary->SetStats(kFALSE);
  hAccPtSummary->Draw("e0");
  
  // save histos
  file = new TFile("acceff_new.root","update");
  hGenPtSummary->Write(0x0, TObject::kOverwrite);
  hRecPtSummary->Write(0x0, TObject::kOverwrite);
  hAccPtSummary->Write(0x0, TObject::kOverwrite);
  cSummary->Write(0x0, TObject::kOverwrite);
  file->Close();

  // acceptance*efficiency for different y bins
  TH1F* hGenYSummary = new TH1F("hGenYSummary","generated J/#psi versus y bins;;a.u.", fYBinLowEdge.GetSize(), -0.5, fYBinLowEdge.GetSize()-0.5);
  TH1F* hRecYSummary = new TH1F("hRecYSummary","reconstructed J/#psi versus y bins;;a.u.", fYBinLowEdge.GetSize(), -0.5, fYBinLowEdge.GetSize()-0.5);
  TH1F* hAccYSummary = new TH1F("hAccYSummary","Acc * Eff versus y bins", fYBinLowEdge.GetSize(), -0.5, fYBinLowEdge.GetSize()-0.5);
  
  // loop over y bins+1 (0 = integrated over y)
  for (Int_t iy = 0; iy < fYBinLowEdge.GetSize(); iy++) {
    
    // acceptance*efficiency versus run
    if (iy > 0) DrawAccEffVsRun(0, iy);
    
    // acceptance*efficiency versus centrality
    if (iy > 0) DrawAccEffVsCent(0, iy);
    
    // print integrated value
    Double_t ae, aee;
    TString title = "---- Integrated acc*eff";
    title += Form(" (%g<pt<%g",fPtBinLowEdge[0],fPtBinLowEdge[fPtBinLowEdge.GetSize()-1]);
    title += (iy == 0) ? Form(" / %g<y<%g):",fYBinLowEdge[0],fYBinLowEdge[fYBinLowEdge.GetSize()-1]) : Form(" / %g<y<%g):",fYBinLowEdge[iy-1],fYBinLowEdge[iy]);
    printf("\n%s\n",title.Data());
    IntegratedAccEff(0, iy, nMatch, 0., -1., ae, aee, kTRUE, hGenYSummary, hRecYSummary, hAccYSummary, iy);
    
    // print integrated value in 0-80%
    IntegratedAccEff(0, iy, nMatch, 0., 80., ae, aee);
    printf("   -  %d matching required - 0-80%%: %f ± %f\n", nMatch, ae, aee);
    
    // print integrated value in 0-80% weighted by nColl (bin=10)
    IntegratedAccEff(0, iy, nMatch, 8, nColl10, nColl10binLowEdge, ae, aee);
    printf("   -  %d matching required - 0-80%% - weighted per nColl (bin=10): %f ± %f\n", nMatch, ae, aee);
    
    // print integrated value in 0-80% weighted by nColl (bin=5)
    IntegratedAccEff(0, iy, nMatch, 16, nColl5, nColl5binLowEdge, ae, aee);
    printf("   -  %d matching required - 0-80%% - weighted per nColl (bin=5): %f ± %f\n", nMatch, ae, aee);
    
    // print integrated value in 0-80% weighted by nJPsi integrated over pt/y
    IntegratedAccEff(0, iy, nMatch, 5, nJPsi[0][0], jPsibinLowEdge[0][0], ae, aee);
    printf("   -  %d matching required - 0-80%% - weighted per nJPsi integrated over pt/y: %f ± %f\n", nMatch, ae, aee);
    
  }
  
  // draw summary histos
  cSummary = new TCanvas("cYSummary", "ySummary", 900, 300);
  cSummary->Divide(3,1);
  cSummary->cd(1);
  hGenYSummary->SetStats(kFALSE);
  hGenYSummary->Draw("e0");
  cSummary->cd(2);
  hRecYSummary->SetStats(kFALSE);
  hRecYSummary->Draw("e0");
  cSummary->cd(3);
  hAccYSummary->SetMarkerStyle(21);
  hAccYSummary->SetMarkerSize(0.6);
  hAccYSummary->SetStats(kFALSE);
  hAccYSummary->Draw("e0");
  
  // save histos
  file = new TFile("acceff_new.root","update");
  hGenYSummary->Write(0x0, TObject::kOverwrite);
  hRecYSummary->Write(0x0, TObject::kOverwrite);
  hAccYSummary->Write(0x0, TObject::kOverwrite);
  cSummary->Write(0x0, TObject::kOverwrite);
  file->Close();
  
  // acceptance*efficiency for different pt/y bins
  TH1F* hGenSummary = new TH1F("hGenSummary","generated J/#psi versus pt/y bins;;a.u.", fPtBinLowEdge.GetSize()*fYBinLowEdge.GetSize(), -0.5, fPtBinLowEdge.GetSize()*fYBinLowEdge.GetSize()-0.5);
  TH1F* hRecSummary = new TH1F("hRecSummary","reconstructed J/#psi versus pt/y bins;;a.u.", fPtBinLowEdge.GetSize()*fYBinLowEdge.GetSize(), -0.5, fPtBinLowEdge.GetSize()*fYBinLowEdge.GetSize()-0.5);
  TH1F* hAccSummary = new TH1F("hAccSummary","Acc * Eff versus pt/y bins", fPtBinLowEdge.GetSize()*fYBinLowEdge.GetSize(), -0.5, fPtBinLowEdge.GetSize()*fYBinLowEdge.GetSize()-0.5);
  
  // loop over pt bins+1 (0 = integrated over pt)
  for (Int_t ipt = 0; ipt < fPtBinLowEdge.GetSize(); ipt++) {
    
    // loop over y bins+1 (0 = integrated over y)
    for (Int_t iy = 0; iy < fYBinLowEdge.GetSize(); iy++) {
      
      // acceptance*efficiency versus run
      if (ipt > 0 && iy > 0) DrawAccEffVsRun(ipt, iy);
      
      // acceptance*efficiency versus centrality
      if (ipt > 0 && iy > 0) DrawAccEffVsCent(ipt, iy);
      
      // print integrated value
      Double_t ae, aee;
      TString title = "---- Integrated acc*eff";
      title += (ipt == 0) ? Form(" (%g<pt<%g",fPtBinLowEdge[0],fPtBinLowEdge[fPtBinLowEdge.GetSize()-1]) : Form(" (%g<pt<%g",fPtBinLowEdge[ipt-1],fPtBinLowEdge[ipt]);
      title += (iy == 0) ? Form(" / %g<y<%g):",fYBinLowEdge[0],fYBinLowEdge[fYBinLowEdge.GetSize()-1]) : Form(" / %g<y<%g):",fYBinLowEdge[iy-1],fYBinLowEdge[iy]);
      printf("\n%s\n",title.Data());
      IntegratedAccEff(ipt, iy, nMatch, 0., -1., ae, aee, kTRUE, hGenSummary, hRecSummary, hAccSummary, ipt*fYBinLowEdge.GetSize()+iy);
      
      // print integrated value in 0-80%
      IntegratedAccEff(ipt, iy, nMatch, 0., 80., ae, aee);
      printf("   -  %d matching required - 0-80%%: %f ± %f\n", nMatch, ae, aee);
      
      // print integrated value in 0-80% weighted by nColl (bin=10)
      IntegratedAccEff(ipt, iy, nMatch, 8, nColl10, nColl10binLowEdge, ae, aee);
      printf("   -  %d matching required - 0-80%% - weighted per nColl (bin=10): %f ± %f\n", nMatch, ae, aee);
      
      // print integrated value in 0-80% weighted by nColl (bin=5)
      IntegratedAccEff(ipt, iy, nMatch, 16, nColl5, nColl5binLowEdge, ae, aee);
      printf("   -  %d matching required - 0-80%% - weighted per nColl (bin=5): %f ± %f\n", nMatch, ae, aee);
      
      // print integrated value in 0-80% weighted by nJPsi in pt/y
      if (ipt < 3 && iy < 3) {
	IntegratedAccEff(ipt, iy, nMatch, 5, nJPsi[ipt][iy], jPsibinLowEdge[ipt][iy], ae, aee, kTRUE);
	printf("   -  %d matching required - 0-80%% - weighted per nJPsi in pt/y: %f ± %f\n", nMatch, ae, aee);
      }
      
      // print integrated value in 0-80% weighted by nJPsi integrated over pt/y
      IntegratedAccEff(ipt, iy, nMatch, 5, nJPsi[0][0], jPsibinLowEdge[0][0], ae, aee);
      printf("   -  %d matching required - 0-80%% - weighted per nJPsi integrated over pt/y: %f ± %f\n", nMatch, ae, aee);
      
    }
    
  }
  
  // draw summary histos
  cSummary = new TCanvas("cSummary", "summary", 900, 300);
  cSummary->Divide(3,1);
  cSummary->cd(1);
  hGenSummary->SetStats(kFALSE);
  hGenSummary->Draw("e0");
  cSummary->cd(2);
  hRecSummary->SetStats(kFALSE);
  hRecSummary->Draw("e0");
  cSummary->cd(3);
  hAccSummary->SetMarkerStyle(21);
  hAccSummary->SetMarkerSize(0.6);
  hAccSummary->SetStats(kFALSE);
  hAccSummary->Draw("e0");
  
  // save histos
  file = new TFile("acceff_new.root","update");
  hGenSummary->Write(0x0, TObject::kOverwrite);
  hRecSummary->Write(0x0, TObject::kOverwrite);
  hAccSummary->Write(0x0, TObject::kOverwrite);
  cSummary->Write(0x0, TObject::kOverwrite);
  file->Close();
  
  // print acceptance*efficiency integrated over pt/y versus reduced centrality
  Double_t ae, aee;
  printf("\n---- acc*eff versus centrality (%g<pt<%g / %g<y<%g):\n",
	 fPtBinLowEdge[0],fPtBinLowEdge[fPtBinLowEdge.GetSize()-1],
	 fYBinLowEdge[0],fYBinLowEdge[fYBinLowEdge.GetSize()-1]);
  IntegratedAccEff(0, 0, nMatch, 0, 90, ae, aee, kTRUE);
  for (Int_t i=0; i<5; i++)
    IntegratedAccEff(0, 0, nMatch, jPsibinLowEdge[0][0][i], jPsibinLowEdge[0][0][i+1], ae, aee, kTRUE);
  for (Int_t i=0; i<9; i++)
    IntegratedAccEff(0, 0, nMatch, nColl10binLowEdge[i], nColl10binLowEdge[i+1], ae, aee, kTRUE);
  
}


//________________________________________________________________________
TH1F* AliAnalysisTaskJPsiAccEffCorr::ComputeAccEff(TH1F &hGen, TH1F &hRec,
						   const Char_t *name, const Char_t *title)
{
  /// Compute acc*eff and binomial errors by hand, i.e. not using TGraphAsymmErrors
  /// Result is identical to divide histograms with option "B", except here error is forced > 1/gen
  
  Int_t nbins = hGen.GetNbinsX();
  TH1F* hAcc = new TH1F(name,title, nbins, hGen.GetXaxis()->GetXmin(), hGen.GetXaxis()->GetXmax());
  for (Int_t i = 1; i <= nbins; i++) {
    Double_t gen = hGen.GetBinContent(i);
    Double_t accEff = (gen>0.) ? hRec.GetBinContent(i)/gen : 0.;
    Double_t accEffErr = (gen>0.) ? TMath::Max(1/gen, TMath::Sqrt(accEff*TMath::Abs(1.-accEff)/gen)) : 0.;
    hAcc->SetBinContent(i, accEff);
    hAcc->SetBinError(i, accEffErr);
  }
  
  return hAcc;
}


//________________________________________________________________________
void AliAnalysisTaskJPsiAccEffCorr::IntegratedAccEff(Int_t ipt, Int_t iy, Int_t nMatch,
						     Float_t centMin, Float_t centMax,
						     Double_t &accEff, Double_t &accEffErr, Bool_t print,
						     TH1F* hGenSum, TH1F* hRecSum, TH1F* hAccSum, Int_t bin)
{
  /// compute AccEff correction integrated over the given centrality range
  
  accEff = 0.;
  accEffErr = 0.;
  TString ptKey = (ipt == 0) ? "any" : Form("%g-%g",fPtBinLowEdge[ipt-1],fPtBinLowEdge[ipt]);
  TString yKey = (iy == 0) ? "any" : Form("%g-%g",fYBinLowEdge[iy-1],fYBinLowEdge[iy]);
  
  // set centrality range according to internal binning
  TString centKey = "";
  Int_t imin = 0, imax = 0;
  if (centMax < centMin) {
    
    centKey = "any";
    
  } else {
    
    centMin -= 1.e-6;
    centMax += 1.e-6;
    imin = 0;
    while (fCentBinLowEdge[imin] < centMin && imin < fCentBinLowEdge.GetSize()) imin++;
    imax = fCentBinLowEdge.GetSize()-1;
    while (fCentBinLowEdge[imax] > centMax && imax > 0) imax--;
    if (imin == fCentBinLowEdge.GetSize() || imax == 0 || imax <= imin) {
      AliError("incorrect centrality integration range");
      return;
    }
    if (TMath::Abs(fCentBinLowEdge[imin]-centMin) > 1.e-4 ||
	TMath::Abs(fCentBinLowEdge[imax]-centMax) > 1.e-4)
      AliWarning(Form("centrality integration range needed to be adjusted:\n %g-%g --> %g-%g",
		      centMin+1.e-6, centMax-1.e-6, fCentBinLowEdge[imin], fCentBinLowEdge[imax]));
    
    // integrate over centrality range
    centKey = Form("%g-%g",fCentBinLowEdge[imin],fCentBinLowEdge[imin+1]);
    for (Int_t icent = imin+1; icent < imax; icent++)
      centKey += Form(",%g-%g",fCentBinLowEdge[icent],fCentBinLowEdge[icent+1]);
    
  }
  
  // compute acc*eff
  TH1D *hgen = fJPsiCounters->Get("match", Form("type:gen/ptbin:%s/ybin:%s/cent:%s",ptKey.Data(),yKey.Data(), centKey.Data()));
  hgen->SetName("hgen");
  TH1D *hrec = fJPsiCounters->Get("match", Form("type:rec/ptbin:%s/ybin:%s/cent:%s",ptKey.Data(),yKey.Data(), centKey.Data()));
  hrec->SetName("hrec");
  TGraphAsymmErrors *gacceff = new TGraphAsymmErrors(hrec, hgen);
  
  // get acc*eff value and error according to matching request and print if required
  TString cent = (centMax < centMin) ? "0-100%" : Form("%g-%g%%", fCentBinLowEdge[imin], fCentBinLowEdge[imax]);
  Int_t iMatch = -1;
  Double_t x;
  if (nMatch < 0 || nMatch == 0) {
    iMatch = 0;
    gacceff->GetPoint(0,x,accEff);
    accEffErr = gacceff->GetErrorYhigh(0);
    if (print) printf("   - no matching required - %s: %f + %f - %f\n", cent.Data(), accEff, accEffErr, gacceff->GetErrorYlow(0));
  }
  if (nMatch < 0 || nMatch == 1) {
    iMatch = 1;
    gacceff->GetPoint(1,x,accEff);
    accEffErr = gacceff->GetErrorYhigh(1);
    if (print) printf("   -  1 matching required - %s: %f + %f - %f\n", cent.Data(), accEff, accEffErr, gacceff->GetErrorYlow(1));
  }
  if (nMatch < 0 || nMatch == 2) {
    iMatch = 2;
    gacceff->GetPoint(2,x,accEff);
    accEffErr = gacceff->GetErrorYhigh(2);
    if (print) printf("   -  2 matching required - %s: %f + %f - %f\n", cent.Data(), accEff, accEffErr, gacceff->GetErrorYlow(2));
  }
  
  // fill summary histos if require
  if (hGenSum && hRecSum && hAccSum) {
    TString label = "";
    label += (ipt == 0) ? Form("%g<pt<%g",fPtBinLowEdge[0],fPtBinLowEdge[fPtBinLowEdge.GetSize()-1]) : Form("%g<pt<%g",fPtBinLowEdge[ipt-1],fPtBinLowEdge[ipt]);
    label += (iy == 0) ? Form(" / %g<y<%g",fYBinLowEdge[0],fYBinLowEdge[fYBinLowEdge.GetSize()-1]) : Form(" / %g<y<%g",fYBinLowEdge[iy-1],fYBinLowEdge[iy]);
    hGenSum->Fill(bin, hgen->GetBinContent(iMatch+1));
    hGenSum->GetXaxis()->SetBinLabel(bin+1, label.Data());
    hRecSum->Fill(bin, hrec->GetBinContent(iMatch+1));
    hRecSum->GetXaxis()->SetBinLabel(bin+1, label.Data());
    hAccSum->SetBinContent(bin+1, accEff);
    hAccSum->SetBinError(bin+1, accEffErr);
    hAccSum->GetXaxis()->SetBinLabel(bin+1, label.Data());
  }
  
  // clean memory
  delete hgen;
  delete hrec;
  delete gacceff;
  
}


//________________________________________________________________________
void AliAnalysisTaskJPsiAccEffCorr::IntegratedAccEff(Int_t ipt, Int_t iy, Int_t nMatch,
						     Int_t nBins, Float_t *nJPsi, Float_t *centbinLowEdge,
						     Double_t &accEff, Double_t &accEffErr, Bool_t print)
{
  /// compute integrated AccEff correction weighted by the number of JPsi in each centrality bin
  
  Double_t nJPsiTot = 0., nJPsiCorrTot = 0., tmpErr = 0.;
  for (Int_t i=0; i<nBins; i++) {
    if (nJPsi[i] > 0.) {
      Double_t accEffi = -1., accEffErri = -1.;
      IntegratedAccEff(ipt, iy, nMatch, centbinLowEdge[i], centbinLowEdge[i+1], accEffi, accEffErri, print);
      if (accEffi <= 0.) ;//AliWarning(Form("acc*eff = 0 in centrality bin %g-%g", centbinLowEdge[i], centbinLowEdge[i+1]));
      else {
	nJPsiTot += nJPsi[i];
	nJPsiCorrTot += nJPsi[i]/accEffi;
	tmpErr += nJPsi[i]*nJPsi[i]*accEffErri*accEffErri/accEffi/accEffi/accEffi/accEffi;
      }
    }
  }
  
  accEff = (nJPsiCorrTot > 0.) ? nJPsiTot/nJPsiCorrTot : 0.;
  accEffErr = (nJPsiCorrTot > 0.) ? accEff*TMath::Sqrt(tmpErr)/nJPsiCorrTot : 0.;
  
}


//________________________________________________________________________
void AliAnalysisTaskJPsiAccEffCorr::DrawAccEffVsRun(Int_t ipt, Int_t iy)
{
  /// Draw acceptance*efficiency versus run for this given pt/y bin
  
  TString ptKey = (ipt == 0) ? "any" : Form("%g-%g",fPtBinLowEdge[ipt-1],fPtBinLowEdge[ipt]);
  TString yKey = (iy == 0) ? "any" : Form("%g-%g",fYBinLowEdge[iy-1],fYBinLowEdge[iy]);
  
  TH1D *hgen = fJPsiCounters->Get("run", Form("type:gen/match:any/ptbin:%s/ybin:%s",ptKey.Data(),yKey.Data()));
  hgen->SetName(Form("hgenVsRun_pt%d_y%d",ipt,iy));
  
  TH1D *hrec = fJPsiCounters->Get("run", Form("type:rec/match:2/ptbin:%s/ybin:%s",ptKey.Data(),yKey.Data()));
  hrec->SetName(Form("hrecVsRun_pt%d_y%d",ipt,iy));
  
  TGraphAsymmErrors *acceff = new TGraphAsymmErrors(hrec, hgen);
  
  acceff->GetXaxis()->Set(hgen->GetNbinsX(), hgen->GetXaxis()->GetXmin(), hgen->GetXaxis()->GetXmax());
  TIter nextRun(hgen->GetXaxis()->GetLabels());
  TObjString *srun = 0x0;
  Int_t irun = 1;
  while ((srun = static_cast<TObjString*>(nextRun()))) acceff->GetXaxis()->SetBinLabel(irun++, srun->GetName());
  acceff->GetXaxis()->SetNdivisions(1,kFALSE);
  
  // draw histos
  TString title = "acceptance * efficiency versus run";
  title += (ipt == 0) ? Form(" (%g<pt<%g",fPtBinLowEdge[0],fPtBinLowEdge[fPtBinLowEdge.GetSize()-1]) : Form(" (%g<pt<%g",fPtBinLowEdge[ipt-1],fPtBinLowEdge[ipt]);
  title += (iy == 0) ? Form(" / %g<y<%g)",fYBinLowEdge[0],fYBinLowEdge[fYBinLowEdge.GetSize()-1]) : Form(" / %g<y<%g)",fYBinLowEdge[iy-1],fYBinLowEdge[iy]);
  new TCanvas(Form("cAccEffVsRun_pt%d_y%d",ipt,iy), title.Data(), 1000, 300);
  acceff->SetNameTitle(Form("accEffVsRun_pt%d_y%d",ipt,iy), title.Data());
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
  
  // save histos
  TFile *file = new TFile("acceff_new.root","update");
  hgen->Write(0x0, TObject::kOverwrite);
  hrec->Write(0x0, TObject::kOverwrite);
  acceff->Write(0x0, TObject::kOverwrite);
  file->Close();
  
  // clean memory
  delete hgen;
  delete hrec;
}


//________________________________________________________________________
void AliAnalysisTaskJPsiAccEffCorr::DrawAccEffVsCent(Int_t ipt, Int_t iy)
{
  /// Draw acceptance*efficiency versus centrality for this given pt/y bin
  
  TString ptKey = (ipt == 0) ? "any" : Form("%g-%g",fPtBinLowEdge[ipt-1],fPtBinLowEdge[ipt]);
  TString yKey = (iy == 0) ? "any" : Form("%g-%g",fYBinLowEdge[iy-1],fYBinLowEdge[iy]);
  
  TH1D *hgen = fJPsiCounters->Get("cent", Form("type:gen/match:any/ptbin:%s/ybin:%s",ptKey.Data(),yKey.Data()));
  hgen->SetName(Form("hgenVsCent_pt%d_y%d",ipt,iy));
  
  TH1D *hrec = fJPsiCounters->Get("cent", Form("type:rec/match:2/ptbin:%s/ybin:%s",ptKey.Data(),yKey.Data()));
  hrec->SetName(Form("hrecVsCent_pt%d_y%d",ipt,iy));
  
  TGraphAsymmErrors *acceff = new TGraphAsymmErrors(hrec, hgen);
  
  acceff->GetXaxis()->Set(hgen->GetNbinsX(), hgen->GetXaxis()->GetXmin(), hgen->GetXaxis()->GetXmax());
  TIter nextCent(hgen->GetXaxis()->GetLabels());
  TObjString *scent = 0x0;
  Int_t icent = 1;
  while ((scent = static_cast<TObjString*>(nextCent()))) acceff->GetXaxis()->SetBinLabel(icent++, scent->GetName());
  acceff->GetXaxis()->SetNdivisions(1,kFALSE);
  
  // draw histos
  TString title = "acceptance * efficiency versus centrality";
  title += (ipt == 0) ? Form(" (%g<pt<%g",fPtBinLowEdge[0],fPtBinLowEdge[fPtBinLowEdge.GetSize()-1]) : Form(" (%g<pt<%g",fPtBinLowEdge[ipt-1],fPtBinLowEdge[ipt]);
  title += (iy == 0) ? Form(" / %g<y<%g)",fYBinLowEdge[0],fYBinLowEdge[fYBinLowEdge.GetSize()-1]) : Form(" / %g<y<%g)",fYBinLowEdge[iy-1],fYBinLowEdge[iy]);
  new TCanvas(Form("cAccEffVsCent_pt%d_y%d",ipt,iy), title.Data(), 900, 300);
  acceff->SetNameTitle(Form("accEffVsCent_pt%d_y%d",ipt,iy), title.Data());
  acceff->SetLineStyle(1);
  acceff->SetLineColor(1); 
  acceff->SetMarkerStyle(20);
  acceff->SetMarkerSize(0.7);
  acceff->SetMarkerColor(2);
  acceff->GetXaxis()->SetTitle("Centrality");
  acceff->GetXaxis()->SetLabelFont(22);
  acceff->GetXaxis()->SetTitleFont(22);
  acceff->GetYaxis()->SetTitle("Acc*Eff");
  acceff->GetYaxis()->SetLabelFont(22);
  acceff->GetYaxis()->SetLabelFont(22);
  acceff->Draw("ap");
  
  // save histos
  TFile *file = new TFile("acceff_new.root","update");
  hgen->Write(0x0, TObject::kOverwrite);
  hrec->Write(0x0, TObject::kOverwrite);
  acceff->Write(0x0, TObject::kOverwrite);
  file->Close();
  
  // clean memory
  delete hgen;
  delete hrec;
}

