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
//#include <TRandom3.h>

// ROOT includes
#include <TMath.h>
#include <TH1F.h>
#include <TParameter.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TClonesArray.h>
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TList.h>
#include <THashList.h>

// STEER includes
#include "AliLog.h"
#include "AliVEvent.h"
#include "AliMCEvent.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliVParticle.h"
#include "AliMultSelection.h"
#include "AliCounterCollection.h"

// ANALYSIS includes
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisDataSlot.h"
#include "AliAnalysisDataContainer.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisMuonUtility.h"
#include "AliAnalysisTaskJPsiAccEffCorr2.h"


//========================================================================
ClassImp(Weights)

//________________________________________________________________________
Weights::Weights() :
TNamed(),
fPtMin(0.),
fPtMax(0.),
fYMin(0.),
fYMax(0.),
fnCentBin(0),
fBinLowEdge(0),
fnSignal(0),
fWeightRec(kFALSE),
fWeights()
{
  /// Dummy constructor
  fWeights.SetOwner();
}

//________________________________________________________________________
Weights::Weights(const Char_t *name, Float_t ptMin, Float_t ptMax, Float_t yMin, Float_t yMax,
		 Int_t nCentBin, Float_t *binLowEdge, Double_t *nSig, Bool_t weightRec) :
TNamed(name, Form("%g<pt<%g / %g<y<%g", ptMin, ptMax, yMin, yMax)),
fPtMin(ptMin),
fPtMax(ptMax),
fYMin(yMin),
fYMax(yMax),
fnCentBin(nCentBin),
fBinLowEdge(nCentBin+1, binLowEdge),
fnSignal(nCentBin, nSig),
fWeightRec(weightRec),
fWeights(fnCentBin)
{
  /// Constructor
  fWeights.SetOwner();
}

//________________________________________________________________________
void Weights::Init(TArrayF &centBinLowEdge)
{
  /// Set the array of weights associated to centrality key words using the given centrality bins
  
  for (Int_t i=0; i<fnCentBin; i++) {
    
    if (fnSignal[i] <= 0.) continue;
    
    // set centrality range according to given binning
    TString centKey = GetCentKey(fBinLowEdge[i], fBinLowEdge[i+1], centBinLowEdge);
    
    // add the weight with the associated centrality key word
    if (!centKey.IsNull()) fWeights.Add(new TParameter<Double_t>(Form("cent:%s",centKey.Data()), fnSignal[i]));
    
  }
  
}

//________________________________________________________________________
TString Weights::GetCentKey(Float_t centMin, Float_t centMax, const TArrayF &centBinLowEdge)
{
  /// Return the centrality key work (or list of key words) for the given range
  
  TString centKey = "";
  Int_t imin = 0, imax = 0;
  if (centMax < centMin) {
    
    centKey = "any";
    
  } else {
    
    centMin -= 1.e-6;
    centMax += 1.e-6;
    imin = 0;
    while (centBinLowEdge[imin] < centMin && imin < centBinLowEdge.GetSize()) imin++;
    imax = centBinLowEdge.GetSize()-1;
    while (centBinLowEdge[imax] > centMax && imax > 0) imax--;
    if (imin == centBinLowEdge.GetSize() || imax == 0 || imax <= imin) {
      AliErrorClass("incorrect centrality integration range");
      return "";
    }
    if (TMath::Abs(centBinLowEdge[imin]-centMin) > 1.e-4 ||
	TMath::Abs(centBinLowEdge[imax]-centMax) > 1.e-4)
      AliWarningClass(Form("centrality integration range needed to be adjusted:\n %g-%g --> %g-%g",
		      centMin+1.e-6, centMax-1.e-6, centBinLowEdge[imin], centBinLowEdge[imax]));
    
    // integrate over centrality range
    centKey = Form("%g-%g",centBinLowEdge[imin],centBinLowEdge[imin+1]);
    for (Int_t icent = imin+1; icent < imax; icent++)
      centKey += Form(",%g-%g",centBinLowEdge[icent],centBinLowEdge[icent+1]);
    
  }
  
  return centKey;
}

//========================================================================
ClassImp(AliAnalysisTaskJPsiAccEffCorr2)

//________________________________________________________________________
AliAnalysisTaskJPsiAccEffCorr2::AliAnalysisTaskJPsiAccEffCorr2() :
AliAnalysisTaskSE(),
fList(0x0),
fEventCounters(0x0),
fJPsiCounters(0x0),
fCentBinLowEdge(0),
fPtBinLowEdge(0),
fYBinLowEdge(0),
fTrigLevel(1),
fNMatch(2),
fMuLowPtCut(-1.),
fSigWeights(0x0),
fRunWeights(0x0)
{
  /// Dummy constructor
}

//________________________________________________________________________
AliAnalysisTaskJPsiAccEffCorr2::AliAnalysisTaskJPsiAccEffCorr2(const char *name) :
AliAnalysisTaskSE(name),
fList(0x0),
fEventCounters(0x0),
fJPsiCounters(0x0),
fCentBinLowEdge(0),
fPtBinLowEdge(0),
fYBinLowEdge(0),
fTrigLevel(1),
fNMatch(2),
fMuLowPtCut(-1.),
fSigWeights(0x0),
fRunWeights(0x0)
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
AliAnalysisTaskJPsiAccEffCorr2::~AliAnalysisTaskJPsiAccEffCorr2()
{
  /// Destructor
  if (!AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
    delete fList;
    delete fEventCounters;
    delete fJPsiCounters;
  }
  delete fSigWeights;
  delete fRunWeights;
}

//___________________________________________________________________________
void AliAnalysisTaskJPsiAccEffCorr2::UserCreateOutputObjects()
{
  /// Create histograms and counters
  
  //gRandom->SetSeed(0);
  
  // initialize histos
  fList = new TObjArray(2000);
  fList->SetOwner();
  
  TH1F* hPtGen = new TH1F("hPtGen","generated J/#psi versus pt;p_{T} (GeV/c);N_{J/#psi} / 0.5 GeV/c", 60, 0., 30.);
  fList->AddAtAndExpand(hPtGen, kPtGen);
  TH1F* hPtRec = new TH1F("hPtRec","reconstructed J/#psi versus pt;p_{T} (GeV/c);N_{J/#psi} / 0.5 GeV/c", 60, 0., 30.);
  fList->AddAtAndExpand(hPtRec, kPtRec);
  
  TH1F* hYGen = new TH1F("hYGen","generated J/#psi versus y;y;N_{J/#psi}", 15, -4., -2.5);
  fList->AddAtAndExpand(hYGen, kYGen);
  TH1F* hYRec = new TH1F("hYRec","reconstructed J/#psi versus y;y;N_{J/#psi}", 15, -4., -2.5);
  fList->AddAtAndExpand(hYRec, kYRec);
  
  TH1F* hPtGenMu = new TH1F("hPtGenMu","generated single-#mu versus pt;p_{T} (GeV/c);N_{single-#mu} / 0.5 GeV/c", 60, 0., 30.);
  fList->AddAtAndExpand(hPtGenMu, kPtGenMu);
  TH1F* hPtRecMu = new TH1F("hPtRecMu","reconstructed single-#mu versus pt;p_{T} (GeV/c);N_{single-#mu} / 0.5 GeV/c", 60, 0., 30.);
  fList->AddAtAndExpand(hPtRecMu, kPtRecMu);
  
  TH1F* hYGenMu = new TH1F("hYGenMu","generated single-#mu versus y;y;N_{single-#mu}", 15, -4., -2.5);
  fList->AddAtAndExpand(hYGenMu, kYGenMu);
  TH1F* hYRecMu = new TH1F("hYRecMu","reconstructed single-#mu versus y;y;N_{single-#mu}", 15, -4., -2.5);
  fList->AddAtAndExpand(hYRecMu, kYRecMu);
  
  TH1F* hDzVtx = new TH1F("hDzVtx","vertex resolution;#DeltaZ (cm)", 200, -1., 1.);
  fList->AddAtAndExpand(hDzVtx, kDzVtx);
  
  TH1F* hDzVtxJPsi = new TH1F("hDzVtxJPsi","vertex resolution for event with J/#psi;#DeltaZ (cm)", 200, -1., 1.);
  fList->AddAtAndExpand(hDzVtxJPsi, kDzVtx2);
  
  TH1F* hMass = new TH1F("hMass","invariant mass distribution;Mass (GeV/c^{2});N_{J/#psi} / 0.025 GeV/c", 200, 0., 5.);
  fList->AddAtAndExpand(hMass, kMass);
  
  TH1F* hDphi = new TH1F("hDphi","J/#psi #phi resolution;#Delta#phi (deg);N_{J/#psi}", 720, -180., 180.);
  fList->AddAtAndExpand(hDphi, kDphi);
  
  TH1F* hCos2Dphi = new TH1F("hCos2Dphi","J/#psi cos(2.#times(#phi_{rec}-#phi_{MC})) distribution;cos(2.#times(#phi_{rec}-#phi_{MC}));N_{J/#psi}", 200, -1., 1.);
  fList->AddAtAndExpand(hCos2Dphi, kCos2Dphi);
  
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
  
  // Post data at least once per task to ensure data synchronisation (required for merging)
  PostData(1, fList);
  PostData(2, fEventCounters);
  PostData(3, fJPsiCounters);
}

//________________________________________________________________________
void AliAnalysisTaskJPsiAccEffCorr2::UserExec(Option_t *)
{
  /// Called for each event
  
  // get event
  AliVEvent *evt = InputEvent();
  if (!dynamic_cast<AliESDEvent*>(evt) && !dynamic_cast<AliAODEvent*>(evt)) return;
  
  // get MC event
  AliMCEvent *mcEvt = MCEvent();
  if (!mcEvt) return;
  
  // get the centrality percentile
  AliMultSelection *multSelection = static_cast<AliMultSelection*>(evt->FindListObject("MultSelection"));
  Float_t centrality = multSelection ? multSelection->GetMultiplicityPercentile("V0M") : -1.;
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
  Double_t mcPhi = 999.;
  
  // ------ MC part ------
  Int_t nMCTracks = mcEvt->GetNumberOfTracks();
  /*
  // remove events with JPsi outside -4.2<y<-2.3
  for (Int_t iMCTrack = 0; iMCTrack < nMCTracks; ++iMCTrack) {
    AliVParticle *mctrack = mcEvt->GetTrack(iMCTrack);
    if (AliAnalysisMuonUtility::IsPrimary(mctrack,mcEvt) && mctrack->PdgCode() == 443 &&
        (mctrack->Y() < -4.2 || mctrack->Y() > -2.3)) return;
  }
  */
  // fill event counter
  fEventCounters->Count(Form("event:any/cent:any/run:%d",fCurrentRunNumber));
  if (!centKey.IsNull()) fEventCounters->Count(Form("event:any/cent:%s/run:%d",centKey.Data(),fCurrentRunNumber));
  
  // vertex resolution
  Double_t zVtx = evt->GetPrimaryVertex()->GetZ();
  Double_t zVtxMC = AliAnalysisMuonUtility::GetMCVertexZ(evt,mcEvt);
  ((TH1F*)fList->UncheckedAt(kDzVtx))->Fill(zVtx-zVtxMC);
  
  for (Int_t iMCTrack = 0; iMCTrack < nMCTracks; ++iMCTrack) {
    
    AliVParticle *mctrack = mcEvt->GetTrack(iMCTrack);
    
    Double_t pT = mctrack->Pt();
    Double_t eta = mctrack->Eta();
    Double_t y = mctrack->Y();
    
    // look for muons
    if (TMath::Abs(mctrack->PdgCode()) == 13 &&
	eta >= fYBinLowEdge[0] && eta <= fYBinLowEdge[fYBinLowEdge.GetSize()-1]) {
      
      ((TH1F*)fList->UncheckedAt(kPtGenMu))->Fill(pT);
      ((TH1F*)fList->UncheckedAt(kYGenMu))->Fill(y);
      
    }
    
    // look for generated particles
    if(AliAnalysisMuonUtility::IsPrimary(mctrack,mcEvt) && mctrack->PdgCode() == 443 &&
       y >= fYBinLowEdge[0] && y < fYBinLowEdge[fYBinLowEdge.GetSize()-1]) {
      
      mcPhi = mctrack->Phi();
      
      ((TH1F*)fList->UncheckedAt(kPtGen))->Fill(pT);
      ((TH1F*)fList->UncheckedAt(kYGen))->Fill(y);
      
      // pt bin
      TString ptKey = "";
      for (Int_t ipt = 0; ipt < fPtBinLowEdge.GetSize()-1; ipt++)
	if (pT >= fPtBinLowEdge[ipt] && pT < fPtBinLowEdge[ipt+1])
	  ptKey = Form("%g-%g",fPtBinLowEdge[ipt],fPtBinLowEdge[ipt+1]);
      if (ptKey.IsNull()) continue;
      
      // y bin
      TString yKey = "";
      for (Int_t iy = 0; iy < fYBinLowEdge.GetSize()-1; iy++)
	if (y >= fYBinLowEdge[iy] && y < fYBinLowEdge[iy+1])
	  yKey = Form("%g-%g",fYBinLowEdge[iy],fYBinLowEdge[iy+1]);
      if (yKey.IsNull()) continue;
      
      for (Int_t itrg = 0; itrg < 3; itrg++) {
	fJPsiCounters->Count(Form("type:gen/match:%s/ptbin:%s/ybin:%s/cent:any/run:%d",
				  matchKey[itrg].Data(),ptKey.Data(),yKey.Data(),fCurrentRunNumber));
	if (!centKey.IsNull()) fJPsiCounters->Count(Form("type:gen/match:%s/ptbin:%s/ybin:%s/cent:%s/run:%d",
							 matchKey[itrg].Data(),ptKey.Data(),yKey.Data(),centKey.Data(),fCurrentRunNumber));
      }
      
    }
    
  }
  
  // ------ Reco part ------
  Bool_t jpsiFound = kFALSE;
  Int_t nTracks = AliAnalysisMuonUtility::GetNTracks(evt);
  for (Int_t iTrack1 = 0; iTrack1 < nTracks; iTrack1++) {
    
    AliVParticle *track1 = AliAnalysisMuonUtility::GetTrack(iTrack1, evt);
    if (!AliAnalysisMuonUtility::IsMuonTrack(track1)) continue;
    
    if (track1->GetLabel() < 0) continue;
    
    TLorentzVector muV1(track1->Px(), track1->Py(), track1->Pz(), track1->E());
    Short_t charge1 = track1->Charge();
    Int_t matchTrig1 = AliAnalysisMuonUtility::GetMatchTrigger(track1);
    Double_t eta1 = track1->Eta();
    Double_t rAbs1 = AliAnalysisMuonUtility::GetRabs(track1);
    //Double_t thetaAbs1 = AliAnalysisMuonUtility::GetThetaAbsDeg(track1);
    Double_t pT1 = track1->Pt();
    
    // fill single muon information
    if (matchTrig1 >= fTrigLevel &&
	eta1 >= fYBinLowEdge[0] && eta1 <= fYBinLowEdge[fYBinLowEdge.GetSize()-1] &&
	rAbs1 >= rAbsMin && rAbs1 <= rAbsMax) {
      
      ((TH1F*)fList->UncheckedAt(kPtRecMu))->Fill(pT1);
      ((TH1F*)fList->UncheckedAt(kYRecMu))->Fill(track1->Y());
      
    }
    
    for (Int_t iTrack2 = iTrack1 + 1; iTrack2 < nTracks; iTrack2++) {
      
      AliVParticle *track2 = AliAnalysisMuonUtility::GetTrack(iTrack2, evt);
      if (!AliAnalysisMuonUtility::IsMuonTrack(track2)) continue;
      
      if (track2->GetLabel() < 0) continue;
      
      TLorentzVector muV2(track2->Px(), track2->Py(), track2->Pz(), track2->E());
      Short_t charge2 = track2->Charge();
      Int_t matchTrig2 = AliAnalysisMuonUtility::GetMatchTrigger(track2);
      Double_t eta2 = track2->Eta();
      Double_t rAbs2 = AliAnalysisMuonUtility::GetRabs(track2);
      //Double_t thetaAbs2 = AliAnalysisMuonUtility::GetThetaAbsDeg(track2);
      Double_t pT2 = track2->Pt();
      
      TLorentzVector dimuV = muV1 + muV2;
      Short_t charge = charge1*charge2;
      Double_t pT = dimuV.Pt();
      Double_t y = dimuV.Rapidity();
      
      // fill dimuon information
      //Double_t ptCut1 = gRandom->Gaus(fMuLowPtCut,0.333);
      //Double_t ptCut2 = gRandom->Gaus(fMuLowPtCut,0.333);
      if(charge < 0 &&
         y >= fYBinLowEdge[0] && y < fYBinLowEdge[fYBinLowEdge.GetSize()-1] &&
         eta1 >= -4. && eta1 <= -2.5 && eta2 >= -4. && eta2 <= -2.5 &&
         //pT1 > ptCut1 && pT2 > ptCut2 &&
         pT1 > fMuLowPtCut && pT2 > fMuLowPtCut &&
         rAbs1 >= rAbsMin && rAbs1 <= rAbsMax && rAbs2 >= rAbsMin && rAbs2 <= rAbsMax
         //thetaAbs1 > 2. && thetaAbs1 < 10. && thetaAbs2 > 2. && thetaAbs2 < 10.
         ) {
        /*
        // remove events with muons generated outside 168.5<Theta<178.5
        Bool_t accOk = kTRUE;
        AliVParticle *mu[2] = {track1, track2};
        for (Int_t imu = 0; imu < 2; imu++) {
          Int_t label = mu[imu]->GetLabel();
          AliVParticle *mctrack = mcEvt->GetTrack(label);
          if (mctrack && mctrack->PdgCode() == 13 && mctrack->GetMother() >= 0 &&
              mcEvt->GetTrack(mctrack->GetMother())->PdgCode() == 443 &&
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
        if(matchTrig1 >= fTrigLevel || matchTrig2 >= fTrigLevel) trigOk[1] = kTRUE;
        if(matchTrig1 >= fTrigLevel && matchTrig2 >= fTrigLevel) trigOk[2] = kTRUE;
        
        if (trigOk[2]) {
          ((TH1F*)fList->UncheckedAt(kPtRec))->Fill(pT);
          ((TH1F*)fList->UncheckedAt(kYRec))->Fill(y);
          ((TH1F*)fList->UncheckedAt(kMass))->Fill(dimuV.M());
          Double_t phi = dimuV.Phi();
          if (phi < 0.) phi += 2.*TMath::Pi();
          Double_t dPhi = (phi-mcPhi);
          if (dPhi < -TMath::Pi()) dPhi += 2.*TMath::Pi();
          else if (dPhi > TMath::Pi()) dPhi -= 2.*TMath::Pi();
          ((TH1F*)fList->UncheckedAt(kDphi))->Fill(dPhi*TMath::RadToDeg());
          if (dPhi > -0.5*TMath::Pi() && dPhi < 0.5*TMath::Pi()) ((TH1F*)fList->UncheckedAt(kCos2Dphi))->Fill(TMath::Cos(2.*dPhi));
          //((TH1F*)fList->UncheckedAt(kPtRecMu))->Fill(pT1);
          //((TH1F*)fList->UncheckedAt(kYRecMu))->Fill(track1->Y());
          //((TH1F*)fList->UncheckedAt(kPtRecMu))->Fill(pT2);
          //((TH1F*)fList->UncheckedAt(kYRecMu))->Fill(track2->Y());
          jpsiFound = kTRUE;
        }
        
        // pt bin
        TString ptKey;
        for (Int_t ipt = 0; ipt < fPtBinLowEdge.GetSize()-1; ipt++)
          if (pT >= fPtBinLowEdge[ipt] && pT < fPtBinLowEdge[ipt+1])
            ptKey = Form("%g-%g",fPtBinLowEdge[ipt],fPtBinLowEdge[ipt+1]);
        if (ptKey.IsNull()) continue;
        
        // y bin
        TString yKey;
        for (Int_t iy = 0; iy < fYBinLowEdge.GetSize()-1; iy++)
          if (y >= fYBinLowEdge[iy] && y < fYBinLowEdge[iy+1])
            yKey = Form("%g-%g",fYBinLowEdge[iy],fYBinLowEdge[iy+1]);
        if (yKey.IsNull()) continue;
        
        for (Int_t itrg = 0; itrg < 3; itrg++) {
          if (!trigOk[itrg]) continue;
          fJPsiCounters->Count(Form("type:rec/match:%s/ptbin:%s/ybin:%s/cent:any/run:%d",
                                    matchKey[itrg].Data(),ptKey.Data(),yKey.Data(),fCurrentRunNumber));
          if (!centKey.IsNull()) fJPsiCounters->Count(Form("type:rec/match:%s/ptbin:%s/ybin:%s/cent:%s/run:%d",
                                                           matchKey[itrg].Data(),ptKey.Data(),yKey.Data(),centKey.Data(),fCurrentRunNumber));
        }
        
      }
      
    }
    
  }
  
  // vertex resolution for events with a valid reconstructed JPsi
  if (jpsiFound) ((TH1F*)fList->UncheckedAt(kDzVtx2))->Fill(zVtx-zVtxMC);
  
  // Post final data. It will be written to a file with option "RECREATE"
  PostData(1, fList);
  PostData(2, fEventCounters);
  PostData(3, fJPsiCounters);
  
}

//________________________________________________________________________
void AliAnalysisTaskJPsiAccEffCorr2::Terminate(Option_t *)
{
  /// draw final results
  
  printf("\n\n\t\t===============\n");
  printf("\t\t||  Results  ||\n");
  printf("\t\t===============\n\n");
  
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
  
  // init the signal weights
  TIter nextSigWeight(fSigWeights);
  Weights *sigWeights = 0x0;
  while ((sigWeights = static_cast<Weights*>(nextSigWeight()))) sigWeights->Init(fCentBinLowEdge);
  
  // clean the run weights
  CleanRunWeights();
  
  // number of events
  printf("Total number of events = %d\n", nEv);
  TH1D *hevVsRun = fEventCounters->Get("run","");
  TH1D *hevVsCent = fEventCounters->Get("cent","");
  TFile* file = new TFile("acceff_new.root","recreate");
  if (hevVsRun) hevVsRun->Write(0x0, TObject::kOverwrite);
  if (hevVsCent) hevVsCent->Write(0x0, TObject::kOverwrite);
  file->Close();
  delete hevVsRun;
  delete hevVsCent;
  
  // compute acc*eff vs. pt (J/Psi)
  TH1F* hPtGen = (TH1F*)fList->UncheckedAt(kPtGen);
  TH1F* hPtRec = (TH1F*)fList->UncheckedAt(kPtRec);
  TH1F* hPtAcc = ComputeAccEff(*hPtGen, *hPtRec, "hPtAcc","Acc * Eff versus pt (J/#psi);p_{T} (GeV/c);Acc*Eff");
  
  // compute acc*eff vs. y (J/Psi)
  TH1F* hYGen = (TH1F*)fList->UncheckedAt(kYGen);
  TH1F* hYRec = (TH1F*)fList->UncheckedAt(kYRec);
  TH1F* hYAcc = ComputeAccEff(*hYGen, *hYRec, "hYAcc","Acc * Eff versus y (J/#psi);y;Acc*Eff");
  
  // compute acc*eff vs. pt (single mu)
  TH1F* hPtGenMu = (TH1F*)fList->UncheckedAt(kPtGenMu);
  TH1F* hPtRecMu = (TH1F*)fList->UncheckedAt(kPtRecMu);
  TH1F* hPtAccMu = ComputeAccEff(*hPtGenMu, *hPtRecMu, "hPtAccMu","Acc * Eff versus pt (single-#mu);p_{T} (GeV/c);Acc*Eff");
  
  // compute acc*eff vs. y (single mu)
  TH1F* hYGenMu = (TH1F*)fList->UncheckedAt(kYGenMu);
  TH1F* hYRecMu = (TH1F*)fList->UncheckedAt(kYRecMu);
  TH1F* hYAccMu = ComputeAccEff(*hYGenMu, *hYRecMu, "hYAccMu","Acc * Eff versus y (single-#mu);y;Acc*Eff");
  
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
  
  // acceptance*efficiency for different pt bins
  TString unit = fRunWeights ? "a.u." : "N_{J/#psi}";
  TH1F* hGenPtSummary = new TH1F("hGenPtSummary",Form("generated J/#psi versus pt bins;;%s",unit.Data()),
				 fPtBinLowEdge.GetSize(), -0.5, fPtBinLowEdge.GetSize()-0.5);
  TH1F* hRecPtSummary = new TH1F("hRecPtSummary",Form("reconstructed J/#psi versus pt bins;;%s",unit.Data()),
				 fPtBinLowEdge.GetSize(), -0.5, fPtBinLowEdge.GetSize()-0.5);
  TH1F* hAccPtSummary = new TH1F("hAccPtSummary","Acc * Eff versus pt bins;;Acc*Eff", fPtBinLowEdge.GetSize(), -0.5,
				 fPtBinLowEdge.GetSize()-0.5);
  Double_t gen[2], rec[2], acc[3];
  
  // acceptance*efficiency for different y bins
  TH1F* hGenYSummary = new TH1F("hGenYSummary",Form("generated J/#psi versus y bins;;%s",unit.Data()),
				fYBinLowEdge.GetSize(), -0.5, fYBinLowEdge.GetSize()-0.5);
  TH1F* hRecYSummary = new TH1F("hRecYSummary",Form("reconstructed J/#psi versus y bins;;%s",unit.Data()),
				fYBinLowEdge.GetSize(), -0.5, fYBinLowEdge.GetSize()-0.5);
  TH1F* hAccYSummary = new TH1F("hAccYSummary","Acc * Eff versus y bins;;Acc*Eff", fYBinLowEdge.GetSize(), -0.5,
				fYBinLowEdge.GetSize()-0.5);
  
  // acceptance*efficiency for different pt/y bins
  TH1F* hGenSummary = new TH1F("hGenSummary",Form("generated J/#psi versus pt/y bins;;%s",unit.Data()),
			       fPtBinLowEdge.GetSize()*fYBinLowEdge.GetSize(), -0.5, fPtBinLowEdge.GetSize()*fYBinLowEdge.GetSize()-0.5);
  TH1F* hRecSummary = new TH1F("hRecSummary",Form("reconstructed J/#psi versus pt/y bins;;%s",unit.Data()),
			       fPtBinLowEdge.GetSize()*fYBinLowEdge.GetSize(), -0.5, fPtBinLowEdge.GetSize()*fYBinLowEdge.GetSize()-0.5);
  TH1F* hAccSummary = new TH1F("hAccSummary","Acc * Eff versus pt/y bins;;Acc*Eff",
			       fPtBinLowEdge.GetSize()*fYBinLowEdge.GetSize(), -0.5, fPtBinLowEdge.GetSize()*fYBinLowEdge.GetSize()-0.5);
  
  // loop over pt bins+1 (0 = integrated over pt)
  for (Int_t ipt = 0; ipt < fPtBinLowEdge.GetSize(); ipt++) {
    
    Float_t ptMin = (ipt == 0) ? fPtBinLowEdge[0] : fPtBinLowEdge[ipt-1];
    Float_t ptMax = (ipt == 0) ? fPtBinLowEdge[fPtBinLowEdge.GetSize()-1] : fPtBinLowEdge[ipt];
    
    // print integrated value
    printf("\n---- Integrated acc*eff (%g<pt<%g / %g<y<%g):\n", ptMin, ptMax, fYBinLowEdge[0], fYBinLowEdge[fYBinLowEdge.GetSize()-1]);
    
    // in 0-100% not weighted over centrality
    IntegratedAccEff(ipt, 0, 0., -1., fNMatch, gen, rec, acc);
    FillHistos(ipt, 0, gen, rec, acc, hGenPtSummary, hRecPtSummary, hAccPtSummary, ipt+1);
    if (ipt == 0) FillHistos(0, 0, gen, rec, acc, hGenYSummary, hRecYSummary, hAccYSummary, 1);
    FillHistos(ipt, 0, gen, rec, acc, hGenSummary, hRecSummary, hAccSummary, ipt*fYBinLowEdge.GetSize()+1);
    
    // in 0-90% not weighted over centrality
    if (fSigWeights) IntegratedAccEff(ipt, 0, 0., 90., fNMatch, gen, rec, acc);
    
    // loop over signal weights
    Weights *weights0 = 0x0;
    nextSigWeight.Reset();
    while ((sigWeights = static_cast<Weights*>(nextSigWeight()))) {
      
      // check if these weights can be used for this pt/y bin
      if (!sigWeights->IsValid(ptMin, ptMax, fYBinLowEdge[0], fYBinLowEdge[fYBinLowEdge.GetSize()-1])) continue;
      
      // get the weights with the biggest pt/y validity range
      if (!weights0 || sigWeights->IsValid(*weights0)) weights0 = sigWeights;
      
      // print acc*eff weighted over centrality
      IntegratedAccEff(ipt, 0, *sigWeights, fNMatch, gen, rec, acc);
      
    }
    
    // acceptance*efficiency versus run
    DrawAccEffVsRun(ipt, 0, weights0, fNMatch);
    
    // acceptance*efficiency versus centrality
    if (fSigWeights) DrawAccEffVsCent(ipt, 0, fNMatch);
    
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

  // loop over y bins+1 (0 = integrated over y)
  for (Int_t iy = 1; iy < fYBinLowEdge.GetSize(); iy++) {
    
    Float_t yMin = (iy == 0) ? fYBinLowEdge[0] : fYBinLowEdge[iy-1];
    Float_t yMax = (iy == 0) ? fYBinLowEdge[fYBinLowEdge.GetSize()-1] : fYBinLowEdge[iy];
    
    // print integrated value
    printf("\n---- Integrated acc*eff (%g<pt<%g / %g<y<%g):\n", fPtBinLowEdge[0],fPtBinLowEdge[fPtBinLowEdge.GetSize()-1], yMin, yMax);
    
    // in 0-100% not weighted over centrality
    IntegratedAccEff(0, iy, 0., -1., fNMatch, gen, rec, acc);
    FillHistos(0, iy, gen, rec, acc, hGenYSummary, hRecYSummary, hAccYSummary, iy+1);
    FillHistos(0, iy, gen, rec, acc, hGenSummary, hRecSummary, hAccSummary, iy+1);
    
    // in 0-90% not weighted over centrality
    if (fSigWeights) IntegratedAccEff(0, iy, 0., 90., fNMatch, gen, rec, acc);
    
    // loop over signal weights
    Weights *weights0 = 0x0;
    nextSigWeight.Reset();
    while ((sigWeights = static_cast<Weights*>(nextSigWeight()))) {
      
      // check if these weights can be used for this pt/y bin
      if (!sigWeights->IsValid(fPtBinLowEdge[0],fPtBinLowEdge[fPtBinLowEdge.GetSize()-1], yMin, yMax)) continue;
      
      // get the weights with the biggest pt/y validity range
      if (!weights0 || sigWeights->IsValid(*weights0)) weights0 = sigWeights;
      
      // print acc*eff weighted over centrality
      IntegratedAccEff(0, iy, *sigWeights, fNMatch, gen, rec, acc);
      
    }
    
    // acceptance*efficiency versus run
    DrawAccEffVsRun(0, iy, weights0, fNMatch);
    
    // acceptance*efficiency versus centrality
    if (fSigWeights) DrawAccEffVsCent(0, iy, fNMatch);
    
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
  /*
  // loop over pt bins+1 (0 = integrated over pt)
  for (Int_t ipt = 1; ipt < fPtBinLowEdge.GetSize(); ipt++) {
    
    Float_t ptMin = (ipt == 0) ? fPtBinLowEdge[0] : fPtBinLowEdge[ipt-1];
    Float_t ptMax = (ipt == 0) ? fPtBinLowEdge[fPtBinLowEdge.GetSize()-1] : fPtBinLowEdge[ipt];
    
    // loop over y bins+1 (0 = integrated over y)
    for (Int_t iy = 1; iy < fYBinLowEdge.GetSize(); iy++) {
      
      Float_t yMin = (iy == 0) ? fYBinLowEdge[0] : fYBinLowEdge[iy-1];
      Float_t yMax = (iy == 0) ? fYBinLowEdge[fYBinLowEdge.GetSize()-1] : fYBinLowEdge[iy];
      
      // print integrated value
      printf("\n---- Integrated acc*eff (%g<pt<%g / %g<y<%g):\n", ptMin, ptMax, yMin, yMax);
      
      // in 0-100% not weighted over centrality
      IntegratedAccEff(ipt, iy, 0., -1., fNMatch, gen, rec, acc);
      FillHistos(ipt, iy, gen, rec, acc, hGenSummary, hRecSummary, hAccSummary, ipt*fYBinLowEdge.GetSize()+iy+1);
      
      // in 0-90% not weighted over centrality
      if (fSigWeights) IntegratedAccEff(ipt, iy, 0., 90., fNMatch, gen, rec, acc);
      
      // loop over signal weights
      Weights *weights0 = 0x0;
      nextSigWeight.Reset();
      while ((sigWeights = static_cast<Weights*>(nextSigWeight()))) {
	
	// check if these weights can be used for this pt/y bin
	if (!sigWeights->IsValid(ptMin, ptMax, yMin, yMax)) continue;
	
	// get the weights with the biggest pt/y validity range
	if (!weights0 || sigWeights->IsValid(*weights0)) weights0 = sigWeights;
	
	// print acc*eff weighted over centrality
	IntegratedAccEff(ipt, iy, *sigWeights, fNMatch, gen, rec, acc);
	
      }
      
      // acceptance*efficiency versus run
      DrawAccEffVsRun(ipt, iy, weights0, fNMatch);
      
      // acceptance*efficiency versus centrality
      if (fSigWeights) DrawAccEffVsCent(ipt, iy, fNMatch);
      
    }
    
  }
  */
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
  
  // print acceptance*efficiency integrated over pt/y versus a reduced number of centrality bins
  if (fSigWeights) {
    
    printf("\n---- acc*eff versus centrality (%g<pt<%g / %g<y<%g):\n",
	   fPtBinLowEdge[0],fPtBinLowEdge[fPtBinLowEdge.GetSize()-1],
	   fYBinLowEdge[0],fYBinLowEdge[fYBinLowEdge.GetSize()-1]);
    
    // loop over signal weights
    nextSigWeight.Reset();
    while ((sigWeights = static_cast<Weights*>(nextSigWeight()))) {
      
      // check if these weights can be used for the integrated pt/y bin
      if (!sigWeights->IsValid(fPtBinLowEdge[0], fPtBinLowEdge[fPtBinLowEdge.GetSize()-1],
			       fYBinLowEdge[0], fYBinLowEdge[fYBinLowEdge.GetSize()-1])) continue;
      
      // loop over centrality bins and print acc*eff
      const TArrayF &centBinLowEdge = sigWeights->GetCentBinLowEdge();
      for (Int_t iCent = 1; iCent < centBinLowEdge.GetSize(); iCent++)
	IntegratedAccEff(0, 0, centBinLowEdge[iCent-1], centBinLowEdge[iCent], fNMatch, gen, rec, acc);
      
    }
    
  }
  
}

//________________________________________________________________________
void AliAnalysisTaskJPsiAccEffCorr2::SetSigWeights(const Char_t *name, Float_t ptMin, Float_t ptMax,
						   Float_t yMin, Float_t yMax,
						   Int_t nCentBin, Float_t *binLowEdge,
						   Double_t *nSig, Bool_t weightRec)
{
  /// Set the number of signal versus centrality in a given pt/y bin
  /// (used to weight the acc*eff correction integrated over centrality for
  /// that pt/y bin or for any pt/y bins included in this one)
  
  if (ptMin > ptMax || yMin > yMax || nCentBin <= 0) {
    AliError("incorrect settings of signal weights");
    return;
  }
  
  if (!fSigWeights) {
    fSigWeights = new TList;
    fSigWeights->SetOwner();
  }
  
  fSigWeights->AddLast(new Weights(name, ptMin, ptMax, yMin, yMax, nCentBin, binLowEdge, nSig, weightRec));
  
}

//________________________________________________________________________
void AliAnalysisTaskJPsiAccEffCorr2::LoadRunWeights(const Char_t *fileName)
{
  /// Set the number of interested events per run
  /// (used to weight the acc*eff correction integrated
  /// over run for any pt/y/centrality bins)
  
  delete fRunWeights;
  fRunWeights = new THashList(1000);
  fRunWeights->SetOwner();
  
  ifstream inFile(fileName);
  if (!inFile.is_open()) {
    AliError(Form("cannot open file %s", fileName));
    exit(0);
  }
  
  TString line;
  while (! inFile.eof() ) {
    
    line.ReadLine(inFile,kTRUE);
    if(line.IsNull()) continue;
    
    TObjArray *param = line.Tokenize(" ");
    if (param->GetEntries() != 2) {
      AliError(Form("bad input line %s", line.Data()));
      continue;
    }
    
    Int_t run = ((TObjString*)param->UncheckedAt(0))->String().Atoi();
    if (run < 0) {
      AliError(Form("invalid run number: %d", run));
      continue;
    }
    
    Float_t weight = ((TObjString*)param->UncheckedAt(1))->String().Atof();
    if (weight <= 0.) {
      AliError(Form("invalid weight: %g", weight));
      continue;
    }
    
    fRunWeights->Add(new TParameter<Double_t>(Form("run:%d",run), weight));
    
    delete param;
  }
  
  inFile.close();
  
}

//________________________________________________________________________
void AliAnalysisTaskJPsiAccEffCorr2::CleanRunWeights()
{
  /// Remove the weights for the runs that have not been analyzed and warn when a weight is missing
  
  if (!fRunWeights) return;
  
  TObjArray *runList = fJPsiCounters->GetKeyWords("run").Tokenize(",");
  
  // loop over run weights
  TIter nextRunWeight(fRunWeights);
  TParameter<Double_t> *runWeight = 0x0;
  while ((runWeight = static_cast<TParameter<Double_t>*>(nextRunWeight()))) {
    
    // get the run number
    TString run = runWeight->GetName();
    run.Remove(0,4);
    
    // remove the weight if the run is not in the list
    if (!runList->FindObject(run.Data())) {
      AliError(Form("run %s has not been analyzed --> remove it from the list of run weights", run.Data()));
      fRunWeights->Remove(runWeight);
    }
    
  }
  
  // loop over analyzed runs
  TIter nextRun(runList);
  TObjString *run = 0x0;
  while ((run = static_cast<TObjString*>(nextRun()))) {
    
    // issue a warning if their is no weight for this run
    if (!fRunWeights->FindObject(Form("run:%s",run->GetName())))
      AliWarning(Form("no weight for run %s --> will be skipped in the run by run weighting", run->GetName()));
    
  }
  
  // clean memory
  delete runList;
  
}

//________________________________________________________________________
TH1F* AliAnalysisTaskJPsiAccEffCorr2::ComputeAccEff(TH1F &hGen, TH1F &hRec,
						   const Char_t *name, const Char_t *title)
{
  /// Compute acc*eff and binomial errors by hand, i.e. not using TGraphAsymmErrors
  /// Result is identical to divide histograms with option "B", except here error is forced > 1/gen
  
  Int_t nbins = hGen.GetNbinsX();
  TH1F* hAcc = new TH1F(name,title, nbins, hGen.GetXaxis()->GetXmin(), hGen.GetXaxis()->GetXmax());
  for (Int_t i = 1; i <= nbins; i++) {
    Double_t gen = hGen.GetBinContent(i);
    Double_t accEff = (gen>0.) ? hRec.GetBinContent(i)/gen : 0.;
    Double_t accEffErr = (gen>0.) ? TMath::Max(1./gen, TMath::Sqrt(accEff*TMath::Abs(1.-accEff)/gen)) : 1.;
    hAcc->SetBinContent(i, accEff);
    hAcc->SetBinError(i, accEffErr);
  }
  
  return hAcc;
}

//________________________________________________________________________
void AliAnalysisTaskJPsiAccEffCorr2::GetAccEff(TString selection, Double_t gen[3][2], Double_t rec[3][2], Double_t acc[3][3])
{
  /// return the number of generated and reconstructed JPsi and the acc*eff correction for the 3 matching requirements
  
  // get generated and reconstructed numbers of JPsi
  TH1D *hgen = fJPsiCounters->Get("match", Form("type:gen/%s", selection.Data()));
  if (hgen) hgen->SetName("hgen");
  TH1D *hrec = fJPsiCounters->Get("match", Form("type:rec/%s", selection.Data()));
  if (hrec) hrec->SetName("hrec");
  
  if (!hgen || !hrec) {
    
    // reset values
    for (Int_t i = 0; i < 3; i++) {
      gen[i][0] = rec[i][0] = acc[i][0] = 0.;
      gen[i][1] = rec[i][1] = acc[i][1] = acc[i][2] = 1.;
    }
    
    // clean memory
    delete hgen;
    delete hrec;
    
  } else {
    
    // compute acc*eff
    TGraphAsymmErrors *gacceff = new TGraphAsymmErrors(hrec, hgen, "cpe0");
    
    // get the results for each matching case
    for (Int_t i = 0; i < 3; i++) {
      
      gen[i][0] = hgen->GetBinContent(i+1);
      gen[i][1] = (gen[i][0] > 0.) ? hgen->GetBinError(i+1) : 1.;
      
      rec[i][0] = hrec->GetBinContent(i+1);
      rec[i][1] = (rec[i][0] > 0.) ? hrec->GetBinError(i+1) : 1.;
      
      // switch to "manual" mode in pathological cases
      if (gen[i][0] <= 0.) {
	
	acc[i][0] = 0.;
	acc[i][1] = acc[i][2] = 1.;
	
      } else if (gen[i][0] < rec[i][0]) {
	
	acc[i][0] = rec[i][0]/gen[i][0];
	acc[i][1] = acc[i][2] = TMath::Max(1./gen[i][0], TMath::Sqrt(acc[i][0]*TMath::Abs(1.-acc[i][0])/gen[i][0]));
	
      } else {
	
	Double_t x = 0.;
	gacceff->GetPoint(i,x,acc[i][0]);
	acc[i][1] = gacceff->GetErrorYlow(i);
	acc[i][2] = gacceff->GetErrorYhigh(i);
	
      }
      
    }
    
    // clean memory
    delete hgen;
    delete hrec;
    delete gacceff;
    
  }
  
}

//________________________________________________________________________
void AliAnalysisTaskJPsiAccEffCorr2::IntegratedAccEff(TString commonKey1, const THashList &weights1, Bool_t weightRec1,
						      Double_t gen[3][2], Double_t rec[3][2], Double_t acc[3][3],
						      const THashList *weights2, Bool_t weightRec2)
{
  /// return the number of generated and reconstructed JPsi and the acc*eff correction for the 3 matching requirements
  
  for (Int_t i = 0; i < 3; i++) {
    for (Int_t j = 0; j < 2; j++) {
      gen[i][j] = 0.;
      rec[i][j] = 0.;
      acc[i][j] = 0.;
    }
    acc[i][2] = 0.;
  }
  
  Double_t geni[3][2], reci[3][2], acci[3][3];
  Double_t tmpErr[3][2] = {{0.,0.},{0.,0.},{0.,0.}};
  
  // loop over bins with weights
  TIter nextw(&weights1);
  TParameter<Double_t> *w = 0x0;
  while ((w = static_cast<TParameter<Double_t>*>(nextw()))) {
    
    // set the common selection key for the second level
    TString commonKey2 = Form("%s/%s", commonKey1.Data(), w->GetName());
    
    // compute acc*eff for the current bin
    if (weights2) IntegratedAccEff(commonKey2, *weights2, weightRec2, geni, reci, acci);
    else GetAccEff(commonKey2, geni, reci, acci);
    
    // integrate
    Double_t w2 = w->GetVal() * w->GetVal();
    for (Int_t iTrig = 0; iTrig < 3; iTrig++) {
      
      if (weightRec1) {
	
	// cannot use a bin where acc*eff = 0 --> skip it
	if (acci[iTrig][0] <= 0.) {
	  AliError(Form("acc*eff = 0 for %d matching in bin %s --> bin skipped", iTrig, w->GetName()));
	  continue;
	}
	
	// weight the number of reconstructed JPsi
	rec[iTrig][0] += w->GetVal();
	rec[iTrig][1] += w2*reci[iTrig][1]*reci[iTrig][1]/reci[iTrig][0]/reci[iTrig][0];
	gen[iTrig][0] += w->GetVal()/acci[iTrig][0];
	gen[iTrig][1] += w2*geni[iTrig][1]*geni[iTrig][1]/reci[iTrig][0]/reci[iTrig][0];
	tmpErr[iTrig][0] += w2*acci[iTrig][1]*acci[iTrig][1]/acci[iTrig][0]/acci[iTrig][0]/acci[iTrig][0]/acci[iTrig][0];
	tmpErr[iTrig][1] += w2*acci[iTrig][2]*acci[iTrig][2]/acci[iTrig][0]/acci[iTrig][0]/acci[iTrig][0]/acci[iTrig][0];
	
      } else {
	
	// cannot use a bin where gen = 0 --> skip it
	if (geni[iTrig][0] <= 0.) {
	  if (iTrig == 0) AliError(Form("gen = 0 in bin %s --> bin skipped", w->GetName()));
	  continue;
	}
	
	// weight the number of generated JPsi
	rec[iTrig][0] += w->GetVal()*acci[iTrig][0];
	rec[iTrig][1] += w2*reci[iTrig][1]*reci[iTrig][1]/geni[iTrig][0]/geni[iTrig][0];
	gen[iTrig][0] += w->GetVal();
	gen[iTrig][1] += w2*geni[iTrig][1]*geni[iTrig][1]/geni[iTrig][0]/geni[iTrig][0];
	tmpErr[iTrig][0] += w2*acci[iTrig][1]*acci[iTrig][1];
	tmpErr[iTrig][1] += w2*acci[iTrig][2]*acci[iTrig][2];
	
      }
      
    }
    
  }
  
  // complete integration
  for (Int_t iTrig = 0; iTrig < 3; iTrig++) {
    gen[iTrig][1] = (gen[iTrig][1] > 0.) ? TMath::Sqrt(gen[iTrig][1]) : 1.;
    rec[iTrig][1] = (rec[iTrig][1] > 0.) ? TMath::Sqrt(rec[iTrig][1]) : 1.;
    acc[iTrig][0] = (gen[iTrig][0] > 0.) ? rec[iTrig][0]/gen[iTrig][0] : 0.;
    Double_t x = weightRec1 ? acc[iTrig][0] : 1.;
    acc[iTrig][1] = (gen[iTrig][0] > 0.) ? x*TMath::Sqrt(tmpErr[iTrig][0])/gen[iTrig][0] : 1.;
    acc[iTrig][2] = (gen[iTrig][0] > 0.) ? x*TMath::Sqrt(tmpErr[iTrig][1])/gen[iTrig][0] : 1.;
  }
  
}

//________________________________________________________________________
void AliAnalysisTaskJPsiAccEffCorr2::IntegratedAccEff(Int_t ipt, Int_t iy, Float_t centMin, Float_t centMax,
						      Int_t nMatch, Double_t gen[2], Double_t rec[2], Double_t acc[3],
						      Bool_t print)
{
  /// compute AccEff correction integrated over the given centrality range
  
  TString ptKey = (ipt == 0) ? "any" : Form("%g-%g",fPtBinLowEdge[ipt-1],fPtBinLowEdge[ipt]);
  TString yKey = (iy == 0) ? "any" : Form("%g-%g",fYBinLowEdge[iy-1],fYBinLowEdge[iy]);
  TString centKey = Weights::GetCentKey(centMin, centMax, fCentBinLowEdge);
  TString selection = Form("ptbin:%s/ybin:%s/cent:%s",ptKey.Data(),yKey.Data(), centKey.Data());
  Int_t iTrig = (nMatch >= 0 && nMatch <= 2) ? nMatch : 2;
  
  // compute acc*eff
  Double_t g[3][2], r[3][2], a[3][3];
  if (fRunWeights) IntegratedAccEff(selection, *fRunWeights, kFALSE, g, r, a);
  else GetAccEff(selection, g, r, a);
  for (Int_t i = 0; i < 2; i++) {
    gen[i] = g[iTrig][i];
    rec[i] = r[iTrig][i];
    acc[i] = a[iTrig][i];
  }
  acc[2] = a[iTrig][2];
  
  // print acc*eff for different matching conditions if required
  if (print) {
    TString cent = (centMax < centMin) ? "0-100%" : Form("%g-%g%%", centMin, centMax);
    if (nMatch < 0 || nMatch == 0) printf("   - 0 matching required - %s: %f + %f - %f\n",
					  cent.Data(), a[0][0], a[0][2], a[0][1]);
    if (nMatch < 0 || nMatch == 1) printf("   - 1 matching required - %s: %f + %f - %f\n",
					  cent.Data(), a[1][0], a[1][2], a[1][1]);
    if (nMatch < 0 || nMatch == 2) printf("   - 2 matching required - %s: %f + %f - %f\n",
					  cent.Data(), a[2][0], a[2][2], a[2][1]);
  }
  
}

//________________________________________________________________________
void AliAnalysisTaskJPsiAccEffCorr2::IntegratedAccEff(Int_t ipt, Int_t iy, Weights &w, Int_t nMatch,
						      Double_t gen[2], Double_t rec[2], Double_t acc[3],
						      Bool_t print)
{
  /// compute AccEff correction weighted over centrality
  
  TString ptKey = (ipt == 0) ? "any" : Form("%g-%g",fPtBinLowEdge[ipt-1],fPtBinLowEdge[ipt]);
  TString yKey = (iy == 0) ? "any" : Form("%g-%g",fYBinLowEdge[iy-1],fYBinLowEdge[iy]);
  TString selection = Form("ptbin:%s/ybin:%s",ptKey.Data(),yKey.Data());
  Int_t iTrig = (nMatch >= 0 && nMatch <= 2) ? nMatch : 2;
  
  // compute acc*eff
  Double_t g[3][2], r[3][2], a[3][3];
  if (fRunWeights) IntegratedAccEff(selection, w.GetWeights(), w.WeightRec(), g, r, a, fRunWeights, kFALSE);
  else IntegratedAccEff(selection, w.GetWeights(), w.WeightRec(), g, r, a);
  for (Int_t i = 0; i < 2; i++) {
    gen[i] = g[iTrig][i];
    rec[i] = r[iTrig][i];
    acc[i] = a[iTrig][i];
  }
  acc[2] = a[iTrig][2];
  
  // print acc*eff for different matching conditions if required
  if (print) {
    Float_t centMin, centMax;
    w.GetCentralityRange(centMin, centMax);
    if (nMatch < 0 || nMatch == 0) printf("   - 0 matching required - %g-%g%% - weighted per %s in %s: %f + %f - %f\n",
					  centMin, centMax, w.GetName(), w.GetTitle(), a[0][0], a[0][2], a[0][1]);
    if (nMatch < 0 || nMatch == 1) printf("   - 1 matching required - %g-%g%% - weighted per %s in %s: %f + %f - %f\n",
					  centMin, centMax, w.GetName(), w.GetTitle(), a[1][0], a[1][2], a[1][1]);
    if (nMatch < 0 || nMatch == 2) printf("   - 2 matching required - %g-%g%% - weighted per %s in %s: %f + %f - %f\n",
					  centMin, centMax, w.GetName(), w.GetTitle(), a[2][0], a[2][2], a[2][1]);
  }
  
}

//________________________________________________________________________
void AliAnalysisTaskJPsiAccEffCorr2::FillHistos(Int_t ipt, Int_t iy, Double_t gen[2], Double_t rec[2], Double_t acc[3],
						TH1F* hGenSum, TH1F* hRecSum, TH1F* hAccSum, Int_t bin)
{
  /// fill summary plots
  
  TString label = "";
  label += (ipt == 0) ? Form("%g<pt<%g",fPtBinLowEdge[0],fPtBinLowEdge[fPtBinLowEdge.GetSize()-1]) : Form("%g<pt<%g",fPtBinLowEdge[ipt-1],fPtBinLowEdge[ipt]);
  label += (iy == 0) ? Form(" / %g<y<%g",fYBinLowEdge[0],fYBinLowEdge[fYBinLowEdge.GetSize()-1]) : Form(" / %g<y<%g",fYBinLowEdge[iy-1],fYBinLowEdge[iy]);
  
  hGenSum->SetBinContent(bin, gen[0]);
  hGenSum->SetBinError(bin, gen[1]);
  hGenSum->GetXaxis()->SetBinLabel(bin, label.Data());
  
  hRecSum->SetBinContent(bin, rec[0]);
  hRecSum->SetBinError(bin, rec[1]);
  hRecSum->GetXaxis()->SetBinLabel(bin, label.Data());
  
  hAccSum->SetBinContent(bin, acc[0]);
  hAccSum->SetBinError(bin, TMath::Max(acc[1],acc[2]));
  hAccSum->GetXaxis()->SetBinLabel(bin, label.Data());
  
}

//________________________________________________________________________
void AliAnalysisTaskJPsiAccEffCorr2::DrawAccEffVsRun(Int_t ipt, Int_t iy, Weights *w, Int_t nMatch)
{
  /// Draw acceptance*efficiency versus run for this given pt/y bin
  
  TString ptKey = (ipt == 0) ? "any" : Form("%g-%g",fPtBinLowEdge[ipt-1],fPtBinLowEdge[ipt]);
  TString yKey = (iy == 0) ? "any" : Form("%g-%g",fYBinLowEdge[iy-1],fYBinLowEdge[iy]);
  TString commonKey = Form("ptbin:%s/ybin:%s",ptKey.Data(),yKey.Data());
  Int_t iTrig = (nMatch >= 0 && nMatch <= 2) ? nMatch : 2;
  
  // get the list of runs
  TObjArray *runs = fJPsiCounters->GetKeyWords("run").Tokenize(",");
  Int_t nRuns = runs->GetEntriesFast();
  
  // create histos
  TString range = (ipt == 0) ? Form("(%g<pt<%g",fPtBinLowEdge[0],fPtBinLowEdge[fPtBinLowEdge.GetSize()-1]) : Form("(%g<pt<%g",fPtBinLowEdge[ipt-1],fPtBinLowEdge[ipt]);
  range += (iy == 0) ? Form(" / %g<y<%g)",fYBinLowEdge[0],fYBinLowEdge[fYBinLowEdge.GetSize()-1]) : Form(" / %g<y<%g)",fYBinLowEdge[iy-1],fYBinLowEdge[iy]);
  TString unit = w ? "a.u." : "N_{J/#psi}";
  TH1F* hgen = new TH1F(Form("hgenVsRun_pt%d_y%d",ipt,iy),Form("generated versus run %s;run;%s",range.Data(), unit.Data()), nRuns, -0.5, nRuns-0.5);
  TH1F* hrec = new TH1F(Form("hrecVsRun_pt%d_y%d",ipt,iy),Form("reconstructed versus run %s;run;%s",range.Data(), unit.Data()), nRuns, -0.5, nRuns-0.5);
  TGraphAsymmErrors *gacceff = new TGraphAsymmErrors(nRuns);
  gacceff->SetNameTitle(Form("accEffVsRun_pt%d_y%d",ipt,iy), Form("acceptance * efficiency versus run %s",range.Data()));
  
  // loop over runs
  Int_t irun = 1;
  TObjString * run = 0x0;
  TIter nextRun(runs);
  while ((run = static_cast<TObjString*>(nextRun()))) {
    
    // set the selection key for this run
    TString selection = Form("%s/run:%s", commonKey.Data(), run->GetName());
    
    // compute acc*eff
    Double_t g[3][2], r[3][2], a[3][3];
    if (w) IntegratedAccEff(selection, w->GetWeights(), w->WeightRec(), g, r, a);
    else GetAccEff(selection, g, r, a);
    
    // fill histos
    hgen->SetBinContent(irun, g[iTrig][0]);
    hgen->SetBinError(irun, g[iTrig][1]);
    hgen->GetXaxis()->SetBinLabel(irun, run->GetName());
    hrec->SetBinContent(irun, r[iTrig][0]);
    hrec->SetBinError(irun, r[iTrig][1]);
    hrec->GetXaxis()->SetBinLabel(irun, run->GetName());
    gacceff->SetPoint(irun-1,irun-1,a[iTrig][0]);
    gacceff->SetPointError(irun-1,0.,0.,a[iTrig][1],a[iTrig][2]);
    irun++;
    
  }
  
  // Set graph labels
  irun = 1;
  nextRun.Reset();
  gacceff->GetXaxis()->Set(nRuns, -0.5, nRuns-0.5);
  gacceff->GetXaxis()->SetNdivisions(1,kFALSE);
  while ((run = static_cast<TObjString*>(nextRun()))) 
    gacceff->GetXaxis()->SetBinLabel(irun++, run->GetName());
  
  // draw histos
  new TCanvas(Form("cAccEffVsRun_pt%d_y%d",ipt,iy), Form("acceptance * efficiency versus run %s",range.Data()), 1000, 300);
  gacceff->SetLineStyle(1);
  gacceff->SetLineColor(1);
  gacceff->SetMarkerStyle(20);
  gacceff->SetMarkerSize(0.7);
  gacceff->SetMarkerColor(2);
  gacceff->GetXaxis()->SetTitle("run");
  gacceff->GetXaxis()->SetLabelFont(22);
  gacceff->GetXaxis()->SetTitleFont(22);
  gacceff->GetYaxis()->SetTitle("Acc*Eff");
  gacceff->GetYaxis()->SetLabelFont(22);
  gacceff->GetYaxis()->SetLabelFont(22);
  gacceff->Draw("ap");
  
  // save histos
  TFile *file = new TFile("acceff_new.root","update");
  hgen->Write(0x0, TObject::kOverwrite);
  hrec->Write(0x0, TObject::kOverwrite);
  gacceff->Write(0x0, TObject::kOverwrite);
  file->Close();
  
  // clean memory
  delete hgen;
  delete hrec;
  delete runs;
  
}

//________________________________________________________________________
void AliAnalysisTaskJPsiAccEffCorr2::DrawAccEffVsCent(Int_t ipt, Int_t iy, Int_t nMatch)
{
  /// Draw acceptance*efficiency versus centrality for this given pt/y bin
  
  // create histos
  Int_t nCents = fCentBinLowEdge.GetSize();
  TString range = (ipt == 0) ? Form("(%g<pt<%g",fPtBinLowEdge[0],fPtBinLowEdge[fPtBinLowEdge.GetSize()-1]) : Form("(%g<pt<%g",fPtBinLowEdge[ipt-1],fPtBinLowEdge[ipt]);
  range += (iy == 0) ? Form(" / %g<y<%g)",fYBinLowEdge[0],fYBinLowEdge[fYBinLowEdge.GetSize()-1]) : Form(" / %g<y<%g)",fYBinLowEdge[iy-1],fYBinLowEdge[iy]);
  TString unit = fRunWeights ? "a.u." : "N_{J/#psi}";
  TH1F* hgen = new TH1F(Form("hgenVsCent_pt%d_y%d",ipt,iy),Form("generated versus centrality %s;centrality;%s",range.Data(),unit.Data()), nCents-1, -0.5, nCents-1.5);
  TH1F* hrec = new TH1F(Form("hrecVsCent_pt%d_y%d",ipt,iy),Form("reconstructed versus centrality %s;centrality;%s",range.Data(),unit.Data()), nCents-1, -0.5, nCents-1.5);
  TGraphAsymmErrors *gacceff = new TGraphAsymmErrors(nCents-1);
  gacceff->SetNameTitle(Form("accEffVsCent_pt%d_y%d",ipt,iy), Form("acceptance * efficiency versus centrality %s",range.Data()));
  
  // loop over centrality
  for (Int_t iCent = 1; iCent < nCents; iCent++) {
    
    // compute acc*eff
    Double_t g[2], r[2], a[3];
    IntegratedAccEff(ipt, iy, fCentBinLowEdge[iCent-1], fCentBinLowEdge[iCent], nMatch, g, r, a, kFALSE);
    
    // fill histos
    TString label = Form("%g-%g",fCentBinLowEdge[iCent-1], fCentBinLowEdge[iCent]);
    hgen->SetBinContent(iCent, g[0]);
    hgen->SetBinError(iCent, g[1]);
    hgen->GetXaxis()->SetBinLabel(iCent, label.Data());
    hrec->SetBinContent(iCent, r[0]);
    hrec->SetBinError(iCent, r[1]);
    hrec->GetXaxis()->SetBinLabel(iCent, label.Data());
    gacceff->SetPoint(iCent-1,iCent-1,a[0]);
    gacceff->SetPointError(iCent-1,0.,0.,a[1],a[2]);
    
  }
  
  // Set graph labels
  gacceff->GetXaxis()->Set(nCents-1, -0.5, nCents-1.5);
  gacceff->GetXaxis()->SetNdivisions(1,kFALSE);
  for (Int_t iCent = 1; iCent < nCents; iCent++)
    gacceff->GetXaxis()->SetBinLabel(iCent, Form("%g-%g",fCentBinLowEdge[iCent-1], fCentBinLowEdge[iCent]));
  
  // draw histos
  new TCanvas(Form("cAccEffVsCent_pt%d_y%d",ipt,iy), Form("acceptance * efficiency versus centrality %s",range.Data()), 1000, 300);
  gacceff->SetLineStyle(1);
  gacceff->SetLineColor(1);
  gacceff->SetMarkerStyle(20);
  gacceff->SetMarkerSize(0.7);
  gacceff->SetMarkerColor(2);
  gacceff->GetXaxis()->SetTitle("centrality");
  gacceff->GetXaxis()->SetLabelFont(22);
  gacceff->GetXaxis()->SetTitleFont(22);
  gacceff->GetYaxis()->SetTitle("Acc*Eff");
  gacceff->GetYaxis()->SetLabelFont(22);
  gacceff->GetYaxis()->SetLabelFont(22);
  gacceff->Draw("ap");
  
  // save histos
  TFile *file = new TFile("acceff_new.root","update");
  hgen->Write(0x0, TObject::kOverwrite);
  hrec->Write(0x0, TObject::kOverwrite);
  gacceff->Write(0x0, TObject::kOverwrite);
  file->Close();
  
  // clean memory
  delete hgen;
  delete hrec;
  
}

