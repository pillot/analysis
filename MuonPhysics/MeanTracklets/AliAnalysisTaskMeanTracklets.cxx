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

// ROOT includes
#include "THnSparse.h"
#include "TMath.h"
#include "TRandom3.h"

// STEER includes
#include "AliVEvent.h"
#include "AliAODEvent.h"
#include "AliVVertex.h"
#include "AliAODTracklets.h"
#include "AliCounterCollection.h"
#include "AliMCEvent.h"
#include "AliGenPythiaEventHeader.h"
#include "AliAODMCParticle.h"

// ANALYSIS includes
#include "AliAnalysisTaskMeanTracklets.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataSlot.h"
#include "AliAnalysisDataContainer.h"
#include "AliInputEventHandler.h"
#include "AliMultSelection.h"
#include "AliAnalysisMuonUtility.h"

ClassImp(AliAnalysisTaskMeanTracklets)

//________________________________________________________________________
AliAnalysisTaskMeanTracklets::AliAnalysisTaskMeanTracklets() :
AliAnalysisTaskSE(),
fEvents(0x0),
fhNtrk(0x0),
fhNtrkCorr(0x0),
fhNtrkCorrVsCuts(0x0),
fpMeanNtrkVsZvtx(0x0),
fpMeanNtrkVsZvtxCorr(0x0),
fpMeanNtrkVsZvtxRef(0x0),
fMeanNtrkRef(-1.),
fUseBinomial(kFALSE),
fRandom(new TRandom3(0)),
fTrigger(""),
fPSTriggerMask(0),
fRejectSD(kFALSE),
fRejectPUFromSPD(kFALSE),
fSelectSPDVtxQA(kTRUE),
fReject0Tracklet(kFALSE)
{
  // Dummy constructor
}

//________________________________________________________________________
AliAnalysisTaskMeanTracklets::AliAnalysisTaskMeanTracklets(const char *name) :
AliAnalysisTaskSE(name),
fEvents(0x0),
fhNtrk(0x0),
fhNtrkCorr(0x0),
fhNtrkCorrVsCuts(0x0),
fpMeanNtrkVsZvtx(0x0),
fpMeanNtrkVsZvtxCorr(0x0),
fpMeanNtrkVsZvtxRef(0x0),
fMeanNtrkRef(-1.),
fUseBinomial(kFALSE),
fRandom(new TRandom3(0)),
fTrigger(""),
fPSTriggerMask(0),
fRejectSD(kFALSE),
fRejectPUFromSPD(kFALSE),
fSelectSPDVtxQA(kTRUE),
fReject0Tracklet(kFALSE)
{
  /// Constructor
  
  // Output slot #1 writes into a AliCounterCollection
  DefineOutput(1,AliCounterCollection::Class());
  // Output slot #2 writes into a THnSparse
  DefineOutput(2,THnSparse::Class());
  // Output slot #3 writes into a THnSparse
  DefineOutput(3,THnSparse::Class());
  // Output slot #4 writes into a THnSparse
  DefineOutput(4,THnSparse::Class());
  // Output slot #5 writes into a TProfile
  DefineOutput(5,TProfile::Class());
  // Output slot #6 writes into a TProfile
  DefineOutput(6,TProfile::Class());
  
}

//________________________________________________________________________
AliAnalysisTaskMeanTracklets::~AliAnalysisTaskMeanTracklets()
{
  /// Destructor
  if (!AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
    delete fEvents;
    delete fhNtrk;
    delete fhNtrkCorr;
    delete fhNtrkCorrVsCuts;
    delete fpMeanNtrkVsZvtx;
    delete fpMeanNtrkVsZvtxCorr;
  }
  delete fpMeanNtrkVsZvtxRef;
  delete fRandom;
}

//___________________________________________________________________________
void AliAnalysisTaskMeanTracklets::UserCreateOutputObjects()
{
  /// Create output objects
  
  printf("\nseed = %u\n\n", fRandom->GetSeed());
  
  // events analyzed
  fEvents = new AliCounterCollection(GetOutputSlot(1)->GetContainer()->GetName());
  fEvents->AddRubric("event", "any/selected");
  fEvents->AddRubric("run", 10000);
  fEvents->Init();
  
  // prepare binning for THnSparse for N SPD tracklets studies
  // 0: N SPD tracklets
  // 1: zVtx
  // 2: run
  // 3: N MC particles
  const Int_t nDims = 4;
  Int_t nBins[nDims] = {350, 80, 300000, 350};
  Double_t xMin[nDims] = {-0.5, -10., 99999.5, -0.5};
  Double_t xMax[nDims] = {349.5, 10., 399999.5, 349.5};

  // prepare binning for THnSparse for cut efficiency studies
  // 0: Corrected (if correction profile is provided) N SPD tracklets
  // 1: N MC particles
  // 2: run
  // 3: flag: physics selection (with or without pile-up rejection)
  // 4: flag: SPD pile-up
  // 5: flag: Vtx QA
  // 6: flag: reconstructed |zVtx| < 10cm
  // 7: flag: MC |zVtx| < 10cm
  // 8: flag: NSD
  // 9: flag: INEL>0
  const Int_t nDims2 = 10;
  Int_t nBins2[nDims2] = {350, 350, 300000, 2, 2, 2, 2, 2, 2, 2};
  Double_t xMin2[nDims2] = {-0.5, -0.5, 99999.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5};
  Double_t xMax2[nDims2] = {349.5, 349.5, 399999.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5};

  // create histogram
  fhNtrk = new THnSparseT<TArrayF>("hNtrk", "N SPD tracklets", nDims, nBins, xMin, xMax);
  fhNtrk->Sumw2();
  fhNtrkCorr = new THnSparseT<TArrayF>("hNtrkCorr", "Corrected N SPD tracklets", nDims, nBins, xMin, xMax);
  fhNtrkCorr->Sumw2();
  fhNtrkCorrVsCuts = new THnSparseT<TArrayF>("hNtrkCorrVsCuts", "Corrected N SPD tracklets versus cuts", nDims2, nBins2, xMin2, xMax2);
  fhNtrkCorrVsCuts->Sumw2();

  fpMeanNtrkVsZvtx = new TProfile("fpMeanNtrkVsZvtx", "<Ntrk> vs Zvtx", nBins[1], xMin[1], xMax[1]);
  fpMeanNtrkVsZvtxCorr = new TProfile("fpMeanNtrkVsZvtxCorr", "corrected <Ntrk> vs Zvtx", nBins[1], xMin[1], xMax[1]);
  
  // Post data at least once per task to ensure data synchronisation (required for merging)
  PostData(1, fEvents);
  PostData(2, fhNtrk);
  PostData(3, fhNtrkCorr);
  PostData(4, fhNtrkCorrVsCuts);
  PostData(5, fpMeanNtrkVsZvtx);
  PostData(6, fpMeanNtrkVsZvtxCorr);
  
}

//________________________________________________________________________
void AliAnalysisTaskMeanTracklets::UserExec(Option_t *)
{
  /// Called for each event
  
  // get the input event
  AliVEvent *evt = InputEvent();
  
  // cast it to AOD and make sure it is actually AOD
  if (!dynamic_cast<AliAODEvent*>(evt)) return;
  
  // fill event counters for all events
  fEvents->Count(Form("event:any/run:%d", fCurrentRunNumber));
  
  // select NSD MC events (only Pythia case supported)
  Int_t nsdFlag = 0;
  AliMCEvent *mcEvt = MCEvent();
  if (mcEvt) {
    AliGenEventHeader* mcGenH = fMCEvent->GenEventHeader();
    if (mcGenH && mcGenH->InheritsFrom(AliGenPythiaEventHeader::Class())) {
      AliGenPythiaEventHeader *hPythia = static_cast<AliGenPythiaEventHeader*>(mcGenH);
      if (hPythia->ProcessType() != 92 && hPythia->ProcessType() != 93) nsdFlag = 1; // non single-diffractive
    }
  }

  // select a specific trigger (select all by default)
  Int_t triggerFlag = 1;
  if (!fTrigger.IsNull()) {
    TString trigger = evt->GetFiredTriggerClasses();
    if (mcEvt && fTrigger == "MB" && !(trigger.Contains("V0L") && trigger.Contains("V0R"))) triggerFlag = 0;
    else if (!trigger.Contains(fTrigger.Data())) triggerFlag = 0;
  }

  // apply physics selection
  Int_t psFlag = (fInputHandler->IsEventSelected() & fPSTriggerMask) ? 1 : 0;

  // reject events with pile-up from SPD
  Int_t spdPUFlag = evt->IsPileupFromSPDInMultBins() ? 1 : 0;

  // apply vertex QA selection
  Int_t vtxQAFlag = 0;
  const AliVVertex* vtx = evt->GetPrimaryVertexSPD();
  if (vtx && vtx->GetNContributors() > 0) {
    Double_t cov[6]={0};
    vtx->GetCovarianceMatrix(cov);
    if (TMath::Sqrt(cov[5]) <= 0.25) vtxQAFlag = 1;
  }

  // get and select on vertex position
  Double_t zVtx = (vtxQAFlag == 1) ? vtx->GetZ() : 0.;
  Int_t zVtxRangeFlag = (zVtx > -10. && zVtx < 10.) ? 1 : 0;

  // get and select on MC vertex position
  Double_t zVtxMC = mcEvt ? AliAnalysisMuonUtility::GetMCVertexZ(evt,mcEvt) : 0.;
  Int_t zVtxMCRangeFlag = (zVtxMC > -10. && zVtxMC < 10.) ? 1 : 0;

  // get number of charged particules
  Int_t nMCpartInEtaRange(0);
  if (mcEvt) {
    Int_t nMCpart = mcEvt->GetNumberOfTracks();
    for (Int_t i = 0; i < nMCpart ; ++i) {
      AliAODMCParticle *p = static_cast<AliAODMCParticle*>(mcEvt->GetTrack(i));
      if (!p->IsPhysicalPrimary()) continue;
      if (p->Charge() == 0) continue;
      if (TMath::Abs(p->Eta()) > 1.) continue;
      ++nMCpartInEtaRange;
    }
  }

  // select INEL>0 MC events
  Int_t inelPosFlag = (nMCpartInEtaRange > 0) ? 1 : 0;

  // get the number of SPD tracklets (WARNING: not reliable without a good SPD vertex)
  AliAODTracklets *tracklets = static_cast<const AliAODEvent*>(evt)->GetTracklets();
  Int_t nTrk = tracklets->GetNumberOfTracklets();
  Int_t nTrkInEtaRange(0);
  for (Int_t i = 0; i < nTrk; i++) {
    Double_t eta = -TMath::Log(TMath::Tan(tracklets->GetTheta(i)/2.));
    if ( eta < -1. || eta > 1. ) continue;
    ++nTrkInEtaRange;
  }
  
  // get the corrected number of SPD tracklets (WARNING: not reliable without a good SPD vertex)
  Int_t nTrkInEtaRangeCorr = fpMeanNtrkVsZvtxRef ? GetCorrectedNtrk(nTrkInEtaRange,zVtx) : GetCorrectedNtrkFromMultSel();

  // fill histograms for cut efficiency studies
  Double_t EventInfoVsCuts[10];
  EventInfoVsCuts[0] = fpMeanNtrkVsZvtxRef ? nTrkInEtaRangeCorr : nTrkInEtaRange;
  EventInfoVsCuts[1] = nMCpartInEtaRange;
  EventInfoVsCuts[2] = fCurrentRunNumber;
  EventInfoVsCuts[3] = psFlag;
  EventInfoVsCuts[4] = spdPUFlag;
  EventInfoVsCuts[5] = vtxQAFlag;
  EventInfoVsCuts[6] = zVtxRangeFlag;
  EventInfoVsCuts[7] = zVtxMCRangeFlag;
  EventInfoVsCuts[8] = nsdFlag;
  EventInfoVsCuts[9] = inelPosFlag;
  fhNtrkCorrVsCuts->Fill(EventInfoVsCuts);

  // apply all events selections
  if (fRejectSD && nsdFlag == 0) return;
  if (triggerFlag == 0) return;
  if (fPSTriggerMask != 0 && psFlag == 0) return;
  if (fRejectPUFromSPD && spdPUFlag == 1) return;
  if (fSelectSPDVtxQA && vtxQAFlag == 0) return;
  if (zVtxRangeFlag == 0) return;
  if (fReject0Tracklet && nTrkInEtaRange == 0) return;
  
  // fill event counters for selected events
  fEvents->Count(Form("event:selected/run:%d", fCurrentRunNumber));

  // fill histograms for N SPD tracklets studies
  Double_t EventInfo[4];
  EventInfo[0] = nTrkInEtaRange;
  EventInfo[1] = zVtx;
  EventInfo[2] = fCurrentRunNumber;
  EventInfo[3] = nMCpartInEtaRange;
  fhNtrk->Fill(EventInfo);
  fpMeanNtrkVsZvtx->Fill(zVtx,nTrkInEtaRange);
  if (nTrkInEtaRangeCorr > 0 || (nTrkInEtaRangeCorr == 0 && !fReject0Tracklet)) {
    EventInfo[0] = nTrkInEtaRangeCorr;
    fhNtrkCorr->Fill(EventInfo);
    fpMeanNtrkVsZvtxCorr->Fill(zVtx,nTrkInEtaRangeCorr);
  }

  // Post data
  PostData(1, fEvents);
  PostData(2, fhNtrk);
  PostData(3, fhNtrkCorr);
  PostData(4, fhNtrkCorrVsCuts);
  PostData(5, fpMeanNtrkVsZvtx);
  PostData(6, fpMeanNtrkVsZvtxCorr);

}

//________________________________________________________________________
Int_t AliAnalysisTaskMeanTracklets::GetCorrectedNtrk(Int_t nTrkInEtaRange, Double_t zVtx)
{
  /// return the number of SPD tracklet corrected for the Zvtx dependence of <Ntrk>.
  /// return -999 in case of invalid correction
  
  Double_t meanNtrk = fpMeanNtrkVsZvtxRef->GetBinContent(fpMeanNtrkVsZvtxRef->FindBin(zVtx));
  if (meanNtrk < 1.e-6) return -999;
  
  if (fUseBinomial && fMeanNtrkRef <= meanNtrk) {
    
    return fRandom->Binomial(nTrkInEtaRange, fMeanNtrkRef/meanNtrk);
    
  } else {
    
    Double_t dN = nTrkInEtaRange*fMeanNtrkRef/meanNtrk - nTrkInEtaRange;
    Int_t sign = (dN > 0.) ? 1 : -1;
    /*
     Int_t nTrkInEtaRangeCorr = -1;
     do {
     nTrkInEtaRangeCorr = nTrkInEtaRange + sign*fRandom->Poisson(TMath::Abs(dN));
     } while (nTrkInEtaRangeCorr < 0);
     
     return nTrkInEtaRangeCorr;
     */
    return TMath::Max(nTrkInEtaRange + sign*fRandom->Poisson(TMath::Abs(dN)), 0);
    
  }
  
}

//________________________________________________________________________
Int_t AliAnalysisTaskMeanTracklets::GetCorrectedNtrkFromMultSel()
{
  /// return the corrected number of SPD tracklet from the multiplicity framework
  /// return -999 if not set
  
  AliMultSelection *multSelection = static_cast<AliMultSelection*>(InputEvent()->FindListObject("MultSelection"));
  if (!multSelection) return -999;
  
  return static_cast<Int_t>(multSelection->GetEstimator("SPDTracklets")->GetValue());
  
}

