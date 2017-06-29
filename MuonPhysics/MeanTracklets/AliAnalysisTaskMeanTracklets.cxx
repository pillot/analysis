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

ClassImp(AliAnalysisTaskMeanTracklets)

//________________________________________________________________________
AliAnalysisTaskMeanTracklets::AliAnalysisTaskMeanTracklets() :
AliAnalysisTaskSE(),
fEvents(0x0),
fhNtrk(0x0),
fhNtrkCorr(0x0),
fpMeanNtrkVsZvtx(0x0),
fpMeanNtrkVsZvtxCorr(0x0),
fpMeanNtrkVsZvtxRef(0x0),
fMeanNtrkRef(-1.),
fRandom(new TRandom3(0)),
fTrigger(""),
fRejectNSD(kFALSE),
fRejectPUFromSPD(kFALSE),
fSelectSPDVtxQA(kTRUE)
{
  // Dummy constructor
}

//________________________________________________________________________
AliAnalysisTaskMeanTracklets::AliAnalysisTaskMeanTracklets(const char *name) :
AliAnalysisTaskSE(name),
fEvents(0x0),
fhNtrk(0x0),
fhNtrkCorr(0x0),
fpMeanNtrkVsZvtx(0x0),
fpMeanNtrkVsZvtxCorr(0x0),
fpMeanNtrkVsZvtxRef(0x0),
fMeanNtrkRef(-1.),
fRandom(new TRandom3(0)),
fTrigger(""),
fRejectNSD(kFALSE),
fRejectPUFromSPD(kFALSE),
fSelectSPDVtxQA(kTRUE)
{
  /// Constructor
  
  // Output slot #1 writes into a AliCounterCollection
  DefineOutput(1,AliCounterCollection::Class());
  // Output slot #2 writes into a THnSparse
  DefineOutput(2,THnSparse::Class());
  // Output slot #3 writes into a THnSparse
  DefineOutput(3,THnSparse::Class());
  // Output slot #4 writes into a TProfile
  DefineOutput(4,TProfile::Class());
  // Output slot #5 writes into a TProfile
  DefineOutput(5,TProfile::Class());
  
}

//________________________________________________________________________
AliAnalysisTaskMeanTracklets::~AliAnalysisTaskMeanTracklets()
{
  /// Destructor
  if (!AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
    delete fEvents;
    delete fhNtrk;
    delete fhNtrkCorr;
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
  
  // prepare binning for THnSparse
  // 1: N SPD tracklets
  // 2: zVtx
  // 3: run
  // 4: N MC particles
  const Int_t nDims = 4;
  Int_t nBins[nDims] = {350, 80, 200000, 350};
  Double_t xMin[nDims] = {-0.5, -10., 99999.5, -0.5};
  Double_t xMax[nDims] = {349.5, 10., 299999.5, 349.5};
  
  // create histogram
  fhNtrk = new THnSparseT<TArrayF>("hNtrk", "N SPD tracklets", nDims, nBins, xMin, xMax);
  fhNtrk->Sumw2();
  fhNtrkCorr = new THnSparseT<TArrayF>("hNtrkCorr", "Corrected N SPD tracklets", nDims, nBins, xMin, xMax);
  fhNtrkCorr->Sumw2();
  
  fpMeanNtrkVsZvtx = new TProfile("fpMeanNtrkVsZvtx", "<Ntrk> vs Zvtx", nBins[1], xMin[1], xMax[1]);
  fpMeanNtrkVsZvtxCorr = new TProfile("fpMeanNtrkVsZvtxCorr", "corrected <Ntrk> vs Zvtx", nBins[1], xMin[1], xMax[1]);
  
  // Post data at least once per task to ensure data synchronisation (required for merging)
  PostData(1, fEvents);
  PostData(2, fhNtrk);
  PostData(3, fhNtrkCorr);
  PostData(4, fpMeanNtrkVsZvtx);
  PostData(5, fpMeanNtrkVsZvtxCorr);
  
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
  AliMCEvent *mcEvt = MCEvent();
  if (fRejectNSD && mcEvt) {
    AliGenEventHeader* mcGenH = fMCEvent->GenEventHeader();
    if (mcGenH && mcGenH->InheritsFrom(AliGenPythiaEventHeader::Class())) {
      AliGenPythiaEventHeader *hPythia = static_cast<AliGenPythiaEventHeader*>(mcGenH);
      if (hPythia->ProcessType() == 92 || hPythia->ProcessType() == 93) return; // single-diffractive
    }
  }
  
  // select a specific trigger
  if (!fTrigger.IsNull()) {
    TString trigger = evt->GetFiredTriggerClasses();
    if (mcEvt && fTrigger == "MB" && !(trigger.Contains("V0L") && trigger.Contains("V0R"))) return;
    else if (!trigger.Contains(fTrigger.Data())) return;
  }
  
  // reject events with pile-up from SPD
  if (fRejectPUFromSPD && evt->IsPileupFromSPDInMultBins()) return;
  
  // apply vertex QA selection
  Bool_t vtxQA = kFALSE;
  const AliVVertex* vtx = evt->GetPrimaryVertexSPD();
  if (vtx && vtx->GetNContributors() > 0) {
    Double_t cov[6]={0};
    vtx->GetCovarianceMatrix(cov);
    if (TMath::Sqrt(cov[5]) <= 0.25) vtxQA = kTRUE;
  }
  if (fSelectSPDVtxQA && !vtxQA) return;
  
  // get and select on vertex position
  Double_t zVtx = vtxQA ? vtx->GetZ() : 0.;
  if (zVtx <= -10. || zVtx >= 10.) return;
  
  // fill event counters for selected events
  fEvents->Count(Form("event:selected/run:%d", fCurrentRunNumber));
  
  // get number of charged particules
  Int_t NMCpartInEtaRange(0);
  if (mcEvt) {
    Int_t nMCpart = mcEvt->GetNumberOfTracks();
    for (Int_t i = 0; i < nMCpart ; ++i) {
      AliAODMCParticle *p = static_cast<AliAODMCParticle*>(mcEvt->GetTrack(i));
      if (!p->IsPhysicalPrimary()) continue;
      if (p->Charge() == 0) continue;
      if (TMath::Abs(p->Eta()) > 1.) continue;
      ++NMCpartInEtaRange;
    }
  }
  
  // get the number of SPD tracklets
  AliAODTracklets *tracklets = static_cast<const AliAODEvent*>(evt)->GetTracklets();
  Int_t nTrk = tracklets->GetNumberOfTracklets();
  Int_t nTrkInEtaRange(0);
  for (Int_t i = 0; i < nTrk; i++) {
    Double_t eta = -TMath::Log(TMath::Tan(tracklets->GetTheta(i)/2.));
    if ( eta < -1. || eta > 1. ) continue;
    ++nTrkInEtaRange;
  }
  
  // get the corrected number of SPD tracklets (no correction without a good SPD vertex)
  Int_t nTrkInEtaRangeCorr = nTrkInEtaRange;
  if (vtxQA) nTrkInEtaRangeCorr = fpMeanNtrkVsZvtxRef ? GetCorrectedNtrk(nTrkInEtaRange,zVtx) : GetCorrectedNtrkFromMultSel();
  
  // fill histograms
  Double_t EventInfo[4];
  EventInfo[0] = nTrkInEtaRange;
  EventInfo[1] = zVtx;
  EventInfo[2] = fCurrentRunNumber;
  EventInfo[3] = NMCpartInEtaRange;
  fhNtrk->Fill(EventInfo);
  fpMeanNtrkVsZvtx->Fill(zVtx,nTrkInEtaRange);
  if (nTrkInEtaRangeCorr >= 0) {
    EventInfo[0] = nTrkInEtaRangeCorr;
    fhNtrkCorr->Fill(EventInfo);
    fpMeanNtrkVsZvtxCorr->Fill(zVtx,nTrkInEtaRangeCorr);
  }
  
  // Post data
  PostData(1, fEvents);
  PostData(2, fhNtrk);
  PostData(3, fhNtrkCorr);
  PostData(4, fpMeanNtrkVsZvtx);
  PostData(5, fpMeanNtrkVsZvtxCorr);
  
}

//________________________________________________________________________
Int_t AliAnalysisTaskMeanTracklets::GetCorrectedNtrk(Int_t nTrkInEtaRange, Double_t zVtx)
{
  /// return the number of SPD tracklet corrected for the Zvtx dependence of <Ntrk>.
  /// return -999 in case of invalid correction
  
  Double_t meanNtrk = fpMeanNtrkVsZvtxRef->GetBinContent(fpMeanNtrkVsZvtxRef->FindBin(zVtx));
  if (meanNtrk < 1.e-6) return -999;
  
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

//________________________________________________________________________
Int_t AliAnalysisTaskMeanTracklets::GetCorrectedNtrkFromMultSel()
{
  /// return the corrected number of SPD tracklet from the multiplicity framework
  /// return -999 if not set
  
  AliMultSelection *multSelection = static_cast<AliMultSelection*>(InputEvent()->FindListObject("MultSelection"));
  if (!multSelection) return -999;
  
  return static_cast<Int_t>(multSelection->GetEstimator("SPDTracklets")->GetValue());
  
}

