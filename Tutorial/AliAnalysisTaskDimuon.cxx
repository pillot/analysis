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
#include "TCanvas.h"

// STEER includes
#include "AliVEvent.h"
#include "AliAODEvent.h"
#include "AliVParticle.h"
#include "AliCounterCollection.h"

// ANALYSIS includes
#include "AliAnalysisTaskDimuon.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataSlot.h"
#include "AliAnalysisDataContainer.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisMuonUtility.h"
#include "AliMultSelection.h"

ClassImp(AliAnalysisTaskDimuon)

//________________________________________________________________________
AliAnalysisTaskDimuon::AliAnalysisTaskDimuon() :
AliAnalysisTaskSE(),
fEvents(0x0),
fhOS(0x0),
fhLS(0x0),
fMuonTrackCuts(0x0)
{
  // Dummy constructor
}

//________________________________________________________________________
AliAnalysisTaskDimuon::AliAnalysisTaskDimuon(const char *name) :
AliAnalysisTaskSE(name),
fEvents(0x0),
fhOS(0x0),
fhLS(0x0),
fMuonTrackCuts(0x0)
{
  /// Constructor
  
  // Output slot #1 writes into a AliCounterCollection
  DefineOutput(1,AliCounterCollection::Class());
  // Output slot #2 writes into a THnSparse
  DefineOutput(2,THnSparse::Class());
  // Output slot #3 writes into a THnSparse
  DefineOutput(3,THnSparse::Class());
  
}

//________________________________________________________________________
AliAnalysisTaskDimuon::~AliAnalysisTaskDimuon()
{
  /// Destructor
  if (!AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
    delete fEvents;
    delete fhOS;
    delete fhLS;
  }
  delete fMuonTrackCuts;
}

//___________________________________________________________________________
void AliAnalysisTaskDimuon::NotifyRun()
{
  /// Load the trackCuts
  
  if (!fMuonTrackCuts) AliFatal("You must specify the requested selections (AliMuonTrackCut obj is missing)");
  fMuonTrackCuts->SetRun(fInputHandler);
  fMuonTrackCuts->Print();
  
}

//___________________________________________________________________________
void AliAnalysisTaskDimuon::UserCreateOutputObjects()
{
  /// Create output objects
  
  // events analyzed
  fEvents = new AliCounterCollection(GetOutputSlot(1)->GetContainer()->GetName());
  fEvents->AddRubric("event", "MUL/MLL/any");
  fEvents->AddRubric("cent", "0-10/10-20/20-30/30-40/40-50/50-60/60-70/70-80/80-90/90-100/other");
  fEvents->Init();
  
  // prepare binning for THnSparse
  // 1: centrality
  // 2: mass
  // 3: pT
  const Int_t nDims = 3;
  Int_t nBins[nDims] = {10, 350, 200};
  Double_t xMin[nDims] = {0., 0., 0.};
  Double_t xMax[nDims] = {100., 7., 20.};
  
  // create histogram
  fhOS = new THnSparseT<TArrayF>("hOS", "OS dimuons", nDims, nBins, xMin, xMax);
  fhLS = new THnSparseT<TArrayF>("hLS", "LS dimuons", nDims, nBins, xMin, xMax);
  
  // Post data at least once per task to ensure data synchronisation (required for merging)
  PostData(1, fEvents);
  PostData(2, fhOS);
  PostData(3, fhLS);
  
}

//________________________________________________________________________
void AliAnalysisTaskDimuon::UserExec(Option_t *)
{
  /// Called for each event
  
  // get the input event
  AliVEvent *evt = InputEvent();
  
  // cast it to AOD and make sure it is actually AOD
  if (!dynamic_cast<AliAODEvent*>(evt)) return;
  
  // get physics selected trigger
  Bool_t trigOS = (fInputHandler->IsEventSelected() & AliVEvent::kMuonUnlikeLowPt7);
  Bool_t trigLS = (fInputHandler->IsEventSelected() & AliVEvent::kMuonLikeLowPt7);
  
  // get centrality percentile
  Double_t trackInfo[3];
  AliMultSelection *multSelection = static_cast<AliMultSelection*>(evt->FindListObject("MultSelection"));
  trackInfo[0] = multSelection ? multSelection->GetMultiplicityPercentile("V0M") : -1.;
  TString centKey = (trackInfo[0] >= 0 && trackInfo[0] < 100) ? Form("cent:%d-%d", ((Int_t)trackInfo[0])/10*10, (((Int_t)trackInfo[0])/10+1)*10) : "cent:other";
  
  // fill event counters
  fEvents->Count(Form("event:any/%s", centKey.Data()));
  if (trigOS) fEvents->Count(Form("event:MUL/%s", centKey.Data()));
  if (trigLS) fEvents->Count(Form("event:MLL/%s", centKey.Data()));
  
  Int_t nTracks = AliAnalysisMuonUtility::GetNTracks(evt);
  for (Int_t iTrack1 = 0; iTrack1 < nTracks; ++iTrack1) {
    
    AliVParticle *track1 = AliAnalysisMuonUtility::GetTrack(iTrack1, evt);
    if (!fMuonTrackCuts->IsSelected(track1)) continue;
    
    for (Int_t iTrack2 = iTrack1 + 1; iTrack2 < nTracks; iTrack2++) {
      
      AliVParticle *track2 = AliAnalysisMuonUtility::GetTrack(iTrack2, evt);
      if (!fMuonTrackCuts->IsSelected(track2)) continue;
      
      // get dimuon mass/pT
      TLorentzVector dimuV = AliAnalysisMuonUtility::GetTrackPair(track1, track2);
      trackInfo[1] = dimuV.M();
      trackInfo[2] = dimuV.Pt();
      Short_t charge = track1->Charge()*track2->Charge();
      
      
      // fill OS or LS histo
      if (trigOS && charge < 0) fhOS->Fill(trackInfo);
      else if (trigLS && charge > 0) fhLS->Fill(trackInfo);
      
    }
    
  }
  
  // Post data
  PostData(1, fEvents);
  PostData(2, fhOS);
  PostData(3, fhLS);
  
}

//________________________________________________________________________
void AliAnalysisTaskDimuon::Terminate(Option_t *)
{
  /// Work on final output objects (after merging)
  
  // initialized internal data members
  // (needed because this instance of the task is not necessarily the one which produced and filled the output object)
  fhOS = dynamic_cast<THnSparse*>(GetOutputData(2));
  fhLS = dynamic_cast<THnSparse*>(GetOutputData(3));
  if (!fhOS || !fhLS) return;
  
  // draw histogram
  TCanvas *c = new TCanvas("c", "c", 900, 600);
  c->Divide(3,2);
  c->cd(1);
  fhOS->Projection(0,"e")->Draw();
  c->cd(2);
  gPad->SetLogy();
  fhOS->Projection(1,"e")->Draw();
  c->cd(3);
  gPad->SetLogy();
  fhOS->Projection(2,"e")->Draw();
  c->cd(4);
  fhLS->Projection(0,"e")->Draw();
  c->cd(5);
  gPad->SetLogy();
  fhLS->Projection(1,"e")->Draw();
  c->cd(6);
  gPad->SetLogy();
  fhLS->Projection(2,"e")->Draw();
  
}

