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
#include "TH1F.h"
#include "TCanvas.h"

// STEER includes
#include "AliVEvent.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"

// ANALYSIS includes
#include "AliAnalysisTaskTest.h"

ClassImp(AliAnalysisTaskTest)

//________________________________________________________________________
AliAnalysisTaskTest::AliAnalysisTaskTest() :
AliAnalysisTaskSE(),
fhPt(0x0)
{
  // Dummy constructor
}

//________________________________________________________________________
AliAnalysisTaskTest::AliAnalysisTaskTest(const char *name) :
AliAnalysisTaskSE(name),
fhPt(0x0)
{
  /// Constructor
  
  // Output slot #1 writes into a TList container
  DefineOutput(1,TH1F::Class());
  
}

//________________________________________________________________________
AliAnalysisTaskTest::~AliAnalysisTaskTest()
{
  /// Destructor
}

//___________________________________________________________________________
void AliAnalysisTaskTest::UserCreateOutputObjects()
{
  /// Create output objects
  
  // create histogram
  fhPt = new TH1F("hPt", "transverse momentum distribution;p_{t} (GeV/c)", 300, 0., 30.);
  
  // Post data at least once per task to ensure data synchronisation (required for merging)
  PostData(1, fhPt);
  
}

//________________________________________________________________________
void AliAnalysisTaskTest::UserExec(Option_t *)
{
  /// Called for each event
  
  // get the input event
  AliVEvent *evt = InputEvent();
  
  // cast it to AOD and make sure it is actually AOD
  AliAODEvent *aodEvt = dynamic_cast<AliAODEvent*>(evt);
  if (!aodEvt) return;
  
  // loop over AOD tracks
  Int_t nTracks = aodEvt->GetNTracks();
  for (Int_t iTrack = 0; iTrack < nTracks; ++iTrack) {
    
    // get the track
    AliAODTrack *track = aodEvt->GetTrack(iTrack);
    
    // select muon tracks
    if (!track->IsMuonTrack()) continue;
    
    // fill histogram
    fhPt->Fill(track->Pt());
    
  }
  
  // Post data
  PostData(1, fhPt);
  
}

//________________________________________________________________________
void AliAnalysisTaskTest::Terminate(Option_t *)
{
  /// Work on final output objects (after merging)
  
  // initialized internal data members
  // (needed because this instance of the task is not necessarily the one which produced and filled the output object)
  fhPt = dynamic_cast<TH1F*>(GetOutputData(1));
  if (!fhPt) return;
  
  // draw histogram
  TCanvas *c = new TCanvas();
  gPad->SetLogy();
  fhPt->Draw();
  
}

