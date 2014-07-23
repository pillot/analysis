/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id$ */

//-----------------------------------------------------------------------------
/// \class AliAnalysisTaskMTRSign
/// Check the deviation in MTR
///
/// \author Diego Stocco
//-----------------------------------------------------------------------------

#define AliAnalysisTaskMTRSign_cxx

#include "AliAnalysisTaskMTRSign.h"

// ROOT includes
#include "TROOT.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TObjString.h"
#include "TAxis.h"

// STEER includes
#include "AliESDMuonTrack.h"
#include "AliVParticle.h"
#include "AliLog.h"
#include "AliInputEventHandler.h"

// ANALYSIS includes
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisDataSlot.h"
#include "AliAnalysisDataContainer.h"

// CORRFW includes
#include "AliCFGridSparse.h"

// PWG includes
#include "AliMergeableCollection.h"
#include "AliAnalysisMuonUtility.h"
#include "AliMuonEventCuts.h"
#include "AliMuonTrackCuts.h"

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskMTRSign) // Class implementation in ROOT context
/// \endcond


//________________________________________________________________________
AliAnalysisTaskMTRSign::AliAnalysisTaskMTRSign() :
  AliAnalysisTaskSE(),
  fMuonEventCuts(0x0),
  fMuonTrackCuts(0x0),
  fUseMCLabel(kFALSE),
  fMergeableCollection(0x0)
{
  /// Default ctor.
}

//________________________________________________________________________
AliAnalysisTaskMTRSign::AliAnalysisTaskMTRSign(const char *name) :
  AliAnalysisTaskSE(name),
  fMuonEventCuts(new AliMuonEventCuts("muEvtCuts","muEvtCuts")),
  fMuonTrackCuts(new AliMuonTrackCuts("muTrackCuts","muTrackCuts")),
  fUseMCLabel(kFALSE),
  fMergeableCollection(0x0)
{
  //
  /// Constructor.
  //
  DefineOutput(1, AliMergeableCollection::Class());
}


//________________________________________________________________________
AliAnalysisTaskMTRSign::~AliAnalysisTaskMTRSign()
{
  //
  /// Destructor
  //

  delete fMuonEventCuts;
  delete fMuonTrackCuts;
  
  // For proof: do not delete output containers
  if ( ! AliAnalysisManager::GetAnalysisManager() || ! AliAnalysisManager::GetAnalysisManager()->IsProofMode() ) {
    delete fMergeableCollection;
  }
}

//________________________________________________________________________
void AliAnalysisTaskMTRSign::NotifyRun()
{
  /// Notofy run
  fMuonTrackCuts->SetRun(fInputHandler);
}

//________________________________________________________________________
AliCFGridSparse* AliAnalysisTaskMTRSign::GetCFGridSparse ( TString identifier )
{
  //
  /// Get CF grid sparse object
  /// (create collection if necessary)
  //
  
  TString sparseName = "MTRSign";
  AliCFGridSparse* gridSparse = static_cast<AliCFGridSparse*>(fMergeableCollection->GetObject(identifier.Data(),sparseName.Data()));
  if ( ! gridSparse ) {
    Int_t nPtBins = 150;
    Double_t ptMin = 0., ptMax = 15.;
    TString ptTitle("p_{T}"), ptUnits("GeV/c");
    
    Int_t nPBins = 150.;
    Double_t pMin = 0., pMax = 150.;
    TString pTitle("p"), pUnits("GeV/c");
    
    Int_t nChargeBins = 2;
    Double_t chargeMin = -2., chargeMax = 2.;
    TString chargeTitle("charge"), chargeUnits("e");
    
    Int_t nMatchTrigBins = 4;
    Double_t matchTrigMin = -0.5, matchTrigMax = 3.5;
    TString matchTrigTitle("matchTrig"), matchTrigUnits("");
    
    Int_t nTrigSignBins = 3;
    Double_t trigSignMin = -1.5, trigSignMax = 1.5;
    TString trigSignTitle("trigSign"), trigSignUnits("a.u.");
    
    Int_t nLoCircuitBins = 234;
    Double_t loCircuitMin = 0.5, loCircuitMax = 234.5;
    TString loCircuitTitle("LoCircuit"), loCircuitUnits("");
    
    Int_t nbins[kNvars] = {nPtBins, nPBins, nChargeBins, nMatchTrigBins, nTrigSignBins, nLoCircuitBins};
    Double_t xmin[kNvars] = {ptMin, pMin, chargeMin, matchTrigMin, trigSignMin, loCircuitMin};
    Double_t xmax[kNvars] = {ptMax, pMax, chargeMax, matchTrigMax, trigSignMax, loCircuitMax};
    TString axisTitle[kNvars] = {ptTitle, pTitle, chargeTitle, matchTrigTitle, trigSignTitle, loCircuitTitle};
    TString axisUnits[kNvars] = {ptUnits, pUnits, chargeUnits, matchTrigUnits, trigSignUnits, loCircuitUnits};
    
    gridSparse = new AliCFGridSparse(sparseName.Data(),"MTR sign",kNvars,nbins);
    
    for ( Int_t idim = 0; idim<kNvars; idim++){
      TString currStr = Form("%s (%s)", axisTitle[idim].Data(), axisUnits[idim].Data());
      currStr.ReplaceAll("()","");
      
      gridSparse->SetVarTitle(idim, currStr.Data());
      gridSparse->SetBinLimits(idim, xmin[idim], xmax[idim]);
    }
    
    fMergeableCollection->Adopt(identifier.Data(),gridSparse);
  }
  
  return gridSparse;
}

//___________________________________________________________________________
void AliAnalysisTaskMTRSign::FinishTaskOutput()
{
  //
  /// Remove empty histograms to reduce the number of histos to be merged
  //
  
  fMergeableCollection->PruneEmptyObjects();
}

//___________________________________________________________________________
void AliAnalysisTaskMTRSign::UserCreateOutputObjects()
{
  fMergeableCollection = new AliMergeableCollection(GetOutputSlot(1)->GetContainer()->GetName());
  
  fMuonEventCuts->Print();
  fMuonTrackCuts->Print("mask");
  
  PostData(1,fMergeableCollection);
}

//________________________________________________________________________
void AliAnalysisTaskMTRSign::UserExec(Option_t * /*option*/)
{
  //
  /// Fill output objects
  //
  
//  if ( InputEvent()->IsA() != AliESDEvent::Class() ) {
//    AliError("The task runs only on AliESDs.root!");
//    return;
//  }
  
  if ( ! fMuonEventCuts->IsSelected(fInputHandler) ) return;
  
  TString physSel = ( fInputHandler->IsEventSelected() & AliVEvent::kAny ) ? "PhysSelPass" : "PhysSelReject";
  Double_t centrality = fMuonEventCuts->GetCentrality(InputEvent());
  Int_t centralityBin = fMuonEventCuts->GetCentralityClasses()->FindBin(centrality);
  TString centralityBinLabel = fMuonEventCuts->GetCentralityClasses()->GetBinLabel(centralityBin);

  Double_t containerInput[kNvars];
  
  const TObjArray* selectTrigClasses = fMuonEventCuts->GetSelectedTrigClassesInEvent(InputEvent());
  
  Int_t nTracks = AliAnalysisMuonUtility::GetNTracks(InputEvent());
  
  for (Int_t itrack = 0; itrack < nTracks; itrack++) {
    AliVParticle* track = AliAnalysisMuonUtility::GetTrack(itrack,InputEvent());
    
    if ( ! fMuonTrackCuts->IsSelected(track) ) continue;
    if (fUseMCLabel && track->GetLabel() < 0) continue;

    containerInput[kHvarPt]        = track->Pt();
    containerInput[kHvarP]         = track->P();
    containerInput[kHvarCharge]    = track->Charge();
    containerInput[kHvarMatchTrig] = AliAnalysisMuonUtility::GetMatchTrigger(track);
    containerInput[kHvarTrigSign]  = TriggerDevSign(track);
//    containerInput[kHvarTrigSign]  = AliAnalysisMuonUtility::GetTrigDevSign(track); // REMEMBER TO CHANGE when commit done
    containerInput[kHvarLoCircuit] = AliAnalysisMuonUtility::GetLoCircuit(track);
    
//    if ( containerInput[kHvarCharge] == containerInput[kHvarTrigSign] ) printf("BADSIGN: %s  %lli\n",CurrentFileName(),Entry()); // REMEMBER TO CUT
//    Int_t devSignRef = TriggerDevSign(track), devSignNew = AliAnalysisMuonUtility::GetTrigDevSign(track); // REMEMBER TO CUT
//    if ( devSignRef != devSignNew ) printf("\n     AAAARGHHH!!!!! %i  %i\n\n",devSignRef,devSignNew); // REMEMBER TO CUT
//    printf("Entry %lli DevSign %i\n",Entry(),(Int_t)containerInput[kHvarTrigSign]); // REMEMBER TO CUT

    TIter next(selectTrigClasses);
    TObjString* trigClassName = 0x0;
    while ( ( trigClassName = static_cast<TObjString*>(next()) ) ) {
      if ( selectTrigClasses->GetEntries() == 3 && ! trigClassName->GetString().Contains("&") ) continue; // REMEMBER TO CUT
      if ( fMuonTrackCuts->TrackPtCutMatchTrigClass(track, fMuonEventCuts->GetTrigClassPtCutLevel(trigClassName->GetString())) ) {
        GetCFGridSparse(Form("/%s/%s/%s/",physSel.Data(),trigClassName->GetString().Data(),centralityBinLabel.Data()))->Fill(containerInput);
      }
    } // loop on triggers
  } // loop on tracks
  
  PostData(1,fMergeableCollection);
}

// REMEMBER TO COMMENT function when commit done
//________________________________________________________________________
Int_t AliAnalysisTaskMTRSign::TriggerDevSign ( AliVParticle *track ) const
{
  /// get the sign (+-1) of track deviation in the trigger (0 = unknown)
  
  AliESDMuonTrack *esdTrack = static_cast<AliESDMuonTrack*>(track);
  
  Int_t deviation = esdTrack->LoDev();
  if ( deviation > 15 ) return 1;
  else if ( deviation < 15 ) return -1;
  else return 0;
}

//________________________________________________________________________
void AliAnalysisTaskMTRSign::Terminate(Option_t *)
{
  //
  /// Draw some histograms at the end.
  //

  fMergeableCollection = dynamic_cast<AliMergeableCollection*>(GetOutputData(1));
  if ( ! fMergeableCollection ) return;
  
//  fMergeableCollection->Print("*");
  
  
  TObjArray* listOfKeys = fMergeableCollection->SortAllIdentifiers();
  TIter next(listOfKeys);
  
  TObjString* objString = 0x0;
  while ( ( objString = static_cast<TObjString*>(next()) ) ) {
    TString identifier = objString->GetString();
    AliCFGridSparse* gridSparse = GetCFGridSparse(identifier);
    identifier.ReplaceAll("/","");
    TString histoName = Form("%s_signCorr", identifier.Data());
    TH1* histo = gridSparse->Project(kHvarCharge,kHvarTrigSign);
    if ( ! histo ) continue;
    histo->SetName(histoName.Data());

    TString currName = Form("%s_%s", GetName(), histoName.Data());
    TCanvas* can = new TCanvas(currName.Data(),currName.Data(),200,10,600,600);
    can->SetLogz();
    
    histo->Draw("TEXTCOLZ");
  } // loop on grid sparse
  
  delete listOfKeys;
  
}
