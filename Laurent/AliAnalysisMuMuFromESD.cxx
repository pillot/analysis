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

// $Id$

#include "AliAnalysisMuMuFromESD.h"

#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliHistogramCollection.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliESDMuonTrack.h"
#include "AliVParticle.h"
#include "Riostream.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLorentzVector.h"
#include "TMap.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TPaveText.h"
#include <cassert>
#include "AliPhysicsSelection.h"
#include "AliPhysicsSelectionTask.h"
#include "AliBackgroundSelection.h"
#include "AliCounterCollection.h"
#include "THashList.h"

#include "TGrid.h"

ClassImp(AliAnalysisMuMuFromESD)

//_____________________________________________________________________________
AliAnalysisMuMuFromESD::AliAnalysisMuMuFromESD() 
: AliAnalysisMuMu()
{
  // default ctor
}

//_____________________________________________________________________________
AliAnalysisMuMuFromESD::AliAnalysisMuMuFromESD(TList* triggerClasses) 
:AliAnalysisMuMu(triggerClasses,kTRUE)
{
  // Constructor
  
  if ( fAA ) 
  {
    AddSingleCut("MATCHLOWETARABS",kAll | kMatchedLow | kEta | kRabs );
    
    AddSingleCut("MATCHLOWETARABSDCA",kAll | kMatchedLow | kEta | kRabs | kDCA );
    
    AddPairCut("MATCHLOWETARABSY",
               kAll | kMatchedLow | kEta | kRabs,
               kAll | kMatchedLow | kEta | kRabs | kPairRapidity
               );
    
    AddPairCut("MATCHLOWETARABS",
               kAll | kMatchedLow | kEta | kRabs);    
  }
  else
  {
    AddSingleCut("ALL",kAll);
    AddPairCut("ALL",kAll);
    AddSingleCut("MATCHRABS",kAll|kMatched|kRabs);
    AddPairCut("MATCHRABS",kAll|kMatched|kRabs);      
  }
  
  if ( fAA ) 
  {
    AddGlobalEventSelection("SD2");
  }
  else
  {
    AddGlobalEventSelection("ALL");    
  }
  
  AddGlobalEventSelection("PS");
  
//  fBranchNames = "ESD:AliESDRun.,AliESDHeader.,MuonTracks";
}

//_____________________________________________________________________________
AliAnalysisMuMuFromESD::~AliAnalysisMuMuFromESD()
{
}

//_____________________________________________________________________________
void AliAnalysisMuMuFromESD::MuUserExec(Option_t*)
{
  // Main loop
  // Called for each event
  
  AliESDEvent* esd = dynamic_cast<AliESDEvent*>(InputEvent());
  
  if (esd)
  {      

    if ( fAA ) 
    {
    // consider only events with OSM2 fired
    
      UInt_t sd2 = (1<<12);
      UInt_t trigger = esd->GetHeader()->GetL0TriggerInputs();
      
      Bool_t ok = ( ( trigger & sd2 ) == sd2);


      if (!ok) return;

    }
        
    TString runNumber(RunNumber(*esd));
  
    TString triggers = esd->GetFiredTriggerClasses();
    TObjArray* a = triggers.Tokenize(" ");
    TIter next(a);
    TObjString* tname;
    
    AliInputEventHandler* inputEventHandler = static_cast<AliInputEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    
    Bool_t isPhysicsSelected =  ( inputEventHandler->IsEventSelected() & AliVEvent::kMB )
    || ( inputEventHandler->IsEventSelected() & AliVEvent::kMUON );

    TString centrality(DefaultCentralityName());
    

    if ( fAA ) 
    {
      // FIXME : should get it from centrality task...
      
      // Double_t cent(-100);
      //    Double_t v0mult = vzero->GetMTotV0A()+vzero->GetMTotV0C();
      
      //    for ( std::vector<double>::size_type i = 0 ; i < fCentralityLimits.size()-1 && cent < -10 ; ++i )
      //    {
      //      if ( v0mult > fCentralityAmplitude[i] ) 
      //      {
      //        cent = fCentralityLimits[i];
      //      }
      //    }
      //    
      // if ( cent > -1 ) 
      // {
      //   centrality = CentralityName(cent);
      // }
      
      
    }

    TString eventtype("all");
    TString et;
    
    if ( fAA ) eventtype="sd2";
    et = eventtype;
    et.ToUpper();
    
    while ( ( tname = static_cast<TObjString*>(next()) ) )
    {
      if ( !fTriggerClasses->FindObject(tname->GetName()) ) 
      {
        // skip the trigger classes we're not supposed to analyse

        continue;
      }

      fEventCounters->Count(Form("event:%s/trigger:%s/run:%d", eventtype.Data(),tname->GetName(), runNumber.Atoi()));
        
      AssertHistogramCollection(et.Data(),tname->String().Data());
        
      FillHistos(et.Data(),tname->String().Data(),centrality.Data(),*esd);    
      
      if ( isPhysicsSelected )
      {
        
        fEventCounters->Count(Form("event:ps/trigger:%s/run:%d", tname->GetName(), runNumber.Atoi()));
        
        AssertHistogramCollection("PS",tname->String().Data());
        
        FillHistos("PS",tname->String().Data(),centrality.Data(),*esd);            
      }
    }
    
    delete a;
  }
}

//_____________________________________________________________________________
const char*
AliAnalysisMuMuFromESD::RunNumber(const AliESDEvent& esd) const
{
  return Form("%09d",esd.GetRunNumber());
}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMuFromESD::TrackMatchCut(const AliESDMuonTrack& track) const
{
  return track.GetMatchTrigger() > 0;
}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMuFromESD::TrackMatchLowCut(const AliESDMuonTrack& track) const
{
  return track.GetMatchTrigger() > 1;
}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMuFromESD::TrackMatchHighCut(const AliESDMuonTrack& track) const
{
  return track.GetMatchTrigger() > 2;
}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMuFromESD::TrackRabsCut(const AliESDMuonTrack& track) const
{
  /// Cut between 2 and 9 degrees
  
  const Double_t kAbsZEnd = 505; // cm
  const Double_t angleMin = 2.0*TMath::DegToRad();
  const Double_t angleMax = 9.0*TMath::DegToRad();
  
  Double_t angle = TMath::ATan(track.GetRAtAbsorberEnd()/kAbsZEnd);
  
  return ( angle > angleMin && angle < angleMax );
}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMuFromESD::TrackPtCut(const AliESDMuonTrack& track) const
{
  return track.Pt() > 1.0;
}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMuFromESD::TrackChi2(const AliESDMuonTrack& track) const
{
  return track.GetNormalizedChi2() < 4.0;
}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMuFromESD::TrackEtaCut(const AliESDMuonTrack& track) const
{
  return track.Eta() < -2.5 && track.Eta() > -4;
}

//_____________________________________________________________________________
UInt_t AliAnalysisMuMuFromESD::GetTrackMask(const AliESDMuonTrack& track) const
{
  UInt_t m(kAll);
  
  if ( TrackRabsCut(track) ) m |= kRabs;
  
  if ( TrackPtCut(track) ) m |= kPt;  
  
  if ( TrackMatchCut(track) ) m |= kMatched;
  
  if ( TrackMatchLowCut(track) ) m |= kMatchedLow;
  
  if ( TrackMatchHighCut(track) ) m |= kMatchedHigh;
  
  if ( TrackEtaCut(track) ) m |= kEta;
  
  if ( TrackChi2(track) ) m |= kChi2;
  
  return m;
  
}

//_____________________________________________________________________________
void AliAnalysisMuMuFromESD::FillHistosForTrack(const char* physics,
                                                const char* triggerClassName, 
                                                const char* centrality,
                                                const AliESDMuonTrack& track,
                                                const char* /*runNumber*/)
{
  
  if ( track.ContainTrackerData() )
  {      
    TLorentzVector p;
    
    track.LorentzP(p);
    
    TString charge("Plus");
    
    if ( track.Charge() < 0 ) charge = "Minus";

    UInt_t mask = GetTrackMask(track);
    
    const Double_t kAbsZEnd = 505; // cm      
    Double_t theta = TMath::ATan(track.GetRAtAbsorberEnd()/kAbsZEnd);

    Double_t xdca = track.GetNonBendingCoorAtDCA();
    Double_t ydca = track.GetBendingCoorAtDCA();
    Double_t dca = TMath::Sqrt(xdca*xdca+ydca*ydca);
    
    TIter next(fSingleTrackCutNames);
    TObjString* str;
    
    while ( ( str = static_cast<TObjString*>(next()) ) )
    {
      Bool_t test = ( ( str->GetUniqueID() & mask ) == str->GetUniqueID() );

      if ( test ) 
      {
        
        Histo(physics,triggerClassName,centrality,str->String().Data(),Form("PtEtaMu%s",charge.Data()))->Fill(p.Eta(),p.Pt());
        
        Histo(physics,triggerClassName,centrality,str->String().Data(),Form("PEtaMu%s",charge.Data()))->Fill(p.Eta(),p.P());
        
        
        Histo(physics,triggerClassName,centrality,str->String().Data(),Form("XYdcaMu%s",charge.Data()))->Fill(ydca,xdca);
        
        Histo(physics,triggerClassName,centrality,str->String().Data(),Form("Chi2Mu%s",charge.Data()))->Fill(track.GetNormalizedChi2());
        
        if ( theta >= Deg2() && theta < Deg3() )         
        {
          //        Histo(physics,triggerClassName,centrality,str->String().Data(),Form("PDCAP23Mu%s",charge.Data()))->Fill(p.P(),p.P()*dca);
          Histo(physics,triggerClassName,centrality,str->String().Data(),Form("dcaP23Mu%s",charge.Data()))->Fill(p.P(),dca);
          if ( p.Pt() > 2 )
          {
            Histo(physics,triggerClassName,centrality,str->String().Data(),Form("dcaPwPtCut23Mu%s",charge.Data()))->Fill(p.P(),dca);
          }
        }
        else if ( theta >= Deg3() && theta < Deg10() )
        {
          //        Histo(physics,triggerClassName,centrality,str->String().Data(),Form("PDCAP310Mu%s",charge.Data()))->Fill(p.P(),p.P()*dca);        
          Histo(physics,triggerClassName,centrality,str->String().Data(),Form("dcaP310Mu%s",charge.Data()))->Fill(p.P(),dca);
          if ( p.Pt() > 2 )
          {
            Histo(physics,triggerClassName,centrality,str->String().Data(),Form("dcaPwPtCut310Mu%s",charge.Data()))->Fill(p.P(),dca);
          }
        }        
      }
    }
  }
}

//_____________________________________________________________________________
void AliAnalysisMuMuFromESD::FillHistos(const char* physics, 
                                        const char* triggerClassName, 
                                        const char* centrality,
                                        const AliESDEvent& esd)
{
  const AliESDVertex* vertex = esd.GetPrimaryVertexSPD();
  
  if ( vertex ) 
  {
    Histo(physics,triggerClassName,centrality,"Xvertex")->Fill(vertex->GetX());
    Histo(physics,triggerClassName,centrality,"Yvertex")->Fill(vertex->GetY());
    Histo(physics,triggerClassName,centrality,"Zvertex")->Fill(vertex->GetZ());
  }
  
  AliESDVZERO* vzero = esd.GetVZEROData();
  
  if ( vzero ) 
  {
    Histo(physics,triggerClassName,centrality,"V0Time")->Fill(vzero->GetV0ATime(),vzero->GetV0CTime());
    
    Double_t v0mult = vzero->GetMTotV0A()+vzero->GetMTotV0C();
    
    Histo(physics,triggerClassName,centrality,"V0Amplitude")->Fill(v0mult);    
  }
  
  // Track loop
  
  TString runNumber(RunNumber(esd));
  
  for (Int_t i = 0; i < esd.GetNumberOfMuonTracks(); ++i) 
  {
    AliESDMuonTrack* tracki = esd.GetMuonTrack(i);

    UInt_t maski = GetTrackMask(*tracki);
    
    FillHistosForTrack(physics,triggerClassName,centrality,*tracki,runNumber.Data());
    
    if ( tracki->ContainTrackerData() ) 
    {
      TLorentzVector pi;

      tracki->LorentzP(pi);
    
      for (Int_t j = i+1; j < esd.GetNumberOfMuonTracks(); ++j) 
      {
        AliESDMuonTrack* trackj = esd.GetMuonTrack(j);
        
        if ( trackj->ContainTrackerData() ) 
        {
          UInt_t maskj = GetTrackMask(*trackj);

          TLorentzVector pj;
          trackj->LorentzP(pj);

          pj += pi;
        
          TIter next(fPairTrackCutNames);
          TObjString* str;
          
          while ( ( str = static_cast<TObjString*>(next()) ) )
          {
            UInt_t singleTrackMask(0);
            UInt_t pairMask(0);
            
            DecodePairCutMask(str->GetUniqueID(),singleTrackMask,pairMask);
            
            Bool_t testi = ( ( maski & singleTrackMask ) == singleTrackMask ) ;
            Bool_t testj = ( ( maskj & singleTrackMask ) == singleTrackMask ) ;
            Bool_t testij(kTRUE);
            
            if ( pairMask > 0 )
            {
              testij = ( (maski & pairMask) == pairMask ) &&
              ( (maskj & pairMask) == pairMask );
            }
            
            if ( ( testi || testj ) && testij )
            {
              if ( tracki->Charge() != trackj->Charge() )
              {
                Histo(physics,triggerClassName,centrality,str->String().Data(),"MinvUSPt")->Fill(pj.Pt(),pj.M());            
              }
              else if ( tracki->Charge() > 0 && trackj->Charge() > 0 )
              {
                Histo(physics,triggerClassName,centrality,str->String().Data(),"MinvPPPt")->Fill(pj.Pt(),pj.M());                        
              }
              else
              {
                Histo(physics,triggerClassName,centrality,str->String().Data(),"MinvMMPt")->Fill(pj.Pt(),pj.M());                        
              }
            }
          }
        }
      }
    }
  } //track loop
  
  
  
}

//==============================================================================
//
//   AddTask ....
//
//==============================================================================

#include "AliAnalysisMuMuFromESD.h"
#include "AliAnalysisManager.h"
#include "TString.h"
#include "AliAnalysisDataContainer.h"

AliAnalysisTask* AddTaskMuMuFromESD(const char* outputname, TList* triggerClassesToConsider)
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskMuMu", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskMuMu", "This task requires an input event handler");
    return NULL;
  }
  TString inputDataType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  
  // Configure analysis
  //===========================================================================
  
  
  AliAnalysisMuMuFromESD* task = new AliAnalysisMuMuFromESD(triggerClassesToConsider);
  mgr->AddTask(task);  
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();  
  AliAnalysisDataContainer *coutput = 
  mgr->CreateContainer("chist", TList::Class(),
                       AliAnalysisManager::kOutputContainer,outputname);
  
  // Connect input/output
  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 1, coutput);
  
  return task;
}
