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

#include "AliAnalysisMuMuFromAOD.h"

#include "AliAODDimuon.h"
#include "AliAODEvent.h"
#include "AliAODHeader.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliAODTrack.h"
#include "AliAnalysisManager.h"
#include "AliCentrality.h"
#include "AliCodeTimer.h"
#include "AliCounterCollection.h"
#include "AliHistogramCollection.h"
#include "AliLog.h"
#include "AliMCEvent.h"
#include "AliVParticle.h"
#include "Riostream.h"
#include "TCanvas.h"
#include "TDatabasePDG.h"
#include "TDirectory.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THashList.h"
#include "TLorentzVector.h"
#include "TMap.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TPaveText.h"
#include <cassert>
#include <set>

ClassImp(AliAnalysisMuMuFromAOD)

//_____________________________________________________________________________
AliAnalysisMuMuFromAOD::AliAnalysisMuMuFromAOD() 
: AliAnalysisMuMu(), fVertex(0), fGlobalEventSelectionName("ALL")
{
  // default ctor
}

//_____________________________________________________________________________
AliAnalysisMuMuFromAOD::AliAnalysisMuMuFromAOD(TList* triggerClasses) 
: AliAnalysisMuMu(triggerClasses,kFALSE), fVertex(0), fGlobalEventSelectionName("ALL")
{
  // Constructor
  // Must define here the single track and track pair cuts,
  // and also the global event selection "names"
  
  if ( fAA ) 
  {
    AddSingleCut("MATCHLOWETARABS",kAll | kMatchedLow | kEta | kRabs );
    
    AddSingleCut("MATCHLOWETARABSDCA",kAll | kMatchedLow | kEta | kRabs | kDCA );
    
    AddSingleCut("MATCHLOWETARABSCHI2",kAll | kMatchedLow | kEta | kRabs | kChi2 );
    
    AddSingleCut("MATCHLOWETARABSDCACHI2",kAll | kMatchedLow | kEta | kRabs | kDCA | kChi2 );
    
    AddPairCut("MATCHLOWETARABSY",
               kAll | kMatchedLow | kEta | kRabs,
               kAll | kMatchedLow | kEta | kRabs | kPairRapidity
               );
    
    AddPairCut("MATCHLOWETARABS",
               kAll | kMatchedLow | kEta | kRabs);

    fGlobalEventSelectionName="SD2";
  }
  else
  {
    AddSingleCut("MATCHRABS",kAll|kMatched|kRabs);
    AddPairCut("MATCHRABS",kAll|kMatched|kRabs);      
  }
  
  AddGlobalEventSelection(fGlobalEventSelectionName.Data());
}

//_____________________________________________________________________________
AliAnalysisMuMuFromAOD::~AliAnalysisMuMuFromAOD()
{
}

//_____________________________________________________________________________
void AliAnalysisMuMuFromAOD::DumpMC(const AliAODEvent& aod)
{
  // dump mc part (if present)
  
  //  MC header
  AliAODMCHeader *mcHeader = static_cast<AliAODMCHeader*>(aod.FindListObject(AliAODMCHeader::StdBranchName()));
  if(mcHeader) 
  {
    AliInfo(Form("================================ Generator %s x %e y %e z %e",
                 mcHeader->GetGeneratorName(),
                 mcHeader->GetVtxX(),
                 mcHeader->GetVtxY(),
                 mcHeader->GetVtxZ()));
  }
  
  // MC particles
  TClonesArray* mcParticles = static_cast<TClonesArray*>(aod.FindListObject(AliAODMCParticle::StdBranchName()));
  if (mcParticles) 
  {
    AliInfo(Form("%d mcparticles in that event",mcParticles->GetEntries()));
    
    // loop on muon tracks are find back their decay chain

    TIter next(aod.GetTracks());
    AliAODTrack* t;
//    Int_t nmu(0);
//    
//    AliInfo("Muon tracks");
//    
//    while ( ( t = static_cast<AliAODTrack*>(next()) ) )
//    {
//      if ( t->IsMuonTrack() ) 
//      {
//        AliInfo(Form("mutrack #%2d Pt %e Eta %e Charge %d Label %d",
//                     nmu,t->Pt(),t->Eta(),t->Charge(),t->GetLabel()));
//        ++nmu;
//      }
//    }
//    
//    AliInfo("MC particles...");
    TIter nextMC(mcParticles);
    AliAODMCParticle* p;
//    Int_t nmc(0);
//    
//    while ( ( p = static_cast<AliAODMCParticle*>(nextMC()) ) )
//    {
//      AliInfo(Form("mcpart #%4d",nmc));
//      p->Print();
//      ++nmc;
//    }
//    next.Reset();
    
    AliInfo("Muon track ancestors...");
        
    while ( ( t = static_cast<AliAODTrack*>(next()) ) )
    {
      if ( t->IsMuonTrack() ) 
      {
        AliInfo(Form("Track label %d",t->GetLabel()));
        Int_t label = t->GetLabel();
        
        while ( label >= 0 ) 
        {
          p = static_cast<AliAODMCParticle*>(mcParticles->At(label));
          if (p)
          {
            p->Print();
            label = p->GetMother();
          }
          else
          {
            label = -1;
          }
        }
      }
    }

  }

}

//_____________________________________________________________________________
void AliAnalysisMuMuFromAOD::MuUserExec(Option_t*)
{
  // Main loop
  // Called for each event
    
  AliAODEvent* aod = dynamic_cast<AliAODEvent*>(InputEvent());
  if ( !aod ) return;
  
//  AliAODHeader* header = aod->GetHeader();

//  AliInfo(Form("nmu=%d",header->GetNumberOfMuons()));
  
//  if ( header->GetNumberOfMuons()==0 ) return;

//  TIter nextV(aod->GetVertices());
//  AliAODVertex* v;
//  
//  while ( ( v = static_cast<AliAODVertex*>(nextV())) )
//  {
//    AliInfo(Form("vertex %e %e %e %s %s %d",v->GetX(),v->GetY(),v->GetZ(),
//                 v->GetTypeName(static_cast<AliAODVertex::AODVtx_t>(v->GetType())),
//                 v->GetTitle(),
//                 v->GetNContributors()
//                 ));
//  }

  
  if ( fAA ) 
  {
    // consider only events with OSM2 fired
    UInt_t trigger = aod->GetHeader()->GetL0TriggerInputs();
    UInt_t sd2 = (1<<12);  
    Bool_t ok = ( ( trigger & sd2 ) == sd2);  
    if (!ok) return;
  }
  
//  DumpMC(*aod);  
  
  TObjArray* a = aod->GetFiredTriggerClasses().Tokenize(" ");
  
  TIter next(a);
  TObjString* tname;
  
  AliCentrality* acent = aod->GetCentrality();
  
  Float_t fcent(-1);
  
  if (acent) fcent = acent->GetCentralityPercentile("V0M");
  
  Double_t cent(-100);
  
  if ( fcent > 0 )
  {
    for ( std::vector<double>::size_type i = 0 ; i < fCentralityLimits.size() && cent < 0 ; ++i )
    {
      if ( fcent < fCentralityLimits[i] ) 
      {
        cent = fCentralityLimits[i];
      }
    }
  }
  
  TString centrality(DefaultCentralityName());
  
  if ( cent > -1 ) 
  {
    centrality = CentralityName(cent);
  }
  
//  AliDebug(1,Form("fcent=%7.2f cent=%7.2f %s ndimu=%d",fcent,cent,centrality.Data(),di->GetEntries()));
  
//  AliInfo(Form("fcent=%7.2f cent=%7.2f %s nmu=%d ispileupspd=%d",
//               fcent,cent,centrality.Data(),header->GetNumberOfMuons(),
//               aod->IsPileupFromSPD(3,0.8)));

//  AliInfo(Form("nmu=%d run number=%d",header->GetNumberOfMuons(),fCurrentRunNumber));
  
  fVertex = aod->GetPrimaryVertex();
  
  TString eventtype(fGlobalEventSelectionName);
  eventtype.ToLower();
  Bool_t validTrigFound = kFALSE;
  
  while ( ( tname = static_cast<TObjString*>(next()) ) )
  {
//    AliInfo(Form("tname=%s",tname->String().Data()));
    
 //   TObjArray* tpart = tname->String().Tokenize("-");
//    
//    TString p1 = (static_cast<TObjString*>(tpart->At(0)))->String();
//    TString p2 = (static_cast<TObjString*>(tpart->At(1)))->String();
//                  
//    TString triggerShortName(Form("%s-%s",p1.Data(),p2.Data()));
//    
//    AliInfo(Form("tname=%s triggerShortName=%s",tname->GetName(),triggerShortName.Data()));
//    delete tpart;
 
    if ( !fTriggerClasses->FindObject(tname->GetName()) ) 
    {
      continue;
      // skip the trigger classes we're not supposed to analyse
    }
    
    validTrigFound = kTRUE;
    
    fEventCounters->Count(Form("event:%s/trigger:%s/run:%d", eventtype.Data(), tname->GetName(), fCurrentRunNumber));
    
    AssertHistogramCollection(fGlobalEventSelectionName.Data(),tname->String().Data());
    
    FillHistos(fGlobalEventSelectionName.Data(),tname->String().Data(),centrality.Data(),*aod);    
    
//    Histo(GLOBALEVENTSELECTIONNAME,tname->String().Data(),centrality.Data(),"Tracklets")->Fill(trackletMult);    
  }
  
  if (validTrigFound) {
    AssertHistogramCollection(fGlobalEventSelectionName.Data(),"ANY");
    FillHistos(fGlobalEventSelectionName.Data(),"ANY",centrality.Data(),*aod);    
  }
  
  delete a;
}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMuFromAOD::TrackChi2(const AliAODTrack& track) const
{
  return track.Chi2perNDF() < 3.5;
}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMuFromAOD::TrackEtaCut(const AliAODTrack& track) const
{
  return track.Eta() < -2.5 && track.Eta() > -4;
}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMuFromAOD::TrackMatchCut(const AliAODTrack& track) const
{
  return track.GetMatchTrigger() > 0;
}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMuFromAOD::TrackMatchLowCut(const AliAODTrack& track) const
{
  return track.GetMatchTrigger() > 1;
}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMuFromAOD::TrackMatchHighCut(const AliAODTrack& track) const
{
  return track.GetMatchTrigger() > 2;
}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMuFromAOD::TrackRabsCut(const AliAODTrack& track) const
{
  /// Cut between 2 and 10 degrees
  
  Double_t angle = TMath::ATan(track.GetRAtAbsorberEnd()/AbsZEnd());
  
  return ( angle > Deg2() && angle < Deg10() );
}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMuFromAOD::TrackPtCut(const AliAODTrack& track) const
{
  return track.Pt() > 1.0;
}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMuFromAOD::TrackDCACut(const AliAODTrack& track) const
{
  if (!fVertex) return kFALSE;
  
  TLorentzVector p(track.Px(),track.Py(),track.Pz(),
                   TMath::Sqrt(MuonMass2()+track.P()*track.P()));

  Double_t theta = TMath::ATan(track.GetRAtAbsorberEnd()/AbsZEnd());
  
  Double_t meanPcorr = ( theta < ( 180. - 3. ) * TMath::DegToRad() ) ? 1.2 : 1.49;
  Double_t pTotMean = p.P() - meanPcorr;
  
  TVector3 eventVertex(fVertex->GetX(),fVertex->GetY(),fVertex->GetZ());
  
  TVector3 trackDcaAtVz(track.XAtDCA(), track.YAtDCA(), fVertex->GetZ());
  TVector3 meanDca(-0.46, -0.92, 0.); // LHC10h1
  TVector3 dcaAtVz = trackDcaAtVz - eventVertex - meanDca;
  Double_t correctedDca = dcaAtVz.Mag(); // it should also be equal to dcaAtVz.Pt().
  Double_t cutVariable = pTotMean * correctedDca;
  
  Double_t cutValue = (theta > Deg3()) ? 63. : 120.;
  cutValue = TMath::Sqrt(cutValue*cutValue + 0.4*0.4*p.P()*p.P());
  
  return ( cutVariable < 5*cutValue );
}

//_____________________________________________________________________________
UInt_t AliAnalysisMuMuFromAOD::GetTrackMask(const AliAODTrack& track) const
{
  UInt_t m(kAll);
  
  if ( TrackRabsCut(track) ) m |= kRabs;
  
  if ( TrackPtCut(track) ) m |= kPt;  
  
  if ( TrackMatchCut(track) ) m |= kMatched;

  if ( TrackMatchLowCut(track) ) m |= kMatchedLow;

  if ( TrackMatchHighCut(track) ) m |= kMatchedHigh;

  if ( TrackEtaCut(track) ) m |= kEta;
  
  if ( TrackChi2(track) ) m |= kChi2;
  
  if ( TrackDCACut(track) ) m |= kDCA;
  
  return m;
}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMuFromAOD::PairRapidityCut(const AliAODTrack& t1, const AliAODTrack& t2) const
{
  TLorentzVector p1(t1.Px(),t1.Py(),t1.Pz(),TMath::Sqrt(MuonMass2()+t1.P()*t1.P()));
  TLorentzVector p2(t2.Px(),t2.Py(),t2.Pz(),TMath::Sqrt(MuonMass2()+t2.P()*t2.P()));
  
  TLorentzVector total(p1+p2);
  
  Double_t y = total.Rapidity();
  
  Bool_t ok = ( y < -2.5 && y > -4.0 );
  
  return ok;
}

//_____________________________________________________________________________
void AliAnalysisMuMuFromAOD::GetPairMask(const AliAODTrack& t1, const AliAODTrack& t2,
                                         UInt_t& mask1, UInt_t& mask2,
                                         UInt_t& mask12) const
{
  mask1 = GetTrackMask(t1);
  mask2 = GetTrackMask(t2);
    
  mask12 = mask1 | mask2;
  
  if ( PairRapidityCut(t1,t2) ) mask12 |= kPairRapidity;
}

//_____________________________________________________________________________
void AliAnalysisMuMuFromAOD::FillHistosForTrack(const char* physics,
                                                const char* triggerClassName, 
                                                const char* centrality,
                                                const AliAODTrack& track)
{
  TLorentzVector p(track.Px(),track.Py(),track.Pz(),
                   TMath::Sqrt(MuonMass2()+track.P()*track.P()));
  
  
  TString charge("Plus");
  
  if ( track.Charge() < 0 ) charge = "Minus";
    
  UInt_t mask = GetTrackMask(track);
  
//  AliDebug(1,Form("mask=%x pt=%e rabs=%e matched=%d",mask,track.Pt(),track.GetRAtAbsorberEnd(),track.GetMatchTrigger()));
  
  Double_t xdca = track.XAtDCA();
  Double_t ydca = track.YAtDCA();
  Double_t dca = TMath::Sqrt(xdca*xdca+ydca*ydca);
  Double_t theta = TMath::ATan(track.GetRAtAbsorberEnd()/AbsZEnd());
  
  Double_t meanPcorr = ( theta < ( 180. - 3. ) * TMath::DegToRad() ) ? 1.2 : 1.49;
  Double_t pTotMean = p.P() - meanPcorr;
  TVector3 eventVertex(fVertex->GetX(),fVertex->GetY(),fVertex->GetZ());
  TVector3 trackDcaAtVz(track.XAtDCA(), track.YAtDCA(), fVertex->GetZ());
  TVector3 meanDca(-0.46, -0.92, 0.); // LHC10h1
  TVector3 dcaAtVz = trackDcaAtVz - eventVertex - meanDca;
  Double_t correctedDca = dcaAtVz.Mag(); // it should also be equal to dcaAtVz.Pt().
  Double_t correctedPDca = pTotMean * correctedDca;
  
  Int_t nChHit = 0;
  for (Int_t i=0; i<10; i++) if (track.HitsMuonChamber(i)) nChHit++;
  
  TIter next(fSingleTrackCutNames);
  TObjString* str;
  
  while ( ( str = static_cast<TObjString*>(next()) ) )
  {
    Bool_t test = ( ( str->GetUniqueID() & mask ) == str->GetUniqueID() );
//    AliDebug(1,Form("  str %s uid %x mask %x test %d",str->String().Data(),str->GetUniqueID(),mask,test));
    
    if ( test ) 
    {
      Histo(physics,triggerClassName,centrality,str->String().Data(),Form("PtEtaMu%s",charge.Data()))->Fill(p.Eta(),p.Pt());
    
      Histo(physics,triggerClassName,centrality,str->String().Data(),"Pt")->Fill(p.Pt());
      
      Histo(physics,triggerClassName,centrality,str->String().Data(),Form("PEtaMu%s",charge.Data()))->Fill(p.Eta(),p.P());
    
    
      Histo(physics,triggerClassName,centrality,str->String().Data(),Form("XYdcaMu%s",charge.Data()))->Fill(ydca,xdca);
      
      Histo(physics,triggerClassName,centrality,str->String().Data(),Form("Chi2Mu%s",charge.Data()))->Fill(track.Chi2perNDF());
      
      Histo(physics,triggerClassName,centrality,str->String().Data(),"Chi2")->Fill(track.Chi2perNDF());
      
      Histo(physics,triggerClassName,centrality,str->String().Data(),"NChHit")->Fill(nChHit);
      
      if (p.Pt() > 2) {
	Histo(physics,triggerClassName,centrality,str->String().Data(),"Chi2Pt2GeV")->Fill(track.Chi2perNDF());
	
	Histo(physics,triggerClassName,centrality,str->String().Data(),"NChHitPt2GeV")->Fill(nChHit);
      }
      
      if ( theta >= Deg2() && theta < Deg3() )         
      {
//        Histo(physics,triggerClassName,centrality,str->String().Data(),Form("PDCAP23Mu%s",charge.Data()))->Fill(p.P(),p.P()*dca);
        Histo(physics,triggerClassName,centrality,str->String().Data(),Form("dcaP23Mu%s",charge.Data()))->Fill(p.P(),dca);
        
        Histo(physics,triggerClassName,centrality,str->String().Data(),"pdca23")->Fill(correctedPDca);
	
	if ( p.Pt() > 2 )
        {
          Histo(physics,triggerClassName,centrality,str->String().Data(),Form("dcaPwPtCut23Mu%s",charge.Data()))->Fill(p.P(),dca);
	  
	  Histo(physics,triggerClassName,centrality,str->String().Data(),"pdca23Pt2GeV")->Fill(correctedPDca);
        }
	
      }
      else if ( theta >= Deg3() && theta < Deg10() )
      {
//        Histo(physics,triggerClassName,centrality,str->String().Data(),Form("PDCAP310Mu%s",charge.Data()))->Fill(p.P(),p.P()*dca);        
        Histo(physics,triggerClassName,centrality,str->String().Data(),Form("dcaP310Mu%s",charge.Data()))->Fill(p.P(),dca);
        
        Histo(physics,triggerClassName,centrality,str->String().Data(),"pdca310")->Fill(correctedPDca);
	
	if ( p.Pt() > 2 )
        {
          Histo(physics,triggerClassName,centrality,str->String().Data(),Form("dcaPwPtCut310Mu%s",charge.Data()))->Fill(p.P(),dca);
	  
	  Histo(physics,triggerClassName,centrality,str->String().Data(),"pdca310Pt2GeV")->Fill(correctedPDca);
        }
      }
    }
  }
  
  
}

//_____________________________________________________________________________
void AliAnalysisMuMuFromAOD::FillHistos(const char* physics, const char* triggerClassName, 
                                        const char* centrality, const AliAODEvent& aod)
{
  if (!aod.GetNumberOfTracks()) 
  {
    return;
  }

  Histo(physics,triggerClassName,centrality,"Nevents")->Fill(1.0);

  if ( fVertex ) 
  {
    Histo(physics,triggerClassName,centrality,"Xvertex")->Fill(fVertex->GetX());
    Histo(physics,triggerClassName,centrality,"Yvertex")->Fill(fVertex->GetY());
    Histo(physics,triggerClassName,centrality,"Zvertex")->Fill(fVertex->GetZ());
  }
  
  // Track loop
  
    for (Int_t i = 0; i < aod.GetNumberOfTracks(); ++i) 
  {
    AliAODTrack* tracki = aod.GetTrack(i);

    if (!tracki->IsMuonTrack()) continue;

    FillHistosForTrack(physics,triggerClassName,centrality,*tracki);
    
    TLorentzVector pi(tracki->Px(),tracki->Py(),tracki->Pz(),
                      TMath::Sqrt(MuonMass2()+tracki->P()*tracki->P()));
    
    for (Int_t j = i+1; j < aod.GetNumberOfTracks(); ++j) 
    {
      AliAODTrack* trackj = aod.GetTrack(j);
      
      if (!trackj->IsMuonTrack()) continue;
      
      TLorentzVector pj(trackj->Px(),trackj->Py(),trackj->Pz(),
                       TMath::Sqrt(MuonMass2()+trackj->P()*trackj->P()));
      
      pj += pi;

      TIter next(fPairTrackCutNames);
      TObjString* str;
      
      UInt_t maski(0),maskj(0),maskij(0);
      
      GetPairMask(*trackj,*tracki,maski,maskj,maskij);
      
      while ( ( str = static_cast<TObjString*>(next()) ) )
      {
        UInt_t singleTrackMask(0);
        UInt_t pairMask(0);
        
        DecodePairCutMask(str->GetUniqueID(),singleTrackMask,pairMask);

        Bool_t testi = ( ( maski & singleTrackMask ) == singleTrackMask ) ;
        Bool_t testj = ( ( maskj & singleTrackMask ) == singleTrackMask ) ;
        Bool_t testij(kTRUE);
        
        if (pairMask>0) testij = ( ( maskij & pairMask ) == pairMask ) ;
        
        if ( testi && testj && testij )
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
  } //track loop
}

//==============================================================================
//
//   AddTask ....
//
//==============================================================================

#include "AliAnalysisMuMuFromAOD.h"
#include "AliAnalysisManager.h"
#include "TString.h"
#include "AliAnalysisDataContainer.h"
#include "AliVEventHandler.h"

AliAnalysisTask* AddTaskMuMuFromAOD(const char* outputname, TList* triggerClassesToConsider)
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
  
  AliAnalysisMuMuFromAOD* task = new AliAnalysisMuMuFromAOD(triggerClassesToConsider);
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
