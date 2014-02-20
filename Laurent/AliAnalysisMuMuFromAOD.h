#ifndef AliAnalysisMuMuFromAOD_H
#define AliAnalysisMuMuFromAOD_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup
/// \class AliAnalysisMuMuFromAOD
/// \brief
/// 
/// \author Laurent Aphecetche

#ifndef ALIANALYSISMUMU_H
#  include "AliAnalysisMuMu.h"
#endif

class TH1;
class AliAODTrack;
class AliAODEvent;
class AliAODVertex;

class AliAnalysisMuMuFromAOD : public AliAnalysisMuMu
{
public:
  AliAnalysisMuMuFromAOD();
  AliAnalysisMuMuFromAOD(TList* triggerClassesToConsider);
  virtual ~AliAnalysisMuMuFromAOD();

protected:
  virtual void MuUserExec(Option_t *option);
  void DumpMC(const AliAODEvent& aod);

private:
  
  UInt_t GetTrackMask(const AliAODTrack& track) const;
  
  void FillHistosForTrack(const char* physics, const char* triggerClassName, 
                          const char* centrality, const AliAODTrack& track);
  
  void FillHistos(const char* physics, const char* triggerClassName, const char* centrality, const AliAODEvent& event);

private:
  Bool_t TrackMatchCut(const AliAODTrack& track) const;
  Bool_t TrackMatchLowCut(const AliAODTrack& track) const;
  Bool_t TrackMatchHighCut(const AliAODTrack& track) const;
  Bool_t TrackRabsCut(const AliAODTrack& track) const;
  Bool_t TrackPtCut(const AliAODTrack& track) const;
  Bool_t TrackEtaCut(const AliAODTrack& track) const;
  Bool_t TrackChi2(const AliAODTrack& track) const;
  Bool_t TrackDCACut(const AliAODTrack& track) const;
  
  Bool_t PairRapidityCut(const AliAODTrack& t1, const AliAODTrack& t2) const;
  void GetPairMask(const AliAODTrack& t1, const AliAODTrack& t2, UInt_t& mask1, UInt_t& mask2, UInt_t& mask12) const;
  
private:
  
  AliAODVertex* fVertex; //! current event vertex
  TString fGlobalEventSelectionName; // global event selection name
  
  ClassDef(AliAnalysisMuMuFromAOD,7) // Analysis of mu-mu pairs from AOD
};


AliAnalysisTask* AddTaskMuMuFromAOD(const char* outputname, TList* triggerClassesToConsider=0x0);

#endif
