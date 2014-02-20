#ifndef AliAnalysisMuMuFromESD_H
#define AliAnalysisMuMuFromESD_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup
/// \class AliAnalysisMuMuFromESD
/// \brief
/// 
/// \author Laurent Aphecetche

#ifndef ALIANALYSISMUMU_H
#  include "AliAnalysisMuMu.h"
#endif
#ifndef ALIESDVZERO_H
#  include "AliESDVZERO.h"
#endif

class AliESDEvent;
class TH1;
class TH2;
class AliESDMuonTrack;

class AliAnalysisMuMuFromESD : public AliAnalysisMuMu
{
public:
  AliAnalysisMuMuFromESD();
  AliAnalysisMuMuFromESD(TList* triggerClassesToConsider);
  virtual ~AliAnalysisMuMuFromESD();

protected:
  virtual void MuUserExec(Option_t *option);

private:
  
  UInt_t GetTrackMask(const AliESDMuonTrack& track) const;

  void FillHistogramCollection(const char* physics, const char* triggerClassName);
  
  void FillHistosForTrack(const char* physics, const char* triggerClassName, const char* centrality, const AliESDMuonTrack& track, const char* runNumber);
  
  void FillHistos(const char* physics, const char* triggerClassName, const char* centrality, const AliESDEvent& esd);

private:
  const char* RunNumber(const AliESDEvent& esd) const;
  Bool_t TrackMatchCut(const AliESDMuonTrack& track) const;
  Bool_t TrackMatchLowCut(const AliESDMuonTrack& track) const;
  Bool_t TrackMatchHighCut(const AliESDMuonTrack& track) const;
  Bool_t TrackRabsCut(const AliESDMuonTrack& track) const;
  Bool_t TrackPtCut(const AliESDMuonTrack& track) const;
  Bool_t TrackChi2(const AliESDMuonTrack& track) const;
  Bool_t TrackEtaCut(const AliESDMuonTrack& track) const;
  
  ClassDef(AliAnalysisMuMuFromESD,4) // 
};


AliAnalysisTask* AddTaskMuMuFromESD(const char* outputname="test.MuMu.ESD.1.root", TList* triggerClassesToConsider=0x0);

#endif
