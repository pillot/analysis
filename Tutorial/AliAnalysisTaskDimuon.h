#ifndef ALIANALYSISTASKDIMUON_H
#define ALIANALYSISTASKDIMUON_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/// \ingroup muon
/// \class AliAnalysisTaskDimuon
/// \brief basic task to analyse dimuons
//Author: Philippe Pillot - SUBATECH Nantes

#include "AliAnalysisTaskSE.h"
#include "AliMuonTrackCuts.h"

class THnSparse;
class AliCounterCollection;

class AliAnalysisTaskDimuon : public AliAnalysisTaskSE {
public:
  
  AliAnalysisTaskDimuon();
  AliAnalysisTaskDimuon(const char *name);
  virtual ~AliAnalysisTaskDimuon();
  
  virtual void   NotifyRun();
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *);
  virtual void   Terminate(Option_t *);
  
  void SetMuonTrackCuts(AliMuonTrackCuts &trackCuts);
  void SetDefaultMuonTrackCuts(Bool_t isMC);
  
private:
  
  /// Not implemented
  AliAnalysisTaskDimuon(const AliAnalysisTaskDimuon& rhs);
  /// Not implemented
  AliAnalysisTaskDimuon& operator = (const AliAnalysisTaskDimuon& rhs);
  
private:
  
  AliCounterCollection* fEvents; //!< number of analyzed events
  THnSparse *fhOS; //!< output histogram for opposite-sign dimuons
  THnSparse *fhLS; //!< output histogram for Like-sign dimuons
  
  AliMuonTrackCuts *fMuonTrackCuts; ///< cuts to select tracks to be considered
  
  ClassDef(AliAnalysisTaskDimuon, 1);
};

//________________________________________________________________________
inline void AliAnalysisTaskDimuon::SetMuonTrackCuts(AliMuonTrackCuts &trackCuts)
{
  /// set cuts to select tracks to be considered
  delete fMuonTrackCuts;
  fMuonTrackCuts = new AliMuonTrackCuts(trackCuts);
}

//________________________________________________________________________
inline void AliAnalysisTaskDimuon::SetDefaultMuonTrackCuts(Bool_t isMC)
{
  /// set default cuts to select tracks to be considered
  delete fMuonTrackCuts;
  fMuonTrackCuts = new AliMuonTrackCuts("stdCuts", "stdCuts");
  fMuonTrackCuts->SetAllowDefaultParams();
  fMuonTrackCuts->SetFilterMask(AliMuonTrackCuts::kMuMatchLpt | AliMuonTrackCuts::kMuEta |
                                AliMuonTrackCuts::kMuThetaAbs | AliMuonTrackCuts::kMuPdca);
  fMuonTrackCuts->SetIsMC(isMC);
}

#endif

