/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#ifndef ALIANALYSISTASKEX02_H
#define ALIANALYSISTASKEX02_H

#ifndef ALIANALYSISTASKSE_H
#include "AliAnalysisTaskSE.h"
#endif

class TH1F;
class TArrayF;
class TList;
class AliCounterCollection;
class AliMuonTrackCuts;


class AliAnalysisTaskEx02 : public AliAnalysisTaskSE {
  
public:
  AliAnalysisTaskEx02();
  AliAnalysisTaskEx02(const char *name);
  virtual ~AliAnalysisTaskEx02();
  
  virtual void     NotifyRun();
  virtual void     UserCreateOutputObjects();
  virtual void     UserExec(Option_t *option);
  virtual void     UserExecMix(Option_t *);
  
  void SetMuonTrackCuts(AliMuonTrackCuts* trackCuts) { fMuonTrackCuts = static_cast<AliMuonTrackCuts*>(trackCuts->Clone()); }
  void SetPtLowEdge(Int_t nPtBins, Float_t *ptLowEdge) { fPtBin[0].Set(nPtBins, ptLowEdge); }
  void SetPtUpEdge(Int_t nPtBins, Float_t *ptUpEdge) { fPtBin[1].Set(nPtBins, ptUpEdge); }
  void SetYLowEdge(Int_t nYBins, Float_t *yLowEdge) { fYBin[0].Set(nYBins, yLowEdge); }
  void SetYUpEdge(Int_t nYBins, Float_t *yUpEdge) { fYBin[1].Set(nYBins, yUpEdge); }
  
protected:
  Double_t MuonMass2() const;
  
private:
  Bool_t TrackPtCut(const AliAODTrack& track) const;
  Bool_t PairRapidityCut(const TLorentzVector& pair) const;
  
  TList                  *fOutput;            // Output list
  TH1F                   *fHistCent;		  // Check centrality
  AliCounterCollection   *fEventCounters;     //! event counters
  AliMuonTrackCuts       *fMuonTrackCuts;     // track cuts
  TArrayF				 fPtBin[2];
  TArrayF				 fYBin[2];
  
  AliAnalysisTaskEx02(const AliAnalysisTaskEx02 &); // not implemented
  AliAnalysisTaskEx02 &operator=(const AliAnalysisTaskEx02 &); // not implemented
  
  ClassDef(AliAnalysisTaskEx02, 1);
};

#endif
