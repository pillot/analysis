#ifndef ALIANALYSISTASKJPSIACCEFFCORR_H
#define ALIANALYSISTASKJPSIACCEFFCORR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/// \ingroup muondep
/// \class AliAnalysisTaskJPsiAccEffCorr
/// \brief task to extrac the JPsi acc*eff corrections
//Author: Philippe Pillot - SUBATECH Nantes

#include "TArrayF.h"
#include "AliAnalysisTaskSE.h"
#include "AliLog.h"

class TH1F;
class TObjArray;
class AliCounterCollection;

class AliAnalysisTaskJPsiAccEffCorr : public AliAnalysisTaskSE {
public:
  
  AliAnalysisTaskJPsiAccEffCorr();
  AliAnalysisTaskJPsiAccEffCorr(const char *name);
  virtual ~AliAnalysisTaskJPsiAccEffCorr();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *);
  virtual void   Terminate(Option_t *);
  
  /// set centrality binning
  void SetCentBins(Int_t nBins, Float_t *binLowEdge) {fCentBinLowEdge.Set(nBins+1, binLowEdge);}
  
  /// set pt binning
  void SetPtBins(Int_t nBins, Float_t *binLowEdge) {fPtBinLowEdge.Set(nBins+1, binLowEdge);}
  
  /// set y binning
  void SetYBins(Int_t nBins, Float_t *binLowEdge) {fYBinLowEdge.Set(nBins+1, binLowEdge);}
  
  /// set the trigger level to be matched with (0=no, 1=all, 2=low, 3=high)
  void SetTrigLevel(Int_t level = 1) {fTrigLevel = (level>=0 && level<4) ? level : 0;}
  
  /// set the muon pt cut value
  void SetMuLowPtCut(Double_t cut) {fMuLowPtCut = cut;}
  
private:
  
  /// Not implemented
  AliAnalysisTaskJPsiAccEffCorr(const AliAnalysisTaskJPsiAccEffCorr& rhs);
  /// Not implemented
  AliAnalysisTaskJPsiAccEffCorr& operator = (const AliAnalysisTaskJPsiAccEffCorr& rhs);
  
  /// Compute acc*eff and binomial errors by hand, i.e. not using TGraphAsymmErrors
  TH1F* ComputeAccEff(TH1F &hGen, TH1F &hRec, const Char_t *name, const Char_t *title);
  
  // compute AccEff correction integrated over the given centrality range
  void IntegratedAccEff(Int_t ipt, Int_t iy, Int_t nMatch, Float_t centMin, Float_t centMax,
			Double_t &accEff, Double_t &accEffErr, Bool_t print = kFALSE,
			TH1F* hGenSum = 0x0, TH1F* hRecSum = 0x0, TH1F* hAccSum = 0x0, Int_t bin = -1);
  
  // compute integrated AccEff correction weighted by the number of JPsi in each centrality bin
  void IntegratedAccEff(Int_t ipt, Int_t iy, Int_t nMatch, Int_t nBins, Float_t *nJPsi, Float_t *centbinLowEdge,
			Double_t &accEff, Double_t &accEffErr, Bool_t print = kFALSE);
  
  // Draw acceptance*efficiency versus run for this given pt/y bin
  void DrawAccEffVsRun(Int_t ipt, Int_t iy);
  
  // Draw acceptance*efficiency versus centrality for this given pt/y bin
  void DrawAccEffVsCent(Int_t ipt, Int_t iy);
  
private:
  
  enum eList {
    kPtGen   = 0, ///< pT distribution of generated JPsi
    kPtRec   = 1, ///< pT distribution of reconstructed JPsi
    kYGen    = 2, ///< y distribution of generated JPsi
    kYRec    = 3, ///< y distribution of reconstructed JPsi
    kPtGenMu = 4, ///< pT distribution of generated muon
    kPtRecMu = 5, ///< pT distribution of reconstructed muon
    kYGenMu  = 6, ///< y distribution of generated muon
    kYRecMu  = 7, ///< y distribution of reconstructed muon
    kDzVtx   = 8, ///< vertex resolution
    kDzVtx2  = 9, ///< vertex resolution for event with JPsi
    kMass    = 10 ///< invariant mass distribution
  };
  
  TObjArray*  fList; //!< List of output object
  AliCounterCollection* fEventCounters; //!< event statistics
  AliCounterCollection* fJPsiCounters;  //!< JPsi statistics
  
  TArrayF  fCentBinLowEdge; ///< centrality bin low-edge values
  TArrayF  fPtBinLowEdge;   ///< pT bin low-edge values
  TArrayF  fYBinLowEdge;    ///< y bin low-edge values
  Int_t    fTrigLevel;      ///< trigger level to be matched with (1=all, 2=low, 3=high)
  Double_t fMuLowPtCut;     ///< muon pt cut value
  
  ClassDef(AliAnalysisTaskJPsiAccEffCorr, 1);
};

#endif

