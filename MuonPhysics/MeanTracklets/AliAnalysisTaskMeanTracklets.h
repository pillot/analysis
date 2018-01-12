#ifndef ALIANALYSISTASKMEANTRACKLETS_H
#define ALIANALYSISTASKMEANTRACKLETS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/// \ingroup muon
/// \class AliAnalysisTaskMeanTracklets
/// \brief task to study the mean number of SPD tracklets and correct for the z dependence
//Author: Philippe Pillot - SUBATECH Nantes

#include "TProfile.h"
#include "AliAnalysisTaskSE.h"

class THnSparse;
class TRandom3;
class AliCounterCollection;


class AliAnalysisTaskMeanTracklets : public AliAnalysisTaskSE {
public:
  
  AliAnalysisTaskMeanTracklets();
  AliAnalysisTaskMeanTracklets(const char *name);
  virtual ~AliAnalysisTaskMeanTracklets();
  
  virtual void   NotifyRun(){}
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *);
  virtual void   Terminate(Option_t *){}
  
  void SetMeanNtrkVsZvtxRef(TProfile &meanNtrkVsZvtxRef, Double_t meanNtrkRef = -1.);
  
  /// use a binomial instead of poissonian distribution to correct Ntrk when fMeanNtrkRef < <Ntrk>
  void UseBinomial(Bool_t flag = kTRUE) {fUseBinomial = flag;}
  
  /// select a specific trigger class
  void SelectTrigger(TString trigger) {fTrigger = trigger;}
  
  /// reject NSD MC events
  void RejectNSD() {fRejectNSD = kTRUE;}
  
  /// reject pile-up from SPD
  void RejectPUFromSPD() {fRejectPUFromSPD = kTRUE;}
  
  /// disable SPD vertex QA selection
  void DisableSPDVtxQA() {fSelectSPDVtxQA = kFALSE;}
  
  /// reject events with 0 (corrected) tracklets in -1 < eta < 1
  void Reject0Tracklet() {fReject0Tracklet = kTRUE;}
  
private:
  
  /// Not implemented
  AliAnalysisTaskMeanTracklets(const AliAnalysisTaskMeanTracklets& rhs);
  /// Not implemented
  AliAnalysisTaskMeanTracklets& operator = (const AliAnalysisTaskMeanTracklets& rhs);
  
  Int_t GetCorrectedNtrk(Int_t NtrkInEtaRange, Double_t zVtx);
  Int_t GetCorrectedNtrkFromMultSel();
  
private:
  
  AliCounterCollection* fEvents;  //!< number of analyzed events
  THnSparse *fhNtrk;              //!< output histogram for number of SPD tracklets
  THnSparse *fhNtrkCorr;          //!< output histogram for corrected number of SPD tracklets
  TProfile *fpMeanNtrkVsZvtx;     //!< <Ntrk> vs Z profile before correction
  TProfile *fpMeanNtrkVsZvtxCorr; //!< <Ntrk> vs Z profile after correction
  
  TProfile *fpMeanNtrkVsZvtxRef;  /// <Ntrk> vs Z profile used to correct Ntrk
  Double_t fMeanNtrkRef;          /// <Ntrk> value used as a reference
  Bool_t   fUseBinomial;          /// use a binomial distribution to correct Ntrk when fMeanNtrkRef < <Ntrk>
  TRandom3 *fRandom;              //!< random number generator
  
  TString fTrigger;               /// select a specific trigger class
  Bool_t fRejectNSD;              /// reject NSD MC events
  Bool_t fRejectPUFromSPD;        /// reject pile-up from SPD
  Bool_t fSelectSPDVtxQA;         /// select events with a good SPD vertex
  Bool_t fReject0Tracklet;        /// reject events with 0 (corrected) tracklets in -1 < eta < 1
  
  ClassDef(AliAnalysisTaskMeanTracklets, 2);
};

//________________________________________________________________________
inline void AliAnalysisTaskMeanTracklets::SetMeanNtrkVsZvtxRef(TProfile &meanNtrkVsZvtxRef, Double_t meanNtrkRef)
{
  /// set the <Ntrk> vs Z profile used to correct Ntrk
  /// set the reference <Ntrk> to the given value or to the minimum if not provided
  delete fpMeanNtrkVsZvtxRef;
  fpMeanNtrkVsZvtxRef = new TProfile(meanNtrkVsZvtxRef);
  fpMeanNtrkVsZvtxRef->SetDirectory(0);
  fMeanNtrkRef = (meanNtrkRef > 0.) ? meanNtrkRef : fpMeanNtrkVsZvtxRef->GetMinimum();
  if (fMeanNtrkRef < 1.e-6) AliFatal(Form("Invalid <Ntrk> reference value: %f",fMeanNtrkRef));
}

#endif

