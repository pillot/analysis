/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#ifndef ALIANALYSISTASKJPSI_H
#define ALIANALYSISTASKJPSI_H

#include <TNamed.h>
#include "AliAnalysisTaskSE.h"

class TH1F;
class TList;
class TArrayF;
class TString;
class TObjString;
class AliCounterCollection;
class AliMuonTrackCuts;


class MyList : public TNamed {
  
public:
  MyList(const char *name = "myList", const char *title = "my list");
  virtual ~MyList();
  
  virtual void Add(TObject* obj);
  virtual Long64_t Merge(TCollection* list);
  virtual void Print(Option_t* option = "") const;
  
  /// return the internal list
  virtual const TList* GetList() const {return fList;}
  
private:
  MyList(const MyList&);
  MyList &operator=(const MyList&);
  
  TList *fList; ///< internal list
  
  ClassDef(MyList, 1);
};


class AliAnalysisTaskJPsi : public AliAnalysisTaskSE {
  
public:
  AliAnalysisTaskJPsi();
  AliAnalysisTaskJPsi(const char *name);
  virtual ~AliAnalysisTaskJPsi();
  
  virtual void NotifyRun();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *);
  
  void SetMuonTrackCuts(AliMuonTrackCuts &trackCuts) { delete fMuonTrackCuts; fMuonTrackCuts = new AliMuonTrackCuts(trackCuts); }
  void SelectTrgSign(Bool_t flag = kTRUE) { fSelectTrgSign = flag; }
  void SelectSameTrgSignFake(Bool_t flag = kTRUE) { fSelectSameTrgSignFake = flag; }
  void UseMCLabel(Bool_t flag = kTRUE) { fUseMCLabel = flag; }
  void SetPtLowEdge(Int_t nPtBins, Float_t *ptLowEdge) { fPtBin[0].Set(nPtBins, ptLowEdge); }
  void SetPtUpEdge(Int_t nPtBins, Float_t *ptUpEdge) { fPtBin[1].Set(nPtBins, ptUpEdge); }
  void SetYLowEdge(Int_t nYBins, Float_t *yLowEdge) { fYBin[0].Set(nYBins, yLowEdge); }
  void SetYUpEdge(Int_t nYBins, Float_t *yUpEdge) { fYBin[1].Set(nYBins, yUpEdge); }
  void RecordEvWithTrgIssues(Bool_t flag = kTRUE) { fRecordEvWithTrgIssues = flag; }
  
protected:
  Double_t MuonMass2() const;
  
private:
  Bool_t PairRapidityCut(const TLorentzVector& pair) const;
  Int_t  TriggerDevSign(AliVParticle *track) const;
  TObjString* FindFile(TString &file, MyList *badFiles) const;
  
  TList                *fOutput;        //! Output list
  TH1F                 *fHistCent;      //! Check centrality
  AliCounterCollection *fEventCounters; //! event counters
  
  AliMuonTrackCuts     *fMuonTrackCuts;        // track cuts
  Bool_t               fTrackCutsSet;          //! set track cuts only once
  Bool_t               fSelectTrgSign;         // Select pairs according to the trigger sign (OS or LS)
  Bool_t               fSelectSameTrgSignFake; // Select pairs matching the same trigger track and producing LS & !OS
  Bool_t               fUseMCLabel;            // Select tracks using MC label
  TArrayF              fPtBin[2];              // pT bins
  TArrayF              fYBin[2];               // y bins
  
  Bool_t               fRecordEvWithTrgIssues;    // record events with trigger issues in the list below
  MyList               *fTrgClassMissTrgL0Ev;     //! events with l0 trigger input missing according to trigger class
  MyList               *fTrgL0MissTrgClassEv;     //! events with trigger class missing according to l0 trigger input
  MyList               *fTrgClassMissTrgOffEv;    //! events with trigger track pair deviation sign missing according to trigger class
  MyList               *fTrgOffMissTrgClassEv;    //! events with trigger class missing according to trigger track pair deviation sign
  MyList               *fTrgOffMissTrgL0Ev;       //! events with l0 trigger input missing according to trigger track pair deviation sign
  MyList               *fTrgL0MissTrgOffEv;       //! events with trigger track pair deviation sign missing according to l0 trigger input
  MyList               *fOSTrkOSTrgMLLOnlyEv;     //! events with selected OS pair(s) matching OS trigger track pair in MLL&!MUL class
  MyList               *fOSTrkLSTrgMLLOnlyEv;     //! events with selected OS pair(s) matching LS trigger track pair in MLL&!MUL class
  MyList               *fOSTrkLSTrgMULEv;         //! events with selected OS pair(s) matching LS trigger track pair in MUL class
  MyList               *fOSTrkOSTrgFakeMLLOnlyEv; //! events with selected OS pair(s) matching same trigger track (->?S) in MLL&!MUL class
  MyList               *fOSTrkLSTrgFakeMLLOnlyEv; //! events with selected OS pair(s) matching same trigger track (->LS) in MLL&!MUL class
  MyList               *fOSTrkLSTrgFakeMULEv;     //! events with selected OS pair(s) matching same trigger track (->LS) in MUL class
  
  AliAnalysisTaskJPsi(const AliAnalysisTaskJPsi &); // not implemented
  AliAnalysisTaskJPsi &operator=(const AliAnalysisTaskJPsi &); // not implemented
  
  ClassDef(AliAnalysisTaskJPsi, 1);
};

#endif
