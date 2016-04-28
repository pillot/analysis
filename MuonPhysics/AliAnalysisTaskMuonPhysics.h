#ifndef ALIANALYSISTASKMUONPHYSICS_H
#define ALIANALYSISTASKMUONPHYSICS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/// \ingroup muondep
/// \class AliAnalysisTaskMuonPhysics
/// \brief task to extract physical quantities
//Author: Philippe Pillot - SUBATECH Nantes

#include "AliAnalysisTaskSE.h"
#include "AliMuonTrackCuts.h"

class TObjArray;
class AliVVertex;
class AliVParticle;
class AliCounterCollection;

class AliAnalysisTaskMuonPhysics : public AliAnalysisTaskSE {
public:
  
  AliAnalysisTaskMuonPhysics();
  AliAnalysisTaskMuonPhysics(const char *name);
  virtual ~AliAnalysisTaskMuonPhysics();
  
  virtual void   UserCreateOutputObjects();
  virtual void   NotifyRun();
  virtual void   UserExec(Option_t *);
  virtual void   Terminate(Option_t *);
  
  // set standard cuts to select tracks to be considered
  void SetMuonTrackCuts(AliMuonTrackCuts &trackCuts);
  
  // set standard cuts to select tracks not to be considered
  void SetMuonTrackCuts2(AliMuonTrackCuts &trackCuts);
  
  /// set the muon low pT cut to select tracks to be considered
  void SetMuonPtCut(Double_t cut) {fPtCut = cut;}
  
  /// set the flag to select tracks using MC label
  void UseMCLabel(Bool_t flag = kTRUE) { fUseMCLabel = flag; }
  
  /// Select tracks marked as bad (or no longer matched) by the reffing task
  void SelectBadTracks(Bool_t flag = kTRUE) {fSelectBadTracks = flag;}
  
  /// Select tracks in the given centrality range
  void SelectCentrality(Double_t min, Double_t max) {fCentMin = min; fCentMax = max;}
  
  /// Fill counters versus run (the size of the counters may explode!)
  void VersusRun(Bool_t ev, Bool_t trk, Bool_t trg) {fEvVsRun = ev; fTrkVsRun = trk; fTrgVsRun = trg;}
  
private:
  
  /// Not implemented
  AliAnalysisTaskMuonPhysics(const AliAnalysisTaskMuonPhysics& rhs);
  /// Not implemented
  AliAnalysisTaskMuonPhysics& operator = (const AliAnalysisTaskMuonPhysics& rhs);
  
  Bool_t IsSelected(AliVParticle& track, Bool_t isESD, Bool_t fillCounters, TString &centKey, Bool_t selectBadTracks);
  
  Bool_t GetVtxStatus(const AliVVertex &vtx) const;
  
private:
  
  enum eList {
    kPt                  = 0,  ///< Pt distribution of single muons
    kRapidity            = 1,  ///< rapidity distribution of single muons
    kDCA                 = 2,  ///< DCA distribution of single muons
    kChi2                = 3,  ///< normalized chi2 distribution of single muons
    kNClustersPerTrack   = 4,  ///< number of clusters per track
    kNChamberHitPerTrack = 5,  ///< number of chamber hit per track
    kMass                = 6,  ///< invariant mass distribution of opposite sign dimuons
    kPUncorrected        = 7,  ///< uncorrected momentum distribution of single muons
    kRAbs                = 8,  ///< track position at the end of the absorber
    kDCAX                = 9,  ///< DCA distribution of single muons in non-bending direction
    kDCAY                = 10, ///< DCA distribution of single muons in non-bending direction
    kNTracks             = 11, ///< number of tracks
    kNTrigAll            = 12, ///< number of trigger tracks
    kNTrigLow            = 13, ///< number of low-pt trigger tracks
    kNTrigHigh           = 14, ///< number of high-pt trigger tracks
    kNMatchTracks        = 15, ///< number of tracks matched with trigger
    kChi2Trig            = 16, ///< matching chi2 distribution of single muons
    kPDCA23              = 17, ///< pDCA distribution of single muons in [2,3] deg
    kPDCA310             = 18, ///< pDCA distribution of single muons in ]3,10] deg
    kPtMuPlus            = 19, ///< mu+ Pt distribution of single muons
    kPtMuMinus           = 20, ///< mu- Pt distribution of single muons
    kSPDzVtx             = 21, ///< z position of the SPD vertex
    kT0zVtx              = 22, ///< z position of the T0 vertex
    kT0SPDDeltaZVtx      = 23, ///< delta_z of SPD-T0 vertices
    kT0VsSPDZVtx         = 24, ///< z position of T0 versus SPD vertices
    kPtDimu              = 25, ///< Pt distribution of dimuons
    kInvBendingP         = 26, ///< inverse bending momentum distribution of single muon
    kDCA23VsP            = 27, ///< DCA distribution of single muons versus p in [2,3] deg
    kDCA310VsP           = 28, ///< DCA distribution of single muons versus p in ]3,10] deg
    kCent                = 29, ///< centrality distribution
    kMult                = 30, ///< tracklet multiplicity
    kMultSelect          = 31, ///< tracklet multiplicity for events with selected tracks
    kPosX                = 32, ///< track x position at vertex
    kPosY                = 33, ///< track y position at vertex
    kPosZ                = 34  ///< track z position at vertex
  };
  
  TObjArray*  fList; //!< List of output object
  
  AliCounterCollection* fTrackCounters; //!< track statistics
  AliCounterCollection* fEventCounters; //!< event statistics
  AliCounterCollection* fTrigCounters;  //!< trigger track statistics
  
  AliMuonTrackCuts* fMuonTrackCuts; ///< cuts to select tracks to be considered
  AliMuonTrackCuts* fMuonTrackCuts2;///< cuts to select tracks not to be considered
  Double_t fPtCut;                  ///< minimum pT cut to select tracks to be considered
  Bool_t   fUseMCLabel;             ///< Select tracks using MC label
  Bool_t   fSelectBadTracks;        ///< Select bad (or no longer matched) tracks
  Double_t fCentMin;                ///< Select centrality > fCentMin
  Double_t fCentMax;                ///< Select centrality <= fCentMax
  Bool_t   fEvVsRun;                ///< Fill event counters versus run
  Bool_t   fTrkVsRun;               ///< Fill track counters versus run
  Bool_t   fTrgVsRun;               ///< Fill trigger counters versus run
  
  ClassDef(AliAnalysisTaskMuonPhysics, 3);
};

//________________________________________________________________________
inline void AliAnalysisTaskMuonPhysics::SetMuonTrackCuts(AliMuonTrackCuts &trackCuts)
{
  /// set standard cuts to select tracks to be considered
  delete fMuonTrackCuts;
  fMuonTrackCuts = new AliMuonTrackCuts(trackCuts);
}

//________________________________________________________________________
inline void AliAnalysisTaskMuonPhysics::SetMuonTrackCuts2(AliMuonTrackCuts &trackCuts)
{
  /// set standard cuts to select tracks to be considered
  delete fMuonTrackCuts2;
  fMuonTrackCuts2 = new AliMuonTrackCuts(trackCuts);
}

#endif

