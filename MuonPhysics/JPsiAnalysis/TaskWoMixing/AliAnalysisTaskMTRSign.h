#ifndef ALIANALYSISTASKMTRSIGN_H
#define ALIANALYSISTASKMTRSIGN_H

//
// AliAnalysisTaskMTRSign
// Analysis task to study trigger sign
//
//  Author: Diego Stocco
//

#include "AliAnalysisTaskSE.h"

class TString;
class AliCFGridSparse;
class AliMuonEventCuts;
class AliMuonTrackCuts;
class AliMergeableCollection;
class AliVParticle;


class AliAnalysisTaskMTRSign : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskMTRSign();
  AliAnalysisTaskMTRSign(const char *name);
  virtual ~AliAnalysisTaskMTRSign();

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   NotifyRun();
  virtual void   Terminate(Option_t *option);
  virtual void   FinishTaskOutput();
  
  void UseMCLabel(Bool_t flag = kTRUE) { fUseMCLabel = flag; }
  
  /// Get muon event cuts
  AliMuonEventCuts* GetMuonEventCuts() { return fMuonEventCuts; }
  /// Get muon track cuts
  AliMuonTrackCuts* GetMuonTrackCuts() { return fMuonTrackCuts; }
  
  Int_t TriggerDevSign ( AliVParticle* track ) const; // REMEMBER TO CUT when commit done
  
  enum {
    kHvarPt,         ///< Pt at vertex
    kHvarP,          ///< Total momentum
    kHvarCharge,     ///< Particle charge (from tracker)
    kHvarMatchTrig,  ///< Tracking trigger matching
    kHvarTrigSign,   ///< Deviation sign from trigger
    kHvarLoCircuit,  ///< Local board id
    kNvars           ///< THnSparse dimensions
  };

 private:

  AliAnalysisTaskMTRSign(const AliAnalysisTaskMTRSign&);
  AliAnalysisTaskMTRSign& operator=(const AliAnalysisTaskMTRSign&);
  
  AliCFGridSparse* GetCFGridSparse ( TString identifier );
  
  AliMuonEventCuts* fMuonEventCuts; ///< Muon event cuts
  AliMuonTrackCuts* fMuonTrackCuts; ///< Muon track cuts
  Bool_t            fUseMCLabel;    ///< Select tracks using MC label
  AliMergeableCollection* fMergeableCollection; //!< collection of mergeable objects

  ClassDef(AliAnalysisTaskMTRSign, 1); // MTR sign analysis
};

#endif
