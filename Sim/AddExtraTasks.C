Bool_t AddExtraTasks()
{
  /// Config macro to add tasks to e.g. AOD train
  
  // check input files
  if (gSystem->AccessPathName("AddTaskMuonPhysics.C") ||
      gSystem->AccessPathName("AliAnalysisTaskMuonPhysics.cxx") ||
      gSystem->AccessPathName("AliAnalysisTaskMuonPhysics.h")) {
    Error("AddExtraTasks","input files missing");
    return kFALSE;
  }
  
  // track selection
  AliMuonTrackCuts *trackCuts = 0x0;
  trackCuts = new AliMuonTrackCuts("stdCuts", "stdCuts");
  trackCuts->SetAllowDefaultParams();
  trackCuts->SetFilterMask(AliMuonTrackCuts::kMuMatchLpt | AliMuonTrackCuts::kMuEta | AliMuonTrackCuts::kMuThetaAbs | AliMuonTrackCuts::kMuPdca);
  trackCuts->SetIsMC(kTRUE);
  
  // add AliAnalysisTaskMuonPhysics
  gROOT->LoadMacro("AliAnalysisTaskMuonPhysics.cxx+");
  gROOT->LoadMacro("AddTaskMuonPhysics.C");
  AliAnalysisTaskMuonPhysics* physics = AddTaskMuonPhysics("phys");
  if (!physics) return kFALSE;
  physics->SetMuonTrackCuts(*trackCuts);
  physics->SetMuonPtCut(0.5);
  physics->UseMCLabel(kTRUE);
  physics->VersusRun(kTRUE, kTRUE, kFALSE);
  
  return kTRUE;

}
