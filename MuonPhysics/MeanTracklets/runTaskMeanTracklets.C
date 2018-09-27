/*
 *  runTaskDimuon.C
 *  aliroot
 *
 *  Created by Philippe Pillot on 10/04/17.
 *  Copyright 2017 SUBATECH. All rights reserved.
 *
 */

//______________________________________________________________________________
void runTaskMeanTracklets(TString smode = "local", TString inputFileName = "AliAOD.Muons.root",
                          Bool_t applyPS = kTRUE, Bool_t applyPileupCuts = kTRUE,
                          Bool_t isMC = kTRUE, TString refInput = "Results_PFwPU.root")
{
  /// Run the baseline task to test the framework
  
  // --- general analysis setup ---
  TString rootVersion = "";
  TString alirootVersion = "";
  TString aliphysicsVersion = "vAN-20180907-1";
  TString extraLibs="";
  TString extraIncs="include";
  TString extraTasks="AliAnalysisTaskMeanTracklets";
  TString extraPkgs="";
  TList pathList; pathList.SetOwner();
  pathList.Add(new TObjString("$WORK/Macros/MuonPhysics/MeanTracklets"));
  TList fileList; fileList.SetOwner();
  fileList.Add(new TObjString("runTaskMeanTracklets.C"));
  fileList.Add(new TObjString("AddTaskMeanTracklets.C"));
  fileList.Add(new TObjString("AliAnalysisTaskMeanTracklets.cxx"));
  fileList.Add(new TObjString("AliAnalysisTaskMeanTracklets.h"));
  
  // --- grid specific setup ---
  TString dataDir = "/alice/sim/2016/LHC16h8b";
  TString dataPattern = "AOD/*AliAOD.root";
  TString runFormat = "%d";
  TString outDir = "Sim/LHC15n/LHC16h8b/AOD/MeanTracklets";
  TString analysisMacroName = "MeanTracklets";
  Int_t ttl = 30000;
  Int_t maxFilesPerJob = 10;
  Int_t maxMergeFiles = 10;
  Int_t maxMergeStages = 1;
  
  // --- saf3 specific setup ---
  Bool_t splitDataset = kFALSE;
  
  gROOT->LoadMacro("$HOME/Work/Alice/Macros/Facilities/runTaskFacilities.C");
  
  // --- prepare the analysis environment ---
  Int_t mode = PrepareAnalysis(smode, inputFileName, extraLibs, extraIncs, extraTasks, extraPkgs, pathList, fileList);
  if (!refInput.IsNull()) {
    CopyInputFileLocally(refInput.Data(), "ReferenceResults.root", 'a');
    fileList.Add(new TObjString("ReferenceResults.root"));
  }
  
  // --- run the analysis (saf3 is a special case as the analysis is launched on the server) ---
  if (mode == kSAF3Connect) {
    
    if (!RunAnalysisOnSAF3(fileList, aliphysicsVersion, inputFileName, splitDataset)) return;
    
  } else {
    
    CreateAnalysisTrain(applyPS, applyPileupCuts, isMC, refInput);
    
    if (smode == "saf3" && splitDataset) AliAnalysisManager::GetAnalysisManager()->SetSkipTerminate(kTRUE);
    
    RunAnalysis(smode, inputFileName, rootVersion, alirootVersion, aliphysicsVersion, extraLibs, extraIncs, extraTasks, extraPkgs, dataDir, dataPattern, outDir, analysisMacroName, runFormat, ttl, maxFilesPerJob, maxMergeFiles, maxMergeStages);
    
  }
  
}

//______________________________________________________________________________
void CreateAnalysisTrain(Bool_t applyPS, Bool_t applyPileupCuts, Bool_t isMC, TString refInput)
{
  /// create the analysis train and configure it
  
  // analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("MeanTrackletsAnalysis");
  //mgr->SetDebugLevel(3);
  
  // AOD handler
  AliInputEventHandler* aodH = new AliAODInputHandler;
  mgr->SetInputEventHandler(aodH);
  
  if (applyPS) {
    // event selection
    gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
    AliPhysicsSelectionTask* physicsSelection = AddTaskPhysicsSelection(isMC,applyPileupCuts);
    if(!physicsSelection) {
      Error("CreateAnalysisTrain","AliPhysicsSelectionTask not created!");
      return;
    }
  }
  /*
  // multiplicity/centrality selection
  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
  AliMultSelectionTask *mult = AddTaskMultSelection(kFALSE);
  if (applyPS) mult->SelectCollisionCandidates(AliVEvent::kINT7inMUON | AliVEvent::kINT7);
//  if (applyPS) mult->SelectCollisionCandidates(AliVEvent::kMuonUnlikeLowPt7);
  */
  // MeanTracklets task
  gROOT->LoadMacro("AddTaskMeanTracklets.C");
  AliAnalysisTaskMeanTracklets* meanTracklets = AddTaskMeanTracklets();
  if(!meanTracklets) {
    Error("CreateAnalysisTrain","AliAnalysisTaskMeanTracklets not created!");
    return;
  }
  if (applyPS) meanTracklets->ApplyPhysicsSelection(AliVEvent::kINT7inMUON | AliVEvent::kINT7);
//  if (applyPS) meanTracklets->ApplyPhysicsSelection(AliVEvent::kMuonUnlikeLowPt7);
//  meanTracklets->SelectTrigger("MB");
//  meanTracklets->SelectTrigger("CINT7-B-NOPF-");
//  meanTracklets->SelectTrigger("CMUL7-B-NOPF-MUFAST");
//  if (isMC) meanTracklets->RejectSD();
//  meanTracklets->RejectPUFromSPD();
//  meanTracklets->DisableSPDVtxQA();
//  meanTracklets->Reject0Tracklet();
  if (!refInput.IsNull() && !SetMeanNtrkVsZvtxRef(meanTracklets)) return;
//  meanTracklets->UseBinomial();
  
}

//______________________________________________________________________________
Bool_t SetMeanNtrkVsZvtxRef(TObject* meanTracklets)
{
  /// set the reference <Ntrk> vs Z profile
  
  TFile *file = new TFile("ReferenceResults.root", "read");
  if (!file || !file->IsOpen()) {
    printf("cannot open file ReferenceResults.root\n");
    return kFALSE;
  }
  
  TProfile *p = static_cast<TProfile*>(file->FindObjectAny("fpMeanNtrkVsZvtx"));
  if (!p) {
    Error("CreateAnalysisTrain","cannot find the reference <Ntrk> vs Z profile!");
    return kFALSE;
  }
  
  AliAnalysisTaskMeanTracklets *mT = static_cast<AliAnalysisTaskMeanTracklets*>(meanTracklets);
//  mT->SetMeanNtrkVsZvtxRef(*p);
  mT->SetMeanNtrkVsZvtxRef(*p, p->GetMaximum());
  
  file->Close();
  delete file;
  
  return kTRUE;
  
}
