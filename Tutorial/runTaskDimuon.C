/*
 *  runTaskDimuon.C
 *  aliroot
 *
 *  Created by Philippe Pillot on 10/04/17.
 *  Copyright 2017 SUBATECH. All rights reserved.
 *
 */

//______________________________________________________________________________
void runTaskDimuon(TString smode = "local", TString inputFileName = "AliAOD.Muons.root")
{
  /// Run the baseline task to test the framework
  
  // --- general analysis setup ---
  TString rootVersion = "";
  TString alirootVersion = "";
  TString aliphysicsVersion = "vAN-20170409-1";
  TString extraLibs="";
  TString extraIncs="include";
  TString extraTasks="AliAnalysisTaskDimuon";
  TString extraPkgs="";
  TList pathList; pathList.SetOwner();
  pathList.Add(new TObjString("$WORK/Macros/Tutorial"));
  TList fileList; fileList.SetOwner();
  fileList.Add(new TObjString("runTaskDimuon.C"));
  fileList.Add(new TObjString("AddTaskDimuon.C"));
  fileList.Add(new TObjString("AliAnalysisTaskDimuon.cxx"));
  fileList.Add(new TObjString("AliAnalysisTaskDimuon.h"));
  
  // --- grid specific setup ---
  TString dataDir = "/alice/cern.ch/user/p/ppillot/Sim/LHC15o/muTuneCMSL7/results";
  TString dataPattern = "*AliESDs.root";
  TString runFormat = "%d";
  TString outDir = "Sim/LHC15o/muTuneCMSL7/Eff";
  TString analysisMacroName = "Dimuon";
  Int_t ttl = 30000;
  Int_t maxFilesPerJob = 100;
  Int_t maxMergeFiles = 10;
  Int_t maxMergeStages = 2;
  
  // --- saf3 specific setup ---
  Bool_t splitDataset = kFALSE;
  
  gROOT->LoadMacro("$HOME/Work/Alice/Macros/Facilities/runTaskFacilities.C");
  
  // --- prepare the analysis environment ---
  Int_t mode = PrepareAnalysis(smode, inputFileName, extraLibs, extraIncs, extraTasks, extraPkgs, pathList, fileList);
  
  // --- run the analysis (saf3 is a special case as the analysis is launched on the server) ---
  if (mode == kSAF3Connect) {
    
    if (!RunAnalysisOnSAF3(fileList, aliphysicsVersion, inputFileName, splitDataset)) return;
    
  } else {
    
    CreateAnalysisTrain();
    
    if (smode == "saf3" && splitDataset) AliAnalysisManager::GetAnalysisManager()->SetSkipTerminate(kTRUE);
    
    RunAnalysis(smode, inputFileName, rootVersion, alirootVersion, aliphysicsVersion, extraLibs, extraIncs, extraTasks, extraPkgs, dataDir, dataPattern, outDir, analysisMacroName, runFormat, ttl, maxFilesPerJob, maxMergeFiles, maxMergeStages);
    
  }
  
}

//______________________________________________________________________________
void CreateAnalysisTrain()
{
  /// create the analysis train and configure it
  
  // analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("DimuonAnalysis");
  //mgr->SetDebugLevel(3);
  
  // AOD handler
  AliInputEventHandler* aodH = new AliAODInputHandler;
  mgr->SetInputEventHandler(aodH);
  /*
  // multiplicity/centrality selection
  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
  AliMultSelectionTask *mult = AddTaskMultSelection(kFALSE);
  mult->SelectCollisionCandidates(AliVEvent::kMuonUnlikeLowPt7 | AliVEvent::kMuonLikeLowPt7);
  */
  
  // dimuon task
  gROOT->LoadMacro("AddTaskDimuon.C");
  AliAnalysisTaskDimuon* dimu = AddTaskDimuon();
  if(!dimu) {
    Error("CreateAnalysisTrain","AliAnalysisTaskDimuon not created!");
    return;
  }
  dimu->SelectCollisionCandidates(AliVEvent::kMuonUnlikeLowPt7 | AliVEvent::kMuonLikeLowPt7);
  dimu->SetDefaultMuonTrackCuts(kFALSE);
  
}

