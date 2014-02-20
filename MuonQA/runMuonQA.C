/*
 *  runMuonQA.C
 *  aliroot
 *
 *  Created by Philippe Pillot on 21/06/10.
 *  Copyright 2010 Subatech. All rights reserved.
 *
 */

TString rootVersion = "v5-34-08";
TString alirootVersion = "v5-05-27-AN";
TString dataDir = "/alice/cern.ch/user/p/ppillot/Sim/LHC12h/muon_calo/mu/tunepA/results";
TString dataPattern = "*AliESDs.root";
//TString dataDir = "/alice/data/2010/LHC10h";
//TString dataPattern = "pass2/*AliESDs.root";
TString runFormat = "%d";
//TString runFormat = "%09d";
TString outDir = "Sim/LHC12h/muon_calo/mu/tunepA/MuonQA";
//TString outDir = "PbPb2.76TeV/LHC10h/pass2/Embedding/MuonQA";
Int_t ttl = 30000;
Int_t maxFilesPerJob = 20;
Int_t maxMergeFiles = 10;
Int_t maxMergeStages = 2;

//______________________________________________________________________________
void runMuonQA(TString smode="test", TString inputFileName="runList.txt", Bool_t simu = kTRUE,
	       Bool_t selectPhysics = kFALSE, Bool_t selectTrigger = kFALSE, TString sTriggerMask = "mb",
	       Bool_t selectMatched = kTRUE, Bool_t applyAccCut = kTRUE, Short_t selectCharge = 0)
{
  /// run the MUON QA task
  
  gROOT->LoadMacro("$ALICE/Macros/Facilities/runTaskFacilities.C");
  
  // --- Check runing mode ---
  Int_t mode = GetMode(smode, inputFileName);
  if(mode < 0) {
    Error("runMuonQA","Please provide either an ESD root file a collection of ESDs or a dataset.");
    return;
  }
  
  // --- copy files needed for this analysis ---
  TList pathList; pathList.SetOwner();
  pathList.Add(new TObjString("$WORK/Macros/MuonQA"));
  pathList.Add(new TObjString("$WORK/Macros/MuonPhysics"));
  pathList.Add(new TObjString("$ALICE_ROOT/PWGPP/MUON/lite"));
  pathList.Add(new TObjString("$ALICE_ROOT/PWGPP/PilotTrain"));
  TList fileList; fileList.SetOwner();
  fileList.Add(new TObjString("runMuonQA.C"));
  fileList.Add(new TObjString("AddTaskMuonQA.C"));
  fileList.Add(new TObjString("AliAnalysisTaskMuonQA.cxx"));
  fileList.Add(new TObjString("AliAnalysisTaskMuonQA.h"));
  fileList.Add(new TObjString("AddTaskMuonPhysics.C"));
  fileList.Add(new TObjString("AliAnalysisTaskMuonPhysics.cxx"));
  fileList.Add(new TObjString("AliAnalysisTaskMuonPhysics.h"));
  CopyFileLocally(pathList, fileList);
  
  // --- prepare environment ---
  TString extraLibs="CORRFW:PWGmuon";
  TString extraIncs="include:PWG/muon";
  TString extraTasks="AliAnalysisTaskMuonQA:AliAnalysisTaskMuonPhysics";
  LoadAlirootLocally(extraLibs, extraIncs, extraTasks);
  AliAnalysisGrid *alienHandler = 0x0;
  if (mode == kProof || mode == kProofLite) LoadAlirootOnProof(smode, rootVersion, alirootVersion, extraLibs, extraIncs, extraTasks, kTRUE);
  else if (mode == kGrid || mode == kTerminate) {
    TString analysisMacroName = "MuonQA";
    alienHandler = static_cast<AliAnalysisGrid*>(CreateAlienHandler(smode, rootVersion, alirootVersion, inputFileName, dataDir, dataPattern, outDir, extraLibs, extraIncs, extraTasks, analysisMacroName, runFormat, ttl, maxFilesPerJob, maxMergeFiles, maxMergeStages));
    if (!alienHandler) return;
  }
  
  // --- Create the analysis train ---
  CreateAnalysisTrain(simu, selectPhysics, selectTrigger, sTriggerMask, selectMatched, applyAccCut, selectCharge, alienHandler);
  
  // --- Create input object ---
  TObject* inputObj = CreateInputObject(mode, inputFileName);
  
  // --- start analysis ---
  StartAnalysis(mode, inputObj);
  
}

//______________________________________________________________________________
void CreateAnalysisTrain(Bool_t simu, Bool_t selectPhysics, Bool_t selectTrigger, TString sTriggerMask,
			 Bool_t selectMatched, Bool_t applyAccCut, Short_t selectCharge, TObject* alienHandler)
{
  /// create the analysis train and configure it
  
  // analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("MuonQAAnalysis");
  
  // Connect plugin to the analysis manager if any
  if (alienHandler) mgr->SetGridHandler(static_cast<AliAnalysisGrid*>(alienHandler));
  
  // ESD input
  AliESDInputHandler* esdH = new AliESDInputHandler();
  esdH->SetReadFriends(kFALSE);
  esdH->SetInactiveBranches("*");
  esdH->SetActiveBranches("MuonTracks MuonClusters MuonPads AliESDRun. AliESDHeader. AliMultiplicity. AliESDFMD. AliESDVZERO. SPDVertex. PrimaryVertex. AliESDZDC. AliESDTZERO.");
  mgr->SetInputEventHandler(esdH);
  
  // event selection
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
  AliPhysicsSelectionTask* physicsSelection = AddTaskPhysicsSelection(simu);
  if(!physicsSelection) {
    Error("RunMuonQA","AliPhysicsSelectionTask not created!");
    return;
  }
  
  // track selection
  AliMuonTrackCuts trackCuts("stdCuts", "stdCuts");
  trackCuts.SetAllowDefaultParams();
  trackCuts.SetFilterMask(AliMuonTrackCuts::kMuMatchApt | AliMuonTrackCuts::kMuEta |
			  AliMuonTrackCuts::kMuThetaAbs);
  if (simu) trackCuts.SetIsMC();
  
  // determine the trigger mask
  sTriggerMask.ToLower();
  UInt_t triggerMask = AliVEvent::kAny;
  if (sTriggerMask == "muon") triggerMask = AliVEvent::kMUON;
  else if (sTriggerMask == "mb") triggerMask = AliVEvent::kMB;
  else if (sTriggerMask == "highmult") triggerMask = AliVEvent::kHighMult;
  
  // Muon QA analysis
  gROOT->LoadMacro("AddTaskMuonQA.C");
  AliAnalysisTaskMuonQA* muonQA = AddTaskMuonQA(selectPhysics, selectMatched, applyAccCut, selectTrigger, triggerMask, selectCharge);
  if(!muonQA) {
    Error("RunMuonQA","AliAnalysisTaskMuonQA not created!");
    return;
  }
  
  // Physics results
  gROOT->LoadMacro("AddTaskMuonPhysics.C");
  AliAnalysisTaskMuonPhysics* physics = AddTaskMuonPhysics("0");
  if(!physics) {
    Error("RunMuonQA","AliAnalysisTaskMuonPhysics not created!");
    return;
  }
  if (selectPhysics) physics->SelectCollisionCandidates(AliVEvent::kAny);
  physics->SetMuonTrackCuts(trackCuts);
}

