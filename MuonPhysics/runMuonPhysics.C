/*
 *  runMuonPhysics.C
 *  aliroot
 *
 *  Created by Philippe Pillot on 05/11/10.
 *  Copyright 2010 SUBATECH. All rights reserved.
 *
 */

//TString rootVersion = "v5-34-05";
//TString alirootVersion = "v5-04-65-AN";
TString rootVersion = "v5-34-08";
TString alirootVersion = "v5-05-29-AN";
TString dataDir = "/alice/data/2013/LHC13f";
TString dataPattern = "ESDs/muon_pass2/*AliESDs.root";
TString runFormat = "%09d";
TString outDir = "Data/LHC13f/muon_pass2/Phys_MB";
Int_t ttl = 30000;
Int_t maxFilesPerJob = 20;
Int_t maxMergeFiles = 10;
Int_t maxMergeStages = 2;

//______________________________________________________________________________
void runMuonPhysics(TString smode = "local", TString inputFileName = "AliESDs.root", Bool_t mc = kFALSE)
{
  /// Extract some physical quantities
  
  gROOT->LoadMacro("$WORK/Macros/Facilities/runTaskFacilities.C");
  
  // --- Check runing mode ---
  Int_t mode = GetMode(smode, inputFileName);
  if(mode < 0){
    Error("runMuonPhysics","Please provide either an ESD root file a collection of ESDs or a dataset.");
    return;
  }
  
  // --- copy files needed for this analysis ---
  TList pathList; pathList.SetOwner();
  pathList.Add(new TObjString("$WORK/Macros/MuonPhysics"));
  TList fileList; fileList.SetOwner();
  fileList.Add(new TObjString("runMuonPhysics.C"));
  fileList.Add(new TObjString("AddTaskMuonPhysics.C"));
  fileList.Add(new TObjString("AliAnalysisTaskMuonPhysics.cxx"));
  fileList.Add(new TObjString("AliAnalysisTaskMuonPhysics.h"));
  CopyFileLocally(pathList, fileList);
  
  // --- prepare environment ---
  TString extraLibs="CORRFW:PWGmuon";
  TString extraIncs="include/PWG/muon";
  TString extraTasks="AliAnalysisTaskMuonPhysics";
//  TString extraLibs="RAWDatabase:RAWDatarec:CDB:STEER:MUONcore:MUONmapping:MUONcalib:MUONgeometry:MUONtrigger:MUONraw:MUONbase:MUONshuttle:MUONrec:RAWDatasim:MUONsim:MUONevaluation:CORRFW:PWGmuon";
//  TString extraIncs="include:MUON:MUON/mapping:PWG/muon";
//  TString extraTasks="AliAnalysisTaskMuonPhysics";
  LoadAlirootLocally(extraLibs, extraIncs, extraTasks);
  AliAnalysisGrid *alienHandler = 0x0;
  if (mode == kProof || mode == kProofLite) LoadAlirootOnProof(smode, rootVersion, alirootVersion, extraLibs, extraIncs, extraTasks, kTRUE);
  else if (mode == kGrid || mode == kTerminate) {
    TString analysisMacroName = "Physics";
    alienHandler = static_cast<AliAnalysisGrid*>(CreateAlienHandler(smode, rootVersion, alirootVersion, inputFileName, dataDir, dataPattern, outDir, extraLibs, extraIncs, extraTasks, analysisMacroName, runFormat, ttl, maxFilesPerJob, maxMergeFiles, maxMergeStages));
    if (!alienHandler) return;
  }
  
  // --- Create the analysis train ---
  CreateAnalysisTrain(mc, alienHandler);
  
  // --- Create input object ---
  TObject* inputObj = CreateInputObject(mode, inputFileName);
  
  // --- start analysis ---
  StartAnalysis(mode, inputObj);
  
}

//______________________________________________________________________________
void CreateAnalysisTrain(Bool_t mc, TObject* alienHandler)
{
  /// create the analysis train and configure it
  
  // analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("MuonPhysicsAnalysis");
  
  // Connect plugin to the analysis manager if any
  if (alienHandler) mgr->SetGridHandler(static_cast<AliAnalysisGrid*>(alienHandler));
  //mgr->SetDebugLevel(3);
  
  // ESD input
  AliESDInputHandler* esdH = new AliESDInputHandler();
  esdH->SetReadFriends(kFALSE);
  esdH->SetInactiveBranches("*");
  esdH->SetActiveBranches("MuonTracks MuonClusters MuonPads AliESDRun. AliESDHeader. AliMultiplicity. AliESDFMD. AliESDVZERO. SPDVertex. PrimaryVertex. AliESDZDC. AliESDTZERO.");
  mgr->SetInputEventHandler(esdH);
  
  // event selection
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
  AliPhysicsSelectionTask* physicsSelection = AddTaskPhysicsSelection(mc);
  if(!physicsSelection) {
    Error("CreateAnalysisTrain","AliPhysicsSelectionTask not created!");
    return;
  }
  
  // centrality selection
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskCentrality.C");
  AliCentralitySelectionTask *taskCentrality = AddTaskCentrality();
  if(!taskCentrality) {
    Error("CreateAnalysisTrain","AliCentralitySelectionTask not created!");
    return;
  }
  if (mc) taskCentrality->SetMCInput();
  
  // track selection
  AliMuonTrackCuts trackCuts("stdCuts", "stdCuts");
  trackCuts.SetAllowDefaultParams();
//  trackCuts.SetFilterMask(0);
  trackCuts.SetFilterMask(AliMuonTrackCuts::kMuMatchLpt | AliMuonTrackCuts::kMuEta |
			  AliMuonTrackCuts::kMuThetaAbs | AliMuonTrackCuts::kMuPdca);
//  trackCuts.SetPassNumber(2);
  
  // physics task
  gROOT->LoadMacro("AddTaskMuonPhysics.C");
  AliAnalysisTaskMuonPhysics* physics1 = AddTaskMuonPhysics("");
  if(!physics1) {
    Error("CreateAnalysisTrain","AliAnalysisTaskMuonPhysics not created!");
    return;
  }
  physics1->SelectCollisionCandidates(AliVEvent::kAnyINT);
//  physics1->SelectCollisionCandidates(AliVEvent::kINT8 |
//				      AliVEvent::kMuonSingleLowPt8 | AliVEvent::kMuonSingleHighPt8 |
//				      AliVEvent::kMuonLikeLowPt8 | AliVEvent::kMuonUnlikeLowPt8);
//  physics1->SelectCollisionCandidates(AliVEvent::kMuonUnlikeLowPt0);
//  physics1->SelectCollisionCandidates(AliVEvent::kMUL7 | AliVEvent::kMUU7 |
//				      AliVEvent::kMuonLikeLowPt8 | AliVEvent::kMuonUnlikeLowPt8);
//  physics1->SelectCollisionCandidates(AliVEvent::kMuonUnlikeLowPt8 | AliVEvent::kMUU7);
//  physics1->SelectCentrality(0., 90.);
  physics1->SetMuonTrackCuts(trackCuts);
  
}

