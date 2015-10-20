/*
 *  runMuonEfficiency.C
 *  aliroot
 *
 *  Created by Philippe Pillot on 16/12/10.
 *  Copyright 2010 SUBATECH. All rights reserved.
 *
 */

TString rootVersion = "";
TString alirootVersion = "";
TString aliphysicsVersion = "vAN-20150930-1";
TString dataDir = "/alice/data/2015/LHC15g";
TString dataPattern = "muon_calo_pass1/*AliESDs.root";
TString runFormat = "%09d";
TString outDir = "Data/LHC15g/muon_calo_pass1/Eff";
Int_t ttl = 30000;
Int_t maxFilesPerJob = 20;
Int_t maxMergeFiles = 10;
Int_t maxMergeStages = 2;

TString alignStorage = "";
//TString alignStorage = "alien://folder=/alice/simulation/2008/v4-15-Release/Residual";
//TString alignStorage = "alien://folder=/alice/cern.ch/user/p/ppillot/OCDB_PbPbSim";

//______________________________________________________________________________
void runMuonEfficiency(TString smode = "local", TString inputFileName = "AliESDs.root",
		       Bool_t applyPhysSel = kTRUE, Bool_t mc = kFALSE, Bool_t embedding = kFALSE)
{
  /// Study the MUON performances
  
  gROOT->LoadMacro("$WORK/Macros/Facilities/runTaskFacilities.C");
  
  // --- Check runing mode ---
  Int_t mode = GetMode(smode, inputFileName);
  if(mode < 0) {
    Error("runMuonEfficiency","Please provide either an ESD root file a collection of ESDs or a dataset.");
    return;
  }
  
  // --- copy files needed for this analysis ---
  TList pathList; pathList.SetOwner();
  pathList.Add(new TObjString("$WORK/Macros/MuonEfficiency"));
  pathList.Add(new TObjString("$WORK/Macros/MuonPhysics"));
  pathList.Add(new TObjString("$DEVPHYS/PWGPP/MUON/dep"));
  pathList.Add(new TObjString("$DEVPHYS/PWGPP/PilotTrain"));
  pathList.Add(new TObjString("$DEVPHYS/../inst/PARfiles"));
  TList fileList; fileList.SetOwner();
  fileList.Add(new TObjString("runMuonEfficiency.C"));
  fileList.Add(new TObjString("AddTaskMUONTrackingEfficiency.C"));
  fileList.Add(new TObjString("PWGPPMUONdep.par"));
  fileList.Add(new TObjString("AddTaskMuonQA.C"));
  fileList.Add(new TObjString("PWGPPMUONlite.par"));
  fileList.Add(new TObjString("AddTaskMuonPhysics.C"));
  fileList.Add(new TObjString("AliAnalysisTaskMuonPhysics.cxx"));
  fileList.Add(new TObjString("AliAnalysisTaskMuonPhysics.h"));
  CopyFileLocally(pathList, fileList);
  
  // --- prepare environment ---
  TString extraLibs="PWGmuon";
  TString extraIncs="include";
  TString extraTasks="AliAnalysisTaskMuonPhysics";
  TString extraPkgs="PWGPPMUONdep:PWGPPMUONlite";
  LoadAlirootLocally(extraLibs, extraIncs, extraTasks, extraPkgs);
  AliAnalysisGrid *alienHandler = 0x0;
  if (mode == kProof || mode == kProofLite) LoadAlirootOnProof(smode, rootVersion, aliphysicsVersion, extraLibs, extraIncs, extraTasks, extraPkgs, kTRUE);
  else if (mode == kGrid || mode == kTerminate) {
    TString analysisMacroName = "Eff";
    alienHandler = static_cast<AliAnalysisGrid*>(CreateAlienHandler(smode, alirootVersion, aliphysicsVersion, inputFileName, dataDir, dataPattern, outDir, extraLibs, extraIncs, extraTasks, extraPkgs, analysisMacroName, runFormat, ttl, maxFilesPerJob, maxMergeFiles, maxMergeStages));
    if (!alienHandler) return;
  }
  
  // --- Create the analysis train ---
  CreateAnalysisTrain(applyPhysSel, mc, embedding, alienHandler);
  
  // --- Create input object ---
  TObject* inputObj = CreateInputObject(mode, inputFileName);
  
  // --- start analysis ---
  StartAnalysis(mode, inputObj);
  
}

//______________________________________________________________________________
void CreateAnalysisTrain(Bool_t applyPhysSel, Bool_t mc, Bool_t embedding, TObject* alienHandler)
{
  /// create the analysis train and configure it
  
  // analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("MuonEfficiencyAnalysis");
  
  // Connect plugin to the analysis manager if any
  if (alienHandler) mgr->SetGridHandler(static_cast<AliAnalysisGrid*>(alienHandler));
  
  // ESD input
  AliESDInputHandler* esdH = new AliESDInputHandler();
  esdH->SetReadFriends(kFALSE);
  esdH->SetInactiveBranches("*");
  esdH->SetActiveBranches("MuonTracks MuonClusters MuonPads AliESDRun. AliESDHeader. AliMultiplicity. AliESDFMD. AliESDVZERO. SPDVertex. PrimaryVertex. AliESDZDC. AliESDTZERO.");
  mgr->SetInputEventHandler(esdH);
  
  // event selection
  UInt_t offlineTriggerMask;
  if (applyPhysSel) {
    gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
    AliPhysicsSelectionTask* physicsSelection = AddTaskPhysicsSelection(mc && !embedding);
    if(!physicsSelection) {
      Error("CreateAnalysisTrain","AliPhysicsSelectionTask not created!");
      return;
    }
    //offlineTriggerMask = AliVEvent::kAny;
    offlineTriggerMask = AliVEvent::kMUU7;
    //offlineTriggerMask = AliVEvent::kMUU7 | AliVEvent::kMuonUnlikeLowPt8;
  }
  /*
  // centrality selection
  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskCentrality.C");
  AliCentralitySelectionTask *taskCentrality = AddTaskCentrality();
  if(!taskCentrality) {
    Error("CreateAnalysisTrain","AliCentralitySelectionTask not created!");
    return;
  }
  if (applyPhysSel) taskCentrality->SelectCollisionCandidates(offlineTriggerMask);
  if (mc && !embedding) taskCentrality->SetMCInput();
  */
  // track selection
  AliMuonTrackCuts trackCuts("stdCuts", "stdCuts");
  trackCuts.SetAllowDefaultParams();
//  trackCuts.SetFilterMask(0);
//  trackCuts.SetCustomParamFromRun(169099, "pass2_muon");
//  trackCuts.CustomParam()->SetChi2NormCut(3.5);
  trackCuts.SetFilterMask(AliMuonTrackCuts::kMuMatchLpt | AliMuonTrackCuts::kMuEta |
			  AliMuonTrackCuts::kMuThetaAbs | AliMuonTrackCuts::kMuPdca);
  trackCuts.SetIsMC(mc && !embedding);
  
  // Muon efficiency analysis
  gROOT->LoadMacro("AddTaskMUONTrackingEfficiency.C");
  AliAnalysisTaskMuonTrackingEff* muonEfficiency = AddTaskMUONTrackingEfficiency(trackCuts,"");
  if(!muonEfficiency) {
    Error("CreateAnalysisTrain","AliAnalysisTaskMuonTrackingEff not created!");
    return;
  }
  if (applyPhysSel) muonEfficiency->SelectCollisionCandidates(offlineTriggerMask);
  if (!alignStorage.IsNull()) muonEfficiency->SetAlignStorage(alignStorage.Data());
  //muonEfficiency->SetRecoParamStorage("alien://folder=/alice/cern.ch/user/p/ppillot/OCDB2012_newReco");
  muonEfficiency->SetMuonPtCut(1.);
  muonEfficiency->UseMCLabel(kFALSE);
  muonEfficiency->EnableDisplay(kFALSE);
  
  // QA task
  gROOT->LoadMacro("AddTaskMuonQA.C");
  AliAnalysisTaskMuonQA* muonQA = AddTaskMuonQA(kFALSE, kFALSE, kFALSE, 0);
  if(!muonQA) {
    Error("CreateAnalysisTrain","AliAnalysisTaskMuonQA not created!");
    return;
  }
  if (applyPhysSel) muonQA->SelectCollisionCandidates(offlineTriggerMask);
  muonQA->SetTrackCuts(&trackCuts);

  // Physics results
  gROOT->LoadMacro("AddTaskMuonPhysics.C");
  AliAnalysisTaskMuonPhysics* physics = AddTaskMuonPhysics("0");
  if(!physics) {
    Error("CreateAnalysisTrain","AliAnalysisTaskMuonPhysics not created!");
    return;
  }
  if (applyPhysSel) physics->SelectCollisionCandidates(offlineTriggerMask);
  physics->SetMuonTrackCuts(trackCuts);
  physics->VersusRun(kTRUE);
  
}

