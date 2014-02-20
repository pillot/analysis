/*
 *  runMuonPerformance.C
 *  aliroot_dev
 *
 *  Created by philippe pillot on 25/11/12.
 *  Copyright 2012 SUBATECH. All rights reserved.
 *
 */

TString rootVersion = "v5-34-02-1";
TString alirootVersion = "v5-04-05-AN";
TString dataDir = "/alice/cern.ch/user/p/ppillot/Sim/LHC12h/muon_calo/JPsipp7/MisIdealVtx0/results";
TString dataPattern = "*AliESDs.root";
TString runFormat = "%d";
TString outDir = "Sim/LHC12h/muon_calo/JPsipp7/MisIdealVtx0";
Int_t ttl = 30000;
Int_t maxFilesPerJob = 50;
Int_t maxMergeFiles = 10;
Int_t maxMergeStages = 3;

//______________________________________________________________________________
void runMuonPerformance(TString smode = "local", TString inputFileName = "AliESDs.root",
			Bool_t correctClusterResForSystematics = kTRUE, Bool_t fitClusterResiduals = kTRUE,
			Bool_t applyPhysicsSelection = kFALSE, Bool_t embedding = kFALSE)
{
  /// Study the MUON performances
  
  gROOT->LoadMacro("$WORK/Macros/Facilities/runTaskFacilities.C");
  
  // --- Check runing mode ---
  Int_t mode = GetMode(smode, inputFileName);
  if(mode < 0) {
    Error("runMuonPerformance","Please provide either an ESD root file a collection of ESDs or a dataset.");
    return;
  }
  
  // --- prepare environment ---
  TString extraLibs="RAWDatabase:RAWDatarec:CDB:STEER:MUONcore:MUONmapping:MUONcalib:MUONgeometry:MUONtrigger:MUONraw:MUONbase:MUONshuttle:MUONrec:RAWDatasim:MUONsim:MUONevaluation:CORRFW:PWGPPMUONdep:PWGmuon:PWGmuondep";
  TString extraIncs="";
  TString extraTasks="";
  LoadAlirootLocally(extraLibs, extraIncs, extraTasks);
  AliAnalysisGrid *alienHandler = 0x0;
  if (mode == kProof) LoadAlirootOnProof(smode, rootVersion, alirootVersion, extraLibs, extraIncs, extraTasks, kTRUE);
  else if (mode == kGrid || mode == kTerminate) {
    TString analysisMacroName = "Perf";
    alienHandler = static_cast<AliAnalysisGrid*>(CreateAlienHandler(smode, rootVersion, alirootVersion, inputFileName, dataDir, dataPattern, outDir, extraLibs, extraIncs, extraTasks, analysisMacroName, runFormat, ttl, maxFilesPerJob, maxMergeFiles, maxMergeStages));
    if (!alienHandler) return;
  }
  
  // --- Create the analysis train ---
  CreateAnalysisTrain(correctClusterResForSystematics, fitClusterResiduals, applyPhysicsSelection, embedding, alienHandler);
  
  // --- Create input object ---
  TObject* inputObj = CreateInputObject(mode, inputFileName);
  
  // --- start analysis ---
  StartAnalysis(mode, inputObj);
  
  if (smode == "offline") {
    TString homeDir = "/alice/cern.ch/user/p/ppillot";
    TString targetDir = Form("%s/%s", homeDir.Data(), outDir.Data());
    gSystem->Exec(Form("alien_cp file:%s.root alien://%s/%s.root", analysisMacroName.Data(), targetDir.Data(), analysisMacroName.Data()));
    gSystem->Exec(Form("alien_cp file:%s_merge.C alien://%s/%s_merge.C", analysisMacroName.Data(), targetDir.Data(), analysisMacroName.Data()));
    gSystem->Exec(Form("alien_cp file:%s_validation_merge.sh alien://%s/%s_validation_merge.sh", analysisMacroName.Data(), targetDir.Data(), analysisMacroName.Data()));
    gSystem->Exec(Form("alien_cp file:%s_merge.jdl alien://%s/results/%s_merge.jdl", analysisMacroName.Data(), targetDir.Data(), analysisMacroName.Data()));
    gSystem->Exec(Form("alien_cp file:%s_merge_final.jdl alien://%s/results/%s_merge_final.jdl", analysisMacroName.Data(), targetDir.Data(), analysisMacroName.Data()));
    gSystem->Exec(Form("alien_cp file:%s_merge.sh alien://%s/bin/%s_merge.sh", analysisMacroName.Data(), homeDir.Data(), analysisMacroName.Data()));
  }
  
}

//______________________________________________________________________________
void CreateAnalysisTrain(Bool_t correctClusterResForSystematics, Bool_t fitClusterResiduals,
			 Bool_t applyPhysicsSelection, Bool_t embedding, TObject* alienHandler)
{
  /// create the analysis train and configure it
  
  // analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("MuonPerformanceAnalysis");
  
  // Connect plugin to the analysis manager if any
  if (alienHandler) mgr->SetGridHandler(static_cast<AliAnalysisGrid*>(alienHandler));
  
  // ESD input
  AliESDInputHandler* esdH = new AliESDInputHandler();
  esdH->SetReadFriends(kFALSE);
  esdH->SetInactiveBranches("*");
  esdH->SetActiveBranches("MuonTracks MuonClusters MuonPads AliESDRun. AliESDHeader. AliMultiplicity. AliESDFMD. AliESDVZERO. SPDVertex. PrimaryVertex. AliESDZDC.");
  mgr->SetInputEventHandler(esdH);
  
  // Monte Carlo handler
  AliMCEventHandler* mcH = new AliMCEventHandler();
  mgr->SetMCtruthEventHandler(mcH);
  
  // event selection
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
  AliPhysicsSelectionTask* physicsSelection = embedding ? AddTaskPhysicsSelection() : AddTaskPhysicsSelection(kTRUE);
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
  if (!embedding) taskCentrality->SetMCInput();
  if (applyPhysicsSelection) taskCentrality->SelectCollisionCandidates(AliVEvent::kAny);
  
  // Add MC labels
  gROOT->LoadMacro("$ALICE_ROOT/PWG/muondep/AddTaskESDMCLabelAddition.C");
  AliAnalysisTaskESDMCLabelAddition *esdMCLabeltask = AddTaskESDMCLabelAddition();
  if(!esdMCLabeltask) {
    Error("CreateAnalysisTrain","AliAnalysisTaskESDMCLabelAddition not created!");
    return;
  }
  if (applyPhysicsSelection) esdMCLabeltask->SelectCollisionCandidates(AliVEvent::kAny);
  
  // Muon performance analysis
  gROOT->LoadMacro("$ALICE_ROOT/PWGPP/MUON/dep/AddTaskMuonPerformance.C");
  AliAnalysisTaskMuonPerformance* muonPerformance = AddTaskMuonPerformance(correctClusterResForSystematics, fitClusterResiduals);
  if(!muonPerformance) {
    Error("CreateAnalysisTrain","AliAnalysisTaskMuonPerformance not created!");
    return;
  }
  if (applyPhysicsSelection) muonPerformance->SelectCollisionCandidates(AliVEvent::kAny);
  muonPerformance->UseMCKinematics(kTRUE);
  muonPerformance->SetMCTrigLevelFromMatchTrk(kTRUE);
  
  // Muon efficiency analysis
  gROOT->LoadMacro("$ALICE_ROOT/PWGPP/MUON/dep/AddTaskMUONTrackingEfficiency.C");
  AliAnalysisTaskMuonTrackingEff* muonEfficiency = AddTaskMUONTrackingEfficiency(kTRUE, kTRUE);
  if(!muonEfficiency) {
    Error("CreateAnalysisTrain","AliAnalysisTaskMuonTrackingEff not created!");
    return;
  }
  if (applyPhysicsSelection) muonEfficiency->SelectCollisionCandidates(AliVEvent::kAny);
  muonEfficiency->UseMCLabel(kTRUE);
  
}

