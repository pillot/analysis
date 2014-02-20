/*
 *  runMuonEfficiency.C
 *  aliroot
 *
 *  Created by Philippe Pillot on 16/12/10.
 *  Copyright 2010 SUBATECH. All rights reserved.
 *
 */

TString rootVersion = "v5-34-08-4";
TString alirootVersion = "v5-05-63-AN";
TString dataDir = "/alice/cern.ch/user/l/lvalenci/MySimulation/v5-02-Rev-37/LHC11hPass2/OldRecoParamsALIGNMENT/Output";
TString dataPattern = "*AliESDs.root";
TString runFormat = "%d";
TString outDir = "Sim/LHC11h/muLizardo/Eff";
Int_t ttl = 30000;
Int_t maxFilesPerJob = 10;
Int_t maxMergeFiles = 10;
Int_t maxMergeStages = 2;

//TString alignStorage = "alien://folder=/alice/simulation/2008/v4-15-Release/Residual";
TString alignStorage = "alien://folder=/alice/cern.ch/user/p/ppillot/OCDB_PbPbSim";

//______________________________________________________________________________
void runMuonEfficiency(TString smode = "local", TString inputFileName = "AliESDs.root",
		       Bool_t applyPhysSel = kFALSE, Bool_t mc = kTRUE, Bool_t embedding = kFALSE)
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
  pathList.Add(new TObjString("$DEV/aliroot/PWGPP/MUON/dep"));
  TList fileList; fileList.SetOwner();
  fileList.Add(new TObjString("runMuonEfficiency.C"));
  fileList.Add(new TObjString("AddTaskMUONTrackingEfficiency.C"));
  fileList.Add(new TObjString("AliAnalysisTaskMuonTrackingEff.cxx"));
  fileList.Add(new TObjString("AliAnalysisTaskMuonTrackingEff.h"));
  fileList.Add(new TObjString("AddTaskMUONTrackingEfficiency_old.C"));
  fileList.Add(new TObjString("AliAnalysisTaskMuonTrackingEff_old.cxx"));
  fileList.Add(new TObjString("AliAnalysisTaskMuonTrackingEff_old.h"));
  fileList.Add(new TObjString("AddTaskMuonPhysics.C"));
  fileList.Add(new TObjString("AliAnalysisTaskMuonPhysics.cxx"));
  fileList.Add(new TObjString("AliAnalysisTaskMuonPhysics.h"));
  CopyFileLocally(pathList, fileList);
  
  // --- prepare environment ---
  TString extraLibs="RAWDatabase:CDB:STEER:MUONcore:MUONmapping:MUONcalib:MUONgeometry:MUONtrigger:MUONraw:MUONbase:MUONrec:CORRFW:PWGmuon";
  TString extraIncs="include:MUON:MUON/mapping";
  TString extraTasks="AliAnalysisTaskMuonTrackingEff:AliAnalysisTaskMuonTrackingEff_old:AliAnalysisTaskMuonPhysics";
  LoadAlirootLocally(extraLibs, extraIncs, extraTasks);
  AliAnalysisGrid *alienHandler = 0x0;
  if (mode == kProof || mode == kProofLite) LoadAlirootOnProof(smode, rootVersion, alirootVersion, extraLibs, extraIncs, extraTasks, kTRUE);
  else if (mode == kGrid || mode == kTerminate) {
    TString analysisMacroName = "Eff";
    alienHandler = static_cast<AliAnalysisGrid*>(CreateAlienHandler(smode, rootVersion, alirootVersion, inputFileName, dataDir, dataPattern, outDir, extraLibs, extraIncs, extraTasks, analysisMacroName, runFormat, ttl, maxFilesPerJob, maxMergeFiles, maxMergeStages));
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
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
    AliPhysicsSelectionTask* physicsSelection = AddTaskPhysicsSelection(mc && !embedding);
    if(!physicsSelection) {
      Error("CreateAnalysisTrain","AliPhysicsSelectionTask not created!");
      return;
    }
    //offlineTriggerMask = AliVEvent::kAny;
    //offlineTriggerMask = AliVEvent::kMUU7;
    offlineTriggerMask = AliVEvent::kMUU7 | AliVEvent::kMuonUnlikeLowPt8;
  }
  /*
  // centrality selection
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskCentrality.C");
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
			  AliMuonTrackCuts::kMuThetaAbs);
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
  muonEfficiency->UseMCLabel(kTRUE);
  muonEfficiency->EnableDisplay(kTRUE);
  
  // Muon efficiency analysis (old without cut but with MC label)
  gROOT->LoadMacro("AddTaskMUONTrackingEfficiency_old.C");
  AliAnalysisTaskMuonTrackingEff_old* muonEfficiency_old_wocut_wMClabel = AddTaskMUONTrackingEfficiency_old(kFALSE, kFALSE,"old_wocut_wMClabel");
  if(!muonEfficiency_old_wocut_wMClabel) {
    Error("CreateAnalysisTrain","AliAnalysisTaskMuonTrackingEff_old not created!");
    return;
  }
  if (applyPhysSel) muonEfficiency_old_wocut_wMClabel->SelectCollisionCandidates(offlineTriggerMask);
  muonEfficiency_old_wocut_wMClabel->UseMCLabel(kTRUE);
  
  // Muon efficiency analysis (old with cut)
  AliAnalysisTaskMuonTrackingEff_old* muonEfficiency_old_wcut = AddTaskMUONTrackingEfficiency_old(kTRUE, kTRUE,"old_wcut");
  if(!muonEfficiency_old_wcut) {
    Error("CreateAnalysisTrain","AliAnalysisTaskMuonTrackingEff_old not created!");
    return;
  }
  if (applyPhysSel) muonEfficiency_old_wcut->SelectCollisionCandidates(offlineTriggerMask);
  muonEfficiency_old_wcut->PtCut(1.);
  
  // Muon efficiency analysis (old with cut and MC label)
  AliAnalysisTaskMuonTrackingEff_old* muonEfficiency_old_wcut_wMClabel = AddTaskMUONTrackingEfficiency_old(kTRUE, kTRUE,"old_wcut_wMClabel");
  if(!muonEfficiency_old_wcut_wMClabel) {
    Error("CreateAnalysisTrain","AliAnalysisTaskMuonTrackingEff_old not created!");
    return;
  }
  if (applyPhysSel) muonEfficiency_old_wcut_wMClabel->SelectCollisionCandidates(offlineTriggerMask);
  muonEfficiency_old_wcut_wMClabel->PtCut(1.);
  muonEfficiency_old_wcut_wMClabel->UseMCLabel(kTRUE);
  /*
  // Muon efficiency analysis -- with pDCA cut
  trackCuts.SetFilterMask(AliMuonTrackCuts::kMuMatchLpt | AliMuonTrackCuts::kMuEta |
			  AliMuonTrackCuts::kMuThetaAbs | AliMuonTrackCuts::kMuPdca);
  AliAnalysisTaskMuonTrackingEff* muonEfficiency2 = AddTaskMUONTrackingEfficiency(trackCuts,"pDCA");
  if(!muonEfficiency2) {
    Error("CreateAnalysisTrain","AliAnalysisTaskMuonTrackingEff not created!");
    return;
  }
  if (applyPhysSel) muonEfficiency2->SelectCollisionCandidates(offlineTriggerMask);
  
  // Muon efficiency analysis -- with chi2 cut
  trackCuts.SetFilterMask(AliMuonTrackCuts::kMuMatchLpt | AliMuonTrackCuts::kMuEta |
			  AliMuonTrackCuts::kMuThetaAbs | AliMuonTrackCuts::kMuTrackChiSquare);
  AliAnalysisTaskMuonTrackingEff* muonEfficiency3 = AddTaskMUONTrackingEfficiency(trackCuts,"chi2");
  if(!muonEfficiency3) {
    Error("CreateAnalysisTrain","AliAnalysisTaskMuonTrackingEff not created!");
    return;
  }
  if (applyPhysSel) muonEfficiency3->SelectCollisionCandidates(offlineTriggerMask);
  
  // Muon efficiency analysis -- with pDCA and chi2 cut
  trackCuts.SetFilterMask(AliMuonTrackCuts::kMuMatchLpt | AliMuonTrackCuts::kMuEta |
			  AliMuonTrackCuts::kMuThetaAbs | AliMuonTrackCuts::kMuPdca |
			  AliMuonTrackCuts::kMuTrackChiSquare);
  AliAnalysisTaskMuonTrackingEff* muonEfficiency4 = AddTaskMUONTrackingEfficiency(trackCuts,"pDCAchi2");
  if(!muonEfficiency4) {
    Error("CreateAnalysisTrain","AliAnalysisTaskMuonTrackingEff not created!");
    return;
  }
  if (applyPhysSel) muonEfficiency4->SelectCollisionCandidates(offlineTriggerMask);
  */
  // Physics results
  gROOT->LoadMacro("AddTaskMuonPhysics.C");
  AliAnalysisTaskMuonPhysics* physics = AddTaskMuonPhysics();
  if(!physics) {
    Error("CreateAnalysisTrain","AliAnalysisTaskMuonPhysics not created!");
    return;
  }
  if (applyPhysSel) physics->SelectCollisionCandidates(offlineTriggerMask);
  physics->SetMuonTrackCuts(trackCuts);
  
}

