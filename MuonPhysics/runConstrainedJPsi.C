/*
 *  runConstrainedJPsi.C
 *  aliroot
 *
 *  Created by Philippe Pillot on 07/06/12.
 *  Copyright 2012 SUBATECH. All rights reserved.
 *
 */


TString rootVersion = "v5-30-03-1";
TString alirootVersion = "v5-01-Rev-10";
TString dataDir = "/alice/sim/2011/LHC11f3";
TString dataPattern = "*AliESDsSignal.root";
TString runFormat = "%d";
TString outDir = "Sim/LHC10h/EmbeddingPbPb276/Phys/sig";
Int_t ttl = 30000;
Int_t maxFilesPerJob = 20;
Int_t maxMergeFiles = 10;
Int_t maxMergeStages = 2;

//______________________________________________________________________________
void runConstrainedJPsi(TString smode = "local", TString inputFileName = "AliESDs.root", Bool_t mc = kTRUE)
{
  /// Compute JPsi kinematics constrained by the JPsi PDG mass
  
  gROOT->LoadMacro("$WORK/Macros/Facilities/runTaskFacilities.C");
  
  // --- Check runing mode ---
  Int_t mode = GetMode(smode, inputFileName);
  if(mode < 0){
    Error("runConstrainedJPsi","Please provide either an ESD root file a collection of ESDs or a dataset.");
    return;
  }
  
  // --- copy files needed for this analysis ---
  TList pathList; pathList.SetOwner();
  pathList.Add(new TObjString("$WORK/Macros/MuonPhysics"));
  TList fileList; fileList.SetOwner();
  fileList.Add(new TObjString("runConstrainedJPsi.C"));
  fileList.Add(new TObjString("AddTaskConstrainedJPsi.C"));
  fileList.Add(new TObjString("AliAnalysisTaskConstrainedJPsi.cxx"));
  fileList.Add(new TObjString("AliAnalysisTaskConstrainedJPsi.h"));
  CopyFileLocally(pathList, fileList);
  
  // --- prepare environment ---
  TString extraLibs="RAWDatabase:RAWDatarec:CDB:STEER:MUONcore:MUONmapping:MUONcalib:MUONgeometry:MUONtrigger:MUONraw:MUONbase:MUONshuttle:MUONrec:RAWDatasim:MUONsim:MUONevaluation:CORRFW:PWGmuon";
  TString extraIncs="include:MUON:MUON/mapping";
  TString extraTasks="AliAnalysisTaskConstrainedJPsi";
  LoadAlirootLocally(extraLibs, extraIncs, extraTasks);
  AliAnalysisGrid *alienHandler = 0x0;
  if (mode == kProof || mode == kProofLite) LoadAlirootOnProof(smode, alirootVersion, extraLibs, extraIncs, extraTasks, kTRUE);
  else if (mode == kGrid || mode == kTerminate) {
    TString analysisMacroName = "ConstrainedJPsi";
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
  AliAnalysisManager *mgr = new AliAnalysisManager("ConstrainedJPsiAnalysis");
  
  // Connect plugin to the analysis manager if any
  if (alienHandler) mgr->SetGridHandler(static_cast<AliAnalysisGrid*>(alienHandler));
  
  // ESD input
  AliESDInputHandler* esdH = new AliESDInputHandler();
  esdH->SetReadFriends(kFALSE);
  esdH->SetInactiveBranches("*");
  esdH->SetActiveBranches("MuonTracks MuonClusters MuonPads AliESDRun. AliESDHeader. AliMultiplicity. AliESDFMD. AliESDVZERO. SPDVertex. PrimaryVertex. AliESDZDC.");
  mgr->SetInputEventHandler(esdH);
  
  // Monte Carlo handler
  if (mc) {
    AliMCEventHandler* mcH = new AliMCEventHandler();
    mgr->SetMCtruthEventHandler(mcH);
  }
  /*
  // event selection
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
  AliPhysicsSelectionTask* physicsSelection = AddTaskPhysicsSelection(mc);
  if(!physicsSelection) {
    Error("CreateAnalysisTrain","AliPhysicsSelectionTask not created!");
    return;
  }
  */
  // track selection
  AliMuonPairCuts pairCuts("stdCuts", "stdCuts");
  pairCuts.SetIsMC();
  pairCuts.SetFilterMask(AliMuonPairCuts::kBothMuEta|AliMuonPairCuts::kBothMuThetaAbs|
			 AliMuonPairCuts::kBothMuPdca|AliMuonPairCuts::kOneMuMatchLpt|
			 AliMuonPairCuts::kDimuRapidity|AliMuonPairCuts::kDimuUnlikeSign);
  pairCuts.ApplySharpPtCutInMatching();
  
  // physics task
  gROOT->LoadMacro("AddTaskConstrainedJPsi.C");
  AliAnalysisTaskConstrainedJPsi* task = AddTaskConstrainedJPsi();
  if(!task) {
    Error("CreateAnalysisTrain","AliAnalysisTaskConstrainedJPsi not created!");
    return;
  }
  task->SetDefaultStorage(Form("local://%s/..",gSystem->pwd()));
  task->SetMuonPairCuts(pairCuts);
  task->SetChResCorrFactor(0.5, 0.25);
  //task->SetVtxZRes(5.);

  
}

