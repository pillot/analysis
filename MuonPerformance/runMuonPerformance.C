/*
 *  runMuonPerformance.C
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
TString outDir = "Sim/LHC11h/muLizardo/Perf";
Int_t ttl = 30000;
Int_t maxFilesPerJob = 10;
Int_t maxMergeFiles = 10;
Int_t maxMergeStages = 2;

//TString alignStorage = "alien://folder=/alice/simulation/2008/v4-15-Release/Residual";
TString alignStorage = "alien://folder=/alice/cern.ch/user/p/ppillot/OCDB_PbPbSim";

//______________________________________________________________________________
void runMuonPerformance(TString smode = "local", TString inputFileName = "AliESDs.root",
			Bool_t correctClusterResForSystematics = kTRUE, Bool_t fitClusterResiduals = kTRUE,
			Bool_t applyPhysicsSelection = kFALSE, Bool_t embedding = kFALSE,
                        Bool_t refit = kFALSE, Bool_t runFakeAnalysis = kFALSE)
{
  /// Study the MUON performances
  
  gROOT->LoadMacro("$WORK/Macros/Facilities/runTaskFacilities.C");
  
  // --- Check runing mode ---
  Int_t mode = GetMode(smode, inputFileName);
  if(mode < 0) {
    Error("runMuonPerformance","Please provide either an ESD root file a collection of ESDs or a dataset.");
    return;
  }
  
  // --- copy files needed for this analysis ---
  TList pathList; pathList.SetOwner();
  pathList.Add(new TObjString("$WORK/Macros/MuonPerformance"));
  pathList.Add(new TObjString("$DEV/aliroot/PWG/muondep"));
  pathList.Add(new TObjString("$DEV/aliroot/PWGPP/MUON/dep"));
  pathList.Add(new TObjString("$WORK/Macros/MuonPhysics"));
//  pathList.Add(new TObjString("$WORK/Macros/Sim"));
  TList fileList; fileList.SetOwner();
  fileList.Add(new TObjString("runMuonPerformance.C"));
  if (refit) {
    fileList.Add(new TObjString("AddTaskMuonRefit.C"));
    fileList.Add(new TObjString("AliAnalysisTaskMuonRefit.cxx"));
    fileList.Add(new TObjString("AliAnalysisTaskMuonRefit.h"));
  }
  if (runFakeAnalysis) {
    fileList.Add(new TObjString("AddTaskMuonFakes.C"));
    fileList.Add(new TObjString("AliAnalysisTaskMuonFakes.cxx"));
    fileList.Add(new TObjString("AliAnalysisTaskMuonFakes.h"));
  }
  fileList.Add(new TObjString("AddTaskESDMCLabelAddition.C"));
  fileList.Add(new TObjString("AliAnalysisTaskESDMCLabelAddition.cxx"));
  fileList.Add(new TObjString("AliAnalysisTaskESDMCLabelAddition.h"));
  fileList.Add(new TObjString("AddTaskMuonPerformance.C"));
  fileList.Add(new TObjString("AliAnalysisTaskMuonPerformance.cxx"));
  fileList.Add(new TObjString("AliAnalysisTaskMuonPerformance.h"));
  fileList.Add(new TObjString("AddTaskMUONTrackingEfficiency.C"));
  fileList.Add(new TObjString("AliAnalysisTaskMuonTrackingEff.cxx"));
  fileList.Add(new TObjString("AliAnalysisTaskMuonTrackingEff.h"));
  fileList.Add(new TObjString("AddTaskMUONTrackingEfficiency_old.C"));
  fileList.Add(new TObjString("AliAnalysisTaskMuonTrackingEff_old.cxx"));
  fileList.Add(new TObjString("AliAnalysisTaskMuonTrackingEff_old.h"));
  fileList.Add(new TObjString("AddTaskMuonPhysics.C"));
  fileList.Add(new TObjString("AliAnalysisTaskMuonPhysics.cxx"));
  fileList.Add(new TObjString("AliAnalysisTaskMuonPhysics.h"));
//  fileList.Add(new TObjString("AliMCMuonEventHandler.cxx"));
//  fileList.Add(new TObjString("AliMCMuonEventHandler.h"));
  CopyFileLocally(pathList, fileList);
  
  // --- prepare environment ---
  TString extraLibs="RAWDatabase:RAWDatarec:CDB:STEER:MUONcore:MUONmapping:MUONcalib:MUONgeometry:MUONtrigger:MUONraw:MUONbase:MUONshuttle:MUONrec:RAWDatasim:MUONsim:MUONevaluation:CORRFW:PWGmuon";
  TString extraIncs="include:MUON:MUON/mapping:PWG/muon";
  TString extraTasks = "AliAnalysisTaskESDMCLabelAddition:AliAnalysisTaskMuonPerformance:AliAnalysisTaskMuonTrackingEff:AliAnalysisTaskMuonTrackingEff_old:AliAnalysisTaskMuonPhysics";
  if (refit) extraTasks += ":AliAnalysisTaskMuonRefit";
  if (runFakeAnalysis) extraTasks += ":AliAnalysisTaskMuonFakes";
//  TString extraTasks="AliAnalysisTaskMuonPerformance:AliAnalysisTaskMuonTrackingEff:AliAnalysisTaskMuonPhysics:AliMCMuonEventHandler";
  LoadAlirootLocally(extraLibs, extraIncs, extraTasks);
  AliAnalysisGrid *alienHandler = 0x0;
  if (mode == kProof) LoadAlirootOnProof(smode, rootVersion, alirootVersion, extraLibs, extraIncs, extraTasks, kTRUE);
  else if (mode == kGrid || mode == kTerminate) {
    TString analysisMacroName = "Perf";
    alienHandler = static_cast<AliAnalysisGrid*>(CreateAlienHandler(smode, rootVersion, alirootVersion, inputFileName, dataDir, dataPattern, outDir, extraLibs, extraIncs, extraTasks, analysisMacroName, runFormat, ttl, maxFilesPerJob, maxMergeFiles, maxMergeStages));
    if (!alienHandler) return;
  }
  
  // --- Create the analysis train ---
  CreateAnalysisTrain(correctClusterResForSystematics, fitClusterResiduals, applyPhysicsSelection, embedding, refit, runFakeAnalysis, alienHandler);
  
  // --- Create input object ---
  TObject* inputObj = CreateInputObject(mode, inputFileName);
  
  // --- start analysis ---
  StartAnalysis(mode, inputObj);
  
  // --- save summary canvases ---
  if (runFakeAnalysis) {
    
    TString outFileName = AliAnalysisManager::GetCommonFileName();
    TFile *outFile = (TFile*)gROOT->GetListOfFiles()->FindObject(outFileName.Data());
    if (outFile) outFile->ReOpen("UPDATE");
    else outFile = TFile::Open(outFileName.Data(),"UPDATE");
    if (outFile && outFile->IsOpen()) {
      
      TObjArray *listOfTasks = AliAnalysisManager::GetAnalysisManager()->GetTasks();
      TIter nextTask(listOfTasks);
      AliAnalysisTask *task = 0x0;
      while ((task = static_cast<AliAnalysisTask*>(nextTask()))) {
        
        TString extention = task->GetName();
        if (extention.Contains("MUONFakes")) extention.ReplaceAll("MUONFakes","");
        else continue;
        
        AliAnalysisTaskMuonFakes *fakes = static_cast<AliAnalysisTaskMuonFakes*>(task);
        if (fakes->GetCanvases()) {
          outFile->Cd(Form("%s:/MUON_Fakes_%s",outFileName.Data(),extention.Data()));
          fakes->GetCanvases()->Write(0x0, TObject::kOverwrite);
        }
        
      }
      
      outFile->Close();
    }
    
  }
  
}

//______________________________________________________________________________
void CreateAnalysisTrain(Bool_t correctClusterResForSystematics, Bool_t fitClusterResiduals,
			 Bool_t applyPhysicsSelection, Bool_t embedding, Bool_t refit,
                         Bool_t runFakeAnalysis, TObject* alienHandler)
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
  esdH->SetActiveBranches("MuonTracks MuonClusters MuonPads AliESDRun. AliESDHeader. AliMultiplicity. AliESDFMD. AliESDVZERO. SPDVertex. PrimaryVertex. AliESDZDC. AliESDTZERO.");
  mgr->SetInputEventHandler(esdH);
  
  // Monte Carlo handler
  AliMCEventHandler* mcH = new AliMCEventHandler();
//  AliMCEventHandler* mcH = new AliMCMuonEventHandler();
  mgr->SetMCtruthEventHandler(mcH);
  
  // event selection
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
  AliPhysicsSelectionTask* physicsSelection = embedding ? AddTaskPhysicsSelection() : AddTaskPhysicsSelection(kTRUE);
  if(!physicsSelection) {
    Error("CreateAnalysisTrain","AliPhysicsSelectionTask not created!");
    return;
  }
  UInt_t offlineTriggerMask = AliVEvent::kMB;
  
  // centrality selection
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskCentrality.C");
  AliCentralitySelectionTask *taskCentrality = AddTaskCentrality();
  if(!taskCentrality) {
    Error("CreateAnalysisTrain","AliCentralitySelectionTask not created!");
    return;
  }
  if (!embedding) taskCentrality->SetMCInput();
  
  // track selection
  AliMuonTrackCuts trackCuts("stdCuts", "stdCuts");
  trackCuts.SetAllowDefaultParams();
  //trackCuts.SetPassNumber(2);
  trackCuts.SetFilterMask(AliMuonTrackCuts::kMuMatchLpt | AliMuonTrackCuts::kMuEta |
			  AliMuonTrackCuts::kMuThetaAbs);
  if (!embedding) trackCuts.SetIsMC(kTRUE);
  
  // track refitting
  if (refit) {
    gROOT->LoadMacro("AddTaskMuonRefit.C");
    AliAnalysisTaskMuonRefit* refit = AddTaskMuonRefit(-1., -1., kTRUE, -1., -1.);
    if(!refit) {
      Error("CreateAnalysisTrain","AliAnalysisTaskMuonRefit not created!");
      return;
    }
    if (applyPhysicsSelection) refit->SelectCollisionCandidates(offlineTriggerMask);
    if (!alignStorage.IsNull()) refit->SetAlignStorage(alignStorage.Data());
    refit->RemoveMonoCathodClusters(kTRUE, kFALSE);
  }
  
  // Fakes analysis
  if (runFakeAnalysis) {
    gROOT->LoadMacro("AddTaskMuonFakes.C");
    AliAnalysisTaskMuonFakes* muonFakes = AddTaskMuonFakes(kFALSE, kTRUE, kTRUE, kTRUE, "CombiOCDB");
    if(!muonFakes) {
      Error("CreateAnalysisTrain","AliAnalysisTaskMuonFakes not created!");
      return;
    }
    if (applyPhysicsSelection) muonFakes->SelectCollisionCandidates(offlineTriggerMask);
    muonFakes->DisableDetailedCounters();
    muonFakes->SetMuonTrackCuts(trackCuts);
    //  muonFakes->ShowProgressBar();
  }
  
  // Add MC labels
  gROOT->LoadMacro("AddTaskESDMCLabelAddition.C");
  AliAnalysisTaskESDMCLabelAddition *esdMCLabeltask = AddTaskESDMCLabelAddition();
  if(!esdMCLabeltask) {
    Error("CreateAnalysisTrain","AliAnalysisTaskESDMCLabelAddition not created!");
    return;
  }
  if (applyPhysicsSelection) esdMCLabeltask->SelectCollisionCandidates(offlineTriggerMask);
  if (!alignStorage.IsNull()) esdMCLabeltask->SetAlignStorage(alignStorage.Data());
  
  // Muon performance analysis
  gROOT->LoadMacro("AddTaskMuonPerformance.C");
  AliAnalysisTaskMuonPerformance* muonPerformance = AddTaskMuonPerformance(correctClusterResForSystematics, fitClusterResiduals);
  if(!muonPerformance) {
    Error("CreateAnalysisTrain","AliAnalysisTaskMuonPerformance not created!");
    return;
  }
  if (applyPhysicsSelection) muonPerformance->SelectCollisionCandidates(offlineTriggerMask);
//  muonPerformance->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  if (!alignStorage.IsNull()) muonPerformance->SetAlignStorage(alignStorage.Data());
  //muonPerformance->SetPBins(15, 0., 100.);
  //muonPerformance->EnforceTrackingCriteria(kTRUE);
  muonPerformance->UseMCKinematics(kTRUE);
  muonPerformance->SetMCTrigLevelFromMatchTrk(kTRUE);
  
  // Muon efficiency analysis (old without cut but with MC label)
  gROOT->LoadMacro("AddTaskMUONTrackingEfficiency_old.C");
  AliAnalysisTaskMuonTrackingEff_old* muonEfficiency_old_wocut_wMClabel = AddTaskMUONTrackingEfficiency_old(kFALSE, kFALSE,"old_wocut_wMClabel");
  if(!muonEfficiency_old_wocut_wMClabel) {
    Error("CreateAnalysisTrain","AliAnalysisTaskMuonTrackingEff_old not created!");
    return;
  }
  if (applyPhysicsSelection) muonEfficiency_old_wocut_wMClabel->SelectCollisionCandidates(offlineTriggerMask);
  //muonEfficiency_old_wocut_wMClabel->PtCut(1.);
  muonEfficiency_old_wocut_wMClabel->UseMCLabel(kTRUE);
  
  // Muon efficiency analysis (old with cut)
  AliAnalysisTaskMuonTrackingEff_old* muonEfficiency_old_wcut = AddTaskMUONTrackingEfficiency_old(kTRUE, kTRUE,"old_wcut");
  if(!muonEfficiency_old_wcut) {
    Error("CreateAnalysisTrain","AliAnalysisTaskMuonTrackingEff_old not created!");
    return;
  }
  if (applyPhysicsSelection) muonEfficiency_old_wcut->SelectCollisionCandidates(offlineTriggerMask);
  muonEfficiency_old_wcut->PtCut(1.);
  
  // Muon efficiency analysis (old with cut and MC label)
  AliAnalysisTaskMuonTrackingEff_old* muonEfficiency_old_wcut_wMClabel = AddTaskMUONTrackingEfficiency_old(kTRUE, kTRUE,"old_wcut_wMClabel");
  if(!muonEfficiency_old_wcut_wMClabel) {
    Error("CreateAnalysisTrain","AliAnalysisTaskMuonTrackingEff_old not created!");
    return;
  }
  if (applyPhysicsSelection) muonEfficiency_old_wcut_wMClabel->SelectCollisionCandidates(offlineTriggerMask);
  muonEfficiency_old_wcut_wMClabel->PtCut(1.);
  muonEfficiency_old_wcut_wMClabel->UseMCLabel(kTRUE);
  
  // Muon efficiency analysis
  gROOT->LoadMacro("AddTaskMUONTrackingEfficiency.C");
  AliAnalysisTaskMuonTrackingEff* muonEfficiency = AddTaskMUONTrackingEfficiency(trackCuts,"");
  if (applyPhysicsSelection) muonEfficiency->SelectCollisionCandidates(offlineTriggerMask);
  if (!alignStorage.IsNull()) muonEfficiency->SetAlignStorage(alignStorage.Data());
  //muonEfficiency->SetRecoParamStorage("alien://folder=/alice/cern.ch/user/p/ppillot/OCDB2012_newReco");
  muonEfficiency->SetMuonPtCut(1.);
  muonEfficiency->UseMCLabel(kTRUE);
  muonEfficiency->EnableDisplay(kFALSE);
  
  // Physics results
  gROOT->LoadMacro("AddTaskMuonPhysics.C");
  AliAnalysisTaskMuonPhysics* physics = AddTaskMuonPhysics();
  if(!physics) {
    Error("CreateAnalysisTrain","AliAnalysisTaskMuonPhysics not created!");
    return;
  }
  if (applyPhysicsSelection) physics->SelectCollisionCandidates(offlineTriggerMask);
  physics->SetMuonTrackCuts(trackCuts);
}

