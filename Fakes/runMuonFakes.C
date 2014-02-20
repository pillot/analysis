/*
 *  runMuonFakes.C
 *  aliroot
 *
 *  Created by Philippe Pillot on 31/03/10.
 *  Copyright 2010 SUBATECH. All rights reserved.
 *
 */

TString rootVersion = "v5-33-02a";
TString alirootVersion = "v5-03-15-AN";
TString dataDir = "/alice/sim/LHC11a10a_bis";
TString dataPattern = "*AliESDs.root";
TString runFormat = "%d";
TString outDir = "Sim/LHC11a10a_bis/Fakes/train2";
Int_t ttl = 60000;
Int_t maxFilesPerJob = 100;
Int_t maxMergeFiles = 10;
Int_t maxMergeStages = 3;

//______________________________________________________________________________
void runMuonFakes(TString smode = "local", TString inputFileName = "AliESDs.root",
		  Bool_t useMCLabels = kFALSE, Bool_t combineMCId = kTRUE,
		  Bool_t matchTrig = kTRUE, Bool_t applyAccCut = kTRUE)
{
  /// Study the MUON performances
  
  gROOT->LoadMacro("$WORK/Macros/Facilities/runTaskFacilities.C");
  
  // --- Check runing mode ---
  Int_t mode = GetMode(smode, inputFileName);
  if(mode < 0) {
    Error("runMuonFakes","Please provide either an ESD root file a collection of ESDs or a dataset.");
    return;
  }
  
  // --- copy files needed for this analysis ---
  TList pathList; pathList.SetOwner();
  pathList.Add(new TObjString("$WORK/Macros/Fakes"));
  pathList.Add(new TObjString("$WORK/aliroot/PWGPP/MUON/dep"));
  pathList.Add(new TObjString("$WORK/Macros/MuonPhysics"));
  TList fileList; fileList.SetOwner();
  fileList.Add(new TObjString("runMuonFakes.C"));
  fileList.Add(new TObjString("AddTaskMuonFakes.C"));
  fileList.Add(new TObjString("AliAnalysisTaskMuonFakes.cxx"));
  fileList.Add(new TObjString("AliAnalysisTaskMuonFakes.h"));
  fileList.Add(new TObjString("AddTaskMuonPhysics.C"));
  fileList.Add(new TObjString("AliAnalysisTaskMuonPhysics.cxx"));
  fileList.Add(new TObjString("AliAnalysisTaskMuonPhysics.h"));
  CopyFileLocally(pathList, fileList);
  
  // --- prepare environment ---
  TString extraLibs="RAWDatabase:RAWDatarec:CDB:STEER:MUONcore:MUONmapping:MUONcalib:MUONgeometry:MUONtrigger:MUONraw:MUONbase:MUONshuttle:MUONrec:RAWDatasim:MUONsim:MUONevaluation:CORRFW:PWGmuon";
  TString extraIncs="include:MUON:MUON/mapping:PWG/muon";
//  TString extraTasks="AliAnalysisTaskMuonFakes";
  TString extraTasks="AliAnalysisTaskMuonFakes:AliAnalysisTaskMuonPhysics";
  LoadAlirootLocally(extraLibs, extraIncs, extraTasks);
  AliAnalysisGrid *alienHandler = 0x0;
  if (mode == kProof) LoadAlirootOnProof(smode, alirootVersion, extraLibs, extraIncs, extraTasks, kTRUE);
  else if (mode == kGrid || mode == kTerminate) {
    TString analysisMacroName = "MuonFakesAnalysis";
    alienHandler = static_cast<AliAnalysisGrid*>(CreateAlienHandler(smode, rootVersion, alirootVersion, inputFileName, dataDir, dataPattern, outDir, extraLibs, extraIncs, extraTasks, analysisMacroName, runFormat, ttl, maxFilesPerJob, maxMergeFiles, maxMergeStages));
    if (!alienHandler) return;
  }
  
  // --- Create the analysis train ---
  CreateAnalysisTrain(useMCLabels, combineMCId, matchTrig, applyAccCut, alienHandler);
  
  // --- Create input object ---
  TObject* inputObj = CreateInputObject(mode, inputFileName);
  
  // --- start analysis ---
  StartAnalysis(mode, inputObj);
  
  // --- save summary canvases ---
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

//______________________________________________________________________________
void CreateAnalysisTrain(Bool_t useMCLabels, Bool_t combineMCId, Bool_t matchTrig, Bool_t applyAccCut, TObject* alienHandler)
{
  /// create the analysis train and configure it
  
  // analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("MuonFakesAnalysis");
  
  // Connect plugin to the analysis manager if any
  if (alienHandler) mgr->SetGridHandler(static_cast<AliAnalysisGrid*>(alienHandler));
  
  // ESD input
  AliESDInputHandler* esdH = new AliESDInputHandler();
  esdH->SetReadFriends(kFALSE);
  esdH->SetInactiveBranches("*");
  esdH->SetActiveBranches("MuonTracks MuonClusters MuonPads AliESDRun. AliESDHeader. AliMultiplicity. AliESDFMD. AliESDVZERO. SPDVertex. PrimaryVertex. AliESDZDC.");
  mgr->SetInputEventHandler(esdH);
  
  // Monte Carlo handler
  AliMCEventHandler* mcHandler = new AliMCEventHandler();
  mgr->SetMCtruthEventHandler(mcHandler);
  
  // event selection
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
  AliPhysicsSelectionTask* physicsSelection = AddTaskPhysicsSelection(kTRUE);
  if (!physicsSelection) {
    Error("CreateAnalysisTrain","AliPhysicsSelectionTask not created!");
    return;
  }
  
  // centrality selection
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskCentrality.C");
  AliCentralitySelectionTask* centralityTask = AddTaskCentrality();
  if (!centralityTask) {
    Error("CreateAnalysisTrain","AliCentralitySelectionTask not created!");
    return;
  }
  centralityTask->SetMCInput();
  
  // ---------------- single task ----------------
  
  // Fakes analysis
  gROOT->LoadMacro("AddTaskMuonFakes.C");
  TString extension = useMCLabels ? "Label" : "PosOCDB";
  if (combineMCId) extension = "CombiOCDB";
  AliAnalysisTaskMuonFakes* muonFakes = AddTaskMuonFakes(useMCLabels, combineMCId, matchTrig, applyAccCut, extension);
  if(!muonFakes) {
    Error("CreateAnalysisTrain","AliAnalysisTaskMuonFakes not created!");
    return;
  }
  //muonFakes->ShowProgressBar();
  //muonFakes->RecoParamLocation("local://$ALICE_ROOT/OCDB");
  //muonFakes->DecayAsFake();
  
  // track selection
  AliMuonTrackCuts trackCuts("stdCuts", "stdCuts");
  trackCuts.SetAllowDefaultParams();
  trackCuts.SetFilterMask(AliMuonTrackCuts::kMuMatchLpt | AliMuonTrackCuts::kMuEta |
			  AliMuonTrackCuts::kMuThetaAbs | AliMuonTrackCuts::kMuPdca);
  trackCuts.SetIsMC(kTRUE);
  
  // Physics results
  gROOT->LoadMacro("AddTaskMuonPhysics.C");
  AliAnalysisTaskMuonPhysics* physics = AddTaskMuonPhysics();
  if(!physics) {
    Error("CreateAnalysisTrain","AliAnalysisTaskMuonPhysics not created!");
    return;
  }
  physics->SetMuonTrackCuts(trackCuts);
  
  // ---------------- train ----------------
  /*
  // Fakes analysis 1 (with label)
  gROOT->LoadMacro("AddTaskMuonFakes.C");
  AliAnalysisTaskMuonFakes* muonFakes = AddTaskMuonFakes(kTRUE, kFALSE, matchTrig, applyAccCut, "Label");
  if(!muonFakes) {
    Error("CreateAnalysisTrain","AliAnalysisTaskMuonFakes not created!");
    return;
  }
  //muonFakes->ShowProgressBar();
  
  // Fakes analysis 2 (position with sigma cut from OCDB)
  muonFakes = AddTaskMuonFakes(kFALSE, kFALSE, matchTrig, applyAccCut, "PosOCDB");
  if(!muonFakes) {
    Error("CreateAnalysisTrain","AliAnalysisTaskMuonFakes not created!");
    return;
  }
  
  // Fakes analysis 3 (combining label and position with sigma cut from OCDB)
  muonFakes = AddTaskMuonFakes(kFALSE, kTRUE, matchTrig, applyAccCut, "CombiOCDB");
  if(!muonFakes) {
    Error("CreateAnalysisTrain","AliAnalysisTaskMuonFakes not created!");
    return;
  }
  
  // Fakes analysis 4 (position with sigma cut = 6)
  muonFakes = AddTaskMuonFakes(kFALSE, kFALSE, matchTrig, applyAccCut, "Pos6Sigma");
  if(!muonFakes) {
    Error("CreateAnalysisTrain","AliAnalysisTaskMuonFakes not created!");
    return;
  }
  muonFakes->SetExternalSigmaCut(6.);
  
  // Fakes analysis 5 (combining label and position with sigma cut = 6)
  muonFakes = AddTaskMuonFakes(kFALSE, kTRUE, matchTrig, applyAccCut, "Combi6Sigma");
  if(!muonFakes) {
    Error("CreateAnalysisTrain","AliAnalysisTaskMuonFakes not created!");
    return;
  }
  muonFakes->SetExternalSigmaCut(6.);
  
  // Fakes analysis 6 (position with sigma cut = 8)
  muonFakes = AddTaskMuonFakes(kFALSE, kFALSE, matchTrig, applyAccCut, "Pos8Sigma");
  if(!muonFakes) {
    Error("CreateAnalysisTrain","AliAnalysisTaskMuonFakes not created!");
    return;
  }
  muonFakes->SetExternalSigmaCut(8.);
  
  // Fakes analysis 7 (combining label and position with sigma cut = 8)
  muonFakes = AddTaskMuonFakes(kFALSE, kTRUE, matchTrig, applyAccCut, "Combi8Sigma");
  if(!muonFakes) {
    Error("CreateAnalysisTrain","AliAnalysisTaskMuonFakes not created!");
    return;
  }
  muonFakes->SetExternalSigmaCut(8.);
  
  // Fakes analysis 8 (position with sigma cut = 10)
  muonFakes = AddTaskMuonFakes(kFALSE, kFALSE, matchTrig, applyAccCut, "Pos10Sigma");
  if(!muonFakes) {
    Error("CreateAnalysisTrain","AliAnalysisTaskMuonFakes not created!");
    return;
  }
  muonFakes->SetExternalSigmaCut(10.);
  
  // Fakes analysis 9 (combining label and position with sigma cut = 10)
  muonFakes = AddTaskMuonFakes(kFALSE, kTRUE, matchTrig, applyAccCut, "Combi10Sigma");
  if(!muonFakes) {
    Error("CreateAnalysisTrain","AliAnalysisTaskMuonFakes not created!");
    return;
  }
  muonFakes->SetExternalSigmaCut(10.);
  
  // Fakes analysis 10 (combining label and position with sigma cut from OCDB + pt > 2 GeV/c)
  muonFakes = AddTaskMuonFakes(kFALSE, kTRUE, matchTrig, applyAccCut, "CombiOCDBpT2GeV");
  if(!muonFakes) {
    Error("CreateAnalysisTrain","AliAnalysisTaskMuonFakes not created!");
    return;
  }
  muonFakes->PtCut(2.);
  
  // Fakes analysis 11 (combining label and position with sigma cut = 10 + pt > 2 GeV/c)
  muonFakes = AddTaskMuonFakes(kFALSE, kTRUE, matchTrig, applyAccCut, "Combi10SigmapT2GeV");
  if(!muonFakes) {
    Error("CreateAnalysisTrain","AliAnalysisTaskMuonFakes not created!");
    return;
  }
  muonFakes->SetExternalSigmaCut(10.);
  muonFakes->PtCut(2.);
  
  // Fakes analysis 12 (combining label and position with sigma cut from OCDB + chi2 < 0.5)
  muonFakes = AddTaskMuonFakes(kFALSE, kTRUE, matchTrig, applyAccCut, "CombiOCDBchi205");
  if(!muonFakes) {
    Error("CreateAnalysisTrain","AliAnalysisTaskMuonFakes not created!");
    return;
  }
  muonFakes->Chi2Cut(0.5);
  
  // Fakes analysis 13 (combining label and position with sigma cut = 10 + chi2 < 0.5)
  muonFakes = AddTaskMuonFakes(kFALSE, kTRUE, matchTrig, applyAccCut, "Combi10Sigmachi205");
  if(!muonFakes) {
    Error("CreateAnalysisTrain","AliAnalysisTaskMuonFakes not created!");
    return;
  }
  muonFakes->SetExternalSigmaCut(10.);
  muonFakes->Chi2Cut(0.5);
  
  // Fakes analysis 14 (combining label and position with sigma cut from OCDB + pDCA)
  AliMuonTrackCuts trackCuts("stdCuts", "stdCuts");
  //trackCuts.SetAllowDefaultParams();
  trackCuts.SetIsMC();
  trackCuts.SetFilterMask(AliMuonTrackCuts::kMuPdca);
  muonFakes = AddTaskMuonFakes(kFALSE, kTRUE, matchTrig, applyAccCut, "CombiOCDBpDCA");
  if(!muonFakes) {
    Error("CreateAnalysisTrain","AliAnalysisTaskMuonFakes not created!");
    return;
  }
  muonFakes->SetMuonTrackCuts(trackCuts);
  
  // Fakes analysis 15 (combining label and position with sigma cut = 10 + pDCA)
  muonFakes = AddTaskMuonFakes(kFALSE, kTRUE, matchTrig, applyAccCut, "Combi10SigmapDCA");
  if(!muonFakes) {
    Error("CreateAnalysisTrain","AliAnalysisTaskMuonFakes not created!");
    return;
  }
  muonFakes->SetExternalSigmaCut(10.);
  muonFakes->SetMuonTrackCuts(trackCuts);
  */
}

