//
//  runAODtrainCustom.C
//  aliphysics-dev
//
//  Created by philippe pillot on 12/02/2016.
//  Copyright © 2016 Philippe Pillot. All rights reserved.
//


//______________________________________________________________________________
void runAODtrainCustom(TString smode = "local", TString inputFileName = "AliESDs.root")
{
  /// run the AOD train using the plugin and getting the wagons from AODtrainCustom.C
  
  gROOT->LoadMacro("$HOME/Work/Alice/Macros/Facilities/runTaskFacilities.C");
  
  // load the AODtrainCustom.C macro to get the setting and know what else to load
  CopyInputFileLocally("$WORK/Macros/Sim/AODtrainCustom.C", "AODtrainCustom.C");
  gROOT->LoadMacro("AODtrainCustom.C");
  
  // --- general analysis setup ---
  TString rootVersion = "";
  TString alirootVersion = "";
  TString aliphysicsVersion = "vAN-20160212-1";
  TString extraLibs="";
  TString extraIncs="include";
  TString extraTasks="";
  TString extraPkgs="";
  Bool_t needPWGmuondep = kFALSE;
  Bool_t needPWGPPMUONlite = kFALSE;
  Bool_t needPWGPPMUONdep = kFALSE;
  TList pathList; pathList.SetOwner();
  pathList.Add(new TObjString("$WORK/Macros/Sim/AOD"));
  TList fileList; fileList.SetOwner();
  fileList.Add(new TObjString("runAODtrainCustom.C"));
  if (iMUONCDBConnect > 1) {
    fileList.Add(new TObjString("AddTaskMuonCDBConnect.C"));
    needPWGmuondep = kTRUE;
  }
  if (iESDMCLabelAddition > 1) {
    fileList.Add(new TObjString("AddTaskESDMCLabelAddition.C"));
    needPWGmuondep = kTRUE;
  }
  if (iMUONRefit > 1) {
    fileList.Add(new TObjString("AddTaskMuonRefit.C"));
    needPWGmuondep = kTRUE;
  }
  if (iMUONRefitVtx > 1) {
    fileList.Add(new TObjString("AddTaskMuonRefitVtx.C"));
    needPWGmuondep = kTRUE;
  }
  if (iMUONQA > 1) {
    fileList.Add(new TObjString("AddTaskMuonQA.C"));
    needPWGPPMUONlite = kTRUE;
  }
  if (useMC && useTR && iMUONPerformance > 1) {
    fileList.Add(new TObjString("AddTaskMuonPerformance.C"));
    needPWGPPMUONdep = kTRUE;
  }
  if (iMUONEfficiency > 1) {
    fileList.Add(new TObjString("AddTaskMUONTrackingEfficiency.C"));
    needPWGPPMUONdep = kTRUE;
  }
  if (needPWGmuondep || needPWGPPMUONlite || needPWGPPMUONdep) {
    pathList.Add(new TObjString("$ALICE_PHYSICS/PARfiles"));
  }
  if (needPWGmuondep) {
    extraPkgs += extraPkgs.IsNull() ? "PWGmuondep" : ":PWGmuondep";
    pathList.Add(new TObjString("$DEVPHYS/PWG/muondep"));
    fileList.Add(new TObjString("PWGmuondep.par"));
  }
  if (needPWGPPMUONlite) {
    extraPkgs += extraPkgs.IsNull() ? "PWGPPMUONlite" : ":PWGPPMUONlite";
    pathList.Add(new TObjString("$DEVPHYS/PWGPP/MUON/lite"));
    fileList.Add(new TObjString("PWGPPMUONlite.par"));
  }
  if (needPWGPPMUONdep) {
    extraPkgs += extraPkgs.IsNull() ? "PWGPPMUONdep" : ":PWGPPMUONdep";
    pathList.Add(new TObjString("$DEVPHYS/PWGPP/MUON/dep"));
    fileList.Add(new TObjString("PWGPPMUONdep.par"));
  }
  if (iMUONPhysics > 0) {
    extraTasks = "AliAnalysisTaskMuonPhysics";
    pathList.Add(new TObjString("$WORK/Macros/MuonPhysics"));
    fileList.Add(new TObjString("AddTaskMuonPhysics.C"));
    fileList.Add(new TObjString("AliAnalysisTaskMuonPhysics.cxx"));
    fileList.Add(new TObjString("AliAnalysisTaskMuonPhysics.h"));
  }
  
  // --- grid specific setup ---
  TString dataDir = "/alice/cern.ch/user/p/ppillot/Sim/LHC15n/muTune2VtxDCA/results";
  TString dataPattern = "*AliESDs.root";
  TString runFormat = "%d";
  TString outDir = "Sim/LHC15n/muTune2VtxDCA/SpectroShift";
  TString analysisMacroName = "AOD";
  Int_t ttl = 30000;
  Int_t maxFilesPerJob = 10;
  Int_t maxMergeFiles = 10;
  Int_t maxMergeStages = 2;
  TString RegisterExcludes = "AliAOD.root pyxsec_hists.root";
  
  // --- prepare the analysis environment ---
  Int_t mode = PrepareAnalysis(smode, inputFileName, extraLibs, extraIncs, extraTasks, extraPkgs, pathList, fileList);
  
  // --- run the analysis (saf3 is a special case as the analysis is launched on the server) ---
  if (mode == kSAF3Connect) {
    
    if (!RunAnalysisOnSAF3(fileList, aliphysicsVersion, inputFileName)) return;
    
  } else {
    
    CreateAnalysisTrain();
    
    RunAnalysis(smode, inputFileName, rootVersion, alirootVersion, aliphysicsVersion, extraLibs, extraIncs, extraTasks, extraPkgs, dataDir, dataPattern, outDir, analysisMacroName, runFormat, ttl, maxFilesPerJob, maxMergeFiles, maxMergeStages, RegisterExcludes);
    
  }
  
}

//______________________________________________________________________________
void CreateAnalysisTrain()
{
  /// create the analysis train and configure it
  
  // analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("AODtrainCustom");
  
  // ESD input
  AliESDInputHandler* esdH = new AliESDInputHandler();
  esdH->SetReadFriends(kFALSE);
  mgr->SetInputEventHandler(esdH);
  
  // Monte Carlo handler
  AliMCEventHandler* mcH = new AliMCEventHandler();
  mgr->SetMCtruthEventHandler(mcH);
  mcH->SetReadTR(kTRUE);
  
  // AOD output handler (AOD output container created automatically when setting an AOD handler)
  AliAODHandler* aodH = new AliAODHandler();
  aodH->SetOutputFileName("AliAOD.root");
  mgr->SetOutputEventHandler(aodH);
  
  // configure the train from AODtrainCustom.C
  AddAnalysisTasks(0);
  
}

