/*
 *  run.C
 *  aliroot_dev
 *
 *  Created by philippe pillot on 28/07/11.
 *  Copyright 2011 SUBATECH. All rights reserved.
 *
 */


TString rootVersion = "v5-28-00e";
TString alirootVersion = "v4-21-31-AN";
TString dataDir = "/alice/data/2010/LHC10h";
TString dataPattern = "pass2/*AliESDs.root";
TString runFormat = "%09d";
TString outDir = "PbPb2.76TeV/LHC10h/pass2/MuonEfficiency/pT2GeV2";
Int_t ttl = 30000;
Int_t maxFilesPerJob = 100;
Int_t maxMergeFiles = 10;
Int_t maxMergeStages = 2;

//______________________________________________________________________________
void run(TString smode = "local", TString inputFileName = "AliESDs.root",
	 TString outputFileName = "AnalysisResults.root", char overwrite = '\0')
{
  /// Study the MUON performances
  
  gROOT->LoadMacro("$ALICE/Macros/Facilities/runTaskFacilities.C");
  
  // --- Check runing mode ---
  Int_t mode = GetMode(smode, inputFileName);
  if(mode < 0) {
    Error("run","Please provide either an ESD root file a collection of ESDs or a dataset.");
    return;
  }
  
  // --- copy files needed for this analysis ---
  TList pathList; pathList.SetOwner();
  pathList.Add(new TObjString("$ALICE/Macros/Laurent"));
  pathList.Add(new TObjString("$ALICE_ROOT/PWG3/muon"));
  TList fileList; fileList.SetOwner();
  fileList.Add(new TObjString("run.C"));
  fileList.Add(new TObjString("AliHistogramCollection.cxx"));
  fileList.Add(new TObjString("AliHistogramCollection.h"));
  fileList.Add(new TObjString("AliAnalysisMuMu.cxx"));
  fileList.Add(new TObjString("AliAnalysisMuMu.h"));
  fileList.Add(new TObjString("AliAnalysisMuMuFromAOD.cxx"));
  fileList.Add(new TObjString("AliAnalysisMuMuFromAOD.h"));
  CopyFileLocally(pathList, fileList, overwrite);
  
  // --- prepare environment ---
  TString extraLibs="CORRFW:PWG3base";
  TString extraIncs="include:PWG3/base";
  TString extraTasks="AliHistogramCollection:AliAnalysisMuMu:AliAnalysisMuMuFromAOD";
  LoadAlirootLocally(extraLibs, extraIncs, extraTasks);
  AliAnalysisGrid *alienHandler = 0x0;
  if (mode == kProof) LoadAlirootOnProof(smode, alirootVersion, extraLibs, extraIncs, extraTasks, kTRUE);
  else if (mode == kGrid || mode == kTerminate) {
    TString analysisMacroName = "AODAnalysis";
    alienHandler = static_cast<AliAnalysisGrid*>(CreateAlienHandler(smode, rootVersion, alirootVersion, inputFileName, dataDir, dataPattern, outDir, extraLibs, extraIncs, extraTasks, analysisMacroName, runFormat, ttl, maxFilesPerJob, maxMergeFiles, maxMergeStages));
    if (!alienHandler) return;
  }
  
  // --- Create the analysis train ---
  CreateAnalysisTrain(alienHandler, outputFileName);
  
  // --- Create input object ---
  TObject* inputObj = CreateInputObject(mode, inputFileName);
  
  // --- start analysis ---
  StartAnalysis(mode, inputObj);
  
}

//______________________________________________________________________________
void CreateAnalysisTrain(TObject* alienHandler, TString outputFileName)
{
  /// create the analysis train and configure it
  
  // analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("AODAnalysis");
  
  // Connect plugin to the analysis manager if any
  if (alienHandler) mgr->SetGridHandler(static_cast<AliAnalysisGrid*>(alienHandler));
  
  // AOD input
  AliInputEventHandler* input = new AliAODInputHandler;
  mgr->SetInputEventHandler(input);
  
  // AOD analysis
  TList triggers;
  triggers.SetOwner(kTRUE);
  triggers.Add(new TObjString("ANY"));
  triggers.Add(new TObjString("CMBAC-B-NOPF-ALL"));
  triggers.Add(new TObjString("CMBACS2-B-NOPF-ALL"));
  triggers.Add(new TObjString("CMBACS2-B-NOPF-ALLNOTRD"));
  AddTaskMuMuFromAOD(outputFileName.Data(),&triggers);
  
}

