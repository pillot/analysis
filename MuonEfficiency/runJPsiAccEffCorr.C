/*
 *  runJPsiAccEffCorr.C
 *  aliroot_dev
 *
 *  Created by philippe pillot on 15/11/11.
 *  Copyright 2011 SUBATECH. All rights reserved.
 *
 */


TString rootVersion = "v5-34-02-1";
TString alirootVersion = "v5-04-10-AN";
TString dataDir = "/alice/data/2010/LHC10h";
TString dataPattern = "pass2/*AliESDs.root";
TString runFormat = "%09d";
TString outDir = "Data/LHC10h/pass2/Eff/pDCAChi2";
Int_t ttl = 30000;
Int_t maxFilesPerJob = 100;
Int_t maxMergeFiles = 10;
Int_t maxMergeStages = 2;

//______________________________________________________________________________
void runJPsiAccEffCorr(TString smode = "local", TString inputFileName = "AliAOD.root")
{
  /// Study the MUON performances
  
  gROOT->LoadMacro("$WORK/Macros/Facilities/runTaskFacilities.C");
  
  // --- Check runing mode ---
  Int_t mode = GetMode(smode, inputFileName);
  if(mode < 0) {
    Error("runJPsiAccEffCorr","Please provide either an ESD root file a collection of ESDs or a dataset.");
    return;
  }
  
  // --- copy files needed for this analysis ---
  TList pathList; pathList.SetOwner();
  pathList.Add(new TObjString("$WORK/Macros/MuonEfficiency"));
  TList fileList; fileList.SetOwner();
  fileList.Add(new TObjString("runJPsiAccEffCorr.C"));
  fileList.Add(new TObjString("AliAnalysisTaskJPsiAccEffCorr.cxx"));
  fileList.Add(new TObjString("AliAnalysisTaskJPsiAccEffCorr.h"));
  CopyFileLocally(pathList, fileList);
  
  // --- prepare environment ---
  TString extraLibs="";
  TString extraIncs="";
//  TString extraLibs="PWG3base";
//  TString extraIncs="PWG3/base";
  TString extraTasks="AliAnalysisTaskJPsiAccEffCorr";
  LoadAlirootLocally(extraLibs, extraIncs, extraTasks);
  AliAnalysisGrid *alienHandler = 0x0;
  if (mode == kProof || mode == kProofLite) LoadAlirootOnProof(smode, rootVersion, alirootVersion, extraLibs, extraIncs, extraTasks, kTRUE);
  else if (mode == kGrid || mode == kTerminate) {
    TString analysisMacroName = "AccEff";
    alienHandler = static_cast<AliAnalysisGrid*>(CreateAlienHandler(smode, rootVersion, alirootVersion, inputFileName, dataDir, dataPattern, outDir, extraLibs, extraIncs, extraTasks, analysisMacroName, runFormat, ttl, maxFilesPerJob, maxMergeFiles, maxMergeStages));
    if (!alienHandler) return;
  }
  
  // --- Create the analysis train ---
  CreateAnalysisTrain(alienHandler);
  
  // --- Create input object ---
  TObject* inputObj = CreateInputObject(mode, inputFileName);
  
  // --- start analysis ---
  StartAnalysis(mode, inputObj);
  
}

//______________________________________________________________________________
void CreateAnalysisTrain(TObject* alienHandler)
{
  /// create the analysis train and configure it
  
  // analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("JPsiAccEffCorrAnalysis");
  
  // Connect plugin to the analysis manager if any
  if (alienHandler) mgr->SetGridHandler(static_cast<AliAnalysisGrid*>(alienHandler));
  
  // AOD handler
  AliInputEventHandler* aodH = new AliAODInputHandler;
  mgr->SetInputEventHandler(aodH);
  
  // Acc*Eff results
  AliAnalysisTaskJPsiAccEffCorr *jPsiAccEffCorr = new AliAnalysisTaskJPsiAccEffCorr("JPsiAccEffCorr");
//  jPsiAccEffCorr->SelectCollisionCandidates(AliVEvent::kMB | AliVEvent::kMUON);
//  jPsiAccEffCorr->SelectCollisionCandidates(AliVEvent::kMB);
//  jPsiAccEffCorr->SelectCollisionCandidates(AliVEvent::kAny);
  jPsiAccEffCorr->SetTrigLevel(1);
//  jPsiAccEffCorr->SetMuLowPtCut(1.15);
//  Float_t pTBinLowEdge[] = {0., 1., 2., 3., 4., 5., 6., 8.};
//  Float_t pTBinLowEdge[] = {0., 2., 5., 8.};
//  jPsiAccEffCorr->SetPtBins((Int_t)(sizeof(pTBinLowEdge)/sizeof(Float_t))-1, pTBinLowEdge);
//  Float_t yBinLowEdge[] = {-4., -3.75, -3.5, -3.25, -3., -2.75, -2.5};
//  Float_t yBinLowEdge[] = {-4., -3.25, -2.5};
//  jPsiAccEffCorr->SetYBins((Int_t)(sizeof(yBinLowEdge)/sizeof(Float_t))-1, yBinLowEdge);
  
  mgr->AddTask(jPsiAccEffCorr);
  mgr->ConnectInput(jPsiAccEffCorr, 0, mgr->GetCommonInputContainer());
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  outputfile += ":JPsiAccEffCorr";
  AliAnalysisDataContainer *histo = mgr->CreateContainer("Histograms", TObjArray::Class(), AliAnalysisManager::kOutputContainer, outputfile);
  mgr->ConnectOutput(jPsiAccEffCorr, 1, histo);
  AliAnalysisDataContainer *eventStat = mgr->CreateContainer("EventCounters", AliCounterCollection::Class(), AliAnalysisManager::kOutputContainer, outputfile);
  mgr->ConnectOutput(jPsiAccEffCorr, 2, eventStat);
  AliAnalysisDataContainer *jPsiStat = mgr->CreateContainer("JPsiCounters", AliCounterCollection::Class(), AliAnalysisManager::kOutputContainer, outputfile);
  mgr->ConnectOutput(jPsiAccEffCorr, 3, jPsiStat);
}

