//
//  runJPsi.C
//  aliroot_dev
//
//  Created by philippe pillot on 13/05/2014.
//  Copyright (c) 2014 Philippe Pillot. All rights reserved.
//

//TString rootVersion = "v5-34-05";
//TString alirootVersion = "v5-04-65-AN";
TString rootVersion = "v5-34-08";
TString alirootVersion = "vAN-20140630";
TString dataDir = "/alice/data/2011/LHC11h";
TString dataPattern = "ESDs/pass2_muon/*AliESDs.root";
TString runFormat = "%09d";
TString outDir = "Data/LHC11h/pass2_muon/JPsipDCA";
Int_t ttl = 30000;
Int_t maxFilesPerJob = 100;
Int_t maxMergeFiles = 10;
Int_t maxMergeStages = 2;

//______________________________________________________________________________
void runJPsi(TString smode = "local", TString inputFileName = "AliESDs.root",
             Bool_t applyPhysSel = kTRUE, Bool_t mc = kFALSE, Bool_t embedding = kFALSE)
{
  /// Extract some physical quantities
  
  gROOT->LoadMacro("$WORK/Macros/Facilities/runTaskFacilities.C");
  
  // --- Check runing mode ---
  Int_t mode = GetMode(smode, inputFileName);
  if(mode < 0) {
    Error("runJPsi","Please provide either an ESD root file, a collection of ESDs or a dataset.");
    return;
  }
  
  // get data type
  TString dataType = "unknown";
  if (mode == kGrid || mode == kTerminate) {
    if (dataPattern.Contains("AOD")) dataType = "AOD";
    else dataType = "ESD";
  } else {
    dataType = GetDataType(inputFileName);
    if (dataType == "unknown") {
      Error("runJPsi","Unable to determine the data type from the input file name or its content in case of a txt file.");
      return;
    }
  }
  
  // --- copy files needed for this analysis ---
  TList pathList; pathList.SetOwner();
  pathList.Add(new TObjString("$WORK/Macros/MuonPhysics/JPsiAnalysis/TaskWoMixing"));
  pathList.Add(new TObjString("$DEV/aliroot/PWG/muondep"));
  TList fileList; fileList.SetOwner();
  fileList.Add(new TObjString("runJPsi.C"));
  fileList.Add(new TObjString("AddTaskJPsi.C"));
  fileList.Add(new TObjString("AliAnalysisTaskJPsi.cxx"));
  fileList.Add(new TObjString("AliAnalysisTaskJPsi.h"));
  if (dataType == "ESD") {
    fileList.Add(new TObjString("AddTaskMTRSign.C"));
    fileList.Add(new TObjString("AliAnalysisTaskMTRSign.cxx"));
    fileList.Add(new TObjString("AliAnalysisTaskMTRSign.h"));
    fileList.Add(new TObjString("AddTaskMuonCDBConnect.C"));
    fileList.Add(new TObjString("AliAnalysisTaskMuonCDBConnect.cxx"));
    fileList.Add(new TObjString("AliAnalysisTaskMuonCDBConnect.h"));
    fileList.Add(new TObjString("AddTaskMuonRefit.C"));
    fileList.Add(new TObjString("AliAnalysisTaskMuonRefit.cxx"));
    fileList.Add(new TObjString("AliAnalysisTaskMuonRefit.h"));
    if (mc) {
      fileList.Add(new TObjString("AddTaskESDMCLabelAddition.C"));
      fileList.Add(new TObjString("AliAnalysisTaskESDMCLabelAddition.cxx"));
      fileList.Add(new TObjString("AliAnalysisTaskESDMCLabelAddition.h"));
    }
  }
  CopyFileLocally(pathList, fileList);
  
  // --- prepare environment ---
  TString extraLibs, extraIncs, extraTasks;
  if (dataType == "AOD") {
    extraLibs = "CORRFW:PWGmuon";
    extraIncs = "include:PWG/muon";
    extraTasks = "AliAnalysisTaskJPsi";
  } else {
    extraLibs = "RAWDatabase:RAWDatarec:CDB:STEER:MUONcore:MUONmapping:MUONcalib:MUONgeometry:MUONtrigger:MUONraw:MUONbase:MUONshuttle:MUONrec:RAWDatasim:MUONsim:MUONevaluation:CORRFW:PWGmuon";
    extraIncs = "include:MUON:MUON/mapping:PWG/muon";
    if (mc) extraTasks = "AliAnalysisTaskMuonCDBConnect:AliAnalysisTaskMuonRefit:AliAnalysisTaskESDMCLabelAddition:AliAnalysisTaskJPsi:AliAnalysisTaskMTRSign";
    else extraTasks = "AliAnalysisTaskMuonCDBConnect:AliAnalysisTaskMuonRefit:AliAnalysisTaskJPsi:AliAnalysisTaskMTRSign";
  }
  LoadAlirootLocally(extraLibs, extraIncs, extraTasks);
  AliAnalysisGrid *alienHandler = 0x0;
  if (mode == kProof || mode == kProofLite) LoadAlirootOnProof(smode, rootVersion, alirootVersion, extraLibs, extraIncs, extraTasks, kTRUE);
  else if (mode == kGrid || mode == kTerminate) {
    TString analysisMacroName = "JPsi";
    alienHandler = static_cast<AliAnalysisGrid*>(CreateAlienHandler(smode, rootVersion, alirootVersion, inputFileName, dataDir, dataPattern, outDir, extraLibs, extraIncs, extraTasks, analysisMacroName, runFormat, ttl, maxFilesPerJob, maxMergeFiles, maxMergeStages));
    if (!alienHandler) return;
  }
  
  // --- Create the analysis train ---
  CreateAnalysisTrain(dataType, applyPhysSel, mc, embedding, alienHandler);
  
  // --- Create input object ---
  TObject* inputObj = CreateInputObject(mode, inputFileName);
  
  // --- start analysis ---
  StartAnalysis(mode, inputObj);
  
}

//______________________________________________________________________________
void CreateAnalysisTrain(TString dataType, Bool_t applyPhysSel, Bool_t mc, Bool_t embedding, TObject *alienHandler)
{
  /// create the analysis train and configure it
  
  // analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("JPsiAnalysis");
  
  // Connect plugin to the analysis manager if any
  if (alienHandler) mgr->SetGridHandler(static_cast<AliAnalysisGrid*>(alienHandler));
  //mgr->SetDebugLevel(3);
  
  if (dataType == "AOD") {
    
    // AOD input
    AliInputEventHandler* aodH = new AliAODInputHandler;
    mgr->SetInputEventHandler(aodH);
    
  } else {
    
    // ESD input
    AliESDInputHandler* esdH = new AliESDInputHandler();
    esdH->SetReadFriends(kFALSE);
    esdH->SetInactiveBranches("*");
    esdH->SetActiveBranches("MuonTracks MuonClusters MuonPads AliESDRun. AliESDHeader. AliMultiplicity. AliESDFMD. AliESDVZERO. SPDVertex. PrimaryVertex. AliESDZDC. AliESDTZERO.");
    mgr->SetInputEventHandler(esdH);
    
    // Monte Carlo input
    if (mc) {
      AliMCEventHandler* mcH = new AliMCEventHandler();
      mgr->SetMCtruthEventHandler(mcH);
    }
    
  }
  
  // CDB connection
  if (dataType == "ESD") {
    gROOT->LoadMacro("AddTaskMuonCDBConnect.C");
    AliAnalysisTaskMuonCDBConnect *taskCDBConnect = AddTaskMuonCDBConnect();
    if(!taskCDBConnect) {
      Error("CreateAnalysisTrain","AliAnalysisTaskMuonCDBConnect not created!");
      return;
    }
    if (mc && !embedding) taskCDBConnect->SetAlignStorage("alien://folder=/alice/simulation/2008/v4-15-Release/Residual");
    taskCDBConnect->LoadGeometry();
    taskCDBConnect->LoadMagField();
  }
  
  // event selection
  if (applyPhysSel && dataType == "ESD") {
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
    AliPhysicsSelectionTask* physicsSelection = AddTaskPhysicsSelection(mc && !embedding);
    if(!physicsSelection) {
      Error("CreateAnalysisTrain","AliPhysicsSelectionTask not created!");
      return;
    }
  }
  UInt_t offlineTriggerMask = AliVEvent::kAny;
  if (mc && embedding) offlineTriggerMask = AliVEvent::kMB;
  else if (!mc) offlineTriggerMask = AliVEvent::kMuonUnlikePB | AliVEvent::kMuonLikePB;
  
  // centrality selection
  if (dataType == "ESD") {
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskCentrality.C");
    AliCentralitySelectionTask *taskCentrality = AddTaskCentrality();
    if(!taskCentrality) {
      Error("CreateAnalysisTrain","AliCentralitySelectionTask not created!");
      return;
    }
    AliInputEventHandler* hdl = static_cast<AliInputEventHandler*>(mgr->GetInputEventHandler());
    hdl->SetNeedField(kFALSE);
    if (applyPhysSel) taskCentrality->SelectCollisionCandidates(offlineTriggerMask);
    if (mc && !embedding) taskCentrality->SetMCInput();
  }
  
  // track selection
  AliMuonTrackCuts trackCuts("stdCuts", "stdCuts");
  trackCuts.SetAllowDefaultParams();
  trackCuts.SetFilterMask(AliMuonTrackCuts::kMuMatchLpt | AliMuonTrackCuts::kMuEta |
			  AliMuonTrackCuts::kMuThetaAbs | AliMuonTrackCuts::kMuPdca);
  trackCuts.SetIsMC(mc);
  //  trackCuts.SetPassNumber(2);
  
  // track refitting
  if (dataType == "ESD") {
    gROOT->LoadMacro("AddTaskMuonRefit.C");
    AliAnalysisTaskMuonRefit *refit = AddTaskMuonRefit(-1., -1., kTRUE, -1., -1.);
    if(!refit) {
      Error("CreateAnalysisTrain","AliAnalysisTaskMuonRefit not created!");
      return;
    }
    if (applyPhysSel) refit->SelectCollisionCandidates(offlineTriggerMask);
    if (mc && !embedding) refit->SetAlignStorage("alien://folder=/alice/simulation/2008/v4-15-Release/Residual");
    refit->RemoveMonoCathodClusters(kTRUE, kFALSE);
  }
  
  // setup the train to run on MC or on data
  if (mc) {
    
    // JPsi task 1
    gROOT->LoadMacro("AddTaskJPsi.C");
    AliAnalysisTaskJPsi* jpsi1 = AddTaskJPsi("MB");
    if(!jpsi1) {
      Error("CreateAnalysisTrain","AliAnalysisTaskJPsi not created!");
      return;
    }
    if (applyPhysSel) jpsi1->SelectCollisionCandidates(offlineTriggerMask);
    jpsi1->SetMuonTrackCuts(trackCuts);
    jpsi1->UseMCLabel();
    SetKinematicsRange(jpsi1);
    
    // JPsi task 2
    AliAnalysisTaskJPsi* jpsi2 = AddTaskJPsi("MB_TrgSign");
    if (applyPhysSel) jpsi2->SelectCollisionCandidates(offlineTriggerMask);
    jpsi2->SetMuonTrackCuts(trackCuts);
    jpsi2->SelectTrgSign();
    jpsi2->UseMCLabel();
    SetKinematicsRange(jpsi2);
    
    // JPsi task 15
    AliAnalysisTaskJPsi* jpsi15 = AddTaskJPsi("MB_TrgFake");
    if (applyPhysSel) jpsi15->SelectCollisionCandidates(offlineTriggerMask);
    jpsi15->SetMuonTrackCuts(trackCuts);
    jpsi15->SelectSameTrgSignFake();
    jpsi15->UseMCLabel();
    SetKinematicsRange(jpsi15);
    
    if (dataType == "ESD") {
      
      // MTR sign task
      gROOT->LoadMacro("AddTaskMTRSign.C");
      AliAnalysisTaskMTRSign* mtrSign = AddTaskMTRSign(mc);
      if(!mtrSign) {
        Error("CreateAnalysisTrain","AliAnalysisTaskMTRSign not created!");
        return;
      }
      if (applyPhysSel) mtrSign->SelectCollisionCandidates(offlineTriggerMask);
      mtrSign->GetMuonEventCuts()->SetFilterMask(AliMuonEventCuts::kSelectedCentrality | AliMuonEventCuts::kSelectedTrig);
      mtrSign->GetMuonEventCuts()->SetTrigClassPatterns("CPBI2_B1-B-NOPF-ALLNOTRD");
      Double_t centralityBins[] = {0., 10., 20., 30., 40., 50., 60., 70., 80., 90.};
      Int_t nCentralityBins = sizeof(centralityBins)/sizeof(centralityBins[0])-1;
      mtrSign->GetMuonEventCuts()->SetCentralityClasses(nCentralityBins, centralityBins);
      *(mtrSign->GetMuonTrackCuts()) = trackCuts;
      mtrSign->UseMCLabel();
      
      // recompute MC labels
      gROOT->LoadMacro("AddTaskESDMCLabelAddition.C");
      AliAnalysisTaskESDMCLabelAddition *esdmclabel = AddTaskESDMCLabelAddition();
      if(!esdmclabel) {
        Error("CreateAnalysisTrain","AliAnalysisTaskESDMCLabelAddition not created!");
        return;
      }
      if (applyPhysSel) esdmclabel->SelectCollisionCandidates(offlineTriggerMask);
      
      // JPsi task 3
      AliAnalysisTaskJPsi* jpsi3 = AddTaskJPsi("MB_relabeled");
      if (applyPhysSel) jpsi3->SelectCollisionCandidates(offlineTriggerMask);
      jpsi3->SetMuonTrackCuts(trackCuts);
      jpsi3->UseMCLabel();
      SetKinematicsRange(jpsi3);
      
      // JPsi task 4
      AliAnalysisTaskJPsi* jpsi4 = AddTaskJPsi("MB_TrgSign_relabeled");
      if (applyPhysSel) jpsi4->SelectCollisionCandidates(offlineTriggerMask);
      jpsi4->SetMuonTrackCuts(trackCuts);
      jpsi4->SelectTrgSign();
      jpsi4->UseMCLabel();
      SetKinematicsRange(jpsi4);
      
      // recompute MC labels discarding the decays
      AliAnalysisTaskESDMCLabelAddition *esdmclabel2 = AddTaskESDMCLabelAddition();
      if (applyPhysSel) esdmclabel2->SelectCollisionCandidates(offlineTriggerMask);
      esdmclabel2->DecayAsFake();
      
      // JPsi task 5
      AliAnalysisTaskJPsi* jpsi5 = AddTaskJPsi("MB_relabeled_woDecays");
      if (applyPhysSel) jpsi5->SelectCollisionCandidates(offlineTriggerMask);
      jpsi5->SetMuonTrackCuts(trackCuts);
      jpsi5->UseMCLabel();
      SetKinematicsRange(jpsi5);
      
      // JPsi task 6
      AliAnalysisTaskJPsi* jpsi6 = AddTaskJPsi("MB_TrgSign_relabeled_woDecays");
      if (applyPhysSel) jpsi6->SelectCollisionCandidates(offlineTriggerMask);
      jpsi6->SetMuonTrackCuts(trackCuts);
      jpsi6->SelectTrgSign();
      jpsi6->UseMCLabel();
      SetKinematicsRange(jpsi6);
      
      // recompute MC labels with 10-sigma cut
      AliAnalysisTaskESDMCLabelAddition *esdmclabel3 = AddTaskESDMCLabelAddition();
      if (applyPhysSel) esdmclabel3->SelectCollisionCandidates(offlineTriggerMask);
      esdmclabel3->SetExternalTrkSigmaCut(10.);
      
      // JPsi task 7
      AliAnalysisTaskJPsi* jpsi7 = AddTaskJPsi("MB_relabeled_10sigma");
      if (applyPhysSel) jpsi7->SelectCollisionCandidates(offlineTriggerMask);
      jpsi7->SetMuonTrackCuts(trackCuts);
      jpsi7->UseMCLabel();
      SetKinematicsRange(jpsi7);
      
      // JPsi task 8
      AliAnalysisTaskJPsi* jpsi8 = AddTaskJPsi("MB_TrgSign_relabeled_10sigma");
      if (applyPhysSel) jpsi8->SelectCollisionCandidates(offlineTriggerMask);
      jpsi8->SetMuonTrackCuts(trackCuts);
      jpsi8->SelectTrgSign();
      jpsi8->UseMCLabel();
      SetKinematicsRange(jpsi8);
      
      // recompute MC labels with 10-sigma cut and discarding the decays
      AliAnalysisTaskESDMCLabelAddition *esdmclabel4 = AddTaskESDMCLabelAddition();
      if (applyPhysSel) esdmclabel4->SelectCollisionCandidates(offlineTriggerMask);
      esdmclabel4->SetExternalTrkSigmaCut(10.);
      esdmclabel4->DecayAsFake();
      
      // JPsi task 9
      AliAnalysisTaskJPsi* jpsi9 = AddTaskJPsi("MB_relabeled_10sigma_woDecays");
      if (applyPhysSel) jpsi9->SelectCollisionCandidates(offlineTriggerMask);
      jpsi9->SetMuonTrackCuts(trackCuts);
      jpsi9->UseMCLabel();
      SetKinematicsRange(jpsi9);
      
      // JPsi task 10
      AliAnalysisTaskJPsi* jpsi10 = AddTaskJPsi("MB_TrgSign_relabeled_10sigma_woDecays");
      if (applyPhysSel) jpsi10->SelectCollisionCandidates(offlineTriggerMask);
      jpsi10->SetMuonTrackCuts(trackCuts);
      jpsi10->SelectTrgSign();
      jpsi10->UseMCLabel();
      SetKinematicsRange(jpsi10);
      
      // recompute MC labels with 25-sigma cut
      AliAnalysisTaskESDMCLabelAddition *esdmclabel5 = AddTaskESDMCLabelAddition();
      if (applyPhysSel) esdmclabel5->SelectCollisionCandidates(offlineTriggerMask);
      esdmclabel5->SetExternalTrkSigmaCut(25.);
      
      // JPsi task 11
      AliAnalysisTaskJPsi* jpsi11 = AddTaskJPsi("MB_relabeled_25sigma");
      if (applyPhysSel) jpsi11->SelectCollisionCandidates(offlineTriggerMask);
      jpsi11->SetMuonTrackCuts(trackCuts);
      jpsi11->UseMCLabel();
      SetKinematicsRange(jpsi11);
      
      // JPsi task 12
      AliAnalysisTaskJPsi* jpsi12 = AddTaskJPsi("MB_TrgSign_relabeled_25sigma");
      if (applyPhysSel) jpsi12->SelectCollisionCandidates(offlineTriggerMask);
      jpsi12->SetMuonTrackCuts(trackCuts);
      jpsi12->SelectTrgSign();
      jpsi12->UseMCLabel();
      SetKinematicsRange(jpsi12);
      
      // recompute MC labels with 25-sigma cut and discarding the decays
      AliAnalysisTaskESDMCLabelAddition *esdmclabel6 = AddTaskESDMCLabelAddition();
      if (applyPhysSel) esdmclabel6->SelectCollisionCandidates(offlineTriggerMask);
      esdmclabel6->SetExternalTrkSigmaCut(25.);
      esdmclabel6->DecayAsFake();
      
      // JPsi task 13
      AliAnalysisTaskJPsi* jpsi13 = AddTaskJPsi("MB_relabeled_25sigma_woDecays");
      if (applyPhysSel) jpsi13->SelectCollisionCandidates(offlineTriggerMask);
      jpsi13->SetMuonTrackCuts(trackCuts);
      jpsi13->UseMCLabel();
      SetKinematicsRange(jpsi13);
      
      // JPsi task 14
      AliAnalysisTaskJPsi* jpsi14 = AddTaskJPsi("MB_TrgSign_relabeled_25sigma_woDecays");
      if (applyPhysSel) jpsi14->SelectCollisionCandidates(offlineTriggerMask);
      jpsi14->SetMuonTrackCuts(trackCuts);
      jpsi14->SelectTrgSign();
      jpsi14->UseMCLabel();
      SetKinematicsRange(jpsi14);
      
    }
    
  } else {
    
    // JPsi task 1
    gROOT->LoadMacro("AddTaskJPsi.C");
    AliAnalysisTaskJPsi* jpsi1 = AddTaskJPsi("MUL");
    if(!jpsi1) {
      Error("CreateAnalysisTrain","AliAnalysisTaskJPsi not created!");
      return;
    }
    if (applyPhysSel) jpsi1->SelectCollisionCandidates(AliVEvent::kMuonUnlikePB);
    jpsi1->SetMuonTrackCuts(trackCuts);
    SetKinematicsRange(jpsi1);
    
    // JPsi task 2
    AliAnalysisTaskJPsi* jpsi2 = AddTaskJPsi("MULorMLL");
    if (applyPhysSel) jpsi2->SelectCollisionCandidates(AliVEvent::kMuonUnlikePB | AliVEvent::kMuonLikePB);
    jpsi2->SetMuonTrackCuts(trackCuts);
    SetKinematicsRange(jpsi2);
    jpsi2->RecordEvWithTrgIssues();
    
    // JPsi task 3
    AliAnalysisTaskJPsi* jpsi3 = AddTaskJPsi("MUL_TrgSign");
    if (applyPhysSel) jpsi3->SelectCollisionCandidates(AliVEvent::kMuonUnlikePB);
    jpsi3->SetMuonTrackCuts(trackCuts);
    jpsi3->SelectTrgSign();
    SetKinematicsRange(jpsi3);
    
    // JPsi task 4
    AliAnalysisTaskJPsi* jpsi4 = AddTaskJPsi("MULorMLL_TrgSign");
    if (applyPhysSel) jpsi4->SelectCollisionCandidates(AliVEvent::kMuonUnlikePB | AliVEvent::kMuonLikePB);
    jpsi4->SetMuonTrackCuts(trackCuts);
    jpsi4->SelectTrgSign();
    SetKinematicsRange(jpsi4);
    
    // JPsi task 5
    AliAnalysisTaskJPsi* jpsi5 = AddTaskJPsi("MULorMLL_TrgFake");
    if (applyPhysSel) jpsi5->SelectCollisionCandidates(AliVEvent::kMuonUnlikePB | AliVEvent::kMuonLikePB);
    jpsi5->SetMuonTrackCuts(trackCuts);
    jpsi5->SelectSameTrgSignFake();
    SetKinematicsRange(jpsi5);
    
    if (dataType == "ESD") {
      
      // MTR sign task
      gROOT->LoadMacro("AddTaskMTRSign.C");
      AliAnalysisTaskMTRSign* mtrSign = AddTaskMTRSign(mc);
      if(!mtrSign) {
        Error("CreateAnalysisTrain","AliAnalysisTaskMTRSign not created!");
        return;
      }
      if (applyPhysSel) mtrSign->SelectCollisionCandidates(AliVEvent::kMuonUnlikePB | AliVEvent::kMuonLikePB);
      mtrSign->GetMuonEventCuts()->SetFilterMask(AliMuonEventCuts::kSelectedCentrality | AliMuonEventCuts::kSelectedTrig);
      mtrSign->GetMuonEventCuts()->SetTrigClassPatterns("CPBI1MUL-B-NOPF-MUON,CPBI1MLL-B-NOPF-MUON,CPBI1MUL-B-NOPF-MUON&CPBI1MLL-B-NOPF-MUON");
      Double_t centralityBins[] = {0., 10., 20., 30., 40., 50., 60., 70., 80., 90.};
      Int_t nCentralityBins = sizeof(centralityBins)/sizeof(centralityBins[0])-1;
      mtrSign->GetMuonEventCuts()->SetCentralityClasses(nCentralityBins, centralityBins);
      *(mtrSign->GetMuonTrackCuts()) = trackCuts;
      
    }
    
  }
  
}

//______________________________________________________________________________
void SetKinematicsRange(TObject *jpsiTask)
{
  /// set pT and y ranges in the JPsi analysis task
  
  // define pT bins
  const Int_t nPtBins = 19;
  Float_t dPtLowEdge[nPtBins] = {0., 0., 2., 5., 0., 1., 2., 3., 4., 5., 6., 8., 2., 0., 3., 0.3, 0.3, 0.3, 0.3};
  Float_t dPtUpEdge[nPtBins] = {8., 2., 5., 8., 1., 2., 3., 4., 5., 6., 8., 20., 8., 3., 8., 1., 8., 2., 3.};
  
  // define y bins
  const Int_t nYBins = 10;
  Float_t dYLowEdge[nYBins] = { -2.5, -2.5, -3., -3.5, -2.5, -2.75, -3., -3.25, -3.5, -3.75};
  Float_t dYUpEdge[nYBins] = { -4., -3., -3.5, -4., -2.75, -3., -3.25, -3.5, -3.75, -4.};
  
  // set ranges
  AliAnalysisTaskJPsi *jpsi = static_cast<AliAnalysisTaskJPsi*>(jpsiTask);
  jpsi->SetPtLowEdge(nPtBins, dPtLowEdge);
  jpsi->SetPtUpEdge(nPtBins, dPtUpEdge);
  jpsi->SetYLowEdge(nYBins, dYLowEdge);
  jpsi->SetYUpEdge(nYBins, dYUpEdge);
  
}

