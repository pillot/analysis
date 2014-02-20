/*
 *  runMuonPhysics.C
 *  aliroot
 *
 *  Created by Philippe Pillot on 11/12/12.
 *  Copyright 2010 SUBATECH. All rights reserved.
 *
 */

TString rootVersion = "v5-34-05";
TString alirootVersion = "v5-04-65-AN";
//TString dataDir = "/alice/cern.ch/user/p/ppillot/Data/LHC11h/pass2_muon/AODs/results";
TString dataDir = "/alice/data/2012/LHC12h";
//TString dataPattern = "*AliAOD.Muons.root";
TString dataPattern = "ESDs/muon_calo_pass2/AOD/*AliAOD.Muons.root";
TString runFormat = "%09d";
//TString outDir = "Data/LHC11h/pass2_muon/AODs/PhysTest";
TString outDir = "Data/LHC12h/muon_calo_pass2/AODs/Phys/any_muon";
Int_t ttl = 30000;
Int_t maxFilesPerJob = 20;
Int_t maxMergeFiles = 10;
Int_t maxMergeStages = 2;

//______________________________________________________________________________
void runMuonPhysicsAOD(TString smode = "local", TString inputFileName = "AliAOD.Muons.root")
{
  /// Extract some physical quantities
  
  gROOT->LoadMacro("$WORK/Macros/Facilities/runTaskFacilities.C");
  
  // --- Check runing mode ---
  Int_t mode = GetMode(smode, inputFileName);
  if(mode < 0){
    Error("runAODMuonPhysics","Please provide either an ESD or AOD root file a collection of ESDs or AODs or a dataset.");
    return;
  }
  
  // --- copy files needed for this analysis ---
  TList pathList; pathList.SetOwner();
  pathList.Add(new TObjString("$WORK/Macros/MuonPhysics"));
  TList fileList; fileList.SetOwner();
  fileList.Add(new TObjString("runMuonPhysicsAOD.C"));
  fileList.Add(new TObjString("AddTaskMuonPhysics.C"));
  fileList.Add(new TObjString("AliAnalysisTaskMuonPhysics.cxx"));
  fileList.Add(new TObjString("AliAnalysisTaskMuonPhysics.h"));
  CopyFileLocally(pathList, fileList);
  
  // --- prepare environment ---
  TString extraLibs="CORRFW:PWGmuon";
  TString extraIncs="include:PWG/muon";
  TString extraTasks="AliAnalysisTaskMuonPhysics";
  LoadAlirootLocally(extraLibs, extraIncs, extraTasks);
  AliAnalysisGrid *alienHandler = 0x0;
  if (mode == kProof || mode == kProofLite) LoadAlirootOnProof(smode, rootVersion, alirootVersion, extraLibs, extraIncs, extraTasks, kTRUE);
  else if (mode == kGrid || mode == kTerminate) {
    TString analysisMacroName = "Physics";
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
  AliAnalysisManager *mgr = new AliAnalysisManager("MuonPhysicsAnalysis");
  
  // Connect plugin to the analysis manager if any
  if (alienHandler) mgr->SetGridHandler(static_cast<AliAnalysisGrid*>(alienHandler));
  //mgr->SetDebugLevel(3);
  
  // AOD handler
  AliInputEventHandler* aodH = new AliAODInputHandler;
  mgr->SetInputEventHandler(aodH);
  
  // track selection
  AliMuonTrackCuts trackCuts("stdCuts", "stdCuts");
  trackCuts.SetAllowDefaultParams();
  trackCuts.SetFilterMask(0);
//  trackCuts.SetFilterMask(AliMuonTrackCuts::kMuMatchLpt | AliMuonTrackCuts::kMuEta |
//			  AliMuonTrackCuts::kMuThetaAbs | AliMuonTrackCuts::kMuPdca);
//  trackCuts.SetIsMC(kTRUE);
//  trackCuts.SetCustomParamFromRun(169859, "pass2_muon");
//  trackCuts.CustomParam()->SetMeanDCA(0., 0., 0.);
//  trackCuts.CustomParam()->SetSigmaPdca(80., 54.);
//  trackCuts.CustomParam()->SetRelPResolution(0.0004);
//  trackCuts.CustomParam()->SetSlopeResolution(0.0005);
//  trackCuts.CustomParam()->SetNSigmaPdca(5.);
//  trackCuts.Print();
/*  
  // track rejection
  AliMuonTrackCuts trackCuts2("stdCuts2", "stdCuts2");
  trackCuts2.SetFilterMask(AliMuonTrackCuts::kMuPdca);
  trackCuts2.SetCustomParamFromRun(169859, "pass2_muon");
//  trackCuts2.CustomParam()->SetMeanDCA(0., 0., 0.);
  trackCuts2.CustomParam()->SetSigmaPdca(80, 54.);
//  trackCuts2.CustomParam()->SetNSigmaPdca(5.);
//  trackCuts2.CustomParam()->SetRelPResolution(0.0005);
//  trackCuts2.CustomParam()->SetSlopeResolution(0.0004);
//  trackCuts2.CustomParam()->SetChi2NormCut(3.5);
*/  
  // physics task
  gROOT->LoadMacro("AddTaskMuonPhysics.C");
  AliAnalysisTaskMuonPhysics* physics = AddTaskMuonPhysics("");
  if(!physics) {
    Error("CreateAnalysisTrain","AliAnalysisTaskMuonPhysics not created!");
    return;
  }
  physics->SelectCollisionCandidates(AliVEvent::kAnyINT);
//  physics->SelectCollisionCandidates(AliVEvent::kMUU7 | AliVEvent::kMuonUnlikeLowPt8);
//  physics->SelectCollisionCandidates(AliVEvent::kMUSH7 | AliVEvent::kMuonSingleHighPt8);
//  physics->SelectCollisionCandidates(AliVEvent::kMUU7);
//  physics->SelectCollisionCandidates(AliVEvent::kMUL7 | AliVEvent::kMUU7 |
//				     AliVEvent::kMuonLikeLowPt8 | AliVEvent::kMuonUnlikeLowPt8);
//  physics->SelectCentrality(0., 90.);
  physics->SetMuonTrackCuts(trackCuts);
//  physics->SetMuonTrackCuts2(trackCuts2);
//  physics->VersusRun(kTRUE);
  
  // track selection 2
  AliMuonTrackCuts trackCuts2("stdCuts2", "stdCuts2");
  trackCuts2.SetAllowDefaultParams();
  trackCuts2.SetFilterMask(AliMuonTrackCuts::kMuMatchApt | AliMuonTrackCuts::kMuEta |
  			  AliMuonTrackCuts::kMuThetaAbs);
  
  // physics task 2
  AliAnalysisTaskMuonPhysics* physics2 = AddTaskMuonPhysics("stdCuts");
  if(!physics2) {
    Error("CreateAnalysisTrain","AliAnalysisTaskMuonPhysics not created!");
    return;
  }
  physics2->SelectCollisionCandidates(AliVEvent::kAnyINT);
  physics2->SetMuonTrackCuts(trackCuts2);
  
  // track selection 3
  AliMuonTrackCuts trackCuts3("stdCuts3", "stdCuts3");
  trackCuts3.SetAllowDefaultParams();
  trackCuts3.SetFilterMask(AliMuonTrackCuts::kMuMatchApt | AliMuonTrackCuts::kMuEta |
  			  AliMuonTrackCuts::kMuThetaAbs | AliMuonTrackCuts::kMuPdca);
  
  // physics task 3
  AliAnalysisTaskMuonPhysics* physics3 = AddTaskMuonPhysics("stdCuts_pDCA");
  if(!physics3) {
    Error("CreateAnalysisTrain","AliAnalysisTaskMuonPhysics not created!");
    return;
  }
  physics3->SelectCollisionCandidates(AliVEvent::kAnyINT);
  physics3->SetMuonTrackCuts(trackCuts3);
  
  // track selection 4
  AliMuonTrackCuts trackCuts4("stdCuts4", "stdCuts4");
  trackCuts4.SetAllowDefaultParams();
  trackCuts4.SetFilterMask(AliMuonTrackCuts::kMuEta | AliMuonTrackCuts::kMuThetaAbs);
  
  // physics task 4
  AliAnalysisTaskMuonPhysics* physics4 = AddTaskMuonPhysics("accCuts");
  if(!physics4) {
    Error("CreateAnalysisTrain","AliAnalysisTaskMuonPhysics not created!");
    return;
  }
  physics4->SelectCollisionCandidates(AliVEvent::kAnyINT);
  physics4->SetMuonTrackCuts(trackCuts4);
  
  // track selection 5
  AliMuonTrackCuts trackCuts5("stdCuts5", "stdCuts5");
  trackCuts5.SetAllowDefaultParams();
  trackCuts5.SetFilterMask(AliMuonTrackCuts::kMuEta | AliMuonTrackCuts::kMuThetaAbs |
                           AliMuonTrackCuts::kMuPdca);
  
  // physics task 5
  AliAnalysisTaskMuonPhysics* physics5 = AddTaskMuonPhysics("accCuts_pDCA");
  if(!physics5) {
    Error("CreateAnalysisTrain","AliAnalysisTaskMuonPhysics not created!");
    return;
  }
  physics5->SelectCollisionCandidates(AliVEvent::kAnyINT);
  physics5->SetMuonTrackCuts(trackCuts5);
  
}

