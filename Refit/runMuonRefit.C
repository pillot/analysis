/*
 *  runMuonRefit.C
 *  aliroot
 *
 *  Created by Philippe Pillot on 02/11/10.
 *  Copyright 2010 SUBATECH. All rights reserved.
 *
 */


//______________________________________________________________________________
void runMuonRefit(TString smode = "local", TString inputFileName = "AliESDs.root",
		  Double_t resNB = -1., Double_t resB = -1., Double_t sigmaCut = -1.,
		  Double_t sigmaCutTrig = -1., Bool_t selectPhysics = kTRUE, Bool_t selectCentrality = kFALSE)
{
  /// Refit the tracks with new recoParam and/or new alignment
  /// check the effect on track quality, efficiency and physical quantities
  
  // --- general analysis setup ---
  TString rootVersion = "";
  TString alirootVersion = "";
  TString aliphysicsVersion = "vAN-20160524-1";
  TString extraLibs="";
  TString extraIncs="include";
  TString extraTasks="AliAnalysisTaskMuonPhysics";
  TString extraPkgs="PWGmuondep";
  TList pathList; pathList.SetOwner();
  pathList.Add(new TObjString("$WORK/Macros/MuonPhysics"));
  pathList.Add(new TObjString("$WORK/Macros/Refit"));
  pathList.Add(new TObjString("$ALICE_PHYSICS/PARfiles"));
  TList fileList; fileList.SetOwner();
  fileList.Add(new TObjString("runMuonRefit.C"));
  fileList.Add(new TObjString("PWGmuondep.par"));
  fileList.Add(new TObjString("AddTaskMuonPhysics.C"));
  fileList.Add(new TObjString("AliAnalysisTaskMuonPhysics.cxx"));
  fileList.Add(new TObjString("AliAnalysisTaskMuonPhysics.h"));
  
  // --- grid specific setup ---
  TString dataDir = "/alice/data/2015/LHC15j";
  TString dataPattern = "muon_calo_pass2/*AliESDs.root";
  TString runFormat = "%09d";
  TString outDir = "Data/LHC15j/muon_calo_pass2/Refit_AlignV6/CMUU7CINT7";
  TString analysisMacroName = "Refit";
  Int_t ttl = 30000;
  Int_t maxFilesPerJob = 100;
  Int_t maxMergeFiles = 10;
  Int_t maxMergeStages = 2;
  
  // --- saf3 specific setup ---
  Bool_t splitDataset = kFALSE;
  
  gROOT->LoadMacro("$HOME/Work/Alice/Macros/Facilities/runTaskFacilities.C");
  
  // --- prepare the analysis environment ---
  Int_t mode = PrepareAnalysis(smode, inputFileName, extraLibs, extraIncs, extraTasks, extraPkgs, pathList, fileList);
  
  // --- run the analysis (saf3 is a special case as the analysis is launched on the server) ---
  if (mode == kSAF3Connect) {
    
    RunAnalysisOnSAF3(fileList, aliphysicsVersion, inputFileName, splitDataset);
    
  } else {
    
    if (!CreateAnalysisTrain(resNB, resB, sigmaCut, sigmaCutTrig, selectPhysics, selectCentrality)) return;
    
    if (smode == "saf3" && splitDataset) AliAnalysisManager::GetAnalysisManager()->SetSkipTerminate(kTRUE);
    
    RunAnalysis(smode, inputFileName, rootVersion, alirootVersion, aliphysicsVersion, extraLibs, extraIncs, extraTasks, extraPkgs, dataDir, dataPattern, outDir, analysisMacroName, runFormat, ttl, maxFilesPerJob, maxMergeFiles, maxMergeStages);
    
  }
  
}

//______________________________________________________________________________
Bool_t CreateAnalysisTrain(Double_t resNB, Double_t resB, Double_t sigmaCut, Double_t sigmaCutTrig,
			   Bool_t selectPhysics, Bool_t selectCentrality)
{
  /// create the analysis train and configure it
  
  // analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("MuonRefitAnalysis");
  //mgr->SetAutoBranchLoading(kFALSE);
  //mgr->SetCacheSize(0);
  //mgr->SetDebugLevel(3);
  
  // ESD input
  AliESDInputHandler* esdH = new AliESDInputHandler();
  esdH->SetReadFriends(kFALSE);
  esdH->SetInactiveBranches("*");
  esdH->SetActiveBranches("MuonTracks MuonClusters MuonPads AliESDRun. AliESDHeader. AliMultiplicity. AliESDFMD. AliESDVZERO. SPDVertex. PrimaryVertex. AliESDZDC. AliESDTZERO.");
  mgr->SetInputEventHandler(esdH);
  
  // CDB connection
  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/muondep/AddTaskMuonCDBConnect.C");
  AliAnalysisTaskMuonCDBConnect *taskCDBConnect = AddTaskMuonCDBConnect();
  if(!taskCDBConnect) {
    Error("CreateAnalysisTrain","AliAnalysisTaskMuonCDBConnect not created!");
    return;
  }
//  taskCDBConnect->SetDefaultStorage("local:///cvmfs/alice-ocdb.cern.ch/calibration/data/2015/OCDB");
//  taskCDBConnect->SetRecoParamStorage("alien://folder=/alice/cern.ch/user/p/ppillot/OCDB2015_PbPb");
//  taskCDBConnect->SetAlignStorage("", 6);
  taskCDBConnect->LoadGeometry();
  taskCDBConnect->LoadMagField();
  
  // event selection
  if (selectPhysics) {
    gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
    AliPhysicsSelectionTask* physicsSelection = AddTaskPhysicsSelection();
    if(!physicsSelection) {
      Error("CreateAnalysisTrain","AliPhysicsSelectionTask not created!");
      return kFALSE;
    }
  }
//  UInt_t offlineTriggerMask = AliVEvent::kMUU7 | AliVEvent::kMUSPB;
  UInt_t offlineTriggerMask = AliVEvent::kMUU7;
//  UInt_t offlineTriggerMask = AliVEvent::kAny;
  
  // multiplicity/centrality selection
  if (selectCentrality) {
    gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
    AliMultSelectionTask *taskCentrality = AddTaskMultSelection(kFALSE);
    if(!taskCentrality) {
      Error("CreateAnalysisTrain","AliMultSelectionTask not created!");
      return;
    }
    //taskCentrality->SetAlternateOADBforEstimators("LHC15o");
    if (selectPhysics) taskCentrality->SelectCollisionCandidates(offlineTriggerMask);
  }
  
  // track selection
  AliMuonTrackCuts trackCuts("stdCuts", "stdCuts");
  trackCuts.SetAllowDefaultParams();
  trackCuts.SetFilterMask(AliMuonTrackCuts::kMuMatchLpt | AliMuonTrackCuts::kMuEta |
			  AliMuonTrackCuts::kMuThetaAbs | AliMuonTrackCuts::kMuPdca);
//  trackCuts.SetFilterMask(0);
  //trackCuts.SetIsMC();
  
  // first task (physics before refit)
  gROOT->LoadMacro("AddTaskMuonPhysics.C");
  AliAnalysisTaskMuonPhysics* physics1 = AddTaskMuonPhysics("before");
  if(!physics1) {
    Error("CreateAnalysisTrain","AliAnalysisTaskMuonPhysics not created!");
    return kFALSE;
  }
  if (selectPhysics) physics1->SelectCollisionCandidates(offlineTriggerMask);
//  if (selectPhysics) physics1->SelectCollisionCandidates(AliVEvent::kMUU7);
  if (selectCentrality) physics1->SelectCentrality(0., 90.);
  physics1->SetMuonTrackCuts(trackCuts);
  /*
  // first task (physics before refit)
  AliAnalysisTaskMuonPhysics* physics12 = AddTaskMuonPhysics("before_LowPt7");
  if(!physics12) {
    Error("CreateAnalysisTrain","AliAnalysisTaskMuonPhysics not created!");
    return kFALSE;
  }
//  if (selectPhysics) physics12->SelectCollisionCandidates(offlineTriggerMask);
  if (selectPhysics) physics12->SelectCollisionCandidates(AliVEvent::kMuonLikeLowPt7 | AliVEvent::kMuonUnlikeLowPt7 | AliVEvent::kMuonSingleLowPt7);
  if (selectCentrality) physics12->SelectCentrality(0., 90.);
  physics12->SetMuonTrackCuts(trackCuts);
  */
  // second task (refit)
  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/muondep/AddTaskMuonRefit.C");
  AliAnalysisTaskMuonRefit* refit = AddTaskMuonRefit(resNB, resB, kTRUE, sigmaCut, sigmaCutTrig);
  if(!refit) {
    Error("CreateAnalysisTrain","AliAnalysisTaskMuonRefit not created!");
    return kFALSE;
  }
  if (selectPhysics) refit->SelectCollisionCandidates(offlineTriggerMask);
  //refit->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
//  refit->SetFieldPath("$(ALICE_ROOT)/data/maps/mfchebKGI_sym.root");
//  refit->SetFieldPath("alien:///alice/cern.ch/user/p/ppillot/FieldRS/mfchebKGI_symMW.root");
//  Double_t valNB[10] = {0.45, 0.45, 0.65, 0.65, 0.65, 0.65, 0.55, 0.55, 0.55, 0.55};
//  Double_t valB[10] = {0.15, 0.15, 0.3, 0.3, 0.3, 0.3, 0.15, 0.15, 0.15, 0.15};
//  Double_t valNB[10] = {0.14, 0.14, 0.2, 0.2, 0.2, 0.2, 0.17, 0.17, 0.17, 0.17};
//  Double_t valB[10] = {0.1, 0.1, 0.2, 0.2, 0.2, 0.2, 0.1, 0.1, 0.1, 0.1};
//  refit->ResetClusterResolution(valNB, valB);
//  refit->SetAlignStorage("", 5);
//  refit->ReAlign("", 6, -1, "");
//  refit->ReAlign("", -1, -1, "alien://folder=/alice/cern.ch/user/h/hupereir/CDB/LHC15_realign_all_4_dca");
  refit->RemoveMonoCathodClusters(kTRUE, kFALSE);
  refit->TagBadTracks();
  
  // third task (physics after refit)
  AliAnalysisTaskMuonPhysics* physics2 = AddTaskMuonPhysics("after");
  if(!physics2) {
    Error("CreateAnalysisTrain","AliAnalysisTaskMuonPhysics not created!");
    return kFALSE;
  }
  if (selectPhysics) physics2->SelectCollisionCandidates(offlineTriggerMask);
//  if (selectPhysics) physics2->SelectCollisionCandidates(AliVEvent::kMUU7);
  if (selectCentrality) physics2->SelectCentrality(0., 90.);
  physics2->SetMuonTrackCuts(trackCuts);
  /*
  // third task (physics after refit)
  AliAnalysisTaskMuonPhysics* physics22 = AddTaskMuonPhysics("after_LowPt7");
  if(!physics22) {
    Error("CreateAnalysisTrain","AliAnalysisTaskMuonPhysics not created!");
    return kFALSE;
  }
//  if (selectPhysics) physics22->SelectCollisionCandidates(offlineTriggerMask);
  if (selectPhysics) physics22->SelectCollisionCandidates(AliVEvent::kMuonLikeLowPt7 | AliVEvent::kMuonUnlikeLowPt7 | AliVEvent::kMuonSingleLowPt7);
  if (selectCentrality) physics22->SelectCentrality(0., 90.);
  physics22->SetMuonTrackCuts(trackCuts);
  */
  // fourth task (physics after refit)
  AliAnalysisTaskMuonPhysics* physics3 = AddTaskMuonPhysics("bad");
  if(!physics3) {
    Error("CreateAnalysisTrain","AliAnalysisTaskMuonPhysics not created!");
    return kFALSE;
  }
  if (selectPhysics) physics3->SelectCollisionCandidates(offlineTriggerMask);
//  if (selectPhysics) physics3->SelectCollisionCandidates(AliVEvent::kMUU7);
  if (selectCentrality) physics3->SelectCentrality(0., 90.);
  physics3->SetMuonTrackCuts(trackCuts);
  physics3->SelectBadTracks();
  /*
  // fourth task (physics after refit)
  AliAnalysisTaskMuonPhysics* physics32 = AddTaskMuonPhysics("bad_LowPt7");
  if(!physics32) {
    Error("CreateAnalysisTrain","AliAnalysisTaskMuonPhysics not created!");
    return kFALSE;
  }
  //  if (selectPhysics) physics32->SelectCollisionCandidates(offlineTriggerMask);
  if (selectPhysics) physics32->SelectCollisionCandidates(AliVEvent::kMuonLikeLowPt7 | AliVEvent::kMuonUnlikeLowPt7 | AliVEvent::kMuonSingleLowPt7);
  if (selectCentrality) physics32->SelectCentrality(0., 90.);
  physics32->SetMuonTrackCuts(trackCuts);
  physics32->SelectBadTracks();
  */
  return kTRUE;
}

