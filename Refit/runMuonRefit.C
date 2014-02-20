/*
 *  runMuonRefit.C
 *  aliroot
 *
 *  Created by Philippe Pillot on 02/11/10.
 *  Copyright 2010 SUBATECH. All rights reserved.
 *
 */


TString rootVersion = "v5-34-05";
TString alirootVersion = "v5-04-58-AN";
TString dataDir = "/alice/data/2013/LHC13f";
TString dataPattern = "muon_pass2/*AliESDs.root";
TString runFormat = "%09d";
TString outDir = "Data/LHC13f/muon_pass2/Refit/trg6sigma";
Int_t ttl = 30000;
Int_t maxFilesPerJob = 100;
Int_t maxMergeFiles = 10;
Int_t maxMergeStages = 1;

//______________________________________________________________________________
void runMuonRefit(TString smode = "local", TString inputFileName = "AliESDs.root",
		  Double_t resNB = -1., Double_t resB = -1., Double_t sigmaCut = -1.,
		  Double_t sigmaCutTrig = 6, Bool_t selectPhysics = kTRUE, Bool_t selectCentrality = kFALSE)
{
  /// Refit the tracks with new recoParam and/or new alignment
  /// check the effect on track quality, efficiency and physical quantities
  
  gROOT->LoadMacro("$WORK/Macros/Facilities/runTaskFacilities.C");
  
  // --- Check runing mode ---
  Int_t mode = GetMode(smode, inputFileName);
  if(mode < 0){
    Error("runMuonRefit","Please provide either an ESD root file a collection of ESDs or a dataset.");
    return;
  }
  
  // --- copy files needed for this analysis ---
  TList pathList; pathList.SetOwner();
  pathList.Add(new TObjString("$WORK/Macros/MuonPhysics"));
  pathList.Add(new TObjString("$WORK/Macros/Refit"));
  pathList.Add(new TObjString("$ALICE_ROOT/PWG/muondep"));
  TList fileList; fileList.SetOwner();
  fileList.Add(new TObjString("runMuonRefit.C"));
  fileList.Add(new TObjString("AddTaskMuonPhysics.C"));
  fileList.Add(new TObjString("AliAnalysisTaskMuonPhysics.cxx"));
  fileList.Add(new TObjString("AliAnalysisTaskMuonPhysics.h"));
  fileList.Add(new TObjString("AliAnalysisTaskMuonRefit.cxx"));
  fileList.Add(new TObjString("AliAnalysisTaskMuonRefit.h"));
  CopyFileLocally(pathList, fileList);
  
  // --- prepare environment ---
  TString extraLibs="RAWDatabase:CDB:STEER:MUONcore:MUONmapping:MUONcalib:MUONgeometry:MUONtrigger:MUONraw:MUONbase:MUONrec:CORRFW:PWGmuon";
  TString extraIncs="include:MUON:MUON/mapping:PWG/muon";
  TString extraTasks="AliAnalysisTaskMuonPhysics:AliAnalysisTaskMuonRefit";
  LoadAlirootLocally(extraLibs, extraIncs, extraTasks);
  AliAnalysisGrid *alienHandler = 0x0;
  if (mode == kProof || mode == kProofLite) LoadAlirootOnProof(smode, rootVersion, alirootVersion, extraLibs, extraIncs, extraTasks, kTRUE);
  else if (mode == kGrid || mode == kTerminate) {
    TString analysisMacroName = "Refit";
    alienHandler = static_cast<AliAnalysisGrid*>(CreateAlienHandler(smode, rootVersion, alirootVersion, inputFileName, dataDir, dataPattern, outDir, extraLibs, extraIncs, extraTasks, analysisMacroName, runFormat, ttl, maxFilesPerJob, maxMergeFiles, maxMergeStages));
    if (!alienHandler) return;
  }
  
  // --- Create the analysis train ---
  if (!CreateAnalysisTrain(resNB, resB, sigmaCut, sigmaCutTrig, selectPhysics, selectCentrality, alienHandler)) return;
  
  // --- Create input object ---
  TObject* inputObj = CreateInputObject(mode, inputFileName);
  
  // --- start analysis ---
  StartAnalysis(mode, inputObj);
  
}

//______________________________________________________________________________
Bool_t CreateAnalysisTrain(Double_t resNB, Double_t resB, Double_t sigmaCut, Double_t sigmaCutTrig,
			   Bool_t selectPhysics, Bool_t selectCentrality, TObject* alienHandler = 0x0)
{
  /// create the analysis train and configure it
  
  // analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("MuonRefitAnalysis");
  //mgr->SetAutoBranchLoading(kFALSE);
  //mgr->SetCacheSize(0);
  
  // Connect plugin to the analysis manager if any
  if (alienHandler) mgr->SetGridHandler(static_cast<AliAnalysisGrid*>(alienHandler));
  
  // ESD input
  AliESDInputHandler* esdH = new AliESDInputHandler();
  esdH->SetReadFriends(kFALSE);
  esdH->SetInactiveBranches("*");
  esdH->SetActiveBranches("MuonTracks MuonClusters MuonPads AliESDRun. AliESDHeader. AliMultiplicity. AliESDFMD. AliESDVZERO. SPDVertex. PrimaryVertex. AliESDZDC. AliESDTZERO.");
  mgr->SetInputEventHandler(esdH);
  
  // event selection
  if (selectPhysics) {
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
    AliPhysicsSelectionTask* physicsSelection = AddTaskPhysicsSelection();
    if(!physicsSelection) {
      Error("CreateAnalysisTrain","AliPhysicsSelectionTask not created!");
      return kFALSE;
    }
  }
  
  // centrality selection
  if (selectCentrality) {
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskCentrality.C");
    AliCentralitySelectionTask *taskCentrality = AddTaskCentrality();
    if(!taskCentrality) {
      Error("CreateAnalysisTrain","AliCentralitySelectionTask not created!");
      return kFALSE;
    }
//    if (selectPhysics) taskCentrality->SelectCollisionCandidates(AliVEvent::kMUL7 | AliVEvent::kMUU7 |
//								 AliVEvent::kMuonLikeLowPt8 | AliVEvent::kMuonUnlikeLowPt8 |
//								 AliVEvent::kMUS7 | AliVEvent::kMuonSingleLowPt8);
//      if (selectPhysics) taskCentrality->SelectCollisionCandidates(AliVEvent::kMUS7 | AliVEvent::kMuonSingleLowPt8);
//      if (selectPhysics) taskCentrality->SelectCollisionCandidates(AliVEvent::kMUSH7 | AliVEvent::kMUSHPB | AliVEvent::kMuonSingleHighPt8);
//    if (selectPhysics) taskCentrality->SelectCollisionCandidates(AliVEvent::kAny);
    if (selectPhysics) taskCentrality->SelectCollisionCandidates(AliVEvent::kMUU7);
  }
  
  // track selection
  AliMuonTrackCuts trackCuts("stdCuts", "stdCuts");
  trackCuts.SetAllowDefaultParams();
  trackCuts.SetFilterMask(AliMuonTrackCuts::kMuMatchLpt | AliMuonTrackCuts::kMuEta |
			  AliMuonTrackCuts::kMuThetaAbs | AliMuonTrackCuts::kMuPdca);
  //trackCuts.SetIsMC();
  
  // first task (physics before refit)
  gROOT->LoadMacro("AddTaskMuonPhysics.C");
  AliAnalysisTaskMuonPhysics* physics1 = AddTaskMuonPhysics("before");
  if(!physics1) {
    Error("CreateAnalysisTrain","AliAnalysisTaskMuonPhysics not created!");
    return kFALSE;
  }
//  if (selectPhysics) physics1->SelectCollisionCandidates(AliVEvent::kCMUS5 | AliVEvent::kMUL7 | AliVEvent::kMUU7 | AliVEvent::kMUS7);
//  if (selectPhysics) physics1->SelectCollisionCandidates(AliVEvent::kMUL7 | AliVEvent::kMUU7 |
//							 AliVEvent::kMuonLikeLowPt8 | AliVEvent::kMuonUnlikeLowPt8 |
//							 AliVEvent::kMUS7 | AliVEvent::kMuonSingleLowPt8);
//  if (selectPhysics) physics1->SelectCollisionCandidates(AliVEvent::kMUS7 | AliVEvent::kMuonSingleLowPt8);
//  if (selectPhysics) physics1->SelectCollisionCandidates(AliVEvent::kMUSH7 | AliVEvent::kMUSHPB | AliVEvent::kMuonSingleHighPt8);
//  if (selectPhysics) physics1->SelectCollisionCandidates(AliVEvent::kCMUS5 | AliVEvent::kMUS7);
//  if (selectPhysics) physics1->SelectCollisionCandidates(AliVEvent::kMUU7);
//  if (selectPhysics) physics1->SelectCollisionCandidates(AliVEvent::kAny);
  if (selectPhysics) physics1->SelectCollisionCandidates(AliVEvent::kMUU7);
  if (selectCentrality) physics1->SelectCentrality(0., 90.);
  physics1->SetMuonTrackCuts(trackCuts);
  
  // second task (refit)
  gROOT->LoadMacro("$ALICE_ROOT/PWG/muondep/AddTaskMuonRefit.C");
  AliAnalysisTaskMuonRefit* refit = AddTaskMuonRefit(resNB, resB, kTRUE, sigmaCut, sigmaCutTrig);
  if(!refit) {
    Error("CreateAnalysisTrain","AliAnalysisTaskMuonRefit not created!");
    return kFALSE;
  }
//  if (selectPhysics) refit->SelectCollisionCandidates(AliVEvent::kCMUS5 | AliVEvent::kMUL7 | AliVEvent::kMUU7 | AliVEvent::kMUS7);
//  if (selectPhysics) refit->SelectCollisionCandidates(AliVEvent::kMUL7 | AliVEvent::kMUU7 |
//						      AliVEvent::kMuonLikeLowPt8 | AliVEvent::kMuonUnlikeLowPt8 |
//						      AliVEvent::kMUS7 | AliVEvent::kMuonSingleLowPt8);
//  if (selectPhysics) refit->SelectCollisionCandidates(AliVEvent::kMUS7 | AliVEvent::kMuonSingleLowPt8);
//  if (selectPhysics) refit->SelectCollisionCandidates(AliVEvent::kMUSH7 | AliVEvent::kMUSHPB | AliVEvent::kMuonSingleHighPt8);
//  if (selectPhysics) refit->SelectCollisionCandidates(AliVEvent::kCMUS5 | AliVEvent::kMUS7);
//  if (selectPhysics) refit->SelectCollisionCandidates(AliVEvent::kMUU7);
//  if (selectPhysics) refit->SelectCollisionCandidates(AliVEvent::kAny);
  if (selectPhysics) refit->SelectCollisionCandidates(AliVEvent::kMUU7);
  if (selectCentrality) refit->SelectCentrality(0., 90.);
  //refit->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
//  refit->SetFieldPath("$(ALICE_ROOT)/data/maps/mfchebKGI_sym.root");
//  refit->SetFieldPath("alien:///alice/cern.ch/user/p/ppillot/FieldRS/mfchebKGI_symMW.root");
//  refit->ReAlign("alien://folder=/alice/cern.ch/user/p/ppillot/OCDB2011_Align1", "alien://folder=/alice/cern.ch/user/p/ppillot/OCDB2011_Align2");
//  refit->ReAlign("alien://folder=/alice/cern.ch/user/p/ppillot/OCDB2011_Align1", "alien://folder=/alice/cern.ch/user/s/shahoian/CorG4Fresmx05y015");
//  refit->ReAlign("alien://folder=/alice/cern.ch/user/p/ppillot/OCDB2011_Align1", "alien://folder=/alice/cern.ch/user/s/shahoian/CorG4Gresmx05y015");
//  refit->ReAlign("alien://folder=/alice/cern.ch/user/p/ppillot/OCDB2011_Align1", "alien://folder=/alice/cern.ch/user/j/jcastill/CorG4Fresmx05y015");
//  refit->ReAlign("alien://folder=/alice/cern.ch/user/p/ppillot/OCDB2011_Align1", "alien://folder=/alice/cern.ch/user/j/jcastill/pbpb11wrk/CorG4Fresmx05y015_pp2PbPb");
//  refit->ReAlign("", "alien://folder=/alice/cern.ch/user/h/hupereir/CDB/LHC12_ReAlign_new_0");
//  refit->ReAlign("alien://folder=/alice/cern.ch/user/p/ppillot/OCDB2012", "alien://folder=/alice/cern.ch/user/h/hupereir/CDB/LHC12_ReAlign_new_1");
//  refit->ReAlign("alien://folder=/alice/cern.ch/user/p/ppillot/OCDB2012", "");
//  refit->ReAlign("alien://folder=/alice/simulation/2008/v4-15-Release/Residual", "alien://folder=/alice/simulation/2008/v4-15-Release/Residual");
//    refit->SetAlignStorage("alien://folder=/alice/simulation/2008/v4-15-Release/Residual");
//  Double_t valNB[10] = {0.45, 0.45, 0.65, 0.65, 0.65, 0.65, 0.55, 0.55, 0.55, 0.55};
//  Double_t valB[10] = {0.15, 0.15, 0.3, 0.3, 0.3, 0.3, 0.15, 0.15, 0.15, 0.15};
//  Double_t valNB[10] = {0.14, 0.14, 0.2, 0.2, 0.2, 0.2, 0.17, 0.17, 0.17, 0.17};
//  Double_t valB[10] = {0.1, 0.1, 0.2, 0.2, 0.2, 0.2, 0.1, 0.1, 0.1, 0.1};
//  refit->ResetClusterResolution(valNB, valB);
//  refit->RemoveMonoCathodClusters(kTRUE, kFALSE);
  refit->TagBadTracks();
  
  // third task (physics after refit)
  AliAnalysisTaskMuonPhysics* physics2 = AddTaskMuonPhysics("after");
  if(!physics2) {
    Error("CreateAnalysisTrain","AliAnalysisTaskMuonPhysics not created!");
    return kFALSE;
  }
//  if (selectPhysics) physics2->SelectCollisionCandidates(AliVEvent::kCMUS5 | AliVEvent::kMUL7 | AliVEvent::kMUU7 | AliVEvent::kMUS7);
//  if (selectPhysics) physics2->SelectCollisionCandidates(AliVEvent::kMUL7 | AliVEvent::kMUU7 |
//							 AliVEvent::kMuonLikeLowPt8 | AliVEvent::kMuonUnlikeLowPt8 |
//							 AliVEvent::kMUS7 | AliVEvent::kMuonSingleLowPt8);
//  if (selectPhysics) physics2->SelectCollisionCandidates(AliVEvent::kMUS7 | AliVEvent::kMuonSingleLowPt8);
//  if (selectPhysics) physics2->SelectCollisionCandidates(AliVEvent::kMUSH7 | AliVEvent::kMUSHPB | AliVEvent::kMuonSingleHighPt8);
//  if (selectPhysics) physics2->SelectCollisionCandidates(AliVEvent::kCMUS5 | AliVEvent::kMUS7);
//  if (selectPhysics) physics2->SelectCollisionCandidates(AliVEvent::kMUU7);
//  if (selectPhysics) physics2->SelectCollisionCandidates(AliVEvent::kAny);
  if (selectPhysics) physics2->SelectCollisionCandidates(AliVEvent::kMUU7);
  if (selectCentrality) physics2->SelectCentrality(0., 90.);
  physics2->SetMuonTrackCuts(trackCuts);
  
  // third task (physics after refit)
  AliAnalysisTaskMuonPhysics* physics3 = AddTaskMuonPhysics("bad");
  if(!physics3) {
    Error("CreateAnalysisTrain","AliAnalysisTaskMuonPhysics not created!");
    return kFALSE;
  }
//  if (selectPhysics) physics3->SelectCollisionCandidates(AliVEvent::kCMUS5 | AliVEvent::kMUL7 | AliVEvent::kMUU7 | AliVEvent::kMUS7);
//  if (selectPhysics) physics3->SelectCollisionCandidates(AliVEvent::kMUL7 | AliVEvent::kMUU7 |
//							 AliVEvent::kMuonLikeLowPt8 | AliVEvent::kMuonUnlikeLowPt8 |
//							 AliVEvent::kMUS7 | AliVEvent::kMuonSingleLowPt8);
//  if (selectPhysics) physics3->SelectCollisionCandidates(AliVEvent::kMUS7 | AliVEvent::kMuonSingleLowPt8);
//  if (selectPhysics) physics3->SelectCollisionCandidates(AliVEvent::kMUSH7 | AliVEvent::kMUSHPB | AliVEvent::kMuonSingleHighPt8);
//  if (selectPhysics) physics3->SelectCollisionCandidates(AliVEvent::kCMUS5 | AliVEvent::kMUS7);
//  if (selectPhysics) physics3->SelectCollisionCandidates(AliVEvent::kMUU7);
//  if (selectPhysics) physics3->SelectCollisionCandidates(AliVEvent::kAny);
  if (selectPhysics) physics3->SelectCollisionCandidates(AliVEvent::kMUU7);
  if (selectCentrality) physics3->SelectCentrality(0., 90.);
  physics3->SetMuonTrackCuts(trackCuts);
  physics3->SelectBadTracks();
  
  return kTRUE;
}

