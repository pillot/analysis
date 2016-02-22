/*
 *  runJPsiAccEffCorr.C
 *  aliroot_dev
 *
 *  Created by philippe pillot on 21/12/12.
 *  Copyright 2011 SUBATECH. All rights reserved.
 *
 */


//______________________________________________________________________________
void runJPsiAccEffCorr2(TString smode = "local", TString inputFileName = "AliAOD.root",
                        Bool_t applyPhysicsSelection = kTRUE, Bool_t embedding = kTRUE)
{
  /// Compute the JPsi acc*eff correction
  
  // --- general analysis setup ---
  TString rootVersion = "";
  TString alirootVersion = "";
  TString aliphysicsVersion = "v5-08-00-01-1";
  TString extraLibs="";
  TString extraIncs="include";
  TString extraTasks="AliAnalysisTaskJPsiAccEffCorr2";
  TString extraPkgs="";
  TList pathList; pathList.SetOwner();
  pathList.Add(new TObjString("$WORK/Macros/MuonEfficiency"));
  TList fileList; fileList.SetOwner();
  fileList.Add(new TObjString("runJPsiAccEffCorr2.C"));
  fileList.Add(new TObjString("AddTaskJPsiAccEffCorr2.C"));
  fileList.Add(new TObjString("AliAnalysisTaskJPsiAccEffCorr2.cxx"));
  fileList.Add(new TObjString("AliAnalysisTaskJPsiAccEffCorr2.h"));
  
  // --- grid specific setup ---
  TString dataDir = "/alice/sim/2016/LHC16b1";
  TString dataPattern = "p80/AOD/*AliAOD.root";
  TString runFormat = "%d";
  TString outDir = "Sim/LHC15o/EmbedJPsi/AccEff/AOD";
  TString analysisMacroName = "AccEff";
  Int_t ttl = 30000;
  Int_t maxFilesPerJob = 20;
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
    
    // get data type
    TString dataType = GetDataType(mode, inputFileName, dataPattern);
    CreateAnalysisTrain(dataType, applyPhysicsSelection, embedding);
    
    if (smode == "saf3" && splitDataset) AliAnalysisManager::GetAnalysisManager()->SetSkipTerminate(kTRUE);
    
    RunAnalysis(smode, inputFileName, rootVersion, alirootVersion, aliphysicsVersion, extraLibs, extraIncs, extraTasks, extraPkgs, dataDir, dataPattern, outDir, analysisMacroName, runFormat, ttl, maxFilesPerJob, maxMergeFiles, maxMergeStages);
    
  }
  
}

//______________________________________________________________________________
void CreateAnalysisTrain(TString dataType, Bool_t applyPhysicsSelection, Bool_t embedding)
{
  /// create the analysis train and configure it
  
  // analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("JPsiAccEffCorrAnalysis");
  
  // ESD or AOD handler
  if (dataType == "ESD") {
    AliESDInputHandler* esdH = new AliESDInputHandler();
    esdH->SetReadFriends(kFALSE);
    mgr->SetInputEventHandler(esdH);
    AliMCEventHandler* mcH = new AliMCEventHandler();
    mgr->SetMCtruthEventHandler(mcH);
  } else if (dataType == "AOD") {
    AliInputEventHandler* aodH = new AliAODInputHandler;
    mgr->SetInputEventHandler(aodH);
  } else {
    Error("CreateAnalysisTrain","Unknown data type. Cannot define input handler!");
    return;
  }
  
  UInt_t offlineTriggerMask = AliVEvent::kMUSPB;
  if (dataType == "ESD") {
    
    // event selection
    if (applyPhysicsSelection) {
      gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
      AliPhysicsSelectionTask* physicsSelection = embedding ? AddTaskPhysicsSelection() : AddTaskPhysicsSelection(kTRUE);
      if(!physicsSelection) {
        Error("CreateAnalysisTrain","AliPhysicsSelectionTask not created!");
        return;
      }
    }
    
    // multiplicity/centrality selection
    gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
    AliMultSelectionTask *mult = AddTaskMultSelection(kFALSE);
    if(!mult) {
      Error("CreateAnalysisTrain","AliMultSelectionTask not created!");
      return;
    }
    mult->SetAlternateOADBforEstimators("LHC15o");
    if (applyPhysicsSelection) mult->SelectCollisionCandidates(offlineTriggerMask);
    
  }
  
  // Acc*Eff results
  gROOT->LoadMacro("AddTaskJPsiAccEffCorr2.C");
  AliAnalysisTaskJPsiAccEffCorr2 *jPsiAccEffCorr = AddTaskJPsiAccEffCorr2();
  if(!jPsiAccEffCorr) {
    Error("CreateAnalysisTrain","AliAnalysisTaskJPsiAccEffCorr2 not created!");
    return;
  }
  if (applyPhysicsSelection) jPsiAccEffCorr->SelectCollisionCandidates(offlineTriggerMask);
  jPsiAccEffCorr->SetTrigLevel(2);
  jPsiAccEffCorr->SetNMatch(2);
  Float_t centBinLowEdge[] = {0.,10.,20.,30.,40.,50.,60.,70.,80.,90.};
  jPsiAccEffCorr->SetCentBins((Int_t)(sizeof(centBinLowEdge)/sizeof(Float_t))-1, centBinLowEdge);
  Float_t pTBinLowEdge[] = {0.,1.,2.,3.,4.,5.,6.,8.};
  jPsiAccEffCorr->SetPtBins((Int_t)(sizeof(pTBinLowEdge)/sizeof(Float_t))-1, pTBinLowEdge);
  Float_t yBinLowEdge[] = {-4.,-3.5,-3.,-2.5};
  jPsiAccEffCorr->SetYBins((Int_t)(sizeof(yBinLowEdge)/sizeof(Float_t))-1, yBinLowEdge);
  jPsiAccEffCorr->LoadRunWeights("../../runWeightCMUL7.txt");
  SetSigWeights(jPsiAccEffCorr);
  
}

//______________________________________________________________________________
void SetSigWeights(TObject* jPsiAccEffCorr)
{
  /// set the number of measured JPsi per pt/y bin or the <Ncoll> used to weight the acc*eff versus centrality
  
  // Ncoll per centrality bin (width=10) to weight the acc*eff calculation
  Float_t nColl10CentBinLowEdge[] = {0., 10., 20., 30., 40., 50., 60., 70., 80., 90.};
  Double_t nColl10[] = {1636., 1001., 601., 344., 183., 90., 40., 16., 6.};
  (static_cast<AliAnalysisTaskJPsiAccEffCorr2*>(jPsiAccEffCorr))->SetSigWeights("nColl (width=10)", 0., 8., -4., -2.5, (Int_t)(sizeof(nColl10)/sizeof(Double_t)), nColl10CentBinLowEdge, nColl10, kFALSE);
  
  // NJpsi per centrality bin (width=10) integrated over pt and y to weight the acc*eff calculation
  Float_t centBinLowEdge00[] = {0., 10., 20., 30., 40., 50., 60., 70., 80., 90.};
  Double_t nJpsi00[] = {118139., 73280., 46010., 25294., 15034., 7884., 4023., 2062., 910.};
  (static_cast<AliAnalysisTaskJPsiAccEffCorr2*>(jPsiAccEffCorr))->SetSigWeights("nJPsi", 0., 8., -4., -2.5, (Int_t)(sizeof(nJpsi00)/sizeof(Double_t)), centBinLowEdge00, nJpsi00, kTRUE);
  /*
  // Ncoll per centrality bin (width=10) to weight the acc*eff calculation
  Float_t nColl10CentBinLowEdge[9] = {0., 10., 20., 30., 40., 50., 60., 70., 80.};
  Double_t nColl10[8] = {1502.7, 923.26, 558.68, 321.20, 171.67, 85.13, 38.51, 15.78};
  (static_cast<AliAnalysisTaskJPsiAccEffCorr2*>(jPsiAccEffCorr))->SetSigWeights("nColl (width=10)", 0., 1.e9, -4., -2.5, (Int_t)(sizeof(nColl10)/sizeof(Double_t)), nColl10CentBinLowEdge, nColl10, kFALSE);
  
  // Ncoll per centrality bin (width=5) to weight the acc*eff calculation
  Float_t nColl5CentBinLowEdge[17] = {0., 5., 10., 15., 20., 25., 30., 35., 40., 45., 50., 55., 60., 65., 70., 75., 80.};
  Double_t nColl5[16] = {1686.87, 1319.89, 1031.9, 807.90, 627.99, 483.95, 369.13, 274.03, 199.30, 143.45, 100.54, 68.82, 46.09, 29.70, 18.80, 11.95};
  (static_cast<AliAnalysisTaskJPsiAccEffCorr2*>(jPsiAccEffCorr))->SetSigWeights("nColl (width=5)", 0., 1.e9, -4., -2.5, (Int_t)(sizeof(nColl5)/sizeof(Double_t)), nColl5CentBinLowEdge, nColl5, kFALSE);
  
  // integrated pt/y
  Float_t centBinLowEdge00[6] = {0., 10., 20., 30., 50., 80.};
  Double_t nJpsi00[5] = {973., 573., 406., 306., 119.};
  (static_cast<AliAnalysisTaskJPsiAccEffCorr2*>(jPsiAccEffCorr))->SetSigWeights("nJPsi", 0., 1.e9, -4., -2.5, (Int_t)(sizeof(nJpsi00)/sizeof(Double_t)), centBinLowEdge00, nJpsi00, kTRUE);
  
  // integrated over pt / -4 < y < -3.25
  Float_t centBinLowEdge01[6] = {0., 10., 20., 30., 50., 80.};
  Double_t nJpsi01[5] = {377., 257., 142., 126., 45.};
  (static_cast<AliAnalysisTaskJPsiAccEffCorr2*>(jPsiAccEffCorr))->SetSigWeights("nJPsi", 0., 1.e9, -4., -3.25, (Int_t)(sizeof(nJpsi01)/sizeof(Double_t)), centBinLowEdge01, nJpsi01, kTRUE);
  
  // integrated over pt / -3.25 < y < -4
  Float_t centBinLowEdge02[6] = {0., 10., 20., 30., 50., 80.};
  Double_t nJpsi02[5] = {586., 327., 260., 182., 79.};
  (static_cast<AliAnalysisTaskJPsiAccEffCorr2*>(jPsiAccEffCorr))->SetSigWeights("nJPsi", 0., 1.e9, -3.25, -2.5, (Int_t)(sizeof(nJpsi02)/sizeof(Double_t)), centBinLowEdge02, nJpsi02, kTRUE);
  
  // 0 < pt < 3 / integrated over y
  Float_t centBinLowEdge10[6] = {0., 10., 20., 30., 50., 80.};
  Double_t nJpsi10[5] = {785., 450., 323., 215., 81.};
  (static_cast<AliAnalysisTaskJPsiAccEffCorr2*>(jPsiAccEffCorr))->SetSigWeights("nJPsi", 0., 3., -4., -2.5, (Int_t)(sizeof(nJpsi10)/sizeof(Double_t)), centBinLowEdge10, nJpsi10, kTRUE);
  
  // 3 < pt < inf / integrated over y
  Float_t centBinLowEdge20[4] = {0., 20., 40., 80.};
  Double_t nJpsi20[3] = {340., 151., 56.};
  (static_cast<AliAnalysisTaskJPsiAccEffCorr2*>(jPsiAccEffCorr))->SetSigWeights("nJPsi", 3., 1.e9, -4., -2.5, (Int_t)(sizeof(nJpsi20)/sizeof(Double_t)), centBinLowEdge20, nJpsi20, kTRUE);
  */
  
}

