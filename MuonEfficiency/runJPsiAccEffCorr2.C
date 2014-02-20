/*
 *  runJPsiAccEffCorr.C
 *  aliroot_dev
 *
 *  Created by philippe pillot on 21/12/12.
 *  Copyright 2011 SUBATECH. All rights reserved.
 *
 */


TString rootVersion = "v5-34-05";
TString alirootVersion = "v5-04-46-AN";
TString dataDir = "/alice/data/2010/LHC10h";
TString dataPattern = "pass2/*AliESDs.root";
TString runFormat = "%09d";
TString outDir = "Data/LHC10h/pass2/Eff/pDCAChi2";
Int_t ttl = 30000;
Int_t maxFilesPerJob = 100;
Int_t maxMergeFiles = 10;
Int_t maxMergeStages = 2;

//______________________________________________________________________________
void runJPsiAccEffCorr2(TString smode = "local", TString inputFileName = "AliAOD.root")
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
  fileList.Add(new TObjString("runJPsiAccEffCorr2.C"));
  fileList.Add(new TObjString("AliAnalysisTaskJPsiAccEffCorr2.cxx"));
  fileList.Add(new TObjString("AliAnalysisTaskJPsiAccEffCorr2.h"));
  CopyFileLocally(pathList, fileList);
  
  // --- prepare environment ---
  TString extraLibs="";
  TString extraIncs="";
  TString extraTasks="AliAnalysisTaskJPsiAccEffCorr2";
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
  AliAnalysisTaskJPsiAccEffCorr2 *jPsiAccEffCorr = new AliAnalysisTaskJPsiAccEffCorr2("JPsiAccEffCorr");
  //  jPsiAccEffCorr->SelectCollisionCandidates(AliVEvent::kMB | AliVEvent::kMUON);
  //  jPsiAccEffCorr->SelectCollisionCandidates(AliVEvent::kMB);
  //  jPsiAccEffCorr->SelectCollisionCandidates(AliVEvent::kAny);
  jPsiAccEffCorr->SetTrigLevel(2);
  jPsiAccEffCorr->SetNMatch(2);
  //jPsiAccEffCorr->SetMuLowPtCut(1.13);
  Float_t centBinLowEdge[] = {-999., 999.};
  jPsiAccEffCorr->SetCentBins((Int_t)(sizeof(centBinLowEdge)/sizeof(Float_t))-1, centBinLowEdge);
  //Float_t pTBinLowEdge[] = {0., 0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5., 6., 7., 8., 9., 11., 13., 15., 1000000.};
  //Float_t pTBinLowEdge[] = {0., 1., 2., 3., 4., 5., 6., 7., 9., 11., 15., 1000000.}; // Igor's ranges
  Float_t pTBinLowEdge[] = {0., 1., 2., 3., 4., 5., 6., 8., 10., 15., 1000000.}; // Roberta's ranges
  //  Float_t pTBinLowEdge[] = {0., 2., 5., 8.};
  jPsiAccEffCorr->SetPtBins((Int_t)(sizeof(pTBinLowEdge)/sizeof(Float_t))-1, pTBinLowEdge);
  Float_t yBinLowEdge[] = {-4., -3.75, -3.5, -3.25, -3., -2.75, -2.5};
  //Float_t yBinLowEdge[] = {-4., -3.81, -3.62, -3.43}; // LHC13de
  //Float_t yBinLowEdge[] = {-3.07, -2.88, -2.69, -2.5}; // LHC13f
  //  Float_t yBinLowEdge[] = {-4., -3.25, -2.5};
  jPsiAccEffCorr->SetYBins((Int_t)(sizeof(yBinLowEdge)/sizeof(Float_t))-1, yBinLowEdge);
//  jPsiAccEffCorr->LoadRunWeights("../../runWeights.txt");
//  SetSigWeights(jPsiAccEffCorr);
  
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

//______________________________________________________________________________
void SetSigWeights(TObject* jPsiAccEffCorr)
{
  /// set the number of measured JPsi per pt/y bin or the <Ncoll> used to weight the acc*eff versus centrality
  
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
  
}

