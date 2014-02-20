/*
 *  runJPsiAccEffCorr.C
 *  aliroot_dev
 *
 *  Created by philippe pillot on 06/03/13.
 *  Copyright 2011 SUBATECH. All rights reserved.
 *
 */


TString rootVersion = "v5-34-05";
TString alirootVersion = "v5-04-38-AN";
TString dataDir = "/alice/data/2010/LHC10h";
TString dataPattern = "pass2/*AliESDs.root";
TString runFormat = "%09d";
TString outDir = "Data/LHC10h/pass2/Eff/pDCAChi2";
Int_t ttl = 30000;
Int_t maxFilesPerJob = 100;
Int_t maxMergeFiles = 10;
Int_t maxMergeStages = 2;

// generator parameters used in the simulation
/*
// tune0 LHC13de
Double_t oldPtParam[6] = {371.909, 0.84614, 0.560486, 9.34831, 0.000474983, -0.853963};
Bool_t oldFixPtParam[6] = {kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE};
Double_t newPtParam[6] = {371.909, 0.84614, 0.560486, 9.34831, 0.000474983, -0.853963};
Bool_t newFixPtParam[6] = {kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE};
*//*
// tune1 LHC13de
Double_t oldPtParam[6] = {371.665, 0.845642, 0.56192, 9.34859, 0.000474519, -0.851091};
Bool_t oldFixPtParam[6] = {kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE};
Double_t newPtParam[6] = {371.665, 0.845642, 0.56192, 9.34859, 0.000474519, -0.851091};
Bool_t newFixPtParam[6] = {kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE};
*//*
// tune0 LHC13f
Double_t oldPtParam[6] = {522.811, 0.997725, 0.705636, 8.52259, 0., -1.};
Bool_t oldFixPtParam[6] = {kFALSE, kFALSE, kFALSE, kFALSE, kTRUE, kTRUE};
Double_t newPtParam[6] = {522.811, 0.997725, 0.705636, 8.52259, 0.0001, -1.};
Bool_t newFixPtParam[6] = {kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE};
*/
// tune1 LHC13f
Double_t oldPtParam[6] = {455.614, 0.942071, 0.706755, 8.69856, 0.000168775, -0.925487};
Bool_t oldFixPtParam[6] = {kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE};
Double_t newPtParam[6] = {455.614, 0.942071, 0.706755, 8.69856, 0.000168775, -0.925487};
Bool_t newFixPtParam[6] = {kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE};

Double_t ptRange[2] = {0.8, 999.};

/*
// tune0 LHC13de
Double_t oldYParam[8] = {0.539134, 1, 0, 0.0499378, 0, -0.00450342, 0, 2};
Bool_t oldFixYParam[8] = {kFALSE, kTRUE, kTRUE, kFALSE, kTRUE, kFALSE, kTRUE, kTRUE};
Double_t newYParam[8] = {0.539134, 1, 0, 0.0499378, 0, -0.00450342, 0, 2};
Bool_t newFixYParam[8] = {kFALSE, kTRUE, kTRUE, kFALSE, kTRUE, kFALSE, kTRUE, kTRUE};
*//*
// tune1 LHC13de
Double_t oldYParam[8] = {0.777922, 1, 0, -0.0184202, 0, -0.00107081, 0, 2};
Bool_t oldFixYParam[8] = {kFALSE, kTRUE, kTRUE, kFALSE, kTRUE, kFALSE, kTRUE, kTRUE};
Double_t newYParam[8] = {0.777922, 1, 0, -0.0184202, 0, -0.00107081, 0, 2};
Bool_t newFixYParam[8] = {kFALSE, kTRUE, kTRUE, kFALSE, kTRUE, kFALSE, kTRUE, kTRUE};
*//*
// tune0 LHC13f
Double_t oldYParam[8] = {1.75646, 1., 8.70262e-05, -0.129939, -0.0190949, 0., 0., 2.};
Bool_t oldFixYParam[8] = {kFALSE, kTRUE, kFALSE, kFALSE, kFALSE, kTRUE, kTRUE, kTRUE};
Double_t newYParam[8] = {1.5712, 1., 0., -0.0893785, 0., 0.00228603, 0., 2.};
Bool_t newFixYParam[8] = {kFALSE, kTRUE, kTRUE, kFALSE, kTRUE, kFALSE, kTRUE, kTRUE};
//Double_t newYParam[8] = {1.8216, 0., 0., 0., 0., 0., 1., 2.0016};
//Bool_t newFixYParam[8] = {kFALSE, kTRUE, kTRUE, kTRUE, kTRUE, kTRUE, kTRUE, kFALSE};
*/
// tune1 LHC13f
Double_t oldYParam[8] = {1.29511, 1., 0., -0.0767846, 0., 0.00176313, 0., 2.};
Bool_t oldFixYParam[8] = {kFALSE, kTRUE, kTRUE, kFALSE, kTRUE, kFALSE, kTRUE, kTRUE};
Double_t newYParam[8] = {1.29511, 1., 0., -0.0767846, 0., 0.00176313, 0., 2.};
Bool_t newFixYParam[8] = {kFALSE, kTRUE, kTRUE, kFALSE, kTRUE, kFALSE, kTRUE, kTRUE};

Double_t yRange[2] = {-4.3, -2.3};


Bool_t isMC = kTRUE;
Bool_t applyPhysicsSelection = kFALSE;

//______________________________________________________________________________
void runGenTuner(TString smode = "local", TString inputFileName = "AliAOD.root",
		 Int_t iStep = -1, char overwrite = '\0')
{
  /// Study the MUON performances
  
  gROOT->LoadMacro("$WORK/Macros/Facilities/runTaskFacilities.C");
  
  // --- Check runing mode ---
  Int_t mode = GetMode(smode, inputFileName);
  if(mode < 0) {
    Error("runGenTuner","Please provide either an AOD root file a collection of AODs or a dataset.");
    return;
  }
  
  // --- copy files needed for this analysis ---
  TList pathList; pathList.SetOwner();
  pathList.Add(new TObjString("$WORK/Macros/MuonEfficiency/GenTuner"));
  TList fileList; fileList.SetOwner();
  fileList.Add(new TObjString("runGenTuner.C"));
  fileList.Add(new TObjString("AddTaskGenTuner.C"));
  fileList.Add(new TObjString("AliAnalysisTaskGenTuner.cxx"));
  fileList.Add(new TObjString("AliAnalysisTaskGenTuner.h"));
  CopyFileLocally(pathList, fileList, overwrite);
  
  // --- prepare environment ---
  TString extraLibs="CORRFW:PWGmuon";
  TString extraIncs="include";
  TString extraTasks="AliAnalysisTaskGenTuner";
  LoadAlirootLocally(extraLibs, extraIncs, extraTasks);
  AliAnalysisGrid *alienHandler = 0x0;
  if (mode == kProof || mode == kProofLite) LoadAlirootOnProof(smode, rootVersion, alirootVersion, extraLibs, extraIncs, extraTasks, kTRUE);
  else if (mode == kGrid || mode == kTerminate) {
    TString analysisMacroName = "GenTuner";
    alienHandler = static_cast<AliAnalysisGrid*>(CreateAlienHandler(smode, rootVersion, alirootVersion, inputFileName, dataDir, dataPattern, outDir, extraLibs, extraIncs, extraTasks, analysisMacroName, runFormat, ttl, maxFilesPerJob, maxMergeFiles, maxMergeStages));
    if (!alienHandler) return;
  }
  
  // --- Create the analysis train ---
  AliAnalysisTaskGenTuner *genTuner = static_cast<AliAnalysisTaskGenTuner*>(CreateAnalysisTrain(alienHandler, iStep));
  if (!genTuner) return;
  
  // --- Create input object ---
  TObject* inputObj = CreateInputObject(mode, inputFileName);
  
  // --- start analysis ---
  StartAnalysis(mode, inputObj);
  
  // --- save fitting functions ---
  TString outFileName = AliAnalysisManager::GetCommonFileName();
  TFile *outFile = (TFile*)gROOT->GetListOfFiles()->FindObject(outFileName.Data());
  if (outFile) outFile->ReOpen("UPDATE");
  else outFile = TFile::Open(outFileName.Data(),"UPDATE");
  if (outFile && outFile->IsOpen()) {
    outFile->Cd(Form("%s:/MUON_GenTuner",outFileName.Data()));
    if (genTuner->GetOldPtFunc()) genTuner->GetOldPtFunc()->Write(0x0, TObject::kOverwrite);
    if (genTuner->GetOldPtFuncMC()) genTuner->GetOldPtFuncMC()->Write(0x0, TObject::kOverwrite);
    if (genTuner->GetNewPtFunc()) genTuner->GetNewPtFunc()->Write(0x0, TObject::kOverwrite);
    if (genTuner->GetOldYFunc()) genTuner->GetOldYFunc()->Write(0x0, TObject::kOverwrite);
    if (genTuner->GetOldYFuncMC()) genTuner->GetOldYFuncMC()->Write(0x0, TObject::kOverwrite);
    if (genTuner->GetNewYFunc()) genTuner->GetNewYFunc()->Write(0x0, TObject::kOverwrite);
    if (genTuner->GetResults()) genTuner->GetResults()->Write(0x0, TObject::kOverwrite);
    if (genTuner->GetRatios()) genTuner->GetRatios()->Write(0x0, TObject::kOverwrite);
    outFile->Close();
  }
  
  // save results of current step if running in a loop
  if (iStep > -1) gSystem->Exec(Form("cp -f AnalysisResults.root Results_step%d.root", iStep));
  
}

//______________________________________________________________________________
TObject* CreateAnalysisTrain(TObject* alienHandler, Int_t iStep)
{
  /// create the analysis train and configure it
  
  // analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("GenTunerAnalysis");
  
  // Connect plugin to the analysis manager if any
  if (alienHandler) mgr->SetGridHandler(static_cast<AliAnalysisGrid*>(alienHandler));
  
  // AOD handler
  AliInputEventHandler* aodH = new AliAODInputHandler;
  mgr->SetInputEventHandler(aodH);
  
  // track selection
  AliMuonTrackCuts trackCuts("stdCuts", "stdCuts");
  trackCuts.SetAllowDefaultParams();
  trackCuts.SetFilterMask(AliMuonTrackCuts::kMuMatchLpt | AliMuonTrackCuts::kMuEta |
			  AliMuonTrackCuts::kMuThetaAbs | AliMuonTrackCuts::kMuPdca);
//  trackCuts.SetFilterMask(AliMuonTrackCuts::kMuMatchHpt | AliMuonTrackCuts::kMuEta |
//			  AliMuonTrackCuts::kMuThetaAbs | AliMuonTrackCuts::kMuPdca);
  if (isMC) trackCuts.SetIsMC(kTRUE);
  
  // generator tuner
  gROOT->LoadMacro("AddTaskGenTuner.C");
  AliAnalysisTaskGenTuner* genTuner = AddTaskGenTuner();
  if(!genTuner) {
    Error("CreateAnalysisTrain","AliAnalysisTaskGenTuner not created!");
    return 0x0;
  }
  if (applyPhysicsSelection) genTuner->SelectCollisionCandidates(AliVEvent::kMUU7);
//  if (applyPhysicsSelection) genTuner->SelectCollisionCandidates(AliVEvent::kMUSH7);
  //genTuner->SelectCentrality(0., 90.);
  genTuner->SetMuonTrackCuts(trackCuts);
  genTuner->SetMuonPtCut(1.);
//  if (isMC) genTuner->SetDataFile("/Users/pillot/Work/Alice/Work/Data/2013/LHC13d/muon_pass2/AOD/GenTuner/pT1GeV/AnalysisResults.root");
//  if (isMC) genTuner->SetDataFile("/Users/pillot/Work/Alice/Work/Data/2013/LHC13d/muon_pass2/AOD/GenTuner/pT2GeV/AnalysisResults.root");
//  if (isMC) genTuner->SetDataFile("/Users/pillot/Work/Alice/Work/Data/2013/LHC13e/muon_pass2/AOD/GenTuner/pT1GeV/AnalysisResults.root");
//  if (isMC) genTuner->SetDataFile("/Users/pillot/Work/Alice/Work/Data/2013/LHC13e/muon_pass2/AOD/GenTuner/pT2GeV/AnalysisResults.root");
//  if (isMC) genTuner->SetDataFile("/Users/pillot/Work/Alice/Work/Data/2013/LHC13de/muon_pass2/AOD/GenTuner/pT1GeV/AnalysisResults.root");
//  if (isMC) genTuner->SetDataFile("/Users/pillot/Work/Alice/Work/Data/2013/LHC13de/muon_pass2/AOD/GenTuner/pT2GeV/AnalysisResults.root");
  if (isMC) genTuner->SetDataFile("/Users/pillot/Work/Alice/Work/Data/2013/LHC13f/muon_calo/AOD127/GenTuner/pT1GeV/AnalysisResults.root");
//  if (isMC) genTuner->SetDataFile("/Users/pillot/Work/Alice/Work/Data/2013/LHC13f/muon_calo/AOD127/GenTuner/pT2GeV/AnalysisResults.root");
//  if (isMC) genTuner->SetDataFile("/Users/pillot/Work/Alice/Work/Data/2013/LHC13f/muon_calo/AOD127/GenTuner/pT4GeV/AnalysisResults.root");
//  if (isMC) genTuner->SetDataFile("/Users/pillot/Work/Alice/Work/Data/2013/LHC13f/muon_calo/AOD127/GenTuner/pT6GeV/AnalysisResults.root");
//  if (isMC) genTuner->SetDataFile("/Users/pillot/Work/Alice/Work/Data/2013/LHC13f/muon_calo/AOD127/GenTuner/pT1GeV_y2.5-3/AnalysisResults.root");
//  if (isMC) genTuner->SetDataFile("/Users/pillot/Work/Alice/Work/Data/2013/LHC13f/muon_calo/AOD127/GenTuner/pT1GeV_y3-4/AnalysisResults.root");
  
  if (iStep == 0) {
    
    // set the generator parameters used in simulation
    genTuner->SetPtParam(oldPtParam, oldFixPtParam, newPtParam, newFixPtParam, ptRange[0], ptRange[1]);
    genTuner->SetYParam(oldYParam, oldFixYParam, newYParam, newFixYParam, yRange[0], yRange[1]);
    
  } else if (iStep > 0) {
    /*
    // get the original generator parameters from first step if any
    TFile *inFile = TFile::Open("Results_step0.root","READ");
    if (inFile && inFile->IsOpen()) {
      TF1 *fOldPtFunc = static_cast<TF1*>(inFile->FindObjectAny("fPtFunc"));
      TF1 *fOldYFunc = static_cast<TF1*>(inFile->FindObjectAny("fYFunc"));
      if (fOldPtFunc && fOldYFunc) {
	fOldPtFunc->GetParameters(oldPtParam);
	fOldYFunc->GetParameters(oldYParam);
      }
      inFile->Close();
    }
    */
    // get the new generator parameters from previous step if any and configure the tuner
    TString inFileName = Form("Results_step%d.root",iStep-1);
    inFile = TFile::Open(inFileName.Data(),"READ");
    if (inFile && inFile->IsOpen()) {
      TF1 *fNewPtFunc = static_cast<TF1*>(inFile->FindObjectAny("fPtFuncNew"));
      TF1 *fNewYFunc = static_cast<TF1*>(inFile->FindObjectAny("fYFuncNew"));
      if (fNewPtFunc && fNewYFunc) {
	genTuner->SetPtParam(oldPtParam, newFixPtParam, fNewPtFunc->GetParameters(), newFixPtParam, fNewPtFunc->GetXmin(), fNewPtFunc->GetXmax());
	genTuner->SetYParam(oldYParam, newFixYParam, fNewYFunc->GetParameters(), newFixYParam, fNewYFunc->GetXmin(), fNewYFunc->GetXmax());
      }
      inFile->Close();
    }
    
    // enable the weighing
    genTuner->Weight(kTRUE);
    
  }
  
  return genTuner;
  
}

