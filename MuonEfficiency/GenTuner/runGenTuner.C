/*
 *  runJPsiAccEffCorr.C
 *  aliroot_dev
 *
 *  Created by philippe pillot on 06/03/13.
 *  Copyright 2011 SUBATECH. All rights reserved.
 *
 */


TString rootVersion = "v5-34-30-1";
TString alirootVersion = "";
TString aliphysicsVersion = "vAN-20150908";
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
*//*
// tune1 LHC13f
Double_t oldPtParam[6] = {455.614, 0.942071, 0.706755, 8.69856, 0.000168775, -0.925487};
Bool_t oldFixPtParam[6] = {kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE};
Double_t newPtParam[6] = {455.614, 0.942071, 0.706755, 8.69856, 0.000168775, -0.925487};
Bool_t newFixPtParam[6] = {kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE};
*/
// tune0 LHC15g
TString oldPtFormula = "[0]/TMath::Power([1]+TMath::Power(x,[2]),[3])";
Double_t oldPtParam[4] = {4.05962, 1.0, 2.46187, 2.08644};
Bool_t oldFixPtParam[4] = {kFALSE, kFALSE, kFALSE, kFALSE};
TString newPtFormula = "[0] * (1. / TMath::Power([1] + TMath::Power(x,[2]), [3]) + [4] * TMath::Exp([5]*x))";
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
*//*
// tune1 LHC13f
Double_t oldYParam[8] = {1.29511, 1., 0., -0.0767846, 0., 0.00176313, 0., 2.};
Bool_t oldFixYParam[8] = {kFALSE, kTRUE, kTRUE, kFALSE, kTRUE, kFALSE, kTRUE, kTRUE};
Double_t newYParam[8] = {1.29511, 1., 0., -0.0767846, 0., 0.00176313, 0., 2.};
Bool_t newFixYParam[8] = {kFALSE, kTRUE, kTRUE, kFALSE, kTRUE, kFALSE, kTRUE, kTRUE};
*/
// tune0 LHC15g
TString oldYFormula = "[0]*(1.+[1]*x+[2]*x*x+[3]*x*x*x)";
Double_t oldYParam[4] = {0.729545, 0.53837, 0.141776, 0.0130173};
Bool_t oldFixYParam[4] = {kFALSE, kFALSE, kFALSE, kFALSE};
//TString newYFormula = "[0] * ([1] * (1. + [2]*x + [3]*x*x + [4]*x*x*x + [5]*x*x*x*x) + [6]*TMath::Exp(-0.5*x*x/[7]/[7]))";
TString newYFormula = "[0] * (1. + [1]*x*x + [2]*x*x*x*x)";
Double_t newYParam[3] = {1.29511, -0.0767846, 0.00176313};
Bool_t newFixYParam[3] = {kFALSE, kFALSE, kFALSE};

Double_t yRange[2] = {-4.2, -2.3};


Bool_t isMC = kTRUE;
Bool_t applyPhysicsSelection = kFALSE;


void UpdateParametersAndRanges(Int_t iStep);


//______________________________________________________________________________
void runGenTuner(TString smode = "local", TString inputFileName = "AliAOD.root",
		 Int_t iStep = -1, char overwrite = '\0')
{
  /// Tune single muon kinematics distribution
  
  gROOT->LoadMacro("$HOME/Work/Alice/Macros/Facilities/runTaskFacilities.C");
  
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
  CopyInputFileLocally("/Users/pillot/Work/Alice/Data/2015/LHC15g/muon_calo_pass1/GenTuner/CMSL7/AnalysisResults.root", "ReferenceResults.root", overwrite);
  fileList.Add(new TObjString("ReferenceResults.root"));
  
  // --- saf3 case ---
  if (mode == kSAF3Connect) {
    
    // run on SAF3
    if (!RunAnalysisOnSAF3(fileList, aliphysicsVersion, inputFileName)) return;
    
    // draw the results locally
    outFile = TFile::Open(Form("Results_step%d.root", iStep),"READ");
    if (outFile && outFile->IsOpen()) {
      outFile->FindObjectAny("cRes")->Draw();
      outFile->FindObjectAny("cRat")->Draw();
      outFile->Close();
    }
    
    // do not try to re-run locally!
    return;
    
  }
  
  // --- prepare environment ---
  TString extraLibs="PWGmuon";
  TString extraIncs="include";
  TString extraTasks="AliAnalysisTaskGenTuner";
  TString extraPkgs="";
  LoadAlirootLocally(extraLibs, extraIncs, extraTasks, extraPkgs);
  AliAnalysisGrid *alienHandler = 0x0;
  if (mode == kProof || mode == kProofLite) LoadAlirootOnProof(smode, rootVersion, aliphysicsVersion, extraLibs, extraIncs, extraTasks, extraPkgs, kTRUE);
  else if (mode == kGrid || mode == kTerminate) {
    TString analysisMacroName = "GenTuner";
    alienHandler = static_cast<AliAnalysisGrid*>(CreateAlienHandler(smode, alirootVersion, aliphysicsVersion, inputFileName, dataDir, dataPattern, outDir, extraLibs, extraIncs, extraTasks, extraPkgs, analysisMacroName, runFormat, ttl, maxFilesPerJob, maxMergeFiles, maxMergeStages));
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
    if (genTuner->GetCurrentPtFunc()) genTuner->GetCurrentPtFunc()->Write(0x0, TObject::kOverwrite);
    if (genTuner->GetCurrentPtFuncMC()) genTuner->GetCurrentPtFuncMC()->Write(0x0, TObject::kOverwrite);
    if (genTuner->GetNewPtFunc()) genTuner->GetNewPtFunc()->Write(0x0, TObject::kOverwrite);
    if (genTuner->GetCurrentYFunc()) genTuner->GetCurrentYFunc()->Write(0x0, TObject::kOverwrite);
    if (genTuner->GetCurrentYFuncMC()) genTuner->GetCurrentYFuncMC()->Write(0x0, TObject::kOverwrite);
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
//  if (applyPhysicsSelection) genTuner->SelectCollisionCandidates(AliVEvent::kMUU7);
  if (applyPhysicsSelection) genTuner->SelectCollisionCandidates(AliVEvent::kMUS7);
//  if (applyPhysicsSelection) genTuner->SelectCollisionCandidates(AliVEvent::kMUSH7);
  //genTuner->SelectCentrality(0., 90.);
  genTuner->SetMuonTrackCuts(trackCuts);
  genTuner->SetMuonPtCut(1.);
  genTuner->SetMuonGenPtCut(0.8);
  
  if (isMC) {
    
    genTuner->SetDataFile("ReferenceResults.root");
    
    // update the parameters and the fitting ranges from the previous step if any
    UpdateParametersAndRanges(iStep);
    
    // set the original function and parameters used in simulation
    genTuner->SetOriginPtFunc(oldPtFormula.Data(), oldPtParam, oldFixPtParam, ptRange[0], ptRange[1]);
    genTuner->SetOriginYFunc(oldYFormula.Data(), oldYParam, oldFixYParam, yRange[0], yRange[1]);
    
    // set the new function and initial parameters
    genTuner->SetNewPtFunc(newPtFormula.Data(), newPtParam, newFixPtParam, ptRange[0], ptRange[1]);
    genTuner->SetNewYFunc(newYFormula.Data(), newYParam, newFixYParam, yRange[0], yRange[1]);
    
    // enable the weighing
    if (iStep > 0) genTuner->Weight(kTRUE);
    
  }
  
  return genTuner;
  
}

//______________________________________________________________________________
void UpdateParametersAndRanges(Int_t iStep)
{
  /// update the parameters and the fitting ranges from the previous step
  
  if (iStep <= 0) return;
  
  TString inFileName = Form("Results_step%d.root",iStep-1);
  inFile = TFile::Open(inFileName.Data(),"READ");
  if (!inFile || !inFile->IsOpen()) {
    printf("cannot open file from previous step\n");
    exit(1);
  }
  
  TF1 *fNewPtFunc = static_cast<TF1*>(inFile->FindObjectAny("fPtFuncNew"));
  TF1 *fNewYFunc = static_cast<TF1*>(inFile->FindObjectAny("fYFuncNew"));
  if (!fNewPtFunc || !fNewYFunc) {
    printf("previous functions not found\n");
    exit(1);
  }
  
  if ((fNewPtFunc->GetNpar() != (Int_t)(sizeof(newPtParam)/sizeof(Double_t))) ||
      (fNewYFunc->GetNpar() != (Int_t)(sizeof(newYParam)/sizeof(Double_t)))) {
    printf("mismatch between the number of parameters in the previous step and in this macro\n");
    exit(1);
  }
  
  for (Int_t i = 0; i < fNewPtFunc->GetNpar(); ++i) newPtParam[i] = fNewPtFunc->GetParameter(i);
  ptRange[0] = fNewPtFunc->GetXmin();
  ptRange[1] = fNewPtFunc->GetXmax();
  
  for (Int_t i = 0; i < fNewYFunc->GetNpar(); ++i) newYParam[i] = fNewYFunc->GetParameter(i);
  yRange[0] = fNewYFunc->GetXmin();
  yRange[1] = fNewYFunc->GetXmax();
  
  inFile->Close();
  
}

