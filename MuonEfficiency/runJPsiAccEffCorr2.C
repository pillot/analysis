/*
 *  runJPsiAccEffCorr.C
 *  aliroot_dev
 *
 *  Created by philippe pillot on 21/12/12.
 *  Copyright 2011 SUBATECH. All rights reserved.
 *
 */


TString runWeight = "runWeightCMUL7.txt";
//TString runWeight = "";

//______________________________________________________________________________
void runJPsiAccEffCorr2(TString smode = "local", TString inputFileName = "AliAOD.root",
                        Int_t ipT = -1, Int_t iy = -1,
                        Bool_t applyPhysicsSelection = kFALSE, Bool_t embedding = kFALSE)
{
  /// Compute the JPsi acc*eff correction
  
  // --- general analysis setup ---
  TString rootVersion = "";
  TString alirootVersion = "";
  TString aliphysicsVersion = "vAN-20160601-1";
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
  TString dataDir = "/alice/cern.ch/user/p/ppillot/Sim/LHC15n/JPsiTune1/VtxShift/results";
  TString dataPattern = "*AliAOD.Muons.root";
  TString runFormat = "%d";
  TString outDir = "Sim/LHC15n/JPsiTune1/VtxShift/AccEff";
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
  if (!runWeight.IsNull()) fileList.Add(new TObjString(runWeight.Data()));
  
  // --- run the analysis (saf3 is a special case as the analysis is launched on the server) ---
  if (mode == kSAF3Connect) {
    
    RunAnalysisOnSAF3(fileList, aliphysicsVersion, inputFileName, splitDataset);
    
  } else {
    
    // get data type
    TString dataType = GetDataType(mode, inputFileName, dataPattern);
    CreateAnalysisTrain(dataType, ipT, iy, applyPhysicsSelection, embedding);
    
    if (smode == "saf3" && splitDataset) AliAnalysisManager::GetAnalysisManager()->SetSkipTerminate(kTRUE);
    
    RunAnalysis(smode, inputFileName, rootVersion, alirootVersion, aliphysicsVersion, extraLibs, extraIncs, extraTasks, extraPkgs, dataDir, dataPattern, outDir, analysisMacroName, runFormat, ttl, maxFilesPerJob, maxMergeFiles, maxMergeStages);
    
  }
  
}

//______________________________________________________________________________
void CreateAnalysisTrain(TString dataType, Int_t ipT, Int_t iy, Bool_t applyPhysicsSelection, Bool_t embedding)
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
  
  // track selection
  AliMuonTrackCuts trackCuts("stdCuts", "stdCuts");
  trackCuts.SetAllowDefaultParams();
  //  trackCuts.SetFilterMask(0);
  trackCuts.SetFilterMask(AliMuonTrackCuts::kMuMatchLpt | AliMuonTrackCuts::kMuEta |
                          AliMuonTrackCuts::kMuThetaAbs);
  
  // Acc*Eff results
  gROOT->LoadMacro("AddTaskJPsiAccEffCorr2.C");
  AliAnalysisTaskJPsiAccEffCorr2 *jPsiAccEffCorr = AddTaskJPsiAccEffCorr2();
  if(!jPsiAccEffCorr) {
    Error("CreateAnalysisTrain","AliAnalysisTaskJPsiAccEffCorr2 not created!");
    return;
  }
  if (applyPhysicsSelection) jPsiAccEffCorr->SelectCollisionCandidates(offlineTriggerMask);
  jPsiAccEffCorr->SetMuonTrackCuts(trackCuts);
  jPsiAccEffCorr->SetNMatch(2);
  jPsiAccEffCorr->UseMCLabel(kTRUE);
//  jPsiAccEffCorr->SetMuLowPtCut(1.);
//  Float_t centBinLowEdge[] = {0.,10.,20.,30.,40.,50.,60.,70.,80.,90.};
  Float_t centBinLowEdge[] = {-9999.,9999.};
  jPsiAccEffCorr->SetCentBins((Int_t)(sizeof(centBinLowEdge)/sizeof(Float_t))-1, centBinLowEdge);
//  Float_t pTBinLowEdge[] = {0.,1.,2.,3.,4.,5.,6.,8.};
  Float_t pTBinLowEdge[] = {0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,12.};
  jPsiAccEffCorr->SetPtBins((Int_t)(sizeof(pTBinLowEdge)/sizeof(Float_t))-1, pTBinLowEdge);
//  Float_t yBinLowEdge[] = {-4.,-3.5,-3.,-2.5};
  Float_t yBinLowEdge[] = {-4.,-3.75,-3.5,-3.25,-3.,-2.75,-2.5};
  jPsiAccEffCorr->SetYBins((Int_t)(sizeof(yBinLowEdge)/sizeof(Float_t))-1, yBinLowEdge);
  if (!runWeight.IsNull()) jPsiAccEffCorr->LoadRunWeights(runWeight.Data());
  SetGenWeights(jPsiAccEffCorr, ipT, iy);
//  SetSigWeights(jPsiAccEffCorr);
  
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

//______________________________________________________________________________
void SetGenWeights(TObject* jPsiAccEffCorr, Int_t ipT, Int_t iy)
{
  /// set the old and new pT/y functions to reweight the generated JPsi
  
  AliAnalysisTaskJPsiAccEffCorr2* accEffCorr = static_cast<AliAnalysisTaskJPsiAccEffCorr2*>(jPsiAccEffCorr);
  
  Double_t ptRange[2] = {0., 999.};
  Double_t yRange[2] = {-4.2, -2.3};
  
  if (ipT < 0 && iy < 0) {
    
    TString oldPtFormula = "[0] * x / TMath::Power([1] + TMath::Power(x,[2]), [3])";
    Double_t oldPtParam[4] = {4654.3, 12.8133, 1.9647, 3.66641};
    accEffCorr->SetOriginPtFunc(oldPtFormula.Data(), oldPtParam, ptRange[0], ptRange[1]);
    
    TString newPtFormula = "[0] * x / TMath::Power([1] + TMath::Power(x,[2]), [3])";
//    Double_t newPtParam[4] = {4654.3, 12.8133, 1.9647, 3.66641};
    Double_t newPtParam[4]  = {5.21831e+04, 1.46939e+01,1.93309, 3.93941}; // fit cross section
//    Double_t newPtParam[4] = {5164.53, 13.868, 1.98408, 3.6322}; // from genTuner
    // AliGenMUONlib "pp E" tuning, where E is the energy set in param[1] (in GeV)
    //  TString newPtFormula = "[0] * x / TMath::Power(1.+0.363*TMath::Power((x/(1.04*TMath::Power([1],0.101))),2.),3.9)";
    //  Double_t newPtParam[2] = {1., 5030.};
    accEffCorr->SetNewPtFunc(newPtFormula.Data(), newPtParam, ptRange[0], ptRange[1]);
    
    TString oldYFormula = "[0] * (1. + [1]*x*x)";
    Double_t oldYParam[2] = {1.18296, -0.0405994};
    accEffCorr->SetOriginYFunc(oldYFormula.Data(), oldYParam, yRange[0], yRange[1]);
    
    TString newYFormula = "[0] * (1. + [1]*x*x)";
//    Double_t newYParam[2] = {1.18296, -0.0405994};
    Double_t newYParam[2] = {6.36959, -3.99165e-02}; // fit cross section
//    Double_t newYParam[2] = {1.17479, -0.040254}; // from genTuner
    // AliGenMUONlib "pp E" tuning, where E is the energy set in param[1] (in GeV)
    //  TString newYFormula = "[0] * TMath::Exp(-0.5*TMath::Power(x/TMath::Log([1]/3.097)/0.4,2.))";
    //  Double_t newYParam[2] = {1., 5030.};
    accEffCorr->SetNewYFunc(newYFormula.Data(), newYParam, yRange[0], yRange[1]);
    
//  } else if (ipT >= 0 && ipT < 7 && iy >= 0 && iy < 13) {
//  } else if (ipT >= 0 && ipT < 6 && iy >= 0 && iy < 13) {
  } else if (ipT >= 0 && ipT < 4 && iy >= 0 && iy < 12) {
    /*
    // pp 13 TeV
    TString ptFormula = "[0] * x / TMath::Power(1. + TMath::Power(x/[1],[2]), [3])";
    const Double_t ptParam[7][4] = {
      {1., 4.75208, 1.69247, 4.49224},
      {1., 4.608, 1.6918, 4.2423},
      {1., 5.0297, 1.6365, 4.7689},
      {1., 4.4274, 1.7741, 4.026},
      {1., 4.7547, 1.6953, 4.4997},
      {1., 4.5144, 1.7925, 4.2223},
      {1., 5.448, 1.5637, 5.5774}};
    *//*
    // pp 7 TeV (LHCb: https://arxiv.org/pdf/1103.0423v2.pdf ) my fits
    TString ptFormula = "[0] * x / TMath::Power([1] + TMath::Power(x,[2]), [3])";
    const Double_t ptParam[6][4] = {
      {135901421.678099, 13.976580, 1.885244, 3.861901},
      {34559125.841439, 15.561190, 1.944794, 3.586307},
      {49631156.645487, 14.642333, 1.893818, 3.793124},
      {28053036.112401, 13.785175, 1.926424, 3.688917},
      {100872181.540582, 13.718217, 1.830482, 4.247041},
      {51184062.999060, 14.769584, 1.953139, 4.019600}};
    *//*
    // pp 7 TeV (LHCb: https://arxiv.org/pdf/1103.0423v2.pdf ) my fits 2
    TString ptFormula = "[0] * x / TMath::Power(1. + TMath::Power(x/[1],[2]), [3])";
    const Double_t ptParam[6][4] = {
      {5607.654784, 4.522230, 1.727098, 4.589162},
      {2637.774943, 5.287302, 1.525239, 5.544391},
      {2074.243503, 4.586102, 1.732606, 4.498168},
      {1956.414017, 4.490100, 1.722058, 4.595126},
      {1582.221713, 4.655218, 1.706322, 4.949379},
      {1173.498991, 5.108379, 1.663448, 5.772206}};
    */
    // pp 7 TeV (LHCb: https://arxiv.org/pdf/1103.0423v2.pdf ) Jana's fits
    TString ptFormula = "[0] * x / TMath::Power((1 + [1]*x*x), [2])";
    const Double_t ptParam[4][3] = {
      {1720.43, 0.07818, 3.27866},
      {1946.4, 0.0768768, 3.21477},
      {1773.8, 0.0777934, 3.30499},
      {1432.46, 0.076854, 3.43851}};
    /*
    TString oldPtFormula = "[0] * x / TMath::Power([1] + TMath::Power(x,[2]), [3])";
    Double_t oldPtParam[4] = {4654.3, 12.8133, 1.9647, 3.66641};
    accEffCorr->SetOriginPtFunc(oldPtFormula.Data(), oldPtParam, ptRange[0], ptRange[1]);
    */
    accEffCorr->SetOriginPtFunc(ptFormula.Data(), ptParam[0], ptRange[0], ptRange[1]);
    accEffCorr->SetNewPtFunc(ptFormula.Data(), ptParam[ipT], ptRange[0], ptRange[1]);
    /*
    // pp 13 TeV
    TString yFormula = "[0] * TMath::Exp(-0.5*x*x/[1]/[1])";
    const Double_t yParam[13][2] = {
      {1., 2.98887},
      {1.6018, 3.5004},
      {3.0763, 3.4912},
      {2.7625, 3.2579},
      {1.832, 3.1353},
      {1.1356, 2.9421},
      {0.78027, 2.5091},
      {0.38906, 2.7312},
      {0.27124, 2.3781},
      {0.10791, 2.5896},
      {0.050845, 2.2532},
      {0.026694, 1.815},
      {0.003918, 2.249}};
    *//*
    // pp 7 TeV (LHCb: https://arxiv.org/pdf/1103.0423v2.pdf ) my fits
    TString yFormula = "[0] * (1. + [1]*x*x)";
    const Double_t yParam[13][2] = {
      {8043.370402, -0.037174},
      {1302.471164, -0.035378},
      {2288.167393, -0.034248},
      {1881.161207, -0.036614},
      {1152.417371, -0.038817},
      {633.980533, -0.039936},
      {346.591513, -0.042033},
      {182.095618, -0.043719},
      {101.616307, -0.045046},
      {58.948019, -0.046820},
      {33.824028, -0.048227},
      {20.465781, -0.048441},
      {13.876576, -0.054479}};
    *//*
    // pp 7 TeV (LHCb: https://arxiv.org/pdf/1103.0423v2.pdf ) my fits 2
    TString yFormula = "[0] * TMath::Exp(-0.5*(x-[1])*(x-[1])/[2]/[2])";
    const Double_t yParam[13][3] = {
      {9787.590971, 0.000000, 2.739344},
      {1462.099491, 0.000000, 2.961196},
      {2510.616034, 0.000000, 3.094858},
      {2175.536695, 0.000000, 2.862361},
      {1403.080318, 0.000000, 2.675295},
      {817.962662, 0.000000, 2.538551},
      {482.813895, 0.000000, 2.364329},
      {263.394556, 0.000000, 2.277106},
      {163.866381, 0.000000, 2.140196},
      {90.951878, 0.000000, 2.119798},
      {58.197163, 0.000000, 2.005567},
      {31.018428, 0.000000, 2.102733},
      {23.571972, 0.000000, 1.898878}};
    */
    // pp 7 TeV (LHCb: https://arxiv.org/pdf/1103.0423v2.pdf ) Jana's fits
    TString yFormula = "[0] * TMath::Exp(-0.5*(x-[1])*(x-[1])/[2]/[2])";
    const Double_t yParam[12][3] = {
      {6189.19, -2.15031, 1.61276},
      {2394.34, 1.78261, 3.37687},
      {1679.64, -2.51158, 1.49752},
      {1392.31, -2.3965, 1.468},
      {884.81, -2.18453, 1.5023},
      {495.167, -2.10634, 1.49666},
      {278.361, -1.98704, 1.46456},
      {144.209, -2.0503, 1.36351},
      {81.6065, -2.01818, 1.32101},
      {49.4501, -1.83735, 1.356},
      {26.2173, -2.11165, 1.16468},
      {16.078, -2.06902, 1.17026}};
    /*
    TString oldYFormula = "[0] * (1. + [1]*x*x)";
    Double_t oldYParam[2] = {1.18296, -0.0405994};
    accEffCorr->SetOriginYFunc(oldYFormula.Data(), oldYParam, yRange[0], yRange[1]);
    */
    accEffCorr->SetOriginYFunc(yFormula.Data(), yParam[0], yRange[0], yRange[1]);
    accEffCorr->SetNewYFunc(yFormula.Data(), yParam[iy], yRange[0], yRange[1]);
    
  }
  
}

