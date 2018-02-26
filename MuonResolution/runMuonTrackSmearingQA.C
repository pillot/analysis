/*
 *  runMuonTrackSmearingQA.C
 *
 *  Created by Philippe Pillot on 08/11/17.
 *  Copyright 2017 SUBATECH
 *  This software is made available under the terms of the GNU GPL 3.0
 *
 */

#if !defined(__runMuonTrackSmearingQA__)
#define __runMuonTrackSmearingQA__

#if __has_include("/Users/pillot/Work/Alice/Macros/Facilities/runTaskFacilities.C")
#include "/Users/pillot/Work/Alice/Macros/Facilities/runTaskFacilities.C"
#else
#include "runTaskFacilities.C"
#endif

//______________________________________________________________________________
void runMuonTrackSmearingQA(TString smode = "local", TString inputFileName = "AliAOD.Muons.root")
{
  /// Extract some physical quantities
  
  // --- general analysis setup ---
  TString rootVersion = "";
  TString alirootVersion = "";
  TString aliphysicsVersion = "vAN-20180105-1";
  TString extraLibs="";
  TString extraIncs="include";
  TString extraTasks="";
  TString extraPkgs="";
//  TString extraPkgs="PWGmuon";
  TList pathList; pathList.SetOwner();
  pathList.Add(new TObjString("$WORK/Macros/MuonResolution"));
//  pathList.Add(new TObjString("$WORK2/devphys/AliPhysics/PWG/muon"));
  //pathList.Add(new TObjString("$ALICE_PHYSICS/PARfiles"));
  TList fileList; fileList.SetOwner();
  fileList.Add(new TObjString("runMuonTrackSmearingQA.C"));
//  fileList.Add(new TObjString("AddTaskMuonTrackSmearingQA.C"));
  //fileList.Add(new TObjString("PWGmuon.par"));

  // --- grid specific setup ---
  TString dataDir = "/alice/data/2017/LHC17n";
  TString dataPattern = "muon_calo_pass1/AOD/*AliAOD.Muons.root";
  TString runFormat = "%09d";
  TString outDir = "Data/LHC17n/muon_calo_pass1/AODAlignv2/PhysNoPS";
  TString analysisMacroName = "MuonSmearingQA";
  Int_t ttl = 30000;
  Int_t maxFilesPerJob = 20;
  Int_t maxMergeFiles = 10;
  Int_t maxMergeStages = 2;
  
  // --- saf3 specific setup ---
  Bool_t splitDataset = kFALSE;
  
  // --- prepare the analysis environment ---
  Int_t mode = PrepareAnalysis(smode, inputFileName, extraLibs, extraIncs, extraTasks, extraPkgs, pathList, fileList);
  
  // --- run the analysis (saf3 is a special case as the analysis is launched on the server) ---
  if (mode == kSAF3Connect) {
    
    RunAnalysisOnSAF3(fileList, aliphysicsVersion, inputFileName, splitDataset);
    
    // draw the results locally
    TFile *outFile = TFile::Open("AnalysisResults.root","READ");
    if (outFile && outFile->IsOpen()) {
      outFile->FindObjectAny("cGen")->Draw();
      outFile->FindObjectAny("cRec")->Draw();
      outFile->FindObjectAny("cRat")->Draw();
      outFile->FindObjectAny("cResVsP")->Draw();
      outFile->Close();
    }
    
  } else {
    
    gSystem->Exec(TString::Format("cp %s __runTask__.C", __FILE__));
    gROOT->LoadMacro("__runTask__.C");
    gSystem->Exec("rm __runTask__.C");
    TString dataType = GetDataType(mode, inputFileName, dataPattern);
    TObject *muonSmearingQA = reinterpret_cast<TObject*>(gROOT->ProcessLineSync(Form("CreateAnalysisTrain(\"%s\")",dataType.Data())));
    if (!muonSmearingQA) return;

    Bool_t terminate = kTRUE;
    if ((smode == "saf3" && splitDataset) || (mode == kGrid && smode != "terminate")) {
      AliAnalysisManager::GetAnalysisManager()->SetSkipTerminate(kTRUE);
      terminate = kFALSE;
    }
    
    RunAnalysis(smode, inputFileName, rootVersion, alirootVersion, aliphysicsVersion, extraLibs, extraIncs, extraTasks, extraPkgs, dataDir, dataPattern, outDir, analysisMacroName, runFormat, ttl, maxFilesPerJob, maxMergeFiles, maxMergeStages);
    
    if (terminate) gROOT->ProcessLineSync(TString::Format("Terminate(reinterpret_cast<TObject*>(%p))",muonSmearingQA));
    
  }
  
}

#else

//______________________________________________________________________________
TObject* CreateAnalysisTrain(TString dataType)
{
  /// create the analysis train and configure it

  // analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("MuonSmearingAnalysis");
  
  // Debug mode
  //mgr->SetDebugLevel(3);
  
  // ESD or AOD handler
  if (dataType == "ESD") {
    AliESDInputHandler* esdH = new AliESDInputHandler();
    esdH->SetReadFriends(kFALSE);
    mgr->SetInputEventHandler(esdH);
    AliMCEventHandler* mcH = new AliMCEventHandler();
    mgr->SetMCtruthEventHandler(mcH);
  } else if (dataType == "AOD") {
    AliInputEventHandler* aodH = new AliAODInputHandler();
    mgr->SetInputEventHandler(aodH);
  } else {
    Error("CreateAnalysisTrain","Unknown data type. Cannot define input handler!");
    return;
  }
  
  // track selection
  AliMuonTrackCuts trackCuts("stdCuts", "stdCuts");
  trackCuts.SetAllowDefaultParams();
  trackCuts.SetFilterMask(AliMuonTrackCuts::kMuMatchHpt | AliMuonTrackCuts::kMuEta |
			  AliMuonTrackCuts::kMuThetaAbs);
  trackCuts.SetIsMC(kTRUE);
  
  // muon track smearing task
  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/muon/AddTaskMuonTrackSmearing.C");
  AliTaskMuonTrackSmearing* muonSmearing = reinterpret_cast<AliTaskMuonTrackSmearing*>(gROOT->ProcessLineSync(Form("AddTaskMuonTrackSmearing(%d)",AliMuonTrackSmearing::kCrystalBall)));
  if(!muonSmearing) {
    Error("CreateAnalysisTrain","AliTaskMuonTrackSmearing not created!");
    return;
  }
//  muonSmearing->GetMuonTrackSmearing().SetSigmaTrk(0.004);
  muonSmearing->GetMuonTrackSmearing().SetSigmaTrkCut(6.);
//  muonSmearing->GetMuonTrackSmearing().SetNSigmaShift(1.2);
  // CB Data2
//  muonSmearing->GetMuonTrackSmearing().SetCrystalBallParams(2.017293,1.890778,1.514588,1.938707,1.331549,1.941949);
  // CB Data3 /Users/pillot/Work/Alice/Data/2015/LHC15n/muon_calo_pass1/244540/MuonResolution/AlignHugoV4_20GeV/results.root
  muonSmearing->GetMuonTrackSmearing().SetSigmaxChSt1(0.000430);
  muonSmearing->GetMuonTrackSmearing().SetSigmayChSt1(0.000095);
  muonSmearing->GetMuonTrackSmearing().SetSigmayCh(0.000125);
  muonSmearing->GetMuonTrackSmearing().SetCrystalBallParams(1.907786,1.846464,1.587050,1.750238,1.331549,1.941949);
  // gauss MC
//  muonSmearing->GetMuonTrackSmearing().SetSigmaxChSt1(0.000506);
//  muonSmearing->GetMuonTrackSmearing().SetSigmayChSt1(0.000065);
//  muonSmearing->GetMuonTrackSmearing().SetSigmayCh(0.000112);
  // CB MC
//  muonSmearing->GetMuonTrackSmearing().SetSigmaxChSt1(0.000519);
//  muonSmearing->GetMuonTrackSmearing().SetSigmayChSt1(0.000064);
//  muonSmearing->GetMuonTrackSmearing().SetSigmayCh(0.000105);
//  muonSmearing->GetMuonTrackSmearing().SetCrystalBallParams(3.671484,4.590817,2.453060,1.715218,2.145035,1.782320);
  // /Users/pillot/Work/Alice/Sim/LHC16r/muon_calo_pass1/mu/TuneCMSH1/pT5GeV/MuonEfficiency/AnalysisResults.root
//  muonSmearing->GetMuonTrackSmearing().SetSigmayCh(0.000125);
//  muonSmearing->GetMuonTrackSmearing().SetCrystalBallParams(3.671484,4.590817,2.453060,1.715218,2.498486,1.364013);
  // CB MC Station 3 /Users/pillot/Work/Alice/Sim/LHC15o/muon_calo_pass1/mu/TuneCMSH7/MuonEfficiency/AnalysisResults.root
//  muonSmearing->GetMuonTrackSmearing().SetSigmaxChSt1(0.000512);
//  muonSmearing->GetMuonTrackSmearing().SetSigmayChSt1(0.000065);
//  muonSmearing->GetMuonTrackSmearing().SetSigmayCh(0.000126);
//  muonSmearing->GetMuonTrackSmearing().SetCrystalBallParams(3.13205,0.813714,2.554298,1.479729,2.356118,1.695968);
  // CB Data Station 3 /Users/pillot/Work/Alice/Data/2015/LHC15n/muon_calo_pass1/244540/MuonResolution/AlignHugoV4_20GeV/results.root
//  muonSmearing->GetMuonTrackSmearing().SetSigmaxChSt1(0.000430);
//  muonSmearing->GetMuonTrackSmearing().SetSigmayChSt1(0.000095);
//  muonSmearing->GetMuonTrackSmearing().SetSigmayCh(0.000195);
//  muonSmearing->GetMuonTrackSmearing().SetCrystalBallParams(1.907786,1.846464,1.587050,1.750238,1.421882,2.139755);
//  muonSmearing->GetMuonTrackSmearing().SetCrystalBallParams(1.907786,1.846464,1.587050,1.750238,1.221338,2.599482); // step0
  //  CB MC measured /Users/pillot/Work/Alice/Sim/LHC15o/muon_calo_pass1/EmbedV2b/MuonResolution/20GeV/chamberResolution_step2.root
//  muonSmearing->GetMuonTrackSmearing().SetSigmaxChSt1(0.000656);
//  muonSmearing->GetMuonTrackSmearing().SetSigmayChSt1(0.000082);
//  muonSmearing->GetMuonTrackSmearing().SetSigmayCh(0.000122);
//  muonSmearing->GetMuonTrackSmearing().SetCrystalBallParams(2.817940,1.475378,1.551441,2.679521,1.756194,2.128461);
  //  CB MC measured 0-10% /Users/pillot/Work/Alice/Sim/LHC15o/muon_calo_pass1/EmbedV2b/MuonResolution/20GeV/chamberResolution_step2.root
//  muonSmearing->GetMuonTrackSmearing().SetSigmaxChSt1(0.000671);
//  muonSmearing->GetMuonTrackSmearing().SetSigmayChSt1(0.000085);
//  muonSmearing->GetMuonTrackSmearing().SetSigmayCh(0.000124);
//  muonSmearing->GetMuonTrackSmearing().SetCrystalBallParams(2.635614,1.436419,1.563039,2.210707,1.690065,1.936893);
  //  CB MC measured 80-90% /Users/pillot/Work/Alice/Sim/LHC15o/muon_calo_pass1/EmbedV2b/MuonResolution/20GeV/chamberResolution_step2.root
//  muonSmearing->GetMuonTrackSmearing().SetSigmaxChSt1(0.000649);
//  muonSmearing->GetMuonTrackSmearing().SetSigmayChSt1(0.000080);
//  muonSmearing->GetMuonTrackSmearing().SetSigmayCh(0.000121);
//  muonSmearing->GetMuonTrackSmearing().SetCrystalBallParams(2.907202,1.564498,1.477585,3.277859,1.772329,2.276131);
  //  CB MC measured Station 3 /Users/pillot/Work/Alice/Sim/LHC15o/muon_calo_pass1/EmbedV2b/MuonResolution/20GeV/chamberResolution_step2.root
//  muonSmearing->GetMuonTrackSmearing().SetSigmaxChSt1(0.000656);
//  muonSmearing->GetMuonTrackSmearing().SetSigmayChSt1(0.000082);
//  muonSmearing->GetMuonTrackSmearing().SetSigmayCh(0.000158);
//  muonSmearing->GetMuonTrackSmearing().SetCrystalBallParams(2.817940,1.475378,1.551441,2.679521,1.923223,1.951717);
  //  CB MC measured Station 3 /Users/pillot/Work/Alice/Sim/LHC15o/muon_calo_pass1/EmbedV2b/MuonResolution/20GeV/chamberResolution_step0.root
//  muonSmearing->GetMuonTrackSmearing().SetSigmaxChSt1(0.000695);
//  muonSmearing->GetMuonTrackSmearing().SetSigmayChSt1(0.000109);
//  muonSmearing->GetMuonTrackSmearing().SetSigmayCh(0.000162);
//  muonSmearing->GetMuonTrackSmearing().SetCrystalBallParams(2.74320,1.61658,1.532941,2.487253,1.773199,2.291951);

  // muon track smearing QA task
  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/muon/AddTaskMuonTrackSmearingQA.C");
//  gROOT->LoadMacro("AddTaskMuonTrackSmearingQA.C");
  AliTaskMuonTrackSmearingQA* muonSmearingQA = reinterpret_cast<AliTaskMuonTrackSmearingQA*>(gROOT->ProcessLineSync("AddTaskMuonTrackSmearingQA()"));
  if(!muonSmearingQA) {
    Error("CreateAnalysisTrain","AliTaskMuonTrackSmearingQA not created!");
    return;
  }
  muonSmearingQA->SetMuonTrackCuts(trackCuts);
  
  return muonSmearingQA;

}

//______________________________________________________________________________
void Terminate(TObject *o)
{
  /// save canvases
  
  AliTaskMuonTrackSmearingQA *muonSmearingQA = reinterpret_cast<AliTaskMuonTrackSmearingQA*>(o);
  
  TString outFileName = AliAnalysisManager::GetCommonFileName();
  TFile *outFile = (TFile*)gROOT->GetListOfFiles()->FindObject(outFileName.Data());
  if (outFile) outFile->ReOpen("UPDATE");
  else outFile = TFile::Open(outFileName.Data(),"UPDATE");
  if (outFile && outFile->IsOpen()) {
    outFile->Cd(Form("%s:/MUON_TrackSmearingQA",outFileName.Data()));
    if (muonSmearingQA->GetGen()) muonSmearingQA->GetGen()->Write(0x0, TObject::kOverwrite);
    if (muonSmearingQA->GetRec()) muonSmearingQA->GetRec()->Write(0x0, TObject::kOverwrite);
    if (muonSmearingQA->GetRat()) muonSmearingQA->GetRat()->Write(0x0, TObject::kOverwrite);
    if (muonSmearingQA->GetResVsP()) muonSmearingQA->GetResVsP()->Write(0x0, TObject::kOverwrite);
    outFile->Close();
  }
  
}

#endif

