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
  TString aliphysicsVersion = "vAN-20171030-1";
  TString extraLibs="";
  TString extraIncs="include";
  TString extraTasks="";
  TString extraPkgs="PWGmuon";
  TList pathList; pathList.SetOwner();
  pathList.Add(new TObjString("$WORK/Macros/MuonResolution"));
  pathList.Add(new TObjString("$WORK2/devphys/AliPhysics/PWG/muon"));
  pathList.Add(new TObjString("$ALICE_PHYSICS/PARfiles"));
  TList fileList; fileList.SetOwner();
  fileList.Add(new TObjString("runMuonTrackSmearingQA.C"));
  fileList.Add(new TObjString("AddTaskMuonTrackSmearingQA.C"));
  fileList.Add(new TObjString("PWGmuon.par"));

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
  trackCuts.SetFilterMask(/*AliMuonTrackCuts::kMuMatchLpt | */AliMuonTrackCuts::kMuEta |
			  AliMuonTrackCuts::kMuThetaAbs);
  trackCuts.SetIsMC(kTRUE);
  
  // muon track smearing task
  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/muon/AddTaskMuonTrackSmearing.C");
  AliTaskMuonTrackSmearing* muonSmearing = reinterpret_cast<AliTaskMuonTrackSmearing*>(gROOT->ProcessLineSync("AddTaskMuonTrackSmearing()"));
  if(!muonSmearing) {
    Error("CreateAnalysisTrain","AliTaskMuonTrackSmearing not created!");
    return;
  }
//  muonSmearing->GetMuonTrackSmearing().SetSigmaxChSt1(0.000519);
//  muonSmearing->GetMuonTrackSmearing().SetSigmayChSt1(0.000064);
//  muonSmearing->GetMuonTrackSmearing().SetSigmayCh(0.000105);
//  muonSmearing->GetMuonTrackSmearing().SetCrystalBallParams(3.671484,4.590817,2.453060,1.715218,2.145035,1.782320);

  // muon track smearing QA task
//  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/muon/AddTaskMuonTrackSmearingQA.C");
  gROOT->LoadMacro("AddTaskMuonTrackSmearingQA.C");
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

