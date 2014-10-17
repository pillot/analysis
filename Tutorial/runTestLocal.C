//
//  runTestLocal.C
//  aliroot_dev
//
//  Created by philippe pillot on 17/10/2014.
//  Copyright (c) 2014 Philippe Pillot. All rights reserved.
//

//______________________________________________________________________________
void runTestLocal(TString inputFileName = "AliAOD.Muons.root")
{
  /// run the test task locally
  
  // --- prepare environment ---
  gSystem->Load("libVMC");
  gSystem->Load("libTree.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libMinuit.so");
  gSystem->Load("libXMLParser.so");
  gSystem->Load("libGui.so");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gROOT->LoadMacro("AliAnalysisTaskTest.cxx++g");
  
  // --- Create the analysis train ---
  CreateAnalysisTrain();
  
  // --- Create input object ---
  TChain *chain = new TChain("aodTree");
  if(!chain->Add(inputFileName.Data(), -1)) {
    Error("runTestLocal","Problem creating the input chain");
    return;
  }
  
  // --- start analysis ---
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (mgr->InitAnalysis()) {
    mgr->PrintStatus();
    mgr->StartAnalysis("local", chain);
  }
  
}

//______________________________________________________________________________
void CreateAnalysisTrain()
{
  /// create the analysis train and configure it
  
  // analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("TestAnalysis");
  
  // AOD handler
  AliInputEventHandler* aodH = new AliAODInputHandler();
  mgr->SetInputEventHandler(aodH);
  
  // create the test task and configure it
  gROOT->LoadMacro("AddTaskTest.C");
  AliAnalysisTaskTest* task = AddTaskTest();
  if(!task) {
    Error("CreateAnalysisTrain","AliAnalysisTaskTest not created");
    return;
  }
  
}

