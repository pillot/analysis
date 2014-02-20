/*
 *  runESDCheck.C
 *  aliroot
 *
 *  Created by Philippe Pillot on 14/04/10.
 *  Copyright 2010 SUBATECH. All rights reserved.
 *
 */

void runESDCheck(TString runMode="test", TString runListName="runList.txt",
		 Bool_t physicsSelect = kTRUE, Short_t selectCharge = 0)
{
  // Load common libraries
  gSystem->Load("libTree.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");    // The plugin is here
  gSystem->Load("libGui");
  gSystem->Load("libMinuit");
  gSystem->Load("libProofPlayer");
  gSystem->Load("libXMLParser");
  gSystem->Load("libRAWDatabase");
  gSystem->Load("libCDB");
  gSystem->Load("libSTEER");
  gSystem->Load("libMUONcore");
  gSystem->Load("libMUONmapping");
  gSystem->Load("libMUONcalib");
  gSystem->Load("libMUONgeometry");
  gSystem->Load("libMUONtrigger");
  gSystem->Load("libMUONraw");
  gSystem->Load("libMUONbase");
  gSystem->Load("libMUONrec");
  
  // Use AliRoot includes to compile our task
  gROOT->ProcessLine(".include $ALICE_ROOT/include");
  gROOT->ProcessLine(".include $ALICE_ROOT/MUON");
  gROOT->ProcessLine(".include $ALICE_ROOT/MUON/mapping");
  
  // copy files needed for this analysis
  TString path1 ("/Users/philippe/Work/Alice/Work/Data/Macro");
  TString path2 ("/Users/philippe/Work/Alice/Work/Data/Macro/ESDCheck");
  TObjArray fileList(100);
  fileList.SetOwner();
  fileList.AddLast(new TObjString("CreateAlienHandler.C"));
  fileList.AddLast(new TObjString("AliCounterCollection.cxx"));
  fileList.AddLast(new TObjString("AliCounterCollection.h"));
  fileList.AddLast(new TObjString("AliAnalysisTaskESDCheck.cxx"));
  fileList.AddLast(new TObjString("AliAnalysisTaskESDCheck.h"));
  fileList.AddLast(new TObjString("ConfigureCuts.C"));
  TIter nextFile(&fileList);
  TObjString *file;
  char overwrite = '';
  while ((file = static_cast<TObjString*>(nextFile()))) {
    if (overwrite != 'a') {
      overwrite = '';
      if (!gSystem->AccessPathName(file->GetName())) {
	cout<<Form("file %s exist in current directory. Overwrite? [y=yes, n=no, a=all, k=keep all] ",file->GetName())<<flush;
	while (overwrite != 'y' && overwrite != 'n' && overwrite != 'a' && overwrite != 'k') cin>>overwrite;
      } else overwrite = 'y';
    }
    if (overwrite == 'y' || overwrite == 'a') {
      if (!gSystem->AccessPathName(Form("%s/%s", path1.Data(), file->GetName())))
	gSystem->Exec(Form("cp %s/%s %s", path1.Data(), file->GetName(),file->GetName()));
      else
	gSystem->Exec(Form("cp %s/%s %s", path2.Data(), file->GetName(),file->GetName()));
    } else if (overwrite == 'k') break;
  }
  
  // Create and configure the alien handler plugin
  gROOT->LoadMacro("CreateAlienHandler.C");
  AliAnalysisGrid *alienHandler = CreateAlienHandler(runMode.Data(), runListName.Data());
  if (!alienHandler) return;
  
  // Create the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("ESDCheckAnalysis");
  
  // Connect plugin to the analysis manager
  mgr->SetGridHandler(alienHandler);
  
  // ESD input handler
  AliESDInputHandler* esdH = new AliESDInputHandler();
  esdH->SetReadFriends(kFALSE);
  mgr->SetInputEventHandler(esdH);
  
  // event selection
  if (physicsSelect) {
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
    AliPhysicsSelectionTask* physicsSelection = AddTaskPhysicsSelection();
    physicsSelection->GetPhysicsSelection()->SetUseMuonTriggers();
    //physicsSelection->GetPhysicsSelection()->SetAnalyzeMC();
  }
  
  // ESD scan analysis
  gROOT->LoadMacro("AliCounterCollection.cxx++g");
  gROOT->LoadMacro("AliAnalysisTaskESDCheck.cxx++g");
  AliAnalysisTaskESDCheck* task = new AliAnalysisTaskESDCheck("ESDCheck");
  task->SelectCharge(selectCharge);
  mgr->AddTask(task);
  
  // Create containers for output
  TString outFileName = "ESDCheck.root";
  AliAnalysisDataContainer *cout_histo1 = mgr->CreateContainer("general", TObjArray::Class(), AliAnalysisManager::kOutputContainer, outFileName.Data());
  AliAnalysisDataContainer *cout_histo2 = mgr->CreateContainer("expert", TObjArray::Class(), AliAnalysisManager::kOutputContainer, outFileName.Data());
  AliAnalysisDataContainer *cout_trackStat = mgr->CreateContainer("trackCounters", AliCounterCollection::Class(), AliAnalysisManager::kOutputContainer, outFileName.Data());
  AliAnalysisDataContainer *cout_eventStat = mgr->CreateContainer("eventCounters", AliCounterCollection::Class(), AliAnalysisManager::kOutputContainer, outFileName.Data());
//  AliAnalysisDataContainer *cout_test = mgr->CreateContainer("test", TH1F::Class(), AliAnalysisManager::kParamContainer, outFileName.Data());
  
  // Connect input/output
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, cout_histo1);
  mgr->ConnectOutput(task, 2, cout_histo2);
  mgr->ConnectOutput(task, 3, cout_trackStat);
  mgr->ConnectOutput(task, 4, cout_eventStat);
//  mgr->ConnectOutput(task, 5, cout_test);
  
  // Enable debug printouts
  //mgr->SetDebugLevel(2);
  
  if (!mgr->InitAnalysis()) return;
  
  mgr->PrintStatus();
  // Start analysis in grid.
  mgr->StartAnalysis("grid");
}

