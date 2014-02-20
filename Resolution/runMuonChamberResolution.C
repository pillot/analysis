/*
 *  runESDCheck.C
 *  aliroot
 *
 *  Created by Philippe Pillot on 14/04/10.
 *  Copyright 2010 SUBATECH. All rights reserved.
 *
 */

void runMuonChamberResolution(TString runMode="test", TString runListName="runList.txt", Int_t extrapMode = 1,
			      Bool_t correctForSystematics = kTRUE, Double_t minMomentum = 0., Bool_t matchTrig = kFALSE)
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
  TString path ("/Users/philippe/Work/Alice/Work/Data/Macro/Resolution");
  TObjArray fileList(100);
  fileList.SetOwner();
  fileList.AddLast(new TObjString("runMuonChamberResolution.C"));
  fileList.AddLast(new TObjString("CreateAlienHandler.C"));
  fileList.AddLast(new TObjString("AliAnalysisTaskMuonChamberResolution.cxx"));
  fileList.AddLast(new TObjString("AliAnalysisTaskMuonChamberResolution.h"));
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
    if (overwrite == 'y' || overwrite == 'a') gSystem->Exec(Form("cp %s/%s %s", path.Data(), file->GetName(),file->GetName()));
    else if (overwrite == 'k') break;
  }
  
  // Create and configure the alien handler plugin
  gROOT->LoadMacro("CreateAlienHandler.C");
  AliAnalysisGrid *alienHandler = CreateAlienHandler(runMode.Data(), runListName.Data());
  if (!alienHandler) return;
  
  // Create the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("MuonChamberResolutionAnalysis");
  
  // Connect plugin to the analysis manager
  mgr->SetGridHandler(alienHandler);
  
  // ESD input handler
  AliESDInputHandler* esdH = new AliESDInputHandler();
  esdH->SetReadFriends(kFALSE);
  mgr->SetInputEventHandler(esdH);
  
  // Muon chamber resolution analysis
  gROOT->LoadMacro("AliAnalysisTaskMuonChamberResolution.cxx++g");
  AliAnalysisTaskMuonChamberResolution* task = new AliAnalysisTaskMuonChamberResolution("MuonChamberResolution");
  //Double_t clusterResNB[10] = {0.075, 0.060, 0.124, 0.120, 0.095, 0.103, 0.125, 0.140, 0.169, 0.180};
  //Double_t clusterResB[10] = {0.0573, 0.0430, 0.0938, 0.0919, 0.1451, 0.1133, 0.1773, 0.1899, 0.1156, 0.1487};
  //task->SetStartingResolution(clusterResNB, clusterResB);
  task->SetMinMomentum(minMomentum);
  task->MatchTrigger(matchTrig);
  task->SetExtrapMode(extrapMode);
  task->CorrectForSystematics(correctForSystematics);
  mgr->AddTask(task);
  
  // Create containers for output
  TString outFileName = "chamberResolution.root";
  AliAnalysisDataContainer *cout_histo1 = mgr->CreateContainer("Residuals", TObjArray::Class(), AliAnalysisManager::kOutputContainer, outFileName.Data());
  AliAnalysisDataContainer *cout_histo2 = mgr->CreateContainer("ResidualsVsP", TObjArray::Class(), AliAnalysisManager::kOutputContainer, outFileName.Data());
  AliAnalysisDataContainer *cout_histo3 = mgr->CreateContainer("Summary", TObjArray::Class(), AliAnalysisManager::kParamContainer, outFileName.Data());
  
  // Connect input/output
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, cout_histo1);
  mgr->ConnectOutput(task, 2, cout_histo2);
  mgr->ConnectOutput(task, 3, cout_histo3);
  
  // Enable debug printouts
  //mgr->SetDebugLevel(2);
  
  if (!mgr->InitAnalysis()) return;
  
  mgr->PrintStatus();
  // Start analysis in grid.
  mgr->StartAnalysis("grid");
  
  // save the summary canvases
  if (task->GetCanvases()) {
    TFile* outFile = TFile::Open("chamberResolution.root","UPDATE");
    if (outFile && outFile->IsOpen()) {
      task->GetCanvases()->Write();
      outFile->Close();
    }
  }
  
}

