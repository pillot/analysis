/*
 *  runMuonRefitFakes.C
 *  aliroot
 *
 *  Created by Philippe Pillot on 08/11/10.
 *  Copyright 2010 SUBATECH. All rights reserved.
 *
 */

enum {kLocal, kProof, kGrid};

//______________________________________________________________________________
void runMuonRefitFakes(TString smode = "local", TString inputFileName = "AliESDs.root",
		       TString alirootVersion = "VO_ALICE@AliRoot::v4-21-03-AN", Bool_t useMCLabels = kFALSE,
		       Double_t resNB = 0.2, Double_t resB = 0.2, Double_t sigmaCut = 4., Double_t sigmaCutTrig = 4.)
{
  /// Refit the tracks with new recoParam and/or new alignment
  /// check the effect on track quality, efficiency and physical quantities
  /// and how it affects the fake tracks
  
  // --- Check runing mode ---
  Int_t mode = GetMode(smode, inputFileName);
  if(mode < 0){
    Error("runMuonFakes","Please provide either an ESD root file a collection of ESDs or a dataset.");
    return;
  }
  
  // --- copy files needed for this analysis ---
  CopyFileLocally();
  
  // --- prepare environment ---
  TString extraLibs="RAWDatabase:RAWDatarec:CDB:STEER:MUONcore:MUONmapping:MUONcalib:MUONgeometry:MUONtrigger:MUONraw:MUONbase:MUONshuttle:MUONrec:RAWDatasim:MUONsim:MUONevaluation:PWG3base";
  LoadAlirootLocally(extraLibs);
  if (mode == kProof) LoadAlirootOnProof(smode, alirootVersion, extraLibs);
  
  // --- Create and configure the alien handler plugin if needed ---
  AliAnalysisGrid *alienHandler = 0x0;
  if (mode == kGrid) {
    gROOT->LoadMacro("CreateAlienHandler.C");
    alienHandler = CreateAlienHandler(smode.Data(), inputFileName.Data());
    if (!alienHandler) return;
  }
  
  // --- Create the analysis train ---
  CreateAnalysisTrain(useMCLabels, resNB, resB, sigmaCut, sigmaCutTrig, static_cast<Long_t>(alienHandler));
  
  // --- Create input object ---
  TObject* inputObj = 0x0;
  if (mode == kProof) inputObj = new TObjString(inputFileName);
  else if (mode != kGrid) inputObj = CreateChainFromFile(inputFileName.Data());
  
  // --- start analysis ---
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (mgr->InitAnalysis()) {
    mgr->PrintStatus();
    if (mode == kGrid) mgr->StartAnalysis("grid");
    else if (mode == kProof && inputObj) mgr->StartAnalysis("proof", Form("%s",static_cast<TObjString*>(inputObj)->GetName()));
    else if (inputObj) {
      mgr->SetUseProgressBar(kTRUE);
      mgr->StartAnalysis("local", static_cast<TChain*>(inputObj));
    }
  }
  
  // --- save summary canvases ---
  AliAnalysisTaskMuonFakes *fakes = static_cast<AliAnalysisTaskMuonFakes*>(mgr->GetTasks()->FindObject("MUONFakes"));
  if (fakes && fakes->GetCanvases()) {
    TFile* outFile = TFile::Open(AliAnalysisManager::GetCommonFileName(),"UPDATE");
    if (outFile && outFile->IsOpen()) {
      outFile->Cd("MUON_Fakes");
      fakes->GetCanvases()->Write();
      outFile->Close();
    }
  }
  
}

//______________________________________________________________________________
void CopyFileLocally()
{
  /// Copy files needed for this analysis
  
  TString path1("/Users/philippe/Work/Alice/Work/Data/Macro/Fakes");
  TString path2("/Users/philippe/Work/Alice/Work/Data/Macro/MuonPhysics");
  TString path3("$ALICE_ROOT/PWG3/muondep");
  TObjArray fileList(100);
  fileList.SetOwner();
  fileList.AddLast(new TObjString("runMuonRefitFakes.C"));
  fileList.AddLast(new TObjString("CreateAlienHandler.C"));
  fileList.AddLast(new TObjString("AliAnalysisTaskMuonPhysics.cxx"));
  fileList.AddLast(new TObjString("AliAnalysisTaskMuonPhysics.h"));
  fileList.AddLast(new TObjString("AliAnalysisTaskMuonRefit.cxx"));
  fileList.AddLast(new TObjString("AliAnalysisTaskMuonRefit.h"));
  fileList.AddLast(new TObjString("AliAnalysisTaskMuonFakes.cxx"));
  fileList.AddLast(new TObjString("AliAnalysisTaskMuonFakes.h"));
  
  TIter nextFile(&fileList);
  TObjString *file;
  char overwrite = '\0';
  while ((file = static_cast<TObjString*>(nextFile()))) {
    
    if (overwrite != 'a') {
      
      overwrite = '\0';
      
      if (!gSystem->AccessPathName(file->GetName())) {
	
	while (overwrite != 'y' && overwrite != 'n' && overwrite != 'a' && overwrite != 'k') {
	  cout<<Form("file %s exist in current directory. Overwrite? [y=yes, n=no, a=all, k=keep all] ",file->GetName())<<flush;
	  cin>>overwrite;
	}
	
      } else overwrite = 'y';
      
    }
    
    if (overwrite == 'y' || overwrite == 'a') {
      
      if (!gSystem->AccessPathName(Form("%s/%s", path1.Data(), file->GetName())))
	gSystem->Exec(Form("cp %s/%s %s", path1.Data(), file->GetName(), file->GetName()));
      else if (!gSystem->AccessPathName(Form("%s/%s", path2.Data(), file->GetName())))
	gSystem->Exec(Form("cp %s/%s %s", path2.Data(), file->GetName(), file->GetName()));
      else
	gSystem->Exec(Form("cp %s/%s %s", path3.Data(), file->GetName(), file->GetName()));
      
    } else if (overwrite == 'k') break;
    
  }
  
}

//______________________________________________________________________________
void LoadAlirootLocally(TString& extraLibs)
{
  /// Load libraries locally
  
  // Load common libraries
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
  
  // Load additional libraries
  gSystem->Load("libProofPlayer");
  TObjArray* libs = extraLibs.Tokenize(":");
  for (Int_t i = 0; i < libs->GetEntriesFast(); i++)
    gSystem->Load(Form("lib%s",static_cast<TObjString*>(libs->UncheckedAt(i))->GetName()));
  delete libs;
  
  // Use AliRoot includes and compile our task
  gROOT->ProcessLine(".include $ALICE_ROOT/PWG3/base");
  gROOT->ProcessLine(".include $ALICE_ROOT/include");
  gROOT->ProcessLine(".include $ALICE_ROOT/MUON");
  gROOT->ProcessLine(".include $ALICE_ROOT/MUON/mapping");
  gROOT->LoadMacro("AliAnalysisTaskMuonPhysics.cxx++g");
  gROOT->LoadMacro("AliAnalysisTaskMuonRefit.cxx++g");
  gROOT->LoadMacro("AliAnalysisTaskMuonFakes.cxx++g");
  
}

//______________________________________________________________________________
void LoadAlirootOnProof(TString& aaf, TString alirootVersion, TString& extraLibs)
{
  /// Load aliroot packages and set environment on Proof
  
  // set general environment and close previous session
  gEnv->SetValue("XSec.GSI.DelegProxy","2");
  
  // connect
  TString location = (aaf == "caf") ? "alice-caf.cern.ch" : "nansafmaster.in2p3.fr";
  //TString location = (aaf == "caf") ? "alice-caf.cern.ch" : "localhost:1093";
  TString nWorkers = (aaf == "caf") ? "workers=40" : "";
  if (gSystem->Getenv("alien_API_USER") == NULL) TProof::Open(location.Data(), nWorkers.Data());
  else TProof::Open(Form("%s@%s",gSystem->Getenv("alien_API_USER"), location.Data()), nWorkers.Data());
  if (!gProof) return;
  
  // set environment and load libraries on workers
  TList* list = new TList();
  list->Add(new TNamed("ALIROOT_MODE", ""));
  list->Add(new TNamed("ALIROOT_EXTRA_LIBS", extraLibs.Data()));
  list->Add(new TNamed("ALIROOT_EXTRA_INCLUDES", "MUON:MUON/mapping:PWG3/base"));
  gProof->EnablePackage(alirootVersion.Data(), list, kTRUE);
  
  // compile task on workers
  gProof->Load("AliAnalysisTaskMuonPhysics.cxx++g", kTRUE);
  gProof->Load("AliAnalysisTaskMuonRefit.cxx++g", kTRUE);
  gProof->Load("AliAnalysisTaskMuonFakes.cxx++g", kTRUE);
  
}

//______________________________________________________________________________
void CreateAnalysisTrain(Bool_t useMCLabels, Double_t resNB, Double_t resB,
			 Double_t sigmaCut, Double_t sigmaCutTrig, Long_t alienHandler)
{
  /// create the analysis train and configure it
  
  // analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("MuonRefitAnalysis");
  
  // Connect plugin to the analysis manager if any
  if (alienHandler) mgr->SetGridHandler(static_cast<AliAnalysisGrid*>(alienHandler));
  
  // ESD input
  AliESDInputHandler* esdH = new AliESDInputHandler();
  esdH->SetReadFriends(kFALSE);
  esdH->SetInactiveBranches(" FMD PHOS  EMCAL  Pmd Trd V0s TPC "
			    "Cascades Kinks CaloClusters ACORDE RawData HLT TZERO ZDC"
			    " Cells ACORDE Pileup");
  mgr->SetInputEventHandler(esdH);
  
  // Monte Carlo handler
  AliMCEventHandler* mcHandler = new AliMCEventHandler();
  mgr->SetMCtruthEventHandler(mcHandler);
  
  // event selection
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
  AliPhysicsSelectionTask* physicsSelection = AddTaskPhysicsSelection(kTRUE);
  if (!physicsSelection) {
    Error("CreateAnalysisTrain","AliPhysicsSelectionTask not created!");
    return;
  }
  
  // first task (physics before refit)
  AliAnalysisTaskMuonPhysics *physics1 = new AliAnalysisTaskMuonPhysics("MuonPhysics1");
  mgr->AddTask(physics1);
  mgr->ConnectInput(physics1, 0, mgr->GetCommonInputContainer());
  AliAnalysisDataContainer *histo1 = mgr->CreateContainer("Histograms_before", TObjArray::Class(), AliAnalysisManager::kOutputContainer, AliAnalysisManager::GetCommonFileName());
  mgr->ConnectOutput(physics1, 1, histo1);
  AliAnalysisDataContainer *trackStat1 = mgr->CreateContainer("trackCounters_before", AliCounterCollection::Class(), AliAnalysisManager::kOutputContainer, AliAnalysisManager::GetCommonFileName());
  mgr->ConnectOutput(physics1, 2, trackStat1);
  
  // second task (refit)
  gROOT->LoadMacro("$ALICE_ROOT/PWG3/muondep/AddTaskMuonRefit.C");
  AliAnalysisTaskMuonRefit* refit = AddTaskMuonRefit(resNB, resB, kTRUE, sigmaCut, sigmaCutTrig);
  if(!refit) {
    Error("CreateAnalysisTrain","AliAnalysisTaskMuonRefit not created!");
    return 0x0;
  }
  
  // third task (physics after refit)
  AliAnalysisTaskMuonPhysics *physics2 = new AliAnalysisTaskMuonPhysics("MuonPhysics2");
  mgr->AddTask(physics2);
  mgr->ConnectInput(physics2, 0, mgr->GetCommonInputContainer());
  AliAnalysisDataContainer *histo2 = mgr->CreateContainer("Histograms_after", TObjArray::Class(), AliAnalysisManager::kOutputContainer, AliAnalysisManager::GetCommonFileName());
  mgr->ConnectOutput(physics2, 1, histo2);
  AliAnalysisDataContainer *trackStat2 = mgr->CreateContainer("trackCounters_after", AliCounterCollection::Class(), AliAnalysisManager::kOutputContainer, AliAnalysisManager::GetCommonFileName());
  mgr->ConnectOutput(physics2, 2, trackStat2);
  
  // forth task (fakes study)
  gROOT->LoadMacro("$ALICE_ROOT/PWG3/muondep/AddTaskMuonFakes.C");
  AliAnalysisTaskMuonFakes* muonFakes = AddTaskMuonFakes(useMCLabels);
  if(!muonFakes) {
    Error("CreateAnalysisTrain","AliAnalysisTaskMuonFakes not created!");
    return;
  }
  
}

//______________________________________________________________________________
Int_t GetMode(TString smode, TString input)
{
  if (smode == "local" && input.EndsWith(".root")) return kLocal;    
  else if (smode == "caf" || smode == "saf") return kProof;
  else if ((smode == "test" || smode == "offline" || smode == "submit" ||
	    smode == "full" || smode == "terminate") && input.EndsWith(".txt")) return kGrid;
  return -1;
}

//______________________________________________________________________________
TChain* CreateChainFromFile(const char *rootfile)
{
  // Create a chain using the root file.
  TChain* chain = new TChain("esdTree");
  chain->Add(rootfile);
  if (!chain->GetNtrees()) return NULL;
  chain->ls();
  return chain;
}

