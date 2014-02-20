/*
 *  runLocalMuonQA.C
 *  aliroot
 *
 *  Created by Philippe Pillot on 29/09/10.
 *  Copyright 2010 SUBATECH. All rights reserved.
 *
 */

Bool_t taskInPWG3 = kFALSE;

enum {kLocal, kInteractif_xml, kInteractif_ESDList};

void runLocalMuonQA(TString inputFileName = "AliESDs.root", Bool_t selectPhysics = kTRUE,
	       Bool_t selectTrigger = kFALSE, Short_t selectCharge = 0)
{
  TStopwatch timer;
  timer.Start();
  
  // Check runing mode
  Int_t mode = GetMode(inputFileName);
  if(mode < 0){
    Error("RunMuonQA","Please provide either an ESD root file or a collection of ESDs.");
    return;
  }
  
  // copy files needed for this analysis
  TString path1("/Users/philippe/Work/Alice/Work/Data/Macro/MuonQA");
  TString path2("$ALICE_ROOT/PWG3/muon");
  TString path3("$ALICE_ROOT/PWG3/base");
  TObjArray fileList(100);
  fileList.SetOwner();
  fileList.AddLast(new TObjString("runMuonQA.C"));
  fileList.AddLast(new TObjString("CreateAlienHandler.C"));
  fileList.AddLast(new TObjString("ConfigureCuts.C"));
  if (!taskInPWG3) {
    fileList.AddLast(new TObjString("AliCounterCollection.cxx"));
    fileList.AddLast(new TObjString("AliCounterCollection.h"));
    fileList.AddLast(new TObjString("AliAnalysisTaskMuonQA.cxx"));
    fileList.AddLast(new TObjString("AliAnalysisTaskMuonQA.h"));
    fileList.AddLast(new TObjString("AddTaskMuonQA.C"));
  }
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
      else if (!taskInPWG3) {
	if (!gSystem->AccessPathName(Form("%s/%s", gSystem->ExpandPathName(path2.Data()), file->GetName())))
	  gSystem->Exec(Form("cp %s/%s %s", path2.Data(), file->GetName(),file->GetName()));
	else
	  gSystem->Exec(Form("cp %s/%s %s", path3.Data(), file->GetName(),file->GetName()));
      }
    } else if (overwrite == 'k') break;
  }
  
  // Load common libraries
  gSystem->Load("libTree");
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");
  gSystem->Load("libPhysics");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  
  // load PWG3muon library or compile the task if not in there
  if (taskInPWG3) {
    gSystem->Load("libCORRFW");
    gSystem->Load("libPWG3base");
    gSystem->Load("libPWG3muon");
  } else {
    gROOT->LoadMacro("AliCounterCollection.cxx++g");
    gROOT->LoadMacro("AliAnalysisTaskMuonQA.cxx++g");
  }
  
  // Create input chain
  TChain* chain = CreateChain(inputFileName);
  if (!chain) return;
  
  // Create the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("MuonQAAnalysis");
  
  // ESD input handler
  AliESDInputHandler* esdH = new AliESDInputHandler();
  esdH->SetReadFriends(kFALSE);
  mgr->SetInputEventHandler(esdH);
  
  // event selection
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
  AliPhysicsSelectionTask* physicsSelection = AddTaskPhysicsSelection();
  if(!physicsSelection) {
    Error("RunMuonQA","AliPhysicsSelectionTask not created!");
    return;
  }
  
  // Muon QA analysis
  gROOT->LoadMacro("$ALICE_ROOT/PWG3/muon/AddTaskMuonQA.C");
  AliAnalysisTaskMuonQA* muonQA = AddTaskMuonQA(selectPhysics, selectTrigger, selectCharge);
  if(!muonQA) {
    Error("RunMuonQA","AliAnalysisTaskMuonQA not created!");
    return;
  }
  
  // Enable debug printouts
  //mgr->SetDebugLevel(2);
  
  // start local analysis
  if (mgr->InitAnalysis()) {
    mgr->PrintStatus();
    mgr->StartAnalysis("local", chain);
  }
  
  timer.Stop();
  timer.Print();
}

//______________________________________________________________________________
Int_t GetMode(TString inputFileName)
{
  if ( inputFileName.EndsWith(".xml") ) return kInteractif_xml;
  else if ( inputFileName.EndsWith(".txt") ) return kInteractif_ESDList;
  else if ( inputFileName.EndsWith(".root") ) return kLocal;
  return -1;
}

//______________________________________________________________________________
TChain* CreateChainFromCollection(const char *xmlfile)
{
  // Create a chain from the collection of tags.
  TAlienCollection* coll = TAlienCollection::Open(xmlfile);
  if (!coll) {
    ::Error("CreateChainFromTags", "Cannot create an AliEn collection from %s", xmlfile);
    return NULL;
  }
  
  TGridResult* tagResult = coll->GetGridResult("",kFALSE,kFALSE);
  AliTagAnalysis *tagAna = new AliTagAnalysis("ESD");
  tagAna->ChainGridTags(tagResult);
  
  AliRunTagCuts      *runCuts = new AliRunTagCuts();
  AliLHCTagCuts      *lhcCuts = new AliLHCTagCuts();
  AliDetectorTagCuts *detCuts = new AliDetectorTagCuts();
  AliEventTagCuts    *evCuts  = new AliEventTagCuts();
  
  // Check if the cuts configuration file was provided
  if (!gSystem->AccessPathName("ConfigureCuts.C")) {
    gROOT->LoadMacro("ConfigureCuts.C");
    ConfigureCuts(runCuts, lhcCuts, detCuts, evCuts);
  }
  
  TChain *chain = tagAna->QueryTags(runCuts, lhcCuts, detCuts, evCuts);
  if (!chain || !chain->GetNtrees()) return NULL;
  chain->ls();
  return chain;
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

//______________________________________________________________________________
TChain* CreateChainFromESDList(const char *esdList)
{
  // Create a chain using tags from the run list.
  TChain* chain = new TChain("esdTree");
  ifstream inFile(esdList);
  TString inFileName;
  if (inFile.is_open()) {
    while (! inFile.eof() ) {
      inFileName.ReadLine(inFile,kFALSE);
      if(!inFileName.EndsWith(".root")) continue;
      chain->Add(inFileName.Data());
    }
  }
  inFile.close();
  if (!chain->GetNtrees()) return NULL;
  chain->ls();
  return chain;
}

//______________________________________________________________________________
TChain* CreateChain(TString inputFileName)
{
  printf("*******************************\n");
  printf("*** Getting the Chain       ***\n");
  printf("*******************************\n");
  Int_t mode = GetMode(inputFileName);
  if(mode == kInteractif_xml) return CreateChainFromCollection(inputFileName.Data());
  else if (mode == kInteractif_ESDList) return CreateChainFromESDList(inputFileName.Data());
  else if (mode == kLocal) return CreateChainFromFile(inputFileName.Data());
  else return NULL;
}

