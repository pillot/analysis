/*
 *  runLocalESDCheck.C
 *  aliroot
 *
 *  Created by Philippe Pillot on 05/04/10.
 *  Copyright 2010 SUBATECH. All rights reserved.
 *
 */

enum {kLocal, kInteractif_xml, kInteractif_ESDList};

void runLocalESDCheck(TString inputFileName = "AliESDs.root", Bool_t fullPrintout = kTRUE,
		      Bool_t physicsSelect = kTRUE, Short_t selectCharge = 0)
{
  TStopwatch timer;
  timer.Start();
  
  // Check runing mode
  Int_t mode = GetMode(inputFileName);
  if(mode < 0){
    printf("Fatal: Please provide either an ESD root file or a collection of ESDs.\nExit!");
    return;
  }
  
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
  gROOT->ProcessLine(".include /Users/philippe/Work/Alice/Work/Data/Macro");
  
  // OCDB access
  if(mode==kLocal) AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  else {
    if (!TGrid::Connect("alien://")) return;
    if(mode==kInteractif_ESDList || !gSystem->AccessPathName("ConfigureCuts.C"))
      AliCDBManager::Instance()->SetDefaultStorage("raw://");
  }
  
  // Create input chain
  TChain* chain = CreateChain(inputFileName);
  if (!chain) return;
  
  // Create the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("ESDCheckAnalysis");
  
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
  gROOT->LoadMacro("/Users/philippe/Work/Alice/Work/Data/Macro/AliCounterCollection.cxx++g");
  gROOT->LoadMacro("/Users/philippe/Work/Alice/Work/Data/Macro/ESDCheck/AliAnalysisTaskESDCheck.cxx++g");
  AliAnalysisTaskESDCheck* task = new AliAnalysisTaskESDCheck("ESDCheck");
  task->SetFullPrintout(fullPrintout);
  task->SelectCharge(selectCharge);
  mgr->AddTask(task);
  
  // Create containers for output
  TString outFileName = "ESDCheck.root";
  AliAnalysisDataContainer *cout_histo1 = mgr->CreateContainer("general", TObjArray::Class(), AliAnalysisManager::kOutputContainer, outFileName.Data());
  AliAnalysisDataContainer *cout_histo2 = mgr->CreateContainer("expert", TObjArray::Class(), AliAnalysisManager::kOutputContainer, outFileName.Data());
  AliAnalysisDataContainer *cout_trackStat = mgr->CreateContainer("trackCounters", AliCounterCollection::Class(), AliAnalysisManager::kOutputContainer, outFileName.Data());
  AliAnalysisDataContainer *cout_eventStat = mgr->CreateContainer("eventCounters", AliCounterCollection::Class(), AliAnalysisManager::kOutputContainer, outFileName.Data());
  
  // Connect input/output
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, cout_histo1);
  mgr->ConnectOutput(task, 2, cout_histo2);
  mgr->ConnectOutput(task, 3, cout_trackStat);
  mgr->ConnectOutput(task, 4, cout_eventStat);
  
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

