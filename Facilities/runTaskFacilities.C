/*
 *  runTaskFacilities.C
 *  aliroot
 *
 *  Created by Philippe Pillot on 16/12/10.
 *  Copyright 2010 SUBATECH. All rights reserved.
 *
 */

enum {kLocal, kProof, kProofLite, kGrid, kTerminate};

//______________________________________________________________________________
Int_t GetMode(TString smode, TString input)
{
  // Get runing mode
  if (smode == "local" && input.EndsWith(".root")) return kLocal;    
  else if (smode == "caf" || smode == "saf") return kProof;
  else if (smode == "prooflite" && input.EndsWith(".txt")) return kProofLite;
  else if ((smode == "test" || smode == "offline" || smode == "submit" || smode == "full" ||
	    smode == "merge" || smode == "terminate") && input.EndsWith(".txt")) return kGrid;
  else if (smode == "terminateonly") return kTerminate;
  return -1;
}

//______________________________________________________________________________
TString GetDataType(TString input)
{
  // Get the data type (ESD or AOD)
  if (input.Contains(".root")) {
    if (input.Contains("AOD")) return "AOD";
    else if (input.Contains("ESD")) return "ESD";
  } else if (input.EndsWith(".txt")) {
    if (gSystem->Exec(Form("grep -q AOD %s", input.Data())) == 0) return "AOD";
    else if (gSystem->Exec(Form("grep -q ESD %s", input.Data())) == 0) return "ESD";
  }
  return "unknown";
}

//______________________________________________________________________________
void CopyFileLocally(TList &pathList, TList &fileList, char overwrite = '\0')
{
  /// Copy files needed for this analysis
  
  TIter nextFile(&fileList);
  TObjString *file;
  while ((file = static_cast<TObjString*>(nextFile()))) {
    
    if (overwrite != 'a' && overwrite != 'k') {
      
      overwrite = '\0';
      
      if (!gSystem->AccessPathName(file->GetName())) {
	
	while (overwrite != 'y' && overwrite != 'n' && overwrite != 'a' && overwrite != 'k') {
	  cout<<Form("file %s exist in current directory. Overwrite? [y=yes, n=no, a=all, k=keep all] ",file->GetName())<<flush;
	  cin>>overwrite;
	}
	
      } else overwrite = 'y';
      
    }
    
    if (overwrite == 'y' || overwrite == 'a' || (overwrite == 'k' && gSystem->AccessPathName(file->GetName()))) {
      
      TIter nextPath(&pathList);
      TObjString *path;
      Bool_t copied = kFALSE;
      while ((path = static_cast<TObjString*>(nextPath()))) {
	
	if (!gSystem->AccessPathName(Form("%s/%s", gSystem->ExpandPathName(path->GetName()), file->GetName()))) {
	  gSystem->Exec(Form("cp %s/%s %s", path->GetName(), file->GetName(), file->GetName()));
	  copied = kTRUE;
	  break;
	}
	
      }
      
      if (!copied) cout<<Form("file %s not found in any directory\n",file->GetName())<<flush;
      
    }
    
  }
  
}

//______________________________________________________________________________
void LoadAlirootLocally(TString& extraLibs, TString& extraIncs, TString& extraTasks)
{
  /// Set environment locally
  
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
  
  // Add additional includes
  TObjArray* incs = extraIncs.Tokenize(":");
  for (Int_t i = 0; i < incs->GetEntriesFast(); i++)
    gROOT->ProcessLine(Form(".include $ALICE_ROOT/%s",static_cast<TObjString*>(incs->UncheckedAt(i))->GetName()));
  delete incs;
  
  // Compile additional tasks
  TObjArray* tasks = extraTasks.Tokenize(":");
  for (Int_t i = 0; i < tasks->GetEntriesFast(); i++)
    gROOT->LoadMacro(Form("%s.cxx+",static_cast<TObjString*>(tasks->UncheckedAt(i))->GetName()));
  delete tasks;
  
}

//______________________________________________________________________________
void LoadAlirootOnProof(TString& aaf, TString rootVersion, TString alirootVersion, TString& extraLibs,
			TString& extraIncs, TString& extraTasks, Bool_t notOnClient = kFALSE, TString alirootMode = "")
{
  /// Load aliroot packages and set environment on Proof
  
  // set general environment
  gEnv->SetValue("XSec.GSI.DelegProxy","2");
  
  // connect
  if (aaf == "prooflite") TProof::Open("");
  //  if (aaf == "prooflite") TProof::Open("workers=2");
  else {
    TString location, nWorkers;
    if (aaf == "caf") {
      location = "alice-caf.cern.ch";
      nWorkers = ""; // "workers=40";
    } else if (aaf == "saf") {
      location = "nansafmaster2.in2p3.fr"; // "localhost:1093"
      nWorkers = ""; // "workers=8x";
    } else return;
    TString user = (gSystem->Getenv("alien_API_USER") == NULL) ? "" : Form("%s@",gSystem->Getenv("alien_API_USER"));
    TProof::Mgr(Form("%s%s",user.Data(), location.Data()))->SetROOTVersion(Form("VO_ALICE@ROOT::%s",rootVersion.Data()));
    TProof::Open(Form("%s%s/?N",user.Data(), location.Data()), nWorkers.Data());
  }
  if (!gProof) return;
  
  // set environment and load libraries on workers
  TList* list = new TList();
  list->Add(new TNamed("ALIROOT_MODE", alirootMode.Data()));
  list->Add(new TNamed("ALIROOT_EXTRA_LIBS", extraLibs.Data()));
  list->Add(new TNamed("ALIROOT_EXTRA_INCLUDES", extraIncs.Data()));
  if (aaf == "prooflite") {
    gProof->UploadPackage("$ALICE_ROOT/ANALYSIS/macros/AliRootProofLite.par");
    gProof->EnablePackage("$ALICE_ROOT/ANALYSIS/macros/AliRootProofLite.par", list);
  } else gProof->EnablePackage(Form("VO_ALICE@AliRoot::%s",alirootVersion.Data()), list, notOnClient);
  
  // compile additional tasks on workers
  TObjArray* tasks = extraTasks.Tokenize(":");
  for (Int_t i = 0; i < tasks->GetEntriesFast(); i++)
    gProof->Load(Form("%s.cxx++g",static_cast<TObjString*>(tasks->UncheckedAt(i))->GetName()), notOnClient);
  delete tasks;
  
}

//______________________________________________________________________________
TObject* CreateAlienHandler(TString runMode, TString& rootVersion, TString& alirootVersion, TString& runListName,
			    TString &dataDir, TString &dataPattern, TString &outDir, TString& extraLibs,
			    TString& extraIncs, TString& extraTasks, TString& analysisMacroName,
			    TString runFormat = "%09d", Int_t ttl = 30000, Int_t maxFilesPerJob = 100,
			    Int_t maxMergeFiles = 10, Int_t maxMergeStages = 1)
{
  // Configure the alien plugin
  AliAnalysisAlien *plugin = new AliAnalysisAlien();
  
  // If the run mode is merge, run in mode terminate to merge via jdl
  // If the run mode is terminate, disable the mergin via jdl
  Bool_t merge = kTRUE;
  if (runMode.Contains("terminate")) merge = kFALSE;
  else if (runMode == "merge") runMode = "terminate";
  
  // Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
  plugin->SetRunMode(runMode.Data());
  
  // Set the number of input files in test mode
  plugin->SetNtestFiles(1);
  
  // Set versions of used packages
  plugin->SetAPIVersion("V1.1x");
  //plugin->SetROOTVersion(rootVersion.Data());
  plugin->SetAliROOTVersion(alirootVersion.Data());
  
  // Declare input data to be processed
  plugin->SetGridDataDir(dataDir.Data());
  plugin->SetDataPattern(dataPattern.Data());
  ifstream inFile(runListName.Data());
  TString currRun;
  if (inFile.is_open())
  {
    while (! inFile.eof() )
    {
      currRun.ReadLine(inFile,kTRUE); // Read line
      if(currRun.IsNull()) continue;
      plugin->AddRunNumber(Form(runFormat.Data(), currRun.Atoi()));
    }
  }
  inFile.close();
  
  // Define alien work directory where all files will be copied. Relative to alien $HOME
  plugin->SetGridWorkingDir(outDir.Data());
  
  // Declare alien output directory. Relative to working directory
  plugin->SetGridOutputDir("results");
  
  // Set the ouput directory of each masterjob to the run number
  plugin->SetOutputToRunNo();
  
  // Declare the analysis source file names separated by blancs. To be compiled runtime using ACLiC on the worker nodes
  TString AddTasks;
  TObjArray* tasks = extraTasks.Tokenize(":");
  for (Int_t i = 0; i < tasks->GetEntriesFast(); i++)
    AddTasks += Form("%s.cxx ",static_cast<TObjString*>(tasks->UncheckedAt(i))->GetName());
  plugin->SetAnalysisSource(AddTasks.Data());
  
  // Declare all libraries (other than the default ones for the framework. These will be
  // loaded by the generated analysis macro. Add all extra files (task .cxx/.h) here
  TString AddLibs;
  TObjArray* libs = extraLibs.Tokenize(":");
  for (Int_t i = 0; i < libs->GetEntriesFast(); i++)
    AddLibs += Form("lib%s.so ",static_cast<TObjString*>(libs->UncheckedAt(i))->GetName());
  delete libs;
  for (Int_t i = 0; i < tasks->GetEntriesFast(); i++) {
    AddLibs += Form("%s.cxx ",static_cast<TObjString*>(tasks->UncheckedAt(i))->GetName());
    AddLibs += Form("%s.h ",static_cast<TObjString*>(tasks->UncheckedAt(i))->GetName());
  }
  delete tasks;
  plugin->SetAdditionalLibs(AddLibs.Data());
  plugin->SetAdditionalRootLibs("libGui.so libProofPlayer.so libXMLParser.so");
  
  // Optionally add include paths
  TObjArray* incs = extraIncs.Tokenize(":");
  for (Int_t i = 0; i < incs->GetEntriesFast(); i++)
    plugin->AddIncludePath(Form("-I$ALICE_ROOT/%s",static_cast<TObjString*>(incs->UncheckedAt(i))->GetName()));
  delete incs;
  plugin->AddIncludePath("-I.");
  
  // Optionally set a name for the generated analysis macro (default MyAnalysis.C)
  plugin->SetAnalysisMacro(Form("%s.C",analysisMacroName.Data()));
  plugin->SetExecutable(Form("%s.sh",analysisMacroName.Data()));
  
  // Optionally set time to live (default 30000 sec)
  plugin->SetTTL(ttl);
  
  // Optionally set input format (default xml-single)
  plugin->SetInputFormat("xml-single");
  
  // Optionally modify the name of the generated JDL (default analysis.jdl)
  plugin->SetJDLName(Form("%s.jdl",analysisMacroName.Data()));
  
  // Optionally modify job price (default 1)
  plugin->SetPrice(1);
  
  // Optionally modify split mode (default 'se')
  plugin->SetSplitMode("se");
  
  // Optionally modify the maximum number of files per job
  plugin->SetSplitMaxInputFileNumber(maxFilesPerJob);
  
  // Merge via JDL
  if (merge) plugin->SetMergeViaJDL(kTRUE);
  
  // Optionally set the maximum number of files merged together in one stage
  plugin->SetMaxMergeFiles(maxMergeFiles);
  
  // Optionally set the maximum number of merging stages
  plugin->SetMaxMergeStages(maxMergeStages);
  
  // Exclude given output file(s) from merging
  plugin->SetMergeExcludes("EventStat_temp.root");
  
  // Save the log files
  plugin->SetKeepLogs();
  
  return plugin;
}

//______________________________________________________________________________
TChain* CreateChainFromFile(const char *rootfile)
{
  // Create a chain using the root file
  TChain* chain = (strstr(rootfile,"AOD")) ? new TChain("aodTree") : new TChain("esdTree");
  chain->Add(rootfile);
  if (!chain->GetNtrees()) return NULL;
  chain->ls();
  return chain;
}

//______________________________________________________________________________
TObject* CreateInputObject(Int_t mode, TString &inputName)
{
  // Create input object
  if (mode == kProof) return new TObjString(inputName);
  else if (mode == kProofLite) {
    TFileCollection *coll = new TFileCollection();
    coll->AddFromFile(inputName.Data());
    gProof->RegisterDataSet("test_collection", coll, "OV");
    return coll;
  }
  else if (mode != kGrid && mode != kTerminate) return CreateChainFromFile(inputName.Data());
  return NULL;
}

//______________________________________________________________________________
void StartAnalysis(Int_t mode, TObject* inputObj)
{
  // Start analysis
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (mgr->InitAnalysis()) {
    mgr->PrintStatus();
    if (mode == kGrid) mgr->StartAnalysis("grid");
    else if (mode == kTerminate) mgr->StartAnalysis("grid terminate");
    else if (mode == kProof && inputObj) mgr->StartAnalysis("proof", Form("%s",static_cast<TObjString*>(inputObj)->GetName()));
    else if (mode == kProofLite) mgr->StartAnalysis("proof", "test_collection");
    else if (mode == kLocal && inputObj) {
      //mgr->SetUseProgressBar(kTRUE);
      mgr->StartAnalysis("local", static_cast<TChain*>(inputObj));
    }
  }
}

