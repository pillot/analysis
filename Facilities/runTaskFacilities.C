/*
 *  runTaskFacilities.C
 *  aliroot
 *
 *  Created by Philippe Pillot on 16/12/10.
 *  Copyright 2010 SUBATECH. All rights reserved.
 *
 */

enum {kLocal, kProof, kProofLite, kGrid, kTerminate, kSAF3Connect};

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
  else if (smode == "saf3")
    return (gSystem->GetFromPipe("hostname") == "nansafmaster3.in2p3.fr") ? kProof : kSAF3Connect;
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
  
  if (gSystem->GetFromPipe("hostname") == "nansafmaster3.in2p3.fr") return; // files already there
  
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
void CopyInputFileLocally(TString inFile, TString outFileName, char overwrite = '\0')
{
  /// Copy an input file needed for this analysis eventually changing its name
  
  if (gSystem->GetFromPipe("hostname") == "nansafmaster3.in2p3.fr") return; // files already there
  
  if (overwrite != 'a' && overwrite != 'k') {
    
    overwrite = '\0';
    
    if (!gSystem->AccessPathName(outFileName.Data())) {
      
      while (overwrite != 'y' && overwrite != 'n') {
        cout<<Form("input file %s exist in current directory. Overwrite? [y=yes, n=no] ",outFileName.Data())<<flush;
        cin>>overwrite;
      }
      
    } else overwrite = 'y';
    
  }
  
  if (overwrite == 'y' || overwrite == 'a' || (overwrite == 'k' && gSystem->AccessPathName(outFileName.Data()))) {
    
    if (!gSystem->AccessPathName(Form("%s", gSystem->ExpandPathName(inFile.Data()))))
      gSystem->Exec(Form("cp -p %s %s", inFile.Data(), outFileName.Data()));
    else cout<<Form("file %s not found\n",inFile.Data())<<flush;
    
  }
  
}

//______________________________________________________________________________
Bool_t CopyFileOnSAF3(TList &fileList, TString dataset)
{
  /// Copy files needed for this analysis from current to saf3 directory
  
  TString saf3dir = gSystem->ExpandPathName("$HOME/saf3");
  if (gSystem->AccessPathName(Form("%s/.vaf", saf3dir.Data()))) {
    cout<<"saf3 folder is not mounted"<<endl;
    cout<<"please retry as it can take some time to mount it"<<endl;
    return kFALSE;
  }
  
  TString remoteLocation = gSystem->pwd();
  remoteLocation.ReplaceAll(gSystem->Getenv("HOME"),saf3dir.Data());
  if (gSystem->AccessPathName(remoteLocation.Data())) gSystem->Exec(Form("mkdir -p %s", remoteLocation.Data()));
  
  TIter nextFile(&fileList);
  TObjString *file;
  while ((file = static_cast<TObjString*>(nextFile()))) {
    
    if (gSystem->AccessPathName(file->GetName())) {
      cout<<Form("file %s not found in current directory\n",file->GetName())<<flush;
      return kFALSE;
    }
    
    TString remoteFile = remoteLocation+"/"+file->GetName();
    gSystem->Exec(Form("cp -p %s %s", file->GetName(), remoteFile.Data()));
    
  }
  
  gSystem->Exec(Form("cp runAnalysis.sh %s/runAnalysis.sh", remoteLocation.Data()));
  
  if (gSystem->AccessPathName(Form("%s/Work/Alice/Macros/Facilities", saf3dir.Data())))
    gSystem->Exec(Form("mkdir -p %s/Work/Alice/Macros/Facilities", saf3dir.Data()));
  gSystem->Exec("cp -p $HOME/Work/Alice/Macros/Facilities/runTaskFacilities.C $HOME/saf3/Work/Alice/Macros/Facilities/runTaskFacilities.C");
  
  if (dataset.EndsWith(".txt")) {
    gSystem->Exec(Form("cat %s | awk {'print $1\";Mode=cache\"}' > datasetSaf3.txt", dataset.Data()));
    gSystem->Exec(Form("cp datasetSaf3.txt %s/datasetSaf3.txt", remoteLocation.Data()));
  }
  
  return kTRUE;
  
}

//______________________________________________________________________________
void LoadAlirootLocally(TString& extraLibs, TString& extraIncs, TString& extraTasks, TString& extraPkgs)
{
  /// Set environment locally
  
  // Load additional libraries
  TObjArray* libs = extraLibs.Tokenize(":");
  for (Int_t i = 0; i < libs->GetEntriesFast(); i++)
    gSystem->Load(Form("lib%s",static_cast<TObjString*>(libs->UncheckedAt(i))->GetName()));
  delete libs;
  
  // Add additional includes
  TObjArray* incs = extraIncs.Tokenize(":");
  for (Int_t i = 0; i < incs->GetEntriesFast(); i++) {
    gROOT->ProcessLine(Form(".include $ALICE_ROOT/%s",static_cast<TObjString*>(incs->UncheckedAt(i))->GetName()));
    gROOT->ProcessLine(Form(".include $ALICE_PHYSICS/%s",static_cast<TObjString*>(incs->UncheckedAt(i))->GetName()));
  }
  delete incs;
  
  // Optionally add packages
  TObjArray* pkgs = extraPkgs.Tokenize(":");
  for (Int_t i = 0; i < pkgs->GetEntriesFast(); i++) {
    AliAnalysisAlien::SetupPar(Form("%s.par",static_cast<TObjString*>(pkgs->UncheckedAt(i))->GetName()));
  }
  delete pkgs;
  
  // Compile additional tasks
  TObjArray* tasks = extraTasks.Tokenize(":");
  for (Int_t i = 0; i < tasks->GetEntriesFast(); i++)
    gROOT->LoadMacro(Form("%s.cxx+",static_cast<TObjString*>(tasks->UncheckedAt(i))->GetName()));
  delete tasks;
  
}

//______________________________________________________________________________
void LoadAlirootOnProof(TString& aaf, TString rootVersion, TString aliphysicsVersion, TString& extraLibs,
			TString& extraIncs, TString& extraTasks, TString& extraPkgs,
                        Bool_t notOnClient = kFALSE, TString alirootMode = "base")
{
  /// Load aliroot packages and set environment on Proof
  
  // set general environment
  gEnv->SetValue("XSec.GSI.DelegProxy","2");
  
  // connect
  if (aaf == "prooflite") TProof::Open("");
  //  if (aaf == "prooflite") TProof::Open("workers=2");
  else if (aaf == "saf3") TProof::Open("pod://");
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
//  list->Add(new TNamed("ALIROOT_ENABLE_ALIEN", "1"));
  if (aaf == "prooflite") {
    gProof->UploadPackage("$ALICE_ROOT/ANALYSIS/macros/AliRootProofLite.par");
    gProof->EnablePackage("$ALICE_ROOT/ANALYSIS/macros/AliRootProofLite.par", list);
  } else if (aaf == "saf3") {
    TString home = gSystem->Getenv("HOME");
    gProof->UploadPackage(Form("%s/AliceVaf.par", home.Data()));
    gProof->EnablePackage(Form("%s/AliceVaf.par", home.Data()), list);
  } else gProof->EnablePackage(Form("VO_ALICE@AliPhysics::%s",aliphysicsVersion.Data()), list, notOnClient);
  
  // Optionally add packages
  TObjArray* pkgs = extraPkgs.Tokenize(":");
  for (Int_t i = 0; i < pkgs->GetEntriesFast(); i++) {
    gProof->UploadPackage(Form("%s.par",static_cast<TObjString*>(pkgs->UncheckedAt(i))->GetName()));
    gProof->EnablePackage(Form("%s.par",static_cast<TObjString*>(pkgs->UncheckedAt(i))->GetName()), notOnClient);
  }
  delete pkgs;
  
  // compile additional tasks on workers
  TObjArray* tasks = extraTasks.Tokenize(":");
  for (Int_t i = 0; i < tasks->GetEntriesFast(); i++)
    gProof->Load(Form("%s.cxx+g",static_cast<TObjString*>(tasks->UncheckedAt(i))->GetName()), notOnClient);
  delete tasks;
  
}

//______________________________________________________________________________
TObject* CreateAlienHandler(TString runMode, TString& alirootVersion, TString& aliphysicsVersion, TString& runListName,
			    TString &dataDir, TString &dataPattern, TString &outDir, TString& extraLibs,
			    TString& extraIncs, TString& extraTasks, TString& extraPkgs, TString& analysisMacroName,
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
  if (!alirootVersion.IsNull()) plugin->SetAliROOTVersion(alirootVersion.Data());
  if (!aliphysicsVersion.IsNull()) plugin->SetAliPhysicsVersion(aliphysicsVersion.Data());
  
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
  for (Int_t i = 0; i < incs->GetEntriesFast(); i++) {
    plugin->AddIncludePath(Form("-I$ALICE_ROOT/%s",static_cast<TObjString*>(incs->UncheckedAt(i))->GetName()));
    plugin->AddIncludePath(Form("-I$ALICE_PHYSICS/%s",static_cast<TObjString*>(incs->UncheckedAt(i))->GetName()));
  }
  delete incs;
  plugin->AddIncludePath("-I.");
  
  // Optionally add packages
  TObjArray* pkgs = extraPkgs.Tokenize(":");
  for (Int_t i = 0; i < pkgs->GetEntriesFast(); i++) {
    plugin->EnablePackage(Form("%s.par",static_cast<TObjString*>(pkgs->UncheckedAt(i))->GetName()));
  }
  delete pkgs;
  
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
  if(!chain->Add(rootfile, -1)) return NULL;
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

//______________________________________________________________________________
Bool_t RunAnalysisOnSAF3(TList &fileList, TString aliphysicsVersion, TString dataset)
{
  /// Run analysis on SAF3
  
  // --- mount nansafmaster3 ---
  TString saf3dir = gSystem->ExpandPathName("$HOME/saf3");
  if (gSystem->AccessPathName(saf3dir.Data())) gSystem->Exec(Form("mkdir %s", saf3dir.Data()));
  if (gSystem->AccessPathName(Form("%s/.vaf", saf3dir.Data()))) {
    Int_t ret = gSystem->Exec(Form("sshfs -o ssh_command=\"gsissh -p1975\" nansafmaster3.in2p3.fr: %s", saf3dir.Data()));
    if (ret != 0) {
      cout<<"mounting of saf3 folder failed"<<endl;
      return kFALSE;
    }
  }
  
  // --- create the executable to run on SAF3 ---
  CreateSAF3Executable(dataset);
  
  // --- copy files needed for this analysis ---
  if (!CopyFileOnSAF3(fileList, dataset)) {
    cout << "cp problem" << endl;
    return kFALSE;
  }
  
  // --- change the AliPhysics version on SAF3 ---
  gSystem->Exec(Form("sed -i '' 's/VafAliPhysicsVersion.*/VafAliPhysicsVersion=\"%s\"/g' $HOME/saf3/.vaf/vaf.conf", aliphysicsVersion.Data()));
  
  // --- enter SAF3 and run analysis ---
  TString analysisLocation = gSystem->pwd();
  analysisLocation.ReplaceAll(Form("%s/", gSystem->Getenv("HOME")), "");
  gSystem->Exec(Form("gsissh -p 1975 -t -Y nansafmaster3.in2p3.fr 'cd %s; ~/saf3-enter \"\" \"./runAnalysis.sh; exit\"'", analysisLocation.Data()));
  
  // --- copy analysis results (assuming analysis run smootly) ---
  gSystem->Exec(Form("cp -p %s/%s/*.root .", saf3dir.Data(), analysisLocation.Data()));
  
  return kTRUE;
  
}

//______________________________________________________________________________
void CreateSAF3Executable(TString dataset)
{
  /// Create the executable to run on SAF3
  ofstream outFile("runAnalysis.sh");
  outFile << "#!/bin/bash" << endl;
  outFile << "vafctl start" << endl;
  Int_t nWorkers = 88;
  outFile << "nWorkers=" << nWorkers << endl;
  outFile << "let \"nWorkers -= `pod-info -n`\"" << endl;
  outFile << "echo \"requesting $nWorkers additional workers\"" << endl;
  outFile << "vafreq $nWorkers" << endl;
  outFile << "vafwait " << nWorkers << endl;
  outFile << "root -b -q ";
  TString macro = gSystem->GetFromPipe("tail -n 1 $HOME/.root_hist | sed 's/(.*)//g;s/^\.x\ //g;s:^.*/::g'");
  TString arg = gSystem->GetFromPipe("tail -n 1 $HOME/.root_hist | sed 's/.*(/(/g'");
  if (dataset.EndsWith(".txt")) arg.ReplaceAll(dataset.Data(), "datasetSaf3.txt");
  outFile << "'" << macro.Data() << arg.Data() << "'" << endl;
  outFile << "vafctl stop" << endl;
  outFile.close();
  gSystem->Exec("chmod u+x runAnalysis.sh");
}

