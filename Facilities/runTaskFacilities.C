/*
 *  runTaskFacilities.C
 *  aliroot
 *
 *  Created by Philippe Pillot on 16/12/10.
 *  Copyright 2010 SUBATECH. All rights reserved.
 *
 */

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TObject.h>
#include <TString.h>
#include <TObjString.h>
#include <TSystem.h>
#include <TChain.h>
#include <TProof.h>
#include <TList.h>
#include <THashList.h>
#include <TGrid.h>
#include <TEnv.h>
#include <TROOT.h>
#include <TList.h>
#include <TFile.h>
#include <TFileInfo.h>
#include <TFileCollection.h>
#include "AliAnalysisManager.h"
#include "AliAnalysisAlien.h"
#endif

enum {kLocal, kProof, kProofLite, kGrid, kTerminate, kSAF3Connect};

//______________________________________________________________________________
Int_t GetMode(TString smode, TString input)
{
  // Get runing mode
  Int_t mode = -1;
  if (smode == "local" && (input.EndsWith(".root") || input.EndsWith(".txt"))) mode = kLocal;
  else if (smode == "prooflite" && input.EndsWith(".txt")) mode = kProofLite;
  else if ((smode == "test" || smode == "offline" || smode == "submit" || smode == "full" ||
	    smode == "merge" || smode == "terminate") && input.EndsWith(".txt")) mode = kGrid;
  else if (smode == "terminateonly") mode = kTerminate;
  else if (smode == "caf" || smode == "saf" || smode == "saf3") {
    if (input.EndsWith(".root")) {
      TFile *inFile = TFile::Open(input.Data(),"READ");
      if (inFile && inFile->IsOpen()) {
        TFileCollection *coll = dynamic_cast<TFileCollection*>(inFile->FindObjectAny("dataset"));
        if (coll) {
          mode = kProof;
          delete coll;
        }
        inFile->Close();
      }
    } else if (!input.IsNull()) mode = kProof;
    if (mode == kProof && smode == "saf3" && gSystem->GetFromPipe("hostname") != "nansafmaster3.in2p3.fr")
      mode = kSAF3Connect;
  } else if (smode == "vaf" && !input.IsNull() && !input.EndsWith(".root")) {
    if (gSystem->GetFromPipe("hostname") != "alivaf-002.cern.ch") mode = kSAF3Connect;
    else mode = kProof;
  }
  return mode;
}

//______________________________________________________________________________
TString GetDataType(Int_t mode, TString input, TString dataPattern)
{
  // Get the data type (ESD or AOD)
  if (mode == kGrid || mode == kTerminate) {
    if (dataPattern.Contains("AOD")) return "AOD";
    else return "ESD";
  } else {
    if (input.Contains(".root")) {
      if (mode == kProof) {
        TString dataType = "unknown";
        TFile *inFile = TFile::Open(input.Data(),"READ");
        if (inFile && inFile->IsOpen()) {
          TFileCollection *coll = dynamic_cast<TFileCollection*>(inFile->FindObjectAny("dataset"));
          if (coll && coll->GetList()) {
            TFileInfo *fileInfo = static_cast<TFileInfo*>(coll->GetList()->At(0));
            if (fileInfo) {
              TString fileName = fileInfo->GetCurrentUrl()->GetFile();
              if (fileName.Contains("AOD")) dataType = "AOD";
              else if (fileName.Contains("ESD")) dataType = "ESD";
            }
            delete coll;
          }
          inFile->Close();
        }
        return dataType.Data();
      } else if (input.Contains("AOD")) return "AOD";
      else if (input.Contains("ESD")) return "ESD";
    } else if (input.EndsWith(".txt")) {
      if (gSystem->Exec(Form("grep -q AOD %s", input.Data())) == 0) return "AOD";
      else if (gSystem->Exec(Form("grep -q ESD %s", input.Data())) == 0) return "ESD";
    }
  }
  return "unknown";
}

//______________________________________________________________________________
void CopyFileLocally(TList &pathList, TList &fileList, char overwrite = '\0')
{
  /// Copy files needed for this analysis
  
  if (gSystem->GetFromPipe("hostname") == "nansafmaster3.in2p3.fr") return; // files already there
  else if (gSystem->GetFromPipe("hostname") == "alivaf-002.cern.ch") return; // files already there
  
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
  else if (gSystem->GetFromPipe("hostname") == "alivaf-002.cern.ch") return; // files already there
  
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
Bool_t CopyFileOnAAF(TString aaf, TList &fileList, TString dataset)
{
  /// Copy files needed for this analysis from current to aaf directory
  
  TString aafdir = gSystem->ExpandPathName(Form("$HOME/%s",aaf.Data()));
  if (gSystem->AccessPathName(Form("%s/.vaf", aafdir.Data()))) {
    cout<<aaf.Data()<<" folder is not mounted"<<endl;
    cout<<"please retry as it can take some time to mount it"<<endl;
    return kFALSE;
  }
  
  TString remoteLocation = gSystem->pwd();
  remoteLocation.ReplaceAll(gSystem->Getenv("HOME"),aafdir.Data());
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
  
  if (gSystem->AccessPathName(Form("%s/Work/Alice/Macros/Facilities", aafdir.Data())))
    gSystem->Exec(Form("mkdir -p %s/Work/Alice/Macros/Facilities", aafdir.Data()));
  gSystem->Exec(Form("cp -p $HOME/Work/Alice/Macros/Facilities/runTaskFacilities.C %s/Work/Alice/Macros/Facilities/runTaskFacilities.C", aafdir.Data()));
  gSystem->Exec(Form("cp -p $HOME/Work/Alice/Macros/Facilities/mergeLocally.C %s/Work/Alice/Macros/Facilities/mergeLocally.C", aafdir.Data()));
  
  if (dataset.EndsWith(".txt")) {
    if (aaf == "saf3") gSystem->Exec(Form("cat %s | awk {'print $1\";Mode=cache\"}' > datasetAAF.txt", dataset.Data()));
    else gSystem->Exec(Form("cat %s | awk {'print $1\";Mode=remote;\"}' > datasetAAF.txt", dataset.Data()));
    gSystem->Exec(Form("cp datasetAAF.txt %s/datasetAAF.txt", remoteLocation.Data()));
  } else if (dataset.EndsWith(".root"))
    gSystem->Exec(Form("cp %s %s/%s", dataset.Data(), remoteLocation.Data(), gSystem->BaseName(dataset.Data())));
  
  return kTRUE;
  
}

//______________________________________________________________________________
void LoadAlirootLocally(TString& extraLibs, TString& extraIncs, TString& extraTasks, TString& extraPkgs)
{
  /// Set environment locally
  
  // Load additional libraries
  if (!extraLibs.IsNull()) {
    TObjArray* libs = extraLibs.Tokenize(":");
    for (Int_t i = 0; i < libs->GetEntriesFast(); i++)
      gSystem->Load(Form("lib%s",static_cast<TObjString*>(libs->UncheckedAt(i))->GetName()));
    delete libs;
  }
  
  // Add additional includes
  if (!extraIncs.IsNull()) {
    TObjArray* incs = extraIncs.Tokenize(":");
    for (Int_t i = 0; i < incs->GetEntriesFast(); i++) {
      gROOT->ProcessLine(Form(".include $ALICE_ROOT/%s",static_cast<TObjString*>(incs->UncheckedAt(i))->GetName()));
      gROOT->ProcessLine(Form(".include $ALICE_PHYSICS/%s",static_cast<TObjString*>(incs->UncheckedAt(i))->GetName()));
    }
    delete incs;
  }
  
  // Optionally add packages (need to change DYLD_LIBRARY_PATH on mac to load only the correct library)
/*  TString dyld_library_path = gSystem->Getenv("DYLD_LIBRARY_PATH");
  if (!dyld_library_path.BeginsWith(".")) {
    gSystem->Exec("export DYLD_LIBRARY_PATH=.:$DYLD_LIBRARY_PATH");
    printf("DYLD_LIBRARY_PATH = %s\n",gSystem->Getenv("DYLD_LIBRARY_PATH"));
  }*/
  if (!extraPkgs.IsNull()) {
    TObjArray* pkgs = extraPkgs.Tokenize(":");
    for (Int_t i = 0; i < pkgs->GetEntriesFast(); i++)
      AliAnalysisAlien::SetupPar(Form("%s.par",static_cast<TObjString*>(pkgs->UncheckedAt(i))->GetName()));
    delete pkgs;
  }
  
  // Compile additional tasks
  if (!extraTasks.IsNull()) {
    TObjArray* tasks = extraTasks.Tokenize(":");
    for (Int_t i = 0; i < tasks->GetEntriesFast(); i++)
      gROOT->LoadMacro(Form("%s.cxx+",static_cast<TObjString*>(tasks->UncheckedAt(i))->GetName()));
    delete tasks;
  }
  
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
  else if (aaf == "saf3" || aaf == "vaf") TProof::Open("pod://");
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
  list->Add(new TNamed("ALIROOT_ENABLE_ALIEN", "1"));
  if (aaf == "prooflite") {
    gProof->UploadPackage("$ALICE_ROOT/ANALYSIS/macros/AliRootProofLite.par");
    gProof->EnablePackage("$ALICE_ROOT/ANALYSIS/macros/AliRootProofLite.par", list);
  } else if (aaf == "saf3") {
    TString home = gSystem->Getenv("HOME");
    gProof->UploadPackage(Form("%s/AliceVaf.par", home.Data()));
    gProof->EnablePackage(Form("%s/AliceVaf.par", home.Data()), list);
  } else if (aaf == "vaf") {
    TFile::Cp("http://alibrary.web.cern.ch/alibrary/vaf/AliceVaf.par", "AliceVaf.par");
    gProof->UploadPackage("AliceVaf.par");
    gProof->EnablePackage("AliceVaf.par", list);
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
    gProof->Load(Form("%s.cxx++g",static_cast<TObjString*>(tasks->UncheckedAt(i))->GetName()), notOnClient);
  delete tasks;
  
}

//______________________________________________________________________________
Int_t PrepareAnalysis(TString smode, TString input, TString &extraLibs, TString &extraIncs, TString &extraTasks,
                      TString &extraPkgs, TList &pathList, TList &fileList, char overwrite = '\0')
{
  /// Prepare the analysis environment
  
  // --- Check runing mode ---
  Int_t mode = GetMode(smode, input);
  if(mode < 0){
    Error("runTaskFacility","invalid mode or incompatible input");
    exit(1);
  }
  
  // --- copy files needed for this analysis ---
  CopyFileLocally(pathList, fileList, overwrite);
  
  // --- prepare environment locally (not needed if runing on aaf) ---
  if (mode != kSAF3Connect) LoadAlirootLocally(extraLibs, extraIncs, extraTasks, extraPkgs);
  
  return mode;
  
}

//______________________________________________________________________________
TObject* CreateAlienHandler(TString runMode, TString& alirootVersion, TString& aliphysicsVersion, TString& runListName,
			    TString &dataDir, TString &dataPattern, TString &outDir, TString& extraLibs,
			    TString& extraIncs, TString& extraTasks, TString& extraPkgs, TString& analysisMacroName,
			    TString runFormat = "%09d", Int_t ttl = 30000, Int_t maxFilesPerJob = 100,
			    Int_t maxMergeFiles = 10, Int_t maxMergeStages = 1, TString RegisterExcludes = "")
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
  plugin->SetNtestFiles(2);
  
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
  
  // Exclude given output file(s) from registration/merging
  plugin->SetRegisterExcludes(RegisterExcludes.Data());
  
  // Save the log files
  plugin->SetKeepLogs();
  
  return plugin;
}

//______________________________________________________________________________
TChain* CreateChainFromFile(const char *file)
{
  // Create a chain using the root file or txt file containing root file
  
  TString dataType = GetDataType(kLocal, file, "");
  if (dataType == "unknown") return NULL;
  
  TChain* chain = (dataType == "AOD") ? new TChain("aodTree") : new TChain("esdTree");
  
  if (strstr(file,".root")) {
    
    if(chain->Add(file, -1) == 0) return NULL;
    
  } else {
    
    char line[1024];
    ifstream in(file);
    while (in.getline(line,1024,'\n')) if(chain->Add(line, -1) == 0) return NULL;
    
  }
  
  return chain;
  
}

//______________________________________________________________________________
TObject* CreateInputObject(Int_t mode, TString &inputName)
{
  // Create input object
  if (mode == kProof) {
    if (inputName.EndsWith(".root")) {
      TFile *inFile = TFile::Open(inputName.Data(),"READ");
      if (!inFile || !inFile->IsOpen()) return NULL;
      TFileCollection *coll = dynamic_cast<TFileCollection*>(inFile->FindObjectAny("dataset"));
      inFile->Close();
      return coll;
    } else return new TObjString(inputName);
  } else if (mode == kProofLite) {
    TFileCollection *coll = new TFileCollection();
    coll->AddFromFile(inputName.Data());
    gProof->RegisterDataSet("test_collection", coll, "OV");
    return coll;
  } else if (mode == kLocal) return CreateChainFromFile(inputName.Data());
  return NULL;
}

//______________________________________________________________________________
void CreateAAFExecutable(TString aaf, TString dataset, Bool_t splitDataset, Bool_t overwriteDataset = kFALSE)
{
  /// Create the executable to run on AAF
  ofstream outFile("runAnalysis.sh");
  outFile << "#!/bin/bash" << endl;
  outFile << "ln -s -f $HOME/Work/Alice/Macros/Facilities/runTaskFacilities.C" << endl;
  outFile << "vafctl start" << endl;
  Int_t nWorkers = (aaf == "saf3") ? 88 : 48;
  outFile << "nWorkers=" << nWorkers << endl;
  outFile << "let \"nWorkers -= `pod-info -n`\"" << endl;
  outFile << "echo \"requesting $nWorkers additional workers\"" << endl;
  outFile << "vafreq $nWorkers" << endl;
  outFile << "vafwait " << nWorkers - 8 << endl;
  TString macro = gSystem->GetFromPipe("tail -n 1 $HOME/.root_hist | sed 's/(.*)//g;s/^.x //g;s:^.*/::g'");
  TString arg = gSystem->GetFromPipe("tail -n 1 $HOME/.root_hist | sed 's/.*(/(/g'");
  arg.ReplaceAll("'", "'\"'\"'");
  if (dataset.EndsWith(".txt")) {
    if (splitDataset) {
      if (overwriteDataset) outFile << "rm -rf results" << endl;
      else {
        outFile << "if [[ -d \"results\" ]]; then" << endl;
        outFile << "  answer=\"\"" << endl;
        outFile << "  while ! [ \"$answer\" = \"d\" -o \"$answer\" = \"r\" ]" << endl;
        outFile << "  do" << endl;
        outFile << "    echo \"results already exist: delete or resume? [d/r] \"" << endl;
        outFile << "    read answer" << endl;
        outFile << "  done" << endl;
        outFile << "  if [ $answer = \"d\" ]; then rm -rf results; fi" << endl;
        outFile << "fi" << endl;
      }
      outFile << "if [[ ! -d \"results\" ]]; then mkdir results; fi" << endl;
      outFile << "for ds in `cat datasetAAF.txt`; do" << endl;
      outFile << "  run=`echo $ds | egrep -o '[1-9][0-9][0-9][0-9][0-9][0-9]'`" << endl;
      outFile << "  if [[ -d \"results/$run\" ]]; then continue; fi" << endl;
      outFile << "  mkdir results/$run" << endl;
      outFile << "  echo $ds > results/$run/ds.txt" << endl;
      arg.ReplaceAll(dataset.Data(), "'$ds'");
      outFile << "  root -b -q '" << macro.Data() << arg.Data() << "'" << endl;
      outFile << "  mv AnalysisResults.root results/$run/." << endl;
      outFile << "done" << endl;
      outFile << "ls `pwd`/results/*/AnalysisResults.root > files2merge.txt" << endl;
      outFile << "root -b -q '$HOME/Work/Alice/Macros/Facilities/mergeLocally.C+(\"files2merge.txt\", kTRUE, kFALSE)'" << endl;
      arg.ReplaceAll(aaf.Data(), "terminateonly");
    } else arg.ReplaceAll(dataset.Data(), "datasetAAF.txt");
  } else if (dataset.EndsWith(".root")) arg.ReplaceAll(dataset.Data(), gSystem->BaseName(dataset.Data()));
  outFile << "root -b -q '" << macro.Data() << arg.Data() << "'" << endl;
  outFile << "vafctl stop" << endl;
  outFile.close();
  gSystem->Exec("chmod u+x runAnalysis.sh");
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
    else if (mode == kProof && inputObj) {
      if (inputObj->IsA() == TFileCollection::Class())
        mgr->StartAnalysis("proof", static_cast<TFileCollection*>(inputObj));
      else if (inputObj->IsA() == TObjString::Class())
        mgr->StartAnalysis("proof", Form("%s",static_cast<TObjString*>(inputObj)->GetName()));
      else {
        cout<<"E-StartAnalysis: invalid input object"<<endl;
        return;
      }
    } else if (mode == kProofLite) mgr->StartAnalysis("proof", "test_collection");
    else if (mode == kLocal && inputObj) {
      //mgr->SetUseProgressBar(kTRUE);
      if (inputObj->IsA() == TChain::Class())
        mgr->StartAnalysis("local", static_cast<TChain*>(inputObj));
      else {
        cout<<"E-StartAnalysis: invalid input object"<<endl;
        return;
      }
    }
  }
}

//______________________________________________________________________________
Bool_t RunAnalysis(TString smode, TString input, TString& rootVersion, TString& alirootVersion, TString& aliphysicsVersion,
                   TString &extraLibs, TString &extraIncs, TString &extraTasks, TString &extraPkgs,
                   TString &dataDir, TString &dataPattern, TString &outDir, TString &analysisMacroName,
                   TString runFormat, Int_t ttl, Int_t maxFilesPerJob, Int_t maxMergeFiles, Int_t maxMergeStages,
                   TString RegisterExcludes = "")
{
  /// Run the analysis locally, on proof or on the grid
  
  // --- Get runing mode ---
  Int_t mode = GetMode(smode, input);
  if (mode < 0) return kFALSE;
  
  // --- prepare proof or grid environment ---
  if (mode == kProof || mode == kProofLite) LoadAlirootOnProof(smode, rootVersion, aliphysicsVersion, extraLibs, extraIncs, extraTasks, extraPkgs, kTRUE);
  else if (mode == kGrid || mode == kTerminate) {
    AliAnalysisGrid *alienHandler = static_cast<AliAnalysisGrid*>(CreateAlienHandler(smode, alirootVersion, aliphysicsVersion, input, dataDir, dataPattern, outDir, extraLibs, extraIncs, extraTasks, extraPkgs, analysisMacroName, runFormat, ttl, maxFilesPerJob, maxMergeFiles, maxMergeStages, RegisterExcludes));
    if (alienHandler) AliAnalysisManager::GetAnalysisManager()->SetGridHandler(alienHandler);
    else return kFALSE;
  }
  
  // --- Create input object ---
  TObject* inputObj = CreateInputObject(mode, input);
  
  // --- start analysis ---
  StartAnalysis(mode, inputObj);
  
  return kTRUE;
  
}

//______________________________________________________________________________
Bool_t RunAnalysisOnSAF3(TList &fileList, TString aliphysicsVersion, TString dataset,
                         Bool_t splitDataset, Bool_t overwriteDataset = kFALSE)
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
  CreateAAFExecutable("saf3", dataset, splitDataset, overwriteDataset);
  
  // --- copy files needed for this analysis ---
  if (!CopyFileOnAAF("saf3", fileList, dataset)) {
    cout << "cp problem" << endl;
    return kFALSE;
  }
  
  // --- change the AliPhysics version on SAF3 ---
  gSystem->Exec(Form("sed -i '' 's/VafAliPhysicsVersion.*/VafAliPhysicsVersion=\"%s\"/g' $HOME/saf3/.vaf/vaf.conf", aliphysicsVersion.Data()));
  
  // --- enter SAF3 and run analysis ---
  TString analysisLocation = gSystem->pwd();
  analysisLocation.ReplaceAll(Form("%s/", gSystem->Getenv("HOME")), "");
  gSystem->Exec(Form("gsissh -p 1975 -t -Y nansafmaster3.in2p3.fr 'cd %s; ~/saf3-enter \"\" \"./runAnalysis.sh 2>&1 | tee runAnalysis.log; exit\"'", analysisLocation.Data()));
  
  // --- copy analysis results (assuming analysis run smootly) ---
  gSystem->Exec(Form("cp -p %s/%s/*.root .", saf3dir.Data(), analysisLocation.Data()));
  
  return kTRUE;
  
}

//______________________________________________________________________________
Bool_t RunAnalysisOnVAF(TList &fileList, TString aliphysicsVersion, TString dataset,
                        Bool_t splitDataset, Bool_t overwriteDataset = kFALSE)
{
  /// Run analysis on VAF
  
  // --- need a tunnel to alivaf-002 ---
  TString user = (gSystem->Getenv("alien_API_USER") == NULL) ? "" : Form("%s@",gSystem->Getenv("alien_API_USER"));
  if (gSystem->Exec("nc -z localhost 5501") != 0) {
    gSystem->Exec(Form("xterm -e 'ssh %slxplus.cern.ch -L 5501:alivaf-002:22' &", user.Data()));
    while (gSystem->Exec("nc -z localhost 5501") != 0);
  }
  
  // --- mount alivaf-002 ---
  TString vafdir = gSystem->ExpandPathName("$HOME/vaf");
  if (gSystem->AccessPathName(vafdir.Data())) gSystem->Exec(Form("mkdir %s", vafdir.Data()));
  if (gSystem->AccessPathName(Form("%s/.vaf", vafdir.Data()))) {
    Int_t ret = gSystem->Exec(Form("sshfs -o ssh_command=\"gsissh -p5501\" %slocalhost: %s", user.Data(), vafdir.Data()));
    if (ret != 0) {
      cout<<"mounting of vaf folder failed"<<endl;
      return kFALSE;
    }
  }
  
  // --- create the executable to run on VAF ---
  CreateAAFExecutable("vaf", dataset, splitDataset, overwriteDataset);
  
  // --- copy files needed for this analysis ---
  if (!CopyFileOnAAF("vaf", fileList, dataset)) {
    cout << "cp problem" << endl;
    return kFALSE;
  }
  
  // --- change the AliPhysics version on VAF ---
  gSystem->Exec(Form("sed -i '' 's/VafAliPhysicsVersion.*/VafAliPhysicsVersion=\"%s\"/g' $HOME/vaf/.vaf/vaf.conf", aliphysicsVersion.Data()));
  
  // --- enter VAF and run analysis ---
  TString analysisLocation = gSystem->pwd();
  analysisLocation.ReplaceAll(Form("%s/", gSystem->Getenv("HOME")), "");
  gSystem->Exec(Form("gsissh -p 5501 -t -Y %slocalhost 'cd %s; ~/vaf-enter \"\" \"./runAnalysis.sh 2>&1 | tee runAnalysis.log; exit\"'", user.Data(), analysisLocation.Data()));
  
  // --- copy analysis results (assuming analysis run smootly) ---
  gSystem->Exec(Form("cp -p %s/%s/*.root .", vafdir.Data(), analysisLocation.Data()));
  
  return kTRUE;
  
}

