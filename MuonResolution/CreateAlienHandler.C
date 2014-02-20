/*
 *  CreateAlienHandler.C
 *  aliroot
 *
 *  Created by Philippe Pillot on 27/06/10.
 *  Copyright 2010 SUBATECH. All rights reserved.
 *
 */

AliAnalysisGrid* CreateAlienHandler(TString runMode, TString runListName, Bool_t taskInPWG3)
{
  // Check if user has a valid token, otherwise make one. This has limitations.
  // One can always follow the standard procedure of calling alien-token-init then
  //   source /tmp/gclient_env_$UID in the current shell.
  if (!AliAnalysisGrid::CreateToken()) return NULL;
  
  AliAnalysisAlien *plugin = new AliAnalysisAlien();
  
  // Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
  plugin->SetRunMode(runMode.Data());  // VERY IMPORTANT - DECRIBED BELOW
  
  // Set versions of used packages
  plugin->SetAPIVersion("V1.1x");
  plugin->SetROOTVersion("v5-26-00b-6");
  plugin->SetAliROOTVersion("v4-19-18-AN");
  
  // Declare input data to be processed.
  // Method 1: Create automatically XML collections using alien 'find' command.
  // Define production directory LFN
//  plugin->SetGridDataDir("/alice/data/2010/LHC10d/");
  plugin->SetGridDataDir("/alice/cern.ch/user/p/ppillot/pp7TeV/LHC10c/pass2_GMS");
  // Set data search pattern
//  plugin->SetDataPattern("pass1/*ESD.tag.root");
  plugin->SetDataPattern("*ESD.tag.root");
  // ...then add run numbers to be considered
  ifstream inFile(runListName.Data());
  TString currRun;
  if (inFile.is_open())
  {
    while (! inFile.eof() )
    {
      currRun.ReadLine(inFile,kTRUE); // Read line
      if(currRun.IsNull()) continue;
      plugin->AddRunNumber(Form("%09d", currRun.Atoi()));
    }
  }
  inFile.close();
  
  // Define alien work directory where all files will be copied. Relative to alien $HOME.
  plugin->SetGridWorkingDir("pp7TeV/LHC10c/pass2_GMS/MuonResolution");
  
  // Declare alien output directory. Relative to working directory.
  plugin->SetGridOutputDir("results");
  
  // Set the ouput directory of each masterjob to the run number
  plugin->SetOutputToRunNo();
  
  // Declare the analysis source files names separated by blancs. To be compiled runtime
  // using ACLiC on the worker nodes.
  if (!taskInPWG3) plugin->SetAnalysisSource("AliAnalysisTaskMuonResolution.cxx");
  
  // Declare all libraries (other than the default ones for the framework. These will be
  // loaded by the generated analysis macro. Add all extra files (task .cxx/.h) here.
  if (taskInPWG3) plugin->SetAdditionalLibs("libRAWDatabase.so libCDB.so libSTEER.so libMUONcore.so libMUONmapping.so libMUONcalib.so libMUONgeometry.so libMUONtrigger.so libMUONraw.so libMUONbase.so libMUONrec.so libPWG3base.so libPWG3muondep.so");
  else plugin->SetAdditionalLibs("libRAWDatabase.so libCDB.so libSTEER.so libMUONcore.so libMUONmapping.so libMUONcalib.so libMUONgeometry.so libMUONtrigger.so libMUONraw.so libMUONbase.so libMUONrec.so AliAnalysisTaskMuonResolution.h AliAnalysisTaskMuonResolution.cxx");
  plugin->SetAdditionalRootLibs("libGui.so libProofPlayer.so libXMLParser.so");
  
  // Optionally add include paths
  plugin->AddIncludePath("-I$ALICE_ROOT/MUON");
  plugin->AddIncludePath("-I$ALICE_ROOT/MUON/mapping");
  plugin->AddIncludePath("-I.");
  
  // Optionally set a name for the generated analysis macro (default MyAnalysis.C)
  plugin->SetAnalysisMacro("MuonResolution.C");
  plugin->SetExecutable("MuonResolution.sh");
  
  // Optionally set maximum number of input files/subjob (default 100, put 0 to ignore)
  //plugin->SetSplitMaxInputFileNumber(0);
  
  // Optionally set number of failed jobs that will trigger killing waiting sub-jobs.
  plugin->SetMaxInitFailed(5);
  
  // Optionally resubmit threshold.
  //plugin->SetMasterResubmitThreshold(90);
  
  // Optionally set time to live (default 30000 sec)
  plugin->SetTTL(30000);
  
  // Optionally set input format (default xml-single)
  plugin->SetInputFormat("xml-single");
  
  // Optionally modify the name of the generated JDL (default analysis.jdl)
  plugin->SetJDLName("MuonResolution.jdl");
  
  // Optionally modify job price (default 1)
  plugin->SetPrice(1);
  
  // Optionally modify split mode (default 'se')
  plugin->SetSplitMode("se");
  
  return plugin;
}

