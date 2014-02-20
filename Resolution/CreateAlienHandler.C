/*
 *  CreateAlienHandler.C
 *  aliroot
 *
 *  Created by Philippe Pillot on 05/04/10.
 *  Copyright 2010 SUBATECH. All rights reserved.
 *
 */

AliAnalysisGrid* CreateAlienHandler(TString runMode, TString runListName)
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
  plugin->SetROOTVersion("v5-26-00b-5");
  plugin->SetAliROOTVersion("v4-19-13-AN");
  
  // Set the number of files used in test mode
  //plugin->SetNtestFiles(1000000);
  
  // Set the number of runs per master job
  //plugin->SetNrunsPerMaster(1000000);
  
  // Declare input data to be processed.
  // Method 1: Create automatically XML collections using alien 'find' command.
  // Define production directory LFN
//  plugin->SetGridDataDir("/alice/data/2009/LHC09c/");
//  plugin->SetGridDataDir("/alice/data/2010/LHC10c/");
  plugin->SetGridDataDir("/alice/data/2010/LHC10d/");
  // Set data search pattern
  plugin->SetDataPattern("pass1/*ESD.tag.root");
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
  // Method 2: Declare existing data files (raw collections, xml collections, root file)
  // If no path mentioned data is supposed to be in the work directory (see SetGridWorkingDir())
  // XML collections added via this method can be combined with the first method if
  // the content is compatible (using or not tags)
  //   plugin->AddDataFile("tag.xml");
  //   plugin->AddDataFile("/alice/data/2008/LHC08c/000057657/raw/Run57657.Merged.RAW.tag.root");
  
  // Define alien work directory where all files will be copied. Relative to alien $HOME.
//  plugin->SetGridWorkingDir("pp7TeV/LHC10c/Pos/Resolution");
  plugin->SetGridWorkingDir("pp7TeV/LHC10d/Resolution");
  
  // Declare alien output directory. Relative to working directory.
  plugin->SetGridOutputDir("results");
  
  // Set the ouput directory of each masterjob to the run number
  plugin->SetOutputToRunNo();
  
  // Declare the analysis source files names separated by blancs. To be compiled runtime
  // using ACLiC on the worker nodes.
  plugin->SetAnalysisSource("AliAnalysisTaskMuonChamberResolution.cxx");
  
  // Declare all libraries (other than the default ones for the framework. These will be
  // loaded by the generated analysis macro. Add all extra files (task .cxx/.h) here.
  plugin->SetAdditionalLibs("libRAWDatabase.so libCDB.so libSTEER.so libMUONcore.so libMUONmapping.so libMUONcalib.so libMUONgeometry.so libMUONtrigger.so libMUONraw.so libMUONbase.so libMUONrec.so AliAnalysisTaskMuonChamberResolution.h AliAnalysisTaskMuonChamberResolution.cxx");
  plugin->SetAdditionalRootLibs("libGui.so libProofPlayer.so libXMLParser.so");
  
  // Optionally add include paths
  plugin->AddIncludePath("-I$ALICE_ROOT/MUON");
  plugin->AddIncludePath("-I$ALICE_ROOT/MUON/mapping");
  plugin->AddIncludePath("-I.");
  
  // Declare the output file names separated by blancs.
  // (can be like: file.root or file.root@ALICE::Niham::File)
//  plugin->SetOutputFiles("ESDCheck.root EventStat_temp.root");
  plugin->SetOutputFiles("chamberResolution.root");
  
  // Optionally define the files to be archived.
  //   plugin->SetOutputArchive("log_archive.zip:stdout,stderr@ALICE::NIHAM::File root_archive.zip:*.root@ALICE::NIHAM::File");
  plugin->SetOutputArchive("log_archive.zip:stdout,stderr");
  //plugin->SetOutputArchive();
  
  // Optionally set a name for the generated analysis macro (default MyAnalysis.C)
  plugin->SetAnalysisMacro("MuonChamberResolution.C");
  plugin->SetExecutable("MuonChamberResolution.sh");
  
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
  plugin->SetJDLName("MuonChamberResolution.jdl");
  
  // Optionally modify job price (default 1)
  plugin->SetPrice(1);
  
  // Optionally modify split mode (default 'se')
  plugin->SetSplitMode("se");
  
  // Optionally set preferred SE
  plugin->SetPreferedSE("ALICE::Subatech::SE,disk=2");
  
  return plugin;
}

