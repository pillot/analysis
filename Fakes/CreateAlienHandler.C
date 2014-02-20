/*
 *  CreateAlienHandler.C
 *  aliroot
 *
 *  Created by Philippe Pillot on 01/04/10.
 *  Copyright 2010 SUBATECH. All rights reserved.
 *
 */

AliAnalysisGrid* CreateAlienHandler(TString runMode, TString runListName)
{
  
  AliAnalysisAlien *plugin = new AliAnalysisAlien();
  
  // Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
  plugin->SetRunMode(runMode.Data());  // VERY IMPORTANT - DECRIBED BELOW
  
  // Set versions of used packages
  plugin->SetAPIVersion("V1.1x");
  plugin->SetROOTVersion("v5-27-06b");
  plugin->SetAliROOTVersion("v4-21-09-AN");
  
  // Declare input data to be processed.
  // Method 1: Create automatically XML collections using alien 'find' command.
  // Define production directory LFN
  plugin->SetGridDataDir("/alice/cern.ch/user/b/bianchil/RealisticSignalProductionLHC10/Sept2010/out/");
  // Set data search pattern
  plugin->SetDataPattern("*AliESDs.root");
//  plugin->SetDataPattern("*ESD.tag.root");
  // ...then add run numbers to be considered
  ifstream inFile(runListName.Data());
  TString currRun;
  if (inFile.is_open())
  {
    while (! inFile.eof() )
    {
      currRun.ReadLine(inFile,kTRUE); // Read line
      if(currRun.IsNull()) continue;
      plugin->AddRunNumber(Form("%d", currRun.Atoi()));
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
  plugin->SetGridWorkingDir("Sim/JPsi_Livio_LHC10d2b_Nov2010/Fakes2");
  
  // Declare alien output directory. Relative to working directory.
  plugin->SetGridOutputDir("results");
  
  // Set the ouput directory of each masterjob to the run number
  plugin->SetOutputToRunNo();
  
  // Declare the analysis source files names separated by blancs. To be compiled runtime
  // using ACLiC on the worker nodes.
  plugin->SetAnalysisSource("AliAnalysisTaskMuonFakes.cxx");
  
  // Declare all libraries (other than the default ones for the framework. These will be
  // loaded by the generated analysis macro. Add all extra files (task .cxx/.h) here.
  plugin->SetAdditionalLibs("libRAWDatabase.so libRAWDatarec.so libCDB.so libSTEER.so libMUONcore.so libMUONmapping.so libMUONcalib.so libMUONgeometry.so libMUONtrigger.so libMUONraw.so libMUONbase.so libMUONshuttle.so libMUONrec.so libRAWDatasim.so libMUONsim.so libMUONevaluation.so libPWG3base.so AliAnalysisTaskMuonFakes.h AliAnalysisTaskMuonFakes.cxx");
  plugin->SetAdditionalRootLibs("libGui.so libProofPlayer.so libXMLParser.so");
  
  // Optionally add include paths
  plugin->AddIncludePath("-I$ALICE_ROOT/MUON");
  plugin->AddIncludePath("-I$ALICE_ROOT/MUON/mapping");
  plugin->AddIncludePath("-I$ALICE_ROOT/PWG3/base");
  plugin->AddIncludePath("-I.");
  
  // Optionally set a name for the generated analysis macro (default MyAnalysis.C)
  plugin->SetAnalysisMacro("MuonFakesAnalysis.C");
  plugin->SetExecutable("MuonFakesAnalysis.sh");
  
  // Optionally set maximum number of input files/subjob (default 100, put 0 to ignore)
  //plugin->SetSplitMaxInputFileNumber(0);
  
  // Optionally set number of failed jobs that will trigger killing waiting sub-jobs.
  //plugin->SetMaxInitFailed(5);
  
  // Optionally resubmit threshold.
  //plugin->SetMasterResubmitThreshold(90);
  
  // Optionally set time to live (default 30000 sec)
  plugin->SetTTL(30000);
  
  // Optionally set input format (default xml-single)
  plugin->SetInputFormat("xml-single");
  
  // Optionally modify the name of the generated JDL (default analysis.jdl)
  plugin->SetJDLName("MuonFakesAnalysis.jdl");
  
  // Optionally modify job price (default 1)
  plugin->SetPrice(1);
  
  // Optionally modify split mode (default 'se')
  plugin->SetSplitMode("se");
  
  // Merge via JDL
  plugin->SetMergeViaJDL(kTRUE);
  
  // Optionally set preferred SE
  //plugin->SetPreferedSE("ALICE::Subatech::SE");
  
  return plugin;
}

