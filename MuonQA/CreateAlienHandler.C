/*
 *  CreateAlienHandler.C
 *  aliroot
 *
 *  Created by Philippe Pillot on 21/06/10.
 *  Copyright 2010 Subatech. All rights reserved.
 *
 */

AliAnalysisGrid* CreateAlienHandler(TString runMode, TString runListName, Bool_t taskInPWG3)
{
  
  AliAnalysisAlien *plugin = new AliAnalysisAlien();
  
  // Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
  plugin->SetRunMode(runMode.Data());  // VERY IMPORTANT - DECRIBED BELOW
  
  // Set versions of used packages
  plugin->SetAPIVersion("V1.1x");
  plugin->SetROOTVersion("v5-27-06b");
  plugin->SetAliROOTVersion("v4-21-14-AN");
  
  // Declare input data to be processed:
  // Define production directory LFN
  plugin->SetGridDataDir("/alice/data/2010/LHC10h");
  //plugin->SetGridDataDir("/alice/cern.ch/user/p/ppillot/PbPb2.76TeV/LHC10h/pass1_4plus");
  // Set data search pattern
  plugin->SetDataPattern("pass1_5plus/*AliESDs.root");
  //plugin->SetDataPattern("*AliESDs.root");
  // ...then add run numbers to be considered
  ifstream inFile(runListName.Data());
  TString currRun;
  if (inFile.is_open())
  {
    while (! inFile.eof() )
    {
      currRun.ReadLine(inFile,kTRUE); // Read line
      if(currRun.IsNull()) continue;
      //plugin->AddRunNumber(Form("%d", currRun.Atoi()));
      plugin->AddRunNumber(Form("%09d", currRun.Atoi()));
    }
  }
  inFile.close();
  
  // Define alien work directory where all files will be copied. Relative to alien $HOME.
  plugin->SetGridWorkingDir("PbPb2.76TeV/LHC10h/pass1_5plus/MuonQA_selected");
  
  // Declare alien output directory. Relative to working directory.
  plugin->SetGridOutputDir("results");
  
  // Set the ouput directory of each masterjob to the run number
  plugin->SetOutputToRunNo();
  
  // Declare the analysis source files names separated by blancs. To be compiled runtime
  // using ACLiC on the worker nodes.
  if (!taskInPWG3) plugin->SetAnalysisSource("AliCounterCollection.cxx AliAnalysisTaskMuonQA.cxx");
  
  // enable specific packages
  //plugin->EnablePackage("PWG3muon.par");
  
  // Declare all libraries (other than the default ones for the framework. These will be
  // loaded by the generated analysis macro. Add all extra files (task .cxx/.h) here.
  //if (taskInPWG3) plugin->SetAdditionalLibs("libPWG3base.so");
  if (taskInPWG3) plugin->SetAdditionalLibs("libPWG3base.so libPWG3muon.so");
  else plugin->SetAdditionalLibs("AliAnalysisTaskMuonQA.h AliAnalysisTaskMuonQA.cxx AliCounterCollection.h AliCounterCollection.cxx");
  
  // Optionally add include paths
  plugin->AddIncludePath("-I.");
  
  // Optionally set a name for the generated analysis macro (default MyAnalysis.C)
  plugin->SetAnalysisMacro("MuonQAAnalysis.C");
  plugin->SetExecutable("MuonQAAnalysis.sh");
  
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
  plugin->SetJDLName("MuonQAAnalysis.jdl");
  
  // Optionally modify job price (default 1)
  plugin->SetPrice(1);
  
  // Optionally modify split mode (default 'se')
  plugin->SetSplitMode("se");
  
  // Merge via JDL
  plugin->SetMergeViaJDL(kTRUE);
  
  // Optionally set the maximum number of files merged together in one stage
  plugin->SetMaxMergeFiles(10);
  
  // Optionally set the maximum number of merging stages
  plugin->SetMaxMergeStages(3);
  
  // List of files not to be merged
  //plugin->SetMergeExcludes("EventStat_temp.root");
  
  // Do not overwrite existing collections and merged outputs
  //plugin->SetOverwriteMode(kFALSE);
  
  return plugin;
}

