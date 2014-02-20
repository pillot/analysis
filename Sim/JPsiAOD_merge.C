/*
 *  JPsiAOD_merge.C
 *  aliroot_dev
 *
 *  Created by philippe pillot on 22/08/11.
 *  Copyright 2011 SUBATECH. All rights reserved.
 *
 */


void JPsiAOD_merge(const char *dir, Int_t stage=0)
{
  // Macro to replace the one produced by the plugin
  // Simply change the list of outputFiles to merge to "AliAOD.root"
  
  TStopwatch timer;
  timer.Start();
  
  // Reset existing include path and add current directory first in the search
  gSystem->SetIncludePath("-I.");
  // Load analysis framework libraries
  gSystem->Load("libANALYSIS");
  gSystem->Load("libOADB");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCORRFW");
  
  // include path
  TString intPath = gInterpreter->GetIncludePath();
  TObjArray *listpaths = intPath.Tokenize(" ");
  TIter nextpath(listpaths);
  TObjString *pname;
  while ((pname=(TObjString*)nextpath())) {
    TString current = pname->GetName();
    if (current.Contains("AliRoot") || current.Contains("ALICE_ROOT")) continue;
    gSystem->AddIncludePath(current);
  }
  if (listpaths) delete listpaths;
  gROOT->ProcessLine(".include $ALICE_ROOT/include");
  printf("Include path: %s\n", gSystem->GetIncludePath());
  
  // Add aditional AliRoot libraries
  gSystem->Load("libPWG3base.so");
  gSystem->Load("libPWG3muon.so");
  
  // Analysis source to be compiled at runtime (if any)
  
  // Set temporary merging directory to current one
  gSystem->Setenv("TMPDIR", gSystem->pwd());
  
  // Set temporary compilation directory to current one
  gSystem->SetBuildDir(gSystem->pwd(), kTRUE);
  
  // Connect to AliEn
  if (!TGrid::Connect("alien://")) return;
  TString outputDir = dir;
  TString outputFiles = "AliAOD.root";
  TString mergeExcludes = "";
  TObjArray *list = outputFiles.Tokenize(",");
  TIter *iter = new TIter(list);
  TObjString *str;
  TString outputFile;
  Bool_t merged = kTRUE;
  while((str=(TObjString*)iter->Next())) {
    outputFile = str->GetString();
    if (outputFile.Contains("*")) continue;
    Int_t index = outputFile.Index("@");
    if (index > 0) outputFile.Remove(index);
    // Skip already merged outputs
    if (!gSystem->AccessPathName(outputFile)) {
      printf("Output file <%s> found. Not merging again.",outputFile.Data());
      continue;
    }
    if (mergeExcludes.Contains(outputFile.Data())) continue;
    merged = AliAnalysisAlien::MergeOutput(outputFile, outputDir, 10, stage);
    if (!merged) {
      printf("ERROR: Cannot merge %s\n", outputFile.Data());
      return;
    }
  }
  // all outputs merged, validate
  ofstream out;
  out.open("outputs_valid", ios::out);
  out.close();
  // read the analysis manager from file
  if (!outputDir.Contains("Stage")) return;
  AliAnalysisManager *mgr = AliAnalysisAlien::LoadAnalysisManager("JPsiAOD.root");
  if (!mgr) return;
  mgr->SetRunFromPath(mgr->GetRunFromAlienPath(dir));
  mgr->SetSkipTerminate(kFALSE);
  mgr->PrintStatus();
  AliLog::SetGlobalLogLevel(AliLog::kError);
  TTree *tree = NULL;
  mgr->StartAnalysis("gridterminate", tree);
}

