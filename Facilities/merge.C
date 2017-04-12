/*
 *  merge.C
 *  aliroot_dev
 *
 *  Created by philippe pillot on 03/06/11.
 *  Copyright 2011 SUBATECH. All rights reserved.
 *
 */

TString home = "/alice/cern.ch/user/p/ppillot";

void merge(TString gridLocation, TString fileName = "AnalysisResults.root", TString runList = "runList.txt")
{
  /// merge the files named "fileName" located in "alien://home/gridLocation/*/run" for each run in the "runList"
  
  // load potentially needed libraries
  gROOT->LoadMacro("$HOME/Work/Alice/Macros/Facilities/runTaskFacilities.C");
  TString extraLibs = "";
  TString extraIncs = "";
  TString extraTasks = "";
  TString extraPkgs="";
//  TString extraLibs = "CORRFW:PWGmuon";
//  TString extraIncs = "include:PWG/muon";
//  TString extraTasks = "AliAnalysisTaskJPsi";
  LoadAlirootLocally(extraLibs, extraIncs, extraTasks, extraPkgs);
  
  // connect to alien
  if (!TGrid::Connect("alien://")) {
    Error("merge","cannot connect to grid");
    return;
  }
  
  // open the run list
  ifstream inFile(runList.Data());
  if (!inFile.is_open()) {
    Error("merge",Form("unable to open file %s", runList.Data()));
    return;
  }
  
  // create TFileMerger
  TFileMerger *fm = new TFileMerger(kFALSE);
  fm->SetFastMethod(kTRUE);
  
  // loop over runs
  while (!inFile.eof()) {
    
    // get the current run number
    TString currRun;
    currRun.ReadLine(inFile, kTRUE);
    if (currRun.IsNull()) continue;
    if (!currRun.IsDigit()) {
      Error("merge","invalid run number: %s", currRun.Data());
      delete fm;
      return;
    }
    
    // get the corresponding file to merge
    TString command = Form("find %s/%s/ *%s/%s", home.Data(), gridLocation.Data(), currRun.Data(), fileName.Data());
    TGridResult *res = gGrid->Command(command);
    if (!res || res->GetSize() == 0) {
      Error("merge",Form("no result for the command: %s", command.Data()));
      continue;
    }
    
    // should be only 1 file to merge per run
    if (res->GetSize() != 1) {
      Error("merge",Form("%d results for the command: %s", res->GetSize(), command.Data()));
      delete res;
      continue;
    }
    
    // get the turl of the current file
    TMap *map = static_cast<TMap*>(res->First());
    TObjString *objs = dynamic_cast<TObjString*>(map->GetValue("turl"));
    if (!objs || !objs->GetString().Length()) {
      Error("merge","turl not found for the run %s... skipping", currRun.Data());
      delete res;
      continue;
    }
    
    // add the file for merging
    fm->AddFile(objs->GetName());
    
    // clean memory
    delete res;
    
  }
  
  // close the run list
  inFile.close();
  
  // nothing found - skip this output
  if (!fm->GetMergeList() || !fm->GetMergeList()->GetSize()) {
    Error("merge",Form("no <%s> files found", fileName.Data()));
    delete fm;
    return;
  }
  
  // open output file
  fm->OutputFile(fileName);
  
  // merge
  printf("merging %d files...\n", fm->GetMergeList()->GetSize());
  if (!fm->Merge()) {
    Error("merge", "could not merge all <%s> files", fileName.Data());
    delete fm;
    return;
  }
  
  // clean memory
  delete fm;
  
}

