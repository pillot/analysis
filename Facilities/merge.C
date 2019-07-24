/*
 *  merge.C
 *  aliroot_dev
 *
 *  Created by philippe pillot on 03/06/11.
 *  Copyright 2011 SUBATECH. All rights reserved.
 *
 */

#if !defined(__CINT__) || defined(__MAKECINT__)

#include <Riostream.h>

// ROOT includes
#include "TError.h"
#include "TROOT.h"
#include "TString.h"
#include "TObjString.h"
#include "TGrid.h"
#include "TGridResult.h"
#include "TMap.h"
#include "TFileMerger.h"

#endif

TString home = "/alice/cern.ch/user/p/ppillot";

void AddFiles(TString gridLocation, TString run, TString fileName, TFileMerger *fm);
Bool_t FileExists(const char *lfn);

//---------------------------------------------------------------------------------
void merge(TString gridLocation, TString fileName = "AnalysisResults.root", TString runList = "runList.txt")
{
  /// merge the files named "fileName" located in "alien://home/gridLocation/run" for each run in the "runList"
  /// if a file is found in the directory run/. then take this one, otherwise take all files in run/*
  /// if runList is empty, then merge the files in gridLocation/
  ///
  /// The same can be achieved by doing:
  /// TGrid::Connect("alien://")
  /// AliAnalysisAlien::MergeOutput("fileName","home/gridLocation",10,0)
  /// if the last argument is 0, the merging is done in several steps (merging 10 files at a time here)
  /// if the last argument is >0, the merging is done in one go
  /// works also with a .xml collection or a .txt file containing alien:///home/gridLocation directories

  // load potentially needed libraries
  gROOT->LoadMacro("$HOME/Work/Alice/Macros/Facilities/runTaskFacilities.C");
	TString extraLibs="";
	TString extraPkgs="";
  gROOT->ProcessLineFast(Form("LoadAlirootLocally(\"%s\", \"\", \"\", \"%s\")",extraLibs.Data(),extraPkgs.Data()));

  // connect to alien
  if (!TGrid::Connect("alien://")) {
    Error("merge","cannot connect to grid");
    return;
  }

  // create TFileMerger
  TFileMerger *fm = new TFileMerger(kFALSE);
  fm->SetFastMethod(kTRUE);

  if (runList.IsNull()) {

    // add the file(s) in gridLocation/ if runList is empty
    AddFiles(gridLocation, "", fileName, fm);

  } else {

    // ... or runList is a text file of run numbers so open it
    ifstream inFile(runList.Data());
    if (!inFile.is_open()) {
      Error("merge", "unable to open file %s", runList.Data());
      return;
    }

    // loop over runs and add the corresponding file(s)
    TString currRun;
    while (!inFile.eof()) {
      currRun.ReadLine(inFile, kTRUE);
      if (currRun.IsNull()) continue;
      if (!currRun.IsDigit()) {
        Error("merge", "invalid run number: %s", currRun.Data());
        continue;
      }
      AddFiles(gridLocation, currRun, fileName, fm);
    }

    inFile.close();
  }

  // nothing found - skip this output
  if (!fm->GetMergeList() || !fm->GetMergeList()->GetSize()) {
    Error("merge", "no <%s> files found", fileName.Data());
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

//---------------------------------------------------------------------------------
void AddFiles(TString gridLocation, TString run, TString fileName, TFileMerger *fm)
{
  /// add file(s) in home/gridLocation/run/(*)

  // if there is a file "fileName" in gridLocation/run/. then add it
  TString file = Form("%s/%s/%s/%s", home.Data(), gridLocation.Data(), run.Data(), fileName.Data());
  if (FileExists(file.Data())) {
    fm->AddFile(TString::Format("alien://%s", file.Data()));
    return;
  }

  // otherwise add all files "fileName" found in gridLocation/run/*
  TString command = Form("find %s/%s/%s */%s", home.Data(), gridLocation.Data(), run.Data(), fileName.Data());
  TGridResult *res = gGrid->Command(command);
  if (!res || res->GetSize() == 0) {
    Error("merge", "no file to merge in %s/%s", gridLocation.Data(), run.Data());
    delete res;
    return;
  }

  TIter resIt(res);
  TMap *map = nullptr;
  while ((map = static_cast<TMap*>(resIt.Next()))) {

    // get the turl of the current file
    TObjString *objs = dynamic_cast<TObjString*>(map->GetValue("turl"));
    if (!objs || !objs->GetString().Length()) {
      Error("merge","turl not found in %s/%s... skipping", gridLocation.Data(), run.Data());
      continue;
    }

    // add the file for merging
    fm->AddFile(objs->GetName());
  }

  delete res;
}

//---------------------------------------------------------------------------------
Bool_t FileExists(const char *lfn)
{
  /// Returns true if file exists.
  if (!gGrid) return kFALSE;
  TGridResult *res = gGrid->Ls(lfn);
  if (!res) return kFALSE;
  TMap *map = dynamic_cast<TMap*>(res->At(0));
  if (!map) {
    delete res;
    return kFALSE;
  }   
  TObjString *objs = dynamic_cast<TObjString*>(map->GetValue("name"));
  if (!objs || !objs->GetString().Length()) {
    delete res;
    return kFALSE;
  }
  delete res;   
  return kTRUE;
}
