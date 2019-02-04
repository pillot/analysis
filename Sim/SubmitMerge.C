/*
 *  Submit.C
 *  aliroot
 *
 *  Created by Philippe Pillot on 09/03/12.
 *  Copyright 2010 SUBATECH. All rights reserved.
 *
 */

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TFile.h>
#include <Riostream.h>
#include <TSystem.h>
#include <TString.h>
#include <TGrid.h>
#include <TGridResult.h>
#include <TMap.h>
#include <TObjString.h>
#endif

Int_t splitLevel = 20;

Bool_t FileExists(const char *lfn);
Bool_t DirectoryExists(const char *dirname);

//______________________________________________________________________________
void SubmitMerge(const char* inDir, const char* outDir, const char* resDir = "results",
                 const char* runList = "runList.txt", TString runFormat = "%d", Int_t stage = 0, Bool_t submit = kFALSE)
{
  /// Submit multiple merging jobs with the format "submit AOD_merge(_final).jdl run# (stage#)".
  /// Also produce the xml collection before sending jobs
  /// Example:
  /// - inDir = "/alice/sim/2012/LHC12a10_bis" (where to find the data to merge)
  ///         = 0x0 --> inDir = homeDir/outDir/resDir
  /// - outDir = "Sim/LHC11h/embedding/AODs" (where to store merged results)
  /// - runList.txt must contains the list of run number
  /// - stage=0 --> final merging / stage>0 --> intermediate merging i
  ///
  /// *** THIS MACRO MUST BE COMPILED <-- CINT PROBLEMS... ***
  
  TString homeDir = "/alice/cern.ch/user/p/ppillot";
  
  if (!TGrid::Connect("alien://")) {
    printf("cannot connect to grid\n");
    return;
  }
  
  TString outputDir = Form("%s/%s", homeDir.Data(), outDir);
  TString inputDir = inDir ? inDir : Form("%s/%s", outputDir.Data(),resDir);
  
  if (!DirectoryExists(outputDir.Data())) {
    printf("directory %s does not exist\n", outputDir.Data());
    return;
  }
  
  gGrid->Cd(outputDir.Data());
  
  //TString jdl = (stage > 0) ? "AOD_merge.jdl" : "AOD_merge_final.jdl";
  TString jdl = (stage > 0) ? "merge.jdl" : "merge_final.jdl";
  if (!FileExists(jdl.Data())) {
    printf("file %s does not exist in %s\n", jdl.Data(), outputDir.Data());
    return;
  }
  
  ifstream inFile(runList);
  if (!inFile.is_open()) {
    printf("cannot open file %s\n", runList);
    return;
  }
  
  TString currRun;
  TString reply = "";
  gSystem->Exec("rm -f __failed__");
  Bool_t failedRun = kFALSE;
  while (!inFile.eof())
  {
    currRun.ReadLine(inFile,kTRUE);
    if(currRun.IsNull()) continue;
    
    Int_t run = currRun.Atoi();
    TString srun = Form(runFormat.Data(), run);
    printf("\n --- processing run %s ---\n", srun.Data());
    
    TString runDir = Form("%s/%s/%s", outputDir.Data(), resDir, srun.Data());
    if (!DirectoryExists(runDir.Data())) {
      printf(" - creating output directory %s\n", runDir.Data());
      gSystem->Exec(Form("alien_mkdir -p %s", runDir.Data()));
    }
    
    if (FileExists(Form("%s/root_archive.zip", runDir.Data()))) {
      printf(" ! final merging already done\n");
      continue;
    }
    
    Int_t n = 0, lastStage = 0;
    gSystem->Exec(Form("alien_ls -F %s | grep Stage_.*/ > __stage__", runDir.Data()));
    ifstream f("__stage__");
    std::string dummy;
    while (std::getline(f, dummy)) n++;
    f.close();
    while (n > 0) if (gSystem->Exec(Form("grep Stage_%d/ __stage__ 2>&1 >/dev/null", ++lastStage)) == 0) n--;
    gSystem->Exec("rm -f __stage__");
    if (stage > 0 && stage != lastStage+1) {
      printf(" ! lastest merging stage = %d. Next must be stage %d or final stage\n", lastStage, lastStage+1);
      continue;
    }
    
//    TString wn = (stage > 0) ? Form("Stage_%d.xml", stage) : "wn.xml";
    TString wn = (stage > 0) ? Form("Stage_%d.xml", stage) : "Stage_final.xml";
    TString find = (lastStage == 0) ?
      Form("alien_find -x %s %s/%s *root_archive.zip", wn.Data(), inputDir.Data(), srun.Data()) :
      Form("alien_find -x %s %s/%s/Stage_%d *root_archive.zip", wn.Data(), inputDir.Data(), srun.Data(), lastStage);
//    TString find = (lastStage == 0) ?
//      Form("alien_find -x %s %s/%s/ESDs/muon_calo_pass2 12%s*/root_archive.zip", wn.Data(), inputDir.Data(), srun.Data(), srun.Data()) :
//      Form("alien_find -x %s %s/%s/vpass1 12%s*/root_archive.zip", wn.Data(), inputDir.Data(), srun.Data(), srun.Data()) :
//    TString find = (lastStage == 0) ?
//      Form("alien_find -x %s %s/%s *Merged.QA.Data.root", wn.Data(), inputDir.Data(), srun.Data()) :
//      Form("alien_find -x %s %s/%s/Stage_%d *Merged.QA.Data.root", wn.Data(), inputDir.Data(), srun.Data(), lastStage);
    gSystem->Exec(Form("%s 1> %s 2>/dev/null", find.Data(), wn.Data()));
    gSystem->Exec(Form("grep -c /event %s > __nfiles__", wn.Data()));
    ifstream f2("__nfiles__");
    TString nFiles;
    nFiles.ReadLine(f2,kTRUE);
    f2.close();
    gSystem->Exec("rm -f __nfiles__");
    printf(" - number of files to merge = %d\n", nFiles.Atoi());
    if (nFiles.Atoi() == 0) {
      printf(" ! collection of files to merge is empty\n");
      gSystem->Exec(Form("rm -f %s", wn.Data()));
      continue;
    } else if (stage > 0 && nFiles.Atoi() <= splitLevel && !reply.BeginsWith("y")) {
      if (!reply.BeginsWith("n")) {
	printf(" ! number of files to merge <= split level (%d). Continue? [Y/n] ", splitLevel);
	fflush(stdout);
	reply.Gets(stdin,kTRUE);
	reply.ToLower();
      }
      if (reply.BeginsWith("n")) {
	gSystem->Exec(Form("rm -f %s", wn.Data()));
	continue;
      } else reply = "y";
    }
    
    if (submit) {
      TString dirwn = Form("%s/%s", runDir.Data(), wn.Data());
      if (FileExists(dirwn.Data())) gGrid->Rm(dirwn.Data());
      gSystem->Exec(Form("alien_cp file:%s alien://%s", wn.Data(), dirwn.Data()));
      gSystem->Exec(Form("rm -f %s", wn.Data()));
    }
    
    TString query;
    if (stage > 0) query = Form("submit %s %s %d", jdl.Data(), srun.Data(), stage);
    else query = Form("submit %s %s", jdl.Data(), srun.Data());
    printf(" - %s ...", query.Data());
    fflush(stdout);
    
    if (!submit) {
      printf(" dry run\n");
      continue;
    }
    
    Bool_t done = kFALSE;
    TGridResult *res = gGrid->Command(query);
    if (res) {
      TString cjobId1 = res->GetKey(0,"jobId");
      if (!cjobId1.IsDec()) {
	printf(" FAILED\n");
	gGrid->Stdout();
	gGrid->Stderr();
      } else {
	printf(" DONE\n   --> the job Id is: %s \n", cjobId1.Data());
	done = kTRUE;
      }
      delete res;
    } else printf(" FAILED\n");
    
    if (!done) {
      gSystem->Exec(Form("echo %s >> __failed__", srun.Data()));
      failedRun = kTRUE;
    }
    
  }
  printf("\n");

  inFile.close();
  
  if (failedRun) {
    printf("\n--------------------\n");
    printf("list of failed runs:\n");
    gSystem->Exec("cat __failed__");
    gSystem->Exec("rm -f __failed__");
    printf("\n");
  }
  
}

//______________________________________________________________________________
Bool_t FileExists(const char *lfn)
{
  // Returns true if file exists.
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

//______________________________________________________________________________
Bool_t DirectoryExists(const char *dirname)
{
  // Returns true if directory exists. Can be also a path.
  if (!gGrid) return kFALSE;
  // Check if dirname is a path
  TString dirstripped = dirname;
  dirstripped = dirstripped.Strip();
  dirstripped = dirstripped.Strip(TString::kTrailing, '/');
  TString dir = gSystem->BaseName(dirstripped);
  dir += "/";
  TString path = gSystem->DirName(dirstripped);
  TGridResult *res = gGrid->Ls(path, "-F");
  if (!res) return kFALSE;
  TIter next(res);
  TMap *map;
  TObject *obj;
  while ((map=dynamic_cast<TMap*>(next()))) {
    obj = map->GetValue("name");
    if (!obj) break;
    if (dir == obj->GetName()) {
      delete res;
      return kTRUE;
    }
  }
  delete res;
  return kFALSE;
}      

