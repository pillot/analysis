/*
 *  Submit.C
 *  aliroot
 *
 *  Created by Philippe Pillot on 17/07/10.
 *  Copyright 2010 SUBATECH. All rights reserved.
 *
 */

Bool_t FileExists(const char *lfn);
Bool_t DirectoryExists(const char *dirname);

//______________________________________________________________________________
void Submit(const char* outputDir, const char* jdl, const char* runList = "runList.txt")
{
  /// Submit multiple jobs with the format "submit jdl 000run#.xml 000run#".
  /// Example:
  /// - outputDir = "/alice/cern.ch/user/p/ppillot/pp7TeV/LHC10d/MuonQA/pass1/results/"
  /// - jdl = "MuonQAAnalysis.jdl"
  
  if (!TGrid::Connect("alien://")) {
    Error("Submit","cannot connect to grid");
    return;
  }
  
  if (!DirectoryExists(outputDir)) {
    Error("Submit", Form("directory %s does not exist", outputDir));
    return;
  }
  
  gGrid->Cd(outputDir);
  
  if (!FileExists(jdl)) {
    Error("Submit", Form("file %s does not exist in %s", jdl, outputDir));
    return;
  }
  
  ifstream inFile(runList);
  if (!inFile.is_open()) {
    Error("Submit", Form("cannot open file %s", runList));
    return;
  }
  
  TString currRun;
  while (! inFile.eof() ) {
    
    currRun.ReadLine(inFile,kTRUE);
    
    if(currRun.IsNull()) continue;
    
    TString query = Form("submit %s %09d.xml %09d", jdl, currRun.Atoi(), currRun.Atoi());
    printf("\n%s ...", query.Data());
    fflush(stdout);
    
    TGridResult *res = gGrid->Command(query);
    
    if (res) {
      
      TString cjobId1 = res->GetKey(0,"jobId");
      
      if (!cjobId1.Length()) {
	printf(" FAILED\n");
	gGrid->Stdout();
	gGrid->Stderr();
      } else printf(" DONE\n   --> the job Id is: %s \n", cjobId1.Data());
      
      delete res;
      
    } else printf(" FAILED\n");
    
  }
  
  inFile.close();
  printf("\n");
  
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

