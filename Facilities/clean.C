/*
 *  clean.C
 *  aliroot_dev
 *
 *  Created by philippe pillot on 16/08/11.
 *  Copyright 2011 SUBATECH. All rights reserved.
 *
 */

void clean(TString gridLocation, TString runList, TString fileName)
{
  /// remove file "fileName" in every "home/gridLocation/run/*" directories
  
  TString home = "/alice/cern.ch/user/p/ppillot";
  
  // connect to alien
  if (!TGrid::Connect("alien://")) {
    Error("clean","cannot connect to grid");
    return;
  }
  
  // open the run list
  ifstream inFile(runList.Data());
  if (!inFile.is_open()) {
    Error("clean",Form("unable to open file %s", runList.Data()));
    return;
  }
  
  // loop over runs
  while (!inFile.eof()) {
    
    // get the current run number
    TString currRun;
    currRun.ReadLine(inFile, kTRUE);
    if (currRun.IsNull()) continue;
    if (!currRun.IsDigit()) {
      Error("clean","invalid run number: %s", currRun.Data());
      return;
    }
    
    // get the corresponding list of files to remove
    TString command = Form("find %s/%s/%s */%s", home.Data(), gridLocation.Data(), currRun.Data(), fileName.Data());
    TGridResult *res = gGrid->Command(command);
    if (!res || res->GetSize() == 0) {
      Error("clean",Form("no result for the command: %s", command.Data()));
      continue;
    }
    
    // loop over the list of files
    TIter next(res);
    TMap *map = 0x0;
    while ((map = static_cast<TMap*>(next()))) {
      
      // get the turl of the current file
      TObjString *objs = dynamic_cast<TObjString*>(map->GetValue("turl"));
      if (!objs || !objs->GetString().Length()) {
	Error("clean","turl not found for the run %s... SKIPPING", currRun.Data());
	delete res;
	continue;
      }
      
      // remove the file
      TString file = objs->GetName();
      file.ReplaceAll("alien://","");
      printf("removing file %s...\n", file.Data());
      //gSystem->Exec(Form("alien_rm %s", file.Data()));
      
    }
    
    // clean memory
    delete res;
    
  }
  
  // close the run list
  inFile.close();
  
}

