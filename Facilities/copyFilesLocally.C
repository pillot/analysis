/*
 *  copyFilesLocally.C
 *  aliroot_dev
 *
 *  Created by philippe pillot on 02/04/11.
 *  Copyright 2011 Subatech. All rights reserved.
 *
 */

void copyFilesLocally(TString gridLocation, TString fileName = "AnalysisResults.root", TString runList = "runList.txt")
{
  /// copy the file named "fileName" located in "alien://~/gridLocation/*/#run"
  /// for each #run in the "runList" to the local directory ./runs/#run/
  
//  TString home = "/alice/sim/2012/LHC12a10_bis/167818/p80d";
  TString home = "/alice/cern.ch/user/p/ppillot";
//  TString home = "/alice/cern.ch/user/a/alardeux";
//  TString home = "";
  
  // connect to alien
  if (!TGrid::Connect("alien://")) {
    Error("copyFilesLocally","cannot connect to grid");
    return;
  }
  
  // open the run list
  ifstream inFile(runList.Data());
  if (!inFile.is_open()) {
    Error("copyFilesLocally",Form("unable to open file %s", runList.Data()));
    return;
  }
  
  // loop over runs
  Int_t nMissingFiles = 0;
  Int_t nCopyFailed = 0;
  Bool_t failedRun = kFALSE;
  gSystem->Exec("rm -f __failed__");
  while (!inFile.eof()) {
    
    // get the current run number
    TString currRun;
    currRun.ReadLine(inFile, kTRUE);
    if (currRun.IsNull()) continue;
    if (!currRun.IsDigit()) {
      Error("copyFilesLocally","invalid run number: %s", currRun.Data());
      return;
    }
    
    // copy the file unless it is already there
    TString localDir = Form("runs/%s", currRun.Data());
    if (!gSystem->AccessPathName(Form("%s/%s", localDir.Data(), fileName.Data()))) {
      printf("file %s/%s already exist\n", localDir.Data(), fileName.Data());
      continue;
    }
    
    // get the corresponding file to copy
    TString command = Form("find %s/%s/ *%s/%s", home.Data(), gridLocation.Data(), currRun.Data(), fileName.Data());
//    TString command = Form("find %s/%s/ *%s/ESDs/muon_calo_pass2/%s", home.Data(), gridLocation.Data(), currRun.Data(), fileName.Data());
    TGridResult *res = gGrid->Command(command);
    if (!res || res->GetSize() == 0) {
      Error("copyFilesLocally",Form("no result for the command: %s", command.Data()));
      nMissingFiles++;
      continue;
    }
    
    // should be only 1 file to copy per run
    if (res->GetSize() != 1) {
      Error("copyFilesLocally",Form("%d results for the command: %s", res->GetSize(), command.Data()));
      nMissingFiles++;
      delete res;
      continue;
    }
    
    // get the turl of the current file
    TMap *map = static_cast<TMap*>(res->First());
    TObjString *objs = dynamic_cast<TObjString*>(map->GetValue("turl"));
    if (!objs || !objs->GetString().Length()) {
      Error("copyFilesLocally","turl not found for the run %s... SKIPPING", currRun.Data());
      nMissingFiles++;
      delete res;
      continue;
    }
    
    // create the local directory for that run if needed
    if (gSystem->AccessPathName(localDir.Data()))
      gSystem->Exec(Form("mkdir -p %s", localDir.Data()));
    
    // copy the file
    printf("copying file %s...\n", objs->GetName());
    if (fileName.EndsWith(".root")) {
      if (!TFile::Cp(objs->GetName(), Form("%s/%s",localDir.Data(),fileName.Data())))  {
	gSystem->Exec(Form("echo %s >> __failed__", currRun.Data()));
	nCopyFailed++;
	failedRun = kTRUE;
      }
    } else gSystem->Exec(Form("alien_cp %s file:%s/%s", objs->GetName(), localDir.Data(), fileName.Data()));
    
    // clean memory
    delete res;
    
  }
  
  printf("\n--------------------\n");
  printf("number of missing file = %d\n",nMissingFiles);
  printf("number of failed copy = %d\n",nCopyFailed);
  if (failedRun) {
    printf("list of failed runs:\n");
    gSystem->Exec("cat __failed__");
    gSystem->Exec("rm -f __failed__");
    printf("\n");
  }
  
  // close the run list
  inFile.close();
  
}

