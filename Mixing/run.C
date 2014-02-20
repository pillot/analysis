/*
 *  run.C
 *  aliroot_dev
 *
 *  Created by philippe pillot on 01/05/11.
 *  Copyright 2011 SUBATECH. All rights reserved.
 *
 */

TString indir = ".";
TString outdir = "./mixed";

void run(TString runList = "runList.txt", Bool_t createPools = kTRUE, Bool_t mixEvents = kTRUE)
{
  /// create Pools and mix events for each run and each centrality bin
  
  const Int_t nCent = 11;
  TString centBinName[nCent] = {"010", "1020", "2040", "4060", "6080", "2030", "3040", "4050", "5060", "6070", "7080"};
  
  if (createPools) gROOT->LoadMacro("$ALICE/Macros/Mixing/CreatePools.C++");
  if (mixEvents) gROOT->LoadMacro("$ALICE/Macros/Mixing/MixEventsCentr.C++");
  
  // open the run list
  ifstream inFile(runList.Data());
  if (!inFile.is_open()) {
    Error("run",Form("unable to open file %s", runList.Data()));
    return;
  }
  
  // loop over runs
  while (!inFile.eof()) {
    
    // get the current run number
    TString currRun;
    currRun.ReadLine(inFile, kTRUE);
    if (currRun.IsNull()) continue;
    if (!currRun.IsDigit()) {
      Error("run","invalid run number: %s", currRun.Data());
      return;
    }
    printf("\n\n-------------------------\n");
    printf("  Processing run %d\n", currRun.Atoi());
    printf("-------------------------\n");
    
    // create outdir if not exist
    TString rundir = Form("%s/%d",outdir.Data(),currRun.Atoi());
    if (gSystem->AccessPathName(rundir.Data())) gSystem->Exec(Form("mkdir -p %s", rundir.Data()));
    
    // create the pools if required
    if (createPools) {
      printf("\n- Creating pools:\n");
      CreatePools(indir.Data(),rundir.Data(),currRun.Atoi());
    }
    
    // mix events if required
    if (mixEvents) {
      for (Int_t i = 0; i < nCent; i++) {
	printf("\n- mix events in centrality bin %s:\n", centBinName[i].Data());
	MixEventsCentr(rundir.Data(),currRun.Atoi(),centBinName[i].Data(),kFALSE);
      }
    }
    
  }
  
  // close the run list
  inFile.close();
  
}

