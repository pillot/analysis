/*
*  cacheDatasets.C
*  aliroot_dev
*
*  Created by philippe pillot on 19/02/15.
*  Copyright 2015 Subatech. All rights reserved.
*
*/

void cacheOneDataset(TString query, TString outFileName);

//-----------------------------------------------------------------------------
void cacheDatasets(TString datasets, TString outFileName = "", Bool_t setRootVersion = kTRUE)
{
  /// update the cache for the given dataset(s)
  
  gEnv->SetValue("XSec.GSI.DelegProxy","2");
  if (setRootVersion) TProof::Mgr("ppillot@nansafmaster2.in2p3.fr")->SetROOTVersion("la-root");
  
  do {
    TProof::Open("ppillot@nansafmaster2.in2p3.fr/?N", "masteronly");
    if (!gProof) gSystem->Exec("sleep 30");
  } while (!gProof);
  
  cout << endl;
  
  if (datasets.EndsWith(".txt")) {
    
    ifstream inFile(datasets.Data());
    if (!inFile.is_open()) {
      Error("cacheDatasets", "unable to open file %s", datasets.Data());
      return;
    }
    
    while (!inFile.eof()) {
      
      TString query;
      query.ReadLine(inFile, kTRUE);
      if (query.IsNull()) continue;
      
      cacheOneDataset(query, outFileName);
      
    }
    
    inFile.close();
    
  } else cacheOneDataset(datasets, outFileName);
  
}

//-----------------------------------------------------------------------------
void cacheOneDataset(TString query, TString outFileName)
{
  /// update the cache for the given dataset
  
  std::ostream *outFile = (!outFileName.IsNull()) ? new ofstream(outFileName.Data(), ios_base::app) : 0x0;
  std::ostream &out = outFile ? *outFile : std::cout;
  
  TFileCollection *fileList = 0;
  if ((fileList = gProof->GetDataSet(Form("%s", query.Data())))) {
    
    out << TString::Format("%s\n--> %6lld files, %5.1f %% staged\n\n", query.Data(), fileList->GetNFiles(), fileList->GetStagedPercentage());
    
    delete fileList;
    
  } else out << TString::Format("%s\nnot found\n\n", query.Data());
  
  delete outFile;
  
}

