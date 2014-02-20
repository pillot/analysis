void rec(const char *filename="raw.root")
{
  /////////////////////////////////////////////////////////////////////////////////////////
  //
  // Reconstruction script for 2012 RAW data
  //
  /////////////////////////////////////////////////////////////////////////////////////////
  
  
  // Set the CDB storage location
  AliCDBManager * man = AliCDBManager::Instance();
  //  man->SetDefaultStorage("raw://");
  man->SetDefaultStorage("alien://folder=/alice/data/2013/OCDB");
  //man->SetSpecificStorage("MUON/Align/Data","alien://folder=/alice/cern.ch/user/j/jcastill/MATFtestCDBSurvey");
  
  // Reconstruction settings
  AliReconstruction MuonRec;
  
  MuonRec.SetRunReconstruction("MUON ITS VZERO");
  //MuonRec.SetRunReconstruction("MUON");
  
  MuonRec.SetOption("MUON","SAVEDIGITS");
  
  //MuonRec.SetRunQA("MUON:ALL");
  MuonRec.SetRunQA(":");
  
  MuonRec.SetQARefDefaultStorage("local://$ALICE_ROOT/QAref") ;
  
  // AliReconstruction settings
  MuonRec.SetWriteESDfriend(kFALSE);
  MuonRec.SetInput(filename);
  
  // switch off cleanESD
  MuonRec.SetCleanESD(kFALSE);
  
  //Ignore SetStopOnError
  MuonRec.SetStopOnError(kFALSE);
  
  AliLog::Flush();
  MuonRec.Run();
}

