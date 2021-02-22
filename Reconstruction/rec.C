void rec(const char *filename="raw.root")
{
  /////////////////////////////////////////////////////////////////////////////////////////
  //
  // Reconstruction script for 2012 RAW data
  //
  /////////////////////////////////////////////////////////////////////////////////////////
  
  AliReconstruction MuonRec;

  // Set the CDB storage location
  AliCDBManager * man = AliCDBManager::Instance();
  if (gSystem->AccessPathName("OCDB.root", kFileExists)==0) {        
    MuonRec.SetCDBSnapshotMode("OCDB.root");
  } else {
    man->SetDefaultStorage("raw://");
    //man->SetDefaultStorage("local:///Users/pillot/Work/Alice/Data/2015/CDBMirror/alice/data/2015/OCDB");
    //man->SetSpecificStorage("MUON/Align/Data","alien://folder=/alice/cern.ch/user/h/hupereir/CDB/LHC15_realign_all_4");
    //man->SetSpecificStorage("MUON/Align/Data","alien://folder=/alice/data/2015/OCDB");
    //man->SetSpecificStorage("MUON/Calib/RecoParam","alien://folder=/alice/cern.ch/user/p/ppillot/OCDB2015_PbPb");
  }

  // Reconstruction settings
  //MuonRec.SetRunReconstruction("MUON ITS VZERO");
  //MuonRec.SetRunReconstruction("MUON");
  MuonRec.SetRunReconstruction("ALL -TPC -TRD -TOF -HLT -PMD -FMD -EMCAL -PHOS -HMPID -ACORDE");
  
  MuonRec.SetOption("MUON","SAVEDIGITS");
  
  MuonRec.SetRunQA("Global MUON:ALL") ;
  //MuonRec.SetRunQA("MUON:ALL");
  //MuonRec.SetRunQA(":");
  MuonRec.SetQARefDefaultStorage("local://$ALICE_ROOT/QAref") ;
  
  // AliReconstruction settings
  MuonRec.SetWriteESDfriend(kFALSE);
  MuonRec.SetWriteAlignmentData();
  MuonRec.SetUseTrackingErrorsForAlignment("ITS");
  
  MuonRec.SetInput(filename);
  //MuonRec.SetInput(Form("%s?Trigger=C0MSL-ABCE-NOPF-ALLNOTRD",filename));
  //AliRawReaderChain::SetSearchPath("/alice/data/2018/LHC18q");
  //MuonRec.SetInput("raw://run295584?Trigger=CMUL7-B-NOPF-MUFAST");
  
  // switch off cleanESD
  MuonRec.SetCleanESD(kFALSE);
  
  //Ignore SetStopOnError
  MuonRec.SetStopOnError(kFALSE);
  
  AliLog::Flush();
  MuonRec.Run();
}

