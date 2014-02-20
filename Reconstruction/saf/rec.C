void rec(TString input, TString output, TString cdbSnapshot = "CDBMirror", Int_t run = -1)
{

  // Set the CDB storage location
  AliCDBManager * man = AliCDBManager::Instance();
  //man->SetDefaultStorage("raw://");
  man->SetDefaultStorage("local:///Users/pillot/Work/Alice/Work/Data/2012/LHC12h/raw/190147/saf/CDBMirror/alice/data/2012/OCDB");
  //man->SetDefaultStorage("local:///Users/pillot/Work/Alice/Work/Data/2012/LHC12h/raw/190206/saf/CDBMirror/alice/data/2012/OCDB");
  //man->SetDefaultStorage("alien://folder=/alice/data/2012/OCDB");
  
  // only SPD reconstruction
  man->SetSpecificStorage("ITS/Calib/RecoParam","local:///Users/pillot/Work/Alice/Work/Data/2012/LHC12h/raw/190147/OCDB_RecoParams");
  //man->SetSpecificStorage("ITS/Calib/RecoParam","local:///Users/pillot/Work/Alice/Work/Data/2012/LHC12h/raw/190206/OCDB_RecoParams");
  
  // new MUON parameters
  //man->SetSpecificStorage("MUON/Calib/RecoParam","local:///Users/pillot/Work/Alice/Work/Data/2012/LHC12h/raw/190147/OCDB_RecoParams");
  //man->SetSpecificStorage("MUON/Calib/RecoParam","local:///Users/pillot/Work/Alice/Work/Data/2012/LHC12h/raw/190147/OCDB_ManuOcc1");
  //man->SetSpecificStorage("MUON/Calib/RecoParam","local:///Users/pillot/Work/Alice/Work/Data/2012/LHC12h/raw/190147/OCDB_ManuOcc1.5");
  //man->SetSpecificStorage("MUON/Calib/RecoParam","local:///Users/pillot/Work/Alice/Work/Data/2012/LHC12h/raw/190147/OCDB_MoreTrkCand_ManuOcc1");
  //man->SetSpecificStorage("MUON/Calib/RecoParam","local:///Users/pillot/Work/Alice/Work/Data/2012/LHC12h/raw/190147/OCDB_MoreTrkCand_ManuOcc1.5");
  //man->SetSpecificStorage("MUON/Calib/RecoParam","local:///Users/pillot/Work/Alice/Work/Data/2012/LHC12h/raw/190147/OCDB_NoLimit");
  //man->SetSpecificStorage("MUON/Calib/RecoParam","local:///Users/pillot/Work/Alice/Work/Data/2012/LHC12h/raw/190147/OCDB_NoLimit_ManuOcc1.5");
  man->SetSpecificStorage("MUON/Calib/RecoParam","local:///Users/pillot/Work/Alice/Work/Data/2012/LHC12h/raw/190147/OCDB_MoreTrkCand_NoLimit");
  //man->SetSpecificStorage("MUON/Calib/RecoParam","local:///Users/pillot/Work/Alice/Work/Data/2012/LHC12h/raw/190147/OCDB_MoreTrkCand_NoLimit_ManuOcc1");
  //man->SetSpecificStorage("MUON/Calib/RecoParam","local:///Users/pillot/Work/Alice/Work/Data/2012/LHC12h/raw/190147/OCDB_MoreTrkCand_NoLimit_ManuOcc1.5");
  //man->SetSpecificStorage("MUON/Calib/RecoParam","local:///Users/pillot/Work/Alice/Work/Data/2012/LHC12h/raw/190147/OCDB_MoreTrkCand_Limit50000");
  //man->SetSpecificStorage("MUON/Calib/RecoParam","local:///Users/pillot/Work/Alice/Work/Data/2012/LHC12h/raw/190147/OCDB_MoreTrkCand_Limit50000_ManuOcc1.5");
  //man->SetSpecificStorage("MUON/Calib/RecoParam","local:///Users/pillot/Work/Alice/Work/Data/2012/LHC12h/raw/190206/OCDB_noHVcut");
  
  // here we check if local CDB snapshot is requested
  if (!cdbSnapshot.IsNull() && run >= 0) {
    TString locCDB = (cdbSnapshot.BeginsWith("/")) ? cdbSnapshot : Form("%s/%s", gSystem->pwd(), cdbSnapshot.Data());
    gROOT->LoadMacro("UpdateCDBSnapshot.C+");
    UpdateCDBSnapshot(run, locCDB.Data());
  }
  
  // Reconstruction settings
  AliReconstruction rec;

  // QA options
  rec.SetRunQA("MUON:ALL");
  rec.SetQAWriteExpert(AliQAv1::kMUON);
//  rec.SetRunQA(":");
  rec.SetRunGlobalQA(kFALSE);
  rec.SetQARefDefaultStorage("local://$ALICE_ROOT/QAref") ;

  // AliReconstruction settings
  rec.SetWriteESDfriend(kFALSE);
  if (input.EndsWith(".root") && !input.Contains("collection")) rec.SetInput(input.Data());
  else if (input.BeginsWith("raw://run")) rec.SetInput(input.Data());
  else rec.SetInput(Form("collection://%s", input.Data()));
  //rec.SetEventRange(0,1257);
  rec.SetRunReconstruction("MUON ITS");
  //rec.SetRunMultFinder(kFALSE);

  // switch off cleanESD
  rec.SetCleanESD(kFALSE);

  //Ignore SetStopOnError
  rec.SetStopOnError(kFALSE);
  
  //  rec.SetOutput(Form("root://localhost/aaf_data/alishift/run%d/root_archive.zip#AliESDs.root:AliESDs.root,AliESDfriends.root,ITS.RecPoints.root,TPC.RecPoints.root,TRD.RecPoints.root,TOF.RecPoints.root,galice.root@dataset://run%d",runNumber,runNumber));
  //rec.SetOutput(Form("/data/acrtest/run%d/root_archive.zip#AliESDs.root:AliESDs.root,AliESDfriends.root,ITS.RecPoints.root,TPC.RecPoints.root,TRD.RecPoints.root,TOF.RecPoints.root,galice.root@dataset://run%d",runNumber,runNumber));
  //rec.SetOutput(Form("root_archive.zip#AliESDs.root:AliESDs.root,AliESDfriends.root@dataset://run%d",runNumber));
  // rec.SetOutput(Form("root_archive.zip#AliESDs.root:AliESDs.root,AliESDfriends.root,*.RecPoints.root,galice.root@dataset://run%d",runNumber));
  if (!output.IsNull()) rec.SetOutput(Form("AliESDs.root@dataset://%s", output.Data()));
  
  rec.Run();
  
}
