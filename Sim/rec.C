void rec() {
  
  AliReconstruction reco;
  reco.SetRunReconstruction("MUON");
//  reco.SetRunReconstruction("MUON ITS");
  reco.SetRunQA("MUON:ALL");
  reco.SetQAWriteExpert(AliQAv1::kMUON);
  
// switch off cleanESD
  reco.SetCleanESD(kFALSE);
  
  // Raw OCDB
  reco.SetDefaultStorage("alien://folder=/alice/data/2011/OCDB");
//  reco.SetDefaultStorage("alien://Folder=/alice/simulation/2008/v4-15-Release/Ideal");
  
  // GRP from local OCDB
  reco.SetSpecificStorage("GRP/GRP/Data",Form("local://%s",gSystem->pwd()));
  
  // MUON Tracker
//  reco.SetSpecificStorage("MUON/Calib/RejectList","alien://folder=/alice/cern.ch/user/l/laphecet/OCDB");
//  reco.SetSpecificStorage("MUON/Calib/RecoParam","alien://folder=/alice/cern.ch/user/p/ppillot/OCDB2012_newReco");
//  reco.SetSpecificStorage("MUON/Align/Data","alien://folder=/alice/simulation/2008/v4-15-Release/Residual");
//  reco.SetSpecificStorage("MUON/Align/Data","alien://folder=/alice/simulation/2008/v4-15-Release/Ideal");
//  reco.SetSpecificStorage("MUON/Align/Data","alien://folder=/alice/data/2011/OCDB");
  reco.SetSpecificStorage("MUON/Align/Data","alien://folder=/alice/cern.ch/user/p/ppillot/OCDB_PbPbSim");
  /*
  AliMUONRecoParam *muonRecoParam = AliMUONRecoParam::GetLowFluxParam();
  for (Int_t iCh=0; iCh<10; iCh++) {
    muonRecoParam->SetDefaultNonBendingReso(iCh,0.2);
    muonRecoParam->SetDefaultBendingReso(iCh,0.2);
  }
  reco.SetRecoParam("MUON",muonRecoParam);
  */
  // RecoParam for vertexerMC
  //reco.SetSpecificStorage("ITS/Calib/RecoParam","alien://folder=/alice/cern.ch/user/p/ppillot/OCDB_PbPbSim");
  
  TStopwatch timer;
  timer.Start();
  reco.Run();
  timer.Stop();
  timer.Print();
}
