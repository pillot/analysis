void sim(Int_t nev=1000) {
  
  AliSimulation simulator;
  simulator.SetTriggerConfig("MUON");
  simulator.SetMakeDigits("MUON");
//  simulator.SetMakeDigits("MUON ITS");
  simulator.SetMakeSDigits("MUON");
  simulator.SetMakeDigitsFromHits("");
//  simulator.SetMakeDigitsFromHits("ITS");
  simulator.SetRunQA("MUON:ALL");
  simulator.SetRunHLT("");
  
  // raw OCDB
  simulator.SetDefaultStorage("alien://folder=/alice/data/2011/OCDB");
  //simulator.SetDefaultStorage("alien://Folder=/alice/simulation/2008/v4-15-Release/Ideal");
  //simulator.SetSpecificStorage("GRP/GRP/Data", "alien://Folder=/alice/data/2011/OCDB");
  
  // MUON Tracker
//  simulator.SetSpecificStorage("MUON/Calib/Gains","alien://folder=/alice/simulation/2008/v4-15-Release/Ideal");
//  simulator.SetSpecificStorage("MUON/Calib/RejectList","alien://folder=/alice/cern.ch/user/l/laphecet/OCDB");
  simulator.SetSpecificStorage("MUON/Align/Data","alien://folder=/alice/simulation/2008/v4-15-Release/Ideal");
//  simulator.SetSpecificStorage("MUON/Align/Data","alien://folder=/alice/simulation/2008/v4-15-Release/Full");
//  simulator.SetSpecificStorage("MUON/Align/Data","alien://folder=/alice/cern.ch/user/j/jcastill/LHC10hMisAlignCDB");
//  simulator.SetSpecificStorage("MUON/Align/Data","alien://folder=/alice/cern.ch/user/j/jcastill/LHC11hMisAlignCDB3");
//  simulator.SetSpecificStorage("MUON/Align/Data","alien://folder=/alice/cern.ch/user/j/jcastill/ReAligni00pbpb11CDB2");
//  simulator.SetSpecificStorage("MUON/Align/Data","alien://folder=/alice/cern.ch/user/j/jcastill/pbpb11wrk/LHC11hMisAlignCDB4");
  
  // Vertex and Mag.field from OCDB
//  simulator.UseVertexFromCDB();
  simulator.UseMagFieldFromGRP();
  
  // The rest
  TStopwatch timer;
  timer.Start();
  simulator.Run(nev);
  timer.Stop();
  timer.Print();
}
