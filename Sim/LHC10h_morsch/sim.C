void sim(Int_t nev=10) {
  AliSimulation simulator;
  simulator.SetWriteRawData("ALL","raw.root",kTRUE);
  simulator.SetMakeSDigits("TRD TOF PHOS HMPID EMCAL MUON FMD ZDC PMD T0 VZERO");
  simulator.SetMakeDigits("ITS TPC TRD TOF PHOS HMPID EMCAL MUON FMD ZDC PMD T0 VZERO");
  simulator.SetMakeDigitsFromHits("ITS TPC");
//
//
// Ideal OCDB
  simulator.SetDefaultStorage("alien://Folder=/alice/simulation/2008/v4-15-Release/Ideal/");
//
// Read GRP Data from RAW
  simulator.SetSpecificStorage("GRP/GRP/Data",                  "alien://Folder=/alice/data/2010/OCDB");
//
//
// Mean vertex from RAW OCDB 
  simulator.SetSpecificStorage("GRP/Calib/MeanVertexSPD",       "alien://folder=/alice/data/2010/OCDB"); 
  simulator.SetSpecificStorage("GRP/Calib/MeanVertex",          "alien://folder=/alice/data/2010/OCDB");
//
// Clock phase from RAW OCDB 
  simulator.SetSpecificStorage("GRP/Calib/LHCClockPhase",       "alien://folder=/alice/data/2010/OCDB");
//
// ITS
//    SDD from RAW OCDB
  simulator.SetSpecificStorage("ITS/Calib/CalibSDD",            "alien://Folder=/alice/data/2010/OCDB");
//    SSD
  simulator.SetSpecificStorage("ITS/Calib/NoiseSSD",            "alien://Folder=/alice/data/2010/OCDB");
  simulator.SetSpecificStorage("ITS/Calib/BadChannelsSSD",      "alien://Folder=/alice/data/2010/OCDB"); 
//
// TRD from RAW OCDB
  simulator.SetSpecificStorage("TRD/Calib/ChamberStatus",       "alien://folder=/alice/data/2010/OCDB");
  simulator.SetSpecificStorage("TRD/Calib/PadStatus",           "alien://folder=/alice/data/2010/OCDB");
//
// TOF from RAW OCDB
  simulator.SetSpecificStorage("TOF/Calib/Status",              "alien://folder=/alice/data/2010/OCDB");
//
// VZERO
  simulator.SetSpecificStorage("VZERO/Calib/Data",              "alien://Folder=/alice/data/2010/OCDB");
  simulator.SetSpecificStorage("VZERO/Trigger/Data",            "alien://folder=/alice/data/2010/OCDB");
  simulator.SetSpecificStorage("VZERO/Calib/RecoParam",         "alien://folder=/alice/data/2010/OCDB");
  simulator.SetSpecificStorage("VZERO/Calib/TimeSlewing",       "alien://folder=/alice/data/2010/OCDB");
  simulator.SetSpecificStorage("VZERO/Calib/TimeDelays",        "alien://folder=/alice/data/2010/OCDB");
//
// FMD from RAW OCDB
  simulator.SetSpecificStorage("FMD/Calib/Pedestal",            "alien://folder=/alice/data/2010/OCDB");
  simulator.SetSpecificStorage("FMD/Calib/PulseGain",           "alien://folder=/alice/data/2010/OCDB");
  simulator.SetSpecificStorage("FMD/Calib/Dead",                "alien://folder=/alice/data/2010/OCDB");
  simulator.SetSpecificStorage("FMD/Calib/AltroMap",            "alien://folder=/alice/data/2010/OCDB");
//
// EMCAL from RAW OCDB
  simulator.SetSpecificStorage("EMCAL/Calib/Data",             "alien://Folder=/alice/data/2010/OCDB");
//
// MUON 
// Pedestals  
  simulator.SetSpecificStorage("MUON/Calib/Pedestals",          "alien://folder=/alice/data/2010/OCDB");
// MappingData
  simulator.SetSpecificStorage("MUON/Calib/MappingData",        "alien://folder=/alice/data/2010/OCDB");
// Trigger LuT 
  simulator.SetSpecificStorage("MUON/Calib/TriggerLut",         "alien://folder=/alice/data/2010/OCDB");
// Trigger efficiency 
  simulator.SetSpecificStorage("MUON/Calib/TriggerEfficiency",  "alien://folder=/alice/simulation/2008/v4-15-Release/Full");
//
// ZDC
  simulator.SetSpecificStorage("ZDC/Calib/Pedestals",           "alien://folder=/alice/data/2010/OCDB");
//
//
// Vertex and Mag.field from OCDB

  simulator.UseVertexFromCDB();
  simulator.UseMagFieldFromGRP();

//
// The rest
//

  TStopwatch timer;
  timer.Start();
  simulator.Run(nev);
  timer.Stop();
  timer.Print();
}
