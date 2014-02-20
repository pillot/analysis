void rec() {  
  AliReconstruction reco;
  reco.SetRunReconstruction("ITS TPC TRD TOF PHOS HMPID EMCAL MUON FMD ZDC PMD T0 VZERO");

  reco.SetWriteESDfriend();
  reco.SetWriteAlignmentData();

  reco.SetDefaultStorage("alien://Folder=/alice/simulation/2008/v4-15-Release/Residual/");
//
// Vertex from RAW OCDB
  reco.SetSpecificStorage("GRP/Calib/MeanVertexTPC",   "alien://folder=/alice/data/2010/OCDB"); 
  reco.SetSpecificStorage("GRP/Calib/MeanVertex",      "alien://folder=/alice/data/2010/OCDB");
  reco.SetSpecificStorage("GRP/Calib/MeanVertexSPD",   "alien://folder=/alice/data/2010/OCDB");
  reco.SetSpecificStorage("GRP/Calib/RecoParam",       "alien://folder=/alice/data/2010/OCDB"); 
//
// Clock phase from RAW OCDB 
  reco.SetSpecificStorage("GRP/Calib/LHCClockPhase",   "alien://folder=/alice/data/2010/OCDB");
//
// ITS
  reco.SetSpecificStorage("ITS/Calib/RecoParam",       "alien://folder=/alice/data/2010/OCDB"); 
  reco.SetSpecificStorage("ITS/Calib/SPDDead",         "alien://folder=/alice/data/2010/OCDB");
  reco.SetSpecificStorage("TRIGGER/SPD/PITConditions", "alien://folder=/alice/data/2010/OCDB");
// SDD
  reco.SetSpecificStorage("ITS/Calib/CalibSDD",        "alien://Folder=/alice/data/2010/OCDB"); 
// SSD
  reco.SetSpecificStorage("ITS/Calib/NoiseSSD",        "alien://Folder=/alice/data/2010/OCDB");
  reco.SetSpecificStorage("ITS/Calib/BadChannelsSSD",  "alien://Folder=/alice/data/2010/OCDB"); 

// TPC
// TPC from RAW OCDB
  reco.SetSpecificStorage("TPC/Calib/PadGainFactor",   "alien://folder=/alice/data/2010/OCDB");
//
// TRD from RAW OCDB
  reco.SetSpecificStorage("TRD/Calib/ChamberStatus",   "alien://folder=/alice/data/2010/OCDB");
  reco.SetSpecificStorage("TRD/Calib/PadStatus",       "alien://folder=/alice/data/2010/OCDB");
//
// TOF from RAW OCDB
  reco.SetSpecificStorage("TOF/Calib/Status",          "alien://folder=/alice/data/2010/OCDB");
//
// EMCAL from RAW OCDB
  reco.SetSpecificStorage("EMCAL/Calib/Data",          "alien://Folder=/alice/data/2010/OCDB");
  reco.SetSpecificStorage("EMCAL/Calib/Pedestals",     "alien://Folder=/alice/data/2010/OCDB");
  reco.SetSpecificStorage("EMCAL/Calib/RecoParam",     "alien://Folder=/alice/data/2010/OCDB");

//
// PHOS from RAW OCDB
  reco.SetSpecificStorage("PHOS/Calib/EmcBadChannels", "alien://Folder=/alice/data/2010/OCDB");
  reco.SetSpecificStorage("PHOS/Calib/RecoParam",      "alien://Folder=/alice/data/2010/OCDB");
  reco.SetSpecificStorage("PHOS/Calib/Mapping",        "alien://Folder=/alice/data/2010/OCDB");
//
// VZERO
  reco.SetSpecificStorage("VZERO/Calib/Data",          "alien://folder=/alice/data/2010/OCDB");
  reco.SetSpecificStorage("VZERO/Trigger/Data",        "alien://folder=/alice/data/2010/OCDB");
  reco.SetSpecificStorage("VZERO/Calib/RecoParam",     "alien://folder=/alice/data/2010/OCDB");
  reco.SetSpecificStorage("VZERO/Calib/TimeSlewing",   "alien://folder=/alice/data/2010/OCDB");
  reco.SetSpecificStorage("VZERO/Calib/TimeDelays",    "alien://folder=/alice/data/2010/OCDB");
//
// FMD from RAW OCDB
  reco.SetSpecificStorage("FMD/Calib/Pedestal",        "alien://folder=/alice/data/2010/OCDB");
  reco.SetSpecificStorage("FMD/Calib/PulseGain",       "alien://folder=/alice/data/2010/OCDB");
  reco.SetSpecificStorage("FMD/Calib/Dead",            "alien://folder=/alice/data/2010/OCDB");
  reco.SetSpecificStorage("FMD/Calib/AltroMap",        "alien://folder=/alice/data/2010/OCDB");
//
// MUON
// Config
  reco.SetSpecificStorage("MUON/Calib/Config",         "alien://folder=/alice/data/2010/OCDB");
// MappingData
  reco.SetSpecificStorage("MUON/Calib/MappingData",    "alien://folder=/alice/data/2010/OCDB");
// MappingRunData
  reco.SetSpecificStorage("MUON/Calib/MappingRunData", "alien://folder=/alice/data/2010/OCDB");
// Neighbours
  reco.SetSpecificStorage("MUON/Calib/Neighbours",     "alien://folder=/alice/data/2010/OCDB");
// OccupancyMap
  reco.SetSpecificStorage("MUON/Calib/OccupancyMap",   "alien://folder=/alice/data/2010/OCDB");
// Pedestals
  reco.SetSpecificStorage("MUON/Calib/Pedestals",      "alien://folder=/alice/data/2010/OCDB");
// RecParam
  reco.SetSpecificStorage("MUON/Calib/RecoParam",      "alien://folder=/alice/data/2010/OCDB");
// RejectList
  reco.SetSpecificStorage("MUON/Calib/RejectList",     "alien://folder=/alice/data/2010/OCDB");
// Alignment (2nd alignment in Residual, not needed Residual is the default)
// reco.SetSpecificStorage("MUON/Align/Data",          "alien://folder=/alice/simulation/2008/v4-15-Release/Residual");
// Trigger LuT 
  reco.SetSpecificStorage("MUON/Calib/TriggerLut",     "alien://folder=/alice/data/2010/OCDB");
// Trigger DCS (for QA)
  reco.SetSpecificStorage("MUON/Calib/TriggerDCS",     "alien://folder=/alice/data/2010/OCDB");
//
// ZDC
  reco.SetSpecificStorage("ZDC/Calib/Pedestals",       "alien://folder=/alice/data/2010/OCDB");
  // No write access to the OCDB => specific storage

  reco.SetSpecificStorage("GRP/GRP/Data",
                          Form("local://%s",gSystem->pwd()));

  TStopwatch timer;
  timer.Start();
  reco.Run();
  timer.Stop();
  timer.Print();
}
