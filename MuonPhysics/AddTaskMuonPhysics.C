AliAnalysisTaskMuonPhysics* AddTaskMuonPhysics(TString extension = "")
{
  /// Add AliAnalysisTaskMuonPhysics to the train (Philippe Pillot)
  
  // Get the pointer to the existing analysis manager via the static access method.
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr) { 
    Error("AddTaskMuonPhysics","AliAnalysisManager not set!");
    return NULL;
  }
  
  // This task run on ESDs or AODs
  TString type = mgr->GetInputEventHandler()->GetDataType();
  if (!type.Contains("ESD") && !type.Contains("AOD")) {
    Error("AddTaskMuonPhysics", "ESD input handler needed!");
    return NULL;
  }
  
  // Create and configure task
  TString suffix = (!extension.IsNull()) ? Form("_%s",extension.Data()) : "";
  AliAnalysisTaskMuonPhysics *task = new AliAnalysisTaskMuonPhysics(Form("MuonPhysics%s",suffix.Data()));
  if (!task) {
    Error("AddTaskMuonPhysics", "Muon physics task cannot be created!");
    return NULL;
  }
  
  // Add task to analysis manager
  mgr->AddTask(task);
  
  // Connect input container
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  
  // Define output file directory
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  if ( outputfile.IsNull() ) {
    Error("AddTaskMuonPhysics", "Common output file is not defined!");
    return NULL;
  }
  outputfile += ":MUON_Physics";
  
  // Create and connect output containers
  AliAnalysisDataContainer *histo = mgr->CreateContainer(Form("Histograms%s",suffix.Data()), TObjArray::Class(), AliAnalysisManager::kOutputContainer, outputfile);
  AliAnalysisDataContainer *trackStat = mgr->CreateContainer(Form("trackCounters%s",suffix.Data()), AliCounterCollection::Class(), AliAnalysisManager::kOutputContainer, outputfile);
  AliAnalysisDataContainer *eventStat = mgr->CreateContainer(Form("eventCounters%s",suffix.Data()), AliCounterCollection::Class(), AliAnalysisManager::kOutputContainer, outputfile);
  AliAnalysisDataContainer *trigStat = mgr->CreateContainer(Form("trigCounters%s",suffix.Data()), AliCounterCollection::Class(), AliAnalysisManager::kOutputContainer, outputfile);
  mgr->ConnectOutput(task, 1, histo);
  mgr->ConnectOutput(task, 2, trackStat);
  mgr->ConnectOutput(task, 3, eventStat);
  mgr->ConnectOutput(task, 4, trigStat);
  
  return task;
}

