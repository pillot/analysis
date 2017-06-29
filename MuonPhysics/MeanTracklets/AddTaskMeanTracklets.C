AliAnalysisTaskMeanTracklets* AddTaskMeanTracklets()
{
  /// Add AliAnalysisTaskMeanTracklets to the train (Philippe Pillot)
  
  // Get the pointer to the existing analysis manager via the static access method.
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr) {
    Error("AddTaskMeanTracklets","AliAnalysisManager not set!");
    return NULL;
  }
  
  // This task run on AODs
  TString type = mgr->GetInputEventHandler()->GetDataType();
  if (!type.Contains("AOD")) {
    Error("AddTaskMeanTracklets", "AOD input handler needed!");
    return NULL;
  }
  
  // Create and configure task
  AliAnalysisTaskMeanTracklets *task = new AliAnalysisTaskMeanTracklets("MeanTracklets");
  if (!task) {
    Error("AddTaskMeanTracklets", "MeanTracklets task cannot be created!");
    return NULL;
  }
  
  // Add task to analysis manager
  mgr->AddTask(task);
  
  // Connect input container
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  
  // Define output file directory
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  if ( outputfile.IsNull() ) {
    Error("AddTaskMeanTracklets", "Common output file is not defined!");
    return NULL;
  }
  
  // Create and connect output containers
  AliAnalysisDataContainer *events = mgr->CreateContainer("events", AliCounterCollection::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data());
  AliAnalysisDataContainer *hNtrk = mgr->CreateContainer("hNtrk", THnSparse::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data());
  AliAnalysisDataContainer *hNtrkCorr = mgr->CreateContainer("hNtrkCorr", THnSparse::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data());
  AliAnalysisDataContainer *pMeanNtrkVsZvtx = mgr->CreateContainer("pMeanNtrkVsZvtx", TProfile::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data());
  AliAnalysisDataContainer *pMeanNtrkVsZvtxCorr = mgr->CreateContainer("pMeanNtrkVsZvtxCorr", TProfile::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data());
  mgr->ConnectOutput(task, 1, events);
  mgr->ConnectOutput(task, 2, hNtrk);
  mgr->ConnectOutput(task, 3, hNtrkCorr);
  mgr->ConnectOutput(task, 4, pMeanNtrkVsZvtx);
  mgr->ConnectOutput(task, 5, pMeanNtrkVsZvtxCorr);
  
  return task;
}

