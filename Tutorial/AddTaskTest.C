AliAnalysisTaskTest* AddTaskTest()
{
  /// Add AliAnalysisTaskTest to the train (Philippe Pillot)
  
  // Get the pointer to the existing analysis manager via the static access method.
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr) {
    Error("AddTaskTest","AliAnalysisManager not set!");
    return NULL;
  }
  
  // This task run on AODs
  TString type = mgr->GetInputEventHandler()->GetDataType();
  if (!type.Contains("AOD")) {
    Error("AddTaskTest", "AOD input handler needed!");
    return NULL;
  }
  
  // Create and configure task
  AliAnalysisTaskTest *task = new AliAnalysisTaskTest("Test");
  if (!task) {
    Error("AddTaskTest", "Test task cannot be created!");
    return NULL;
  }
  
  // Add task to analysis manager
  mgr->AddTask(task);
  
  // Connect input container
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  
  // Define output file directory
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  if ( outputfile.IsNull() ) {
    Error("AddTaskTest", "Common output file is not defined!");
    return NULL;
  }
  
  // Create and connect output containers
  AliAnalysisDataContainer *histo = mgr->CreateContainer("hPt", TH1F::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data());
  mgr->ConnectOutput(task, 1, histo);
  
  return task;
}

