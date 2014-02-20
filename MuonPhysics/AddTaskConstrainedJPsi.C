AliAnalysisTaskConstrainedJPsi* AddTaskConstrainedJPsi()
{
  /// Add AliAnalysisTaskConstrainedJPsi to the train (Philippe Pillot)
  
  // Get the pointer to the existing analysis manager via the static access method.
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr) { 
    Error("AddTaskConstrainedJPsi","AliAnalysisManager not set!");
    return NULL;
  }
  
  // This task run on ESDs
  TString type = mgr->GetInputEventHandler()->GetDataType();
  if (!type.Contains("ESD")) {
    Error("AddTaskConstrainedJPsi", "ESD input handler needed!");
    return NULL;
  }
  
  // Create and configure task
  AliAnalysisTaskConstrainedJPsi *task = new AliAnalysisTaskConstrainedJPsi("ConstrainedJPsi");
  if (!task) {
    Error("AddTaskConstrainedJPsi", "Constrained JPsi task cannot be created!");
    return NULL;
  }
  
  // Add task to analysis manager
  mgr->AddTask(task);
  
  // Connect input container
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  
  // Define output file directory
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  if ( outputfile.IsNull() ) {
    Error("AddTaskConstrainedJPsi", "Common output file is not defined!");
    return NULL;
  }
  outputfile += ":MUON_ConstrainedJPsi";
  
  // Create and connect output containers
  AliAnalysisDataContainer *data = mgr->CreateContainer("Data", TObjArray::Class(), AliAnalysisManager::kOutputContainer, outputfile);
  AliAnalysisDataContainer *constData = mgr->CreateContainer("ConstData", TObjArray::Class(), AliAnalysisManager::kOutputContainer, outputfile);
  AliAnalysisDataContainer *resK = mgr->CreateContainer("ResK", TObjArray::Class(), AliAnalysisManager::kOutputContainer, outputfile);
  AliAnalysisDataContainer *constResK = mgr->CreateContainer("ConstResK", TObjArray::Class(), AliAnalysisManager::kOutputContainer, outputfile);
  AliAnalysisDataContainer *res = mgr->CreateContainer("Res", TObjArray::Class(), AliAnalysisManager::kOutputContainer, outputfile);
  AliAnalysisDataContainer *constRes = mgr->CreateContainer("ConstRes", TObjArray::Class(), AliAnalysisManager::kOutputContainer, outputfile);
  mgr->ConnectOutput(task, 1, data);
  mgr->ConnectOutput(task, 2, constData);
  mgr->ConnectOutput(task, 3, resK);
  mgr->ConnectOutput(task, 4, constResK);
  mgr->ConnectOutput(task, 5, res);
  mgr->ConnectOutput(task, 6, constRes);
  
  return task;
}

