AliAnalysisTaskJPsi* AddTaskJPsi(TString extension = "")
{
  // Get the pointer to the existing analysis manager via the static access method.
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr) {
    Error("AddTaskJPsi","AliAnalysisManager not set!");
    return NULL;
  }
  
  // Create and configure task
  TString suffix = (!extension.IsNull()) ? Form("_%s",extension.Data()) : "";
  AliAnalysisTaskJPsi *task = new AliAnalysisTaskJPsi(Form("JPsi%s",suffix.Data()));
  if (!task) {
    Error("AddTaskJPsi", "JPsi task cannot be created!");
    return NULL;
  }
  
  // Add task to analysis manager
  mgr->AddTask(task);
  
  // Connect input container
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  
  // Create and connect output containers
  AliAnalysisDataContainer *output = mgr->CreateContainer(Form("cOut%s",suffix.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, "Output.root");
  mgr->ConnectOutput(task, 1, output);
  
  return task;
}

