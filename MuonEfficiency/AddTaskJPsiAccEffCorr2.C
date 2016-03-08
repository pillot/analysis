AliAnalysisTaskJPsiAccEffCorr2* AddTaskJPsiAccEffCorr2(TString extension = "")
{
  /// Add AliAnalysisTaskJPsiAccEffCorr2 to the train (Philippe Pillot)
  
  // Get the pointer to the existing analysis manager via the static access method.
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr) {
    Error("AddTaskJPsiAccEffCorr2","AliAnalysisManager not set!");
    return NULL;
  }
  
  // This task run on AODs
  TString type = mgr->GetInputEventHandler()->GetDataType();
  if (!type.Contains("ESD") && !type.Contains("AOD")) {
    Error("AddTaskJPsiAccEffCorr2", "ESD or AOD input handler needed!");
    return NULL;
  }
  
  // Create and configure task
  TString suffix = (!extension.IsNull()) ? Form("_%s",extension.Data()) : "";
  AliAnalysisTaskJPsiAccEffCorr2 *task = new AliAnalysisTaskJPsiAccEffCorr2("JPsiAccEffCorr");
  if (!task) {
    Error("AddTaskJPsiAccEffCorr2", "JPsi acc*eff task cannot be created!");
    return NULL;
  }
  
  // Add task to analysis manager
  mgr->AddTask(task);
  
  // Connect input container
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  
  // Define output file directory
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  if ( outputfile.IsNull() ) {
    Error("AddTaskJPsiAccEffCorr2", "Common output file is not defined!");
    return NULL;
  }
  outputfile += ":MUON_JPsiAccEff";
  
  // Create and connect output containers
  AliAnalysisDataContainer *histo = mgr->CreateContainer(Form("Histograms%s",suffix.Data()), TObjArray::Class(), AliAnalysisManager::kOutputContainer, outputfile);
  AliAnalysisDataContainer *eventStat = mgr->CreateContainer(Form("EventCounters%s",suffix.Data()), AliCounterCollection::Class(), AliAnalysisManager::kOutputContainer, outputfile);
  AliAnalysisDataContainer *jPsiStat = mgr->CreateContainer(Form("JPsiCounters%s",suffix.Data()), AliCounterCollection::Class(), AliAnalysisManager::kOutputContainer, outputfile);
  AliAnalysisDataContainer *massVspT = mgr->CreateContainer(Form("MassVspT%s",suffix.Data()), TObjArray::Class(), AliAnalysisManager::kOutputContainer, outputfile);
  AliAnalysisDataContainer *massVsy = mgr->CreateContainer(Form("MassVsy%s",suffix.Data()), TObjArray::Class(), AliAnalysisManager::kOutputContainer, outputfile);
  mgr->ConnectOutput(task, 1, histo);
  mgr->ConnectOutput(task, 2, eventStat);
  mgr->ConnectOutput(task, 3, jPsiStat);
  mgr->ConnectOutput(task, 4, massVspT);
  mgr->ConnectOutput(task, 5, massVsy);
  
  return task;
}

