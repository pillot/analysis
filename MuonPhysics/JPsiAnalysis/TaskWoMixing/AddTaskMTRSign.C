#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TString.h"
#include "TObjArray.h"

#include "AliLog.h"
#include "AliVEventHandler.h"

#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"

#include "AliAnalysisTaskMTRSign.h"
#endif

AliAnalysisTaskMTRSign* AddTaskMTRSign ( Bool_t isMC, TString changeName = "" )
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddtaskMTRSign", "No analysis manager to connect to.");
    return NULL;
  }

//  TString type = mgr->GetInputEventHandler()->GetDataType();
//  if (!type.Contains("ESD") && !type.Contains("AOD") ) {
//    ::Error("AddtaskMTRSign", "MTRSign task needs the manager to have ESD or AOD input handler.");
//    return NULL;
//  }

  // Create container
  TString outputfile = "Output.root"; // mgr->GetCommonFileName();
  if ( ! outputfile.IsNull() ) outputfile += ":PWG_MTRSign" + changeName;
  else outputfile = "MTRSignAnalysis" + changeName + ".root";

  TString containerName = "MTRSignOut" + changeName;
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(containerName.Data(),AliMergeableCollection::Class(),AliAnalysisManager::kOutputContainer,outputfile);

  // Create task
  TString taskName = "MTRSignTask" + changeName;
  AliAnalysisTaskMTRSign *task = new AliAnalysisTaskMTRSign(taskName.Data());
  if ( isMC ) task->GetMuonTrackCuts()->SetIsMC();
  task->GetMuonEventCuts()->SkipTestsNonInFilterMask();
  mgr->AddTask(task);

   // Connect containers
  mgr->ConnectInput  (task,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput (task,  1, coutput1);

  return task;
}
