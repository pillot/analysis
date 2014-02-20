#ifndef __CINT__
#include <AliAnalysisManager.h>
#include <AliMultiInputEventHandler.h>
#include <ANALYSIS/EventMixing/AliMixEventPool.h>
#include <ANALYSIS/EventMixing/AliMixEventCutObj.h>
#include <PWGLF/RESONANCES/AliRsnAnalysisTask.h>
#include <PWGLF/RESONANCES/AliRsnMiniAnalysisTask.h>
#endif

void AddMixingHandler ( AliMultiInputEventHandler* multiInputHandler,TString format = "esd", Bool_t useMC = kFALSE, TString opts = "" ) {

  const Int_t bufferSize = 1;
  const Int_t mixNum = 20;  // 6
  if ( !multiInputHandler ) return;

    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    AliMixInputEventHandler *mixHandler = new AliMixInputEventHandler ( bufferSize, mixNum );
    mixHandler->SetInputHandlerForMixing ( dynamic_cast<AliMultiInputEventHandler*> ( mgr->GetInputEventHandler() ) );
    AliMixEventPool *evPool = new AliMixEventPool();

    AliMixEventCutObj *centrality = new AliMixEventCutObj(AliMixEventCutObj::kCentrality, 0., 90., 10., "V0M");
    evPool->AddCut(centrality);  
  
    // adds event pool (comment it and u will have default mixing)
    mixHandler->SetEventPool(evPool);

    mixHandler->SelectCollisionCandidates(AliVEvent::kMUSPB);
    mixHandler->DoMixIfNotEnoughEvents(kTRUE);

    multiInputHandler->AddInputEventHandler(mixHandler);

    // adds mixing info task
    gROOT->LoadMacro("AddAnalysisTaskMixInfo.C");
    AddAnalysisTaskMixInfo (opts );
    
}
