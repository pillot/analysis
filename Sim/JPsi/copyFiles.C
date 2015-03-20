/*
 *  copyFiles.C
 *  aliphysics/dev/src_dev
 *
 *  Created by philippe pillot on 18/08/11.
 *  Copyright 2011 SUBATECH. All rights reserved.
 *
 */

//______________________________________________________________________________
void copyFiles(TString targetDir)
{
  /// copy all files necessary to run the simulation into "targetDir/" directory
  
  if (targetDir.IsNull()) {
    Error("copyFiles","you must provide the location where to copy the files");
    return;
  }
  
  TString localDir = "/Users/pillot/Work/Alice/Macros/Sim";
  
  if (targetDir.BeginsWith("alien://")) { // grid copy
    
    // connect to alien
    if (!TGrid::Connect("alien://")) {
      Error("copyFiles","cannot connect to grid");
      return;
    }
    
    if (!DirectoryExists(targetDir)) {
      Error("copyFiles", Form("directory %s does not exist", targetDir.Data()));
      return;
    }
    
    gSystem->Exec(Form("alien_cp file:%s/rec.C %s/rec.C", localDir.Data(), targetDir.Data()));
    gSystem->Exec(Form("alien_cp file:%s/sim.C %s/sim.C", localDir.Data(), targetDir.Data()));
    gSystem->Exec(Form("alien_cp file:%s/tag.C %s/tag.C", localDir.Data(), targetDir.Data()));
    gSystem->Exec(Form("alien_cp file:%s/validation.sh %s/validation.sh", localDir.Data(), targetDir.Data()));
    gSystem->Exec(Form("alien_cp file:%s/CheckESD.C %s/CheckESD.C", localDir.Data(), targetDir.Data()));
    gSystem->Exec(Form("alien_cp file:%s/CheckAOD.C %s/CheckAOD.C", localDir.Data(), targetDir.Data()));
    gSystem->Exec(Form("alien_cp file:%s/JPsi/Config.C %s/Config.C", localDir.Data(), targetDir.Data()));
    gSystem->Exec(Form("alien_cp file:%s/JPsi/run.jdl %s/run.jdl", localDir.Data(), targetDir.Data()));
    gSystem->Exec(Form("alien_cp file:%s/JPsi/simrun.C %s/simrun.C", localDir.Data(), targetDir.Data()));
    //  gSystem->Exec(Form("alien_cp file:%s/AODtrain.C %s/AODtrain.C", localDir.Data(), targetDir.Data()));
    gSystem->Exec(Form("alien_cp file:%s/AODtrainCustom.C %s/AODtrainCustom.C", localDir.Data(), targetDir.Data()));
    gSystem->Exec(Form("alien_cp file:$WORK/aliphysics/dev/src/PWG/muondep/AddTaskMuonRefit.C %s/AddTaskMuonRefit.C", targetDir.Data()));
    gSystem->Exec(Form("alien_cp file:$WORK/aliphysics/dev/src/PWG/muondep/AliAnalysisTaskMuonRefit.cxx %s/AliAnalysisTaskMuonRefit.cxx", targetDir.Data()));
    gSystem->Exec(Form("alien_cp file:$WORK/aliphysics/dev/src/PWG/muondep/AliAnalysisTaskMuonRefit.h %s/AliAnalysisTaskMuonRefit.h", targetDir.Data()));
    gSystem->Exec(Form("alien_cp file:$WORK/aliphysics/dev/src/PWG/muondep/AddTaskESDMCLabelAddition.C %s/AddTaskESDMCLabelAddition.C", targetDir.Data()));
    gSystem->Exec(Form("alien_cp file:$WORK/aliphysics/dev/src/PWG/muondep/AliAnalysisTaskESDMCLabelAddition.cxx %s/AliAnalysisTaskESDMCLabelAddition.cxx", targetDir.Data()));
    gSystem->Exec(Form("alien_cp file:$WORK/aliphysics/dev/src/PWG/muondep/AliAnalysisTaskESDMCLabelAddition.h %s/AliAnalysisTaskESDMCLabelAddition.h", targetDir.Data()));
    gSystem->Exec(Form("alien_cp file:$WORK/aliphysics/dev/src/PWGPP/MUON/dep/AddTaskMuonPerformance.C %s/AddTaskMuonPerformance.C", targetDir.Data()));
    gSystem->Exec(Form("alien_cp file:$WORK/aliphysics/dev/src/PWGPP/MUON/dep/AliAnalysisTaskMuonPerformance.cxx %s/AliAnalysisTaskMuonPerformance.cxx", targetDir.Data()));
    gSystem->Exec(Form("alien_cp file:$WORK/aliphysics/dev/src/PWGPP/MUON/dep/AliAnalysisTaskMuonPerformance.h %s/AliAnalysisTaskMuonPerformance.h", targetDir.Data()));
    gSystem->Exec(Form("alien_cp file:$WORK/aliphysics/dev/src/PWGPP/MUON/dep/AddTaskMUONTrackingEfficiency.C %s/AddTaskMUONTrackingEfficiency.C", targetDir.Data()));
    gSystem->Exec(Form("alien_cp file:$WORK/aliphysics/dev/src/PWGPP/MUON/dep/AliAnalysisTaskMuonTrackingEff.cxx %s/AliAnalysisTaskMuonTrackingEff.cxx", targetDir.Data()));
    gSystem->Exec(Form("alien_cp file:$WORK/aliphysics/dev/src/PWGPP/MUON/dep/AliAnalysisTaskMuonTrackingEff.h %s/AliAnalysisTaskMuonTrackingEff.h", targetDir.Data()));
    
  } else { // local copy
    
    if (gSystem->AccessPathName(targetDir.Data())) {
      Error("copyFiles", Form("directory %s does not exist", targetDir.Data()));
      return;
    }
    
    gSystem->Exec(Form("cp %s/rec.C %s/rec.C", localDir.Data(), targetDir.Data()));
    gSystem->Exec(Form("cp %s/sim.C %s/sim.C", localDir.Data(), targetDir.Data()));
    gSystem->Exec(Form("cp %s/tag.C %s/tag.C", localDir.Data(), targetDir.Data()));
    gSystem->Exec(Form("cp %s/validation.sh %s/validation.sh", localDir.Data(), targetDir.Data()));
    gSystem->Exec(Form("cp %s/CheckESD.C %s/CheckESD.C", localDir.Data(), targetDir.Data()));
    gSystem->Exec(Form("cp %s/CheckAOD.C %s/CheckAOD.C", localDir.Data(), targetDir.Data()));
    gSystem->Exec(Form("cp %s/JPsi/Config.C %s/Config.C", localDir.Data(), targetDir.Data()));
    gSystem->Exec(Form("cp %s/JPsi/run.jdl %s/run.jdl", localDir.Data(), targetDir.Data()));
    gSystem->Exec(Form("cp %s/JPsi/simrun.C %s/simrun.C", localDir.Data(), targetDir.Data()));
    //  gSystem->Exec(Form("cp %s/AODtrain.C %s/AODtrain.C", localDir.Data(), targetDir.Data()));
    gSystem->Exec(Form("cp %s/AODtrainCustom.C %s/AODtrainCustom.C", localDir.Data(), targetDir.Data()));
    gSystem->Exec(Form("cp $WORK/aliphysics/dev/src/PWG/muondep/AddTaskMuonRefit.C %s/AddTaskMuonRefit.C", targetDir.Data()));
    gSystem->Exec(Form("cp $WORK/aliphysics/dev/src/PWG/muondep/AliAnalysisTaskMuonRefit.cxx %s/AliAnalysisTaskMuonRefit.cxx", targetDir.Data()));
    gSystem->Exec(Form("cp $WORK/aliphysics/dev/src/PWG/muondep/AliAnalysisTaskMuonRefit.h %s/AliAnalysisTaskMuonRefit.h", targetDir.Data()));
    gSystem->Exec(Form("cp $WORK/aliphysics/dev/src/PWG/muondep/AddTaskESDMCLabelAddition.C %s/AddTaskESDMCLabelAddition.C", targetDir.Data()));
    gSystem->Exec(Form("cp $WORK/aliphysics/dev/src/PWG/muondep/AliAnalysisTaskESDMCLabelAddition.cxx %s/AliAnalysisTaskESDMCLabelAddition.cxx", targetDir.Data()));
    gSystem->Exec(Form("cp $WORK/aliphysics/dev/src/PWG/muondep/AliAnalysisTaskESDMCLabelAddition.h %s/AliAnalysisTaskESDMCLabelAddition.h", targetDir.Data()));
    gSystem->Exec(Form("cp $WORK/aliphysics/dev/src/PWGPP/MUON/dep/AddTaskMuonPerformance.C %s/AddTaskMuonPerformance.C", targetDir.Data()));
    gSystem->Exec(Form("cp $WORK/aliphysics/dev/src/PWGPP/MUON/dep/AliAnalysisTaskMuonPerformance.cxx %s/AliAnalysisTaskMuonPerformance.cxx", targetDir.Data()));
    gSystem->Exec(Form("cp $WORK/aliphysics/dev/src/PWGPP/MUON/dep/AliAnalysisTaskMuonPerformance.h %s/AliAnalysisTaskMuonPerformance.h", targetDir.Data()));
    gSystem->Exec(Form("cp $WORK/aliphysics/dev/src/PWGPP/MUON/dep/AddTaskMUONTrackingEfficiency.C %s/AddTaskMUONTrackingEfficiency.C", targetDir.Data()));
    gSystem->Exec(Form("cp $WORK/aliphysics/dev/src/PWGPP/MUON/dep/AliAnalysisTaskMuonTrackingEff.cxx %s/AliAnalysisTaskMuonTrackingEff.cxx", targetDir.Data()));
    gSystem->Exec(Form("cp $WORK/aliphysics/dev/src/PWGPP/MUON/dep/AliAnalysisTaskMuonTrackingEff.h %s/AliAnalysisTaskMuonTrackingEff.h", targetDir.Data()));
    
  }
  
}

//______________________________________________________________________________
Bool_t DirectoryExists(const char *dirname)
{
  // Returns true if directory exists. Can be also a path.
  if (!gGrid) return kFALSE;
  // Check if dirname is a path
  TString dirstripped = dirname;
  dirstripped.ReplaceAll("alien://","");
  dirstripped = dirstripped.Strip();
  dirstripped = dirstripped.Strip(TString::kTrailing, '/');
  TString dir = gSystem->BaseName(dirstripped);
  dir += "/";
  TString path = gSystem->DirName(dirstripped);
  TGridResult *res = gGrid->Ls(path, "-F");
  if (!res) return kFALSE;
  TIter next(res);
  TMap *map;
  TObject *obj;
  while ((map=dynamic_cast<TMap*>(next()))) {
    obj = map->GetValue("name");
    if (!obj) break;
    if (dir == obj->GetName()) {
      delete res;
      return kTRUE;
    }
  }
  delete res;
  return kFALSE;
}      

