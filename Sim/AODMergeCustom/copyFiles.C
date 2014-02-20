/*
 *  copyFiles.C
 *  aliroot_dev
 *
 *  Created by philippe pillot on 18/08/11.
 *  Copyright 2011 SUBATECH. All rights reserved.
 *
 */

//______________________________________________________________________________
void copyFiles(TString gridLocation)
{
  /// copy all files necessary to run the simulation into "homeDir/gridLocation/" directory
  
  if (gridLocation.IsNull()) {
    Error("copyFiles","you must provide the grid location where to copy the files");
    return;
  }
  
  TString localDir = "$WORK/Macros/Sim";
  TString targetDir = "/alice/cern.ch/user/p/ppillot/";
  targetDir += gridLocation.Data();
  
  // connect to alien
  if (!TGrid::Connect("alien://")) {
    Error("copyFiles","cannot connect to grid");
    return;
  }
  
  if (!DirectoryExists(targetDir)) {
    Error("copyFiles", Form("directory %s does not exist", targetDir.Data()));
    return;
  }
  
  gSystem->Exec(Form("alien_cp file:%s/AODtrainCustom.C alien://%s/AODtrainCustom.C",
                     localDir.Data(), targetDir.Data()));
  gSystem->Exec(Form("alien_cp file:%s/AODMergeCustom/AOD_merge.jdl alien://%s/AOD_merge.jdl",
                     localDir.Data(), targetDir.Data()));
  gSystem->Exec(Form("alien_cp file:%s/AODMergeCustom/AOD_merge_final.jdl alien://%s/AOD_merge_final.jdl",
                     localDir.Data(), targetDir.Data()));
  gSystem->Exec(Form("alien_cp file:%s/AODMerge/validation_merge.sh alien://%s/validation_merge.sh",
                     localDir.Data(), targetDir.Data()));
/*  gSystem->Exec(Form("alien_cp file:$DEV/aliroot/PWG/muondep/AddTaskMuonRefit.C alien://%s/AddTaskMuonRefit.C",
                     targetDir.Data()));
  gSystem->Exec(Form("alien_cp file:$DEV/aliroot/PWG/muondep/AliAnalysisTaskMuonRefit.cxx alien://%s/AliAnalysisTaskMuonRefit.cxx",
                     targetDir.Data()));
  gSystem->Exec(Form("alien_cp file:$DEV/aliroot/PWG/muondep/AliAnalysisTaskMuonRefit.h alien://%s/AliAnalysisTaskMuonRefit.h",
                     targetDir.Data()));
*/  gSystem->Exec(Form("alien_cp file:$DEV/aliroot/PWG/muondep/AddTaskESDMCLabelAddition.C alien://%s/AddTaskESDMCLabelAddition.C",
                     targetDir.Data()));
  gSystem->Exec(Form("alien_cp file:$DEV/aliroot/PWG/muondep/AliAnalysisTaskESDMCLabelAddition.cxx alien://%s/AliAnalysisTaskESDMCLabelAddition.cxx",
                     targetDir.Data()));
  gSystem->Exec(Form("alien_cp file:$DEV/aliroot/PWG/muondep/AliAnalysisTaskESDMCLabelAddition.h alien://%s/AliAnalysisTaskESDMCLabelAddition.h",
                     targetDir.Data()));
  gSystem->Exec(Form("alien_cp file:$DEV/aliroot/PWGPP/MUON/dep/AddTaskMuonPerformance.C alien://%s/AddTaskMuonPerformance.C",
                     targetDir.Data()));
  gSystem->Exec(Form("alien_cp file:$DEV/aliroot/PWGPP/MUON/dep/AliAnalysisTaskMuonPerformance.cxx alien://%s/AliAnalysisTaskMuonPerformance.cxx",
                     targetDir.Data()));
  gSystem->Exec(Form("alien_cp file:$DEV/aliroot/PWGPP/MUON/dep/AliAnalysisTaskMuonPerformance.h alien://%s/AliAnalysisTaskMuonPerformance.h",
                     targetDir.Data()));
  gSystem->Exec(Form("alien_cp file:$DEV/aliroot/PWGPP/MUON/dep/AddTaskMUONTrackingEfficiency_old.C alien://%s/AddTaskMUONTrackingEfficiency_old.C",
                     targetDir.Data()));
  gSystem->Exec(Form("alien_cp file:$DEV/aliroot/PWGPP/MUON/dep/AliAnalysisTaskMuonTrackingEff_old.cxx alien://%s/AliAnalysisTaskMuonTrackingEff_old.cxx",
                     targetDir.Data()));
  gSystem->Exec(Form("alien_cp file:$DEV/aliroot/PWGPP/MUON/dep/AliAnalysisTaskMuonTrackingEff_old.h alien://%s/AliAnalysisTaskMuonTrackingEff_old.h",
                     targetDir.Data()));
  gSystem->Exec(Form("alien_cp file:$DEV/aliroot/PWGPP/MUON/dep/AddTaskMUONTrackingEfficiency.C alien://%s/AddTaskMUONTrackingEfficiency.C",
                     targetDir.Data()));
  gSystem->Exec(Form("alien_cp file:$DEV/aliroot/PWGPP/MUON/dep/AliAnalysisTaskMuonTrackingEff.cxx alien://%s/AliAnalysisTaskMuonTrackingEff.cxx",
                     targetDir.Data()));
  gSystem->Exec(Form("alien_cp file:$DEV/aliroot/PWGPP/MUON/dep/AliAnalysisTaskMuonTrackingEff.h alien://%s/AliAnalysisTaskMuonTrackingEff.h",
                     targetDir.Data()));
  gSystem->Exec(Form("alien_cp file:$WORK/Macros/MuonPhysics/AddTaskMuonPhysics.C alien://%s/AddTaskMuonPhysics.C",
                     targetDir.Data()));
  gSystem->Exec(Form("alien_cp file:$WORK/Macros/MuonPhysics/AliAnalysisTaskMuonPhysics.cxx alien://%s/AliAnalysisTaskMuonPhysics.cxx",
                     targetDir.Data()));
  gSystem->Exec(Form("alien_cp file:$WORK/Macros/MuonPhysics/AliAnalysisTaskMuonPhysics.h alien://%s/AliAnalysisTaskMuonPhysics.h",
                     targetDir.Data()));
  
}

//______________________________________________________________________________
Bool_t DirectoryExists(const char *dirname)
{
  // Returns true if directory exists. Can be also a path.
  if (!gGrid) return kFALSE;
  // Check if dirname is a path
  TString dirstripped = dirname;
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

