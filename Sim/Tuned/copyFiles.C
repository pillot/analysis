/*
 *  copyFiles.C
 *  aliroot_dev
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
  
  TString localDir = "$WORK/Macros/Sim";
  TString cp = "cp ";
  
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
    
    cp = "alien_cp file:";
    
  } else { // local copy
    
    if (gSystem->AccessPathName(targetDir.Data())) {
      Error("copyFiles", Form("directory %s does not exist", targetDir.Data()));
      return;
    }
    
  }
  
  gSystem->Exec(Form("%s%s/rec.C %s/rec.C",
                     cp.Data(), localDir.Data(), targetDir.Data()));
  gSystem->Exec(Form("%s%s/sim.C %s/sim.C",
                     cp.Data(), localDir.Data(), targetDir.Data()));
  gSystem->Exec(Form("%s%s/tag.C %s/tag.C",
                     cp.Data(), localDir.Data(), targetDir.Data()));
  gSystem->Exec(Form("%s%s/validation.sh %s/validation.sh",
                     cp.Data(), localDir.Data(), targetDir.Data()));
  gSystem->Exec(Form("%s%s/CheckESD.C %s/CheckESD.C",
                     cp.Data(), localDir.Data(), targetDir.Data()));
  gSystem->Exec(Form("%s%s/CheckAOD.C %s/CheckAOD.C",
                     cp.Data(), localDir.Data(), targetDir.Data()));
  gSystem->Exec(Form("%s%s/Tuned/Config.C %s/Config.C",
                     cp.Data(), localDir.Data(), targetDir.Data()));
  gSystem->Exec(Form("%s%s/Tuned/MuonGenerator.C %s/MuonGenerator.C",
                     cp.Data(), localDir.Data(), targetDir.Data()));
  gSystem->Exec(Form("%s%s/Tuned/run.jdl %s/run.jdl",
                     cp.Data(), localDir.Data(), targetDir.Data()));
  gSystem->Exec(Form("%s%s/Tuned/simrun.C %s/simrun.C",
                     cp.Data(), localDir.Data(), targetDir.Data()));
  //gSystem->Exec(Form("%s%s/AODtrain.C %s/AODtrain.C",
  //                   cp.Data(), localDir.Data(), targetDir.Data()));
  gSystem->Exec(Form("%s%s/AODtrainCustom.C %s/AODtrainCustom.C",
                     cp.Data(), localDir.Data(), targetDir.Data()));
  
/*  gSystem->Exec(Form("%s$DEVPHYS/PWG/muondep/AddTaskMuonRefit.C %s/AddTaskMuonRefit.C",
                     cp.Data(), targetDir.Data()));
*/  gSystem->Exec(Form("%s$DEVPHYS/PWG/muondep/AddTaskESDMCLabelAddition.C %s/AddTaskESDMCLabelAddition.C",
                     cp.Data(), targetDir.Data()));
  gSystem->Exec(Form("%s$DEVPHYS/../build/PWG/muondep/PWGmuondep.par %s/PWGmuondep.par",
                     cp.Data(), targetDir.Data()));
  
  gSystem->Exec(Form("%s$DEVPHYS/PWGPP/MUON/dep/AddTaskMuonPerformance.C %s/AddTaskMuonPerformance.C",
                     cp.Data(), targetDir.Data()));
  gSystem->Exec(Form("%s$DEVPHYS/PWGPP/MUON/dep/AddTaskMUONTrackingEfficiency.C %s/AddTaskMUONTrackingEfficiency.C",
                     cp.Data(), targetDir.Data()));
  gSystem->Exec(Form("%s$DEVPHYS/../build/PWGPP/MUON/dep/PWGPPMUONdep.par %s/PWGPPMUONdep.par",
                     cp.Data(), targetDir.Data()));
  
  gSystem->Exec(Form("%s$DEVPHYS/PWGPP/PilotTrain/AddTaskMuonQA.C %s/AddTaskMuonQA.C",
                     cp.Data(), targetDir.Data()));
  gSystem->Exec(Form("%s$DEVPHYS/../build/PWGPP/MUON/lite/PWGPPMUONlite.par %s/PWGPPMUONlite.par",
                     cp.Data(), targetDir.Data()));
  
  gSystem->Exec(Form("%s$WORK/Macros/MuonPhysics/AddTaskMuonPhysics.C %s/AddTaskMuonPhysics.C",
                     cp.Data(), targetDir.Data()));
  gSystem->Exec(Form("%s$WORK/Macros/MuonPhysics/AliAnalysisTaskMuonPhysics.cxx %s/AliAnalysisTaskMuonPhysics.cxx",
                     cp.Data(), targetDir.Data()));
  gSystem->Exec(Form("%s$WORK/Macros/MuonPhysics/AliAnalysisTaskMuonPhysics.h %s/AliAnalysisTaskMuonPhysics.h",
                     cp.Data(), targetDir.Data()));
  
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

