#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include "TSystem.h"
#include "TString.h"
#include "TObjString.h"
#include "TFile.h"
#include "TMath.h"
#include "AliMuonAccEffSubmitter.h"
#endif

enum eParams {
  kLHC13deMUL1,
  kLHC16nMSL0
};

//_________________________________________________
Bool_t SetupGenParams(AliMuonAccEffSubmitter& sub, eParams settings)
{
  /// set the generation parameters
  
  switch (settings) {
      
    case kLHC13deMUL1: // LHC13de tuning 1 for CMUL7
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLEBEN_PTMIN","0.3");
      
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLEBEN_PT_P0","371.665");
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLEBEN_PT_P1","0.845642");
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLEBEN_PT_P2","0.56192");
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLEBEN_PT_P3","9.34859");
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLEBEN_PT_P4","0.000474519");
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLEBEN_PT_P5","-0.851091");
      
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLEBEN_Y_P0","0.777922");
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLEBEN_Y_P1","-0.0184202");
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLEBEN_Y_P2","0.00107081");
      break;
      
    case kLHC16nMSL0: // LHC16n tuning 0 for CMSL7
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLEBEN_PTMIN","0.3");
      
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLEBEN_PT_P0","135.137");
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLEBEN_PT_P1","0.555323");
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLEBEN_PT_P2","0.578374");
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLEBEN_PT_P3","10.1345");
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLEBEN_PT_P4","0.000232233");
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLEBEN_PT_P5","-0.924726");
      
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLEBEN_Y_P0","1.95551");
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLEBEN_Y_P1","-0.104761");
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLEBEN_Y_P2","0.00311324");
      break;
      
    default:
      printf("Error: unknown generator settings\n");
      return kFALSE;
      break;
  }
  
  return kTRUE;
  
}

//_________________________________________________
Float_t GetVtxSigmaXY(Float_t energy)
{
  /// Vertex transverse dispersion
  
  Float_t betast  = 3.5;                      // beta* [m]
  Float_t eps     = 3.75e-6;                   // emittance [m]
  Float_t gamma   = energy / 2.0 / 0.938272;  // relativistic gamma [1]
  return TMath::Sqrt(eps * betast / gamma) / TMath::Sqrt(2.) * 100.;  // [cm]
  
}

//_________________________________________________
void AddMuonPhysicsTask(AliMuonAccEffSubmitter &sub)
{
  /// add the task AliAnalysisTaskMuonPhysics
  
  // copy the task locally (this to avoid a call to AliMuonGridSubmitter::CleanLocal to erase the original!)
  gSystem->Exec("cp $WORK/Macros/Sim/AddExtraTasks.C .");
  gSystem->Exec("cp $WORK/Macros/MuonPhysics/AddTaskMuonPhysics.C .");
  gSystem->Exec("cp $WORK/Macros/MuonPhysics/AliAnalysisTaskMuonPhysics.cxx .");
  gSystem->Exec("cp $WORK/Macros/MuonPhysics/AliAnalysisTaskMuonPhysics.h .");
  
  // add the files to the list for submission
  sub.AddToLocalFileList("AddExtraTasks.C", kTRUE);
  sub.AddToLocalFileList("AddTaskMuonPhysics.C", kTRUE);
  sub.AddToLocalFileList("AliAnalysisTaskMuonPhysics.cxx", kTRUE);
  sub.AddToLocalFileList("AliAnalysisTaskMuonPhysics.h", kTRUE);
  
  // tell the AOD train to configure the task
  sub.SetVar("VAR_EXTRATASKS_CONFIGMACRO","\"AddExtraTasks.C\"");
  
}

//_________________________________________________
void SubmitSingleMu(TString runMode, TString runList, TString outDir)
{
  /// example of outDir = "/alice/cern.ch/user/p/ppillot/Sim/LHC15n/muTune0CMSL7"
  
  gSystem->Load("libpythia6");     // needed when testing the compilation of GenParamCustomSingleBen
  AliMuonAccEffSubmitter sub("GenParamCustomSingleBen",kFALSE);
  sub.ShouldOverwriteFiles(true);
  
  // generator config
  if (!SetupGenParams(sub, kLHC13deMUL1)) return;
  sub.SetVar("VAR_USE_ITS_RECO","0");
  sub.SetVar("VAR_USE_MC_VERTEX","1");
  Float_t sigmaxy = GetVtxSigmaXY(13000.);
  sub.SetVar("VAR_VERTEX_SIGMA_X",Form("%f",sigmaxy));
  sub.SetVar("VAR_VERTEX_SIGMA_Y",Form("%f",sigmaxy));

  // OCDB config
  sub.UseOCDBSnapshots(kFALSE);
  sub.SetVar("VAR_USE_RAW_ALIGN","0");
  sub.SetVar("VAR_SIM_ALIGNDATA","\"alien://folder=/alice/simulation/2008/v4-15-Release/Ideal\"");
  sub.SetVar("VAR_REC_ALIGNDATA","\"alien://folder=/alice/simulation/2008/v4-15-Release/Residual\"");
  
  // AODtrain config
  sub.SetVar("VAR_MUONMCMODE","1");
  sub.SetVar("VAR_EFFTASK_PTMIN","0.5");
  AddMuonPhysicsTask(sub);
  sub.UseAODMerging(kTRUE);
  sub.SetVar("VAR_AOD_MERGE_FILES","\"Merged.QA.Data.root,AnalysisResults.root\"");
  
  // jdl config
  sub.SetAliPhysicsVersion("VO_ALICE@AliPhysics::vAN-20161114-1");
  sub.SetCustomOutFiles("log_archive.zip:stderr,stdout@disk=1", "root_archive.zip:AliAOD.Muons.root,Merged.QA.Data.root,AnalysisResults.root@disk=2");

  // input/output config
  gSystem->ExpandPathName(runList);
  sub.SetRunList(runList.Data());
  sub.SetRemoteDir(outDir.Data());
  sub.SetMergedDir(outDir.Data());
  
  // run config
  runMode.ToUpper();
  if ( runMode == "LOCALTEST" ) sub.MakeNofEventsFixed(100);
  else sub.MakeNofEventsPropToTriggerCount("CMUL7-B-NOPF-MUFAST",2.);
  sub.SetMaxEventsPerChunk(5000);
  sub.Print();
  
  // produce or merge
  if (runMode.Contains("MERGE")) {
    TObjArray* mergeModeStep = runMode.Tokenize(":");
    if (mergeModeStep->GetEntriesFast() != 2) {
      printf("Error: merging format = MODE:STEP where MODE = MERGE or TESTMERGE and STEP is the step number (0=final)\n");
      return;
    }
    Int_t stage = static_cast<TObjString*>(mergeModeStep->UncheckedAt(1))->String().Atoi();
    Bool_t dryRun = (static_cast<TObjString*>(mergeModeStep->UncheckedAt(0))->String().Contains("TEST"));
    sub.Merge(stage, dryRun);
  } else {
    sub.Run(runMode.Data());
  }
  
}

