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
  kLHC13fMUL1,
  kLHC16nMSL0,
  kLHC16rMUL1,
  kLHC16rMUL2,
  kLHC16rMSH1,
  kLHC16sMUL1,
  kLHC16sMUL2,
  kLHC16sMSH1
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
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLEBEN_Y_P2","-0.00107081");
      break;
      
    case kLHC13fMUL1: // LHC13f tuning 1 for CMUL7
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLEBEN_PTMIN","0.3");
      
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLEBEN_PT_P0","455.614");
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLEBEN_PT_P1","0.942071");
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLEBEN_PT_P2","0.706755");
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLEBEN_PT_P3","8.69856");
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLEBEN_PT_P4","0.000168775");
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLEBEN_PT_P5","-0.925487");
      
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLEBEN_Y_P0","1.29511");
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLEBEN_Y_P1","-0.0767846");
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLEBEN_Y_P2","0.00176313");
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
      
    case kLHC16rMUL1: // LHC16r tuning 1 for CMUL7 (with pT > 0.5 GeV/c)
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLE_PTMIN","0.3");
      
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLE_PT_P0","86.4396");
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLE_PT_P1","0.952196");
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLE_PT_P2","0.596445");
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLE_PT_P3","8.55775");
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLE_PT_P4","0.000535784");
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLE_PT_P5","-0.794826");
      
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLE_Y_P0","-7.09298");
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLE_Y_P1","1.0469");
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLE_Y_P2","0.324374");
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLE_Y_P3","0.032018");
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLE_Y_P4","0.");
      break;
      
    case kLHC16rMUL2: // LHC16r tuning 2 for CMUL7 (with pT > 1 GeV/c)
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLE_PTMIN","0.8");
      
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLE_PT_P0","409.198");
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLE_PT_P1","0.945457");
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLE_PT_P2","0.577475");
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLE_PT_P3","8.87429");
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLE_PT_P4","0.00050164");
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLE_PT_P5","-0.809412");
      
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLE_Y_P0","-23.9284");
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLE_Y_P1","1.27419");
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLE_Y_P2","0.586072");
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLE_Y_P3","0.118576");
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLE_Y_P4","0.00895485");
      break;
      
    case kLHC16rMSH1: // LHC16r tuning 1 for CMSH7 (with pT > 5 GeV/c)
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLE_PTMIN","4.7");
      
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLE_PT_P0","17167.7");
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLE_PT_P1","0.76396");
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLE_PT_P2","0.686429");
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLE_PT_P3","7.84528");
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLE_PT_P4","0.000306233");
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLE_PT_P5","-0.877303");
      
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLE_Y_P0","1.45681");
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLE_Y_P1","0.");
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLE_Y_P2","-0.083698");
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLE_Y_P3","0.");
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLE_Y_P4","0.00201026");
      break;
      
    case kLHC16sMUL1: // LHC16s tuning 1 for CMUL7  (with pT > 0.5 GeV/c)
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLE_PTMIN","0.3");

      sub.SetVar("VAR_GENPARAMCUSTOMSINGLE_PT_P0","77.7982");
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLE_PT_P1","0.97131");
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLE_PT_P2","0.686076");
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLE_PT_P3","8.40771");
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLE_PT_P4","0.000262219");
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLE_PT_P5","-0.883664");
      
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLE_Y_P0","-6.87096");
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLE_Y_P1","1.10931");
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLE_Y_P2","0.357909");
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLE_Y_P3","0.0364195");
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLE_Y_P4","0.");
      break;
      
    case kLHC16sMUL2: // LHC16s tuning 2 for CMUL7 (with pT > 1 GeV/c)
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLE_PTMIN","0.8");
      
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLE_PT_P0","483.244");
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLE_PT_P1","0.977042");
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLE_PT_P2","0.67239");
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLE_PT_P3","8.66863");
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLE_PT_P4","0.000266572");
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLE_PT_P5","-0.903221");
      
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLE_Y_P0","-15.6829");
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLE_Y_P1","1.2897");
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLE_Y_P2","0.58274");
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLE_Y_P3","0.114462");
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLE_Y_P4","0.00835864");
      break;
      
    case kLHC16sMSH1: // LHC16s tuning 1 for CMSH7 (with pT > 5 GeV/c)
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLE_PTMIN","4.7");
      
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLE_PT_P0","64161.1");
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLE_PT_P1","1.96091");
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLE_PT_P2","1.12885");
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLE_PT_P3","5.52326");
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLE_PT_P4","0.0684062");
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLE_PT_P5","-2.39476");
      
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLE_Y_P0","2.07495");
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLE_Y_P1","0.");
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLE_Y_P2","-0.102246");
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLE_Y_P3","0.");
      sub.SetVar("VAR_GENPARAMCUSTOMSINGLE_Y_P4","0.00276377");
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
  /// example of outDir = "/alice/cern.ch/user/p/ppillot/Sim/LHC16r/muTune0"
  
  gSystem->Load("libpythia6");     // needed when testing the compilation of GenParamCustomSingleBen
  AliMuonAccEffSubmitter sub("GenParamCustomSingle",kFALSE);
  sub.ShouldOverwriteFiles(true);
  
  // generator config
  if (!SetupGenParams(sub, kLHC16sMSH1)) return;
  sub.SetVar("VAR_USE_ITS_RECO","0");
  sub.SetVar("VAR_USE_MC_VERTEX","1");
  Float_t sigmaxy = GetVtxSigmaXY(8000.);
  sub.SetVar("VAR_VERTEX_SIGMA_X",Form("%f",sigmaxy));
  sub.SetVar("VAR_VERTEX_SIGMA_Y",Form("%f",sigmaxy));

  // OCDB config
  sub.UseOCDBSnapshots(kFALSE);
  sub.SetVar("VAR_OCDB_PATH","\"raw://\"");
  sub.SetVar("VAR_USE_RAW_ALIGN","1");
//  sub.SetVar("VAR_SIM_ALIGNDATA","\"alien://folder=/alice/simulation/2008/v4-15-Release/Ideal\"");
//  sub.SetVar("VAR_REC_ALIGNDATA","\"alien://folder=/alice/simulation/2008/v4-15-Release/Residual\"");
//  sub.SetVar("VAR_SIM_ALIGNDATA","\"alien://folder=/alice/cern.ch/user/j/jcastill/pp16wrk/LHC16_mcp1vsrealv2_tr_MisAlignCDB\"");
  sub.SetVar("VAR_SIM_ALIGNDATA","\"alien://folder=/alice/simulation/2008/v4-15-Release/Full\"");
  sub.SetVar("VAR_REC_ALIGNDATA","\"alien://folder=/alice/data/2016/OCDB\"");
  
  // AODtrain config
  sub.SetVar("VAR_MUONMCMODE","1");
//  sub.SetVar("VAR_EFFTASK_PTMIN","0.5");
  AddMuonPhysicsTask(sub);
  sub.UseAODMerging(kTRUE);
  sub.SetVar("VAR_AOD_MERGE_FILES","\"AnalysisResults.root\"");
  
  // jdl config
  sub.SetAliPhysicsVersion("VO_ALICE@AliPhysics::vAN-20170102-1");
  sub.SetCustomOutFiles("log_archive.zip:stderr,stdout@disk=1", "root_archive.zip:AliAOD.Muons.root,AnalysisResults.root@disk=2");
//  sub.SetCustomOutFiles("log_archive.zip:stderr,stdout@disk=1", "root_archive.zip:AnalysisResults.root@disk=2");

  // input/output config
  gSystem->ExpandPathName(runList);
  sub.SetRunList(runList.Data());
  sub.SetRemoteDir(outDir.Data());
  sub.SetMergedDir(outDir.Data());
  
  // run config
  runMode.ToUpper();
  if ( runMode == "LOCALTEST" ) sub.MakeNofEventsFixed(100);
  else sub.MakeNofEventsPropToTriggerCount("CMUL7-B-NOPF-MUFAST",0.03);
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

