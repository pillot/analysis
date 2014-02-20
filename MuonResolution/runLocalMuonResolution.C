/*
 *  runLocalMuonResolution.C
 *  aliroot
 *
 *  Created by Philippe Pillot on 25/06/10.
 *  Copyright 2010 SUBATECH. All rights reserved.
 *
 */

void LoadAlirootLocally(TString& extraLibs);

//______________________________________________________________________________
void runLocalMuonResolution(TString smode = "saf", TString inputFileName = "/MUON/laphecet/SIM_MuTune0_LHC13f_muon_calo_refit_ESD_000196474",
			    TString rootVersion = "v5-34-05", TString alirootVersion = "v5-04-65-AN", Int_t nSteps = 5,
			    Bool_t selectPhysics = kFALSE, Bool_t selectTrigger = kFALSE, Bool_t matchTrig = kTRUE,
			    Bool_t applyAccCut = kTRUE, Double_t minMomentum = 0., Bool_t correctForSystematics = kTRUE,
			    Int_t extrapMode = 1, Bool_t shiftHalfCh = kFALSE, Bool_t shiftDE = kFALSE, Int_t nevents = 1234567890)
{
  /// Compute the cluster resolution by studying cluster-track residual, deconvoluting from track resolution
  /// if extrapMode == 0: extrapolate from the closest cluster
  /// if extrapMode == 1: extrapolate from the previous cluster except between stations 2-3-4
  /// if correctForSystematics == kTRUE: the systematic shifts of the residuals is included in the resolution
  
  // copy files needed for this analysis
  CopyFileLocally();
  
  // Load libraries locally
  TString extraLibs = "RAWDatabase:CDB:STEER:MUONcore:MUONmapping:MUONcalib:MUONgeometry:MUONtrigger:MUONraw:MUONbase:MUONrec:CORRFW:PWGmuon";
  LoadAlirootLocally(extraLibs);
  
  // compile analysis macro locally
  gROOT->LoadMacro("MuonResolution.C++g");
  MuonResolution(smode, inputFileName, rootVersion, alirootVersion, nSteps, selectPhysics, selectTrigger, matchTrig,
		 applyAccCut, minMomentum, correctForSystematics, extrapMode, shiftHalfCh, shiftDE, nevents, extraLibs);
  
}

//______________________________________________________________________________
void CopyFileLocally()
{
  /// Copy files needed for this analysis
  
  TString path1("$WORK/Macros/MuonResolution");
  TString path2("$WORK/aliroot/PWGPP/MUON/dep");
  TObjArray fileList(100);
  fileList.SetOwner();
  fileList.AddLast(new TObjString("runLocalMuonResolution.C"));
  fileList.AddLast(new TObjString("MuonResolution.C"));
  fileList.AddLast(new TObjString("AliAnalysisTaskMuonResolution.cxx"));
  fileList.AddLast(new TObjString("AliAnalysisTaskMuonResolution.h"));
  fileList.AddLast(new TObjString("AddTaskMuonResolution.C"));
  
  TIter nextFile(&fileList);
  TObjString *file;
  char overwrite = '\0';
  while ((file = static_cast<TObjString*>(nextFile()))) {
    
    if (overwrite != 'a') {
      
      overwrite = '\0';
      
      if (!gSystem->AccessPathName(file->GetName())) {
	
	while (overwrite != 'y' && overwrite != 'n' && overwrite != 'a' && overwrite != 'k') {
	  cout<<Form("file %s exist in current directory. Overwrite? [y=yes, n=no, a=all, k=keep all] ",file->GetName())<<flush;
	  cin>>overwrite;
	}
	
      } else overwrite = 'y';
      
    }
    
    if (overwrite == 'y' || overwrite == 'a') {
      
      if (!gSystem->AccessPathName(Form("%s/%s", gSystem->ExpandPathName(path1.Data()), file->GetName())))
	gSystem->Exec(Form("cp %s/%s %s", path1.Data(), file->GetName(), file->GetName()));
      else
	gSystem->Exec(Form("cp %s/%s %s", path2.Data(), file->GetName(), file->GetName()));
      
    } else if (overwrite == 'k') break;
    
  }
  
}

//______________________________________________________________________________
void LoadAlirootLocally(TString& extraLibs)
{
  /// Load libraries locally
  
  // Load common libraries
  gSystem->Load("libVMC");
  gSystem->Load("libTree.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libMinuit.so");
  gSystem->Load("libXMLParser.so");
  gSystem->Load("libGui.so");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  
  // Load additional libraries
  gSystem->Load("libProofPlayer");
  TObjArray* libs = extraLibs.Tokenize(":");
  for (Int_t i = 0; i < libs->GetEntriesFast(); i++)
    gSystem->Load(Form("lib%s",static_cast<TObjString*>(libs->UncheckedAt(i))->GetName()));
  delete libs;
  
  // Load lib for final mchview display
  gSystem->Load("libMUONgraphics");
  
  // Use AliRoot includes and compile our task
  gROOT->ProcessLine(".include $ALICE_ROOT/ANALYSIS/macros");
  gROOT->ProcessLine(".include $ALICE_ROOT/include");
  gROOT->ProcessLine(".include $ALICE_ROOT/MUON");
  gROOT->ProcessLine(".include $ALICE_ROOT/MUON/mapping");
  gROOT->ProcessLine(".include $ALICE_ROOT/PWG/muon");
  gROOT->LoadMacro("AliAnalysisTaskMuonResolution.cxx++g");
  
}

