/*
 *  runLocalMuonRefitResolution.C
 *  aliroot
 *
 *  Created by Philippe Pillot on 11/11/10.
 *  Copyright 2010 SUBATECH. All rights reserved.
 *
 */

Bool_t taskInPWG3 = kFALSE;

void LoadAlirootLocally(TString& extraLibs);

//______________________________________________________________________________
void runLocalMuonRefitResolution(TString smode = "local", TString inputFileName = "AliESDs.root", Int_t nSteps = 5,
				 Bool_t selectPhysics = kFALSE, Bool_t matchTrig = kTRUE, Double_t minMomentum = 0.,
				 Bool_t correctForSystematics = kTRUE, Int_t extrapMode = 1, Int_t nevents = 1234567890)
{
  /// Compute the cluster resolution by studying cluster-track residual, deconvoluting from track resolution
  /// if extrapMode == 0: extrapolate from the closest cluster
  /// if extrapMode == 1: extrapolate from the previous cluster except between stations 2-3-4
  /// if correctForSystematics == kTRUE: the systematic shifts of the residuals is included in the resolution
  
  // copy files needed for this analysis
  TString path1("/Users/philippe/Work/Alice/Work/Data/Macro/MuonResolution");
  TString path2("/Users/philippe/Work/Alice/Work/Data/Macro/MuonPhysics");
  TString path3("$ALICE_ROOT/PWG3/muondep");
  TObjArray fileList(100);
  fileList.SetOwner();
  fileList.AddLast(new TObjString("runLocalMuonRefitResolution.C"));
  fileList.AddLast(new TObjString("MuonRefitResolution.C"));
  fileList.AddLast(new TObjString("AliAnalysisTaskMuonResolution.cxx"));
  fileList.AddLast(new TObjString("AliAnalysisTaskMuonResolution.h"));
  fileList.AddLast(new TObjString("AddTaskMuonResolution.C"));
  fileList.AddLast(new TObjString("ConfigureCuts.C"));
  fileList.AddLast(new TObjString("AliAnalysisTaskMuonRefit.cxx"));
  fileList.AddLast(new TObjString("AliAnalysisTaskMuonRefit.h"));
  TIter nextFile(&fileList);
  TObjString *file;
  char overwrite = '\0';
  while ((file = static_cast<TObjString*>(nextFile()))) {
    if (overwrite != 'a') {
      overwrite = '\0';
      if (!gSystem->AccessPathName(file->GetName())) {
	cout<<Form("file %s exist in current directory. Overwrite? [y=yes, n=no, a=all, k=keep all] ",file->GetName())<<flush;
	while (overwrite != 'y' && overwrite != 'n' && overwrite != 'a' && overwrite != 'k') cin>>overwrite;
      } else overwrite = 'y';
    }
    if (overwrite == 'y' || overwrite == 'a') {
      if (!gSystem->AccessPathName(Form("%s/%s", path1.Data(), file->GetName())))
	gSystem->Exec(Form("cp %s/%s %s", path1.Data(), file->GetName(), file->GetName()));
      else if (!gSystem->AccessPathName(Form("%s/%s", path2.Data(), file->GetName())))
	gSystem->Exec(Form("cp %s/%s %s", path2.Data(), file->GetName(), file->GetName()));
      else if (!taskInPWG3)
	gSystem->Exec(Form("cp %s/%s %s", path3.Data(), file->GetName(), file->GetName()));
    } else if (overwrite == 'k') break;
  }
  
  // Load libraries locally
  TString extraLibs;
  if (taskInPWG3) extraLibs="RAWDatabase:CDB:STEER:MUONcore:MUONmapping:MUONcalib:MUONgeometry:MUONtrigger:MUONraw:MUONbase:MUONrec:CORRFW:PWG3muondep";
  else extraLibs="RAWDatabase:CDB:STEER:MUONcore:MUONmapping:MUONcalib:MUONgeometry:MUONtrigger:MUONraw:MUONbase:MUONrec";
  LoadAlirootLocally(extraLibs);
  
  // compile analysis macro locally
  gROOT->LoadMacro("MuonRefitResolution.C++g");
  MuonRefitResolution(smode, inputFileName, nSteps, selectPhysics, matchTrig, minMomentum, correctForSystematics, extrapMode, nevents, extraLibs);
  
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
  if (taskInPWG3) gROOT->ProcessLine(".include $ALICE_ROOT/PWG3/muondep");
  else {
    gROOT->ProcessLine(".include $ALICE_ROOT/include");
    gROOT->ProcessLine(".include $ALICE_ROOT/MUON");
    gROOT->ProcessLine(".include $ALICE_ROOT/MUON/mapping");
    gROOT->LoadMacro("AliAnalysisTaskMuonResolution.cxx++g");
    gROOT->LoadMacro("AliAnalysisTaskMuonRefit.cxx++g");
  }
  
}

