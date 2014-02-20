/*
 *  runLocalMuonResolution.C
 *  aliroot
 *
 *  Created by Philippe Pillot on 25/06/10.
 *  Copyright 2010 SUBATECH. All rights reserved.
 *
 */

#if !defined(__CINT__) || defined(__MAKECINT__)
// ROOT includes
#include <fstream>
#include <TString.h>
#include <TStopwatch.h>
#include <TMultiGraph.h>
#include <TSystem.h>
#include <TChain.h>
#include <TGraphErrors.h>
#include <TProof.h>
#include <TList.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TGrid.h>
#include <TEnv.h>
#include <TROOT.h>
#include <TAlienCollection.h>
#include <TGridCollection.h>
#include <TGridResult.h>

// STEER includes
#include "AliLog.h"
#include "AliCDBManager.h"
#include "AliAnalysisManager.h"
#include "AliESDInputHandler.h"
#include "AliTagAnalysis.h"
#include "AliRunTagCuts.h"
#include "AliLHCTagCuts.h"
#include "AliDetectorTagCuts.h"
#include "AliEventTagCuts.h"

#endif

Bool_t taskInPWG3 = kFALSE;
TString aliroot="VO_ALICE@AliRoot::v4-20-11-AN";

enum {kLocal, kInteractif_xml, kInteractif_ESDList, kProof};

void LoadAlirootLocaly(TString& extraLibs);
void run(Int_t mode, Int_t iStep, Int_t nevents, TObject* input, Int_t extrapMode, Bool_t correctForSystematics,
	 Double_t minMomentum, Bool_t matchTrig, Bool_t selectPhysics, Double_t clusterResNB[10], Double_t clusterResB[10]);
Bool_t GetChamberResolution(Int_t iStep, Double_t clusterResNB[10], Double_t clusterResB[10], Double_t clusterResNBErr[10], Double_t clusterResBErr[10]);
Int_t GetMode(TString smode, TString input);
TChain* CreateChainFromCollection(const char *xmlfile);
TChain* CreateChainFromFile(const char *rootfile);
TChain* CreateChainFromESDList(const char *esdList);
TChain* CreateChain(Int_t mode, TString input);

//______________________________________________________________________________
void runLocalMuonResolution(TString smode = "local", TString inputFileName = "AliESDs.root", Int_t nSteps = 10,
			    Bool_t selectPhysics = kFALSE, Bool_t matchTrig = kFALSE, Double_t minMomentum = 0.,
			    Bool_t correctForSystematics = kTRUE, Int_t extrapMode = 1, Int_t nevents = 1234567890)
{
  /// Compute the cluster resolution by studying cluster-track residual, deconvoluting from track resolution
  /// if extrapMode == 0: extrapolate from the closest cluster
  /// if extrapMode == 1: extrapolate from the previous cluster except between stations 2-3-4
  /// if correctForSystematics == kTRUE: the systematic shifts of the residuals is included in the resolution
  
  Int_t gNEvents = nevents;
  Int_t gExtrapMode = extrapMode;
  Bool_t gCorrectForSystematics = correctForSystematics;
  Double_t gMinMomentum = minMomentum;
  Bool_t gMatchTrig = matchTrig;
  Bool_t gSelectPhysics = selectPhysics;
  
  // timer start...
  TStopwatch* localTimer = new TStopwatch;
  
  // check parameters
  nSteps = TMath::Max(nSteps,1);
  if (gExtrapMode != 0 && gExtrapMode != 1) {
    Error("runLocalMuonResolution","incorrect extrapolation mode!");
    return;
  }
  
  // Check runing mode
  Int_t gMode = GetMode(smode, inputFileName);
  if(gMode < 0){
    Error("runLocalMuonResolution","Please provide either an ESD root file a collection of ESDs or a dataset.");
    return;
  }
  
  // copy files needed for this analysis
  TString path1("/Users/philippe/Work/Alice/Work/Data/Macro/MuonResolution");
  TString path2("$ALICE_ROOT/PWG3/muondep");
  TObjArray fileList(100);
  fileList.SetOwner();
  fileList.AddLast(new TObjString("runLocalMuonResolution.C"));
  fileList.AddLast(new TObjString("AliAnalysisTaskMuonResolution.cxx"));
  fileList.AddLast(new TObjString("AliAnalysisTaskMuonResolution.h"));
  fileList.AddLast(new TObjString("AddTaskMuonResolution.C"));
  fileList.AddLast(new TObjString("ConfigureCuts.C"));
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
	gSystem->Exec(Form("cp %s/%s %s", path1.Data(), file->GetName(),file->GetName()));
      else if (!taskInPWG3)
	gSystem->Exec(Form("cp %s/%s %s", path2.Data(), file->GetName(),file->GetName()));
    } else if (overwrite == 'k') break;
  }
  
  // Load libraries locally
  TString extraLibs;
  if (taskInPWG3) extraLibs="RAWDatabase:CDB:STEER:MUONcore:MUONmapping:MUONcalib:MUONgeometry:MUONtrigger:MUONraw:MUONbase:MUONrec:CORRFW:PWG3muondep";
  else extraLibs="RAWDatabase:CDB:STEER:MUONcore:MUONmapping:MUONcalib:MUONgeometry:MUONtrigger:MUONraw:MUONbase:MUONrec";
  if (gMode != kProof) LoadAlirootLocaly(extraLibs);
  
  // check for old output file to removed
  char remove = '\0';
  if (!gSystem->Exec("ls chamberResolution_step*[0-9].root")) {
    cout<<"above files must be removed from the current directory. Delete? [y=yes, n=no] "<<flush;
    while (remove != 'y' && remove != 'n') cin>>remove;
    if (remove == 'y') gSystem->Exec("rm -f chamberResolution_step*[0-9].root");
    else {
      Error("runLocalMuonResolution","cannot proceed with these files there otherwise results will be mixed up!");
      return;
    }
  }
  
  // OCDB access
  if(gMode == kLocal) AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  else if (gMode != kProof) {
    if (!TGrid::Connect("alien://")) return;
    if(gMode == kInteractif_ESDList || !gSystem->AccessPathName("ConfigureCuts.C"))
      AliCDBManager::Instance()->SetDefaultStorage("raw://");
  }
  
  // Create input object
  TObject* inputObj = 0x0;
  if (gMode == kProof) inputObj = new TObjString(inputFileName);
  else inputObj = CreateChain(gMode, inputFileName);
  if (!inputObj) return;
  
  // set starting chamber resolution (if -1 they will be loaded from recoParam in the task)
  Double_t gClusterResNB[10];
  Double_t gClusterResB[10];
  Double_t gClusterResNBErr[10];
  Double_t gClusterResBErr[10];
  for (Int_t i = 0; i < 10; i++) {
    gClusterResNB[i] = -1.;
    gClusterResB[i] = -1.;
    gClusterResNBErr[i] = 0.;
    gClusterResBErr[i] = 0.;
  }
  
  // output graphs
  TMultiGraph* mgClusterResXVsStep = new TMultiGraph("mgClusterResXVsStep","cluster X-resolution versus step;step;#sigma_{X} (cm)");
  TMultiGraph* mgClusterResYVsStep = new TMultiGraph("mgClusterResYVsStep","cluster Y-resolution versus step;step;#sigma_{Y} (cm)");
  TGraphErrors* gClusterResXVsStep[10];
  TGraphErrors* gClusterResYVsStep[10];
  for (Int_t i = 0; i < 10; i++) {
    gClusterResXVsStep[i] = new TGraphErrors(nSteps+1);
    gClusterResXVsStep[i]->SetName(Form("gResX_ch%d",i+1));
    gClusterResXVsStep[i]->SetMarkerStyle(kFullDotMedium);
    gClusterResXVsStep[i]->SetMarkerColor(i+1+i/9);
    mgClusterResXVsStep->Add(gClusterResXVsStep[i],"lp");
    
    gClusterResYVsStep[i] = new TGraphErrors(nSteps+1);
    gClusterResYVsStep[i]->SetName(Form("gResY_ch%d",i+1));
    gClusterResYVsStep[i]->SetMarkerStyle(kFullDotMedium);
    gClusterResYVsStep[i]->SetMarkerColor(i+1+i/9);
    mgClusterResYVsStep->Add(gClusterResYVsStep[i],"lp");
  }
  
  // loop over step
  for (Int_t gStep = 0; gStep < nSteps; gStep++) {
    cout<<"step "<<gStep+1<<"/"<<nSteps<<endl;
    
    // Connect to proof if needed and prepare environment
    if (gMode == kProof) {
      
      // set general environment and close previous session
      if (gStep == 0) gEnv->SetValue("XSec.GSI.DelegProxy","2");
      else gProof->Close("s");
      
      // connect
      if (gSystem->Getenv("alien_API_USER") == NULL) TProof::Open("alice-caf.cern.ch","workers=40");
      else TProof::Open(Form("%s@alice-caf.cern.ch",gSystem->Getenv("alien_API_USER")),"workers=40");
      if (!gProof) return;
      
      // set environment and load libraries on workers
      TList* list = new TList();
      list->Add(new TNamed("ALIROOT_MODE", ""));
      list->Add(new TNamed("ALIROOT_EXTRA_LIBS", extraLibs.Data()));
      if (!taskInPWG3) list->Add(new TNamed("ALIROOT_EXTRA_INCLUDES", "MUON:MUON/mapping"));
      gProof->EnablePackage(aliroot.Data(), list);
      
      // compile task on workers
      if (!taskInPWG3) gProof->Load("AliAnalysisTaskMuonResolution.cxx++g");
      
      // prepare OCDB access on workers
      gProof->Exec("AliCDBManager::Instance()->SetDefaultStorage(\"raw://\")");
      
    }
    
    cout<<gMode<<endl;
    cout<<gStep<<endl;
    cout<<gNEvents<<endl;
    cout<<inputObj<<endl;
    cout<<gExtrapMode<<endl;
    cout<<gCorrectForSystematics<<endl;
    cout<<gMinMomentum<<endl;
    cout<<gMatchTrig<<endl;
    cout<<gSelectPhysics<<endl;
    
    // run one step
    run(gMode, gStep, gNEvents, inputObj, gExtrapMode, gCorrectForSystematics, gMinMomentum, gMatchTrig, gSelectPhysics, gClusterResNB, gClusterResB);
    
    // fill graph with starting resolutions from the task at first step
    if (gStep == 0) for (Int_t i = 0; i < 10; i++) {
      gClusterResXVsStep[i]->SetPoint(0, 0, gClusterResNB[i]);
      gClusterResXVsStep[i]->SetPointError(0, 0., gClusterResNBErr[i]);
      gClusterResYVsStep[i]->SetPoint(0, 0, gClusterResB[i]);
      gClusterResYVsStep[i]->SetPointError(0, 0., gClusterResBErr[i]);
    }
    
    // read the chamber resolution from the output file
    if (!GetChamberResolution(gStep, gClusterResNB, gClusterResB, gClusterResNBErr, gClusterResBErr)) return;
    
    // fill graphs with computed resolutions
    for (Int_t i = 0; i < 10; i++) {
      gClusterResXVsStep[i]->SetPoint(gStep+1, gStep+1, gClusterResNB[i]);
      gClusterResXVsStep[i]->SetPointError(gStep+1, 0., gClusterResNBErr[i]);
      gClusterResYVsStep[i]->SetPoint(gStep+1, gStep+1, gClusterResB[i]);
      gClusterResYVsStep[i]->SetPointError(gStep+1, 0., gClusterResBErr[i]);
    }
    
  }
  
  // copy final results in results.root file
  gSystem->Exec(Form("cp chamberResolution_step%d.root results.root", nSteps-1));
  
  // display convergence
  TCanvas* convergence = new TCanvas("convergence","convergence");
  convergence->Divide(1,2);
  convergence->cd(1);
  mgClusterResXVsStep->Draw("ap");
  convergence->cd(2);
  mgClusterResYVsStep->Draw("ap");
  
  // save convergence plots
  TFile* outFile = TFile::Open("results.root","UPDATE");
  if (!outFile || !outFile->IsOpen()) return;
  outFile->cd();
  mgClusterResXVsStep->Write();
  mgClusterResYVsStep->Write();
  convergence->Write();
  outFile->Close();
  
  // print results
  printf("\nchamber resolution:\n");
  printf(" - non-bending:");
  for (Int_t i = 0; i < 10; i++) printf((i==0)?" %5.3f":", %5.3f",gClusterResNB[i]);
  printf("\n -     bending:");
  for (Int_t i = 0; i < 10; i++) printf((i==0)?" %6.4f":", %6.4f",gClusterResB[i]);
  printf("\n\n");
  
  // ...timer stop
  localTimer->Stop();
  localTimer->Print();
  
}

//______________________________________________________________________________
void LoadAlirootLocaly(TString& extraLibs)
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
  
  // Use AliRoot includes and compile our task
  if (!taskInPWG3) {
    gROOT->ProcessLine(".include $ALICE_ROOT/include");
    gROOT->ProcessLine(".include $ALICE_ROOT/MUON");
    gROOT->ProcessLine(".include $ALICE_ROOT/MUON/mapping");
    gROOT->LoadMacro("AliAnalysisTaskMuonResolution.cxx++g");
  }
  
}

//______________________________________________________________________________
void run(Int_t mode, Int_t iStep, Int_t nevents, TObject* input, Int_t extrapMode, Bool_t correctForSystematics,
	 Double_t minMomentum, Bool_t matchTrig, Bool_t selectPhysics, Double_t clusterResNB[10], Double_t clusterResB[10])
{
  /// launch the analysis with these parameters
  
  // Create the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("MuonResolutionAnalysis");
  
  // ESD input handler
  AliESDInputHandler* esdH = new AliESDInputHandler();
  esdH->SetReadFriends(kFALSE);
  mgr->SetInputEventHandler(esdH);
  
  // event selection
  if (selectPhysics) {
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
    AliPhysicsSelectionTask* physicsSelection = AddTaskPhysicsSelection();
    if(!physicsSelection) {
      Error("RunLocalMuonResolution","AliPhysicsSelectionTask not created!");
      return;
    }
    physicsSelection->GetPhysicsSelection()->SetUseMuonTriggers();
  }
  
  // Muon Resolution analysis
  if (taskInPWG3) gROOT->LoadMacro("$ALICE_ROOT/PWG3/muondep/AddTaskMuonResolution.C");
  else gROOT->LoadMacro("AddTaskMuonResolution.C");
  AliAnalysisManager::SetCommonFileName(Form("chamberResolution_step%d.root", iStep));
  AliAnalysisTaskMuonResolution* muonResolution = AddTaskMuonResolution(selectPhysics, matchTrig, minMomentum, correctForSystematics, extrapMode);
  if(!muonResolution) {
    Error("RunLocalMuonResolution","AliAnalysisTaskMuonResolution not created!");
    return;
  }
//  muonResolution->ReAlign("", "alien://folder=/alice/cern.ch/user/p/ppillot/OCDB");
//  muonResolution->ReAlign(0x0, "");
  muonResolution->SetStartingResolution(clusterResNB, clusterResB);
  
  // Enable debug printouts
  //mgr->SetDebugLevel(2);
  
  // start local analysis
  if (mgr->InitAnalysis()) {
    mgr->PrintStatus();
    if (mode == kProof) mgr->StartAnalysis("proof", Form("%s",static_cast<TObjString*>(input)->GetName()), nevents);
    else mgr->StartAnalysis("local", static_cast<TChain*>(input), nevents);
  }
  
  // save the summary canvases
  if (muonResolution->GetCanvases()) {
    TFile* outFile = TFile::Open(Form("chamberResolution_step%d.root", iStep),"UPDATE");
    if (outFile && outFile->IsOpen()) {
      muonResolution->GetCanvases()->Write();
      outFile->Close();
    }
  }
  
  // return starting chamber resolution from the task
  muonResolution->GetStartingResolution(clusterResNB, clusterResB);
  
  // clean memory
  delete mgr;
  TObject::SetObjectStat(kFALSE);
}

//______________________________________________________________________________
Bool_t GetChamberResolution(Int_t iStep, Double_t clusterResNB[10], Double_t clusterResB[10], Double_t clusterResNBErr[10], Double_t clusterResBErr[10])
{
  /// read the chamber resolution from the output file
  
  TFile* outFile = TFile::Open(Form("chamberResolution_step%d.root", iStep),"READ");
  
  if (!outFile || !outFile->IsOpen()) {
    Error("GetChamberResolution","output file does not exist!");
    return kFALSE;
  }
  
  TObjArray* summary = static_cast<TObjArray*>(outFile->FindObjectAny("ChamberRes"));
  TGraphErrors* gCombinedResidualXPerChSigma = (summary) ? static_cast<TGraphErrors*>(summary->FindObject("gCombinedResidualXPerChSigma")) : 0x0;
  TGraphErrors* gCombinedResidualYPerChSigma = (summary) ? static_cast<TGraphErrors*>(summary->FindObject("gCombinedResidualYPerChSigma")) : 0x0;
  
  if (!gCombinedResidualXPerChSigma || !gCombinedResidualYPerChSigma) {
    Error("GetChamberResolution","resolution graphs do not exist!");
    return kFALSE;
  }
  
  Double_t dummy;
  for (Int_t i = 0; i < 10; i++) {
    gCombinedResidualXPerChSigma->GetPoint(i, dummy, clusterResNB[i]);
    gCombinedResidualYPerChSigma->GetPoint(i, dummy, clusterResB[i]);
    clusterResNBErr[i] = gCombinedResidualXPerChSigma->GetErrorY(i);
    clusterResBErr[i] = gCombinedResidualYPerChSigma->GetErrorY(i);
  }
  
  outFile->Close();
  
  return kTRUE;
}

//______________________________________________________________________________
Int_t GetMode(TString smode, TString input)
{
  if (smode == "local") {
    if ( input.EndsWith(".xml") ) return kInteractif_xml;
    else if ( input.EndsWith(".txt") ) return kInteractif_ESDList;
    else if ( input.EndsWith(".root") ) return kLocal;    
  } else if (smode == "proof") return kProof;
  return -1;
}

//______________________________________________________________________________
TChain* CreateChainFromCollection(const char *xmlfile)
{
  // Create a chain from the collection of tags.
  TGridCollection* coll = TAlienCollection::Open(xmlfile);
  if (!coll) {
    Error("CreateChainFromCollection", "Cannot create the AliEn collection");
    return NULL;
  }
  
  TGridResult* tagResult = coll->GetGridResult("",kFALSE,kFALSE);
  AliTagAnalysis *tagAna = new AliTagAnalysis("ESD");
  tagAna->ChainGridTags(tagResult);
  
  AliRunTagCuts      *runCuts = new AliRunTagCuts();
  AliLHCTagCuts      *lhcCuts = new AliLHCTagCuts();
  AliDetectorTagCuts *detCuts = new AliDetectorTagCuts();
  AliEventTagCuts    *evCuts  = new AliEventTagCuts();
  
  // Check if the cuts configuration file was provided
  if (!gSystem->AccessPathName("ConfigureCuts.C")) {
    gROOT->LoadMacro("ConfigureCuts.C");
    ConfigureCuts(runCuts, lhcCuts, detCuts, evCuts);
  }
  
  TChain *chain = tagAna->QueryTags(runCuts, lhcCuts, detCuts, evCuts);
  if (!chain || !chain->GetNtrees()) return NULL;
  chain->ls();
  return chain;
}

//______________________________________________________________________________
TChain* CreateChainFromFile(const char *rootfile)
{
  // Create a chain using the root file.
  TChain* chain = new TChain("esdTree");
  chain->Add(rootfile);
  if (!chain->GetNtrees()) return NULL;
  chain->ls();
  return chain;
}

//______________________________________________________________________________
TChain* CreateChainFromESDList(const char *esdList)
{
  // Create a chain using tags from the run list.
  TChain* chain = new TChain("esdTree");
  ifstream inFile(esdList);
  TString inFileName;
  if (inFile.is_open()) {
    while (! inFile.eof() ) {
      inFileName.ReadLine(inFile,kFALSE);
      if(!inFileName.EndsWith(".root")) continue;
      chain->Add(inFileName.Data());
    }
  }
  inFile.close();
  if (!chain->GetNtrees()) return NULL;
  chain->ls();
  return chain;
}

//______________________________________________________________________________
TChain* CreateChain(Int_t mode, TString input)
{
  printf("*******************************\n");
  printf("*** Getting the Chain       ***\n");
  printf("*******************************\n");
  if(mode == kInteractif_xml) return CreateChainFromCollection(input.Data());
  else if (mode == kInteractif_ESDList) return CreateChainFromESDList(input.Data());
  else if (mode == kLocal) return CreateChainFromFile(input.Data());
  else return NULL;
}

