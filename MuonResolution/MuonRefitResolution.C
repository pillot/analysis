/*
 *  MuonRefitResolution.C
 *  aliroot
 *
 *  Created by Philippe Pillot on 11/11/10.
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
#include "TAxis.h"
#include "THashList.h"
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
#include "AliPhysicsSelectionTask.h"
#include "AliPhysicsSelection.h"
#include "AliBackgroundSelection.h"
#include "AliAnalysisDataContainer.h"
#include "AliAnalysisTaskMuonResolution.h"
#include "AliAnalysisTaskMuonRefit.h"

// MUON includes
#include "AliMpCDB.h"
#include "AliMpDetElement.h"
#include "AliMpDDLStore.h"
#include "AliMUONCalibParamND.h"
#include "AliMUON2DMap.h"
#include "AliMUONTrackerData.h"
#include "AliMUONPainterDataRegistry.h"
#include "AliMUONTrackerDataWrapper.h"

#include "AddTaskPhysicsSelection.C"
#include "AddTaskMuonResolution.C"
#include "ConfigureCuts.C"

#endif

TString aliroot="VO_ALICE@AliRoot::v4-21-03-AN";

enum {kLocal, kInteractif_xml, kInteractif_ESDList, kProof};

Bool_t  run(Int_t mode, Int_t iStep, Int_t nevents, TObject* inputObj, Int_t extrapMode, Bool_t correctForSystematics,
	    Double_t minMomentum, Bool_t matchTrig, Bool_t selectPhysics, Double_t clusterResNB[10], Double_t clusterResB[10]);
Bool_t  GetChamberResolution(Int_t iStep, Double_t clusterResNB[10], Double_t clusterResB[10],
			     Double_t clusterResNBErr[10], Double_t clusterResBErr[10]);
void    AddMCHViews(TFile* file);
AliMUONTrackerData* ConvertGraph(TGraphErrors& g, const char* name);
Int_t   GetMode(TString smode, TString input);
TChain* CreateChainFromCollection(const char *xmlfile);
TChain* CreateChainFromFile(const char *rootfile);
TChain* CreateChainFromESDList(const char *esdList);
TChain* CreateChain(Int_t mode, TString input);

//______________________________________________________________________________
void MuonRefitResolution(TString smode, TString inputFileName, Int_t nSteps, Bool_t selectPhysics, Bool_t matchTrig, Double_t minMomentum,
			 Bool_t correctForSystematics, Int_t extrapMode, Int_t nevents, TString extraLibs)
{
  /// Compute the cluster resolution by studying cluster-track residual, deconvoluting from track resolution
  /// if extrapMode == 0: extrapolate from the closest cluster
  /// if extrapMode == 1: extrapolate from the previous cluster except between stations 2-3-4
  /// if correctForSystematics == kTRUE: the systematic shifts of the residuals is included in the resolution
  
  // timer start...
  TStopwatch* localTimer = new TStopwatch;
  
  // check parameters
  nSteps = TMath::Max(nSteps,1);
  if (extrapMode != 0 && extrapMode != 1) {
    Error("runLocalMuonResolution","incorrect extrapolation mode!");
    return;
  }
  
  // Check runing mode
  Int_t mode = GetMode(smode, inputFileName);
  if(mode < 0){
    Error("runLocalMuonResolution","Please provide either an ESD root file a collection of ESDs or a dataset.");
    return;
  }
  
  // Create input object
  TObject* inputObj = 0x0;
  if (mode == kProof) inputObj = new TObjString(inputFileName);
  else inputObj = CreateChain(mode, inputFileName);
  if (!inputObj) return;
  
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
  //  if(mode == kLocal) AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  if(mode == kLocal) AliCDBManager::Instance()->SetDefaultStorage("raw://");
  else if (mode != kProof) {
    if (!TGrid::Connect("alien://")) return;
    if(mode == kInteractif_ESDList || !gSystem->AccessPathName("ConfigureCuts.C"))
      AliCDBManager::Instance()->SetDefaultStorage("raw://");
  }
  
  // set starting chamber resolution (if -1 they will be loaded from recoParam in the task)
  Double_t clusterResNB[10] ={-1., -1., -1., -1., -1., -1., -1., -1., -1., -1.};
  Double_t clusterResB[10] ={-1., -1., -1., -1., -1., -1., -1., -1., -1., -1.};
  
  //  Double_t clusterResNB[10] ={0.06, 0.06, 0.06, 0.06, 0.06, 0.06, 0.06, 0.06, 0.06, 0.06};
  //  Double_t clusterResB[10] ={0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02};
  
  //Double_t clusterResNB[10] ={0.061, 0.059, 0.101, 0.103, 0.092, 0.108, 0.137, 0.132, 0.249, 0.145};
  //Double_t clusterResB[10] ={0.0392, 0.0335, 0.0718, 0.0701, 0.1127, 0.0881, 0.1112, 0.1160, 0.1338, 0.1515};
  
  // simu ideal realRecoParam (RMS)
  //Double_t clusterResNB[10] ={0.062, 0.061, 0.073, 0.073, 0.073, 0.072, 0.072, 0.071, 0.071, 0.073};
  //Double_t clusterResB[10] ={0.0070, 0.0064, 0.0104, 0.0111, 0.0139, 0.0138, 0.0125, 0.0116, 0.0109, 0.0133};
  
  // simu ideal realRecoParam (fit)
  //Double_t clusterResNB[10] ={0.067, 0.068, 0.086, 0.086, 0.086, 0.085, 0.085, 0.085, 0.084, 0.083};
  //Double_t clusterResB[10] ={0.0050, 0.0047, 0.0071, 0.0070, 0.0101, 0.0101, 0.0080, 0.0074, 0.0069, 0.0081};
  
  Double_t clusterResNBErr[10] ={0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
  Double_t clusterResBErr[10] ={0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
  
  // output graphs
  TMultiGraph* mgClusterResXVsStep = new TMultiGraph("mgClusterResXVsStep","cluster X-resolution versus step;step;#sigma_{X} (cm)");
  TMultiGraph* mgClusterResYVsStep = new TMultiGraph("mgClusterResYVsStep","cluster Y-resolution versus step;step;#sigma_{Y} (cm)");
  TGraphErrors* clusterResXVsStep[10];
  TGraphErrors* clusterResYVsStep[10];
  for (Int_t i = 0; i < 10; i++) {
    clusterResXVsStep[i] = new TGraphErrors(nSteps+1);
    clusterResXVsStep[i]->SetName(Form("gResX_ch%d",i+1));
    clusterResXVsStep[i]->SetMarkerStyle(kFullDotMedium);
    clusterResXVsStep[i]->SetMarkerColor(i+1+i/9);
    mgClusterResXVsStep->Add(clusterResXVsStep[i],"lp");
    
    clusterResYVsStep[i] = new TGraphErrors(nSteps+1);
    clusterResYVsStep[i]->SetName(Form("gResY_ch%d",i+1));
    clusterResYVsStep[i]->SetMarkerStyle(kFullDotMedium);
    clusterResYVsStep[i]->SetMarkerColor(i+1+i/9);
    mgClusterResYVsStep->Add(clusterResYVsStep[i],"lp");
  }
  
  // loop over step
  for (Int_t iStep = 0; iStep < nSteps; iStep++) {
    cout<<"step "<<iStep+1<<"/"<<nSteps<<endl;
    
    // Connect to proof if needed and prepare environment
    if (mode == kProof) {
      
      // set general environment and close previous session
      if (iStep == 0) gEnv->SetValue("XSec.GSI.DelegProxy","2");
      else gProof->Close("s");
      
      // connect
//      TString location = (smode == "caf") ? "alice-caf.cern.ch" : "nansafmaster.in2p3.fr";
      TString location = (smode == "caf") ? "alice-caf.cern.ch" : "localhost:1093";
      TString nWorkers = (smode == "caf") ? "workers=40" : "";
      if (gSystem->Getenv("alien_API_USER") == NULL) TProof::Open(location.Data(), nWorkers.Data());
      else TProof::Open(Form("%s@%s",gSystem->Getenv("alien_API_USER"), location.Data()), nWorkers.Data());
      if (!gProof) return;
      
      // set environment and load libraries on workers
      TList* list = new TList();
      list->Add(new TNamed("ALIROOT_MODE", ""));
      list->Add(new TNamed("ALIROOT_EXTRA_LIBS", extraLibs.Data()));
      if (!gSystem->AccessPathName("AliAnalysisTaskMuonResolution.cxx"))
	list->Add(new TNamed("ALIROOT_EXTRA_INCLUDES", "MUON:MUON/mapping"));
      gProof->EnablePackage(aliroot.Data(), list, kTRUE);
      
      // compile task on workers
      gProof->Load("AliAnalysisTaskMuonResolution.cxx++g", kTRUE);
      gProof->Load("AliAnalysisTaskMuonRefit.cxx++g", kTRUE);
      
      // prepare OCDB access on workers
      gProof->Exec("AliCDBManager::Instance()->SetDefaultStorage(\"raw://\")");
      
    }
    
    // run one step
    if (!run(mode, iStep, nevents, inputObj, extrapMode, correctForSystematics, minMomentum, matchTrig, selectPhysics, clusterResNB, clusterResB)) return;
    
    // fill graph with starting resolutions from the task at first step
    if (iStep == 0) for (Int_t i = 0; i < 10; i++) {
      clusterResXVsStep[i]->SetPoint(0, 0, clusterResNB[i]);
      clusterResXVsStep[i]->SetPointError(0, 0., clusterResNBErr[i]);
      clusterResYVsStep[i]->SetPoint(0, 0, clusterResB[i]);
      clusterResYVsStep[i]->SetPointError(0, 0., clusterResBErr[i]);
    }
    
    // read the chamber resolution from the output file
    if (!GetChamberResolution(iStep, clusterResNB, clusterResB, clusterResNBErr, clusterResBErr)) return;
    
    // fill graphs with computed resolutions
    for (Int_t i = 0; i < 10; i++) {
      clusterResXVsStep[i]->SetPoint(iStep+1, iStep+1, clusterResNB[i]);
      clusterResXVsStep[i]->SetPointError(iStep+1, 0., clusterResNBErr[i]);
      clusterResYVsStep[i]->SetPoint(iStep+1, iStep+1, clusterResB[i]);
      clusterResYVsStep[i]->SetPointError(iStep+1, 0., clusterResBErr[i]);
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
  
  // ...timer stop
  localTimer->Stop();
  localTimer->Print();
  
}

//______________________________________________________________________________
Bool_t run(Int_t mode, Int_t iStep, Int_t nevents, TObject* inputObj, Int_t extrapMode, Bool_t correctForSystematics,
	   Double_t minMomentum, Bool_t matchTrig, Bool_t selectPhysics, Double_t clusterResNB[10], Double_t clusterResB[10])
{
  /// launch the analysis with these parameters
  
  // Create the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("MuonResolutionAnalysis");
  
  // ESD input handler
  AliESDInputHandler* esdH = new AliESDInputHandler();
  esdH->SetReadFriends(kFALSE);
  esdH->SetInactiveBranches(" FMD PHOS  EMCAL  Pmd Trd V0s TPC "
			    "Cascades Kinks CaloClusters ACORDE RawData HLT TZERO ZDC"
			    " Cells ACORDE Pileup");
  mgr->SetInputEventHandler(esdH);
  
  // event selection
  if (selectPhysics) {
    AliPhysicsSelectionTask* physicsSelection = AddTaskPhysicsSelection();
    if(!physicsSelection) {
      Error("MuonResolution","AliPhysicsSelectionTask not created!");
      return kFALSE;
    }
  }
  
  // second task (refit)
  AliAnalysisTaskMuonRefit *refit = new AliAnalysisTaskMuonRefit("MuonRefit");
  refit->SetDefaultStorage("raw://");
  //refit->ReAlign("alien://folder=/alice/cern.ch/user/p/ppillot/OCDB", "");
  Double_t clResNB[10];
  Double_t clResB[10];
  for (Int_t i=0; i<10; i++) { clResNB[i] = 0.2; clResB[i] = 0.2; }
  refit->ResetClusterResolution(clResNB, clResB);
  refit->ImproveTracks(kTRUE, 4.);
  mgr->AddTask(refit);
  mgr->ConnectInput(refit, 0, mgr->GetCommonInputContainer());
  
  // Muon Resolution analysis
  TString outputFileName = Form("chamberResolution_step%d.root", iStep);
  AliAnalysisManager::SetCommonFileName(outputFileName.Data());
  AliAnalysisTaskMuonResolution* muonResolution = AddTaskMuonResolution(selectPhysics, matchTrig, minMomentum, correctForSystematics, extrapMode);
  if(!muonResolution) {
    Error("MuonResolution","AliAnalysisTaskMuonResolution not created!");
    return kFALSE;
  }
  if (mode != kProof) muonResolution->ShowProgressBar();
  muonResolution->PrintClusterRes(kTRUE, kTRUE);
  muonResolution->ApplyAccCut();
  //muonResolution->SelectTrigger();
  muonResolution->FitResiduals();
  //muonResolution->ReAlign("alien://folder=/alice/cern.ch/user/p/ppillot/OCDB", "alien://folder=/alice/cern.ch/user/p/ppillot/OCDB_GMS");
  //muonResolution->ReAlign("alien://folder=/alice/cern.ch/user/p/ppillot/OCDB", "");
  //  muonResolution->ReAlign(0x0, "");
  muonResolution->SetStartingResolution(clusterResNB, clusterResB);
  
  // Enable debug printouts
  //mgr->SetDebugLevel(2);
  
  // start local analysis
  if (mgr->InitAnalysis()) {
    mgr->PrintStatus();
    if (mode == kProof) mgr->StartAnalysis("proof", Form("%s",static_cast<TObjString*>(inputObj)->GetName()), nevents);
    else mgr->StartAnalysis("local", static_cast<TChain*>(inputObj), nevents);
  }
  
  // save the summary canvases and mchview display
  if (muonResolution->GetCanvases()) {
    TFile* outFile = TFile::Open(outputFileName.Data(),"UPDATE");
    if (outFile && outFile->IsOpen()) {
      outFile->cd();
      muonResolution->GetCanvases()->Write();
      AddMCHViews(outFile);
      outFile->Close();
    }
  }
  
  // return starting chamber resolution from the task
  muonResolution->GetStartingResolution(clusterResNB, clusterResB);
  
  // clean memory
  delete mgr;
  TObject::SetObjectStat(kFALSE);
  
  return kTRUE;
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
void AddMCHViews(TFile* file)
{
  /// Get from the file the graphs containing data per DE, convert them into mchview objects and save them
  
  if (  ! AliMpDDLStore::Instance(false) )
  {
    Warning("AddMCHViews","mapping was not loaded. Loading it from $ALICE_ROOT/OCDB");
    AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
    AliCDBManager::Instance()->SetRun(999999999);
  }
  
  AliMpCDB::LoadAll();
  
  TObjArray* summary = static_cast<TObjArray*>(file->FindObjectAny("ChamberRes"));
  if (!summary) {
    Error("AddMCHViews","resolution graphs do not exist!");
    return;
  }
  
  TGraphErrors* g = 0x0;
  g = static_cast<TGraphErrors*>(summary->FindObject("gCombinedResidualXPerDESigma"));
  if (g) {
    file->cd();
    ConvertGraph(*g, "resoX")->Write();
  }
  
  g = static_cast<TGraphErrors*>(summary->FindObject("gCombinedResidualYPerDESigma"));
  if (g) {
    file->cd();
    ConvertGraph(*g, "resoY")->Write();
  }
  
  g = static_cast<TGraphErrors*>(summary->FindObject("gResidualXPerDEMean_ClusterOut"));
  if (g) {
    file->cd();
    ConvertGraph(*g, "shiftX")->Write();
  }
  
  g = static_cast<TGraphErrors*>(summary->FindObject("gResidualYPerDEMean_ClusterOut"));
  if (g) {
    file->cd();
    ConvertGraph(*g, "shiftY")->Write();
  }
}

//______________________________________________________________________________
AliMUONTrackerData* ConvertGraph(TGraphErrors& g, const char* name)
{
  /// Convert graph containing data per DE into mchview object
  
  AliMUON2DMap deValues(kFALSE);
  
  for ( Int_t i = 0 ; i < g.GetN(); ++i ) 
  {
    double y = g.GetY()[i];
    double ey = g.GetEY()[i];
    int detElemId;
    sscanf(g.GetXaxis()->GetBinLabel(i+1),"%d",&detElemId);
    
    AliMpDetElement* de = AliMpDDLStore::Instance()->GetDetElement(detElemId);
    
    AliMUONVCalibParam* param = new AliMUONCalibParamND(5, 1, detElemId, 0);
    
    Double_t sumn = 1000.0;
    Double_t sumw = sumn*y;
    Double_t sumw2 = (sumn-1)*ey*ey+sumw*sumw/sumn;
    
    param->SetValueAsDouble(0,0,sumw);
    param->SetValueAsDouble(0,1,sumw2);
    param->SetValueAsDouble(0,2,sumn);
    param->SetValueAsDouble(0,3,de->NofChannels());
    param->SetValueAsDouble(0,4,1);
    
    deValues.Add(param);
  }
  
  AliMUONTrackerData* data = new AliMUONTrackerData(name,name,deValues,1);
  data->SetDimensionName(0,name);
  
  return data;
}

//______________________________________________________________________________
Int_t GetMode(TString smode, TString input)
{
  if (smode == "local") {
    if ( input.EndsWith(".xml") ) return kInteractif_xml;
    else if ( input.EndsWith(".txt") ) return kInteractif_ESDList;
    else if ( input.EndsWith(".root") ) return kLocal;    
  } else if (smode == "caf" || smode == "saf") return kProof;
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
  ConfigureCuts(runCuts, lhcCuts, detCuts, evCuts);
  
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

