/*
 *  runLocalESDCheck.C
 *  aliroot
 *
 *  Created by Philippe Pillot on 05/04/10.
 *  Copyright 2010 SUBATECH. All rights reserved.
 *
 */

enum {kLocal, kInteractif_xml, kInteractif_ESDList, kProof};

// needed in Proof mode, certainly because of the interpretor
TStopwatch* localTimer;
TMultiGraph* mgClusterResXVsStep;
TMultiGraph* mgClusterResYVsStep;
Double_t clusterResNB[10];
Double_t clusterResB[10];
//Double_t clusterResNB[10] = {0.075, 0.060, 0.124, 0.120, 0.095, 0.103, 0.125, 0.140, 0.169, 0.180};
//Double_t clusterResB[10] = {0.0573, 0.0430, 0.0938, 0.0919, 0.1451, 0.1133, 0.1773, 0.1899, 0.1156, 0.1487};
Double_t clusterResNBErr[10];
Double_t clusterResBErr[10];

void runLocalMuonChamberResolution(TString smode = "local", TString input = "AliESDs.root", Int_t nSteps = 10,
				   Int_t extrapMode = 1, Bool_t correctForSystematics = kTRUE, Double_t minMomentum = 0.,
				   Bool_t matchTrig = kFALSE, Int_t nevents = 1234567890)
{
  /// Compute the cluster resolution by studying cluster-track residual, deconvoluting from track resolution
  /// if extrapMode == 0: extrapolate from the closest cluster
  /// if extrapMode == 1: extrapolate from the previous cluster except between stations 2-3-4
  /// if correctForSystematics == kTRUE: the systematic shifts of the residuals is included in the resolution
  
  // timer start...
  localTimer = new TStopwatch;
  
  // check parameters
  nSteps = TMath::Max(nSteps,1);
  if (extrapMode != 0 && extrapMode != 1) {
    printf("Fatal: wrong extrapolation mode... Exiting!\n");
    return;
  }
  
  // Check runing mode
  Int_t mode = GetMode(smode, input);
  if(mode < 0){
    printf("Fatal: Please provide either an ESD root file or a collection of ESDs... Exiting!\n");
    return;
  }
  
  // Load common libraries
  gSystem->Load("libVMC");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  
  // Load additional libraries
  gSystem->Load("libProofPlayer");
  TString extraLibs="Physics:Minuit:XMLParser:Gui:RAWDatabase:CDB:STEER:MUONcore:MUONmapping:MUONcalib:MUONgeometry:MUONtrigger:MUONraw:MUONbase:MUONrec";
  TObjArray* libs = extraLibs.Tokenize(":");
  for (Int_t i = 0; i < libs->GetEntriesFast(); i++)
    gSystem->Load(Form("lib%s",static_cast<TObjString*>(libs->UncheckedAt(i))->GetName()));
    
  // Use AliRoot includes to compile our task
  gROOT->ProcessLine(".include $ALICE_ROOT/include");
  gROOT->ProcessLine(".include $ALICE_ROOT/MUON");
  gROOT->ProcessLine(".include $ALICE_ROOT/MUON/mapping");
  
  // copy files needed for this analysis
  TString path ("/Users/philippe/Work/Alice/Work/Data/Macro/Resolution");
  TObjArray fileList(100);
  fileList.SetOwner();
  fileList.AddLast(new TObjString("runLocalMuonChamberResolution.C"));
  fileList.AddLast(new TObjString("AliAnalysisTaskMuonChamberResolution.cxx"));
  fileList.AddLast(new TObjString("AliAnalysisTaskMuonChamberResolution.h"));
  fileList.AddLast(new TObjString("ConfigureCuts.C"));
  TIter nextFile(&fileList);
  TObjString *file;
  char overwrite = '';
  while ((file = static_cast<TObjString*>(nextFile()))) {
    if (overwrite != 'a') {
      overwrite = '';
      if (!gSystem->AccessPathName(file->GetName())) {
	cout<<Form("file %s exist in current directory. Overwrite? [y=yes, n=no, a=all, k=keep all] ",file->GetName())<<flush;
	while (overwrite != 'y' && overwrite != 'n' && overwrite != 'a' && overwrite != 'k') cin>>overwrite;
      } else overwrite = 'y';
    }
    if (overwrite == 'y' || overwrite == 'a') gSystem->Exec(Form("cp %s/%s %s", path.Data(), file->GetName(),file->GetName()));
    else if (overwrite == 'k') break;
  }
  
  // check for old output file to removed
  char remove = '';
  if (!gSystem->Exec("ls chamberResolution_step*[0-9].root")) {
    cout<<"above files must be removed from the current directory. Delete? [y=yes, n=no] "<<flush;
    while (remove != 'y' && remove != 'n') cin>>remove;
    if (remove == 'y') gSystem->Exec("rm -f chamberResolution_step*[0-9].root");
    else {
      cout<<"cannot proceed with these files there otherwise results will be mixed up... Exiting!"<<endl;
      return;
    }
  }
  
  // compile task
  if (mode != kProof) gROOT->LoadMacro("AliAnalysisTaskMuonChamberResolution.cxx++g");
  
  // OCDB access
  if(mode == kLocal) AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  else if (mode != kProof) {
    if (!TGrid::Connect("alien://")) return;
    if(mode == kInteractif_ESDList || !gSystem->AccessPathName("ConfigureCuts.C"))
      AliCDBManager::Instance()->SetDefaultStorage("raw://");
  }
  
  // Create input object
  TObject* inputObj = 0x0;
  if (mode == kProof) inputObj = new TObjString(input);
  else inputObj = CreateChain(mode, input);
  if (!inputObj) return;
  
  // set starting chamber resolution (if -1 they will be loaded from recoParam in the task)
  for (Int_t i = 0; i < AliMUONConstants::NTrackingCh(); i++) {
    clusterResNB[i] = -1.;
    clusterResB[i] = -1.;
    clusterResNBErr[i] = 0.;
    clusterResBErr[i] = 0.;
  }
  
  // output graphs
  mgClusterResXVsStep = new TMultiGraph("mgClusterResXVsStep","cluster X-resolution versus step;step;#sigma_{X} (cm)");
  mgClusterResYVsStep = new TMultiGraph("mgClusterResYVsStep","cluster Y-resolution versus step;step;#sigma_{Y} (cm)");
  TGraphErrors* gClusterResXVsStep[10];
  TGraphErrors* gClusterResYVsStep[10];
  for (Int_t i = 0; i < AliMUONConstants::NTrackingCh(); i++) {
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
  for (Int_t iStep = 0; iStep < nSteps; iStep++) {
    cout<<"step "<<iStep+1<<"/"<<nSteps<<endl;
    
    // Connect to proof if needed and prepare environment
    if (mode == kProof) {
      
      // set general environment and close previous session
      if (iStep == 0) {
	gEnv->SetValue("XSec.GSI.DelegProxy","2");
	TProof::AddEnvVar("ALIROOT_EXTRA_LIBS", extraLibs.Data());
      } else gProof->Close("s");
      
      // connect
      TProof::Open("ppillot@alice-caf.cern.ch","workers=40");
//      TProof::Open("ppillot@alice-caf.cern.ch");
      if (!gProof) return;
      
      // set environment and compile task on workers
      gProof->EnablePackage("VO_ALICE@AliRoot::v4-19-15-AN");
      gProof->AddIncludePath("$ALICE_ROOT/MUON");
      gProof->AddIncludePath("$ALICE_ROOT/MUON/mapping");
      
      // compile task on workers
      gProof->Load("AliAnalysisTaskMuonChamberResolution.cxx++g");
      
      // prepare OCDB access on workers
      gProof->Exec("AliCDBManager::Instance()->SetDefaultStorage(\"raw://\")");
      
    }
    
    // run one step
    run(mode, iStep, nevents, inputObj, extrapMode, correctForSystematics, minMomentum, matchTrig, clusterResNB, clusterResB);
    
    // fill graph with starting resolutions from the task at first step
    if (iStep == 0) for (Int_t i = 0; i < AliMUONConstants::NTrackingCh(); i++) {
      gClusterResXVsStep[i]->SetPoint(0, 0, clusterResNB[i]);
      gClusterResXVsStep[i]->SetPointError(0, 0., clusterResNBErr[i]);
      gClusterResYVsStep[i]->SetPoint(0, 0, clusterResB[i]);
      gClusterResYVsStep[i]->SetPointError(0, 0., clusterResBErr[i]);
    }
    
    // read the chamber resolution from the output file
    if (!GetChamberResolution(iStep, clusterResNB, clusterResB, clusterResNBErr, clusterResBErr)) return;
    
    // fill graphs with computed resolutions
    for (Int_t i = 0; i < AliMUONConstants::NTrackingCh(); i++) {
      gClusterResXVsStep[i]->SetPoint(iStep+1, iStep+1, clusterResNB[i]);
      gClusterResXVsStep[i]->SetPointError(iStep+1, 0., clusterResNBErr[i]);
      gClusterResYVsStep[i]->SetPoint(iStep+1, iStep+1, clusterResB[i]);
      gClusterResYVsStep[i]->SetPointError(iStep+1, 0., clusterResBErr[i]);
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
  for (Int_t i = 0; i < AliMUONConstants::NTrackingCh(); i++) printf((i==0)?" %5.3f":", %5.3f",clusterResNB[i]);
  printf("\n -     bending:");
  for (Int_t i = 0; i < AliMUONConstants::NTrackingCh(); i++) printf((i==0)?" %6.4f":", %6.4f",clusterResB[i]);
  printf("\n\n");
  
  // ...timer stop
  localTimer->Stop();
  localTimer->Print();
  
}

//______________________________________________________________________________
void run(Int_t mode, Int_t iStep, Int_t nevents, TObject* input, Int_t extrapMode, Bool_t correctForSystematics,
	 Double_t minMomentum, Bool_t matchTrig, Double_t clusterResNB[10], Double_t clusterResB[10])
{
  /// launch the analysis with these parameters
  
  // Create the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("MuonChamberResolutionAnalysis");
  
  // ESD input handler
  AliESDInputHandler* esdH = new AliESDInputHandler();
  esdH->SetReadFriends(kFALSE);
  mgr->SetInputEventHandler(esdH);
  
  // ESD scan analysis
  AliAnalysisTaskMuonChamberResolution* task = new AliAnalysisTaskMuonChamberResolution("MuonChamberResolution");
  task->SetStartingResolution(clusterResNB, clusterResB);
  task->SetMinMomentum(minMomentum);
  task->MatchTrigger(matchTrig);
  task->SetExtrapMode(extrapMode);
  task->CorrectForSystematics(correctForSystematics);
  task->ReAlign("", "alien://folder=/alice/cern.ch/user/p/ppillot/OCDB");
//  task->ReAlign();
//  task->ReAlign(NULL, "");
  mgr->AddTask(task);
  
  // Create containers for output
  TString outFileName = Form("chamberResolution_step%d.root", iStep);
  AliAnalysisDataContainer *cout_histo1 = mgr->CreateContainer("Residuals", TObjArray::Class(), AliAnalysisManager::kOutputContainer, outFileName.Data());
  AliAnalysisDataContainer *cout_histo2 = mgr->CreateContainer("ResidualsVsP", TObjArray::Class(), AliAnalysisManager::kOutputContainer, outFileName.Data());
  AliAnalysisDataContainer *cout_histo3 = mgr->CreateContainer("Summary", TObjArray::Class(), AliAnalysisManager::kParamContainer, outFileName.Data());
  
  // Connect input/output
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, cout_histo1);
  mgr->ConnectOutput(task, 2, cout_histo2);
  mgr->ConnectOutput(task, 3, cout_histo3);
  
  // Enable debug printouts
  //mgr->SetDebugLevel(2);
  
  // start local analysis
  if (mgr->InitAnalysis()) {
    mgr->PrintStatus();
    if (mode == kProof) mgr->StartAnalysis("proof", Form("%s",static_cast<TObjString*>(input)->GetName()), nevents);
    else mgr->StartAnalysis("local", static_cast<TChain*>(input), nevents);
  }
  
  // save the summary canvases
  if (task->GetCanvases()) {
    TFile* outFile = TFile::Open(Form("chamberResolution_step%d.root", iStep),"UPDATE");
    if (outFile && outFile->IsOpen()) {
      task->GetCanvases()->Write();
      outFile->Close();
    }
  }
  
  // return starting chamber resolution from the task
  task->GetStartingResolution(clusterResNB, clusterResB);
  
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
    printf("error: output file does not exist\n");
    return kFALSE;
  }
  
  TObjArray* summary = static_cast<TObjArray*>(outFile->FindObjectAny("Summary"));
  TGraphErrors* gCombinedResidualXPerChSigma = (summary) ? static_cast<TGraphErrors*>(summary->FindObject("gCombinedResidualXPerChSigma")) : 0x0;
  TGraphErrors* gCombinedResidualYPerChSigma = (summary) ? static_cast<TGraphErrors*>(summary->FindObject("gCombinedResidualYPerChSigma")) : 0x0;
  
  if (!gCombinedResidualXPerChSigma || !gCombinedResidualYPerChSigma) {
    printf("error: resolution graphs do not exist\n");
    return kFALSE;
  }
  
  Double_t dummy;
  for (Int_t i = 0; i < AliMUONConstants::NTrackingCh(); i++) {
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
  TAlienCollection* coll = TAlienCollection::Open(xmlfile);
  if (!coll) {
    ::Error("CreateChainFromTags", "Cannot create an AliEn collection from %s", xmlfile);
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

