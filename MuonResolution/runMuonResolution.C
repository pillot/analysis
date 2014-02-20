/*
 *  runMuonResolution.C
 *  aliroot
 *
 *  Created by Philippe Pillot on 27/06/10.
 *  Copyright 2010 SUBATECH. All rights reserved.
 *
 */

TString rootVersion = "v5-28-00e";
TString alirootVersion = "v4-21-31-AN";
TString dataDir = "/alice/sim/LHC11a10a_bis";
TString dataPattern = "*AliESDs.root";
TString runFormat = "%d";
TString outDir = "Sim/LHC11a10a_bis/Fakes/selected";
Int_t ttl = 60000;
Int_t maxFilesPerJob = 100;
Int_t maxMergeFiles = 10;
Int_t maxMergeStages = 2;

//______________________________________________________________________________
void runMuonResolution(TString smode = "local", TString inputFileName = "AliESDs.root",
		       Bool_t selectPhysics = kTRUE, Bool_t selectTrigger = kFALSE, Bool_t matchTrig = kTRUE,
		       Bool_t applyAccCut = kTRUE, Double_t minMomentum = 0., Bool_t correctForSystematics = kTRUE,
		       Int_t extrapMode = 1, Int_t nevents = 1234567890)
{
  
  gStyle->SetOptFit(1);
  gROOT->LoadMacro("/Users/pillot/Work/Alice/Work/Macros/Facilities/runTaskFacilities.C");
  
  // --- Check runing mode ---
  Int_t mode = GetMode(smode, inputFileName);
  if(mode < 0) {
    Error("runMuonFakes","Please provide either an ESD root file a collection of ESDs or a dataset.");
    return;
  }
  
  // --- copy files needed for this analysis ---
  TList pathList; pathList.SetOwner();
  pathList.Add(new TObjString("/Users/pillot/Work/Alice/Work/Macros/MuonResolution"));
  pathList.Add(new TObjString("/Users/pillot/Work/Alice/Work/aliroot/PWG3/muondep"));
  TList fileList; fileList.SetOwner();
  fileList.Add(new TObjString("runMuonResolution.C"));
  fileList.Add(new TObjString("AddTaskMuonResolution.C"));
  fileList.Add(new TObjString("AliAnalysisTaskMuonResolution.cxx"));
  fileList.Add(new TObjString("AliAnalysisTaskMuonResolution.h"));
  CopyFileLocally(pathList, fileList);
  
  // --- prepare environment ---
  TString extraLibs = "RAWDatabase:CDB:STEER:MUONcore:MUONmapping:MUONcalib:MUONgeometry:MUONtrigger:MUONraw:MUONbase:MUONrec";
  TString extraIncs="include:MUON:MUON/mapping";
  TString extraTasks="AliAnalysisTaskMuonResolution";
  LoadAlirootLocally(extraLibs, extraIncs, extraTasks);
  AliAnalysisGrid *alienHandler = 0x0;
  if (mode == kProof) LoadAlirootOnProof(smode, alirootVersion, extraLibs, extraIncs, extraTasks, kTRUE);
  else if (mode == kGrid || mode == kTerminate) {
    TString analysisMacroName = "MuonResolutionAnalysis";
    alienHandler = static_cast<AliAnalysisGrid*>(CreateAlienHandler(smode, rootVersion, alirootVersion, inputFileName, dataDir, dataPattern, outDir, extraLibs, extraIncs, extraTasks, analysisMacroName, runFormat, ttl, maxFilesPerJob, maxMergeFiles, maxMergeStages));
    if (!alienHandler) return;
  }
  
  // --- Create the analysis train ---
  CreateAnalysisTrain(selectPhysics, selectTrigger, matchTrig, applyAccCut, minMomentum, correctForSystematics, extrapMode, alienHandler);
  
  // --- Create input object ---
  TObject* inputObj = CreateInputObject(mode, inputFileName);
  
  // --- start analysis ---
  StartAnalysis(mode, inputObj);
  
  // --- save summary canvases ---
//  if (mode != kTerminate) {
    AliAnalysisTaskMuonResolution *muonResolution = static_cast<AliAnalysisTaskMuonResolution*>(AliAnalysisManager::GetAnalysisManager()->GetTasks()->FindObject("MuonResolution"));
    if (muonResolution && muonResolution->GetCanvases()) {
      TFile* outFile = TFile::Open(AliAnalysisManager::GetCommonFileName(),"UPDATE");
      if (outFile && outFile->IsOpen()) {
	muonResolution->GetCanvases()->Write();
	AddMCHViews(outFile);
	outFile->Close();
	delete outFile;
      }
    }
//  }
  
}

//______________________________________________________________________________
void CreateAnalysisTrain(Bool_t selectPhysics, Bool_t selectTrigger, Bool_t matchTrig, Bool_t applyAccCut,
			 Double_t minMomentum, Bool_t correctForSystematics, Int_t extrapMode, TObject* alienHandler)
{
  /// create the analysis train and configure it
  
  // analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("MuonResolutionAnalysis");
  
  // Connect plugin to the analysis manager if any
  if (alienHandler) mgr->SetGridHandler(static_cast<AliAnalysisGrid*>(alienHandler));
  
  // ESD input
  AliESDInputHandler* esdH = new AliESDInputHandler();
  esdH->SetReadFriends(kFALSE);
  esdH->SetInactiveBranches("*");
  esdH->SetActiveBranches("MuonTracks AliESDRun. AliESDHeader. AliMultiplicity. AliESDFMD. AliESDVZERO. SPDVertex. PrimaryVertex. AliESDZDC.");
  mgr->SetInputEventHandler(esdH);
  
  // event selection
  if (selectPhysics) {
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
    AliPhysicsSelectionTask* physicsSelection = AddTaskPhysicsSelection();
    if (!physicsSelection) {
      Error("CreateAnalysisTrain","AliPhysicsSelectionTask not created!");
      return 0x0;
    }
  }
  
  // centrality selection
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskCentrality.C");
  AliCentralitySelectionTask* centralityTask = AddTaskCentrality();
  if (!centralityTask) {
    Error("CreateAnalysisTrain","AliCentralitySelectionTask not created!");
    return 0x0;
  }
  
  // Muon Resolution analysis
  gROOT->LoadMacro("AddTaskMuonResolution.C");
  AliAnalysisTaskMuonResolution *muonResolution = AddTaskMuonResolution(selectPhysics, selectTrigger, matchTrig, applyAccCut, minMomentum, correctForSystematics, extrapMode);
  if (!muonResolution) {
    Error("CreateAnalysisTrain","AliAnalysisTaskMuonResolution not created!");
    return 0x0;
  }
  /*if (mode == kLocal) muonResolution->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
   else */muonResolution->SetDefaultStorage("alien://folder=/alice/data/2011/OCDB");
  muonResolution->PrintClusterRes(kTRUE, kTRUE);
  Double_t clusterResNB[10] = {0.050, 0.049, 0.075, 0.074, 0.079, 0.081, 0.075, 0.070, 0.065, 0.068};
  Double_t clusterResB[10] = {0.0156, 0.0091, 0.0566, 0.0554, 0.0265, 0.0245, 0.0300, 0.0234, 0.0477, 0.0382};
  muonResolution->SetStartingResolution(clusterResNB, clusterResB);
  //muonResolution->RemoveMonoCathodClusters(kTRUE, kFALSE);
  //  muonResolution->FitResiduals(kFALSE);
  //  muonResolution->ReAlign("alien://folder=/alice/cern.ch/user/p/ppillot/OCDB2011","");
  muonResolution->ReAlign("", "alien://folder=/alice/cern.ch/user/p/ppillot/OCDB2011_Align2");
  //  muonResolution->ReAlign("", "alien://folder=/alice/cern.ch/user/j/jcastill/MATFtestCDBVanik");
  //  muonResolution->ReAlign("", "alien://folder=/alice/cern.ch/user/j/jcastill/MATFtestCDBJavier");
  //  muonResolution->ReAlign("alien://folder=/alice/cern.ch/user/p/ppillot/OCDB2011", "alien://folder=/alice/cern.ch/user/p/ppillot/OCDB2011_Align2bis");
  
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
    AliMUONTrackerData* data = static_cast<AliMUONTrackerData*>(ConvertGraph(*g, "resoX"));
    data->Write();
    delete data;
  }
  
  g = static_cast<TGraphErrors*>(summary->FindObject("gCombinedResidualYPerDESigma"));
  if (g) {
    file->cd();
    AliMUONTrackerData* data = static_cast<AliMUONTrackerData*>(ConvertGraph(*g, "resoY"));
    data->Write();
    delete data;
  }
  
  g = static_cast<TGraphErrors*>(summary->FindObject("gResidualXPerDEMean_ClusterOut"));
  if (g) {
    file->cd();
    AliMUONTrackerData* data = static_cast<AliMUONTrackerData*>(ConvertGraph(*g, "shiftX"));
    data->Write();
    delete data;
  }
  
  g = static_cast<TGraphErrors*>(summary->FindObject("gResidualYPerDEMean_ClusterOut"));
  if (g) {
    file->cd();
    AliMUONTrackerData* data = static_cast<AliMUONTrackerData*>(ConvertGraph(*g, "shiftY"));
    data->Write();
    delete data;
  }
}

//______________________________________________________________________________
TObject* ConvertGraph(TGraphErrors& g, const char* name)
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

