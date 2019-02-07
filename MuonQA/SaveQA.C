//
//  SetPath.C
//  
//
//  Created by Javier Martin-Blanco on 29/03/13.
//
//

#include <stdio.h>
#include <TArrayI.h>
#include <TString.h>
#include <TObjString.h>
#include <TObjArray.h>
#include <Riostream.h>
#include <fstream>
#include "TSystem.h"
#include "TString.h"
#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TGrid.h"
#include "AliAnalysisMuonUtility.h"

TString MultDirData = "HighMultiplicity";
TString MultDirMC = "HighMultiplicity";

//----------------------------------------------------------------------------
void SaveQA(TString inputDataFile, TString inputMCFile)
{
  
  ifstream inDataFile(inputDataFile.Data()); // open input file
  if (!inDataFile.is_open())
  {
    printf("cannot open file %s\n",inputDataFile.Data());
    return;
  }
  
  ifstream inMCFile(inputMCFile.Data()); // open input file
  if (!inMCFile.is_open())
  {
    printf("cannot open file %s\n",inputMCFile.Data());
    return;
  }
  
  if (TString(gSystem->GetFromPipe(TString::Format("grep -c alien:// %s", inputDataFile.Data()))).Atoi() > 0 ||
      TString(gSystem->GetFromPipe(TString::Format("grep -c alien:// %s", inputMCFile.Data()))).Atoi() > 0) {
    TGrid::Connect("alien://");
  }

  if (gSystem->AccessPathName("displays")) gSystem->Exec("mkdir displays");
  
  TString currDataName; // Current file path in Data
  TString currMCName; // Current file path in Data
  TString failedRuns;
  
  TCanvas cTmp;
  gStyle->SetOptStat(0);
  cTmp.Divide(2,1,0,0);
  
  while (!inDataFile.eof()) //loop over elements in input file
  {
    currDataName.ReadLine(inDataFile,kTRUE); // Read line in input Data file
    if(currDataName.IsNull()) continue;
    
    currMCName.ReadLine(inMCFile,kTRUE); // Read line in input MC file
    if(currMCName.IsNull()) continue;
    
    Int_t runNumber(AliAnalysisMuonUtility::GetRunNumber(currDataName.Data()));
    Int_t runNumberMC(AliAnalysisMuonUtility::GetRunNumber(currMCName.Data()));
    if (runNumberMC != runNumber) {
      cout << "ERROR: mismatch found between the 2 files... Aborting." << endl;
      return;
    }
    
    if (!gSystem->AccessPathName(Form("displays/%d/ESDclusterMapChamber1.png",runNumberMC))) {
      cout<<"run "<<runNumberMC<<" already processed. Remove files to reprocess..." << endl;
      continue;
    }
    
    cout<<"\rprocessing run "<<runNumberMC<<" ...\r"<<flush;
    
    TFile* currDataFile = TFile::Open(currDataName.Data(),"read");
    TFile* currMCFile = TFile::Open(currMCName.Data(),"read");
    if (!currDataFile || !currDataFile->IsOpen() || !currMCFile || !currMCFile->IsOpen()) {
      failedRuns += Form("%d ", runNumberMC);
      continue;
    }
    
    if (gSystem->AccessPathName(Form("displays/%d",runNumberMC)))
      gSystem->Exec(Form("mkdir displays/%d", runNumberMC));
    
    TObjArray *objsData = static_cast<TObjArray*>(currDataFile->Get("MUON_QA/expert"));
    TObjArray *objsMC = static_cast<TObjArray*>(currMCFile->Get("MUON_QA/expert"));
    
    for ( Int_t iCh = 1 ; iCh < 11 ; iCh++ )
    {
      TH1* clusterMapData = objsData ?
      static_cast<TH1*>(objsData->FindObject(Form("hClusterHitMapInCh%d",iCh))) :
      static_cast<TH1*>(currDataFile->Get(Form("MUON/ESDs/%s/Expert/%s_hESDClusterHitMap%d",
					       MultDirData.Data(),MultDirData.Data(),iCh)));
      TH1* clusterMapMC = objsMC ?
      static_cast<TH1*>(objsMC->FindObject(Form("hClusterHitMapInCh%d",iCh))) :
      static_cast<TH1*>(currMCFile->Get(Form("MUON/ESDs/%s/Expert/%s_hESDClusterHitMap%d",
					     MultDirMC.Data(),MultDirMC.Data(),iCh)));
      if (!clusterMapData || !clusterMapMC) {
	cout << "histograms not found for run " << runNumber << endl;
        failedRuns += Form("%d ", runNumberMC);
	break;
      }
      
      cTmp.cd(1);
      clusterMapData->Draw();
      
      cTmp.cd(2);
      clusterMapMC->Draw();

      cTmp.Print(Form("displays/%d/ESDclusterMapChamber%d.png",runNumber,iCh), "png");

    }
    
    currDataFile->Close();
    currMCFile->Close();
  }
  
  inDataFile.close();
  inMCFile.close();
  
  if (!failedRuns.IsNull()) {
    cout << "processing failed for run(s) " << failedRuns.Data() << endl;
  } else cout << "done                     " << endl;

}

//----------------------------------------------------------------------------
void CountTracks(TString inputFile)
{
  /// return the number of tracker track per run in the reconstruction QA
  
  ifstream inFile(inputFile.Data()); // open input file
  if (!inFile.is_open())
  {
    printf("cannot open file %s\n",inputFile.Data());
    return;
  }
  
  if (TString(gSystem->GetFromPipe(TString::Format("grep -c alien:// %s", inputFile.Data()))).Atoi() > 0) {
    TGrid::Connect("alien://");
  }
  
  TString currName;
  Int_t nRuns = 0;
  TArrayI runs(10000);
  Int_t nTracksTot = 0;
  TArrayI nTracks(10000);
  TList badFiles;
  badFiles.SetOwner(kTRUE);
  
  while (!inFile.eof()) //loop over elements in input file
  {
    currName.ReadLine(inFile,kTRUE); // Read line in input Data file
    if(currName.IsNull()) continue;
    
    Int_t runNumber(AliAnalysisMuonUtility::GetRunNumber(currName.Data()));
    
    cout<<"\rprocessing run "<<runNumber<<" ...\r"<<flush;
    
    TFile* currFile = TFile::Open(currName.Data(),"read");
    if (!currFile || !currFile->IsOpen()) {
      badFiles.AddLast(new TObjString(currName.Data()));
      continue;
    }
    
    TH1* hnClusters = 0x0;
    TObjArray *objs = static_cast<TObjArray*>(currFile->Get("MUON_QA/general1"));
    if (objs) hnClusters = static_cast<TH1*>(objs->FindObject("hNClustersPerTrack"));
    else hnClusters = static_cast<TH1*>(currFile->Get(Form("MUON/ESDs/%s/%s_hESDnClustersPerTrack",
                                                           MultDirData.Data(),MultDirData.Data())));
    
    if (!hnClusters) {
      currFile->Close();
      badFiles.AddLast(new TObjString(currName.Data()));
      continue;
    }
    
    runs[nRuns] = runNumber;
    nTracks[nRuns] = (Int_t)hnClusters->GetEntries();
    nTracksTot += nTracks[nRuns];
    nRuns++;
    
    currFile->Close();
  }
  
  cout << "done                     " << endl;
  inFile.close();
  
  Float_t percentTot = 0.;
  for (Int_t i = 0; i < nRuns; i++) {
    Float_t percent = 100.*((Float_t)nTracks[i])/((Float_t)nTracksTot);
    printf("%d %9d   %05.2f%%\n", runs[i], nTracks[i], percent);
    percentTot += percent;
  }
  printf("\ntotal number of tracks in the %d good runs = %d (%06.2f%%)\n", nRuns, nTracksTot, percentTot);
  
  if (badFiles.GetEntries() > 0) {
    printf("\nlist of bad files:\n");
    badFiles.Print();
  }
  
}

