#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TFile.h>
#include <TTree.h>
#include <TCollection.h>
#include <Riostream.h>
#include <TROOT.h>
#include <TObjectTable.h>

// STEER includes
#include "AliESDEvent.h"
#include "AliCDBManager.h"

// MUON includes
#include "AliMUONCDB.h"
#include "AliMUONRecoParam.h"
#include "AliMUONESDInterface.h"
#include "AliESDMuonTrack.h"
#include "AliMUONTrack.h"
#include "AliMUONVCluster.h"
#include "AliMUONVDigit.h"
#endif


TTree* GetESDTree(TFile *esdFile);


//-----------------------------------------------------------------------
void CheckESDInterface(Int_t iEv = 0)
{
  
  // open the ESD file and tree
  TFile* esdFile = TFile::Open("AliESDs.root");
  TTree* esdTree = GetESDTree(esdFile);
  
  // connect ESD event to the ESD tree
  AliESDEvent* esd = new AliESDEvent();
  esd->ReadFromTree(esdTree);
  
  // get the ESD of the given event
  if (iEv >= (Int_t)esdTree->GetEntries()) iEv = 0;
  if (esdTree->GetEvent(iEv) <= 0) {
    Error("CheckESDInterface", Form("no ESD object found for event %d",iEv));
    return;
  }
  
  // load necessary data from OCDB
  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  AliCDBManager::Instance()->SetRun(esd->GetRunNumber());
  if (!AliMUONCDB::LoadField()) return;
  if (!AliMUONCDB::LoadMapping(kTRUE)) return;
  AliMUONRecoParam* recoParam = AliMUONCDB::LoadRecoParam();
  if (!recoParam) return;
  
  // reset tracker for track restoring initial track parameters at cluster
  AliMUONESDInterface::ResetTracker(recoParam);
  
  // ESD interface
  AliMUONESDInterface *esdInterface = new AliMUONESDInterface();
  esdInterface->LoadEvent(*esd);
  esdInterface->LoadEvent(*esd);
  esdInterface->LoadEvent(*esd);
  esdInterface->LoadEvent(*esd);
  esdInterface->LoadEvent(*esd);
  esdInterface->LoadEvent(*esd);
  esdInterface->LoadEvent(*esd);
  esdInterface->LoadEvent(*esd);
  esdInterface->LoadEvent(*esd);
  esdInterface->LoadEvent(*esd);
  esdInterface->LoadEvent(*esd);
  esdInterface->LoadEvent(*esd);
  
  // check #1
  cout<<"\t----------------------check #1---------------------\t"<<endl;
  Int_t nTracks = 0;
  UInt_t trackId = 0;
  AliMUONTrack *track;
  TIter nextTrack(esdInterface->CreateTrackIterator());
  while ((track = (AliMUONTrack*)nextTrack())) {
    track->Print();
    trackId = track->GetUniqueID();
    nTracks++;
  }
  cout<< "printed tracks: "<<nTracks<<"/"<<esdInterface->GetNTracks()<<endl;
  esdInterface->FindTrack(trackId)->Print();
  cout<<"\t---------------------------------------------------\t"<<endl<<endl;
  
  // check #2
  cout<<"\t----------------------check #2---------------------\t"<<endl;
  Int_t nClusters = 0;
  UInt_t clusterId = 0;
  AliMUONVCluster *cluster;
  TIter nextCluster(esdInterface->CreateClusterIterator());
  while ((cluster = (AliMUONVCluster*)nextCluster())) {
    cluster->Print();
    clusterId = cluster->GetUniqueID();
    nClusters++;
  }
  cout<< "printed clusters: "<<nClusters<<"/"<<esdInterface->GetNClusters()<<endl;
  esdInterface->FindCluster(clusterId)->Print();
  cout<<"\t---------------------------------------------------\t"<<endl<<endl;
  
  // check #3
  cout<<"\t----------------------check #3---------------------\t"<<endl;
  nClusters = 0;
  TIter nextCluster2(esdInterface->CreateClusterIterator(trackId));
  while ((cluster = (AliMUONVCluster*)nextCluster2())) {
    cluster->Print();
    clusterId = cluster->GetUniqueID();
    nClusters++;
  }
  cout<< "printed clusters: "<<nClusters<<"/"<<esdInterface->GetNClusters(trackId)<<endl;
  esdInterface->FindCluster(trackId, clusterId)->Print();
  cout<<"\t---------------------------------------------------\t"<<endl<<endl;
  
  // check #4
  cout<<"\t----------------------check #4---------------------\t"<<endl;
  Int_t nDigits = 0;
  UInt_t digitId = 0;
  AliMUONVDigit *digit;
  TIter nextDigit(esdInterface->CreateDigitIterator());
  while ((digit = (AliMUONVDigit*)nextDigit())) {
    digit->Print();
    digitId = digit->GetUniqueID();
    nDigits++;
  }
  cout<< "printed clusters: "<<nDigits<<"/"<<esdInterface->GetNDigits()<<endl;
  esdInterface->FindDigit(digitId)->Print();
  cout<<"\t---------------------------------------------------\t"<<endl<<endl;
  
  // check #5
  cout<<"\t----------------------check #5---------------------\t"<<endl;
  nDigits = 0;
  TIter nextDigit2(esdInterface->CreateDigitIterator(trackId));
  while ((digit = (AliMUONVDigit*)nextDigit2())) {
    digit->Print();
    digitId = digit->GetUniqueID();
    nDigits++;
  }
  cout<< "printed clusters: "<<nDigits<<"/"<<esdInterface->GetNDigits(trackId)<<endl;
  esdInterface->FindDigit(digitId)->Print();
  cout<<"\t---------------------------------------------------\t"<<endl<<endl;
  
  // check #6
  cout<<"\t----------------------check #6---------------------\t"<<endl;
  nDigits = 0;
  TIter nextDigit3(esdInterface->CreateDigitIterator(trackId, clusterId));
  while ((digit = (AliMUONVDigit*)nextDigit3())) {
    digit->Print();
    nDigits++;
  }
  cout<< "printed clusters: "<<nDigits<<"/"<<esdInterface->GetNDigits(trackId, clusterId)<<endl;
  esdInterface->FindDigit(digitId)->Print();
  cout<<"\t---------------------------------------------------\t"<<endl<<endl;
  
  // check #7
  cout<<"\t----------------------check #7---------------------\t"<<endl;
  nDigits = 0;
  TIter nextDigit4(esdInterface->CreateDigitIteratorInCluster(clusterId));
  while ((digit = (AliMUONVDigit*)nextDigit4())) {
    digit->Print();
    nDigits++;
  }
  cout<< "printed clusters: "<<nDigits<<"/"<<esdInterface->GetNDigitsInCluster(clusterId)<<endl;
  esdInterface->FindDigit(esdInterface->FindCluster(clusterId)->GetDigitId(nDigits-1))->Print();
  cout<<"\t---------------------------------------------------\t"<<endl<<endl;
  
  delete esdInterface;
  gObjectTable->Print();
}

//-----------------------------------------------------------------------
TTree* GetESDTree(TFile *esdFile)
{
  /// Check that the file is properly open
  /// Return pointer to the ESD Tree
  
  if (!esdFile || !esdFile->IsOpen()) {
    Error("GetESDTree", "opening ESD file failed");
    exit(-1);
  }
  
  TTree* tree = (TTree*) esdFile->Get("esdTree");
  if (!tree) {
    Error("GetESDTree", "no ESD tree found");
    exit(-1);
  }
  
  return tree;
  
}

