#include <chrono>

#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>

#include "AliCDBManager.h"
#include "AliGeomManager.h"

#include "AliMpArea.h"

#include "AliMUONCDB.h"
#include "AliMUONGeometryTransformer.h"
#include "AliMUONRecoParam.h"
#include "AliMUONPreClusterFinder.h"
#include "AliMUONSimpleClusterServer.h"
#include "AliMUONVDigitStore.h"
#include "AliMUONClusterStoreV2.h"

//------------------------------------------------------------------
void RunAliRootClusterizer(const char* digitFileName, const char* clusterFileName, int event = -1)
{
  /// run the Aliroot clustering

  // load geometry and reconstruction parameters
//  AliCDBManager::Instance()->SetDefaultStorage("raw://");
  AliCDBManager::Instance()->SetDefaultStorage("local://./OCDB");
  AliCDBManager::Instance()->SetRun(196099);
  AliGeomManager::LoadGeometry();
  if (!AliGeomManager::GetGeometry() || !AliGeomManager::ApplyAlignObjsFromCDB("MUON")) return;
  AliMUONGeometryTransformer transformer;
  transformer.LoadGeometryData(); // also load the mapping segmentation
  AliMUONRecoParam* recoParam = AliMUONCDB::LoadRecoParam();
  if (!recoParam) return;

  // prepare the precluster finder
  gROOT->LoadMacro("$WORK/Macros/PreClustering/AliMUONPreClusterFinderV4.cxx+");
  AliMUONVClusterFinder* clusterFinder = reinterpret_cast<AliMUONVClusterFinder*>(gROOT->ProcessLineSync("new AliMUONClusterFinderMLEM(false, new AliMUONPreClusterFinderV4())"));
  AliMUONVClusterServer* clusterServer = new AliMUONSimpleClusterServer(clusterFinder, transformer);

  // prepare to read digits
  TFile* digitFile = TFile::Open(digitFileName);
  if (!digitFile) return;
  TTree* treeD = static_cast<TTree*>(digitFile->Get("TreeD"));
  if (!treeD) return;
  AliMUONVDigitStore* digitStore = AliMUONVDigitStore::Create(*treeD);
  if (!digitStore->Connect(*treeD)) return;

  // prepare to write clusters
  AliMUONVClusterStore* clusterStore = new AliMUONClusterStoreV2();
  TFile* clusterFile = TFile::Open(clusterFileName,"RECREATE");
  TTree* treeR = new TTree("TreeR","Cluster Container");
  clusterStore->Connect(*treeR,kTRUE);

  // loop over events
  AliMpArea area; // invalid area to clusterize everything
  int nEvents = treeD->GetEntries();
  std::chrono::duration<double> clusteringTime{};
  for (int iEvent = 0; iEvent < nEvents; ++iEvent) {

    if (event >= 0 && iEvent != event) {
      continue;
    }

    // load the digits
    treeD->GetEntry(iEvent);
    TIter next(digitStore->CreateIterator());
    clusterServer->UseDigits(next);

    // clusterize every chambers
    auto tStart = std::chrono::high_resolution_clock::now();
    for (int iCh = 0; iCh < 10; ++iCh) {
      clusterServer->Clusterize(iCh, *clusterStore, area, recoParam);
    }
    auto tEnd = std::chrono::high_resolution_clock::now();
    clusteringTime += tEnd - tStart;

    // store the clusters
    treeR->Fill();
    
    // clear the stores
    clusterStore->Clear();
    digitStore->Clear();
  }

  // write the clusters
  clusterFile->cd();
  treeR->Write();
  clusterFile->Close();

  // print timer
  cout << "clustering duration = " << clusteringTime.count() << " s" << endl;
}
