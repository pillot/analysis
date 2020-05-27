#include <fstream>
#include <vector>

#include <TString.h>
#include <TFile.h>
#include <TTree.h>

#include "AliCDBManager.h"
#include "AliGeomManager.h"

#include "AliMUONGeometryTransformer.h"
#include "AliMUONVCluster.h"
#include "AliMUONClusterStoreV2.h"

#include "MCHBase/ClusterBlock.h"
#include "MCHBase/Digit.h"

using namespace o2::mch;

AliMUONGeometryTransformer transformer;

void StoreClusters(std::vector<ClusterStruct>& clusters, std::vector<Digit>& digits, AliMUONClusterStoreV2* clusterStore);

//------------------------------------------------------------------
void ConvertO2Clusters(TString inFileName, TString outFileName = "clusters.root")
{
  /// Convert O2 clusters to AliRoot clusters
  /// input binary file with the following format:
  ///
  /// Number of clusters
  /// Number of digits
  /// All Clusters
  /// All Digits

  // load geometry to change cluster position from local to global coordinates
  //  AliCDBManager::Instance()->SetDefaultStorage("raw://");
  AliCDBManager::Instance()->SetDefaultStorage("local://./OCDB");
  AliCDBManager::Instance()->SetRun(196099);
  AliGeomManager::LoadGeometry();
  if (!AliGeomManager::GetGeometry() || !AliGeomManager::ApplyAlignObjsFromCDB("MUON")) return;
  transformer.LoadGeometryData(); // also load the mapping segmentation

  // open input file
  ifstream inFile(inFileName,ios::binary);
  if (!inFile.is_open()) return;

  // prepare storage of clusters in an AliMUONVCluster format
  TFile* fileOut = new TFile(outFileName.Data(),"RECREATE");
  if (!fileOut) return;
  TTree* clusterTree = new TTree("TreeR","Clusters");
  if (!clusterTree) return;
  AliMUONClusterStoreV2* clusterStore = new AliMUONClusterStoreV2();
  if (!clusterStore) return;
  if (!clusterStore->Connect(*clusterTree, kTRUE)) return;

  int nClusters(0);
  int nDigits(0);
  std::vector<ClusterStruct> clusters{};
  std::vector<Digit> digits{};

  // loop over events, each event starting with the number of clusters
  // (testing inFile.eof() does not work because the eof bit is not set until trying to read beyond the eof)
  while (inFile.read(reinterpret_cast<char*>(&nClusters),sizeof(int))) {
    
    // get the total number of digits
    inFile.read(reinterpret_cast<char*>(&nDigits),sizeof(int));

    // get the clusters
    clusters.resize(nClusters);
    inFile.read(reinterpret_cast<char*>(clusters.data()), nClusters * sizeof(ClusterStruct));

    // get the digits
    digits.resize(nDigits);
    inFile.read(reinterpret_cast<char*>(digits.data()), nDigits * sizeof(Digit));

    // convert clusters in AliRoot clusters
    StoreClusters(clusters, digits, clusterStore);

    clusterTree->Fill();
    clusterStore->Clear();
  }

  fileOut->cd();
  clusterTree->Write();
  fileOut->Close();
  inFile.close();
}

//------------------------------------------------------------------
void StoreClusters(std::vector<ClusterStruct>& clusters, std::vector<Digit>& digits, AliMUONClusterStoreV2* clusterStore)
{
  /// convert the clusters in AliRoot clusters and store them in the cluster store

  int clIndex = clusterStore->GetSize();
  std::vector<uint32_t> digitIds(100);

  for (const auto& o2Cluster : clusters) {

    int deId = digits[o2Cluster.firstDigit].getDetID();
    int chId = deId / 100 - 1;

    // store a new cluster (its ID has to be unique to add it to the new store)
    AliMUONVCluster* cluster = clusterStore->Add(chId, deId, clIndex);
    ++clIndex;

    // store the position in the global coordinate system
    double xg(0.), yg(0.), zg(0.);
    transformer.Local2Global(deId, o2Cluster.x, o2Cluster.y, 0., xg, yg, zg);
    cluster->SetXYZ(xg, yg, zg);

    // add the list of digit Ids
    digitIds.clear();
    for (uint32_t iDigit = 0; iDigit <= o2Cluster.nDigits; ++iDigit) {
      digitIds.push_back(static_cast<uint32_t>(digits[o2Cluster.firstDigit + iDigit].getPadID()));
    }
    cluster->SetDigitsId(o2Cluster.nDigits, digitIds.data());
  }
}
