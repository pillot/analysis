#include <fstream>
#include <vector>

#include <TError.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TString.h>
#include <TFile.h>
#include <TTree.h>

#include "AliCDBManager.h"
#include "AliGeomManager.h"

#include "AliMUONGeometryTransformer.h"
#include "AliMUONVCluster.h"
#include "AliMUONClusterStoreV2.h"

#include "MCHGeometryTransformer/Transformations.h"
#include "MathUtils/Cartesian.h"

#include "DataFormatsMCH/Cluster.h"
#include "DataFormatsMCH/Digit.h"

using namespace o2::mch;
using namespace std;

AliMUONGeometryTransformer alirootTransformer;
geo::TransformationCreator o2Transformer;

void StoreClusters(std::vector<Cluster>& clusters, std::vector<Digit>& digits, AliMUONClusterStoreV2* clusterStore, bool impl4);
uint32_t PadId2DigitId(int deId, int padId, bool impl4);

//------------------------------------------------------------------
void ConvertO2Clusters(TString inFileName, TString outFileName = "clusters.root", bool impl4 = true,
                       bool localtoGobal = true, TString geoFileName = "",
                       TString ocdb = "local://./OCDB", int run = 169099)
{
  /// convert O2 clusters to AliRoot clusters
  /// change clusters position from local to global coordinate system if requested,
  /// using provided geometry (with AliRoot if *.root or O2 if *.json) or loading it from OCDB (location or snapshot)
  /// the OCDB is also needed when using an AliRoot geometry file to get the mapping segmentation
  ///
  /// input binary file must have the following format:
  /// number of clusters
  /// number of digits
  /// all Clusters
  /// all Digits

  // load the digitId converter linked with the requested mapping implementation
  gSystem->Load(impl4 ? "libO2MCHMappingImpl4" : "libO2MCHMappingImpl3");
  gROOT->LoadMacro("/Users/PILLOT/Work/Alice/Macros/PreClustering/ConvertDigitId.C++");

  // load geometry to change cluster position from local to global coordinates, if needed
  if (localtoGobal) {
    if (geoFileName.IsNull()) {
      if (ocdb.EndsWith(".root")) {
        AliCDBManager::Instance()->SetDefaultStorage("local:///dev/null");
        AliCDBManager::Instance()->SetSnapshotMode(ocdb.Data());
      } else {
        AliCDBManager::Instance()->SetDefaultStorage(ocdb.Data());
      }
      AliCDBManager::Instance()->SetRun(run);
      AliGeomManager::LoadGeometry();
      if (!AliGeomManager::GetGeometry() || !AliGeomManager::ApplyAlignObjsFromCDB("MUON")) {
        Error("ConvertO2Clusters", "unable to load geometry and apply alignment from OCDB");
        return;
      }
      alirootTransformer.LoadGeometryData(); // also load the mapping segmentation
    } else if (gSystem->AccessPathName(geoFileName.Data(), kFileExists) != 0) {
      Error("ConvertO2Clusters", "geometry file %s not found", geoFileName.Data());
      return;
    } else if (geoFileName.EndsWith(".root")) {
      if (ocdb.EndsWith(".root")) {
        AliCDBManager::Instance()->SetDefaultStorage("local:///dev/null");
        AliCDBManager::Instance()->SetSnapshotMode(ocdb.Data());
      } else {
        AliCDBManager::Instance()->SetDefaultStorage(ocdb.Data());
      }
      AliCDBManager::Instance()->SetRun(run);
      AliGeomManager::LoadGeometry(geoFileName.Data());
      alirootTransformer.LoadGeometryData(); // also load the mapping segmentation
    } else if (geoFileName.EndsWith(".json")) {
      std::ifstream inGeoFile(geoFileName.Data());
      o2Transformer = geo::transformationFromJSON(inGeoFile);
    } else {
      Error("ConvertO2Clusters", "invalid geometry file");
      return;
    }
  }

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
  std::vector<Cluster> clusters{};
  std::vector<Digit> digits{};

  // loop over events, each event starting with the number of clusters
  // (testing inFile.eof() does not work because the eof bit is not set until trying to read beyond the eof)
  while (inFile.read(reinterpret_cast<char*>(&nClusters),sizeof(int))) {
    
    // get the total number of digits
    inFile.read(reinterpret_cast<char*>(&nDigits),sizeof(int));

    // get the clusters
    clusters.resize(nClusters);
    inFile.read(reinterpret_cast<char*>(clusters.data()), nClusters * sizeof(Cluster));

    // get the digits
    digits.resize(nDigits);
    inFile.read(reinterpret_cast<char*>(digits.data()), nDigits * sizeof(Digit));

    // convert clusters in AliRoot clusters
    StoreClusters(clusters, digits, clusterStore, impl4);

    clusterTree->Fill();
    clusterStore->Clear();
  }

  fileOut->cd();
  clusterTree->Write();
  fileOut->Close();
  inFile.close();
}

//------------------------------------------------------------------
void StoreClusters(std::vector<Cluster>& clusters, std::vector<Digit>& digits, AliMUONClusterStoreV2* clusterStore, bool impl4)
{
  /// convert the clusters in AliRoot clusters and store them in the cluster store
  /// change clusters position from local to global coordinate system if requested

  std::vector<uint32_t> digitIds(100);

  for (const auto& o2Cluster : clusters) {

    // store a new cluster (its ID has to be unique to add it to the new store)
    int deId = o2Cluster.getDEId();
    AliMUONVCluster* cluster = clusterStore->Add(o2Cluster.getChamberId(), deId, o2Cluster.getClusterIndex());

    // change the coordinate system if needed
    if (alirootTransformer.GetNofModuleTransformers() > 0) {
      double xg(0.), yg(0.), zg(0.);
      alirootTransformer.Local2Global(deId, o2Cluster.x, o2Cluster.y, 0., xg, yg, zg);
      cluster->SetXYZ(xg, yg, zg);
    } else if (o2Transformer) {
      o2::math_utils::Point3D<double> lpos{o2Cluster.x, o2Cluster.y, 0.};
      auto local2Global = o2Transformer(deId);
      auto gpos = local2Global(lpos);
      cluster->SetXYZ(gpos.x(), gpos.y(), gpos.z());
    } else {
      cluster->SetXYZ(o2Cluster.x, o2Cluster.y, o2Cluster.z);
    }

    // store cluster resolution
    cluster->SetErrXY(o2Cluster.ex, o2Cluster.ey);

    // add the list of digit Ids, if any
    if (!digits.empty()) {
      digitIds.clear();
      for (uint32_t iDigit = 0; iDigit < o2Cluster.nDigits; ++iDigit) {
        uint32_t digitId = PadId2DigitId(deId, digits[o2Cluster.firstDigit + iDigit].getPadID(), impl4);
        if (digitId > 0) {
          digitIds.push_back(digitId);
        }
      }
      cluster->SetDigitsId(o2Cluster.nDigits, digitIds.data());
    }
  }
}
