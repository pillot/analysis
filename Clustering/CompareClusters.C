#include <stdio.h>
#include <list>
#include <map>
#include <unordered_map>
#include <utility>
#include <tuple>

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>

#include "AliMUONVCluster.h"
#include "AliMUONVClusterStore.h"

struct preCluster {
  preCluster(AliMUONVCluster* cluster, std::list<u_int32_t>& digits) : clusters{cluster}, digitIds(std::move(digits)) {}

  std::list<AliMUONVCluster*> clusters{}; // link to the corresponding clusters
  std::list<u_int32_t> digitIds{};        // list of digit Ids
  bool compareMe = true;                  // kFALSE if already associate to an identical precluster
};

const int nDEs = 156;
std::unordered_map<int, int> deIndices(0);
TH1F hDiff("hDiff", "distance between clusters", 100001, -0.000005, 1.000005);
double kPrecision = 0.;

int LoadClusters(AliMUONVClusterStore* clusterStore, std::list<preCluster>* preclusters, bool perPrecluster);
std::tuple<int, int, int> Compare(std::list<preCluster>* preclusters1, std::list<preCluster>* preclusters2);
std::tuple<int, int, int> Compare(std::list<AliMUONVCluster*>& clusters1, std::list<AliMUONVCluster*>& clusters2);
preCluster* FindPrecluster(std::list<u_int32_t>& digitIds, std::list<preCluster>& preclusters);

//------------------------------------------------------------------
void CompareClusters(const char *clusterFileName1, const char *clusterFileName2, double precision = 1.e-3,
                     bool perPrecluster = true)
{
  /// Compare the clusters DE per DE between the two input files
  /// - precision: print the difference between clusters if it higher than the given precision
  /// - perPrecluster: compare clusters from identical preclusters or from all preclusters

  kPrecision = precision;

  // input clusters 1
  TFile* clusterFile1 = new TFile(clusterFileName1);
  if (!clusterFile1) return;
  TTree* treeR1 = static_cast<TTree*>(clusterFile1->Get("TreeR"));
  if (!treeR1) return;
  AliMUONVClusterStore *clusterStore1 = AliMUONVClusterStore::Create(*treeR1);
  clusterStore1->Connect(*treeR1);
  
  // input clusters 2
  TFile* clusterFile2 = new TFile(clusterFileName2);
  if (!clusterFile2) return;
  TTree* treeR2 = static_cast<TTree*>(clusterFile2->Get("TreeR"));
  if (!treeR2) return;
  AliMUONVClusterStore *clusterStore2 = AliMUONVClusterStore::Create(*treeR2);
  clusterStore2->Connect(*treeR2);
  
  std::list<preCluster> preclusters1[nDEs];
  std::list<preCluster> preclusters2[nDEs];
  deIndices.reserve(nDEs);
  int integratedNCl1(0);
  int integratedNCl2(0);
  int integratedNDiff1(0);
  int integratedNDiff2(0);
  int integratedNDiff12(0);

  // number of events to process
  int64_t nEvents = treeR1->GetEntries();
  if (treeR2->GetEntries() != nEvents) {
    printf("Warning: not the same number of events in the two cluster trees --> try with the smallest one\n");
    nEvents = TMath::Min(nEvents, treeR2->GetEntries());
  }

  for (int64_t iEv = 0; iEv < nEvents; ++iEv) {

    printf("Event %lld:\n", iEv);
    
    treeR1->GetEntry(iEv);
    treeR2->GetEntry(iEv);

    int nClusters1 = LoadClusters(clusterStore1, preclusters1, perPrecluster);
    int nClusters2 = LoadClusters(clusterStore2, preclusters2, perPrecluster);
    printf("\tTotal number of clusters = %d / %d\n", nClusters1, nClusters2);
    integratedNCl1 += nClusters1;
    integratedNCl2 += nClusters2;

    auto nDiffClusters = Compare(preclusters1, preclusters2);
    printf("Total number of isolated clusters = %d / %d\n\n", std::get<0>(nDiffClusters), std::get<1>(nDiffClusters));
    integratedNDiff1 += std::get<0>(nDiffClusters);
    integratedNDiff2 += std::get<1>(nDiffClusters);
    integratedNDiff12 += std::get<2>(nDiffClusters);

    for (int iDE = 0; iDE < nDEs; ++iDE) {
      preclusters1[iDE].clear();
      preclusters2[iDE].clear();
    }
    clusterStore1->Clear();
    clusterStore2->Clear();
  }

  printf("Integrated number of clusters = %d / %d\n", integratedNCl1, integratedNCl2);
  printf("Integrated number of isolated clusters = %d / %d\n", integratedNDiff1, integratedNDiff2);
  printf("Integrated number of associated clusters with different position = %d\n\n", integratedNDiff12);

  new TCanvas;
  hDiff.Draw();
}

//------------------------------------------------------------------
int LoadClusters(AliMUONVClusterStore* clusterStore, std::list<preCluster>* preclusters, bool perPrecluster)
{
  /// fill the precluster structures with reconstructed clusters per DE
  /// if perPrecluster = false, attach all clusters to the same dummy precluster without digit

  int nClusters(0);
  AliMUONVCluster* cluster(nullptr);
  TIter nextCluster(clusterStore->CreateIterator());
  while ((cluster = static_cast<AliMUONVCluster*>(nextCluster()))) {

    // get the DE index
    int deId = cluster->GetDetElemId();
    int iDE = deIndices.size();
    auto itDE = deIndices.find(deId);
    if (itDE == deIndices.end()) {
      deIndices.emplace(deId, iDE);
    } else {
      iDE = itDE->second;
    }

    // get the list of digits if the comparison is done per precluster
    std::list<u_int32_t> digitIds;
    if (perPrecluster) {
      for (int iDigit = 0; iDigit < cluster->GetNDigits(); ++iDigit) {
        digitIds.push_back(cluster->GetDigitId(iDigit));
      }
      digitIds.sort();
    }

    // find the precluster with the same digits, or create it, and attached the cluster to it
    auto* precluster = FindPrecluster(digitIds, preclusters[iDE]);
    if (precluster == nullptr) {
      preclusters[iDE].emplace_back(cluster, digitIds);
    } else {
      precluster->clusters.push_back(cluster);
    }

    ++nClusters;
  }
  
  return nClusters;
}

//------------------------------------------------------------------
std::tuple<int, int, int> Compare(std::list<preCluster>* preclusters1, std::list<preCluster>* preclusters2)
{
  /// For every DE: compare every clusters attached to every preclusters between the two lists

  int nDiffCluster1(0);
  int nDiffCluster2(0);
  int nDiffCluster12(0);

  for (size_t iDE = 0; iDE < deIndices.size(); iDE++) {

    // find identical preclusters in both lists and compare their attached clusters
    for (auto& precluster1 : preclusters1[iDE]) {
      auto* precluster2 = FindPrecluster(precluster1.digitIds, preclusters2[iDE]);
      if (precluster2 == nullptr) {
        nDiffCluster1 += precluster1.clusters.size();
        printf("DE: %d, precluster %u containing %lu clusters is missing in file 2\n",
               precluster1.clusters.front()->GetDetElemId(),
               precluster1.clusters.front()->GetDigitId(0),
               precluster1.clusters.size());
      } else {
        precluster2->compareMe = false;
        auto nDiffClusters = Compare(precluster1.clusters, precluster2->clusters);
        nDiffCluster1 += std::get<0>(nDiffClusters);
        nDiffCluster2 += std::get<1>(nDiffClusters);
        nDiffCluster12 += std::get<2>(nDiffClusters);
      }
    }

    // check for remaining clusters in the second list
    for (auto& precluster2 : preclusters2[iDE]) {
      if (precluster2.compareMe) {
        nDiffCluster2 += precluster2.clusters.size();
        printf("DE: %d, precluster %u containing %lu clusters is missing in file 1\n",
               precluster2.clusters.front()->GetDetElemId(),
               precluster2.clusters.front()->GetDigitId(0),
               precluster2.clusters.size());
      }
    }
  }

  return std::make_tuple(nDiffCluster1, nDiffCluster2, nDiffCluster12);
}

//------------------------------------------------------------------
std::tuple<int, int, int> Compare(std::list<AliMUONVCluster*>& clusters1, std::list<AliMUONVCluster*>& clusters2)
{
  /// compare every clusters between the two lists

  // make every combinations of clusters between both lists ordered per distance between clusters
  std::multimap<double, std::pair<AliMUONVCluster*, AliMUONVCluster*>> clusterPairs{};
  for (auto cluster1 : clusters1) {
    for (auto cluster2 : clusters2) {
      double dx = cluster2->GetX() - cluster1->GetX();
      double dy = cluster2->GetY() - cluster1->GetY();
      clusterPairs.emplace(sqrt(dx * dx + dy * dy), std::make_pair(cluster1, cluster2));
    }
  }

  // find the pair of clusters closest to each others (use BIT(20) to tag clusters already associated to another)
  int nDiffCluster12(0);
  for (auto& clusterPair : clusterPairs) {
    if (!clusterPair.second.first->TestBit(BIT(20)) && !clusterPair.second.second->TestBit(BIT(20))) {
      clusterPair.second.first->SetBit(BIT(20));
      clusterPair.second.second->SetBit(BIT(20));
      hDiff.Fill(clusterPair.first);
      if (clusterPair.first > kPrecision) {
        ++nDiffCluster12;
        printf("DE: %d, precluster: %u / %u: distance between clusters = %f\n",
               clusterPair.second.first->GetDetElemId(),
               clusterPair.second.first->GetDigitId(0),
               clusterPair.second.second->GetDigitId(0),
               clusterPair.first);
      }
    }
  }

  // count the number of clusters in first list not assocated with another
  int nDiffCluster1(0);
  for (auto cluster1 : clusters1) {
    if (!cluster1->TestBit(BIT(20))) {
      ++nDiffCluster1;
    }
  }
  if (nDiffCluster1 > 0) {
    printf("DE: %d, precluster: %u: %d clusters missing in file 2\n",
           clusters1.front()->GetDetElemId(),
           clusters1.front()->GetDigitId(0),
           nDiffCluster1);
  }

  // count the number of clusters in second list not assocated with another
  int nDiffCluster2(0);
  for (auto cluster2 : clusters2) {
    if (!cluster2->TestBit(BIT(20))) {
      ++nDiffCluster2;
    }
  }
  if (nDiffCluster2 > 0) {
    printf("DE: %d, precluster: %u: %d clusters missing in file 1\n",
           clusters2.front()->GetDetElemId(),
           clusters2.front()->GetDigitId(0),
           nDiffCluster2);
  }

  return std::make_tuple(nDiffCluster1, nDiffCluster2, nDiffCluster12);
}

//------------------------------------------------------------------
preCluster* FindPrecluster(std::list<u_int32_t>& digitIds, std::list<preCluster>& preclusters)
{
  /// find the precluster with the exact same list of digits if any

  for (auto& precluster : preclusters) {
    if (precluster.compareMe && precluster.digitIds == digitIds) {
      return &precluster;
    }
  }

  return nullptr;
}
