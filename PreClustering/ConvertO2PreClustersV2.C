#include <fstream>
#include <vector>

#include <TROOT.h>
#include <TSystem.h>
#include <TString.h>
#include <TFile.h>
#include <TTree.h>

#include "AliMUONVCluster.h"
#include "AliMUONClusterStoreV2.h"

#include "MCHBase/PreCluster.h"
#include "DataFormatsMCH/Digit.h"

using namespace o2::mch;

void StorePreclusters(std::vector<PreCluster>& preClusters, std::vector<Digit>& digits, AliMUONClusterStoreV2* clusterStore, bool impl4);
uint32_t PadId2DigitId(int deId, int padId, bool impl4);

//------------------------------------------------------------------
void ConvertO2PreClustersV2(TString inFileName, TString outFileName = "preclusters.v2.root", bool impl4 = true)
{
  /// Convert O2 preclusters to AliRoot clusters
  /// input binary file with the following format:
  ///
  /// Number of preclusters
  /// Number of digits
  /// All PreClusters
  /// All Digits

  // load the digitId converter linked with the requested mapping implementation
  gSystem->Load(impl4 ? "libO2MCHMappingImpl4" : "libO2MCHMappingImpl3");
  gROOT->LoadMacro("/Users/PILLOT/Work/Alice/Macros/PreClustering/ConvertDigitId.C++");

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

  int nPreClusters(0);
  int nDigits(0);
  std::vector<PreCluster> preClusters{};
  std::vector<Digit> digits{};

  // loop over events, each event starting with the number of preclusters
  // (testing inFile.eof() does not work because the eof bit is not set until trying to read beyond the eof)
  while (inFile.read(reinterpret_cast<char*>(&nPreClusters),sizeof(int))) {
    
    // get the total number of digits
    inFile.read(reinterpret_cast<char*>(&nDigits),sizeof(int));

    // get the preclusters
    preClusters.resize(nPreClusters);
    inFile.read(reinterpret_cast<char*>(preClusters.data()), nPreClusters * sizeof(PreCluster));

    // get the digits
    digits.resize(nDigits);
    inFile.read(reinterpret_cast<char*>(digits.data()), nDigits * sizeof(Digit));

    // convert preclusters in AliRoot clusters
    StorePreclusters(preClusters, digits, clusterStore, impl4);

    clusterTree->Fill();
    clusterStore->Clear();
  }

  fileOut->cd();
  clusterTree->Write();
  fileOut->Close();
  inFile.close();
}

//------------------------------------------------------------------
void StorePreclusters(std::vector<PreCluster>& preClusters, std::vector<Digit>& digits, AliMUONClusterStoreV2* clusterStore, bool impl4)
{
  /// convert the preclusters in AliRoot clusters and store them in the cluster store

  int clIndex = clusterStore->GetSize();
  std::vector<uint32_t> digitIds(100);

  for (const auto& preCluster : preClusters) {

    int deId = digits[preCluster.firstDigit].getDetID();
    int chId = deId / 100 - 1;

    // store a new cluster (its ID has to be unique to add it to the new store)
    AliMUONVCluster* cluster = clusterStore->Add(chId, deId, clIndex);
    ++clIndex;

    // add the list of digit Ids
    digitIds.clear();
    for (auto iDigit = preCluster.firstDigit; iDigit <= preCluster.lastDigit(); ++iDigit) {
      uint32_t digitId = PadId2DigitId(deId, digits[iDigit].getPadID(), impl4);
      if (digitId > 0) {
        digitIds.push_back(digitId);
      }
    }
    cluster->SetDigitsId(preCluster.nDigits, digitIds.data());
  }
}
