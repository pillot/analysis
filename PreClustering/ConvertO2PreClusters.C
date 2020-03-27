#include <fstream>
#include <vector>

#include <TString.h>
#include <TFile.h>
#include <TTree.h>

#include "AliMUONVCluster.h"
#include "AliMUONClusterStoreV2.h"

struct DataBlockHeader {
   uint16_t fType;        // The type of the data block. Must contain a value defined by DataBlockType.
   uint16_t fRecordWidth; // The number of bytes each record uses.
   uint32_t fNrecords;    // Number of records in this data block.
 };

  struct DigitBlock {
   DataBlockHeader header; // Common data block header
 };

  struct DigitStruct {
   uint32_t uid;   // Digit ID in the current mapping (from OCDB)
   uint16_t index; // Digit index in the new mapping (produced internally)
   uint16_t adc;   // ADC value of signal
 };
 
#include "PreClusterBlock.cxx"

void StorePreclusters(int deId, const PreClusterBlock& preClusterBlock, AliMUONClusterStoreV2* clusterStore);

//------------------------------------------------------------------
void ConvertO2PreClusters(TString inFileName, TString outFileName = "preclusters.root")
{
  /// Convert O2 preclusters (old format) to AliRoot clusters
  /// input binary file with the following format:
  ///
  /// Number of DE with preclusters
  /// DE ID
  /// PreClusterBlock
  /// DE ID
  /// PreClusterBlock
  /// ...
  /// Number of DE with preclusters
  /// ...

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

  int nDEWithPreClusters(0);
  int deId(0);
  uint32_t size(0);
  std::vector<char> buffer{};
  PreClusterBlock preClusterBlock{};

  // loop over events, each event starting with the number of DE with preclusters
  // (testing inFile.eof() does not work because the eof bit is not set until trying to read beyond the eof)
  while (inFile.read(reinterpret_cast<char*>(&nDEWithPreClusters),sizeof(int))) {
    
    for (int iDE = 0, nDEs = nDEWithPreClusters; iDE < nDEs; ++iDE) {
      
      // get the DE ID
      inFile.read(reinterpret_cast<char*>(&deId),sizeof(int));

      // get the size of the PreClusterBlock
      inFile.read(reinterpret_cast<char*>(&size),sizeof(int));

      // get the preclusters
      buffer.reserve(size);
      inFile.read(buffer.data(), size);
      if (preClusterBlock.reset(buffer.data(), size, false) < 0) return;
      if (preClusterBlock.getCurrentSize() != size) return;

      // convert them in AliRoot clusters
      StorePreclusters(deId, preClusterBlock, clusterStore);
    }

    clusterTree->Fill();
    clusterStore->Clear();
  }

  fileOut->cd();
  clusterTree->Write();
  fileOut->Close();
  inFile.close();
}

//------------------------------------------------------------------
void StorePreclusters(int deId, const PreClusterBlock& preClusterBlock, AliMUONClusterStoreV2* clusterStore)
{
  /// convert the preclusters in AliRoot clusters and store them in the cluster store

  AliMUONVCluster *cluster(0x0);
  int clIndex = clusterStore->GetSize();
  int chId(deId/100-1);
  const auto& preclusters(preClusterBlock.getPreClusters());
  std::vector<uint32_t> digitIds(100);

  for (const auto& precluster : preclusters) {

    // store a new cluster (its ID has to be unique to add it to the new store)
    cluster = clusterStore->Add(chId, deId, clIndex);
    ++clIndex;

    // add the list of digit Ids
    digitIds.clear();
    for (uint16_t iDigit = 0; iDigit < precluster.nDigits; ++iDigit) {
      digitIds.push_back(precluster.digits[iDigit].uid);
    }
    cluster->SetDigitsId(precluster.nDigits, digitIds.data());
  }
}
