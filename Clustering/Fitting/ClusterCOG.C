#include <cmath>
#include <tuple>
#include <utility>
#include <vector>

#include <TCanvas.h>
#include <TFile.h>
#include <TStyle.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

#include "DataFormatsMCH/Cluster.h"
#include "DataFormatsMCH/Digit.h"
#include "MCHBase/TrackBlock.h"

#include "CCDBUtils.h"
#include "ClusterUtils.h"
#include "DataUtils.h"
#include "PreClusterUtils.h"

using o2::mch::Cluster;
using o2::mch::Digit;
using o2::mch::TrackParamStruct;

//_________________________________________________________________________________________________
void ClusterCOG(int run, const char* inFile = "clusters.root", const char* outFile = "newclusters.root")
{
  /// compute the center of gravity of the digits attached to the selected clusters
  /// store the new clusters together with the corresponding input data in outFile
  /// require the MCH mapping to be loaded: gSystem->Load("libO2MCHMappingImpl4")

  /// load CCDB objects
  InitFromCCDB(run, true, true, false);

  // load input data
  auto [dataFileIn, dataReader] = LoadData(inFile, "data");
  TTreeReaderValue<TrackParamStruct> trackParam(*dataReader, "trackParameters");
  TTreeReaderValue<Cluster> cluster(*dataReader, "clusters");
  TTreeReaderValue<std::vector<Digit>> digits(*dataReader, "digits");

  // setup the output
  TFile dataFileOut(outFile, "recreate");
  TTree* dataTreeOut = dataReader->GetTree()->CloneTree(0);
  dataTreeOut->SetTitle("tree with input and output data");
  Cluster newCluster;
  dataTreeOut->Branch("newClusters", &newCluster);

  std::vector<TH1*> preClusterInfo{};
  CreatePreClusterInfo(preClusterInfo);
  std::vector<TH1*> preClusterInfoSt[3] = {{}, {}, {}};
  CreatePreClusterInfo(preClusterInfoSt[0], "St1");
  CreatePreClusterInfo(preClusterInfoSt[1], "St2");
  CreatePreClusterInfo(preClusterInfoSt[2], "St345");

  int nClusters = dataReader->GetEntries(false);
  int iCluster(0);
  while (dataReader->Next()) {
    if (++iCluster % 1000 == 0) {
      std::cout << "\rprocessing cluster " << iCluster << " / " << nClusters << "..." << std::flush;
    }

    // those 2 DE have lower HV for the run 529691
    // if (cluster->getDEId() == 202 || cluster->getDEId() == 300) {
    //   continue;
    // }

    // check precluster characteristics (charge, size, ...)
    const auto [sizeX, sizeY] = GetSize(*digits);
    const auto [chargeNB, chargeB] = GetCharge(*digits);
    double chargeAsymm = (chargeNB - chargeB) / (chargeNB + chargeB);

    // apply further precluster selections (e.g. to have enough constraints for the fit)
    if (cluster->nDigits < 4 || sizeY < 2 || std::abs(chargeAsymm) > 0.5) {
      continue;
    }

    // fill selected precluster characteristics
    FillPreClusterInfo(chargeNB, chargeB, sizeX, sizeY, preClusterInfo);
    int iSt = (cluster->getChamberId() < 4) ? cluster->getChamberId() / 2 : 2;
    FillPreClusterInfo(chargeNB, chargeB, sizeX, sizeY, preClusterInfoSt[iSt]);

    // refit the cluster
    const auto [x, y] = (cluster->getChamberId() < 4) ? GetCOG2(*digits) : GetCOG(*digits);
    newCluster = MakeCluster(cluster->uid, x, y);

    // fill output tree
    *trackParam; // this is just to trigger the reading/copying of track parameters
    dataTreeOut->Fill();
  }
  cout << "\r\033[Kprocessing completed" << endl;

  gStyle->SetOptStat(1);

  auto c = DrawPreClusterInfo(preClusterInfo);
  auto cSt1 = DrawPreClusterInfo(preClusterInfoSt[0], "St1");
  auto cSt2 = DrawPreClusterInfo(preClusterInfoSt[1], "St2");
  auto cSt345 = DrawPreClusterInfo(preClusterInfoSt[2], "St345");

  dataFileOut.Write();
  c->Write();
  cSt1->Write();
  cSt2->Write();
  cSt345->Write();
  dataFileOut.Close();
  dataFileIn->Close();
}
