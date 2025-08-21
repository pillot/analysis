#include <algorithm>
#include <array>
#include <cmath>
#include <chrono>
#include <string>
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
#include "FitUtils.h"
#include "PreClusterUtils.h"

using o2::mch::Cluster;
using o2::mch::Digit;
using o2::mch::TrackParamStruct;

static constexpr double pi = 3.14159265358979323846;

//_________________________________________________________________________________________________
void ClusterFit(int run, bool fitAsymm = true, std::string errorMode = "MLS", double errorAlpha = 1.,
                double k3x = 0.3, double k3y = 0.3, bool correctCharge = false,
                std::array<int, 6> fix = {0, 0, 1, 1, 0, 0},
                std::string inFile = "clusters.root", std::string outFile = "newclusters.root")
{
  /// fit the digits attached to the selected clusters with Mathieson functions
  /// fitAsymm = true --> fit bending and non-bending cluster charge separately
  /// errorMode = "LS", "MLS", "MC" or "const", combined with errorAlpha (see FitUtils.C)
  /// store the new clusters together with the corresponding input data in outFile
  /// require the MCH mapping to be loaded: gSystem->Load("libO2MCHMappingImpl4")

  if (errorMode != "LS" && errorMode != "MLS" && errorMode != "MC" && errorMode != "const") {
    LOGP(error, "unknown error mode. Must be \"LS\", \"MLS\", \"MC\" or \"const\"");
    exit(-1);
  }

  // load CCDB objects
  InitFromCCDB(run, true, true, false);

  // load input data
  auto [dataFileIn, dataReader] = LoadData(inFile.c_str(), "data");
  TTreeReaderValue<TrackParamStruct> trackParam(*dataReader, "trackParameters");
  TTreeReaderValue<int> trackTime(*dataReader, "trackTime");
  TTreeReaderValue<Cluster> cluster(*dataReader, "clusters");
  TTreeReaderValue<std::vector<Digit>> digits(*dataReader, "digits");

  // setup the output
  TFile dataFileOut(outFile.c_str(), "recreate");
  TTree* dataTreeOut = dataReader->GetTree()->CloneTree(0);
  dataTreeOut->SetTitle("tree with input and output data");
  Cluster newCluster;
  dataTreeOut->Branch("newClusters", &newCluster);
  std::array<double, 6> fitParameters;
  dataTreeOut->Branch("fitParameters", &fitParameters);
  double pvalue;
  dataTreeOut->Branch("pvalue", &pvalue);
  double chi2;
  dataTreeOut->Branch("chi2", &chi2);

  std::vector<TH1*> preClusterInfo{};
  CreatePreClusterInfo(preClusterInfo);
  std::vector<TH1*> preClusterInfoSt[3] = {{}, {}, {}};
  CreatePreClusterInfo(preClusterInfoSt[0], "St1");
  CreatePreClusterInfo(preClusterInfoSt[1], "St2");
  CreatePreClusterInfo(preClusterInfoSt[2], "St345");

  int nClusters = dataReader->GetEntries(false);
  int iCluster = 0;
  int selected = 0;
  int fitted = 0;
  auto tStart = std::chrono::high_resolution_clock::now();
  while (dataReader->Next()) {
    if (++iCluster % 10000 == 0) {
      std::cout << "\rprocessing cluster " << iCluster << " / " << nClusters << "..." << std::flush;
    }

    // those 2 DE have lower HV for the run 529691
    if (run == 529691 && (cluster->getDEId() == 202 || cluster->getDEId() == 300)) {
      continue;
    }

    // cut on track angle at chamber
    if (std::abs(std::atan2(trackParam->py, -trackParam->pz)) / pi * 180. > 10.) {
      continue;
    }

    // cut on digit time
    std::vector<Digit> selectedDigits(*digits);
    selectedDigits.erase(
      std::remove_if(selectedDigits.begin(), selectedDigits.end(), [&trackTime](const auto& digit) {
        return std::abs(digit.getTime() + 1.5 - *trackTime) > 10.;
      }),
      selectedDigits.end());
    if (selectedDigits.empty()) {
      continue;
    }

    // reject mono-cathode preclusters after digit selection
    if (IsMonoCathode(selectedDigits)) {
      continue;
    }

    // reject composite preclusters
    if (IsComposite(selectedDigits, true)) {
      continue;
    }

    // cut on precluster charge asymmetry
    auto [chargeNB, chargeB] = GetCharge(selectedDigits, run < 300000);
    double chargeAsymm = (chargeNB - chargeB) / (chargeNB + chargeB);
    if (std::abs(chargeAsymm) > 0.5) {
      continue;
    }

    // check if precluster pass the fit selection
    if (!IsFittable(selectedDigits)) {
      continue;
    }

    ++selected;

    // init fit parameters (x, y, k3x, k3y, qBtot, qNBtot)
    auto local = GlobalToLocal(cluster->getDEId(), cluster->x, cluster->y, cluster->z, run < 300000);
    // correct charge total bending and charge total nonbending
    if (correctCharge) {
      auto [chargeFracNB, chargeFracB] = GetChargeFraction(selectedDigits, local.x(), local.y());
      chargeNB /= chargeFracNB;
      chargeB /= chargeFracB;
    }
    std::array<double, 6> param = {local.x(), local.y(), k3x, k3y, chargeB, chargeNB};

    // do the fit
    auto result = Fit(selectedDigits, param, fix, fitAsymm, errorMode, errorAlpha);
    if (result.Status() != 0) {
      continue;
    }
    chi2 = result.Chi2();
    pvalue = result.Prob();
    ++fitted;

    // fill characteristics of selected preclusters successfully fitted
    const auto [sizeX, sizeY] = GetSize(selectedDigits);
    FillPreClusterInfo(chargeNB, chargeB, sizeX, sizeY, preClusterInfo);
    int iSt = (cluster->getChamberId() < 4) ? cluster->getChamberId() / 2 : 2;
    FillPreClusterInfo(chargeNB, chargeB, sizeX, sizeY, preClusterInfoSt[iSt]);

    // fill output tree
    for (int i = 0; i < 5; ++i) {
      fitParameters[i] = result.Parameter(i);
    }
    fitParameters[5] = fitAsymm ? result.Parameter(5) : result.Parameter(4);
    newCluster = MakeCluster(cluster->uid, result.Parameter(0), result.Parameter(1));
    dataTreeOut->Fill();
  }

  auto tEnd = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> timer = tEnd - tStart;
  cout << "\r\033[Kprocessing completed. Duration = " << timer.count() << " s" << endl;
  cout << "selected clusters = " << selected << " / " << nClusters << endl;
  cout << "successfully fitted = " << fitted << " / " << selected << endl;

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
