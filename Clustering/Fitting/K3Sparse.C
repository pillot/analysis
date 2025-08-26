#include <cmath>
#include <numeric>
#include <vector>

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <THnSparse.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <fmt/format.h>
#include "TROOT.h"
#include "TTreeReaderArray.h"

#include "CCDBUtils.h"
#include "ClusterUtils.h"
#include "CommonUtils/ConfigurableParam.h"
#include "DataFormatsMCH/Cluster.h"
#include "DataFormatsMCH/Digit.h"
#include "DataUtils.h"
#include "MCHBase/TrackBlock.h"
#include "PlotsUtils.h"
#include "PreClusterUtils.h"
#include "ResolutionUtils.h"

using o2::mch::Cluster;
using o2::mch::Digit;
using o2::mch::TrackParamStruct;

static constexpr double pi = 3.14159265358979323846;

//_________________________________________________________________________________________________
void K3Sparse(int run, const char* inFile = "clusters.root", const char* outFile = "residuals_sparse.root", int correctADCfit = 5, bool uniqueResponse = true)
{
  /// require the MCH mapping to be loaded:
  /// gSystem->Load("libO2MCHGeometryTransformer"),
  /// gSystem->Load("libO2MCHMappingImpl4"), gSystem->Load("libO2MCHTracking")

  /// load CCDB objects
  InitFromCCDB(run, true, true, false);

  if (correctADCfit != 0) {
    std::cout << "-- WARNING -- : Fit ADC selection is activated " << std::endl;
    std::cout << "SELECTION : " << std::abs(correctADCfit) << "< FitADC"
              << std::endl;
  }
  // load histograms
  LoadHist();

  //________________________________________________________________________________________________
  // load input data and loop
  //________________________________________________________________________________________________
  std::cout << "loading data ..." << std::endl;

  auto [dataFileIn, dataReader] = LoadData(inFile, "data");
  TTreeReaderValue<TrackParamStruct> trackParam(*dataReader, "trackParameters");
  TTreeReaderValue<int> trackTime(*dataReader, "trackTime");
  TTreeReaderValue<Cluster> cluster(*dataReader, "clusters");
  TTreeReaderValue<std::vector<Digit>> digits(*dataReader, "digits");
  std::unique_ptr<TTreeReaderArray<double>> fitParameters{};

  if (!dataReader->GetTree()->FindBranch("fitParameters")) {
    LOGP(error, "unable to load branch \"fitParameters\" from {}", inFile);
    exit(-1);
  }
  fitParameters =
    std::make_unique<TTreeReaderArray<double>>(*dataReader, "fitParameters");

  std::unique_ptr<TTreeReaderValue<double>> pvalue{};
  if (!dataReader->GetTree()->FindBranch("pvalue")) {
    LOGP(error, "unable to load branch \"pvalue\" from {}", inFile);
    exit(-1);
  }
  pvalue = std::make_unique<TTreeReaderValue<double>>(*dataReader, "pvalue");

  std::unique_ptr<TTreeReaderValue<double>> chi2{};
  if (!dataReader->GetTree()->FindBranch("chi2")) {
    LOGP(error, "unable to load branch \"chi2\" from {}", inFile);
    exit(-1);
  }
  chi2 = std::make_unique<TTreeReaderValue<double>>(*dataReader, "chi2");

  int nClusters = dataReader->GetEntries(false);
  int iCluster(0);

  // multi dimensional histogram : contains as dim : {p-value, k3x, k3y,
  // ADC_cluster, p , phi, Asymm, Wire, Bending}
  THnSparseD* hPreClusterInfoMULTIK3[3];
  hPreClusterInfoMULTIK3[0] = CreatePreClusterInfoMULTIK3("St1");
  hPreClusterInfoMULTIK3[1] = CreatePreClusterInfoMULTIK3("St2");
  hPreClusterInfoMULTIK3[2] = CreatePreClusterInfoMULTIK3("St345");

  auto tStart = std::chrono::high_resolution_clock::now();
  std::cout << "looping over data ..." << std::endl;

  std::vector<double> delta_k3x, delta_k3y;
  // loop precluster
  while (dataReader->Next()) {
    if (++iCluster % 10000 == 0) {
      std::cout << "\rprocessing cluster " << iCluster << " / " << nClusters
                << "..." << std::flush;
    }
    //___________________SELECTION__________________________
    // those 2 DE have lower HV for the run 529691
    if (run == 529691 &&
        (cluster->getDEId() == 202 || cluster->getDEId() == 300)) {
      continue;
    }

    // cut on track angle at chamber
    double Track_angle =
      std::abs(std::atan2(trackParam->py, -trackParam->pz)) / pi * 180.;
    double Track_momentum =
      sqrt(trackParam->px * trackParam->px + trackParam->py * trackParam->py +
           trackParam->pz * trackParam->pz);
    if (Track_angle > 10.) {
      continue;
    }
    // cut on digit time
    std::vector<Digit> selectedDigits(*digits);
    selectedDigits.erase(
      std::remove_if(selectedDigits.begin(), selectedDigits.end(),
                     [&trackTime](const auto& digit) {
                       return std::abs(digit.getTime() + 1.5 - *trackTime) >
                              10.;
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

    // check if precluster pass the fit selection if needed (N° pads, size, ...)
    if (!IsFittable(selectedDigits)) {
      continue;
    }

    // cut on precluster charge asymmetry
    auto [chargeNB, chargeB] = GetCharge(selectedDigits, run < 300000);
    double chargeAsymm = (chargeNB - chargeB) / (chargeNB + chargeB);
    double charge = sqrt(chargeNB * chargeB);

    if (std::abs(chargeAsymm) > 0.5) {
      continue;
    }

    std::vector<double> parameters;
    for (int i = 0; i < 6; ++i) {
      parameters.push_back((*fitParameters)[i]);
    }

    bool skip = false;
    // determine if we have a logical issue (i.e. using an unique response while requiring the response to vary)
    delta_k3x.push_back(parameters[2]); // vector that keep growing in size
    delta_k3y.push_back(parameters[3]);
    double avg_k3x = std::accumulate(delta_k3x.begin(), delta_k3x.end(), 0.0) / delta_k3x.size(); // average calculation
    double avg_k3y = std::accumulate(delta_k3y.begin(), delta_k3y.end(), 0.0) / delta_k3y.size();
    bool err_k3x = std::all_of(delta_k3x.begin(), delta_k3x.end(), [avg_k3x](double v) { return std::abs(v - avg_k3x) < 1e-5; }); // difference between all elements of the vector and the mean value calculated before
    bool err_k3y = std::all_of(delta_k3y.begin(), delta_k3y.end(), [avg_k3y](double v) { return std::abs(v - avg_k3y) < 1e-5; });

    if (uniqueResponse && (!err_k3x || !err_k3y)) {
      std::cerr << "Unique response is asked but k3 is not constant..." << std::endl;
      exit(-1);
    }
    ////cut on ADCfit
    for (auto digit : selectedDigits) {
      if (ADCFit(digit, parameters, uniqueResponse) < std::abs(correctADCfit)) {
        skip = true;
      }
    }
    if (skip) {
      continue;
    } // skip cluster who have at least one ADCfit < cut
    if (!skip && (correctADCfit < 0)) {
      continue;
    } // skip cluster who have all ADCfit > cut

    int iSt = (cluster->getChamberId() < 4) ? cluster->getChamberId() / 2 : 2;
    float dx_new = DistanceToClosestWire(cluster->getDEId(),
                                         parameters[0]); // use local X

    parameters.push_back(charge);         // parameters[6]
    parameters.push_back(chargeAsymm);    // parameters[7]
    parameters.push_back(dx_new);         // parameters[8]
    parameters.push_back(**pvalue);       // parameters[9]
    parameters.push_back(Track_angle);    // parameters[10]
    parameters.push_back(Track_momentum); // parameters[11]

    FillK3Info(parameters, hPreClusterInfoMULTIK3[iSt]);

    // histograms that can't be in the THnSparse
    auto [nPadsNB, nPadsB] = GetNPads(selectedDigits);
    h2chi2_ndf[iSt]->Fill((nPadsNB + nPadsB - 4), **chi2);
    hprob[iSt]->Fill(**pvalue);
    hk3x[iSt]->Fill(parameters[2]);
    hk3y[iSt]->Fill(parameters[3]);
  }

  dataFileIn->Close();

  std::cout << "Creating fit status plots ..." << std::endl;
  gStyle->SetOptStat(1);
  plot2D(h2chi2_ndf, "c_chi2_ndf", "chi2 vs ndf");
  plot1D(hprob, "c_prob", "p-value");
  plot1D(hk3x, "c_k3x", "K_{3X}");
  plot1D(hk3y, "c_k3y", "K_{3Y}");

  // output
  std::cout << "Saving plots ..." << std::endl;
  TFile fOut(outFile, "recreate");
  for (THnSparseD* const& h : hPreClusterInfoMULTIK3) {
    if (h)
      h->Write();
  }
  std::vector<std::string> canvasNames = {
    "c_chi2_ndf",
    "c_prob",
    "c_k3x",
    "c_k3y",
  };
  for (const auto& name : canvasNames) {
    if (TCanvas* c = (TCanvas*)gROOT->FindObject(name.c_str())) {
      c->Write();
    } else {
      std::cerr << "Warning: Canvas " << name << " not found." << std::endl;
    }
  }

  fOut.Close();

  auto tEnd = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> timer = tEnd - tStart;
  cout << "\r\033[Kprocessing completed. Duration = " << timer.count() << " s"
       << endl;
}
