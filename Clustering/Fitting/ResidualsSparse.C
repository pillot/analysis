#include <cmath>
#include <vector>

#include <fmt/format.h>

#include <TFile.h>
#include <TROOT.h>
#include <TTreeReaderArray.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <THnSparse.h>
#include <TH2D.h>
#include <TH1D.h>

#include "CommonUtils/ConfigurableParam.h"
#include "DataFormatsMCH/Cluster.h"
#include "DataFormatsMCH/Digit.h"
#include "MCHBase/TrackBlock.h"

#include "CCDBUtils.h"
#include "ClusterUtils.h"
#include "DataUtils.h"
#include "ResolutionUtils.h"
#include "PreClusterUtils.h"
#include "PlotsUtils.h"

using o2::mch::Cluster;
using o2::mch::Digit;
using o2::mch::TrackParamStruct;

static constexpr double pi = 3.14159265358979323846;
//_________________________________________________________________________________________________
void SetupMathieson(const double sqrtk3x_1, const double sqrtk3y_1, const double sqrtk3x_2345, const double sqrtk3y_2345)
{
  std::string K3X_1 = std::to_string(sqrtk3x_1);
  std::string K3Y_1 = std::to_string(sqrtk3y_1);
  std::string K3X_2345 = std::to_string(sqrtk3x_2345);
  std::string K3Y_2345 = std::to_string(sqrtk3y_2345);

  o2::conf::ConfigurableParam::setValue("MCHResponse.mathiesonSqrtKx3St1", K3X_1);
  o2::conf::ConfigurableParam::setValue("MCHResponse.mathiesonSqrtKy3St1", K3Y_1);

  o2::conf::ConfigurableParam::setValue("MCHResponse.mathiesonSqrtKx3St2345", K3X_2345);
  o2::conf::ConfigurableParam::setValue("MCHResponse.mathiesonSqrtKy3St2345", K3Y_2345);
}

//_________________________________________________________________________________________________
void ResidualsSparse(int run, const char* inFile = "clusters.root", const char* outFile = "residuals_sparse.root", int correctADCfit = 5)
{
  // require the MCH mapping to be loaded:
  // gSystem->Load("libO2MCHGeometryTransformer"),  gSystem->Load("libO2MCHMappingImpl4"), gSystem->Load("libO2MCHTracking")

  // load CCDB objects
  InitFromCCDB(run, true, true, false);

  if (correctADCfit != 0) {
    std::cout << "-- WARNING -- : Fit ADC selection is activated " << std::endl;
    std::cout << "SELECTION : " << std::abs(correctADCfit) << "< FitADC" << std::endl;
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
  fitParameters = std::make_unique<TTreeReaderArray<double>>(*dataReader, "fitParameters");

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

  // multi dimensional histogram : contains as dim : {p-value, residuals, ADC_fit, ADC_mes, ADC_cluster, nSamples, Asymm, Wire, Bending}
  THnSparseD* hPreClusterInfoMULTI[3];
  hPreClusterInfoMULTI[0] = CreatePreClusterInfoMULTI("St1");
  hPreClusterInfoMULTI[1] = CreatePreClusterInfoMULTI("St2");
  hPreClusterInfoMULTI[2] = CreatePreClusterInfoMULTI("St345");

  auto tStart = std::chrono::high_resolution_clock::now();
  std::cout << "looping over data ..." << std::endl;

  // loop precluster
  while (dataReader->Next()) {

    if (++iCluster % 10000 == 0) {
      std::cout << "\rprocessing cluster " << iCluster << " / " << nClusters << "..." << std::flush;
    }
    //___________________SELECTION__________________________
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
      std::remove_if(selectedDigits.begin(), selectedDigits.end(), [&trackTime](const auto& digit) { return std::abs(digit.getTime() + 1.5 - *trackTime) > 10.; }),
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

    // cut on ADCfit
    bool skip = false;
    // setup the mathieson
    auto sqrtK3x = sqrt(parameters[2]);
    auto sqrtK3y = sqrt(parameters[3]);
    SetupMathieson(sqrtK3x, sqrtK3y, sqrtK3x, sqrtK3y);
    for (auto digit : selectedDigits) {
      if (ADCFit(digit, parameters) < std::abs(correctADCfit)) {
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
    float dx_new = DistanceToClosestWire(cluster->getDEId(), parameters[0]); // use local X

    parameters.push_back(charge);      // parameters[6]
    parameters.push_back(chargeAsymm); // parameters[7]
    parameters.push_back(dx_new);      // parameters[8]
    parameters.push_back(**pvalue);    // parameters[9]

    for (auto digit : selectedDigits) {
      FillResolutionInfo(digit, parameters, hPreClusterInfoMULTI[iSt]);
    }

    // histograms that can't be in the THnSparse
    hAsymm[2 * iSt]->Fill(chargeAsymm);
    hAsymm[2 * iSt + 1]->Fill((parameters[5] - parameters[4]) / (parameters[4] + parameters[5]));
    auto [nPadsNB, nPadsB] = GetNPads(selectedDigits);
    h2chi2_ndf[iSt]->Fill((nPadsNB + nPadsB - 4), **chi2);
    hprob[iSt]->Fill(**pvalue);
  }

  dataFileIn->Close();

  std::cout << "Creating fit status plots ..." << std::endl;
  gStyle->SetOptStat(1);
  plot2D(h2chi2_ndf, "c_chi2_ndf", "chi2 vs ndf");
  plot1D(hprob, "c_prob", "p-value");
  plotSAME(hAsymm, "c_asymm", "Asymmetry");

  // output
  std::cout << "Saving plots ..." << std::endl;
  TFile fOut(outFile, "recreate");
  for (THnSparseD* const& h : hPreClusterInfoMULTI) {
    if (h)
      h->Write();
  }
  std::vector<std::string> canvasNames = {
    "c_chi2_ndf",
    "c_prob",
    "c_asymm"
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
  cout << "\r\033[Kprocessing completed. Duration = " << timer.count() << " s" << endl;
}
