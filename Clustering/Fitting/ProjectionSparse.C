#include <cmath>
#include <vector>
#include <string>
#include <optional>
#include <fmt/format.h>

#include <TCanvas.h>
#include <TFile.h>
#include <TROOT.h>
#include <THnSparse.h>
#include <TH2D.h>
#include <TH1D.h>

#include "ResolutionUtils.h"
#include "PlotsUtils.h"

//_________________________________________________________________________________________________
// require the MCH mapping to be loaded:
// gSystem->Load("libO2MCHGeometryTransformer"),  gSystem->Load("libO2MCHMappingImpl4"), gSystem->Load("libO2MCHTracking")

void ProjectionSparse(
  const std::string& inFile = "residuals_sparse.root",
  const std::string& outFile = "projection_sparse.root",
  const bool auto_bin = true,                   // adjust bin in the resolution extraction
  const std::pair<int, int>& projYX = { 1, 2 }, // default is Residual vs ADC from fitted pad
  const std::pair<std::optional<double>, std::optional<double>>& asym = { std::nullopt, std::nullopt },
  const std::pair<std::optional<double>, std::optional<double>>& chargetot = { std::nullopt, std::nullopt },
  const std::pair<std::optional<double>, std::optional<double>>& pvalue = { std::nullopt, std::nullopt },
  const std::pair<std::optional<double>, std::optional<double>>& nsamples = { std::nullopt, std::nullopt },
  const std::string& wire = "")
{

  auto warning = [](const std::string& label, const auto& range) {
    if (range.first && range.second) {
      std::cout << "-- WARNING -- : " << label << " selection is activated\n";
      std::cout << "SELECTION : " << *range.first << " < " << label << " < " << *range.second << "\n";
    }
  };

  warning("Asymmetry", asym);
  warning("Total Charge (ADC)", chargetot);
  warning("P-Value", pvalue);
  warning("NofSamples", nsamples);

  if (!wire.empty()) {
    std::cout << "-- WARNING -- : WIRE selection is activated\n";
    std::cout << "SELECTION : " << wire << std::endl;
  }

  std::cout << "loading data ..." << std::endl;
  TFile f(inFile.c_str(), "read");
  if (f.IsZombie()) {
    std::cerr << "Error: Cannot open file " << inFile << std::endl;
    return;
  }

  // 6 TList => odd : NBending, even : Bending (ordered by station number)
  std::vector<TH2D*> Residual2D;
  std::string sStation[3] = { "St1", "St2", "St345" };
  TList* ListResidual[6];

  for (int i = 0; i < 6; ++i) {
    ListResidual[i] = new TList();
    std::string cathode = (i % 2 == 0) ? "Bend" : "NBend";
    auto lName = fmt::format("Residual_{}_{}", sStation[i / 2], cathode);
    ListResidual[i]->SetName(lName.c_str());
  }

  auto tStart = std::chrono::high_resolution_clock::now();
  std::cout << "looping over the THnSparses ..." << std::endl;

  for (int i = 0; i < 6; i++) {

    auto sName = fmt::format("MultiResolutionPreCluster{}", sStation[i / 2]);

    // multi dimensional histogram whose axes are : {pvalue, residuals, ADC_fit, ADC_mes, ADC_cluster, nSamples, Asymm, Wire, Cathode}
    auto hSparse = dynamic_cast<THnSparse*>(f.Get(sName.c_str()));
    if (!hSparse) {
      std::cerr << "Warning: Could not find THnSparse " << sName << std::endl;
      continue;
    }

    // Apply any range cuts if needed :

    // P-Value
    if (pvalue.first && pvalue.second) {
      hSparse->GetAxis(0)->SetRangeUser(*pvalue.first, *pvalue.second);
    }

    // Asymmetry
    if (asym.first && asym.second) {
      hSparse->GetAxis(6)->SetRangeUser(*asym.first, *asym.second);
    }

    // Total Charge (ADC)
    if (chargetot.first && chargetot.second) {
      hSparse->GetAxis(4)->SetRangeUser(*chargetot.first, *chargetot.second);
    }

    // nSamples
    if (nsamples.first && nsamples.second) {
      hSparse->GetAxis(5)->SetRangeUser(*nsamples.first, *nsamples.second);
    }

    // Wire
    TAxis* axis7 = hSparse->GetAxis(7);
    if (wire == "top") {
      axis7->SetRange(1, 1);
    } else if (wire == "between") {
      axis7->SetRange(3, 3);
    } else if (wire == "crossover") {
      axis7->SetRange(2, 2);
    }

    // Cathode
    TAxis* axis8 = hSparse->GetAxis(8);
    if (i % 2 == 0) {
      axis8->SetRange(2, 2); // Bending bin
    } else {
      axis8->SetRange(1, 1); // NonBending bin
    }

    // 2D projection
    TH2D* h2D = dynamic_cast<TH2D*>(hSparse->Projection(projYX.first, projYX.second));
    std::string cathode = (i % 2 == 0) ? "Bend" : "NBend";
    auto hName = fmt::format("h2D_{}_{}", sStation[i / 2], cathode);
    auto hTitle = fmt::format("h2D_{}_{}", sStation[i / 2], cathode);
    h2D->SetName(hName.c_str());
    h2D->SetTitle(hTitle.c_str());

    Residual2D.push_back(h2D);
    int statistic = 2000; // minimum entries to construct a bin
    Resolution(ListResidual[i], h2D, statistic, auto_bin);
  }
  std::cout << "Saving plots ..." << std::endl;
  gStyle->SetOptStat(1);
  plot2D(Residual2D, "c_2Dresiduals", "residuals vs ADC", true);

  TFile fOut(outFile.c_str(), "recreate");
  if (fOut.IsZombie()) {
    std::cerr << "Error: Cannot open output file " << outFile << std::endl;
    return;
  }

  for (int i = 0; i < 6; ++i) {
    fOut.WriteTObject(ListResidual[i], ListResidual[i]->GetName());
  }

  if (TCanvas* c = dynamic_cast<TCanvas*>(gROOT->FindObject("c_2Dresiduals"))) {
    c->Write();
  } else {
    std::cerr << "Warning: Canvas not found." << std::endl;
  }

  fOut.Close();

  auto tEnd = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> timer = tEnd - tStart;
  cout << "\r\033[Kprocessing completed. Duration = " << timer.count() << " s" << endl;
}
