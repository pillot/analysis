#include <cmath>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <TCanvas.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TGraphAsymmErrors.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TLegend.h>
#include <TMath.h>
#include <TStyle.h>

#include <fmt/format.h>

#include "Framework/Logger.h"

TH3* GetAsymm3D(TFile& f, std::string hName);
std::vector<TGraphAsymmErrors*> GetAsymmDispersionVsCharge(TH2* hAsymm2D, const char* extension, int minEntriesPerBin);

//_________________________________________________________________________________________________
void FitPreClustersAsymm(std::string inFileName = "displays.root", std::string outFileName = "asymm.root",
                         int minEntriesPerBin = 1000)
{
  /// fit the cluster charge asymmetry versus charge for different distances to closest wire

  gStyle->SetOptStat(1);

  TFile inFile(inFileName.c_str(), "read");

  static const char* sSt[] = {"St1", "St2", "St345", ""};

  std::vector<TH3*> hAsymm3D{};
  for (auto st : sSt) {
    hAsymm3D.push_back(GetAsymm3D(inFile, fmt::format("hChargeAsymm3D{}", st).c_str()));
  }

  std::vector<TCanvas*> cAsymm{};
  std::vector<std::pair<double, double>> yRange{};
  cAsymm.emplace_back(new TCanvas("cAsymmMean", "asymmetry mean vs charge", 10, 10, 1200, 600));
  yRange.emplace_back(-0.1, 0.1);
  cAsymm.emplace_back(new TCanvas("cAsymmRMS", "asymmetry RMS vs charge", 10, 10, 1200, 600));
  yRange.emplace_back(0., 0.4);
  cAsymm.emplace_back(new TCanvas("cAsymmMeanG", "asymmetry gaus mean vs charge", 10, 10, 1200, 600));
  yRange.emplace_back(-0.1, 0.1);
  cAsymm.emplace_back(new TCanvas("cAsymmSigmaG", "asymmetry gaus sigma vs charge", 10, 10, 1200, 600));
  yRange.emplace_back(0., 0.5);
  cAsymm.emplace_back(new TCanvas("cAsymmMeanDG", "asymmetry double gaus mean vs charge", 10, 10, 1200, 600));
  yRange.emplace_back(-0.1, 0.1);
  cAsymm.emplace_back(new TCanvas("cAsymmDeltaDG", "asymmetry double gaus delta vs charge", 10, 10, 1200, 600));
  yRange.emplace_back(0., 0.4);
  cAsymm.emplace_back(new TCanvas("cAsymmSigmaDG", "asymmetry double gaus sigma vs charge", 10, 10, 1200, 600));
  yRange.emplace_back(0., 0.4);
  cAsymm.emplace_back(new TCanvas("cAsymmFracTG", "asymmetry triple gaus repartition vs charge", 10, 10, 1200, 600));
  yRange.emplace_back(-0.5, 1.5);
  cAsymm.emplace_back(new TCanvas("cAsymmMeanTG", "asymmetry triple gaus mean vs charge", 10, 10, 1200, 600));
  yRange.emplace_back(-0.1, 0.1);
  cAsymm.emplace_back(new TCanvas("cAsymmDeltaTG", "asymmetry triple gaus delta vs charge", 10, 10, 1200, 600));
  yRange.emplace_back(0., 0.4);
  cAsymm.emplace_back(new TCanvas("cAsymmSigmaTG", "asymmetry triple gaus sigma vs charge", 10, 10, 1200, 600));
  yRange.emplace_back(0., 0.4);
  for (auto c : cAsymm) {
    c->Divide(2, 2);
  }

  auto draw = [](auto g1, auto g2, auto g3, auto yRange) {
    g1->Draw("ap");
    g1->GetXaxis()->SetRangeUser(0., 20000.);
    g1->GetYaxis()->SetRangeUser(yRange.first, yRange.second);
    g2->SetMarkerColor(2);
    g2->SetLineColor(2);
    g2->Draw("p");
    g3->SetMarkerColor(4);
    g3->SetLineColor(4);
    g3->Draw("p");
  };

  for (size_t iSt = 0; iSt < hAsymm3D.size(); ++iSt) {

    std::unique_ptr<TH2> hAsymm2D(static_cast<TH2*>(hAsymm3D[iSt]->Project3D("zy")));
    auto gAsymm = GetAsymmDispersionVsCharge(hAsymm2D.get(), sSt[iSt], minEntriesPerBin);

    hAsymm3D[iSt]->GetXaxis()->SetRange(12, 14);
    std::unique_ptr<TH2> hAsymm2DClose(static_cast<TH2*>(hAsymm3D[iSt]->Project3D("close_zy")));
    auto extension = fmt::format("{}Close", sSt[iSt]);
    auto gAsymmClose = GetAsymmDispersionVsCharge(hAsymm2DClose.get(), extension.c_str(), minEntriesPerBin);

    hAsymm3D[iSt]->GetXaxis()->SetRange(1, 4);
    std::unique_ptr<TH2> hAsymm2DFar(static_cast<TH2*>(hAsymm3D[iSt]->Project3D("far1_zy")));
    hAsymm3D[iSt]->GetXaxis()->SetRange(22, 25);
    std::unique_ptr<TH2> hAsymm2DFar2(static_cast<TH2*>(hAsymm3D[iSt]->Project3D("far2_zy")));
    hAsymm2DFar->Add(hAsymm2DFar2.get());
    extension = fmt::format("{}Far", sSt[iSt]);
    auto gAsymmFar = GetAsymmDispersionVsCharge(hAsymm2DFar.get(), extension.c_str(), minEntriesPerBin);

    TLegend* l = new TLegend(0.15, 0.75, 0.35, 0.9);
    l->SetFillStyle(0);
    l->SetBorderSize(0);
    l->AddEntry(gAsymm[0], "any distance", "l");
    l->AddEntry(gAsymmClose[0], "|dx| < 0.015 cm", "l");
    l->AddEntry(gAsymmFar[0], "|dx| > 0.085 cm", "l");

    for (size_t i = 0; i < cAsymm.size(); ++i) {
      cAsymm[i]->cd(iSt + 1);
      draw(gAsymm[i], gAsymmClose[i], gAsymmFar[i], yRange[i]);
      l->Draw("same");
    }
  }

  inFile.Close();

  TFile outFile(outFileName.c_str(), "recreate");
  for (auto c : cAsymm) {
    c->Write();
  }
  outFile.Close();
}

//_________________________________________________________________________________________________
TH3* GetAsymm3D(TFile& f, std::string hName)
{
  /// return the requested 3D histogram

  auto h = f.FindObjectAny(hName.c_str());
  if (h == nullptr) {
    LOGP(error, "histo {} not found", hName);
    exit(-1);
  }

  return static_cast<TH3*>(h);
}

//_________________________________________________________________________________________________
double DoubleGaus(double* x, double* par)
{
  /// double gaussian with same norm and width, placed at -delta and delta from mean
  /// par[0] = normalization
  /// par[1] = mean
  /// par[2] = delta
  /// par[3] = sigma

  return par[0] * 0.5 *
         (TMath::Gaus(x[0], par[1] - par[2], par[3], true) +
          TMath::Gaus(x[0], par[1] + par[2], par[3], true));
}

//_________________________________________________________________________________________________
double TripleGaus(double* x, double* par)
{
  /// triple gaussian with same width, one at mean and two shifted by ± delta with a different norm
  /// par[0] = global normalization
  /// par[1] = fraction in the two gaus at mean ± delta
  /// par[2] = mean
  /// par[3] = delta
  /// par[4] = sigma

  auto frac = std::abs(par[1]);
  return par[0] * (frac * 0.5 *
                     (TMath::Gaus(x[0], par[2] - par[3], par[4], true) +
                      TMath::Gaus(x[0], par[2] + par[3], par[4], true)) +
                   (1. - frac) * TMath::Gaus(x[0], par[2], par[4], true));
}

//_________________________________________________________________________________________________
std::vector<TGraphAsymmErrors*> GetAsymmDispersionVsCharge(TH2* hAsymm2D, const char* extension, int minEntriesPerBin)
{
  /// extract the charge asymmetry dispersion versus charge from the 2D histogram

  std::unique_ptr<TH1> hCharge(hAsymm2D->ProjectionX());
  std::vector<TGraphAsymmErrors*> gAsymm{};

  auto addGraph = [&gAsymm](std::string name, std::string title) -> TGraphAsymmErrors* {
    auto g = gAsymm.emplace_back(new TGraphAsymmErrors());
    g->SetNameTitle(name.c_str(), title.c_str());
    return g;
  };

  auto gAsymmMeanvsCharge = addGraph(fmt::format("gAsymmMeanvsCharge{}", extension).c_str(),
                                     fmt::format("{} asymmetry mean vs charge", extension).c_str());
  auto gAsymmRMSvsCharge = addGraph(fmt::format("gAsymmRMSvsCharge{}", extension).c_str(),
                                    fmt::format("{} asymmetry RMS vs charge", extension).c_str());

  auto gAsymmMeanGvsCharge = addGraph(fmt::format("gAsymmMeanGvsCharge{}", extension).c_str(),
                                      fmt::format("{} asymmetry gaus mean vs charge", extension).c_str());
  auto gAsymmSigmaGvsCharge = addGraph(fmt::format("gAsymmSigmaGvsCharge{}", extension).c_str(),
                                       fmt::format("{} asymmetry gaus sigma vs charge", extension).c_str());

  auto gAsymmMeanDGvsCharge = addGraph(fmt::format("gAsymmMeanDGvsCharge{}", extension).c_str(),
                                       fmt::format("{} asymmetry double gaus mean vs charge", extension).c_str());
  auto gAsymmDeltaDGvsCharge = addGraph(fmt::format("gAsymmDeltaDGvsCharge{}", extension).c_str(),
                                        fmt::format("{} asymmetry double gaus delta vs charge", extension).c_str());
  auto gAsymmSigmaDGvsCharge = addGraph(fmt::format("gAsymmSigmaDGvsCharge{}", extension).c_str(),
                                        fmt::format("{} asymmetry double gaus sigma vs charge", extension).c_str());

  auto gAsymmFracTGvsCharge = addGraph(fmt::format("gAsymmFracTGvsCharge{}", extension).c_str(),
                                       fmt::format("{} asymmetry triple gaus repartition vs charge", extension).c_str());
  auto gAsymmMeanTGvsCharge = addGraph(fmt::format("gAsymmMeanTGvsCharge{}", extension).c_str(),
                                       fmt::format("{} asymmetry triple gaus mean vs charge", extension).c_str());
  auto gAsymmDeltaTGvsCharge = addGraph(fmt::format("gAsymmDeltaTGvsCharge{}", extension).c_str(),
                                        fmt::format("{} asymmetry triple gaus delta vs charge", extension).c_str());
  auto gAsymmSigmaTGvsCharge = addGraph(fmt::format("gAsymmSigmaTGvsCharge{}", extension).c_str(),
                                        fmt::format("{} asymmetry triple gaus sigma vs charge", extension).c_str());

  auto addDispersionValues = [&](int firstBin, int lastBin) {
    hCharge->GetXaxis()->SetRange(firstBin, lastBin);
    auto chargeMean = hCharge->GetMean();
    auto chargeMeanErrLow = chargeMean - hCharge->GetBinLowEdge(firstBin);
    auto chargeMeanErrHigh = hCharge->GetBinLowEdge(lastBin + 1) - chargeMean;

    auto addPoint = [&](TGraphAsymmErrors* g, double y, double ey) {
      g->AddPoint(chargeMean, y);
      g->SetPointError(g->GetN() - 1, chargeMeanErrLow, chargeMeanErrHigh, ey, ey);
    };

    std::unique_ptr<TH1> hAsymm(hAsymm2D->ProjectionY("tmp", firstBin, lastBin));
    addPoint(gAsymmMeanvsCharge, hAsymm->GetMean(), hAsymm->GetMeanError());
    addPoint(gAsymmRMSvsCharge, hAsymm->GetStdDev(), hAsymm->GetStdDevError());

    static TF1* fGaus = new TF1("Gaus", "gausn", -0.5, 0.5);
    auto norm = hAsymm->GetEntries() * hAsymm->GetBinWidth(1);
    fGaus->SetParameters(norm, 0., 0.05);
    auto fitResultG = hAsymm->Fit(fGaus, "RSNQ");
    addPoint(gAsymmMeanGvsCharge, fitResultG->Parameter(1), fitResultG->ParError(1));
    addPoint(gAsymmSigmaGvsCharge, std::abs(fitResultG->Parameter(2)), fitResultG->ParError(2));

    static TF1* fDoubleGaus = new TF1("DoubleGaus", DoubleGaus, -0.5, 0.5, 4);
    fDoubleGaus->SetParameters(norm, 0., 0., 0.1);
    auto fitResultDG = hAsymm->Fit(fDoubleGaus, "RSNQ");
    addPoint(gAsymmMeanDGvsCharge, fitResultDG->Parameter(1), fitResultDG->ParError(1));
    addPoint(gAsymmDeltaDGvsCharge, std::abs(fitResultDG->Parameter(2)), fitResultDG->ParError(2));
    addPoint(gAsymmSigmaDGvsCharge, std::abs(fitResultDG->Parameter(3)), fitResultDG->ParError(3));

    static TF1* fTripleGaus = new TF1("TripleGaus", TripleGaus, -0.5, 0.5, 5);
    fTripleGaus->SetParameters(norm, 0.3, 0., 0.15, 0.05);
    fTripleGaus->SetParLimits(1, -1., 1.);
    auto fitResultTG = hAsymm->Fit(fTripleGaus, "RSNQ");
    addPoint(gAsymmFracTGvsCharge, std::abs(fitResultTG->Parameter(1)), fitResultTG->ParError(1));
    addPoint(gAsymmMeanTGvsCharge, fitResultTG->Parameter(2), fitResultTG->ParError(2));
    addPoint(gAsymmDeltaTGvsCharge, std::abs(fitResultTG->Parameter(3)), fitResultTG->ParError(3));
    addPoint(gAsymmSigmaTGvsCharge, std::abs(fitResultTG->Parameter(4)), fitResultTG->ParError(4));
  };

  int lastBin = 0;
  auto remainingEntries = static_cast<int>(hCharge->Integral());
  while (remainingEntries > 2 * minEntriesPerBin) {

    int firstBin = lastBin + 1;
    int nEntries = 0;
    do {
      nEntries += hCharge->GetBinContent(++lastBin);
    } while (nEntries < minEntriesPerBin);
    remainingEntries -= nEntries;

    addDispersionValues(firstBin, lastBin);
  }

  if (remainingEntries > 0 && lastBin < hCharge->GetNbinsX()) {
    addDispersionValues(lastBin + 1, hCharge->GetNbinsX());
  }

  return gAsymm;
}
