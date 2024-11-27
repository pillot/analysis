#include <cmath>
#include <utility>
#include <vector>

#include <gsl/span>

#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>

#include "DataFormatsMCH/Digit.h"
#include "MCHMappingInterface/Segmentation.h"

using o2::mch::Digit;

/*
 * This file contains utility functions to determine the
 * characteristics of a precluster and produce control plots.
 * To be used, they require the MCH mapping to be loaded:
 * gSystem->Load("libO2MCHMappingImpl4")
 */

//_________________________________________________________________________________________________
bool IsMonoCathode(const gsl::span<const Digit> digits)
{
  /// return true if all the digits are on the same cathode

  const auto& segmentation = o2::mch::mapping::segmentation(digits[0].getDetID());

  bool hasBendingPad(false);
  bool hasNonBendingPad(false);

  for (const auto& digit : digits) {
    if (segmentation.isBendingPad(digit.getPadID())) {
      if (hasNonBendingPad) {
        return false;
      }
      hasBendingPad = true;
    } else {
      if (hasBendingPad) {
        return false;
      }
      hasNonBendingPad = true;
    }
  }

  return true;
}

//_________________________________________________________________________________________________
std::pair<int, int> GetSize(const gsl::span<const Digit> digits)
{
  /// return the size of the precluster in pad unit
  /// the size in x (y) direction is given by the number of lines of pads
  /// on the non-bending (bending) plane in that direction
  /// note: the pad size is constant in the x (y) direction on the non-bending (bending) plane

  const auto& segmentation = o2::mch::mapping::segmentation(digits[0].getDetID());

  double dimX[2] = {1.e6, -1.e6};
  double padSizeX = -1.;
  double dimY[2] = {1.e6, -1.e6};
  double padSizeY = -1.;

  for (const auto& digit : digits) {
    if (segmentation.isBendingPad(digit.getPadID())) {
      double padY = segmentation.padPositionY(digit.getPadID());
      dimY[0] = std::min(dimY[0], padY);
      dimY[1] = std::max(dimY[1], padY);
      padSizeY = segmentation.padSizeY(digit.getPadID());
    } else {
      double padX = segmentation.padPositionX(digit.getPadID());
      dimX[0] = std::min(dimX[0], padX);
      dimX[1] = std::max(dimX[1], padX);
      padSizeX = segmentation.padSizeX(digit.getPadID());
    }
  }

  int sizeX = (padSizeX > 0.) ? std::lround((dimX[1] - dimX[0]) / padSizeX) + 1 : 0;
  int sizeY = (padSizeY > 0.) ? std::lround((dimY[1] - dimY[0]) / padSizeY) + 1 : 0;

  return std::make_pair(sizeX, sizeY);
}

//_________________________________________________________________________________________________
std::pair<double, double> GetCharge(const gsl::span<const Digit> digits, bool run2 = false)
{
  /// return the total charge of the digits on both cathodes

  // function to reinterpret digit ADC as calibrated charge in run2
  static constexpr auto adcToCharge = [](uint32_t adc) {
    float charge(0.);
    std::memcpy(&charge, &adc, sizeof(adc));
    return static_cast<double>(charge);
  };

  const auto& segmentation = o2::mch::mapping::segmentation(digits[0].getDetID());

  std::pair<double, double> charge{0., 0.};

  for (const auto& digit : digits) {
    if (segmentation.isBendingPad(digit.getPadID())) {
      charge.second += run2 ? adcToCharge(digit.getADC()) : digit.getADC();
    } else {
      charge.first += run2 ? adcToCharge(digit.getADC()) : digit.getADC();
    }
  }

  return charge;
}

//_________________________________________________________________________________________________
std::pair<double, double> GetCOG(const gsl::span<const Digit> digits)
{
  /// return a the center of gravity of the digits
  /// the weight of each digit is given by its ADC charge
  /// for bi-cathode clusters, x (y) position is given by digits in the non-bending (bending) plane
  /// note: the pad size is constant in the x (y) direction on the non-bending (bending) plane

  const auto& segmentation = o2::mch::mapping::segmentation(digits[0].getDetID());

  double x[2] = {0., 0.};
  double y[2] = {0., 0.};
  double q[2] = {0., 0.};

  for (const auto& digit : digits) {
    int i = segmentation.isBendingPad(digit.getPadID()) ? 1 : 0;
    x[i] += digit.getADC() * segmentation.padPositionX(digit.getPadID());
    y[i] += digit.getADC() * segmentation.padPositionY(digit.getPadID());
    q[i] += digit.getADC();
  }

  double xCl = (q[0] > 0.) ? x[0] / q[0] : x[1] / q[1];
  double yCl = (q[1] > 0.) ? y[1] / q[1] : y[0] / q[0];

  return std::make_pair(xCl, yCl);
}

//_________________________________________________________________________________________________
std::pair<double, double> GetCOG2(const gsl::span<const Digit> digits)
{
  /// return the center of gravity of the digits
  /// the weight of each digit is given by its ADC charge / (its size / sqrt(12))^2
  /// note: the pad size is constant in the x (y) direction on the non-bending (bending) plane
  /// the purpose of weighting by the pad size is to combine the digits from both planes

  const auto& segmentation = o2::mch::mapping::segmentation(digits[0].getDetID());

  double x(0.);
  double y(0.);
  double wx(0.);
  double wy(0.);

  for (const auto& digit : digits) {
    double wxi = digit.getADC() / std::pow(segmentation.padSizeX(digit.getPadID()), 2.);
    double wyi = digit.getADC() / std::pow(segmentation.padSizeY(digit.getPadID()), 2.);
    x += wxi * segmentation.padPositionX(digit.getPadID());
    y += wyi * segmentation.padPositionY(digit.getPadID());
    wx += wxi;
    wy += wyi;
  }

  return std::make_pair(x / wx, y / wy);
}

//_________________________________________________________________________________________________
void CreatePreClusterInfo(std::vector<TH1*>& histos, const char* extension = "")
{
  /// create histograms of precluster characteristics

  histos.emplace_back(new TH1F(Form("hCharge%s", extension),
                               "cluster charge;charge (ADC)", 5000, -0.25, 49999.75));
  histos.emplace_back(new TH1F(Form("hChargeAsymm%s", extension),
                               "cluster charge asymmetry;asymmetry", 201, -1.005, 1.005));
  histos.emplace_back(new TH2F(Form("hChargeAsymm2D%s", extension),
                               "cluster charge asymmetry vs charge;charge (ADC);asymmetry",
                               200, -0.25, 19999.75, 201, -1.005, 1.005));
  histos.emplace_back(new TH1F(Form("hDimX%s", extension),
                               "cluster size X;size X (#pads)", 21, -0.5, 20.5));
  histos.emplace_back(new TH1F(Form("hDimY%s", extension),
                               "cluster size Y;size Y (#pads)", 31, -0.5, 30.5));
  histos.emplace_back(new TH2F(Form("hDimXY%s", extension),
                               "cluster size Y vs X;size X (#pads);size Y (#pads)",
                               21, -0.5, 20.5, 31, -0.5, 30.5));

  for (auto h : histos) {
    h->SetDirectory(0);
  }
}

//_________________________________________________________________________________________________
void FillPreClusterInfo(double charge, double chargeAsymm, int sizeX, int sizeY, std::vector<TH1*>& histos)
{
  /// fill histograms of precluster characteristics

  histos[0]->Fill(charge);
  histos[1]->Fill(chargeAsymm);
  histos[2]->Fill(charge, chargeAsymm);
  histos[3]->Fill(sizeX);
  histos[4]->Fill(sizeY);
  histos[5]->Fill(sizeX, sizeY);
}

//_________________________________________________________________________________________________
TCanvas* DrawPreClusterInfo(std::vector<TH1*>& histos, const char* extension = "")
{
  /// draw histograms of precluster characteristics

  static int logy[6] = {1, 1, 0, 1, 1, 0};
  static int logz[6] = {0, 0, 1, 0, 0, 1};
  static const char* opt[6] = {"", "", "colz", "", "", "colz"};

  TCanvas* c = new TCanvas(Form("preClusterInfo%s", extension),
                           Form("precluster characteristics %s", extension), 10, 10, 900, 600);
  c->Divide(3, 2);
  int i(0);
  for (const auto& h : histos) {
    c->cd(i + 1);
    gPad->SetLogy(logy[i]);
    gPad->SetLogz(logz[i]);
    h->Draw(opt[i]);
    ++i;
  }

  return c;
}
