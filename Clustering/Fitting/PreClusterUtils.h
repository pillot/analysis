#include <cmath>
#include <memory>
#include <tuple>
#include <utility>
#include <vector>

#include <fmt/format.h>

#include <gsl/span>

#include <TCanvas.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TGraphAsymmErrors.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TMath.h>

#include "DataFormatsMCH/Digit.h"
#include "MCHMappingInterface/Segmentation.h"

#include "DigitUtils.h"

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
  histos.emplace_back(new TH1F(Form("hChargeB%s", extension),
                               "cluster charge bending;charge (ADC)", 5000, -0.25, 49999.75));
  histos.emplace_back(new TH1F(Form("hChargeNB%s", extension),
                               "cluster charge non bending;charge (ADC)", 5000, -0.25, 49999.75));
  histos.emplace_back(new TH1F(Form("hChargeAsymm%s", extension),
                               "cluster charge asymmetry;(NB-B)/(NB+B)", 201, -1.005, 1.005));
  histos.emplace_back(new TH2F(Form("hChargeAsymm2D%s", extension),
                               "cluster charge asymmetry vs charge;(NB+B)/2 (ADC);(NB-B)/(NB+B)",
                               200, -0.25, 19999.75, 201, -1.005, 1.005));
  histos.back()->GetYaxis()->SetTitleOffset(1.5);
  histos.emplace_back(new TH2F(Form("hChargeCorr%s", extension),
                               "cluster charge bending vs non-bending;charge NB (ADC);charge B (ADC)",
                               200, -0.25, 19999.75, 200, -0.25, 19999.75));
  histos.back()->GetYaxis()->SetTitleOffset(1.5);
  histos.emplace_back(new TH1F(Form("hChargeAsymm2%s", extension),
                               "cluster charge asymmetry2;1/2 * ln(NB/B)", 201, -1.005, 1.005));
  histos.emplace_back(new TH2F(Form("hChargeAsymm22D%s", extension),
                               "cluster charge asymmetry2 vs charge2;#sqrt{NB*B} (ADC);1/2 * ln(NB/B)",
                               200, -0.25, 19999.75, 201, -1.005, 1.005));
  histos.back()->GetYaxis()->SetTitleOffset(1.5);
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
void FillPreClusterInfo(double chargeNB, double chargeB, int sizeX, int sizeY, std::vector<TH1*>& histos)
{
  /// fill histograms of precluster characteristics

  double charge = 0.5 * (chargeNB + chargeB);
  double charge2 = std::sqrt(chargeB * chargeNB);
  double chargeAsymm = (chargeNB - chargeB) / (chargeNB + chargeB);
  double chargeRatio = 0.5 * std::log(chargeNB / chargeB);

  histos[0]->Fill(charge);
  histos[1]->Fill(chargeB);
  histos[2]->Fill(chargeNB);
  histos[3]->Fill(chargeAsymm);
  histos[4]->Fill(charge, chargeAsymm);
  histos[5]->Fill(chargeNB, chargeB);
  histos[6]->Fill(chargeRatio);
  histos[7]->Fill(charge2, chargeRatio);
  histos[8]->Fill(sizeX);
  histos[9]->Fill(sizeY);
  histos[10]->Fill(sizeX, sizeY);
}

//_________________________________________________________________________________________________
TCanvas* DrawPreClusterInfo(std::vector<TH1*>& histos, const char* extension = "")
{
  /// draw histograms of precluster characteristics

  static int logy[] = {1, 1, 0, 0, 1, 0, 1, 1, 0};
  static int logz[] = {0, 0, 1, 1, 0, 1, 0, 0, 1};
  static const char* opt[] = {"", "", "colz", "colz", "", "colz", "", "", "colz"};

  TCanvas* c = new TCanvas(Form("preClusterInfo%s", extension),
                           Form("precluster characteristics %s", extension), 10, 10, 900, 900);
  c->Divide(3, 3);

  TLegend* l = new TLegend(0.45, 0.7, 0.85, 0.9);
  l->SetFillStyle(0);
  l->SetBorderSize(0);

  c->cd(1);
  gPad->SetLogy();
  histos[0]->SetLineColor(1);
  histos[0]->Draw();
  l->AddEntry(histos[0], "(NB+B)/2", "l");
  histos[1]->SetStats(false);
  histos[1]->SetLineColor(2);
  histos[1]->Draw("sames");
  l->AddEntry(histos[1], "Bending", "l");
  histos[2]->SetStats(false);
  histos[2]->SetLineColor(4);
  histos[2]->Draw("sames");
  l->AddEntry(histos[2], "Non Bending", "l");
  l->Draw("same");

  int i(0);
  for (int i = 1; i < 9; ++i) {
    c->cd(i + 1);
    gPad->SetLogy(logy[i]);
    gPad->SetLogz(logz[i]);
    histos[i + 2]->Draw(opt[i]);
  }

  return c;
}

//_________________________________________________________________________________________________
static const std::vector<double> chargeLimits{200., 400., 700., 1200., 2200., 4000.};

//_________________________________________________________________________________________________
void CreatePreClusterInfoVsWire(std::vector<TH1*>& histos, const char* extension = "")
{
  /// create histograms of precluster characteristics versus distance to closest wire

  histos.emplace_back(new TH2F(Form("hChargeAsymmVsWire%s", extension),
                               "cluster charge asymmetry vs distance to wire;#Delta_x (cm);(NB-B)/(NB+B)",
                               25, -0.125, 0.125, 101, -1.01, 1.01));
  histos.back()->GetYaxis()->SetTitleOffset(1.5);

  for (size_t i = 0; i <= chargeLimits.size(); ++i) {
    auto chargeRange = (i == 0)                    ? fmt::format("charge < {}", chargeLimits[i])
                       : (i < chargeLimits.size()) ? fmt::format("{} <= charge < {}", chargeLimits[i - 1], chargeLimits[i])
                                                   : fmt::format("charge >= {}", chargeLimits[i - 1]);
    histos.emplace_back(new TH2F(Form("hChargeAsymmVsWire%zu%s", i, extension),
                                 Form("cluster charge asymmetry vs distance to wire (%s);#Delta_x (cm);(NB-B)/(NB+B)", chargeRange.c_str()),
                                 25, -0.125, 0.125, 101, -1.01, 1.01));
    histos.back()->GetYaxis()->SetTitleOffset(1.5);
  }

  for (auto h : histos) {
    h->SetDirectory(0);
  }
}

//_________________________________________________________________________________________________
void FillPreClusterInfoVsWire(double chargeNB, double chargeB, float dx, std::vector<TH1*>& histos)
{
  /// fill histograms of precluster characteristics versus distance to closest wire

  double charge = 0.5 * (chargeNB + chargeB);
  double chargeAsymm = (chargeNB - chargeB) / (chargeNB + chargeB);
  auto i = std::distance(chargeLimits.begin(), std::upper_bound(chargeLimits.begin(), chargeLimits.end(), charge));

  histos[0]->Fill(dx, chargeAsymm);
  histos[i + 1]->Fill(dx, chargeAsymm);
}

//_________________________________________________________________________________________________
TCanvas* DrawPreClusterInfoVsWire(std::vector<TH1*>& histos, const char* extension = "")
{
  /// draw histograms of precluster characteristics versus distance to closest wire

  TCanvas* c = new TCanvas(Form("preClusterInfoVsWire%s", extension),
                           Form("precluster characteristics versus distance to wire%s", extension),
                           10, 10, 1200, 600);
  c->Divide((chargeLimits.size() + 3) / 2, 2);

  int i = 0;
  for (auto h : histos) {
    c->cd(++i);
    gPad->SetLogz();
    h->Draw("colz");
  }

  return c;
}

//_________________________________________________________________________________________________
TH3* CreatePreClusterInfo3D(const char* extension = "")
{
  /// create 3D histogram of precluster characteristics

  auto h = new TH3F(Form("hChargeAsymm3D%s", extension),
                    "cluster charge asymmetry vs charge vs distance to wire;#Delta_x (cm);#sqrt{NB*B} (ADC);1/2 * ln(NB/B)",
                    25, -0.125, 0.125, 200, -0.25, 19999.75, 201, -1.005, 1.005);
  h->GetZaxis()->SetTitleOffset(1.5);
  h->SetDirectory(0);

  return h;
}

//_________________________________________________________________________________________________
void FillPreClusterInfo3D(double chargeNB, double chargeB, float dx, TH3* h)
{
  /// fill 3D histogram of precluster characteristics

  double charge = std::sqrt(chargeB * chargeNB);
  double chargeAsymm = 0.5 * std::log(chargeNB / chargeB);

  h->Fill(dx, charge, chargeAsymm);
}

//_________________________________________________________________________________________________
double DoubleGaus(double* x, double* par)
{
  /// double gaussian with same norm and width, placed at -mean and mean
  /// par[0] = Normalization
  /// par[1] = mean
  /// par[2] = sigma

  return par[0] * (TMath::Gaus(x[0], -par[1], par[2], true) + TMath::Gaus(x[0], par[1], par[2], true));
}

//_________________________________________________________________________________________________
std::tuple<TGraphAsymmErrors*, TGraphAsymmErrors*, TGraphAsymmErrors*>
  GetAsymmDispersionVsCharge(TH2* hAsymm2D, const char* extension = "")
{
  /// extract the charge asymmetry dispersion versus charge from the 2D histogram

  static const int minEntriesPerBin = 1000;
  static TF1* fDoubleGaus = new TF1("DoubleGaus", DoubleGaus, -0.5, 0.5, 3);

  std::unique_ptr<TH1> hCharge(hAsymm2D->ProjectionX());
  auto gAsymmRMSvsCharge = new TGraphAsymmErrors();
  gAsymmRMSvsCharge->SetNameTitle(fmt::format("gChargeAsymmRMS{}", extension).c_str(),
                                  fmt::format("{} cluster charge asymmetry RMS vs charge", extension).c_str());
  auto gAsymmDeltavsCharge = new TGraphAsymmErrors();
  gAsymmDeltavsCharge->SetNameTitle(fmt::format("gAsymmDeltavsCharge{}", extension).c_str(),
                                    fmt::format("{} cluster charge asymmetry delta vs charge", extension).c_str());
  auto gAsymmSigmavsCharge = new TGraphAsymmErrors();
  gAsymmSigmavsCharge->SetNameTitle(fmt::format("gAsymmSigmavsCharge{}", extension).c_str(),
                                    fmt::format("{} cluster charge asymmetry sigma vs charge", extension).c_str());

  auto addDispersionValues = [&](int firstBin, int lastBin) {
    hCharge->GetXaxis()->SetRange(firstBin, lastBin);
    auto chargeMean = hCharge->GetMean();
    auto chargeMeanErrLow = chargeMean - hCharge->GetBinLowEdge(firstBin);
    auto chargeMeanErrHigh = hCharge->GetBinLowEdge(lastBin + 1) - chargeMean;

    std::unique_ptr<TH1> hAsymm(hAsymm2D->ProjectionY("tmp", firstBin, lastBin));
    auto asymmStdDevError = hAsymm->GetStdDevError();
    gAsymmRMSvsCharge->AddPoint(chargeMean, hAsymm->GetStdDev());
    gAsymmRMSvsCharge->SetPointError(gAsymmRMSvsCharge->GetN() - 1, chargeMeanErrLow, chargeMeanErrHigh,
                                     asymmStdDevError, asymmStdDevError);

    fDoubleGaus->SetParameters(hAsymm->GetEntries() * hAsymm->GetBinWidth(1) / 2., 0., 0.1);
    auto fitResult = hAsymm->Fit(fDoubleGaus, "LRSNQ");
    auto asymmDeltaError = std::min(std::abs(fitResult->Parameter(1)), fitResult->ParError(1));
    gAsymmDeltavsCharge->AddPoint(chargeMean, std::abs(fitResult->Parameter(1)));
    gAsymmDeltavsCharge->SetPointError(gAsymmDeltavsCharge->GetN() - 1, chargeMeanErrLow, chargeMeanErrHigh,
                                       asymmDeltaError, asymmDeltaError);
    auto asymmSigmaError = fitResult->ParError(2);
    gAsymmSigmavsCharge->AddPoint(chargeMean, fitResult->Parameter(2));
    gAsymmSigmavsCharge->SetPointError(gAsymmSigmavsCharge->GetN() - 1, chargeMeanErrLow, chargeMeanErrHigh,
                                       asymmSigmaError, asymmSigmaError);
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

  if (lastBin < hCharge->GetNbinsX()) {
    addDispersionValues(lastBin + 1, hCharge->GetNbinsX());
  }

  return std::make_tuple(gAsymmRMSvsCharge, gAsymmDeltavsCharge, gAsymmSigmavsCharge);
}

//_________________________________________________________________________________________________
std::tuple<TCanvas*, TCanvas*, TCanvas*> DrawPreClusterInfo3D(TH3* hAsymm3D[4])
{
  /// draw the charge asymmetry dispersion versus charge for different distances to wire

  TCanvas* c1 = new TCanvas("cAsymmRMS", "cluster charge asymmetry RMS vs charge", 10, 10, 1200, 600);
  c1->Divide(2, 2);
  TCanvas* c2 = new TCanvas("cAsymmDelta", "cluster charge asymmetry delta vs charge", 10, 10, 1200, 600);
  c2->Divide(2, 2);
  TCanvas* c3 = new TCanvas("cAsymmSigma", "cluster charge asymmetry sigma vs charge", 10, 10, 1200, 600);
  c3->Divide(2, 2);

  auto draw = [](auto g1, auto g2, auto g3) {
    g1->Draw("ap");
    g1->GetXaxis()->SetRangeUser(0., 20000.);
    g1->GetYaxis()->SetRangeUser(0., 0.4);
    g2->SetMarkerColor(2);
    g2->SetLineColor(2);
    g2->Draw("p");
    g3->SetMarkerColor(4);
    g3->SetLineColor(4);
    g3->Draw("p");
  };

  static const char* sSt[] = {"St1", "St2", "St345", ""};
  for (int iSt = 0; iSt < 4; ++iSt) {

    std::unique_ptr<TH2> hAsymm2D(static_cast<TH2*>(hAsymm3D[iSt]->Project3D("zy")));
    auto [gRMS, gDelta, gSigma] = GetAsymmDispersionVsCharge(hAsymm2D.get(), sSt[iSt]);

    hAsymm3D[iSt]->GetXaxis()->SetRange(12, 14);
    std::unique_ptr<TH2> hAsymm2DClose(static_cast<TH2*>(hAsymm3D[iSt]->Project3D("close_zy")));
    auto [gRMSClose, gDeltaClose, gSigmaClose] = GetAsymmDispersionVsCharge(hAsymm2DClose.get(), sSt[iSt]);

    hAsymm3D[iSt]->GetXaxis()->SetRange(1, 4);
    std::unique_ptr<TH2> hAsymm2DFar(static_cast<TH2*>(hAsymm3D[iSt]->Project3D("far1_zy")));
    hAsymm3D[iSt]->GetXaxis()->SetRange(22, 25);
    std::unique_ptr<TH2> hAsymm2DFar2(static_cast<TH2*>(hAsymm3D[iSt]->Project3D("far2_zy")));
    hAsymm2DFar->Add(hAsymm2DFar2.get());
    auto [gRMSFar, gDeltaFar, gSigmaFar] = GetAsymmDispersionVsCharge(hAsymm2DFar.get(), sSt[iSt]);

    TLegend* l = new TLegend(0.15, 0.75, 0.35, 0.9);
    l->SetFillStyle(0);
    l->SetBorderSize(0);
    l->AddEntry(gRMS, "any distance", "l");
    l->AddEntry(gRMSClose, "|dx| < 0.015 cm", "l");
    l->AddEntry(gRMSFar, "|dx| > 0.085 cm", "l");

    c1->cd(iSt + 1);
    draw(gRMS, gRMSClose, gRMSFar);
    l->Draw("same");

    c2->cd(iSt + 1);
    draw(gDelta, gDeltaClose, gDeltaFar);
    l->Clone()->Draw("same");

    c3->cd(iSt + 1);
    draw(gSigma, gSigmaClose, gSigmaFar);
    l->Clone()->Draw("same");
  }

  return std::make_tuple(c1, c2, c3);
}
