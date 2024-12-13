#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <fmt/format.h>

#include <TCanvas.h>
#include <TColor.h>
#include <TDialogCanvas.h>
#include <TBox.h>
#include <TButton.h>
#include <TFile.h>
#include <TGeoManager.h>
#include <TH1F.h>
#include <TMarker.h>
#include <TPad.h>
#include <TPaletteAxis.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TText.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

#include "DataFormatsMCH/Cluster.h"
#include "DataFormatsMCH/Digit.h"
#include "DetectorsBase/GeometryManager.h"
#include "MathUtils/Cartesian.h"
#include "MCHBase/TrackBlock.h"
#include "MCHGeometryTransformer/Transformations.h"
#include "MCHMappingInterface/Segmentation.h"

#include "CCDBUtils.h"
#include "DataUtils.h"
#include "DigitUtils.h"
#include "PreClusterUtils.h"

using o2::mch::Cluster;
using o2::mch::Digit;
using o2::mch::TrackParamStruct;

static constexpr double pi = 3.14159265358979323846;
TPad* mainPad = nullptr;
TPad* bendingPad = nullptr;
TPad* nonBendingPad = nullptr;

void Next(int run, const char* inFile, bool backward = false);
bool IsSelected(int run, const TrackParamStruct& trackParam, int trackTime,
                const Cluster& cluster, std::vector<Digit>& selectedDigits);
void Display(int entry, int nEntries, const std::vector<Digit>& digits, const Cluster& cluster, bool run2);
void Display(double x, double y, std::string text, std::string name = "txt");
void DigitView(const std::vector<Digit>& digits, std::vector<TBox*> pads[2], std::vector<double> charges[2],
               double& xMin, double& yMin, double& xMax, double& yMax,
               double& minCharge, double& maxCharge, bool run2 = false);
TMarker* ClusterView(const Cluster& cluster, bool run2 = false);

//_________________________________________________________________________________________________
void VisualizePreClusters(int run, const char* inFile = "clusters.root")
{
  /// visualize preclusters
  /// require the MCH mapping to be loaded: gSystem->Load("libO2MCHMappingImpl4")

  gStyle->SetPalette(kSolar);
  TColor::InvertPalette();

  if (run < 300000) {
    o2::base::GeometryManager::loadGeometry("O2geometry.root");
  } else {
    InitFromCCDB(run, false, true, false);
  }

  TDialogCanvas* dialog = new TDialogCanvas("dialog", "precluster display", 500, 500);
  mainPad = new TPad("display", "", 0.03, 0.1, 0.97, 0.97);
  mainPad->SetRightMargin(0.15);
  mainPad->Draw();
  std::string command = fmt::format("Next({}, \"{}\", true)", run, inFile);
  TButton* previous = new TButton("previous", command.c_str(), 0.03, 0.03, 0.33, 0.08);
  previous->Draw();
  command = fmt::format("Next({}, \"{}\")", run, inFile);
  TButton* next = new TButton("next", command.c_str(), 0.35, 0.03, 0.65, 0.08);
  next->Draw();
  TButton* quit = new TButton("quit", "exit(0)", 0.67, 0.03, 0.97, 0.08);
  quit->Draw();

  bendingPad = new TCanvas("cBending", "precluster bending display", 510, 0, 470, 455);
  bendingPad->SetRightMargin(0.15);

  nonBendingPad = new TCanvas("cNonBending", "precluster non-bending display", 990, 0, 470, 455);
  nonBendingPad->SetRightMargin(0.15);

  Next(run, inFile);
}

//_________________________________________________________________________________________________
void Next(int run, const char* inFile, bool backward)
{
  /// display next precluster

  static auto [dataFileIn, dataReader] = LoadData(inFile, "data");
  static TTreeReaderValue<TrackParamStruct> trackParam(*dataReader, "trackParameters");
  static TTreeReaderValue<int> trackTime(*dataReader, "trackTime");
  static TTreeReaderValue<Cluster> cluster(*dataReader, "clusters");
  static TTreeReaderValue<std::vector<Digit>> digits(*dataReader, "digits");
  static int currentDisplay = -1;

  while (true) {

    int nextEntry = backward ? dataReader->GetCurrentEntry() - 1 : dataReader->GetCurrentEntry() + 1;
    if (nextEntry < 0) {
      mainPad->cd();
      Display(0.19, 0.15, "this is the first selected precluster");
    } else if (nextEntry >= dataReader->GetEntries()) {
      mainPad->cd();
      Display(0.19, 0.15, "this is the last selected precluster");
    } else if (dataReader->SetEntry(nextEntry) == TTreeReader::kEntryValid) {
      std::vector<Digit> selectedDigits(*digits);
      if (!IsSelected(run, *trackParam, *trackTime, *cluster, selectedDigits)) {
        continue;
      }
      Display(nextEntry, dataReader->GetEntries(), selectedDigits, *cluster, run < 300000);
      currentDisplay = nextEntry;
    } else {
      std::cout << "Error reading entry " << nextEntry << "... Exiting." << std::endl;
      exit(-1);
    }

    if (currentDisplay >= 0 && dataReader->GetCurrentEntry() != currentDisplay) {
      dataReader->SetEntry(currentDisplay);
    }

    break;
  }
}

//_________________________________________________________________________________________________
bool IsSelected(int run, const TrackParamStruct& trackParam, int trackTime,
                const Cluster& cluster, std::vector<Digit>& selectedDigits)
{
  /// apply precluster selections

  // those 2 DE have lower HV for the run 529691
  if (run == 529691 && (cluster.getDEId() == 202 || cluster.getDEId() == 300)) {
    return false;
  }

  // selection stations 3, 4 and 5
  // if (cluster.getChamberId() < 4) {
  //   return false;
  // }

  // cut on track angle at chamber
  if (std::abs(std::atan2(trackParam.py, -trackParam.pz)) / pi * 180. > 10.) {
    return false;
  }

  // cut on digit time
  selectedDigits.erase(
    std::remove_if(selectedDigits.begin(), selectedDigits.end(), [trackTime](const auto& digit) {
      return std::abs(digit.getTime() + 1.5 - trackTime) > 10.;
    }),
    selectedDigits.end());
  if (selectedDigits.empty()) {
    return false;
  }

  // cut on cluster charge asymmetry
  const auto [chargeNB, chargeB] = GetCharge(selectedDigits, run < 300000);
  double chargeAsymm = (chargeNB - chargeB) / (chargeNB + chargeB);
  if (std::abs(chargeAsymm) > 0.5) {
    return false;
  }

  // cut on cluster charge
  // double charge = 0.5 * (chargeNB + chargeB);
  // if (charge < 4000.) {
  //   return false;
  // }

  // cut on cluster size asymmetry
  // const auto [sizeX, sizeY] = GetSize(selectedDigits);
  // if (sizeY > sizeX + 3 || sizeX > sizeY + 2) {
  //   return false;
  // }

  return true;
}

//_________________________________________________________________________________________________
void Display(int entry, int nEntries, const std::vector<Digit>& digits, const Cluster& cluster, bool run2)
{
  /// display precluster

  std::vector<TBox*> pads[2];
  std::vector<double> charges[2];
  double xMin, yMin, xMax, yMax, minCharge, maxCharge;
  DigitView(digits, pads, charges, xMin, yMin, xMax, yMax, minCharge, maxCharge, run2);

  std::string title = fmt::format("precluster {}/{}, DE {}", entry + 1, nEntries, digits[0].getDetID());
  double xMargin = 0.2 * (xMax - xMin);
  double yMargin = 0.2 * (yMax - yMin);
  for (auto pad : {mainPad, bendingPad, nonBendingPad}) {
    pad->cd();
    TH1F* hFrame = gPad->DrawFrame(xMin - xMargin, yMin - yMargin, xMax + xMargin, yMax + yMargin, title.c_str());
    hFrame->GetXaxis()->SetLabelSize(0.04);
    hFrame->GetYaxis()->SetLabelSize(0.04);
  }

  TPaletteAxis* palette = new TPaletteAxis(
    xMax + xMargin, yMin - yMargin, xMax + xMargin + 0.05 * (xMax - xMin + 2 * xMargin), yMax + yMargin,
    minCharge, maxCharge);
  palette->SetBit(TObject::kCanDelete);
  palette->SetNdivisions(10);
  palette->SetLabelFont(40);
  palette->SetLabelSize(0.04);
  mainPad->cd();
  palette->Draw();
  for (auto pad : {bendingPad, nonBendingPad}) {
    pad->cd();
    auto* clone = palette->Clone();
    clone->SetBit(TObject::kCanDelete);
    clone->Draw();
  }

  mainPad->cd();
  for (auto pad : pads[1]) {
    pad->Draw();
  }
  for (auto pad : pads[0]) {
    pad->Draw();
  }

  bendingPad->cd();
  for (size_t i = 0; i < pads[1].size(); ++i) {
    auto* clone = static_cast<TBox*>(pads[1][i]->Clone());
    clone->SetBit(TObject::kCanDelete);
    clone->SetToolTipText(fmt::format("{:g}", charges[1][i]).c_str(), 10);
    clone->Draw();
  }

  nonBendingPad->cd();
  for (size_t i = 0; i < pads[0].size(); ++i) {
    auto* clone = static_cast<TBox*>(pads[0][i]->Clone());
    clone->SetFillColorAlpha(pads[0][i]->GetFillColor(), 1.);
    clone->SetToolTipText(fmt::format("{:g}", charges[0][i]).c_str(), 10);
    clone->SetBit(TObject::kCanDelete);
    clone->Draw();
  }

  auto point = ClusterView(cluster, run2);
  mainPad->cd();
  point->Draw();
  for (auto pad : {bendingPad, nonBendingPad}) {
    pad->cd();
    auto* clone = point->Clone();
    clone->SetBit(TObject::kCanDelete);
    clone->Draw();
  }

  const auto [chargeNB, chargeB] = GetCharge(digits, run2);
  double chargeAsymm = (chargeNB - chargeB) / (chargeNB + chargeB);
  mainPad->cd();
  Display(0.23, 0.85, fmt::format("charge asymmetry = {:+.4f}", chargeAsymm).c_str(), "asymm");
  bendingPad->cd();
  Display(0.28, 0.85, fmt::format("bending charge = {:g}", chargeB).c_str(), "chargeB");
  nonBendingPad->cd();
  Display(0.24, 0.85, fmt::format("non-bending charge = {:g}", chargeNB).c_str(), "chargeNB");

  bendingPad->Update();
  nonBendingPad->Update();
}

//_________________________________________________________________________________________________
void Display(double x, double y, std::string text, std::string name)
{
  /// display text

  gPad->RecursiveRemove(gPad->FindObject(name.c_str()));

  TText* txt = new TText(x, y, text.c_str());
  txt->SetBit(TObject::kCanDelete);
  txt->SetName(name.c_str());
  txt->SetNDC();
  txt->SetTextFont(40);
  txt->SetTextSize(0.045);
  txt->Draw();
}

//_________________________________________________________________________________________________
void DigitView(const std::vector<Digit>& digits, std::vector<TBox*> pads[2], std::vector<double> charges[2],
               double& xMin, double& yMin, double& xMax, double& yMax,
               double& minCharge, double& maxCharge, bool run2)
{
  /// make the visual objects for digits

  const auto& segmentation = o2::mch::mapping::segmentation(digits[0].getDetID());

  minCharge = 1.e6;
  maxCharge = -1.e6;

  for (const auto& digit : digits) {
    auto charge = run2 ? adcToCharge(digit.getADC()) : digit.getADC();
    minCharge = std::min(minCharge, charge);
    maxCharge = std::max(maxCharge, charge);
  }

  double chargeToColor = (gStyle->GetNumberOfColors() - 1.) / (maxCharge - minCharge);

  xMin = 1.e6;
  yMin = 1.e6;
  xMax = -1.e6;
  yMax = -1.e6;

  for (const auto& digit : digits) {

    double x = segmentation.padPositionX(digit.getPadID());
    double y = segmentation.padPositionY(digit.getPadID());
    double hx = 0.5 * segmentation.padSizeX(digit.getPadID());
    double hy = 0.5 * segmentation.padSizeY(digit.getPadID());
    bool isBending = segmentation.isBendingPad(digit.getPadID());
    int i = isBending ? 1 : 0;
    float transparency = isBending ? 1. : 0.6;
    auto charge = run2 ? adcToCharge(digit.getADC()) : digit.getADC();
    int color = std::round((charge - minCharge) * chargeToColor);

    pads[i].emplace_back(new TBox(x - hx, y - hy, x + hx, y + hy));
    pads[i].back()->SetFillColorAlpha(gStyle->GetColorPalette(color), transparency);
    pads[i].back()->SetToolTipText(fmt::format("{:g}", charge).c_str(), 10);
    pads[i].back()->SetBit(TObject::kCanDelete);
    charges[i].emplace_back(charge);

    xMin = std::min(xMin, x - hx);
    yMin = std::min(yMin, y - hy);
    xMax = std::max(xMax, x + hx);
    yMax = std::max(yMax, y + hy);
  }
}

//_________________________________________________________________________________________________
TMarker* ClusterView(const Cluster& cluster, bool run2)
{
  /// make the visual object for cluster

  static o2::mch::geo::TransformationCreator transformation;
  if (!transformation) {
    if (run2) {
      std::ifstream geoFile("AlignedGeometry.json");
      if (!geoFile.is_open()) {
        std::cout << "cannot open geometry file AlignedGeometry.json" << std::endl;
        exit(-1);
      }
      transformation = o2::mch::geo::transformationFromJSON(geoFile);
    } else {
      transformation = o2::mch::geo::transformationFromTGeoManager(*gGeoManager);
    }
  }

  auto de = cluster.getDEId();
  o2::math_utils::Point3D<float> global{cluster.x, cluster.y, cluster.z};
  auto local = transformation(de) ^ global;

  TMarker* point = new TMarker(local.x(), local.y(), kFullDotLarge);
  point->SetBit(TObject::kCanDelete);

  return point;
}
