#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include <fmt/format.h>

#include <TColor.h>
#include <TDialogCanvas.h>
#include <TBox.h>
#include <TButton.h>
#include <TFile.h>
#include <TGeoManager.h>
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

void Next(int run, const char* inFile, bool backward = false);
bool IsSelected(int run, const TrackParamStruct& trackParam, int trackTime,
                const Cluster& cluster, std::vector<Digit>& selectedDigits);
void Display(int entry, int nEntries, const std::vector<Digit>& digits, const Cluster& cluster, bool run2);
void Display(std::string text);
void DigitView(const std::vector<Digit>& digits, std::vector<TBox*> pads[2],
               double& xMin, double& yMin, double& xMax, double& yMax, bool run2 = false);
TMarker* ClusterView(const Cluster& cluster);

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
  TPad* pad = new TPad("display", "", 0.05, 0.15, 0.95, 0.95);
  pad->Draw();
  std::string command = fmt::format("Next({}, \"{}\", true)", run, inFile);
  TButton* previous = new TButton("previous", command.c_str(), 0.05, 0.05, 0.32, 0.1);
  previous->Draw();
  command = fmt::format("Next({}, \"{}\")", run, inFile);
  TButton* next = new TButton("next", command.c_str(), 0.365, 0.05, 0.635, 0.1);
  next->Draw();
  TButton* quit = new TButton("quit", "exit(0)", 0.68, 0.05, 0.95, 0.1);
  quit->Draw();

  pad->cd();
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
      Display("this is the first selected precluster");
    } else if (nextEntry >= dataReader->GetEntries()) {
      Display("this is the last selected precluster");
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
  double xMin, yMin, xMax, yMax;
  DigitView(digits, pads, xMin, yMin, xMax, yMax, run2);

  std::string title = fmt::format("precluster {}/{}, DE {}", entry, nEntries, digits[0].getDetID());
  gPad->DrawFrame(xMin - 1., yMin - 1., xMax + 1., yMax + 1., title.c_str());

  TPaletteAxis* palette = new TPaletteAxis(xMax + 1., yMin - 1., xMax + 1. + 0.05 * (xMax - xMin + 2), yMax + 1., 0., 1.);
  palette->SetBit(TObject::kCanDelete);
  palette->Draw();

  for (auto pad : pads[1]) {
    pad->Draw();
  }

  for (auto pad : pads[0]) {
    pad->Draw();
  }

  auto point = ClusterView(cluster);
  point->Draw();
}

//_________________________________________________________________________________________________
void Display(std::string text)
{
  /// display text

  gPad->RecursiveRemove(gPad->FindObject("txt"));

  TText* txt = new TText(0.2, 0.15, text.c_str());
  txt->SetBit(TObject::kCanDelete);
  txt->SetName("txt");
  txt->SetNDC();
  txt->SetTextFont(40);
  txt->Draw();
}

//_________________________________________________________________________________________________
void DigitView(const std::vector<Digit>& digits, std::vector<TBox*> pads[2],
               double& xMin, double& yMin, double& xMax, double& yMax, bool run2)
{
  /// make the visual objects for digits

  const auto& segmentation = o2::mch::mapping::segmentation(digits[0].getDetID());

  double minCharge = 1.e6;
  double maxCharge = -1.e6;
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
    pads[i].back()->SetBit(TObject::kCanDelete);

    xMin = std::min(xMin, x - hx);
    yMin = std::min(yMin, y - hy);
    xMax = std::max(xMax, x + hx);
    yMax = std::max(yMax, y + hy);
  }
}

//_________________________________________________________________________________________________
TMarker* ClusterView(const Cluster& cluster)
{
  /// make the visual object for cluster

  static auto transformation = o2::mch::geo::transformationFromTGeoManager(*gGeoManager);

  auto de = cluster.getDEId();
  o2::math_utils::Point3D<float> global{cluster.x, cluster.y, cluster.z};
  auto local = transformation(de) ^ (global);

  TMarker* point = new TMarker(local.x(), local.y(), kFullDotLarge);
  point->SetBit(TObject::kCanDelete);

  return point;
}
