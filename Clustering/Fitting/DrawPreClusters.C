#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>

#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

#include "DataFormatsMCH/Cluster.h"
#include "DataFormatsMCH/Digit.h"
#include "DetectorsBase/GeometryManager.h"
#include "MCHBase/TrackBlock.h"

#include "CCDBUtils.h"
#include "ClusterUtils.h"
#include "DataUtils.h"
#include "DigitUtils.h"
#include "PreClusterUtils.h"

using o2::mch::Cluster;
using o2::mch::Digit;
using o2::mch::TrackParamStruct;

static constexpr double pi = 3.14159265358979323846;

//_________________________________________________________________________________________________
void DrawPreClusters(int run, bool applyTrackSelection = false, bool applyClusterSelection = false,
                     bool applyTimeSelection = false, const char* inFile = "clusters.root")
{
  /// draw precluster and associated digits informations
  /// require the MCH mapping to be loaded: gSystem->Load("libO2MCHMappingImpl4")

  if (run < 300000) {
    o2::base::GeometryManager::loadGeometry("O2geometry.root");
  } else {
    InitFromCCDB(run, false, true, false);
  }

  auto [dataFileIn, dataReader] = LoadData(inFile, "data");
  TTreeReaderValue<TrackParamStruct> trackParam(*dataReader, "trackParameters");
  TTreeReaderValue<int> trackTime(*dataReader, "trackTime");
  TTreeReaderValue<Cluster> cluster(*dataReader, "clusters");
  TTreeReaderValue<std::vector<Digit>> digits(*dataReader, "digits");

  std::vector<TH1*> preClusterInfo{};
  CreatePreClusterInfo(preClusterInfo);
  std::vector<TH1*> preClusterInfoSt[3] = {{}, {}, {}};
  CreatePreClusterInfo(preClusterInfoSt[0], "St1");
  CreatePreClusterInfo(preClusterInfoSt[1], "St2");
  CreatePreClusterInfo(preClusterInfoSt[2], "St345");

  std::vector<TH1*> preClusterInfoVsWireSt[3] = {{}, {}, {}};
  CreatePreClusterInfoVsWire(preClusterInfoVsWireSt[0], "St1");
  CreatePreClusterInfoVsWire(preClusterInfoVsWireSt[1], "St2");
  CreatePreClusterInfoVsWire(preClusterInfoVsWireSt[2], "St345");

  std::vector<TH1*> digitTimeInfo{};
  CreateDigitTimeInfo(digitTimeInfo);
  std::vector<TH1*> digitChargeInfo{};
  CreateDigitChargeInfo(digitChargeInfo);
  std::vector<TH1*> digitChargeInfoSt[3] = {{}, {}, {}};
  CreateDigitChargeInfo(digitChargeInfoSt[0], "St1");
  CreateDigitChargeInfo(digitChargeInfoSt[1], "St2");
  CreateDigitChargeInfo(digitChargeInfoSt[2], "St345");

  while (dataReader->Next()) {

    // those 2 DE have lower HV for the run 529691
    if (run == 529691 && (cluster->getDEId() == 202 || cluster->getDEId() == 300)) {
      continue;
    }

    // cut on track angle at chamber
    if (applyTrackSelection && std::abs(std::atan2(trackParam->py, -trackParam->pz)) / pi * 180. > 10.) {
      continue;
    }

    // cut on track sign
    // if (applyTrackSelection && trackParam->sign > 0) {
    //   continue;
    // }

    // cut on digit time
    std::vector<Digit> selectedDigits(*digits);
    if (applyTimeSelection) {
      selectedDigits.erase(
        std::remove_if(selectedDigits.begin(), selectedDigits.end(), [&trackTime](const auto& digit) {
          return std::abs(digit.getTime() + 1.5 - *trackTime) > 10.;
        }),
        selectedDigits.end());
      if (selectedDigits.empty()) {
        continue;
      }
    }

    // reject mono-cathode clusters after digit selection
    const auto [sizeX, sizeY] = GetSize(selectedDigits);
    if (sizeX == 0 || sizeY == 0) {
      continue;
    }

    float dx = DistanceToClosestWire(cluster->getDEId(), cluster->x, cluster->y, cluster->z, run < 300000);
    // float dx = DistanceToClosestWire(cluster->getDEId(), trackParam->x, trackParam->y, trackParam->z, run < 300000);
    const auto [chargeNB, chargeB] = GetCharge(selectedDigits, run < 300000);
    double chargeAsymm = (chargeNB - chargeB) / (chargeNB + chargeB);

    // cut on distance to closest wire
    // if (applyClusterSelection && std::abs(dx) > 0.015) {
    //   continue;
    // }

    // cut on cluster charge
    // double charge = 0.5 * (chargeNB + chargeB);
    // if (applyClusterSelection && charge < 4000.) {
    //   continue;
    // }

    // cut on cluster charge asymmetry
    if (applyClusterSelection && std::abs(chargeAsymm) > 0.5) {
      continue;
    }
    // if (applyClusterSelection && (chargeAsymm >= -0.2 || chargeAsymm < -0.3)) {
    //   continue;
    // }

    // cut on cluster size
    // if (applyClusterSelection && sizeX < 2) {
    //   continue;
    // }

    // cut on cluster size asymmetry
    // if (applyClusterSelection && (sizeY > sizeX + 3 || sizeX > sizeY + 2)) {
    //   continue;
    // }

    // double maxChargeB = -1;
    // Digit* digitMaxB = nullptr;
    // double maxChargeNB = -1;
    // Digit* digitMaxNB = nullptr;
    // for (auto& digit : selectedDigits) {
    //   const auto& segmentation = o2::mch::mapping::segmentation(digit.getDetID());
    //   if (segmentation.isBendingPad(digit.getPadID())) {
    //     if (digit.getADC() > maxChargeB) {
    //       digitMaxB = &digit;
    //       maxChargeB = digit.getADC();
    //     }
    //   } else {
    //     if (digit.getADC() > maxChargeNB) {
    //       digitMaxNB = &digit;
    //       maxChargeNB = digit.getADC();
    //     }
    //   }
    // }

    FillPreClusterInfo(chargeNB, chargeB, sizeX, sizeY, preClusterInfo);
    int iSt = (cluster->getChamberId() < 4) ? cluster->getChamberId() / 2 : 2;
    FillPreClusterInfo(chargeNB, chargeB, sizeX, sizeY, preClusterInfoSt[iSt]);
    FillPreClusterInfoVsWire(chargeNB, chargeB, dx, preClusterInfoVsWireSt[iSt]);

    for (const auto& digit : selectedDigits) {
      // if (!(&digit == digitMaxB || &digit == digitMaxNB)) {
      //   continue;
      // }
      FillDigitTimeInfo(digit, *trackTime, digitTimeInfo);
      FillDigitChargeInfo(digit, digitChargeInfo, chargeAsymm, run < 300000);
      FillDigitChargeInfo(digit, digitChargeInfoSt[iSt], chargeAsymm, run < 300000);
    }
  }

  gStyle->SetOptStat(1);

  auto c = DrawPreClusterInfo(preClusterInfo);
  auto cSt1 = DrawPreClusterInfo(preClusterInfoSt[0], "St1");
  auto cSt2 = DrawPreClusterInfo(preClusterInfoSt[1], "St2");
  auto cSt345 = DrawPreClusterInfo(preClusterInfoSt[2], "St345");

  auto cwSt1 = DrawPreClusterInfoVsWire(preClusterInfoVsWireSt[0], "St1");
  auto cwSt2 = DrawPreClusterInfoVsWire(preClusterInfoVsWireSt[1], "St2");
  auto cwSt345 = DrawPreClusterInfoVsWire(preClusterInfoVsWireSt[2], "St345");

  auto ct = DrawDigitTimeInfo(digitTimeInfo);
  auto cc = DrawDigitChargeInfo(digitChargeInfo);
  auto ccSt1 = DrawDigitChargeInfo(digitChargeInfoSt[0], "St1");
  auto ccSt2 = DrawDigitChargeInfo(digitChargeInfoSt[1], "St2");
  auto ccSt345 = DrawDigitChargeInfo(digitChargeInfoSt[2], "St345");

  dataFileIn->Close();
}
