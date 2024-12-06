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
#include "MCHBase/TrackBlock.h"

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

    const auto [sizeX, sizeY] = GetSize(selectedDigits);
    const auto [chargeNB, chargeB] = GetCharge(selectedDigits, run < 300000);
    double charge = 0.5 * (chargeNB + chargeB);
    double chargeAsymm = (chargeNB - chargeB) / (chargeNB + chargeB);

    // cut on cluster charge
    // if (applyClusterSelection && charge < 4000.) {
    //   continue;
    // }

    // cut on cluster charge asymmetry
    if (applyClusterSelection && std::abs(chargeAsymm) > 0.5) {
      continue;
    }

    // cut on cluster size asymmetry
    // if (applyClusterSelection && (sizeY > sizeX + 3 || sizeX > sizeY + 2)) {
    //   continue;
    // }

    FillPreClusterInfo(charge, chargeAsymm, sizeX, sizeY, preClusterInfo);
    int iSt = (cluster->getChamberId() < 4) ? cluster->getChamberId() / 2 : 2;
    FillPreClusterInfo(charge, chargeAsymm, sizeX, sizeY, preClusterInfoSt[iSt]);

    for (const auto& digit : selectedDigits) {
      FillDigitTimeInfo(digit, *trackTime, digitTimeInfo);
      FillDigitChargeInfo(digit, digitChargeInfo, run < 300000);
      FillDigitChargeInfo(digit, digitChargeInfoSt[iSt], run < 300000);
    }
  }

  gStyle->SetOptStat(1);

  auto c = DrawPreClusterInfo(preClusterInfo);
  auto cSt1 = DrawPreClusterInfo(preClusterInfoSt[0], "St1");
  auto cSt2 = DrawPreClusterInfo(preClusterInfoSt[1], "St2");
  auto cSt345 = DrawPreClusterInfo(preClusterInfoSt[2], "St345");

  auto ct = DrawDigitTimeInfo(digitTimeInfo);
  auto cc = DrawDigitChargeInfo(digitChargeInfo);
  auto ccSt1 = DrawDigitChargeInfo(digitChargeInfoSt[0], "St1");
  auto ccSt2 = DrawDigitChargeInfo(digitChargeInfoSt[1], "St2");
  auto ccSt345 = DrawDigitChargeInfo(digitChargeInfoSt[2], "St345");

  dataFileIn->Close();
}
