#include <algorithm>
#include <array>
#include <cmath>
#include <random>
#include <string>
#include <vector>

#include "DataFormatsMCH/Digit.h"
#include "MCHMappingInterface/Segmentation.h"
#include "MCHSimulation/Response.h"

#include "PreClusterUtils.h"

using o2::mch::Digit;
using o2::mch::Response;

// generate random numbers
std::mt19937 mRandom{std::random_device{}()};

//_________________________________________________________________________________________________
void AddNoise(double& charge, uint32_t nSamples, std::string mode)
{
  auto salpha = mode.substr(mode.find('_') + 1);
  std::replace(salpha.begin(), salpha.end(), 'p', '.');
  double alpha = std::stod(salpha);

  if (mode.starts_with("MC_")) {

    // for MC th. noise
    static std::normal_distribution mNoise{0., 0.5};
    charge += mNoise(mRandom) * (std::sqrt(nSamples) + alpha);

  } else if (mode.starts_with("sADC_")) {

    // sqrt of ADC noise
    static std::normal_distribution mNoise{0., 1.};
    charge += mNoise(mRandom) * std::sqrt(charge) * alpha;
  }
}

//_________________________________________________________________________________________________
void GenerateAsymm(double& chargeB, double& chargeNB, std::string mode)
{
  double charge = std::sqrt(chargeB * chargeNB);

  if (mode == "none") {

    chargeB = charge;
    chargeNB = charge;

  } else if (mode.starts_with("gaus_")) {

    mode.erase(0, 5);
    std::replace(mode.begin(), mode.end(), 'p', '.');
    auto gamma = std::stod(mode);

    static std::normal_distribution mAsym{0., 0.055};
    double Y = gamma * mAsym(mRandom);
    chargeB = charge * exp(Y);
    chargeNB = charge * exp(-Y);

  } else if (mode == "tripleGaus") {
    // to be done
  }
}

//_________________________________________________________________________________________________
bool IsAboveThreshold(double charge, std::string mode)
{
  if (mode == "gaus") {

    static std::normal_distribution mThreshold{22.2, 2.8};
    return charge > mThreshold(mRandom);

  } else if (mode == "uniform") {

    return charge > 22.2;
  }

  return charge > 0.;
}

//_________________________________________________________________________________________________
void TMC(std::vector<Digit>& digits, int32_t time, int deId, std::array<double, 6> param,
         std::string asymm, std::string noise, std::string threshold)
{
  static const Response response[] = {{o2::mch::Station::Type1}, {o2::mch::Station::Type2345}};
  const auto& mSegmentation = o2::mch::mapping::segmentation(deId);
  int iSt = (deId < 300) ? 0 : 1;

  // generate charge asymmetry if needed
  if (asymm != "copy") {
    GenerateAsymm(param[4], param[5], asymm);
  }

  // borders of charge integration area
  auto dxy = response[iSt].getSigmaIntegration() * response[iSt].getChargeSpread();
  auto xMin = param[0] - dxy;
  auto xMax = param[0] + dxy;
  auto yMin = param[1] - dxy;
  auto yMax = param[1] + dxy;

  mSegmentation.forEachPadInArea(xMin, yMin, xMax, yMax, [&](int padid) {
    auto dx = mSegmentation.padSizeX(padid) * 0.5;
    auto dy = mSegmentation.padSizeY(padid) * 0.5;
    auto xPad = mSegmentation.padPositionX(padid) - param[0];
    auto yPad = mSegmentation.padPositionY(padid) - param[1];
    double q = response[iSt].chargePadfraction(xPad - dx, xPad + dx, yPad - dy, yPad + dy);
    if (response[iSt].isAboveThreshold(q)) {
      q *= mSegmentation.isBendingPad(padid) ? param[4] : param[5];
      auto nSamples = response[iSt].nSamples(q);
      if (noise != "none") {
        AddNoise(q, nSamples, noise);
      }
      if (IsAboveThreshold(q, threshold)) {
        digits.emplace_back(deId, padid, std::round(q), time - 2, std::min(nSamples, 0x3FFU), false);
      }
    }
  });
}
