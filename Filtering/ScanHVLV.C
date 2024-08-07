#include <algorithm>
#include <ctime>
#include <fstream>
#include <iterator>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#include <fmt/format.h>

#include "TCanvas.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TLine.h"
#include "TMultiGraph.h"
#include "TStyle.h"

#include "CCDB/BasicCCDBManager.h"
#include "CCDB/CcdbApi.h"
#include "DetectorsDCS/DataPointIdentifier.h"
#include "DetectorsDCS/DataPointValue.h"
#include "MCHConditions/DCSAliases.h"
#include "MCHStatus/HVStatusCreator.h"
#include "MCHStatus/StatusMapCreatorParam.h"

using namespace o2;
using DPID = dcs::DataPointIdentifier;
using DPVAL = dcs::DataPointValue;
using DPMAP = std::unordered_map<DPID, std::vector<DPVAL>>;
using DPMAP2 = std::map<std::string, std::map<uint64_t, double>>;
using RBMAP = std::map<int, std::pair<uint64_t, uint64_t>>;
using DPBMAP = std::map<uint64_t, uint64_t>;
using ISSUE = std::tuple<uint64_t, uint64_t, double, double, std::string>;
using ISSUELIST = std::vector<ISSUE>;
using ISSUEMAP = std::map<std::string, ISSUELIST>;

double hvRange[2] = {-10., 1700.};
double hvLimits[10] = {1550., 1550., 1600., 1600., 1600., 1600., 1600., 1600., 1600., 1600.};
const char* hvFormat = "%7.2f";

double lvRange[3] = {-1., 4., 8.};   // min, max FeeAnalog/FeeDigital, max Solar
double lvLimits[3] = {1.5, 1.5, 6.}; // FeeAnalog, FeeDigital, Solar
const char* lvFormat = "%4.2f";

bool ContainsAKey(std::string data, const std::set<std::string>& Keys);
bool isValid(std::string alias);
void ScanWhat(std::string what, std::string& path, bool& scanHV, bool& scanAll, std::set<std::string>& aliases);
std::set<int> GetRuns(std::string runList);
RBMAP GetRunBoundaries(ccdb::CcdbApi const& api, std::string runList);
void CheckRunBoundaries(const RBMAP& runBoundaries);
void PrintRunBoundaries(const RBMAP& runBoundaries);
void DrawRunBoudaries(const RBMAP& runBoundaries, TCanvas* c);
DPBMAP GetDPBoundaries(ccdb::CcdbApi const& api, std::string what, uint64_t tStart, uint64_t tStop);
void CheckDPBoundaries(const DPBMAP& dpBoundaries, bool scanHV, uint64_t tStart, uint64_t tStop);
void PrintDPBoundaries(const DPBMAP& dpBoundaries, bool scanHV);
double GetLVLimit(std::string alias);
void DrawLimit(double limit, TCanvas* c);
uint64_t MSToS(uint64_t ts);
std::string GetTime(uint64_t ts);
std::string GetDuration(uint64_t tStart, uint64_t tStop);
double GetValue(DPVAL dp);
std::string GetDE(std::string alias);
void FillDataPoints(const std::vector<DPVAL>& dps, std::map<uint64_t, double>& dps2, uint64_t tMin, uint64_t tMax, int warningLevel);
void SelectDataPoints(DPMAP2 dpsMapsPerCh[10], uint64_t tStart, uint64_t tStop);
void PrintDataPoints(const DPMAP2 dpsMapsPerCh[10], bool scanHV, bool all);
TGraph* MapToGraph(std::string alias, const std::map<uint64_t, double>& dps);
TCanvas* DrawDataPoints(TMultiGraph* mg, double min, double max);
void FindIssues(const std::map<uint64_t, double>& dps, double limit, ISSUELIST& issues);
void FillO2Issues(const std::vector<mch::HVStatusCreator::TimeRange>& o2issues, ISSUELIST& issues, uint64_t tMin, uint64_t tMax);
void SelectIssues(ISSUEMAP issuesPerCh[10], const RBMAP& runBoundaries, uint64_t minDuration);
void SelectO2Issues(ISSUEMAP issuesPerCh[10], const RBMAP& runBoundaries);
std::string FindAffectedRuns(const RBMAP& runBoundaries, uint64_t tStart, uint64_t tStop);
bool eraseIssue(const ISSUE& issue, ISSUELIST& issues);
void PrintIssues(const ISSUEMAP issuesPerCh[10], const ISSUEMAP o2IssuesPerCh[10], bool scanHV);

//----------------------------------------------------------------------------
void ScanHVLV(std::string runList, std::string what = "HV", uint64_t minDuration = 0,
              int warningLevel = 1, int printLevel = 1)
{
  /// scan every HV or LV channels to look for issues
  /// runList = an ASCII file with one run per line
  ///         = a comma separated run list
  ///         = a single run
  /// what = "HV" to scan all HV channels
  ///      = "LV" to scan all LV channels
  ///      = a comma separated list of (part of) aliases to consider, which must be all of the
  ///        same type, i.e contain either Quad/Slat (type = HV), or Group/an/di/Sol (type = LV)
  /// minDuration = minimum duration (ms) of HV/LV issues to consider
  /// warningLevel == 1: check data points w.r.t. HV/LV file validity range with ±2s tolerance
  /// warningLevel >= 2: check data points w.r.t. HV/LV file validity range
  /// printLevel >= 1: print time stamps of runs and HV/LV files
  /// printLevel >= 2: print the first and last data points of each selected channel
  /// printLevel >= 3: print all the data points of each selected channel

  gStyle->SetPalette(kVisibleSpectrum);

  // setup O2 algorithm for searching HV issues
  conf::ConfigurableParam::setValue("MCHStatusMap.hvMinDuration", std::to_string(minDuration));
  conf::ConfigurableParam::setValue("MCHStatusMap.timeMargin", "0");

  // determine what is scanned
  std::string path{};
  bool scanHV = false;
  bool scanAll = false;
  std::set<std::string> aliases{};
  ScanWhat(what, path, scanHV, scanAll, aliases);

  ccdb::CcdbApi api;
  api.init("http://alice-ccdb.cern.ch");

  // get the SOR/EOR of every runs from the list, ordered in run number
  auto runBoundaries = GetRunBoundaries(api, runList);
  if (printLevel > 0) {
    PrintRunBoundaries(runBoundaries);
  }
  CheckRunBoundaries(runBoundaries);

  // extract the time boundaries for each HV/LV file in the full time range
  auto dpBoundaries = GetDPBoundaries(api, path.c_str(), runBoundaries.begin()->second.first,
                                      runBoundaries.rbegin()->second.second);
  if (printLevel > 0) {
    PrintDPBoundaries(dpBoundaries, scanHV);
  }
  CheckDPBoundaries(dpBoundaries, scanHV, runBoundaries.begin()->second.first,
                    runBoundaries.rbegin()->second.second);

  // loop over the HV/LV files, fill the lists of data points per chamber and find issues using O2 algorithm
  DPMAP2 dpsMapsPerCh[10];
  mch::HVStatusCreator hvStatusCreator{};
  ISSUEMAP o2issuesPerCh[10];
  std::map<std::string, std::string> metadata;
  for (auto boundaries : dpBoundaries) {

    auto* dpMap = api.retrieveFromTFileAny<DPMAP>(path.c_str(), metadata, boundaries.first);

    // fill the lists of data points per chamber for requested aliases
    for (const auto& [dpid, dps] : *dpMap) {
      std::string alias(dpid.get_alias());
      if (!mch::dcs::isValid(alias)) {
        printf("error: invalid DCS alias: %s\n", alias.c_str());
        exit(1);
      }
      if ((scanAll || ContainsAKey(alias, aliases)) && (!scanHV || alias.find(".iMon") == alias.npos)) {
        int chamber = mch::dcs::toInt(mch::dcs::aliasToChamber(alias));
        FillDataPoints(dps, dpsMapsPerCh[chamber][alias], boundaries.first, boundaries.second, warningLevel);
      }
    }

    // find issues for requested aliases using O2 algorithm (only for HV)
    if (scanHV) {
      hvStatusCreator.findBadHVs(*dpMap);
      for (const auto& [alias, issues] : hvStatusCreator.getBadHVs()) {
        if (scanAll || ContainsAKey(alias, aliases)) {
          int chamber = mch::dcs::toInt(mch::dcs::aliasToChamber(alias));
          FillO2Issues(issues, o2issuesPerCh[chamber][alias], boundaries.first, boundaries.second);
        }
      }
    }
  }
  if (printLevel > 1) {
    PrintDataPoints(dpsMapsPerCh, scanHV, printLevel > 2);
  }

  // select the data points in the time range
  SelectDataPoints(dpsMapsPerCh, runBoundaries.begin()->second.first, runBoundaries.rbegin()->second.second);
  if (printLevel > 1) {
    PrintDataPoints(dpsMapsPerCh, scanHV, printLevel > 2);
  }

  // create and fill the graphs, and find HV/LV issues
  ISSUEMAP issuesPerCh[10];
  TMultiGraph* mg[10];
  std::set<double> limits;
  for (int ch = 0; ch < 10; ++ch) {
    mg[ch] = new TMultiGraph;
    mg[ch]->SetNameTitle(fmt::format("ch{}", ch).c_str(),
                         fmt::format("chamber {};time;{} (V)", ch + 1, scanHV ? "HV" : "LV").c_str());
    for (const auto& [alias, dps] : dpsMapsPerCh[ch]) {
      mg[ch]->Add(MapToGraph(alias, dps), "lp");
      auto limit = scanHV ? hvLimits[ch] : GetLVLimit(alias);
      limits.emplace(limit);
      FindIssues(dps, limit, issuesPerCh[ch][alias]);
    }
  }

  // select HV/LV issues of a minimum duration (ms) occurring during runs
  SelectIssues(issuesPerCh, runBoundaries, minDuration);
  SelectO2Issues(o2issuesPerCh, runBoundaries);
  PrintIssues(issuesPerCh, o2issuesPerCh, scanHV);

  // display
  for (int ch = 0; ch < 10; ++ch) {
    TCanvas* c;
    if (scanHV) {
      c = DrawDataPoints(mg[ch], hvRange[0], hvRange[1]);
      DrawLimit(hvLimits[ch], c);
    } else {
      auto lvMax = (what.find("LV") != what.npos || what.find("Sol") != what.npos) ? lvRange[2] : lvRange[1];
      c = DrawDataPoints(mg[ch], lvRange[0], lvMax);
      for (auto limit : limits) {
        DrawLimit(limit, c);
      }
    }
    DrawRunBoudaries(runBoundaries, c);
  }
}

//----------------------------------------------------------------------------
bool ContainsAKey(std::string data, const std::set<std::string>& Keys)
{
  /// check if the data contains one of the keys

  auto itKey = std::find_if(Keys.begin(), Keys.end(), [&data](const auto& key) {
    return data.find(key) != data.npos;
  });

  return itKey != Keys.end();
}

//----------------------------------------------------------------------------
bool isValid(std::string alias)
{
  /// check if the alias is a valid (part of a) DCS alias

  static const std::vector<std::string> aliases =
    mch::dcs::aliases({mch::dcs::MeasurementType::HV_V,
                       mch::dcs::MeasurementType::LV_V_FEE_ANALOG,
                       mch::dcs::MeasurementType::LV_V_FEE_DIGITAL,
                       mch::dcs::MeasurementType::LV_V_SOLAR});

  auto itAlias = std::find_if(aliases.begin(), aliases.end(), [&alias](const auto& a) {
    return a.find(alias) != a.npos;
  });

  return itAlias != aliases.end();
}

//----------------------------------------------------------------------------
void ScanWhat(std::string what, std::string& path, bool& scanHV, bool& scanAll, std::set<std::string>& aliases)
{
  /// get what to scan and where

  static const std::set<std::string> hvKeys{"HV", "Quad", "Slat"};
  static const std::set<std::string> lvKeys{"LV", "Group", "an", "di", "Sol"};

  // HV or LV ?
  path = "";
  scanHV = false;
  if (ContainsAKey(what, hvKeys)) {
    path = "MCH/Calib/HV";
    scanHV = true;
  }
  if (ContainsAKey(what, lvKeys)) {
    if (scanHV) {
      printf("error: cannot scan HV and LV channels at the same time\n");
      exit(1);
    }
    path = "MCH/Calib/LV";
  }
  if (path.empty()) {
    printf("error: no valid HV or LV channel to scan\n");
    exit(1);
  }

  // everything or specific aliases ?
  if (what.find(scanHV ? "HV" : "LV") != what.npos) {
    scanAll = true;
    aliases.clear();
  } else {
    scanAll = false;
    std::istringstream input(what);
    for (std::string alias; std::getline(input, alias, ',');) {
      if (isValid(alias)) {
        aliases.insert(alias);
      } else {
        printf("error: \"%s\" invalid (part of) HV or LV alias\n", alias.c_str());
        exit(1);
      }
    }
  }
}

//----------------------------------------------------------------------------
std::set<int> GetRuns(std::string runList)
{
  /// read the runList from an ASCII file, or a comma separated run list, or a single run

  std::set<int> runs{};

  auto isNumber = [](std::string val) { return !val.empty() && val.find_first_not_of("0123456789") == val.npos; };

  if (isNumber(runList)) {

    runs.insert(std::stoi(runList));

  } else if (runList.find(",") != runList.npos) {

    std::istringstream input(runList);
    for (std::string run; std::getline(input, run, ',');) {
      if (isNumber(run)) {
        runs.insert(std::stoi(run));
      }
    }

  } else {

    std::ifstream input(runList);
    if (input.is_open()) {
      for (std::string run; std::getline(input, run);) {
        if (isNumber(run)) {
          runs.insert(std::stoi(run));
        }
      }
    }
  }

  return runs;
}

//----------------------------------------------------------------------------
RBMAP GetRunBoundaries(ccdb::CcdbApi const& api, std::string runList)
{
  /// return the SOR / EOR time stamps for every runs in the list

  RBMAP runBoundaries{};

  auto runs = GetRuns(runList);

  for (auto run : runs) {
    auto boundaries = ccdb::CCDBManagerInstance::getRunDuration(api, run);
    runBoundaries.emplace(run, boundaries);
  }

  return runBoundaries;
}

//----------------------------------------------------------------------------
void CheckRunBoundaries(const RBMAP& runBoundaries)
{
  /// check the consistency of the run time boundaries

  if (runBoundaries.empty()) {
    printf("error: no run found from the list\n");
    exit(1);
  }

  bool error = false;
  int previousRun = 0;
  uint64_t endOfPreviousRun = 0;

  for (const auto& [run, boundaries] : runBoundaries) {
    if (boundaries.second <= boundaries.first) {
      printf("error: run %d EOR <= SOR: %llu - %llu (%s - %s)\n",
             run, boundaries.first, boundaries.second,
             GetTime(boundaries.first).c_str(), GetTime(boundaries.second).c_str());
      error = true;
    }
    if (boundaries.first <= endOfPreviousRun) {
      printf("error: SOR run %d <= EOR run %d: %llu (%s) <= %llu (%s)\n",
             run, previousRun, boundaries.first, GetTime(boundaries.first).c_str(),
             endOfPreviousRun, GetTime(endOfPreviousRun).c_str());
      error = true;
    }
    previousRun = run;
    endOfPreviousRun = boundaries.second;
  }

  if (error) {
    exit(1);
  }
}

//----------------------------------------------------------------------------
void PrintRunBoundaries(const RBMAP& runBoundaries)
{
  /// print the list of runs with their time boundaries

  printf("\nlist of runs with their boundaries:\n");
  printf("------------------------------------\n");

  for (const auto& [run, boundaries] : runBoundaries) {
    printf("%d: %llu - %llu (%s - %s)\n", run, boundaries.first, boundaries.second,
           GetTime(boundaries.first).c_str(), GetTime(boundaries.second).c_str());
  }

  printf("------------------------------------\n");
}

//----------------------------------------------------------------------------
void DrawRunBoudaries(const RBMAP& runBoundaries, TCanvas* c)
{
  /// draw the run time boundaries

  c->cd();

  for (const auto& [run, boundaries] : runBoundaries) {

    TLine* startRunLine = new TLine(MSToS(boundaries.first), c->GetUymin(), MSToS(boundaries.first), c->GetUymax());
    startRunLine->SetUniqueID(run);
    startRunLine->SetLineColor(4);
    startRunLine->SetLineWidth(1);
    startRunLine->Draw();

    TLine* endRunLine = new TLine(MSToS(boundaries.second), c->GetUymin(), MSToS(boundaries.second), c->GetUymax());
    endRunLine->SetUniqueID(run);
    endRunLine->SetLineColor(2);
    endRunLine->SetLineWidth(1);
    endRunLine->Draw();
  }
}

//----------------------------------------------------------------------------
DPBMAP GetDPBoundaries(ccdb::CcdbApi const& api, std::string what, uint64_t tStart, uint64_t tStop)
{
  /// get the time boundaries of every HV/LV files found in the time range

  // add extra margin (ms) of ± 1 min to the creation time, which occurs every 30 min
  static const uint64_t timeMarging[2] = {60000, 1860000};

  std::istringstream fileInfo(api.list(what.c_str(), false, "text/plain",
                                       tStop + timeMarging[1], tStart - timeMarging[0]));

  DPBMAP dpBoundaries{};
  std::string dummy{};
  uint64_t begin = 0;
  uint64_t end = 0;

  for (std::string line; std::getline(fileInfo, line);) {
    if (line.find("Validity:") == 0) {
      std::istringstream in(line);
      in >> dummy >> begin >> dummy >> end;
      dpBoundaries.emplace(begin, end);
    }
  }

  if (dpBoundaries.empty()) {
    printf("\e[0;31merror: no file found in %s in time range %llu - %llu (%s - %s) --> use the default one\e[0m\n",
           what.c_str(), tStart, tStop, GetTime(tStart).c_str(), GetTime(tStop).c_str());
    dpBoundaries.emplace(1, 9999999999999);
  }

  return dpBoundaries;
}

//----------------------------------------------------------------------------
void CheckDPBoundaries(const DPBMAP& dpBoundaries, bool scanHV, uint64_t tStart, uint64_t tStop)
{
  /// check the consistency of HV/LV file time boundaries

  bool error = false;

  if (dpBoundaries.begin()->first > tStart) {
    printf("error: the beginning of the time range is not covered: %llu > %llu (%s > %s)\n",
           dpBoundaries.begin()->first, tStart,
           GetTime(dpBoundaries.begin()->first).c_str(), GetTime(tStart).c_str());
    error = true;
  }
  if (dpBoundaries.rbegin()->second < tStop) {
    printf("error: the end of the time range is not covered: %llu < %llu (%s < %s)\n",
           dpBoundaries.rbegin()->second, tStop,
           GetTime(dpBoundaries.rbegin()->second).c_str(), GetTime(tStop).c_str());
    error = true;
  }

  uint64_t previousTStop = dpBoundaries.begin()->first;
  for (auto [tStart, tStop] : dpBoundaries) {
    if (tStop <= tStart) {
      printf("error: EOF <= SOF: %llu - %llu (%s - %s)\n",
             tStart, tStop, GetTime(tStart).c_str(), GetTime(tStop).c_str());
      error = true;
    }
    if (tStart != previousTStop) {
      printf("error: end of %s file != start of next %s file: %llu (%s) != %llu (%s))\n",
             scanHV ? "HV" : "LV", scanHV ? "HV" : "LV",
             previousTStop, GetTime(previousTStop).c_str(), tStart, GetTime(tStart).c_str());
      error = true;
    }
    previousTStop = tStop;
  }

  if (error) {
    exit(1);
  }
}

//----------------------------------------------------------------------------
void PrintDPBoundaries(const DPBMAP& dpBoundaries, bool scanHV)
{
  /// print the time boundaries of every HV/LV files found in the full time range

  printf("\nlist of %s file time boundaries:\n", scanHV ? "HV" : "LV");
  printf("------------------------------------\n");

  for (auto [tStart, tStop] : dpBoundaries) {
    printf("%llu - %llu (%s - %s)\n", tStart, tStop, GetTime(tStart).c_str(), GetTime(tStop).c_str());
  }

  printf("------------------------------------\n");
}

//----------------------------------------------------------------------------
double GetLVLimit(std::string alias)
{
  /// return the LV limit for that channel

  if (alias.find("an") != alias.npos) {
    return lvLimits[0];
  } else if (alias.find("di") != alias.npos) {
    return lvLimits[1];
  }
  return lvLimits[2];
}

//----------------------------------------------------------------------------
void DrawLimit(double limit, TCanvas* c)
{
  /// draw the HV/LV limit for the displayed chamber

  c->cd();

  TLine* l = new TLine(c->GetUxmin(), limit, c->GetUxmax(), limit);
  l->SetLineColor(1);
  l->SetLineWidth(1);
  l->SetLineStyle(2);
  l->Draw();
}

//----------------------------------------------------------------------------
uint64_t MSToS(uint64_t ts)
{
  /// convert the time stamp from ms to s

  return (ts + 500) / 1000;
}

//----------------------------------------------------------------------------
std::string GetTime(uint64_t ts)
{
  /// convert the time stamp (ms) to local time

  time_t t = MSToS(ts);

  std::string time = std::ctime(&t);
  time.pop_back(); // remove trailing \n

  return time;
}

//----------------------------------------------------------------------------
std::string GetDuration(uint64_t tStart, uint64_t tStop)
{
  /// get the duration (dd hh:mm:ss) between the two time stamps (ms)

  auto dt = MSToS(tStop - tStart);
  auto s = dt % 60;
  auto m = (dt / 60) % 60;
  auto h = (dt / 3600) % 24;
  auto d = dt / 86400;

  return fmt::format("{:02}d {:02}:{:02}:{:02}", d, h, m, s);
}

//----------------------------------------------------------------------------
double GetValue(DPVAL dp)
{
  /// return the value of this data point

  union Converter {
    uint64_t raw_data;
    double value;
  } converter;

  converter.raw_data = dp.payload_pt1;

  return converter.value;
}

//----------------------------------------------------------------------------
std::string GetDE(std::string alias)
{
  /// for DCS HV alias: return the corresponding DE (and sector)
  /// for DCS LV alias: return an empty string

  auto de = mch::dcs::aliasToDetElemId(alias);

  if (de) {
    return (mch::dcs::isQuadrant(mch::dcs::aliasToChamber(alias)))
             ? fmt::format("DE{}-{}", *de, mch::dcs::aliasToNumber(alias) % 10)
             : fmt::format("DE{}", *de);
  }

  return "";
}

//----------------------------------------------------------------------------
void FillDataPoints(const std::vector<DPVAL>& dps, std::map<uint64_t, double>& dps2,
                    uint64_t tMin, uint64_t tMax, int warningLevel)
{
  /// fill the map of data points

  if (dps.empty()) {
    printf("error: the file does not contain any data point\n");
    exit(1);
  }

  auto itDP = dps.begin();
  auto ts = itDP->get_epoch_time();
  std::string header = "warning:";
  std::string color = (ts < tMin - 2000 || ts > tMin + 2000) ? "\e[0;31m" : "\e[0;34m";
  bool printWarning = warningLevel > 1 || (warningLevel == 1 && color == "\e[0;31m");

  // check if the first data point is a copy of the last one from previous file
  if (!dps2.empty()) {
    auto previousTS = dps2.rbegin()->first;
    if (ts != previousTS || GetValue(*itDP) != dps2.rbegin()->second) {
      if (ts <= previousTS) {
        printf("error: wrong data point order (%llu <= %llu)\n", ts, previousTS);
        exit(1);
      }
      if (printWarning) {
        printf("%s%s missing the previous data point (dt = %s%llu ms)", color.c_str(), header.c_str(),
               (previousTS < tMin) ? "-" : "+", (previousTS < tMin) ? tMin - previousTS : previousTS - tMin);
        if (ts <= tMin) {
          printf(" but get one at dt = -%llu ms\e[0m\n", tMin - ts);
        } else {
          printf("\e[0m\n");
        }
        header = "        ";
      }
    }
  }

  // add the first data point (should be before the start of validity of the file)
  if (ts >= tMax) {
    printf("error: first data point exceeding file validity range (dt = +%llu ms)\n", ts - tMax);
    exit(1);
  } else if (ts > tMin && printWarning) {
    printf("%s%s missing data point prior file start of validity (dt = +%llu ms)\e[0m\n",
           color.c_str(), header.c_str(), ts - tMin);
    header = "        ";
  }
  dps2.emplace(ts, GetValue(*itDP));

  // add other data points (should be within the validity range of the file)
  auto previousTS = ts;
  for (++itDP; itDP < dps.end(); ++itDP) {
    ts = itDP->get_epoch_time();
    if (ts <= previousTS) {
      printf("error: wrong data point order (%llu <= %llu)\n", ts, previousTS);
      exit(1);
    }
    if (ts < tMin && (warningLevel > 1 || (warningLevel == 1 && ts < tMin - 2000))) {
      printf("%s%s data point outside of file validity range (dt = -%llu ms)\e[0m\n",
             (ts < tMin - 2000) ? "\e[0;31m" : "\e[0;34m", header.c_str(), tMin - ts);
    } else if (ts >= tMax && warningLevel >= 1) {
      printf("\e[0;31m%s data point outside of file validity range (dt = +%llu ms)\e[0m\n",
             header.c_str(), ts - tMax);
    }
    dps2.emplace(ts, GetValue(*itDP));
    previousTS = ts;
  }
}

//----------------------------------------------------------------------------
void SelectDataPoints(DPMAP2 dpsMapsPerCh[10], uint64_t tStart, uint64_t tStop)
{
  /// remove the data points outside of the given time range and, if needed,
  /// add a data point at the boundaries with value equal to the preceding one

  for (int ch = 0; ch < 10; ++ch) {
    for (auto& [alias, dps] : dpsMapsPerCh[ch]) {

      // get the first data point in the time range, remove the previous ones
      // and add a data point with value equal to the preceding one if it exits
      // or to this one otherwise
      auto itFirst = dps.lower_bound(tStart);
      if (itFirst != dps.begin()) {
        double previousVal = std::prev(itFirst)->second;
        for (auto it = dps.begin(); it != itFirst;) {
          it = dps.erase(it);
        }
        dps.emplace(tStart, previousVal);
      } else if (itFirst->first != tStart) {
        if (itFirst->first > tStop) {
          printf("error (%s): all data points are posterior to the end of the time range\n", alias.c_str());
        } else {
          printf("error (%s): first data point is posterior to the beginning of the time range\n", alias.c_str());
        }
        dps.emplace(tStart, itFirst->second);
      }

      // get the first data point exceeding the time range, remove it and the next ones
      // and add a data point with value equal to the preceding one if needed
      auto itLast = dps.upper_bound(tStop);
      double previousVal = std::prev(itLast)->second;
      for (auto it = itLast; it != dps.end();) {
        it = dps.erase(it);
      }
      dps.emplace(tStop, previousVal);
    }
  }
}

//----------------------------------------------------------------------------
void PrintDataPoints(const DPMAP2 dpsMapsPerCh[10], bool scanHV, bool all)
{
  /// print all the registered data points

  const auto format1 = fmt::format("  %llu (%s): {} V\n", scanHV ? hvFormat : lvFormat);
  const auto format2 = fmt::format(": %llu (%s): {} V -- %llu (%s): {} V\n",
                                   scanHV ? hvFormat : lvFormat, scanHV ? hvFormat : lvFormat);

  for (int ch = 0; ch < 10; ++ch) {

    printf("\n------------ chamber %d ------------\n", ch + 1);

    for (const auto& [alias, dps] : dpsMapsPerCh[ch]) {

      printf("- %s: %lu values", alias.c_str(), dps.size());

      if (all) {

        printf("\n");
        for (const auto& [ts, val] : dps) {
          printf(format1.c_str(), ts, GetTime(ts).c_str(), val);
        }

      } else if (!dps.empty()) {

        const auto firstdt = dps.begin();
        const auto lastdt = dps.rbegin();
        printf(format2.c_str(),
               firstdt->first, GetTime(firstdt->first).c_str(), firstdt->second,
               lastdt->first, GetTime(lastdt->first).c_str(), lastdt->second);

      } else {
        printf("\n");
      }
    }
  }
}

//----------------------------------------------------------------------------
TGraph* MapToGraph(std::string alias, const std::map<uint64_t, double>& dps)
{
  /// create a graph for the DCS channel and add the data points

  TGraph* g = new TGraph(dps.size());

  auto pos = alias.find(".");
  auto shortAlias = alias.substr(0, pos);
  auto de = GetDE(alias);
  auto title = de.empty() ? fmt::format("{}", shortAlias.c_str())
                          : fmt::format("{} ({})", de.c_str(), shortAlias.c_str());
  g->SetNameTitle(alias.c_str(), title.c_str());

  int i(0);
  for (auto [ts, val] : dps) {
    g->SetPoint(i, MSToS(ts), val);
    ++i;
  }

  g->SetMarkerSize(1.5);
  g->SetMarkerStyle(2);
  g->SetLineStyle(2);

  return g;
}

//----------------------------------------------------------------------------
TCanvas* DrawDataPoints(TMultiGraph* mg, double min, double max)
{
  /// display the data points of the given chamber

  TCanvas* c = new TCanvas(mg->GetName(), mg->GetHistogram()->GetTitle(), 1500, 900);

  mg->Draw("A plc pmc");
  mg->SetMinimum(min);
  mg->SetMaximum(max);
  mg->GetXaxis()->SetTimeDisplay(1);
  mg->GetXaxis()->SetTimeFormat("%d/%m %H:%M");
  mg->GetXaxis()->SetTimeOffset(0, "local");
  mg->GetXaxis()->SetNdivisions(21010);

  c->BuildLegend();
  c->Update();

  return c;
}

//----------------------------------------------------------------------------
void FindIssues(const std::map<uint64_t, double>& dps, double limit, ISSUELIST& issues)
{
  /// return the list of HV/LV issues (time range, min value, mean value) for each DCS channel

  uint64_t tStart(0);
  double min(0.);
  double mean(0.);
  uint64_t prevTS(0);
  double prevVal(-1.);

  for (auto [ts, val] : dps) {

    if (val < limit) {

      if (tStart == 0) {

        // start a new issue...
        tStart = ts;
        min = val;
        mean = 0.;
        prevTS = ts;
        prevVal = val;

      } else {

        // ... or complement the current one
        min = std::min(min, val);
        mean += prevVal * (ts - prevTS);
        prevTS = ts;
        prevVal = val;
      }

    } else if (tStart > 0) {

      // complete the current issue, if any, and register it
      mean += prevVal * (ts - prevTS);
      mean /= (ts - tStart);
      issues.emplace_back(tStart, ts, min, mean, "");
      tStart = 0;
    }
  }

  // complete the last issue, if any and its duration is != 0, and register it
  if (tStart > 0 && prevTS != tStart) {
    mean /= (prevTS - tStart);
    issues.emplace_back(tStart, prevTS, min, mean, "");
  }
}

//----------------------------------------------------------------------------
void FillO2Issues(const std::vector<mch::HVStatusCreator::TimeRange>& o2issues, ISSUELIST& issues,
                  uint64_t tMin, uint64_t tMax)
{
  /// fill the list of issues from O2 (extend the previous one and/or create new ones)

  // the list must not be empty
  if (o2issues.empty()) {
    printf("error: O2 returns an empty list of issues\n");
    exit(1);
  }

  for (auto itIssue = o2issues.begin(); itIssue != o2issues.end(); ++itIssue) {

    // exclude issues fully outside of the DP file boudaries
    if (itIssue->end <= tMin || itIssue->begin >= tMax) {
      printf("\e[0;35mwarning: skipping O2 issue outside of file boundaries (%llu - %llu)\e[0m\n",
             itIssue->begin, itIssue->end);
      continue;
    }

    // only the first issue could in principle extend before the start of the DP file, to O
    if (itIssue->begin < tMin - mch::StatusMapCreatorParam::Instance().timeMargin &&
        (itIssue != o2issues.begin() || itIssue->begin != 0)) {
      printf("\e[0;35mwarning: O2 returns an issue with uncommon start time (%llu < %llu)\e[0m\n",
             itIssue->begin, tMin - mch::StatusMapCreatorParam::Instance().timeMargin);
    }

    // only the last issue could in principle extend beyond the end of the DP file, to infinity
    if (itIssue->end >= tMax + mch::StatusMapCreatorParam::Instance().timeMargin &&
        (itIssue != std::prev(o2issues.end()) || itIssue->end != std::numeric_limits<uint64_t>::max())) {
      printf("\e[0;35mwarning: O2 returns an issue with uncommon end time (%llu >= %llu)\e[0m\n",
             itIssue->end, tMax + mch::StatusMapCreatorParam::Instance().timeMargin);
    }

    // extend the last issue in case of continuity accross the DP files or add a new one,
    // restricting their time range within the DP file boundaries
    if (itIssue->begin <= tMin && !issues.empty() && std::get<1>(issues.back()) == tMin) {
      std::get<1>(issues.back()) = std::min(itIssue->end, tMax);
    } else {
      issues.emplace_back(std::max(itIssue->begin, tMin), std::min(itIssue->end, tMax), 0., 0., "");
    }
  }
}

//----------------------------------------------------------------------------
void SelectIssues(ISSUEMAP issuesPerCh[10], const RBMAP& runBoundaries, uint64_t minDuration)
{
  /// select HV/LV issues of a minimum duration (ms) occurring during runs

  for (int ch = 0; ch < 10; ++ch) {
    for (auto& issues : issuesPerCh[ch]) {
      for (auto itIssue = issues.second.begin(); itIssue != issues.second.end();) {

        auto tStart = std::get<0>(*itIssue);
        auto tStop = std::get<1>(*itIssue);

        if (tStop - tStart < minDuration) {

          itIssue = issues.second.erase(itIssue);

        } else {

          auto runs = FindAffectedRuns(runBoundaries, tStart, tStop);

          if (runs.empty()) {

            itIssue = issues.second.erase(itIssue);

          } else {

            std::get<4>(*itIssue) = runs;
            ++itIssue;
          }
        }
      }
    }
  }
}

//----------------------------------------------------------------------------
void SelectO2Issues(ISSUEMAP issuesPerCh[10], const RBMAP& runBoundaries)
{
  /// select HV issues from O2 algorithm occurring during runs
  /// and restrict the range of issues to the run range

  for (int ch = 0; ch < 10; ++ch) {
    for (auto& issues : issuesPerCh[ch]) {
      for (auto itIssue = issues.second.begin(); itIssue != issues.second.end();) {

        auto& tStart = std::get<0>(*itIssue);
        auto& tStop = std::get<1>(*itIssue);

        auto runs = FindAffectedRuns(runBoundaries, tStart, tStop);

        if (runs.empty()) {

          itIssue = issues.second.erase(itIssue);

        } else {

          tStart = std::max(tStart, runBoundaries.begin()->second.first);
          tStop = std::min(tStop, runBoundaries.rbegin()->second.second);
          std::get<4>(*itIssue) = runs;
          ++itIssue;
        }
      }
    }
  }
}

//----------------------------------------------------------------------------
std::string FindAffectedRuns(const RBMAP& runBoundaries, uint64_t tStart, uint64_t tStop)
{
  /// return the list of affected runs in this time range

  std::string runs;

  for (const auto& [run, boundaries] : runBoundaries) {

    if (boundaries.second <= tStart) {
      continue;
    } else if (boundaries.first >= tStop) {
      break;
    }

    runs += fmt::format("{},", run);
  }

  if (!runs.empty()) {
    runs.pop_back();
  }

  return runs;
}

//----------------------------------------------------------------------------
bool eraseIssue(const ISSUE& issue, ISSUELIST& issues)
{
  /// find an issue with the same time range and associated run list and erase it
  /// return true in case of success

  auto itIssue = std::find_if(issues.begin(), issues.end(), [&issue](const auto& i) {
    return (std::get<0>(i) == std::get<0>(issue) &&
            std::get<1>(i) == std::get<1>(issue) &&
            std::get<4>(i) == std::get<4>(issue));
  });

  if (itIssue != issues.end()) {
    issues.erase(itIssue);
    return true;
  }

  return false;
}

//----------------------------------------------------------------------------
void PrintIssues(const ISSUEMAP issuesPerCh[10], const ISSUEMAP o2IssuesPerCh[10], bool scanHV)
{
  /// print all HV/LV issues

  // copy the issues so that we can modify them (i.e. add empty lists or delete issues after printing)
  ISSUEMAP issuesPerChCopy[10];
  ISSUEMAP o2IssuesPerChCopy[10];
  for (int ch = 0; ch < 10; ++ch) {
    issuesPerChCopy[ch] = issuesPerCh[ch];
    o2IssuesPerChCopy[ch] = o2IssuesPerCh[ch];
  }

  // make sure that all alias keys in the map o2IssuesPerChCopy are also in issuesPerChCopy in order to
  // simplify the loop over all issues from both algorithms and fix the order in which they are printed
  for (int ch = 0; ch < 10; ++ch) {
    for (const auto& [alias, o2Issues] : o2IssuesPerChCopy[ch]) {
      if (!o2Issues.empty()) {
        issuesPerChCopy[ch].try_emplace(alias, ISSUELIST{});
      }
    }
  }

  auto printHeader = [](std::string alias) {
    auto de = GetDE(alias);
    if (de.empty()) {
      printf("Problem found for %s:\n", alias.c_str());
    } else {
      printf("Problem found for %s (%s):\n", alias.c_str(), de.c_str());
    }
  };

  const auto format = fmt::format("%llu - %llu: %s (duration = %s, min = {} V, mean = {} V) --> run(s) %s\n",
                                  scanHV ? hvFormat : lvFormat, scanHV ? hvFormat : lvFormat);

  auto printIssue = [&format](ISSUE issue, std::string color) {
    const auto& [tStart, tStop, min, mean, runs] = issue;
    printf("%s", color.c_str());
    printf(format.c_str(), tStart, tStop,
           GetTime(tStart).c_str(), GetDuration(tStart, tStop).c_str(), min, mean, runs.c_str());
    printf("\e[0m");
  };

  if (scanHV) {
    printf("\n------ list of issues from \e[0;31mthis macro only\e[0m, \e[0;35mO2 only\e[0m, or \e[0;32mboth\e[0m ------\n");
  } else {
    printf("\n------ list of issues ------\n");
  }

  bool foundIssues = false;

  for (int ch = 0; ch < 10; ++ch) {
    for (const auto& [alias, issues] : issuesPerChCopy[ch]) {

      auto& o2Issues = o2IssuesPerChCopy[ch][alias];

      if (!issues.empty() || !o2Issues.empty()) {

        foundIssues = true;
        printHeader(alias);

        // print all issues found by this macro
        for (const auto& issue : issues) {
          // change color if the issue is not found by the O2 algorithm (only for HV)
          std::string color = (scanHV && !eraseIssue(issue, o2Issues)) ? "\e[0;31m" : "\e[0;32m";
          printIssue(issue, color);
        }

        // print other issues found by the O2 algorithm
        for (const auto& issue : o2Issues) {
          printIssue(issue, "\e[0;35m");
        }

        printf("----------------------------\n");
      }
    }
  }

  if (!foundIssues) {
    printf("----------------------------\n");
  }
}
