#include <algorithm>
#include <ctime>
#include <fstream>
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

using namespace o2;
using DPID = dcs::DataPointIdentifier;
using DPVAL = dcs::DataPointValue;
using DPMAP = std::unordered_map<DPID, std::vector<DPVAL>>;
using DPMAP2 = std::map<std::string, std::map<uint64_t, double>>;
using RBMAP = std::map<int, std::pair<uint64_t, uint64_t>>;
using HVBMAP = std::map<uint64_t, uint64_t>;
using BADHVLIST = std::vector<std::tuple<uint64_t, uint64_t, double, double, std::string>>;
using BADHVMAP = std::map<std::string, BADHVLIST>;

double yRange[2] = {-1., 1700.};
double hvLimits[10] = {1550., 1550., 1590., 1590., 1590., 1590., 1590., 1590., 1590., 1590.};

std::set<int> GetRuns(std::string runList);
RBMAP GetRunBoundaries(ccdb::CcdbApi const& api, std::string runList);
void PrintRunBoundaries(const RBMAP& runBoundaries);
void DrawRunBoudaries(const RBMAP& runBoundaries, TCanvas* c);
HVBMAP GetHVBoundaries(ccdb::CcdbApi const& api, uint64_t tStart, uint64_t tStop);
void PrintHVBoundaries(const HVBMAP& hvBoundaries);
void DrawHVLimits(int ch, TCanvas* c);
uint64_t MSToS(uint64_t ts);
std::string GetTime(uint64_t ts);
std::string GetDuration(uint64_t tStart, uint64_t tStop);
double GetHV(DPVAL dp);
std::string GetDE(std::string alias);
void SelectDataPoints(DPMAP2 dpsMapsPerCh[10], uint64_t tStart, uint64_t tStop);
void PrintDataPoints(const DPMAP2 dpsMapsPerCh[10], bool all);
TGraph* MapToGraph(std::string alias, const std::map<uint64_t, double>& dps);
TCanvas* DrawDataPoints(TMultiGraph* mg);
void FindHVIssues(const std::map<uint64_t, double>& dps, double hvLimit, BADHVLIST& hvIssues);
void SelectHVIssues(BADHVMAP hvIssuesPerCh[10], const RBMAP& runBoundaries, uint64_t minDuration);
std::string FindAffectedRuns(const RBMAP& runBoundaries, uint64_t tStart, uint64_t tStop);
void PrintHVIssues(const BADHVMAP hvIssuesPerCh[10]);

//----------------------------------------------------------------------------
void ScanHV(std::string runList, std::string aliases = "", int printLevel = 1)
{
  /// scan the HV of every sectors to check for issues
  /// runList can be an ASCII file with one run per line, or a comma separated run list, or a single run
  /// aliases is the list of DCS channels to be considered (empty = all)
  /// printLevel >= 1: print time stamps of runs and HV files
  /// printLevel >= 2: print the first and last data points of each selected channel
  /// printLevel >= 3: print all the data points of each selected channel

  gStyle->SetPalette(kVisibleSpectrum);

  ccdb::CcdbApi api;
  api.init("http://alice-ccdb.cern.ch");
 
  // get the SOR/EOR of every runs from the list, ordered in run number
  auto runBoundaries = GetRunBoundaries(api, runList);
  if (runBoundaries.empty()) {
    printf("no run found from the list\n");
    return;
  }
  if (printLevel > 0) {
    PrintRunBoundaries(runBoundaries);
  }

  // extract the time boundaries for each HV file in the full time range
  auto hvBoundaries = GetHVBoundaries(api, runBoundaries.begin()->second.first, runBoundaries.rbegin()->second.second);
  if (printLevel > 0) {
    PrintHVBoundaries(hvBoundaries);
  }

  // loop over the HV files and fill the lists of data points per chamber
  DPMAP2 dpsMapsPerCh[10];
  std::map<std::string, std::string> metadata;
  for (auto boundaries : hvBoundaries) {
    auto* dpMap = api.retrieveFromTFileAny<DPMAP>("MCH/Calib/HV", metadata, boundaries.first);
    for (const auto& [dpid, dps] : *dpMap) {
      std::string alias(dpid.get_alias());
      if ((aliases.empty() || aliases.find(alias) != aliases.npos) && alias.find(".vMon") != alias.npos) {
        int chamber = mch::dcs::toInt(mch::dcs::aliasToChamber(alias));
        auto& dps2 = dpsMapsPerCh[chamber][alias];
        for (const auto& dp : dps) {
          dps2.emplace(dp.get_epoch_time(), GetHV(dp));
        }
      }
    }
  }
  if (printLevel > 1) {
    PrintDataPoints(dpsMapsPerCh, printLevel > 2);
  }

  // select the data points in the time range
  SelectDataPoints(dpsMapsPerCh, runBoundaries.begin()->second.first, runBoundaries.rbegin()->second.second);
  if (printLevel > 1) {
    PrintDataPoints(dpsMapsPerCh, printLevel > 2);
  }

  // create and fill the graphs, and find HV issues
  BADHVMAP hvIssuesPerCh[10];
  TMultiGraph* mg[10];
  for (int ch = 0; ch < 10; ++ch) {
    mg[ch] = new TMultiGraph;
    mg[ch]->SetNameTitle(fmt::format("ch{}", ch).c_str(), fmt::format("chamber {};time;HV (V)", ch + 1).c_str());
    for (const auto& [alias, dps] : dpsMapsPerCh[ch]) {
      mg[ch]->Add(MapToGraph(alias, dps), "lp");
      FindHVIssues(dps, hvLimits[ch], hvIssuesPerCh[ch][alias]);
    }
  }

  // select HV issues of a minimum duration (ms) occurring during runs
  SelectHVIssues(hvIssuesPerCh, runBoundaries, 0);
  PrintHVIssues(hvIssuesPerCh);

  // display
  for (int ch = 0; ch < 10; ++ch) {
    TCanvas* c = DrawDataPoints(mg[ch]);
    DrawRunBoudaries(runBoundaries, c);
    DrawHVLimits(ch, c);
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
void PrintRunBoundaries(const RBMAP& runBoundaries)
{
  /// print the list of runs with their time boundaries

  printf("\nlist of runs with their boundaries:\n");
  printf("------------------------------------\n");

  for (const auto& [run, boundaries] : runBoundaries) {
    printf("%d: %lld - %lld (%s - %s)\n", run, boundaries.first, boundaries.second,
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

    TLine* startRunLine = new TLine(MSToS(boundaries.first), yRange[0], MSToS(boundaries.first), yRange[1]);
    startRunLine->SetUniqueID(run);
    startRunLine->SetLineColor(4);
    startRunLine->SetLineWidth(1);
    startRunLine->Draw();

    TLine* endRunLine = new TLine(MSToS(boundaries.second), yRange[0], MSToS(boundaries.second), yRange[1]);
    endRunLine->SetUniqueID(run);
    endRunLine->SetLineColor(2);
    endRunLine->SetLineWidth(1);
    endRunLine->Draw();
  }
}

//----------------------------------------------------------------------------
HVBMAP GetHVBoundaries(ccdb::CcdbApi const& api, uint64_t tStart, uint64_t tStop)
{
  /// get the time boundaries of every HV files found in the time range

  // add extra margin (ms) of ± 1 min to the creation time, which occurs every 30 min
  static const uint64_t timeMarging[2] = {60000, 1860000};

  HVBMAP hvBoundaries{};

  std::istringstream fileInfo(api.list("MCH/Calib/HV", false, "text/plain",
                                       tStop + timeMarging[1], tStart - timeMarging[0]));

  for (std::string line; std::getline(fileInfo, line);) {
    if (line.find("Validity:") == 0) {
      hvBoundaries.emplace(std::stoull(line.substr(10, 13)), std::stoull(line.substr(26, 13)));
    }
  }

  return hvBoundaries;
}

//----------------------------------------------------------------------------
void PrintHVBoundaries(const HVBMAP& hvBoundaries)
{
  /// print the time boundaries of every HV files found in the full time range

  printf("\nlist of HV file time boundaries:\n");
  printf("------------------------------------\n");

  for (auto [tStart, tStop] : hvBoundaries) {
    printf("%lld - %lld (%s - %s)\n", tStart, tStop, GetTime(tStart).c_str(), GetTime(tStop).c_str());
  }

  printf("------------------------------------\n");
}

//----------------------------------------------------------------------------
void DrawHVLimits(int ch, TCanvas* c)
{
  /// draw the HV limits

  c->cd();

  TLine* hvLimit = new TLine(c->GetUxmin(), hvLimits[ch], c->GetUxmax(), hvLimits[ch]);
  hvLimit->SetLineColor(1);
  hvLimit->SetLineWidth(1);
  hvLimit->SetLineStyle(2);
  hvLimit->Draw();
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
double GetHV(DPVAL dp)
{
  /// return the HV value of this data point

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
  /// get the DE (and sector) corresponding to the DCS alias

  auto de = mch::dcs::aliasToDetElemId(alias);

  return (mch::dcs::isQuadrant(mch::dcs::aliasToChamber(alias)))
           ? fmt::format("DE{}-{}", *de, mch::dcs::aliasToNumber(alias) % 10)
           : fmt::format("DE{}", *de);
}

//----------------------------------------------------------------------------
void SelectDataPoints(DPMAP2 dpsMapsPerCh[10], uint64_t tStart, uint64_t tStop)
{
  /// remove the data points outside of the given time range and, if needed,
  /// add a data point at the boundaries with HV equal to the preceding value

  for (int ch = 0; ch < 10; ++ch) {
    for (auto& [alias, dps] : dpsMapsPerCh[ch]) {

      // get the first data point in the time range, remove the previous ones
      // and add a data point with HV equal to the preceding value if needed
      auto itFirst = dps.lower_bound(tStart);
      if (itFirst != dps.begin()) {
        double previousHV = std::prev(itFirst)->second;
        for (auto it = dps.begin(); it != itFirst;) {
          it = dps.erase(it);
        }
        dps.emplace(tStart, previousHV);
      } else if (itFirst->first != tStart) {
        printf("error (%s): first data point is posterior to the beginning of the time range\n", alias.c_str());
      }

      // get the first data point exceeding the time range, remove it and the next ones
      // and add a data point with HV equal to the preceding value if needed
      auto itLast = dps.upper_bound(tStop);
      if (itLast != dps.begin()) {
        double previousHV = std::prev(itLast)->second;
        for (auto it = itLast; it != dps.end();) {
          it = dps.erase(it);
        }
        dps.emplace(tStop, previousHV);
      } else {
        printf("error (%s): all data points are posterior to the end of the time range\n", alias.c_str());
        dps.clear();
      }
    }
  }
}

//----------------------------------------------------------------------------
void PrintDataPoints(const DPMAP2 dpsMapsPerCh[10], bool all)
{
  /// print all the registered data points

  for (int ch = 0; ch < 10; ++ch) {

    printf("\n------------ chamber %d ------------\n", ch + 1);

    for (const auto& [alias, dps] : dpsMapsPerCh[ch]) {

      printf("- %s: %lu values", alias.c_str(), dps.size());

      if (all) {

        printf("\n");
        for (const auto& [ts, hv] : dps) {
          printf("  %lld (%s): %7.2f V\n", ts, GetTime(ts).c_str(), hv);
        }

      } else if (!dps.empty()) {

        const auto firstdt = dps.begin();
        const auto lastdt = dps.rbegin();
        printf(": %lld (%s): %7.2f V -- %lld (%s): %7.2f V\n",
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

  auto shortAlias = alias.substr(0, alias.size() - 12);
  auto title = fmt::format("{} ({})", GetDE(alias).c_str(), shortAlias.c_str());
  g->SetNameTitle(alias.c_str(), title.c_str());

  int i(0);
  for (auto [ts, hv] : dps) {
    g->SetPoint(i, MSToS(ts), hv);
    ++i;
  }

  g->SetMarkerSize(1.5);
  g->SetMarkerStyle(2);
  g->SetLineStyle(2);

  return g;
}

//----------------------------------------------------------------------------
TCanvas* DrawDataPoints(TMultiGraph* mg)
{
  /// display the data points of the given chamber

  TCanvas* c = new TCanvas(mg->GetName(), mg->GetHistogram()->GetTitle(), 1500, 900);

  mg->Draw("A plc pmc");
  mg->SetMinimum(yRange[0]);
  mg->SetMaximum(yRange[1]);
  mg->GetXaxis()->SetTimeDisplay(1);
  mg->GetXaxis()->SetTimeFormat("%d/%m %H:%M");
  mg->GetXaxis()->SetTimeOffset(0, "local");
  mg->GetXaxis()->SetNdivisions(21010);

  c->BuildLegend();
  c->Update();

  return c;
}

//----------------------------------------------------------------------------
void FindHVIssues(const std::map<uint64_t, double>& dps, double hvLimit, BADHVLIST& hvIssues)
{
  /// return the list of HV issues (time range, min HV, mean HV) for each DCS channel

  uint64_t tStart(0);
  double minHV(0.);
  double meanHV(0.);
  uint64_t prevTS(0);
  double prevHV(-1.);

  for (auto [ts, hv] : dps) {

    if (hv < hvLimit) {

      if (tStart == 0) {

        // start a new HV issue...
        tStart = ts;
        minHV = hv;
        meanHV = 0.;
        prevTS = ts;
        prevHV = hv;

      } else {

        // ... or complement the current one
        minHV = std::min(minHV, hv);
        meanHV += prevHV * (ts - prevTS);
        prevTS = ts;
        prevHV = hv;
      }

    } else if (tStart > 0) {

      // complete the current HV issue, if any, and register it
      meanHV += prevHV * (ts - prevTS);
      meanHV /= (ts - tStart);
      hvIssues.emplace_back(tStart, ts, minHV, meanHV, "");
      tStart = 0;
    }
  }

  // complete the last HV issue, if any and its duration is != 0, and register it
  if (tStart > 0 && prevTS != tStart) {
    meanHV /= (prevTS - tStart);
    hvIssues.emplace_back(tStart, prevTS, minHV, meanHV, "");
  }
}

//----------------------------------------------------------------------------
void SelectHVIssues(BADHVMAP hvIssuesPerCh[10], const RBMAP& runBoundaries, uint64_t minDuration)
{
  /// select HV issues of a minimum duration (ms) occurring during runs

  for (int ch = 0; ch < 10; ++ch) {
    for (auto& hvIssues : hvIssuesPerCh[ch]) {
      for (auto itIssue = hvIssues.second.begin(); itIssue != hvIssues.second.end();) {

        auto tStart = std::get<0>(*itIssue);
        auto tStop = std::get<1>(*itIssue);

        if (tStop - tStart < minDuration) {

          itIssue = hvIssues.second.erase(itIssue);

        } else {

          auto runs = FindAffectedRuns(runBoundaries, tStart, tStop);

          if (runs.empty()) {

            itIssue = hvIssues.second.erase(itIssue);

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
void PrintHVIssues(const BADHVMAP hvIssuesPerCh[10])
{
  /// print all HV issues

  printf("\n------ list of issues ------\n");

  bool foundIssues = false;
  for (int ch = 0; ch < 10; ++ch) {
    for (const auto& [alias, hvIssues] : hvIssuesPerCh[ch]) {

      if (!hvIssues.empty()) {

        foundIssues = true;
        printf("Problem found for %s (%s):\n", alias.c_str(), GetDE(alias).c_str());

        for (const auto& [tStart, tStop, minHV, meanHV, runs] : hvIssues) {
          printf("- %s (duration = %s, min = %7.2f V, mean = %7.2f V) --> run(s) %s\n",
                 GetTime(tStart).c_str(), GetDuration(tStart, tStop).c_str(), minHV, meanHV, runs.c_str());
        }

        printf("----------------------------\n");
      }
    }
  }

  if (!foundIssues) {
    printf("----------------------------\n");
  }
}
