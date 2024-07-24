#include <fstream>
#include <iostream>
#include <list>
#include <vector>
#include <algorithm>
#include <type_traits>
#include <tuple>
#include <string>

#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TMath.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TMatrixD.h>
#include <TDatabasePDG.h>
#include <Math/Vector4D.h>
#include <TGeoGlobalMagField.h>
// #include <TRandom.h>

#include "AliCDBManager.h"
#include "AliGRPManager.h"
#include "AliGRPObject.h"
#include "AliGeomManager.h"
#include "AliMUONCDB.h"

#include "CCDB/BasicCCDBManager.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "Field/MagneticField.h"

#include "DataFormatsMCH/ROFRecord.h"
#include "DataFormatsMCH/TrackMCH.h"
#include "DataFormatsMCH/Cluster.h"
#include "MCHBase/TrackBlock.h"
#include "MCHTracking/TrackParam.h"
#include "MCHTracking/Track.h"
#include "MCHTracking/TrackFitter.h"
#include "MCHTracking/TrackExtrap.h"
#include "ReconstructionDataFormats/TrackMCHMID.h"

#include "/Users/PILLOT/Work/Alice/Macros/Tracking/TrackMCHv1.h"
#include "/Users/PILLOT/Work/Alice/Macros/Tracking/TrackMCHv2.h"

using namespace std;
using namespace ROOT::Math;
using o2::dataformats::TrackMCHMID;
using o2::mch::Cluster;
using o2::mch::ROFRecord;
using o2::mch::Track;
using o2::mch::TrackExtrap;
using o2::mch::TrackFitter;
using o2::mch::TrackMCH;
using o2::mch::TrackMCHv1;
using o2::mch::TrackMCHv2;
using o2::mch::TrackParam;
using o2::mch::TrackParamStruct;

constexpr double pi() { return 3.14159265358979323846; }
static const double muMass = TDatabasePDG::Instance()->GetParticle("mu-")->Mass();
double chi2Max = 2. * 4. * 4.;
TrackFitter trackFitter{};
bool ocdbLoaded = false;
int run = 0;

//_________________________________________________________________________________________________
struct ClusterStructV1 {
  float x;             ///< cluster position along x
  float y;             ///< cluster position along y
  float z;             ///< cluster position along z
  float ex;            ///< cluster resolution along x
  float ey;            ///< cluster resolution along y
  uint32_t uid;        ///< cluster unique ID
};

//_________________________________________________________________________________________________
struct TrackAtVtxStruct {
  TrackParamStruct paramAtVertex{};
  double dca = 0.;
  double rAbs = 0.;
  int mchTrackIdx = 0;
};

//_________________________________________________________________________________________________
struct VertexStruct {
  double x;
  double y;
  double z;
};

//_________________________________________________________________________________________________
struct TrackStruct {
  PxPyPzMVector pxpypzm{};
  short sign = 0;
  double dca = 0.;
  double rAbs = 0.;
  double chi2 = 0.;
  double matchedChi2 = 0.;
  TrackParamStruct param{};
  TMatrixD cov{5, 5};
  TrackParamStruct paramAtMID{};
  TMatrixD covAtMID{5, 5};
  std::vector<Cluster> clusters{};
  bool matchFound = false;
  bool matchIdentical = false;
  bool connected = false;
  bool needExtrapToVtx = false;

  bool operator==(const TrackStruct& track) const
  {
    /// tracks are considered identical when all their clusters match within chi2Max
    if (this->clusters.size() != track.clusters.size()) {
      return false;
    }
    for (size_t iCl = 0; iCl != this->clusters.size(); ++iCl) {
      auto& cl1 = this->clusters[iCl];
      auto& cl2 = track.clusters[iCl];
      if (cl1.getDEId() != cl2.getDEId()) {
        return false;
      }
      double dx = cl1.getX() - cl2.getX();
      double dy = cl1.getY() - cl2.getY();
      double chi2 = dx * dx / (cl1.getEx2() + cl2.getEx2()) + dy * dy / (cl1.getEy2() + cl2.getEy2());
      if (chi2 > chi2Max) {
        return false;
      }
    }
    return true;
  }

  bool match(const TrackStruct& track) const
  {
    /// Try to match this track with the given track. Matching conditions:
    /// - more than 50% of clusters from one of the two tracks matched with clusters from the other
    /// - at least 1 cluster matched before and 1 cluster matched after the dipole

    size_t nMatchClusters(0);
    bool matchCluster[10] = {false, false, false, false, false, false, false, false, false, false};

    for (const auto& cl1 : this->clusters) {
      for(const auto& cl2 : track.clusters) {
        if (cl1.getDEId() == cl2.getDEId()) {
          double dx = cl1.getX() - cl2.getX();
          double dy = cl1.getY() - cl2.getY();
          double chi2 = dx * dx / (cl1.getEx2() + cl2.getEx2()) + dy * dy / (cl1.getEy2() + cl2.getEy2());
          if (chi2 <= chi2Max) {
            matchCluster[cl1.getChamberId()] = true;
            ++nMatchClusters;
            break;
          }
        }
      }
    }

    return ((matchCluster[0] || matchCluster[1] || matchCluster[2] || matchCluster[3]) &&
            (matchCluster[6] || matchCluster[7] || matchCluster[8] || matchCluster[9]) &&
            (2 * nMatchClusters > this->clusters.size() || 2 * nMatchClusters > track.clusters.size()));
  }

  void printClusterDifferences(const TrackStruct& track) const
  {
    /// print the DE where clusters differ between the 2 tracks
    printf("additional clusters: track1 (%f) : ", getChi2N2());
    for (const auto& cluster1 : clusters) {
      bool found(false);
      for (const auto& cluster2 : track.clusters) {
        if (cluster1.uid == cluster2.uid) {
          found = true;
          break;
        }
      }
      if (found) {
        if (cluster1.getEx() > 1. || cluster1.getEy() > 1.) {
          printf("%d (mono), ", cluster1.getDEId());
        } else {
          printf("%d, ", cluster1.getDEId());
        }
      } else {
        if (cluster1.getEx() > 1. || cluster1.getEy() > 1.) {
          printf("\e[91m%d (mono)\e[0m, ", cluster1.getDEId());
        } else {
          printf("\e[91m%d\e[0m, ", cluster1.getDEId());
        }
      }
    }
    printf("track2 (%f) : ", track.getChi2N2());
    for (const auto& cluster2 : track.clusters) {
      bool found(false);
      for (const auto& cluster1 : clusters) {
        if (cluster1.uid == cluster2.uid) {
          found = true;
          break;
        }
      }
      if (found) {
        if (cluster2.getEx() > 1. || cluster2.getEy() > 1.) {
          printf("%d (mono), ", cluster2.getDEId());
        } else {
          printf("%d, ", cluster2.getDEId());
        }
      } else {
        if (cluster2.getEx() > 1. || cluster2.getEy() > 1.) {
          printf("\e[91m%d (mono)\e[0m, ", cluster2.getDEId());
        } else {
          printf("\e[91m%d\e[0m, ", cluster2.getDEId());
        }
      }
    }
    printf("\n");
  }

  int getNClustersInCommon(const TrackStruct& track) const
  {
    /// return the number of clusters in common between this track and the one given as parameter
    int nClustersInCommon(0);
    for (const auto& cluster1 : clusters) {
      for (const auto& cluster2 : track.clusters) {
        if (cluster1.uid == cluster2.uid) {
          ++nClustersInCommon;
          break;
        }
      }
    }
    return nClustersInCommon;
  }

  bool missClusterOnSt345(const TrackStruct& track) const
  {
    /// return true if the given track misses clusters from this track on station 3, 4 and 5
    for (const auto& cluster1 : clusters) {
      if (cluster1.getChamberId() < 4) {
        continue;
      }
      bool found(false);
      for (const auto& cluster2 : track.clusters) {
        if (cluster2.getChamberId() < 4) {
          continue;
        }
        if (cluster1.uid == cluster2.uid) {
          found = true;
          break;
        }
      }
      if (!found) {
        return true;
      }
    }
    return false;
  }

  bool isValid() const
  {
    /// check if the track passes *all* the tracking criteria, that are:
    /// 1 chamber fired per station and 3 chambers fired in stations 4 & 5
    int nChFired[5] = {0, 0, 0, 0, 0};
    int previousCh(-1);
    for (const auto& cluster : clusters) {
      int chId = cluster.getChamberId();
      if (chId != previousCh) {
        nChFired[chId / 2]++;
        previousCh = chId;
      }
    }
    if (nChFired[0] == 0 || nChFired[1] == 0 || nChFired[2] == 0 || nChFired[3] + nChFired[4] < 3) {
      return false;
    }
    return true;
  }

  bool containsMonoCathodClusters() const
  {
    /// return true if the track contains mono-cathod cluster(s)
    for (const auto& cluster : clusters) {
      if (cluster.getEx() > 1. || cluster.getEy() > 1.) {
        return true;
      }
    }
    return false;
  }

  int getNChamberFired() const
  {
    /// return the number of chamber fired
    int nCh(0);
    int previousCh(-1);
    for (const auto& cluster : clusters) {
      int chId = cluster.getChamberId();
      if (chId != previousCh) {
        nCh++;
        previousCh = chId;
      }
    }
    return nCh;
  }

  int getNdf() const
  {
    /// return the number of degrees of freedom
    return 2 * clusters.size() - 5;
  }

  int getNdf2() const
  {
    /// return the number of degrees of freedom - 1
    return 2 * clusters.size() - 6;
  }

  int getNdf3() const
  {
    /// return the number of degrees of freedom accounting for missing information from mono-cathod clusters
    int ndf(-5);
    for (const auto& cluster : clusters) {
      if (cluster.getEx() < 1.) {
        ++ndf;
      }
      if (cluster.getEy() < 1.) {
        ++ndf;
      }
    }
    return ndf;
  }

  double getChi2N() const
  {
    /// return the normalized chi2
    return chi2 / getNdf();
  }

  double getChi2N2() const
  {
    /// return the normalized chi2
    return chi2 / getNdf2();
  }

  double getChi2N3() const
  {
    /// return the normalized chi2 (accounting for missing info)
    return chi2 / getNdf3();
  }

  bool isBetter1(const TrackStruct& track) const
  {
    /// Return true if this track is better than the one given as parameter
    /// It is better if it has more clusters or a better chi2 in case of equality
    int nCl1 = this->clusters.size();
    int nCl2 = track.clusters.size();
    return ((nCl1 > nCl2) || ((nCl1 == nCl2) && (this->chi2 < track.chi2)));
  }

  bool isBetter2(const TrackStruct& track) const
  {
    /// Return true if this track is better than the one given as parameter
    /// It is better if it has a better normalized chi2
    return (this->getChi2N() < track.getChi2N());
  }

  bool isBetter22(const TrackStruct& track) const
  {
    /// Return true if this track is better than the one given as parameter
    /// It is better if it has a better normalized chi2
    return (this->getChi2N2() < track.getChi2N2());
  }

  bool isBetter23(const TrackStruct& track) const
  {
    /// Return true if this track is better than the one given as parameter
    /// It is better if it has a better normalized chi2 (accounting for missing info)
    return (this->getChi2N3() < track.getChi2N3());
  }

  bool isBetter3(const TrackStruct& track) const
  {
    /// Return true if this track is better than the one given as parameter
    /// It is better if it has a higher probability of chi2
    return (TMath::Prob(this->chi2, this->getNdf()) > TMath::Prob(track.chi2, track.getNdf()));
  }

  bool isBetter4(const TrackStruct& track) const
  {
    /// Return true if this track is better than the one given as parameter
    /// It is better if it has more chambers fired or a better chi2/ndf in case of equality
    int nCh1 = getNChamberFired();
    int nCh2 = track.getNChamberFired();
    return ((nCh1 > nCh2) || ((nCh1 == nCh2) && (getChi2N() < track.getChi2N())));
  }

  bool isBetter42(const TrackStruct& track) const
  {
    /// Return true if this track is better than the one given as parameter
    /// It is better if it has more chambers fired or a better chi2/ndf in case of equality
    int nCh1 = getNChamberFired();
    int nCh2 = track.getNChamberFired();
    return ((nCh1 > nCh2) || ((nCh1 == nCh2) && (getChi2N2() < track.getChi2N2())));
  }

  bool isBetter43(const TrackStruct& track) const
  {
    /// Return true if this track is better than the one given as parameter
    /// It is better if it has more chambers fired or a better chi2/ndf (accounting for missing info) in case of equality
    int nCh1 = getNChamberFired();
    int nCh2 = track.getNChamberFired();
    return ((nCh1 > nCh2) || ((nCh1 == nCh2) && (getChi2N3() < track.getChi2N3())));
  }
};

//_________________________________________________________________________________________________
int ReadNextEvent(ifstream& inFile, std::vector<VertexStruct>& vertices);
int ReadNextEvent(ifstream& inFile, int version, std::list<TrackStruct>& tracks);
void ReadTrack(ifstream& inFile, TrackStruct& track, int version);
template <class T>
void ReadNextEventV5(ifstream& inFile, int& event, std::list<TrackStruct>& tracks);
std::tuple<TFile*, TTreeReader*> LoadData(const char* fileName, const char* treeName);
template <typename T>
void FillTrack(TrackStruct& track, const TrackAtVtxStruct* trackAtVtx, const T* mchTrack, std::vector<Cluster>& clusters);
void UpdateTrack(TrackStruct& track, const Track& tmpTrack);
void RefitTracks(std::list<TrackStruct>& tracks);
void ImproveTracks(std::list<TrackStruct>& tracks);
void RemoveInvalidTracks(std::list<TrackStruct>& tracks);
void RemoveConnectedTracks(std::list<TrackStruct>& tracks, int stMin, int stMax, std::function<bool(const TrackStruct&, const TrackStruct&)> isBetter);
bool LoadCCDB();
bool LoadOCDB();
bool SetMagField();
void ExtrapToVertex(std::list<TrackStruct>& tracks, const std::vector<VertexStruct>& vertices, int event);
void ExtrapToVertex(TrackStruct& track, VertexStruct& vertex);
void ExtrapToMID(TrackParam& param);
void selectTracks(std::list<TrackStruct>& tracks);
bool IsSelected(TrackStruct& track);
const TrackMCHMID* FindMuon(uint32_t iMCHTrack, const std::vector<TrackMCHMID>& muonTracks);
int CompareEvents(std::list<TrackStruct>& tracks1, std::list<TrackStruct>& tracks2, double precision, bool printDiff, bool printAll, std::vector<TH1*>& histos);
bool AreTrackParamCompatible(const TrackParamStruct& param1, const TrackParamStruct& param2, double precision);
bool AreTrackParamCovCompatible(const TMatrixD& cov1, const TMatrixD& cov2, double precision);
int PrintEvent(const std::list<TrackStruct>& tracks);
void PrintTrack(const TrackStruct& track);
void PrintResiduals(const TrackParamStruct& param1, const TrackParamStruct& param2);
void PrintCovResiduals(const TMatrixD& cov1, const TMatrixD& cov2);
void FillResiduals(const TrackParamStruct& param1, const TrackParamStruct& param2, std::vector<TH1*>& histos);
void DrawResiduals(std::vector<TH1*>& histos);
void CreateHistosAtVertex(std::vector<TH1*>& histos, const char* extension);
void FillHistosAtVertex(const std::list<TrackStruct>& tracks, std::vector<TH1*>& histos);
void FillHistosMuAtVertex(const TrackStruct& track, std::vector<TH1*>& histos);
void FillHistosDimuAtVertex(const TrackStruct& track1, const TrackStruct& track2, std::vector<TH1*>& histos);
void DrawHistosAtVertex(std::vector<TH1*> histos[2]);
void FillComparisonsAtVertex(std::list<TrackStruct>& tracks1, std::list<TrackStruct>& tracks2, std::vector<TH1*> histos[5]);
void DrawComparisonsAtVertex(std::vector<TH1*> histos[5]);
void PrintStat(TString what, int n1, int n2);
void PrintStat(int nTracksAll[2], int nTracksMatch[2], std::vector<TH1*> histos[2]);

//_________________________________________________________________________________________________
void CompareTracks(int runNumber, string inFileName1, int versionFile1, string inFileName2, int versionFile2,
                   string vtxFileName = "", double precision = 1.e-4, bool applyTrackSelection = false,
                   bool printDiff = true, bool printAll = false)
{
  /// Compare the tracks stored in the 2 binary files
  /// file version 1: param at 1st cluster + clusters
  /// file version 2: param at 1st cluster + chi2 + clusters
  /// file version 3: param at vertex + dca + rAbs + chi2 + param at 1st cluster + clusters
  /// file version 4: param at vertex + dca + rAbs + chi2 + param at 1st cluster + clusters (v2)
  /// file version 5: nTracksAtVtx + nMCHTracks + nClusters + list of Tracks at vertex (param at vertex + dca + rAbs + mchTrackIdx) + list of MCH tracks (v1) + list of clusters (v2)
  /// file version 6: nTracksAtVtx + nMCHTracks + nClusters + list of Tracks at vertex (param at vertex + dca + rAbs + mchTrackIdx) + list of MCH tracks (v2) + list of clusters (v2)
  /// file version 7: nTracksAtVtx + nMCHTracks + nClusters + list of Tracks at vertex (param at vertex + dca + rAbs + mchTrackIdx) + list of MCH tracks (v3) + list of clusters (v2)

  run = runNumber;

  // get vertices and prepare track extrapolation
  std::vector<VertexStruct> vertices{};
  if (!vtxFileName.empty()) {
    ifstream inFile(vtxFileName,ios::binary);
    if (!inFile.is_open()) {
      cout << "fail opening vertex file" << endl;
      return;
    }
    int event(-1), iEvent(-1);
    while((event = ReadNextEvent(inFile, vertices)) >= 0) {
      if (event != ++iEvent) {
        cout << "event " << iEvent << " missing" << endl;
        return;
      }
    }
  }

  // open files
  ifstream inFile1(inFileName1,ios::binary);
  ifstream inFile2(inFileName2,ios::binary);
  if (!inFile1.is_open() || !inFile2.is_open()) {
    cout << "fail opening files" << endl;
    return;
  }

  int event1(-1), event2(-1);
  std::list<TrackStruct> tracks1{};
  std::list<TrackStruct> tracks2{};
  bool readNextEvent1(true), readNextEvent2(true);
  int nDifferences(0);
  std::vector<TH1*> residualsAtFirstCluster{};
  std::vector<TH1*> comparisonsAtVertex[5] = {{}, {}, {}, {}, {}};
  CreateHistosAtVertex(comparisonsAtVertex[0], "identical");
  CreateHistosAtVertex(comparisonsAtVertex[1], "similar1");
  CreateHistosAtVertex(comparisonsAtVertex[2], "similar2");
  CreateHistosAtVertex(comparisonsAtVertex[3], "additional");
  CreateHistosAtVertex(comparisonsAtVertex[4], "missing");
  std::vector<TH1*> histosAtVertex[2] = {{}, {}};
  CreateHistosAtVertex(histosAtVertex[0], "1");
  CreateHistosAtVertex(histosAtVertex[1], "2");
 
  while (true) {

    if (readNextEvent1) {
      if (versionFile1 < 5) {
        event1 = ReadNextEvent(inFile1, versionFile1, tracks1);
      } else if (versionFile1 == 5) {
        ReadNextEventV5<TrackMCHv1>(inFile1, event1, tracks1);
      } else if (versionFile1 == 6) {
        ReadNextEventV5<TrackMCHv2>(inFile1, event1, tracks1);
      } else {
        ReadNextEventV5<TrackMCH>(inFile1, event1, tracks1);
      }
      //trackFitter.useChamberResolution();
      //ImproveTracks(tracks1);
      //RemoveInvalidTracks(tracks1);
      //RemoveConnectedTracks(tracks1, 2, 4, &TrackStruct::isBetter1);
      //trackFitter.useClusterResolution();
      //RefitTracks(tracks1);
      ExtrapToVertex(tracks1, vertices, event1);
      if (applyTrackSelection) {
        selectTracks(tracks1);
      }
      FillHistosAtVertex(tracks1, histosAtVertex[0]);
    }

    if (readNextEvent2) {
      if (versionFile2 < 5) {
        event2 = ReadNextEvent(inFile2, versionFile2, tracks2);
      } else if (versionFile2 == 5) {
        ReadNextEventV5<TrackMCHv1>(inFile2, event2, tracks2);
      } else if (versionFile2 == 6) {
        ReadNextEventV5<TrackMCHv2>(inFile2, event2, tracks2);
      } else {
        ReadNextEventV5<TrackMCH>(inFile2, event2, tracks2);
      }
      //RemoveInvalidTracks(tracks2);
      //RemoveConnectedTracks(tracks2, 2, 4, &TrackStruct::isBetter42);
      //trackFitter.useClusterResolution();
      //RefitTracks(tracks2);
      ExtrapToVertex(tracks2, vertices, event2);
      if (applyTrackSelection) {
        selectTracks(tracks2);
      }
      FillHistosAtVertex(tracks2, histosAtVertex[1]);
    }

    // reaching end of both files
    if (event1 < 0 && event2 < 0) {
      break;
    }

    if (event1 == event2) {
      // reading the same event --> we can compare tracks
      int nDiff = CompareEvents(tracks1, tracks2, precision, printDiff, printAll, residualsAtFirstCluster);
      if (nDiff > 0) {
        if (printDiff) {
          cout << "--> " << nDiff << " differences found in event " << event1 << endl;
        }
        nDifferences += nDiff;
      }
      FillComparisonsAtVertex(tracks1, tracks2, comparisonsAtVertex);
      readNextEvent1 = true;
      readNextEvent2 = true;
    } else if (event2 < 0 || (event1 >= 0 && event1 < event2)) {
      // event 2 > event 1 or reaching end of file 2
      if (tracks1.size() > 0) {
        if (printDiff) {
          cout << "tracks in event " << event1 << " are missing in file 2" << endl;
          nDifferences += PrintEvent(tracks1);
        } else {
          nDifferences += tracks1.size();
        }
        FillHistosAtVertex(tracks1, comparisonsAtVertex[4]);
      }
      readNextEvent1 = true;
      readNextEvent2 = false;
    } else {
      // event 1 > event 2 or reaching end of file 1
      if (tracks2.size() > 0) {
        if (printDiff) {
          cout << "tracks in event " << event2 << " are missing in file 1" << endl;
          nDifferences += PrintEvent(tracks2);
        } else {
          nDifferences += tracks2.size();
        }
        FillHistosAtVertex(tracks2, comparisonsAtVertex[3]);
      }
      readNextEvent1 = false;
      readNextEvent2 = true;
    }
  }
  
  cout << "number of different tracks = " << nDifferences << endl;
  gStyle->SetOptStat(111111);
  DrawResiduals(residualsAtFirstCluster);
  DrawHistosAtVertex(histosAtVertex);
  DrawComparisonsAtVertex(comparisonsAtVertex);

  inFile1.close();
  inFile2.close();

}

//_________________________________________________________________________________________________
void CompareTracks(int runNumber, string mchFileName1, string muonFileName1, string mchFileName2, string muonFileName2,
                   string vtxFileName = "", double precision = 1.e-4,
                   bool applyTrackSelection = false, bool selectMatched = false,
                   bool printStat = true, bool printDiff = true, bool printAll = false)
{
  /// Compare the tracks stored in the 2 sets of root files

  run = runNumber;

  // get vertices and prepare track extrapolation
  std::vector<VertexStruct> vertices{};
  if (!vtxFileName.empty()) {
    ifstream inFile(vtxFileName, ios::binary);
    if (!inFile.is_open()) {
      cout << "fail opening vertex file" << endl;
      return;
    }
    int event(-1), iEvent(-1);
    while ((event = ReadNextEvent(inFile, vertices)) >= 0) {
      if (event != ++iEvent) {
        cout << "event " << iEvent << " missing" << endl;
        return;
      }
    }
  }

  // load tracks
  auto [fMCH1, mchReader1] = LoadData(mchFileName1.c_str(), "o2sim");
  TTreeReaderValue<std::vector<ROFRecord>> mchROFs1 = {*mchReader1, "trackrofs"};
  TTreeReaderValue<std::vector<TrackMCH>> mchTracks1 = {*mchReader1, "tracks"};
  TTreeReaderValue<std::vector<Cluster>> mchClusters1 = {*mchReader1, "trackclusters"};
  auto [fMUON1, muonReader1] = LoadData(muonFileName1.c_str(), "o2sim");
  TTreeReaderValue<std::vector<TrackMCHMID>> muonTracks1 = {*muonReader1, "tracks"};
  auto [fMCH2, mchReader2] = LoadData(mchFileName2.c_str(), "o2sim");
  TTreeReaderValue<std::vector<ROFRecord>> mchROFs2 = {*mchReader2, "trackrofs"};
  TTreeReaderValue<std::vector<TrackMCH>> mchTracks2 = {*mchReader2, "tracks"};
  TTreeReaderValue<std::vector<Cluster>> mchClusters2 = {*mchReader2, "trackclusters"};
  auto [fMUON2, muonReader2] = LoadData(muonFileName2.c_str(), "o2sim");
  TTreeReaderValue<std::vector<TrackMCHMID>> muonTracks2 = {*muonReader2, "tracks"};
  int nTF = mchReader1->GetEntries(false);
  if (muonReader1->GetEntries(false) != nTF ||
      mchReader2->GetEntries(false) != nTF ||
      muonReader2->GetEntries(false) != nTF) {
    LOG(error) << " not all files contain the same number of TF";
    exit(-1);
  }

  int event1(-1), event2(-1), iTF(-1);
  std::list<TrackStruct> tracks1{};
  std::list<TrackStruct> tracks2{};
  int nDifferences(0);
  int nTracksAll[2] = {0, 0};
  int nTracksMatch[2] = {0, 0};
  std::vector<TH1*> residualsAtFirstCluster{};
  std::vector<TH1*> comparisonsAtVertex[5] = {{}, {}, {}, {}, {}};
  CreateHistosAtVertex(comparisonsAtVertex[0], "identical");
  CreateHistosAtVertex(comparisonsAtVertex[1], "similar1");
  CreateHistosAtVertex(comparisonsAtVertex[2], "similar2");
  CreateHistosAtVertex(comparisonsAtVertex[3], "additional");
  CreateHistosAtVertex(comparisonsAtVertex[4], "missing");
  std::vector<TH1*> histosAtVertex[2] = {{}, {}};
  CreateHistosAtVertex(histosAtVertex[0], "1");
  CreateHistosAtVertex(histosAtVertex[1], "2");

  while (mchReader1->Next() && muonReader1->Next() && mchReader2->Next() && muonReader2->Next()) {

    ++iTF;
    auto nROFs = TMath::Max(mchROFs1->size(), mchROFs2->size());
    int iROF1(-1), iROF2(-1);
    for (size_t iROF = 0; iROF < nROFs; ++iROF) {

      if (iROF < mchROFs1->size()) {
        ++iROF1;
        const auto& mchROF1 = (*mchROFs1)[iROF];
        tracks1.clear();
        for (int iMCHTrack = mchROF1.getFirstIdx(); iMCHTrack <= mchROF1.getLastIdx(); ++iMCHTrack) {
          ++nTracksAll[0];
          auto muon = FindMuon(iMCHTrack, *muonTracks1);
          if (muon) {
            ++nTracksMatch[0];
          } else if (selectMatched) {
            continue;
          }
          tracks1.emplace_back();
          FillTrack<TrackMCH>(tracks1.back(), nullptr, &(*mchTracks1)[iMCHTrack], *mchClusters1);
          tracks1.back().matchedChi2 = muon ? muon->getMatchChi2OverNDF() : -1.;
          tracks1.back().needExtrapToVtx = true;
        }
        // RefitTracks(tracks1);
        ExtrapToVertex(tracks1, vertices, ++event1);
        if (applyTrackSelection) {
          selectTracks(tracks1);
        }
        FillHistosAtVertex(tracks1, histosAtVertex[0]);
      }

      if (iROF < mchROFs2->size()) {
        ++iROF2;
        const auto& mchROF2 = (*mchROFs2)[iROF];
        tracks2.clear();
        for (int iMCHTrack = mchROF2.getFirstIdx(); iMCHTrack <= mchROF2.getLastIdx(); ++iMCHTrack) {
          ++nTracksAll[1];
          auto muon = FindMuon(iMCHTrack, *muonTracks2);
          if (muon) {
            ++nTracksMatch[1];
          } else if (selectMatched) {
            continue;
          }
          tracks2.emplace_back();
          FillTrack<TrackMCH>(tracks2.back(), nullptr, &(*mchTracks2)[iMCHTrack], *mchClusters2);
          tracks2.back().matchedChi2 = muon ? muon->getMatchChi2OverNDF() : -1.;
          tracks2.back().needExtrapToVtx = true;
        }
        // RefitTracks(tracks2);
        ExtrapToVertex(tracks2, vertices, ++event2);
        if (applyTrackSelection) {
          selectTracks(tracks2);
        }
        FillHistosAtVertex(tracks2, histosAtVertex[1]);
      }

      if (iROF1 == iROF2) {
        // reading the same ROF --> we can compare tracks
        int nDiff = CompareEvents(tracks1, tracks2, precision, printDiff, printAll, residualsAtFirstCluster);
        if (nDiff > 0) {
          if (printDiff) {
            cout << "--> " << nDiff << " differences found in TF/ROF " << iTF << "/" << iROF1 << endl;
          }
          nDifferences += nDiff;
        }
        FillComparisonsAtVertex(tracks1, tracks2, comparisonsAtVertex);
      } else if (iROF1 > iROF2) {
        // ROF 1 missing in file 2
        if (printDiff) {
          cout << "tracks in TF/ROF " << iTF << "/" << iROF1 << " are missing in file 2" << endl;
          nDifferences += PrintEvent(tracks1);
        } else {
          nDifferences += tracks1.size();
        }
        FillHistosAtVertex(tracks1, comparisonsAtVertex[4]);
      } else {
        // ROF 2 missing in file 1
        if (printDiff) {
          cout << "tracks in TF/ROF " << iTF << "/" << iROF2 << " are missing in file 1" << endl;
          nDifferences += PrintEvent(tracks2);
        } else {
          nDifferences += tracks2.size();
        }
        FillHistosAtVertex(tracks2, comparisonsAtVertex[3]);
      }
    }
  }

  cout << "number of different tracks = " << nDifferences << endl;
  gStyle->SetOptStat(111111);
  DrawResiduals(residualsAtFirstCluster);
  DrawHistosAtVertex(histosAtVertex);
  DrawComparisonsAtVertex(comparisonsAtVertex);

  if (printStat) {
    PrintStat(nTracksAll, nTracksMatch, histosAtVertex);
  }

  fMCH1->Close();
  fMUON1->Close();
  fMCH2->Close();
  fMUON2->Close();
}

//_________________________________________________________________________________________________
int ReadNextEvent(ifstream& inFile, std::vector<VertexStruct>& vertices)
{
  /// read the next vertex in the input file
  /// return the event number or -1 when reaching end of file

  int event(-1);
  if (!inFile.read(reinterpret_cast<char*>(&event), sizeof(int))) {
    // reaching end of file
    return -1;
  }

  // read the vertex
  vertices.emplace_back();
  inFile.read(reinterpret_cast<char*>(&vertices.back()), sizeof(VertexStruct));

  return event;
}

//_________________________________________________________________________________________________
int ReadNextEvent(ifstream& inFile, int version, std::list<TrackStruct>& tracks)
{
  /// read the next event in the input file

  tracks.clear();

  int event(-1);
  if (!inFile.read(reinterpret_cast<char*>(&event), sizeof(int))) {
    // reaching end of file
    return -1;
  }
 
  int size(0);
  inFile.read(reinterpret_cast<char*>(&size), sizeof(int));
  if (size == 0) {
    // empty event
    return event;
  }

  // read the tracks
  int nTracks(0);
  inFile.read(reinterpret_cast<char*>(&nTracks), sizeof(int));
  for (Int_t iTrack = 0; iTrack < nTracks; ++iTrack) {
    tracks.emplace_back();
    ReadTrack(inFile, tracks.back(), version);
    if (version < 3) {
      tracks.back().needExtrapToVtx = true;
    }
  }

  return event;
}

//_________________________________________________________________________________________________
void ReadTrack(ifstream& inFile, TrackStruct& track, int version)
{
  /// read one track from the input file

  if (version > 2) {
    TrackParamStruct paramAtVtx{};
    inFile.read(reinterpret_cast<char*>(&paramAtVtx), sizeof(TrackParamStruct));
    track.pxpypzm.SetPx(paramAtVtx.px);
    track.pxpypzm.SetPy(paramAtVtx.py);
    track.pxpypzm.SetPz(paramAtVtx.pz);
    track.pxpypzm.SetM(muMass);
    track.sign = paramAtVtx.sign;
    inFile.read(reinterpret_cast<char*>(&(track.dca)), sizeof(double));
    inFile.read(reinterpret_cast<char*>(&(track.rAbs)), sizeof(double));
  }

  inFile.read(reinterpret_cast<char*>(&(track.param)), sizeof(TrackParamStruct));

  if (version > 1) {
    inFile.read(reinterpret_cast<char*>(&(track.chi2)), sizeof(double));
  }

  int nClusters(0);
  inFile.read(reinterpret_cast<char*>(&nClusters), sizeof(int));
  track.clusters.resize(nClusters);
  for (Int_t iCl = 0; iCl < nClusters; ++iCl) {
    if (version < 4) {
      inFile.read(reinterpret_cast<char*>(&track.clusters[iCl]), sizeof(ClusterStructV1));
    } else {
      inFile.read(reinterpret_cast<char*>(&track.clusters[iCl]), sizeof(Cluster));
    }
  }
}

//_________________________________________________________________________________________________
template <class T>
void ReadNextEventV5(ifstream& inFile, int& event, std::list<TrackStruct>& tracks)
{
  /// read the next event in the input file

  tracks.clear();

  int nTracksAtVtx(-1);
  if (!inFile.read(reinterpret_cast<char*>(&nTracksAtVtx), sizeof(int))) {
    event = -1;
    return; // reaching end of file ...
  } else {
    ++event; // ... or reading the next event
  }
  int nMCHTracks(-1);
  inFile.read(reinterpret_cast<char*>(&nMCHTracks), sizeof(int));
  int nClusters(-1);
  inFile.read(reinterpret_cast<char*>(&nClusters), sizeof(int));
  if (nTracksAtVtx == 0 && nMCHTracks == 0) {
    return; // event is empty
  }

  // read the tracks at vertex, the MCH tracks and the attached clusters
  std::vector<TrackAtVtxStruct> tracksAtVtx(nTracksAtVtx);
  inFile.read(reinterpret_cast<char*>(tracksAtVtx.data()), nTracksAtVtx * sizeof(TrackAtVtxStruct));
  std::vector<T> mchTracks(nMCHTracks);
  inFile.read(reinterpret_cast<char*>(mchTracks.data()), nMCHTracks * sizeof(T));
  std::vector<Cluster> clusters(nClusters);
  inFile.read(reinterpret_cast<char*>(clusters.data()), nClusters * sizeof(Cluster));

  if (nTracksAtVtx > 0) {

    // fill the internal track structure based on the provided tracks at vertex
    for (const auto& trackAtVtx : tracksAtVtx) {
      tracks.emplace_back();
      FillTrack(tracks.back(), &trackAtVtx, (nMCHTracks > 0) ? &(mchTracks[trackAtVtx.mchTrackIdx]) : nullptr, clusters);
    }

  } else {

    // fill the internal track structure based on the MCH tracks
    for (const auto& mchTrack : mchTracks) {
      tracks.emplace_back();
      FillTrack(tracks.back(), nullptr, &mchTrack, clusters);
      tracks.back().needExtrapToVtx = true;
    }
  }
}

//_________________________________________________________________________________________________
std::tuple<TFile*, TTreeReader*> LoadData(const char* fileName, const char* treeName)
{
  /// open the input file and get the intput tree

  TFile* f = TFile::Open(fileName, "READ");
  if (!f || f->IsZombie()) {
    LOG(error) << "opening file " << fileName << " failed";
    exit(-1);
  }

  TTreeReader* r = new TTreeReader(treeName, f);
  if (r->IsZombie()) {
    LOG(error) << "tree " << treeName << " not found";
    exit(-1);
  }

  return std::make_tuple(f, r);
}

//_________________________________________________________________________________________________
template <typename T>
void FillTrack(TrackStruct& track, const TrackAtVtxStruct* trackAtVtx, const T* mchTrack, std::vector<Cluster>& clusters)
{
  /// fill the internal track structure from the provided informations

  if (trackAtVtx) {
    track.pxpypzm.SetPx(trackAtVtx->paramAtVertex.px);
    track.pxpypzm.SetPy(trackAtVtx->paramAtVertex.py);
    track.pxpypzm.SetPz(trackAtVtx->paramAtVertex.pz);
    track.pxpypzm.SetM(muMass);
    track.sign = trackAtVtx->paramAtVertex.sign;
    track.dca = trackAtVtx->dca;
    track.rAbs = trackAtVtx->rAbs;
  }

  if (mchTrack) {
    track.chi2 = mchTrack->getChi2();
    track.param.x = mchTrack->getX();
    track.param.y = mchTrack->getY();
    track.param.z = mchTrack->getZ();
    track.param.px = mchTrack->getPx();
    track.param.py = mchTrack->getPy();
    track.param.pz = mchTrack->getPz();
    track.param.sign = mchTrack->getSign();
    for (int i = 0; i < 5; i++) {
      for (int j = 0; j < 5; j++) {
        track.cov(i, j) = mchTrack->getCovariance(i, j);
      }
    }

    if (std::is_same_v<T, TrackMCH>) {
      const auto* t = reinterpret_cast<const TrackMCH*>(mchTrack);
      TrackParam paramAtMID(t->getZAtMID(), t->getParametersAtMID(), t->getCovariancesAtMID());
      track.paramAtMID.x = paramAtMID.getNonBendingCoor();
      track.paramAtMID.y = paramAtMID.getBendingCoor();
      track.paramAtMID.z = paramAtMID.getZ();
      track.paramAtMID.px = paramAtMID.px();
      track.paramAtMID.py = paramAtMID.py();
      track.paramAtMID.pz = paramAtMID.pz();
      track.paramAtMID.sign = paramAtMID.getCharge();
      track.covAtMID = paramAtMID.getCovariances();
    }

    track.clusters.insert(track.clusters.begin(), clusters.begin() + mchTrack->getFirstClusterIdx(),
                          clusters.begin() + mchTrack->getLastClusterIdx() + 1);
  }
}

//_________________________________________________________________________________________________
void UpdateTrack(TrackStruct& track, const Track& tmpTrack)
{
  /// update the track parameters from the temporary track

  const auto& param = tmpTrack.first();
  track.chi2 = param.getTrackChi2();
  track.param.x = param.getNonBendingCoor();
  track.param.y = param.getBendingCoor();
  track.param.z = param.getZ();
  track.param.px = param.px();
  track.param.py = param.py();
  track.param.pz = param.pz();
  track.param.sign = param.getCharge();

  const TMatrixD& cov = param.getCovariances();
  for (int i = 0; i < 5; i++) {
    for (int j = 0; j <= i; j++) {
      track.cov(i, j) = track.cov(j, i) = cov(i, j);
    }
  }

  TrackParam paramAtMID(tmpTrack.last());
  ExtrapToMID(paramAtMID);
  track.paramAtMID.x = paramAtMID.getNonBendingCoor();
  track.paramAtMID.y = paramAtMID.getBendingCoor();
  track.paramAtMID.z = paramAtMID.getZ();
  track.paramAtMID.px = paramAtMID.px();
  track.paramAtMID.py = paramAtMID.py();
  track.paramAtMID.pz = paramAtMID.pz();
  track.paramAtMID.sign = paramAtMID.getCharge();

  const TMatrixD& covAtMID = paramAtMID.getCovariances();
  for (int i = 0; i < 5; i++) {
    for (int j = 0; j <= i; j++) {
      track.covAtMID(i, j) = track.covAtMID(j, i) = covAtMID(i, j);
    }
  }

  // update the list of clusters only if some have been removed
  if (tmpTrack.getNClusters() != (int)track.clusters.size()) {
    std::vector<Cluster> clusters{};
    clusters.reserve(tmpTrack.getNClusters());
    for (auto& param : tmpTrack) {
      clusters.push_back(std::move(*param.getClusterPtr()));
    }
    track.clusters.swap(clusters);
  }

  track.needExtrapToVtx = true;
}

//_________________________________________________________________________________________________
void RefitTracks(std::list<TrackStruct>& tracks)
{
  /// refit the tracks to the attached clusters

  // load OCDB (only once)
  if (!ocdbLoaded) {
    if (!LoadOCDB()) {
      cout << "fail loading OCDB objects for track extrapolation" << endl;
      exit(1);
    }
    ocdbLoaded = true;
  }

  // same extrapolation method as in TrackFinder
  TrackExtrap::useExtrapV2();

  for (auto itTrack = tracks.begin(); itTrack != tracks.end();) {

    if (itTrack->clusters.size() == 0) {
      cout << "track does not contain clusters --> unable to refit" << endl;
      exit(1);
    }

    Track track{};
    for (auto& cluster : itTrack->clusters) {
      // if (cluster.ey < 1.) {
      //   cluster.ey = 0.05;
      // }
      // if (cluster.getChamberId() == 2 || cluster.getChamberId() == 3) {
      //   cluster.y += gRandom->Gaus(0., 0.05);
      //   cluster.ey *= 2.;
      // }
      track.createParamAtCluster(cluster);
    }

    try {
      trackFitter.fit(track);
      UpdateTrack(*itTrack, track);
      ++itTrack;
    } catch (exception const& e) {
      itTrack = tracks.erase(itTrack);
    }
  }
}

//_________________________________________________________________________________________________
void ImproveTracks(std::list<TrackStruct>& tracks)
{
  /// Improve tracks by removing removable clusters with local chi2 higher than the defined cut
  /// Removable clusters are identified by the method Track::tagRemovableClusters()
  /// Recompute track parameters and covariances at the remaining clusters
  /// Remove the track if it cannot be improved or in case of failure

  // load OCDB (only once)
  if (!ocdbLoaded) {
    if (!LoadOCDB()) {
      cout << "fail loading OCDB objects for track extrapolation" << endl;
      exit(1);
    }
    ocdbLoaded = true;
  }

  // same extrapolation method as in TrackFinder
  TrackExtrap::useExtrapV2();

  // Maximum chi2 to keep a cluster (the factor 2 is for the 2 degrees of freedom: x and y)
  static const double maxChi2OfCluster = 2. * 4. * 4.;

  for (auto itTrack = tracks.begin(); itTrack != tracks.end();) {

    if (itTrack->clusters.size() == 0) {
      cout << "track does not contain clusters --> unable to improve" << endl;
      exit(1);
    }

    // create a temporary tracking track and fit it without running the smoother
    Track track{};
    for (const auto& cluster : itTrack->clusters) {
      track.createParamAtCluster(cluster);
    }
    try {
      trackFitter.fit(track, false);
    } catch (exception const& e) {
      itTrack = tracks.erase(itTrack);
      continue;
    }

    bool removeTrack(false);

    // At the first step, only run the smoother
    auto itStartingParam = std::prev(track.rend());

    while (true) {

      // Refit the part of the track affected by the cluster removal, run the smoother, but do not finalize
      try {
        trackFitter.fit(track, true, false, (itStartingParam == track.rbegin()) ? nullptr : &itStartingParam);
      } catch (exception const&) {
        removeTrack = true;
        break;
      }

      // Identify removable clusters
      track.tagRemovableClusters(0x1F, false);

      // Look for the cluster with the worst local chi2
      double worstLocalChi2(-1.);
      auto itWorstParam(track.end());
      for (auto itParam = track.begin(); itParam != track.end(); ++itParam) {
        if (itParam->getLocalChi2() > worstLocalChi2) {
          worstLocalChi2 = itParam->getLocalChi2();
          itWorstParam = itParam;
        }
      }

      // If the worst chi2 is under requirement then the track is improved
      if (worstLocalChi2 < maxChi2OfCluster) {
        break;
      }

      // If the worst cluster is not removable then the track cannot be improved
      if (!itWorstParam->isRemovable()) {
        removeTrack = true;
        break;
      }

      // Remove the worst cluster
      auto itNextParam = track.removeParamAtCluster(itWorstParam);

      // Decide from where to refit the track: from the cluster next the one suppressed or
      // from scratch if the removed cluster was used to compute the tracking seed
      itStartingParam = track.rbegin();
      auto itNextToNextParam = (itNextParam == track.end()) ? itNextParam : std::next(itNextParam);
      while (itNextToNextParam != track.end()) {
        if (itNextToNextParam->getClusterPtr()->getChamberId() != itNextParam->getClusterPtr()->getChamberId()) {
          itStartingParam = std::make_reverse_iterator(++itNextParam);
          break;
        }
        ++itNextToNextParam;
      }
    }

    // Remove the track if it couldn't be improved or update the parameters
    if (removeTrack) {
      itTrack = tracks.erase(itTrack);
    } else {
      for (auto& param : track) {
        param.setParameters(param.getSmoothParameters());
        param.setCovariances(param.getSmoothCovariances());
      }
      UpdateTrack(*itTrack, track);
      ++itTrack;
    }
  }
}

//_________________________________________________________________________________________________
void RemoveInvalidTracks(std::list<TrackStruct>& tracks)
{
  /// Find and remove tracks that do not pass all the tracking criteria,
  /// including 3 chambers fired out of 4 in stations 4 & 5
  for (auto itTrack = tracks.begin(); itTrack != tracks.end();) {
    if (!itTrack->isValid()) {
      itTrack = tracks.erase(itTrack);
    } else {
      ++itTrack;
    }
  }
}

//_________________________________________________________________________________________________
void RemoveConnectedTracks(std::list<TrackStruct>& tracks, int stMin, int stMax,
                           std::function<bool(const TrackStruct&, const TrackStruct&)> isBetter)
{
  /// Find and remove tracks sharing 1 cluster or more in station(s) [stMin, stMax]
  /// For each couple of connected tracks, one removes the worst one as determined by the isBetter function

  if (tracks.size() < 2) {
    return;
  }

  int chMin = 2 * stMin;
  int chMax = 2 * stMax + 1;
  int nPlane = 2 * (chMax - chMin + 1);

  // first loop to fill the array of cluster Ids
  std::vector<uint32_t> ClIds{};
  ClIds.resize(nPlane * tracks.size());
  int iTrack(0);
  for (auto itTrack = tracks.begin(); itTrack != tracks.end(); ++itTrack, ++iTrack) {
    for (auto itCl = itTrack->clusters.rbegin(); itCl != itTrack->clusters.rend(); ++itCl) {
      int ch = itCl->getChamberId();
      if (ch > chMax) {
        continue;
      } else if (ch < chMin) {
        break;
      }
      ClIds[nPlane * iTrack + 2 * (ch - chMin) + itCl->getDEId() % 2] = itCl->uid;
    }
  }

  // second loop to tag the tracks to remove
  int iindex = ClIds.size() - 1;
  for (auto itTrack1 = tracks.rbegin(); itTrack1 != tracks.rend(); ++itTrack1, iindex -= nPlane) {
    int jindex = iindex - nPlane;
    for (auto itTrack2 = std::next(itTrack1); itTrack2 != tracks.rend(); ++itTrack2) {
      for (int iPlane = nPlane; iPlane > 0; --iPlane) {
        if (ClIds[iindex] > 0 && ClIds[iindex] == ClIds[jindex]) {
          if (isBetter(*itTrack2, *itTrack1)) {
            itTrack1->connected = true;
          } else {
            itTrack2->connected = true;
          }
          iindex -= iPlane;
          jindex -= iPlane;
          break;
        }
        --iindex;
        --jindex;
      }
      iindex += nPlane;
    }
  }

  // third loop to remove them. That way all combinations are tested.
  for (auto itTrack = tracks.begin(); itTrack != tracks.end();) {
    if (itTrack->connected) {
      itTrack = tracks.erase(itTrack);
    } else {
      ++itTrack;
    }
  }
}

//_________________________________________________________________________________________________
bool LoadCCDB()
{
  /// access CCDB and prepare track extrapolation to vertex and track fitting

  // load magnetic field and geometry from CCDB
  auto& ccdb = o2::ccdb::BasicCCDBManager::instance();
  auto [tStart, tEnd] = ccdb.getRunDuration(run);
  ccdb.setTimestamp(tEnd);
  auto grp = ccdb.get<o2::parameters::GRPMagField>("GLO/Config/GRPMagField");
  // ccdb.setURL("http://localhost:6060");
  auto geom = ccdb.get<TGeoManager>("GLO/Config/GeometryAligned");

  // and prepare track extrapolation to vertex and track fitting
  o2::base::Propagator::initFieldFromGRP(grp);
  TrackExtrap::setField();
  trackFitter.smoothTracks(true);
  trackFitter.setChamberResolution(0.2, 0.2);

  return true;
}

//_________________________________________________________________________________________________
bool LoadOCDB()
{
  /// access OCDB and prepare track extrapolation to vertex and track fitting

  // return LoadCCDB();

  AliCDBManager* man = AliCDBManager::Instance();
  if (gSystem->AccessPathName("OCDB.root", kFileExists)==0) {
    man->SetDefaultStorage("local:///dev/null");
    man->SetSnapshotMode("OCDB.root");
  } else {
    man->SetDefaultStorage("local://./OCDB");
  }
  man->SetRun(run);

  if (!SetMagField()) {
    //  if (!AliMUONCDB::LoadField()) {
    return false;
  }
  TrackExtrap::setField();

  // prepare the track fitter
  trackFitter.smoothTracks(true);
  trackFitter.setChamberResolution(0.2, 0.2);

  o2::base::GeometryManager::loadGeometry("O2geometry.root");
  if (!gGeoManager) {
    return false;
  }
  //  AliGeomManager::LoadGeometry();
  //  if (!AliGeomManager::GetGeometry() || !AliGeomManager::ApplyAlignObjsFromCDB("MUON")) {
  //    return false;
  //  }

  return true;
}

//_________________________________________________________________________________________________
bool SetMagField()
{
  /// Set the magnetic field using O2 maps and GRP info

  AliGRPManager grpMan;
  if (!grpMan.ReadGRPEntry()) {
    Error("SetMagField", "failed to load GRP Data from OCDB");
    return false;
  }

  const AliGRPObject* grpData = grpMan.GetGRPData();
  if (!grpData) {
    Error("SetMagField", "GRP Data is not loaded");
    return false;
  }

  float l3Current = grpData->GetL3Current((AliGRPObject::Stats)0);
  if (l3Current == AliGRPObject::GetInvalidFloat()) {
    Error("SetMagField", "GRP/GRP/Data entry:  missing value for the L3 current !");
    return false;
  }

  char l3Polarity = grpData->GetL3Polarity();
  if (l3Polarity == AliGRPObject::GetInvalidChar()) {
    Error("SetMagField", "GRP/GRP/Data entry:  missing value for the L3 polarity !");
    return false;
  }

  float diCurrent = grpData->GetDipoleCurrent((AliGRPObject::Stats)0);
  if (diCurrent == AliGRPObject::GetInvalidFloat()) {
    Error("SetMagField", "GRP/GRP/Data entry:  missing value for the dipole current !");
    return false;
  }

  char diPolarity = grpData->GetDipolePolarity();
  if (diPolarity == AliGRPObject::GetInvalidChar()) {
    Error("SetMagField", "GRP/GRP/Data entry:  missing value for the dipole polarity !");
    return false;
  }

  float beamEnergy = grpData->GetBeamEnergy();
  if (beamEnergy == AliGRPObject::GetInvalidFloat()) {
    Error("SetMagField", "GRP/GRP/Data entry:  missing value for the beam energy !");
    return false;
  }

  TString beamType = grpData->GetBeamType();
  if (beamType == AliGRPObject::GetInvalidString()) {
    Error("SetMagField", "GRP/GRP/Data entry:  missing value for the beam type !");
    return false;
  }

  Info("SetMagField", "l3Current = %f, diCurrent = %f", TMath::Abs(l3Current) * (l3Polarity ? -1 : 1),
       TMath::Abs(diCurrent) * (diPolarity ? -1 : 1));

  auto field = o2::field::MagneticField::createFieldMap(TMath::Abs(l3Current) * (l3Polarity ? -1 : 1),
                                                        TMath::Abs(diCurrent) * (diPolarity ? -1 : 1),
                                                        o2::field::MagneticField::kConvLHC, false, beamEnergy, beamType.Data(),
                                                        "$(O2_ROOT)/share/Common/maps/mfchebKGI_sym.root");
  TGeoGlobalMagField::Instance()->SetField(field);
  TGeoGlobalMagField::Instance()->Lock();

  return true;
}

//_________________________________________________________________________________________________
void ExtrapToVertex(std::list<TrackStruct>& tracks, const std::vector<VertexStruct>& vertices, int event)
{
  /// extrapolate to the vertex all tracks that need to be

  // same extrapolation method as in TrackAtVertexSpec
  TrackExtrap::useExtrapV2(false);

  // get the vertex corresponding to this event or use (0.,0.,0.)
  VertexStruct vertex = {0., 0., 0.};
  if (!vertices.empty()) {
    if (event < (int)vertices.size()) {
      vertex = vertices[event];
    } else {
      cout << "missing vertex for event " << event << endl;
    }
  }

  // loop over tracks and extrapolate them to the vertex if needed
  for (auto& track : tracks) {
    if (track.needExtrapToVtx) {
      ExtrapToVertex(track, vertex);
    }
  }
}

//_________________________________________________________________________________________________
void ExtrapToVertex(TrackStruct& track, VertexStruct& vertex)
{
  /// compute the track parameters at vertex, at DCA and at the end of the absorber

  // load OCDB (only once)
  if (!ocdbLoaded) {
    if (!LoadOCDB()) {
      cout << "fail loading OCDB objects for track extrapolation" << endl;
      exit(1);
    }
    ocdbLoaded = true;
  }

  // convert parameters at first cluster in TrackParam format
  TrackParam trackParam;
  trackParam.setNonBendingCoor(track.param.x);
  trackParam.setBendingCoor(track.param.y);
  trackParam.setZ(track.param.z);
  trackParam.setNonBendingSlope(track.param.px / track.param.pz);
  trackParam.setBendingSlope(track.param.py / track.param.pz);
  trackParam.setInverseBendingMomentum(track.param.sign/TMath::Sqrt(track.param.py*track.param.py + track.param.pz*track.param.pz));

  // extrapolate to vertex
  TrackParam trackParamAtVertex(trackParam);
  TrackExtrap::extrapToVertex(trackParamAtVertex, vertex.x, vertex.y, vertex.z, 0., 0.);
  track.pxpypzm.SetPx(trackParamAtVertex.px());
  track.pxpypzm.SetPy(trackParamAtVertex.py());
  track.pxpypzm.SetPz(trackParamAtVertex.pz());
  track.pxpypzm.SetM(muMass);
  track.sign = trackParamAtVertex.getCharge();

  // extrapolate to DCA
  TrackParam trackParamAtDCA(trackParam);
  TrackExtrap::extrapToVertexWithoutBranson(trackParamAtDCA, vertex.z);
  double dcaX = trackParamAtDCA.getNonBendingCoor() - vertex.x;
  double dcaY = trackParamAtDCA.getBendingCoor() - vertex.y;
  track.dca = TMath::Sqrt(dcaX*dcaX + dcaY*dcaY);

  // extrapolate to the end of the absorber
  TrackExtrap::extrapToZ(trackParam, -505.);
  double xAbs = trackParam.getNonBendingCoor();
  double yAbs = trackParam.getBendingCoor();
  track.rAbs = TMath::Sqrt(xAbs*xAbs + yAbs*yAbs);

  track.needExtrapToVtx = false;
}

//_________________________________________________________________________________________________
void ExtrapToMID(TrackParam& param)
{
  /// extrapolate the track parameters to the z position of the first MID chamber

  // load OCDB (only once)
  if (!ocdbLoaded) {
    if (!LoadOCDB()) {
      cout << "fail loading OCDB objects for track extrapolation" << endl;
      exit(1);
    }
    ocdbLoaded = true;
  }

  // same extrapolation method as in TrackFinder
  TrackExtrap::useExtrapV2();

  TrackExtrap::extrapToMID(param);
}

//_________________________________________________________________________________________________
void selectTracks(std::list<TrackStruct>& tracks)
{
  /// remove tracks that do not pass the selection criteria
  for (auto itTrack = tracks.begin(); itTrack != tracks.end();) {
    if (!IsSelected(*itTrack)) {
      itTrack = tracks.erase(itTrack);
    } else {
      ++itTrack;
    }
  }
}

//_________________________________________________________________________________________________
bool IsSelected(TrackStruct& track)
{
  /// apply standard track selections + pDCA

  static const double sigmaPDCA23 = 80.;
  static const double sigmaPDCA310 = 54.;
  static const double nSigmaPDCA = 6.;
  static const double relPRes = 0.0004;
  static const double slopeRes = 0.0005;

  double thetaAbs = TMath::ATan(track.rAbs/505.) * TMath::RadToDeg();
  if (thetaAbs < 2. || thetaAbs > 10.) {
    return false;
  }

  double eta = track.pxpypzm.Eta();
  if (eta < -4. || eta > -2.5) {
    return false;
  }

  double pUncorr = TMath::Sqrt(track.param.px*track.param.px + track.param.py*track.param.py + track.param.pz*track.param.pz);
  double pDCA = pUncorr * track.dca;
  double sigmaPDCA = (thetaAbs < 3) ? sigmaPDCA23 : sigmaPDCA310;
  double pTot = track.pxpypzm.P();
  double nrp = nSigmaPDCA * relPRes * pTot;
  double pResEffect = sigmaPDCA / (1. - nrp / (1. + nrp));
  double slopeResEffect = 535. * slopeRes * pTot;
  double sigmaPDCAWithRes = TMath::Sqrt(pResEffect*pResEffect + slopeResEffect*slopeResEffect);
  if (pDCA > nSigmaPDCA * sigmaPDCAWithRes) {
    return false;
  }

  // if (track.pxpypzm.Pt() < 1.5) {
  //   return false;
  // }

  // if (track.sign > 0) {
  //   return false;
  // }

  return true;
}

//_________________________________________________________________________________________________
const TrackMCHMID* FindMuon(uint32_t iMCHTrack, const std::vector<TrackMCHMID>& muonTracks)
{
  /// find the MCH-MID matched track corresponding to this MCH track
  for (const auto& muon : muonTracks) {
    if (muon.getMCHRef().getIndex() == iMCHTrack) {
      return &muon;
    }
  }
  return nullptr;
}

//_________________________________________________________________________________________________
int CompareEvents(std::list<TrackStruct>& tracks1, std::list<TrackStruct>& tracks2, double precision,
                  bool printDiff, bool printAll, std::vector<TH1*>& histos)
{
  /// compare the tracks between the 2 events

  int nDifferences(0);

  // first look for identical tracks in the 2 events
  for (auto& track1 : tracks1) {
    // find a track in the second event identical to track1 and not already matched
    auto itTrack2 = tracks2.begin();
    do {
      itTrack2 = find(itTrack2, tracks2.end(), track1);
    } while (itTrack2 != tracks2.end() && itTrack2->matchFound && ++itTrack2 != tracks2.end());
    if (itTrack2 != tracks2.end()) {
      track1.matchFound = true;
      itTrack2->matchFound = true;
      track1.matchIdentical = true;
      itTrack2->matchIdentical = true;
      // compare the track parameters
      bool areCompatible = AreTrackParamCompatible(track1.param, itTrack2->param, precision);
      bool areCovCompatible = AreTrackParamCovCompatible(track1.cov, itTrack2->cov, precision);
      if (!areCompatible || !areCovCompatible) {
        ++nDifferences;
      }
      if (printAll || (printDiff && !areCompatible)) {
        PrintResiduals(track1.param, itTrack2->param);
      }
      if (printAll || (printDiff && !areCovCompatible)) {
        PrintCovResiduals(track1.cov, itTrack2->cov);
      }
      FillResiduals(track1.param, itTrack2->param, histos);
      // compare the track parameters at MID (if any)
      if (printDiff && track1.paramAtMID.z < 0. && itTrack2->paramAtMID.z < 0.) {
        if (!AreTrackParamCompatible(track1.paramAtMID, itTrack2->paramAtMID, precision)) {
          PrintResiduals(track1.paramAtMID, itTrack2->paramAtMID);
        }
        if (!AreTrackParamCovCompatible(track1.covAtMID, itTrack2->covAtMID, precision)) {
          PrintCovResiduals(track1.covAtMID, itTrack2->covAtMID);
        }
      }
    }
  }

  // then look for similar tracks in the 2 events
  for (auto& track1 : tracks1) {
    // skip already matched tracks
    if (track1.matchFound) {
      continue;
    }
    // find a track in the second event similar to track1 and not already matched
    for (auto& track2 : tracks2) {
      if (!track2.matchFound && track2.match(track1)) {
        track1.matchFound = true;
        track2.matchFound = true;
        // if (track2.getNClustersInCommon(track1) == track2.clusters.size()) {
        //   track1.matchIdentical = true;
        //   track2.matchIdentical = true;
        // } else {
        //   track1.printClusterDifferences(track2);
        // }
        // compare the track parameters
        bool areCompatible = AreTrackParamCompatible(track1.param, track2.param, precision);
        bool areCovCompatible = AreTrackParamCovCompatible(track1.cov, track2.cov, precision);
        if (!areCompatible || !areCovCompatible) {
          ++nDifferences;
        }
        if (printAll || (printDiff && !areCompatible)) {
          PrintResiduals(track1.param, track2.param);
        }
        if (printAll || (printDiff && !areCovCompatible)) {
          PrintCovResiduals(track1.cov, track2.cov);
        }
        FillResiduals(track1.param, track2.param, histos);
        // compare the track parameters at MID (if any)
        if (printDiff && track1.paramAtMID.z < 0. && track2.paramAtMID.z < 0.) {
          if (!AreTrackParamCompatible(track1.paramAtMID, track2.paramAtMID, precision)) {
            PrintResiduals(track1.paramAtMID, track2.paramAtMID);
          }
          if (!AreTrackParamCovCompatible(track1.covAtMID, track2.covAtMID, precision)) {
            PrintCovResiduals(track1.covAtMID, track2.covAtMID);
          }
        }
        break;
      }
    }
  }

  // then print the missing tracks
  for (const auto& track1 : tracks1) {
    if (!track1.matchFound) {
      if (printDiff) {
        cout << "did not find a track in file 2 matching" << endl;
        PrintTrack(track1);
      }
      ++nDifferences;
    }
  }

  // and finally print the additional tracks
  for (const auto& track2 : tracks2) {
    if (!track2.matchFound) {
      if (printDiff) {
        cout << "did not find a track in file 1 matching" << endl;
        PrintTrack(track2);
      }
      ++nDifferences;
    }
  }

  return nDifferences;
}

//_________________________________________________________________________________________________
bool AreTrackParamCompatible(const TrackParamStruct& param1, const TrackParamStruct& param2, double precision)
{
  /// compare track parameters within precision
  return (abs(param2.x - param1.x) <= precision &&
          abs(param2.y - param1.y) <= precision &&
          abs(param2.z - param1.z) <= precision &&
          abs(param2.px - param1.px) <= precision &&
          abs(param2.py - param1.py) <= precision &&
          abs(param2.pz - param1.pz) <= precision &&
          param2.sign == param1.sign);
}

//_________________________________________________________________________________________________
bool AreTrackParamCovCompatible(const TMatrixD& cov1, const TMatrixD& cov2, double precision)
{
  /// compare track parameters covariances (if any) within precision
  if (cov1.NonZeros() == 0 || cov2.NonZeros() == 0) {
    return true;
  }
  TMatrixD diff(cov2, TMatrixD::kMinus, cov1);
  return (diff <= precision && diff >= -precision);
}

//_________________________________________________________________________________________________
int PrintEvent(const std::list<TrackStruct>& tracks)
{
  /// print all tracks in the events
  for (const auto& track : tracks) {
    PrintTrack(track);
  }
  return tracks.size();
}

//_________________________________________________________________________________________________
void PrintTrack(const TrackStruct& track)
{
  /// print the track parameters
  cout << "{x = " << track.param.x
       << ", y = " << track.param.y
       << ", z = " << track.param.z
       << ", px = " << track.param.px
       << ", py = " << track.param.py
       << ", pz = " << track.param.pz
       << ", sign = " << track.param.sign
       << "}" << endl;
}

//_________________________________________________________________________________________________
void PrintResiduals(const TrackParamStruct& param1, const TrackParamStruct& param2)
{
  /// print param2 - param1
  cout << "{dx = " << param2.x - param1.x
       << ", dy = " << param2.y - param1.y
       << ", dz = " << param2.z - param1.z
       << ", dpx = " << (param2.px - param1.px) << " (" << 100. * (param2.px - param1.px) / param1.px << "\%)"
       << ", dpy = " << (param2.py - param1.py) << " (" << 100. * (param2.py - param1.py) / param1.py << "\%)"
       << ", dpz = " << (param2.pz - param1.pz) << " (" << 100. * (param2.pz - param1.pz) / param1.pz << "\%)"
       << ", dsign = " << param2.sign - param1.sign
       << "}" << endl;
}

//_________________________________________________________________________________________________
void PrintCovResiduals(const TMatrixD& cov1, const TMatrixD& cov2)
{
  /// print cov2 - cov1
  TMatrixD diff(cov2, TMatrixD::kMinus, cov1);
  diff.Print();
}

//_________________________________________________________________________________________________
void FillResiduals(const TrackParamStruct& param1, const TrackParamStruct& param2, std::vector<TH1*>& histos)
{
  /// fill param2 - param1 histos
  if (histos.size() == 0) {
    histos.emplace_back(new TH1F("dx", "dx;dx (cm)", 20001, -1.00005, 1.00005));
    histos.emplace_back(new TH1F("dy", "dy;dy (cm)", 20001, -1.00005, 1.00005));
    histos.emplace_back(new TH1F("dz", "dz;dz (cm)", 20001, -1.00005, 1.00005));
    histos.emplace_back(new TH1F("dpx", "dpx;dpx (GeV/c)", 20001, -1.00005, 1.00005));
    histos.emplace_back(new TH1F("dpy", "dpy;dpy (GeV/c)", 20001, -1.00005, 1.00005));
    histos.emplace_back(new TH1F("dpz", "dpz;dpz (GeV/c)", 20001, -1.00005, 1.00005));
    histos.emplace_back(new TH2F("dpxvspx", "dpxvspx;px (GeV/c);dpx/px (\%)", 2000, 0., 20., 2001, -10.005, 10.005));
    histos.emplace_back(new TH2F("dpyvspy", "dpyvspy;py (GeV/c);dpy/py (\%)", 2000, 0., 20., 2001, -10.005, 10.005));
    histos.emplace_back(new TH2F("dpzvspz", "dpzvspz;pz (GeV/c);dpz/pz (\%)", 2000, 0., 200., 2001, -10.005, 10.005));
    histos.emplace_back(new TH2F("dslopexvsp", "dslopexvsp;p (GeV/c);dslopex", 2000, 0., 200., 2001, -0.0010005, 0.0010005));
    histos.emplace_back(new TH2F("dslopeyvsp", "dslopeyvsp;p (GeV/c);dslopey", 2000, 0., 200., 2001, -0.0010005, 0.0010005));
    histos.emplace_back(new TH2F("dpvsp", "dpvsp;p (GeV/c);dp/p (\%)", 2000, 0., 200., 2001, -10.005, 10.005));
    for (auto h : histos) {
      h->SetDirectory(0);
    }
  }
  double p1 = TMath::Sqrt(param1.px*param1.px + param1.py*param1.py + param1.pz*param1.pz);
  double p2 = TMath::Sqrt(param2.px*param2.px + param2.py*param2.py + param2.pz*param2.pz);
  histos[0]->Fill(param2.x - param1.x);
  histos[1]->Fill(param2.y - param1.y);
  histos[2]->Fill(param2.z - param1.z);
  histos[3]->Fill(param2.px - param1.px);
  histos[4]->Fill(param2.py - param1.py);
  histos[5]->Fill(param2.pz - param1.pz);
  histos[6]->Fill(TMath::Abs(param1.px), 100. * (param2.px - param1.px) / param1.px);
  histos[7]->Fill(TMath::Abs(param1.py), 100. * (param2.py - param1.py) / param1.py);
  histos[8]->Fill(TMath::Abs(param1.pz), 100. * (param2.pz - param1.pz) / param1.pz);
  histos[9]->Fill(p1, param2.px / param2.pz - param1.px / param1.pz);
  histos[10]->Fill(p1, param2.py / param2.pz - param1.py / param1.pz);
  histos[11]->Fill(p1, 100. * (p2 - p1) / p1);
}

//_________________________________________________________________________________________________
void DrawResiduals(std::vector<TH1*>& histos)
{
  /// draw param2 - param1 histos
  int nPadsx = (histos.size()+2)/3;
  TCanvas* c = new TCanvas("residuals", "residuals", 10, 10, nPadsx*300, 900);
  c->Divide(nPadsx,3);
  int i(0);
  for (const auto& h : histos) {
    c->cd((i % 3) * nPadsx + i / 3 + 1);
    if (dynamic_cast<TH1F*>(h) == nullptr) {
      h->Draw("colz");
      gPad->SetLogz();
    } else {
      h->Draw();
      gPad->SetLogy();
    }
    ++i;
  }
}

//_________________________________________________________________________________________________
void CreateHistosAtVertex(std::vector<TH1*>& histos, const char* extension)
{
  /// create single muon and dimuon histograms at vertex
  if (histos.size() == 0) {
    histos.emplace_back(new TH1F(Form("pT%s", extension), "pT;p_{T} (GeV/c)", 300, 0., 30.));
    histos.emplace_back(new TH1F(Form("eta%s", extension), "eta;eta", 200, -4.5, -2.));
    histos.emplace_back(new TH1F(Form("phi%s", extension), "phi;phi", 360, 0., 360.));
    histos.emplace_back(new TH1F(Form("rAbs%s", extension), "rAbs;R_{abs} (cm)", 1000, 0., 100.));
    histos.emplace_back(new TH1F(Form("p%s", extension), "p;p (GeV/c)", 300, 0., 300.));
    histos.emplace_back(new TH1F(Form("dca%s", extension), "DCA;DCA (cm)", 500, 0., 500.));
    histos.emplace_back(new TH1F(Form("pDCA23%s", extension), "pDCA for #theta_{abs} < 3#circ;pDCA (GeV.cm/c)", 2500, 0., 5000.));
    histos.emplace_back(new TH1F(Form("pDCA310%s", extension), "pDCA for #theta_{abs} > 3#circ;pDCA (GeV.cm/c)", 2500, 0., 5000.));
    histos.emplace_back(new TH1F(Form("nClusters%s", extension), "number of clusters per track;n_{clusters}", 20, 0., 20.));
    histos.emplace_back(new TH1F(Form("chi2%s", extension), "normalized #chi^{2};#chi^{2} / ndf", 500, 0., 50.));
    histos.emplace_back(new TH1F(Form("matchChi2%s", extension), "normalized matched #chi^{2};#chi^{2} / ndf", 160, 0., 16.));
    histos.emplace_back(new TH1F(Form("mass%s", extension), "#mu^{+}#mu^{-} invariant mass;mass (GeV/c^{2})", 1600, 0., 20.));
  }
  for (auto h : histos) {
    h->SetDirectory(0);
  }
}

//_________________________________________________________________________________________________
void FillHistosAtVertex(const std::list<TrackStruct>& tracks, std::vector<TH1*>& histos)
{
  /// fill single muon and dimuon histograms at vertex
  for (auto itTrack1 = tracks.begin(); itTrack1 != tracks.end(); ++itTrack1) {
    FillHistosMuAtVertex(*itTrack1, histos);
    for (auto itTrack2 = std::next(itTrack1); itTrack2 != tracks.end(); ++itTrack2) {
      FillHistosDimuAtVertex(*itTrack1, *itTrack2, histos);
    }
  }
}

//_________________________________________________________________________________________________
void FillHistosMuAtVertex(const TrackStruct& track, std::vector<TH1*>& histos)
{
  /// fill single muon histograms at vertex
/*
  if (track.pxpypzm.Pt() < 1.) {
    return;
  }
*/
  double thetaAbs = TMath::ATan(track.rAbs/505.) * TMath::RadToDeg();
  double pUncorr = TMath::Sqrt(track.param.px*track.param.px + track.param.py*track.param.py + track.param.pz*track.param.pz);
  double pDCA = pUncorr * track.dca;
  
  histos[0]->Fill(track.pxpypzm.Pt());
  histos[1]->Fill(track.pxpypzm.Eta());
  histos[2]->Fill(180. + atan2(-track.pxpypzm.Py(), -track.pxpypzm.Px()) / pi() * 180.);
  histos[3]->Fill(track.rAbs);
  histos[4]->Fill(track.pxpypzm.P());
  histos[5]->Fill(track.dca);
  if (thetaAbs < 3) {
    histos[6]->Fill(pDCA);
  } else {
    histos[7]->Fill(pDCA);
  }
  histos[8]->Fill(track.clusters.size());
  histos[9]->Fill(track.getChi2N2());
  histos[10]->Fill(track.matchedChi2);
}

//_________________________________________________________________________________________________
void FillHistosDimuAtVertex(const TrackStruct& track1, const TrackStruct& track2, std::vector<TH1*>& histos)
{
  /// fill dimuon histograms at vertex
/*
  if (track1.pxpypzm.Pt() < 1. || track2.pxpypzm.Pt() < 1.) {
    return;
  }
*/
  if (track1.sign * track2.sign < 0) {
    PxPyPzMVector dimu = track1.pxpypzm + track2.pxpypzm;
    histos[11]->Fill(dimu.M());
  }
}

//_________________________________________________________________________________________________
void DrawHistosAtVertex(std::vector<TH1*> histos[2])
{
  /// Draw histograms at vertex and differences between the 2 inputs
  
  // find the optimal number of pads
  int nPadsx(1), nPadsy(1);
  while ((int)histos[0].size() > nPadsx*nPadsy) {
    if (nPadsx == nPadsy) {
      ++nPadsx;
    } else {
      ++nPadsy;
    }
  }

  // draw histograms
  TCanvas* cHist = new TCanvas("histos", "histos", 10, 10, TMath::Max(nPadsx * 300, 1200), TMath::Max(nPadsy * 300, 900));
  cHist->Divide(nPadsx, nPadsy);
  for (int i = 0; i < (int)histos[0].size(); ++i) {
    cHist->cd(i + 1);
    gPad->SetLogy();
    histos[0][i]->SetStats(false);
    histos[0][i]->SetLineColor(4);
    histos[0][i]->Draw();
    histos[1][i]->SetLineColor(2);
    histos[1][i]->Draw("same");
  }

  // add a legend
  TLegend* lHist = new TLegend(0.5, 0.65, 0.9, 0.8);
  lHist->SetFillStyle(0);
  lHist->SetBorderSize(0);
  lHist->AddEntry(histos[0][0], Form("%g tracks in file 1", histos[0][0]->GetEntries()), "l");
  lHist->AddEntry(histos[1][0], Form("%g tracks in file 2", histos[1][0]->GetEntries()), "l");
  cHist->cd(1);
  lHist->Draw("same");

  // draw differences
  TCanvas* cDiff = new TCanvas("differences", "histos2 - histos1", 10, 10, TMath::Max(nPadsx * 300, 1200), TMath::Max(nPadsy * 300, 900));
  cDiff->Divide(nPadsx, nPadsy);
  for (int i = 0; i < (int)histos[0].size(); ++i) {
    cDiff->cd(i + 1);
    TH1F* hDiff = static_cast<TH1F*>(histos[1][i]->Clone());
    hDiff->Add(histos[0][i], -1.);
    hDiff->SetStats(false);
    hDiff->SetLineColor(2);
    hDiff->Draw();
  }

  // draw ratios
  TCanvas* cRat = new TCanvas("ratios", "histos2 / histos1", 10, 10, TMath::Max(nPadsx * 300, 1200), TMath::Max(nPadsy * 300, 900));
  cRat->Divide(nPadsx, nPadsy);
  for (int i = 0; i < (int)histos[0].size(); ++i) {
    cRat->cd(i + 1);
    TH1F* hRat = static_cast<TH1F*>(histos[1][i]->Clone());
    hRat->Divide(histos[0][i]);
    hRat->SetStats(false);
    hRat->SetLineColor(2);
    hRat->Draw();
  }
}

//_________________________________________________________________________________________________
void FillComparisonsAtVertex(std::list<TrackStruct>& tracks1, std::list<TrackStruct>& tracks2, std::vector<TH1*> histos[5])
{
  /// fill comparison histograms at vertex

  for (auto itTrack21 = tracks2.begin(); itTrack21 != tracks2.end(); ++itTrack21) {

    // fill histograms for identical, similar (from file 2) and additional muons
    if (itTrack21->matchIdentical) {
      FillHistosMuAtVertex(*itTrack21, histos[0]);
    } else if (itTrack21->matchFound) {
      FillHistosMuAtVertex(*itTrack21, histos[2]);
    } else {
      FillHistosMuAtVertex(*itTrack21, histos[3]);
    }

    for (auto itTrack22 = std::next(itTrack21); itTrack22 != tracks2.end(); ++itTrack22) {

      // fill histograms for identical, similar (from file 2) and additional dimuons
      if (itTrack21->matchIdentical && itTrack22->matchIdentical) {
        FillHistosDimuAtVertex(*itTrack21, *itTrack22, histos[0]);
      } else if (itTrack21->matchFound && itTrack22->matchFound) {
        FillHistosDimuAtVertex(*itTrack21, *itTrack22, histos[2]);
      } else {
        FillHistosDimuAtVertex(*itTrack21, *itTrack22, histos[3]);
      }
    }
  }

  for (auto itTrack11 = tracks1.begin(); itTrack11 != tracks1.end(); ++itTrack11) {

    // fill histograms for missing and similar (from file 1) muons
    if (!itTrack11->matchFound) {
      FillHistosMuAtVertex(*itTrack11, histos[4]);
    } else if (!itTrack11->matchIdentical) {
      FillHistosMuAtVertex(*itTrack11, histos[1]);
    }

    for (auto itTrack12 = std::next(itTrack11); itTrack12 != tracks1.end(); ++itTrack12) {

      // fill histograms for missing and similar (from file 1) dimuons
      if (!itTrack11->matchFound || !itTrack12->matchFound) {
        FillHistosDimuAtVertex(*itTrack11, *itTrack12, histos[4]);
      } else if (!itTrack11->matchIdentical || !itTrack12->matchIdentical) {
        FillHistosDimuAtVertex(*itTrack11, *itTrack12, histos[1]);
      }
    }
  }
/*
  // fill histograms for missing dimuons in both files
  for (const auto& track1 : tracks1) {
    if (!track1.matchFound) {
      for (const auto& track2 : tracks2) {
        if (!track2.matchFound) {
          FillHistosDimuAtVertex(track1, track2, histos[4]);
        }
      }
    }
  }*/
}

//_________________________________________________________________________________________________
void DrawComparisonsAtVertex(std::vector<TH1*> histos[5])
{
  /// draw comparison histograms at vertex

  // find the optimal number of pads
  int nPadsx(1), nPadsy(1);
  while ((int)histos[0].size() > nPadsx*nPadsy) {
    if (nPadsx == nPadsy) {
      ++nPadsx;
    } else {
      ++nPadsy;
    }
  }

  // draw histograms
  TCanvas* cHist = new TCanvas("comparisons", "comparisons", 10, 10, TMath::Max(nPadsx * 300, 1200), TMath::Max(nPadsy * 300, 900));
  cHist->Divide(nPadsx, nPadsy);
  for (int i = 0; i < (int)histos[0].size(); ++i) {
    cHist->cd(i + 1);
    gPad->SetLogy();
    histos[0][i]->SetStats(false);
    histos[0][i]->SetLineColor(1);
    histos[0][i]->SetMinimum(0.5);
    histos[0][i]->Draw();
    histos[1][i]->SetLineColor(4);
    histos[1][i]->Draw("same");
    histos[2][i]->SetLineColor(877);
    histos[2][i]->Draw("same");
    histos[3][i]->SetLineColor(3);
    histos[3][i]->Draw("same");
    histos[4][i]->SetLineColor(2);
    histos[4][i]->Draw("same");
  }

  // add a legend
  TLegend* lHist = new TLegend(0.5, 0.5, 0.9, 0.8);
  lHist->SetFillStyle(0);
  lHist->SetBorderSize(0);
  lHist->AddEntry(histos[0][0], Form("%g tracks identical", histos[0][0]->GetEntries()), "l");
  lHist->AddEntry(histos[1][0], Form("%g tracks similar (1)", histos[1][0]->GetEntries()), "l");
  lHist->AddEntry(histos[2][0], Form("%g tracks similar (2)", histos[2][0]->GetEntries()), "l");
  lHist->AddEntry(histos[3][0], Form("%g tracks additional", histos[3][0]->GetEntries()), "l");
  lHist->AddEntry(histos[4][0], Form("%g tracks missing", histos[4][0]->GetEntries()), "l");
  cHist->cd(1);
  lHist->Draw("same");
}

//_________________________________________________________________________________________________
void PrintStat(TString what, int n1, int n2)
{
  /// print given statistics and the relative difference (in %)
  /// the uncertainty on the difference is the normal approximation of the binomial error
  if (n1 == 0) {
    printf("%s | %8d | %8d | %7s ± %4s %%\n", what.Data(), n1, n2, "nan", "nan");
  } else {
    double eff = double(n2) / n1;
    double diff = 100. * (eff - 1.);
    double err = 100. * TMath::Max(1. / n1, TMath::Sqrt(eff * TMath::Abs(1. - eff) / n1));
    printf("%s | %8d | %8d | %7.2f ± %4.2f %%\n", what.Data(), n1, n2, diff, err);
  }
}

//_________________________________________________________________________________________________
void PrintStat(int nTracksAll[2], int nTracksMatch[2], std::vector<TH1*> histos[2])
{
  /// print some statistics

  printf("\n");
  printf("-------------------------------------------------------\n");
  printf("selection      |  file 1  |  file 2  |       diff\n");
  printf("-------------------------------------------------------\n");

  PrintStat("all           ", nTracksAll[0], nTracksAll[1]);

  PrintStat("matched       ", nTracksMatch[0], nTracksMatch[1]);

  PrintStat("selected      ", histos[0][0]->GetEntries(), histos[1][0]->GetEntries());

  double pTRange[6] = {0., 0.5, 1., 2., 4., 1000.};
  for (int i = 0; i < 5; ++i) {
    TString selection = (i == 0) ? TString::Format("pT < %.1f GeV/c", pTRange[1])
                                 : ((i == 4) ? TString::Format("pT > %.1f      ", pTRange[4])
                                             : TString::Format("%.1f < pT < %.1f", pTRange[i], pTRange[i + 1]));
    int n1 = histos[0][0]->Integral(histos[0][0]->GetXaxis()->FindBin(pTRange[i] + 0.01),
                                    histos[0][0]->GetXaxis()->FindBin(pTRange[i + 1] - 0.01));
    int n2 = histos[1][0]->Integral(histos[1][0]->GetXaxis()->FindBin(pTRange[i] + 0.01),
                                    histos[1][0]->GetXaxis()->FindBin(pTRange[i + 1] - 0.01));
    PrintStat(selection, n1, n2);
  }

  double pRange[6] = {0., 5., 10., 20., 40., 10000.};
  for (int i = 0; i < 5; ++i) {
    TString selection = (i == 0) ? TString::Format("p < %02.0f GeV/c  ", pRange[1])
                                 : ((i == 4) ? TString::Format("p > %02.0f        ", pRange[4])
                                             : TString::Format("%02.0f < p < %02.0f   ", pRange[i], pRange[i + 1]));
    int n1 = histos[0][4]->Integral(histos[0][4]->GetXaxis()->FindBin(pRange[i] + 0.01),
                                    histos[0][4]->GetXaxis()->FindBin(pRange[i + 1] - 0.01));
    int n2 = histos[1][4]->Integral(histos[1][4]->GetXaxis()->FindBin(pRange[i] + 0.01),
                                    histos[1][4]->GetXaxis()->FindBin(pRange[i + 1] - 0.01));
    PrintStat(selection, n1, n2);
  }

  printf("-------------------------------------------------------\n");
  printf("\n");
}
