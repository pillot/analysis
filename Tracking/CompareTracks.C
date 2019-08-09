#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>

#include <TMath.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TH1F.h>
#include <TH2F.h>

#include "MCHBase/ClusterBlock.h"
#include "MCHBase/TrackBlock.h"

using namespace std;
using namespace o2::mch;

//_________________________________________________________________________________________________
struct TrackStruct {
  TrackParamStruct param{};
  std::vector<ClusterStruct> clusters{};
  bool matchFound = false;

  bool operator==(const TrackStruct& track) const
  {
    /// tracks are considered identical when they share exactly the same clusters
    if (this->clusters.size() != track.clusters.size()) {
      return false;
    }
    for (size_t iCl = 0; iCl != this->clusters.size(); ++iCl) {
      if (this->clusters[iCl].uid != track.clusters[iCl].uid) {
        return false;
      }
    }
    return true;
  }
};

//_________________________________________________________________________________________________
void ReadNextEvent(ifstream& inFile, int& event, std::vector<TrackStruct>& tracks);
void ReadTrack(ifstream& inFile, TrackStruct& track);
int CompareEvents(std::vector<TrackStruct>& tracks1, std::vector<TrackStruct>& tracks2, double precision, bool printAll, std::vector<TH1*>& histos);
bool AreTrackParamCompatible(const TrackParamStruct& param1, const TrackParamStruct& param2, double precision);
int PrintEvent(const std::vector<TrackStruct>& tracks);
void PrintTrack(const TrackStruct& track);
void PrintResiduals(const TrackParamStruct& param1, const TrackParamStruct& param2);
void FillResiduals(const TrackParamStruct& param1, const TrackParamStruct& param2, std::vector<TH1*>& histos);
void DrawResiduals(std::vector<TH1*>& histos);

//_________________________________________________________________________________________________
void CompareTracks(string inFileName1, string inFileName2, double precision = 1.e-4, bool printAll = false)
{
  /// Compare the tracks stored in the 2 binary files
  
  // open files
  ifstream inFile1(inFileName1,ios::binary);
  ifstream inFile2(inFileName2,ios::binary);
  if (!inFile1.is_open() || !inFile2.is_open()) {
    cout << "fail opening files" << endl;
    return;
  }

  int event1(0), event2(0);
  std::vector<TrackStruct> tracks1{};
  std::vector<TrackStruct> tracks2{};
  bool readNextEvent1(true), readNextEvent2(true);
  int nDifferences(0);
  std::vector<TH1*> histos{};
 
  while (true) {

    if (readNextEvent1) {
      ReadNextEvent(inFile1, event1, tracks1);
    }

    if (readNextEvent2) {
      ReadNextEvent(inFile2, event2, tracks2);
    }

    // reaching end of both files
    if (event1 < 0 && event2 < 0) {
      break;
    }

    if (event1 == event2) {
      // reading the same event --> we can compare tracks
      nDifferences += CompareEvents(tracks1, tracks2, precision, printAll, histos);
      readNextEvent1 = true;
      readNextEvent2 = true;
    } else if (event2 < 0 || (event1 >= 0 && event1 < event2)) {
      // event 2 > event 1 or reaching end of file 2
      if (tracks1.size() > 0) {
        cout << "tracks in event " << event1 << " are missing in file 2" << endl;
        nDifferences += PrintEvent(tracks1);
      }
      readNextEvent1 = true;
      readNextEvent2 = false;
    } else {
      // event 1 > event 2 or reaching end of file 1
      if (tracks2.size() > 0) {
        cout << "tracks in event " << event2 << " are missing in file 1" << endl;
        nDifferences += PrintEvent(tracks2);
      }
      readNextEvent1 = false;
      readNextEvent2 = true;
    }
  }
  
  cout << "number of different tracks = " << nDifferences << endl;
  DrawResiduals(histos);

  inFile1.close();
  inFile2.close();

}

//_________________________________________________________________________________________________
void ReadNextEvent(ifstream& inFile, int& event, std::vector<TrackStruct>& tracks)
{
  /// read the next event in the input file

  tracks.clear();
  
  if (!inFile.read(reinterpret_cast<char*>(&event), sizeof(int))) {
    // reaching end of file
    event = -1;
    return;
  }
  
  int size(0);
  inFile.read(reinterpret_cast<char*>(&size), sizeof(int));
  if (size == 0) {
    // empty event
    return;
  }

  // read the tracks
  int nTracks(0);
  inFile.read(reinterpret_cast<char*>(&nTracks), sizeof(int));
  tracks.reserve(nTracks);
  TrackStruct track{};
  for (Int_t iTrack = 0; iTrack < nTracks; ++iTrack) {
    ReadTrack(inFile, track);
    tracks.push_back(track);
  }
}

//_________________________________________________________________________________________________
void ReadTrack(ifstream& inFile, TrackStruct& track)
{
  /// read one track from the input file
  inFile.read(reinterpret_cast<char*>(&(track.param)), sizeof(TrackParamStruct));
  int nClusters(0);
  inFile.read(reinterpret_cast<char*>(&nClusters), sizeof(int));
  track.clusters.resize(nClusters);
  for (Int_t iCl = 0; iCl < nClusters; ++iCl) {
    inFile.read(reinterpret_cast<char*>(&(track.clusters[iCl])), sizeof(ClusterStruct));
  }
}

//_________________________________________________________________________________________________
int CompareEvents(std::vector<TrackStruct>& tracks1, std::vector<TrackStruct>& tracks2, double precision, bool printAll, std::vector<TH1*>& histos)
{
  /// compare the tracks between the 2 events

  int nDifferences(0);

  for (const auto& track1 : tracks1) {
    // find a track in the second event corresponding to track1 and not already matched
    auto itTrack2 = tracks2.begin();
    do {
      itTrack2 = find(itTrack2, tracks2.end(), track1);
    } while (itTrack2 != tracks2.end() && itTrack2->matchFound && ++itTrack2 != tracks2.end());
    if (itTrack2 != tracks2.end()) {
      itTrack2->matchFound = true;
      // compare the track parameters
      bool areCompatible = AreTrackParamCompatible(track1.param, itTrack2->param, precision);
      if (!areCompatible) {
        ++nDifferences;
      }
      if (printAll || !areCompatible) {
        PrintResiduals(track1.param, itTrack2->param);
      }
      FillResiduals(track1.param, itTrack2->param, histos);
    } else {
      cout << "did not find a track in file 2 matching" << endl;
      PrintTrack(track1);
      ++nDifferences;
    }
  }

  for (const auto& track2 : tracks2) {
    if (!track2.matchFound) {
      cout << "did not find a track in file 1 matching" << endl;
      PrintTrack(track2);
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
int PrintEvent(const std::vector<TrackStruct>& tracks)
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
void FillResiduals(const TrackParamStruct& param1, const TrackParamStruct& param2, std::vector<TH1*>& histos)
{
  /// fill param2 - param1 histos
  if (histos.size() == 0) {
    histos.emplace_back(new TH1F("dx", "dx;dx (cm)", 20000, -1., 1.));
    histos.emplace_back(new TH1F("dy", "dy;dy (cm)", 20000, -1., 1.));
    histos.emplace_back(new TH1F("dz", "dz;dz (cm)", 20000, -1., 1.));
    histos.emplace_back(new TH1F("dpx", "dpx;dpx (GeV/c)", 20000, -1., 1.));
    histos.emplace_back(new TH1F("dpy", "dpy;dpy (GeV/c)", 20000, -1., 1.));
    histos.emplace_back(new TH1F("dpz", "dpz;dpz (GeV/c)", 20000, -1., 1.));
    histos.emplace_back(new TH2F("dpxvspx", "dpxvspx;px (GeV/c);dpx (\%)", 2000, 0., 20., 2000, -10., 10.));
    histos.emplace_back(new TH2F("dpyvspy", "dpyvspy;py (GeV/c);dpy (\%)", 2000, 0., 20., 2000, -10., 10.));
    histos.emplace_back(new TH2F("dpzvspz", "dpzvspz;pz (GeV/c);dpz (\%)", 2000, 0., 200., 2000, -10., 10.));
    histos.emplace_back(new TH2F("dslopexvsp", "dslopexvsp;p (GeV/c);dslopex", 2000, 0., 200., 1000, 0., 0.001));
    histos.emplace_back(new TH2F("dslopeyvsp", "dslopeyvsp;p (GeV/c);dslopey", 2000, 0., 200., 1000, 0., 0.001));
    histos.emplace_back(new TH2F("dpvsp", "dpvsp;p (GeV/c);dp (\%)", 2000, 0., 200., 1000, 0., 10.));
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
  histos[9]->Fill(p1, TMath::Abs(param2.px / param2.pz - param1.px / param1.pz));
  histos[10]->Fill(p1, TMath::Abs(param2.py / param2.pz - param1.py / param1.pz));
  histos[11]->Fill(p1, 100. * TMath::Abs(p2 - p1) / p1);
}

//_________________________________________________________________________________________________
void DrawResiduals(std::vector<TH1*>& histos)
{
  /// draw param2 - param1 histos
  gStyle->SetOptStat(111111);
  int nPadsx = (histos.size()+2)/3;
  TCanvas* c = new TCanvas("residuals", "residuals", 10, 10, nPadsx*300, 900);
  c->Divide(nPadsx,3);
  int i(0);
  for (const auto& h : histos) {
    c->cd((i%3)*nPadsx + i/3 + 1);
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
