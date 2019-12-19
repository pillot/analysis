#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>

#include <TMath.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TDatabasePDG.h>
#include <Math/Vector4D.h>

#include "AliCDBManager.h"
#include "AliGeomManager.h"
#include "AliMUONCDB.h"
#include "AliMUONTrackExtrap.h"
#include "AliMUONTrackParam.h"

#include "MCHBase/ClusterBlock.h"
#include "MCHBase/TrackBlock.h"

using namespace std;
using namespace o2::mch;
using namespace ROOT::Math;

static const double muMass = TDatabasePDG::Instance()->GetParticle("mu-")->Mass();

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
  TrackParamStruct param{};
  std::vector<ClusterStruct> clusters{};
  bool matchFound = false;
  bool matchIdentical = false;

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

  bool match(const TrackStruct& track) const
  {
    /// Try to match this track with the given track. Matching conditions:
    /// - more than 50% of clusters from one of the two tracks matched with clusters from the other
    /// - at least 1 cluster matched before and 1 cluster matched after the dipole

    size_t nMatchClusters(0);
    bool matchCluster[10] = {false, false, false, false, false, false, false, false, false, false};

    for(const auto& cluster1 : this->clusters) {
      for(const auto& cluster2 : track.clusters) {
        if (cluster1.uid == cluster2.uid) {
          matchCluster[cluster1.getChamberId()] = true;
          ++nMatchClusters;
          break;
        }
      }
    }
  
    return ((matchCluster[0] || matchCluster[1] || matchCluster[2] || matchCluster[3]) &&
            (matchCluster[6] || matchCluster[7] || matchCluster[8] || matchCluster[9]) &&
            (2 * nMatchClusters > this->clusters.size() || 2 * nMatchClusters > track.clusters.size()));
  }
};

//_________________________________________________________________________________________________
int ReadNextEvent(ifstream& inFile, std::vector<VertexStruct>& vertices);
int ReadNextEvent(ifstream& inFile, int version, std::vector<VertexStruct>& vertices, bool selectTracks, std::vector<TrackStruct>& tracks);
void ReadTrack(ifstream& inFile, TrackStruct& track, int version);
bool LoadOCDB();
void ExtrapToVertex(TrackStruct& track, VertexStruct& vertex);
bool IsSelected(TrackStruct& track);
int CompareEvents(std::vector<TrackStruct>& tracks1, std::vector<TrackStruct>& tracks2, double precision, bool printAll, std::vector<TH1*>& histos);
bool AreTrackParamCompatible(const TrackParamStruct& param1, const TrackParamStruct& param2, double precision);
int PrintEvent(const std::vector<TrackStruct>& tracks);
void PrintTrack(const TrackStruct& track);
void PrintResiduals(const TrackParamStruct& param1, const TrackParamStruct& param2);
void FillResiduals(const TrackParamStruct& param1, const TrackParamStruct& param2, std::vector<TH1*>& histos);
void DrawResiduals(std::vector<TH1*>& histos);
void CreateHistosAtVertex(std::vector<TH1*>& histos, const char* extension);
void FillHistosAtVertex(const std::vector<TrackStruct>& tracks, std::vector<TH1*>& histos);
void FillHistosMuAtVertex(const TrackStruct& track, std::vector<TH1*>& histos);
void FillHistosDimuAtVertex(const TrackStruct& track1, const TrackStruct& track2, std::vector<TH1*>& histos);
void DrawHistosAtVertex(std::vector<TH1*> histos[2]);
void FillComparisonsAtVertex(std::vector<TrackStruct>& tracks1, std::vector<TrackStruct>& tracks2, std::vector<TH1*> histos[4]);
void DrawComparisonsAtVertex(std::vector<TH1*> histos[4]);

//_________________________________________________________________________________________________
void CompareTracks(string inFileName1, int versionFile1, string inFileName2, int versionFile2,
                   string vtxFileName = "", double precision = 1.e-4, bool selectTracks = false, bool printAll = false)
{
  /// Compare the tracks stored in the 2 binary files
  /// file version 1: param at 1st cluster + clusters
  /// file version 2: param at 1st cluster + chi2 + clusters
  /// file version 3: param at vertex + dca +rAbs + chi2 + param at 1st cluster + clusters

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
    if ((versionFile1 < 3 || versionFile2 < 3) && !LoadOCDB()) {
      cout << "fail loading OCDB objects for track extrapolation" << endl;
      return;
    }
  }

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
  std::vector<TH1*> residualsAtFirstCluster{};
  std::vector<TH1*> comparisonsAtVertex[4] = {{}, {}, {}, {}};
  CreateHistosAtVertex(comparisonsAtVertex[0], "identical");
  CreateHistosAtVertex(comparisonsAtVertex[1], "similar");
  CreateHistosAtVertex(comparisonsAtVertex[2], "additional");
  CreateHistosAtVertex(comparisonsAtVertex[3], "missing");
  std::vector<TH1*> histosAtVertex[2] = {{}, {}};
  CreateHistosAtVertex(histosAtVertex[0], "1");
  CreateHistosAtVertex(histosAtVertex[1], "2");
 
  while (true) {

    if (readNextEvent1) {
      event1 = ReadNextEvent(inFile1, versionFile1, vertices, selectTracks, tracks1);
      FillHistosAtVertex(tracks1, histosAtVertex[0]);
    }

    if (readNextEvent2) {
      event2 = ReadNextEvent(inFile2, versionFile2, vertices, selectTracks, tracks2);
      FillHistosAtVertex(tracks2, histosAtVertex[1]);
    }

    // reaching end of both files
    if (event1 < 0 && event2 < 0) {
      break;
    }

    if (event1 == event2) {
      // reading the same event --> we can compare tracks
      int nDiff = CompareEvents(tracks1, tracks2, precision, printAll, residualsAtFirstCluster);
      if (nDiff > 0) {
        cout << "--> " << nDiff << " differences found in event " << event1 << endl;
        nDifferences += nDiff;
      }
      FillComparisonsAtVertex(tracks1, tracks2, comparisonsAtVertex);
      readNextEvent1 = true;
      readNextEvent2 = true;
    } else if (event2 < 0 || (event1 >= 0 && event1 < event2)) {
      // event 2 > event 1 or reaching end of file 2
      if (tracks1.size() > 0) {
        cout << "tracks in event " << event1 << " are missing in file 2" << endl;
        nDifferences += PrintEvent(tracks1);
        FillHistosAtVertex(tracks1, comparisonsAtVertex[3]);
      }
      readNextEvent1 = true;
      readNextEvent2 = false;
    } else {
      // event 1 > event 2 or reaching end of file 1
      if (tracks2.size() > 0) {
        cout << "tracks in event " << event2 << " are missing in file 1" << endl;
        nDifferences += PrintEvent(tracks2);
        FillHistosAtVertex(tracks2, comparisonsAtVertex[2]);
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
int ReadNextEvent(ifstream& inFile, int version, std::vector<VertexStruct>& vertices, bool selectTracks, std::vector<TrackStruct>& tracks)
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

  // get the vertex corresponding to this event or use (0.,0.,0.)
  VertexStruct vertex = {0., 0., 0.};
  if (!vertices.empty()) {
    if (event < (int)vertices.size()) {
      vertex = vertices[event];
    } else {
      cout << "missing vertex for event " << event << endl;
    }
  }

  // read the tracks
  int nTracks(0);
  inFile.read(reinterpret_cast<char*>(&nTracks), sizeof(int));
  tracks.reserve(nTracks);
  TrackStruct track{};
  for (Int_t iTrack = 0; iTrack < nTracks; ++iTrack) {
    ReadTrack(inFile, track, version);
    if (version < 3) {
      ExtrapToVertex(track, vertex);
    }
    if (!selectTracks || IsSelected(track)) {
      tracks.push_back(track);
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
    inFile.read(reinterpret_cast<char*>(&(track.clusters[iCl])), sizeof(ClusterStruct));
  }
}

//_________________________________________________________________________________________________
bool LoadOCDB()
{
  /// access OCDB and prepare track extrapolation to vertex

  AliCDBManager* man = AliCDBManager::Instance();
  if (gSystem->AccessPathName("OCDB.root", kFileExists)==0) {
    man->SetDefaultStorage("local:///dev/null");
    man->SetSnapshotMode("OCDB.root");
  } else {
    man->SetDefaultStorage("local://./OCDB");
  }
  man->SetRun(295584);

  if (!AliMUONCDB::LoadField()) {
    return false;
  }
  AliMUONTrackExtrap::SetField();

  AliGeomManager::LoadGeometry();
  if (!AliGeomManager::GetGeometry() || !AliGeomManager::ApplyAlignObjsFromCDB("MUON")) {
    return false;
  }

  return true;
}

//_________________________________________________________________________________________________
void ExtrapToVertex(TrackStruct& track, VertexStruct& vertex)
{
  /// compute the track parameters at vertex, at DCA and at the end of the absorber

  if (!AliGeomManager::GetGeometry()) {
    return;
  }

  // convert parameters at first cluster in MUON format
  AliMUONTrackParam trackParam;
  trackParam.SetNonBendingCoor(track.param.x);
  trackParam.SetBendingCoor(track.param.y);
  trackParam.SetZ(track.param.z);
  trackParam.SetNonBendingSlope(track.param.px / track.param.pz);
  trackParam.SetBendingSlope(track.param.py / track.param.pz);
  trackParam.SetInverseBendingMomentum(track.param.sign/TMath::Sqrt(track.param.py*track.param.py + track.param.pz*track.param.pz));

  // extrapolate to vertex
  AliMUONTrackParam trackParamAtVertex(trackParam);
  AliMUONTrackExtrap::ExtrapToVertex(&trackParamAtVertex, vertex.x, vertex.y, vertex.z, 0., 0.);
  track.pxpypzm.SetPx(trackParamAtVertex.Px());
  track.pxpypzm.SetPy(trackParamAtVertex.Py());
  track.pxpypzm.SetPz(trackParamAtVertex.Pz());
  track.pxpypzm.SetM(muMass);
  track.sign = trackParamAtVertex.GetCharge();

  // extrapolate to DCA
  AliMUONTrackParam trackParamAtDCA(trackParam);
  AliMUONTrackExtrap::ExtrapToVertexWithoutBranson(&trackParamAtDCA, vertex.z);
  double dcaX = trackParamAtDCA.GetNonBendingCoor() - vertex.x;
  double dcaY = trackParamAtDCA.GetBendingCoor() - vertex.y;
  track.dca = TMath::Sqrt(dcaX*dcaX + dcaY*dcaY);
 
  // extrapolate to the end of the absorber
  AliMUONTrackExtrap::ExtrapToZ(&trackParam, -505.);
  double xAbs = trackParam.GetNonBendingCoor();
  double yAbs = trackParam.GetBendingCoor();
  track.rAbs = TMath::Sqrt(xAbs*xAbs + yAbs*yAbs);
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

  return true;
}

//_________________________________________________________________________________________________
int CompareEvents(std::vector<TrackStruct>& tracks1, std::vector<TrackStruct>& tracks2, double precision, bool printAll, std::vector<TH1*>& histos)
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
      if (!areCompatible) {
        ++nDifferences;
      }
      if (printAll || !areCompatible) {
        PrintResiduals(track1.param, itTrack2->param);
      }
      FillResiduals(track1.param, itTrack2->param, histos);
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
        // compare the track parameters
        bool areCompatible = AreTrackParamCompatible(track1.param, track2.param, precision);
        if (!areCompatible) {
          ++nDifferences;
        }
        if (printAll || !areCompatible) {
          PrintResiduals(track1.param, track2.param);
        }
        FillResiduals(track1.param, track2.param, histos);
        break;
      }
    }
  }

  // then print the missing tracks
  for (const auto& track1 : tracks1) {
    if (!track1.matchFound) {
      cout << "did not find a track in file 2 matching" << endl;
      PrintTrack(track1);
      ++nDifferences;
    }
  }

  // and finally print the additional tracks
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
    histos.emplace_back(new TH1F(Form("pT%s",extension), "pT;p_{T} (GeV/c)", 300, 0., 30.));
    histos.emplace_back(new TH1F(Form("rapidity%s",extension), "rapidity;rapidity", 200, -4.5, -2.));
    histos.emplace_back(new TH1F(Form("rAbs%s",extension), "rAbs;R_{abs} (cm)", 1000, 0., 100.));
    histos.emplace_back(new TH1F(Form("dca%s",extension), "DCA;DCA (cm)", 500, 0., 500.));
    histos.emplace_back(new TH1F(Form("pDCA23%s",extension), "pDCA for #theta_{abs} < 3#circ;pDCA (GeV.cm/c)", 2500, 0., 5000.));
    histos.emplace_back(new TH1F(Form("pDCA310%s",extension), "pDCA for #theta_{abs} > 3#circ;pDCA (GeV.cm/c)", 2500, 0., 5000.));
    histos.emplace_back(new TH1F(Form("nClusters%s",extension), "number of clusters per track;n_{clusters}", 20, 0., 20.));
    histos.emplace_back(new TH1F(Form("chi2%s",extension), "normalized #chi^{2};#chi^{2} / ndf", 500, 0., 50.));
    histos.emplace_back(new TH1F(Form("mass%s",extension), "#mu^{+}#mu^{-} invariant mass;mass (GeV/c^{2})", 1600, 0., 20.));
  }
}

//_________________________________________________________________________________________________
void FillHistosAtVertex(const std::vector<TrackStruct>& tracks, std::vector<TH1*>& histos)
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

  double thetaAbs = TMath::ATan(track.rAbs/505.) * TMath::RadToDeg();
  double pUncorr = TMath::Sqrt(track.param.px*track.param.px + track.param.py*track.param.py + track.param.pz*track.param.pz);
  double pDCA = pUncorr * track.dca;
  
  histos[0]->Fill(track.pxpypzm.Pt());
  histos[1]->Fill(track.pxpypzm.Rapidity());
  histos[2]->Fill(track.rAbs);
  histos[3]->Fill(track.dca);
  if (thetaAbs < 3) {
    histos[4]->Fill(pDCA);
  } else {
    histos[5]->Fill(pDCA);
  }
  histos[6]->Fill(track.clusters.size());
  histos[7]->Fill(track.chi2 / (2. * track.clusters.size() - 5.));
}

//_________________________________________________________________________________________________
void FillHistosDimuAtVertex(const TrackStruct& track1, const TrackStruct& track2, std::vector<TH1*>& histos)
{
  /// fill dimuon histograms at vertex
  if (track1.sign * track2.sign < 0) {
    PxPyPzMVector dimu = track1.pxpypzm + track2.pxpypzm;
    histos[8]->Fill(dimu.M());
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
    cHist->cd((i / nPadsx) * nPadsx + i % nPadsx + 1);
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
    cDiff->cd((i / nPadsx) * nPadsx + i % nPadsx + 1);
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
    cRat->cd((i / nPadsx) * nPadsx + i % nPadsx + 1);
    TH1F* hRat = static_cast<TH1F*>(histos[1][i]->Clone());
    hRat->Divide(histos[0][i]);
    hRat->SetStats(false);
    hRat->SetLineColor(2);
    hRat->Draw();
  }
}

//_________________________________________________________________________________________________
void FillComparisonsAtVertex(std::vector<TrackStruct>& tracks1, std::vector<TrackStruct>& tracks2, std::vector<TH1*> histos[4])
{
  /// fill comparison histograms at vertex

  for (auto itTrack21 = tracks2.begin(); itTrack21 != tracks2.end(); ++itTrack21) {

    // fill histograms for identical, similar and additional muons
    if (itTrack21->matchIdentical) {
      FillHistosMuAtVertex(*itTrack21, histos[0]);
    } else if (itTrack21->matchFound) {
      FillHistosMuAtVertex(*itTrack21, histos[1]);
    } else {
      FillHistosMuAtVertex(*itTrack21, histos[2]);
    }

    for (auto itTrack22 = std::next(itTrack21); itTrack22 != tracks2.end(); ++itTrack22) {

      // fill histograms for identical, similar and additional dimuons
      if (itTrack21->matchIdentical && itTrack22->matchIdentical) {
        FillHistosDimuAtVertex(*itTrack21, *itTrack22, histos[0]);
      } else if (itTrack21->matchFound && itTrack22->matchFound) {
        FillHistosDimuAtVertex(*itTrack21, *itTrack22, histos[1]);
      } else {
        FillHistosDimuAtVertex(*itTrack21, *itTrack22, histos[2]);
      }
    }

    for (const auto& track1 : tracks1) {
      // fill histograms for missing dimuons
      if (!track1.matchFound) {
        FillHistosDimuAtVertex(*itTrack21, track1, histos[3]);
      }
    }
  }

  for (auto itTrack11 = tracks1.begin(); itTrack11 != tracks1.end(); ++itTrack11) {

    if (itTrack11->matchFound) {
      continue;
    }

    // fill histograms for missing muons
    FillHistosMuAtVertex(*itTrack11, histos[3]);

    for (auto itTrack12 = std::next(itTrack11); itTrack12 != tracks1.end(); ++itTrack12) {

      // fill histograms for missing dimuons
      if (!itTrack12->matchFound) {
        FillHistosDimuAtVertex(*itTrack11, *itTrack12, histos[3]);
      }
    }
  }
}

//_________________________________________________________________________________________________
void DrawComparisonsAtVertex(std::vector<TH1*> histos[4])
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
    cHist->cd((i / nPadsx) * nPadsx + i % nPadsx + 1);
    gPad->SetLogy();
    histos[0][i]->SetStats(false);
    histos[0][i]->SetLineColor(1);
    histos[0][i]->Draw();
    histos[1][i]->SetLineColor(4);
    histos[1][i]->Draw("same");
    histos[2][i]->SetLineColor(3);
    histos[2][i]->Draw("same");
    histos[3][i]->SetLineColor(2);
    histos[3][i]->Draw("same");
  }

  // add a legend
  TLegend* lHist = new TLegend(0.5, 0.5, 0.9, 0.8);
  lHist->SetFillStyle(0);
  lHist->SetBorderSize(0);
  lHist->AddEntry(histos[0][0], Form("%g tracks identical", histos[0][0]->GetEntries()), "l");
  lHist->AddEntry(histos[1][0], Form("%g tracks similar", histos[1][0]->GetEntries()), "l");
  lHist->AddEntry(histos[2][0], Form("%g tracks additional", histos[2][0]->GetEntries()), "l");
  lHist->AddEntry(histos[3][0], Form("%g tracks missing", histos[3][0]->GetEntries()), "l");
  cHist->cd(1);
  lHist->Draw("same");
}
