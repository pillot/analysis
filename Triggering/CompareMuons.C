#include <vector>
#include <tuple>
#include <string>

#include <gsl/span>

#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TMath.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TPaveStats.h>
#include <TH1F.h>
#include <TH2F.h>

#include "CCDB/BasicCCDBManager.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"

#include "DataFormatsMCH/ROFRecord.h"
#include "DataFormatsMCH/TrackMCH.h"
#include "DataFormatsMCH/Cluster.h"
#include "MCHTracking/TrackParam.h"
#include "MCHTracking/TrackExtrap.h"
#include "ReconstructionDataFormats/TrackMCHMID.h"

using namespace ROOT::Math;
using o2::InteractionRecord;
using o2::dataformats::TrackMCHMID;
using o2::mch::Cluster;
using o2::mch::ROFRecord;
using o2::mch::TrackExtrap;
using o2::mch::TrackMCH;
using o2::mch::TrackParam;

struct TrackAtVtx {
  TrackParam param{};
  double dca = 0.;
  double rAbs = 0.;
};

constexpr double pi() { return 3.14159265358979323846; }
double chi2Max = 2. * 4. * 4.;

//_________________________________________________________________________________________________
bool LoadCCDB(int runNumber);
std::tuple<TFile*, TTreeReader*> LoadData(const char* fileName, const char* treeName);
bool ExtrapToVertex(const TrackMCH& track, TrackAtVtx& trackAtVtx);
bool IsSelected(const TrackAtVtx& trackAtVtx, double pUncorr);
std::tuple<size_t, size_t> FindCompatibleMuonTracks(const TrackMCHMID& mu1, const std::vector<TrackMCHMID>& muons);
bool AreIdentical(const gsl::span<const Cluster> clusters1, const gsl::span<const Cluster> clusters2);
bool AreSimilar(const gsl::span<const Cluster> clusters1, const gsl::span<const Cluster> clusters2);
const ROFRecord* FindCompatibleROF(const InteractionRecord& ir, std::vector<ROFRecord>& rofs);
void CreateResidualsAt1stCl(std::vector<TH1*>& histos);
void FillResidualsAt1stCl(const TrackMCH& track1, const TrackMCH& track2, std::vector<TH1*>& histos);
void DrawResidualsAt1stCl(std::vector<TH1*>& histos);
void CreateHistosAtVertex(std::vector<TH1*>& histos, const char* extension);
void FillHistosAtVertex(const TrackMCH& track, const TrackAtVtx& trackAtVtx, double matchedChi2, std::vector<TH1*>& histos);
void DrawHistosAtVertex(std::vector<TH1*> histos[2]);
void DrawComparisonsAtVertex(std::vector<TH1*> histos[5]);
void CreateROFTimeHistos(std::vector<TH1*>& histos, const char* extension);
void FillROFTimeHistos(const ROFRecord& mchROF, const InteractionRecord& midTime, std::vector<TH1*>& histos);
void DrawROFTimeHistos(std::vector<TH1*> histos[2]);

//_________________________________________________________________________________________________
void CompareMuons(int runNumber,
                  string mchFileName1, string muonFileName1,
                  string mchFileName2, string muonFileName2,
                  bool applyTrackSelection = false, bool print = false)
{
  /// Compare the muon tracks stored in the 2 sets of root files and display the differences
  /// Muons are paired first by comparing their MID time then by comparing their clusters' position

  /// access CCDB and prepare track extrapolation to vertex
  LoadCCDB(runNumber);

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

  std::vector<TH1*> residualsAt1stCl{};
  CreateResidualsAt1stCl(residualsAt1stCl);
  std::vector<TH1*> histosAtVertex[2] = {{}, {}};
  CreateHistosAtVertex(histosAtVertex[0], "1");
  CreateHistosAtVertex(histosAtVertex[1], "2");
  std::vector<TH1*> comparisonsAtVertex[5] = {{}, {}, {}, {}, {}};
  CreateHistosAtVertex(comparisonsAtVertex[0], "identical");
  CreateHistosAtVertex(comparisonsAtVertex[1], "similar1");
  CreateHistosAtVertex(comparisonsAtVertex[2], "similar2");
  CreateHistosAtVertex(comparisonsAtVertex[3], "additional");
  CreateHistosAtVertex(comparisonsAtVertex[4], "missing");
  std::vector<TH1*> rofTimeHistos[2] = {{}, {}};
  CreateROFTimeHistos(rofTimeHistos[0], "1");
  CreateROFTimeHistos(rofTimeHistos[1], "2");

  int iTF = -1;
  int nMissingMu1OutOfROF2 = 0;
  int nMissingMu2OutOfROF1 = 0;
  while (mchReader1->Next() && muonReader1->Next() && mchReader2->Next() && muonReader2->Next()) {
    ++iTF;

    std::vector<TrackAtVtx> track2AtVtx(muonTracks2->size());
    std::vector<bool> isTrack2Selected(muonTracks2->size(), true);
    std::vector<bool> track2MatchFound(muonTracks2->size(), false);

    // loop over muons in file2 to extrapolate them to the vertex and select them only once
    for (size_t iMu2 = 0; iMu2 < muonTracks2->size(); ++iMu2) {
      const auto& mu2 = (*muonTracks2)[iMu2];
      const auto& track2 = (*mchTracks2)[mu2.getMCHRef().getIndex()];
      if (!ExtrapToVertex(track2, track2AtVtx[iMu2])) {
        isTrack2Selected[iMu2] = false;
        continue;
      }
      if (applyTrackSelection && !IsSelected(track2AtVtx[iMu2], track2.getP())) {
        isTrack2Selected[iMu2] = false;
        continue;
      }
      FillHistosAtVertex(track2, track2AtVtx[iMu2], mu2.getMatchChi2OverNDF(), histosAtVertex[1]);
    }

    // loop over muons in file1 and try to associate them with muons in file2
    for (size_t iMu1 = 0; iMu1 < muonTracks1->size(); ++iMu1) {
      const auto& mu1 = (*muonTracks1)[iMu1];
      const auto& track1 = (*mchTracks1)[mu1.getMCHRef().getIndex()];
      
      // extrapolate to vertex and apply track selection
      TrackAtVtx trackAtVtx1;
      if (!ExtrapToVertex(track1, trackAtVtx1)) {
        continue;
      }
      if (applyTrackSelection && !IsSelected(trackAtVtx1, track1.getP())) {
        continue;
      }

      FillHistosAtVertex(track1, trackAtVtx1, mu1.getMatchChi2OverNDF(), histosAtVertex[0]);

      // find potentially compatible muons in file2 based on MID information (IR)
      bool track1MatchFound = false;
      auto [iMu2First, iMu2Last] = FindCompatibleMuonTracks(mu1, *muonTracks2);

      if (iMu2First <= iMu2Last) {
        const gsl::span<const Cluster> clusters1(&(*mchClusters1)[track1.getFirstClusterIdx()], track1.getNClusters());
        for (size_t iMu2 = iMu2First; iMu2 <= iMu2Last; ++iMu2) {

          // skip not selected muons and muons already matched with another muon in file1
          if (!isTrack2Selected[iMu2] || track2MatchFound[iMu2]) {
            continue;
          }

          const auto& mu2 = (*muonTracks2)[iMu2];
          const auto& track2 = (*mchTracks2)[mu2.getMCHRef().getIndex()];
          const gsl::span<const Cluster> clusters2(&(*mchClusters2)[track2.getFirstClusterIdx()], track2.getNClusters());

          // check if the 2 muons are identical or similar by comparing their clusters, and record them accordingly
          if (AreIdentical(clusters1, clusters2)) {
            track1MatchFound = true;
            track2MatchFound[iMu2] = true;
            FillResidualsAt1stCl(track1, track2, residualsAt1stCl);
            FillHistosAtVertex(track1, trackAtVtx1, mu1.getMatchChi2OverNDF(), comparisonsAtVertex[0]);
            break;
          } else if (AreSimilar(clusters1, clusters2)) {
            track1MatchFound = true;
            track2MatchFound[iMu2] = true;
            FillHistosAtVertex(track1, trackAtVtx1, mu1.getMatchChi2OverNDF(), comparisonsAtVertex[1]);
            FillHistosAtVertex(track2, track2AtVtx[iMu2], mu2.getMatchChi2OverNDF(), comparisonsAtVertex[2]);
            break;
          }
        }
      }

      // if not match in found, record the muon as missing in file2
      if (!track1MatchFound) {
        FillHistosAtVertex(track1, trackAtVtx1, mu1.getMatchChi2OverNDF(), comparisonsAtVertex[4]);
        if (print) {
          printf("muon in TF %d (%s) not found in file 2.", iTF, mu1.getIR().asString().c_str());
        }

        // find the MCH ROF in file2 compatible with the muon time
        auto rof2 = FindCompatibleROF(mu1.getIR(), *mchROFs2);

        // compare the muon time with the ROF time range, if any
        if (rof2) {
          FillROFTimeHistos(*rof2, mu1.getIR(), rofTimeHistos[1]);
          if (print) {
            printf(" compatible ROF found (%s Width: %d)\n", rof2->getBCData().asString().c_str(), rof2->getBCWidth());
          }
        } else {
          ++nMissingMu1OutOfROF2;
          if (print) {
            printf(" compatible ROF not found\n");
          }
        }
      }
    }

    // loop over selected muons in file2 not found in file1 and record them as additional
    for (size_t iMu2 = 0; iMu2 < muonTracks2->size(); ++iMu2) {
      if (isTrack2Selected[iMu2] && !track2MatchFound[iMu2]) {
        const auto& mu2 = (*muonTracks2)[iMu2];
        const auto& track2 = (*mchTracks2)[mu2.getMCHRef().getIndex()];

        FillHistosAtVertex(track2, track2AtVtx[iMu2], mu2.getMatchChi2OverNDF(), comparisonsAtVertex[3]);
        if (print) {
          printf("muon in TF %d (%s) not found in file 1.", iTF, mu2.getIR().asString().c_str());
        }

        // find the MCH ROF in file1 compatible with the muon time
        auto rof1 = FindCompatibleROF(mu2.getIR(), *mchROFs1);

        // compare the muon time with the ROF time range, if any
        if (rof1) {
          FillROFTimeHistos(*rof1, mu2.getIR(), rofTimeHistos[0]);
          if (print) {
            printf(" compatible ROF found (%s Width: %d)\n", rof1->getBCData().asString().c_str(), rof1->getBCWidth());
          }
        } else {
          ++nMissingMu2OutOfROF1;
          if (print) {
            printf(" compatible ROF not found\n");
          }
        }
      }
    }
  }

  gStyle->SetOptStat(111111);
  DrawResidualsAt1stCl(residualsAt1stCl);
  DrawHistosAtVertex(histosAtVertex);
  DrawComparisonsAtVertex(comparisonsAtVertex);
  DrawROFTimeHistos(rofTimeHistos);

  printf("- number of missing muons outside of any ROF in file1 = %d\n", nMissingMu2OutOfROF1);
  printf("- number of missing muons outside of any ROF in file2 = %d\n", nMissingMu1OutOfROF2);

  fMCH1->Close();
  fMUON1->Close();
  fMCH2->Close();
  fMUON2->Close();
}

//_________________________________________________________________________________________________
bool LoadCCDB(int runNumber)
{
  /// access CCDB and prepare track extrapolation to vertex

  // load magnetic field and geometry from CCDB
  auto ccdb = o2::ccdb::BasicCCDBManager::instance();
  auto [tStart, tEnd] = ccdb.getRunDuration(runNumber);
  ccdb.setTimestamp(tEnd);
  auto grp = ccdb.get<o2::parameters::GRPMagField>("GLO/Config/GRPMagField");
  auto geom = ccdb.get<TGeoManager>("GLO/Config/GeometryAligned");

  // prepare track extrapolation to vertex (0,0,0)
  o2::base::Propagator::initFieldFromGRP(grp);
  TrackExtrap::setField();

  return true;
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
bool ExtrapToVertex(const TrackMCH& track, TrackAtVtx& trackAtVtx)
{
  /// compute the track parameters at vertex, at DCA and at the end of the absorber
  /// return false if the propagation fails

  // extrapolate to vertex
  trackAtVtx.param.setZ(track.getZ());
  trackAtVtx.param.setParameters(track.getParameters());
  if (!TrackExtrap::extrapToVertex(trackAtVtx.param, 0., 0., 0., 0., 0.)) {
    return false;
  }

  // extrapolate to DCA
  TrackParam trackParamAtDCA(track.getZ(), track.getParameters());
  if (!TrackExtrap::extrapToVertexWithoutBranson(trackParamAtDCA, 0.)) {
    return false;
  }
  double dcaX = trackParamAtDCA.getNonBendingCoor();
  double dcaY = trackParamAtDCA.getBendingCoor();
  trackAtVtx.dca = sqrt(dcaX * dcaX + dcaY * dcaY);

  // extrapolate to the end of the absorber
  TrackParam trackParamAtRAbs(track.getZ(), track.getParameters());
  if (!TrackExtrap::extrapToZ(trackParamAtRAbs, -505.)) {
    return false;
  }
  double xAbs = trackParamAtRAbs.getNonBendingCoor();
  double yAbs = trackParamAtRAbs.getBendingCoor();
  trackAtVtx.rAbs = sqrt(xAbs * xAbs + yAbs * yAbs);

  return true;
}

//_________________________________________________________________________________________________
bool IsSelected(const TrackAtVtx& trackAtVtx, double pUncorr)
{
  /// apply standard track selections + pDCA

  static const double sigmaPDCA23 = 80.;
  static const double sigmaPDCA310 = 54.;
  static const double nSigmaPDCA = 6.;
  static const double relPRes = 0.0004;
  static const double slopeRes = 0.0005;

  double thetaAbs = TMath::ATan(trackAtVtx.rAbs / 505.) * TMath::RadToDeg();
  if (thetaAbs < 2. || thetaAbs > 10.) {
    return false;
  }

  double p = trackAtVtx.param.p();
  double eta = 0.5 * log((p + trackAtVtx.param.pz()) / (p - trackAtVtx.param.pz()));
  if (eta < -4. || eta > -2.5) {
    return false;
  }

  double pDCA = pUncorr * trackAtVtx.dca;
  double sigmaPDCA = (thetaAbs < 3) ? sigmaPDCA23 : sigmaPDCA310;
  double nrp = nSigmaPDCA * relPRes * p;
  double pResEffect = sigmaPDCA / (1. - nrp / (1. + nrp));
  double slopeResEffect = 535. * slopeRes * p;
  double sigmaPDCAWithRes = TMath::Sqrt(pResEffect * pResEffect + slopeResEffect * slopeResEffect);
  if (pDCA > nSigmaPDCA * sigmaPDCAWithRes) {
    return false;
  }

  return true;
}

//_________________________________________________________________________________________________
std::tuple<size_t, size_t> FindCompatibleMuonTracks(const TrackMCHMID& mu1, const std::vector<TrackMCHMID>& muons)
{
  /// look for muons compatible in time with mu1
  /// (comparison is done with time and not with MID idx to allow matching
  ///  with another similar MID track produced by the MID reconstruction)
  size_t iFirst = 999999999;
  size_t iLast = 0;
  for (size_t iMu2 = 0; iMu2 < muons.size(); ++iMu2) {
    const auto& mu2 = muons[iMu2];
    if (mu2.getIR() < mu1.getIR()) {
      continue;
    }
    if (mu2.getIR() > mu1.getIR()) {
      break;
    }
    if (iFirst > iLast) {
      iFirst = iMu2;
    }
    iLast = iMu2;
  }
  return std::make_tuple(iFirst, iLast);
}

//_________________________________________________________________________________________________
bool AreIdentical(const gsl::span<const Cluster> clusters1, const gsl::span<const Cluster> clusters2)
{
  /// tracks are considered identical when all their clusters match within chi2Max
  if (clusters1.size() != clusters2.size()) {
    return false;
  }
  for (size_t iCl = 0; iCl != clusters1.size(); ++iCl) {
    auto& cl1 = clusters1[iCl];
    auto& cl2 = clusters2[iCl];
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

//_________________________________________________________________________________________________
bool AreSimilar(const gsl::span<const Cluster> clusters1, const gsl::span<const Cluster> clusters2)
{
  /// tracks are considered similar when:
  /// - more than 50% of clusters from one of the two tracks match within chi2Max with clusters from the other
  /// - at least 1 cluster matches before and 1 cluster matches after the dipole

  size_t nMatchClusters(0);
  bool matchCluster[10] = {false, false, false, false, false, false, false, false, false, false};

  for (const auto& cl1 : clusters1) {
    for (const auto& cl2 : clusters2) {
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
          (2 * nMatchClusters > clusters1.size() || 2 * nMatchClusters > clusters2.size()));
}

//_________________________________________________________________________________________________
const ROFRecord* FindCompatibleROF(const InteractionRecord& ir, std::vector<ROFRecord>& rofs)
{
  /// find the MCH ROF containing the IR, if any
  for (const auto& rof : rofs) {
    auto timeDiff = rof.getBCData().differenceInBC(ir);
    if (timeDiff + rof.getBCWidth() - 1 < 0) {
      continue;
    }
    if (timeDiff > 0) {
      break;
    }
    return &rof;
  }
  return nullptr;
}

//_________________________________________________________________________________________________
void CreateResidualsAt1stCl(std::vector<TH1*>& histos)
{
  /// create histograms for track parameter differences at first cluster
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
  for (auto h : histos) {
    h->SetDirectory(0);
  }
}

//_________________________________________________________________________________________________
void FillResidualsAt1stCl(const TrackMCH& track1, const TrackMCH& track2, std::vector<TH1*>& histos)
{
  /// fill histograms for track parameter differences at first cluster
  double px1 = track1.getPx();
  double py1 = track1.getPy();
  double pz1 = track1.getPz();
  double p1 = track1.getP();
  double px2 = track2.getPx();
  double py2 = track2.getPy();
  double pz2 = track2.getPz();
  double p2 = track2.getP();

  histos[0]->Fill(track2.getX() - track1.getX());
  histos[1]->Fill(track2.getY() - track1.getY());
  histos[2]->Fill(track2.getZ() - track1.getZ());
  histos[3]->Fill(px2 - px1);
  histos[4]->Fill(py2 - py1);
  histos[5]->Fill(pz2 - pz1);
  histos[6]->Fill(TMath::Abs(px1), 100. * (px2 - px1) / px1);
  histos[7]->Fill(TMath::Abs(py1), 100. * (py2 - py1) / py1);
  histos[8]->Fill(TMath::Abs(pz1), 100. * (pz2 - pz1) / pz1);
  histos[9]->Fill(p1, track2.getParameters()[1] - track1.getParameters()[1]);
  histos[10]->Fill(p1, track2.getParameters()[3] - track1.getParameters()[3]);
  histos[11]->Fill(p1, 100. * (p2 - p1) / p1);
}

//_________________________________________________________________________________________________
void DrawResidualsAt1stCl(std::vector<TH1*>& histos)
{
  /// draw histograms for track parameter differences at first cluster
  int nPadsx = (histos.size() + 2) / 3;
  TCanvas* c = new TCanvas("residuals", "residuals", 10, 10, nPadsx * 300, 900);
  c->Divide(nPadsx, 3);
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
  /// create single muon histograms at vertex
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
  }
  for (auto h : histos) {
    h->SetDirectory(0);
  }
}

//_________________________________________________________________________________________________
void FillHistosAtVertex(const TrackMCH& track, const TrackAtVtx& trackAtVtx, double matchedChi2, std::vector<TH1*>& histos)
{
  /// fill single muon histograms at vertex
  double thetaAbs = TMath::ATan(trackAtVtx.rAbs / 505.) * TMath::RadToDeg();
  double pDCA = track.getP() * trackAtVtx.dca;
  double px = trackAtVtx.param.px();
  double py = trackAtVtx.param.py();
  double pz = trackAtVtx.param.pz();
  double p = trackAtVtx.param.p();

  histos[0]->Fill(sqrt(px * px + py * py));
  histos[1]->Fill(0.5 * log((p + pz) / (p - pz)));
  histos[2]->Fill(180. + atan2(-py, -px) / pi() * 180.);
  histos[3]->Fill(trackAtVtx.rAbs);
  histos[4]->Fill(p);
  histos[5]->Fill(trackAtVtx.dca);
  if (thetaAbs < 3) {
    histos[6]->Fill(pDCA);
  } else {
    histos[7]->Fill(pDCA);
  }
  histos[8]->Fill(track.getNClusters());
  histos[9]->Fill(track.getChi2OverNDF());
  histos[10]->Fill(matchedChi2);
}

//_________________________________________________________________________________________________
void DrawHistosAtVertex(std::vector<TH1*> histos[2])
{
  /// Draw single muon histograms at vertex and differences between the 2 inputs

  // find the optimal number of pads
  int nPadsx(1), nPadsy(1);
  while ((int)histos[0].size() > nPadsx * nPadsy) {
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
    hDiff->SetDirectory(0);
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
    hRat->SetDirectory(0);
    hRat->Divide(histos[0][i]);
    hRat->SetStats(false);
    hRat->SetLineColor(2);
    hRat->Draw();
  }
}

//_________________________________________________________________________________________________
void DrawComparisonsAtVertex(std::vector<TH1*> histos[5])
{
  /// draw comparison histograms at vertex

  // find the optimal number of pads
  int nPadsx(1), nPadsy(1);
  while ((int)histos[0].size() > nPadsx * nPadsy) {
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
void CreateROFTimeHistos(std::vector<TH1*>& histos, const char* extension)
{
  /// create histograms to compare the track time with ROF time range
  if (histos.size() == 0) {
    histos.emplace_back(new TH1F(Form("rofMinVsMIDTime%s", extension), Form("ROF boundaries - MID time of missing track in file %s;#Deltat (BC)", extension), 4001, -2000.5, 2000.5));
    histos.emplace_back(new TH1F(Form("rofMaxVsMIDTime%s", extension), Form("ROF boundaries - MID time of missing track in file %s;#Deltat (BC)", extension), 4001, -2000.5, 2000.5));
    histos.emplace_back(new TH1F(Form("rofCenterVsMIDTime%s", extension), Form("ROF center - MID time of missing track in file %s;#Deltat (BC)", extension), 8001, -2000.25, 2000.25));
  }
  for (auto h : histos) {
    h->SetDirectory(0);
  }
}

//_________________________________________________________________________________________________
void FillROFTimeHistos(const ROFRecord& mchROF, const InteractionRecord& midTime, std::vector<TH1*>& histos)
{
  /// fill histograms to compare the track time with ROF time range
  auto timeDiff = mchROF.getBCData().differenceInBC(midTime);
  histos[0]->Fill(timeDiff);
  histos[1]->Fill(timeDiff + mchROF.getBCWidth() - 1);
  histos[2]->Fill(timeDiff + 0.5 * (mchROF.getBCWidth() - 1));
}

//_________________________________________________________________________________________________
void DrawROFTimeHistos(std::vector<TH1*> histos[2])
{
  /// draw histograms to compare the track time with ROF time range
  TCanvas* cHist = new TCanvas("cROFTime", "cROFTime", 10, 10, 800, 800);
  cHist->Divide(2, 2);
  cHist->cd(1);
  gPad->SetLogy();
  histos[0][0]->SetStats(false);
  histos[0][0]->SetLineColor(3);
  histos[0][0]->Draw();
  histos[0][1]->SetLineColor(6);
  histos[0][1]->Draw("same");
  cHist->cd(2);
  gPad->SetLogy();
  histos[1][0]->SetStats(false);
  histos[1][0]->SetLineColor(3);
  histos[1][0]->Draw();
  histos[1][1]->SetLineColor(6);
  histos[1][1]->Draw("same");
  cHist->cd(3);
  gPad->SetLogy();
  histos[0][2]->Draw();
  gPad->Update();
  auto st = static_cast<TPaveStats*>(histos[0][2]->FindObject("stats"));
  st->SetOptStat(10);
  st->SetY1NDC(0.85);
  gPad->Modified();
  gPad->Update();
  cHist->cd(4);
  gPad->SetLogy();
  histos[1][2]->Draw();
  gPad->Update();
  st = static_cast<TPaveStats*>(histos[1][2]->FindObject("stats"));
  st->SetOptStat(10);
  st->SetY1NDC(0.85);
  gPad->Modified();
  gPad->Update();

  // add a legend
  TLegend* lHist = new TLegend(0.65, 0.8, 0.9, 0.9);
  lHist->SetFillStyle(0);
  lHist->SetBorderSize(0);
  lHist->AddEntry(histos[0][0], "lower bound", "l");
  lHist->AddEntry(histos[0][1], "upper bound", "l");
  cHist->cd(1);
  lHist->Draw("same");
}
