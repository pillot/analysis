#include <cmath>
#include <stdexcept>
#include <vector>

#include <gsl/span>

#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>

#include "Steer/MCKinematicsReader.h"
#include "CCDB/BasicCCDBManager.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "SimulationDataFormat/MCTrack.h"
#include "SimulationDataFormat/TrackReference.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCEventHeader.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DetectorsCommonDataFormats/DetID.h"
#include "CommonUtils/NameConf.h"
#include "DataFormatsMCH/ROFRecord.h"
#include "DataFormatsMCH/TrackMCH.h"
#include "DataFormatsMCH/Cluster.h"
#include "MCHBase/TrackBlock.h"
#include "MCHTracking/TrackParam.h"
#include "MCHTracking/TrackExtrap.h"

using namespace std;
using namespace o2;
using namespace o2::base;
using namespace o2::steer;
using namespace o2::parameters;
using namespace o2::dataformats;
using o2::mch::Cluster;
using o2::mch::ROFRecord;
using o2::mch::TrackExtrap;
using o2::mch::TrackMCH;
using o2::mch::TrackParam;
using o2::mch::TrackParamStruct;

struct TrackAtVtxStruct {
  TrackParamStruct paramAtVertex{};
  double dca = 0.;
  double rAbs = 0.;
  int mchTrackIdx = 0;
};

constexpr double pi() { return 3.14159265358979323846; }
void LoadCCDB(int run);
const TrackReference& findTrackRefAt1stCl(const Cluster& cluster, gsl::span<const TrackReference> mcTrackRefs);
void makeMCClusters(gsl::span<const TrackReference> mcTrackRefs, std::vector<Cluster>& mcClusters);
bool extrapToVertex(const TrackMCH& track, double x, double y, double z, TrackAtVtxStruct& trackAtVtx);
bool IsSelected(const TrackMCH& track, const TrackAtVtxStruct& trackAtVtx);
void CreateHistos(std::vector<TH1*>& histos, const char* extension);
void CreateHistosExtra(std::vector<TH1*>& histos);
void FillHistos(const TrackMCH& track, const TrackAtVtxStruct& trackAtVtx, const MCTrack& mcTrack, int nMCClusters,
                std::vector<TH1*>& histReco, std::vector<TH1*>& histMC, std::vector<TH1*>& histRecoExtra);
void DrawHistos(std::vector<TH1*>& histReco, std::vector<TH1*>& histMC, std::vector<TH1*>& histRecoExtra);
void CreateResiduals(std::vector<TH1*>& histos);
void FillResiduals(const TrackAtVtxStruct& trackAtVtx, const MCTrack& mcTrack, std::vector<TH1*>& histos);
void FillResiduals(const TrackReference& mcTrackRef, const MCTrack& mcTrack, std::vector<TH1*>& histos);
void FillResiduals(const TrackMCH& track, const TrackAtVtxStruct& trackAtVtx, std::vector<TH1*>& histos);
void FillResiduals(const TrackMCH& track, const TrackReference& mcTrackRef, std::vector<TH1*>& histos);
void DrawResiduals(std::vector<TH1*>& histos);
void CreateResidualsAt1stCl(std::vector<TH1*>& histos);
void FillResidualsAt1stCl(const TrackMCH& track, const TrackReference& mcTrackRef, std::vector<TH1*>& histos);
void DrawResidualsAt1stCl(std::vector<TH1*>& histos);
void CreateClResiduals(std::vector<TH1*>& histos);
void FillClResiduals(const gsl::span<const Cluster> clusters, std::vector<Cluster>& mcClusters,
                     std::vector<TH1*>& histos);
void DrawClResiduals(std::vector<TH1*>& histos);

//_________________________________________________________________________________________________
void CompareTrackMC(bool applyTrackSelection = false)
{
  /// compare the reconstructed tracks to the simulated ones

  // prepare display histograms
  std::vector<TH1*> histReco{};
  CreateHistos(histReco, "");
  std::vector<TH1*> histMC{};
  CreateHistos(histMC, "MC");
  std::vector<TH1*> histRecoExtra{};
  CreateHistosExtra(histRecoExtra);
  std::vector<TH1*> residuals{};
  CreateResiduals(residuals);
  std::vector<TH1*> clResiduals{};
  CreateClResiduals(clResiduals);
  std::vector<TH1*> residualsAt1stCl{};
  CreateResidualsAt1stCl(residualsAt1stCl);

  // read the MC tracks
  MCKinematicsReader mcReader{};
  if (!mcReader.initFromDigitContext("collisioncontext.root")) {
    throw invalid_argument("initialization of MCKinematicsReader failed");
  }

  // read the reconstructed tracks
  TFile trackFile("mchtracks.root", "READ");
  if (trackFile.IsZombie()) {
    throw invalid_argument("opening of track file failed");
  }
  TTree* trackTree = static_cast<TTree*>(trackFile.Get("o2sim"));
  if (!trackTree) {
    throw invalid_argument("track tree not found");
  }
  TTreeReader dataReader;
  TTreeReaderValue<std::vector<ROFRecord>> rofs = {dataReader, "trackrofs"};
  TTreeReaderValue<std::vector<TrackMCH>> tracks = {dataReader, "tracks"};
  TTreeReaderValue<std::vector<Cluster>> clusters = {dataReader, "trackclusters"};
  TTreeReaderValue<std::vector<MCCompLabel>> labels = {dataReader, "tracklabels"};
  dataReader.SetTree(trackTree);

  // prepare track extrapolation to vertex
  // LoadCCDB(529691);
  const auto grp = GRPObject::loadFrom(NameConf::getGRPFileName());
  Propagator::initFieldFromGRP(grp);
  TrackExtrap::setField();
  GeometryManager::loadGeometry();

  int nFakes(0);
  std::vector<Cluster> mcClusters{};
  TrackAtVtxStruct trackAtVtx{};

  while (dataReader.Next()) {
    for (const auto& rof : *rofs) {
      for (int iTrack = rof.getFirstIdx(); iTrack <= rof.getLastIdx(); ++iTrack) {

        // get the reconstructed track and associated clusters and labels
        const auto& track = (*tracks)[iTrack];
        const gsl::span<const Cluster> trackClusters(&(*clusters)[track.getFirstClusterIdx()], track.getNClusters());
        const auto& l = (*labels)[iTrack];

        // skip fake tracks
        if (l.isFake()) {
          ++nFakes;
          continue;
        }

        // get the MC track, at vertex and at 1st cluster, and make the MC clusters from the associated trackRefs
        const auto* mcTrack = mcReader.getTrack(l);
        // if (!(mcTrack->GetPdgCode() == 13 && mcReader.getTrack(l.getSourceID(), l.getEventID(), mcTrack->getMotherTrackId())->GetPdgCode() == 443)) {
        //   continue;
        // }
        const auto mcTrackRefs = mcReader.getTrackRefs(l.getSourceID(), l.getEventID(), l.getTrackID());
        const auto& mcTrackRefAt1stCl = findTrackRefAt1stCl(trackClusters[0], mcTrackRefs);
        makeMCClusters(mcTrackRefs, mcClusters);

        // propagate the reconstructed track to the MC vertex
        const auto& mcHeader = mcReader.getMCEventHeader(l.getSourceID(), l.getEventID());
        if (!extrapToVertex(track, mcHeader.GetX(), mcHeader.GetY(), mcHeader.GetZ(), trackAtVtx)) {
          cout << "extrapolation to vertex failed" << endl;
          continue;
        }
        if (applyTrackSelection && !IsSelected(track, trackAtVtx)) {
          continue;
        }

        // fill histograms
        FillHistos(track, trackAtVtx, *mcTrack, mcClusters.size(), histReco, histMC, histRecoExtra);
        FillResiduals(trackAtVtx, *mcTrack, residuals);
        // FillResiduals(mcTrackRefs[0], *mcTrack, residuals);
        // FillResiduals(track, trackAtVtx, residuals);
        // FillResiduals(track, mcTrackRefAt1stCl, residuals);
        FillResidualsAt1stCl(track, mcTrackRefAt1stCl, residualsAt1stCl);
        FillClResiduals(trackClusters, mcClusters, clResiduals);
      }
    }
  }

  // draw histograms
  DrawHistos(histReco, histMC, histRecoExtra);
  DrawResiduals(residuals);
  DrawResidualsAt1stCl(residualsAt1stCl);
  DrawClResiduals(clResiduals);

  cout << "number of fake tracks = " << nFakes << endl;
}

//_________________________________________________________________________________________________
void LoadCCDB(int run)
{
  /// load magnetic field and geometry from CCDB
  auto& ccdb = o2::ccdb::BasicCCDBManager::instance();
  auto [tStart, tEnd] = ccdb.getRunDuration(run);
  ccdb.setTimestamp(tEnd);
  auto grp = ccdb.get<o2::parameters::GRPMagField>("GLO/Config/GRPMagField");
  o2::base::Propagator::initFieldFromGRP(grp);
  TrackExtrap::setField();
  auto geom = ccdb.get<TGeoManager>("GLO/Config/GeometryAligned");
}

//_________________________________________________________________________________________________
const TrackReference& findTrackRefAt1stCl(const Cluster& cluster, gsl::span<const TrackReference> mcTrackRefs)
{
  /// find the MC track reference corresponding to the reconstructed track at first cluster
  int deId = cluster.getDEId();
  for (const auto& trackRef : mcTrackRefs) {
    if (trackRef.getDetectorId() == o2::detectors::DetID::MCH && trackRef.getUserId() == deId) {
      return trackRef;
    }
  }
  cout << "MC trackRef at 1st reconstructed cluster not found --> return the 1st MC trackRef" << endl;
  return mcTrackRefs[0];
}

//_________________________________________________________________________________________________
void makeMCClusters(gsl::span<const TrackReference> mcTrackRefs, std::vector<Cluster>& mcClusters)
{
  /// produce the MC clusters by taking the average position of the trackRefs at the entry and exit of each DE
  mcClusters.clear();
  int deId(-1);
  int clusterIdx(0);
  for (const auto& trackRef : mcTrackRefs) {
    if (trackRef.getDetectorId() != o2::detectors::DetID::MCH) {
      deId = -1;
      continue;
    }
    if (trackRef.getUserId() == deId) {
      auto& cluster = mcClusters.back();
      cluster.x = (cluster.x + trackRef.X()) / 2.;
      cluster.y = (cluster.y + trackRef.Y()) / 2.;
      cluster.z = (cluster.z + trackRef.Z()) / 2.;
      deId = -1; // to create a new cluster in case the track re-enter the DE (loop)
    } else {
      deId = trackRef.getUserId();
      mcClusters.push_back({trackRef.X(), trackRef.Y(), trackRef.Z(), 0.f, 0.f,
                            Cluster::buildUniqueId(deId / 100 - 1, deId, clusterIdx++), 0u, 0u});
    }
  }
}

//_________________________________________________________________________________________________
bool extrapToVertex(const TrackMCH& track, double x, double y, double z, TrackAtVtxStruct& trackAtVtx)
{
  /// compute the track parameters at vertex, at DCA and at the end of the absorber
  /// return false if the propagation fails

  // extrapolate to vertex
  TrackParam trackParamAtVertex(track.getZ(), track.getParameters());
  if (!TrackExtrap::extrapToVertex(trackParamAtVertex, x, y, z, 0., 0.)) {
    return false;
  }
  trackAtVtx.paramAtVertex.x = trackParamAtVertex.getNonBendingCoor();
  trackAtVtx.paramAtVertex.y = trackParamAtVertex.getBendingCoor();
  trackAtVtx.paramAtVertex.z = trackParamAtVertex.getZ();
  trackAtVtx.paramAtVertex.px = trackParamAtVertex.px();
  trackAtVtx.paramAtVertex.py = trackParamAtVertex.py();
  trackAtVtx.paramAtVertex.pz = trackParamAtVertex.pz();
  trackAtVtx.paramAtVertex.sign = trackParamAtVertex.getCharge();

  // extrapolate to DCA
  TrackParam trackParamAtDCA(track.getZ(), track.getParameters());
  if (!TrackExtrap::extrapToVertexWithoutBranson(trackParamAtDCA, z)) {
    return false;
  }
  double dcaX = trackParamAtDCA.getNonBendingCoor() - x;
  double dcaY = trackParamAtDCA.getBendingCoor() - y;
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
bool IsSelected(const TrackMCH& track, const TrackAtVtxStruct& trackAtVtx)
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

  double p = sqrt(trackAtVtx.paramAtVertex.px * trackAtVtx.paramAtVertex.px +
                  trackAtVtx.paramAtVertex.py * trackAtVtx.paramAtVertex.py +
                  trackAtVtx.paramAtVertex.pz * trackAtVtx.paramAtVertex.pz);
  double eta = 0.5 * log((p + trackAtVtx.paramAtVertex.pz) / (p - trackAtVtx.paramAtVertex.pz));
  if (eta < -4. || eta > -2.5) {
    return false;
  }

  double pDCA = track.getP() * trackAtVtx.dca;
  double sigmaPDCA = (thetaAbs < 3) ? sigmaPDCA23 : sigmaPDCA310;
  double nrp = nSigmaPDCA * relPRes * p;
  double pResEffect = sigmaPDCA / (1. - nrp / (1. + nrp));
  double slopeResEffect = 535. * slopeRes * p;
  double sigmaPDCAWithRes = TMath::Sqrt(pResEffect * pResEffect + slopeResEffect * slopeResEffect);
  if (pDCA > nSigmaPDCA * sigmaPDCAWithRes) {
    return false;
  }

  // double pT = sqrt(trackAtVtx.paramAtVertex.px * trackAtVtx.paramAtVertex.px +
  //                  trackAtVtx.paramAtVertex.py * trackAtVtx.paramAtVertex.py);
  // if (pT < 1.5) {
  //   return false;
  // }

  return true;
}

//_________________________________________________________________________________________________
void CreateHistos(std::vector<TH1*>& histos, const char* extension)
{
  /// create histograms to compare simulated and reconstructed variables
  histos.emplace_back(new TH1F(Form("pT%s", extension), "pT;p_{T} (GeV/c)", 300, 0., 30.));
  histos.emplace_back(new TH1F(Form("eta%s", extension), "eta;eta", 200, -4.5, -2.));
  histos.emplace_back(new TH1F(Form("phi%s", extension), "phi;phi (deg)", 360, 0., 360.));
  histos.emplace_back(new TH1F(Form("nClusters%s", extension), "number of clusters per track;n_{clusters}", 20, 0., 20.));
}

//_________________________________________________________________________________________________
void CreateHistosExtra(std::vector<TH1*>& histos)
{
  /// create extra histograms to display reconstructed variables
  histos.emplace_back(new TH1F("dca", "DCA;DCA (cm)", 2500, 0., 500.));
  histos.emplace_back(new TH1F("pDCA23", "pDCA for #theta_{abs} < 3#circ;pDCA (GeV.cm/c)", 2500, 0., 5000.));
  histos.emplace_back(new TH1F("pDCA310", "pDCA for #theta_{abs} > 3#circ;pDCA (GeV.cm/c)", 2500, 0., 5000.));
  histos.emplace_back(new TH1F("rAbs", "rAbs;R_{abs} (cm)", 1000, 0., 100.));
  histos.emplace_back(new TH1F("chi2", "normalized #chi^{2};#chi^{2} / ndf", 500, 0., 50.));
}

//_________________________________________________________________________________________________
void FillHistos(const TrackMCH& track, const TrackAtVtxStruct& trackAtVtx, const MCTrack& mcTrack, int nMCClusters,
                std::vector<TH1*>& histReco, std::vector<TH1*>& histMC, std::vector<TH1*>& histRecoExtra)
{
  /// fill histograms of simulated and reconstructed variables

  double thetaAbs = atan(trackAtVtx.rAbs / 505.) * 180. / pi();
  double pDCA = track.getP() * trackAtVtx.dca;
  double pT = sqrt(trackAtVtx.paramAtVertex.px * trackAtVtx.paramAtVertex.px +
                   trackAtVtx.paramAtVertex.py * trackAtVtx.paramAtVertex.py);
  double p = sqrt(trackAtVtx.paramAtVertex.px * trackAtVtx.paramAtVertex.px +
                  trackAtVtx.paramAtVertex.py * trackAtVtx.paramAtVertex.py +
                  trackAtVtx.paramAtVertex.pz * trackAtVtx.paramAtVertex.pz);
  double eta = 0.5 * log((p + trackAtVtx.paramAtVertex.pz) / (p - trackAtVtx.paramAtVertex.pz));
  double phi = 180. + atan2(-trackAtVtx.paramAtVertex.py, -trackAtVtx.paramAtVertex.px) / pi() * 180.;

  histReco[0]->Fill(pT);
  histReco[1]->Fill(eta);
  histReco[2]->Fill(phi);
  histReco[3]->Fill(track.getNClusters());

  histMC[0]->Fill(mcTrack.GetPt());
  histMC[1]->Fill(mcTrack.GetEta());
  histMC[2]->Fill(mcTrack.GetPhi() / pi() * 180.);
  histMC[3]->Fill(nMCClusters);

  histRecoExtra[0]->Fill(trackAtVtx.dca);
  if (thetaAbs < 3) {
    histRecoExtra[1]->Fill(pDCA);
  } else {
    histRecoExtra[2]->Fill(pDCA);
  }
  histRecoExtra[3]->Fill(trackAtVtx.rAbs);
  histRecoExtra[4]->Fill(track.getChi2OverNDF());
}

//_________________________________________________________________________________________________
void DrawHistos(std::vector<TH1*>& histReco, std::vector<TH1*>& histMC, std::vector<TH1*>& histRecoExtra)
{
  /// draw histograms of simulated and reconstructed variables

  TCanvas* cHist = new TCanvas("histos", "histos", 10, 10, 1200, 600);
  cHist->Divide(4, 2);
  for (int i = 0; i < 4; ++i) {
    cHist->cd(i + 1);
    gPad->SetLogy();
    histMC[i]->SetLineColor(2);
    histMC[i]->Draw();
    histReco[i]->Draw("sames");
    cHist->cd(i + 5);
    auto* hRatio = static_cast<TH1*>(histReco[i]->Clone(Form("%s_ratio", histReco[i]->GetName())));
    hRatio->SetTitle("ratio");
    hRatio->Divide(histMC[i]);
    hRatio->SetStats(false);
    hRatio->Draw();
  }

  TCanvas* cHistExtra = new TCanvas("histos2", "histos2", 20, 20, 900, 600);
  cHistExtra->Divide(3, 2);
  for (int i = 0; i < 5; ++i) {
    cHistExtra->cd(i + 1);
    gPad->SetLogy();
    histRecoExtra[i]->Draw();
  }
}

//_________________________________________________________________________________________________
void CreateResiduals(std::vector<TH1*>& histos)
{
  /// create histograms holding residuals between simulated and reconstructed variables
  histos.emplace_back(new TH2F("dslopexvsp", "dslopexvsp;p (GeV/c);dslopex", 300, 0., 300., 1001, -0.05005, 0.05005));
  histos.emplace_back(new TH2F("dslopeyvsp", "dslopeyvsp;p (GeV/c);dslopey", 300, 0., 300., 1001, -0.05005, 0.05005));
  histos.emplace_back(new TH2F("dpvsp", "dpvsp;p (GeV/c);dp (GeV/c)", 300, 0., 300., 501, -50.1, 50.1));
  histos.emplace_back(new TH1F("dpT", "dpT;dp_{T} (GeV/c)", 501, -5.01, 5.01));
  histos.emplace_back(new TH1F("deta", "deta;deta", 1001, -0.5005, 0.5005));
  histos.emplace_back(new TH1F("dphi", "dphi;dphi (rad)", 1001, -0.5005, 0.5005));
}

//_________________________________________________________________________________________________
void FillResiduals(const TrackAtVtxStruct& trackAtVtx, const MCTrack& mcTrack, std::vector<TH1*>& histos)
{
  /// fill residuals between simulated and reconstructed variables at vertex
  double pT = sqrt(trackAtVtx.paramAtVertex.px * trackAtVtx.paramAtVertex.px +
                   trackAtVtx.paramAtVertex.py * trackAtVtx.paramAtVertex.py);
  double p = sqrt(trackAtVtx.paramAtVertex.px * trackAtVtx.paramAtVertex.px +
                  trackAtVtx.paramAtVertex.py * trackAtVtx.paramAtVertex.py +
                  trackAtVtx.paramAtVertex.pz * trackAtVtx.paramAtVertex.pz);
  double eta = 0.5 * log((p + trackAtVtx.paramAtVertex.pz) / (p - trackAtVtx.paramAtVertex.pz));
  double phi = pi() + atan2(-trackAtVtx.paramAtVertex.py, -trackAtVtx.paramAtVertex.px);
  double pMC = mcTrack.GetP();
  histos[0]->Fill(pMC, trackAtVtx.paramAtVertex.px / trackAtVtx.paramAtVertex.pz - mcTrack.Px() / mcTrack.Pz());
  histos[1]->Fill(pMC, trackAtVtx.paramAtVertex.py / trackAtVtx.paramAtVertex.pz - mcTrack.Py() / mcTrack.Pz());
  histos[2]->Fill(pMC, p - pMC);
  histos[3]->Fill(pT - mcTrack.GetPt());
  histos[4]->Fill(eta - mcTrack.GetEta());
  histos[5]->Fill(phi - mcTrack.GetPhi());
}

//_________________________________________________________________________________________________
void FillResiduals(const TrackReference& mcTrackRef, const MCTrack& mcTrack, std::vector<TH1*>& histos)
{
  /// fill residuals between simulated variables at vertex and at first cluster
  double p = mcTrackRef.P();
  double eta = 0.5 * log((p + mcTrackRef.Pz()) / (p - mcTrackRef.Pz()));
  double phi = pi() + atan2(-mcTrackRef.Py(), -mcTrackRef.Px());
  double pMC = mcTrack.GetP();
  histos[0]->Fill(pMC, mcTrackRef.Px() / mcTrackRef.Pz() - mcTrack.Px() / mcTrack.Pz());
  histos[1]->Fill(pMC, mcTrackRef.Py() / mcTrackRef.Pz() - mcTrack.Py() / mcTrack.Pz());
  histos[2]->Fill(pMC, p - pMC);
  histos[3]->Fill(mcTrackRef.Pt() - mcTrack.GetPt());
  histos[4]->Fill(eta - mcTrack.GetEta());
  histos[5]->Fill(phi - mcTrack.GetPhi());
}

//_________________________________________________________________________________________________
void FillResiduals(const TrackMCH& track, const TrackAtVtxStruct& trackAtVtx, std::vector<TH1*>& histos)
{
  /// fill residuals between reconstructed variables at vertex and at first cluster
  double px1 = track.getPx();
  double py1 = track.getPy();
  double pz1 = track.getPz();
  double pT1 = sqrt(px1 * px1 + py1 * py1);
  double p1 = track.getP();
  double eta1 = 0.5 * log((p1 + pz1) / (p1 - pz1));
  double phi1 = pi() + atan2(-py1, -px1);
  double pT = sqrt(trackAtVtx.paramAtVertex.px * trackAtVtx.paramAtVertex.px +
                   trackAtVtx.paramAtVertex.py * trackAtVtx.paramAtVertex.py);
  double p = sqrt(trackAtVtx.paramAtVertex.px * trackAtVtx.paramAtVertex.px +
                  trackAtVtx.paramAtVertex.py * trackAtVtx.paramAtVertex.py +
                  trackAtVtx.paramAtVertex.pz * trackAtVtx.paramAtVertex.pz);
  double eta = 0.5 * log((p + trackAtVtx.paramAtVertex.pz) / (p - trackAtVtx.paramAtVertex.pz));
  double phi = pi() + atan2(-trackAtVtx.paramAtVertex.py, -trackAtVtx.paramAtVertex.px);
  histos[0]->Fill(p, px1 / pz1 - trackAtVtx.paramAtVertex.px / trackAtVtx.paramAtVertex.pz);
  histos[1]->Fill(p, py1 / pz1 - trackAtVtx.paramAtVertex.py / trackAtVtx.paramAtVertex.pz);
  histos[2]->Fill(p, p1 - p);
  histos[3]->Fill(pT1 - pT);
  histos[4]->Fill(eta1 - eta);
  histos[5]->Fill(phi1 - phi);
}

//_________________________________________________________________________________________________
void FillResiduals(const TrackMCH& track, const TrackReference& mcTrackRef, std::vector<TH1*>& histos)
{
  /// fill residuals between simulated and reconstructed variables at first cluster
  double px1 = track.getPx();
  double py1 = track.getPy();
  double pz1 = track.getPz();
  double pT1 = sqrt(px1 * px1 + py1 * py1);
  double p1 = track.getP();
  double eta1 = 0.5 * log((p1 + pz1) / (p1 - pz1));
  double phi1 = pi() + atan2(-py1, -px1);
  double p = mcTrackRef.P();
  double eta = 0.5 * log((p + mcTrackRef.Pz()) / (p - mcTrackRef.Pz()));
  double phi = pi() + atan2(-mcTrackRef.Py(), -mcTrackRef.Px());
  histos[0]->Fill(p, px1 / pz1 - mcTrackRef.Px() / mcTrackRef.Pz());
  histos[1]->Fill(p, py1 / pz1 - mcTrackRef.Py() / mcTrackRef.Pz());
  histos[2]->Fill(p, p1 - p);
  histos[3]->Fill(pT1 - mcTrackRef.Pt());
  histos[4]->Fill(eta1 - eta);
  histos[5]->Fill(phi1 - phi);
}

//_________________________________________________________________________________________________
void DrawResiduals(std::vector<TH1*>& histos)
{
  /// draw residuals between simulated and reconstructed variables
  TCanvas* c = new TCanvas("residuals", "residuals", 30, 30, 900, 600);
  c->Divide(3, 2);
  for (int i = 0; i < 3; ++i) {
    c->cd(i + 1);
    gPad->SetLogz();
    histos[i]->Draw("colz");
  }
  for (int i = 3; i < 6; ++i) {
    c->cd(i + 1);
    gPad->SetLogy();
    histos[i]->Draw();
  }
}

//_________________________________________________________________________________________________
void CreateResidualsAt1stCl(std::vector<TH1*>& histos)
{
  /// create histograms holding residuals between simulated and reconstructed track parameters at first cluster
  histos.emplace_back(new TH1F("dx", "dx;dx (cm)", 20001, -1.00005, 1.00005));
  histos.emplace_back(new TH1F("dy", "dy;dy (cm)", 20001, -1.00005, 1.00005));
  histos.emplace_back(new TH1F("dz", "dz;dz (cm)", 20001, -1.00005, 1.00005));
  histos.emplace_back(new TH1F("dpx", "dpx;dpx (GeV/c)", 20001, -1.00005, 1.00005));
  histos.emplace_back(new TH1F("dpy", "dpy;dpy (GeV/c)", 20001, -1.00005, 1.00005));
  histos.emplace_back(new TH1F("dpz", "dpz;dpz (GeV/c)", 20001, -1.00005, 1.00005));
  histos.emplace_back(new TH2F("dpxvspx", "dpxvspx;px (GeV/c);dpx/px (\%)", 2000, 0., 20., 2001, -10.005, 10.005));
  histos.emplace_back(new TH2F("dpyvspy", "dpyvspy;py (GeV/c);dpy/py (\%)", 2000, 0., 20., 2001, -10.005, 10.005));
  histos.emplace_back(new TH2F("dpzvspz", "dpzvspz;pz (GeV/c);dpz/pz (\%)", 2000, 0., 200., 2001, -10.005, 10.005));
  histos.emplace_back(new TH2F("dslopexvsp1", "dslopexvsp;p (GeV/c);dslopex", 2000, 0., 200., 2001, -0.0010005, 0.0010005));
  histos.emplace_back(new TH2F("dslopeyvsp1", "dslopeyvsp;p (GeV/c);dslopey", 2000, 0., 200., 2001, -0.0010005, 0.0010005));
  histos.emplace_back(new TH2F("dpvsp1", "dpvsp;p (GeV/c);dp/p (\%)", 2000, 0., 200., 2001, -10.005, 10.005));
}

//_________________________________________________________________________________________________
void FillResidualsAt1stCl(const TrackMCH& track, const TrackReference& mcTrackRef, std::vector<TH1*>& histos)
{
  /// fill residuals between simulated and reconstructed track parameters at first cluster
  histos[0]->Fill(track.getX() - mcTrackRef.X());
  histos[1]->Fill(track.getY() - mcTrackRef.Y());
  histos[2]->Fill(track.getZ() - mcTrackRef.Z());
  histos[3]->Fill(track.getPx() - mcTrackRef.Px());
  histos[4]->Fill(track.getPy() - mcTrackRef.Py());
  histos[5]->Fill(track.getPz() - mcTrackRef.Pz());
  histos[6]->Fill(TMath::Abs(mcTrackRef.Px()), 100. * (track.getPx() - mcTrackRef.Px()) / mcTrackRef.Px());
  histos[7]->Fill(TMath::Abs(mcTrackRef.Py()), 100. * (track.getPy() - mcTrackRef.Py()) / mcTrackRef.Py());
  histos[8]->Fill(TMath::Abs(mcTrackRef.Pz()), 100. * (track.getPz() - mcTrackRef.Pz()) / mcTrackRef.Pz());
  histos[9]->Fill(mcTrackRef.P(), track.getParameters()[1] - mcTrackRef.Px() / mcTrackRef.Pz());
  histos[10]->Fill(mcTrackRef.P(), track.getParameters()[3] - mcTrackRef.Py() / mcTrackRef.Pz());
  histos[11]->Fill(mcTrackRef.P(), 100. * (track.getP() - mcTrackRef.P()) / mcTrackRef.P());
}

//_________________________________________________________________________________________________
void DrawResidualsAt1stCl(std::vector<TH1*>& histos)
{
  /// draw residuals between simulated and reconstructed track parameters at first cluster
  int nPadsx = (histos.size() + 2) / 3;
  TCanvas* c = new TCanvas("residualsAt1stCl", "residualsAt1stCl", 10, 10, nPadsx * 300, 900);
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
void CreateClResiduals(std::vector<TH1*>& histos)
{
  /// create histograms holding residuals between simulated and reconstructed clusters
  const double range = 1.0005;
  for (int iSt = 1; iSt <= 5; ++iSt) {
    histos.emplace_back(new TH1F(Form("resX_St%d", iSt),
                                 Form("#DeltaX Station %d;#DeltaX (cm)", iSt), 2001, -range, range));
    histos.emplace_back(new TH1F(Form("resY_St%d", iSt),
                                 Form("#DeltaY Station %d;#DeltaY (cm)", iSt), 2001, -range, range));
  }
  histos.emplace_back(new TH1F("resX", "#DeltaX;#DeltaX (cm)", 2001, -range, range));
  histos.emplace_back(new TH1F("resY", "#DeltaY;#DeltaY (cm)", 2001, -range, range));
}

//_________________________________________________________________________________________________
void FillClResiduals(const gsl::span<const Cluster> clusters, std::vector<Cluster>& mcClusters,
                     std::vector<TH1*>& histos)
{
  /// fill residuals between simulated and reconstructed clusters
  for (const auto& cluster : clusters) {
    for (const auto& mcCluster : mcClusters) {
      if (cluster.getDEId() == mcCluster.getDEId()) {
        histos[cluster.getChamberId() / 2 * 2]->Fill(cluster.x - mcCluster.x);
        histos[cluster.getChamberId() / 2 * 2 + 1]->Fill(cluster.y - mcCluster.y);
        histos[10]->Fill(cluster.x - mcCluster.x);
        histos[11]->Fill(cluster.y - mcCluster.y);
      }
    }
  }
}

//_________________________________________________________________________________________________
void DrawClResiduals(std::vector<TH1*>& histos)
{
  /// draw residuals between simulated and reconstructed clusters
  int nPadsx = (histos.size() + 1) / 2;
  TCanvas* c = new TCanvas("clResiduals", "cluster residuals", 40, 40, nPadsx * 300, 600);
  c->Divide(nPadsx, 2);
  int i(0);
  for (const auto& h : histos) {
    c->cd((i % 2) * nPadsx + i / 2 + 1);
    gPad->SetLogy();
    h->Draw();
    ++i;
  }
}
