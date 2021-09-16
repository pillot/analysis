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
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "SimulationDataFormat/MCTrack.h"
#include "SimulationDataFormat/TrackReference.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include "SimulationDataFormat/MCEventHeader.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsCommonDataFormats/DetID.h"
#include "DetectorsCommonDataFormats/NameConf.h"
#include "DataFormatsMCH/ROFRecord.h"
#include "DataFormatsMCH/TrackMCH.h"
#include "MCHBase/ClusterBlock.h"
#include "MCHBase/TrackBlock.h"
#include "MCHTracking/TrackParam.h"
#include "MCHTracking/TrackExtrap.h"

using namespace std;
using namespace o2;
using namespace o2::base;
using namespace o2::steer;
using namespace o2::parameters;
using namespace o2::dataformats;
using namespace o2::mch;

struct TrackAtVtxStruct {
  TrackParamStruct paramAtVertex{};
  double dca = 0.;
  double rAbs = 0.;
  int mchTrackIdx = 0;
};

constexpr double pi() { return 3.14159265358979323846; }
void makeMCClusters(gsl::span<const TrackReference> mcTrackRefs, std::vector<ClusterStruct>& mcClusters);
bool extrapToVertex(const TrackMCH& track, double x, double y, double z, TrackAtVtxStruct& trackAtVtx);
void CreateHistos(std::vector<TH1*>& histos, const char* extension);
void CreateHistosExtra(std::vector<TH1*>& histos);
void FillHistos(const TrackMCH& track, const TrackAtVtxStruct& trackAtVtx, const MCTrack& mcTrack, int nMCClusters,
                std::vector<TH1*>& histReco, std::vector<TH1*>& histMC, std::vector<TH1*>& histRecoExtra);
void DrawHistos(std::vector<TH1*>& histReco, std::vector<TH1*>& histMC, std::vector<TH1*>& histRecoExtra);
void CreateResiduals(std::vector<TH1*>& histos);
void FillResiduals(const TrackAtVtxStruct& trackAtVtx, const MCTrack& mcTrack, std::vector<TH1*>& histos);
void DrawResiduals(std::vector<TH1*>& histos);
void CreateClResiduals(std::vector<TH1*>& histos);
void FillClResiduals(const gsl::span<const ClusterStruct> clusters, std::vector<ClusterStruct>& mcClusters,
                     std::vector<TH1*>& histos);
void DrawClResiduals(std::vector<TH1*>& histos);

//_________________________________________________________________________________________________
void CompareTrackMC()
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
  TTreeReaderValue<std::vector<ClusterStruct>> clusters = {dataReader, "trackclusters"};
  TTreeReaderValue<MCTruthContainer<MCCompLabel>> labels = {dataReader, "tracklabels"};
  dataReader.SetTree(trackTree);

  // prepare track extrapolation to vertex
  const auto grp = GRPObject::loadFrom(NameConf::getGRPFileName());
  Propagator::initFieldFromGRP(grp);
  TrackExtrap::setField();
  GeometryManager::loadGeometry();

  int nFakes(0);
  int nMultipleLabels(0);
  std::vector<ClusterStruct> mcClusters{};
  TrackAtVtxStruct trackAtVtx{};

  while (dataReader.Next()) {
    for (const auto& rof : *rofs) {
      for (int iTrack = rof.getFirstIdx(); iTrack <= rof.getLastIdx(); ++iTrack) {

        // get the reconstructed track and associated clusters and labels
        const auto& track = (*tracks)[iTrack];
        const gsl::span<const ClusterStruct> trackClusters(&(*clusters)[track.getFirstClusterIdx()], track.getNClusters());
        const auto trackLabels = labels->getLabels(iTrack);

        // skip fake tracks
        if (trackLabels.size() < 1) {
          ++nFakes;
          continue;
        }

        // check if the reconstructed track matches several MC tracks
        if (trackLabels.size() > 1) {
          ++nMultipleLabels;
        }

        // get the MC track and make the MC clusters from the associated trackRefs
        const auto& l = trackLabels[0];
        const auto* mcTrack = mcReader.getTrack(l);
        const auto mcTrackRefs = mcReader.getTrackRefs(l.getSourceID(), l.getEventID(), l.getTrackID());
        makeMCClusters(mcTrackRefs, mcClusters);

        // propagate the reconstructed track to the MC vertex
        const auto& mcHeader = mcReader.getMCEventHeader(l.getSourceID(), l.getEventID());
        if (!extrapToVertex(track, mcHeader.GetX(), mcHeader.GetY(), mcHeader.GetZ(), trackAtVtx)) {
          cout << "extrapolation to vertex failed" << endl;
          continue;
        }

        // fill histograms
        FillHistos(track, trackAtVtx, *mcTrack, mcClusters.size(), histReco, histMC, histRecoExtra);
        FillResiduals(trackAtVtx, *mcTrack, residuals);
        FillClResiduals(trackClusters, mcClusters, clResiduals);
      }
    }
  }

  // draw histograms
  DrawHistos(histReco, histMC, histRecoExtra);
  DrawResiduals(residuals);
  DrawClResiduals(clResiduals);

  cout << "number of fake tracks = " << nFakes << endl;
  cout << "number of tracks matched with several MC tracks = " << nMultipleLabels << endl;
}

//_________________________________________________________________________________________________
void makeMCClusters(gsl::span<const TrackReference> mcTrackRefs, std::vector<ClusterStruct>& mcClusters)
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
                            ClusterStruct::buildUniqueId(deId / 100 - 1, deId, clusterIdx++), 0u, 0u});
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
void CreateHistos(std::vector<TH1*>& histos, const char* extension)
{
  /// create histograms to compare simulated and reconstructed variables
  histos.emplace_back(new TH1F(Form("pT%s", extension), "pT;p_{T} (GeV/c)", 300, 0., 30.));
  histos.emplace_back(new TH1F(Form("eta%s", extension), "eta;eta", 200, -4.5, -2.));
  histos.emplace_back(new TH1F(Form("phi%s", extension), "phi;phi", 360, 0., 360.));
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
  double phi = 180. + atan2(-trackAtVtx.paramAtVertex.px, -trackAtVtx.paramAtVertex.py) / pi() * 180.;

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
}

//_________________________________________________________________________________________________
void FillResiduals(const TrackAtVtxStruct& trackAtVtx, const MCTrack& mcTrack, std::vector<TH1*>& histos)
{
  /// fill residuals between simulated and reconstructed variables
  double p = sqrt(trackAtVtx.paramAtVertex.px * trackAtVtx.paramAtVertex.px +
                  trackAtVtx.paramAtVertex.py * trackAtVtx.paramAtVertex.py +
                  trackAtVtx.paramAtVertex.pz * trackAtVtx.paramAtVertex.pz);
  double pMC = mcTrack.GetP();
  histos[0]->Fill(pMC, trackAtVtx.paramAtVertex.px / trackAtVtx.paramAtVertex.pz - mcTrack.Px() / mcTrack.Pz());
  histos[1]->Fill(pMC, trackAtVtx.paramAtVertex.py / trackAtVtx.paramAtVertex.pz - mcTrack.Py() / mcTrack.Pz());
  histos[2]->Fill(pMC, p - pMC);
}

//_________________________________________________________________________________________________
void DrawResiduals(std::vector<TH1*>& histos)
{
  /// draw residuals between simulated and reconstructed variables
  TCanvas* c = new TCanvas("residuals", "residuals", 30, 30, 900, 300);
  c->Divide(3, 1);
  for (int i = 0; i < 3; ++i) {
    c->cd(i + 1);
    gPad->SetLogz();
    histos[i]->Draw("colz");
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
void FillClResiduals(const gsl::span<const ClusterStruct> clusters, std::vector<ClusterStruct>& mcClusters,
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
