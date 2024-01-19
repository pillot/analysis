#include <cmath>
#include <memory>
#include <tuple>
#include <utility>
#include <vector>

#include <gsl/span>

#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGeoManager.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLegend.h>
#include <TMath.h>
#include <TStyle.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DetectorsBase/Propagator.h"
#include "MathUtils/Cartesian.h"

#include "DataFormatsMCH/Cluster.h"
#include "DataFormatsMCH/Digit.h"
#include "DataFormatsMCH/TrackMCH.h"
#include "MCHGeometryTransformer/Transformations.h"
#include "MCHMappingInterface/Segmentation.h"
#include "MCHTracking/Track.h"
#include "MCHTracking/TrackExtrap.h"
#include "MCHTracking/TrackFitter.h"
#include "MCHTracking/TrackParam.h"
#include "ReconstructionDataFormats/TrackMCHMID.h"

using o2::dataformats::TrackMCHMID;
using o2::mch::Cluster;
using o2::mch::Digit;
using o2::mch::Track;
using o2::mch::TrackExtrap;
using o2::mch::TrackFitter;
using o2::mch::TrackMCH;
using o2::mch::TrackParam;

struct TrackStruct {
  Track track{};
  TrackParam paramAtVtx{};
  double dca = 0.;
  double rAbs = 0.;
};

void InitFromCCDB(int run);
std::tuple<TFile*, TTreeReader*> LoadData(const char* fileName, const char* treeName);
bool ClustersToTrack(const gsl::span<const Cluster> trackClusters, TrackStruct& track);
bool ExtrapToVertex(TrackStruct& track);
bool IsSelected(const TrackStruct& track);
bool IsMonoCathode(const gsl::span<const Digit> digits);
std::pair<double, double> GetCharge(const gsl::span<const Digit> digits);
std::unique_ptr<Cluster> COG(const gsl::span<const Digit> digits);
std::unique_ptr<Cluster> COG2(const gsl::span<const Digit> digits);
std::unique_ptr<Cluster> MakeCluster(int de, float x, float y, float ex, float ey);
void CreateResiduals(std::vector<TH1*>& histos, const char* extension, double range);
void FillResiduals(const TrackParam& param, const Cluster& cluster, std::vector<TH1*>& histos);
void FillResiduals(const Cluster& cluster1, const Cluster& cluster2, std::vector<TH1*>& histos);
void DrawResiduals(std::vector<TH1*>& histos);
void DrawResiduals(std::vector<TH1*>& oldHistos, std::vector<TH1*>& newHistos);
void DrawRatios(std::vector<TH1*>& oldHistos, std::vector<TH1*>& newHistos);
std::pair<double, double> GetSigma(TH1* h, int color);
double CrystalBallSymmetric(double* xx, double* par);

//_________________________________________________________________________________________________
void ClusterFitter(int run, const char* clusterFile, const char* trackFile, const char* muonFile,
                   bool applyTrackSelection = false)
{
  /// select the isolated clusters attached to a (selected) muon track
  /// refit the corresponding preclusters with a Mathieson function
  /// check the effect on the cluster-track residuals (~ resolution)
  /// MCH mapping need to be loaded before: gSystem->Load("libO2MCHMappingImpl4")

  /// load CCDB objects and prepare track fitting, geometry transformation, ...
  InitFromCCDB(run);

  // load clusters and tracks
  auto [fClusters, clusterReader] = LoadData(clusterFile, "o2sim");
  TTreeReaderValue<std::vector<Cluster>> allClusters(*clusterReader, "clusters");
  TTreeReaderValue<std::vector<Digit>> allClusterDigits(*clusterReader, "clusterdigits");
  auto [fTracks, trackReader] = LoadData(trackFile, "o2sim");
  TTreeReaderValue<std::vector<TrackMCH>> allTracks(*trackReader, "tracks");
  TTreeReaderValue<std::vector<Cluster>> allTrackClusters(*trackReader, "trackclusters");
  auto [fMuons, muonReader] = LoadData(muonFile, "o2sim");
  TTreeReaderValue<std::vector<TrackMCHMID>> allMuons(*muonReader, "tracks");

  // make sure the track clusters point to digits in the cluster file
  if (trackReader->GetTree()->FindBranch("trackdigits")) {
    LOG(error) << " track clusters point to digits in the track file. Re-run the MCH reco without the option --digits";
    exit(-1);
  }

  int nTF = clusterReader->GetEntries(false);
  if (trackReader->GetEntries(false) != nTF ||
      muonReader->GetEntries(false) != nTF) {
    LOG(error) << " not all files contain the same number of TF";
    exit(-1);
  }

  TH2* hChargeAsymm = new TH2F("hChargeAsymm", "cluster charge asymmetry;charge (ADC); asymmetry",
                               200, 0., 20000., 201, -1.005, 1.005);
  hChargeAsymm->SetDirectory(0);
  std::vector<TH1*> oldResiduals{};
  CreateResiduals(oldResiduals, "old", 2.);
  std::vector<TH1*> newResiduals{};
  CreateResiduals(newResiduals, "new", 2.);
  std::vector<TH1*> clclResiduals{};
  CreateResiduals(clclResiduals, "clcl", 0.2);

  int iTF(0);
  while (clusterReader->Next() && trackReader->Next() && muonReader->Next()) {
    std::cout << "\rprocessing TF " << ++iTF << " / " << nTF << "..." << std::flush;

    // count the number of clusters associated to each precluster (identified by the first digit index)
    std::unordered_map<uint32_t, int> nClusters{};
    for (const auto& cluster : *allClusters) {
      nClusters[cluster.firstDigit]++;
    }

    for (const auto& muon : *allMuons) {
      const auto& mchTrack = (*allTracks)[muon.getMCHRef().getIndex()];
      const gsl::span<const Cluster> trackClusters(&(*allTrackClusters)[mchTrack.getFirstClusterIdx()], mchTrack.getNClusters());

      // create an internal track and compute its parameters at each attached cluster
      TrackStruct track{};
      if (!ClustersToTrack(trackClusters, track)) {
        continue;
      }

      // propagate the track to the vertex (0,0,0) and apply selections if requested
      if (applyTrackSelection && !(ExtrapToVertex(track) && IsSelected(track))) {
        continue;
      }

      for (const auto& param : track.track) {
        const auto& cluster = *(param.getClusterPtr());
        const gsl::span<const Digit> digits(&(*allClusterDigits)[cluster.firstDigit], cluster.nDigits);

        // those 2 DE have lower HV for the run 529691
        // if (digits[0].getDetID() == 202 || digits[0].getDetID() == 300) {
        //   continue;
        // }

        // select clusters of interest (isolated, bi-cathode, #digits > xxx, ...)
        if (cluster.nDigits < 4 || nClusters[cluster.firstDigit] != 1 || IsMonoCathode(digits)) {
          continue;
        }

        // check the charge asymmetry between the 2 cathodes
        const auto [chargeNB, chargeB] = GetCharge(digits);
        double chargeAsymm = (chargeNB - chargeB) / (chargeNB + chargeB);
        hChargeAsymm->Fill(0.5 * (chargeNB + chargeB), chargeAsymm);
        if (std::abs(chargeAsymm) > 0.5) {
          continue;
        }

        // fill cluster-track residuals
        FillResiduals(param, cluster, oldResiduals);

        // refit the cluster
        auto newCluster = (cluster.getChamberId() < 4) ? COG2(digits) : COG(digits);
        // auto newCluster = FitMathieson(digits);
        newCluster->uid = cluster.uid;
        newCluster->firstDigit = cluster.firstDigit;
        newCluster->nDigits = cluster.nDigits;

        // fill the new cluster-track and cluster-cluster residuals
        FillResiduals(param, *newCluster, newResiduals);
        FillResiduals(cluster, *newCluster, clclResiduals);
      }
    }
  }
  cout << "\r\033[Kprocessing completed" << endl;

  gStyle->SetOptStat(1);

  TCanvas* cChargeAsymm = new TCanvas("cChargeAsymm", "cluster charge asymmetry");
  hChargeAsymm->Draw("colz");
  gPad->SetLogz();

  DrawResiduals(clclResiduals);
  DrawResiduals(oldResiduals, newResiduals);
  DrawRatios(oldResiduals, newResiduals);

  fClusters->Close();
  fTracks->Close();
  fMuons->Close();
}

//_________________________________________________________________________________________________
void InitFromCCDB(int run)
{
  /// load necessary objects from CCDB and prepare track fitting, geometry transformation, ...

  // load magnetic field and geometry from CCDB
  auto& ccdb = o2::ccdb::BasicCCDBManager::instance();
  auto [tStart, tEnd] = ccdb.getRunDuration(run);
  ccdb.setTimestamp(tEnd);
  auto grp = ccdb.get<o2::parameters::GRPMagField>("GLO/Config/GRPMagField");
  auto geom = ccdb.get<TGeoManager>("GLO/Config/GeometryAligned");

  // prepare track fitting and extrapolation to vertex (0,0,0)
  o2::base::Propagator::initFieldFromGRP(grp);
  TrackExtrap::setField();
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
bool ClustersToTrack(const gsl::span<const Cluster> trackClusters, TrackStruct& track)
{
  /// fill the internal track with the associated clusters and fit it
  /// return false if the fit fails

  static TrackFitter trackFitter{};
  trackFitter.smoothTracks(true);
  TrackExtrap::useExtrapV2();

  for (const auto& cluster : trackClusters) {
    track.track.createParamAtCluster(cluster);
  }

  try {
    trackFitter.fit(track.track);
    return true;
  } catch (exception const& e) {
    return false;
  }
}

//_________________________________________________________________________________________________
bool ExtrapToVertex(TrackStruct& track)
{
  /// compute the track parameters at vertex, at DCA and at the end of the absorber

  TrackExtrap::useExtrapV2(false);

  // extrapolate to vertex
  track.paramAtVtx = track.track.first();
  if (!TrackExtrap::extrapToVertex(track.paramAtVtx, 0., 0., 0., 0., 0.)) {
    return false;
  }

  // extrapolate to DCA
  TrackParam paramAtDCA(track.track.first());
  if (!TrackExtrap::extrapToVertexWithoutBranson(paramAtDCA, 0.)) {
    return false;
  }
  double dcaX = paramAtDCA.getNonBendingCoor();
  double dcaY = paramAtDCA.getBendingCoor();
  track.dca = TMath::Sqrt(dcaX * dcaX + dcaY * dcaY);

  // extrapolate to the end of the absorber
  TrackParam paramAtAbs(track.track.first());
  if (!TrackExtrap::extrapToZ(paramAtAbs, -505.)) {
    return false;
  }
  double xAbs = paramAtAbs.getNonBendingCoor();
  double yAbs = paramAtAbs.getBendingCoor();
  track.rAbs = TMath::Sqrt(xAbs * xAbs + yAbs * yAbs);

  return true;
}

//_________________________________________________________________________________________________
bool IsSelected(const TrackStruct& track)
{
  /// apply standard track selections + pDCA + p > 10 GeV/c

  static const double sigmaPDCA23 = 80.;
  static const double sigmaPDCA310 = 54.;
  static const double nSigmaPDCA = 6.;
  static const double relPRes = 0.0004;
  static const double slopeRes = 0.0005;

  double thetaAbs = TMath::ATan(track.rAbs / 505.) * TMath::RadToDeg();
  if (thetaAbs < 2. || thetaAbs > 10.) {
    return false;
  }

  double p = track.paramAtVtx.p();
  if (p < 10.) {
    return false;
  }

  double pz = track.paramAtVtx.pz();
  double eta = 0.5 * log((p + pz) / (p - pz));
  if (eta < -4. || eta > -2.5) {
    return false;
  }

  double pDCA = track.track.first().p() * track.dca;
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
bool IsMonoCathode(const gsl::span<const Digit> digits)
{
  /// return true if all the digits are on the same cathode

  const auto& segmentation = o2::mch::mapping::segmentation(digits[0].getDetID());

  bool hasBendingPad(false);
  bool hasNonBendingPad(false);

  for (const auto& digit : digits) {
    if (segmentation.isBendingPad(digit.getPadID())) {
      if (hasNonBendingPad) {
        return false;
      }
      hasBendingPad = true;
    } else {
      if (hasBendingPad) {
        return false;
      }
      hasNonBendingPad = true;
    }
  }

  return true;
}

//_________________________________________________________________________________________________
std::pair<double, double> GetCharge(const gsl::span<const Digit> digits)
{
  /// return the total charge of the digits on both cathodes

  const auto& segmentation = o2::mch::mapping::segmentation(digits[0].getDetID());

  std::pair<double, double> charge{0., 0.};

  for (const auto& digit : digits) {
    if (segmentation.isBendingPad(digit.getPadID())) {
      charge.second += digit.getADC();
    } else {
      charge.first += digit.getADC();
    }
  }

  return charge;
}

//_________________________________________________________________________________________________
std::unique_ptr<Cluster> COG(const gsl::span<const Digit> digits)
{
  /// return a cluster positioned at the center of gravity of the digits
  /// the weight of each digit is given by its ADC charge
  /// for bi-cathode clusters, x (y) position is given by digits in the non-bending (bending) plane
  /// note: the pad size is constant in the (non-)bending direction on the (non-)bending plane

  const auto& segmentation = o2::mch::mapping::segmentation(digits[0].getDetID());

  double x[2] = {0., 0.};
  double y[2] = {0., 0.};
  double q[2] = {0., 0.};

  for (const auto& digit : digits) {
    int i = segmentation.isBendingPad(digit.getPadID()) ? 1 : 0;
    x[i] += digit.getADC() * segmentation.padPositionX(digit.getPadID());
    y[i] += digit.getADC() * segmentation.padPositionY(digit.getPadID());
    q[i] += digit.getADC();
  }

  double xCl = (q[0] > 0.) ? x[0] / q[0] : x[1] / q[1];
  double yCl = (q[1] > 0.) ? y[1] / q[1] : y[0] / q[0];

  double ex = (q[0] > 0. || digits[0].getDetID() < 500) ? 0.2 : 10.;
  double ey = (q[1] > 0. || digits[0].getDetID() < 500) ? 0.2 : 10.;

  return MakeCluster(digits[0].getDetID(), xCl, yCl, ex, ey);
}

//_________________________________________________________________________________________________
std::unique_ptr<Cluster> COG2(const gsl::span<const Digit> digits)
{
  /// return a cluster positioned at the center of gravity of the digits
  /// the weight of each digit is given by its ADC charge / (its size / sqrt(12))^2
  /// note: the pad size is constant in the (non-)bending direction on the (non-)bending plane
  /// the purpose of weighting by the pad size is to combine the digits from both planes

  const auto& segmentation = o2::mch::mapping::segmentation(digits[0].getDetID());

  double x(0.);
  double y(0.);
  double wx(0.);
  double wy(0.);

  for (const auto& digit : digits) {
    double wxi = digit.getADC() / std::pow(segmentation.padSizeX(digit.getPadID()), 2.);
    double wyi = digit.getADC() / std::pow(segmentation.padSizeY(digit.getPadID()), 2.);
    x += wxi * segmentation.padPositionX(digit.getPadID());
    y += wyi * segmentation.padPositionY(digit.getPadID());
    wx += wxi;
    wy += wyi;
  }

  return MakeCluster(digits[0].getDetID(), x / wx, y / wy, 0.2, 0.2);
}

//_________________________________________________________________________________________________
std::unique_ptr<Cluster> MakeCluster(int de, float x, float y, float ex, float ey)
{
  /// create a cluster in the global coordinate system from the local position on the DE

  static auto transformation = o2::mch::geo::transformationFromTGeoManager(*gGeoManager);

  o2::math_utils::Point3D<float> local{x, y, 0.};
  auto global = transformation(de)(local);

  return std::make_unique<Cluster>(Cluster{global.x(), global.y(), global.z(), ex, ey, 0, 0, 0});
}

//_________________________________________________________________________________________________
void CreateResiduals(std::vector<TH1*>& histos, const char* extension, double range)
{
  /// create histograms of cluster-track residuals

  int nBins = 1001;
  double binWidth = 2. * range / (nBins - 1);
  double min = -range - 0.5 * binWidth;
  double max = range + 0.5 * binWidth;

  if (histos.size() == 0) {
    for (int iSt = 1; iSt <= 5; ++iSt) {
      histos.emplace_back(new TH1F(Form("resX%sSt%d", extension, iSt),
                                   Form("#DeltaX Station %d;#DeltaX (cm)", iSt), nBins, min, max));
      histos.emplace_back(new TH1F(Form("resY%sSt%d", extension, iSt),
                                   Form("#DeltaY Station %d;#DeltaY (cm)", iSt), nBins, min, max));
    }
    histos.emplace_back(new TH1F(Form("resX%s", extension), "#DeltaX;#DeltaX (cm)", nBins, min, max));
    histos.emplace_back(new TH1F(Form("resY%s", extension), "#DeltaY;#DeltaY (cm)", nBins, min, max));
  }

  for (auto h : histos) {
    h->SetDirectory(0);
  }
}

//_________________________________________________________________________________________________
void FillResiduals(const TrackParam& param, const Cluster& cluster, std::vector<TH1*>& histos)
{
  /// fill histograms of cluster-track residuals

  double dx = cluster.getX() - param.getNonBendingCoor();
  double dy = cluster.getY() - param.getBendingCoor();

  histos[cluster.getChamberId() / 2 * 2]->Fill(dx);
  histos[cluster.getChamberId() / 2 * 2 + 1]->Fill(dy);
  histos[10]->Fill(dx);
  histos[11]->Fill(dy);
}

//_________________________________________________________________________________________________
void FillResiduals(const Cluster& cluster1, const Cluster& cluster2, std::vector<TH1*>& histos)
{
  /// fill histograms of cluster-cluster residuals

  double dx = cluster2.getX() - cluster1.getX();
  double dy = cluster2.getY() - cluster1.getY();

  histos[cluster1.getChamberId() / 2 * 2]->Fill(dx);
  histos[cluster1.getChamberId() / 2 * 2 + 1]->Fill(dy);
  histos[10]->Fill(dx);
  histos[11]->Fill(dy);
}

//_________________________________________________________________________________________________
void DrawResiduals(std::vector<TH1*>& histos)
{
  /// draw cluster-cluster residuals

  int nPadsx = (histos.size() + 1) / 2;
  TCanvas* c = new TCanvas("clclResiduals", "cluster-cluster residuals", 10, 10, nPadsx * 300, 600);
  c->Divide(nPadsx, 2);
  int i(0);
  for (const auto& h : histos) {
    c->cd((i % 2) * nPadsx + i / 2 + 1);
    gPad->SetLogy();
    h->Draw();
    ++i;
  }
}

//_________________________________________________________________________________________________
void DrawResiduals(std::vector<TH1*>& oldHistos, std::vector<TH1*>& newHistos)
{
  /// draw cluster-track residuals

  TGraphErrors* g[2][2] = {nullptr};
  const char* dir[2] = {"X", "Y"};
  int color[2] = {4, 2};
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 2; ++j) {
      g[i][j] = new TGraphErrors(6);
      g[i][j]->SetName(Form("sigma%s%d", dir[i], j));
      g[i][j]->SetTitle(Form("#sigma%s per station;station ID;#sigma%s (cm)", dir[i], dir[i]));
      g[i][j]->SetMarkerStyle(kFullDotLarge);
      g[i][j]->SetMarkerColor(color[j]);
      g[i][j]->SetLineColor(color[j]);
    }
  }

  int nPadsx = (oldHistos.size() + 1) / 2;
  TCanvas* c = new TCanvas("cltrResiduals", "cluster-track residuals", 10, 10, nPadsx * 300, 600);
  c->Divide(nPadsx, 2);
  int i(0);
  for (const auto& h : oldHistos) {
    c->cd((i % 2) * nPadsx + i / 2 + 1);
    gPad->SetLogy();
    h->SetLineColor(color[0]);
    h->Draw();
    auto res1 = GetSigma(h, color[0]);
    g[i % 2][0]->SetPoint(i / 2, i / 2 + 1., res1.first);
    g[i % 2][0]->SetPointError(i / 2, 0., res1.second);
    newHistos[i]->SetLineColor(color[1]);
    newHistos[i]->Draw("sames");
    auto res2 = GetSigma(newHistos[i], color[1]);
    g[i % 2][1]->SetPoint(i / 2, i / 2 + 1., res2.first);
    g[i % 2][1]->SetPointError(i / 2, 0., res2.second);
    printf("%s: %f ± %f cm --> %f ± %f cm\n", h->GetName(), res1.first, res1.second, res2.first, res2.second);
    ++i;
  }

  TCanvas* c2 = new TCanvas("sigma", "sigma", 10, 10, 600, 300);
  c2->Divide(2, 1);
  c2->cd(1);
  g[0][0]->Draw("ap");
  g[0][1]->Draw("p");
  c2->cd(2);
  g[1][0]->Draw("ap");
  g[1][1]->Draw("p");

  // add a legend
  TLegend* lHist = new TLegend(0.2, 0.65, 0.4, 0.8);
  lHist->SetFillStyle(0);
  lHist->SetBorderSize(0);
  lHist->AddEntry(oldHistos[0], "old", "l");
  lHist->AddEntry(newHistos[0], "new", "l");
  c->cd(1);
  lHist->Draw("same");
}

//_________________________________________________________________________________________________
void DrawRatios(std::vector<TH1*>& oldHistos, std::vector<TH1*>& newHistos)
{
  /// draw ratios of cluster-track residuals

  int nPadsx = (oldHistos.size() + 1) / 2;
  TCanvas* c = new TCanvas("ratio", "new/old ratio", 10, 10, nPadsx * 300, 600);
  c->Divide(nPadsx, 2);
  int i(0);
  for (const auto& h : newHistos) {
    c->cd((i % 2) * nPadsx + i / 2 + 1);
    TH1* hRat = new TH1F(*static_cast<TH1F*>(h));
    hRat->SetDirectory(0);
    hRat->Divide(oldHistos[i]);
    hRat->SetStats(false);
    hRat->SetLineColor(2);
    hRat->Draw();
    hRat->GetXaxis()->SetRangeUser(-0.5, 0.5);
    ++i;
  }
}

//_________________________________________________________________________________________________
std::pair<double, double> GetSigma(TH1* h, int color)
{
  /// get the dispersion of the histo

  static TF1* fCrystalBall = new TF1("CrystalBall", CrystalBallSymmetric, -2., 2., 5);
  fCrystalBall->SetLineColor(color);

  if (h->GetEntries() < 100.) {
    return std::make_pair(0., 0.);
  }

  double sigmaTrk = 0.2;   // 2 mm
  double sigmaTrkCut = 4.; // 4 sigma

  // first fit
  double xMin = -0.5 * sigmaTrkCut * sigmaTrk;
  double xMax = 0.5 * sigmaTrkCut * sigmaTrk;
  fCrystalBall->SetRange(xMin, xMax);
  fCrystalBall->SetParameters(h->GetEntries(), 0., 0.1, 2., 1.5);
  fCrystalBall->SetParLimits(1, xMin, xMax);
  fCrystalBall->SetParLimits(2, 0., 1.);
  fCrystalBall->FixParameter(3, 1.e6);
  h->Fit(fCrystalBall, "RNQ");

  // second fit
  fCrystalBall->ReleaseParameter(3);
  fCrystalBall->SetParameter(3, 2.);
  fCrystalBall->SetParameter(4, 1.5);
  h->Fit(fCrystalBall, "RQ");

  return std::make_pair(fCrystalBall->GetParameter(2), fCrystalBall->GetParError(2));
}

//_________________________________________________________________________________________________
double CrystalBallSymmetric(double* xx, double* par)
{
  /// Crystal Ball definition

  /// par[0] = Normalization
  /// par[1] = mean
  /// par[2] = sigma
  /// par[3] = alpha = alpha'
  /// par[4] = n = n'

  double tp = std::fabs((xx[0] - par[1]) / par[2]);

  double absAlpha = std::fabs(par[3]);
  double ap = std::pow(par[4] / absAlpha, par[4]) * std::exp(-0.5 * absAlpha * absAlpha);
  double bp = par[4] / absAlpha - absAlpha;

  if (tp < absAlpha) {
    return par[0] * (std::exp(-0.5 * tp * tp)); // gaussian core
  } else {
    return par[0] * (ap / std::pow(bp + tp, par[4])); // left and right tails
  }
}
