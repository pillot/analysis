#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <tuple>
#include <vector>

#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TGeoGlobalMagField.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TLegend.h>

#include "Framework/Logger.h"
#include "DetectorsBase/GeometryManager.h"
#include "Field/MagneticField.h"
#include "MCHTracking/TrackParam.h"
#include "MCHTracking/TrackExtrap.h"
#include "ReconstructionDataFormats/TrackMCHMID.h"
#include "DataFormatsMCH/ROFRecord.h"
#include "DataFormatsMCH/TrackMCH.h"
#include "DataFormatsMID/Track.h"

using namespace o2;

struct VertexStruct {
  double x;
  double y;
  double z;
};

struct TrackAtVtxStruct {
  mch::TrackParam paramAtVertex{};
  double dca = 0.;
  double rAbs = 0.;
};

constexpr double pi() { return 3.14159265358979323846; }
void LoadVertices(std::string vtxFileName, std::vector<VertexStruct>& vertices);
void PrepareTrackExtrapolation(float l3Current, float dipoleCurrent);
std::tuple<TFile*, TTreeReader*> LoadData(const char* fileName, const char* treeName);
const dataformats::TrackMCHMID* FindMuon(uint32_t iMCHTrack, const std::vector<dataformats::TrackMCHMID>& muonTracks);
bool ExtrapToVertex(const mch::TrackMCH& track, VertexStruct& vertex, TrackAtVtxStruct& trackAtVtx);
bool IsSelected(const mch::TrackMCH& track, const TrackAtVtxStruct& trackAtVtx);
void CreateHistosAtVertex(std::vector<TH1*>& histos, const char* extension);
void FillHistosAtVertex(const mch::TrackMCH& track, const TrackAtVtxStruct& trackAtVtx, std::vector<TH1*>& histos);
void DrawHistosAtVertex(std::vector<TH1*> histos[2]);
void DrawComparisonsAtVertex(std::vector<TH1*> histos[4]);
void CreateHistosAtMID(std::vector<TH1*>& histos);
void FillHistosAtMID(double chi21, double matchChi21, double chi22, double matchChi22, std::vector<TH1*>& histos);
void DrawHistosAtMID(std::vector<TH1*>& histos);

//_________________________________________________________________________________________________
void CompareMatching(float l3Current, float dipoleCurrent, std::string vtxFileName, const char* mchFileName,
                     const char* midFileName1, const char* muonFileName1,
                     const char* midFileName2, const char* muonFileName2,
                     bool applyTrackSelection = false, bool printDifferences = true, double precision = 1.e-4)
{
  /// compare the matching of the same MCH tracks with different MID track and/or different matching algorithm

  // get vertices and prepare track extrapolation
  std::vector<VertexStruct> vertices{};
  LoadVertices(vtxFileName, vertices);
  PrepareTrackExtrapolation(l3Current, dipoleCurrent);

  // read the input data
  auto [fMCH, mchReader] = LoadData(mchFileName, "o2sim");
  TTreeReaderValue<std::vector<mch::ROFRecord>> mchROFs = {*mchReader, "trackrofs"};
  TTreeReaderValue<std::vector<mch::TrackMCH>> mchTracks = {*mchReader, "tracks"};
  auto [fMID1, midReader1] = LoadData(midFileName1, "midreco");
  TTreeReaderValue<std::vector<char>> midTracks1 = {*midReader1, "MIDTrack"};
  auto [fMUON1, muonReader1] = LoadData(muonFileName1, "o2sim");
  TTreeReaderValue<std::vector<dataformats::TrackMCHMID>> muonTracks1 = {*muonReader1, "tracks"};
  auto [fMID2, midReader2] = LoadData(midFileName2, "midreco");
  TTreeReaderValue<std::vector<char>> midTracks2 = {*midReader2, "MIDTrack"};
  auto [fMUON2, muonReader2] = LoadData(muonFileName2, "o2sim");
  TTreeReaderValue<std::vector<dataformats::TrackMCHMID>> muonTracks2 = {*muonReader2, "tracks"};

  std::vector<TH1*> comparisonsAtVertex[4] = {{}, {}, {}, {}};
  CreateHistosAtVertex(comparisonsAtVertex[0], "all");
  CreateHistosAtVertex(comparisonsAtVertex[1], "matched");
  CreateHistosAtVertex(comparisonsAtVertex[2], "additional");
  CreateHistosAtVertex(comparisonsAtVertex[3], "missing");
  std::vector<TH1*> histosAtVertex[2] = {{}, {}};
  CreateHistosAtVertex(histosAtVertex[0], "1");
  CreateHistosAtVertex(histosAtVertex[1], "2");
  std::vector<TH1*> histosAtMID{};
  CreateHistosAtMID(histosAtMID);

  VertexStruct vertex{};
  TrackAtVtxStruct mchTrackAtVtx{};

  while (mchReader->Next() && midReader1->Next() && muonReader1->Next() && midReader2->Next() && muonReader2->Next()) {

    for (const auto& mchROF : *mchROFs) {

      // get the vertex corresponding to this event or use (0.,0.,0.)
      if (!vertices.empty()) {
        auto event = mchROF.getBCData().orbit;
        if (event < vertices.size()) {
          vertex = vertices[event];
        } else {
          LOG(ERROR) << "missing vertex for event" << event;
          exit(-1);
        }
      }

      for (int iMCHTrack = mchROF.getFirstIdx(); iMCHTrack <= mchROF.getLastIdx(); ++iMCHTrack) {

        // compute the track parameters at vertex
        const auto& mchTrack = (*mchTracks)[iMCHTrack];
        if (!ExtrapToVertex(mchTrack, vertex, mchTrackAtVtx)) {
          LOG(ERROR) << "track extrapolation to vertex failed";
          continue;
        }

        // apply track selection if requested
        if (applyTrackSelection && !IsSelected(mchTrack, mchTrackAtVtx)) {
          continue;
        }

        // find the corresponding MCH-MID matched tracks
        auto muon1 = FindMuon(iMCHTrack, *muonTracks1);
        auto muon2 = FindMuon(iMCHTrack, *muonTracks2);

        // fill histograms
        double chi21(-1.), matchChi21(-1.), chi22(-1.), matchChi22(-1.);
        FillHistosAtVertex(mchTrack, mchTrackAtVtx, comparisonsAtVertex[0]);
        if (muon1) {
          FillHistosAtVertex(mchTrack, mchTrackAtVtx, histosAtVertex[0]);
          auto midTrack1 = reinterpret_cast<mid::Track*>(midTracks1->data()) + muon1->getMIDRef().getIndex();
          chi21 = midTrack1->getNDF() ? midTrack1->getChi2OverNDF() : 0.;
          matchChi21 = muon1->getMatchChi2OverNDF();
        }
        if (muon2) {
          FillHistosAtVertex(mchTrack, mchTrackAtVtx, histosAtVertex[1]);
          auto midTrack2 = reinterpret_cast<mid::Track*>(midTracks2->data()) + muon2->getMIDRef().getIndex();
          chi22 = midTrack2->getNDF() ? midTrack2->getChi2OverNDF() : 0.;
          matchChi22 = muon2->getMatchChi2OverNDF();
        }
        if (muon1 && muon2) {
          FillHistosAtVertex(mchTrack, mchTrackAtVtx, comparisonsAtVertex[1]);
          if (printDifferences && abs(matchChi22 - matchChi21) > precision) {
            std::cout << "chi2 difference = " << matchChi22 - matchChi21 << std::endl;
          }
        } else if (muon2) {
          FillHistosAtVertex(mchTrack, mchTrackAtVtx, comparisonsAtVertex[2]);
          if (printDifferences) {
            std::cout << "additional matching: chi2 = " << matchChi22 << std::endl;
          }
        } else if (muon1) {
          FillHistosAtVertex(mchTrack, mchTrackAtVtx, comparisonsAtVertex[3]);
          if (printDifferences) {
            std::cout << "missing matching: chi2 = " << matchChi21 << std::endl;
          }
        }
        FillHistosAtMID(chi21, matchChi21, chi22, matchChi22, histosAtMID);
      }
    }
  }

  DrawHistosAtVertex(histosAtVertex);
  DrawComparisonsAtVertex(comparisonsAtVertex);
  DrawHistosAtMID(histosAtMID);

  fMCH->Close();
  fMID1->Close();
  fMUON1->Close();
  fMID2->Close();
  fMUON2->Close();
}

//_________________________________________________________________________________________________
void LoadVertices(std::string vtxFileName, std::vector<VertexStruct>& vertices)
{
  /// get all event vertices

  if (vtxFileName.empty()) {
    return;
  }

  ifstream inFile(vtxFileName, ios::binary);
  if (!inFile.is_open()) {
    LOG(ERROR) << "fail opening vertex file";
    exit(-1);
  }

  int event(-1), expectedEvent(-1);
  while (inFile.read(reinterpret_cast<char*>(&event), sizeof(int))) {

    if (event != ++expectedEvent) {
      LOG(ERROR) << "event " << expectedEvent << " missing";
      exit(-1);
    }

    vertices.emplace_back();
    inFile.read(reinterpret_cast<char*>(&vertices.back()), sizeof(VertexStruct));
  }
}

//_________________________________________________________________________________________________
void PrepareTrackExtrapolation(float l3Current, float dipoleCurrent)
{
  /// prepare magnetic field and geometry for track extrapolation to vertex

  auto field = field::MagneticField::createFieldMap(l3Current, dipoleCurrent, field::MagneticField::kConvLHC, false, 3500.,
                                                    "A-A", "$(O2_ROOT)/share/Common/maps/mfchebKGI_sym.root");
  TGeoGlobalMagField::Instance()->SetField(field);
  TGeoGlobalMagField::Instance()->Lock();
  mch::TrackExtrap::setField();
  mch::TrackExtrap::useExtrapV2();
  
  base::GeometryManager::loadGeometry("O2geometry.root");
  if (!gGeoManager) {
    exit(-1);
  }
}

//_________________________________________________________________________________________________
std::tuple<TFile*, TTreeReader*> LoadData(const char* fileName, const char* treeName)
{
  /// open the input file and get the intput tree

  TFile* f = TFile::Open(fileName, "READ");
  if (f->IsZombie()) {
    LOG(ERROR) << "opening file " << fileName << " failed";
    exit(-1);
  }

  TTreeReader* r = new TTreeReader(treeName, f);
  if (r->IsZombie()) {
    LOG(ERROR) << "tree " << treeName << " not found";
    exit(-1);
  }

  return std::make_tuple(f, r);
}

//_________________________________________________________________________________________________
const dataformats::TrackMCHMID* FindMuon(uint32_t iMCHTrack, const std::vector<dataformats::TrackMCHMID>& muonTracks)
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
bool ExtrapToVertex(const mch::TrackMCH& track, VertexStruct& vertex, TrackAtVtxStruct& trackAtVtx)
{
  /// compute the track parameters at vertex, at DCA and at the end of the absorber
  /// return false if the propagation fails

  // extrapolate to vertex
  trackAtVtx.paramAtVertex.setZ(track.getZ());
  trackAtVtx.paramAtVertex.setParameters(track.getParameters());
  trackAtVtx.paramAtVertex.deleteCovariances();
  if (!mch::TrackExtrap::extrapToVertex(trackAtVtx.paramAtVertex, vertex.x, vertex.y, vertex.z, 0., 0.)) {
    return false;
  }

  // extrapolate to DCA
  mch::TrackParam trackParamAtDCA(track.getZ(), track.getParameters());
  if (!mch::TrackExtrap::extrapToVertexWithoutBranson(trackParamAtDCA, vertex.z)) {
    return false;
  }
  double dcaX = trackParamAtDCA.getNonBendingCoor() - vertex.x;
  double dcaY = trackParamAtDCA.getBendingCoor() - vertex.y;
  trackAtVtx.dca = sqrt(dcaX * dcaX + dcaY * dcaY);

  // extrapolate to the end of the absorber
  mch::TrackParam trackParamAtRAbs(track.getZ(), track.getParameters());
  if (!mch::TrackExtrap::extrapToZ(trackParamAtRAbs, -505.)) {
    return false;
  }
  double xAbs = trackParamAtRAbs.getNonBendingCoor();
  double yAbs = trackParamAtRAbs.getBendingCoor();
  trackAtVtx.rAbs = sqrt(xAbs * xAbs + yAbs * yAbs);

  return true;
}

//_________________________________________________________________________________________________
bool IsSelected(const mch::TrackMCH& track, const TrackAtVtxStruct& trackAtVtx)
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

  double p = trackAtVtx.paramAtVertex.p();
  double eta = 0.5 * log((p + trackAtVtx.paramAtVertex.pz()) / (p - trackAtVtx.paramAtVertex.pz()));
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

  return true;
}

//_________________________________________________________________________________________________
void CreateHistosAtVertex(std::vector<TH1*>& histos, const char* extension)
{
  /// create single muon histograms at vertex

  histos.emplace_back(new TH1F(Form("pT%s", extension), "pT;p_{T} (GeV/c)", 300, 0., 30.));
  histos.emplace_back(new TH1F(Form("eta%s", extension), "eta;eta", 200, -4.5, -2.));
  histos.emplace_back(new TH1F(Form("phi%s", extension), "phi;phi", 360, 0., 360.));
  histos.emplace_back(new TH1F(Form("dca%s", extension), "DCA;DCA (cm)", 500, 0., 500.));
  histos.emplace_back(new TH1F(Form("pDCA23%s", extension), "pDCA for #theta_{abs} < 3#circ;pDCA (GeV.cm/c)", 2500, 0., 5000.));
  histos.emplace_back(new TH1F(Form("pDCA310%s", extension), "pDCA for #theta_{abs} > 3#circ;pDCA (GeV.cm/c)", 2500, 0., 5000.));
  histos.emplace_back(new TH1F(Form("rAbs%s", extension), "rAbs;R_{abs} (cm)", 1000, 0., 100.));
  histos.emplace_back(new TH1F(Form("nClusters%s", extension), "number of clusters per track;n_{clusters}", 20, 0., 20.));
  histos.emplace_back(new TH1F(Form("chi2%s", extension), "normalized #chi^{2};#chi^{2} / ndf", 500, 0., 50.));

  for (auto h : histos) {
    h->SetDirectory(0);
  }
}

//_________________________________________________________________________________________________
void FillHistosAtVertex(const mch::TrackMCH& track, const TrackAtVtxStruct& trackAtVtx, std::vector<TH1*>& histos)
{
  /// fill single muon histograms at vertex

  double thetaAbs = TMath::ATan(trackAtVtx.rAbs / 505.) * TMath::RadToDeg();
  double pDCA = track.getP() * trackAtVtx.dca;
  double pT = sqrt(trackAtVtx.paramAtVertex.px() * trackAtVtx.paramAtVertex.px() +
                   trackAtVtx.paramAtVertex.py() * trackAtVtx.paramAtVertex.py());
  double p = trackAtVtx.paramAtVertex.p();
  double eta = 0.5 * log((p + trackAtVtx.paramAtVertex.pz()) / (p - trackAtVtx.paramAtVertex.pz()));
  double phi = 180. + atan2(-trackAtVtx.paramAtVertex.px(), -trackAtVtx.paramAtVertex.py()) / pi() * 180.;

  histos[0]->Fill(pT);
  histos[1]->Fill(eta);
  histos[2]->Fill(phi);
  histos[3]->Fill(trackAtVtx.dca);
  if (thetaAbs < 3) {
    histos[4]->Fill(pDCA);
  } else {
    histos[5]->Fill(pDCA);
  }
  histos[6]->Fill(trackAtVtx.rAbs);
  histos[7]->Fill(track.getNClusters());
  histos[8]->Fill(track.getChi2OverNDF());
}

//_________________________________________________________________________________________________
void DrawHistosAtVertex(std::vector<TH1*> histos[2])
{
  /// draw histograms at vertex and differences between the 2 inputs

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
    cRat->cd((i / nPadsx) * nPadsx + i % nPadsx + 1);
    TH1F* hRat = static_cast<TH1F*>(histos[1][i]->Clone());
    hRat->SetDirectory(0);
    hRat->Divide(histos[0][i]);
    hRat->SetStats(false);
    hRat->SetLineColor(2);
    hRat->Draw();
  }
}

//_________________________________________________________________________________________________
void DrawComparisonsAtVertex(std::vector<TH1*> histos[4])
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
  lHist->AddEntry(histos[0][0], Form("%g all tracks", histos[0][0]->GetEntries()), "l");
  lHist->AddEntry(histos[1][0], Form("%g matched tracks", histos[1][0]->GetEntries()), "l");
  lHist->AddEntry(histos[2][0], Form("%g matched additional", histos[2][0]->GetEntries()), "l");
  lHist->AddEntry(histos[3][0], Form("%g matched missing", histos[3][0]->GetEntries()), "l");
  cHist->cd(1);
  lHist->Draw("same");
}

//_________________________________________________________________________________________________
void CreateHistosAtMID(std::vector<TH1*>& histos)
{
  /// create single muon histograms at MID

  histos.emplace_back(new TH1F("matchChi2_1all", "normalized matching #chi^{2};#chi^{2} / ndf", 500, 0., 50.));
  histos.emplace_back(new TH1F("matchChi2_2all", "normalized matching #chi^{2};#chi^{2} / ndf", 500, 0., 50.));
  histos.emplace_back(new TH1F("matchChi2_1both", "normalized matching #chi^{2};#chi^{2} / ndf", 500, 0., 50.));
  histos.emplace_back(new TH1F("matchChi2_2both", "normalized matching #chi^{2};#chi^{2} / ndf", 500, 0., 50.));
  histos.emplace_back(new TH1F("matchChi2_1add", "normalized matching #chi^{2};#chi^{2} / ndf", 500, 0., 50.));
  histos.emplace_back(new TH1F("matchChi2_2add", "normalized matching #chi^{2};#chi^{2} / ndf", 500, 0., 50.));
  histos.emplace_back(new TH2F("chi2VsmatchChi2_1", "MID #chi^{2} vs matching #chi^{2} - track 1;match #chi^{2} / ndf;MID #chi^{2} / ndf", 500, 0., 50., 500, 0., 50.));
  histos.emplace_back(new TH2F("chi2VsmatchChi2_2", "MID #chi^{2} vs matching #chi^{2} - track 2;match #chi^{2} / ndf;MID #chi^{2} / ndf", 500, 0., 50., 500, 0., 50.));

  for (auto h : histos) {
    h->SetDirectory(0);
  }
}

//_________________________________________________________________________________________________
void FillHistosAtMID(double chi21, double matchChi21, double chi22, double matchChi22, std::vector<TH1*>& histos)
{
  /// fill single muon histograms at MID

  if (matchChi21 >= 0) {
    histos[0]->Fill(matchChi21);
    histos[6]->Fill(matchChi21, chi21);
  }

  if (matchChi22 >= 0) {
    histos[1]->Fill(matchChi22);
    histos[7]->Fill(matchChi22, chi22);
  }

  if (matchChi21 >= 0 && matchChi22 >= 0) {
    histos[2]->Fill(matchChi21);
    histos[3]->Fill(matchChi22);
  } else if (matchChi21 >= 0) {
    histos[4]->Fill(matchChi21);
  } else if (matchChi22 >= 0) {
    histos[5]->Fill(matchChi22);
  }
}

//_________________________________________________________________________________________________
void DrawHistosAtMID(std::vector<TH1*>& histos)
{
  /// draw single muon histograms at MID

  TCanvas* cHist = new TCanvas("histosMID", "histosMID", 10, 10, 900, 600);
  cHist->Divide(3, 2);
  cHist->cd(1);
  gPad->SetLogy();
  histos[0]->SetStats(false);
  histos[0]->SetLineColor(4);
  histos[0]->Draw();
  histos[1]->SetLineColor(2);
  histos[1]->Draw("same");
  cHist->cd(4);
  TH1F* hRat = static_cast<TH1F*>(histos[1]->Clone());
  hRat->SetDirectory(0);
  hRat->SetTitle("track2 / track1 ratio");
  hRat->Divide(histos[0]);
  hRat->SetStats(false);
  hRat->SetLineColor(2);
  hRat->Draw();
  cHist->cd(2);
  gPad->SetLogy();
  histos[2]->SetStats(false);
  histos[2]->SetLineColor(4);
  histos[2]->Draw();
  histos[3]->SetLineColor(877);
  histos[3]->Draw("same");
  histos[4]->SetLineColor(2);
  histos[4]->Draw("same");
  histos[5]->SetLineColor(3);
  histos[5]->Draw("same");
  cHist->cd(5);
  hRat = static_cast<TH1F*>(histos[3]->Clone());
  hRat->SetDirectory(0);
  hRat->SetTitle("track2 / track1 ratio");
  hRat->Divide(histos[2]);
  hRat->SetStats(false);
  hRat->SetLineColor(2);
  hRat->Draw();
  cHist->cd(3);
  gPad->SetLogz();
  histos[6]->SetStats(false);
  histos[6]->Draw("colz");
  cHist->cd(6);
  gPad->SetLogz();
  histos[7]->SetStats(false);
  histos[7]->Draw("colz");

  TLegend* lHist1 = new TLegend(0.5, 0.75, 0.9, 0.85);
  lHist1->SetFillStyle(0);
  lHist1->SetBorderSize(0);
  lHist1->AddEntry(histos[0], Form("%g tracks in file 1", histos[0]->GetEntries()), "l");
  lHist1->AddEntry(histos[1], Form("%g tracks in file 2", histos[1]->GetEntries()), "l");
  cHist->cd(1);
  lHist1->Draw("same");

  TLegend* lHist2 = new TLegend(0.4, 0.65, 0.9, 0.85);
  lHist2->SetFillStyle(0);
  lHist2->SetBorderSize(0);
  lHist2->AddEntry(histos[2], Form("%g tracks 1 match both", histos[2]->GetEntries()), "l");
  lHist2->AddEntry(histos[3], Form("%g tracks 2 match both", histos[3]->GetEntries()), "l");
  lHist2->AddEntry(histos[4], Form("%g tracks 1 match add", histos[4]->GetEntries()), "l");
  lHist2->AddEntry(histos[5], Form("%g tracks 2 match add", histos[5]->GetEntries()), "l");
  cHist->cd(2);
  lHist2->Draw("same");
}
