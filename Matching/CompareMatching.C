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
void CreateHistosAtVertex(std::vector<TH1*>& histos, const char* extension);
void FillHistosAtVertex(const mch::TrackMCH& track, const TrackAtVtxStruct& trackAtVtx, double matchChi2, std::vector<TH1*>& histos);
void DrawHistosAtVertex(std::vector<TH1*> histos[2]);
void DrawComparisonsAtVertex(std::vector<TH1*> histos[4]);

//_________________________________________________________________________________________________
void CompareMatching(float l3Current, float dipoleCurrent, std::string vtxFileName, const char* mchFileName,
                     const char* midFileName1, const char* muonFileName1,
                     const char* midFileName2, const char* muonFileName2,
                     double precision = 1.e-4)
{
  /// compare the matching of the same MCH tracks with different MID track and/or different matching algorithm
  /// O2 tracking need to be loaded before: gSystem->Load("libO2MCHTracking")

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

        // find the corresponding MCH-MID matched tracks
        auto muon1 = FindMuon(iMCHTrack, *muonTracks1);
        auto muon2 = FindMuon(iMCHTrack, *muonTracks2);

        // fill histograms
        FillHistosAtVertex(mchTrack, mchTrackAtVtx, 0., comparisonsAtVertex[0]);
        if (muon1) {
          FillHistosAtVertex(mchTrack, mchTrackAtVtx, muon1->getMatchChi2OverNDF(), histosAtVertex[0]);
        }
        if (muon2) {
          FillHistosAtVertex(mchTrack, mchTrackAtVtx, muon2->getMatchChi2OverNDF(), histosAtVertex[1]);
        }
        if (muon1 && muon2) {
          FillHistosAtVertex(mchTrack, mchTrackAtVtx, muon2->getMatchChi2OverNDF(), comparisonsAtVertex[1]);
          if (abs(muon2->getMatchChi2OverNDF() - muon1->getMatchChi2OverNDF()) > precision) {
            std::cout << "chi2 difference = " << muon2->getMatchChi2OverNDF() - muon1->getMatchChi2OverNDF() << std::endl;
          }
        } else if (muon2) {
          FillHistosAtVertex(mchTrack, mchTrackAtVtx, muon2->getMatchChi2OverNDF(), comparisonsAtVertex[2]);
          std::cout << "additional matching: chi2 = " << muon2->getMatchChi2OverNDF() << std::endl;
        } else if (muon1) {
          FillHistosAtVertex(mchTrack, mchTrackAtVtx, muon1->getMatchChi2OverNDF(), comparisonsAtVertex[3]);
          std::cout << "missing matching: chi2 = " << muon1->getMatchChi2OverNDF() << std::endl;
        }
      }
    }
  }

  DrawHistosAtVertex(histosAtVertex);
  DrawComparisonsAtVertex(comparisonsAtVertex);

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
void CreateHistosAtVertex(std::vector<TH1*>& histos, const char* extension)
{
  /// create single muon and dimuon histograms at vertex

  histos.emplace_back(new TH1F(Form("pT%s", extension), "pT;p_{T} (GeV/c)", 300, 0., 30.));
  histos.emplace_back(new TH1F(Form("eta%s", extension), "eta;eta", 200, -4.5, -2.));
  histos.emplace_back(new TH1F(Form("phi%s", extension), "phi;phi", 360, 0., 360.));
  histos.emplace_back(new TH1F(Form("rAbs%s", extension), "rAbs;R_{abs} (cm)", 1000, 0., 100.));
  histos.emplace_back(new TH1F(Form("dca%s", extension), "DCA;DCA (cm)", 500, 0., 500.));
  histos.emplace_back(new TH1F(Form("pDCA23%s", extension), "pDCA for #theta_{abs} < 3#circ;pDCA (GeV.cm/c)", 2500, 0., 5000.));
  histos.emplace_back(new TH1F(Form("pDCA310%s", extension), "pDCA for #theta_{abs} > 3#circ;pDCA (GeV.cm/c)", 2500, 0., 5000.));
  histos.emplace_back(new TH1F(Form("nClusters%s", extension), "number of clusters per track;n_{clusters}", 20, 0., 20.));
  histos.emplace_back(new TH1F(Form("chi2%s", extension), "normalized #chi^{2};#chi^{2} / ndf", 500, 0., 50.));
  histos.emplace_back(new TH1F(Form("matchChi2%s", extension), "normalized matching #chi^{2};#chi^{2} / ndf", 500, 0., 50.));

  for (auto h : histos) {
    h->SetDirectory(0);
  }
}

//_________________________________________________________________________________________________
void FillHistosAtVertex(const mch::TrackMCH& track, const TrackAtVtxStruct& trackAtVtx, double matchChi2, std::vector<TH1*>& histos)
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
  histos[3]->Fill(trackAtVtx.rAbs);
  histos[4]->Fill(trackAtVtx.dca);
  if (thetaAbs < 3) {
    histos[5]->Fill(pDCA);
  } else {
    histos[6]->Fill(pDCA);
  }
  histos[7]->Fill(track.getNClusters());
  histos[8]->Fill(track.getChi2OverNDF());
  histos[9]->Fill(matchChi2);
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
