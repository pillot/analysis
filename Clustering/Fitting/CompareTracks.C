#include <cmath>
#include <vector>

#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TMath.h>
#include <TStyle.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

#include "MCHBase/TrackBlock.h"

#include "DataUtils.h"
#include "TrackUtils.h"

using o2::mch::TrackParamStruct;

void CreateHistosAtVertex(std::vector<TH1*>& histos, const char* extension);
void FillHistosAtVertex(const char* inFile, std::vector<TH1*>& histos);
void FillHistosAtVertex(const TrackLite& track, std::vector<TH1*>& histos);
void NormHistosAtVertex(std::vector<TH1*>& histos);
void DrawHistosAtVertex(std::vector<TH1*>& histos1, std::vector<TH1*>& histos2);

//_________________________________________________________________________________________________
void CompareTracks(const char* inFile1, const char* inFile2)
{
  /// compare the track parameters at vertex between the 2 inFile

  std::vector<TH1*> histosAtVertex1{};
  CreateHistosAtVertex(histosAtVertex1, "1");
  FillHistosAtVertex(inFile1, histosAtVertex1);

  std::vector<TH1*> histosAtVertex2{};
  CreateHistosAtVertex(histosAtVertex2, "2");
  FillHistosAtVertex(inFile2, histosAtVertex2);

  gStyle->SetOptStat(1);

  NormHistosAtVertex(histosAtVertex1);
  NormHistosAtVertex(histosAtVertex2);

  DrawHistosAtVertex(histosAtVertex1, histosAtVertex2);
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
void FillHistosAtVertex(const char* inFile, std::vector<TH1*>& histos)
{
  /// fill histograms of track parameters from data in inFile

  auto [dataFileIn, dataReader] = LoadData(inFile, "data");
  TTreeReaderValue<TrackLite> track(*dataReader, "track");

  while (dataReader->Next()) {
    FillHistosAtVertex(*track, histos);
  }

  dataFileIn->Close();
}

//_________________________________________________________________________________________________
void FillHistosAtVertex(const TrackLite& track, std::vector<TH1*>& histos)
{
  /// fill single muon histograms at vertex

  static constexpr double pi = 3.14159265358979323846;

  double thetaAbs = TMath::ATan(track.rAbs / 505.) * TMath::RadToDeg();
  double pDCA = track.pUncorr * track.dca;
  double p = sqrt(track.param.px * track.param.px + track.param.py * track.param.py + track.param.pz * track.param.pz);

  histos[0]->Fill(sqrt(track.param.px * track.param.px + track.param.py * track.param.py));
  histos[1]->Fill(0.5 * log((p + track.param.pz) / (p - track.param.pz)));
  histos[2]->Fill(180. + atan2(-track.param.py, -track.param.px) / pi * 180.);
  histos[3]->Fill(track.rAbs);
  histos[4]->Fill(p);
  histos[5]->Fill(track.dca);
  if (thetaAbs < 3) {
    histos[6]->Fill(pDCA);
  } else {
    histos[7]->Fill(pDCA);
  }
  histos[8]->Fill(track.nClusters);
  histos[9]->Fill(track.chi2);
  histos[10]->Fill(track.matchChi2);
}

//_________________________________________________________________________________________________
void NormHistosAtVertex(std::vector<TH1*>& histos)
{
  /// normalize histograms

  for (auto& h : histos) {
    h->Scale(1. / h->GetEntries());
  }
}

//_________________________________________________________________________________________________
void DrawHistosAtVertex(std::vector<TH1*>& histos1, std::vector<TH1*>& histos2)
{
  /// Draw single muon histograms at vertex and differences between the 2 inputs

  // find the optimal number of pads
  int nPadsx(1), nPadsy(1);
  while ((int)histos1.size() > nPadsx * nPadsy) {
    if (nPadsx == nPadsy) {
      ++nPadsx;
    } else {
      ++nPadsy;
    }
  }

  // draw histograms
  TCanvas* cHist = new TCanvas("histos", "histos", 10, 10, TMath::Max(nPadsx * 300, 1200), TMath::Max(nPadsy * 300, 900));
  cHist->Divide(nPadsx, nPadsy);
  for (int i = 0; i < (int)histos1.size(); ++i) {
    cHist->cd(i + 1);
    gPad->SetLogy();
    histos1[i]->SetStats(false);
    histos1[i]->SetLineColor(4);
    histos1[i]->Draw("hist");
    histos2[i]->SetLineColor(2);
    histos2[i]->Draw("same hist");
  }

  // add a legend
  TLegend* lHist = new TLegend(0.5, 0.65, 0.9, 0.8);
  lHist->SetFillStyle(0);
  lHist->SetBorderSize(0);
  lHist->AddEntry(histos1[0], Form("%g tracks in file 1", histos1[0]->GetEntries()), "l");
  lHist->AddEntry(histos2[0], Form("%g tracks in file 2", histos2[0]->GetEntries()), "l");
  cHist->cd(1);
  lHist->Draw("same");

  // draw differences
  TCanvas* cDiff = new TCanvas("differences", "histos2 - histos1", 10, 10, TMath::Max(nPadsx * 300, 1200), TMath::Max(nPadsy * 300, 900));
  cDiff->Divide(nPadsx, nPadsy);
  for (int i = 0; i < (int)histos1.size(); ++i) {
    cDiff->cd(i + 1);
    TH1F* hDiff = static_cast<TH1F*>(histos2[i]->Clone());
    hDiff->SetDirectory(0);
    hDiff->Add(histos1[i], -1.);
    hDiff->SetStats(false);
    hDiff->SetLineColor(2);
    hDiff->Draw("hist");
  }

  // draw ratios
  TCanvas* cRat = new TCanvas("ratios", "histos2 / histos1", 10, 10, TMath::Max(nPadsx * 300, 1200), TMath::Max(nPadsy * 300, 900));
  cRat->Divide(nPadsx, nPadsy);
  for (int i = 0; i < (int)histos1.size(); ++i) {
    cRat->cd(i + 1);
    TH1F* hRat = static_cast<TH1F*>(histos2[i]->Clone());
    hRat->SetDirectory(0);
    hRat->Divide(histos1[i]);
    hRat->SetStats(false);
    hRat->SetLineColor(2);
    hRat->Draw();
  }
}
