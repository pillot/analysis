#include <cmath>
#include <utility>
#include <vector>

#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

#include "DataFormatsMCH/Cluster.h"
#include "MCHBase/TrackBlock.h"

#include "DataUtils.h"

using o2::mch::Cluster;
using o2::mch::TrackParamStruct;

void CreateResiduals(std::vector<TH1*>& histos, const char* extension, double range);
void FillResiduals(const char* inFile, bool useNewClusters, std::vector<TH1*>& histos);
void FillResiduals(const TrackParamStruct& param, const Cluster& cluster, std::vector<TH1*>& histos);
void FillResiduals(const Cluster& cluster1, const Cluster& cluster2, std::vector<TH1*>& histos);
void NormResiduals(std::vector<TH1*>& histos);
void DrawResiduals(std::vector<TH1*>& histos);
void DrawResiduals(std::vector<TH1*>& oldHistos, std::vector<TH1*>& newHistos);
void DrawRatios(std::vector<TH1*>& oldHistos, std::vector<TH1*>& newHistos, int rebin);
std::pair<double, double> GetSigma(TH1* h, int color);
double CrystalBallSymmetric(double* xx, double* par);

//_________________________________________________________________________________________________
void CompareClusters(const char* inFile = "newclusters.root")
{
  /// compare the original and new clusters from inFile
  /// and their corresponding cluster-track residuals

  auto [dataFileIn, dataReader] = LoadData(inFile, "data");
  TTreeReaderValue<TrackParamStruct> trackParam(*dataReader, "trackParameters");
  TTreeReaderValue<Cluster> oldCluster(*dataReader, "clusters");
  TTreeReaderValue<Cluster> newCluster(*dataReader, "newClusters");

  std::vector<TH1*> oldResiduals{};
  CreateResiduals(oldResiduals, "old", 2.);
  std::vector<TH1*> newResiduals{};
  CreateResiduals(newResiduals, "new", 2.);
  std::vector<TH1*> clclResiduals{};
  CreateResiduals(clclResiduals, "clcl", 0.2);

  while (dataReader->Next()) {
    FillResiduals(*trackParam, *oldCluster, oldResiduals);
    FillResiduals(*trackParam, *newCluster, newResiduals);
    FillResiduals(*oldCluster, *newCluster, clclResiduals);
  }

  gStyle->SetOptStat(1);

  DrawResiduals(clclResiduals);
  DrawResiduals(oldResiduals, newResiduals);
  DrawRatios(oldResiduals, newResiduals, 5);

  dataFileIn->Close();
}

//_________________________________________________________________________________________________
void CompareClusters(const char* inFile1, const char* inFile2, bool compareNewClusters = false)
{
  /// compare the cluster-track residuals between the 2 inFile
  /// if compareNewClusters = true, compare the new clusters

  std::vector<TH1*> residuals1{};
  CreateResiduals(residuals1, "file1", 2.);
  FillResiduals(inFile1, compareNewClusters, residuals1);

  std::vector<TH1*> residuals2{};
  CreateResiduals(residuals2, "file2", 2.);
  FillResiduals(inFile2, compareNewClusters, residuals2);

  gStyle->SetOptStat(1);

  NormResiduals(residuals1);
  NormResiduals(residuals2);

  DrawResiduals(residuals1, residuals2);
  DrawRatios(residuals1, residuals2, 5);
}

//_________________________________________________________________________________________________
void CreateResiduals(std::vector<TH1*>& histos, const char* extension, double range)
{
  /// create histograms of cluster-track residuals

  int nBins = 1000;

  if (histos.size() == 0) {
    for (int iSt = 1; iSt <= 5; ++iSt) {
      histos.emplace_back(new TH1F(Form("resX%sSt%d", extension, iSt),
                                   Form("#DeltaX Station %d;#DeltaX (cm)", iSt), nBins, -range, range));
      histos.emplace_back(new TH1F(Form("resY%sSt%d", extension, iSt),
                                   Form("#DeltaY Station %d;#DeltaY (cm)", iSt), nBins, -range, range));
    }
    histos.emplace_back(new TH1F(Form("resX%s", extension), "#DeltaX;#DeltaX (cm)", nBins, -range, range));
    histos.emplace_back(new TH1F(Form("resY%s", extension), "#DeltaY;#DeltaY (cm)", nBins, -range, range));
  }

  for (auto h : histos) {
    h->SetDirectory(0);
  }
}

//_________________________________________________________________________________________________
void FillResiduals(const char* inFile, bool useNewClusters, std::vector<TH1*>& histos)
{
  /// fill histograms of cluster-track residuals from data in inFile
  /// if useNewClusters = true, use the new clusters

  auto [dataFileIn, dataReader] = LoadData(inFile, "data");
  TTreeReaderValue<TrackParamStruct> trackParam(*dataReader, "trackParameters");
  TTreeReaderValue<Cluster> cluster(*dataReader, useNewClusters ? "newClusters" : "clusters");

  while (dataReader->Next()) {
    FillResiduals(*trackParam, *cluster, histos);
  }

  dataFileIn->Close();
}

//_________________________________________________________________________________________________
void FillResiduals(const TrackParamStruct& param, const Cluster& cluster, std::vector<TH1*>& histos)
{
  /// fill histograms of cluster-track residuals

  double dx = cluster.getX() - param.x;
  double dy = cluster.getY() - param.y;

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
void NormResiduals(std::vector<TH1*>& histos)
{
  /// normalize histograms

  for (auto& h : histos) {
    h->Scale(1. / h->GetEntries());
  }
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
    // h->GetXaxis()->SetRangeUser(-1., 1.);
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
void DrawRatios(std::vector<TH1*>& oldHistos, std::vector<TH1*>& newHistos, int rebin)
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
    hRat->Rebin(rebin);
    auto* h1 = static_cast<TH1*>(oldHistos[i]->Clone());
    h1->Rebin(rebin);
    hRat->Divide(h1);
    delete h1;
    hRat->SetStats(false);
    hRat->SetLineColor(2);
    hRat->Draw();
    hRat->GetXaxis()->SetRangeUser(-0.5, 0.5);
    if (i > 9) {
      hRat->GetYaxis()->SetRangeUser(0.95, 1.05);
    } else {
      hRat->GetYaxis()->SetRangeUser(0.9, 1.1);
    }
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
  fCrystalBall->SetParameter(3, 1.);
  fCrystalBall->SetParameter(4, 2.5);
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
