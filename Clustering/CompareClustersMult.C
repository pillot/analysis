#include <algorithm>
#include <iostream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <fmt/format.h>

#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

#include "DataFormatsMCH/Cluster.h"
#include "DataFormatsMCH/ROFRecord.h"
#include "Framework/Logger.h"

using o2::mch::Cluster;
using o2::mch::ROFRecord;

std::tuple<TFile*, TTreeReader*> LoadData(const char* fileName, const char* treeName);
void CreateMultPlots(std::vector<TH1*>& histos, std::string extension);
std::pair<int, int> FillMultPlots(std::string clusterFileName, std::vector<TH1*>& histos);
void FillMultPlots(const std::vector<ROFRecord>& rofs, const std::vector<Cluster>& clusters, std::vector<TH1*>& histos);
void DrawMultPlots(std::vector<TH1*>& histos1, int nROF1, int nTF1, std::vector<TH1*>& histos2, int nROF2, int nTF2, bool perROF);

//_________________________________________________________________________________________________
void CompareClustersMult(std::string clusterFileName1, std::string clusterFileName2, bool perROF = true)
{
  /// compare the multiplicity of all clusters

  std::vector<TH1*> multPlots1{};
  CreateMultPlots(multPlots1, "File1");
  auto [nROF1, nTF1] = FillMultPlots(clusterFileName1, multPlots1);

  std::vector<TH1*> multPlots2{};
  CreateMultPlots(multPlots2, "File2");
  auto [nROF2, nTF2] = FillMultPlots(clusterFileName2, multPlots2);

  gStyle->SetOptStat(1);

  DrawMultPlots(multPlots1, nROF1, nTF1, multPlots2, nROF2, nTF2, perROF);
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
void CreateMultPlots(std::vector<TH1*>& histos, std::string extension)
{
  /// create multiplicity histograms

  static const int nDE[10] = {4, 4, 4, 4, 18, 18, 26, 26, 26, 26};

  histos.emplace_back(new TH1D(fmt::format("nClustersPerCh{}", extension).c_str(),
                               "number of clusters per chamber;chamber", 10, 0.5, 10.5));

  for (int i = 0; i < 10; ++i) {
    histos.emplace_back(new TH1D(fmt::format("nClustersPerDEinCh{}{}", i + 1, extension).c_str(),
                                 fmt::format("number of clusters per DE in chamber {};DE", i + 1).c_str(),
                                 nDE[i], 100 * (i + 1) - 0.5, 100 * (i + 1) + nDE[i] - 0.5));
  }

  histos.emplace_back(new TH1D(fmt::format("nClusters{}", extension).c_str(),
                               "number of clusters;#clusters", 10001, -0.5, 10000.5));

  for (int i = 0; i < 10; ++i) {
    histos.emplace_back(new TH1D(fmt::format("nClustersinCh{}{}", i + 1, extension).c_str(),
                                 fmt::format("number of clusters in chamber {};#clusters", i + 1).c_str(),
                                 2001, -0.5, 2000.5));
  }

  for (auto h : histos) {
    h->SetDirectory(0);
  }
}

//_________________________________________________________________________________________________
std::pair<int, int> FillMultPlots(std::string clusterFileName, std::vector<TH1*>& histos)
{
  /// load clusters from clusterFileName and fill multiplicity histograms

  auto [fClusters, clusterReader] = LoadData(clusterFileName.c_str(), "o2sim");
  TTreeReaderValue<std::vector<Cluster>> clusters = {*clusterReader, "clusters"};
  TTreeReaderValue<std::vector<ROFRecord>> rofs = {*clusterReader, "clusterrofs"};

  int nTF = clusterReader->GetEntries(false);
  int iTF = 0;
  int nROF = 0;

  while (clusterReader->Next()) {
    if (++iTF % 1000 == 0) {
      std::cout << "\rprocessing TF " << iTF << " / " << nTF << "..." << std::flush;
    }
    nROF += rofs->size();
    FillMultPlots(*rofs, *clusters, histos);
  }
  std::cout << "\r\033[Kprocessing " << nTF << " TF completed" << std::endl;

  fClusters->Close();

  return std::make_pair(nROF, nTF);
}

//_________________________________________________________________________________________________
void FillMultPlots(const std::vector<ROFRecord>& rofs, const std::vector<Cluster>& clusters, std::vector<TH1*>& histos)
{
  /// fill multiplicity histograms

  for (const auto& cluster : clusters) {
    int chId = cluster.getChamberId() + 1;
    histos[0]->Fill(chId);
    histos[chId]->Fill(cluster.getDEId());
  }

  for (const auto& rof : rofs) {
    histos[11]->Fill(rof.getNEntries());
    std::array<int, 10> nClusters{};
    for (int i = rof.getFirstIdx(); i <= rof.getLastIdx(); ++i) {
      nClusters[clusters[i].getChamberId()]++;
    }
    for (int i = 0; i < 10; ++i) {
      histos[12 + i]->Fill(nClusters[i]);
    }
  }
}

//_________________________________________________________________________________________________
void DrawMultPlots(std::vector<TH1*>& histos1, int nROF1, int nTF1, std::vector<TH1*>& histos2, int nROF2, int nTF2, bool perROF)
{
  /// draw multiplicity histograms

  int norm1 = perROF ? nROF1 : nTF1;
  int norm2 = perROF ? nROF2 : nTF2;

  auto drawHistos = [&norm1, &norm2](TH1* h1, TH1* h2) {
    h1->SetStats(false);
    h1->Scale(1. / norm1);
    h1->SetLineColor(4);
    h1->Draw();
    h2->Scale(1. / norm2);
    h2->SetLineColor(2);
    h2->Draw("same");
    double h1minmax[2] = {0., 0.};
    h1->GetMinimumAndMaximum(h1minmax[0], h1minmax[1]);
    double h2minmax[2] = {0., 0.};
    h2->GetMinimumAndMaximum(h2minmax[0], h2minmax[1]);
    h1minmax[0] = std::min(h1minmax[0], h2minmax[0]);
    h1minmax[1] = std::max(h1minmax[1], h2minmax[1]);
    double h1range = h1minmax[1] - h1minmax[0];
    h1->SetMinimum(h1minmax[0] - 0.1 * h1range);
    h1->SetMaximum(h1minmax[1] + 0.1 * h1range);
  };

  auto drawRatio = [](TH1* h1, TH1* h2, std::string title) {
    TH1D* hRat = static_cast<TH1D*>(h2->Clone());
    hRat->SetTitle(title.c_str());
    hRat->Divide(h1);
    hRat->SetStats(false);
    hRat->SetLineColor(2);
    hRat->Draw();
  };

  TLegend* lHist = new TLegend(0.15, 0.7, 0.75, 0.85);
  lHist->SetFillStyle(0);
  lHist->SetBorderSize(0);
  lHist->AddEntry(histos1[0], fmt::format("file 1: {:.2f} ROFs/TF", static_cast<float>(nROF1) / nTF1).c_str(), "l");
  lHist->AddEntry(histos2[0], fmt::format("file 2: {:.2f} ROFs/TF", static_cast<float>(nROF2) / nTF2).c_str(), "l");

  TCanvas* cCh = new TCanvas("cNClustersPerCh", "cNClustersPerCh", 10, 10, 600, 300);
  cCh->Divide(2, 1);
  cCh->cd(1);
  drawHistos(histos1[0], histos2[0]);
  lHist->Clone()->Draw("same");
  cCh->cd(2);
  drawRatio(histos1[0], histos2[0], "h2 / h1");

  TCanvas* cDE = new TCanvas("cNClustersPerDE", "cNClustersPerDE", 10, 10, 1500, 600);
  cDE->Divide(5, 2);
  TCanvas* cDERat = new TCanvas("cNClustersPerDERatio", "cNClustersPerDERatio", 10, 10, 1500, 600);
  cDERat->Divide(5, 2);
  for (int i = 1; i <= 10; ++i) {
    cDE->cd(i);
    drawHistos(histos1[i], histos2[i]);
    cDERat->cd(i);
    drawRatio(histos1[i], histos2[i], fmt::format("h2 / h1 in chamber {}", i));
  }
  cDE->cd(1);
  lHist->Clone()->Draw("same");

  TCanvas* cN = new TCanvas("cNClusters", "cNClusters", 10, 10, 600, 300);
  cN->Divide(2, 1);
  cN->cd(1);
  drawHistos(histos1[11], histos2[11]);
  histos1[11]->SetMinimum();
  gPad->SetLogy();
  lHist->Clone()->Draw("same");
  cN->cd(2);
  drawRatio(histos1[11], histos2[11], "h2 / h1");

  TCanvas* cNCh = new TCanvas("cNClustersPerChDist", "cNClustersPerChDist", 10, 10, 1500, 600);
  cNCh->Divide(5, 2);
  TCanvas* cNChRat = new TCanvas("cNClustersPerChRatio", "cNClustersPerChRatio", 10, 10, 1500, 600);
  cNChRat->Divide(5, 2);
  for (int i = 1; i <= 10; ++i) {
    cNCh->cd(i);
    drawHistos(histos1[11 + i], histos2[11 + i]);
    histos1[11 + i]->SetMinimum();
    gPad->SetLogy();
    cNChRat->cd(i);
    drawRatio(histos1[11 + i], histos2[11 + i], fmt::format("h2 / h1 in chamber {}", i));
  }
  cNCh->cd(1);
  lHist->Clone()->Draw("same");
}
