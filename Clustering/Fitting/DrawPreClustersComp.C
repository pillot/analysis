#include <string>

#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TText.h>

#include <fmt/format.h>

#include "Framework/Logger.h"

template <class T>
T* GetClone(TFile& f, std::string canvasName, std::string objectName, std::string extension);

//_________________________________________________________________________________________________
void DrawPreClustersComp(std::string file1 = "displays.root", std::string file2 = "displays.root",
                         bool norm = true)
{
  /// compare precluster and associated digits informations

  gStyle->SetOptStat(1);
  TH1::AddDirectory(false);

  TFile f[] = {{file1.c_str(), "read"}, {file2.c_str(), "read"}};
  int color[] = {2, 4};
  std::string hDrawOpt[] = {"", "sames"};

  auto l = new TPaveText(0.55, 0.75, 0.75, 0.85, "NDC");
  l->SetFillStyle(0);
  l->SetBorderSize(0);
  for (int i = 0; i < 2; ++i) {
    auto t = l->AddText(fmt::format("file {}", i + 1).c_str());
    t->SetTextColor(color[i]);
  }

  auto hDraw = [norm, color, hDrawOpt](auto h, int i) {
    if (norm) {
      h->Scale(1. / h->GetEntries());
    }
    h->SetStats();
    h->SetLineColor(color[i]);
    h->Draw(hDrawOpt[i].c_str());
  };

  std::string sStation[] = {"St1", "St2", "St345"};
  for (auto sSt : sStation) {

    auto cName = fmt::format("preClusterInfo{}", sSt);
    auto c = new TCanvas(fmt::format("{}Comp", cName).c_str(),
                         fmt::format("precluster characteristics {}", sSt).c_str(), 10, 10, 1200, 600);
    c->Divide(3, 2);
    std::string hNameBase[] = {"hCharge", "hChargeB", "hChargeNB", "hChargeAsymm2", "hDimX", "hDimY"};
    for (int i = 0; i < 6; ++i) {
      auto hName = fmt::format("{}{}", hNameBase[i], sSt);
      for (int j = 0; j < 2; ++j) {
        auto h = GetClone<TH1>(f[j], cName, hName, fmt::format("_{}", j + 1));
        c->cd(i + 1);
        gPad->SetLogy();
        hDraw(h, j);
      }
    }
    c->cd(1);
    l->Clone()->Draw("same");

    cName = fmt::format("digitChargeInfo{}", sSt);
    c = new TCanvas(fmt::format("{}Comp", cName).c_str(),
                    fmt::format("digit characteristics {}", sSt).c_str(), 10, 10, 1200, 600);
    c->Divide(3, 2);
    std::string hNameBase2[] = {"ADC", "ADCB", "ADCNB", "Samples", "SamplesB", "SamplesNB"};
    for (int i = 0; i < 6; ++i) {
      auto hName = fmt::format("{}{}", hNameBase2[i], sSt);
      for (int j = 0; j < 2; ++j) {
        auto h = GetClone<TH1>(f[j], cName, hName, fmt::format("_{}", j + 1));
        c->cd(i + 1);
        gPad->SetLogy();
        hDraw(h, j);
      }
    }
    c->cd(1);
    l->Clone()->Draw("same");
  }

  f[0].Close();
  f[1].Close();
}

//_________________________________________________________________________________________________
template <class T>
T* GetClone(TFile& f, std::string canvasName, std::string objectName, std::string extension)
{
  /// return a clone of the object from the canvas

  auto c = f.FindObjectAny(canvasName.c_str());
  if (c == nullptr) {
    LOGP(error, "canvas {} not found", canvasName);
    exit(-1);
  }

  auto o = c->FindObject(objectName.c_str());
  if (o == nullptr) {
    LOGP(error, "object {} in canvas {} not found", objectName, canvasName);
    exit(-1);
  }

  return static_cast<T*>(o->Clone(fmt::format("{}{}", objectName, extension).c_str()));
}
