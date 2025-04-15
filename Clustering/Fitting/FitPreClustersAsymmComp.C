#include <string>

#include <TAxis.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TText.h>

#include <fmt/format.h>

#include "Framework/Logger.h"

template <class T>
T* GetClone(TFile& f, std::string canvasName, std::string objectName, std::string extension);

//_________________________________________________________________________________________________
void FitPreClustersAsymmComp(std::string file1 = "asymm.root", std::string file2 = "asymm.root")
{
  /// compare fitted cluster charge asymmetry

  gStyle->SetOptStat(1);

  TFile f[] = {{file1.c_str(), "read"}, {file2.c_str(), "read"}};
  int color[] = {2, 4};

  auto l = new TPaveText(0.15, 0.75, 0.35, 0.85, "NDC");
  l->SetFillStyle(0);
  l->SetBorderSize(0);
  for (int i = 0; i < 2; ++i) {
    auto t = l->AddText(fmt::format("file {}", i + 1).c_str());
    t->SetTextColor(color[i]);
  }

  auto gDraw = [color](auto g, int i, std::string opt, double range[2]) {
    g->GetXaxis()->SetRangeUser(0., 20000.);
    g->GetYaxis()->SetRangeUser(range[0], range[1]);
    g->SetMarkerColor(color[i]);
    g->SetLineColor(color[i]);
    g->Draw(opt.c_str());
  };

  std::string sStation[] = {"St1", "St2", "St345"};
  std::string gNameExt[] = {"", "Close", "Far"};
  for (auto sSt : sStation) {

    TCanvas* c = new TCanvas(fmt::format("cAsymm{}Comp", sSt).c_str(),
                             fmt::format("asymmetry vs charge {}", sSt).c_str(), 10, 10, 1200, 600);
    c->Divide(3, 2);
    std::string cNames[] = {"cAsymmMean", "cAsymmRMS"};
    std::string gNameBase[] = {"gAsymmMeanvsCharge", "gAsymmRMSvsCharge"};
    double range[2][2] = {{-0.1, 0.1}, {0., 0.4}};
    for (int i = 0; i < 2; ++i) {
      for (int j = 0; j < 3; ++j) {
        auto gName = fmt::format("{}{}{}", gNameBase[i], sSt, gNameExt[j]);
        std::string opt = "ap";
        for (int k = 0; k < 2; ++k) {
          auto g = GetClone<TGraphAsymmErrors>(f[k], cNames[i], gName, fmt::format("_{}", k + 1));
          if (g != nullptr) {
            c->cd(3 * i + j + 1);
            gDraw(g, k, opt, range[i]);
            opt = "p";
          }
        }
      }
    }
    c->cd(1);
    l->Clone()->Draw("same");

    c = new TCanvas(fmt::format("cAsymm{}GComp", sSt).c_str(),
                    fmt::format("gaus asymmetry vs charge {}", sSt).c_str(), 10, 10, 1200, 600);
    c->Divide(3, 2);
    std::string cNamesG[] = {"cAsymmMeanG", "cAsymmSigmaG"};
    std::string gNameBaseG[] = {"gAsymmMeanGvsCharge", "gAsymmSigmaGvsCharge"};
    double rangeG[][2] = {{-0.1, 0.1}, {0., 0.5}};
    for (int i = 0; i < 2; ++i) {
      for (int j = 0; j < 3; ++j) {
        auto gName = fmt::format("{}{}{}", gNameBaseG[i], sSt, gNameExt[j]);
        std::string opt = "ap";
        for (int k = 0; k < 2; ++k) {
          auto g = GetClone<TGraphAsymmErrors>(f[k], cNamesG[i], gName, fmt::format("_{}", k + 1));
          if (g != nullptr) {
            c->cd(3 * i + j + 1);
            gDraw(g, k, opt, rangeG[i]);
            opt = "p";
          }
        }
      }
    }
    c->cd(1);
    l->Clone()->Draw("same");

    c = new TCanvas(fmt::format("cAsymm{}DGComp", sSt).c_str(),
                    fmt::format("double gaus asymmetry vs charge {}", sSt).c_str(), 10, 10, 1200, 900);
    c->Divide(3, 3);
    std::string cNamesDG[] = {"cAsymmMeanDG", "cAsymmDeltaDG", "cAsymmSigmaDG"};
    std::string gNameBaseDG[] = {"gAsymmMeanDGvsCharge", "gAsymmDeltaDGvsCharge", "gAsymmSigmaDGvsCharge"};
    double rangeDG[][2] = {{-0.1, 0.1}, {0., 0.4}, {0., 0.4}};
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        auto gName = fmt::format("{}{}{}", gNameBaseDG[i], sSt, gNameExt[j]);
        std::string opt = "ap";
        for (int k = 0; k < 2; ++k) {
          auto g = GetClone<TGraphAsymmErrors>(f[k], cNamesDG[i], gName, fmt::format("_{}", k + 1));
          if (g != nullptr) {
            c->cd(3 * i + j + 1);
            gDraw(g, k, opt, rangeDG[i]);
            opt = "p";
          }
        }
      }
    }
    c->cd(1);
    l->Clone()->Draw("same");

    c = new TCanvas(fmt::format("cAsymm{}TGComp", sSt).c_str(),
                    fmt::format("triple gaus asymmetry vs charge {}", sSt).c_str(), 10, 10, 1200, 900);
    c->Divide(3, 4);
    std::string cNamesTG[] = {"cAsymmFracTG", "cAsymmMeanTG", "cAsymmDeltaTG", "cAsymmSigmaTG"};
    std::string gNameBaseTG[] = {"gAsymmFracTGvsCharge", "gAsymmMeanTGvsCharge", "gAsymmDeltaTGvsCharge", "gAsymmSigmaTGvsCharge"};
    double rangeTG[][2] = {{-0.5, 1.5}, {-0.1, 0.1}, {0., 0.4}, {0., 0.4}};
    for (int i = 0; i < 4; ++i) {
      for (int j = 0; j < 3; ++j) {
        auto gName = fmt::format("{}{}{}", gNameBaseTG[i], sSt, gNameExt[j]);
        std::string opt = "ap";
        for (int k = 0; k < 2; ++k) {
          auto g = GetClone<TGraphAsymmErrors>(f[k], cNamesTG[i], gName, fmt::format("_{}", k + 1));
          if (g != nullptr) {
            c->cd(3 * i + j + 1);
            gDraw(g, k, opt, rangeTG[i]);
            opt = "p";
          }
        }
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
    return nullptr;
  }

  auto o = c->FindObject(objectName.c_str());
  if (o == nullptr) {
    LOGP(error, "object {} in canvas {} not found", objectName, canvasName);
    return nullptr;
  }

  return static_cast<T*>(o->Clone(fmt::format("{}{}", objectName, extension).c_str()));
}
