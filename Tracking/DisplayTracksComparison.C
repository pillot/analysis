#include <string>
#include <vector>

#include <gsl/span>

#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TParameter.h>

int LoadHistos(const char* fileName, std::vector<TH1*> histosAtVertex[2], TH1*& hmatchChi2, std::vector<TH1*>& chargeHistos);
void LoadHistosAtVertex(TFile* f, std::vector<TH1*>& histos, const char* extension);
void CompareHistosAtVertex(std::vector<TH1*> histos1, int nTF1, std::vector<TH1*> histos2, int nTF2, const char* extension);
void CompareMatchChi2(TH1* h1, int nTF1, TH1* h2, int nTF2);
void LoadChargeHistos(TFile* f, std::vector<TH1*>& histos, const char* extension);
void CompareChargeHistos(gsl::span<TH1*> histos1, int nTF1, gsl::span<TH1*> histos2, int nTF2, const char* extension);

//_________________________________________________________________________________________________
void DisplayTracksComparison(std::string inFileName1, std::string inFileName2 = "")
{
  /// Compare histograms between the 2 input files

  std::vector<TH1*> histosAtVertex1[2] = {{}, {}};
  TH1* hmatchChi21 = nullptr;
  std::vector<TH1*> chargeHistos1{};
  int nTF1 = LoadHistos(inFileName1.c_str(), histosAtVertex1, hmatchChi21, chargeHistos1);

  std::vector<TH1*> histosAtVertex2[2] = {{}, {}};
  TH1* hmatchChi22 = nullptr;
  std::vector<TH1*> chargeHistos2{};
  int nTF2 = LoadHistos(inFileName2.c_str(), histosAtVertex2, hmatchChi22, chargeHistos2);

  // display histograms
  CompareHistosAtVertex(histosAtVertex1[0], nTF1, histosAtVertex2[0], nTF2, "mch");
  CompareHistosAtVertex(histosAtVertex1[1], nTF1, histosAtVertex2[1], nTF2, "muon");
  CompareMatchChi2(hmatchChi21, nTF1, hmatchChi22, nTF2);
  CompareChargeHistos({&chargeHistos1[0], 2}, nTF1, {&chargeHistos2[0], 2}, nTF2, "AllDigits_mch");
  CompareChargeHistos({&chargeHistos1[3], 2}, nTF1, {&chargeHistos2[3], 2}, nTF2, "AllDigits_muon");
  CompareChargeHistos({&chargeHistos1[6], 2}, nTF1, {&chargeHistos2[6], 2}, nTF2, "DigitsAtClusterPos_mch");
  CompareChargeHistos({&chargeHistos1[9], 2}, nTF1, {&chargeHistos2[9], 2}, nTF2, "DigitsAtClusterPos_muon");
}

//_________________________________________________________________________________________________
int LoadHistos(const char* fileName, std::vector<TH1*> histosAtVertex[2], TH1*& hmatchChi2, std::vector<TH1*>& chargeHistos)
{
  /// load histograms from input file
  /// return the number of TF analyzed to produce them

  TFile* f = TFile::Open(fileName, "READ");
  if (!f || f->IsZombie()) {
    std::cout << "opening file " << fileName << " failed" << std::endl;
    exit(-1);
  }

  auto p = static_cast<TParameter<int>*>(f->FindObjectAny("nTF"));
  if (!p) {
    std::cout << "missing number of TF" << std::endl;
    exit(1);
  }

  LoadHistosAtVertex(f, histosAtVertex[0], "mch");
  LoadHistosAtVertex(f, histosAtVertex[1], "muon");
  hmatchChi2 = static_cast<TH1*>(f->FindObjectAny("matchChi2"));
  LoadChargeHistos(f, chargeHistos, "AllDig");
  LoadChargeHistos(f, chargeHistos, "AllDigMatch");
  LoadChargeHistos(f, chargeHistos, "");
  LoadChargeHistos(f, chargeHistos, "Match");

  return p->GetVal();
}

//_________________________________________________________________________________________________
void LoadHistosAtVertex(TFile* f, std::vector<TH1*>& histos, const char* extension)
{
  /// load single muon histograms at vertex

  histos.emplace_back(static_cast<TH1*>(f->FindObjectAny(Form("pT%s", extension))));
  histos.emplace_back(static_cast<TH1*>(f->FindObjectAny(Form("eta%s", extension))));
  histos.emplace_back(static_cast<TH1*>(f->FindObjectAny(Form("phi%s", extension))));
  histos.emplace_back(static_cast<TH1*>(f->FindObjectAny(Form("dca%s", extension))));
  histos.emplace_back(static_cast<TH1*>(f->FindObjectAny(Form("pDCA23%s", extension))));
  histos.emplace_back(static_cast<TH1*>(f->FindObjectAny(Form("pDCA310%s", extension))));
  histos.emplace_back(static_cast<TH1*>(f->FindObjectAny(Form("rAbs%s", extension))));
  histos.emplace_back(static_cast<TH1*>(f->FindObjectAny(Form("nClusters%s", extension))));
  histos.emplace_back(static_cast<TH1*>(f->FindObjectAny(Form("chi2%s", extension))));

  for (auto h : histos) {
    if (!h) {
      std::cout << "missing histo at vertex" << std::endl;
      exit(1);
    }
  }
}

//_________________________________________________________________________________________________
void CompareHistosAtVertex(std::vector<TH1*> histos1, int nTF1, std::vector<TH1*> histos2, int nTF2, const char* extension)
{
  /// draw histograms at vertex

  double nTFref = (nTF1 > nTF2) ? nTF1 : nTF2;

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
  TCanvas* cHist = new TCanvas(Form("histos_%s", extension), Form("histos_%s", extension),
                               10, 10, TMath::Max(nPadsx * 300, 1200), TMath::Max(nPadsy * 300, 900));
  cHist->Divide(nPadsx, nPadsy);
  for (int i = 0; i < (int)histos1.size(); ++i) {
    cHist->cd((i / nPadsx) * nPadsx + i % nPadsx + 1);
    gPad->SetLogy();
    histos1[i]->SetStats(false);
    histos1[i]->Scale(nTFref / nTF1);
    histos1[i]->SetLineColor(4);
    histos1[i]->Draw("hist");
    histos2[i]->Scale(nTFref / nTF2);
    histos2[i]->SetLineColor(2);
    histos2[i]->Draw("histsame");
  }

  // add a legend
  TLegend* lHist = new TLegend(0.5, 0.65, 0.9, 0.8);
  lHist->SetFillStyle(0);
  lHist->SetBorderSize(0);
  lHist->AddEntry(histos1[0], Form("file 1: %g tracks", histos1[0]->GetEntries()), "l");
  lHist->AddEntry(histos2[0], Form("file 2: %g tracks", histos2[0]->GetEntries()), "l");
  cHist->cd(1);
  lHist->Draw("same");

  // draw ratios
  TCanvas* cRat = new TCanvas(Form("ratios_%s", extension), Form("histos2 / histos1 - %s", extension),
                              10, 10, TMath::Max(nPadsx * 300, 1200), TMath::Max(nPadsy * 300, 900));
  cRat->Divide(nPadsx, nPadsy);
  for (int i = 0; i < (int)histos1.size(); ++i) {
    cRat->cd((i / nPadsx) * nPadsx + i % nPadsx + 1);
    TH1F* hRat = static_cast<TH1F*>(histos2[i]->Clone());
    hRat->Divide(histos1[i]);
    hRat->SetStats(false);
    hRat->SetLineColor(2);
    hRat->Draw();
  }

}

//_________________________________________________________________________________________________
void CompareMatchChi2(TH1* h1, int nTF1, TH1* h2, int nTF2)
{
  /// compare matching chi2

  double nTFref = (nTF1 > nTF2) ? nTF1 : nTF2;

  TCanvas* cHist = new TCanvas("cMatch", "cMatch", 10, 10, 400, 600);
  cHist->Divide(1, 2);
  cHist->cd(1);
  gPad->SetLogy();
  h1->SetStats(false);
  h1->Scale(nTFref / nTF1);
  h1->SetLineColor(4);
  h1->Draw("hist");
  h2->Scale(nTFref / nTF2);
  h2->SetLineColor(2);
  h2->Draw("histsame");
  cHist->cd(2);
  TH1F* hRat = static_cast<TH1F*>(h2->Clone());
  hRat->SetTitle("h2 / h1 ratio");
  hRat->Divide(h1);
  hRat->SetStats(false);
  hRat->SetLineColor(2);
  hRat->Draw();

  TLegend* lHist = new TLegend(0.5, 0.65, 0.9, 0.8);
  lHist->SetFillStyle(0);
  lHist->SetBorderSize(0);
  lHist->AddEntry(h1, Form("file 1: %d TF", nTF1), "l");
  lHist->AddEntry(h2, Form("file 2: %d TF", nTF2), "l");
  cHist->cd(1);
  lHist->Draw("same");
}

//_________________________________________________________________________________________________
void LoadChargeHistos(TFile* f, std::vector<TH1*>& histos, const char* extension)
{
  /// load track charge histograms

  histos.emplace_back(static_cast<TH1*>(f->FindObjectAny(Form("ADC%s", extension))));
  histos.emplace_back(static_cast<TH1*>(f->FindObjectAny(Form("Samples%s", extension))));
  histos.emplace_back(static_cast<TH1*>(f->FindObjectAny(Form("ADCvsSample%s", extension))));

  for (auto h : histos) {
    if (!h) {
      std::cout << "missing charge histo" << std::endl;
      exit(1);
    }
  }
}

//_________________________________________________________________________________________________
void CompareChargeHistos(gsl::span<TH1*> histos1, int nTF1, gsl::span<TH1*> histos2, int nTF2, const char* extension)
{
  /// draw track charge histograms

  double nTFref = (nTF1 > nTF2) ? nTF1 : nTF2;

  TCanvas* cHist = new TCanvas(Form("cCharge%s", extension), Form("cCharge%s", extension), 10, 10, 1200, 800);
  cHist->Divide(3, 2);
  for (int i = 0; i < 2; ++i) {
    cHist->cd(3 * i + 1);
    gPad->SetLogy();
    histos1[i]->SetStats(false);
    histos1[i]->Scale(nTFref / nTF1);
//    histos1[i]->Scale(1. / histos1[i]->GetEntries());
    histos1[i]->SetLineColor(4);
    histos1[i]->Draw("hist");
    histos2[i]->Scale(nTFref / nTF2);
//    histos2[i]->Scale(1. / histos2[i]->GetEntries());
    histos2[i]->SetLineColor(2);
    histos2[i]->Draw("histsame");
    cHist->cd(3 * i + 2);
    TH1F* hRat = static_cast<TH1F*>(histos2[i]->Clone());
    hRat->SetTitle("h2 / h1");
    hRat->Divide(histos1[i]);
    hRat->SetStats(false);
    hRat->SetLineColor(2);
    hRat->Draw("hist");
    cHist->cd(3 * i + 3);
    TH1F* hDiff = static_cast<TH1F*>(histos2[i]->Clone());
    hDiff->SetTitle("h2 - h1");
    hDiff->Add(histos1[i], -1.);
    hDiff->SetStats(false);
    hDiff->SetLineColor(2);
    hDiff->Draw("hist");
  }

  TLegend* lHist = new TLegend(0.5, 0.65, 0.9, 0.8);
  lHist->SetFillStyle(0);
  lHist->SetBorderSize(0);
  lHist->AddEntry(histos1[0], Form("file 1: %d TF", nTF1), "l");
  lHist->AddEntry(histos2[0], Form("file 2: %d TF", nTF2), "l");
  cHist->cd(1);
  lHist->Draw("same");
}
