#ifndef DIGITUTILS_H_
#define DIGITUTILS_H_

#include <cmath>
#include <utility>
#include <vector>

#include <gsl/span>

#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLegend.h>

#include "DataFormatsMCH/Digit.h"
#include "MCHMappingInterface/Segmentation.h"

using o2::mch::Digit;

/*
 * This file contains utility functions to produce digit control plots.
 * To be used, they require the MCH mapping to be loaded:
 * gSystem->Load("libO2MCHMappingImpl4")
 */

//_________________________________________________________________________________________________
double adcToCharge(uint32_t adc)
{
  /// function to reinterpret digit ADC as calibrated charge in run2

  float charge(0.);
  std::memcpy(&charge, &adc, sizeof(adc));

  return static_cast<double>(charge);
}

//_________________________________________________________________________________________________
void CreateDigitTimeInfo(std::vector<TH1*>& histos, const char* extension = "")
{
  /// create digit time histograms

  histos.emplace_back(new TH1F(Form("time%s", extension), "#Deltat vs track time;#Deltat (BC)", 101, -50.5, 50.5));

  for (auto h : histos) {
    h->SetDirectory(0);
  }
}

//_________________________________________________________________________________________________
void FillDigitTimeInfo(const Digit& digit, int trackTime, gsl::span<TH1*> histos)
{
  /// fill digit time histograms

  histos[0]->Fill(digit.getTime() + 1.5 - trackTime);
}

//_________________________________________________________________________________________________
TCanvas* DrawDigitTimeInfo(gsl::span<TH1*> histos, const char* extension = "")
{
  /// draw digit time histograms

  TCanvas* c = new TCanvas(Form("digitTimeInfo%s", extension),
                           Form("digit time characteristics %s", extension), 10, 10, 400, 400);
  gPad->SetLogy();
  histos[0]->Draw();

  return c;
}

//_________________________________________________________________________________________________
void CreateDigitChargeInfo(std::vector<TH1*>& histos, const char* extension = "")
{
  /// create digit charge histograms

  histos.emplace_back(new TH1F(Form("ADC%s", extension), "ADC;ADC", 10001, -0.5, 10000.5));
  histos.emplace_back(new TH1F(Form("ADCB%s", extension), "ADC Bending;ADC", 10001, -0.5, 10000.5));
  histos.emplace_back(new TH1F(Form("ADCNB%s", extension), "ADC Non Bending;ADC", 10001, -0.5, 10000.5));
  histos.emplace_back(new TH1F(Form("Samples%s", extension), "N samples;N samples", 1024, -0.5, 1023.5));
  histos.emplace_back(new TH1F(Form("SamplesB%s", extension), "N samples Bending;N samples", 1024, -0.5, 1023.5));
  histos.emplace_back(new TH1F(Form("SamplesNB%s", extension), "N samples Non Bending;N samples", 1024, -0.5, 1023.5));
  histos.emplace_back(new TH2F(Form("ADCvsSampleB%s", extension), "ADC vs N samples Bending;N samples;ADC", 1024, -0.5, 1023.5, 10001, -0.5, 100009.5));
  histos.emplace_back(new TH2F(Form("ADCvsSampleNB%s", extension), "ADC vs N samples Non Bending;N samples;ADC", 1024, -0.5, 1023.5, 10001, -0.5, 100009.5));

  for (auto h : histos) {
    h->SetDirectory(0);
  }
}

//_________________________________________________________________________________________________
void FillDigitChargeInfo(const Digit& digit, gsl::span<TH1*> histos, bool run2 = false)
{
  /// fill digit charge histograms

  const auto& segmentation = o2::mch::mapping::segmentation(digit.getDetID());

  auto charge = run2 ? adcToCharge(digit.getADC()) : digit.getADC();

  histos[0]->Fill(charge);
  histos[3]->Fill(digit.getNofSamples());
  if (segmentation.isBendingPad(digit.getPadID())) {
    histos[1]->Fill(charge);
    histos[4]->Fill(digit.getNofSamples());
    histos[6]->Fill(digit.getNofSamples(), charge);
  } else {
    histos[2]->Fill(charge);
    histos[5]->Fill(digit.getNofSamples());
    histos[7]->Fill(digit.getNofSamples(), charge);
  }
}

//_________________________________________________________________________________________________
TCanvas* DrawDigitChargeInfo(gsl::span<TH1*> histos, const char* extension = "")
{
  /// draw digit charge histograms

  static int logy[3] = {1, 1, 0};
  static int logz[3] = {0, 0, 1};
  static const char* opt[3] = {"", "", "colz"};

  TCanvas* c = new TCanvas(Form("digitChargeInfo%s", extension),
                           Form("digit charge characteristics %s", extension), 10, 10, 800, 800);
  c->Divide(2, 2);
  c->cd(1);
  gPad->SetLogy();
  histos[0]->SetLineColor(1);
  histos[0]->Draw();
  histos[1]->SetStats(false);
  histos[1]->SetLineColor(2);
  histos[1]->Draw("sames");
  histos[2]->SetStats(false);
  histos[2]->SetLineColor(4);
  histos[2]->Draw("sames");
  c->cd(2);
  gPad->SetLogy();
  histos[3]->SetLineColor(1);
  histos[3]->Draw();
  histos[4]->SetStats(false);
  histos[4]->SetLineColor(2);
  histos[4]->Draw("sames");
  histos[5]->SetStats(false);
  histos[5]->SetLineColor(4);
  histos[5]->Draw("sames");
  c->cd(3);
  gPad->SetLogz();
  histos[6]->Draw("colz");
  c->cd(4);
  gPad->SetLogz();
  histos[7]->Draw("colz");

  TLegend* l = new TLegend(0.5, 0.7, 0.9, 0.9);
  l->SetFillStyle(0);
  l->SetBorderSize(0);
  l->AddEntry(histos[0], "all", "l");
  l->AddEntry(histos[1], "bending", "l");
  l->AddEntry(histos[2], "non bending", "l");
  c->cd(1);
  l->Draw("same");

  return c;
}

#endif // DIGITUTILS_H_