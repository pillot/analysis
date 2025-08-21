#ifndef RESOLUTIONUTILS_H_
#define RESOLUTIONUTILS_H_

#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>
#include <fmt/format.h>
#include <gsl/span>

#include <TF1.h>
#include <THnSparse.h>
#include <TFile.h>
#include "TTreeReaderArray.h"
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

#include "DataFormatsMCH/Digit.h"
#include "MCHMappingInterface/Segmentation.h"
#include "MCHSimulation/Response.h"

using o2::mch::Digit;

//_________________________________________________________________________________________________
double ADCFit(const Digit digit, std::vector<double> parameters)
{
  //return the charge fraction seen by digit on a cathode given the cluster position

  static const o2::mch::Response response[] = {{o2::mch::Station::Type1}, {o2::mch::Station::Type2345}};

  const auto& segmentation = o2::mch::mapping::segmentation(digit.getDetID());
  int iSt = (digit.getDetID() < 300) ? 0 : 1;

  auto padid = digit.getPadID();
  auto dx = segmentation.padSizeX(padid) * 0.5;
  auto dy = segmentation.padSizeY(padid) * 0.5;
  auto xPad = segmentation.padPositionX(padid) - parameters[0];
  auto yPad = segmentation.padPositionY(padid) - parameters[1];
  auto qPad = response[iSt].chargePadfraction(xPad - dx, xPad + dx, yPad - dy, yPad + dy);
  return qPad * (segmentation.isBendingPad(padid) ? parameters[4] : parameters[5]);

}
//_________________________________________________________________________________________________
THnSparseD* CreatePreClusterInfoMULTI(const char* extension = "")
{
  const Int_t nDim = 9;

  Int_t nbins[nDim] = {
    1000,   // pvalue
    1122,  // residuals
    5001,  // ADC_fit
    5001,  // ADC_mes
    5001,  // ADC_cluster
    301,   // nSamples
    280,  // Asymm
    3,  // Wire
    2   // Cathode
  };

  Double_t xmin[nDim] = {
      0.,        // pvalue
    -280.5,    // residuals
    -0.5,     // ADC_fit
    -0.5,     // ADC_mes
    -5,     // ADC_cluster
    -0.5,     // nSamples
    -0.7,  // Asymm
    -0.5,     // Wire
    -2  // NonBending
  };

  Double_t xmax[nDim] = {
    1.,       // pvalue
    280.5,     // residuals
    5000.5,  // ADC_fit
    5000.5,  // ADC_mes
    50005,  // ADC_cluster
    300.5,    // nSamples
    0.7,   // Asymm
    2.5,   // Wire
    2   // Bending
  };

  TString name = Form("MultiResolutionPreCluster%s", extension);
  TString title = "9D Sparse Histograms for Pre-Cluster ";

  THnSparseD* hSparse = new THnSparseD(name, title, nDim, nbins, xmin, xmax);

  const char* axisTitles[nDim] = {
    "pvalue", "residuals", "ADC_fit", "ADC_mes", "ADC_cluster",
    "nSamples", "Asymm", "Wire", "Cathode"
  };

  for (Int_t i = 0; i < nDim; ++i) {
    hSparse->GetAxis(i)->SetTitle(axisTitles[i]);
  }

  return hSparse;
}
//_________________________________________________________________________________________________
void FillResolutionInfo(const Digit digit, std::vector<double> parameters, THnSparseD* h)
{
  //fill THnSparse9D histogram of precluster 

  std::vector<double> pre_parameters(parameters.begin(), parameters.begin() + 6);

  const auto& segmentation = o2::mch::mapping::segmentation(digit.getDetID());
  auto padid = digit.getPadID();

  Double_t position = -1.;
  if( (std::abs(parameters[7]) < 0.015)){ //"top"
      position = 0.;
  }else if(std::abs(parameters[7]) > 0.075){ //"between"
      position = 2.;
  }else if((std::abs(parameters[7]) > 0.015) || (std::abs(parameters[7]) < 0.075)){ //"crossover"
      position = 1.;
  }

  Double_t ADC_fit       = ADCFit(digit, pre_parameters);
  Double_t ADC_mes       = digit.getADC();
  Double_t residuals     = (ADC_mes - ADC_fit);
  Double_t ADC_cluster   = parameters[6];
  Double_t nSamples      = digit.getNofSamples();
  Double_t Asymm         = parameters[7];
  Double_t Wire          = position;
  Double_t Bending       = (segmentation.isBendingPad(padid) ? 1. : -1.);
  Double_t pvalue        = parameters[8];

  Double_t values[9] = {
      pvalue, residuals, ADC_fit, ADC_mes, ADC_cluster,
      nSamples, Asymm, Wire, Bending
  };

  h->Fill(values);
 

}
//_________________________________________________________________________________________________
THnSparseD* CreatePreClusterInfoMULTIK3(const char* extension = "")
{
  const Int_t nDim = 8;

  Int_t nbins[nDim] = {
    1000,   // pvalue
    40,  // k3x
    40,  // k3y
    5001,  // ADC_cluster
    51  ,    // p 
    62  ,   // phi
    280,  // Asymm
    3  // Wire
  };

  Double_t xmin[nDim] = {
    0.,        // pvalue
    0.,  // k3x
    0.,  // k3y
    -5,  // ADC_cluster
    -0.5,    // p 
    -15.5,   // phi
    -0.7,  // Asymm
    -0.5     // Wire
  };

  Double_t xmax[nDim] = {
    1.,       // pvalue
    1.,  // k3x
    1.,  // k3y
    50005,  // ADC_cluster
    50.5,    // p 
    15.5,   // phi
    0.7,   // Asymm
    2.5   // Wire
  };

  TString name = Form("MultiK3PreCluster%s", extension);
  TString title = "10D Sparse Histograms for Pre-Cluster ";

  THnSparseD* hSparse = new THnSparseD(name, title, nDim, nbins, xmin, xmax);

  const char* axisTitles[nDim] = {
    "pvalue", "k3x", "k3y", "ADC_cluster", "p", "phi",
    "Asymm", "Wire"
  };

  for (Int_t i = 0; i < nDim; ++i) {
    hSparse->GetAxis(i)->SetTitle(axisTitles[i]);
  }

  return hSparse;
}
//_________________________________________________________________________________________________
void FillK3Info(std::vector<double> parameters, THnSparseD* h)
{
  //fill THnSparse10D histogram of precluster 
  Double_t position = -1.;
  if( (std::abs(parameters[7]) < 0.015)) { //"top"
      position = 0.;
  } else if(std::abs(parameters[7]) > 0.075) { //"between"
      position = 2.;
  } else if((std::abs(parameters[7]) > 0.015) || (std::abs(parameters[7]) < 0.075)) { //"crossover"
      position = 1.;
  }

  Double_t ADC_cluster   = parameters[6];
  Double_t Asymm         = parameters[7];
  Double_t Wire          = position;
  Double_t pvalue        = parameters[8];
  Double_t k3x           = parameters[2];
  Double_t k3y           = parameters[3];
  Double_t p             = parameters[9];
  Double_t phi           = parameters[10];

  Double_t values[8] = {
      pvalue, k3x, k3y, ADC_cluster, p, phi,
      Asymm, Wire, 
  };

  h->Fill(values);


}

//_________________________________________________________________________________________________
void Resolution(TList* &list, TH2D* hist2D, int statistics, bool auto_bin)
{
  //default digit range value : 20 - 10000 ADC
  int start = 20, end = 10000;
  std::vector<double> wavg_charge, err_charge_l, err_charge_r, sigma, err_sigma;

  TH1D* projX = hist2D->ProjectionX();
  int binStart = projX->FindBin(start);
  int binEnd = projX->FindBin(end);
  int bin_i = binStart;
  double wavg = 0.;

  std::vector<std::pair<int,int>> intervals;

  if(auto_bin){
    for(int bin_j = binStart; bin_j < binEnd; bin_j++) {

        int integral = projX->Integral(bin_i, bin_j);
        wavg +=  projX->GetBinCenter(bin_j) * projX->GetBinContent(bin_j);

        if(integral > statistics){
            wavg_charge.push_back(wavg / integral);
            err_charge_l.push_back((wavg / integral) - projX->GetBinLowEdge(bin_i));
            err_charge_r.push_back(projX->GetBinLowEdge(bin_j + 1) - (wavg / integral));    
            intervals.push_back(std::make_pair(bin_i, bin_j));
            bin_i = bin_j + 1;
            wavg = 0.; //reset wavg to 0
        }
    }
  }
  else
  {
    //hardcoded binning
    std::vector<std::pair<int,int>>  charge_bin;

    for (int i = start; i <= 500; i += 1) {
      charge_bin.push_back({i,i});
    }
    for (int i = charge_bin.back().second + 1; i <= 1000; i += 4) {
      charge_bin.push_back({i,i+3});
    }
    for (int i = charge_bin.back().second + 1; i <= 2500; i += 10) {
      charge_bin.push_back({i,i+9});
    }
    for (int i = charge_bin.back().second + 1; i <= 5000; i += 16) {
      charge_bin.push_back({i,i+15});
    }
    for (int i = charge_bin.back().second + 1; i <= end; i += 26) {
      charge_bin.push_back({i,i+25});
    }


    for(auto interval : charge_bin) {
      int bin_i = projX->FindBin(interval.first);
      int bin_j = projX->FindBin(interval.second);
      int integral = projX->Integral(bin_i, bin_j);

      if(integral > statistics) {
          double wavg = (interval.first + interval.second)/2.;
          wavg_charge.push_back(wavg);
          err_charge_l.push_back(wavg - interval.first + 0.5); //+0.5 to fill the empty space between the bin edges (left)
          err_charge_r.push_back(interval.second - wavg + 0.5); //+0.5 to fill the empty space between the bin edges (right)  
          intervals.push_back(std::make_pair(bin_k, bin_k_1));
        }
    }
  }
  //auto rebin after first gaussian fit
  int index = 0;
  for(auto I : intervals) {

    TH1D* projY = hist2D->ProjectionY(Form("projY_%s_%d_%d", hist2D->GetName(), I.first, I.second), I.first, I.second);

    //to guide the gaussian fit 
    double initial_sigma = sqrt(wavg_charge[index]);
    double min = -1.8 * initial_sigma;
    double max = 1.8 * initial_sigma;

    //method of extraction
    TF1* fit = new TF1("fit", "gaus", min, max);
    fit->SetParameter(0, projY->GetMaximum());
    fit->SetParameter(1, 0.0);
    fit->SetParameter(2, 0.5 * initial_sigma);
    projY->Fit(fit, "RQ");

    double sigma1 = fit->GetParameter(2);
    double mean1  = fit->GetParameter(1);

    //---------- SECOND FIT ----------
    double range2 = 1.5 * sigma1;
    fit->SetRange(mean1 - range2, mean1 + range2);
    fit->SetParameters(fit->GetParameter(0), mean1, sigma1);
    projY->Fit(fit, "RQ");

    double sigma2 = fit->GetParameter(2);
    double mean2  = fit->GetParameter(1);

    //---------- THIRD FIT ----------
    double range3 = 1.5 * sigma2;
    fit->SetRange(mean2 - range3, mean2 + range3);
    fit->SetParameters(fit->GetParameter(0), mean2, sigma2);
    projY->Fit(fit, "RQ");
    
    double par[3];
    fit->GetParameters(par);
    sigma.push_back(fit->GetParameter(2));
    err_sigma.push_back(fit->GetParError(2));
    list->Add(projY);
    delete fit;
    index++;
  }
  delete projX; 
}

#endif