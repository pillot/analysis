#include "Riostream.h"
#include "TH1.h"
#include <vector>
#include "TMath.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TFitResult.h"
#include "TStyle.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TString.h"
#include "TVirtualPad.h"
#include "TPaveStats.h"

//
// Macro to fit the pT cross-section
//

//Input file should have the format :
/*
yrange               2.5-3                     3-3.5              ...
//pTrange  dsigmadpT   stat  sys uncor  dsigmadpT   stat  sys uncor ...
0-1      117.697  2.759   9.044 ...
1-2      242.073  3.866  18.778 ...
2-3      207.527  3.313  15.971 ...
3-4      138.534  2.675  10.836 ...
4-5       81.585  1.919   6.353 ...
5-6       44.827  1.130   3.484 ...
6-8       17.842  0.442   1.453 ...
8-15       2.419  0.071   0.187 ...
 ...
*/

//------------------------------------------------------------------
void FitCrossSectionVspTy(TString input="pt.txt", TString input2="", Int_t errorType = 2)
{
  // Results in the 2 files are summed (in case cross sections are separated between prompt and non prompt)
  // errorType: stat = 0, sys = 1, all = other int
  
  //gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  
  Bool_t sum = !(input2.IsNull());
  
  ifstream in(input.Data());
  ifstream in2(input2.Data());
  
  Double_t a, b;
  std::vector<double> pT;
  std::vector<double> y;
  std::vector<double> csVspT;
  std::vector<double> csVspTerr;
  std::vector<double> csVsy;
  std::vector<double> csVsyerr;
  std::vector<double> csVspTy[100];
  std::vector<double> csVspTyerr[100];
  
  std::cout << std::endl;
  while (!in.eof()) {
    
    TString line;
    line.ReadLine(in,kTRUE);
    
    TString line2;
    if (sum) {
      if (in2.eof()) {
        printf("mismatch between the 2 files\n");
        return;
      }
      line2.ReadLine(in2,kTRUE);
    }
    
    if (line.IsNull() || line.BeginsWith("//")) continue;
    if (sum && (line2.IsNull() || line2.BeginsWith("//"))) continue;
    
    TObjArray *inputs = line.Tokenize(" ");
    TObjArray *inputs2 = sum ? line2.Tokenize(" ") : 0x0;
    
    if (y.size() == 0) {
      
      if (sum && line != line2) {
        printf("mismatch between the 2 files\n");
        return;
      }
      
      if (inputs->GetEntriesFast() < 2) {
        printf("invalid format\n");
        return;
      }
      
      // decode the rapidity bins
      for (Int_t i = 1; i < inputs->GetEntriesFast(); ++i) {
        sscanf(inputs->At(i)->GetName(),"%le-%le",&a,&b);
        y.push_back(a);
        csVsy.push_back(0.);
        csVsyerr.push_back(0.);
      }
      y.push_back(b);
      
    } else {
      
      if (sum && inputs->GetEntriesFast() != inputs2->GetEntriesFast()) {
        printf("mismatch between the 2 files\n");
        return;
      }
      
      if ((ULong_t)inputs->GetEntriesFast() < (1+3*(y.size()-1))) {
        printf("invalid format\n");
        return;
      }
      
      // decode pT bins
      sscanf(inputs->At(0)->GetName(),"%le-%le",&a,&b);
      pT.push_back(a);
      csVspT.push_back(0.);
      csVspTerr.push_back(0.);
      
      for (ULong_t i = 0; i < (y.size()-1); ++i) {
        
        // decode values
        Double_t value = static_cast<TObjString*>(inputs->At(3*i+1))->String().Atof();
        Double_t stat = static_cast<TObjString*>(inputs->At(3*i+2))->String().Atof();
        Double_t sys = static_cast<TObjString*>(inputs->At(3*i+3))->String().Atof();
        if (errorType == 1) stat=0;
        else if (errorType == 0) sys=0;
        Double_t error = TMath::Sqrt(stat*stat + sys*sys);
        
        if (sum) {
          Double_t value2 = static_cast<TObjString*>(inputs2->At(3*i+1))->String().Atof();
          Double_t stat2 = static_cast<TObjString*>(inputs2->At(3*i+2))->String().Atof();
          Double_t sys2 = static_cast<TObjString*>(inputs2->At(3*i+3))->String().Atof();
          if (errorType == 1) stat2=0;
          else if (errorType == 0) sys2=0;
          Double_t error2 = stat2*stat2 + sys2*sys2;
          value += value2;
          error = TMath::Sqrt(error*error + error2);
        }
        
        // cross section vs pT/y in the current y bin
        csVspTy[pT.size()-1].push_back(value);
        csVspTyerr[pT.size()-1].push_back(error);
        
        // cross section vs pT integrated of over y
        if (y[i] > 2.49 && y[i+1] < 4.01) {
          csVspT[pT.size()-1] += value;
          csVspTerr[pT.size()-1] += error*error;
        }
        
        // cross section vs y integrated of over pT
        if (a > -0.1 && b < 12.1) {
          csVsy[i] += value;
          csVsyerr[i] += error*error;
        }
        
      }
      
    }
    
  }
  
  pT.push_back(b);
  
  for (ULong_t i = 0; i < (pT.size()-1); ++i) {
    csVspTerr[i] = TMath::Sqrt(csVspTerr[i]);
    //printf("%4.1lf %4.1lf\n", csVspT[i], csVspTerr[i]);
  }
  
  for (ULong_t i = 0; i < (y.size()-1); ++i) {
    csVsyerr[i] = TMath::Sqrt(csVsyerr[i]);
    //printf("%4.1lf %4.1lf\n", csVsy[i], csVsyerr[i]);
  }
  /*
  for (ULong_t i = 0; i < (pT.size()-1); ++i) {
    for (ULong_t j = 0; j < (y.size()-1); ++j) {
      printf("%4.1lf %4.1lf ", csVspTy[i][j], csVspTyerr[i][j]);
    }
    printf("\n");
  }
  */
  // histo vs pT
  TH1* hpTInt = new TH1F("pT_(2.5-4)","ds/pTdy versus pT (2.5 < y < 4)",pT.size()-1,&pT[0]);
  for (ULong_t i = 0; i < (pT.size()-1); ++i) {
    hpTInt->SetBinContent(i+1,csVspT[i]);
    hpTInt->SetBinError(i+1,csVspTerr[i]);
  }
  TH1* hpT[100];
  for (ULong_t i = 0; i < (y.size()-1); ++i) {
    hpT[i] = new TH1F(Form("pT_(%03.1f-%03.1f)",y[i],y[i+1]),Form("ds/pTdy versus pT (%3.1f < y < %3.1f)",y[i],y[i+1]),pT.size()-1,&pT[0]);
    for (ULong_t j = 0; j < (pT.size()-1); ++j) {
      hpT[i]->SetBinContent(j+1,csVspTy[j][i]);
      hpT[i]->SetBinError(j+1,csVspTyerr[j][i]);
    }
  }
  
  // fit function vs pT
  Int_t fitpTFunc = 1; // switch between fit functions
  TF1* fpT = 0x0;
  if (fitpTFunc == 1) fpT = new TF1("fpT","[0] * x / TMath::Power([1] + TMath::Power(x,[2]), [3])",hpT[0]->GetXaxis()->GetXmin(),hpT[0]->GetXaxis()->GetXmax());
  else fpT = new TF1("fpT","[0] * x / TMath::Power(1. + TMath::Power(x/[1],[2]), [3])",hpT[0]->GetXaxis()->GetXmin(),hpT[0]->GetXaxis()->GetXmax());
  fpT->SetNpx(1000);
  Double_t pTParam[100][4];
  TF1* fpT2;
  Int_t color = 1;
  
  // Draw results vs pT
  TCanvas* cpT = new TCanvas("cpT","cpT",800,800);
  cpT->Divide(3,3);
  cpT->cd(1);
  gPad->SetLogy();
  hpTInt->Draw("histe");
  gPad->Update();
  TPaveStats *stats = (TPaveStats*)gPad->GetPrimitive("stats");
  stats->SetX1NDC(0.45);
  stats->SetY1NDC(0.7);
  stats->SetOptStat(0);
  if (fitpTFunc == 1) fpT->SetParameters(10000., 14., 2., 4.);
  else {
    fpT->SetParameters(1000., 5., 2., 5.);
//    fpT->FixParameter(2, 2.);
  }
  hpTInt->Fit(fpT,"MIESQ","same",hpTInt->GetXaxis()->GetXmin(),hpTInt->GetXaxis()->GetXmax());
  fpT->Clone()->Draw("same");
  pTParam[0][0] = fpT->GetParameter(0);
  pTParam[0][1] = fpT->GetParameter(1);
  pTParam[0][2] = fpT->GetParameter(2);
  pTParam[0][3] = fpT->GetParameter(3);
  fpT2 = static_cast<TF1*>(fpT->Clone());
  fpT2->SetParameter(0,pTParam[0][0]/fpT2->Integral(0.,12.));
  cpT->cd(9);
  gPad->SetLogy();
  fpT2->SetLineColor(color++);
  fpT2->Draw();
  for (ULong_t i = 0; i < (y.size()-1); ++i) {
    cpT->cd(i+2);
    gPad->SetLogy();
    hpT[i]->Draw("histe");
    gPad->Update();
    stats = (TPaveStats*)gPad->GetPrimitive("stats");
    stats->SetX1NDC(0.4);
    stats->SetY1NDC(0.7);
    stats->SetOptStat(0);
    if (fitpTFunc == 1) fpT->SetParameters(10000., 14., 2., 4.);
    else {
      fpT->SetParameters(1000., 5., 2., 5.);
//      fpT->FixParameter(2, 2.);
    }
    hpT[i]->Fit(fpT,"MIESQ","same",hpT[i]->GetXaxis()->GetXmin(),hpT[i]->GetXaxis()->GetXmax());
    fpT->Clone()->Draw("same");
    pTParam[i+1][0] = fpT->GetParameter(0);
    pTParam[i+1][1] = fpT->GetParameter(1);
    pTParam[i+1][2] = fpT->GetParameter(2);
    pTParam[i+1][3] = fpT->GetParameter(3);
    fpT2 = static_cast<TF1*>(fpT->Clone());
    fpT2->SetParameter(0,pTParam[i+1][0]/fpT2->Integral(0.,12.));
    cpT->cd(9);
    if (color == 5 || color == 10) color++;
    fpT2->SetLineColor(color++);
    fpT2->Draw("same");
  }
  
  // histo vs y
  TH1* hyInt = new TH1F("y_(0-12)","ds/pTdy versus y (0 < pT < 12)",y.size()-1,&y[0]);
  for (ULong_t i = 0; i < (y.size()-1); ++i) {
    hyInt->SetBinContent(i+1,csVsy[i]);
    hyInt->SetBinError(i+1,csVsyerr[i]);
  }
  TH1* hy[100];
  for (ULong_t i = 0; i < (pT.size()-1); ++i) {
    hy[i] = new TH1F(Form("y_(%04.1f-%04.1f)",pT[i],pT[i+1]),Form("ds/pTdy versus y (%4.1f < pT < %4.1f)",pT[i],pT[i+1]),y.size()-1,&y[0]);
    for (ULong_t j = 0; j < (y.size()-1); ++j) {
      hy[i]->SetBinContent(j+1,csVspTy[i][j]);
      hy[i]->SetBinError(j+1,csVspTyerr[i][j]);
    }
  }
  
  // fit function vs y
  Int_t fityFunc = 1; // switch between fit functions
  TF1* fy = 0x0;
  if (fityFunc == 1) fy = new TF1("fy","[0] * (1. + [1]*x*x)",hy[0]->GetXaxis()->GetXmin(),hy[0]->GetXaxis()->GetXmax());
  else fy = new TF1("fy","[0] * TMath::Exp(-0.5*(x-[1])*(x-[1])/[2]/[2])",hy[0]->GetXaxis()->GetXmin(),hy[0]->GetXaxis()->GetXmax());
  Double_t yParam[100][3];
  TF1* fy2;
  color = 1;
  
  // Draw results vs y
  TCanvas* cy = new TCanvas("cy","cy",800,800);
  cy->Divide(4,4);
  cy->cd(1);
  hyInt->Draw("histe");
  gPad->Update();
  stats = (TPaveStats*)gPad->GetPrimitive("stats");
  stats->SetX1NDC(0.4);
  stats->SetY1NDC(0.75);
  stats->SetOptStat(0);
  if (fityFunc == 1) {
    fy->SetParameters(1., -0.04);
    hyInt->Fit(fy,"MIESQ","same",hy[0]->GetXaxis()->GetXmin(),hy[0]->GetXaxis()->GetXmax());
  } else {
    fy->SetParameters(1000., 0., 2.);
    fy->FixParameter(1, 0.);
    hyInt->Fit(fy,"MIESQ","same",2.5,4.);
  }
  fy->Clone()->Draw("same");
  yParam[0][0] = fy->GetParameter(0);
  yParam[0][1] = fy->GetParameter(1);
  if (fityFunc == 2) yParam[0][2] = fy->GetParameter(2);
  fy2 = static_cast<TF1*>(fy->Clone());
  fy2->SetParameter(0,yParam[0][0]/fy2->Integral(2.5,4.));
  cy->cd(16);
  fy2->SetLineColor(color++);
  fy2->Draw();
  for (ULong_t i = 0; i < (pT.size()-1); ++i) {
    cy->cd(i+2);
    hy[i]->Draw("histe");
    gPad->Update();
    stats = (TPaveStats*)gPad->GetPrimitive("stats");
    stats->SetX1NDC(0.4);
    stats->SetY1NDC(0.75);
    stats->SetOptStat(0);
    if (fityFunc == 1) {
      fy->SetParameters(1., -0.04);
      hy[i]->Fit(fy,"MIESQ","same",hy[i]->GetXaxis()->GetXmin(),hy[i]->GetXaxis()->GetXmax());
    } else {
      fy->SetParameters(1000., 0., 2.);
      fy->FixParameter(1, 0.);
      hy[i]->Fit(fy,"MIESQ","same",2.5,4.);
    }
    fy->Clone()->Draw("same");
    yParam[i+1][0] = fy->GetParameter(0);
    yParam[i+1][1] = fy->GetParameter(1);
    if (fityFunc == 2) yParam[i+1][2] = fy->GetParameter(2);
    fy2 = static_cast<TF1*>(fy->Clone());
    fy2->SetParameter(0,yParam[i+1][0]/fy2->Integral(2.5,4.));
    cy->cd(16);
    if (color == 5 || color == 10) color++;
    fy2->SetLineColor(color++);
    fy2->Draw("same");
  }
  
  printf("\n");
  for (ULong_t i = 0; i < y.size(); ++i) {
    printf("{%f, %f, %f, %f}\n", pTParam[i][0], pTParam[i][1], pTParam[i][2], pTParam[i][3]);
  }
  printf("\n");
  for (ULong_t i = 0; i < pT.size(); ++i) {
    if (fityFunc == 1) printf("{%f, %f}\n", yParam[i][0], yParam[i][1]);
    else printf("{%f, %f, %f}\n", yParam[i][0], yParam[i][1], yParam[i][2]);
  }
  printf("\n");

}
