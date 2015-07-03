#include "Riostream.h"
#include "TH1.h"
#include <vector>
#include "TMath.h"
#include "TF1.h"
#include "TVirtualFitter.h"
#include "TGraph.h"
#include "TMinuit.h"
#include "TCanvas.h"
#include "TFitResult.h"
#include "TStyle.h"
#include "TMatrixD.h"
#include "AliPWGFunc.h"

//
// Macro to get the <pt> from a pt distribution
//

//Input file should have the format :
//range  dsigmadpT   stat  sys uncor
//0-1      117.697  2.759   9.044
//1-2      242.073  3.866  18.778
//2-3      207.527  3.313  15.971
//3-4      138.534  2.675  10.836
//4-5       81.585  1.919   6.353
//5-6       44.827  1.130   3.484
//6-8       17.842  0.442   1.453
//8-15       2.419  0.071   0.187

// range to integrate the function
Double_t xmin = 0.;
Double_t xmax = 8.;

Double_t mJPsi = 3.096916;


//------------------------------------------------------------------
Double_t MeanPtVsx(Double_t *x, Double_t *p)
{
  /// <pT> with p[1..5] = UA1 function parameters and p[0] = index of the varying parameter to be replaced by x[0]
  
  TF1* fx = new TF1("fx", AliPWGFunc::StaticUA1Func, xmin, xmax, 6);
  
  fx->SetParameters(p);
  fx->SetParameter(0,mJPsi);
  fx->SetParameter(p[0],x[0]);
  
  return fx->Mean(xmin,xmax);
}

//------------------------------------------------------------------
Double_t MeanPt2Vsx(Double_t *x, Double_t *p)
{
  /// <pT2> with p[1..5] = UA1 function parameters and p[0] = index of the varying parameter to be replaced by x[0]
  
  TF1* fx2 = new TF1("fx2", AliPWGFunc::StaticUA1Func, xmin, xmax, 6);
  
  fx2->SetParameters(p);
  fx2->SetParameter(0,mJPsi);
  fx2->SetParameter(p[0],x[0]);
  
  return fx2->Moment(2,xmin,xmax);
}

//------------------------------------------------------------------
TH1* FitPtUA1(const char* input="pt.pA.igor.txt", Bool_t norm = kFALSE, Int_t errorType = 2)
{
  // norm = kTRUE if yields have been divided by pT
  // errorType: stat = 0, sys = 1, all = other int
  // delta: range for <pT> scan to <pT> ± delta*<pT>
  //
  // .include $ALICE_PHYSICS/include
  // .x FitPtUA1.C+(...)
  
  
  ifstream in(input);
  char line[1024];
  
  float a,b;
  float value,stat,sys;
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> yerr;
  
  std::cout << std::endl;
  while (in.getline(line,1023,'\n'))
  {
    
    if (strstr(line, "//")) continue;
    
    sscanf(line,"%e-%e %e %e %e",&a,&b,&value,&stat,&sys);
    
    if (norm) {
      Double_t s = 0.5*(a+b);
      value *= s;
      stat *= s;
      sys *= s;
    }
    
    if (errorType == 1) stat=0;
    else if (errorType == 0) sys=0;
    
    printf("%e-%e %e %e %e\n",a,b,value,stat,sys);
    
    x.push_back(a);
    y.push_back(value);
    
    Double_t error = TMath::Sqrt(stat*stat/value/value+sys*sys/value/value)*value;
    yerr.push_back(error);
    
    printf("%e-%e %g +- %g\n",a,b,value,error);
    
  }
  std::cout << std::endl;
  
  gStyle->SetOptFit(1);
  
  x.push_back(b);
  
  TH1* h = new TH1F("pt","pt",x.size()-1,&x[0]);
  
  for ( Int_t i = 1; i <= h->GetXaxis()->GetNbins(); ++i )
  {
    h->SetBinContent(i,y[i-1]);
    h->SetBinError(i,yerr[i-1]);
  }
  
  
  TCanvas* c1 = new TCanvas(input,input);
  c1->Divide(2,2);
  c1->cd(1);
  
  gPad->SetLogy();
  
  h->Draw("histe");
  
  TF1* f = new TF1 ("f", AliPWGFunc::StaticUA1Func, xmin, xmax, 6);
  f->FixParameter(0,mJPsi);
  f->SetParameter(1,2.);
  f->SetParameter(2,1.);
  f->SetParameter(3,5.);
  f->SetParameter(4,0.7);
  f->SetParameter(5,0.1);
  f->SetParLimits(1,0.00001,100000.);
  f->SetParLimits(2,0.00001,100000.);
  f->SetParLimits(3,0.00001,100000.);
  f->SetParLimits(4,0.00001,100000.);
  f->SetParLimits(5,0.00001,100000.);
  f->SetParNames("mass","p0star","pt0","n","T","norm");
  
  TFitResultPtr r = h->Fit(f,"MIS","",h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax());
  f->Draw("same");
  
  r->Print("V");
  
  std::cout << std::endl;
  for ( Int_t i = 1; i <= h->GetXaxis()->GetNbins(); ++i )
  {
    Double_t x1 = h->GetBinLowEdge(i);
    Double_t x2 = x1+h->GetBinWidth(i);
    
    std::cout << Form("Pt %4.2f - %4.2f Mean %5.3f",x1,x2,f->Mean(x1,x2)) << std::endl;
  }
  
  std::cout << std::endl;
  Double_t meanptval = f->Mean(xmin,xmax);
  std::cout << Form("<pT> = %e (%f,%f)",meanptval,xmin,xmax) << std::endl;
  Double_t meanpt2val = f->Moment(2,xmin,xmax);
  std::cout << Form("<pT2> = %e (%f,%f)",meanpt2val,xmin,xmax) << std::endl;
  std::cout << std::endl;
  
  // --- mean pT from fit function and mean pT error from error propagation ---
  TMatrixD cov(r->GetCovarianceMatrix());
  
  TMatrixD jac(1,6);
  jac(0,0) = 0.;
  TF1* fMeanpTVsx = new TF1("fMeanpTVsx",MeanPtVsx,0.,1.,6);
  fMeanpTVsx->SetParameters(f->GetParameters());
  for (Int_t ip = 1; ip < 6; ++ip) {
    fMeanpTVsx->SetParameter(0, ip);
    fMeanpTVsx->SetRange(f->GetParameter(ip)-f->GetParError(ip),f->GetParameter(ip)+f->GetParError(ip));
    jac(0,ip) = fMeanpTVsx->Derivative(f->GetParameter(ip),0x0,0.005);
    std::cout << Form("d<pT>/dp%d = %e ± %e",ip, jac(0,ip), fMeanpTVsx->DerivativeError()) << std::endl;
  }
  
  TMatrixD tmp(jac,TMatrixD::kMult,cov);
  TMatrixD meanpTErr2(tmp,TMatrixD::kMultTranspose,jac);
  
  // --- mean pT**2 from fit function and mean pT**2 error from error propagation ---
  TMatrixD jac2(1,6);
  jac2(0,0) = 0.;
  TF1* fMeanpT2Vsx = new TF1("fMeanpT2Vsx",MeanPt2Vsx,0.,1.,6);
  fMeanpT2Vsx->SetParameters(f->GetParameters());
  for (Int_t ip = 1; ip < 6; ++ip) {
    fMeanpT2Vsx->SetParameter(0, ip);
    fMeanpT2Vsx->SetRange(f->GetParameter(ip)-f->GetParError(ip),f->GetParameter(ip)+f->GetParError(ip));
    jac2(0,ip) = fMeanpT2Vsx->Derivative(f->GetParameter(ip),0x0,0.005);
    std::cout << Form("d<pT2>/dp%d = %e ± %e",ip, jac2(0,ip), fMeanpT2Vsx->DerivativeError()) << std::endl;
  }
  
  TMatrixD tmp2(jac2,TMatrixD::kMult,cov);
  TMatrixD meanpT2Err2(tmp2,TMatrixD::kMultTranspose,jac2);
  
  std::cout << std::endl;
  std::cout << Form("<pT> = %5.6f +- %5.6f   (with error propagation)",
                    meanptval,TMath::Sqrt(meanpTErr2(0,0))) << std::endl;
  std::cout << Form("<pT2> = %5.6f +- %5.6f   (with error propagation)",
                    meanpt2val,TMath::Sqrt(meanpT2Err2(0,0))) << std::endl;
  std::cout << std::endl;
  
  return h;
}

