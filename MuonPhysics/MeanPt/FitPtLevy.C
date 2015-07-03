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
Double_t MeanPtVsn(Double_t *x, Double_t *p)
{
  /// <pT> with x[0] = n and p[0] = T
  
  static TF1 *fn = new TF1("fn","([0]*x*([1]-1)*([1]-2))/([1]*[2]*([1]*[2]+[3]*([1]-2)))*(1+(sqrt([3]*[3]+x*x)-[3])/([1]*[2]))**(-[1])",xmin,xmax);
  
  fn->SetParameters(1., x[0], p[0], mJPsi);
  
  return fn->Mean(xmin,xmax);
}

//------------------------------------------------------------------
Double_t MeanPtVsT(Double_t *x, Double_t *p)
{
  /// <pT> with x[0] = T and p[0] = n
  
  static TF1 *fT = new TF1("fT","([0]*x*([1]-1)*([1]-2))/([1]*[2]*([1]*[2]+[3]*([1]-2)))*(1+(sqrt([3]*[3]+x*x)-[3])/([1]*[2]))**(-[1])",xmin,xmax);
  
  fT->SetParameters(1., p[0], x[0], mJPsi);
  
  return fT->Mean(xmin,xmax);
}

//------------------------------------------------------------------
Double_t MeanPt2Vsn(Double_t *x, Double_t *p)
{
  /// <pT2> with x[0] = n and p[0] = T
  
  static TF1 *fn2 = new TF1("fn2","([0]*x*([1]-1)*([1]-2))/([1]*[2]*([1]*[2]+[3]*([1]-2)))*(1+(sqrt([3]*[3]+x*x)-[3])/([1]*[2]))**(-[1])",xmin,xmax);
  
  fn2->SetParameters(1., x[0], p[0], mJPsi);
  
  return fn2->Moment(2,xmin,xmax);
}

//------------------------------------------------------------------
Double_t MeanPt2VsT(Double_t *x, Double_t *p)
{
  /// <pT2> with x[0] = T and p[0] = n
  
  static TF1 *fT2 = new TF1("fT2","([0]*x*([1]-1)*([1]-2))/([1]*[2]*([1]*[2]+[3]*([1]-2)))*(1+(sqrt([3]*[3]+x*x)-[3])/([1]*[2]))**(-[1])",xmin,xmax);
  
  fT2->SetParameters(1., p[0], x[0], mJPsi);
  
  return fT2->Moment(2,xmin,xmax);
}

//------------------------------------------------------------------
TH1* FitPtLevy(const char* input="pt.pA.igor.txt", Bool_t norm = kFALSE, Int_t errorType = 2, Double_t delta = 0.1)
{
  // norm = kTRUE if yields have been divided by pT
  // errorType: stat = 0, sys = 1, all = other int
  // delta: range for <pT> scan to <pT> ± delta*<pT>
  
  
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
  
  TF1* f = new TF1("f","([0]*x*([1]-1)*([1]-2))/([1]*[2]*([1]*[2]+[3]*([1]-2)))*(1+(sqrt([3]*[3]+x*x)-[3])/([1]*[2]))**(-[1])",xmin,xmax);
  
  f->SetParNames("C","n","T","M");
  
  f->SetParameters(2.2,16.,0.7,mJPsi);
  f->FixParameter(3,mJPsi);
  
  TFitResultPtr r = h->Fit(f,"MIES","",h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax());
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

  TMatrixD jac(1,4);
  jac(0,0) = 0.;
  TF1* fMeanpTVsn = new TF1("fMeanpTVsn",MeanPtVsn,
                            f->GetParameter(1)-f->GetParError(1),f->GetParameter(1)+f->GetParError(1),1);
  fMeanpTVsn->SetParameter(0,f->GetParameter(2));
  jac(0,1) = fMeanpTVsn->Derivative(f->GetParameter(1),0x0,0.005);
  std::cout << Form("d<pT>/dn = %e ± %e", jac(0,1), fMeanpTVsn->DerivativeError()) << std::endl;
  TF1* fMeanpTVsT = new TF1("fMeanpTVsT",MeanPtVsT,
                            f->GetParameter(2)-f->GetParError(2),f->GetParameter(2)+f->GetParError(2),1);
  fMeanpTVsT->SetParameter(0,f->GetParameter(1));
  jac(0,2) = fMeanpTVsT->Derivative(f->GetParameter(2),0x0,0.005);
  std::cout << Form("d<pT>/dT = %e ± %e", jac(0,2), fMeanpTVsT->DerivativeError()) << std::endl;
  jac(0,3) = 0.;
  
  TMatrixD tmp(jac,TMatrixD::kMult,cov);
  TMatrixD meanpTErr2(tmp,TMatrixD::kMultTranspose,jac);
  
  // --- mean pT**2 from fit function and mean pT**2 error from error propagation ---
  TMatrixD jac2(1,4);
  jac2(0,0) = 0.;
  TF1* fMeanpT2Vsn = new TF1("fMeanpT2Vsn",MeanPt2Vsn,
                             f->GetParameter(1)-f->GetParError(1),f->GetParameter(1)+f->GetParError(1),1);
  fMeanpT2Vsn->SetParameter(0,f->GetParameter(2));
  jac2(0,1) = fMeanpT2Vsn->Derivative(f->GetParameter(1),0x0,0.005);
  std::cout << Form("d<pT2>/dn = %e ± %e", jac2(0,1), fMeanpT2Vsn->DerivativeError()) << std::endl;
  TF1* fMeanpT2VsT = new TF1("fMeanpT2VsT",MeanPt2VsT,
                             f->GetParameter(2)-f->GetParError(2),f->GetParameter(2)+f->GetParError(2),1);
  fMeanpT2VsT->SetParameter(0,f->GetParameter(1));
  jac2(0,2) = fMeanpT2VsT->Derivative(f->GetParameter(2),0x0,0.005);
  std::cout << Form("d<pT2>/dT = %e ± %e", jac2(0,2), fMeanpT2VsT->DerivativeError()) << std::endl;
  jac2(0,3) = 0.;
  
  TMatrixD tmp2(jac2,TMatrixD::kMult,cov);
  TMatrixD meanpT2Err2(tmp2,TMatrixD::kMultTranspose,jac2);
  
  // --- mean pT, mean pT**2 and related errors from contour plot ---
  c1->cd(2);
  
  gMinuit->SetErrorDef(1); // 1 sigmas (1^2)
  TGraph* g12 = static_cast<TGraph*>(gMinuit->Contour(1000,1,2));
  Double_t meanPtmin2 = -1.e9;
  Double_t meanPtmax2 = 1.e9;
  Double_t meanPt2min2 = -1.e9;
  Double_t meanPt2max2 = 1.e9;
  if (g12) {
    
    g12->SetLineColor(2);
    g12->Draw("al");
    Int_t nPoints = g12->GetN();
    
    TH1* hmeanpterror = new TH1F("hmeanpterror","h mean pt error",nPoints/5,meanptval*(1.-delta),meanptval*(1.+delta));
    hmeanpterror->SetLineColor(1);
    TH1* hmeanpt2error = new TH1F("hmeanpt2error","h mean pt**2 error",nPoints/5,meanpt2val*(1.-2.*delta),meanpt2val*(1.+2.*delta));
    hmeanpt2error->SetLineColor(1);
    
    Double_t meanPtmin = 1.e9;
    Double_t meanPtmax = -1.e9;
    Double_t meanPt2min = 1.e9;
    Double_t meanPt2max = -1.e9;
    for ( Int_t i = 1; i <= nPoints; ++i )
    {
      Double_t n,T;
      
      g12->GetPoint(i,n,T);
      f->FixParameter(1,n);
      f->FixParameter(2,T);
      
      Double_t meanpt = f->Mean(xmin,xmax);
      if (meanpt < meanPtmin) meanPtmin = meanpt;
      if (meanpt > meanPtmax) meanPtmax = meanpt;
      hmeanpterror->Fill(meanpt);
      
      Double_t meanpt2 = f->Moment(2,xmin,xmax);
      if (meanpt2 < meanPt2min) meanPt2min = meanpt2;
      if (meanpt2 > meanPt2max) meanPt2max = meanpt2;
      hmeanpt2error->Fill(meanpt2);
    }
    
    c1->cd(3);
    hmeanpterror->Draw("histe");
    c1->cd(4);
    hmeanpt2error->Draw("histe");
    
    /*
     std::cout << Form("Mean Pt %5.6f +- %5.6f   (with contour (RMS))",
     meanptval,hmeanpterror->GetRMS()) << std::endl;
     
     std::cout << Form("<pT> = %5.6f +- %5.6f   (with contour (extreme values))",
     meanptval,0.5*(meanPtmax-meanPtmin)) << std::endl;
     
     std::cout << Form("<pT2> = %5.6f +- %5.6f   (with contour (extreme values))",
     meanpt2val,0.5*(meanPt2max-meanPt2min)) << std::endl;
     */
    
    for ( Int_t i = 1; i <= hmeanpterror->GetNbinsX(); ++i ) {
      if (hmeanpterror->GetBinContent(i) > 20.) {
        meanPtmin2 = hmeanpterror->GetBinCenter(i);
        break;
      }
    }
    
    for ( Int_t i = hmeanpterror->GetNbinsX(); i >= 1; --i ) {
      if (hmeanpterror->GetBinContent(i) > 20.) {
        meanPtmax2 = hmeanpterror->GetBinCenter(i);
        break;
      }
    }
    
    for ( Int_t i = 1; i <= hmeanpt2error->GetNbinsX(); ++i ) {
      if (hmeanpt2error->GetBinContent(i) > 20.) {
        meanPt2min2 = hmeanpt2error->GetBinCenter(i);
        break;
      }
    }
    
    for ( Int_t i = hmeanpt2error->GetNbinsX(); i >= 1; --i ) {
      if (hmeanpt2error->GetBinContent(i) > 20.) {
        meanPt2max2 = hmeanpt2error->GetBinCenter(i);
        break;
      }
    }
    
  }
  
  std::cout << std::endl;
  std::cout << Form("<pT> = %5.6f +- %5.6f   (with error propagation)",
                    meanptval,TMath::Sqrt(meanpTErr2(0,0))) << std::endl;
  
  if (g12) std::cout << Form("<pT> = %5.6f +- %5.6f   (with contour (extreme peaks))",
                             meanptval,0.5*(meanPtmax2-meanPtmin2)) << std::endl;
  
  std::cout << Form("<pT2> = %5.6f +- %5.6f   (with error propagation)",
                    meanpt2val,TMath::Sqrt(meanpT2Err2(0,0))) << std::endl;
  
  if (g12) std::cout << Form("<pT2> = %5.6f +- %5.6f   (with contour (extreme peaks))",
                             meanpt2val,0.5*(meanPt2max2-meanPt2min2)) << std::endl;
  std::cout << std::endl;
  
  return h;
}

