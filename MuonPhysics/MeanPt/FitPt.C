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
//
//0- 1     117.697  2.759   9.044
//1- 2     242.073  3.866  18.778
//2- 3     207.527  3.313  15.971
//3- 4     138.534  2.675  10.836
//4- 5      81.585  1.919   6.353
//5- 6      44.827  1.130   3.484
//6- 8      17.842  0.442   1.453
//8-15       2.419  0.071   0.187

Double_t xmin = 0.;
Double_t xmax = 1.e6;

//------------------------------------------------------------------
Double_t MeanPtVsp0(Double_t *x, Double_t *p)
{
  /// <pT> with x[0] = p0 and p[0] = n
  
  static TF1 *fp0 = new TF1("fp0","x/((1+(x/[0])**2)**[1])",xmin,xmax);
  
  fp0->SetParameters(x[0], p[0]);
  
  return fp0->Mean(xmin,xmax);
}

//------------------------------------------------------------------
Double_t MeanPtVsn(Double_t *x, Double_t *p)
{
  /// <pT> with x[0] = n and p[0] = p0
  
  static TF1 *fn = new TF1("fn","x/((1+(x/[0])**2)**[1])",xmin,xmax);
  
  fn->SetParameters(p[0], x[0]);
  
  return fn->Mean(xmin,xmax);
}

//------------------------------------------------------------------
TH1* FitPt(const char* input="pt.pA.igor.txt")
{
  ifstream in(input);
  char line[1024];
  
  int a,b;
  float value,stat,sys;
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> yerr;
  
  while (in.getline(line,1023,'\n'))
  {
  	sscanf(line,"%d-%2d %e %e %e",&a,&b,&value,&stat,&sys);
//   	stat=0;
//	sys=0;
  	printf("%d-%2d %e %e %e\n",a,b,value,stat,sys);
  	x.push_back(a);
  	y.push_back(value);
  	Double_t error = TMath::Sqrt(stat*stat/value/value+sys*sys/value/value)*value;
  	yerr.push_back(error);
  	  	printf("%d-%2d %5.2f +- %5.2f\n",a,b,value,error);

  }
  
  gStyle->SetOptFit(1);
  
  x.push_back(b);
  
  TH1* h = new TH1F("pt","pt",x.size()-1,&x[0]);
  
  for ( Int_t i = 1; i <= h->GetXaxis()->GetNbins(); ++i ) 
  {
  	h->SetBinContent(i,y[i-1]);
  	h->SetBinError(i,yerr[i-1]);
  }
  

  TCanvas* c1 = new TCanvas(input,input);
  c1->Divide(2,1);
  c1->cd(1);
  
  gPad->SetLogy();
  
  h->Draw("histe");
  
  xmin=h->GetXaxis()->GetXmin();
  xmax=h->GetXaxis()->GetXmax();
  
//  TVirtualFitter::SetDefaultFitter("Minuit");
  
//	C*pT/(1+(pT/p0)^2)n
//   TF1* f = new TF1("f","[0]*x/((1+(x/[1])**2)**[2])",h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax());
//   f->SetParNames("C","p0","n");

  TF1* f = new TF1("f","[0]*x/((1+(x/[1])**2)**[2])",xmin,xmax);

  f->SetParNames("C","p0","n");
  
  f->SetParameter(0,1);
//  f->FixParameter(0,2.49424e+02);
  f->SetParameter(1,0);
//  f->FixParameter(1,4.32592e+00);
  f->SetParameter(2,2);
  
  TFitResultPtr r = h->Fit(f,"MIRES");
  
  r->Print("V");
  
  for ( Int_t i = 1; i <= h->GetXaxis()->GetNbins(); ++i ) 
  {
  	Double_t x1 = h->GetBinLowEdge(i);
  	Double_t x2 = x1+h->GetBinWidth(i);

	std::cout << Form("Pt %2.0f - %2.0f Mean %3.2f",x1,x2,f->Mean(x1,x2)) << std::endl;
  }
  
  Double_t meanptval = f->Mean(xmin,xmax);
  std::cout << Form("Mean Pt %e (%f,%f)",meanptval,xmin,xmax) << std::endl;

//  std::cout << Form("Mean Pt %3.2f (0-8)",f->Mean(0,8)) << std::endl;
//  std::cout << Form("Mean Pt %3.2f (0-20)",f->Mean(0,20)) << std::endl;  	

  // --- mean pT from fit function and mean pT error from error propagation ---
  TMatrixD cov(r->GetCovarianceMatrix());
  
  TMatrixD jac(1,3);
  jac(0,0) = 0.;
  TF1* fMeanpTVspO = new TF1("fMeanpTVspO",MeanPtVsp0,
			     f->GetParameter(1)-f->GetParError(1),f->GetParameter(1)+f->GetParError(1),1);
  fMeanpTVspO->SetParameter(0,f->GetParameter(2));
  jac(0,1) = fMeanpTVspO->Derivative(f->GetParameter(1),0x0,0.005);
  std::cout << Form("d<pT>/dp0 = %e ± %e", jac(0,1), fMeanpTVspO->DerivativeError()) << std::endl;
  TF1* fMeanpTVsn = new TF1("fMeanpTVsn",MeanPtVsn,
			    f->GetParameter(2)-f->GetParError(2),f->GetParameter(2)+f->GetParError(2),1);
  fMeanpTVsn->SetParameter(0,f->GetParameter(1));
  jac(0,2) = fMeanpTVsn->Derivative(f->GetParameter(2),0x0,0.005);
  std::cout << Form("d<pT>/dn = %e ± %e", jac(0,2), fMeanpTVsn->DerivativeError()) << std::endl;
  
  TMatrixD tmp(jac,TMatrixD::kMult,cov);
  TMatrixD meanpTErr2(tmp,TMatrixD::kMultTranspose,jac);
  
  std::cout << Form("Mean Pt %5.6f +- %5.6f   (with error propagation)",
		    meanptval,TMath::Sqrt(meanpTErr2(0,0))) << std::endl;
  
  // --- mean pT and mean pT error from contour plot ---
  c1->cd(2);
  TVirtualPad* p = gPad;
  gPad->Divide(1,2);
  p->cd(1);
  
//   gMinuit->SetErrorDef(4); // 2 sigmas (2^2)
//   TGraph* g12at2sigmas = static_cast<TGraph*>(gMinuit->Contour(100,1,2));
//   g12at2sigmas->SetFillColor(4);
//   g12at2sigmas->Draw("alf"); 

  gMinuit->SetErrorDef(1); // 1 sigmas (1^2)
  TGraph* g12 = static_cast<TGraph*>(gMinuit->Contour(5000,1,2));
  g12->SetFillColor(2);
  g12->Draw("alf"); 

  p->cd(2);
  TH1* hmeanpterror = new TH1F("hmeanpterror","h mean pt error",5000,2.4,3.2);
  hmeanpterror->SetLineColor(1);
  
  Double_t meanPtmin = 1.e9;
  Double_t meanPtmax = -1.e9;
  for ( Int_t i = 1; i <= g12->GetN(); ++i ) 
  {
	  Double_t p0,n;
  
     g12->GetPoint(i,p0,n);
     f->FixParameter(1,p0);
     f->FixParameter(2,n);
     //h->Fit(f);
     Double_t meanpt = f->Mean(xmin,xmax);
    if (meanpt < meanPtmin) meanPtmin = meanpt;
    if (meanpt > meanPtmax) meanPtmax = meanpt;
     hmeanpterror->Fill(meanpt);
  }
  
  hmeanpterror->Draw("histe");
  /*
  std::cout << Form("Mean Pt %5.6f +- %5.6f   (with contour (RMS))",
		    meanptval,hmeanpterror->GetRMS()) << std::endl;

  std::cout << Form("Mean Pt %5.6f +- %5.6f   (with contour (extreme values))",
		    meanptval,0.5*(meanPtmax-meanPtmin)) << std::endl;
  */
  Double_t meanPtmin2 = -1.e9;
  for ( Int_t i = 1; i <= hmeanpterror->GetNbinsX(); ++i ) {
    if (hmeanpterror->GetBinContent(i) > 20.) {
      meanPtmin2 = hmeanpterror->GetBinCenter(i);
      break;
    }
  }
  
  Double_t meanPtmax2 = 1.e9;
  for ( Int_t i = hmeanpterror->GetNbinsX(); i >= 1; --i ) {
    if (hmeanpterror->GetBinContent(i) > 20.) {
      meanPtmax2 = hmeanpterror->GetBinCenter(i);
      break;
    }
  }
  
  std::cout << Form("Mean Pt %5.6f +- %5.6f   (with contour (extreme peaks))",
		    meanptval,0.5*(meanPtmax2-meanPtmin2)) << std::endl;
  
  return h;
}
