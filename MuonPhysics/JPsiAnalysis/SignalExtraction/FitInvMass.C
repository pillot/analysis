/*  FitInvMass.C
 *
 *  Created by Antoine Lardeux on 01/08/12.
 *  Copyright 2012 Subatech. All rights reserved.
 */

#include <TH1.h>
#include <TAxis.h>
#include <TString.h>
#include <stdio.h>
#include <Riostream.h>
#include <TFile.h>
#include <TF1.h>
#include <TList.h>
#include "TObjString.h"
#include <TCanvas.h>
#include <TMath.h>
#include <TStyle.h>
#include <TPaveText.h>
#include <TMatrixTSym.h>
#include <TFitResult.h>
#include <TGraphErrors.h>
#include <TASImage.h>
#include "AliCounterCollection.h"


void Getparam( Double_t *&par, TString cent, TString tails, TString obs, Int_t iobs);
void FitParam( Int_t nparBkg, TF1 *&fit, Int_t i, TString tails, TString name, Double_t *fix090);
void Param090( Int_t nparBkg, TF1 *fit, Int_t i, TString tails, Double_t *&fix090);

void FitSpectrum(Double_t *&fix1CB2090,Double_t *&fix2CB2090,TString dirOut, Bool_t Mixing, AliCounterCollection *eventCounters, Int_t rebin, TString obs, Int_t iobs, TString cent, TH1F *hSig, Float_t binWidthInit, TString fitBkgFct);
void BackgroundMix(Double_t rangeDo,Double_t rangeUp, Double_t *&parBkg, TH1F *&hSpec, TString fitBkgFct);
void BackgroundRaw(Double_t rangeDo,Double_t rangeUp, Double_t *&parBkg, TH1F *&hSpec);
void Fit1CB2(Double_t *&out1CB2, Double_t *fix1CB2090,Double_t rangeDo, Double_t rangeUp, Double_t binWidth, TF1 *&fitFct1CB2, TF1 *&fitBG1CB2, TF1 *&fit1CB2,TH1F *&hSpec, Int_t iobs, TString fixMS, TString tails, Int_t nparBkg, Double_t *par, Double_t *parBkg, TString fitBkgFct);
void Fit2CB2(Double_t *&out2CB2, Double_t *fix2CB2090,Double_t rangeDo, Double_t rangeUp, Double_t binWidth, TF1 *&fitFct2CB2, TF1 *&fitBG2CB2, TF1 *&fit1CB2, TF1 *&fit2CB2,TH1F *&hSpec, Int_t iobs, TString fixMS, TString tails, Int_t nparBkg, Double_t *par, Double_t *parBkg, TString fitBkgFct);
void Display(TCanvas *&c1, Double_t *val, TString cent, AliCounterCollection *eventCounters, Int_t nparBkg, TH1F *&hSpec, TF1 *fitFct, TF1 *fitBG, TF1 *fit1CB2, TF1 *fit2CB2=0x0);


//##############################################################
//     FUNCTIONS
//##############################################################

//---------------------------------------------------------------------------
Double_t BackgroundVWG(Double_t *x,Double_t *par)
{  
  if (x[0] > 2.95 && x[0] < 3.3) TF1::RejectPoint();
  Double_t sigma = par[2]+par[3]*((x[0]-par[1])/par[1]);
  return par[0]*exp(-(x[0]-par[1])*(x[0]-par[1])/(2.*sigma*sigma));  
}

//---------------------------------------------------------------------------
Double_t BackgroundExpo(Double_t *x,Double_t *par)
{
  if (x[0] > 2.95 && x[0] < 3.3) TF1::RejectPoint();
  return exp(par[0]+par[1]*x[0]);  
}

//---------------------------------------------------------------------------
Double_t BackgroundPol1(Double_t *x,Double_t *par)
{
  if (x[0] > 2.95 && x[0] < 3.3) TF1::RejectPoint();
  return par[0]+par[1]*x[0];
}

//---------------------------------------------------------------------------
Double_t VWG(Double_t *x,Double_t *par)
{  
  Double_t sigma = par[2]+par[3]*((x[0]-par[1])/par[1]);
  return par[0]*exp(-(x[0]-par[1])*(x[0]-par[1])/(2.*sigma*sigma));  
}

//---------------------------------------------------------------------------
Double_t Expo(Double_t *x,Double_t *par)
{
  return exp(par[0]+par[1]*x[0]);
}

//---------------------------------------------------------------------------
Double_t Pol1(Double_t *x,Double_t *par)
{
  return par[0]+par[1]*x[0];      
}

//---------------------------------------------------------------------------
Double_t CrystalBallExtended(Double_t *x,Double_t *par)
{
  //par[0] = Normalization
  //par[1] = mean
  //par[2] = sigma
  //par[3] = alpha
  //par[4] = n
  //par[5] = alpha'
  //par[6] = n'  
  
  Double_t t = (x[0]-par[1])/par[2];
  if (par[3] < 0) t = -t;
  
  Double_t absAlpha = fabs((Double_t)par[3]);
  Double_t absAlpha2 = fabs((Double_t)par[5]);
  
  if (t >= -absAlpha && t < absAlpha2) // gaussian core
  {
	return par[0]*(exp(-0.5*t*t));
  }
  
  if (t < -absAlpha) //left tail
  {
	Double_t a =  TMath::Power(par[4]/absAlpha,par[4])*exp(-0.5*absAlpha*absAlpha);
	Double_t b = par[4]/absAlpha - absAlpha;
	return par[0]*(a/TMath::Power(b - t, par[4]));
  }
  
  if (t >= absAlpha2) //right tail
  {
	
	Double_t c =  TMath::Power(par[6]/absAlpha2,par[6])*exp(-0.5*absAlpha2*absAlpha2);
	Double_t d = par[6]/absAlpha2 - absAlpha2;
	return par[0]*(c/TMath::Power(d + t, par[6]));
  }
  
  return 0. ; 
} 

//---------------------------------------------------------------------------
Double_t fitFunctionCB2_VWG(Double_t *x, Double_t *par)
{
  return VWG(x, par) + CrystalBallExtended(x, &par[4]);
}

//---------------------------------------------------------------------------
Double_t fitFunctionCB2_Expo(Double_t *x, Double_t *par)
{
  return Expo(x, par) + CrystalBallExtended(x, &par[2]);
}

//---------------------------------------------------------------------------
Double_t fitFunctionCB2_Pol1(Double_t *x, Double_t *par)
{
  return Pol1(x, par) + CrystalBallExtended(x, &par[2]);
}

//---------------------------------------------------------------------------
Double_t fitFunction2CB2_VWG(Double_t *x, Double_t *par)
{
  Double_t parPrim[7]={par[11],par[5],par[6],par[7],par[8],par[9],par[10]};
  parPrim[1] = parPrim[1] / 3.097 * 3.686;
  
  return VWG(x, par) + CrystalBallExtended(x, &par[4]) + CrystalBallExtended(x, parPrim);
}

//---------------------------------------------------------------------------
Double_t fitFunction2CB2_Expo(Double_t *x, Double_t *par)
{
  Double_t parPrim[7]={par[9],par[3],par[4],par[5],par[6],par[7],par[8]};
  parPrim[1] = parPrim[1] / 3.097 * 3.686;
  
  return Expo(x, par) + CrystalBallExtended(x, &par[2]) + CrystalBallExtended(x, parPrim);
}

//---------------------------------------------------------------------------
Double_t fitFunction2CB2_Pol1(Double_t *x, Double_t *par)
{
  Double_t parPrim[7]={par[9],par[3],par[4],par[5],par[6],par[7],par[8]};
  parPrim[1] = parPrim[1] / 3.097 * 3.686;
  
  return Pol1(x, par) + CrystalBallExtended(x, &par[2]) + CrystalBallExtended(x, parPrim);
}

//---------------------------------------------------------------------------
void SetProperties(TH1F* h) 
{
  h->SetTitle("");
  h->SetStats(0);
  h->SetMarkerStyle(20);
  h->SetMarkerSize(0.8);
  h->SetXTitle("#mu+#mu- invariant mass (GeV/c^{2})");
  h->GetYaxis()->SetTitleOffset(1.72);  // 1.32
}



//##############################################################
//     TOPIC
//##############################################################

// options
Bool_t both = kFALSE;        // case 2 CB2
Bool_t fixSigmaMC = kTRUE;   // fix sigma evolution to MC
Bool_t counter = kTRUE;      // if counter
Bool_t display = kTRUE;
Bool_t ptCut = kFALSE;       // if pt sharp cut on tracks
Bool_t perf = kFALSE;        // draw performance mode

// Fit Range
Float_t RangeDown = 2.2;
Float_t RangeUp = 5.;
Float_t RangeDoMax = 2.4;


//---------------------------------------------------------------------------
void FitInvMass(Bool_t Mixing = kTRUE, 
				Int_t rebin = 1,
				TString dirOut = "SL")
{
  // Open File
  TString nFile = "MixNorm";
  if (Mixing) nFile = "Signal";
  TFile *file = new TFile(Form("../Normalization/%s/%s.root",dirOut.Data(),nFile.Data()), "read");
  if (!file || !file->IsOpen()) {
    printf("cannot open file Signal.root\n");
    return;
  }
  TList *lhisto = static_cast<TList*>(file->FindObjectAny(Form("%s",nFile.Data())));
  
  // get counters
  AliCounterCollection *eventCounters = new AliCounterCollection();
  if (counter) {
	TFile *file1 = new TFile(Form("../Task/Outputs/%s/MixOutput_MergeNoNorme.root",dirOut.Data()), "read");
	if (!file1 || !file1->IsOpen()) {
	  printf("cannot open file with counter \n");
	  return;
	}
	TList *lCounter = static_cast<TList*>(file1->FindObjectAny("cOut"));
	eventCounters = static_cast<AliCounterCollection*>(lCounter->FindObject("eventCounters"));
	if (!eventCounters) {
	  printf("Could not find back counters in output... \n");
	  return;
	}
  }
  
  // param of Mean and Sigma fix to 0-90
  Double_t *fix1CB2090 = new Double_t[6];
  Double_t *fix2CB2090 = new Double_t[6];

   
  //########################################
  // loop on histogramms
  //########################################
  Int_t nhisto = lhisto->GetEntries();
  for (int i=0; i<nhisto; i++) {
	
	// Get histogramm 
	TH1F *hSig = static_cast<TH1F*>(lhisto->At(i));
	if (!hSig) {
	  printf("cannot access to histogram n°%d\n", i);
	  return;
	}
	TString hname = hSig->GetName();
	TObjArray* arr = hname.Tokenize("_");
	TString cent = arr->At(3)->GetName();
	TString obs = arr->At(1)->GetName();
	TString siobs = arr->At(2)->GetName();
	Int_t iobs = siobs.Atoi();
	
	if (!((obs=="pt" && iobs==0 && (cent=="090" || cent=="010" || cent=="1020" || cent=="2030" || cent=="3040" || cent=="4050" || cent=="5060" || cent=="6070" || cent=="7080" || cent=="8090")) ||
		  (obs=="pt" && iobs>0 && iobs<4 && (cent=="090" || cent=="010" || cent=="1020" || cent=="2030" || cent=="3040" || cent=="4050" || cent=="5060" || cent=="6090")) ||
		  (obs=="pt" && iobs>3 && iobs<11 && (cent=="090" || cent=="010" || cent=="1020" || cent=="020" || cent=="2040" || cent=="4060" || cent=="6090" || cent=="4090")) ||
		  (obs=="y" && iobs==0 && (cent=="090" || cent=="6090" || cent=="010" || cent=="1020" || cent=="2030" || cent=="3040" || cent=="4050" || cent=="5060" || cent=="6070" || cent=="7080" || cent=="8090")) ||
		  (obs=="y" && iobs>0 && iobs<4 && (cent=="090" || cent=="010" || cent=="1020" || cent=="2030" || cent=="3040" || cent=="4050" || cent=="5060" || cent=="6090")) ||
		  (obs=="y" && iobs>3 && cent=="090")) ) continue;	

	// Special cases   
	if (iobs==0 && (cent=="7080" || cent=="8090")) rebin=2;  // cent
	else if (obs=="pt" && iobs==3 && (cent=="4050" || cent=="5060" || cent=="6090")) rebin=2;   // 5<pt<8
	else if (obs=="pt" && (cent=="4090" || cent=="2040") && (i==8 || i==9 || i==10)) rebin=2;   // pt bin
	else rebin=1;
	 
	// Compute binWidth
	Int_t nbin = hSig->GetNbinsX();
	Double_t xmin = hSig->GetXaxis()->GetBinLowEdge(1);
	Double_t xmax = hSig->GetXaxis()->GetBinUpEdge(nbin);
	Double_t binWidthInit = (xmax-xmin)/nbin;
  
	if (Mixing) {
	  //########################################
	  // loop over cases of fit fct for bkg
	  //########################################
	  TList* Bkg = new TList();
	  Bkg->Add(new TObjString("expo"));
	  //Bkg->Add(new TObjString("pol"));
	  TIter nextBkg(Bkg);
	  TObjString* fitBkgFct;
	  while ( ( fitBkgFct = static_cast<TObjString*>(nextBkg()) ) ) {
		FitSpectrum(fix1CB2090,fix2CB2090,dirOut,Mixing,eventCounters,rebin,obs,iobs,cent,hSig,binWidthInit,fitBkgFct->String()); // Fit
	  }
	}
	else FitSpectrum(fix1CB2090,fix2CB2090,dirOut,Mixing,eventCounters,rebin,obs,iobs,cent,hSig,binWidthInit,""); // Fit
  }
  
  // close input file
  file->Close();
}

//---------------------------------------------------------------------------
void FitSpectrum(Double_t *&fix1CB2090, Double_t *&fix2CB2090, TString dirOut, Bool_t Mixing, AliCounterCollection *eventCounters, 
				 Int_t rebin, TString obs, Int_t iobs, TString cent, TH1F *hSig, Float_t binWidthInit, TString fitBkgFct)
{
  // Output file for systematic
  TString study="";
  if (Mixing) study="Mixing_";
  TString testsyst = Form("%s/FitResults_%s%s_%i_cent_%s.txt", dirOut.Data(), study.Data(), obs.Data(), iobs, cent.Data());
  ofstream outputTest(testsyst.Data());
  
  // param of CB2 Tails
  Double_t *par;
  
  //########################################
  // loop over param (embedding, simu, pp)
  //########################################
  TList* param = new TList();
  param->Add(new TObjString("embedding2011"));
  //param->Add(new TObjString("simuJpsi2011"));
  TIter nextParam(param);
  TObjString* tails;
  while ((tails = static_cast<TObjString*>(nextParam()))) {
	
	// get param of Tails
	Getparam(par, cent, tails->String(), obs, iobs);
	
	//Rebin
	Int_t initrebin = rebin;
	// Range
	Double_t rangeDo = RangeDown;  
	Double_t rangeUp = RangeUp;
	
	//########################################
	// loop over cases of Mean and Sigma
	//########################################
	TList* fixPar = new TList();
	fixPar->Add(new TObjString("free"));
	fixPar->Add(new TObjString("0-90"));
	fixPar->Add(new TObjString("psigma"));
	fixPar->Add(new TObjString("msigma"));
	TIter nextfix(fixPar);
	TObjString* fixMS;
	while ( ( fixMS = static_cast<TObjString*>(nextfix()) ) ) {
	  
	  Int_t error=-1, errorBoth=-1;
	  if (!both) errorBoth=1;
	  
	  while (error<0 || errorBoth<0) {
		
		printf("############   %s_%i_cent_%s_Rebin_%i_RangeDo_%.2f_MS_%s \n",obs.Data(),iobs,cent.Data(),initrebin,rangeDo,fixMS->String().Data());
		
		//Rebin if necessary
		TH1F *hSpec = (TH1F*) hSig->Clone();
		hSpec->Rebin(initrebin);
		Double_t binWidth = binWidthInit*initrebin;
		
		// set proprties of TH1
		TString title = hSpec->GetTitle();
		if (title!="") SetProperties(hSpec);
		hSpec->SetYTitle(Form("#Evts / (%.0f MeV/c^{2})",binWidth*1000));
		hSpec->SetAxisRange(rangeDo, rangeUp-0.01);
		
		//#########################
		// Fit Background
		//#########################
		Double_t *parBkg;
		Int_t nparBkg = 2;
		if (Mixing) BackgroundMix(rangeDo,rangeUp,parBkg,hSpec,fitBkgFct);
		else {BackgroundRaw(rangeDo,rangeUp,parBkg,hSpec); nparBkg=4;}
		
		
		//#########################
		// Fit Signal
		//#########################
		TF1 *fitFct1CB2, *fitFct2CB2, *fitBG1CB2, *fitBG2CB2, *fit1CB2, *fit2CB2, *fit22CB2;
		Double_t *out1CB2, *out2CB2;
		TCanvas *c1CB2, *c2CB2;
		
		// 1 CB2
		if (error<0) {
		Fit1CB2(out1CB2,fix1CB2090,rangeDo,rangeUp,binWidth,fitFct1CB2,fitBG1CB2,fit1CB2,hSpec,iobs,fixMS->String(),tails->String(),nparBkg,par,parBkg,fitBkgFct);
		if (cent=="090" && fixMS->String()=="free") Param090(nparBkg,fitFct1CB2,iobs,tails->String(),fix1CB2090); // Get param from fit 0-90% in order to fix parameters
		if (display) Display(c1CB2,out1CB2,cent,eventCounters,nparBkg,hSpec,fitFct1CB2,fitBG1CB2,fit1CB2);
		  // Write
		  if (out1CB2[1]!=0. && out1CB2[6]<2.) {
			error=1;
			outputTest << "1CB2 " << tails->String() << " Msigma " << fixMS->String() << " Njpsi " << out1CB2[0] << " errNjpsi " << out1CB2[1] << " Chi2 " << out1CB2[6] << " Range " << rangeDo << "-" << rangeUp << " Rebin " << initrebin << " SoverB " << out1CB2[7]/out1CB2[9] << " Significance " << out1CB2[7]/TMath::Sqrt(out1CB2[7]+out1CB2[9]) << endl;  // " Chi2loal " << 
			TFile *outputfile1CB2 = new TFile(Form("%s/InvMassCB2_%s%s_%i_Range_%.2f_%.0f_Rebin%i.root",dirOut.Data(),study.Data(),obs.Data(),iobs,rangeDo,rangeUp,initrebin),"update");
			c1CB2->Write(Form("1CB2_%s%s_%s_%s_%i_FixMS_%s", study.Data(), tails->String().Data(), cent.Data(), obs.Data(), iobs, fixMS->String().Data()), TObject::kOverwrite);
			outputfile1CB2->Close();
		  }
		}
		// 2 CB2
		if (both && errorBoth<0) {
		  Fit2CB2(out2CB2,fix2CB2090,rangeDo, rangeUp,binWidth,fitFct2CB2,fitBG2CB2,fit2CB2,fit22CB2,hSpec,iobs,fixMS->String(),tails->String(),nparBkg,par,parBkg,fitBkgFct);
		  if (cent=="090" && fixMS->String()=="free") Param090(nparBkg, fitFct2CB2, iobs, tails->String(),fix2CB2090); // Get param from fit 0-90% in order to fix parameters
		  if (display) Display(c2CB2, out2CB2, cent, eventCounters, nparBkg, hSpec, fitFct2CB2, fitBG2CB2, fit2CB2, fit22CB2); 
		  // Write
		  if (out2CB2[1]!=0. && out2CB2[6]<2.) {
			errorBoth=1;
			outputTest << "2CB2 " << tails->String() << " Msigma " << fixMS->String() << " Njpsi " << out2CB2[0] << " errNjpsi " << out2CB2[1] << " Chi2 " << out2CB2[6] << " Range " << rangeDo << "-" << rangeUp << " Rebin " << initrebin << endl;  // " Chi2loal " << 
			outputTest << "22CB2 " << tails->String() << " Msigma " << fixMS->String() << " Njpsi " << out2CB2[13] << " errNjpsi " << out2CB2[14] << " Chi2 " << out2CB2[19] << " Range " << rangeDo << "-" << rangeUp << " Rebin " << initrebin << endl;  // " Chi2loal " << 
			TFile *outputfile2CB2 = new TFile(Form("%s/InvMass2CB2_%s%s_%i_Range_%.2f_%.0f_Rebin%i.root",dirOut.Data(),study.Data(),obs.Data(),iobs,rangeDo,rangeUp,initrebin),"update");
			c2CB2->Write(Form("2CB2_%s%s_%s_%s_%i_FixMS_%s", study.Data(),tails->String().Data(), cent.Data(), obs.Data(), iobs, fixMS->String().Data()), TObject::kOverwrite);
			outputfile2CB2->Close();  
		  }
		}
		
		/*
		// Bin and range condition loop  (before rebin then reduce range)
		if (error<0 || errorBoth<0) {
		  printf("Fit procedure : next step \n");
		  if (initrebin<4) initrebin+=1;
		  else if (rangeDo < RangeDoMax) {initrebin=rebin; rangeDo+=0.05;}
		  else {
			printf("************************************* \n");
			printf("Fit not converge -> Exit \n");
			printf("************************************* \n");
			break;
		  }
		}
		 */
		
		// Bin and range condition loop (before reduce range then rebin)
		if (error<0 || errorBoth<0) {
		  if (initrebin<4) {
			if (rangeDo < RangeDoMax) rangeDo+=0.05;
			else {initrebin+=1; rangeDo = 2.2;}
		  }
		  else {
			if (rangeDo < RangeDoMax) rangeDo+=0.05;
			else {
			  printf("###################################### \n");
			  printf("Fit not converge -> Exit \n");
			  printf("###################################### \n");
			  break;
			}
		  }
		}
		
		
	  } // end while err Njpsi and reduce Range
	} // end loop on Mean and Sigma condition
  } // end loop on Tails
  
  // close txt file
  outputTest.close();
}

//---------------------------------------------------------------------------
void Getparam( Double_t *&par, TString cent, TString tails, TString obs, Int_t i)
{
  par = new Double_t[4];

  if (tails.Contains("simuJpsi2011")) {  // No Centrality dependence in simuJpsi
	if (ptCut){
	  par[0] = 0.972; par[1] = 6.55; par[2] = 1.67; par[3] = 3.6; 
	}
	else {
	  // defaut value
	  par[0] = 0.907; par[1] = 9.990; par[2] = 1.637; par[3] = 4.480;
	  
	  if(obs=="pt"){
		if (i==0) {  // Integer over Pt and Y  
		  par[0] = 0.907; par[1] = 9.990; par[2] = 1.637; par[3] = 4.480;
		}
		else if (i==1) {  //  0<Pt<2 
		  par[0] = 0.913; par[1] = 10.770; par[2] = 1.815; par[3] = 5.614;
		}
		else if (i==2) {  //  2<Pt<5  
		  par[0] = 0.888; par[1] = 10.719; par[2] = 1.534; par[3] = 5.223;
		}
		else if (i==3) {  // 5<Pt<8
		  par[0] = 0.929; par[1] = 6.870; par[2] = 1.319; par[3] = 3.280;
		}
		else if (i==4) {  //  0<Pt<1 
		  par[0] = 0.920;	par[1] = 9.808;	par[2] = 1.922;	par[3] = 5.732;
		}
		else if (i==5) {  //  1<Pt<2 
		  par[0] = 0.907; par[1] = 11.482; par[2] = 1.757; par[3] = 5.690;
		}
		else if (i==6) {  //  2<Pt<3 
		  par[0] = 0.891; par[1] = 11.970; par[2] = 1.620; par[3] = 5.590;
		}
		else if (i==7) {  //  3<Pt<4  
		  par[0] = 0.892; par[1] = 9.989; par[2] = 1.485; par[3] = 5.794;
		}
		else if (i==8) {  //  4<Pt<5 
		  par[0] = 0.863; par[1] = 9.910; par[2] = 1.388; par[3] = 4.619;
		}
		else if (i==9) {  //  5<Pt<6 
		  par[0] = 0.922; par[1] = 6.962; par[2] = 1.336; par[3] = 3.850;
		}
		else if (i==10) {  //  6<Pt<8 
		  par[0] = 0.935; par[1] = 6.869; par[2] = 1.294; par[3] = 2.887;
		}
		else if (i==11 && i==12) {  //  0<Pt<1 
		  par[0] = 0.920;	par[1] = 9.808;	par[2] = 1.922;	par[3] = 5.732;
		}
		else if (i==13 && i==14) {  //  1<Pt<2 
		  par[0] = 0.907; par[1] = 11.482; par[2] = 1.757; par[3] = 5.690;
		}
		else if (i==15 && i==16) {  //  2<Pt<3 
		  par[0] = 0.891; par[1] = 11.970; par[2] = 1.620; par[3] = 5.590;
		}
		else if (i==17 && i==18) {  //  3<Pt<4  
		  par[0] = 0.892; par[1] = 9.989; par[2] = 1.485; par[3] = 5.794;
		}
		else if (i==19 && i==20) {  //  4<Pt<5 
		  par[0] = 0.863; par[1] = 9.910; par[2] = 1.388; par[3] = 4.619;
		}
		else if (i==21 && i==22) {  //  5<Pt<6 
		  par[0] = 0.922; par[1] = 6.962; par[2] = 1.336; par[3] = 3.850;
		}
		else if (i==23 && i==24) {  //  6<Pt<8 
		  par[0] = 0.935; par[1] = 6.869; par[2] = 1.294; par[3] = 2.887;
		}
	  }
	  else if (obs=="y"){
		if (i==0) {	  // Integer over Pt and Y
		  par[0] = 0.907; par[1] = 9.990; par[2] = 1.637; par[3] = 4.480;
		}
		else if (i==1) {  // 2.5<y<3.
		  par[0] = 0.907; par[1] = 9.990; par[2] = 1.637; par[3] = 4.480;
		}
		else if (i==2) {  // // 3.<y<3.5
		  par[0] = 0.907; par[1] = 9.990; par[2] = 1.637; par[3] = 4.480;
		}
		else if (i==3) {  // // 3.5<y<4.
		  par[0] = 0.907; par[1] = 9.990; par[2] = 1.637; par[3] = 4.480;
		}
		else if (i==4) {  //  2.5<y<2.75 
		  par[0] = 0.804;	par[1] = 6.775;	par[2] = 2.592;	par[3] = 4.196;
		}
		else if (i==5) {  //  2.75<y<3 
		  par[0] = 0.874; par[1] = 8.064; par[2] = 2.107; par[3] = 3.958;
		}
		else if (i==6) {  //  3<y<3.25 
		  par[0] = 0.917; par[1] = 8.956; par[2] = 1.749; par[3] = 4.806;
		}
		else if (i==7) {  //  3.25<y<3.5
		  par[0] = 0.922; par[1] = 11.755; par[2] = 1.529; par[3] = 5.735;
		}
		else if (i==8) {  //  3.5<y<3.75 
		  par[0] = 1.015; par[1] = 11.097; par[2] = 1.376; par[3] = 8.862;
		}
		else if (i==9) {  //  3.75<y<4 
		  par[0] = 1.090; par[1] = 7.711; par[2] = 1.301; par[3] = 5.489;
		}
	  }
	}
  }
  
  else if (tails.Contains("embedding2011")) {
	if (ptCut) {
	  par[0] = 0.972; par[1] = 6.55; par[2] = 1.67; par[3] = 3.6; 
	}
	else {
	  // defaut value
	  par[0] = 0.93; par[1] = 7.04; par[2] = 1.64; par[3] = 4.56;
	  
	  if(obs=="pt"){
		
		if (cent=="020"){
		  if (i==4 || i==11 || i==12) {	//  0<Pt<1 
			par[0] = 0.91; par[1] = 9.51; par[2] = 2.01; par[3] = 4.84;
		  }
		  else if (i==5 || i==13 || i==14) { //  1<Pt<2 
			par[0] = 0.94; par[1] = 8.26; par[2] = 1.77; par[3] = 5.86;    		
		  }
		  else if (i==6 || i==15 || i==16) { //  2<Pt<3 
			par[0] = 0.94; par[1] = 6.37; par[2] = 1.57; par[3] = 5.92;  		
		  }
		  else if (i==7 || i==17 || i==18) { //  3<Pt<4  
			par[0] = 0.93; par[1] = 6.86; par[2] = 1.62; par[3] = 4.23;  		
		  }
		  else if (i==8 || i==19 || i==20) { //  4<Pt<5 
			par[0] = 1.00; par[1] = 5.48; par[2] = 1.41; par[3] = 5.49;  		
		  }
		  else if (i==9 || i==21 || i==22) { //  5<Pt<6 
			par[0] = 1.02; par[1] = 4.72; par[2] = 1.38; par[3] = 4.02; 		
		  }
		  else if (i==10 || i==23 || i==24) { //  6<Pt<8 
			par[0] = 0.80; par[1] = 9.22; par[2] = 1.15; par[3] = 3.51; 		
		  }
		}
		else if (cent=="2040"){
		  if (i==4 || i==11 || i==12) { //  0<Pt<1 
			par[0] = 0.92; par[1] = 10.22; par[2] = 2.34; par[3] = 2.67;
		  }
		  else if (i==5 || i==13 || i==14) { //  1<Pt<2 
			par[0] = 0.93; par[1] = 7.52; par[2] = 1.82; par[3] = 4.88;    		
		  }
		  else if (i==6 || i==15 || i==16) { //  2<Pt<3 
			par[0] = 0.90; par[1] = 7.75; par[2] = 1.63; par[3] = 5.06;  		
		  }
		  else if (i==7 || i==17 || i==18) { //  3<Pt<4  
			par[0] = 0.97; par[1] = 5.64; par[2] = 1.57; par[3] = 4.72;  		
		  }
		  else if (i==8 || i==19 || i==20) { //  4<Pt<5 
			par[0] = 0.95; par[1] = 5.66; par[2] = 1.39; par[3] = 5.44;  		
		  }
		  else if (i==9 || i==21 || i==22) { //  5<Pt<6 
			par[0] = 0.93; par[1] = 5.95; par[2] = 1.32; par[3] = 3.90; 		
		  }
		  else if (i==10 || i==23 || i==24) { //  6<Pt<8 
			par[0] = 0.85; par[1] = 5.77; par[2] = 1.13; par[3] = 3.77; 		
		  }
		}
		else if (cent=="4090"){
		  if (i==4 || i==11 || i==12) {	//  0<Pt<1 
			par[0] = 0.95; par[1] = 7.59; par[2] = 1.90; par[3] = 5.94;
		  }
		  else if (i==5 || i==13 || i==14) { //  1<Pt<2 
			par[0] = 0.94; par[1] = 7.51; par[2] = 1.85; par[3] = 4.76;    		
		  }
		  else if (i==6 || i==15 || i==16) { //  2<Pt<3 
			par[0] = 0.96; par[1] = 5.97; par[2] = 1.60; par[3] = 6.01;  		
		  }
		  else if (i==7 || i==17 || i==18) { //  3<Pt<4  
			par[0] = 0.95; par[1] = 5.54; par[2] = 1.48; par[3] = 5.72;  		
		  }
		  else if (i==8 || i==19 || i==20) { //  4<Pt<5 
			par[0] = 0.90; par[1] = 6.42; par[2] = 1.25; par[3] = 7.34;  		
		  }
		  else if (i==9 || i==21 || i==22) { //  5<Pt<6 
			par[0] = 0.90; par[1] = 6.74; par[2] = 1.30; par[3] = 5.05; 		
		  }
		  else if (i==10 || i==23 || i==24) { //  6<Pt<8 
			par[0] = 0.86; par[1] = 6.76; par[2] = 1.11; par[3] = 4.17; 		
		  }
		}
		/*
		else if (cent=="6090"){
		  if (i==0) {
			//  0<Pt<2 
			par[0] = 0.93; par[1] = 7.58; par[2] = 1.84; par[3] = 5.57;
		  }
		  else if (i==1) {
			//  2<Pt<4 
			par[0] = 0.93; par[1] = 6.08; par[2] = 1.56; par[3] = 5.69;    		
		  }
		  else if (i==2) {
			//  4<Pt<8 
			par[0] = 0.87; par[1] = 6.85; par[2] = 1.19; par[3] = 5.32;  		
		  }
		}
		 */
		else { 
		  if (i==0) {
			// Integer over Pt and Y
			par[0] = 0.93; par[1] = 7.04; par[2] = 1.64; par[3] = 4.56;
		  }
		  else if (i==1) {	//  0<Pt<2 
			par[0] = 0.92; par[1] = 8.52; par[2] = 1.90; par[3] = 4.67;
		  }
		  else if (i==2) {	//  2<Pt<5 
			par[0] = 0.94; par[1] = 6.45; par[2] = 1.51; par[3] = 5.78;
		  }
		  else if (i==3) {	//  5<Pt<8 
			par[0] = 0.93; par[1] = 5.79; par[2] = 1.29; par[3] = 3.62;
		  }
		  else if (i==4 || i==11 || i==12) { //  0<Pt<1 
			par[0] = 0.93; par[1] = 8.42; par[2] = 2.03; par[3] = 4.25;
		  }
		  else if (i==5 || i==13 || i==14) { //  1<Pt<2 
			par[0] = 0.94; par[1] = 7.46; par[2] = 1.83; par[3] = 4.79;    		
		  }
		  else if (i==6 || i==15 || i==16) { //  2<Pt<3 
			par[0] = 0.94; par[1] = 6.39; par[2] = 1.59; par[3] = 5.83;  		
		  }
		  else if (i==7 || i==17 || i==18) { //  3<Pt<4  
			par[0] = 0.95; par[1] = 5.80; par[2] = 1.54; par[3] = 4.84;  		
		  }
		  else if (i==8 || i==19 || i==20) { //  4<Pt<5 
			par[0] = 0.93; par[1] = 5.99; par[2] = 1.29; par[3] = 6.78;  		
		  }
		  else if (i==9 || i==21 || i==22) { //  5<Pt<6 
			par[0] = 0.91; par[1] = 6.35; par[2] = 1.31; par[3] = 4.32; 		
		  }
		  else if (i==10 || i==23 || i==24) { //  6<Pt<8 
			par[0] = 0.84; par[1] = 7.16; par[2] = 1.09; par[3] = 4.38; 		
		  }
		}
	  }
	  else if (obs=="y"){
		if (cent=="090") {
		  if (i==0) {	// Integer over Pt and Y 
			par[0] = 0.93; par[1] = 7.04; par[2] = 1.64; par[3] = 4.56;
		  }
		  else if (i==1) {	// 2.5<Pt<3.
			par[0] = 0.94; par[1] = 4.32; par[2] = 2.09; par[3] = 4.53;
		  }
		  else if (i==2) {	// 3.<Pt<3.5
			par[0] = 0.95; par[1] = 7.00; par[2] = 1.65; par[3] = 4.68;
		  }
		  else if (i==3) {	// 3.5<Pt<4.
			par[0] = 1.01; par[1] = 8.74; par[2] = 1.35; par[3] = 5.79;
		  }
		  else if (i==4) {	//  2.5<Pt<2.75 
			par[0] = 0.82;	par[1] = 4.25;	par[2] = 2.56;	par[3] = 2.55;
		  }
		  else if (i==5) {	//  2.75<Pt<3 
			par[0] = 0.96; par[1] = 4.79; par[2] = 1.99; par[3] = 5.44;		
		  }
		  else if (i==6) {	//  3<Pt<3.25 
			par[0] = 0.95; par[1] = 6.12; par[2] = 1.74; par[3] = 4.90;
		  }
		  else if (i==7) {	//  3.25<Pt<3.5
			par[0] = 0.93; par[1] = 9.22; par[2] = 1.54; par[3] = 5.45;
		  }
		  else if (i==8) {	//  3.5<Pt<3.75 
			par[0] = 0.96; par[1] = 10.79; par[2] = 1.35; par[3] = 6.69;
		  }
		  else if (i==9) {	//  3.75<Pt<4 
			par[0] = 1.10; par[1] = 7.52; par[2] = 1.27; par[3] = 5.84;
		  }
		}
	  }
	}
  }
}

//---------------------------------------------------------------------------
void FitParam( Int_t nparBkg, TF1 *&fit, Int_t i, TString tails, TString name, Double_t *fix090)
{
  Double_t beta = 1.;
  
  if (name.Contains("free")) {
	fit->SetParameter(nparBkg+1, 3.13); 
	fit->SetParLimits(nparBkg+1, 3.05,3.15); 
	fit->SetParameter(nparBkg+2, 0.08);
	fit->SetParLimits(nparBkg+2, 0.05, 0.15);
  }
  else {
	if (tails=="simuJpsi2011") {
	  if (fixSigmaMC && i>3 && i<11) {
		Double_t sigmaIntPt = 0.0686;
		Double_t sigma[7] = {0.0666, 0.0676, 0.0688, 0.0697, 0.0686, 0.0709, 0.0715}; 
		beta = sigma[i-4] / sigmaIntPt;
	  }
	  fit->FixParameter(nparBkg+1, fix090[0]);
	  if (name.Contains("0-90")) fit->FixParameter(nparBkg+2, beta*fix090[1]);
	  else if (name.Contains("psigma")) fit->FixParameter(nparBkg+2, beta*(fix090[1]+fix090[2]));
	  else if (name.Contains("msigma")) fit->FixParameter(nparBkg+2, beta*(fix090[1]-fix090[2]));
	}
	else if (tails=="embedding2011") {
	  if (fixSigmaMC && i>3 && i<11) {
		Double_t sigmaIntPt = 0.0683;
		Double_t sigma[7] = {0.067, 0.0677, 0.0682, 0.0694, 0.0696, 0.0679, 0.0662};
		beta = sigma[i-4] / sigmaIntPt;
	  }
	  fit->FixParameter(nparBkg+1, fix090[3]);
	  if (name.Contains("0-90")) fit->FixParameter(nparBkg+2, beta*fix090[4]); 
	  else if (name.Contains("psigma")) fit->FixParameter(nparBkg+2, beta*(fix090[4]+fix090[5]));
	  else if (name.Contains("msigma")) fit->FixParameter(nparBkg+2, beta*(fix090[4]-fix090[5]));
	}
  }
}  

//---------------------------------------------------------------------------
void Param090( Int_t nparBkg, TF1 *fit, Int_t i, TString tails, Double_t *&fix090)
{
  if (tails=="simuJpsi2011") {
	fix090[0] = fit->GetParameter(nparBkg+1); 
	if (fixSigmaMC && i>3 && i<11) {
	  fix090[1] = 0.076;  // sigma from 0<pt<8 cent 0-90%
	  fix090[2] = 0.002;  // err sigma from 0<pt<8 cent 0-90%
	}
	else {
	  fix090[1] = fit->GetParameter(nparBkg+2);
	  fix090[2] = fit->GetParError(nparBkg+2);
	}
  }
  else if (tails=="embedding2011") {
	fix090[3] = fit->GetParameter(nparBkg+1); 
	if (fixSigmaMC && i>3 && i<11) {
	  fix090[4] = 0.076;  // sigma from 0<pt<8 cent 0-90%
	  fix090[5] = 0.002;  // err sigma from 0<pt<8 cent 0-90%
	}
	else {
	  fix090[4] = fit->GetParameter(nparBkg+2);
	  fix090[5] = fit->GetParError(nparBkg+2);
	}
  }
}

//---------------------------------------------------------------------------
void BackgroundMix(Double_t rangeDo, Double_t rangeUp, Double_t *&parBkg, TH1F *&hSpec, TString fitBkgFct)
{
  printf("Fit Background \n");

  parBkg = new Double_t[2];
  TF1 *fitBkg;
  if ( fitBkgFct=="expo") fitBkg = new TF1("fitBkg",BackgroundExpo,rangeDo,rangeUp,2);
  else fitBkg = new TF1("fitBkg",BackgroundPol1,rangeDo,rangeUp,2);
  hSpec->Fit(fitBkg,"BEMQR0");
  fitBkg->GetParameters(parBkg);
  delete fitBkg;
}

//---------------------------------------------------------------------------
void BackgroundRaw(Double_t rangeDo, Double_t rangeUp, Double_t *&parBkg, TH1F *&hSpec)
{
  printf("Fit Background \n");

  parBkg = new Double_t[4];
  TF1 *fitBkg = new TF1("fitBkg",BackgroundVWG,rangeDo,rangeUp,4);
  fitBkg->SetParameter(0, 77363.4);
  fitBkg->SetParameter(1, 1.80243);
  fitBkg->SetParameter(2, 0.586137);
  fitBkg->SetParameter(3, 0.268486);
  hSpec->Fit(fitBkg,"BEMQR0");
  fitBkg->GetParameters(parBkg);
  delete fitBkg;  
}
	  	  
//---------------------------------------------------------------------------
void Fit1CB2(Double_t *&out1CB2, Double_t *fix1CB2090, Double_t rangeDo, Double_t rangeUp, Double_t binWidth, TF1 *&fitFct1CB2, TF1 *&fitBG1CB2, TF1 *&fit1CB2, 
			 TH1F *&hSpec, Int_t iobs, TString fixMS, TString tails, Int_t nparBkg, Double_t *par, Double_t *parBkg, TString fitBkgFct)
{
  printf("Fit Signal 1CB2\n");
  const Int_t npartTot = nparBkg+7;
  
  // set params of fitFct CB2
  if (fitBkgFct=="") fitFct1CB2 = new TF1("fitFct1CB2",fitFunctionCB2_VWG,rangeDo,rangeUp,npartTot);
  else if (fitBkgFct=="expo") fitFct1CB2 = new TF1("fitFct1CB2",fitFunctionCB2_Expo,rangeDo,rangeUp,npartTot);
  else if (fitBkgFct=="pol") fitFct1CB2 = new TF1("fitFct1CB2",fitFunctionCB2_Pol1,rangeDo,rangeUp,npartTot);
  fitFct1CB2->SetLineColor(4);
  fitFct1CB2->SetParameter(0, parBkg[0]);
  fitFct1CB2->SetParameter(1, parBkg[1]);
  if (nparBkg==4) {
	fitFct1CB2->SetParameter(2, parBkg[2]);
    fitFct1CB2->SetParameter(3, parBkg[3]);
  }
  fitFct1CB2->SetParameter(nparBkg, 10.);   // norm[i]
  fitFct1CB2->SetParLimits(nparBkg, 0., 10000.);
  FitParam(nparBkg, fitFct1CB2, iobs, tails, fixMS, fix1CB2090);
  fitFct1CB2->FixParameter(nparBkg+3, par[0]);   
  fitFct1CB2->FixParameter(nparBkg+4, par[1]);   
  fitFct1CB2->FixParameter(nparBkg+5, par[2]);   
  fitFct1CB2->FixParameter(nparBkg+6, par[3]);
  TFitResultPtr r1CB2 = hSpec->Fit(fitFct1CB2,"BEMQRS0");
  
  // Get param of Fit
  Double_t *parFct1CB2 = new Double_t[npartTot];
  fitFct1CB2->GetParameters(parFct1CB2);
  Double_t *par1CB2 = new Double_t[7];
  par1CB2[0]=parFct1CB2[nparBkg];
  par1CB2[1]=parFct1CB2[nparBkg+1];
  par1CB2[2]=parFct1CB2[nparBkg+2];
  par1CB2[3]=parFct1CB2[nparBkg+3];
  par1CB2[4]=parFct1CB2[nparBkg+4];
  par1CB2[5]=parFct1CB2[nparBkg+5];
  par1CB2[6]=parFct1CB2[nparBkg+6];
  Double_t *parBkg1CB2 = new Double_t[nparBkg];
  parBkg1CB2[0]=parFct1CB2[0];
  parBkg1CB2[1]=parFct1CB2[1];
  if (nparBkg==4) {parBkg1CB2[2]=parFct1CB2[2];parBkg1CB2[3]=parFct1CB2[3];}
  
  // set params of crystal ball 2
  fit1CB2 = new TF1("fit1CB2", CrystalBallExtended, rangeDo, rangeUp, 7);
  fit1CB2->SetParameters(par1CB2);
  fit1CB2->SetLineColor(2);  
  TMatrixDSym covMat1CB2 = r1CB2->GetCovarianceMatrix().GetSub(nparBkg, npartTot-1, nparBkg, npartTot-1);
  Double_t *covMatArray1CB2 = new Double_t[49];
  covMatArray1CB2 = covMat1CB2.GetMatrixArray();
  
  // set params of background
  if ( fitBkgFct=="") fitBG1CB2 = new TF1("fitBG1CB2",VWG,rangeDo,rangeUp,nparBkg);
  else if ( fitBkgFct=="expo") fitBG1CB2 = new TF1("fitBG1CB2",Expo,rangeDo,rangeUp,nparBkg);
  else if ( fitBkgFct=="pol") fitBG1CB2 = new TF1("fitBG1CB2",Pol1,rangeDo,rangeUp,nparBkg);
  fitBG1CB2->SetParameters(parBkg1CB2);
  fitBG1CB2->SetLineColor(4);
  fitBG1CB2->SetLineStyle(2);
  TMatrixDSym covMatBG1CB2 = r1CB2->GetCovarianceMatrix().GetSub(0, nparBkg-1, 0, nparBkg-1);
  Double_t *covMatArrayBG1CB2 = new Double_t[nparBkg*nparBkg];
  covMatArrayBG1CB2 = covMatBG1CB2.GetMatrixArray();
  
  // Compute
  out1CB2 = new Double_t[13];
  
  out1CB2[0] = fit1CB2->Integral(rangeDo, rangeUp)/binWidth;														  // nJpsi
  out1CB2[1] = fit1CB2->IntegralError(rangeDo,rangeUp,par1CB2,covMatArray1CB2)/binWidth;							  // err nJpsi
  out1CB2[2] = fitFct1CB2->GetParameter(nparBkg+1);																	  // mass Jpsi
  out1CB2[3] = fitFct1CB2->GetParError(nparBkg+1);																	  // err mass Jpsi
  out1CB2[4] = fitFct1CB2->GetParameter(nparBkg+2);																	  // sigma Jpsi
  out1CB2[5] = (fitFct1CB2->GetParError(nparBkg+2))*1000;															  // err sigma Jpsi
  out1CB2[6] = fitFct1CB2->GetChisquare()/fitFct1CB2->GetNDF();														  // Chi 2 
  
  Double_t massm3sigma1CB2 = out1CB2[2]-3*out1CB2[4];																																							  // mass - 3 sigma
  Double_t massp3sigma1CB2 = out1CB2[2]+3*out1CB2[4];																																							  // mass + 3 sigma
  out1CB2[7] = fit1CB2->Integral(massm3sigma1CB2,massp3sigma1CB2)/binWidth;																																		  // nJpsi 3 sigma
  out1CB2[8] = fit1CB2->IntegralError(massm3sigma1CB2,massp3sigma1CB2,par1CB2,covMatArray1CB2)/binWidth;																										  // err nJpsi 3 sigma
  out1CB2[9] = fitBG1CB2->Integral(massm3sigma1CB2,massp3sigma1CB2)/binWidth;																																	  // nBkg 3 sigma
  out1CB2[10] = fitBG1CB2->IntegralError(massm3sigma1CB2,massp3sigma1CB2,parBkg1CB2,covMatArrayBG1CB2)/binWidth;																								  // err nBkg 3 sigma
  out1CB2[11] = (out1CB2[7]/out1CB2[9])*TMath::Sqrt( (out1CB2[8]/out1CB2[7])*(out1CB2[8]/out1CB2[7])+(out1CB2[10]/out1CB2[9])*(out1CB2[10]/out1CB2[9]) );														  // err S over B
  out1CB2[12] = (out1CB2[7]/TMath::Sqrt(out1CB2[7]+out1CB2[9]))/4/(out1CB2[7]*out1CB2[7]+out1CB2[9]*out1CB2[9])*(out1CB2[7]*out1CB2[7]+out1CB2[9]*out1CB2[9]+2*out1CB2[9]/out1CB2[7]*out1CB2[7]*out1CB2[7]);      // err significance
  //out1CB2[12] = (out1CB2[7]/TMath::Sqrt(out1CB2[7]+out1CB2[9]))*((out1CB2[8]/out1CB2[7])+((out1CB2[8]+out1CB2[10])/(2*(out1CB2[7]+out1CB2[9]))));																  // old err significance

}  
	  
//---------------------------------------------------------------------------
void Fit2CB2(Double_t *&out2CB2, Double_t *fix2CB2090, Double_t rangeDo, Double_t rangeUp, Double_t binWidth, TF1 *&fitFct2CB2, TF1 *&fitBG2CB2, TF1 *&fit1CB2, TF1 *&fit2CB2,
			 TH1F *&hSpec, Int_t iobs, TString fixMS, TString tails, Int_t nparBkg, Double_t *par, Double_t *parBkg, TString fitBkgFct)
{
  printf("Fit Signal 2CB2\n");
  const Int_t npartTot = nparBkg+7+1;
  
  // set params of fitFct CB2
  if ( fitBkgFct=="") fitFct2CB2 = new TF1("fitFct2CB2",fitFunction2CB2_VWG,rangeDo,rangeUp,npartTot);
  else if ( fitBkgFct=="expo") fitFct2CB2 = new TF1("fitFct2CB2",fitFunction2CB2_Expo,rangeDo,rangeUp,npartTot);
  else if ( fitBkgFct=="pol") fitFct2CB2 = new TF1("fitFct2CB2",fitFunction2CB2_Pol1,rangeDo,rangeUp,npartTot);
  fitFct2CB2->SetLineColor(4);
  fitFct2CB2->SetParameter(0, parBkg[0]);
  fitFct2CB2->SetParameter(1, parBkg[1]);
  if (nparBkg==4) {
	fitFct2CB2->SetParameter(2, parBkg[2]);
    fitFct2CB2->SetParameter(3, parBkg[3]);
  }
  fitFct2CB2->SetParameter(nparBkg, 10.);  // norm[i]
  fitFct2CB2->SetParLimits(nparBkg, 0., 10000.);
  FitParam(nparBkg, fitFct2CB2, iobs, tails, fixMS, fix2CB2090);
  fitFct2CB2->FixParameter(nparBkg+3, par[0]);   
  fitFct2CB2->FixParameter(nparBkg+4, par[1]);   
  fitFct2CB2->FixParameter(nparBkg+5, par[2]);   
  fitFct2CB2->FixParameter(nparBkg+6, par[3]);
  fitFct2CB2->SetParameter(npartTot-1, 1.);
  fitFct2CB2->SetParLimits(npartTot-1, 0., 1000.);
  TFitResultPtr r2CB2 = hSpec->Fit(fitFct2CB2,"BEMQRS0");
  
  // Get param of Fit
  Double_t *parFct2CB2 = new Double_t[npartTot];
  fitFct2CB2->GetParameters(parFct2CB2);
  Double_t *par1CB2 = new Double_t[7];
  par1CB2[0]=parFct2CB2[nparBkg];
  par1CB2[1]=parFct2CB2[nparBkg+1];
  par1CB2[2]=parFct2CB2[nparBkg+2];
  par1CB2[3]=parFct2CB2[nparBkg+3];
  par1CB2[4]=parFct2CB2[nparBkg+4];
  par1CB2[5]=parFct2CB2[nparBkg+5];
  par1CB2[6]=parFct2CB2[nparBkg+6];
  Double_t *parBkg2CB2 = new Double_t[nparBkg];
  parBkg2CB2[0]=parFct2CB2[0];
  parBkg2CB2[1]=parFct2CB2[1];
  if (nparBkg==4) {parBkg2CB2[2]=parFct2CB2[2];parBkg2CB2[3]=parFct2CB2[3];}  
  
  // set params of 1st crystal ball 2
  fit1CB2 = new TF1("fit1CB2", CrystalBallExtended, rangeDo, rangeUp, 7);
  fit1CB2->SetParameters(par1CB2);
  fit1CB2->SetLineColor(2);  
  TMatrixDSym covMat1CB2 = r2CB2->GetCovarianceMatrix().GetSub(nparBkg, nparBkg+6, nparBkg, nparBkg+6);
  Double_t *covMatArray1CB2 = new Double_t[49];
  covMatArray1CB2 = covMat1CB2.GetMatrixArray();  
  
  // set params of 2nd crystal ball 2 
  fit2CB2 = new TF1("fit2CB2", CrystalBallExtended, rangeDo, rangeUp, 7);
  Double_t *par2CB2 = new Double_t[7];
  par2CB2[0] = parFct2CB2[npartTot-1];
  par2CB2[1] = parFct2CB2[nparBkg+1] / 3.097 * 3.686;
  par2CB2[2] = parFct2CB2[nparBkg+2];
  par2CB2[3] = parFct2CB2[nparBkg+3];
  par2CB2[4] = parFct2CB2[nparBkg+4];
  par2CB2[5] = parFct2CB2[nparBkg+5];
  par2CB2[6] = parFct2CB2[nparBkg+6];
  fit2CB2->SetParameters(par2CB2);
  fit2CB2->SetLineColor(6);  
  TMatrixDSym covMat2CB2 = r2CB2->GetCovarianceMatrix().GetSub(npartTot-1, npartTot-1, npartTot-1, npartTot-1);
  Double_t* covMatArray2CB2Temp = covMat2CB2.GetMatrixArray();
  Double_t *covMatArray2CB2 = new Double_t[49];
  for (Int_t j=0; j<49; j++) covMatArray2CB2[j] = 0.; 
  covMatArray2CB2[0] = covMatArray2CB2Temp[0];
  
  // set params of background
  if ( fitBkgFct=="") fitBG2CB2 = new TF1("fitBG2CB2",VWG,rangeDo,rangeUp,nparBkg);
  else if ( fitBkgFct=="expo") fitBG2CB2 = new TF1("fitBG2CB2",Expo,rangeDo,rangeUp,nparBkg);
  else if ( fitBkgFct=="pol") fitBG2CB2 = new TF1("fitBG2CB2",Pol1,rangeDo,rangeUp,nparBkg);
  fitBG2CB2->SetParameters(parBkg2CB2);
  fitBG2CB2->SetLineColor(4);
  fitBG2CB2->SetLineStyle(2);
  TMatrixDSym covMatBG2CB2 = r2CB2->GetCovarianceMatrix().GetSub(0, nparBkg-1, 0, nparBkg-1);
  Double_t *covMatArrayBG2CB2 = new Double_t[nparBkg*nparBkg];
  covMatArrayBG2CB2 = covMatBG2CB2.GetMatrixArray();
  
  
  out2CB2 = new Double_t[26];
  // Compute 1CB2  
  out2CB2[0] = fit1CB2->Integral(rangeDo, rangeUp)/binWidth;										// nJpsi
  out2CB2[1] = fit1CB2->IntegralError(rangeDo,rangeUp,par1CB2,covMatArray1CB2)/binWidth;			// err nJpsi
  out2CB2[2] = fitFct2CB2->GetParameter(nparBkg+1);													// mass Jpsi
  out2CB2[3] = fitFct2CB2->GetParError(nparBkg+1);													// err mass Jpsi
  out2CB2[4] = fitFct2CB2->GetParameter(nparBkg+2);													// sigma Jpsi
  out2CB2[5] = (fitFct2CB2->GetParError(nparBkg+2))*1000;											// err sigma Jpsi
  out2CB2[6] = fitFct2CB2->GetChisquare()/fitFct2CB2->GetNDF();										// Chi 2 
  Double_t massm3sigma1CB2 = out2CB2[2]-3*out2CB2[4];																								// mass - 3 sigma
  Double_t massp3sigma1CB2 = out2CB2[2]+3*out2CB2[4];																								// mass + 3 sigma
  out2CB2[7] = fit1CB2->Integral(massm3sigma1CB2,massp3sigma1CB2)/binWidth;																			// nJpsi 3 sigma
  out2CB2[8] = fit1CB2->IntegralError(massm3sigma1CB2,massp3sigma1CB2,par1CB2,covMatArray1CB2)/binWidth;											// err nJpsi 3 sigma
  out2CB2[9] = fitBG2CB2->Integral(massm3sigma1CB2,massp3sigma1CB2)/binWidth;																		// nBkg 3 sigma
  out2CB2[10] = fitBG2CB2->IntegralError(massm3sigma1CB2,massp3sigma1CB2,parBkg2CB2,covMatArrayBG2CB2)/binWidth;									// err nBkg 3 sigma
  out2CB2[11] = (out2CB2[7]/out2CB2[9])*TMath::Sqrt( (out2CB2[8]/out2CB2[7])*(out2CB2[8]/out2CB2[7])+(out2CB2[10]/out2CB2[9])*(out2CB2[10]/out2CB2[9]) );														  // err S over B
  out2CB2[12] = (out2CB2[7]/TMath::Sqrt(out2CB2[7]+out2CB2[9]))/4/(out2CB2[7]*out2CB2[7]+out2CB2[9]*out2CB2[9])*(out2CB2[7]*out2CB2[7]+out2CB2[9]*out2CB2[9]+2*out2CB2[9]/out2CB2[7]*out2CB2[7]*out2CB2[7]);      // err significance
  //Compute 2CB2
  out2CB2[13] = fit2CB2->Integral(rangeDo, rangeUp)/binWidth;
  out2CB2[14] = fit2CB2->IntegralError(rangeDo,rangeUp,par2CB2,covMatArray2CB2)/binWidth; 
  out2CB2[15] = fitFct2CB2->GetParameter(nparBkg+1) / 3.097 * 3.686;
  out2CB2[16] = 0.;
  out2CB2[17] = fitFct2CB2->GetParameter(nparBkg+2);
  out2CB2[18] = 0.;
  out2CB2[19] = fitFct2CB2->GetChisquare()/fitFct2CB2->GetNDF();
  Double_t massm3sigma2CB2 = out2CB2[15]-3*out2CB2[17];
  Double_t massp3sigma2CB2 = out2CB2[15]+3*out2CB2[17];
  out2CB2[20] = fit2CB2->Integral(massm3sigma2CB2,massp3sigma2CB2)/binWidth;
  out2CB2[21] = fit2CB2->IntegralError(massm3sigma2CB2,massp3sigma2CB2,par2CB2,covMatArray2CB2)/binWidth;
  out2CB2[22] = fitBG2CB2->Integral(massm3sigma2CB2,massp3sigma2CB2)/binWidth;
  out2CB2[23] = fitBG2CB2->IntegralError(massm3sigma2CB2,massp3sigma2CB2,parBkg2CB2,covMatArrayBG2CB2)/binWidth;
  out2CB2[24] = (out2CB2[7]/out2CB2[9])*TMath::Sqrt( (out2CB2[8]/out2CB2[7])*(out2CB2[8]/out2CB2[7])+(out2CB2[10]/out2CB2[9])*(out2CB2[10]/out2CB2[9]) );														  // err S over B
  out2CB2[25] = (out2CB2[7]/TMath::Sqrt(out2CB2[7]+out2CB2[9]))/4/(out2CB2[7]*out2CB2[7]+out2CB2[9]*out2CB2[9])*(out2CB2[7]*out2CB2[7]+out2CB2[9]*out2CB2[9]+2*out2CB2[9]/out2CB2[7]*out2CB2[7]*out2CB2[7]);      // err significance
  
}
  
//---------------------------------------------------------------------------
void Display(TCanvas *&c1, Double_t *val, TString cent, AliCounterCollection *eventCounters, Int_t nparBkg, TH1F *&hSpec, TF1 *fitFct, TF1 *fitBG, TF1 *fit1CB2, TF1 *fit2CB2)
{
  // TString PaveText
  TString nbevents = "N_{evts MUL} = ";
  TString centrality = "centrality bin : ";
  TString nbJpsi = "N_{J/#psi} = ";
  TString nbprim = "N_{#psi'} = ";
  TString masse = "Mass = ";
  TString sigma = "#sigma_{J/#psi} = ";
  TString sigmaPrim = "#sigma_{#psi'} = ";
  TString alphaCB = "#alpha_{CB} = ";
  TString nCB = "n_{CB} = ";
  TString alphaUpCB = "#alpha.up_{CB} = ";
  TString nUpCB = "n.up_{CB} = ";
  TString SoverB = "S/B_{#pm3#sigma} = ";
  TString signalPrim = "S/ #sqrt{S+B} #approx ";
  TString signal = "S/ #sqrt{S+B} = ";
  TString chi2 = "#chi^{2}/nDoF = ";
  TString ratio = "N_{#psi'}/N_{J/#psi} = ";
  TDatime dt;
  TString today = Form("%02i/%02i/%4i", dt.GetDay(), dt.GetMonth(), dt.GetYear());
  TString todayFile = Form("%02i_%02i_%4i", dt.GetDay(), dt.GetMonth(), dt.GetYear());
  
  // define canvas
  c1 = new TCanvas("c1", "Inv Mass",1100,800);
  c1->Clear();
  gPad->SetTopMargin(0.03);
  gPad->SetRightMargin(0.03);
  gPad->SetLeftMargin(0.13);
  c1->SetTickx(1);
  c1->SetTicky(1);
  c1->cd();
  
  hSpec->Draw("E");							  // draw points
  fitBG->Draw("same");						  // draw background
  fit1CB2->Draw("same");					  // draw function 1CB2
  if (fit2CB2!=NULL) fit2CB2->Draw("same");	  // draw function 2CB2
  fitFct->Draw("same");						  // draw global function

  // Nevent
  Double_t nbevt=0.;
  if (counter) {
	if (cent=="090") nbevt = eventCounters->GetSum("centrality:010,1020,2030,3040,4050,5060,6070,7080,8090");
	else if (cent=="020") nbevt = eventCounters->GetSum("centrality:010,1020");
	else if (cent=="2040") nbevt = eventCounters->GetSum("centrality:2030,3040");
	else if (cent=="4060") nbevt = eventCounters->GetSum("centrality:4050,5060");
	else if (cent=="6090") nbevt = eventCounters->GetSum("centrality:6070,7080,8090");
	else if (cent=="4090") nbevt = eventCounters->GetSum("centrality:4050,5060,6070,7080,8090");
	else  nbevt = eventCounters->GetSum(Form("centrality:%s",cent.Data()));
  }
						
  TPaveText* t1 = new TPaveText(0.69, 0.50, 0.97, 0.97,"NDC");
  t1->SetFillColor(0);
  t1->SetBorderSize(1);
  t1->SetTextColor(1);
  t1->SetTextFont(42);
  t1->SetTextSize(0.03);
  t1->SetTextAlign();
  t1->AddText(0.,0.95,Form("%s%s", centrality.Data(), cent.Data()));														// cent
  if (counter) t1->AddText(0.,0.,Form("%s%.0f", nbevents.Data(), nbevt));													// nevt
  TText *tnbJpsi = t1->AddText(0.,0.,Form("%s%.0f #pm %.0f", nbJpsi.Data(), val[0], val[1]));								// njpsi
  tnbJpsi->SetTextColor(kRed);
  t1->AddText(0.,0.,Form("%s%.3f #pm %.3f GeV/c^{2}", masse.Data(), val[2], val[3]));										//  mass
  t1->AddText(0.,0.,Form("%s%.0f #pm %.0f MeV/c^{2}", sigma.Data(), val[4]*1000, val[5]));
  t1->AddText(0.,0.,Form("%s%.2f #pm %.2f", SoverB.Data(), val[7]/val[9], val[11]));													//  S over B
  t1->AddText(0.,0.,Form("%s%.1f #pm %.1f", signal.Data(), val[7]/TMath::Sqrt(val[7]+val[9]), val[12]));								// significance
//  sigma
  if (fit2CB2!=NULL) {
	TText *tnbpsi = t1->AddText(0.,0.,Form("%s%.0f #pm %.0f", nbprim.Data(), val[13], val[14]));								
	tnbpsi->SetTextColor(kRed);
	t1->AddText(0.,0.,Form("%s%.3f #pm %.3f GeV/c^{2}", masse.Data(), val[15], val[16]));										
	t1->AddText(0.,0.,Form("%s%.0f #pm %.0f MeV/c^{2}", sigmaPrim.Data(), val[17]*1000, val[18]));
	t1->AddText(0.,0.,Form("%s%.2f #pm %.2f", SoverB.Data(), val[20]/val[22], val[24]));													//  S over B
	t1->AddText(0.,0.,Form("%s%.1f #pm %.1f", signal.Data(), val[20]/TMath::Sqrt(val[20]+val[22]), val[25]));								// significance
  }
  else {
	t1->AddText(0.,0.,Form("%s%.2f #pm %.2f", alphaCB.Data(), fitFct->GetParameter(nparBkg+3), fitFct->GetParError(nparBkg+3)));		// alpha
	t1->AddText(0.,0.,Form("%s%.1f #pm %.1f", nCB.Data(), fitFct->GetParameter(nparBkg+4), fitFct->GetParError(nparBkg+4) ));			// n
	t1->AddText(0.,0.,Form("%s%.2f #pm %.2f", alphaUpCB.Data(), fitFct->GetParameter(nparBkg+5), fitFct->GetParError(nparBkg+5)));		// alpha prim
	t1->AddText(0.,0.,Form("%s%.1f #pm %.1f", nUpCB.Data(), fitFct->GetParameter(nparBkg+6), fitFct->GetParError(nparBkg+6)));			// n prim
  }
  t1->AddText(0.,0.,Form("%s%.2f", chi2.Data(), val[6]));																				//  chi 2
  t1->Draw("same");
  
  TPaveText* t2=new TPaveText(0.15,0.15,0.35,0.35,"NDC");
  t2->SetFillStyle(0);
  t2->SetBorderSize(0);
  t2->AddText(0.,0.,"w/ physics selection");
  t2->AddText(0.,0.,"Trigger CPBI1MUL-B-NOPF-MUON");
  t2->AddText(0.,0.,"Muon cuts :");
  t2->AddText(0.,0.,"#eta, #theta_{abs}, Match Low trigger");
  t2->AddText(0.,0.,"Dimuon cut :");
  t2->AddText(0.,0.,"rapidity");
  t2->SetTextColor(kBlack);
  t2->SetTextFont(42);
  t2->SetTextSize(0.025);
  //t2->Draw("same");
  
  if (perf) {
	// ALICE logo
	TString type = "ALICE Performance";
	TString type2 = "Pb-Pb, #sqrt{s_{NN}} = 2.76 TeV";
	TString type3 = "2.5< y_{J/#psi} < 4";
	Double_t xPad = 0.3;
	Double_t yPad = 0.835;
	
	TVirtualPad* currPad = gPad;
	TPad *myPadLogo = new TPad("myPadLogo", "Pad for ALICE Logo",xPad,yPad,xPad+0.10,yPad+0.10);
	myPadLogo->SetFillColor(0);
	myPadLogo->SetBorderMode(0);
	myPadLogo->SetBorderSize(2);
	myPadLogo->SetFrameBorderMode(0);
	myPadLogo->SetLeftMargin(0.0);
	myPadLogo->SetTopMargin(0.0);
	myPadLogo->SetBottomMargin(0.0);
	myPadLogo->SetRightMargin(0.0);
	if (perf) myPadLogo->Draw("same");
	myPadLogo->cd();
	TASImage *myAliceLogo = new TASImage("alice_logo.png");
	if (perf) myAliceLogo->Draw("same");
	
	currPad->cd();
	Double_t x1 = xPad - 0.10, y1 = yPad - 0.06; //-0.07  -0.06
	Double_t x2 = x1 + 0.3, y2 = y1 + 0.08;   //0.08
	
	TPaveText* t3=new TPaveText(x1,y1,x2,y2,"NDC");
	t3->SetFillStyle(0);
	t3->SetBorderSize(0);
	t3->AddText(0.,0.,Form("%s", type.Data()));
	t3->SetTextColor(kRed);
	t3->SetTextFont(42);
	t3->SetTextSize(0.03);
	t3->Draw("same");
	
	TPaveText* t4=new TPaveText(x1+0.04,y1-0.04,x2-0.04,y2-0.04,"NDC"); // +0.06 - - -
	t4->SetFillStyle(0);
	t4->SetBorderSize(0);
	t4->SetTextColor(kBlack);
	t4->SetTextFont(52);
	t4->AddText(0.,0.,today.Data());
	t4->SetTextSize(0.03);
	t4->Draw("same");
	
	TPaveText* t5=new TPaveText(x1+0.08,y1-0.08,x2-0.08,y2-0.08,"NDC");
	t5->SetFillStyle(0);
	t5->SetBorderSize(0);
	t5->SetTextColor(kBlue);
	t5->SetTextFont(52);
	t5->AddText(0.,1.,Form("%s", type2.Data()));
	t5->SetTextSize(0.03);
	t5->Draw("same");
	
	TPaveText* t6=new TPaveText(x1+0.12,y1-0.12,x2-0.12,y2-0.12,"NDC");
	t6->SetFillStyle(0);
	t6->SetBorderSize(0);
	t6->SetTextColor(kBlue);
	t6->SetTextFont(52);
	t6->AddText(0.,1.,Form("%s", type3.Data()));
	t6->SetTextSize(0.03);
	t6->Draw("same");
	
	delete t3;
	delete t4;
	delete t5;
	delete t6;
  }
  
  
}
