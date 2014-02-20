/*
 *  WeightedAverageAndSys.C
 *  aliroot_dev
 *
 *  Created by philippe pillot on 27/06/12.
 *  Copyright 2012 SUBATECH. All rights reserved.
 *
 */

#if !defined(__CINT__) || defined(__MAKECINT__)

#include <Riostream.h>
#include "TMath.h"
#include "TString.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TList.h"

#endif

void WeightedAverageAndSys(Double_t *val, Double_t *err, Int_t nTests, Double_t accEff, TString data,
			   Double_t &mean, Double_t &stat, Double_t &sys);
void Ratio(Double_t valApT, Double_t errApT, Double_t valLpT, Double_t errLpT, Double_t &rat, Double_t &err);
Bool_t ExtractNJPsiVsPtY(Int_t pt, Int_t y, TString trig, Double_t *&val, Double_t *&err);
Double_t GetAccEffApT(Int_t pt, Int_t y);
Double_t GetAccEffLpT(Int_t pt, Int_t y);

//------------------------------------------------------------------------
void WeightedAverageAndSys(Int_t pt, Int_t y)
{
  /// compute the average of the NJPsi weighted by their statistical uncertainties
  /// as well as the associated statistical and systematic (dispersion) uncertainties
  
  // for ApT trigger
  Double_t *val, *err, meanApT, statApT, sysApT;
  TString trig = "Cent_AllPt";
  if (!ExtractNJPsiVsPtY(pt, y, trig, val, err)) return;
  Double_t accEffApT = GetAccEffApT(pt, y);
  WeightedAverageAndSys(val, err, 12, accEffApT, trig, meanApT, statApT, sysApT);
  
  // for LpT trigger
  Double_t meanLpT, statLpT, sysLpT;
  trig = "Cent_LPt";
  if (!ExtractNJPsiVsPtY(pt, y, trig, val, err)) return;
  Double_t accEffLpT = GetAccEffLpT(pt, y);
  WeightedAverageAndSys(val, err, 12, accEffLpT, trig, meanLpT, statLpT, sysLpT);
  
  // ratio
  Double_t ratio, stat, sys;
  Ratio(meanApT, statApT, meanLpT, statLpT, ratio, stat);
  sys = ratio * TMath::Sqrt(sysApT*sysApT/meanApT/meanApT + sysLpT*sysLpT/meanLpT/meanLpT);
  Double_t accEff = accEffLpT / accEffApT;
  printf("ratio: %f ± %f ± %f\n", ratio/accEff, stat/accEff, sys/accEff);
  
}

//------------------------------------------------------------------------
void WeightedRatio(Int_t pt, Int_t y)
{
  /// compute the average of LpT/ApT ratio weighted by their statistical uncertainties
  /// as well as the associated statistical and systematic (dispersion) uncertainties
  
  Double_t *valApT, *errApT;
  TString trig = "Cent_AllPt";
  if (!ExtractNJPsiVsPtY(pt, y, trig, valApT, errApT)) return;
  Double_t *valLpT, *errLpT;
  trig = "Cent_LPt";
  if (!ExtractNJPsiVsPtY(pt, y, trig, valLpT, errLpT)) return;
  const Int_t nTests = 12;
  
  // loop over tests
  Double_t rat[nTests], err[nTests];
  for (Int_t i=0; i<nTests; i++) {
    
    // compute ratio and associated statistical error
    Ratio(valApT[i], errApT[i], valLpT[i], errLpT[i], rat[i], err[i]);
    
  }
  
  // acc*eff correction
  Double_t accEff = GetAccEffLpT(pt, y) / GetAccEffApT(pt, y);
  
  // weighted ratio
  Double_t mean, stat, sys;
  WeightedAverageAndSys(rat, err, 12, accEff, "ratio", mean, stat, sys);
  
}

//------------------------------------------------------------------------
void WeightedAverageAndSys(Double_t *val, Double_t *err, Int_t nTests, Double_t accEff, TString data,
			   Double_t &mean, Double_t &stat, Double_t &sys)
{
  /// compute the average of the given values weighted by their statistical uncertainties
  /// as well as the associated statistical and systematic (dispersion) uncertainties
  /// all corrected for the Acc*Eff value
  
  // loop over tests
  Double_t wval = 0, wval2 = 0, w2err2 = 0, sumw = 0;
  for(Int_t i=0; i<nTests; i++) {
    
    // weight
    Double_t wstat = (data == "ratio") ? 1. : 1./val[i];
    Double_t w = 1./err[i]/err[i]/wstat;
    sumw += w;
    
    // mean
    wval += w*val[i];
    
    // stat
    w2err2 += w*w*err[i]*err[i];
    
    // rms
    wval2 += w*val[i]*val[i];
    
  }
  
  // results
  mean = wval/sumw;
  stat = TMath::Sqrt(w2err2*nTests)/sumw;
  sys = TMath::Sqrt(wval2/sumw - mean*mean);
  
  printf("%s: %f ± %f ± %f\n", data.Data(), mean/accEff, stat/accEff, sys/accEff);
  
}

//------------------------------------------------------------------------
void Ratio(Double_t valApT, Double_t errApT, Double_t valLpT, Double_t errLpT, Double_t &rat, Double_t &err)
{
  /// compute the average of LpT/ApT ratio and the associated statistical uncertainty
  
  // ratio
  rat = valLpT/valApT;
  
  // error
  Double_t sysA2 = TMath::Max(0., errApT*errApT - valApT);
  Double_t sysL2 = TMath::Max(0., errLpT*errLpT - valLpT);
  err = TMath::Sqrt(rat*(1.-rat)/valApT + sysA2/valApT/valApT + sysL2/valLpT/valLpT);
  
}

//------------------------------------------------------------------------
Bool_t ExtractNJPsiVsPtY(Int_t pt, Int_t y, TString trig, Double_t *&val, Double_t *&err)
{
  /// get the number to JPsi and statistical errors from the various fits in the given pT/y bin
  
  TString Study;
  if (y == 0 && pt == 0) Study="Cent_0pt8_4y25";
  else if (y == 0 && pt == 1) Study="Cent_0pt2_4y25";
  else if (y == 0 && pt == 2) Study="Cent_2pt5_4y25";
  else if (y == 0 && pt == 3) Study="Cent_5pt8_4y25";
  else if (y == 1 && pt == 0) Study="Cent_325y25_0pt8";
  else if (y == 2 && pt == 0) Study="Cent_4y325_0pt8";
  else {
    printf("wrong pT/y bin\n");
    return kFALSE;
  }
  
  TString dir = "/Users/pillot/Work/Alice/Work/Sim/LHC11h/LptOverApt";
  Int_t nTests=12, iTest = 0;
  val = new Double_t[nTests];
  err = new Double_t[nTests];
  
  // Signal shape
  TString SigFct = "1CB2";
  
  // BKG shape for mixing
  TString namebkg = "expo";
  
  // Different tails
  TList* Tails = new TList();
  Tails->Add(new TObjString("embedding2011"));
  Tails->Add(new TObjString("simuJpsi2011"));
  //########################################
  // loop over cases of tails
  //########################################
  TIter nextTails(Tails);
  TObjString* nameTails;
  while (( nameTails = static_cast<TObjString*>(nextTails()) ))
  {
    
    // Different Msigma
    TList* Msigma = new TList();
    Msigma->Add(new TObjString("0-90"));
    Msigma->Add(new TObjString("psigma"));
    Msigma->Add(new TObjString("msigma"));
    //########################################
    // loop over cases of Mass and sigma
    //########################################
    TIter nextMsigma(Msigma);
    TObjString* nameMsigma;
    while (( nameMsigma = static_cast<TObjString*>(nextMsigma()) ))
    {
      
      //###############  Fit
      // open file		
      TString file = Form( "%s/fit/%s/%s/FitResults_Cent_010.txt", dir.Data(), trig.Data(), Study.Data() );
      ifstream inFileRaw( file.Data() );
      TString currTestRaw;
      if (inFileRaw.is_open())
      {
	
	while (! inFileRaw.eof() )
	{
	  currTestRaw.ReadLine(inFileRaw,kTRUE);
	  if(currTestRaw.IsNull()) continue;
	  
	  TObjArray* arr = currTestRaw.Tokenize(" ");
	  
	  // Select Test
	  if (arr->At(0)->GetName() != SigFct) continue;
	  else if (arr->At(1)->GetName() != nameTails->String()) continue;
	  else if (arr->At(3)->GetName() != nameMsigma->String()) continue;
	  
	  // Get values
	  TString strnbJpsi = arr->At(5)->GetName();
	  TString strerrnbJpsi = arr->At(7)->GetName();
	  val[iTest] = strnbJpsi.Atof();
	  err[iTest] = strerrnbJpsi.Atof();
	  //printf("%f ± %f\n", val[iTest], err[iTest]);
	  break;
	  
	}  // loop in file
	
	inFileRaw.close();
	
      } else {
	printf("Input file %s unknow !!!!!! \n", file.Data());
	return kFALSE;
      }
      
      
      //###############  Mix
      // open file
      file = Form( "%s/mix/%s/%s/FitResults_Mixing_Cent_010.txt", dir.Data(), trig.Data(), Study.Data() );
      ifstream inFileMix( file.Data() );
      TString currTestMix;
      if (inFileMix.is_open())
      {
	
	while (! inFileMix.eof() )
	{
	  currTestMix.ReadLine(inFileMix,kTRUE);
	  if(currTestMix.IsNull()) continue;
	  
	  TObjArray* arr = currTestMix.Tokenize(" ");
	  
	  // Select Test
	  if (arr->At(0)->GetName() != SigFct) continue;
	  else if (arr->At(1)->GetName() != namebkg) continue;
	  else if (arr->At(2)->GetName() != nameTails->String()) continue;
	  else if (arr->At(4)->GetName() != nameMsigma->String()) continue;
	  
	  // Get values
	  TString strnbJpsi = arr->At(6)->GetName();
	  TString strerrnbJpsi = arr->At(8)->GetName();
	  val[iTest+6] = strnbJpsi.Atof();
	  err[iTest+6] = strerrnbJpsi.Atof();
	  //printf("%f ± %f\n", val[iTest+6], err[iTest+6]);
	  break;
	  
	}  // loop in file
	
	inFileMix.close();
	
      } else {
	printf("Input file %s unknow !!!!!! \n", file.Data());
	return kFALSE;
      }
      
      iTest++;
      
    } // end loop on Mas and sigma
    
  } // end loop on Tails
  
  return kTRUE;
  
}

//------------------------------------------------------------------------
Double_t GetAccEffApT(Int_t pt, Int_t y)
{
  /// return the Acc*Eff correction for the given pT/y bin
  
  if (y == 0 && pt == 0) return 0.1818;
  else if (y == 0 && pt == 1) return 0.1724;
  else if (y == 0 && pt == 2) return 0.1807;
  else if (y == 0 && pt == 3) return 0.2669;
  else if (y == 1 && pt == 0) return 0.1684;
  else if (y == 2 && pt == 0) return 0.1997;
  else {
    printf("wrong pT/y bin\n");
    return 1.;
  }
  
}

//------------------------------------------------------------------------
Double_t GetAccEffLpT(Int_t pt, Int_t y)
{
  /// return the Acc*Eff correction for the given pT/y bin
  
  if (y == 0 && pt == 0) return 0.1328;
  else if (y == 0 && pt == 1) return 0.1249;
  else if (y == 0 && pt == 2) return 0.1287;
  else if (y == 0 && pt == 3) return 0.2315;
  else if (y == 1 && pt == 0) return 0.1201;
  else if (y == 2 && pt == 0) return 0.1500;
  else {
    printf("wrong pT/y bin\n");
    return 1.;
  }
  
}

