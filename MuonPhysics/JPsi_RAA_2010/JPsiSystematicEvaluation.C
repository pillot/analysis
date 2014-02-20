/*
 *  JPsiSystematicEvaluation.C
 *  aliroot_dev
 *
 *  Created by philippe pillot on 05/05/11.
 *  Copyright 2011 SUBATECH. All rights reserved.
 *
 */

#if !defined(__CINT__) || defined(__MAKECINT__)

#include <Riostream.h>
#include "TCanvas.h"
#include "TH1.h"
#include "TMath.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TString.h"
#include "TASImage.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TROOT.h"

#endif

const Int_t nCentBins = 5;

TString gName[2*nCentBins-2] = {"080", "010", "1020", "2040", "4080", "010Over4080", "1020Over4080", "2040Over4080"};
TString gTitle[2*nCentBins-2] = {"0-80%", "0-10%", "10-20%", "20-40%", "40-80%", "0-10% / 40-80%", "10-20% / 40-80%", "20-40% / 40-80%"};

Double_t nMB[nCentBins]={16932459., 2088200., 2097283., 4258151., 8488825.}; // Christophe pass2
//Double_t nMB[nCentBins]={10097131., 1256765., 1260109., 2531675., 5048582.}; // Christophe pass1 reanalyzed
//Double_t nMB[nCentBins]={7.83e6, 0.98e6, 1.00e6, 1.97e6, 3.89e6}; // Christophe pass1
//Double_t nMB[nCentBins]={8.54e6, 1.06e6, 1.09e6, 2.15e6, 4.24e6};  // Roberta pass1
//Double_t nMB[nCentBins]={14820426., 1832821., 1834146., 3721847., 7431612.}; // + |vtx| < 10 cm

Double_t TAA[nCentBins]={7.04, 23.48, 14.43, 6.86, 1.20}; // mb^-1
Double_t eTAA[nCentBins]={0.27, 0.97, 0.57, 0.28, 0.07};

Double_t TAAOver4080[nCentBins-2]={19.49, 11.98, 5.69};
Double_t eTAAOver4080[nCentBins-2]={1.21, 0.63, 0.20};

//Double_t NPart[nCentBins] = {138.9, 356.5, 260.5, 157.3, 45.5};
//Double_t NPart[nCentBins] = {262.9, 359.5, 262.4, 166.9, 67.1}; // weighted by Gines
Double_t NPart[nCentBins] = {264.0, 360.7, 264.0, 168.2, 69.6}; // using Alberica's macro
Double_t eNPart[nCentBins] = {3.2, 3.6, 4.4, 3.4, 2.1};

// relative systematic errors
Double_t eGen[nCentBins+1] = {0., 0., 0., 0., 0., 0.03};
Double_t eTrkUncorr[nCentBins+1] = {0., 0., 0., 0., 0., 0.04};
Double_t eTrkCorr[nCentBins+1] = {0., 0., 0., 0., 0., 0.02};
Double_t eTrkRes[nCentBins+1] = {0., 0., 0., 0., 0., 0.02};
Double_t eTrkCent[nCentBins+1] = {0.02, 0.04, 0.02, 0.01, 0., 0.};
Double_t eTrg[nCentBins+1] = {0., 0., 0., 0., 0., 0.04};

Double_t AccEff=0.1944;
Double_t eAccEff=0.0005/AccEff;

// interpolated
//Double_t sigmaJPsipp = 4.1; //mub
//Double_t esigmaJPsipp=0.6/sigmaJPsipp;
// measured: sigma_J/psi (2.5<y<4) @ 2.76TeV=3.37 +- 0.13 (stat) +- 0.40 (syst) mub
//Double_t sigmaJPsipp = 3.37; //mub
//Double_t esigmaJPsipp=0.42/sigmaJPsipp;
// measured: sigma_J/psi (2.5<y<4) @ 2.76TeV=3.46 +- 0.13(stat) +- 0.32(syst) +- 0.28(syst. lumi) mub
Double_t sigmaJPsipp = 3.46; //mub
Double_t esigmaJPsipp=0.445/sigmaJPsipp;
Double_t BR=0.059;
Double_t eBR=0.01;

TString label[nCentBins] = {"0-80%", "0-10%", "10-20%", "20-40%", "40-80%"};
Double_t x[nCentBins] = {100., 5., 15., 30., 60.};
Double_t ex[nCentBins] = {1., 5., 5., 10., 20.};

Bool_t fineDisplay = kFALSE;

void ComputeRAA(Double_t *nJPsi, Double_t *enJPsi, Double_t *sysnJPsi, Double_t *renJPsi,
		Double_t *RAA, Double_t *eRAA, Double_t *eRAA2, Double_t &eRAACorr, Bool_t print = kFALSE);
void DisplayRAA(Double_t *renJPsi, Double_t *RAA, Double_t *eRAA, Double_t *eRAA2, Double_t eRAACorr,
		TGraphErrors *&gRAA, TGraphErrors *&gRAASys, Bool_t display = kFALSE, Bool_t vsCentClass = kTRUE);
void DisplayRAAWithPHENIX(Double_t *renJPsi, Double_t *RAA, Double_t *eRAA, Double_t *eRAA2, Double_t eRAACorr,
			  TGraphErrors *&gRAA, TGraphErrors *&gRAASys, Bool_t display = kFALSE, Bool_t vsCentClass = kTRUE);
void DisplayRAAWithPHENIX2(Double_t *renJPsi, Double_t *RAA, Double_t *eRAA, Double_t *eRAA2, Double_t eRAACorr,
			   TGraphErrors *&gRAA, TGraphErrors *&gRAASys, Bool_t display = kFALSE, Bool_t vsCentClass = kTRUE);
void DisplayRAAWithEPS09(Double_t *renJPsi, Double_t *RAA, Double_t *eRAA, Double_t *eRAA2, Double_t eRAACorr,
			 TGraphErrors *&gRAA, TGraphErrors *&gRAASys, Bool_t display, Bool_t vsCentClass);
void ComputeRCP(Double_t *nJPsiOver4080, Double_t *enJPsiOver4080, Double_t *sysnJPsiOver4080,
		Double_t *renJPsiOver4080, Double_t *RCP, Double_t *eRCP, Double_t *eRCP2, Bool_t print = kFALSE);
void DisplayRCP(Double_t *renJPsiOver4080, Double_t *RCP, Double_t *eRCP, Double_t *eRCP2,
		TGraphErrors *&gRCP, TGraphErrors *&gRCPSys, Bool_t display = kFALSE, Bool_t vsCentClass = kTRUE);
void DisplayRCPWithATLAS(Double_t *renJPsiOver4080, Double_t *RCP, Double_t *eRCP, Double_t *eRCP2,
			 TGraphErrors *&gRCP, TGraphErrors *&gRCPSys, Bool_t display, Bool_t vsCentClass);
void DisplayRCPWithATLAS2(Double_t *renJPsi, Double_t *RCP, Double_t *reRAA, Double_t *reRAA2,
			  TGraphErrors *&gRCP, TGraphErrors *&gRCPSys, Bool_t display, Bool_t vsCentClass);
void DisplayRCPWithATLASandEE(Double_t *renJPsiOver4080, Double_t *RCP, Double_t *eRCP, Double_t *eRCP2,
			      TGraphErrors *&gRCP, TGraphErrors *&gRCPSys, Bool_t display, Bool_t vsCentClass);
void ALICEseal(TString type, Double_t xPad, Double_t yPad);

//------------------------------------------------------------------------
void JPsiSystematicEvaluation(Bool_t weightedMeanRMS = kTRUE, Bool_t vsCentClass = kTRUE)
{
  
  Int_t font=42;
  gStyle->SetCanvasColor(10);
  gStyle->SetFrameFillColor(10);
  gStyle->SetStatColor(10);
  gStyle->SetFillColor(10);
  gStyle->SetTitleFillColor(10);
  gStyle->SetTextFont(font);
  gStyle->SetLabelFont(font,"x");
  gStyle->SetTitleFont(font,"x");
  gStyle->SetLabelFont(font,"y");
  gStyle->SetTitleFont(font,"y");
  gROOT->SetStyle("plain");
  gROOT->ForceStyle();
  
  const Int_t nTests = 27;
  
  Double_t nJPsi[nTests][nCentBins] = {
    
    // ------ CB1 ------
    
    // CB tail fixed to simulations (pure J/Ψ) in all cases (α=0.98 n=5.2) M and σ free
    // {2096, 809, 454, 663, 210}, // old values
//    {2162.96, 809.411, 454.392, 663.401, 212.865},
    // CB tail fixed to simulations (pure J/Ψ) in all cases (α=0.98 n=5.2) M and σ fixed to 0-80% fit
    // {2089, 862, 499, 540, 212}, // old values
    {2153.85, 876.35, 506.604, 551.718, 221.364},
    // CB tail fixed to simulations (pure J/Ψ) in all cases (α=0.98 n=5.2) M and σ fixed to 0-80% fit + 1 sigma
    {2273.25, 921.832, 529.306, 592.544, 233.296},
    // CB tail fixed to simulations (pure J/Ψ) in all cases (α=0.98 n=5.2) M and σ fixed to 0-80% fit - 1 sigma
    {2024.26, 827.281, 481.042, 507.819, 204.935},
    // CB tail fixed to simulations (pure J/Ψ) in all cases (α=0.98 n=5.2) M and σ fixed to 80-40% fit
    // {2022, 837, 488, 510, 210}, // old values
    // CB tail fixed to simulations (pure J/Ψ) in all cases (α=0.98 n=5.2) M and σ fixed to simu
    // {2044, 849, 495, 515, 214}, // old values
    
    // CB tail fixed to pp LHC10g in all cases (α=1.15 n=1.59) M and σ free
    // {2298, 894, 499, 724, 221}, // old values
//    {2390.11, 893.837, 498.857, 723.754, 221.304},
    // CB tail fixed to pp LHC10g in all cases (α=1.15 n=1.59) M and σ fix to 0-80%
    // {2291, 949, 550, 595, 234}, // old values
    {2391.67, 971.503, 562.676, 613.426, 240.272},
    // CB tail fixed to pp LHC10g in all cases (α=1.15 n=1.59) M and σ fix to 0-80% + 1 sigma
    {2509.3, 1016.26, 584.146, 654.633, 250.984},
    // CB tail fixed to pp LHC10g in all cases (α=1.15 n=1.59) M and σ fix to 0-80% - 1 sigma
    {2261.47, 922.325, 537.769, 568.482, 229.134},
    
    // CB tail fixed to pp “old” MC in all cases (α=1.15 n=3.6) M and σ free
    // {2090, 808, 447, 656, 201}, // old values
//     {2150.67, 807.882, 446.877, 656.033, 201.033},
    // CB tail fixed to pp “old” MC in all cases (α=1.15 n=3.6) M and σ fix to 0-80%
    // {2097, 864, 498, 544, 214}, // old values
     {2143.07, 871.758, 502.391, 549.856, 215.855},
    // CB tail fixed to pp “old” MC in all cases (α=1.15 n=3.6) M and σ fix to 0-80% + 1 sigma
     {2257.42, 915.522, 523.867, 589.171, 226.088},
    // CB tail fixed to pp “old” MC in all cases (α=1.15 n=3.6) M and σ fix to 0-80% - 1 sigma
     {2017.65, 824.049, 477.891, 507.23, 205.409},
     
    // ------ CB2 ------
     
    // CB2 tails fixed to embedding (pass1, Javier) in all cases M and σ free
    // {2159, 835, 481, 722, 217}, // old values
//    {2229.59, 834.767, 481.277, 721.493, 216.88},
    // CB2 tails fixed to embedding (pass1, Javier) in all cases M and σ fixed to 0-80% fit
    // {2164, 910, 527, 575, 218}, // old values
    {2234.58, 926.819, 536.676, 587.362, 223.299},
    // CB2 tails fixed to embedding (pass1, Javier) in all cases M and σ fixed to 0-80% fit + 1 sigma
    {2378.63, 982.654, 565.203, 636.732, 242.545},
    // CB2 tails fixed to embedding (pass1, Javier) in all cases M and σ fixed to 0-80% fit - 1 sigma
    {2078.83, 866.898, 504.08, 534.939, 209.316},
    // CB2 tails fixed to embedding (pass1, Javier) in all cases M and σ fixed to 80-40% fit
    // {2095, 885, 515, 544, 216}, // old values
    // CB2 tails fixed to embedding (pass1, Javier) in all cases M and σ fixed to 0-80% values from embedding
    // {2154, 913, 531, 563, 227}, // old values
    // CB2 tails fixed to embedding (pass1, Javier) in all cases σ fixed to 0-80% values from embedding
    // {2201, 928, 539, 586, 218}, // old values
    
    // CB2 tails fixed to pure J/Ψ simulation (pass2, Phillipe) M and σ free
    // {2142, 822, 466, 703, 207}, // old values
//    {2210.21, 821.691, 466.224, 703.097, 206.872},
    // CB2 tails fixed to pure J/Ψ simulation (pass2, Phillipe) M and σ fixed to 0-80% fit
    // {2134, 881, 508, 551, 219}, // old values
    {2202.92, 897.858, 517.553, 562.783, 223.775},
    // CB2 tails fixed to pure J/Ψ simulation (pass2, Phillipe) M and σ fixed to 0-80% fit + 1 sigma
    {2347.5, 952.58, 546.047, 611.474, 243.159},
    // CB2 tails fixed to pure J/Ψ simulation (pass2, Phillipe) M and σ fixed to 0-80% fit - 1 sigma
    {2046.76, 838.765, 485.59, 511.051, 209.373},
    // CB2 tails fixed to pure J/Ψ simulation (pass2, Phillipe) M and σ fixed to simu
    // {2110, 877, 509, 534, 223}, // old values
    
    // ------ mixing + CB1 ------
    
    // Background from mixing, CB (α=0.98  n=5.2) + 1 expo. + range [2;4.4] M and σ free
//    {2391.16, 1007.57, 499.691, 672.942, 226.851},
    // Background from mixing, CB (α=0.98  n=5.2) + 1 expo. + range [2;4.4] M and σ fixed to 0-80% fit
    {2395.34, 1017.48, 553.493, 610.733, 234.824},
    // Background from mixing, CB (α=0.98  n=5.2) + 1 expo. + range [2;4.4] M and σ fixed to 0-80% fit + 1 sigma
    {2501.83, 1063.2, 573.14, 644.14, 242.544},
    // Background from mixing, CB (α=0.98  n=5.2) + 1 expo. + range [2;4.4] M and σ fixed to 0-80% fit - 1 sigma
    {2277.64, 967.858, 530.947, 573.82, 225.794},
    
    // Background from mixing, CB (α=1.15  n=1.59) + 1 expo. + range [2;4.4] M and σ free
//    {2611.29, 1100.19, 543.469, 729.544, 252.708},
    // Background from mixing, CB (α=1.15  n=1.59) + 1 expo. + range [2;4.4] M and σ fixed to 0-80% fit
    {2606.13, 1106.41, 603.094, 661.453, 256.887},
    // Background from mixing, CB (α=1.15  n=1.59) + 1 expo. + range [2;4.4] M and σ fixed to 0-80% fit + 1 sigma
    {2710.69, 1150.83, 621.452, 694.894, 262.949},
    // Background from mixing, CB (α=1.15  n=1.59) + 1 expo. + range [2;4.4] M and σ fixed to 0-80% fit - 1 sigma
    {2487.58, 1056.85, 581.125, 623.741, 247.985},
    
    // Background from mixing, CB (α=0.98  n=5.2) + pol1 + range [2;4.4] M and σ free
//    {2308.37, 965.144, 483.555, 655.327, 196.402},
    // Background from mixing, CB (α=0.98  n=5.2) + pol1 + range [2;4.4] M and σ fixed to 0-80% fit
    {2300.89, 978.38, 536.43, 579.709, 212.235},
    // Background from mixing, CB (α=0.98  n=5.2) + pol1 + range [2;4.4] M and σ fixed to 0-80% fit + 1 sigma
    {2417.51, 1027.93, 558.608, 616.223, 219.132},
    // Background from mixing, CB (α=0.98  n=5.2) + pol1 + range [2;4.4] M and σ fixed to 0-80% fit - 1 sigma
    {2173.66, 925.174, 511.39, 540.095, 204.183},
    
    // Background from mixing, CB (α=1.15  n=1.59) + pol1 + range [2;4.4] M and σ free
//    {2533.7, 1061.17, 529.945, 714.527, 218.2},
    // Background from mixing, CB (α=1.15  n=1.59) + pol1 + range [2;4.4] M and σ fixed to 0-80% fit
    {2526.57, 1073.77, 590.995, 635.611, 232.116},
    // Background from mixing, CB (α=1.15  n=1.59) + pol1 + range [2;4.4] M and σ fixed to 0-80% fit + 1 sigma
    {2640.63, 1122.28, 611.893, 672.176, 238.669},
    // Background from mixing, CB (α=1.15  n=1.59) + pol1 + range [2;4.4] M and σ fixed to 0-80% fit - 1 sigma
    {2399.27, 1020.62, 566.599, 595.2, 224.137}
    
    // ------ LS + CB1 ------
    /*
    // Background substracted via likesign, CB (α=0.98 n=5.2) + pol(1) Norm[1.7,2.7] M and σ free
//    {2555, 1056, 536, 723, 158},
    // Background substracted via likesign, CB (α=0.98 n=5.2) + pol(1) Norm[1.7,2.7] M and σ fixed to 0-80% fit
    {2557, 1047, 577, 716, 198},
    // Background substracted via likesign, CB (α=0.98 n=5.2) + pol(0) Norm[1.7,2.7] M and σ free
//    {2599, 1037, 608, 717, 158},
    // Background substracted via likesign, CB (α=0.98 n=5.2) + pol(0) Norm[1.7,2.7] M and σ fixed to 0-80% fit
    {2596, 1043, 612, 719, 194},
    
    // Background substracted via likesign, CB (α=1.15 n=1.59) + pol(1) Norm[1.7,2.7] M and σ free
//    {2716, 1120, 534, 767, 168},
    // Background substracted via likesign, CB (α=1.15 n=1.59) + pol(1) Norm[1.7,2.7] M and σ fixed to 0-80% fit
    {2714, 1107, 608, 760, 208},
    // Background substracted via likesign, CB (α=1.15 n=1.59) + pol(0) Norm[1.7,2.7] M and σ free
//    {2736, 1101, 620, 760, 168},
    // Background substracted via likesign, CB (α=1.15 n=1.59) + pol(0) Norm[1.7,2.7] M and σ fixed to 0-80% fit
    {2739, 1102, 643, 759, 205}
    */
    // ------ tests ------
    
    // test pDCA
    // CB tail fixed to simulations (pure J/Ψ) in all cases (α=0.98  n=5.2) M and σ free + pDCA cut on both muons
    //    {2114.42, 779.412, 447.144, 668.686, 200.728},
    // CB tail fixed to simulations (pure J/Ψ) in all cases (α=0.98  n=5.2) M and σ fixed to 0-80% fit + pDCA cut on both muons
//    {2112.19, 856.625, 493.28, 547.333, 215.66}
    
    // test vertex
    // CB tail fixed to simulations (pure J/Ψ) in all cases (α=0.98  n=5.2) M and σ free + |vtx| < 10 cm
//    {1890.34, 720.969, 399.052, 563.088, 170.047},
    // CB tail fixed to simulations (pure J/Ψ) in all cases (α=0.98  n=5.2) M and σ fixed to 0-80% fit + |vtx| < 10 cm
//    {1885.65, 780.415, 442.798, 463.09, 185.024}
    
    // ------ comparison pass1 / pass2
    /*
    // Christophe's results
    //CB tail fixed to pp LHC10g in all cases (α=1.15 n=1.59 )
    //Double_t NjpsiC[5]={961, 344, 229, 263, 124}; 
    //Double_t e_NjpsiC[5]={127, 142, 94, 57, 18};
    //CB tail fixed to pp LHC10g in all cases (α=1.15 n=1.59 ) fix to 0-80%
    {957, 363, 239, 231, 124},
    //CB tail fixed to pp “old” MC in all cases (α=1.15 n=3.6 )
    //Double_t NjpsiC[5]={1049, 381, 250, 280, 140}; 
    //Double_t e_NjpsiC[5]={228, 150, 108, 59, 20};
    //CB tail fixed to pp “old” MC in all cases (α=1.15 n=3.6 ) fix to 0-80%
    {1046, 399, 259, 252, 136}
    
    // Roberta's results
    //CB tail fixed to pp LHC10g in all cases (α=1.15 n=1.59 ) fix to mass=3.132 and sigma=71 (Christophe's fit in 0-80%)
    {0., 518, 232, 259, 144},
    //CB tail fixed to pp “old” MC in all cases (α=1.15 n=3.6 ) fix to mass=3.134 and sigma=68 (Christophe's fit in 0-80%)
    {0., 440, 201, 221, 124}
    */
    // ------ pass1 reanalyzed with v4-21-20-AN
    /*
    // CB tail fixed to simulations (pure J/Ψ) in all cases  (α=0.98  n=5.2) M and σ free
//    {1271.78, 520.431, 238.685, 382.498, 137.197},
    // CB tail fixed to simulations (pure J/Ψ) in all cases  (α=0.98  n=5.2) M and σ fixed to 0-80%
//    {1270.85, 538.437, 251.891, 332.248, 144.243},
    // CB tail fixed to pp  LHC10g in all cases  (α=1.15  n=1.59) M and σ free
//    {1398.1, 557.292, 266.841, 418.198, 152.542},
    // CB tail fixed to pp  LHC10g in all cases  (α=1.15  n=1.59) M and σ fixed to 0-80%
    {1402.52, 594.752, 277.609, 365.538, 159.352},
    // CB tail fixed to pp  LHC10g in all cases  (α=1.15  n=3.6) M and σ free
//    {1263.9, 513.084, 238.167, 384.334, 136.594},
    // CB tail fixed to pp  LHC10g in all cases  (α=1.15  n=3.6) M and σ fixed to 0-80%
    {1262.8, 532.65, 252.039, 329.889, 144.012}
    */
  };
  
  Double_t enJPsi[nTests][nCentBins] = {
    
    // ------ CB1 ------
    
    // CB tail fixed to simulations (pure J/Ψ) in all cases (α=0.98 n=5.2) M and σ free
    // {298, 128, 79, 78, 58}, // old values
//    {183.501, 129.298, 141.612, 121.683, 60.1488},
    // CB tail fixed to simulations (pure J/Ψ) in all cases (α=0.98 n=5.2) M and σ fixed to 0-80% fit
    // {125, 97, 63, 52, 27}, // old values
    {129.275, 89.6631, 64.5531, 54.2925, 23.633},
    // CB tail fixed to simulations (pure J/Ψ) in all cases (α=0.98 n=5.2) M and σ fixed to 0-80% fit + 1 sigma
    {136.626, 102.23, 68.2121, 58.1256, 25.5611},
    // CB tail fixed to simulations (pure J/Ψ) in all cases (α=0.98 n=5.2) M and σ fixed to 0-80% fit - 1 sigma
    {121.845, 89.8367, 59.7223, 49.9765, 23.7522},
    // CB tail fixed to simulations (pure J/Ψ) in all cases (α=0.98 n=5.2) M and σ fixed to 80-40% fit
    // {122, 91, 62, 52, 25}, // old values
    // CB tail fixed to simulations (pure J/Ψ) in all cases (α=0.98 n=5.2) M and σ fixed to simu
    // {127, 95, 64, 57, 25}, // old values
    
    // CB tail fixed to pp LHC10g in all cases (α=1.15 n=1.59) M and σ free
    // {331, 146, 84, 85, 62}, // old values
//    {332.572, 256.135, 81.0814, 90.6398, 63.0544},
    // CB tail fixed to pp LHC10g in all cases (α=1.15 n=1.59) M and σ fix to 0-80%
    // {138, 104, 70, 61, 27}, // old values
    {143.942, 106.187, 71.861, 60.7605, 28.3302},
    // CB tail fixed to pp LHC10g in all cases (α=1.15 n=1.59) M and σ fix to 0-80% + 1 sigma
    {151.269, 106.719, 75.5043, 63.6141, 27.1848},
    // CB tail fixed to pp LHC10g in all cases (α=1.15 n=1.59) M and σ fix to 0-80% - 1 sigma
    {136.422, 100.55, 68.1512, 55.6069, 26.0801},
    
    // CB tail fixed to pp “old” MC in all cases (α=1.15 n=3.6) M and σ free
    // {301, 133, 140, 121, 56}, // old values
//    {181.75, 241.676, 78.9437, 120.51, 56.2427},
    // CB tail fixed to pp “old” MC in all cases (α=1.15 n=3.6) M and σ fix to 0-80%
    // {126, 94, 63, 55, 24}, // old values
    {128.725, 95.002, 64.2779, 55.6465, 23.6216},
    // CB tail fixed to pp “old” MC in all cases (α=1.15 n=3.6) M and σ fix to 0-80% + 1 sigma
    {135.768, 100.284, 67.783, 53.1789, 24.548},
    // CB tail fixed to pp “old” MC in all cases (α=1.15 n=3.6) M and σ fix to 0-80% - 1 sigma
    {121.547, 89.6193, 60.7035, 52.4904, 23.4889},
    
    // ------ CB2 ------
    
    // CB2 tails fixed to embedding (pass1, Javier) in all cases M and σ free
    // {330, 127, 145, 85, 62}, // old values
//    {204.883, 220.135, 151.647, 91.0779, 62.6237},
    // CB2 tails fixed to embedding (pass1, Javier) in all cases M and σ fixed to 0-80% fit
    // {130, 99, 67, 56, 29}, // old values
    {134.041, 100.705, 68.2956, 58.6176, 26.5908},
    // CB2 tails fixed to embedding (pass1, Javier) in all cases M and σ fixed to 0-80% fit + 1 sigma
    {142.916, 107.581, 72.8167, 60.29, 26.3148},
    // CB2 tails fixed to embedding (pass1, Javier) in all cases M and σ fixed to 0-80% fit - 1 sigma
    {125.074, 93.809, 64.6385, 56.0196, 25.6041},
    // CB2 tails fixed to embedding (pass1, Javier) in all cases M and σ fixed to 80-40% fit
    // {126, 96, 65, 57, 26}, // old values
    // CB2 tails fixed to embedding (pass1, Javier) in all cases M and σ fixed to 0-80% values from embedding
    // {135, 104, 70, 61, 25}, // old values
    // CB2 tails fixed to embedding (pass1, Javier) in all cases σ fixed to 0-80% values from embedding
    // {132, 101, 69, 58, 25}, // old values
    
    // CB2 tails fixed to pure J/Ψ simulation (pass2, Phillipe), no centrality dependence M and σ free
    // {324, 128, 149, 92, 59}, // old values
//    {198.079, 224.597, 146.741, 84.2902, 58.4129},
    // CB2 tails fixed to pure J/Ψ simulation (pass2, Phillipe), no centrality dependence M and σ fixed to 0-80% fit
    // {128, 96, 64, 53, 24}, // old values
    {132.122, 97.5311, 64.93, 56.7741, 27.5679},
    // CB2 tails fixed to pure J/Ψ simulation (pass2, Phillipe) M and σ fixed to 0-80% fit + 1 sigma
    {140.98, 104.177, 70.3913, 59.656, 26.361},
    // CB2 tails fixed to pure J/Ψ simulation (pass2, Phillipe) M and σ fixed to 0-80% fit - 1 sigma
    {123.2, 90.8416, 61.5515, 54.8211, 25.74},
    // CB2 tails fixed to pure J/Ψ simulation (pass2, Phillipe), no centrality dependence M and σ fixed to simu
    // {128, 97, 66, 58, 30}, // old values
    
    // ------ mixing + CB1 ------
    
    // Background from mixing, CB (α=0.98  n=5.2) + 1 expo. + range [2;4.4] M and σ free
//    {176.527, 136.822, 78.1132, 73.4777, 34.031},
    // Background from mixing, CB (α=0.98  n=5.2) + 1 expo. + range [2;4.4] M and σ fixed to 0-80% fit
    {127.806, 94.7901, 63.6377, 52.5848, 24.7286},
    // Background from mixing, CB (α=0.98  n=5.2) + 1 expo. + range [2;4.4] M and σ fixed to 0-80% fit + 1 sigma
    {133.934, 99.7752, 66.5451, 55.0206, 24.9866},
    // Background from mixing, CB (α=0.98  n=5.2) + 1 expo. + range [2;4.4] M and σ fixed to 0-80% fit - 1 sigma
    {121.827, 89.8983, 60.635, 50.0617, 24.717},
    
    // Background from mixing, CB (α=1.15  n=1.59) + 1 expo. + range [2;4.4] M and σ free
//    {185.348, 148.772, 80.2905, 77.6998, 31.9619},
    // Background from mixing, CB (α=1.15  n=1.59) + 1 expo. + range [2;4.4] M and σ fixed to 0-80% fit
    {139.625, 103.536, 69.1915, 57.3357, 26.0085},
    // Background from mixing, CB (α=1.15  n=1.59) + 1 expo. + range [2;4.4] M and σ fixed to 0-80% fit + 1 sigma
    {145.882, 108.581, 71.9887, 59.8247, 25.7328},
    // Background from mixing, CB (α=1.15  n=1.59) + 1 expo. + range [2;4.4] M and σ fixed to 0-80% fit - 1 sigma
    {133.247, 98.5327, 66.2313, 54.7693, 27.0107},
    
    // Background from mixing, CB (α=0.98  n=5.2) + pol1 + range [2;4.4] M and σ free
//    {185.626, 148.194, 81.187, 78.7773, 27.7218},
    // Background from mixing, CB (α=0.98  n=5.2) + pol1 + range [2;4.4] M and σ fixed to 0-80% fit
    {130.026, 95.8811, 64.9775, 53.5687, 24.5214},
    // Background from mixing, CB (α=0.98  n=5.2) + pol1 + range [2;4.4] M and σ fixed to 0-80% fit + 1 sigma
    {136.826, 100.981, 68.3486, 56.3233, 25.5404},
    // Background from mixing, CB (α=0.98  n=5.2) + pol1 + range [2;4.4] M and σ fixed to 0-80% fit - 1 sigma
    {123.093, 90.6531, 61.54, 50.7222, 23.3847},
    
    // Background from mixing, CB (α=1.15  n=1.59) + pol1 + range [2;4.4] M and σ free
//    {196.104, 158.551, 85.9328, 82.3399, 30.5106},
    // Background from mixing, CB (α=1.15  n=1.59) + pol1 + range [2;4.4] M and σ fixed to 0-80% fit
    {142.977, 105.45, 71.4385, 58.9279, 26.9793},
    // Background from mixing, CB (α=1.15  n=1.59) + pol1 + range [2;4.4] M and σ fixed to 0-80% fit + 1 sigma
    {149.684, 110.509, 74.76, 61.5461, 27.9855},
    // Background from mixing, CB (α=1.15  n=1.59) + pol1 + range [2;4.4] M and σ fixed to 0-80% fit - 1 sigma
    {136.019, 100.235, 67.9918, 56.0315, 25.849}
    
    // ------ LS + CB1 ------
    /*
    // Background substracted via likesign, CB (α=0.98 n=5.2) + pol(1) Norm[1.7,2.7] M and σ free
//    {258, 417, 132, 84, 24},
    // Background substracted via likesign, CB (α=0.98 n=5.2) + pol(1) Norm[1.7,2.7] M and σ fixed to 0-80% fit
    {186, 127, 93, 74, 31},
    // Background substracted via likesign, CB (α=0.98 n=5.2) + pol(0) Norm[1.7,2.7] M and σ free
//    {220, 402, 143, 84, 24},
    // Background substracted via likesign, CB (α=0.98 n=5.2) + pol(0) Norm[1.7,2.7] M and σ fixed to 0-80% fit
    {172, 127, 86, 69, 31},
    
    // Background substracted via likesign, CB (α=1.15 n=1.59) + pol(1) Norm[1.7,2.7] M and σ free
//    {270, 175, 122, 86, 30},
    // Background substracted via likesign, CB (α=1.15 n=1.59) + pol(1) Norm[1.7,2.7] M and σ fixed to 0-80% fit
    {200, 134, 99, 79, 32},
    // Background substracted via likesign, CB (α=1.15 n=1.59) + pol(0) Norm[1.7,2.7] M and σ free
//    {220, 171, 173, 85, 168},
    // Background substracted via likesign, CB (α=1.15 n=1.59) + pol(0) Norm[1.7,2.7] M and σ fixed to 0-80% fit
    {180, 134, 91, 73, 33}
     */
    // ------ tests ------
    
    // test pDCA
    // CB tail fixed to simulations (pure J/Ψ) in all cases  (α=0.98  n=5.2) M and σ free + pDCA cut on both muons
//    {286.778, 210.986, 139.563, 121.983, 57.3914},
    // CB tail fixed to simulations (pure J/Ψ) in all cases  (α=0.98  n=5.2) M and σ fixed to 0-80% fit + pDCA cut on both muons
//    {126.615, 92.463, 63.4685, 53.5879, 25.1315}
    
    // test vertex
    // CB tail fixed to simulations (pure J/Ψ) in all cases (α=0.98  n=5.2) M and σ free + |vtx| < 10 cm
//    {262.878, 115.814, 74.8944, 77.2995, 26.0748},
    // CB tail fixed to simulations (pure J/Ψ) in all cases (α=0.98  n=5.2) M and σ fixed to 0-80% fit + |vtx| < 10 cm
//    {121.871, 91.4927, 61.7258, 57.5254, 24.4733}    
    
    // ------ comparison pass1 / pass2
    /*
    // Christophe's results
    //CB tail fixed to pp LHC10g in all cases (α=1.15 n=1.59 )
    //Double_t NjpsiC[5]={961, 344, 229, 263, 124}; 
    //Double_t e_NjpsiC[5]={127, 142, 94, 57, 18};
    //CB tail fixed to pp LHC10g in all cases (α=1.15 n=1.59 ) fix to 0-80%
    {85, 63, 42, 34, 16},
    //CB tail fixed to pp “old” MC in all cases (α=1.15 n=3.6 )
    //Double_t NjpsiC[5]={1049, 381, 250, 280, 140}; 
    //Double_t e_NjpsiC[5]={228, 150, 108, 59, 20};
    //CB tail fixed to pp “old” MC in all cases (α=1.15 n=3.6 ) fix to 0-80%
    {93, 69, 46, 37, 18}
    
    // Roberta's results
    //CB tail fixed to pp LHC10g in all cases (α=1.15 n=1.59 ) fix to mass=3.132 and sigma=71 (Christophe's fit in 0-80%)
    {0., 73, 50, 40, 19},
    //CB tail fixed to pp “old” MC in all cases (α=1.15 n=3.6 ) fix to mass=3.134 and sigma=68 (Christophe's fit in 0-80%)
    {0., 62, 43, 35, 17}
    */
    // ------ pass1 reanalyzed with v4-21-20-AN
    /*
    // CB tail fixed to simulations (pure J/Ψ) in all cases  (α=0.98  n=5.2) M and σ free
//    {133.717, 102.784, 52.399, 64.7477, 40.01},
    // CB tail fixed to simulations (pure J/Ψ) in all cases  (α=0.98  n=5.2) M and σ fixed to 0-80%
//    {97.0021, 70.0608, 47.8053, 39.9424, 18.5321},
    // CB tail fixed to pp  LHC10g in all cases  (α=1.15  n=1.59) M and σ free
//    {143.566, 99.3737, 57.1956, 123.29, 21.8205},
    // CB tail fixed to pp  LHC10g in all cases  (α=1.15  n=1.59) M and σ fixed to 0-80%
    {107.852, 79.9492, 52.9088, 44.4149, 20.4593},
    // CB tail fixed to pp  LHC10g in all cases  (α=1.15  n=3.6) M and σ free
//    {135.846, 192.302, 77.8501, 106.525, 39.3303},
    // CB tail fixed to pp  LHC10g in all cases  (α=1.15  n=3.6) M and σ fixed to 0-80%
    {92.0026, 69.0733, 47.6538, 39.77, 18.4855}
    */
  };
    
  // define histograms
  TGraphErrors *gNJpsi[nCentBins];
  TH1F *hNJpsi[nCentBins];
  TH1F *heNJpsi[nCentBins];
  for(Int_t i=0; i<nCentBins; i++) {
    gNJpsi[i] = new TGraphErrors(nTests);
    gNJpsi[i]->SetName(Form("g%s",gName[i].Data()));
    gNJpsi[i]->SetTitle(Form("g%s",gTitle[i].Data()));
    hNJpsi[i] = new TH1F(Form("h%s",gName[i].Data()), Form("h%s",gTitle[i].Data()), 3000, 0., 3000.);
    heNJpsi[i] = new TH1F(Form("he%s",gName[i].Data()), Form("he%s",gTitle[i].Data()), 500, 0., 500.);
  }
  TGraphErrors *gNJpsiOver4080[nCentBins-2];
  TH1F *hNJpsiOver4080[nCentBins-2];
  TH1F *heNJpsiOver4080[nCentBins-2];
  for(Int_t i=0; i<nCentBins-2; i++) {
    gNJpsiOver4080[i] = new TGraphErrors(nTests);
    gNJpsiOver4080[i]->SetName(Form("g%s",gName[nCentBins+i].Data()));
    gNJpsiOver4080[i]->SetTitle(Form("g%s",gTitle[nCentBins+i].Data()));
    hNJpsiOver4080[i] = new TH1F(Form("h%s",gName[nCentBins+i].Data()), Form("h%s",gTitle[nCentBins+i].Data()), 1000, 0., 10.);
    heNJpsiOver4080[i] = new TH1F(Form("he%s",gName[nCentBins+i].Data()), Form("he%s",gTitle[nCentBins+i].Data()), 500, 0., 5.);
  }
  
  // loop over tests
  Double_t enJPsiMean2[nCentBins] = {0., 0., 0., 0., 0.};
  Double_t enJPsiOver4080Mean2[nCentBins-2] = {0., 0., 0.};
  Double_t noSys[nCentBins] = {0., 0., 0., 0., 0.};
  Double_t eRAACorr;
  Double_t renJPsi[nCentBins], RAA[nCentBins], eRAA[nCentBins], eRAA2[nCentBins];
  TGraphErrors *g = 0x0, *gSys = 0x0;
  TCanvas *c = new TCanvas("compareRAA", "compareRAA");
  gPad->SetFrameBorderMode(0);
  for(Int_t i=0; i<nTests; i++) {
    
    // display current RAA
    ComputeRAA(nJPsi[i], enJPsi[i], noSys, renJPsi, RAA, eRAA, eRAA2, eRAACorr);
    DisplayRAA(renJPsi, RAA, eRAA, eRAA2, eRAACorr, g, gSys);
    g->SetMarkerStyle(21);
    g->SetMarkerColor(i+1);
    g->SetLineColor(i+1);
    c->cd();
    if (i==0) {
      g->GetYaxis()->SetRangeUser(0., 1.2);
      g->Draw("ap");
    }
    else g->Draw("p");
    
    // loop over centrality classes
    for(Int_t j=0; j<nCentBins; j++) {
      
      Double_t wstat = 1./nJPsi[i][j];
      Double_t w = weightedMeanRMS ? 1./enJPsi[i][j]/enJPsi[i][j]/wstat : 1;
      enJPsiMean2[j] += w*w*enJPsi[i][j]*enJPsi[i][j];
      
      // fill histogram
      gNJpsi[j]->SetPoint(i, i, nJPsi[i][j]);
      gNJpsi[j]->SetPointError(i, 0, enJPsi[i][j]);
      hNJpsi[j]->Fill(nJPsi[i][j], w);
      heNJpsi[j]->Fill(enJPsi[i][j], w);
      
    }
    
    // loop over centrality classes
    for(Int_t j=0; j<nCentBins-2; j++) {
      
      Double_t r = nJPsi[i][j+1]/nJPsi[i][nCentBins-1];
      Double_t e = r * TMath::Sqrt(enJPsi[i][j+1]*enJPsi[i][j+1]/nJPsi[i][j+1]/nJPsi[i][j+1] + enJPsi[i][nCentBins-1]*enJPsi[i][nCentBins-1]/nJPsi[i][nCentBins-1]/nJPsi[i][nCentBins-1]);
      Double_t wstat = 1. / r / r / (1./nJPsi[i][j+1] + 1./nJPsi[i][nCentBins-1]);
      Double_t w = weightedMeanRMS ? 1./e/e/wstat : 1;
      enJPsiOver4080Mean2[j] += w*w*e*e;
      
      // fill histogram
      gNJpsiOver4080[j]->SetPoint(i, i, r);
      gNJpsiOver4080[j]->SetPointError(i, 0, e);
      hNJpsiOver4080[j]->Fill(r, w);
      heNJpsiOver4080[j]->Fill(e, w);
      
    }
    
  }
  
  // compute NJpsi means and errors and display histograms
  Double_t nJPsiMean[nCentBins];
  Double_t enJPsiMean[nCentBins];
  Double_t sysnJPsiMean[nCentBins];
  TCanvas *sys1;
  if (fineDisplay) sys1 = new TCanvas("sys1", "systematics of NJpsi", 1200, 900);
  else sys1 = new TCanvas("sys1", "systematics of NJpsi", 800, 900);
  gPad->SetFrameBorderMode(0);
  if (fineDisplay) sys1->Divide(3,nCentBins);
  else sys1->Divide(1,nCentBins);
  for(Int_t i=0; i<nCentBins; i++) {
    
    // compute NJpsi means and errors
    nJPsiMean[i] = hNJpsi[i]->GetMean();
    sysnJPsiMean[i] = hNJpsi[i]->GetRMS();
    //enJPsiMean[i] = heNJpsi[i]->GetMean();
    enJPsiMean[i] = TMath::Sqrt(enJPsiMean2[i]*nTests)/hNJpsi[i]->GetSumOfWeights();
    
    // find largest difference from mean value
    Double_t maxDiff = -1.;
    for(Int_t j=0; j<nTests; j++) {
      if (TMath::Abs(nJPsi[j][i]-nJPsiMean[i]) > maxDiff) maxDiff = TMath::Abs(nJPsi[j][i]-nJPsiMean[i]);
    }
    
    // print results
    printf("NJpsi %s = %4.0f ± %3.0f (%4.1f%%) ± %3.0f (%4.1f%%) --- maxDiff = %3.0f (%3.1f RMS)\n",gName[i].Data(),
	   nJPsiMean[i], enJPsiMean[i], 100.*enJPsiMean[i]/nJPsiMean[i],
	   sysnJPsiMean[i], 100.*sysnJPsiMean[i]/nJPsiMean[i],
	   maxDiff, maxDiff / sysnJPsiMean[i]);
    
    // display histograms
    if (fineDisplay) sys1->cd(3*i+1);
    else sys1->cd(i+1);
    gPad->SetFrameBorderMode(0);
    gPad->SetTopMargin(0.01);
    gPad->SetBottomMargin(0.01);
    gNJpsi[i]->GetXaxis()->Set(nTests,-0.5,nTests-0.5);
    gNJpsi[i]->GetXaxis()->SetLabelColor(10);
    gNJpsi[i]->GetYaxis()->SetTitle("N_{J/#psi}");
    gNJpsi[i]->GetYaxis()->SetTitleSize(0.1);
    gNJpsi[i]->GetYaxis()->SetTitleOffset(0.4);
    gNJpsi[i]->GetYaxis()->SetLabelSize(0.08);
    gNJpsi[i]->SetMarkerStyle(20);
    gNJpsi[i]->Draw("ap");
    TPaveText *t = new TPaveText(0.2, 0.8, 0.5, 0.95,"NDC");
    t->AddText(0.,0.,gTitle[i].Data());
    t->SetFillStyle(0);
    t->SetBorderSize(0);
    t->SetTextFont(42);
    t->Draw();
    TLine *l0 = new TLine(-0.5, nJPsiMean[i], nTests-0.5, nJPsiMean[i]);
    l0->Draw();
    TLine *l1 = new TLine(-0.5, nJPsiMean[i]-sysnJPsiMean[i], nTests-0.5, nJPsiMean[i]-sysnJPsiMean[i]);
    l1->SetLineStyle(3);
    l1->Draw();
    TLine *l2 = new TLine(-0.5, nJPsiMean[i]+sysnJPsiMean[i], nTests-0.5, nJPsiMean[i]+sysnJPsiMean[i]);
    l2->SetLineStyle(3);
    l2->Draw();
    if (fineDisplay) {
      sys1->cd(3*i+2);
      gPad->SetFrameBorderMode(0);
      hNJpsi[i]->Draw();
      sys1->cd(3*i+3);
      gPad->SetFrameBorderMode(0);
      heNJpsi[i]->Draw();
    }
  }
  printf("\n");
  
  // compute NJpsiOver4080 means and errors and display histograms
  Double_t nJPsiOver4080Mean[nCentBins-2];
  Double_t enJPsiOver4080Mean[nCentBins-2];
  Double_t sysnJPsiOver4080Mean[nCentBins-2];
  TCanvas *sys2;
  if (fineDisplay) sys2 = new TCanvas("sys2", "systematics of NJpsi over 40-80%", 1200, 900);
  else sys2 = new TCanvas("sys2", "systematics of NJpsi over 40-80%", 800, 540);
  if (fineDisplay) sys2->Divide(3,nCentBins-2);
  else sys2->Divide(1,nCentBins-2);
  for(Int_t i=0; i<nCentBins-2; i++) {
    
    // compute NJpsi means and errors
    nJPsiOver4080Mean[i] = hNJpsiOver4080[i]->GetMean();
    sysnJPsiOver4080Mean[i] = hNJpsiOver4080[i]->GetRMS();
    //enJPsiOver4080Mean[i] = heNJpsiOver4080[i]->GetMean();
    enJPsiOver4080Mean[i] = TMath::Sqrt(enJPsiOver4080Mean2[i]*nTests)/hNJpsiOver4080[i]->GetSumOfWeights();
    
    // find largest difference from mean value
    Double_t maxDiff = -1.;
    for(Int_t j=0; j<nTests; j++) {
      Double_t r = nJPsi[j][i+1]/nJPsi[j][nCentBins-1];
      if (TMath::Abs(r-nJPsiOver4080Mean[i]) > maxDiff) maxDiff = TMath::Abs(r-nJPsiOver4080Mean[i]);
    }
    
    // print results
    printf("NJpsi %s = %5.3f ± %5.3f (%4.1f%%) ± %5.3f (%4.1f%%) --- maxDiff = %5.3f (%3.1f RMS)\n",gName[nCentBins+i].Data(),
	   nJPsiOver4080Mean[i], enJPsiOver4080Mean[i], 100.*enJPsiOver4080Mean[i]/nJPsiOver4080Mean[i],
	   sysnJPsiOver4080Mean[i], 100.*sysnJPsiOver4080Mean[i]/nJPsiOver4080Mean[i],
	   maxDiff, maxDiff / sysnJPsiOver4080Mean[i]);
    
    // display histograms
    if (fineDisplay) sys2->cd(3*i+1);
    else sys2->cd(i+1);
    gPad->SetFrameBorderMode(0);
    gPad->SetTopMargin(0.01);
    gPad->SetBottomMargin(0.01);
    gNJpsiOver4080[i]->GetXaxis()->Set(nTests,-0.5,nTests-0.5);
    gNJpsiOver4080[i]->GetXaxis()->SetLabelColor(10);
    gNJpsiOver4080[i]->GetYaxis()->SetTitle("N_{J/#psi} ratio");
    gNJpsiOver4080[i]->GetYaxis()->SetTitleSize(0.1);
    gNJpsiOver4080[i]->GetYaxis()->SetTitleOffset(0.4);
    gNJpsiOver4080[i]->GetYaxis()->SetLabelSize(0.08);
    gNJpsiOver4080[i]->SetMarkerStyle(20);
    gNJpsiOver4080[i]->Draw("ap");
    TPaveText *t = new TPaveText(0.2, 0.8, 0.5, 0.95,"NDC");
    t->AddText(0.,0.,gTitle[nCentBins+i].Data());
    t->SetFillStyle(0);
    t->SetBorderSize(0);
    t->SetTextFont(42);
    t->Draw();
    TLine *l0 = new TLine(-0.5, nJPsiOver4080Mean[i], nTests-0.5, nJPsiOver4080Mean[i]);
    l0->Draw();
    TLine *l1 = new TLine(-0.5, nJPsiOver4080Mean[i]-sysnJPsiOver4080Mean[i], nTests-0.5, nJPsiOver4080Mean[i]-sysnJPsiOver4080Mean[i]);
    l1->SetLineStyle(3);
    l1->Draw();
    TLine *l2 = new TLine(-0.5, nJPsiOver4080Mean[i]+sysnJPsiOver4080Mean[i], nTests-0.5, nJPsiOver4080Mean[i]+sysnJPsiOver4080Mean[i]);
    l2->SetLineStyle(3);
    l2->Draw();
    if (fineDisplay) {
      sys2->cd(3*i+2);
      gPad->SetFrameBorderMode(0);
      hNJpsiOver4080[i]->Draw();
      sys2->cd(3*i+3);
      gPad->SetFrameBorderMode(0);
      heNJpsiOver4080[i]->Draw();
    }
  }
  printf("\n");
  
  ComputeRAA(nJPsiMean, enJPsiMean, sysnJPsiMean, renJPsi, RAA, eRAA, eRAA2, eRAACorr, kTRUE);
  
  TGraphErrors *gRAA = 0x0, *gRAASys = 0x0;
  DisplayRAA(renJPsi, RAA, eRAA, eRAA2, eRAACorr, gRAA, gRAASys, kTRUE, vsCentClass);
  
  DisplayRAAWithPHENIX(renJPsi, RAA, eRAA, eRAA2, eRAACorr, gRAA, gRAASys, kTRUE, vsCentClass);
  
  DisplayRAAWithPHENIX2(renJPsi, RAA, eRAA, eRAA2, eRAACorr, gRAA, gRAASys, kTRUE, vsCentClass);
  
  DisplayRAAWithEPS09(renJPsi, RAA, eRAA, eRAA2, eRAACorr, gRAA, gRAASys, kTRUE, vsCentClass);
  
  Double_t renJPsiOver4080[nCentBins-2], RCP[nCentBins], eRCP[nCentBins], eRCP2[nCentBins];
  ComputeRCP(nJPsiOver4080Mean, enJPsiOver4080Mean, sysnJPsiOver4080Mean, renJPsiOver4080, RCP, eRCP, eRCP2, kTRUE);
  
  TGraphErrors *gRCP = 0x0, *gRCPSys = 0x0;
  DisplayRCP(renJPsiOver4080, RCP, eRCP, eRCP2, gRCP, gRCPSys, kTRUE, vsCentClass);
  
  DisplayRCPWithATLAS(renJPsiOver4080, RCP, eRCP, eRCP2, gRCP, gRCPSys, kTRUE, vsCentClass);
  
  Double_t reRAA[nCentBins], reRAA2[nCentBins];
  for(Int_t i=0; i<nCentBins; i++) {
    reRAA[i] = eRAA[i]/RAA[i];
    reRAA2[i] = eRAA2[i]/RAA[i];
  }
  TGraphErrors *gRCP2 = 0x0, *gRCPSys2 = 0x0;
  DisplayRCPWithATLAS2(renJPsi, RCP, reRAA, reRAA2, gRCP2, gRCPSys2, kTRUE, vsCentClass);
  
  DisplayRCPWithATLASandEE(renJPsiOver4080, RCP, eRCP, eRCP2, gRCP, gRCPSys, kTRUE, vsCentClass);
  
}

//------------------------------------------------------------------------
void ComputeRAA(Double_t *nJPsi, Double_t *enJPsi, Double_t *sysnJPsi, Double_t *renJPsi,
		Double_t *RAA, Double_t *eRAA, Double_t *eRAA2, Double_t &eRAACorr, Bool_t print)
{
  
  // Compute relative statistic & systematic errors
  Double_t rsysnJPsi[nCentBins];
  Double_t renMB[nCentBins];
  Double_t reTAA[nCentBins];
  for(Int_t i=0; i<nCentBins; i++) {    
    renJPsi[i]=enJPsi[i]/nJPsi[i];
    rsysnJPsi[i]=sysnJPsi[i]/nJPsi[i];
    renMB[i]=1./TMath::Sqrt(nMB[i]);
    reTAA[i]=eTAA[i]/TAA[i];
  }
  
  for(Int_t i=0; i<nCentBins; i++) { 
    
    RAA[i] = nJPsi[i]/BR/AccEff/nMB[i]/sigmaJPsipp*1000./TAA[i]; // 1000. factor mub to mb 
    
    eRAA[i] = RAA[i] * TMath::Sqrt(2.*2.*rsysnJPsi[i]*rsysnJPsi[i] +/* 0.1*0.1 +*/
				   eTrkCent[i]*eTrkCent[i] +
				   renMB[i]*renMB[i]);
    
    eRAA2[i] = RAA[i] * reTAA[i];
    
    Double_t eRAATot = TMath::Sqrt(eRAA[i]*eRAA[i] + eRAA2[i]*eRAA2[i]);
    
    if (print) printf("RAA %s = %5.3f ± %5.3f ± %5.3f\n", gName[i].Data(), RAA[i], renJPsi[i]*RAA[i], eRAATot);
    
  }
  if (print) printf("\n");
  
  // correlated error for RAA
  eRAACorr = TMath::Sqrt(eGen[nCentBins]*eGen[nCentBins] +
			 eTrkUncorr[nCentBins]*eTrkUncorr[nCentBins] +
			 eTrkCorr[nCentBins]*eTrkCorr[nCentBins] +
			 eTrkRes[nCentBins]*eTrkRes[nCentBins] +
			 eTrg[nCentBins]*eTrg[nCentBins] +
			 eBR*eBR +
			 esigmaJPsipp*esigmaJPsipp);
  
  if (print) printf("correlated RAA error = %4.1f%%\n\n", 100.*eRAACorr);
  
}

//------------------------------------------------------------------------
void DisplayRAA(Double_t *renJPsi, Double_t *RAA, Double_t *eRAA, Double_t *eRAA2, Double_t eRAACorr,
		TGraphErrors *&gRAA, TGraphErrors *&gRAASys, Bool_t display, Bool_t vsCentClass)
{
  
  // fill RAA graphs
  gRAA = new TGraphErrors(nCentBins-1);
  gRAASys = new TGraphErrors(nCentBins-1);
//  TGraphErrors *gRAASysAll = new TGraphErrors(nCentBins-1);
  TGraphErrors *gRAASysAll = new TGraphErrors(1);
  if (vsCentClass) {
    gRAASysAll->SetPoint(0,100.-x[nCentBins-1]-ex[nCentBins-1]+2,1.);
    gRAASysAll->SetPointError(0,1.,eRAACorr);
  } else {
    gRAASysAll->SetPoint(0,5.,1.);
    gRAASysAll->SetPointError(0,5.,eRAACorr);
  }
  for(Int_t i=nCentBins-1; i>=1; i--) {
    if (vsCentClass) {
      gRAA->SetPoint(nCentBins-1-i,100.-x[i],RAA[i]);
      gRAA->SetPointError(nCentBins-1-i,ex[i],renJPsi[i]*RAA[i]);
      gRAASys->SetPoint(nCentBins-1-i,100.-x[i],RAA[i]);
      //gRAASys->SetPointError(nCentBins-1-i,1.,eRAA[i]);
      gRAASys->SetPointError(nCentBins-1-i,1.,TMath::Sqrt(eRAA[i]*eRAA[i] + eRAA2[i]*eRAA2[i]));
/*      gRAASysAll->SetPoint(nCentBins-1-i,100.-x[i],1.);
      //gRAASysAll->SetPointError(nCentBins-1-i,ex[i],TMath::Sqrt(eRAA2[i]*eRAA2[i]/RAA[i]/RAA[i] + eRAACorr*eRAACorr));
      gRAASysAll->SetPointError(nCentBins-1-i,ex[i],eRAACorr);
*/    } else {
      gRAA->SetPoint(nCentBins-1-i,NPart[i],RAA[i]);
      gRAA->SetPointError(nCentBins-1-i,eNPart[i],renJPsi[i]*RAA[i]);
      gRAASys->SetPoint(nCentBins-1-i,NPart[i],RAA[i]);
      //gRAASys->SetPointError(nCentBins-1-i,eNPart[i],eRAA[i]);
      gRAASys->SetPointError(nCentBins-1-i,5.,TMath::Sqrt(eRAA[i]*eRAA[i] + eRAA2[i]*eRAA2[i]));
/*      if (i==nCentBins-1) gRAASysAll->SetPoint(nCentBins-1-i,NPart[i]-eNPart[i],1.);
      else if (i==1) gRAASysAll->SetPoint(nCentBins-1-i,NPart[i]+eNPart[i],1.);
      else gRAASysAll->SetPoint(nCentBins-1-i,NPart[i],1.);
      //gRAASysAll->SetPointError(nCentBins-1-i,eNPart[i],TMath::Sqrt(eRAA2[i]*eRAA2[i]/RAA[i]/RAA[i] + eRAACorr*eRAACorr));
      gRAASysAll->SetPointError(nCentBins-1-i,eNPart[i],eRAACorr);
*/    }
  }
  
  if (vsCentClass) {
    Double_t bin[nCentBins] = {20., 60., 80., 90., 100.};
    gRAA->GetXaxis()->Set(4, bin);
    for(Int_t i=nCentBins-1; i>=1; i--) {
      gRAA->GetXaxis()->SetBinLabel(nCentBins-i, label[i].Data());
    }
  } else {
    gRAA->GetXaxis()->Set(40, 0., 400.);
  }
  
  // plot RAA
  if (display) {
    new TCanvas("RAA","RAA",1400,1000);
    gPad->SetTopMargin(0.05);
    gPad->SetBottomMargin(0.12);
    gPad->SetRightMargin(0.05);
    if (vsCentClass) {
      gRAA->GetXaxis()->SetTitle("centrality");
      gRAA->GetXaxis()->SetLabelSize(0.07);
      gRAA->GetXaxis()->SetTitleOffset(1.1);
      gRAA->GetXaxis()->SetTitleSize(0.05);
      gRAA->GetXaxis()->SetTickLength(0.);
    } else {
//      gRAA->GetXaxis()->SetTitle("<N_{part}>");
      gRAA->GetXaxis()->SetTitle("<N_{part}> weighted by N_{coll}");
      gRAA->GetXaxis()->SetTitleSize(0.05);
      gRAA->GetXaxis()->SetTitleOffset(1.1);
      gRAA->GetXaxis()->SetLabelSize(0.05);
    }
    gRAA->GetYaxis()->SetTitle("R_{AA}");
    gRAA->GetYaxis()->SetTitleSize(0.05);
    gRAA->GetYaxis()->SetTitleOffset(0.9);
    gRAA->GetYaxis()->SetLabelSize(0.05);
    gRAA->GetYaxis()->SetRangeUser(0., 1.2);
    gRAA->SetMarkerStyle(21);
    gRAA->SetMarkerSize(2);
    gRAA->SetMarkerColor(1);
    gRAA->SetLineWidth(2);
    gRAA->SetLineColor(1);
    gRAA->Draw("ap");
    gRAASys->SetLineColor(1);
    gRAASys->SetFillStyle(0);
    gRAASys->Draw("e2");
    gRAASysAll->SetLineWidth(1);
    gRAASysAll->SetFillStyle(3001);
    gRAASysAll->SetFillColor(4);
    gRAASysAll->Draw("e2");
    TLine *l = vsCentClass ? new TLine(20., 1., 100., 1.) : new TLine(0., 1., 400., 1.);
    l->SetLineStyle(3);
    l->Draw();
    TPaveText* t = new TPaveText(0.15, 0.82, 0.9, 0.92,"NDC");
    t->AddText(0.,0.,"inclusive J/#psi in Pb-Pb  #sqrt{s_{NN}} = 2.76 TeV, 2.5<y<4, p_{T}>0");
    t->SetFillStyle(0);
    t->SetBorderSize(0);
    t->SetTextFont(42);
    t->Draw();
    if (vsCentClass) {
      TLine *l1 = new TLine(60, 0, 60, 0.03);
      l1->Draw();
      TLine *l2 = new TLine(80, 0, 80, 0.03);
      l2->Draw();
      TLine *l3 = new TLine(90, 0, 90, 0.03);
      l3->Draw();
    }
    ALICEseal("ALICE preliminary", 0.2, 0.2);
  }
  
}
 
//------------------------------------------------------------------------
void DisplayRAAWithPHENIX(Double_t *renJPsi, Double_t *RAA, Double_t *eRAA, Double_t *eRAA2, Double_t eRAACorr,
			  TGraphErrors *&gRAA, TGraphErrors *&gRAASys, Bool_t display, Bool_t vsCentClass)
{
  
  // PHENIX points
  Double_t cent_PHENIX[17] = {2.5, 7.5, 12.5, 17.5, 22.5, 27.5, 32.5, 37.5, 42.5, 47.5, 52.5, 57.5, 62.5, 67.5, 72.5, 77.5, 86.};
  Double_t ecent_PHENIX[17] = {2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 6.};
  Double_t NPart_PHENIX[17] = {350.8, 301.7, 255.7, 216.4, 182.4, 152.7, 126.8, 104.2, 84.6, 67.7, 53.2, 41., 30.8, 22.6, 16.1, 11.2, 5.6};
  Double_t eNPart_PHENIX[17] = {3.1, 4.7, 5.4, 5.6, 5.7, 5.9, 5.9, 5.8, 5.6, 5.4, 5., 4.5, 3.9, 3.4, 2.8, 2.2, 0.8};
  Double_t RAA_PHENIX[17] = {0.17, 0.16, 0.20, 0.25, 0.25, 0.35, 0.35, 0.41, 0.52, 0.49, 0.54, 0.80, 0.68, 0.72, 0.91, 1.03, 1.20};
  Double_t eRAA_PHENIX[17] = {0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.05, 0.06, 0.06, 0.08, 0.11, 0.10};
  Double_t syspRAA_PHENIX[17] = {0.04, 0.02, 0.03, 0.04, 0.03, 0.04, 0.04, 0.06, 0.07, 0.07, 0.09, 0.14, 0.13, 0.15, 0.21, 0.26, 0.23};
  Double_t sysmRAA_PHENIX[17] = {0.02, 0.02, 0.03, 0.04, 0.03, 0.04, 0.04, 0.06, 0.07, 0.07, 0.09, 0.14, 0.13, 0.15, 0.21, 0.26, 0.23};
  TGraphErrors *gRAA_PHENIX = new TGraphErrors(17);
  TGraphAsymmErrors *gRAA_PHENIXSys = new TGraphAsymmErrors(17);
  for (Int_t i=16; i>=0; i--) {
    if (vsCentClass) {
      gRAA_PHENIX->SetPoint(16-i, 100.-cent_PHENIX[i], RAA_PHENIX[i]);
      gRAA_PHENIX->SetPointError(16-i, ecent_PHENIX[i], eRAA_PHENIX[i]);
      gRAA_PHENIXSys->SetPoint(16-i, 100.-cent_PHENIX[i], RAA_PHENIX[i]);
      gRAA_PHENIXSys->SetPointError(16-i, 1., 1., sysmRAA_PHENIX[i], syspRAA_PHENIX[i]);
    } else {
      gRAA_PHENIX->SetPoint(16-i, NPart_PHENIX[i], RAA_PHENIX[i]);
      gRAA_PHENIX->SetPointError(16-i, eNPart_PHENIX[i], eRAA_PHENIX[i]);
      gRAA_PHENIXSys->SetPoint(16-i, NPart_PHENIX[i], RAA_PHENIX[i]);
      gRAA_PHENIXSys->SetPointError(16-i, 5., 5., sysmRAA_PHENIX[i], syspRAA_PHENIX[i]);
    }
  }
  TGraphErrors *gRAA_PHENIXSysAll = new TGraphErrors(1);
  if (vsCentClass) {
    gRAA_PHENIXSysAll->SetPoint(0,5.,1.);
    gRAA_PHENIXSysAll->SetPointError(0,1.,0.092);
  } else {
    gRAA_PHENIXSysAll->SetPoint(0,383.,1.);
    gRAA_PHENIXSysAll->SetPointError(0,5.,0.092);
  }
  
  // fill RAA graphs
  gRAA = new TGraphErrors(nCentBins-1);
  gRAASys = new TGraphErrors(nCentBins-1);
  TGraphErrors *gRAASysAll = new TGraphErrors(1);
  if (vsCentClass) {
    gRAASysAll->SetPoint(0,2.,1.);
    gRAASysAll->SetPointError(0,1.,eRAACorr);
  } else {
    gRAASysAll->SetPoint(0,395.,1.);
    gRAASysAll->SetPointError(0,5.,eRAACorr);
  }
  for(Int_t i=nCentBins-1; i>=1; i--) {
    if (vsCentClass) {
      gRAA->SetPoint(nCentBins-1-i,100.-x[i],RAA[i]);
      gRAA->SetPointError(nCentBins-1-i,ex[i],renJPsi[i]*RAA[i]);
      gRAASys->SetPoint(nCentBins-1-i,100.-x[i],RAA[i]);
      gRAASys->SetPointError(nCentBins-1-i,1.,TMath::Sqrt(eRAA[i]*eRAA[i] + eRAA2[i]*eRAA2[i]));
    } else {
      gRAA->SetPoint(nCentBins-1-i,NPart[i],RAA[i]);
      gRAA->SetPointError(nCentBins-1-i,eNPart[i],renJPsi[i]*RAA[i]);
      gRAASys->SetPoint(nCentBins-1-i,NPart[i],RAA[i]);
      gRAASys->SetPointError(nCentBins-1-i,5.,TMath::Sqrt(eRAA[i]*eRAA[i] + eRAA2[i]*eRAA2[i]));
    }
  }
  
  // plot RAA
  if (display) {
    new TCanvas("RAAvsPHENIX","RAAvsPHENIX");
    gPad->SetTopMargin(0.05);
    gPad->SetBottomMargin(0.12);
    gPad->SetRightMargin(0.05);
    if (vsCentClass) {
      gRAA_PHENIX->GetXaxis()->SetTitle("1 - centrality");
      gRAA_PHENIX->GetXaxis()->Set(100., 0., 100.);
    } else {
      //      gRAA->GetXaxis()->SetTitle("<N_{part}>");
      gRAA_PHENIX->GetXaxis()->SetTitle("<N_{part}*>");
      gRAA_PHENIX->GetXaxis()->Set(40, 0., 400.);
    }
    gRAA_PHENIX->GetXaxis()->SetTitleSize(0.05);
    gRAA_PHENIX->GetXaxis()->SetTitleOffset(1.1);
    gRAA_PHENIX->GetXaxis()->SetLabelSize(0.05);
    gRAA_PHENIX->GetYaxis()->SetTitle("R_{AA}");
    gRAA_PHENIX->GetYaxis()->SetTitleSize(0.05);
    gRAA_PHENIX->GetYaxis()->SetTitleOffset(0.9);
    gRAA_PHENIX->GetYaxis()->SetLabelSize(0.05);
    gRAA_PHENIX->GetYaxis()->SetRangeUser(0., 1.5);
    gRAA_PHENIX->SetMarkerStyle(24);
    gRAA_PHENIX->SetMarkerColor(4);
    gRAA_PHENIX->SetLineWidth(2);
    gRAA_PHENIX->SetLineColor(4);
    gRAA_PHENIX->Draw("ap");
    gRAA_PHENIXSys->SetLineColor(4);
    gRAA_PHENIXSys->SetFillStyle(0);
    gRAA_PHENIXSys->Draw("e2");
    gRAA_PHENIXSysAll->SetLineWidth(1);
    gRAA_PHENIXSysAll->SetFillStyle(3001);
    gRAA_PHENIXSysAll->SetFillColor(4);
    gRAA_PHENIXSysAll->Draw("e2");
    gRAA->SetMarkerStyle(21);
    gRAA->SetMarkerColor(2);
    gRAA->SetLineWidth(2);
    gRAA->SetLineColor(2);
    gRAA->Draw("p");
    gRAASys->SetLineColor(2);
    gRAASys->SetFillStyle(0);
    gRAASys->Draw("e2");
    gRAASysAll->SetLineWidth(1);
    gRAASysAll->SetFillStyle(3001);
    gRAASysAll->SetFillColor(2);
    gRAASysAll->Draw("e2");
    TLine *l = vsCentClass ? new TLine(0., 1., 100., 1.) : new TLine(0., 1., 400., 1.);
    l->SetLineStyle(3);
    l->Draw();
    TLegend *lg = new TLegend(0.17, 0.81, 0.92, 0.93,"","NDC");
    lg->AddEntry(gRAA,"ALICE (Pb-Pb #sqrt{s_{NN}} = 2.76 TeV), 2.5<y<4, p_{T}>0 (preliminary)","p");
    lg->AddEntry(gRAA_PHENIX,"PHENIX (Au-Au #sqrt{s_{NN}} = 0.2 TeV), 1.2<|y|<2.2, p_{T}>0 (arXiv:1103.6269)","p");
    lg->SetFillStyle(0);
    lg->SetBorderSize(0);
    lg->SetTextFont(42);
    lg->SetMargin(0.05);
    lg->Draw();
    if (!vsCentClass) {
      TPaveText *t3 = new TPaveText(0.11, 0.11, 0.51, 0.21,"NDC");
      t3->AddText("(*) ALICE <N_{part}> is weighted by N_{coll}");
      t3->SetFillStyle(0);
      t3->SetBorderSize(0);
      t3->SetTextFont(42);
      t3->Draw();
    }
  }
  
}

//------------------------------------------------------------------------
void DisplayRAAWithPHENIX2(Double_t *renJPsi, Double_t *RAA, Double_t *eRAA, Double_t *eRAA2, Double_t eRAACorr,
			   TGraphErrors *&gRAA, TGraphErrors *&gRAASys, Bool_t display, Bool_t vsCentClass)
{
  
  // PHENIX points forward rapidity
  Double_t cent_PHENIX[17] = {2.5, 7.5, 12.5, 17.5, 22.5, 27.5, 32.5, 37.5, 42.5, 47.5, 52.5, 57.5, 62.5, 67.5, 72.5, 77.5, 86.};
  Double_t ecent_PHENIX[17] = {2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 6.};
  Double_t NPart_PHENIX[17] = {350.8, 301.7, 255.7, 216.4, 182.4, 152.7, 126.8, 104.2, 84.6, 67.7, 53.2, 41., 30.8, 22.6, 16.1, 11.2, 5.6};
  Double_t eNPart_PHENIX[17] = {3.1, 4.7, 5.4, 5.6, 5.7, 5.9, 5.9, 5.8, 5.6, 5.4, 5., 4.5, 3.9, 3.4, 2.8, 2.2, 0.8};
  Double_t RAA_PHENIX[17] = {0.17, 0.16, 0.20, 0.25, 0.25, 0.35, 0.35, 0.41, 0.52, 0.49, 0.54, 0.80, 0.68, 0.72, 0.91, 1.03, 1.20};
  Double_t eRAA_PHENIX[17] = {0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.05, 0.06, 0.06, 0.08, 0.11, 0.10};
  Double_t syspRAA_PHENIX[17] = {0.04, 0.02, 0.03, 0.04, 0.03, 0.04, 0.04, 0.06, 0.07, 0.07, 0.09, 0.14, 0.13, 0.15, 0.21, 0.26, 0.23};
  Double_t sysmRAA_PHENIX[17] = {0.02, 0.02, 0.03, 0.04, 0.03, 0.04, 0.04, 0.06, 0.07, 0.07, 0.09, 0.14, 0.13, 0.15, 0.21, 0.26, 0.23};
  TGraphErrors *gRAA_PHENIX = new TGraphErrors(17);
  TGraphAsymmErrors *gRAA_PHENIXSys = new TGraphAsymmErrors(17);
  for (Int_t i=16; i>=0; i--) {
    if (vsCentClass) {
      gRAA_PHENIX->SetPoint(16-i, 100.-cent_PHENIX[i], RAA_PHENIX[i]);
      gRAA_PHENIX->SetPointError(16-i, ecent_PHENIX[i], eRAA_PHENIX[i]);
      gRAA_PHENIXSys->SetPoint(16-i, 100.-cent_PHENIX[i], RAA_PHENIX[i]);
      gRAA_PHENIXSys->SetPointError(16-i, 1., 1., sysmRAA_PHENIX[i], syspRAA_PHENIX[i]);
    } else {
      gRAA_PHENIX->SetPoint(16-i, NPart_PHENIX[i], RAA_PHENIX[i]);
      gRAA_PHENIX->SetPointError(16-i, eNPart_PHENIX[i], eRAA_PHENIX[i]);
      gRAA_PHENIXSys->SetPoint(16-i, NPart_PHENIX[i], RAA_PHENIX[i]);
      gRAA_PHENIXSys->SetPointError(16-i, 5., 5., sysmRAA_PHENIX[i], syspRAA_PHENIX[i]);
    }
  }
  TGraphErrors *gRAA_PHENIXSysAll = new TGraphErrors(1);
  if (vsCentClass) {
    gRAA_PHENIXSysAll->SetPoint(0,5.,1.);
    gRAA_PHENIXSysAll->SetPointError(0,1.,0.092);
  } else {
    gRAA_PHENIXSysAll->SetPoint(0,383.,1.);
    gRAA_PHENIXSysAll->SetPointError(0,5.,0.092);
  }
  
  // PHENIX points mid rapidity
  Double_t cent_PHENIX2[8] = {2.5, 7.5, 12.5, 17.5, 25., 35., 50., 76.5};
  Double_t ecent_PHENIX2[8] = {2.5, 2.5, 2.5, 2.5, 5., 5., 10., 16.5};
  Double_t NPart_PHENIX2[8] = {351.4, 299.0, 253.9, 215.3, 166.6, 114.2, 58.4, 14.5};
  Double_t eNPart_PHENIX2[8] = {2.9, 3.8, 4.3, 5.3, 5.4, 4.4, 4., 2.5};
  Double_t RAA_PHENIX2[8] = {0.26, 0.34, 0.36, 0.45, 0.58, 0.58, 0.65, 0.74};
  Double_t eRAA_PHENIX2[8] = {0.05, 0.06, 0.06, 0.07, 0.07, 0.08, 0.07, 0.12};
  Double_t sysRAA_PHENIX2[8] = {0.04, 0.05, 0.05, 0.07, 0.08, 0.08, 0.1, 0.21};
  TGraphErrors *gRAA_PHENIX2 = new TGraphErrors(8);
  TGraphErrors *gRAA_PHENIXSys2 = new TGraphErrors(8);
  for (Int_t i=7; i>=0; i--) {
    if (vsCentClass) {
      gRAA_PHENIX2->SetPoint(7-i, 100.-cent_PHENIX2[i], RAA_PHENIX2[i]);
      gRAA_PHENIX2->SetPointError(7-i, ecent_PHENIX2[i], eRAA_PHENIX2[i]);
      gRAA_PHENIXSys2->SetPoint(7-i, 100.-cent_PHENIX2[i], RAA_PHENIX2[i]);
      gRAA_PHENIXSys2->SetPointError(7-i, 1., sysRAA_PHENIX2[i]);
    } else {
      gRAA_PHENIX2->SetPoint(7-i, NPart_PHENIX2[i], RAA_PHENIX2[i]);
      gRAA_PHENIX2->SetPointError(7-i, eNPart_PHENIX2[i], eRAA_PHENIX2[i]);
      gRAA_PHENIXSys2->SetPoint(7-i, NPart_PHENIX2[i], RAA_PHENIX2[i]);
      gRAA_PHENIXSys2->SetPointError(7-i, 5., sysRAA_PHENIX2[i]);
    }
  }
  TGraphErrors *gRAA_PHENIXSysAll2 = new TGraphErrors(1);
  if (vsCentClass) {
    gRAA_PHENIXSysAll2->SetPoint(0,8.,1.);
    gRAA_PHENIXSysAll2->SetPointError(0,1.,0.12);
  } else {
    gRAA_PHENIXSysAll2->SetPoint(0,371.,1.);
    gRAA_PHENIXSysAll2->SetPointError(0,5.,0.12);
  }
  
  // fill RAA graphs
  gRAA = new TGraphErrors(nCentBins-1);
  gRAASys = new TGraphErrors(nCentBins-1);
  TGraphErrors *gRAASysAll = new TGraphErrors(1);
  if (vsCentClass) {
    gRAASysAll->SetPoint(0,2.,1.);
    gRAASysAll->SetPointError(0,1.,eRAACorr);
  } else {
    gRAASysAll->SetPoint(0,395.,1.);
    gRAASysAll->SetPointError(0,5.,eRAACorr);
  }
  for(Int_t i=nCentBins-1; i>=1; i--) {
    if (vsCentClass) {
      gRAA->SetPoint(nCentBins-1-i,100.-x[i],RAA[i]);
      gRAA->SetPointError(nCentBins-1-i,ex[i],renJPsi[i]*RAA[i]);
      gRAASys->SetPoint(nCentBins-1-i,100.-x[i],RAA[i]);
      gRAASys->SetPointError(nCentBins-1-i,1.,TMath::Sqrt(eRAA[i]*eRAA[i] + eRAA2[i]*eRAA2[i]));
    } else {
      gRAA->SetPoint(nCentBins-1-i,NPart[i],RAA[i]);
      gRAA->SetPointError(nCentBins-1-i,eNPart[i],renJPsi[i]*RAA[i]);
      gRAASys->SetPoint(nCentBins-1-i,NPart[i],RAA[i]);
      gRAASys->SetPointError(nCentBins-1-i,5.,TMath::Sqrt(eRAA[i]*eRAA[i] + eRAA2[i]*eRAA2[i]));
    }
  }
  
  // plot RAA
  if (display) {
    new TCanvas("RAAvsPHENIX2","RAAvsPHENIX2");
    gPad->SetTopMargin(0.05);
    gPad->SetBottomMargin(0.12);
    gPad->SetRightMargin(0.05);
    if (vsCentClass) {
      gRAA_PHENIX->GetXaxis()->SetTitle("1 - centrality");
      gRAA_PHENIX->GetXaxis()->Set(100., 0., 100.);
    } else {
      //      gRAA->GetXaxis()->SetTitle("<N_{part}>");
      gRAA_PHENIX->GetXaxis()->SetTitle("<N_{part}*>");
      gRAA_PHENIX->GetXaxis()->Set(40, 0., 400.);
    }
    gRAA_PHENIX->GetXaxis()->SetTitleSize(0.05);
    gRAA_PHENIX->GetXaxis()->SetTitleOffset(1.1);
    gRAA_PHENIX->GetXaxis()->SetLabelSize(0.05);
    gRAA_PHENIX->GetYaxis()->SetTitle("R_{AA}");
    gRAA_PHENIX->GetYaxis()->SetTitleSize(0.05);
    gRAA_PHENIX->GetYaxis()->SetTitleOffset(0.9);
    gRAA_PHENIX->GetYaxis()->SetLabelSize(0.05);
    gRAA_PHENIX->GetYaxis()->SetRangeUser(0., 1.5);
    gRAA_PHENIX->SetMarkerStyle(24);
    gRAA_PHENIX->SetMarkerColor(4);
    gRAA_PHENIX->SetLineWidth(2);
    gRAA_PHENIX->SetLineColor(4);
    gRAA_PHENIX->Draw("ap");
    gRAA_PHENIXSys->SetLineColor(4);
    gRAA_PHENIXSys->SetFillStyle(0);
    gRAA_PHENIXSys->Draw("e2");
    gRAA_PHENIXSysAll->SetLineWidth(1);
    gRAA_PHENIXSysAll->SetFillStyle(3001);
    gRAA_PHENIXSysAll->SetFillColor(4);
    gRAA_PHENIXSysAll->Draw("e2");
    gRAA_PHENIX2->SetMarkerStyle(27);
    gRAA_PHENIX2->SetMarkerColor(kGreen+2);
    gRAA_PHENIX2->SetLineWidth(2);
    gRAA_PHENIX2->SetLineColor(kGreen+2);
    gRAA_PHENIX2->Draw("p");
    gRAA_PHENIXSys2->SetLineColor(kGreen+2);
    gRAA_PHENIXSys2->SetFillStyle(0);
    gRAA_PHENIXSys2->Draw("e2");
    gRAA_PHENIXSysAll2->SetLineWidth(1);
    gRAA_PHENIXSysAll2->SetFillStyle(3001);
    gRAA_PHENIXSysAll2->SetFillColor(kGreen+2);
    gRAA_PHENIXSysAll2->Draw("e2");
    gRAA->SetMarkerStyle(21);
    gRAA->SetMarkerColor(2);
    gRAA->SetLineWidth(2);
    gRAA->SetLineColor(2);
    gRAA->Draw("p");
    gRAASys->SetLineColor(2);
    gRAASys->SetFillStyle(0);
    gRAASys->Draw("e2");
    gRAASysAll->SetLineWidth(1);
    gRAASysAll->SetFillStyle(3001);
    gRAASysAll->SetFillColor(2);
    gRAASysAll->Draw("e2");
    TLine *l = vsCentClass ? new TLine(0., 1., 100., 1.) : new TLine(0., 1., 400., 1.);
    l->SetLineStyle(3);
    l->Draw();
    TLegend *lg = new TLegend(0.17, 0.75, 0.92, 0.93,"","NDC");
    lg->AddEntry(gRAA,"ALICE (Pb-Pb #sqrt{s_{NN}} = 2.76 TeV), 2.5<y<4, p_{T}>0 (preliminary)","p");
    lg->AddEntry(gRAA_PHENIX,"PHENIX (Au-Au #sqrt{s_{NN}} = 0.2 TeV), 1.2<|y|<2.2, p_{T}>0 (arXiv:1103.6269)","p");
    lg->AddEntry(gRAA_PHENIX2,"PHENIX (Au-Au #sqrt{s_{NN}} = 0.2 TeV), |y|<0.35, p_{T}>0 (nucl-ex/0611020)","p");
    lg->SetFillStyle(0);
    lg->SetBorderSize(0);
    lg->SetTextFont(42);
    lg->SetMargin(0.05);
    lg->Draw();
    if (!vsCentClass) {
      TPaveText *t3 = new TPaveText(0.11, 0.11, 0.51, 0.21,"NDC");
      t3->AddText("(*) ALICE <N_{part}> is weighted by N_{coll}");
      t3->SetFillStyle(0);
      t3->SetBorderSize(0);
      t3->SetTextFont(42);
      t3->Draw();
    }
  }
  
}

//------------------------------------------------------------------------
void DisplayRAAWithEPS09(Double_t *renJPsi, Double_t *RAA, Double_t *eRAA, Double_t *eRAA2, Double_t eRAACorr,
			 TGraphErrors *&gRAA, TGraphErrors *&gRAASys, Bool_t display, Bool_t vsCentClass)
{
  
  // Shadowing from Ramona
  Double_t ShadN[4]={0.79,0.70,0.65,0.63};
  Double_t ShadLowN[4]={0.97,0.95,0.94,0.94};
  Double_t ShadHighN[4]={0.59,0.43,0.34,0.30};
  Double_t xshad[4]={40.,70.,85.,95.};
  TGraph *gShadN = new TGraph(4,xshad,ShadN);
  TGraph *gShadLowN = new TGraph(4,xshad,ShadLowN);
  TGraph *gShadHighN = new TGraph(4,xshad,ShadHighN);
  
  // fill RAA graphs
  gRAA = new TGraphErrors(nCentBins-1);
  gRAASys = new TGraphErrors(nCentBins-1);
  TGraphErrors *gRAASysAll = new TGraphErrors(1);
  if (vsCentClass) {
    gRAASysAll->SetPoint(0,2.,1.);
    gRAASysAll->SetPointError(0,1.,eRAACorr);
  } else {
    gRAASysAll->SetPoint(0,395.,1.);
    gRAASysAll->SetPointError(0,5.,eRAACorr);
  }
  for(Int_t i=nCentBins-1; i>=1; i--) {
    if (vsCentClass) {
      gRAA->SetPoint(nCentBins-1-i,100.-x[i],RAA[i]);
      gRAA->SetPointError(nCentBins-1-i,ex[i],renJPsi[i]*RAA[i]);
      gRAASys->SetPoint(nCentBins-1-i,100.-x[i],RAA[i]);
      gRAASys->SetPointError(nCentBins-1-i,1.,TMath::Sqrt(eRAA[i]*eRAA[i] + eRAA2[i]*eRAA2[i]));
    } else {
      gRAA->SetPoint(nCentBins-1-i,NPart[i],RAA[i]);
      gRAA->SetPointError(nCentBins-1-i,eNPart[i],renJPsi[i]*RAA[i]);
      gRAASys->SetPoint(nCentBins-1-i,NPart[i],RAA[i]);
      gRAASys->SetPointError(nCentBins-1-i,5.,TMath::Sqrt(eRAA[i]*eRAA[i] + eRAA2[i]*eRAA2[i]));
    }
  }
  
  if (vsCentClass) {
    Double_t bin[nCentBins] = {20., 60., 80., 90., 100.};
    gShadN->GetXaxis()->Set(4, bin);
    for(Int_t i=nCentBins-1; i>=1; i--) {
      gShadN->GetXaxis()->SetBinLabel(nCentBins-i, label[i].Data());
    }
  } else {
    gShadN->GetXaxis()->Set(40, 0., 400.);
  }
  
  // plot RAA
  if (display) {
    new TCanvas("RAAvsEPS","RAAvsEPS");
    gPad->SetTopMargin(0.05);
    gPad->SetBottomMargin(0.12);
    gPad->SetRightMargin(0.05);
    if (vsCentClass) {
      gShadN->GetXaxis()->SetTitle("centrality");
      gShadN->GetXaxis()->SetTickLength(0.);
    } else {
      //      gRAA->GetXaxis()->SetTitle("<N_{part}>");
      gShadN->GetXaxis()->SetTitle("<N_{part}*>");
      gShadN->GetXaxis()->Set(40, 0., 400.);
    }
    gShadN->GetXaxis()->SetTitleSize(0.05);
    gShadN->GetXaxis()->SetTitleOffset(1.1);
    gShadN->GetXaxis()->SetLabelSize(0.07);
    gShadN->GetYaxis()->SetTitle("R_{AA}");
    gShadN->GetYaxis()->SetTitleSize(0.05);
    gShadN->GetYaxis()->SetTitleOffset(0.9);
    gShadN->GetYaxis()->SetLabelSize(0.05);
    gShadN->GetYaxis()->SetRangeUser(0., 1.4);
    gShadN->SetMarkerColor(kRed);
    gShadN->SetMarkerStyle(20);
    gShadN->SetLineColor(kRed);
    gShadN->SetLineWidth(2);
    gShadN->Draw("al");
    gShadLowN->SetMarkerColor(kRed);
    gShadLowN->SetMarkerStyle(20);
    gShadLowN->SetLineColor(kRed);
    gShadLowN->SetLineWidth(2);
    gShadLowN->SetLineStyle(2);
    gShadLowN->Draw("l");
    gShadHighN->SetMarkerColor(kRed);
    gShadHighN->SetMarkerStyle(20);
    gShadHighN->SetLineColor(kRed);
    gShadHighN->SetLineWidth(2);
    gShadHighN->SetLineStyle(2);
    gShadHighN->Draw("l");
    gRAA->SetMarkerStyle(21);
    gRAA->SetMarkerColor(4);
    gRAA->SetLineWidth(2);
    gRAA->SetLineColor(4);
    gRAA->Draw("p");
    gRAASys->SetLineColor(4);
    gRAASys->SetFillStyle(0);
    gRAASys->Draw("e2");
    gRAASysAll->SetLineWidth(1);
    gRAASysAll->SetFillStyle(3001);
    gRAASysAll->SetFillColor(4);
    gRAASysAll->Draw("e2");
    TLine *l = vsCentClass ? new TLine(20., 1., 100., 1.) : new TLine(0., 1., 400., 1.);
    l->SetLineStyle(3);
    l->Draw();
    TLegend *lg = new TLegend(0.15, 0.77, 0.88, 0.90,"","NDC");
    lg->AddEntry(gRAA,"ALICE (Pb-Pb #sqrt{s_{NN}} = 2.76 TeV), 2.5<y<4, p_{T}>0 (preliminary)","p");
    lg->AddEntry(gShadN,"EPS09 (R. Vogt , priv. comm.)","l");
    lg->SetFillStyle(0);
    lg->SetBorderSize(0);
    lg->SetTextFont(42);
    lg->SetMargin(0.1);
    lg->Draw();
    if (!vsCentClass) {
      TPaveText *t3 = new TPaveText(0.11, 0.11, 0.51, 0.21,"NDC");
      t3->AddText("(*) ALICE <N_{part}> is weighted by N_{coll}");
      t3->SetFillStyle(0);
      t3->SetBorderSize(0);
      t3->SetTextFont(42);
      t3->Draw();
    } else {
      TLine *l1 = new TLine(60, 0, 60, 0.03);
      l1->Draw();
      TLine *l2 = new TLine(80, 0, 80, 0.03);
      l2->Draw();
      TLine *l3 = new TLine(90, 0, 90, 0.03);
      l3->Draw();
    }
  }
  
}

//------------------------------------------------------------------------
void ComputeRCP(Double_t *nJPsiOver4080, Double_t *enJPsiOver4080, Double_t *sysnJPsiOver4080,
		Double_t *renJPsiOver4080, Double_t *RCP, Double_t *eRCP, Double_t *eRCP2, Bool_t print)
{
  
  // Compute relative statistic & systematic errors
  Double_t renMB[nCentBins];
  for(Int_t i=0; i<nCentBins; i++) {    
    renMB[i]=1./TMath::Sqrt(nMB[i]);
  }
  Double_t rsysnJPsiOver4080[nCentBins-2];
  Double_t reTAAOver4080[nCentBins-2];
  for(Int_t i=0; i<nCentBins-2; i++) {    
    renJPsiOver4080[i]=enJPsiOver4080[i]/nJPsiOver4080[i];
    rsysnJPsiOver4080[i]=sysnJPsiOver4080[i]/nJPsiOver4080[i];
    reTAAOver4080[i]=eTAAOver4080[i]/TAAOver4080[i];
  }
  
  for(Int_t i=0; i<nCentBins-2; i++) {    
    
    RCP[i] = nJPsiOver4080[i]*nMB[nCentBins-1]/nMB[i+1]/TAAOver4080[i];
    
    eRCP[i] = RCP[i] * TMath::Sqrt(2.*2.*rsysnJPsiOver4080[i]*rsysnJPsiOver4080[i] +/* 0.1*0.1 +*/
				   renMB[i+1]*renMB[i+1] +
				   renMB[nCentBins-1]*renMB[nCentBins-1] +
				   eTrkCent[i+1]*eTrkCent[i+1] +
				   eTrkCent[nCentBins-1]*eTrkCent[nCentBins-1]);
    
    eRCP2[i] = RCP[i] * reTAAOver4080[i];
    
    Double_t eRCPTot = TMath::Sqrt(eRCP[i]*eRCP[i] + eRCP2[i]*eRCP2[i]);
    
    if (print) printf("RCP %s = %4.2f ± %4.2f ± %4.2f\n", gName[nCentBins+i].Data(), RCP[i], renJPsiOver4080[i]*RCP[i], eRCPTot);
    
  }
  if (print) printf("\n");
  
}

//------------------------------------------------------------------------
void DisplayRCP(Double_t *renJPsiOver4080, Double_t *RCP, Double_t *eRCP, Double_t *eRCP2,
		TGraphErrors *&gRCP, TGraphErrors *&gRCPSys, Bool_t display, Bool_t vsCentClass)
{
  
  // fill RCP graphs
  gRCP = new TGraphErrors(nCentBins-1);
  gRCPSys = new TGraphErrors(nCentBins-1);
  //  TGraphErrors *gRCPSys2 = new TGraphErrors(nCentBins-1);
  if (vsCentClass) {
    gRCP->SetPoint(0,100.-x[nCentBins-1],1.);
    gRCP->SetPointError(0,ex[nCentBins-1],0.);
    gRCPSys->SetPoint(0,100.-x[nCentBins-1],1.);
    gRCPSys->SetPointError(0,1.,0.);
    //    gRCPSys2->SetPoint(0,100.-x[nCentBins-1],1.);
    //    gRCPSys2->SetPointError(0,ex[nCentBins-1],0.);
  } else {
    gRCP->SetPoint(0,NPart[nCentBins-1],1.);
    gRCP->SetPointError(0,eNPart[nCentBins-1],0.);
    gRCPSys->SetPoint(0,NPart[nCentBins-1],1.);
    gRCPSys->SetPointError(0,5.,0.);
    //    gRCPSys2->SetPoint(0,NPart[nCentBins-1]-eNPart[nCentBins-1],1.);
    //    gRCPSys2->SetPointError(0,eNPart[nCentBins-1],0.);
  }
  for(Int_t i=nCentBins-2; i>=1; i--) {
    if (vsCentClass) {
      gRCP->SetPoint(nCentBins-1-i,100.-x[i],RCP[i-1]);
      gRCP->SetPointError(nCentBins-1-i,ex[i],renJPsiOver4080[i-1]*RCP[i-1]);
      gRCPSys->SetPoint(nCentBins-1-i,100.-x[i],RCP[i-1]);
      //      gRCPSys->SetPointError(nCentBins-1-i,1.,eRCP[i-1]);
      gRCPSys->SetPointError(nCentBins-1-i,1.,TMath::Sqrt(eRCP[i-1]*eRCP[i-1] + eRCP2[i-1]*eRCP2[i-1]));
      //      gRCPSys2->SetPoint(nCentBins-1-i,100.-x[i],1.);
      //      gRCPSys2->SetPointError(nCentBins-1-i,ex[i],eRCP2[i-1]/RCP[i-1]);
    } else {
      gRCP->SetPoint(nCentBins-1-i,NPart[i],RCP[i-1]);
      gRCP->SetPointError(nCentBins-1-i,eNPart[i],renJPsiOver4080[i-1]*RCP[i-1]);
      gRCPSys->SetPoint(nCentBins-1-i,NPart[i],RCP[i-1]);
      //      gRCPSys->SetPointError(nCentBins-1-i,eNPart[i],eRCP[i-1]);
      gRCPSys->SetPointError(nCentBins-1-i,5.,TMath::Sqrt(eRCP[i-1]*eRCP[i-1] + eRCP2[i-1]*eRCP2[i-1]));
      //      if (i==1) gRCPSys2->SetPoint(nCentBins-1-i,NPart[i]+eNPart[i],1.);
      //      else gRCPSys2->SetPoint(nCentBins-1-i,NPart[i],1.);
      //      gRCPSys2->SetPointError(nCentBins-1-i,eNPart[i],eRCP2[i-1]/RCP[i-1]);
    }
  }
  
  // plot RCP
  if (vsCentClass) {
    Double_t bin[nCentBins] = {20., 60., 80., 90., 100.};
    gRCP->GetXaxis()->Set(4, bin);
    for(Int_t i=nCentBins-1; i>=1; i--) {
      gRCP->GetXaxis()->SetBinLabel(nCentBins-i, label[i].Data());
    }
  } else {
    gRCP->GetXaxis()->Set(40., 0., 400.);
  }
  
  if (display) {
    new TCanvas("RCP","RCP",1400,1000);
    gPad->SetTopMargin(0.05);
    gPad->SetBottomMargin(0.12);
    gPad->SetRightMargin(0.05);
    if (vsCentClass) {
      gRCP->GetXaxis()->SetTitle("centrality");
      gRCP->GetXaxis()->SetLabelSize(0.07);
      gRCP->GetXaxis()->SetTitleOffset(1.1);
      gRCP->GetXaxis()->SetTitleSize(0.05);
      gRCP->GetXaxis()->SetTickLength(0.);
    } else {
      //      gRCP->GetXaxis()->SetTitle("<N_{part}>");
      gRCP->GetXaxis()->SetTitle("<N_{part}> weighted by N_{coll}");
      gRCP->GetXaxis()->SetTitleSize(0.05);
      gRCP->GetXaxis()->SetTitleOffset(1.1);
      gRCP->GetXaxis()->SetLabelSize(0.05);
    }
    gRCP->GetYaxis()->SetTitle("R_{CP} normalized to 40-80%");
    gRCP->GetYaxis()->SetTitleSize(0.05);
    gRCP->GetYaxis()->SetTitleOffset(0.9);
    gRCP->GetYaxis()->SetLabelSize(0.05);
    gRCP->GetYaxis()->SetRangeUser(0., 1.2);
    gRCP->SetMarkerStyle(21);
    gRCP->SetMarkerSize(2);
    gRCP->SetMarkerColor(1);
    gRCP->SetLineWidth(2);
    gRCP->SetLineColor(1);
    gRCP->Draw("ap");
    gRCPSys->SetLineColor(1);
    gRCPSys->SetFillStyle(0);
    gRCPSys->Draw("e2");
    /*    gRCPSys2->SetLineWidth(0);
     gRCPSys2->SetFillStyle(3001);
     gRCPSys2->SetFillColor(4);
     if (vsCentClass) gRCPSys2->Draw("e2");
     else gRCPSys2->Draw("e3");
     */    TLine *l = vsCentClass ? new TLine(20., 1., 100., 1.) : new TLine(0., 1., 400., 1.);
    l->SetLineStyle(3);
    l->Draw();
    TPaveText* t = new TPaveText(0.15, 0.2, 0.9, 0.3,"NDC");
    t->AddText(0.,0.,"inclusive J/#psi in Pb-Pb  #sqrt{s_{NN}} = 2.76 TeV, 2.5<y<4, p_{T}>0");
    t->SetFillStyle(0);
    t->SetBorderSize(0);
    t->SetTextFont(42);
    t->Draw();
    if (vsCentClass) {
      TLine *l1 = new TLine(60, 0, 60, 0.03);
      l1->Draw();
      TLine *l2 = new TLine(80, 0, 80, 0.03);
      l2->Draw();
      TLine *l3 = new TLine(90, 0, 90, 0.03);
      l3->Draw();
    }
    ALICEseal("ALICE preliminary", 0.3, 0.4);
  }
  
}

//------------------------------------------------------------------------
void DisplayRCPWithATLAS(Double_t *renJPsiOver4080, Double_t *RCP, Double_t *eRCP, Double_t *eRCP2,
			 TGraphErrors *&gRCP, TGraphErrors *&gRCPSys, Bool_t display, Bool_t vsCentClass)
{
  
  // ATLAS points
  Double_t RCP_ATLAS[nCentBins-1] = {0.46412, 0.61796, 0.71267, 1.};
  Double_t eRCP_ATLAS[nCentBins-1] = {0.052, 0.066, 0.066, 0.115};
  Double_t sysRCP_ATLAS[nCentBins-1] = {0.066, 0.087, 0.087, 0.122};
  
  //Centrality	Rcp	Stat	Syst.
  //40-100%		1	0.11	0.13
  //20-40%		0.72	0.06	0.09
  //10-20%		0.62	0.07	0.09
  //0-10%		0.47	0.05	0.07
  
  // fill RCP graphs
  gRCP = new TGraphErrors(nCentBins-1);
  gRCPSys = new TGraphErrors(nCentBins-1);
//  TGraphErrors *gRCPSys2 = new TGraphErrors(nCentBins-1);
  TGraphErrors *gRCP_ATLAS = new TGraphErrors(nCentBins-1);
  TGraphErrors *gRCPSys_ATLAS = new TGraphErrors(nCentBins-1);
  if (vsCentClass) {
    gRCP->SetPoint(0,100.-x[nCentBins-1]+0.5,1.);
    gRCP->SetPointError(0,ex[nCentBins-1],0.);
    gRCPSys->SetPoint(0,100.-x[nCentBins-1]+0.5,1.);
    gRCPSys->SetPointError(0,1.,0.);
    gRCP_ATLAS->SetPoint(0,100.-x[nCentBins-1]-0.5,RCP_ATLAS[nCentBins-2]);
    gRCP_ATLAS->SetPointError(0,ex[nCentBins-1],eRCP_ATLAS[nCentBins-2]);
    gRCPSys_ATLAS->SetPoint(0,100.-x[nCentBins-1]-0.5,RCP_ATLAS[nCentBins-2]);
    gRCPSys_ATLAS->SetPointError(0,1.,sysRCP_ATLAS[nCentBins-2]);
  } else {
    gRCP->SetPoint(0,NPart[nCentBins-1]+2.5,1.);
    gRCP->SetPointError(0,eNPart[nCentBins-1],0.);
    gRCPSys->SetPoint(0,NPart[nCentBins-1]+2.5,1.);
    gRCPSys->SetPointError(0,5.,0.);
    gRCP_ATLAS->SetPoint(0,NPart[nCentBins-1]-2.5,RCP_ATLAS[nCentBins-2]);
    gRCP_ATLAS->SetPointError(0,eNPart[nCentBins-1],eRCP_ATLAS[nCentBins-2]);
    gRCPSys_ATLAS->SetPoint(0,NPart[nCentBins-1]-2.5,RCP_ATLAS[nCentBins-2]);
    gRCPSys_ATLAS->SetPointError(0,5.,sysRCP_ATLAS[nCentBins-2]);
  }
  for(Int_t i=nCentBins-2; i>=1; i--) {
    if (vsCentClass) {
      gRCP->SetPoint(nCentBins-1-i,100.-x[i],RCP[i-1]);
      gRCP->SetPointError(nCentBins-1-i,ex[i],renJPsiOver4080[i-1]*RCP[i-1]);
      gRCPSys->SetPoint(nCentBins-1-i,100.-x[i],RCP[i-1]);
      gRCPSys->SetPointError(nCentBins-1-i,1.,TMath::Sqrt(eRCP[i-1]*eRCP[i-1] + eRCP2[i-1]*eRCP2[i-1]));
      gRCP_ATLAS->SetPoint(nCentBins-1-i,100.-x[i],RCP_ATLAS[i-1]);
      gRCP_ATLAS->SetPointError(nCentBins-1-i,ex[i],eRCP_ATLAS[i-1]);
      gRCPSys_ATLAS->SetPoint(nCentBins-1-i,100.-x[i],RCP_ATLAS[i-1]);
      gRCPSys_ATLAS->SetPointError(nCentBins-1-i,1.,sysRCP_ATLAS[i-1]);
    } else {
      gRCP->SetPoint(nCentBins-1-i,NPart[i],RCP[i-1]);
      gRCP->SetPointError(nCentBins-1-i,eNPart[i],renJPsiOver4080[i-1]*RCP[i-1]);
      gRCPSys->SetPoint(nCentBins-1-i,NPart[i],RCP[i-1]);
      gRCPSys->SetPointError(nCentBins-1-i,5.,TMath::Sqrt(eRCP[i-1]*eRCP[i-1] + eRCP2[i-1]*eRCP2[i-1]));
      gRCP_ATLAS->SetPoint(nCentBins-1-i,NPart[i],RCP_ATLAS[i-1]);
      gRCP_ATLAS->SetPointError(nCentBins-1-i,eNPart[i],eRCP_ATLAS[i-1]);
      gRCPSys_ATLAS->SetPoint(nCentBins-1-i,NPart[i],RCP_ATLAS[i-1]);
      gRCPSys_ATLAS->SetPointError(nCentBins-1-i,5.,sysRCP_ATLAS[i-1]);
    }
  }
  
  // plot RCP
  if (vsCentClass) {
    Double_t bin[nCentBins] = {20., 60., 80., 90., 100.};
    gRCP_ATLAS->GetXaxis()->Set(4, bin);
    for(Int_t i=nCentBins-1; i>=1; i--) {
      gRCP_ATLAS->GetXaxis()->SetBinLabel(nCentBins-i, label[i].Data());
    }
  } else {
    gRCP_ATLAS->GetXaxis()->Set(40., 0., 400.);
  }
  
  if (display) {
    new TCanvas("RCPvsATLAS","RCPvsATLAS");
    gPad->SetTopMargin(0.05);
    gPad->SetBottomMargin(0.12);
    gPad->SetRightMargin(0.05);
    if (vsCentClass) {
      gRCP_ATLAS->GetXaxis()->SetTitle("centrality");
      gRCP_ATLAS->GetXaxis()->SetLabelSize(0.07);
      gRCP_ATLAS->GetXaxis()->SetTitleOffset(1.1);
      gRCP_ATLAS->GetXaxis()->SetTitleSize(0.05);
      gRCP_ATLAS->GetXaxis()->SetTickLength(0.);
    } else {
//      gRCP_ATLAS->GetXaxis()->SetTitle("<N_{part}>");
      gRCP_ATLAS->GetXaxis()->SetTitle("<N_{part}> weighted by N_{coll}");
      gRCP_ATLAS->GetXaxis()->SetTitleSize(0.05);
      gRCP_ATLAS->GetXaxis()->SetTitleOffset(1.1);
      gRCP_ATLAS->GetXaxis()->SetLabelSize(0.05);
    }
    gRCP_ATLAS->GetYaxis()->SetTitle("R_{CP} normalized to 40-80%");
    gRCP_ATLAS->GetYaxis()->SetTitleSize(0.05);
    gRCP_ATLAS->GetYaxis()->SetTitleOffset(0.9);
    gRCP_ATLAS->GetYaxis()->SetLabelSize(0.05);
    gRCP_ATLAS->GetYaxis()->SetRangeUser(0., 1.2);
    gRCP_ATLAS->SetMarkerStyle(24);
    gRCP_ATLAS->SetMarkerColor(kGreen+2);
    gRCP_ATLAS->SetLineWidth(2);
    gRCP_ATLAS->SetLineColor(kGreen+2);
    gRCP_ATLAS->Draw("ap");
    gRCPSys_ATLAS->SetLineColor(kGreen+2);
    gRCPSys_ATLAS->SetFillStyle(0);
    gRCPSys_ATLAS->Draw("e2");
    gRCP->SetMarkerStyle(21);
    gRCP->SetMarkerColor(2);
    gRCP->SetLineWidth(2);
    gRCP->SetLineColor(2);
    gRCP->Draw("p");
    gRCPSys->SetLineColor(2);
    gRCPSys->SetFillStyle(0);
    gRCPSys->Draw("e2");
    TLine *l = vsCentClass ? new TLine(20., 1., 100., 1.) : new TLine(0., 1., 400., 1.);
    l->SetLineStyle(3);
    l->Draw();
    TLegend *lg = new TLegend(0.15, 0.2, 0.8, 0.4,"Pb-Pb  #sqrt{s_{NN}} = 2.76 TeV","NDC");
    lg->AddEntry(gRCP,"ALICE, 2.5<y<4, p_{T}>0 (preliminary)","p");
    lg->AddEntry(gRCP_ATLAS,"ATLAS, |y|<2.5, p_{T}>6.5 GeV/c (arXiv:1012.5419)","p");
    lg->SetFillStyle(0);
    lg->SetBorderSize(0);
    lg->SetTextFont(42);
    lg->SetMargin(0.1);
    lg->Draw();
    if (vsCentClass) {
      TLine *l1 = new TLine(60, 0, 60, 0.03);
      l1->Draw();
      TLine *l2 = new TLine(80, 0, 80, 0.03);
      l2->Draw();
      TLine *l3 = new TLine(90, 0, 90, 0.03);
      l3->Draw();
    }
  }
  
}

//------------------------------------------------------------------------
void DisplayRCPWithATLAS2(Double_t *renJPsi, Double_t *RCP, Double_t *reRAA, Double_t *reRAA2,
			 TGraphErrors *&gRCP, TGraphErrors *&gRCPSys, Bool_t display, Bool_t vsCentClass)
{
  
  // ATLAS points
  Double_t RCP_ATLAS[nCentBins-1] = {0.46412, 0.61796, 0.71267, 1.};
  Double_t eRCP_ATLAS[nCentBins-1] = {0.052, 0.066, 0.066, 0.115};
  Double_t sysRCP_ATLAS[nCentBins-1] = {0.066, 0.087, 0.087, 0.122};
  
  // fill RCP graphs
  gRCP = new TGraphErrors(nCentBins-1);
  gRCPSys = new TGraphErrors(nCentBins-1);
  //  TGraphErrors *gRCPSys2 = new TGraphErrors(nCentBins-1);
  TGraphErrors *gRCP_ATLAS = new TGraphErrors(nCentBins-1);
  TGraphErrors *gRCPSys_ATLAS = new TGraphErrors(nCentBins-1);
  if (vsCentClass) {
    gRCP->SetPoint(0,100.-x[nCentBins-1]+0.5,1.);
    gRCP->SetPointError(0,ex[nCentBins-1],renJPsi[nCentBins-1]);
    gRCPSys->SetPoint(0,100.-x[nCentBins-1]+0.5,1.);
    gRCPSys->SetPointError(0,1.,TMath::Sqrt(reRAA[nCentBins-1]*reRAA[nCentBins-1] + reRAA2[nCentBins-1]*reRAA2[nCentBins-1]));
    gRCP_ATLAS->SetPoint(0,100.-x[nCentBins-1]-0.5,RCP_ATLAS[nCentBins-2]);
    gRCP_ATLAS->SetPointError(0,ex[nCentBins-1],eRCP_ATLAS[nCentBins-2]);
    gRCPSys_ATLAS->SetPoint(0,100.-x[nCentBins-1]-0.5,RCP_ATLAS[nCentBins-2]);
    gRCPSys_ATLAS->SetPointError(0,1.,sysRCP_ATLAS[nCentBins-2]);
  } else {
    gRCP->SetPoint(0,NPart[nCentBins-1]+2.5,1.);
    gRCP->SetPointError(0,eNPart[nCentBins-1],renJPsi[nCentBins-1]);
    gRCPSys->SetPoint(0,NPart[nCentBins-1]+2.5,1.);
    gRCPSys->SetPointError(0,5.,TMath::Sqrt(reRAA[nCentBins-1]*reRAA[nCentBins-1] + reRAA2[nCentBins-1]*reRAA2[nCentBins-1]));
    gRCP_ATLAS->SetPoint(0,NPart[nCentBins-1]-2.5,RCP_ATLAS[nCentBins-2]);
    gRCP_ATLAS->SetPointError(0,eNPart[nCentBins-1],eRCP_ATLAS[nCentBins-2]);
    gRCPSys_ATLAS->SetPoint(0,NPart[nCentBins-1]-2.5,RCP_ATLAS[nCentBins-2]);
    gRCPSys_ATLAS->SetPointError(0,5.,sysRCP_ATLAS[nCentBins-2]);
  }
  for(Int_t i=nCentBins-2; i>=1; i--) {
    if (vsCentClass) {
      gRCP->SetPoint(nCentBins-1-i,100.-x[i],RCP[i-1]);
      gRCP->SetPointError(nCentBins-1-i,ex[i],renJPsi[i]*RCP[i-1]);
      gRCPSys->SetPoint(nCentBins-1-i,100.-x[i],RCP[i-1]);
      gRCPSys->SetPointError(nCentBins-1-i,1.,RCP[i-1]*TMath::Sqrt(reRAA[i]*reRAA[i] + reRAA2[i]*reRAA2[i]));
      gRCP_ATLAS->SetPoint(nCentBins-1-i,100.-x[i],RCP_ATLAS[i-1]);
      gRCP_ATLAS->SetPointError(nCentBins-1-i,ex[i],eRCP_ATLAS[i-1]);
      gRCPSys_ATLAS->SetPoint(nCentBins-1-i,100.-x[i],RCP_ATLAS[i-1]);
      gRCPSys_ATLAS->SetPointError(nCentBins-1-i,1.,sysRCP_ATLAS[i-1]);
    } else {
      gRCP->SetPoint(nCentBins-1-i,NPart[i],RCP[i-1]);
      gRCP->SetPointError(nCentBins-1-i,eNPart[i],renJPsi[i]*RCP[i-1]);
      gRCPSys->SetPoint(nCentBins-1-i,NPart[i],RCP[i-1]);
      gRCPSys->SetPointError(nCentBins-1-i,5.,RCP[i-1]*TMath::Sqrt(reRAA[i]*reRAA[i] + reRAA2[i]*reRAA2[i]));
      gRCP_ATLAS->SetPoint(nCentBins-1-i,NPart[i],RCP_ATLAS[i-1]);
      gRCP_ATLAS->SetPointError(nCentBins-1-i,eNPart[i],eRCP_ATLAS[i-1]);
      gRCPSys_ATLAS->SetPoint(nCentBins-1-i,NPart[i],RCP_ATLAS[i-1]);
      gRCPSys_ATLAS->SetPointError(nCentBins-1-i,5.,sysRCP_ATLAS[i-1]);
    }
  }
  
  // plot RCP
  if (vsCentClass) {
    Double_t bin[nCentBins] = {20., 60., 80., 90., 100.};
    gRCP_ATLAS->GetXaxis()->Set(4, bin);
    for(Int_t i=nCentBins-1; i>=1; i--) {
      gRCP_ATLAS->GetXaxis()->SetBinLabel(nCentBins-i, label[i].Data());
    }
  } else {
    gRCP_ATLAS->GetXaxis()->Set(40., 0., 400.);
  }
  
  if (display) {
    new TCanvas("RCPvsATLAS2","RCPvsATLAS2");
    gPad->SetTopMargin(0.05);
    gPad->SetBottomMargin(0.12);
    gPad->SetRightMargin(0.05);
    if (vsCentClass) {
      gRCP_ATLAS->GetXaxis()->SetTitle("centrality");
      gRCP_ATLAS->GetXaxis()->SetLabelSize(0.07);
      gRCP_ATLAS->GetXaxis()->SetTitleOffset(1.1);
      gRCP_ATLAS->GetXaxis()->SetTitleSize(0.05);
      gRCP_ATLAS->GetXaxis()->SetTickLength(0.);
    } else {
      //      gRCP_ATLAS->GetXaxis()->SetTitle("<N_{part}>");
      gRCP_ATLAS->GetXaxis()->SetTitle("<N_{part}> weighted by N_{coll}");
      gRCP_ATLAS->GetXaxis()->SetTitleSize(0.05);
      gRCP_ATLAS->GetXaxis()->SetTitleOffset(1.1);
      gRCP_ATLAS->GetXaxis()->SetLabelSize(0.05);
    }
    gRCP_ATLAS->GetYaxis()->SetTitle("R_{CP} normalized to 40-80%");
    gRCP_ATLAS->GetYaxis()->SetTitleSize(0.05);
    gRCP_ATLAS->GetYaxis()->SetTitleOffset(0.9);
    gRCP_ATLAS->GetYaxis()->SetLabelSize(0.05);
    gRCP_ATLAS->GetYaxis()->SetRangeUser(0., 1.2);
    gRCP_ATLAS->SetMarkerStyle(24);
    gRCP_ATLAS->SetMarkerColor(kGreen+2);
    gRCP_ATLAS->SetLineWidth(2);
    gRCP_ATLAS->SetLineColor(kGreen+2);
    gRCP_ATLAS->Draw("ap");
    gRCPSys_ATLAS->SetLineColor(kGreen+2);
    gRCPSys_ATLAS->SetFillStyle(0);
    gRCPSys_ATLAS->Draw("e2");
    gRCP->SetMarkerStyle(21);
    gRCP->SetMarkerColor(2);
    gRCP->SetLineWidth(2);
    gRCP->SetLineColor(2);
    gRCP->Draw("p");
    gRCPSys->SetLineColor(2);
    gRCPSys->SetFillStyle(0);
    gRCPSys->Draw("e2");
    TLine *l = vsCentClass ? new TLine(20., 1., 100., 1.) : new TLine(0., 1., 400., 1.);
    l->SetLineStyle(3);
    l->Draw();
    TLegend *lg = new TLegend(0.15, 0.2, 0.8, 0.4,"Pb-Pb  #sqrt{s_{NN}} = 2.76 TeV","NDC");
    lg->AddEntry(gRCP,"ALICE, 2.5<y<4, p_{T}>0 (preliminary)","p");
    lg->AddEntry(gRCP_ATLAS,"ATLAS, |y|<2.5, p_{T}>6.5 GeV/c (arXiv:1012.5419)","p");
    lg->SetFillStyle(0);
    lg->SetBorderSize(0);
    lg->SetTextFont(42);
    lg->SetMargin(0.1);
    lg->Draw();
    if (vsCentClass) {
      TLine *l1 = new TLine(60, 0, 60, 0.03);
      l1->Draw();
      TLine *l2 = new TLine(80, 0, 80, 0.03);
      l2->Draw();
      TLine *l3 = new TLine(90, 0, 90, 0.03);
      l3->Draw();
    }
  }
  
}

//------------------------------------------------------------------------
void DisplayRCPWithATLASandEE(Double_t *renJPsiOver4080, Double_t *RCP, Double_t *eRCP, Double_t *eRCP2,
			      TGraphErrors *&gRCP, TGraphErrors *&gRCPSys, Bool_t display, Bool_t vsCentClass)
{
  
  // ATLAS points
  Double_t RCP_ATLAS[nCentBins-1] = {0.46412, 0.61796, 0.71267, 1.};
  Double_t eRCP_ATLAS[nCentBins-1] = {0.052, 0.066, 0.066, 0.115};
  Double_t sysRCP_ATLAS[nCentBins-1] = {0.066, 0.087, 0.087, 0.122};
  
  // ALICE point
  Double_t RCP_ALICEee[2] = {1., 0.83};
  Double_t eRCP_ALICEee[2] = {0., 0.24};
  Double_t sysRCP_ALICEee[2] = {0., 0.18};
  Double_t cent[2] = {40., 80.};
  Double_t ecent[2] = {20., 20.};
  Double_t syscent[2] = {1., 1.};
  TGraphErrors *gRCP_ALICEee = new TGraphErrors(4,cent,RCP_ALICEee,ecent,eRCP_ALICEee);
  TGraphErrors *gRCP_ALICEeeSys = new TGraphErrors(4,cent,RCP_ALICEee,syscent,sysRCP_ALICEee);
  
  // fill RCP graphs
  gRCP = new TGraphErrors(nCentBins-1);
  gRCPSys = new TGraphErrors(nCentBins-1);
  //  TGraphErrors *gRCPSys2 = new TGraphErrors(nCentBins-1);
  TGraphErrors *gRCP_ATLAS = new TGraphErrors(nCentBins-1);
  TGraphErrors *gRCPSys_ATLAS = new TGraphErrors(nCentBins-1);
  if (vsCentClass) {
    gRCP->SetPoint(0,100.-x[nCentBins-1]+0.5,1.);
    gRCP->SetPointError(0,ex[nCentBins-1],0.);
    gRCPSys->SetPoint(0,100.-x[nCentBins-1]+0.5,1.);
    gRCPSys->SetPointError(0,1.,0.);
    gRCP_ATLAS->SetPoint(0,100.-x[nCentBins-1]-0.5,RCP_ATLAS[nCentBins-2]);
    gRCP_ATLAS->SetPointError(0,ex[nCentBins-1],eRCP_ATLAS[nCentBins-2]);
    gRCPSys_ATLAS->SetPoint(0,100.-x[nCentBins-1]-0.5,RCP_ATLAS[nCentBins-2]);
    gRCPSys_ATLAS->SetPointError(0,1.,sysRCP_ATLAS[nCentBins-2]);
  } else {
    gRCP->SetPoint(0,NPart[nCentBins-1]+2.5,1.);
    gRCP->SetPointError(0,eNPart[nCentBins-1],0.);
    gRCPSys->SetPoint(0,NPart[nCentBins-1]+2.5,1.);
    gRCPSys->SetPointError(0,5.,0.);
    gRCP_ATLAS->SetPoint(0,NPart[nCentBins-1]-2.5,RCP_ATLAS[nCentBins-2]);
    gRCP_ATLAS->SetPointError(0,eNPart[nCentBins-1],eRCP_ATLAS[nCentBins-2]);
    gRCPSys_ATLAS->SetPoint(0,NPart[nCentBins-1]-2.5,RCP_ATLAS[nCentBins-2]);
    gRCPSys_ATLAS->SetPointError(0,5.,sysRCP_ATLAS[nCentBins-2]);
  }
  for(Int_t i=nCentBins-2; i>=1; i--) {
    if (vsCentClass) {
      gRCP->SetPoint(nCentBins-1-i,100.-x[i],RCP[i-1]);
      gRCP->SetPointError(nCentBins-1-i,ex[i],renJPsiOver4080[i-1]*RCP[i-1]);
      gRCPSys->SetPoint(nCentBins-1-i,100.-x[i],RCP[i-1]);
      gRCPSys->SetPointError(nCentBins-1-i,1.,TMath::Sqrt(eRCP[i-1]*eRCP[i-1] + eRCP2[i-1]*eRCP2[i-1]));
      gRCP_ATLAS->SetPoint(nCentBins-1-i,100.-x[i],RCP_ATLAS[i-1]);
      gRCP_ATLAS->SetPointError(nCentBins-1-i,ex[i],eRCP_ATLAS[i-1]);
      gRCPSys_ATLAS->SetPoint(nCentBins-1-i,100.-x[i],RCP_ATLAS[i-1]);
      gRCPSys_ATLAS->SetPointError(nCentBins-1-i,1.,sysRCP_ATLAS[i-1]);
    } else {
      gRCP->SetPoint(nCentBins-1-i,NPart[i],RCP[i-1]);
      gRCP->SetPointError(nCentBins-1-i,eNPart[i],renJPsiOver4080[i-1]*RCP[i-1]);
      gRCPSys->SetPoint(nCentBins-1-i,NPart[i],RCP[i-1]);
      gRCPSys->SetPointError(nCentBins-1-i,5.,TMath::Sqrt(eRCP[i-1]*eRCP[i-1] + eRCP2[i-1]*eRCP2[i-1]));
      gRCP_ATLAS->SetPoint(nCentBins-1-i,NPart[i],RCP_ATLAS[i-1]);
      gRCP_ATLAS->SetPointError(nCentBins-1-i,eNPart[i],eRCP_ATLAS[i-1]);
      gRCPSys_ATLAS->SetPoint(nCentBins-1-i,NPart[i],RCP_ATLAS[i-1]);
      gRCPSys_ATLAS->SetPointError(nCentBins-1-i,5.,sysRCP_ATLAS[i-1]);
    }
  }
  
  // plot RCP
  if (vsCentClass) {
    Double_t bin[nCentBins] = {20., 60., 80., 90., 100.};
    gRCP_ATLAS->GetXaxis()->Set(4, bin);
    for(Int_t i=nCentBins-1; i>=1; i--) {
      gRCP_ATLAS->GetXaxis()->SetBinLabel(nCentBins-i, label[i].Data());
    }
  } else {
    gRCP_ATLAS->GetXaxis()->Set(40., 0., 400.);
  }
  
  if (display) {
    new TCanvas("RCPvsATLASandALICEee","RCPvsATLASandALICEee");
    gPad->SetTopMargin(0.05);
    gPad->SetBottomMargin(0.12);
    gPad->SetRightMargin(0.05);
    if (vsCentClass) {
      gRCP_ATLAS->GetXaxis()->SetTitle("centrality");
      gRCP_ATLAS->GetXaxis()->SetLabelSize(0.07);
      gRCP_ATLAS->GetXaxis()->SetTitleOffset(1.1);
      gRCP_ATLAS->GetXaxis()->SetTitleSize(0.05);
      gRCP_ATLAS->GetXaxis()->SetTickLength(0.);
    } else {
      //      gRCP_ATLAS->GetXaxis()->SetTitle("<N_{part}>");
      gRCP_ATLAS->GetXaxis()->SetTitle("<N_{part}> weighted by N_{coll}");
      gRCP_ATLAS->GetXaxis()->SetTitleSize(0.05);
      gRCP_ATLAS->GetXaxis()->SetTitleOffset(1.1);
      gRCP_ATLAS->GetXaxis()->SetLabelSize(0.05);
    }
    gRCP_ATLAS->GetYaxis()->SetTitle("R_{CP} normalized to 40-80%");
    gRCP_ATLAS->GetYaxis()->SetTitleSize(0.05);
    gRCP_ATLAS->GetYaxis()->SetTitleOffset(0.9);
    gRCP_ATLAS->GetYaxis()->SetLabelSize(0.05);
    gRCP_ATLAS->GetYaxis()->SetRangeUser(0., 1.2);
    gRCP_ATLAS->SetMarkerStyle(24);
    gRCP_ATLAS->SetMarkerColor(kGreen+2);
    gRCP_ATLAS->SetLineWidth(2);
    gRCP_ATLAS->SetLineColor(kGreen+2);
    gRCP_ATLAS->Draw("ap");
    gRCPSys_ATLAS->SetLineColor(kGreen+2);
    gRCPSys_ATLAS->SetFillStyle(0);
    gRCPSys_ATLAS->Draw("e2");
    gRCP->SetMarkerStyle(21);
    gRCP->SetMarkerColor(2);
    gRCP->SetLineWidth(2);
    gRCP->SetLineColor(2);
    gRCP->Draw("p");
    gRCPSys->SetLineColor(2);
    gRCPSys->SetFillStyle(0);
    gRCPSys->Draw("e2");
    gRCP_ALICEee->SetMarkerStyle(20);
    gRCP_ALICEee->SetMarkerColor(4);
    gRCP_ALICEee->SetLineWidth(2);
    gRCP_ALICEee->SetLineColor(4);
    gRCP_ALICEee->Draw("p");
    gRCP_ALICEeeSys->SetLineColor(4);
    gRCP_ALICEeeSys->SetFillStyle(0);
    gRCP_ALICEeeSys->Draw("e2");
    TLine *l = vsCentClass ? new TLine(20., 1., 100., 1.) : new TLine(0., 1., 400., 1.);
    l->SetLineStyle(3);
    l->Draw();
    TLegend *lg = new TLegend(0.15, 0.15, 0.8, 0.4,"Pb-Pb  #sqrt{s_{NN}} = 2.76 TeV","NDC");
    lg->AddEntry(gRCP,"ALICE, 2.5<y<4, p_{T}>0 (preliminary)","p");
    lg->AddEntry(gRCP_ALICEee,"ALICE, |y|<0.8, p_{T}>0 (preliminary)","p");
    lg->AddEntry(gRCP_ATLAS,"ATLAS, |y|<2.5, p_{T}>6.5 GeV/c (arXiv:1012.5419)","p");
    lg->SetFillStyle(0);
    lg->SetBorderSize(0);
    lg->SetTextFont(42);
    lg->SetMargin(0.1);
    lg->Draw();
    if (vsCentClass) {
      TLine *l1 = new TLine(60, 0, 60, 0.03);
      l1->Draw();
      TLine *l2 = new TLine(80, 0, 80, 0.03);
      l2->Draw();
      TLine *l3 = new TLine(90, 0, 90, 0.03);
      l3->Draw();
    }
  }
  
}

//------------------------------------------------------------------------
void ALICEseal(TString type, Double_t xPad, Double_t yPad)
{
  TVirtualPad* currPad = gPad;
  TPad *myPadLogo = new TPad("myPadLogo", "Pad for ALICE Logo",xPad,yPad,xPad+0.17,yPad+0.17);
  myPadLogo->SetFillColor(0);
  myPadLogo->SetBorderMode(0);
  myPadLogo->SetBorderSize(2);
  myPadLogo->SetFrameBorderMode(0);
  myPadLogo->SetLeftMargin(0.0);
  myPadLogo->SetTopMargin(0.0);
  myPadLogo->SetBottomMargin(0.0);
  myPadLogo->SetRightMargin(0.0);
  myPadLogo->Draw();
  myPadLogo->cd();
  TASImage *myAliceLogo = new TASImage("/Users/pillot/Pictures/alice_logo.png");
  myAliceLogo->Draw();
  currPad->cd();
  Double_t x1 = xPad - 0.07, y1 = yPad - 0.06;
  Double_t x2 = x1 + 0.25, y2 = y1 + 0.08;
  TPaveText* t1=new TPaveText(x1,y1,x2,y2,"NDC");
  t1->SetFillStyle(0);
  t1->SetBorderSize(0);
  t1->AddText(0.,0.,Form("%s", type.Data()));
  t1->SetTextColor(kRed);
  t1->SetTextFont(42);
  t1->Draw();
  TPaveText* t2=new TPaveText(x1+0.06,y1-0.06,x2-0.06,y2-0.06,"NDC");
  t2->SetFillStyle(0);
  t2->SetBorderSize(0);
  t2->SetTextColor(kRed);
  t2->SetTextFont(52);
  //TDatime dt;
  //TString today = Form("%02i/%02i/%4i", dt.GetDay(), dt.GetMonth(), dt.GetYear());
  //t2->AddText(0.,0.,today.Data());
  //t2->Draw();
}

