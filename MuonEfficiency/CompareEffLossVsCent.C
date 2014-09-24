//
//  CompareEffLossVsCent.C
//  aliroot_dev
//
//  Created by philippe pillot on 19/09/2014.
//  Copyright (c) 2014 Philippe Pillot. All rights reserved.
//

#if !defined(__CINT__) || defined(__MAKECINT__)

#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TH1F.h>
#include <TMath.h>
#include <TFitResult.h>
#include <TMatrixD.h>

#endif


// data block ([dpT1, dy1], [dpT2, dy2], ...):
//  0 = [0-8, 2.5-4]
//  1 = [0-0.3, 2.5-4], [0.3-8, 2.5-4]
//  2 = [0-0.3, 2.5-4], [0.3-1, 2.5-4], [1-8, 2.5-4]
//  3 = [0-1, 2.5-4], [1-2, 2.5-4], [2-3, 2.5-4], [3-4, 2.5-4], [4-5, 2.5-4], [5-6, 2.5-4], [6-8, 2.5-4]
//  4 = [0-0.5, 2.5-4], [0.5-1, 2.5-4], [1-1.5, 2.5-4], [1.5-2, 2.5-4], [2-2.5, 2.5-4], [2.5-3, 2.5-4], [3-3.5, 2.5-4], [3.5-4, 2.5-4], [4-4.5, 2.5-4], [4.5-5, 2.5-4], [5-5.5, 2.5-4], [5.5-6, 2.5-4], [6-8, 2.5-4]
//  5 = [0-2, 2.5-4], [2-5, 2.5-4], [5-8, 2.5-4]
//  6 = [0-8, 2.5-2.75], [0-8, 2.75-3], [0-8, 3-3.25], [0-8, 3.25-3.5], [0-8, 3.5-3.75], [0-8, 3.75-4]
//  7 = [0-8, 2.5-3], [0-8, 3-3.5], [0-8, 3.5-4]
//  8 = all
//  9 = 3+6
//  10 = 4+6
const Int_t nCases[11] = {1, 2, 3, 7, 13, 3, 6, 3, 36, 13, 19};
const UInt_t iBlock = 8;

// fix parameter in Eloss = a*(100-cent)^b
// -1 = fixed to 2.,2.5,3.,3.5,4.,4.5,5.,5.5
// -2 = free
Double_t b = -1.;

Double_t EffVsCent(Double_t *x, Double_t *p);
Double_t ELoss(Double_t *x, Double_t *p);

//______________________________________________________________________
void CompareEffLossVsCent()
{
  /// compare efficiency loss versus centrality between different pT/y bins
  
  // draw dNdEta
  TF1 *fMultVsCent = new TF1("fMultVsCent","[0]+[1]*TMath::Power(100.-x,[2])",0.,80.);
  Double_t centmult[4][9] = {
    {2.5,7.5,15.,25.,35.,45.,55.,65.,75.},
    {2.5,2.5,5.,5.,5.,5.,5.,5.,5.},
    {1601,1294,966,649,426,261,149,76,35},
    {60,49,37,23,15,9,6,4,2}
  };
  TGraphErrors *gm = new TGraphErrors(9, centmult[0], centmult[2], centmult[1], centmult[3]);
  gm->SetNameTitle("gm","multiplicity");
  TCanvas *cm = new TCanvas();
  gm->Draw("ap");
  fMultVsCent->SetParameters(20.,0.001,2.);
  gm->Fit(fMultVsCent,"IR");
  
  TString *sCase = new TString[nCases[iBlock]];
  Int_t *cCase = new Int_t[nCases[iBlock]];
  TGraphErrors **g = new TGraphErrors*[nCases[iBlock]];
  
  //###############################
  // efficiency values ( https://twiki.cern.ch/twiki/bin/viewauth/ALICE/JPsiPbPb2011Raa )
  Int_t iCase = -1;
  if (iBlock == 0 || iBlock == 8) {
    // 0-8GeV, 2.5-4
    ++iCase;
    sCase[iCase] = "0-8GeV, 2.5-4";
    cCase[iCase] = 1;
    const Int_t nBins = 9;
    Double_t centeff[4][nBins] = {
      {5.,15.,25.,35.,45.,55.,65.,75.,85.},
      {5.,5.,5.,5.,5.,5.,5.,5.,5.},
      {0.117331,0.121002,0.121340,0.123950,0.124594,0.125687,0.126275,0.126707,0.126709},
      {0.000676,0.000690,0.000684,0.000696,0.000693,0.000701,0.000711,0.000706,0.000704}
    };
/*    const Int_t nBins = 5; // estimated for photo-produced JPsi
    Double_t centeff[4][nBins] = {
      {5.,20.,40.,60.,80.},
      {5.,10.,10.,10.,10.},
      {0.184,0.191,0.197,0.199,0.2},
      {0.00552,0.002865,0.00197,0.000995,0.0002}
    };
*/    g[iCase] = new TGraphErrors(nBins, centeff[0], centeff[2], centeff[1], centeff[3]);
    g[iCase]->SetNameTitle(Form("g%d",iCase+1),sCase[iCase].Data());
    g[iCase]->SetLineColor(cCase[iCase]);
  }
  //-------------------------------
  if (iBlock == 1 || iBlock == 2 || iBlock == 8) {
    // 0-0.3GeV, 2.5-4
    ++iCase;
    sCase[iCase] = "0-0.3GeV, 2.5-4";
    cCase[iCase] = kBlue;
    const Int_t nBins = 5;
    Double_t centeff[4][nBins] = {
      {5.,20.,40.,60.,80.},
      {5.,10.,10.,10.,10.},
      {0.123295,0.143935,0.141405,0.139400,0.143132},
      {0.005607,0.004327,0.004518,0.004513,0.004460}
    };
    g[iCase] = new TGraphErrors(nBins, centeff[0], centeff[2], centeff[1], centeff[3]);
    g[iCase]->SetNameTitle(Form("g%d",iCase+1),sCase[iCase].Data());
    g[iCase]->SetLineColor(cCase[iCase]);
  }
  if (iBlock == 1 || iBlock == 8) {
    // 0.3-8GeV, 2.5-4
    ++iCase;
    sCase[iCase] = "0.3-8GeV, 2.5-4";
    cCase[iCase] = kBlue+1;
    const Int_t nBins = 9;
    Double_t centeff[4][nBins] = {
      {5.,15.,25.,35.,45.,55.,65.,75.,85.},
      {5.,5.,5.,5.,5.,5.,5.,5.,5.},
      {0.115476,0.118914,0.119074,0.121894,0.122481,0.123544,0.124184,0.124430,0.124724},
      {0.000672,0.000685,0.000678,0.000691,0.000687,0.000696,0.000705,0.000701,0.000699}
    };
    g[iCase] = new TGraphErrors(nBins, centeff[0], centeff[2], centeff[1], centeff[3]);
    g[iCase]->SetNameTitle(Form("g%d",iCase+1),sCase[iCase].Data());
    g[iCase]->SetLineColor(cCase[iCase]);
  }
  if (iBlock == 2 || iBlock == 8) {
    // 0.3-1GeV, 2.5-4
    ++iCase;
    sCase[iCase] = "0.3-1GeV, 2.5-4";
    cCase[iCase] = kBlue+2;
    const Int_t nBins = 5;
    Double_t centeff[4][nBins] = {
      {5.,20.,40.,60.,80.},
      {5.,10.,10.,10.,10.},
      {0.121276,0.122526,0.127127,0.130286,0.129799},
      {0.001872,0.001350,0.001391,0.001420,0.001455}
    };
    g[iCase] = new TGraphErrors(nBins, centeff[0], centeff[2], centeff[1], centeff[3]);
    g[iCase]->SetNameTitle(Form("g%d",iCase+1),sCase[iCase].Data());
    g[iCase]->SetLineColor(cCase[iCase]);
  }
  if (iBlock == 2 || iBlock == 8) {
    // 1-8GeV, 2.5-4
    ++iCase;
    sCase[iCase] = "1-8GeV, 2.5-4";
    cCase[iCase] = kBlue+3;
    const Int_t nBins = 5;
    Double_t centeff[4][nBins] = {
      {5.,20.,40.,60.,80.},
      {5.,10.,10.,10.,10.},
      {0.114304,0.118150,0.124166,0.126457,0.127444},
      {0.000717,0.000529,0.000556,0.000569,0.000576}
    };
    g[iCase] = new TGraphErrors(nBins, centeff[0], centeff[2], centeff[1], centeff[3]);
    g[iCase]->SetNameTitle(Form("g%d",iCase+1),sCase[iCase].Data());
    g[iCase]->SetLineColor(cCase[iCase]);
  }
  //-------------------------------
  if (iBlock == 3 || iBlock == 8 || iBlock == 9) {
    // 0-1GeV, 2.5-4
    ++iCase;
    sCase[iCase] = "0-1GeV, 2.5-4";
    cCase[iCase] = kGreen;
    const Int_t nBins = 3;
    Double_t centeff[4][nBins] = {
      {10.,30.,65.},
      {10.,10.,25.},
      {0.122080,0.125035,0.129318},
      {0.001282,0.001284,0.001052}
    };
    g[iCase] = new TGraphErrors(nBins, centeff[0], centeff[2], centeff[1], centeff[3]);
    g[iCase]->SetNameTitle(Form("g%d",iCase+1),sCase[iCase].Data());
    g[iCase]->SetLineColor(cCase[iCase]);
  }
  if (iBlock == 3 || iBlock == 8 || iBlock == 9) {
    // 1-2GeV, 2.5-4
    ++iCase;
    sCase[iCase] = "1-2GeV, 2.5-4";
    cCase[iCase] = kGreen-1;
    const Int_t nBins = 3;
    Double_t centeff[4][nBins] = {
      {10.,30.,65.},
      {10.,10.,25.},
      {0.101688,0.105794,0.107163},
      {0.000823,0.000844,0.000673}
    };
    g[iCase] = new TGraphErrors(nBins, centeff[0], centeff[2], centeff[1], centeff[3]);
    g[iCase]->SetNameTitle(Form("g%d",iCase+1),sCase[iCase].Data());
    g[iCase]->SetLineColor(cCase[iCase]);
  }
  if (iBlock == 3 || iBlock == 8 || iBlock == 9) {
    // 2-3GeV, 2.5-4
    ++iCase;
    sCase[iCase] = "2-3GeV, 2.5-4";
    cCase[iCase] = kGreen-2;
    const Int_t nBins = 3;
    Double_t centeff[4][nBins] = {
      {10.,30.,65.},
      {10.,10.,25.},
      {0.105007,0.108130,0.111975},
      {0.000928,0.000958,0.000773}
    };
    g[iCase] = new TGraphErrors(nBins, centeff[0], centeff[2], centeff[1], centeff[3]);
    g[iCase]->SetNameTitle(Form("g%d",iCase+1),sCase[iCase].Data());
    g[iCase]->SetLineColor(cCase[iCase]);
  }
  if (iBlock == 3 || iBlock == 8 || iBlock == 9) {
    // 3-4GeV, 2.5-4
    ++iCase;
    sCase[iCase] = "3-4GeV, 2.5-4";
    cCase[iCase] = kGreen-3;
    const Int_t nBins = 3;
    Double_t centeff[4][nBins] = {
      {10.,30.,65.},
      {10.,10.,25.},
      {0.128315,0.129528,0.133568},
      {0.001356,0.001394,0.001111}
    };
    g[iCase] = new TGraphErrors(nBins, centeff[0], centeff[2], centeff[1], centeff[3]);
    g[iCase]->SetNameTitle(Form("g%d",iCase+1),sCase[iCase].Data());
    g[iCase]->SetLineColor(cCase[iCase]);
  }
  if (iBlock == 3 || iBlock == 8 || iBlock == 9) {
    // 4-5GeV, 2.5-4
    ++iCase;
    sCase[iCase] = "4-5GeV, 2.5-4";
    cCase[iCase] = kGreen-4;
    const Int_t nBins = 3;
    Double_t centeff[4][nBins] = {
      {10.,30.,65.},
      {10.,10.,25.},
      {0.166024,0.172906,0.173167},
      {0.002173,0.002226,0.001789}
    };
    g[iCase] = new TGraphErrors(nBins, centeff[0], centeff[2], centeff[1], centeff[3]);
    g[iCase]->SetNameTitle(Form("g%d",iCase+1),sCase[iCase].Data());
    g[iCase]->SetLineColor(cCase[iCase]);
  }
  if (iBlock == 3 || iBlock == 8 || iBlock == 9) {
    // 5-6GeV, 2.5-4
    ++iCase;
    sCase[iCase] = "5-6GeV, 2.5-4";
    cCase[iCase] = kGreen-5;
    const Int_t nBins = 3;
    Double_t centeff[4][nBins] = {
      {10.,30.,65.},
      {10.,10.,25.},
      {0.209163,0.221530,0.224495},
      {0.003445,0.003609,0.002894}
    };
    g[iCase] = new TGraphErrors(nBins, centeff[0], centeff[2], centeff[1], centeff[3]);
    g[iCase]->SetNameTitle(Form("g%d",iCase+1),sCase[iCase].Data());
    g[iCase]->SetLineColor(cCase[iCase]);
  }
  //-------------------------------
  if (iBlock == 4 || iBlock == 8 || iBlock == 10) {
    // 0-0.5GeV, 2.5-4
    ++iCase;
    sCase[iCase] = "0-0.5GeV, 2.5-4";
    cCase[iCase] = kCyan;
    const Int_t nBins = 3;
    Double_t centeff[4][nBins] = {
      {10.,30.,65.},
      {10.,10.,25.},
      {0.128311,0.132861,0.136867},
      {0.002481,0.002567,0.002134}
    };
    g[iCase] = new TGraphErrors(nBins, centeff[0], centeff[2], centeff[1], centeff[3]);
    g[iCase]->SetNameTitle(Form("g%d",iCase+1),sCase[iCase].Data());
    g[iCase]->SetLineColor(cCase[iCase]);
  }
  if (iBlock == 4 || iBlock == 8 || iBlock == 10) {
    // 0.5-1GeV, 2.5-4
    ++iCase;
    sCase[iCase] = "0.5-1GeV, 2.5-4";
    cCase[iCase] = kCyan-1;
    const Int_t nBins = 3;
    Double_t centeff[4][nBins] = {
      {10.,30.,65.},
      {10.,10.,25.},
      {0.119606,0.122016,0.126470},
      {0.001495,0.001479,0.001205}
    };
    g[iCase] = new TGraphErrors(nBins, centeff[0], centeff[2], centeff[1], centeff[3]);
    g[iCase]->SetNameTitle(Form("g%d",iCase+1),sCase[iCase].Data());
    g[iCase]->SetLineColor(cCase[iCase]);
  }
  if (iBlock == 4 || iBlock == 8 || iBlock == 10) {
    // 1-1.5GeV, 2.5-4
    ++iCase;
    sCase[iCase] = "1-1.5GeV, 2.5-4";
    cCase[iCase] = kCyan-2;
    const Int_t nBins = 3;
    Double_t centeff[4][nBins] = {
      {10.,30.,65.},
      {10.,10.,25.},
      {0.103977,0.106959,0.109063},
      {0.001183,0.001203,0.000963}
    };
    g[iCase] = new TGraphErrors(nBins, centeff[0], centeff[2], centeff[1], centeff[3]);
    g[iCase]->SetNameTitle(Form("g%d",iCase+1),sCase[iCase].Data());
    g[iCase]->SetLineColor(cCase[iCase]);
  }
  if (iBlock == 4 || iBlock == 8 || iBlock == 10) {
    // 1.5-2GeV, 2.5-4
    ++iCase;
    sCase[iCase] = "1.5-2GeV, 2.5-4";
    cCase[iCase] = kCyan-3;
    const Int_t nBins = 3;
    Double_t centeff[4][nBins] = {
      {10.,30.,65.},
      {10.,10.,25.},
      {0.099167,0.104443,0.105303},
      {0.001142,0.001182,0.000939}
    };
    g[iCase] = new TGraphErrors(nBins, centeff[0], centeff[2], centeff[1], centeff[3]);
    g[iCase]->SetNameTitle(Form("g%d",iCase+1),sCase[iCase].Data());
    g[iCase]->SetLineColor(cCase[iCase]);
  }
  if (iBlock == 4 || iBlock == 8 || iBlock == 10) {
    // 2-2.5GeV, 2.5-4
    ++iCase;
    sCase[iCase] = "2-2.5GeV, 2.5-4";
    cCase[iCase] = kCyan-4;
    const Int_t nBins = 3;
    Double_t centeff[4][nBins] = {
      {10.,30.,65.},
      {10.,10.,25.},
      {0.102011,0.104693,0.109320},
      {0.001227,0.001270,0.001034}
    };
    g[iCase] = new TGraphErrors(nBins, centeff[0], centeff[2], centeff[1], centeff[3]);
    g[iCase]->SetNameTitle(Form("g%d",iCase+1),sCase[iCase].Data());
    g[iCase]->SetLineColor(cCase[iCase]);
  }
  if (iBlock == 4 || iBlock == 8 || iBlock == 10) {
    // 2.5-3GeV, 2.5-4
    ++iCase;
    sCase[iCase] = "2.5-3GeV, 2.5-4";
    cCase[iCase] = kCyan-5;
    const Int_t nBins = 3;
    Double_t centeff[4][nBins] = {
      {10.,30.,65.},
      {10.,10.,25.},
      {0.109254,0.112770,0.115379},
      {0.001424,0.001462,0.001164}
    };
    g[iCase] = new TGraphErrors(nBins, centeff[0], centeff[2], centeff[1], centeff[3]);
    g[iCase]->SetNameTitle(Form("g%d",iCase+1),sCase[iCase].Data());
    g[iCase]->SetLineColor(cCase[iCase]);
  }
  if (iBlock == 4 || iBlock == 8 || iBlock == 10) {
    // 3-3.5GeV, 2.5-4
    ++iCase;
    sCase[iCase] = "3-3.5GeV, 2.5-4";
    cCase[iCase] = kCyan-6;
    const Int_t nBins = 3;
    Double_t centeff[4][nBins] = {
      {10.,30.,65.},
      {10.,10.,25.},
      {0.123955,0.124426,0.126293},
      {0.001748,0.001774,0.001419}
    };
    g[iCase] = new TGraphErrors(nBins, centeff[0], centeff[2], centeff[1], centeff[3]);
    g[iCase]->SetNameTitle(Form("g%d",iCase+1),sCase[iCase].Data());
    g[iCase]->SetLineColor(cCase[iCase]);
  }
  if (iBlock == 4 || iBlock == 8 || iBlock == 10) {
    // 3.5-4GeV, 2.5-4
    ++iCase;
    sCase[iCase] = "3.5-4GeV, 2.5-4";
    cCase[iCase] = kCyan-7;
    const Int_t nBins = 3;
    Double_t centeff[4][nBins] = {
      {10.,30.,65.},
      {10.,10.,25.},
      {0.135567,0.137134,0.144129},
      {0.002162,0.002248,0.001781}
    };
    g[iCase] = new TGraphErrors(nBins, centeff[0], centeff[2], centeff[1], centeff[3]);
    g[iCase]->SetNameTitle(Form("g%d",iCase+1),sCase[iCase].Data());
    g[iCase]->SetLineColor(cCase[iCase]);
  }
  if (iBlock == 4 || iBlock == 8 || iBlock == 10) {
    // 4-4.5GeV, 2.5-4
    ++iCase;
    sCase[iCase] = "4-4.5GeV, 2.5-4";
    cCase[iCase] = kCyan-8;
    const Int_t nBins = 3;
    Double_t centeff[4][nBins] = {
      {10.,30.,65.},
      {10.,10.,25.},
      {0.155313,0.165941,0.165783},
      {0.002736,0.002839,0.002320}
    };
    g[iCase] = new TGraphErrors(nBins, centeff[0], centeff[2], centeff[1], centeff[3]);
    g[iCase]->SetNameTitle(Form("g%d",iCase+1),sCase[iCase].Data());
    g[iCase]->SetLineColor(cCase[iCase]);
  }
  if (iBlock == 4 || iBlock == 8 || iBlock == 10) {
    // 4.5-5GeV, 2.5-4
    ++iCase;
    sCase[iCase] = "4.5-5GeV, 2.5-4";
    cCase[iCase] = kCyan-9;
    const Int_t nBins = 3;
    Double_t centeff[4][nBins] = {
      {10.,30.,65.},
      {10.,10.,25.},
      {0.182936,0.183470,0.184306},
      {0.003568,0.003583,0.002819}
    };
    g[iCase] = new TGraphErrors(nBins, centeff[0], centeff[2], centeff[1], centeff[3]);
    g[iCase]->SetNameTitle(Form("g%d",iCase+1),sCase[iCase].Data());
    g[iCase]->SetLineColor(cCase[iCase]);
  }
  if (iBlock == 4 || iBlock == 8 || iBlock == 10) {
    // 5-5.5GeV, 2.5-4
    ++iCase;
    sCase[iCase] = "5-5.5GeV, 2.5-4";
    cCase[iCase] = kCyan-10;
    const Int_t nBins = 3;
    Double_t centeff[4][nBins] = {
      {10.,30.,65.},
      {10.,10.,25.},
      {0.202938,0.210695,0.214104},
      {0.004449,0.004632,0.003646}
    };
    g[iCase] = new TGraphErrors(nBins, centeff[0], centeff[2], centeff[1], centeff[3]);
    g[iCase]->SetNameTitle(Form("g%d",iCase+1),sCase[iCase].Data());
    g[iCase]->SetLineColor(cCase[iCase]);
  }
  if (iBlock == 4 || iBlock == 8 || iBlock == 10) {
    // 5.5-6GeV, 2.5-4
    ++iCase;
    sCase[iCase] = "5.5-6GeV, 2.5-4";
    cCase[iCase] = kCyan+1;
    const Int_t nBins = 3;
    Double_t centeff[4][nBins] = {
      {10.,30.,65.},
      {10.,10.,25.},
      {0.217722,0.237544,0.240261},
      {0.005428,0.005749,0.004732}
    };
    g[iCase] = new TGraphErrors(nBins, centeff[0], centeff[2], centeff[1], centeff[3]);
    g[iCase]->SetNameTitle(Form("g%d",iCase+1),sCase[iCase].Data());
    g[iCase]->SetLineColor(cCase[iCase]);
  }
  if (iBlock == 3 || iBlock == 4 || iBlock == 8 || iBlock == 9 || iBlock == 10) {
    // 6-8GeV, 2.5-4
    ++iCase;
    sCase[iCase] = "6-8GeV, 2.5-4";
    cCase[iCase] = kCyan+2;
    const Int_t nBins = 3;
    Double_t centeff[4][nBins] = {
      {10.,30.,65.},
      {10.,10.,25.},
      {0.257179,0.266027,0.275767},
      {0.004375,0.004470,0.003607}
    };
    g[iCase] = new TGraphErrors(nBins, centeff[0], centeff[2], centeff[1], centeff[3]);
    g[iCase]->SetNameTitle(Form("g%d",iCase+1),sCase[iCase].Data());
    g[iCase]->SetLineColor(cCase[iCase]);
  }
  //-------------------------------
  if (iBlock == 5 || iBlock == 8) {
    // 0-2GeV, 2.5-4
    ++iCase;
    sCase[iCase] = "0-2GeV, 2.5-4";
    cCase[iCase] = kBlue-1;
    const Int_t nBins = 7;
    Double_t centeff[4][nBins] = {
      {5.,15.,25.,35.,45.,55.,75.},
      {5.,5.,5.,5.,5.,5.,15.},
      {0.104442,0.109119,0.109090,0.111445,0.110886,0.113665,0.114081},
      {0.000934,0.000965,0.000940,0.000965,0.000953,0.000970,0.000569}
    };
    g[iCase] = new TGraphErrors(nBins, centeff[0], centeff[2], centeff[1], centeff[3]);
    g[iCase]->SetNameTitle(Form("g%d",iCase+1),sCase[iCase].Data());
    g[iCase]->SetLineColor(cCase[iCase]);
  }
  if (iBlock == 5 || iBlock == 8) {
    // 2-5GeV, 2.5-4
    ++iCase;
    sCase[iCase] = "2-5GeV, 2.5-4";
    cCase[iCase] = kBlue-2;
    const Int_t nBins = 7;
    Double_t centeff[4][nBins] = {
      {5.,15.,25.,35.,45.,55.,75.},
      {5.,5.,5.,5.,5.,5.,15.},
      {0.117998,0.120297,0.120532,0.123277,0.124889,0.123965,0.125386},
      {0.000990,0.001001,0.001004,0.001011,0.001018,0.001014,0.000592}
    };
    g[iCase] = new TGraphErrors(nBins, centeff[0], centeff[2], centeff[1], centeff[3]);
    g[iCase]->SetNameTitle(Form("g%d",iCase+1),sCase[iCase].Data());
    g[iCase]->SetLineColor(cCase[iCase]);
  }
  if (iBlock == 5 || iBlock == 8) {
    // 5-8GeV, 2.5-4
    ++iCase;
    sCase[iCase] = "5-8GeV, 2.5-4";
    cCase[iCase] = kBlue-3;
    const Int_t nBins = 7;
    Double_t centeff[4][nBins] = {
      {5.,15.,25.,35.,45.,55.,75.},
      {5.,5.,5.,5.,5.,5.,15.},
      {0.220664,0.227089,0.232206,0.234008,0.237939,0.240561,0.239267},
      {0.003652,0.003686,0.003713,0.003790,0.003766,0.003902,0.002192}
    };
    g[iCase] = new TGraphErrors(nBins, centeff[0], centeff[2], centeff[1], centeff[3]);
    g[iCase]->SetNameTitle(Form("g%d",iCase+1),sCase[iCase].Data());
    g[iCase]->SetLineColor(cCase[iCase]);
  }
  //-------------------------------
  if (iBlock == 6 || iBlock == 8 || iBlock == 9 || iBlock == 10) {
    // 0-8GeV, 2.5-2.75
    ++iCase;
    sCase[iCase] = "0-8GeV, 2.5-2.75";
    cCase[iCase] = kMagenta;
    const Int_t nBins = 3;
    Double_t centeff[4][nBins] = {
      {10.,30.,65.},
      {10.,10.,25.},
      {0.031810,0.033525,0.034324},
      {0.000578,0.000605,0.000487}
    };
    g[iCase] = new TGraphErrors(nBins, centeff[0], centeff[2], centeff[1], centeff[3]);
    g[iCase]->SetNameTitle(Form("g%d",iCase+1),sCase[iCase].Data());
    g[iCase]->SetLineColor(cCase[iCase]);
  }
  if (iBlock == 6 || iBlock == 8 || iBlock == 9 || iBlock == 10) {
    // 0-8GeV, 2.75-3
    ++iCase;
    sCase[iCase] = "0-8GeV, 2.75-3";
    cCase[iCase] = kMagenta-1;
    const Int_t nBins = 3;
    Double_t centeff[4][nBins] = {
      {10.,30.,65.},
      {10.,10.,25.},
      {0.115548,0.122216,0.124454},
      {0.001091,0.001154,0.000946}
    };
    g[iCase] = new TGraphErrors(nBins, centeff[0], centeff[2], centeff[1], centeff[3]);
    g[iCase]->SetNameTitle(Form("g%d",iCase+1),sCase[iCase].Data());
    g[iCase]->SetLineColor(cCase[iCase]);
  }
  if (iBlock == 6 || iBlock == 8 || iBlock == 9 || iBlock == 10) {
    // 0-8GeV, 3-3.25
    ++iCase;
    sCase[iCase] = "0-8GeV, 3-3.25";
    cCase[iCase] = kMagenta-2;
    const Int_t nBins = 3;
    Double_t centeff[4][nBins] = {
      {10.,30.,65.},
      {10.,10.,25.},
      {0.174683,0.181970,0.189288},
      {0.001375,0.001422,0.001149}
    };
    g[iCase] = new TGraphErrors(nBins, centeff[0], centeff[2], centeff[1], centeff[3]);
    g[iCase]->SetNameTitle(Form("g%d",iCase+1),sCase[iCase].Data());
    g[iCase]->SetLineColor(cCase[iCase]);
  }
  if (iBlock == 6 || iBlock == 8 || iBlock == 9 || iBlock == 10) {
    // 0-8GeV, 3.25-3.5
    ++iCase;
    sCase[iCase] = "0-8GeV, 3.25-3.5";
    cCase[iCase] = kMagenta-3;
    const Int_t nBins = 3;
    Double_t centeff[4][nBins] = {
      {10.,30.,65.},
      {10.,10.,25.},
      {0.182387,0.192443,0.198140},
      {0.001439,0.001508,0.001215}
    };
    g[iCase] = new TGraphErrors(nBins, centeff[0], centeff[2], centeff[1], centeff[3]);
    g[iCase]->SetNameTitle(Form("g%d",iCase+1),sCase[iCase].Data());
    g[iCase]->SetLineColor(cCase[iCase]);
  }
  if (iBlock == 6 || iBlock == 8 || iBlock == 9 || iBlock == 10) {
    // 0-8GeV, 3.5-3.75
    ++iCase;
    sCase[iCase] = "0-8GeV, 3.5-3.75";
    cCase[iCase] = kMagenta-4;
    const Int_t nBins = 3;
    Double_t centeff[4][nBins] = {
      {10.,30.,65.},
      {10.,10.,25.},
      {0.144259,0.153308,0.159697},
      {0.001390,0.001467,0.001185}
    };
    g[iCase] = new TGraphErrors(nBins, centeff[0], centeff[2], centeff[1], centeff[3]);
    g[iCase]->SetNameTitle(Form("g%d",iCase+1),sCase[iCase].Data());
    g[iCase]->SetLineColor(cCase[iCase]);
  }
  if (iBlock == 6 || iBlock == 8 || iBlock == 9 || iBlock == 10) {
    // 0-8GeV, 3.75-4
    ++iCase;
    sCase[iCase] = "0-8GeV, 3.75-4";
    cCase[iCase] = kMagenta-5;
    const Int_t nBins = 3;
    Double_t centeff[4][nBins] = {
      {10.,30.,65.},
      {10.,10.,25.},
      {0.057257,0.062132,0.063068},
      {0.000985,0.001036,0.000831}
    };
    g[iCase] = new TGraphErrors(nBins, centeff[0], centeff[2], centeff[1], centeff[3]);
    g[iCase]->SetNameTitle(Form("g%d",iCase+1),sCase[iCase].Data());
    g[iCase]->SetLineColor(cCase[iCase]);
  }
  //-------------------------------
  if (iBlock == 7 || iBlock == 8) {
    // 0-8GeV, 2.5-3
    ++iCase;
    sCase[iCase] = "0-8GeV, 2.5-3";
    cCase[iCase] = kRed;
    const Int_t nBins = 7;
    Double_t centeff[4][nBins] = {
      {5.,15.,25.,35.,45.,55.,75.},
      {5.,5.,5.,5.,5.,5.,15.},
      {0.074456,0.076946,0.076984,0.077465,0.077522,0.077969,0.078285},
      {0.000881,0.000902,0.000890,0.000909,0.000915,0.000898,0.000527}
    };
    g[iCase] = new TGraphErrors(nBins, centeff[0], centeff[2], centeff[1], centeff[3]);
    g[iCase]->SetNameTitle(Form("g%d",iCase+1),sCase[iCase].Data());
    g[iCase]->SetLineColor(cCase[iCase]);
  }
  if (iBlock == 7 || iBlock == 8) {
    // 0-8GeV, 3-3.5
    ++iCase;
    sCase[iCase] = "0-8GeV, 3-3.5";
    cCase[iCase] = kRed+1;
    const Int_t nBins = 7;
    Double_t centeff[4][nBins] = {
      {5.,15.,25.,35.,45.,55.,75.},
      {5.,5.,5.,5.,5.,5.,15.},
      {0.180188,0.185859,0.185834,0.188193,0.190911,0.194141,0.193746},
      {0.001404,0.001418,0.001408,0.001421,0.001419,0.001448,0.000833}
    };
    g[iCase] = new TGraphErrors(nBins, centeff[0], centeff[2], centeff[1], centeff[3]);
    g[iCase]->SetNameTitle(Form("g%d",iCase+1),sCase[iCase].Data());
    g[iCase]->SetLineColor(cCase[iCase]);
  }
  if (iBlock == 7 || iBlock == 8) {
    // 0-8GeV, 3.5-4
    ++iCase;
    sCase[iCase] = "0-8GeV, 3.5-4";
    cCase[iCase] = kRed+2;
    const Int_t nBins = 7;
    Double_t centeff[4][nBins] = {
      {5.,15.,25.,35.,45.,55.,75.},
      {5.,5.,5.,5.,5.,5.,15.},
      {0.102125,0.104557,0.105520,0.111232,0.110738,0.110521,0.113163},
      {0.001203,0.001240,0.001226,0.001265,0.001246,0.001276,0.000751}
    };
    g[iCase] = new TGraphErrors(nBins, centeff[0], centeff[2], centeff[1], centeff[3]);
    g[iCase]->SetNameTitle(Form("g%d",iCase+1),sCase[iCase].Data());
    g[iCase]->SetLineColor(cCase[iCase]);
  }
  //###############################
  
  // fit function
  TF1 *fEffVsCent = new TF1("fEffVsCent",EffVsCent,0.,90.,3);
  TF1 *fELoss = new TF1("fELoss",ELoss,0.,90.,2);
  
  // display eff vs cent and fit them
  const Int_t nb = 8;
  Double_t p[nb][nCases[iBlock]][2];
  TMatrixD *cov[nb][nCases[iBlock]];
  TCanvas *c = new TCanvas();
  for (Int_t i = 0; i < nCases[iBlock]; ++i) {
    g[i]->Draw((i==0)?"ap":"p");
    fEffVsCent->SetLineColor(cCase[i]);
    TFitResultPtr r;
    if (b < -1.) {
      Double_t dummy, eff;
      g[i]->GetPoint(g[i]->GetN()-1,dummy,eff);
      fEffVsCent->SetParameters(eff,0.5,2.);
      fEffVsCent->SetParLimits(1,0.,10.);
      fEffVsCent->FixParameter(2,2.);
      g[i]->Fit(fEffVsCent,"RISNQ");
      fEffVsCent->ReleaseParameter(2);
      fEffVsCent->SetParLimits(2,1.,10.);
      r = g[i]->Fit(fEffVsCent,"RIS");
      p[0][i][0] = fEffVsCent->GetParameter(1);
      p[0][i][1] = fEffVsCent->GetParameter(2);
      cov[0][i] = new TMatrixD(r->GetCovarianceMatrix().GetSub(1,2,1,2));
    } else if (b < 0.) {
      Double_t dummy, eff;
      g[i]->GetPoint(g[i]->GetN()-1,dummy,eff);
      for (Int_t ib = 0; ib < nb; ++ib) {
        fEffVsCent->SetParameters(eff,0.5,2.);
        fEffVsCent->SetParLimits(1,0.,10.);
        fEffVsCent->FixParameter(2,2.+0.5*ib);
        r = g[i]->Fit(fEffVsCent,"RIS+");
        p[ib][i][0] = fEffVsCent->GetParameter(1);
        p[ib][i][1] = fEffVsCent->GetParameter(2);
        cov[ib][i] = new TMatrixD(r->GetCovarianceMatrix().GetSub(1,2,1,2));
      }
    } else {
      Double_t dummy, eff;
      g[i]->GetPoint(g[i]->GetN()-1,dummy,eff);
      fEffVsCent->SetParameters(eff,0.5,2.);
      fEffVsCent->SetParLimits(1,0.,10.);
      fEffVsCent->FixParameter(2,b);
      r = g[i]->Fit(fEffVsCent,"RIS");
      p[0][i][0] = fEffVsCent->GetParameter(1);
      p[0][i][1] = fEffVsCent->GetParameter(2);
      cov[0][i] = new TMatrixD(r->GetCovarianceMatrix().GetSub(1,2,1,2));
    }
  }
  g[0]->SetMinimum(0.05);
  g[0]->SetMaximum(0.3);
  
  // compute and display eff loss vs cent
  const Int_t nBins = 5;
  Double_t cent[2][nBins] = {
    {5.,20.,40.,60.,80.},
    {5.,10.,10.,10.,10.}
  };
  TGraphErrors *gl[nb][nCases[iBlock]];
  TCanvas *cl = new TCanvas();
  for (Int_t i = 0; i < nCases[iBlock]; ++i) {
    for (Int_t ib = 0; (b > -2. && b < 0. && ib < nb) || ib < 1; ++ib) {
      Double_t effLoss[2][nBins];
      for (Int_t j = 0; j < nBins; ++j) {
        effLoss[0][j] = fELoss->Integral(cent[0][j]-cent[1][j],cent[0][j]+cent[1][j],p[ib][i])/(2.*cent[1][j]);
        effLoss[1][j] = fELoss->IntegralError(cent[0][j]-cent[1][j],cent[0][j]+cent[1][j],p[ib][i],cov[ib][i]->GetMatrixArray())/(2.*cent[1][j]);
      }
      gl[ib][i] = new TGraphErrors(nBins, cent[0], effLoss[0], cent[1], effLoss[1]);
      gl[ib][i]->SetNameTitle((ib < 1)?Form("gl%d",i+1):Form("gl%d_%d",i+1,ib),sCase[i].Data());
      gl[ib][i]->SetLineColor(cCase[i]);
      gl[ib][i]->Draw((i == 0 && ib == 0)?"ap":"p");
    }
  }
  gl[0][0]->SetMinimum(0.);
  gl[0][0]->SetMaximum(0.2);
  
  // display Eloss vs case in most central from data
  TH1F *hEffLoss = new TH1F("hELoss","eff loss in most central;eff loss",10,0.,0.2);
  hEffLoss->Sumw2();
  TH1F *hEffLossVsCase = new TH1F("hELossVsCase","eff loss vs case in most central;case;eff loss",nCases[iBlock],0.5,nCases[iBlock]+0.5);
  hEffLossVsCase->Sumw2();
  for (Int_t i = 0; i < nCases[iBlock]; ++i) {
    Double_t dummy, effc, effcErr, effp, effpErr, effLoss, effLossErr;
    g[i]->GetPoint(0,dummy,effc);
    effcErr = g[i]->GetErrorY(0);
    g[i]->GetPoint(g[i]->GetN()-1,dummy,effp);
    effpErr = g[i]->GetErrorY(g[i]->GetN()-1);
    effLoss = 1. - effc/effp;
    effLossErr = (1.-effLoss)*TMath::Sqrt(effcErr*effcErr/effc/effc + effpErr*effpErr/effp/effp);
    hEffLoss->Fill(effLoss,1./effLossErr/effLossErr);
    hEffLossVsCase->SetBinContent(i+1,effLoss);
    hEffLossVsCase->SetBinError(i+1,effLossErr);
  }
  TCanvas *cEffLoss = new TCanvas("cELoss","eff loss in most central",10,10,800,400);
  cEffLoss->Divide(2,1);
  cEffLoss->cd(1);
  hEffLossVsCase->Draw();
  cEffLoss->cd(2);
  hEffLoss->Draw();
  
  // display Eloss vs case from fit
  TH1F *hEffLossFit[nBins];
  TH1F *hEffLossFitVsCase[nBins][nb];
  for (Int_t j = 0; j < nBins; ++j) {
    hEffLossFit[j] = new TH1F(Form("hELossFit%d%d",(Int_t)(cent[0][j]-cent[1][j]),(Int_t)(cent[0][j]+cent[1][j])),Form("eff loss from fit in %d-%d%%;eff loss",(Int_t)(cent[0][j]-cent[1][j]),(Int_t)(cent[0][j]+cent[1][j])),10,0.,0.2-0.06*j);
    hEffLossFit[j]->Sumw2();
    for (Int_t ib = 0; (b > -2 && b < 0. && ib < nb) || ib < 1; ++ib) {
      hEffLossFitVsCase[j][ib] = new TH1F(Form("hELossFitVsCase%d%d_%d",(Int_t)(cent[0][j]-cent[1][j]),(Int_t)(cent[0][j]+cent[1][j]),ib+1),Form("eff loss from fit vs case in %d-%d%%;case;eff loss",(Int_t)(cent[0][j]-cent[1][j]),(Int_t)(cent[0][j]+cent[1][j])),nCases[iBlock],0.5,nCases[iBlock]+0.5);
      hEffLossFitVsCase[j][ib]->Sumw2();
      hEffLossFitVsCase[j][ib]->SetLineColor(ib+1);
      Double_t dummy, effLoss, effLossErr;
      for (Int_t i = 0; i < nCases[iBlock]; ++i) {
        gl[ib][i]->GetPoint(j,dummy,effLoss);
        effLossErr = TMath::Max(gl[ib][i]->GetErrorY(j),0.00001);
        hEffLossFit[j]->Fill(effLoss,1./effLossErr/effLossErr);
//        hEffLossFit[j]->Fill(effLoss);
        hEffLossFitVsCase[j][ib]->SetBinContent(i+1,effLoss);
        hEffLossFitVsCase[j][ib]->SetBinError(i+1,effLossErr);
      }
    }
  }
  TCanvas *cEffLossFit = new TCanvas("cELossFit","eff loss from fit",10,10,1200,500);
  cEffLossFit->Divide(5,2);
  for (Int_t j = 0; j < nBins; ++j) {
    cEffLossFit->cd(j+1);
    hEffLossFitVsCase[j][0]->SetMinimum(0.);
    for (Int_t ib = 0; (b > -2 && b < 0. && ib < nb) || ib < 1; ++ib) hEffLossFitVsCase[j][ib]->Draw((ib < 1)?"":"same");
    cEffLossFit->cd(j+6);
    hEffLossFit[j]->Draw();
  }
  
}

//______________________________________________________________________
Double_t EffVsCent(Double_t *x, Double_t *p)
{
  /// efficiency evolution versus centrality
  Double_t eff100 = p[0];
  Double_t a = p[1];
  Double_t b = p[2];
  return eff100 * (1. - a*TMath::Power(1.-x[0]/100.,b));
}

//______________________________________________________________________
Double_t ELoss(Double_t *x, Double_t *p)
{
  /// efficiency loss versus centrality
  Double_t a = p[0];
  Double_t b = p[1];
  return a*TMath::Power(1.-x[0]/100.,b);
}



