//
//  CompareJPsi.C
//  aliroot_dev
//
//  Created by philippe pillot on 09/07/2014.
//  Copyright (c) 2014 Philippe Pillot. All rights reserved.
//

Int_t size0, size1, size2;

//------------
TString label1 = "MULorMLL";
const Int_t nCent = 10;
TString centBinName[nCent] = {"010", "1020", "2030", "3040", "4050", "5060", "6070", "7080", "8090", "090"};
Double_t nMULMLLvsCent[nCent] = {
  16280.4,
  10391.7,
  6571.66,
  3640.80,
  2103.33,
  1253.36,
  661.988,
  272.671,
  139.588,
  41092.8
};
const Int_t nCent2 = 7;
TString centpTBinName[nCent2] = {"010", "1020", "2030", "3040", "4050", "5060", "6090"};
const Int_t npTCent = 3;
Float_t dpTCentLowEdge[npTCent] = {0., 2., 5.};
Float_t dpTCentUpEdge[npTCent] = {2., 5., 8.};
Double_t nMULMLLvsCentpT[npTCent][nCent2] = {{
  9497.41,
  5894.44,
  3426.81,
  1860.05,
  942.431,
  613.993,
  589.485
},{
  5658.82,
  3846.60,
  2665.49,
  1419.30,
  971.597,
  531.115,
  406.913
},{
  737.320,
  469.771,
  420.855,
  327.058,
  179.485,
  98.8554,
  74.4495
}};
const Int_t nCent3 = 5;
TString pTCentBinName[nCent3] = {"020", "2040", "4090", "040", "090"};
const Int_t npT = 7;
Float_t dpTLowEdge[npT] = {0., 1., 2., 3., 4., 5., 6.};
Float_t dpTUpEdge[npT] = {1., 2., 3., 4., 5., 6., 8.};
Double_t nMULMLLvspTCent[nCent3][npT] = {{
  6860.71,
  8779.06,
  5420.92,
  2704.56,
  1393.19,
  755.666,
  466.016
},{
  2166.04,
  3176.42,
  2156.06,
  1253.84,
  654.479,
  405.100,
  326.147
},{
  910.984,
  1257.23,
  944.857,
  617.405,
  375.280,
  183.967,
  156.992
},{
  9018.99,
  11938.2,
  7575.45,
  3971.67,
  2046.77,
  1150.02,
  787.476
},{
  9909.05,
  13203.3,
  8514.03,
  4575.99,
  2424.80,
  1328.42,
  947.906
}};
const Int_t nCent4 = 7;
TString centyBinName[nCent4] = {"010", "1020", "2030", "3040", "4050", "5060", "6090"};
const Int_t nyCent = 3;
Float_t dyCentLowEdge[nyCent] = {2.5, 3., 3.5};
Float_t dyCentUpEdge[nyCent] = {3., 3.5, 4.};
Double_t nMULMLLvsCenty[nyCent][nCent4] = {{
  4776.87,
  3467.73,
  1866.68,
  926.026,
  609.347,
  361.456,
  293.445
},{
  8291.93,
  5173.62,
  3358.74,
  1966.56,
  1039.94,
  650.841,
  523.483
},{
  3217.76,
  1716.09,
  1327.12,
  721.111,
  440.992,
  243.100,
  247.507
}};
const Int_t nCent5 = 4;
TString yCentBinName[nCent5] = {"020", "2040", "4090", "090"};
const Int_t ny = 6;
Float_t dyLowEdge[ny] = {2.5, 2.75, 3., 3.25, 3.5, 3.75};
Float_t dyUpEdge[ny] = {2.75, 3., 3.25, 3.5, 3.75, 4.};
Double_t nMULMLLvsyCent[nCent5][ny] = {{
  1960.10,
  6288.71,
  7579.61,
  5909.01,
  3901.25,
  1023.08
},{
  624.713,
  2194.74,
  2703.82,
  2639.86,
  1546.76,
  494.736
},{
  296.810,
  964.806,
  1115.74,
  1103.95,
  710.298,
  215.089
},{
  2870.52,
  9443.88,
  11376.8,
  9632.29,
  6184.74,
  1755.75
}};

//------------
TString label2 = "MULorMLL_TrgSign";
Double_t nMULMLLTrgvsCent[nCent] = {
  16161.8,
  10248.4,
  6436.96,
  3632.41,
  2090.47,
  1231.62,
  656.758,
  267.461,
  137.795,
  40827.0
};
Double_t nMULMLLTrgvsCentpT[npTCent][nCent2] = {{
  9430.03,
  5835.19,
  3352.28,
  1849.65,
  943.205,
  601.370,
  586.665
},{
  5612.85,
  3777.84,
  2614.18,
  1423.67,
  963.598,
  527.657,
  402.418
},{
  740.558,
  457.465,
  414.647,
  326.501,
  169.878,
  96.7483,
  70.8339
}};
Double_t nMULMLLTrgvspTCent[nCent3][npT] = {{
  6865.93,
  8651.77,
  5345.26,
  2671.71,
  1386.13,
  752.327,
  465.500
},{
  2138.95,
  3122.39,
  2133.64,
  1236.26,
  646.991,
  395.370,
  324.105
},{
  909.218,
  1240.99,
  929.541,
  612.294,
  375.936,
  175.052,
  153.691
},{
  8994.88,
  11757.2,
  7480.64,
  3923.09,
  2029.05,
  1137.93,
  785.095
},{
  9879.98,
  13009.8,
  8408.78,
  4522.12,
  2404.31,
  1308.67,
  939.232
}};
Double_t nMULMLLTrgvsCenty[nyCent][nCent4] = {{
  4787.20,
  3420.85,
  1862.44,
  917.159,
  609.272,
  354.907,
  291.467
},{
  8199.88,
  5110.07,
  3277.38,
  1974.67,
  1035.46,
  638.825,
  519.525
},{
  3185.13,
  1676.87,
  1282.79,
  713.634,
  433.575,
  243.776,
  242.107
}};
Double_t nMULMLLTrgvsyCent[nCent5][ny] = {{
  1983.01,
  6231.08,
  7544.09,
  5789.81,
  3818.43,
  1018.83
},{
  626.176,
  2180.96,
  2664.45,
  2598.11,
  1517.88,
  471.495
},{
  293.490,
  961.403,
  1108.70,
  1090.22,
  702.002,
  211.168
},{
  2876.09,
  9379.72,
  11288.1,
  9465.57,
  6072.00,
  1722.81
}};

void PrintStat(TH1F *h1, TH1F *h2, TString label);
void DrawStat(Double_t n1, Double_t n2, TH1F *d, Int_t bin, TString label);

void CompareJPsi()
{
  /// compare JPsi results
  
  size1 = TMath::Max(9,label1.Length());
  size2 = TMath::Max(9,label2.Length());
  
  gROOT->SetStyle("Plain");
  gROOT->ForceStyle();
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  
  TFile* outfile = TFile::Open(Form("DiffJPsi_%s_%s.root",label1.Data(),label2.Data()),"RECREATE");
  
  // ###### plots versus centrality ######
  TH1F *hDiffVsCent = new TH1F("hDiffVsCent","diff vs cent",nCent,0.,nCent);
  TString label0 = "centrality";
  size0 = TMath::Max(9,label0.Length());
  printf("\nstat versus centrality\n");
  printf(Form("%%%ds %%%ds %%%ds       diff\n",size0,size1,size2),label0.Data(),label1.Data(),label2.Data());
  for (Int_t i=0; i<nCent; i++) {
    PrintStat(nMULMLLvsCent[i], nMULMLLTrgvsCent[i], centBinName[i].Data());
    DrawStat(nMULMLLvsCent[i], nMULMLLTrgvsCent[i], hDiffVsCent, i+1, centBinName[i].Data());
  }
  TCanvas* cDiffVsCent = new TCanvas("cDiffVsCent", "diff vs cent", 800, 400);
  hDiffVsCent->Draw();
  outfile->cd();
  hDiffVsCent->Write();
  
  // ###### plots versus centrality in pT bins ######
  for (Int_t ipT = 0; ipT < npTCent; ++ipT) {
    TH1F *hDiffVsCentpT = new TH1F(Form("hDiffVsCentpT%d",ipT+1),Form("diff vs cent in %3.1f < pT < %3.1f",dpTCentLowEdge[ipT],dpTCentUpEdge[ipT]),nCent2,0.,nCent2);
    label0 = "centrality";
    size0 = TMath::Max(9,label0.Length());
    printf("\nstat versus centrality in %3.1f < pT < %3.1f\n",dpTCentLowEdge[ipT],dpTCentUpEdge[ipT]);
    printf(Form("%%%ds %%%ds %%%ds       diff\n",size0,size1,size2),label0.Data(),label1.Data(),label2.Data());
    for (Int_t i=0; i<nCent2; i++) {
      PrintStat(nMULMLLvsCentpT[ipT][i], nMULMLLTrgvsCentpT[ipT][i], centpTBinName[i].Data());
      DrawStat(nMULMLLvsCentpT[ipT][i], nMULMLLTrgvsCentpT[ipT][i], hDiffVsCentpT, i+1, centpTBinName[i].Data());
    }
    TCanvas* cDiffVsCentpT = new TCanvas(Form("cDiffVsCentpT%d",ipT), Form("diff vs cent in %3.1f < pT < %3.1f",dpTCentLowEdge[ipT],dpTCentUpEdge[ipT]), 800, 400);
    hDiffVsCentpT->Draw();
    outfile->cd();
    hDiffVsCentpT->Write();
  }
  
  // ###### plots versus pT in centrality bins ######
  for (Int_t iCent = 0; iCent < nCent3; ++iCent) {
    TH1F *hDiffVspTCent = new TH1F(Form("hDiffVspTCent%d",iCent),Form("diff vs pT in cent %s%%",pTCentBinName[iCent].Data()),npT,0.,npT);
    label0 = "pT";
    size0 = TMath::Max(9,label0.Length());
    printf("\nstat versus pT in cent %s%%\n",pTCentBinName[iCent].Data());
    printf(Form("%%%ds %%%ds %%%ds       diff\n",size0,size1,size2),label0.Data(),label1.Data(),label2.Data());
    for (Int_t ipT = 0; ipT < npT; ++ipT) {
      PrintStat(nMULMLLvspTCent[iCent][ipT], nMULMLLTrgvspTCent[iCent][ipT], Form("%3.1f-%3.1f",dpTLowEdge[ipT],dpTUpEdge[ipT]));
      DrawStat(nMULMLLvspTCent[iCent][ipT], nMULMLLTrgvspTCent[iCent][ipT], hDiffVspTCent, ipT+1, Form("%3.1f-%3.1f",dpTLowEdge[ipT],dpTUpEdge[ipT]));
    }
    TCanvas* cDiffVspTCent = new TCanvas(Form("cDiffVspTCent%d",iCent), Form("diff vs pT in cent %s%%",pTCentBinName[iCent].Data()), 800, 400);
    hDiffVspTCent->Draw();
    outfile->cd();
    hDiffVspTCent->Write();
  }
  
  // ###### plots versus centrality in y bins ######
  for (Int_t iy = 0; iy < nyCent; ++iy) {
    TH1F *hDiffVsCenty = new TH1F(Form("hDiffVsCenty%d",iy+1),Form("diff vs cent in %4.2f < y < %4.2f",dyCentLowEdge[iy],dyCentUpEdge[iy]),nCent4,0.,nCent4);
    label0 = "centrality";
    size0 = TMath::Max(9,label0.Length());
    printf("\nstat versus centrality in %4.2f < y < %4.2f\n",dyCentLowEdge[iy],dyCentUpEdge[iy]);
    printf(Form("%%%ds %%%ds %%%ds       diff\n",size0,size1,size2),label0.Data(),label1.Data(),label2.Data());
    for (Int_t i=0; i<nCent4; i++) {
      PrintStat(nMULMLLvsCenty[iy][i], nMULMLLTrgvsCenty[iy][i], centyBinName[i].Data());
      DrawStat(nMULMLLvsCenty[iy][i], nMULMLLTrgvsCenty[iy][i], hDiffVsCenty, i+1, centyBinName[i].Data());
    }
    TCanvas* cDiffVsCenty = new TCanvas(Form("cDiffVsCenty%d",iy), Form("diff vs cent in %4.2f < y < %4.2f",dyCentLowEdge[iy],dyCentUpEdge[iy]), 800, 400);
    hDiffVsCenty->Draw();
    outfile->cd();
    hDiffVsCenty->Write();
  }
  
  // ###### plots versus y in centrality bins ######
  for (Int_t iCent = 0; iCent < nCent5; ++iCent) {
    TH1F *hDiffVsyCent = new TH1F(Form("hDiffVsyCent%d",(iCent<3)?iCent:4),Form("diff vs y in cent %s%%",yCentBinName[iCent].Data()),ny,0.,ny);
    label0 = "y";
    size0 = TMath::Max(9,label0.Length());
    printf("\nstat versus y in cent %s%%\n",yCentBinName[iCent].Data());
    printf(Form("%%%ds %%%ds %%%ds       diff\n",size0,size1,size2),label0.Data(),label1.Data(),label2.Data());
    for (Int_t iy = 0; iy < ny; ++iy) {
      PrintStat(nMULMLLvsyCent[iCent][iy], nMULMLLTrgvsyCent[iCent][iy], Form("%4.2f-%4.2f",dyLowEdge[iy],dyUpEdge[iy]));
      DrawStat(nMULMLLvsyCent[iCent][iy], nMULMLLTrgvsyCent[iCent][iy], hDiffVsyCent, iy+1, Form("%4.2f-%4.2f",dyLowEdge[iy],dyUpEdge[iy]));
    }
    TCanvas* cDiffVsyCent = new TCanvas(Form("cDiffVsyCent%d",iCent), Form("diff vs y in cent %s%%",yCentBinName[iCent].Data()), 800, 400);
    hDiffVsyCent->Draw();
    outfile->cd();
    hDiffVsyCent->Write();
  }
  
}

void PrintStat(Double_t n1, Double_t n2, TString label)
{
  /// print number of JPsi
  
  Double_t diff = 0.;
  if (n1 > 0.) diff = (n2-n1)/n1;
  Double_t err = 1.;
  if (TMath::Abs(diff) > 0. && TMath::Abs(diff) < 1.) err = TMath::Max(1./n1,TMath::Sqrt(TMath::Abs(diff*(1.-TMath::Abs(diff)))/n1));
  else if (TMath::Abs(diff) > 1.) err = TMath::Abs(diff);
  printf(Form("%%%ds %%%dd %%%dd   %%4.2f ± %%4.2f %%%%\n",size0,size1,size2), label.Data(), (Int_t)n1, (Int_t)n2, 100.*diff, 100.*err);
  
}

void DrawStat(Double_t n1, Double_t n2, TH1F *d, Int_t bin, TString label)
{
  /// print number of JPsi
  
  Double_t diff = 0.;
  if (n1 > 0.) diff = (n2-n1)/n1;
  Double_t err = 1.;
  if (TMath::Abs(diff) > 0. && TMath::Abs(diff) < 1.) err = TMath::Max(1./n1,TMath::Sqrt(TMath::Abs(diff*(1.-TMath::Abs(diff)))/n1));
  else if (TMath::Abs(diff) > 1.) err = TMath::Abs(diff);
  d->SetBinContent(bin, 100.*diff);
  d->SetBinError(bin, 100.*err);
  d->GetXaxis()->SetBinLabel(bin, label.Data());
  
}
