/*
 *  JpsiCalculationinPbPbPass2.C
 *  aliroot_dev
 *
 *  Created by philippe pillot on 04/05/11.
 *  Copyright 2011 SUBATECH. All rights reserved.
 *
 */


void ALICEseal(TString type, Double_t xPad, Double_t yPad);

void JpsiCalculationinPbPbPass2()
{
  
  gStyle->SetFillColor(0);
  /*
  // Calculation RAA Analyse Christophe i=0 0-80, i=1 0-10, i=2 10-20, i=3 20-40, i=4 40-60, i=5 60-80
  const Int_t nCentBins = 6;
  Double_t NjpsiC[nCentBins]={2089., 862., 499., 540., 170., 34.}; 
  Double_t e_NjpsiC[nCentBins]={125., 97., 63., 52., 26., 8.};
  Double_t NmbC[nCentBins]={16932459., 2088200., 2097283., 4258151., 4248453., 4240372.};
  Double_t e_NmbC[nCentBins];
  
  Double_t TAA[nCentBins]={7.04, 23.48, 14.43, 6.86, 2.00, 0.42}; // mb^-1  
  Double_t e_TAA[nCentBins]={0.27, 0.97, 0.57, 0.28, 0.11, 0.03}; //
  
  Double_t TAArat[nCentBins]={1., 56.02, 34.43, 16.36, 4.78, 1.};
  Double_t e_TAArat[nCentBins]={0., 4.19, 2.21, 0.86, 0.18, 0.};
  
  // Systematic error
  Double_t e_bck[nCentBins+1] = {0.085, 0.085, 0.085, 0.085, 0.085, 0.085, 0.};
  Double_t e_signal[nCentBins+1] = {0.06, 0.06, 0.06, 0.06, 0.06, 0.06, 0.06};
  Double_t e_gen[nCentBins+1] = {0., 0., 0., 0., 0., 0., 0.02};
  Double_t e_trkUncorr[nCentBins+1] = {0., 0., 0., 0., 0., 0., 0.015};
  Double_t e_trkCorr[nCentBins+1] = {0., 0., 0., 0., 0., 0., 0.035};
  Double_t e_trkRes[nCentBins+1] = {0., 0., 0., 0., 0., 0., 0.01};
  Double_t e_trkCent[nCentBins+1] = {0.01, 0.02, 0.01, 0.005, 0., 0., 0.};
  Double_t e_trg[nCentBins+1] = {0., 0., 0., 0., 0., 0., 0.04};
   
  TString label[nCentBins] = {"0-80%", "0-10%", "10-20%", "20-40%", "40-60%", "60-80%"};
  Double_t x[nCentBins] = {100., 5., 15., 30., 50., 70.};
  Double_t ex[nCentBins] = {1., 5., 5., 10., 10., 10.};
  */
  
  // Calculation RAA Analyse Christophe i=0 0-80, i=1 0-10, i=2 10-20, i=3 20-40, i=4 40-80
  const Int_t nCentBins = 5;
  //Double_t NjpsiC[nCentBins]={2096., 809., 454., 663., 210.}; // CB free
  //Double_t e_NjpsiC[nCentBins]={298., 128., 79., 78., 58.};
  Double_t NjpsiC[nCentBins]={2089., 862., 499., 540., 212.}; // CB fix 0-80
  Double_t e_NjpsiC[nCentBins]={125., 97., 63., 52., 27.};
  //Double_t NjpsiC[nCentBins]={2333., 966., 561., 593., 209.};  // CB fix 0-80 + mixing pol3
  //Double_t e_NjpsiC[nCentBins]={125., 92., 62., 51., 23.};
  //Double_t NjpsiC[nCentBins]={2386., 989., 572., 600., 249.};  // CB fix 0-80 + mixing pol1
  //Double_t e_NjpsiC[nCentBins]={128., 95., 64., 53., 23.};
  Double_t NmbC[nCentBins]={16932459., 2088200., 2097283., 4258151., 8488825.};
  Double_t e_NmbC[nCentBins];
  
  Double_t TAA[nCentBins]={7.04, 23.48, 14.43, 6.86, 1.20}; // mb^-1  
  Double_t e_TAA[nCentBins]={0.27, 0.97, 0.57, 0.28, 0.07}; //
  
  Double_t TAArat[nCentBins]={1., 19.49, 11.98, 5.69, 1.};
  Double_t e_TAArat[nCentBins]={0., 1.21, 0.63, 0.20, 0.};
  
  Double_t NPart[nCentBins] = {138.9, 356.5, 260.5, 157.3, 45.5};
  Double_t e_NPart[nCentBins] = {3.2, 3.6, 4.4, 3.4, 2.1};
  
  
  // Systematic error
//  Double_t e_bck[nCentBins+1] = {0.085, 0.085, 0.085, 0.085, 0.085, 0.};
//  Double_t e_signal[nCentBins+1] = {0.06, 0.06, 0.06, 0.06, 0.06, 0.06};
  Double_t e_bck[nCentBins+1] = {0., 0., 0., 0., 0., 0.};
  Double_t e_signal[nCentBins+1] = {0., 0., 0., 0., 0., 0.};
  Double_t e_gen[nCentBins+1] = {0., 0., 0., 0., 0., 0.02};
  Double_t e_trkUncorr[nCentBins+1] = {0., 0., 0., 0., 0., 0.03};
  Double_t e_trkCorr[nCentBins+1] = {0., 0., 0., 0., 0., 0.02};
  Double_t e_trkRes[nCentBins+1] = {0., 0., 0., 0., 0., 0.02};
  Double_t e_trkCent[nCentBins+1] = {0.01, 0.02, 0.01, 0.005, 0., 0.};
  Double_t e_trg[nCentBins+1] = {0., 0., 0., 0., 0., 0.04};
  
  Bool_t vsCentClass = kTRUE;
  TString label[nCentBins] = {"0-80%", "0-10%", "10-20%", "20-40%", "40-80%"};
  Double_t x[nCentBins] = {100., 5., 15., 30., 60.};
  Double_t ex[nCentBins] = {1., 5., 5., 10., 20.};
  
  
  for(Int_t i=0; i<nCentBins; i++) {    
    e_NjpsiC[i]=e_NjpsiC[i]/NjpsiC[i];
    e_NmbC[i]=1./TMath::Sqrt(NmbC[i]);
    e_TAA[i]=e_TAA[i]/TAA[i];
    e_TAArat[i]=e_TAArat[i]/TAArat[i];
  }
  
  Double_t aeffC=0.1944;
  Double_t e_aeffC=0.0005/aeffC;  
  
  Double_t sigmajpsipp = 4.1; //mub
  Double_t e_sigmajpsipp=0.6/sigmajpsipp;
  Double_t BR=0.059;
  Double_t e_BR=0.01;
  
  Double_t JpsiYieldPbPbC[nCentBins], RAAC[nCentBins], e_RAAC[nCentBins], RCPC[nCentBins], e_RCPC[nCentBins];
  
  //Double_t RAAC[nCentBins]={0.023, 0.014, 0.022, 0.066, 0.198}; // difference pass2 - pass1
  //Double_t e_RAAC[nCentBins]={0.046, 0.043, 0.131, 0.086, 0.063};
  
  for(i=0; i<nCentBins; i++) { 
    
    JpsiYieldPbPbC[i] = (NjpsiC[i])/(BR*aeffC*NmbC[i]);
    
    RAAC[i]= JpsiYieldPbPbC[i]/sigmajpsipp*1000./TAA[i]; // 1000. factor mub to mb 
    
    e_RAAC[i] = TMath::Sqrt(e_signal[i]*e_signal[i] +
			    e_bck[i]*e_bck[i] +
			    e_trkCent[i]*e_trkCent[i] +
			    e_NmbC[i]*e_NmbC[i] +
			    e_TAA[i]*e_TAA[i]);
    
    cout << "Christophe RAA["<<i<<"] is " << RAAC[i] << "+-" << e_NjpsiC[i]*RAAC[i] << "+-" << e_RAAC[i]*RAAC[i] << endl;
    
  }
  
  cout << endl;
  
  for(i=0; i<nCentBins; i++) { 
    
    RCPC[i]= JpsiYieldPbPbC[i]/JpsiYieldPbPbC[nCentBins-1]/TAArat[i];
    
    e_RCPC[i] = TMath::Sqrt(e_signal[i]*e_signal[i] +
			    e_bck[i]*e_bck[i] +
			    e_trkCent[i]*e_trkCent[i] +
			    e_NmbC[i]*e_NmbC[i] +
			    e_TAArat[i]*e_TAArat[i]);
    
    cout << "Christophe RCP["<<i<<"] is " << RCPC[i] << "+-" << e_NjpsiC[i]*RCPC[i] << "+-" << e_RCPC[i]*RCPC[i] << endl;
    
  }
  
  // correlated error for RAA
  Double_t e_RAARCorr = TMath::Sqrt(e_signal[nCentBins]*e_signal[nCentBins] +
				    e_gen[nCentBins]*e_gen[nCentBins] +
				    e_trkUncorr[nCentBins]*e_trkUncorr[nCentBins] +
				    e_trkCorr[nCentBins]*e_trkCorr[nCentBins] +
				    e_trkRes[nCentBins]*e_trkRes[nCentBins] +
				    e_trg[nCentBins]*e_trg[nCentBins] +
				    e_BR*e_BR +
				    e_sigmajpsipp[nCentBins]*e_sigmajpsipp[nCentBins]);
  
  cout << "correlated RAA error = " << e_RAARCorr*100. << "%" << endl;
  
  // plot RAA
  TGraphErrors *gRAAC = new TGraphErrors(nCentBins);
  TGraphErrors *gRAACSys = new TGraphErrors(nCentBins);
  for(i=nCentBins-1; i>=1; i--) {
    if (vsCentClass) {
      gRAAC->SetPoint(nCentBins-1-i,100.-x[i],RAAC[i]);
      gRAAC->SetPointError(nCentBins-1-i,ex[i],e_NjpsiC[i]*RAAC[i]);
      gRAACSys->SetPoint(nCentBins-1-i,100.-x[i],RAAC[i]);
      gRAACSys->SetPointError(nCentBins-1-i,ex[i],e_RAAC[i]*RAAC[i]);
    } else {
      gRAAC->SetPoint(nCentBins-1-i,NPart[i],RAAC[i]);
      gRAAC->SetPointError(nCentBins-1-i,e_NPart[i],e_NjpsiC[i]*RAAC[i]);
      gRAACSys->SetPoint(nCentBins-1-i,NPart[i],RAAC[i]);
      gRAACSys->SetPointError(nCentBins-1-i,e_NPart[i],e_RAAC[i]*RAAC[i]);
    }
  }
  
  if (vsCentClass) {
    gRAAC->GetXaxis()->Set(16, 21., 101.);
    for(i=nCentBins-1; i>=1; i--) {
      gRAAC->GetXaxis()->SetBinLabel(gRAAC->GetXaxis()->FindBin(100-x[i]), label[i].Data());
    }
  } else {
    gRAAC->GetXaxis()->Set(40, 0., 400.);
  }
  
  new TCanvas("RAAC","RAAC");
  if (vsCentClass) {
    gRAAC->GetXaxis()->SetLabelSize(0.06);
  } else {
    gRAAC->GetXaxis()->SetTitle("N_{part}");
    gRAAC->GetXaxis()->SetTitleSize(0.05);
    gRAAC->GetXaxis()->SetTitleOffset(0.9);
    gRAAC->GetXaxis()->SetLabelSize(0.05);
  }
  gRAAC->GetYaxis()->SetTitle("RAA");
  gRAAC->GetYaxis()->SetTitleSize(0.05);
  gRAAC->GetYaxis()->SetTitleOffset(0.9);
  gRAAC->GetYaxis()->SetLabelSize(0.05);
  gRAAC->GetYaxis()->SetRangeUser(0.2, 1.);
  gRAAC->SetMarkerStyle(21);
  gRAAC->SetMarkerColor(1);
  gRAAC->SetLineColor(1);
  gRAAC->Draw("ap");
  gRAACSys->SetLineColor(1);
  gRAACSys->SetFillStyle(0);
  gRAACSys->Draw("e2");
  ALICEseal("proposed as preliminary", 0.6, 0.6);
  
  // plot RCP
  TGraphErrors *gRCPR = new TGraphErrors(nCentBins-1);
  TGraphErrors *gRCPRSys = new TGraphErrors(nCentBins-1);
  TGraphErrors *gRCPC = new TGraphErrors(nCentBins-1);
  TGraphErrors *gRCPCSys = new TGraphErrors(nCentBins-1);
  for(i=nCentBins-1; i>=1; i--) {
    if (vsCentClass) {
      gRCPC->SetPoint(nCentBins-1-i,100.-x[i],RCPC[i]);
      gRCPC->SetPointError(nCentBins-1-i,ex[i],e_NjpsiC[i]*RCPC[i]);
      gRCPCSys->SetPoint(nCentBins-1-i,100.-x[i],RCPC[i]);
      gRCPCSys->SetPointError(nCentBins-1-i,ex[i],e_RCPC[i]*RCPC[i]);
    } else {
      gRCPC->SetPoint(nCentBins-1-i,NPart[i],RCPC[i]);
      gRCPC->SetPointError(nCentBins-1-i,e_NPart[i],e_NjpsiC[i]*RCPC[i]);
      gRCPCSys->SetPoint(nCentBins-1-i,NPart[i],RCPC[i]);
      gRCPCSys->SetPointError(nCentBins-1-i,e_NPart[i],e_RCPC[i]*RCPC[i]);
    }
  }
  
  if (vsCentClass) {
    gRCPC->GetXaxis()->Set(16., 21., 101.);
    for(i=nCentBins-1; i>=1; i--) {
      gRCPC->GetXaxis()->SetBinLabel(gRCPC->GetXaxis()->FindBin(100-x[i]), label[i].Data());
    }
  } else {
    gRCPC->GetXaxis()->Set(40., 0., 400.);
  }
  
  new TCanvas("RCPC","RCPC");
  if (vsCentClass) {
    gRCPC->GetXaxis()->SetLabelSize(0.06);
  } else {
    gRCPC->GetXaxis()->SetTitle("N_{part}");
    gRCPC->GetXaxis()->SetTitleSize(0.05);
    gRCPC->GetXaxis()->SetTitleOffset(0.9);
    gRCPC->GetXaxis()->SetLabelSize(0.05);
  }
  gRCPC->GetYaxis()->SetTitle("RCP");
  gRCPC->GetYaxis()->SetTitleSize(0.05);
  gRCPC->GetYaxis()->SetTitleOffset(0.9);
  gRCPC->GetYaxis()->SetLabelSize(0.05);
  gRCPC->GetYaxis()->SetRangeUser(0.5, 1.3);
  gRCPC->SetMarkerStyle(21);
  gRCPC->SetMarkerColor(1);
  gRCPC->SetLineColor(1);
  gRCPC->Draw("ap");
  gRCPCSys->SetLineColor(1);
  gRCPCSys->SetFillStyle(0);
  gRCPCSys->Draw("e2");
  ALICEseal("proposed as preliminary", 0.2, 0.2);
  
  
  // cout << "Statistical error: " << e_Njpsi*100. << "%" << endl;
  // cout << "Aeffi error:" << e_aeffi*100. << "%" << endl;
  // cout << "Signal extraction error:" << e_signalExtraction*100. << "%" << endl;
  // cout << "Bck subtraction error:" << e_bcksubtraction*100. << "%" << endl;
  // cout << "pp reference error:" << e_sigmajpsipp*100. << "%" << endl;
  // cout << "BR error:" << e_BR*100. << "%" << endl;
  
  /*
   // NJpsi from Christophe with different fitting technics
   // gauss all free
   Double_t gausFree[4] = {104,263,219,346};
   Double_t egausFree[4] = {15,49,48,105};
   
   // CB fixed tails to pp data
   Double_t cbFixedpp[4] = {136,282,245,377};
   Double_t ecbfixedpp[4] = {18,58,91,63};
   
   // CB extended, tails/sigma  fixed to embedding
   Double_t cb2FixedEmbed[4] = {129,254,239,373};
   Double_t ecb2FixedEmbed[4] = {17,37,43,62};
   
   // likeSign bckgd subracted, free gauss
   Double_t lsBkgGausFree[4] = {127,304,222,303};
   Double_t elsBkgGausFree[4] = {19,55,56,60};
   
   // plot
   TGraphErrors *ggausFree = new TGraphErrors(4);
   TGraphErrors *gcbFixedpp = new TGraphErrors(4);
   TGraphErrors *gcb2FixedEmbed = new TGraphErrors(4);
   TGraphErrors *glsBkgGausFree = new TGraphErrors(4);
   for(i=0; i<4; i++) {
   ggausFree->SetPoint(i,x[i+1],gausFree[3-i]);
   ggausFree->SetPointError(i,ex[i+1],egausFree[3-i]);
   gcbFixedpp->SetPoint(i,x[i+1],cbFixedpp[3-i]);
   gcbFixedpp->SetPointError(i,ex[i+1],ecbfixedpp[3-i]);
   gcb2FixedEmbed->SetPoint(i,x[i+1],cb2FixedEmbed[3-i]);
   gcb2FixedEmbed->SetPointError(i,ex[i+1],ecb2FixedEmbed[3-i]);
   glsBkgGausFree->SetPoint(i,x[i+1],lsBkgGausFree[3-i]);
   glsBkgGausFree->SetPointError(i,ex[i+1],elsBkgGausFree[3-i]);
   }
   ggausFree->GetXaxis()->Set(16., 1., 81.);
   for(i=0; i<4; i++) {
   ggausFree->GetXaxis()->SetBinLabel(ggausFree->GetXaxis()->FindBin(x[i+1]), label[i+1].Data());
   }
   TLegend *l = new TLegend(0.2,0.7,0.9,0.9);
   l->SetBorderSize(0);
   l->AddEntry(ggausFree, "gaus all free", "lp");
   l->AddEntry(gcbFixedpp, "CB fixed tails to pp", "lp");
   l->AddEntry(gcb2FixedEmbed, "CB2, tails/sigma fixed to embedding", "lp");
   l->AddEntry(glsBkgGausFree, "likeSign bkg subracted, free gaus", "lp");
   new TCanvas("Christophe","Christophe");
   ggausFree->GetYaxis()->SetTitle("raw N_{J/#psi}");
   ggausFree->SetMarkerColor(2);
   ggausFree->SetLineColor(2);
   ggausFree->Draw("ap");
   gcbFixedpp->SetMarkerColor(1);
   gcbFixedpp->SetLineColor(1);
   gcbFixedpp->Draw("p");
   gcb2FixedEmbed->SetMarkerColor(8);
   gcb2FixedEmbed->SetLineColor(8);
   gcb2FixedEmbed->Draw("p");
   glsBkgGausFree->SetMarkerColor(4);
   glsBkgGausFree->SetLineColor(4);
   glsBkgGausFree->Draw("p");
   l->Draw("same");
   */
}

//______________________________________________________
void ALICEseal(TString type, Double_t xPad, Double_t yPad)
{
  TVirtualPad* currPad = gPad;
  TPad *myPadLogo = new TPad("myPadLogo", "Pad for ALICE Logo",xPad,yPad,xPad+0.17,yPad+0.17);
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
  Double_t x2 = x1 + 0.3, y2 = y1 + 0.08;
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

