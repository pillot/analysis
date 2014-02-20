void JpsiCalculationinPbPb()
{
  
  gStyle->SetFillColor(0);
  /*
// Calculation RAA Analyse Roberta i=0 0-80, i=1 0-10, i=2 10-20, i=3 20-40, i=4 40-80
  Double_t NjpsiR[5]={1113., 498., 217., 271., 126.}; 
  Double_t e_NjpsiR[5]={108.,82., 48., 45, 17.};
  Double_t NmbR[5]={8.54e6, 1.06e6, 1.09e6, 2.15e6, 4.24e6}; 
  Double_t e_NmbR[5];
  Double_t aeffR=0.187;
  Double_t e_aeffR=0.001/aeffR;
  */
  ////// pass2 results //////
  Double_t NmbR[5]={16932459., 2088200., 2097283., 4258151., 8488825.};
  Double_t e_NmbR[5];
  Double_t aeffR=0.1944;
  Double_t e_aeffR=0.001/aeffR;
  //CB tail fixed to pp LHC10g in all cases (α=1.15 n=1.59)
  //Double_t NjpsiR[5]={2298, 894, 499, 724, 221}; 
  //Double_t e_NjpsiR[5]={331, 146, 84, 85, 62};
  //CB tail fixed to pp LHC10g in all cases (α=1.15 n=1.59) fix to 0-80%
  Double_t NjpsiR[5]={2291, 949, 550, 595, 234}; 
  Double_t e_NjpsiR[5]={138, 104, 70, 61, 27};
  //CB tail fixed to pp “old” MC in all cases (α=1.15 n=3.6)
  //Double_t NjpsiR[5]={2090, 808, 447, 656, 201}; 
  //Double_t e_NjpsiR[5]={301, 133, 140, 121, 56};
  //CB tail fixed to pp “old” MC in all cases (α=1.15 n=3.6) fix to 0-80%
  //Double_t NjpsiR[5]={2097, 864, 498, 544, 214}; 
  //Double_t e_NjpsiR[5]={126, 94, 63, 55, 24};
  
// Calculation RAA Analyse Christophe i=0 0-80, i=1 0-10, i=2 10-20, i=3 20-40, i=4 40-80
  //Double_t NjpsiC[5]={993., 377., 245., 282., 136.}; 
  //Double_t e_NjpsiC[5]={129., 63., 91., 58., 18.};
  Double_t NmbC[5]={7.83e6, 0.98e6, 1.00e6, 1.97e6, 3.89e6}; 
  Double_t e_NmbC[5];
  Double_t aeffC=0.188;
  Double_t e_aeffC=0.001/aeffC;  
  
  ////// pass1 results //////
  //CB tail fixed to pp LHC10g in all cases (α=1.15 n=1.59 )
  //Double_t NjpsiC[5]={961, 344, 229, 263, 124}; 
  //Double_t e_NjpsiC[5]={127, 142, 94, 57, 18};
  //CB tail fixed to pp LHC10g in all cases (α=1.15 n=1.59 ) fix to 0-80%
  Double_t NjpsiC[5]={957, 363, 239, 231, 124}; 
  Double_t e_NjpsiC[5]={85, 63, 42, 34, 16};
  //CB tail fixed to pp “old” MC in all cases (α=1.15 n=3.6 )
  //Double_t NjpsiC[5]={1049, 381, 250, 280, 140}; 
  //Double_t e_NjpsiC[5]={228, 150, 108, 59, 20};
  //CB tail fixed to pp “old” MC in all cases (α=1.15 n=3.6 ) fix to 0-80%
  //Double_t NjpsiC[5]={1046, 399, 259, 252, 136}; 
  //Double_t e_NjpsiC[5]={93, 69, 46, 37, 18};
  

  Double_t TAA[5]={7.04, 23.48, 14.43, 6.86, 1.20}; // mb^-1  
  Double_t e_TAA[5]={0.27, 0.97, 0.57, 0.28, 0.07}; //

  Double_t TAArat[5]={1., 19.49, 11.98, 5.69, 1.};
  Double_t e_TAArat[5]={0., 1.21, 0.63, 0.20, 0.};
  
  
  for(Int_t i=0; i<5; i++) {    
    e_NjpsiR[i]=e_NjpsiR[i]/NjpsiR[i];
    e_NmbR[i]=1./TMath::Sqrt(NmbR[i]);  
    e_NjpsiC[i]=e_NjpsiC[i]/NjpsiC[i];
    e_NmbC[i]=1./TMath::Sqrt(NmbC[i]);
    e_TAA[i]=e_TAA[i]/TAA[i];
    e_TAArat[i]=e_TAArat[i]/TAArat[i];
  }
  
  
  // Systematic error
  Double_t e_bck[6] = {0.085, 0.085, 0.085, 0.085, 0.085, 0.};
  Double_t e_signal[6] = {0.06, 0.06, 0.06, 0.06, 0.06, 0.06};
  Double_t e_gen[6] = {0., 0., 0., 0., 0., 0.02};
  Double_t e_trkUncorr[6] = {0., 0., 0., 0., 0., 0.05};
  Double_t e_trkCorr[6] = {0., 0., 0., 0., 0., 0.035};
  Double_t e_trkRes[6] = {0., 0., 0., 0., 0., 0.01};
  Double_t e_trkCent[6] = {0.01, 0.02, 0.01, 0.005, 0., 0.};
  Double_t e_trg[6] = {0., 0., 0., 0., 0., 0.04};
  
  
  Double_t sigmajpsipp = 4.1; //mub
  Double_t e_sigmajpsipp=0.6/sigmajpsipp;
  Double_t BR=0.059;
  Double_t e_BR=0.01;

  Double_t JpsiYieldPbPbR[5], RAAR[5], e_RAAR[5], RCPR[5], e_RCPR[5];
  Double_t JpsiYieldPbPbC[5], RAAC[5], e_RAAC[5], RCPC[5], e_RCPC[5];

  for(i=0; i<5; i++) { 
    
    JpsiYieldPbPbR[i] = (NjpsiR[i])/(BR*aeffR*NmbR[i]);
    
    JpsiYieldPbPbC[i] = (NjpsiC[i])/(BR*aeffC*NmbC[i]);
    
    RAAR[i]= JpsiYieldPbPbR[i]/sigmajpsipp*1000./TAA[i]; // 1000. factor mub to mb 
    
    e_RAAR[i] = TMath::Sqrt(e_signal[i]*e_signal[i] +
			    e_bck[i]*e_bck[i] +
			    e_trkCent[i]*e_trkCent[i] +
			    e_NmbR[i]*e_NmbR[i] +
			    e_TAA[i]*e_TAA[i]);
    
    RAAC[i]= JpsiYieldPbPbC[i]/sigmajpsipp*1000./TAA[i]; // 1000. factor mub to mb 
    
    e_RAAC[i] = TMath::Sqrt(e_signal[i]*e_signal[i] +
			    e_bck[i]*e_bck[i] +
			    e_trkCent[i]*e_trkCent[i] +
			    e_NmbC[i]*e_NmbC[i] +
			    e_TAA[i]*e_TAA[i]);
    
    cout << "Roberta RAA["<<i<<"] is " << RAAR[i] << "+-" << e_NjpsiR[i]*RAAR[i] << "+-" << e_RAAR[i]*RAAR[i] << " and ";
    cout << "Christophe RAA["<<i<<"] is " << RAAC[i] << "+-" << e_NjpsiC[i]*RAAC[i] << "+-" << e_RAAC[i]*RAAC[i] << endl;

  }
  
  cout << endl;
  
  for(i=0; i<5; i++) { 
    
    RCPR[i]= JpsiYieldPbPbR[i]/JpsiYieldPbPbR[4]/TAArat[i];
    
    e_RCPR[i] = TMath::Sqrt(e_signal[i]*e_signal[i] +
			    e_bck[i]*e_bck[i] +
			    e_trkCent[i]*e_trkCent[i] +
			    e_NmbR[i]*e_NmbR[i] +
			    e_TAArat[i]*e_TAArat[i]);
    
    RCPC[i]= JpsiYieldPbPbC[i]/JpsiYieldPbPbC[4]/TAArat[i];
    
    e_RCPC[i] = TMath::Sqrt(e_signal[i]*e_signal[i] +
			    e_bck[i]*e_bck[i] +
			    e_trkCent[i]*e_trkCent[i] +
			    e_NmbC[i]*e_NmbC[i] +
			    e_TAArat[i]*e_TAArat[i]);
    
    cout << "Roberta RCP["<<i<<"] is " << RCPR[i] << "+-" << e_NjpsiR[i]*RCPR[i] << "+-" << e_RCPR[i]*RCPR[i] << " and ";
    cout << "Christophe RCP["<<i<<"] is " << RCPC[i] << "+-" << e_NjpsiC[i]*RCPC[i] << "+-" << e_RCPC[i]*RCPC[i] << endl;
    
  }
  
  // correlated error for RAA
  Double_t e_RAARCorr = TMath::Sqrt(e_signal[5]*e_signal[5] +
				    e_gen[5]*e_gen[5] +
				    e_trkUncorr[5]*e_trkUncorr[5] +
				    e_trkCorr[5]*e_trkCorr[5] +
				    e_trkRes[5]*e_trkRes[5] +
				    e_trg[5]*e_trg[5] +
				    e_BR*e_BR +
				    e_sigmajpsipp[5]*e_sigmajpsipp[5]);  
  
  cout << "correlated RAA error = " << e_RAARCorr*100. << "%" << endl;
  
  // plot RAA
  TString label[5] = {"0-80%", "0-10%", "10-20%", "20-40%", "40-80%"};
  Double_t x[5] = {100., 5., 15., 30., 60.};
  Double_t ex[5] = {1., 5., 5., 10., 20.};
  TGraphErrors *gRAAR = new TGraphErrors(5);
  TGraphErrors *gRAARSys = new TGraphErrors(5);
  TGraphErrors *gRAAC = new TGraphErrors(5);
  TGraphErrors *gRAACSys = new TGraphErrors(5);
  for(i=4; i>=1; i--) {
    gRAAR->SetPoint(4-i,100.-x[i],RAAR[i]);
    gRAAR->SetPointError(4-i,ex[i],e_NjpsiR[i]*RAAR[i]);
    gRAARSys->SetPoint(4-i,100.-x[i],RAAR[i]);
    gRAARSys->SetPointError(4-i,ex[i],e_RAAR[i]*RAAR[i]);
    gRAAC->SetPoint(4-i,100.-x[i],RAAC[i]);
    gRAAC->SetPointError(4-i,ex[i],e_NjpsiC[i]*RAAC[i]);
    gRAACSys->SetPoint(4-i,100.-x[i],RAAC[i]);
    gRAACSys->SetPointError(4-i,ex[i],e_RAAC[i]*RAAC[i]);
  }
  
  gRAAR->GetXaxis()->Set(16., 21., 101.);
  gRAAC->GetXaxis()->Set(16., 21., 101.);
  for(i=4; i>=1; i--) {
    gRAAR->GetXaxis()->SetBinLabel(gRAAR->GetXaxis()->FindBin(100-x[i]), label[i].Data());
    gRAAC->GetXaxis()->SetBinLabel(gRAAC->GetXaxis()->FindBin(100-x[i]), label[i].Data());
  }
  
  new TCanvas("RAAR1","RAAR");
  gRAAR->GetYaxis()->SetTitle("RAA");
  gRAAR->GetYaxis()->SetRangeUser(0.2, 1.);
  gRAAR->SetMarkerStyle(20);
  gRAAR->SetMarkerColor(2);
  gRAAR->SetLineColor(2);
  gRAAR->Draw("ap");
  gRAARSys->SetFillStyle(0);
  gRAARSys->SetLineColor(2);
  gRAARSys->Draw("e2");
  gRAAC->GetYaxis()->SetTitle("RAA");
  gRAAC->GetYaxis()->SetRangeUser(0.2, 1.);
  gRAAC->SetMarkerStyle(21);
  gRAAC->SetMarkerColor(4);
  gRAAC->SetLineColor(4);
  gRAAC->Draw("p");
  gRAACSys->SetLineColor(4);
  gRAACSys->SetFillStyle(0);
  gRAACSys->Draw("e2");
  
  new TCanvas("RAAC1","RAAC");
  gRAAC->GetYaxis()->SetTitle("RAA");
  gRAAC->GetYaxis()->SetRangeUser(0.2, 1.);
  gRAAC->SetMarkerStyle(21);
  gRAAC->SetMarkerColor(4);
  gRAAC->SetLineColor(4);
  gRAAC->Draw("ap");
  gRAACSys->SetLineColor(4);
  gRAACSys->SetFillStyle(0);
  gRAACSys->Draw("e2");
  
  // plot RCP
  TGraphErrors *gRCPR = new TGraphErrors(4);
  TGraphErrors *gRCPRSys = new TGraphErrors(4);
  TGraphErrors *gRCPC = new TGraphErrors(4);
  TGraphErrors *gRCPCSys = new TGraphErrors(4);
  for(i=4; i>=1; i--) {
    gRCPR->SetPoint(4-i,100.-x[i],RCPR[i]);
    if (i == 4) gRCPR->SetPointError(4-i,ex[i],0.);
    else gRCPR->SetPointError(4-i,ex[i],TMath::Sqrt(e_NjpsiR[i]*e_NjpsiR[i]+e_NjpsiR[4]*e_NjpsiR[4])*RCPR[i]);
    gRCPRSys->SetPoint(4-i,100.-x[i],RCPR[i]);
    gRCPRSys->SetPointError(4-i,ex[i],e_RCPR[i]*RCPR[i]);
    gRCPC->SetPoint(4-i,100.-x[i],RCPC[i]);
    if (i == 4) gRCPC->SetPointError(4-i,ex[i],0.);
    else gRCPC->SetPointError(4-i,ex[i],TMath::Sqrt(e_NjpsiC[i]*e_NjpsiC[i]+e_NjpsiC[4]*e_NjpsiC[4])*RCPC[i]);
    gRCPCSys->SetPoint(4-i,100.-x[i],RCPC[i]);
    gRCPCSys->SetPointError(4-i,ex[i],e_RCPC[i]*RCPC[i]);
  }
  
  gRCPR->GetXaxis()->Set(16., 21., 101.);
  gRCPC->GetXaxis()->Set(16., 21., 101.);
  for(i=4; i>=1; i--) {
    gRCPR->GetXaxis()->SetBinLabel(gRCPR->GetXaxis()->FindBin(100-x[i]), label[i].Data());
    gRCPC->GetXaxis()->SetBinLabel(gRCPC->GetXaxis()->FindBin(100-x[i]), label[i].Data());
  }
  
  new TCanvas("RCPR1","RCPR");
  gRCPR->GetYaxis()->SetTitle("RCP");
  gRCPR->GetYaxis()->SetRangeUser(0., 1.2);
  gRCPR->SetMarkerStyle(20);
  gRCPR->SetMarkerColor(2);
  gRCPR->SetLineColor(2);
  gRCPR->Draw("ap");
  gRCPRSys->SetFillStyle(0);
  gRCPRSys->SetLineColor(2);
  gRCPRSys->Draw("e2");
  gRCPC->GetYaxis()->SetTitle("RCP");
  gRCPC->GetYaxis()->SetRangeUser(0., 1.2);
  gRCPC->SetMarkerStyle(21);
  gRCPC->SetMarkerColor(4);
  gRCPC->SetLineColor(4);
  gRCPC->Draw("p");
  gRCPCSys->SetLineColor(4);
  gRCPCSys->SetFillStyle(0);
  gRCPCSys->Draw("e2");
  
  new TCanvas("RCPC1","RCPC");
  gRCPC->GetYaxis()->SetTitle("RCP");
  gRCPC->GetYaxis()->SetRangeUser(0., 1.2);
  gRCPC->SetMarkerStyle(21);
  gRCPC->SetMarkerColor(4);
  gRCPC->SetLineColor(4);
  gRCPC->Draw("ap");
  gRCPCSys->SetLineColor(4);
  gRCPCSys->SetFillStyle(0);
  gRCPCSys->Draw("e2");
  
  
  // cout << "Statistical error: " << e_Njpsi*100. << "%" << endl;
  // cout << "Aeffi error:" << e_aeffi*100. << "%" << endl;
  // cout << "Signal extraction error:" << e_signalExtraction*100. << "%" << endl;
  // cout << "Bck subtraction error:" << e_bcksubtraction*100. << "%" << endl;
  // cout << "pp reference error:" << e_sigmajpsipp*100. << "%" << endl;
  // cout << "BR error:" << e_BR*100. << "%" << endl;
  
  
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
  
}
