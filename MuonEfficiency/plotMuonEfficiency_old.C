#include <TH1.h>
#include <TH2.h>
#include <TAxis.h>
#include <TString.h>
#include <Riostream.h>
#include <TFile.h>
#include <TList.h>
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include <TMath.h>
#include <TArrayD.h>
#include <TStyle.h>

void plotMuonEfficiencyVsCent();
void plotMuonEfficiencyVsRun(TString runList = "");
void integratedEfficiency(Double_t minCent = 0., Double_t maxCent = 80.);
void computeTrackingEfficiency(TArrayD &chambersEff, TArrayD *chambersEffErr, Bool_t print = kFALSE);

//---------------------------------------------------------------------------
void plotMuonEfficiency_old(TString runList = "")
{
  /// compute the tracking efficiency versus centrality and versus run if runList != ""
  /// !!! to be compiled because of CINT !!!
  
  plotMuonEfficiencyVsCent();
  plotMuonEfficiencyVsRun(runList);
  integratedEfficiency(0., 80.);
  
}

//---------------------------------------------------------------------------
void plotMuonEfficiencyVsCent()
{
  /// compute the tracking efficiency versus centrality
  /// !!! to be compiled because of CINT !!!
  
  gStyle->SetFillColor(0);
  
  // output graphs
  TGraphAsymmErrors *trackingEffVsCent = new TGraphAsymmErrors;
  
  // get input hists
  TFile *file = new TFile("AnalysisResults.root", "read");
  if (!file || !file->IsOpen()) {
    printf("cannot open file ./AnalysisResults.root\n");
    return;
  }
  TList *listTT = static_cast<TList*>(file->FindObjectAny("TotalTracksPerChamber"));
  TList *listTD = static_cast<TList*>(file->FindObjectAny("TracksDetectedPerChamber"));
  TH2F *TT = static_cast<TH2F*>(listTT->At(10));
  TH2F *TD = static_cast<TH2F*>(listTD->At(10));
  
  // loop over centrality bins
  TArrayD chambersEff(11);
  TArrayD chambersEffErr[2];
  chambersEffErr[0].Set(11);
  chambersEffErr[1].Set(11);
  for (Int_t centBin = 1; centBin <= TT->GetYaxis()->GetNbins(); centBin++)
  {
    cout << endl << "Centrality " << TT->GetYaxis()->GetBinLowEdge(centBin) << "-" << TT->GetYaxis()->GetBinUpEdge(centBin) << "%:" << endl;
    
    // project efficiency hists
    TH1D *TTdraw1 = TT->ProjectionX("TTdraw1", centBin, centBin);
    TH1D *TDdraw1 = TD->ProjectionX("TDdraw1", centBin, centBin);
    TGraphAsymmErrors *efficiency = 0x0;
    if (TTdraw1->GetEntries() > 0) efficiency = new TGraphAsymmErrors(TDdraw1, TTdraw1);
    delete TTdraw1;
    delete TDdraw1;
    if (!efficiency) {
      cout << "empty bin" << endl;
      trackingEffVsCent->SetPoint(centBin-1,TT->GetYaxis()->GetBinCenter(centBin),0.);
      trackingEffVsCent->SetPointError(centBin-1,0.,0.,0.,0.);
      continue;
    }
    
    // get the individual chamber efficiency
    for (Int_t i = 0; i < 10; i++) {
      chambersEff[i+1] = efficiency->GetY()[i];
      chambersEffErr[0][i+1] = efficiency->GetErrorYlow(i);
      chambersEffErr[1][i+1] = efficiency->GetErrorYhigh(i);
    }
    delete efficiency;
    
    // compute the overall tracking efficiency
    computeTrackingEfficiency(chambersEff, chambersEffErr, kTRUE);
    
    // fill graph
    trackingEffVsCent->SetPoint(centBin-1,TT->GetYaxis()->GetBinCenter(centBin),chambersEff[0]);
    trackingEffVsCent->SetPointError(centBin-1,0.,0.,chambersEffErr[0][0],chambersEffErr[1][0]);
    
  }
  
  // close input file
  file->Close();
  
  // display
  new TCanvas("cTrackingEffVsCent", "Measured tracking efficiency versus centrality",1000,400);
  trackingEffVsCent->SetName("trackingEffVsCent");
  trackingEffVsCent->SetTitle("Measured tracking efficiency versus centrality");
  trackingEffVsCent->SetLineStyle(1);
  trackingEffVsCent->SetLineColor(1); 
  trackingEffVsCent->SetMarkerStyle(20);
  trackingEffVsCent->SetMarkerSize(0.7);
  trackingEffVsCent->SetMarkerColor(2);
  trackingEffVsCent->GetXaxis()->SetTitle("Centrality %");
  trackingEffVsCent->GetXaxis()->SetLabelFont(22);
  trackingEffVsCent->GetXaxis()->SetTitleFont(22);
  trackingEffVsCent->GetYaxis()->SetTitle("Efficiency");
  trackingEffVsCent->GetYaxis()->SetLabelFont(22);
  trackingEffVsCent->GetYaxis()->SetLabelFont(22);
  trackingEffVsCent->SetMinimum(0.6);
  trackingEffVsCent->SetMaximum(0.95);
  trackingEffVsCent->Draw("ap");
  
  // save output
  file = new TFile("efficiency.root","update");
  trackingEffVsCent->Write(0x0, TObject::kOverwrite);
  file->Close();
  
}

//---------------------------------------------------------------------------
void plotMuonEfficiencyVsRun(TString runList)
{
  /// compute the tracking efficiency versus run
  /// !!! to be compiled because of CINT !!!
  
  gStyle->SetFillColor(0);
  
  // output graphs
  TGraphAsymmErrors *trackingEffVsRun = new TGraphAsymmErrors;
  
  TArrayD chambersEff(11);
  TArrayD chambersEffErr[2];
  chambersEffErr[0].Set(11);
  chambersEffErr[1].Set(11);
  
  ifstream inFile(runList.Data());
  if (!inFile.is_open()) {
    printf("cannot open file %s\n",runList.Data());
    return;
  }
  
  // loop over runs
  Int_t irun = -1;
  TString currRun;
  TList runs;
  runs.SetOwner();
  while (!inFile.eof()) {
    
    // get current run number
    currRun.ReadLine(inFile,kTRUE);
    if(currRun.IsNull()) continue;
    runs.AddLast(new TObjString(currRun));
    irun++;
    
    // get input hists
    TFile *file = new TFile(Form("runs/%d/AnalysisResults.root",currRun.Atoi()), "read");
    if (!file || !file->IsOpen()) {
      printf("cannot open file runs/%d/AnalysisResults.root\n",currRun.Atoi());
      trackingEffVsRun->SetPoint(irun,irun,-1.);
      trackingEffVsRun->SetPointError(irun,0.,0.,0.,0.);
      delete file;
      continue;
    }
    TList *listTT = static_cast<TList*>(file->FindObjectAny("TotalTracksPerChamber"));
    TList *listTD = static_cast<TList*>(file->FindObjectAny("TracksDetectedPerChamber"));
    TH2F *TT = static_cast<TH2F*>(listTT->At(10));
    TH2F *TD = static_cast<TH2F*>(listTD->At(10));
    
    // project efficiency hists integrated over the centrality bin 0-80%
    Int_t lowBin = TT->GetYaxis()->FindBin(0.1);
    Int_t upBin = TT->GetYaxis()->FindBin(79.9);
    TH1D *TTdraw1 = TT->ProjectionX("TTdraw1", lowBin, upBin);
    TH1D *TDdraw1 = TD->ProjectionX("TDdraw1", lowBin, upBin);
    TGraphAsymmErrors *efficiency = 0x0;
    if (TTdraw1->GetEntries() > 0) efficiency = new TGraphAsymmErrors(TDdraw1, TTdraw1);
    delete TTdraw1;
    delete TDdraw1;
    if (!efficiency) {
      cout << "run" << currRun.Atoi() << ": empty bin" << endl;
      trackingEffVsRun->SetPoint(irun,irun,-1.);
      trackingEffVsRun->SetPointError(irun,0.,0.,0.,0.);
      file->Close();
      continue;
    }
    
    // get the individual chamber efficiency
    for (Int_t i = 0; i < 10; i++) {
      chambersEff[i+1] = efficiency->GetY()[i];
      chambersEffErr[0][i+1] = efficiency->GetErrorYlow(i);
      chambersEffErr[1][i+1] = efficiency->GetErrorYhigh(i);
    }
    delete efficiency;
    
    // compute the overall tracking efficiency
    computeTrackingEfficiency(chambersEff, chambersEffErr);
    
    // fill graph
    trackingEffVsRun->SetPoint(irun,irun,chambersEff[0]);
    trackingEffVsRun->SetPointError(irun,0.,0.,chambersEffErr[0][0],chambersEffErr[1][0]);
    
    file->Close();
  }
  inFile.close();
  
  // set bin labels
  trackingEffVsRun->GetXaxis()->Set(irun+1, -0.5, irun+0.5);
  TIter nextRun(&runs);
  TObjString *srun = 0x0;
  irun = 1;
  while ((srun = static_cast<TObjString*>(nextRun())))
    trackingEffVsRun->GetXaxis()->SetBinLabel(irun++, srun->GetName());
  
  // display
  new TCanvas("cTrackingEffVsRun", "Measured tracking efficiency versus run",1000,400);
  trackingEffVsRun->SetName("trackingEffVsRun");
  trackingEffVsRun->SetTitle("Measured tracking efficiency versus run");
  trackingEffVsRun->SetLineStyle(1);
  trackingEffVsRun->SetLineColor(1); 
  trackingEffVsRun->SetMarkerStyle(20);
  trackingEffVsRun->SetMarkerSize(0.7);
  trackingEffVsRun->SetMarkerColor(2);
  trackingEffVsRun->GetXaxis()->SetTitle("Run #");
  trackingEffVsRun->GetXaxis()->SetLabelFont(22);
  trackingEffVsRun->GetXaxis()->SetTitleFont(22);
  trackingEffVsRun->GetYaxis()->SetTitle("Efficiency");
  trackingEffVsRun->GetYaxis()->SetLabelFont(22);
  trackingEffVsRun->GetYaxis()->SetLabelFont(22);
  trackingEffVsRun->Draw("ap");
  
  // save output
  TFile* file = new TFile("efficiency.root","update");
  trackingEffVsRun->Write(0x0, TObject::kOverwrite);
  file->Close();
  
}

//---------------------------------------------------------------------------
void integratedEfficiency(Double_t minCent, Double_t maxCent)
{
  /// compute the integrated tracking efficiency in minCent-maxCent% centrality range
  /// !!! to be compiled because of CINT !!!
  
  // adjust centrality range
  maxCent = TMath::Max(maxCent-1.e-12, minCent);
  
  // get input hists
  TFile *file = new TFile("AnalysisResults.root", "read");
  if (!file || !file->IsOpen()) {
    printf("cannot open file ./AnalysisResults.root\n");
    return;
  }
  TList *listTT = static_cast<TList*>(file->FindObjectAny("TotalTracksPerChamber"));
  TList *listTD = static_cast<TList*>(file->FindObjectAny("TracksDetectedPerChamber"));
  TH2F *TT = static_cast<TH2F*>(listTT->At(10));
  TH2F *TD = static_cast<TH2F*>(listTD->At(10));
  
  // get the centrality range for integration
  Int_t lowBin = TT->GetYaxis()->FindBin(minCent);
  Int_t upBin = TT->GetYaxis()->FindBin(maxCent);
  cout << endl << "Integrated efficiency in "<< TT->GetYaxis()->GetBinLowEdge(lowBin) << "-" << TT->GetYaxis()->GetBinUpEdge(upBin) << "%:" << endl;
  
  // project efficiency hists integrated over this centrality range
  TH1D *TTdraw1 = TT->ProjectionX("TTdraw1", lowBin, upBin);
  TH1D *TDdraw1 = TD->ProjectionX("TDdraw1", lowBin, upBin);
  TGraphAsymmErrors *efficiency = 0x0;
  if (TTdraw1->GetEntries() > 0) efficiency = new TGraphAsymmErrors(TDdraw1, TTdraw1);
  delete TTdraw1;
  delete TDdraw1;
  if (!efficiency) {
    cout << "integrated efficiency: empty bin" << endl;
    return;
  }
  
  // get the individual chamber efficiency
  TArrayD chambersEff(11);
  TArrayD chambersEffErr[2];
  chambersEffErr[0].Set(11);
  chambersEffErr[1].Set(11);
  for (Int_t i = 0; i < 10; i++) {
    chambersEff[i+1] = efficiency->GetY()[i];
    chambersEffErr[0][i+1] = efficiency->GetErrorYlow(i);
    chambersEffErr[1][i+1] = efficiency->GetErrorYhigh(i);
  }
  delete efficiency;
  
  // compute the overall tracking efficiency
  computeTrackingEfficiency(chambersEff, chambersEffErr, kTRUE);
  
  // close input file
  file->Close();
  
}

//---------------------------------------------------------------------------
void computeTrackingEfficiency(TArrayD &chambersEff, TArrayD *chambersEffErr, Bool_t print)
{
  /// compute the tracking efficiency from the individual chamber efficiencies
  
  Double_t st1Eff = 1.0 - (1.0-chambersEff[1])*(1.0-chambersEff[2]);
  Double_t st2Eff = 1.0 - (1.0-chambersEff[3])*(1.0-chambersEff[4]);
  Double_t st3Eff = 1.0 - (1.0-chambersEff[5])*(1.0-chambersEff[6]);
  Double_t st45Eff = chambersEff[7] * chambersEff[8] * chambersEff[9] * chambersEff[10] 
  +(1.0 - chambersEff[7]) * chambersEff[8] * chambersEff[9] * chambersEff[10]
  + chambersEff[7] * (1.0 - chambersEff[8]) * chambersEff[9] * chambersEff[10]
  + chambersEff[7] * chambersEff[8] * (1.0 - chambersEff[9]) * chambersEff[10]
  + chambersEff[7] * chambersEff[8] * chambersEff[9] * (1.0 - chambersEff[10]);
  
  chambersEff[0] = st1Eff * st2Eff * st3Eff * st45Eff;
  
  Double_t st1EffErr[2], st2EffErr[2], st3EffErr[2], st45EffErr[2];
  for (Int_t i = 0; i < 2; i++) {
    st1EffErr[i] = TMath::Sqrt((1.-chambersEff[2])*(1.-chambersEff[2])*chambersEffErr[i][1]*chambersEffErr[i][1] + (1.-chambersEff[1])*(1.-chambersEff[1])*chambersEffErr[i][2]*chambersEffErr[i][2]);
    st2EffErr[i] = TMath::Sqrt((1.-chambersEff[3])*(1.-chambersEff[3])*chambersEffErr[i][4]*chambersEffErr[i][4] + (1.-chambersEff[4])*(1.-chambersEff[4])*chambersEffErr[i][3]*chambersEffErr[i][3]);
    st3EffErr[i] = TMath::Sqrt((1.-chambersEff[5])*(1.-chambersEff[5])*chambersEffErr[i][6]*chambersEffErr[i][6] + (1.-chambersEff[6])*(1.-chambersEff[6])*chambersEffErr[i][5]*chambersEffErr[i][5]);
    Double_t x7 = (1.-chambersEff[8])*chambersEff[9]*chambersEff[10] + chambersEff[8]*(1.-chambersEff[9])*chambersEff[10] + chambersEff[8]*chambersEff[9]*(1.-chambersEff[10]);
    Double_t x8 = (1.-chambersEff[7])*chambersEff[9]*chambersEff[10] + chambersEff[7]*(1.-chambersEff[9])*chambersEff[10] + chambersEff[7]*chambersEff[9]*(1.-chambersEff[10]);
    Double_t x9 = (1.-chambersEff[7])*chambersEff[8]*chambersEff[10] + chambersEff[7]*(1.-chambersEff[8])*chambersEff[10] + chambersEff[7]*chambersEff[8]*(1.-chambersEff[10]);
    Double_t x10 = (1.-chambersEff[7])*chambersEff[8]*chambersEff[9] + chambersEff[7]*(1.-chambersEff[8])*chambersEff[9] + chambersEff[7]*chambersEff[8]*(1.-chambersEff[9]);
    st45EffErr[i] = TMath::Sqrt(x7*x7*chambersEffErr[i][7]*chambersEffErr[i][7] + x8*x8*chambersEffErr[i][8]*chambersEffErr[i][8] + x9*x9*chambersEffErr[i][9]*chambersEffErr[i][9] + x10*x10*chambersEffErr[i][10]*chambersEffErr[i][10]);
    
    chambersEffErr[i][0] = chambersEff[0]*TMath::Sqrt(st1EffErr[i]*st1EffErr[i]/st1Eff/st1Eff + st2EffErr[i]*st2EffErr[i]/st2Eff/st2Eff + st3EffErr[i]*st3EffErr[i]/st3Eff/st3Eff + st45EffErr[i]*st45EffErr[i]/st45Eff/st45Eff);
  }
  
  if (print) {
    for (Int_t i = 1; i <= 10; i++) {
      cout << "Efficiency chamber " << i << " : " << chambersEff[i] << " + " << chambersEffErr[1][i] << " - " << chambersEffErr[0][i] << endl;
    }
    cout << "Station 1 = " << st1Eff << " + " << st1EffErr[1] << " - " << st1EffErr[0] << endl;
    cout << "Station 2 = " << st2Eff << " + " << st2EffErr[1] << " - " << st2EffErr[0] << endl;
    cout << "Station 3 = " << st3Eff << " + " << st3EffErr[1] << " - " << st3EffErr[0] << endl;
    cout << "Station 45 = " << st45Eff << " + " << st45EffErr[1] << " - " << st45EffErr[0] << endl;
    cout << "Total tracking efficiency : " << chambersEff[0] << " + " << chambersEffErr[1][0] << " - " << chambersEffErr[0][0] << endl;
  }
  
}
