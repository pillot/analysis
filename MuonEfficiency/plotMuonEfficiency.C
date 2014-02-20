#include <TH1.h>
#include <TH2.h>
#include <THnSparse.h>
#include <TAxis.h>
#include <TString.h>
#include <TObjString.h>
#include <Riostream.h>
#include <TFile.h>
#include <TList.h>
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include <TMath.h>
#include <TArrayD.h>
#include <TStyle.h>

Bool_t request3outOf4 = kTRUE;

Double_t ptMin = 0.;
//TString suffix = "_old_wocut_wMClabel";
TString suffix = "_old_wcut";
//TString suffix = "_old_wcut_wMClabel";
//TString suffix = "";

void plotMuonEfficiencyVsX(TString var, Double_t minCent = -999., Double_t maxCent = 999., Bool_t print = kFALSE);
void plotMuonEfficiencyVsXY(TString xVar, TString yVar, Double_t minCent = -999., Double_t maxCent = 999.);
void plotMuonEfficiencyVsRun(TString runList = "", Double_t minCent = 0., Double_t maxCent = 80.);
void integratedEfficiency(Double_t minCent = 0., Double_t maxCent = 80.);
void computeTrackingEfficiency(TArrayD &chambersEff, TArrayD *chambersEffErr, Bool_t print = kFALSE);

//---------------------------------------------------------------------------
void plotMuonEfficiency(TString runList = "")
{
  /// compute the tracking efficiency versus centrality and versus run if runList != ""
  /// !!! to be compiled because of CINT !!!
  
  plotMuonEfficiencyVsX("centrality", -999., 999., kTRUE);
  plotMuonEfficiencyVsX("pt");
  plotMuonEfficiencyVsX("y");
  plotMuonEfficiencyVsX("phi");
  //plotMuonEfficiencyVsX("charge");
  plotMuonEfficiencyVsXY("pt","centrality");
  plotMuonEfficiencyVsXY("y","centrality");
  plotMuonEfficiencyVsXY("pt","y");
  plotMuonEfficiencyVsXY("y","phi");
  plotMuonEfficiencyVsRun(runList, -999., 999.);
  integratedEfficiency(-999., 999.);
  
}

//---------------------------------------------------------------------------
void plotMuonEfficiencyVsX(TString var, Double_t minCent, Double_t maxCent, Bool_t print)
{
  /// plot the tracking efficiency vs X integrated over minCent-maxCent% centrality range
  /// !!! to be compiled because of CINT !!!
  
  Int_t xDim = -1;
  if (var == "centrality") xDim = 1;
  else if (var == "pt") xDim = 2;
  else if (var == "y") xDim = 3;
  else if (var == "phi") xDim = 4;
  else if (var == "charge") xDim = 5;
  else {
    printf("incorrect variable. Choices are centrality, pt, y, phi and charge.\n");
    return;
  }
  
  // adjust centrality range
  maxCent = TMath::Max(maxCent-1.e-12, minCent);
  
  // output graph
  TGraphAsymmErrors *effVsX = new TGraphAsymmErrors;
  
  // get input hists
  TFile *file = new TFile("AnalysisResults.root", "read");
  if (!file || !file->IsOpen()) {
    printf("cannot open file ./AnalysisResults.root\n");
    return;
  }
  TList *listTT = static_cast<TList*>(file->FindObjectAny(Form("TotalTracksPerChamber%s",suffix.Data())));
  TList *listTD = static_cast<TList*>(file->FindObjectAny(Form("TracksDetectedPerChamber%s",suffix.Data())));
  THnSparse *TT = static_cast<THnSparse*>(listTT->At(10));
  THnSparse *TD = static_cast<THnSparse*>(listTD->At(10));
  
  // set the centrality range for integration
  Int_t lowBin = TT->GetAxis(1)->FindBin(minCent);
  Int_t upBin = TT->GetAxis(1)->FindBin(maxCent);
  TT->GetAxis(1)->SetRange(lowBin, upBin);
  TD->GetAxis(1)->SetRange(lowBin, upBin);
  
  // set the pt range for integration
  lowBin = TT->GetAxis(2)->FindBin(ptMin);
  upBin = TT->GetAxis(2)->GetNbins()+1;
  TT->GetAxis(2)->SetRange(lowBin, upBin);
  TD->GetAxis(2)->SetRange(lowBin, upBin);
  
  // loop over X bins
  TArrayD chambersEff(11);
  TArrayD chambersEffErr[2];
  chambersEffErr[0].Set(11);
  chambersEffErr[1].Set(11);
  for (Int_t ix = 1; ix <= TT->GetAxis(xDim)->GetNbins(); ix++) {
    
    if (print) cout << endl << var.Data() << " " << TT->GetAxis(xDim)->GetBinLowEdge(ix) << "-" << TT->GetAxis(xDim)->GetBinUpEdge(ix) << ":" << endl;
    
    // project efficiency hists
    TT->GetAxis(xDim)->SetRange(ix, ix);
    TH1D *TTdraw1 = TT->Projection(0,"e");
    TD->GetAxis(xDim)->SetRange(ix, ix);
    TH1D *TDdraw1 = TD->Projection(0,"e");
    /*
    new TCanvas;
    TTdraw1->DrawClone();
    TDdraw1->SetLineColor(2);
    TDdraw1->DrawClone("same");
    */
    TGraphAsymmErrors *efficiency = 0x0;
    if (TTdraw1->GetEntries() > 0) efficiency = new TGraphAsymmErrors(TDdraw1, TTdraw1, "cpe0");
    delete TTdraw1;
    delete TDdraw1;
    if (!efficiency) {
      effVsX->SetPoint(ix-1,TT->GetAxis(xDim)->GetBinCenter(ix),0.);
      effVsX->SetPointError(ix-1,0.,0.,0.,0.);
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
    computeTrackingEfficiency(chambersEff, chambersEffErr, print);
    
    // fill graph
    effVsX->SetPoint(ix-1,TT->GetAxis(xDim)->GetBinCenter(ix),chambersEff[0]);
    effVsX->SetPointError(ix-1,0.,0.,chambersEffErr[0][0],chambersEffErr[1][0]);
    
  }
  
  // close input file
  file->Close();
  
  // display
  new TCanvas(Form("cTrackingEffVs%s",var.Data()), Form("Measured tracking efficiency versus %s",var.Data()),1000,400);
  effVsX->SetName(Form("trackingEffVs%s",var.Data()));
  effVsX->SetTitle(Form("Measured tracking efficiency versus %s",var.Data()));
  effVsX->SetLineStyle(1);
  effVsX->SetLineColor(1); 
  effVsX->SetMarkerStyle(20);
  effVsX->SetMarkerSize(0.7);
  effVsX->SetMarkerColor(2);
  effVsX->GetXaxis()->SetTitle(var.Data());
  effVsX->GetXaxis()->SetLabelFont(22);
  effVsX->GetXaxis()->SetTitleFont(22);
  effVsX->GetYaxis()->SetTitle("Efficiency");
  effVsX->GetYaxis()->SetLabelFont(22);
  effVsX->GetYaxis()->SetTitleFont(22);
  effVsX->SetMinimum(0.6);
  effVsX->SetMaximum(0.95);
  effVsX->Draw("ap");
  
  // save output
  file = new TFile("efficiency.root","update");
  effVsX->Write(0x0, TObject::kOverwrite);
  file->Close();
  
}

//---------------------------------------------------------------------------
void plotMuonEfficiencyVsXY(TString xVar, TString yVar, Double_t minCent, Double_t maxCent)
{
  /// plot the tracking efficiency vs X,Y integrated over minCent-maxCent% centrality range
  /// !!! to be compiled because of CINT !!!
  
  Int_t xDim = -1;
  if (xVar == "centrality") xDim = 1;
  else if (xVar == "pt") xDim = 2;
  else if (xVar == "y") xDim = 3;
  else if (xVar == "phi") xDim = 4;
  else if (xVar == "charge") xDim = 5;
  else {
    printf("incorrect variable. Choices are centrality, pt, y, phi and charge.\n");
    return;
  }
  Int_t yDim = -1;
  if (yVar == "centrality") yDim = 1;
  else if (yVar == "pt") yDim = 2;
  else if (yVar == "y") yDim = 3;
  else if (yVar == "phi") yDim = 4;
  else if (yVar == "charge") yDim = 5;
  else {
    printf("incorrect variable. Choices are centrality, pt, y, phi and charge.\n");
    return;
  }
  
  // adjust centrality range
  maxCent = TMath::Max(maxCent-1.e-12, minCent);
  
  // get input hists
  TFile *file = new TFile("AnalysisResults.root", "read");
  if (!file || !file->IsOpen()) {
    printf("cannot open file ./AnalysisResults.root\n");
    return;
  }
  TList *listTT = static_cast<TList*>(file->FindObjectAny(Form("TotalTracksPerChamber%s",suffix.Data())));
  TList *listTD = static_cast<TList*>(file->FindObjectAny(Form("TracksDetectedPerChamber%s",suffix.Data())));
  THnSparse *TT = static_cast<THnSparse*>(listTT->At(10));
  THnSparse *TD = static_cast<THnSparse*>(listTD->At(10));
  
  // set the centrality range for integration
  Int_t lowBin = TT->GetAxis(1)->FindBin(minCent);
  Int_t upBin = TT->GetAxis(1)->FindBin(maxCent);
  TT->GetAxis(1)->SetRange(lowBin, upBin);
  TD->GetAxis(1)->SetRange(lowBin, upBin);
  
  // set the pt range for integration
  lowBin = TT->GetAxis(2)->FindBin(ptMin);
  upBin = TT->GetAxis(2)->GetNbins()+1;
  TT->GetAxis(2)->SetRange(lowBin, upBin);
  TD->GetAxis(2)->SetRange(lowBin, upBin);
  
  // output map
  Int_t nxBins = TT->GetAxis(xDim)->GetNbins();
  Int_t nyBins = TT->GetAxis(yDim)->GetNbins();
  TH2F *effVsXY = new TH2F(Form("trackingEffVs%s-%s",xVar.Data(),yVar.Data()),
			   Form("Measured tracking efficiency versus %s / %s",xVar.Data(),yVar.Data()),
			   nxBins, TT->GetAxis(xDim)->GetBinLowEdge(1), TT->GetAxis(xDim)->GetBinUpEdge(nxBins),
			   nyBins, TT->GetAxis(yDim)->GetBinLowEdge(1), TT->GetAxis(yDim)->GetBinUpEdge(nyBins));
  effVsXY->SetDirectory(0);
  
  // loop over X/Y bins
  TArrayD chambersEff(11);
  TArrayD chambersEffErr[2];
  chambersEffErr[0].Set(11);
  chambersEffErr[1].Set(11);
  for (Int_t ix = 1; ix <= nxBins; ix++) {
    
    // set x range
    TT->GetAxis(xDim)->SetRange(ix, ix);
    TD->GetAxis(xDim)->SetRange(ix, ix);
    
    for (Int_t iy = 1; iy <= nyBins; iy++) {
      
      // project efficiency hists
      TT->GetAxis(yDim)->SetRange(iy, iy);
      TH1D *TTdraw1 = TT->Projection(0,"e");
      TD->GetAxis(yDim)->SetRange(iy, iy);
      TH1D *TDdraw1 = TD->Projection(0,"e");
      
      TGraphAsymmErrors *efficiency = 0x0;
      if (TTdraw1->GetEntries() > 0) efficiency = new TGraphAsymmErrors(TDdraw1, TTdraw1, "cpe0");
      delete TTdraw1;
      delete TDdraw1;
      if (!efficiency) {
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
      if (chambersEff[0] > 1.) chambersEff[0] = 1.;
      if (chambersEffErr[0][0] > 1.) chambersEffErr[0][0] = 1.;
      
      // fill histo
      if (chambersEff[0] > 0.) {
	effVsXY->Fill(TT->GetAxis(xDim)->GetBinCenter(ix),TT->GetAxis(yDim)->GetBinCenter(iy),chambersEff[0]);
	effVsXY->SetBinError(ix,iy,chambersEffErr[0][0]);
      }
      
    }
    
  }
  
  // close input file
  file->Close();
  
  // display
  new TCanvas(Form("cTrackingEffVs%s-%s",xVar.Data(),yVar.Data()), Form("Measured tracking efficiency versus %s/%s",xVar.Data(),yVar.Data()),700,600);
  effVsXY->GetXaxis()->SetTitle(xVar.Data());
  effVsXY->GetXaxis()->SetLabelFont(22);
  effVsXY->GetXaxis()->SetTitleFont(22);
  effVsXY->GetYaxis()->SetTitle(yVar.Data());
  effVsXY->GetYaxis()->SetLabelFont(22);
  effVsXY->GetYaxis()->SetTitleFont(22);
  effVsXY->GetZaxis()->SetTitle("Efficiency");
  effVsXY->GetZaxis()->SetLabelFont(22);
  effVsXY->GetZaxis()->SetTitleFont(22);
  //effVsXY->SetMinimum(0.6);
  //effVsXY->SetMaximum(0.95);
  effVsXY->Draw("surf1");
  
  // save output
  file = new TFile("efficiency.root","update");
  effVsXY->Write(0x0, TObject::kOverwrite);
  file->Close();
  
}

//---------------------------------------------------------------------------
void plotMuonEfficiencyVsRun(TString runList, Double_t minCent, Double_t maxCent)
{
  /// compute the tracking efficiency versus run
  /// !!! to be compiled because of CINT !!!
  
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
    TList *listTT = static_cast<TList*>(file->FindObjectAny(Form("TotalTracksPerChamber%s",suffix.Data())));
    TList *listTD = static_cast<TList*>(file->FindObjectAny(Form("TracksDetectedPerChamber%s",suffix.Data())));
    THnSparse *TT = static_cast<THnSparse*>(listTT->At(10));
    THnSparse *TD = static_cast<THnSparse*>(listTD->At(10));
    
    // set the centrality range for integration
    Int_t lowBin = TT->GetAxis(1)->FindBin(minCent);
    Int_t upBin = TT->GetAxis(1)->FindBin(maxCent);
    TT->GetAxis(1)->SetRange(lowBin, upBin);
    TD->GetAxis(1)->SetRange(lowBin, upBin);
    
    // set the pt range for integration
    lowBin = TT->GetAxis(2)->FindBin(ptMin);
    upBin = TT->GetAxis(2)->GetNbins()+1;
    TT->GetAxis(2)->SetRange(lowBin, upBin);
    TD->GetAxis(2)->SetRange(lowBin, upBin);
    
    // project efficiency hists
    TH1D *TTdraw1 = TT->Projection(0,"e");
    TH1D *TDdraw1 = TD->Projection(0,"e");
    
    TGraphAsymmErrors *efficiency = 0x0;
    if (TTdraw1->GetEntries() > 0) efficiency = new TGraphAsymmErrors(TDdraw1, TTdraw1, "cpe0");
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
//    cout << "run" << currRun.Atoi() << ":" << endl;
//    computeTrackingEfficiency(chambersEff, chambersEffErr, kTRUE);
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
  trackingEffVsRun->GetYaxis()->SetTitleFont(22);
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
  TList *listTT = static_cast<TList*>(file->FindObjectAny(Form("TotalTracksPerChamber%s",suffix.Data())));
  TList *listTD = static_cast<TList*>(file->FindObjectAny(Form("TracksDetectedPerChamber%s",suffix.Data())));
  THnSparse *TT = static_cast<THnSparse*>(listTT->At(10));
  THnSparse *TD = static_cast<THnSparse*>(listTD->At(10));
  
  // get the centrality range for integration
  Int_t lowBin = TT->GetAxis(1)->FindBin(minCent);
  Int_t upBin = TT->GetAxis(1)->FindBin(maxCent);
  cout << endl << "Integrated efficiency in "<< TT->GetAxis(1)->GetBinLowEdge(lowBin) << "-" << TT->GetAxis(1)->GetBinUpEdge(upBin) << "%:" << endl;
  
  // set the centrality range for integration
  TT->GetAxis(1)->SetRange(lowBin, upBin);
  TD->GetAxis(1)->SetRange(lowBin, upBin);
  
  // set the pt range for integration
  lowBin = TT->GetAxis(2)->FindBin(ptMin);
  upBin = TT->GetAxis(2)->GetNbins()+1;
  TT->GetAxis(2)->SetRange(lowBin, upBin);
  TD->GetAxis(2)->SetRange(lowBin, upBin);
  
  // project efficiency hists
  TH1D *TTdraw1 = TT->Projection(0,"e");
  TH1D *TDdraw1 = TD->Projection(0,"e");
  TGraphAsymmErrors *efficiency = 0x0;
  if (TTdraw1->GetEntries() > 0) efficiency = new TGraphAsymmErrors(TDdraw1, TTdraw1, "cpe0");
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
  Double_t st4Eff = 1.0 - (1.0-chambersEff[7])*(1.0-chambersEff[8]);
  Double_t st5Eff = 1.0 - (1.0-chambersEff[9])*(1.0-chambersEff[10]);
  Double_t st45Eff = chambersEff[7] * chambersEff[8] * chambersEff[9] * chambersEff[10] 
  +(1.0 - chambersEff[7]) * chambersEff[8] * chambersEff[9] * chambersEff[10]
  + chambersEff[7] * (1.0 - chambersEff[8]) * chambersEff[9] * chambersEff[10]
  + chambersEff[7] * chambersEff[8] * (1.0 - chambersEff[9]) * chambersEff[10]
  + chambersEff[7] * chambersEff[8] * chambersEff[9] * (1.0 - chambersEff[10]);
  
  if (request3outOf4) chambersEff[0] = st1Eff * st2Eff * st3Eff * st45Eff;
  else chambersEff[0] = st1Eff * st2Eff * st3Eff * st4Eff * st5Eff;
  
  Double_t st1EffErr[2], st2EffErr[2], st3EffErr[2], st4EffErr[2], st5EffErr[2], st45EffErr[2];
  for (Int_t i = 0; i < 2; i++) {
    st1EffErr[i] = TMath::Sqrt((1.-chambersEff[1])*(1.-chambersEff[1])*chambersEffErr[i][2]*chambersEffErr[i][2] + (1.-chambersEff[2])*(1.-chambersEff[2])*chambersEffErr[i][1]*chambersEffErr[i][1]);
    st2EffErr[i] = TMath::Sqrt((1.-chambersEff[3])*(1.-chambersEff[3])*chambersEffErr[i][4]*chambersEffErr[i][4] + (1.-chambersEff[4])*(1.-chambersEff[4])*chambersEffErr[i][3]*chambersEffErr[i][3]);
    st3EffErr[i] = TMath::Sqrt((1.-chambersEff[5])*(1.-chambersEff[5])*chambersEffErr[i][6]*chambersEffErr[i][6] + (1.-chambersEff[6])*(1.-chambersEff[6])*chambersEffErr[i][5]*chambersEffErr[i][5]);
    st4EffErr[i] = TMath::Sqrt((1.-chambersEff[7])*(1.-chambersEff[7])*chambersEffErr[i][8]*chambersEffErr[i][8] + (1.-chambersEff[8])*(1.-chambersEff[8])*chambersEffErr[i][7]*chambersEffErr[i][7]);
    st5EffErr[i] = TMath::Sqrt((1.-chambersEff[9])*(1.-chambersEff[9])*chambersEffErr[i][10]*chambersEffErr[i][10] + (1.-chambersEff[10])*(1.-chambersEff[10])*chambersEffErr[i][9]*chambersEffErr[i][9]);
    Double_t x7 = (1.-chambersEff[8])*chambersEff[9]*chambersEff[10] + chambersEff[8]*(1.-chambersEff[9])*chambersEff[10] + chambersEff[8]*chambersEff[9]*(1.-chambersEff[10]);
    Double_t x8 = (1.-chambersEff[7])*chambersEff[9]*chambersEff[10] + chambersEff[7]*(1.-chambersEff[9])*chambersEff[10] + chambersEff[7]*chambersEff[9]*(1.-chambersEff[10]);
    Double_t x9 = (1.-chambersEff[7])*chambersEff[8]*chambersEff[10] + chambersEff[7]*(1.-chambersEff[8])*chambersEff[10] + chambersEff[7]*chambersEff[8]*(1.-chambersEff[10]);
    Double_t x10 = (1.-chambersEff[7])*chambersEff[8]*chambersEff[9] + chambersEff[7]*(1.-chambersEff[8])*chambersEff[9] + chambersEff[7]*chambersEff[8]*(1.-chambersEff[9]);
    st45EffErr[i] = TMath::Sqrt(x7*x7*chambersEffErr[i][7]*chambersEffErr[i][7] + x8*x8*chambersEffErr[i][8]*chambersEffErr[i][8] + x9*x9*chambersEffErr[i][9]*chambersEffErr[i][9] + x10*x10*chambersEffErr[i][10]*chambersEffErr[i][10]);
    
    if (request3outOf4) chambersEffErr[i][0] = chambersEff[0]*TMath::Sqrt(st1EffErr[i]*st1EffErr[i]/st1Eff/st1Eff + st2EffErr[i]*st2EffErr[i]/st2Eff/st2Eff + st3EffErr[i]*st3EffErr[i]/st3Eff/st3Eff + st45EffErr[i]*st45EffErr[i]/st45Eff/st45Eff);
    else chambersEffErr[i][0] = chambersEff[0]*TMath::Sqrt(st1EffErr[i]*st1EffErr[i]/st1Eff/st1Eff + st2EffErr[i]*st2EffErr[i]/st2Eff/st2Eff + st3EffErr[i]*st3EffErr[i]/st3Eff/st3Eff + st4EffErr[i]*st4EffErr[i]/st4Eff/st4Eff + st5EffErr[i]*st5EffErr[i]/st5Eff/st5Eff);
  }
  
  if (print) {
    for (Int_t i = 1; i <= 10; i++) {
      cout << "Efficiency chamber " << i << " : " << chambersEff[i] << " + " << chambersEffErr[1][i] << " - " << chambersEffErr[0][i] << endl;
    }
    cout << "Station 1 = " << st1Eff << " + " << st1EffErr[1] << " - " << st1EffErr[0] << endl;
    cout << "Station 2 = " << st2Eff << " + " << st2EffErr[1] << " - " << st2EffErr[0] << endl;
    cout << "Station 3 = " << st3Eff << " + " << st3EffErr[1] << " - " << st3EffErr[0] << endl;
    cout << "Station 4 = " << st4Eff << " + " << st4EffErr[1] << " - " << st4EffErr[0] << endl;
    cout << "Station 5 = " << st5Eff << " + " << st5EffErr[1] << " - " << st5EffErr[0] << endl;
    cout << "Station 45 = " << st45Eff << " + " << st45EffErr[1] << " - " << st45EffErr[0] << endl;
    cout << "Total tracking efficiency : " << chambersEff[0] << " + " << chambersEffErr[1][0] << " - " << chambersEffErr[0][0] << endl;
  }
  
}
