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
#include <THashList.h>
#include <TParameter.h>
//--------------------------------//
#include "TGrid.h"
#include "AliMpDEIterator.h"
#include "AliMUONCDB.h"
#include "AliCDBManager.h"
#include "AliMpDEManager.h"
#include <TGaxis.h>
#include <TLegend.h>
//--------------------------------//

// TODO
  // 1) efficiency estimator and error calculation at chamber and DE level:
  //  2 options:
  //    - using bayesian method with uniform prior
  //    - using Feldman Cousins frequentist method
  //  if n ≠ 0: use above methods
  //  if n = 0: eff = -1 ± 0
  // 2) efficiency and error propagation at station and spectrometer level:
  //  if eff = -1 for one or several ch/DE:
  //    - assume eff_ch = 1 ± 0 to compute eff_up and err_up with std error propagation at nth order
  //    - assume eff_ch = 0 ± 0 to compute eff_low and err_low with std error propagation at nth order
  //    - eff_spectro = eff_up + err_up - (sqrt((eff_up-eff_low)^2 + err_low^2))  (ou - (eff_up-eff_low + err_low))
  //  otherwise: std efficiency and error propagation at nth order

//Possible improvements:
  //Set properly the cent and pt range
  //>Load AnalysisResults.root and tracks THnsparse graphs just once
  //>Not display histos
  //>Condense plotMuonEfficiencyVsRun
  //>Set output file name and path as imput
  //>Reorganize outputs in "folders" to see them easily in the TBrowser
  //>Set input file name and path (AnalysisResults.root) as input
  //>Change the way of setting paths
  //>Write efficiencyDEperDE so it can be used in plotMuonEfficiencyVsRun
  //>Create functions to set the cuts in pt, and centrality as in plotMuonEfficiencyVsX and plotMuonEfficiencyVsXY, replace them there and add them also to plotMuonEfficiencyVsRun and plotMuonEfficiencyDEperDE
  //> Save the computetrakingEfficiency result in a txt file
  //>Change the print purpose from printing the efficincy to display histos and let print the eff values just when using the integratedEfficiency function

//Questions:
  //> In efficiencyDEperDe what is better, to fill the efficiency graph point by point, or remove the xError by hand? Is bad to cast too much?
  //> In the efficiency task the bin of DE are not integers, should we change it to make easier the filling of efficiency histos?

Double_t ptMin = 0.;
Int_t Charge = 0; // 0 selects + and -, -1 and +1 selects - or + muons respectively
Bool_t moreTrackCandidates = kTRUE;
THashList *runWeights = 0x0;

void plotMuonEfficiencyVsX(TString var, Double_t minCent, Double_t maxCent, Bool_t print, TString fileNameData = "AnalysisResults.root", TString fileNameSave = "efficiency_new.root");
void plotMuonEfficiencyVsX(TString var, Double_t minCent, Double_t maxCent, TString runList, TString fileNameWeights, Bool_t print, TString fileNameData = "AnalysisResults.root", TString fileNameSave = "efficiency_new.root");
void plotMuonEfficiencyVsXY(TString xVar, TString yVar, Double_t minCent, Double_t maxCent, TString fileNameData = "AnalysisResults.root", TString fileNameSave = "efficiency_new.root", Bool_t rap = kFALSE );
void plotMuonEfficiencyVsRun(Double_t minCent, Double_t maxCent, TString runList = "", TString fileNameSave = "efficiency_new.root",TString alienPath =""); //Modifyed by Javier Martin
void plotMuonEfficiencyDEperDE(Double_t minCent, Double_t maxCent, TString fileNameData = "AnalysisResults.root", TString fileNameSave = "efficiency_new.root"); //Added by Javier Martin
void integratedEfficiency(Double_t minCent, Double_t maxCent, TString fileNameData = "AnalysisResults.root");
void integratedEfficiency(TString fileNameWeights, TString fileNameSave = "efficiency_new.root");
void computeTrackingEfficiency(TArrayD &chambersEff, TArrayD *chambersEffErr, Bool_t print = kFALSE);
void LoadRunWeights(TString fileName);

//---------------------------------------------------------------------------
void RealEfficiency(TString runList = "runList.txt",
		    TString fileNameWeights = "",
		    TString fileNameData ="AnalysisResults.root",
		    TString fileNameSave = "efficiency_new.root",
		    TString alienPath ="")//Change paths to the desired one, the alien path is relative to the alice/cern.ch/user one (token needed)
{
  /// compute the tracking efficiency versus centrality and versus run if runList != ""
  /// !!! to be compiled because of CINT !!!
  /*
  .x $ALICE_ROOT/MUON/rootlogon.C
  .x $WORK/Macros/MuonEfficiency/RealEfficiency.C+
  */
  /*
  plotMuonEfficiencyVsX("centrality",-999.,999.,kFALSE,fileNameData,fileNameSave);
  plotMuonEfficiencyVsX("pt",-999.,999.,kFALSE,fileNameData,fileNameSave);
  plotMuonEfficiencyVsX("y",-999.,999.,kFALSE,fileNameData,fileNameSave);
  plotMuonEfficiencyVsX("phi",-999.,999.,kFALSE,fileNameData,fileNameSave);
  plotMuonEfficiencyVsX("charge",-999.,999.,kFALSE,fileNameData,fileNameSave);
  */
  plotMuonEfficiencyVsX("centrality", -999., 999., runList, fileNameWeights, kFALSE, fileNameData, fileNameSave);
  plotMuonEfficiencyVsX("pt", -999., 999., runList, fileNameWeights, kFALSE, fileNameData, fileNameSave);
  plotMuonEfficiencyVsX("y", -999., 999., runList, fileNameWeights, kFALSE, fileNameData, fileNameSave);
  plotMuonEfficiencyVsX("phi", -999., 999., runList, fileNameWeights, kFALSE, fileNameData, fileNameSave);
  plotMuonEfficiencyVsX("charge", -999., 999., runList, fileNameWeights, kFALSE, fileNameData, fileNameSave);
  /*
  plotMuonEfficiencyVsXY("pt","centrality",-999.,999.,fileNameData,fileNameSave);
  plotMuonEfficiencyVsXY("y","centrality",-999.,999.,fileNameData,fileNameSave);
  plotMuonEfficiencyVsXY("pt","y",-999.,999.,fileNameData,fileNameSave);
  plotMuonEfficiencyVsXY("phi","y",-999.,999.,fileNameData,fileNameSave,kTRUE);
  plotMuonEfficiencyDEperDE(-999.,999.,fileNameData,fileNameSave);
  plotMuonEfficiencyVsRun(-999.,999.,runList,fileNameSave,alienPath);
  */
  integratedEfficiency(-999.,999.,fileNameData);
  if (!fileNameWeights.IsNull()) integratedEfficiency(fileNameWeights,fileNameSave);

  
}


//---------------------------------------------------------------------------
void setCentPtCh(THnSparse& SparseData, Double_t minCent = 0., Double_t maxCent = 100., Int_t charge = Charge, Double_t minPt = ptMin, Double_t maxPt = -1) //Added by Javier Martin
{
  //Sets the Centrality and Pt range for integration. If maxPt = -1, it sets maxPt as the maximum Pt value on the THnSparse
  
//  if (minCent < 0 || maxCent < 0 )
//  {
//    printf("Centrality must be positive, but you cuoted minCent = %f and maxCent = %f\n ",minCent,maxCent);
//    return;
//  }
  
  if ( maxCent < minCent)
  {
    printf("Wrong centrality range, maxCent must be higher than minCent\n setCentPtCh usage: minCent, maxCent, Charge, minPt, maxPt\n");
    return;
  }
  
  if (maxPt != -1 && (minPt < 0 || maxPt < 0) )
  {
    printf("Pt must be positive, but you cuoted minPt = %f and maxPt = %f\n ",minPt,maxPt);
    return;
  }
  
  if (maxPt != -1 && maxPt < minPt)
  {
    printf("Wrong Pt range, maxPt must be higher than minPt\n setCentPtCh usage: minCent, maxCent, Charge, minPt, maxPt\n");
    return;
  }
  
  if ( (charge != 0 && charge != 1 && charge != -1) )
  {
    printf("Selected charge must be 0, 1 or -1\n setCentPtCh usage: minCent, maxCent, Charge ,minPt, maxPt\n");
    return;
  }
  
  // set the centrality range for integration
  Int_t lowBin = SparseData.GetAxis(1)->FindBin(minCent);
  Int_t upBin = SparseData.GetAxis(1)->FindBin(maxCent);
  SparseData.GetAxis(1)->SetRange(lowBin, upBin);

  // set the pt range for integration
  lowBin = SparseData.GetAxis(2)->FindBin(minPt);
  if (maxPt == -1 ){ upBin = SparseData.GetAxis(2)->GetNbins()+1; }
  else { upBin = SparseData.GetAxis(2)->FindBin(maxPt);}
  SparseData.GetAxis(2)->SetRange(lowBin, upBin);
  
  // set the charge range
  if (charge != 0)
  {
    lowBin = SparseData.GetAxis(5)->FindBin(charge);
    SparseData.GetAxis(5)->SetRange(lowBin, lowBin);
  }
}

//---------------------------------------------------------------------------
TGraphAsymmErrors* CreateGraph(const char* name, const char* title, int value=-1) //Added by Javier Martin
{
  TGraphAsymmErrors* g = new TGraphAsymmErrors();
  
  g->SetName(name);
  g->SetTitle(title);
  
  if ( value >= 0 )
  {
    g->SetTitle(Form(title,value));
    g->SetName(Form(name,value));
  }
  return g;
}

//---------------------------------------------------------------------------
void Zero(TObjArray& array, Int_t irun, Int_t firstIndex=0, Int_t lastIndex=-1) //Added by Javier Martin
{
  TGraphAsymmErrors* g;
  
  for ( Int_t i = firstIndex; i <= ( lastIndex >= 0 ? lastIndex : array.GetLast()) ; ++ i )
  {
    g = static_cast<TGraphAsymmErrors*>(array.At(i));
    g->SetPoint(irun,irun,-1);
    g->SetPointError(irun,0.,0.,0.,0.);
  }
}

//---------------------------------------------------------------------------
void BeautifyGraphs(TObjArray& array, const char* xAxisName, const char* yAxisName) //Added by Javier Martin
{
  TGraphAsymmErrors* g;
  
  for ( Int_t i = 0 ; i <= array.GetLast() ; ++ i )
  {
    g = static_cast<TGraphAsymmErrors*>(array.At(i));
    
    g->SetLineStyle(1);
    g->SetLineColor(1);
    g->SetMarkerStyle(20);
    g->SetMarkerSize(0.7);
    g->SetMarkerColor(2);
    g->GetXaxis()->SetTitle(xAxisName);
    g->GetXaxis()->CenterTitle(kTRUE);
    g->GetXaxis()->SetLabelFont(22);
    g->GetXaxis()->SetTitleFont(22);
    g->GetYaxis()->SetTitle(yAxisName);
    g->GetYaxis()->CenterTitle(kTRUE);
    g->GetYaxis()->SetLabelFont(22);
    g->GetYaxis()->SetTitleFont(22);
  }
  
}

//---------------------------------------------------------------------------
void SetRunLabel(TObjArray& array, Int_t irun, const TList& runs) //Added by Javier Martin
{
  TGraphAsymmErrors* g;
  
  for ( Int_t i = 0 ; i <= array.GetLast() ; ++ i )
  {
    g = static_cast<TGraphAsymmErrors*>(array.At(i));
    
    g->GetXaxis()->Set(irun+1, -0.5, irun+0.5);
    
    TIter nextRun(&runs);
    TObjString *srun = 0x0;
    Int_t iirun = 1;
    while ((srun = static_cast<TObjString*>(nextRun())))
    {
      g->GetXaxis()->SetBinLabel(iirun, srun->GetName());
      iirun++;
    }
    
  }
}


//---------------------------------------------------------------------------
void plotMuonEfficiencyVsX(TString var, Double_t minCent, Double_t maxCent, Bool_t print, TString fileNameData, TString fileNameSave)
{
  /// plot the tracking efficiency vs X integrated over minCent-maxCent% centrality range
  /// !!! to be compiled because of CINT !!!
  
  printf("plotting efficiency versus %s...\n", var.Data());
  
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
  TFile *file = new TFile(fileNameData.Data(), "read");
  if (!file || !file->IsOpen()) {
    printf("cannot open file %s \n",fileNameData.Data());
    return;
  }
  TList *listTT = static_cast<TList*>(file->FindObjectAny("TotalTracksPerChamber"));
  TList *listTD = static_cast<TList*>(file->FindObjectAny("TracksDetectedPerChamber"));
  THnSparse *TT = static_cast<THnSparse*>(listTT->At(10));
  THnSparse *TD = static_cast<THnSparse*>(listTD->At(10));
  
  // set the centrality and Pt range for integration
  setCentPtCh(*TT,minCent,maxCent);  //Added by Javier Martin
  setCentPtCh(*TD,minCent,maxCent);
  
  // loop over X bins
  TArrayD chambersEff(11);
  TArrayD chambersEffErr[2];
  chambersEffErr[0].Set(11);
  chambersEffErr[1].Set(11);
  Double_t minEff = 1.;
  for (Int_t ix = 1; ix <= TT->GetAxis(xDim)->GetNbins(); ix++) {
    
    if (print) cout << endl << var.Data() << " " << TT->GetAxis(xDim)->GetBinLowEdge(ix) << "-" << TT->GetAxis(xDim)->GetBinUpEdge(ix) << ":" << endl;
    
    // project efficiency hists
    TT->GetAxis(xDim)->SetRange(ix, ix); // In this way we take the bin ix for the XDim axis
    TH1D *TTdraw1 = TT->Projection(0,"e"); // Now we project over the chamber axis, integrating the others, so we get a TH1 with Traces VS chamber for the ix bin
    TD->GetAxis(xDim)->SetRange(ix, ix);
    TH1D *TDdraw1 = TD->Projection(0,"e");
    TGraphAsymmErrors *efficiency = 0x0;
    if (TTdraw1->GetEntries() > 0) efficiency = new TGraphAsymmErrors(TDdraw1, TTdraw1,"cpe0");
    if (!efficiency) {
      effVsX->SetPoint(ix-1,TT->GetAxis(xDim)->GetBinCenter(ix),1.);
      effVsX->SetPointError(ix-1,0.,0.,1.,0.);
      delete TTdraw1;
      delete TDdraw1;
      continue;
    }
    
    // get the individual chamber efficiency
    for (Int_t i = 0; i < 10; i++) {
      if (TTdraw1->GetBinContent(i+1) > 0) {
	chambersEff[i+1] = efficiency->GetY()[i];
	chambersEffErr[0][i+1] = efficiency->GetErrorYlow(i);
	chambersEffErr[1][i+1] = efficiency->GetErrorYhigh(i);
      } else {
	chambersEff[i+1] = -1.;
	chambersEffErr[0][i+1] = 0.;
	chambersEffErr[1][i+1] = 0.;
      }
    }
    delete TTdraw1;
    delete TDdraw1;
    delete efficiency;
    
    // compute the overall tracking efficiency
    computeTrackingEfficiency(chambersEff, chambersEffErr, print);
    
    // fill graph
    effVsX->SetPoint(ix-1,TT->GetAxis(xDim)->GetBinCenter(ix),chambersEff[0]);
    effVsX->SetPointError(ix-1,0.,0.,chambersEffErr[0][0],chambersEffErr[1][0]);
    if (chambersEff[0] < minEff) minEff = chambersEff[0];
    
  }
  
  // close input file
  file->Close();
  
  // display
 // new TCanvas(Form("cTrackingEffVs%s",var.Data()), Form("Measured tracking efficiency versus %s",var.Data()),1000,400);
  effVsX->SetName(Form("trackingEffVs%s",var.Data()));
  effVsX->SetTitle(Form("Measured tracking efficiency versus %s",var.Data()));
  effVsX->SetLineStyle(1);
  effVsX->SetLineColor(1); 
  effVsX->SetMarkerStyle(20);
  effVsX->SetMarkerSize(0.7);
  effVsX->SetMarkerColor(2);
  effVsX->GetXaxis()->SetTitle(var.Data());
  effVsX->GetXaxis()->CenterTitle(kTRUE);
  effVsX->GetXaxis()->SetLabelFont(22);
  effVsX->GetXaxis()->SetTitleFont(22);
  effVsX->GetYaxis()->SetTitle("Efficiency");
  effVsX->GetYaxis()->SetLabelFont(22);
  effVsX->GetYaxis()->SetTitleFont(22);
  effVsX->SetMinimum(minEff);
//  effVsX->SetMaximum(0.95);
 // effVsX->Draw("ap");
  
  // save output
  file = new TFile(fileNameSave.Data(),"update");
  effVsX->Write(0x0, TObject::kOverwrite);
  //effVsX->Write();
  file->Close();
  
}

//---------------------------------------------------------------------------
void plotMuonEfficiencyVsX(TString var, Double_t minCent, Double_t maxCent, TString runList, TString fileNameWeights, Bool_t print, TString fileNameData, TString fileNameSave)
{
  /// plot the tracking efficiency vs X integrated over minCent-maxCent% centrality range, for each run and integrated
  /// !!! to be compiled because of CINT !!!
  
  if (!fileNameWeights.IsNull()) LoadRunWeights(fileNameWeights);
  
  TGraphAsymmErrors *intEffVsX = runWeights ? new TGraphAsymmErrors : 0x0;
  Int_t n = -1;
  TArrayD x;
  TArrayD rec;
  TArrayD gen;
  TArrayD effErr[2];
  
  ifstream inFile(runList.Data());
  if (!inFile.is_open()) {
    printf("cannot open file %s\n",runList.Data());
    return;
  }
  
  TString currRun;
  while (!inFile.eof()) {
    
    // get current run number
    currRun.ReadLine(inFile,kTRUE);
    if(currRun.IsNull() || !currRun.IsDec()) continue;
    Int_t run = currRun.Atoi();
    
    printf("run %d: ", run);
    
    // compute efficiency vs var
    TString dataFile = Form("runs/%d/%s", currRun.Atoi(), fileNameData.Data());
    TString outFile = Form("runs/%d/%s", currRun.Atoi(), fileNameSave.Data());
    plotMuonEfficiencyVsX(var, minCent, maxCent, kTRUE, dataFile, outFile);
    
    // get run weight
    if (!runWeights) continue;
    TParameter<Double_t> *weight = static_cast<TParameter<Double_t>*>(runWeights->FindObject(currRun.Data()));
    if (!weight) {
      printf("weight not found for run %s\n", currRun.Data());
      continue;
    }
    Double_t w = weight->GetVal();
    Double_t w2 = w*w;
    
    // get result
    TFile *file = new TFile(outFile.Data(), "read");
    if (!file || !file->IsOpen()) {
      printf("cannot open file\n");
      continue;
    }
    TGraphAsymmErrors *effVsX = static_cast<TGraphAsymmErrors*>(file->FindObjectAny(Form("trackingEffVs%s",var.Data())));
    if (!effVsX) {
      printf("trackingEffVs%s object not found\n", var.Data());
      continue;
    }
    
    // prepare the arrays if not already done
    if (n < 0) {
      n = effVsX->GetN();
      x.Set(n, effVsX->GetX());
      rec.Set(n);
      gen.Set(n);
      effErr[0].Set(n);
      effErr[1].Set(n);
    } else if (n != effVsX->GetN()) {
      printf("number of points in graph trackingEffVs%s for run %d is different than from previous runs\n", var.Data(), run);
      continue;
    }
    
    // loop over bins
    for (Int_t ix = 0; ix < n; ix++) {
      
      Double_t ieff = effVsX->GetY()[ix];
      
      Double_t ieffErr[2] = {effVsX->GetErrorYlow(ix), effVsX->GetErrorYhigh(ix)};
      
      rec[ix] += w*ieff;
      gen[ix] += w;
      effErr[0][ix] += w2*ieffErr[0]*ieffErr[0];
      effErr[1][ix] += w2*ieffErr[1]*ieffErr[1];
      
    }
    
    // close input file
    file->Close();
    
  } 
  inFile.close();
  
  // fill output graph
  for (Int_t ix = 0; ix < n; ix++) {
    
    if (gen[ix] > 0.) {
      
      intEffVsX->SetPoint(ix, x[ix], rec[ix]/gen[ix]);
      intEffVsX->SetPointError(ix, 0., 0., TMath::Sqrt(effErr[0][ix])/gen[ix], TMath::Sqrt(effErr[1][ix])/gen[ix]);
      
    } else {
      
      intEffVsX->SetPoint(ix, x[ix], 1.);
      intEffVsX->SetPointError(ix, 0., 0., 1., 0.);
      
    }
    
  }
  
  // display
  intEffVsX->SetName(Form("integratedTrackingEffVs%s",var.Data()));
  intEffVsX->SetTitle(Form("Integrated tracking efficiency versus %s",var.Data()));
  intEffVsX->SetLineStyle(1);
  intEffVsX->SetLineColor(1); 
  intEffVsX->SetMarkerStyle(20);
  intEffVsX->SetMarkerSize(0.7);
  intEffVsX->SetMarkerColor(2);
  intEffVsX->GetXaxis()->SetTitle(var.Data());
  intEffVsX->GetXaxis()->CenterTitle(kTRUE);
  intEffVsX->GetXaxis()->SetLabelFont(22);
  intEffVsX->GetXaxis()->SetTitleFont(22);
  intEffVsX->GetYaxis()->SetTitle("Efficiency");
  intEffVsX->GetYaxis()->SetLabelFont(22);
  intEffVsX->GetYaxis()->SetTitleFont(22);
  
  // save output
  TFile *file = new TFile(fileNameSave.Data(),"update");
  intEffVsX->Write(0x0, TObject::kOverwrite);
  file->Close();
  
}

//---------------------------------------------------------------------------
void plotMuonEfficiencyVsXY(TString xVar, TString yVar, Double_t minCent, Double_t maxCent, TString fileNameData, TString fileNameSave, Bool_t rap )
{
  /// plot the tracking efficiency vs X,Y integrated over minCent-maxCent% centrality range
  /// !!! to be compiled because of CINT !!!
  
  printf("plotting efficiency versus %s/%s...\n", xVar.Data(), yVar.Data());
  
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
  TFile *file = new TFile(fileNameData.Data(), "read");
  if (!file || !file->IsOpen()) {
    printf("cannot open file %s\n",fileNameData.Data());
    return;
  }
  TList *listTT = static_cast<TList*>(file->FindObjectAny("TotalTracksPerChamber"));
  TList *listTD = static_cast<TList*>(file->FindObjectAny("TracksDetectedPerChamber"));
  THnSparse *TT = static_cast<THnSparse*>(listTT->At(10));
  THnSparse *TD = static_cast<THnSparse*>(listTD->At(10));
  
  
  // set the centrality and Pt range for integration
  setCentPtCh(*TT,minCent,maxCent);  //Added by Javier Martin
  setCentPtCh(*TD,minCent,maxCent);

//  // set the centrality range for integration
//  Int_t lowBin = TT->GetAxis(1)->FindBin(minCent);
//  Int_t upBin = TT->GetAxis(1)->FindBin(maxCent);
//  TT->GetAxis(1)->SetRange(lowBin, upBin);
//  TD->GetAxis(1)->SetRange(lowBin, upBin);
//  
//  // set the pt range for integration
//  lowBin = TT->GetAxis(2)->FindBin(ptMin);
//  upBin = TT->GetAxis(2)->GetNbins();
//  TT->GetAxis(2)->SetRange(lowBin, upBin);
//  TD->GetAxis(2)->SetRange(lowBin, upBin);
  
  // output map
  Int_t nxBins = TT->GetAxis(xDim)->GetNbins();
  Int_t nyBins = TT->GetAxis(yDim)->GetNbins();
  TH2F *effVsXY = new TH2F(Form("trackingEffVs%s-%s",xVar.Data(),yVar.Data()),
			   Form("Measured tracking efficiency versus %s and %s",xVar.Data(),yVar.Data()),
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
      if (TTdraw1->GetEntries() > 0) efficiency = new TGraphAsymmErrors(TDdraw1, TTdraw1,"cpe0");
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
  //new TCanvas(Form("cTrackingEffVs%s-%s",xVar.Data(),yVar.Data()), Form("Measured tracking efficiency versus %s and %s",xVar.Data(),yVar.Data()),700,600);
  effVsXY->GetXaxis()->SetTitle(xVar.Data());
  effVsXY->GetXaxis()->CenterTitle(kTRUE);
  effVsXY->GetXaxis()->SetLabelFont(22);
  effVsXY->GetXaxis()->SetTitleFont(22);
  effVsXY->GetYaxis()->SetTitle(yVar.Data());
  effVsXY->GetYaxis()->CenterTitle(kTRUE);
  effVsXY->GetYaxis()->SetLabelFont(22);
  effVsXY->GetYaxis()->SetTitleFont(22);
  effVsXY->GetZaxis()->SetTitle("Efficiency");
  effVsXY->GetZaxis()->SetLabelFont(22);
  effVsXY->GetZaxis()->SetTitleFont(22);
  //effVsXY->SetMinimum(0.6);
  //effVsXY->SetMaximum(0.95);
  //effVsXY->Draw("surf1");
  
  // save output
  file = new TFile(fileNameSave.Data(),"update");
  effVsXY->Write(0x0, TObject::kOverwrite);
  
  //------------------------------------------------
  if (yDim == 3 && rap)
  {
    TH2F* effVsXYrap = new TH2F();
    TString rapName = Form("trackingEffVs%s-%sRapBins",xVar.Data(),yVar.Data());
    TString rapTitle = Form("Measured tracking efficiency versus %s and %s",xVar.Data(),yVar.Data());
    effVsXYrap->SetTitle(rapTitle.Data());
    effVsXYrap->SetName(rapName.Data());
    
    Double_t xBinEdge[nxBins+1];
    Double_t yBinEdge[nyBins+1];
    
    for (Int_t ybin = 0 ; ybin <= nyBins ; ybin++)
    {
      yBinEdge[ybin] = 2*TMath::ATan(TMath::Exp((effVsXY->GetYaxis()->GetBinLowEdge(ybin+1))));
    }
    for (Int_t xbin = 0 ; xbin <= nxBins ; xbin++)
    {
      xBinEdge[xbin] = effVsXY->GetXaxis()->GetBinLowEdge(xbin+1);
    }
    
    effVsXYrap->SetBins(nxBins,xBinEdge,nyBins,yBinEdge);
    
    for (Int_t xbin = 1 ; xbin <= nxBins ; xbin++)
    {
      for (Int_t ybin = 1 ; ybin <= nyBins ; ybin++)
      {
        effVsXYrap->SetBinContent(xbin,ybin,effVsXY->GetBinContent(xbin,ybin));
      }
    }
    effVsXYrap->Write(0x0, TObject::kOverwrite);

  }
  //----------------------------------------------------
  file->Close();
  
}
//----------------------------------------------------------------------------
void loadOCDB()
{
  //  if (!gGrid)
  //  {
  //    TGrid::Connect("alien://");
  //  }
  
  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  AliCDBManager::Instance()->SetRun(0);
  AliMUONCDB::LoadMapping();
  
  return;
  
}
//---------------------------------------------------------------------------
void  plotMuonEfficiencyVsRun(Double_t minCent, Double_t maxCent, TString runList, TString fileNameSave, TString alienPath)
{
  /// compute the tracking efficiency versus run
  /// !!! to be compiled because of CINT !!!
    
  ///Modify this function because maybe I can put in the same loop the chmbers eff VS and the DE eff VS run
  
  printf("plotting efficiency versus run...\n");
  
  //------------------Added by Javier Martin----------------------------//
  ifstream inFile(runList.Data());
  if (!inFile.is_open()) {
    printf("cannot open file %s\n",runList.Data());
    return;
  }
  
  // output graphs
  TObjArray chamberVSrunGraphs; // 10 graphs, one for each chamber, and 1 for the global efficiency
  TObjArray deVSrunGraphs; // 1 graph per DE
  TObjArray deStationVSrunGraphs; // 1 graph per DE and station
  
  chamberVSrunGraphs.Add(CreateGraph("trackingEffVsRun","Measured tracking efficiency versus run")); // Global efficiency
  
  for ( Int_t ich = 1 ; ich <11 ; ich++ ) // Create Chamber efficiencies vs run
  { 
    chamberVSrunGraphs.Add(CreateGraph("EffCh%dVSrun","Measured Efficiency for Chamber %d VS run",ich));
  }
    
  //Load the mapping for the DE graphs
  loadOCDB();
  AliMpDEIterator deit;
    
  for ( Int_t ich = 0; ich < 10; ++ich ) // Loop over all the DE to initialize and name the graphs for DE vs run
  {
      deit.First(ich);
      
      while ( !deit.IsDone() )
      {
        deVSrunGraphs.Add(CreateGraph("EffDE%dVSrun","Measured Efficiency for DE %d VS run ",deit.CurrentDEId()));
        
        deit.Next();
      }
  }
  
  Int_t station = 1; // Create Stations vs run graphs
  for (Int_t ich = 0 ; ich < 9 ; ich = ich + 2) // Loop over the first chambers in each station
  {
    for ( Int_t iDE = 0 ; iDE < AliMpDEManager::GetNofDEInChamber(ich) ; iDE++)
    {
      if (ich < 6 || moreTrackCandidates) //Stations 1, 2 and 3
      {
        deStationVSrunGraphs.Add(CreateGraph(Form("EffStat%dDE%sVSrun",station,"%d"),"Measured Efficiency for DE %d VS run ",iDE));
      }
      else if (ich < 7)
      {
        deStationVSrunGraphs.Add(CreateGraph("EffStat4-5dDE%dVSrun","Measured Efficiency for DE %d VS run ",iDE));
      }
    }
    station++;
  }
  //-------------------------------------------------------------//
  
  TArrayD chambersEff(11);
  TArrayD chambersEffErr[2];
  chambersEffErr[0].Set(11);
  chambersEffErr[1].Set(11);

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
    //TFile* file = TFile::Open(Form("alien:///alice/cern.ch/user/%s/%09d/AnalysisResults.root",alienPath.Data(),currRun.Atoi()), "read"); //Use this to take the run per run results directly from alien
    TFile* file = (alienPath.IsNull()) ?
      TFile::Open(Form("runs/%d/AnalysisResults.root",currRun.Atoi()), "read") :
      TFile::Open(Form("/alice/cern.ch/user/%s/%09d/AnalysisResults.root",alienPath.Data(),currRun.Atoi()), "read");
    if (!file || !file->IsOpen()) {
      printf("cannot open file runs/%d/AnalysisResults.root\n",currRun.Atoi());
      
      //------Added by Javier Martin-----//
      
      Zero(chamberVSrunGraphs,irun);
      Zero(deVSrunGraphs,irun);
      
      //---------------------------------//
      delete file;
      continue;
    }
    
    TList *listTT = static_cast<TList*>(file->FindObjectAny("TotalTracksPerChamber"));
    TList *listTD = static_cast<TList*>(file->FindObjectAny("TracksDetectedPerChamber"));
    
    //----------------------Added by Javier Martin-------------------//
    // Efficiency per DE vs run
    for (Int_t ich = 0 ; ich <10 ; ich++)
    {
      //Get input hists per DE
      THnSparse *TTDE = static_cast<THnSparse*>(listTT->At(ich)); // Sparse Tracks in ich
      THnSparse *TDDE = static_cast<THnSparse*>(listTD->At(ich));
      
      // set the centrality and Pt range for integration
      setCentPtCh(*TTDE,minCent,maxCent);
      setCentPtCh(*TDDE,minCent,maxCent);
      
      // project efficiency hists
      TH1D *TTDEdraw1 = TTDE->Projection(0,"e");  //Tracks in ich VS DE
      TH1D *TDDEdraw1 = TDDE->Projection(0,"e");
      
      //Compute de efficiency per DE
      TGraphAsymmErrors *efficiencyDE = 0x0; //Efficiency in ich VS DE
      Int_t nDEich = AliMpDEManager::GetNofDEInChamber(ich); //Number of DE in ich
      if (TTDEdraw1->GetEntries() > 0){ efficiencyDE = new TGraphAsymmErrors(TDDEdraw1, TTDEdraw1,"cpe0");}
      delete TTDEdraw1;
      delete TDDEdraw1;
      
      //Get the "global" (1-156) id of the 1st DE in the current chamber
      Int_t shiftGlobalDE = 0;
      
      for (Int_t i = 0 ; i < ich ; i++)
      {
        shiftGlobalDE += AliMpDEManager::GetNofDEInChamber(i);
      }
      
      
      if (!efficiencyDE)
      {
        cout << "run" << currRun.Atoi() << ": empty bin for DEs :" << endl;
        
        deit.First(ich);
        while ( !deit.IsDone() )
        {
          std::cout << deit.CurrentDEId() << std::endl;
          deit.Next();
        }
        
        Zero(deVSrunGraphs,irun,shiftGlobalDE,(shiftGlobalDE + nDEich));
      }
      else
      {
        //Loop over DE in the Chamber ich to fill eff histos
        for (Int_t iDE = 0 ; iDE < nDEich ; iDE++)
        {
          //static_cast<TGraphAsymmErrors*>(deGraphs.At(shiftGlobalDE))->SetPoint(irun,irun,efficiencyDE->GetY()[iDE]);
          //static_cast<TGraphAsymmErrors*>(deGraphs.At(shiftGlobalDE))->SetPointError(irun,0.,0.,efficiencyDE->GetErrorYlow(iDE),efficiencyDE->GetErrorYhigh(iDE));
          static_cast<TGraphAsymmErrors*>(deVSrunGraphs.FindObject(Form("EffDE%dVSrun",100*(ich+1)+iDE)))->SetPoint(irun,irun,efficiencyDE->GetY()[iDE]); //Use FindObject to get the graph by name
          static_cast<TGraphAsymmErrors*>(deVSrunGraphs.FindObject(Form("EffDE%dVSrun",100*(ich+1)+iDE)))->SetPointError(irun,0.,0.,efficiencyDE->GetErrorYlow(iDE),efficiencyDE->GetErrorYhigh(iDE));
          
          //shiftGlobalDE++; //Not necessary if using FindObject to get the graph to fill
        }
      }
      delete efficiencyDE;
    
    }
    
    // Efficiency per station DE vs run
    station = 1;
    for (Int_t ich = 0 ; ich < 9 ; ich = ich + 2) // Loop over the first chambers in each station
    {
      Int_t nDEich = AliMpDEManager::GetNofDEInChamber(ich); //Number of DE in ich
      
      if (ich < 6 || moreTrackCandidates) //Stations 1, 2 and 3
      {
        for (Int_t iDE = 0 ; iDE < nDEich ; iDE++)
        {

          TGraphAsymmErrors* effDE1Graph = static_cast<TGraphAsymmErrors*>(deVSrunGraphs.FindObject(Form("EffDE%dVSrun",100*(ich+1)+iDE)));
          TGraphAsymmErrors* effDE2Graph = static_cast<TGraphAsymmErrors*>(deVSrunGraphs.FindObject(Form("EffDE%dVSrun",100*(ich+2)+iDE)));

          Double_t effDE1 = effDE1Graph->GetY()[irun];
          Double_t effDE2 = effDE2Graph->GetY()[irun];
          Double_t effDE1ErrLow = effDE1Graph->GetErrorYlow(irun);
          Double_t effDE2ErrLow = effDE2Graph->GetErrorYlow(irun);
          Double_t effDE1ErrHigh = effDE1Graph->GetErrorYhigh(irun);
          Double_t effDE2ErrHigh = effDE2Graph->GetErrorYhigh(irun);
          
          Double_t effStatDE = 1.0 - (1.0-effDE1)*(1.0-effDE2);
          Double_t effStatDEErrLow = TMath::Sqrt((1.-effDE1)*(1.-effDE1)*effDE2ErrLow*effDE2ErrLow + (1.-effDE2)*(1.-effDE2)*effDE1ErrLow*effDE1ErrLow);
          Double_t effStatDEErrHigh = TMath::Sqrt((1.-effDE1)*(1.-effDE1)*effDE2ErrHigh*effDE2ErrHigh + (1.-effDE2)*(1.-effDE2)*effDE1ErrHigh*effDE1ErrHigh);
   
          static_cast<TGraphAsymmErrors*>(deStationVSrunGraphs.FindObject(Form("EffStat%dDE%dVSrun",station,iDE)))->SetPoint(irun,irun,effStatDE);
          static_cast<TGraphAsymmErrors*>(deStationVSrunGraphs.FindObject(Form("EffStat%dDE%dVSrun",station,iDE)))->SetPointError(irun,0,0,effStatDEErrLow,effStatDEErrHigh);

        }
      }
      
      else if (ich < 7)
      {
        // get the chambers 7, 8, 9 and 10 efficiency
        
        // Create arrays for the eff of each DE
        TArrayD chambersEffPerDE[26];
        TArrayD chambersEffPerDEErrLow[26];
        TArrayD chambersEffPerDEErrHigh[26];
        
        for (Int_t i = 0 ; i < nDEich ; i++)
        {
          chambersEffPerDE[i].Set(4); // 4 chambers
          chambersEffPerDEErrLow[i].Set(4);
          chambersEffPerDEErrHigh[i].Set(4);
        }
        
        for (Int_t i = 0; i < 4; i++) // Loop over chambers in the last station
        {

          for (Int_t iDE = 0 ; iDE < nDEich ; iDE++) // Loop over DE to get the eff per Chamber and DE
          {
            TGraphAsymmErrors* effCh = static_cast<TGraphAsymmErrors*>(deVSrunGraphs.FindObject(Form("EffDE%dVSrun",100*(ich+1+i)+iDE)));
            
            chambersEffPerDE[iDE][i] = effCh->GetY()[irun];
            chambersEffPerDEErrLow[iDE][i] = effCh->GetErrorYlow(irun);
            chambersEffPerDEErrHigh[iDE][i] = effCh->GetErrorYhigh(irun);
            
          }
          
        }
        
        for (Int_t iDE = 0 ; iDE < nDEich ; iDE++) // Loop over DE to get the station eff per DE
        {
          Double_t effStatDE = chambersEffPerDE[iDE][0] * chambersEffPerDE[iDE][1] * chambersEffPerDE[iDE][2] * chambersEffPerDE[iDE][3]
          +(1.0 - chambersEffPerDE[iDE][0]) * chambersEffPerDE[iDE][1] * chambersEffPerDE[iDE][2] * chambersEffPerDE[iDE][3]
          + chambersEffPerDE[iDE][0] * (1.0 - chambersEffPerDE[iDE][1]) * chambersEffPerDE[iDE][2] * chambersEffPerDE[iDE][3]
          + chambersEffPerDE[iDE][0] * chambersEffPerDE[iDE][1] * (1.0 - chambersEffPerDE[iDE][2]) * chambersEffPerDE[iDE][3]
          + chambersEffPerDE[iDE][0] * chambersEffPerDE[iDE][1] * chambersEffPerDE[iDE][2] * (1.0 - chambersEffPerDE[iDE][3]);
          
          Double_t x7 = (1.-chambersEffPerDE[iDE][1])*chambersEffPerDE[iDE][2]*chambersEffPerDE[iDE][3] + chambersEffPerDE[iDE][1]*(1.-chambersEffPerDE[iDE][2])*chambersEffPerDE[iDE][3] + chambersEffPerDE[iDE][1]*chambersEffPerDE[iDE][2]*(1.-chambersEffPerDE[iDE][3]);
          Double_t x8 = (1.-chambersEffPerDE[iDE][0])*chambersEffPerDE[iDE][2]*chambersEffPerDE[iDE][3] + chambersEffPerDE[iDE][0]*(1.-chambersEffPerDE[iDE][2])*chambersEffPerDE[iDE][3] + chambersEffPerDE[iDE][0]*chambersEffPerDE[iDE][2]*(1.-chambersEffPerDE[iDE][3]);
          Double_t x9 = (1.-chambersEffPerDE[iDE][0])*chambersEffPerDE[iDE][1]*chambersEffPerDE[iDE][3] + chambersEffPerDE[iDE][0]*(1.-chambersEffPerDE[iDE][1])*chambersEffPerDE[iDE][3] + chambersEffPerDE[iDE][0]*chambersEffPerDE[iDE][1]*(1.-chambersEffPerDE[iDE][3]);
          Double_t x10 = (1.-chambersEffPerDE[iDE][0])*chambersEffPerDE[iDE][1]*chambersEffPerDE[iDE][2] + chambersEffPerDE[iDE][0]*(1.-chambersEffPerDE[iDE][1])*chambersEffPerDE[iDE][2] + chambersEffPerDE[iDE][0]*chambersEffPerDE[iDE][1]*(1.-chambersEffPerDE[iDE][2]);
          Double_t effStatDEErrLow = TMath::Sqrt(x7*x7*chambersEffPerDEErrLow[iDE][0]*chambersEffPerDEErrLow[iDE][0] + x8*x8*chambersEffPerDEErrLow[iDE][1]*chambersEffPerDEErrLow[iDE][1] + x9*x9*chambersEffPerDEErrLow[iDE][2]*chambersEffPerDEErrLow[iDE][2] + x10*x10*chambersEffPerDEErrLow[iDE][3]*chambersEffPerDEErrLow[iDE][3]);
          Double_t effStatDEErrHigh = TMath::Sqrt(x7*x7*chambersEffPerDEErrHigh[iDE][0]*chambersEffPerDEErrHigh[iDE][0] + x8*x8*chambersEffPerDEErrHigh[iDE][1]*chambersEffPerDEErrHigh[iDE][1] + x9*x9*chambersEffPerDEErrHigh[iDE][2]*chambersEffPerDEErrHigh[iDE][2] + x10*x10*chambersEffPerDEErrHigh[iDE][3]*chambersEffPerDEErrHigh[iDE][3]);
          
          static_cast<TGraphAsymmErrors*>(deStationVSrunGraphs.FindObject(Form("EffStat4-5dDE%dVSrun",iDE)))->SetPoint(irun,irun,effStatDE);
          static_cast<TGraphAsymmErrors*>(deStationVSrunGraphs.FindObject(Form("EffStat4-5dDE%dVSrun",iDE)))->SetPointError(irun,0,0,effStatDEErrLow,effStatDEErrHigh);
        }
                                          
      }
      station++;
    }
    //-----------------------------------------------------//
    
    //Get input hists per chamber
    THnSparse *TT = static_cast<THnSparse*>(listTT->At(10));
    THnSparse *TD = static_cast<THnSparse*>(listTD->At(10));
    
    // set the centrality and Pt range for integration
    setCentPtCh(*TT,minCent,maxCent);  //Added by Javier Martin
    setCentPtCh(*TD,minCent,maxCent);

    // project efficiency hists
    TH1D *TTdraw1 = TT->Projection(0,"e");
    TH1D *TDdraw1 = TD->Projection(0,"e");
    TGraphAsymmErrors *efficiency = 0x0;
    if (TTdraw1->GetEntries() > 0) efficiency = new TGraphAsymmErrors(TDdraw1, TTdraw1,"cpe0");
    if (!efficiency) {
      cout << "run" << currRun.Atoi() << ": empty bin" << endl;
      //------Added by Javier Martin-----//
      Zero(chamberVSrunGraphs,irun);
      //---------------------------------//
      file->Close();
      delete TTdraw1;
      delete TDdraw1;
      continue;
    }
    
    // get the individual chamber efficiency
    for (Int_t i = 0; i < 10; i++) {
      if (TTdraw1->GetBinContent(i+1) > 0) {
        chambersEff[i+1] = efficiency->GetY()[i];
        chambersEffErr[0][i+1] = efficiency->GetErrorYlow(i);
        chambersEffErr[1][i+1] = efficiency->GetErrorYhigh(i);
      } else {
        chambersEff[i+1] = -1.;
        chambersEffErr[0][i+1] = 0.;
        chambersEffErr[1][i+1] = 0.;
      }
    }
    delete TTdraw1;
    delete TDdraw1;
    delete efficiency;
      
    // compute the overall tracking efficiency
//    cout << "run" << currRun.Atoi() << ":" << endl;
//    computeTrackingEfficiency(chambersEff, chambersEffErr, kTRUE);
    computeTrackingEfficiency(chambersEff, chambersEffErr);
//
    // fill graph
    //------Added by Javier Martin-----//
    for (Int_t ich = 0 ; ich <11 ; ich++) {
        static_cast<TGraphAsymmErrors*>(chamberVSrunGraphs.At(ich))->SetPoint(irun,irun,chambersEff[ich]);
        static_cast<TGraphAsymmErrors*>(chamberVSrunGraphs.At(ich))->SetPointError(irun,0.,0.,chambersEffErr[0][ich],chambersEffErr[1][ich]);
    }
    //---------------------------------//
    
    file->Close();
    
    //Compute efficiency per DE in a Station and per run
    
  } 
  inFile.close();
  
  //Beautify graphs
  BeautifyGraphs(deVSrunGraphs,"Run Number","Efficiency");
  BeautifyGraphs(chamberVSrunGraphs,"Run Number","Efficiency");
  BeautifyGraphs(deStationVSrunGraphs,"Run Number","Efficiency");
  
  // set bin labels
  SetRunLabel(deVSrunGraphs,irun,runs);
  SetRunLabel(chamberVSrunGraphs,irun,runs);
  SetRunLabel(deStationVSrunGraphs,irun,runs);
  
  // save output
  TFile* file = new TFile(fileNameSave.Data(),"update");
  //------Added by Javier Martin-----//
  
  chamberVSrunGraphs.Write("ChambersEffVSrun", TObject::kOverwrite | TObject::kSingleKey);
  deVSrunGraphs.Write("DEEffVSrun", TObject::kOverwrite | TObject::kSingleKey);
  deStationVSrunGraphs.Write("StationsEffPerDEVSrun", TObject::kOverwrite | TObject::kSingleKey);

  //---------------------------------//
  file->Close();
  
}

//---------------------------------------------------------------------------
void plotMuonEfficiencyDEperDE(Double_t minCent, Double_t maxCent, TString fileNameData, TString fileNameSave) //Added by Javier Martin
{
  
  printf("plotting efficiency DE per DE...\n");
  
  // get input hists
  TFile *file = new TFile(fileNameData.Data(), "read");
  if (!file || !file->IsOpen()) {
    printf("cannot open file %s \n",fileNameData.Data());
    return;
  }
  //Load mapping
  loadOCDB();
  
  TList *listTT = static_cast<TList*>(file->FindObjectAny("TotalTracksPerChamber"));
  TList *listTD = static_cast<TList*>(file->FindObjectAny("TracksDetectedPerChamber"));
  
  //-------Chambers efficiency per DE
  
  TObjArray chamberGraphs;
  
  for (Int_t ich = 0 ; ich <10 ; ich++) // ich = chamber ID - 1
  {
    //Create output graphs
    chamberGraphs.Add(CreateGraph("Ch%dEffVSDE","Measured Efficiency for Chamber %d",ich+1));
    
    //Get input hists per DE
    THnSparse *TTDE = static_cast<THnSparse*>(listTT->At(ich)); // Sparse Tracks in ich
    THnSparse *TDDE = static_cast<THnSparse*>(listTD->At(ich));
    
    // project efficiency hists integrated over the centrality bins
    setCentPtCh(*TTDE,minCent,maxCent);
    setCentPtCh(*TDDE,minCent,maxCent);
    TH1D *TTDEdraw1 = TTDE->Projection(0,"e");  //Tracks in ich VS DE
    TH1D *TDDEdraw1 = TDDE->Projection(0,"e");
    
    TGraphAsymmErrors *efficiency = 0x0;
    if (TTDEdraw1->GetEntries() > 0) efficiency = new TGraphAsymmErrors(TDDEdraw1, TTDEdraw1,"cpe0");
    delete TTDEdraw1;
    delete TDDEdraw1;
    
    Int_t nDEich = efficiency->GetN(); // Number of DE (# points in the eff graph) in the ich chamber
    
    for (Int_t iDE = 0 ; iDE < nDEich ; iDE++)
    {
      static_cast<TGraphAsymmErrors*>(chamberGraphs.At(ich))->SetPoint(iDE,iDE,efficiency->GetY()[iDE]);
      static_cast<TGraphAsymmErrors*>(chamberGraphs.At(ich))->SetPointError(iDE,0,0,efficiency->GetErrorYlow(iDE),efficiency->GetErrorYhigh(iDE));
    }
    delete efficiency;
    
    static_cast<TGraphAsymmErrors*>(chamberGraphs.At(ich))->GetXaxis()->Set(nDEich+1, -0.5, nDEich-0.5); //Set the length of the x axis graph for an appropiate display of DE
    static_cast<TGraphAsymmErrors*>(chamberGraphs.At(ich))->GetXaxis()->SetNdivisions(nDEich);
  }
  file->Close();
  
  //----- Stations efficiency per DE
  
  TObjArray stationGraphs;
  Int_t station = 1;
  
  for (Int_t ich = 0 ; ich < 9 ; ich = ich + 2) // Loop over the first chambers in each station 
  {
    if (ich < 6 || moreTrackCandidates) //Stations 1, 2 and 3
    {
      stationGraphs.Add(CreateGraph("St%dEffVSDE","Measured Efficiency for Station %d",station));
      
      TGraphAsymmErrors* effChiGraph = static_cast<TGraphAsymmErrors*>(chamberGraphs.At(ich));
      TGraphAsymmErrors* effChjGraph = static_cast<TGraphAsymmErrors*>(chamberGraphs.At(ich+1));
      
      Int_t nDEich = effChiGraph->GetN(); // Number of DE (# points in the eff graph) in the ich chamber
      
      for (Int_t iDE = 0 ; iDE < nDEich ; iDE++)
      {
        Double_t effDEChi = effChiGraph->GetY()[iDE];
        Double_t effDEChj = effChjGraph->GetY()[iDE];
        Double_t effDEChiErrLow = effChiGraph->GetErrorYlow(iDE);
        Double_t effDEChjErrLow = effChjGraph->GetErrorYlow(iDE);
        Double_t effDEChiErrHigh = effChiGraph->GetErrorYhigh(iDE);
        Double_t effDEChjErrHigh = effChjGraph->GetErrorYhigh(iDE);
        
        Double_t effStatDE = 1.0 - (1.0-effDEChi)*(1.0-effDEChj);
        Double_t effStatDEErrLow = TMath::Sqrt((1.-effDEChi)*(1.-effDEChi)*effDEChjErrLow*effDEChjErrLow + (1.-effDEChj)*(1.-effDEChj)*effDEChiErrLow*effDEChiErrLow);
        Double_t effStatDEErrHigh = TMath::Sqrt((1.-effDEChi)*(1.-effDEChi)*effDEChjErrHigh*effDEChjErrHigh + (1.-effDEChj)*(1.-effDEChj)*effDEChiErrHigh*effDEChiErrHigh);
        
        static_cast<TGraphAsymmErrors*>(stationGraphs.At(station-1))->SetPoint(iDE,iDE,effStatDE);
        static_cast<TGraphAsymmErrors*>(stationGraphs.At(station-1))->SetPointError(iDE,0,0,effStatDEErrLow,effStatDEErrHigh);
    
      }
      
      static_cast<TGraphAsymmErrors*>(stationGraphs.At(station-1))->GetXaxis()->Set(nDEich+1, -0.5, nDEich-0.5); //Set the length of the x axis for an appropiate display of DE
      static_cast<TGraphAsymmErrors*>(stationGraphs.At(station-1))->GetXaxis()->SetNdivisions(nDEich);
      
      station++;
    }
    
    else if (ich < 7) // Stations 4 and 5
    {
      stationGraphs.Add(CreateGraph("St4-5EffVSDE","Measured Efficiency for Station 4-5"));
      
      // get the chambers 7, 8, 9 and 10 efficiency
      
      // Create arrays for the eff of each DE
      Int_t nDEichArr = AliMpDEManager::GetNofDEInChamber(9);; // 26 DE in each chamber
      TArrayD chambersEff[26];
      TArrayD chambersEffErrLow[26];
      TArrayD chambersEffErrHigh[26];
      
      for (Int_t i = 0 ; i < nDEichArr ; i++)
      {
        chambersEff[i].Set(4); // 4 chambers
        chambersEffErrLow[i].Set(4);
        chambersEffErrHigh[i].Set(4);
      }
      
      for (Int_t i = 0; i < 4; i++) // Loop over chambers
      {
        TGraphAsymmErrors* effCh = static_cast<TGraphAsymmErrors*>(chamberGraphs.At(ich+i));
        Int_t nDEich = effCh->GetN(); // Number of DE (# points in the eff graph) in the ich chamber
        
        if( nDEich != nDEichArr )
        {
          printf("Error, number of DE not correct in chamberGraphs");
          return;
        }
        
        for (Int_t iDE = 0 ; iDE < nDEich ; iDE++) // Loop over DE to get the eff per Chamber and DE
        {
        chambersEff[iDE][i] = effCh->GetY()[iDE];
        chambersEffErrLow[iDE][i] = effCh->GetErrorYlow(iDE);
        chambersEffErrHigh[iDE][i] = effCh->GetErrorYhigh(iDE);
        }
        
      }
      
      for (Int_t iDE = 0 ; iDE < nDEichArr ; iDE++) // Loop over DE to get the station eff per DE
      {
        Double_t effStatDE = chambersEff[iDE][0] * chambersEff[iDE][1] * chambersEff[iDE][2] * chambersEff[iDE][3]
        +(1.0 - chambersEff[iDE][0]) * chambersEff[iDE][1] * chambersEff[iDE][2] * chambersEff[iDE][3]
        + chambersEff[iDE][0] * (1.0 - chambersEff[iDE][1]) * chambersEff[iDE][2] * chambersEff[iDE][3]
        + chambersEff[iDE][0] * chambersEff[iDE][1] * (1.0 - chambersEff[iDE][2]) * chambersEff[iDE][3]
        + chambersEff[iDE][0] * chambersEff[iDE][1] * chambersEff[iDE][2] * (1.0 - chambersEff[iDE][3]);
        
        Double_t x7 = (1.-chambersEff[iDE][1])*chambersEff[iDE][2]*chambersEff[iDE][3] + chambersEff[iDE][1]*(1.-chambersEff[iDE][2])*chambersEff[iDE][3] + chambersEff[iDE][1]*chambersEff[iDE][2]*(1.-chambersEff[iDE][3]);
        Double_t x8 = (1.-chambersEff[iDE][0])*chambersEff[iDE][2]*chambersEff[iDE][3] + chambersEff[iDE][0]*(1.-chambersEff[iDE][2])*chambersEff[iDE][3] + chambersEff[iDE][0]*chambersEff[iDE][2]*(1.-chambersEff[iDE][3]);
        Double_t x9 = (1.-chambersEff[iDE][0])*chambersEff[iDE][1]*chambersEff[iDE][3] + chambersEff[iDE][0]*(1.-chambersEff[iDE][1])*chambersEff[iDE][3] + chambersEff[iDE][0]*chambersEff[iDE][1]*(1.-chambersEff[iDE][3]);
        Double_t x10 = (1.-chambersEff[iDE][0])*chambersEff[iDE][1]*chambersEff[iDE][2] + chambersEff[iDE][0]*(1.-chambersEff[iDE][1])*chambersEff[iDE][2] + chambersEff[iDE][0]*chambersEff[iDE][1]*(1.-chambersEff[iDE][2]);
        Double_t effStatDEErrLow = TMath::Sqrt(x7*x7*chambersEffErrLow[iDE][0]*chambersEffErrLow[iDE][0] + x8*x8*chambersEffErrLow[iDE][1]*chambersEffErrLow[iDE][1] + x9*x9*chambersEffErrLow[iDE][2]*chambersEffErrLow[iDE][2] + x10*x10*chambersEffErrLow[iDE][3]*chambersEffErrLow[iDE][3]);
        Double_t effStatDEErrHigh = TMath::Sqrt(x7*x7*chambersEffErrHigh[iDE][0]*chambersEffErrHigh[iDE][0] + x8*x8*chambersEffErrHigh[iDE][1]*chambersEffErrHigh[iDE][1] + x9*x9*chambersEffErrHigh[iDE][2]*chambersEffErrHigh[iDE][2] + x10*x10*chambersEffErrHigh[iDE][3]*chambersEffErrHigh[iDE][3]);
        
        static_cast<TGraphAsymmErrors*>(stationGraphs.At(station-1))->SetPoint(iDE,iDE,effStatDE);
        static_cast<TGraphAsymmErrors*>(stationGraphs.At(station-1))->SetPointError(iDE,0,0,effStatDEErrLow,effStatDEErrHigh);
      }
      
      static_cast<TGraphAsymmErrors*>(stationGraphs.FindObject("St4-5EffVSDE"))->GetXaxis()->Set(nDEichArr+1, -0.5, nDEichArr-0.5); //Set the length of the x axis for an appropiate display of DE
      static_cast<TGraphAsymmErrors*>(stationGraphs.FindObject("St4-5EffVSDE"))->GetXaxis()->SetNdivisions(nDEichArr);
      
    }
    
  }
  
  
  
  //Beautify graphs
  BeautifyGraphs(chamberGraphs,"Detection Element","Efficiency");
  BeautifyGraphs(stationGraphs,"Detection Element","Efficiency");
 
  
  //Save Output
  TFile* savefile = new TFile(fileNameSave.Data(),"update");
  chamberGraphs.Write("ChambersEffVSDE", TObject::kOverwrite | TObject::kSingleKey);
  stationGraphs.Write("StationsEffVSDE", TObject::kOverwrite | TObject::kSingleKey);
  savefile->Close();
  
}

//---------------------------------------------------------------------------
void integratedEfficiency(Double_t minCent, Double_t maxCent, TString fileNameData)
{
  /// compute the integrated tracking efficiency in minCent-maxCent% centrality range
  /// !!! to be compiled because of CINT !!!
  
  // adjust centrality range
  maxCent = TMath::Max(maxCent-1.e-12, minCent);
  
  // get input hists
  TFile *file = new TFile(fileNameData.Data(), "read");
  if (!file || !file->IsOpen()) {
    printf("cannot open file\n");
    return;
  }
  TList *listTT = static_cast<TList*>(file->FindObjectAny("TotalTracksPerChamber"));
  TList *listTD = static_cast<TList*>(file->FindObjectAny("TracksDetectedPerChamber"));
  THnSparse *TT = static_cast<THnSparse*>(listTT->At(10));
  THnSparse *TD = static_cast<THnSparse*>(listTD->At(10));
  
  // get the centrality range for integration
  Int_t lowBin = TT->GetAxis(1)->FindBin(minCent);
  Int_t upBin = TT->GetAxis(1)->FindBin(maxCent);
  cout << endl << "Integrated efficiency in "<< TT->GetAxis(1)->GetBinLowEdge(lowBin) << "-" << TT->GetAxis(1)->GetBinUpEdge(upBin) << "%:" << endl;
  setCentPtCh(*TT,minCent,maxCent);
  setCentPtCh(*TD,minCent,maxCent);
  
  // project efficiency hists integrated over this centrality range
  TH1D *TTdraw1 = TT->Projection(0,"e");
  TH1D *TDdraw1 = TD->Projection(0,"e");
  TGraphAsymmErrors *efficiency = 0x0;
  if (TTdraw1->GetEntries() > 0) efficiency = new TGraphAsymmErrors(TDdraw1, TTdraw1,"cpe0");
  if (!efficiency) {
    cout << "integrated efficiency: empty bin" << endl;
    delete TTdraw1;
    delete TDdraw1;
    return;
  }
  
  // get the individual chamber efficiency
  TArrayD chambersEff(11);
  TArrayD chambersEffErr[2];
  chambersEffErr[0].Set(11);
  chambersEffErr[1].Set(11);
  for (Int_t i = 0; i < 10; i++) {
    if (TTdraw1->GetBinContent(i+1) > 0) {
      chambersEff[i+1] = efficiency->GetY()[i];
      chambersEffErr[0][i+1] = efficiency->GetErrorYlow(i);
      chambersEffErr[1][i+1] = efficiency->GetErrorYhigh(i);
    } else {
      chambersEff[i+1] = -1.;
      chambersEffErr[0][i+1] = 0.;
      chambersEffErr[1][i+1] = 0.;
    }
  }
  delete TTdraw1;
  delete TDdraw1;
  delete efficiency;
  
  // compute the overall tracking efficiency
  computeTrackingEfficiency(chambersEff, chambersEffErr, kTRUE);
  
  // close input file
  file->Close();
  
}

//---------------------------------------------------------------------------
void integratedEfficiency(TString fileNameWeights, TString fileNameSave)
{
  /// compute the integrated tracking efficiency from run per run chamber efficiencies
  /// !!! to be compiled because of CINT !!!
  
  LoadRunWeights(fileNameWeights);
  if (!runWeights) return;
  
  // get input hists
  TFile *file = new TFile(fileNameSave.Data(), "read");
  if (!file || !file->IsOpen()) {
    printf("cannot open file\n");
    return;
  }
  TObjArray *chamberVSrunGraphs = static_cast<TObjArray*>(file->FindObjectAny("ChambersEffVSrun"));
  if (!chamberVSrunGraphs) {
    printf("ChambersEffVSrun object not found\n");
    return;
  }
  
  TGraphAsymmErrors *efficiency[11];
  Double_t rec[11], gen[11], eff[11], effErr[2][11];
  for (Int_t ich = 0; ich < 11; ich++) {
    efficiency[ich] = static_cast<TGraphAsymmErrors*>(chamberVSrunGraphs->At(ich));
    rec[ich] = 0.;
    gen[ich] = 0.;
    eff[ich] = 0.;
    effErr[0][ich] = 0.;
    effErr[1][ich] = 0.;
  }
  
  // loop over runs
  Int_t nRuns = efficiency[0]->GetN();
  for (Int_t iRun = 0; iRun < nRuns; iRun++) {
    
    // get run weight
    TString sRun = efficiency[0]->GetXaxis()->GetBinLabel(iRun+1);
    TParameter<Double_t> *weight = static_cast<TParameter<Double_t>*>(runWeights->FindObject(sRun.Data()));
    if (!weight) {
      printf("weight not found for run %s\n", sRun.Data());
      continue;
    }
    Double_t w = weight->GetVal();
    Double_t w2 = w*w;
    
    // loop over chamber
    for (Int_t ich = 0; ich < 11; ich++) {
      
      Double_t effCh = efficiency[ich]->GetY()[iRun];
      Double_t errCh[2] = {efficiency[ich]->GetErrorYlow(iRun), efficiency[ich]->GetErrorYhigh(iRun)};
      if (effCh < 0.) {
	printf("no efficiency measurement for chamber %d in run %s --> use 1 +0 -1\n", ich, sRun.Data());
	effCh = 1.;
	errCh[0] = 1.;
	errCh[1] = 0.;
      }      
      
      rec[ich] += w*effCh;
      gen[ich] += w;
      effErr[0][ich] += w2*errCh[0]*errCh[0];
      effErr[1][ich] += w2*errCh[1]*errCh[1];
      
    }
    
  }
  
  // compute efficiencies
  cout << endl << "Integrated efficiency weighted run per run:" << endl;
  for (Int_t ich = 0 ; ich <11 ; ich++) {
    
    if (gen[ich] > 0.) {
      
      eff[ich] = rec[ich]/gen[ich];
      effErr[0][ich] = TMath::Sqrt(effErr[0][ich])/gen[ich];
      effErr[1][ich] = TMath::Sqrt(effErr[1][ich])/gen[ich];
      
    } else {
      
      eff[ich] = 1.;
      effErr[0][ich] = 1.;
      effErr[1][ich] = 0.;
      
    }
    
  }
  
  // print results
  for (Int_t i = 1; i <= 10; i++) {
    cout << "Efficiency chamber " << i << " : " << eff[i] << " + " << effErr[1][i] << " - " << effErr[0][i] << endl;
  }
  cout << "Total tracking efficiency : " << eff[0] << " + " << effErr[1][0] << " - " << effErr[0][0] << endl;
  
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
  
  if (moreTrackCandidates) chambersEff[0] = st1Eff * st2Eff * st3Eff * st4Eff * st5Eff;
  else chambersEff[0] = st1Eff * st2Eff * st3Eff * st45Eff;
  
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
    
    if (moreTrackCandidates) chambersEffErr[i][0] = chambersEff[0]*TMath::Sqrt(st1EffErr[i]*st1EffErr[i]/st1Eff/st1Eff + st2EffErr[i]*st2EffErr[i]/st2Eff/st2Eff + st3EffErr[i]*st3EffErr[i]/st3Eff/st3Eff + st4EffErr[i]*st4EffErr[i]/st4Eff/st4Eff + st5EffErr[i]*st5EffErr[i]/st5Eff/st5Eff);
    else chambersEffErr[i][0] = chambersEff[0]*TMath::Sqrt(st1EffErr[i]*st1EffErr[i]/st1Eff/st1Eff + st2EffErr[i]*st2EffErr[i]/st2Eff/st2Eff + st3EffErr[i]*st3EffErr[i]/st3Eff/st3Eff + st45EffErr[i]*st45EffErr[i]/st45Eff/st45Eff);
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

//---------------------------------------------------------------------------
void LoadRunWeights(TString fileName)
{
  /// Set the number of interested events per run
  /// (used to weight the acc*eff correction integrated
  /// over run for any pt/y/centrality bins)
  
  if (runWeights) return;
  
  ifstream inFile(fileName.Data());
  if (!inFile.is_open()) {
    printf("cannot open file %s", fileName.Data());
    return;
  }
  
  runWeights = new THashList(1000);
  runWeights->SetOwner();
  
  TString line;
  while (! inFile.eof() ) {
    
    line.ReadLine(inFile,kTRUE);
    if(line.IsNull()) continue;
    
    TObjArray *param = line.Tokenize(" ");
    if (param->GetEntries() != 2) {
      printf("bad input line %s", line.Data());
      continue;
    }
    
    Int_t run = ((TObjString*)param->UncheckedAt(0))->String().Atoi();
    if (run < 0) {
      printf("invalid run number: %d", run);
      continue;
    }
    
    Float_t weight = ((TObjString*)param->UncheckedAt(1))->String().Atof();
    if (weight <= 0.) {
      printf("invalid weight: %g", weight);
      continue;
    }
    
    runWeights->Add(new TParameter<Double_t>(Form("%d",run), weight));
    
    delete param;
  }
  
  inFile.close();
  
}

