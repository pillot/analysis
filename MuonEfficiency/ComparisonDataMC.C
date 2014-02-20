//
//  CoparisonDataMC.C
//  
//
//  Created by Javier Martin-Blanco on 25/03/13.
//
//

#include <stdio.h>
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
#include <TGaxis.h>
#include <TLegend.h>

#include "AliMpDEIterator.h"
#include "AliMUONCDB.h"
#include "AliCDBManager.h"
#include "AliMpDEManager.h"

#include <TCanvas.h>
//#include <TTree.h>
//#include <TChain.h>


#include <TFitResult.h>
#include <TH2.h>
#include <TProfile.h>
#include <TLegend.h>
#include <TList.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TString.h>

//This task is intended to run over the results of RealEfficiency.C and it is to compare the results for Data and MC simulation

TGraphAsymmErrors* CreateRatioGraph(const char* name, const char* title, TGraphAsymmErrors& GraphData, TGraphAsymmErrors& GraphSim) // Change this function as in RealEfficiency.C (value)
{
  if (GraphData.GetN() != GraphSim.GetN() )
  {
    printf("Error, dividing graphs of different entries number\n");
    return 0x0;
  }
  
  // compute efficiency ratio data/sim
  Int_t nBins = GraphData.GetN();
  TGraphAsymmErrors *ratio = new TGraphAsymmErrors(nBins);
  Double_t x,effD,effS,effDErrh,effDErrl,effSErrh,effSErrl,rat,ratErrh,ratErrl;
  
  // Set ratio's Xaxis lenght, name and title
  
  ratio->SetName(name);
  ratio->SetTitle(title);
  
  // Loop over bins
  for (Int_t i = 0; i < nBins; i++)
  {
    // Get points and errors from individual efficiencies
    GraphData.GetPoint(i,x,effD);
    GraphSim.GetPoint(i,x,effS);
    effDErrh = GraphData.GetErrorYhigh(i);
    effDErrl = GraphData.GetErrorYlow(i);
    effSErrh = GraphSim.GetErrorYhigh(i);
    effSErrl = GraphSim.GetErrorYlow(i);
    
    // Compute the ratio and assym errors
    if (effD > 0. && effS > 0.)
    {
      rat = effD/effS;
      ratErrh = rat*TMath::Sqrt(effDErrh*effDErrh/effD*effD + effSErrl*effSErrl/effS*effS);
      ratErrl = rat*TMath::Sqrt(effDErrl*effDErrl/effD*effD + effSErrh*effSErrh/effS*effS);
    }
    if (effD == 0 && effS == 0)
    {
      rat = 1.;
      ratErrh = 0.;
      ratErrl = 0.;
    }
    if (effD == 0 && effS > 0.)
    {
      rat = 0.;
      ratErrh = 0.;
      ratErrl = 0.;
    }
    if (effD > 0. && effS == 0)
    {
      rat = 2.;
      ratErrh = 0.;
      ratErrl = 0.;
    }
    // Fill the ratio Graph
    ratio->SetPoint(i,x,rat);
    ratio->SetPointError(i,0.,0.,ratErrl,ratErrh);
    
//     //------------------------------------------------
//    if ( !strncmp(name,"RatioEffVsRun",15) )
//    {
//      cout << "run" << GraphData.GetXaxis()->GetBinLabel(i+1) << ":" << endl;
//      cout << "Ratio tracking efficiency : " << rat << " + " << ratErrh << " - " << ratErrl << endl;
//     
//    }
//     //------------------------------------------------
  }

  ratio->GetXaxis()->SetTitle(GraphData.GetXaxis()->GetTitle());
  ratio->GetYaxis()->SetTitle("EffData/EffSim");

  if ( !strncmp(GraphData.GetXaxis()->GetTitle(),"Run Number",10) )
  {
    ratio->GetXaxis()->Set(nBins, -0.5, nBins-0.5);
  
    
    for (Int_t i = 1; i <= nBins; i++)
    {
      // set bin labels
      ratio->GetXaxis()->SetBinLabel(i, GraphData.GetXaxis()->GetBinLabel(i));
      
    }
  }
  
  if ( !strncmp(GraphData.GetXaxis()->GetTitle(),"Detection Element",10) )
  {
    ratio->GetXaxis()->Set(nBins+1, -0.5, nBins-0.5); //Set the length of the x axis for an appropiate display of DE
    ratio->GetXaxis()->SetNdivisions(nBins);    
  }
  
  
  
  return ratio;

}

//-------------------------------------------------------------------------------------------------------------------------------------------
void BeautifyGraphs(TObjArray& array, const char* yAxisName) // Is differrent from the one in RealEfficiency.C since now the xaxis name is set in CreateRatioGraph
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
    g->GetXaxis()->CenterTitle(kTRUE);
    g->GetXaxis()->SetLabelFont(22);
    g->GetXaxis()->SetTitleFont(22);
    g->GetYaxis()->SetTitle(yAxisName);
    g->GetYaxis()->CenterTitle(kTRUE);
    g->GetYaxis()->SetLabelFont(22);
    g->GetYaxis()->SetTitleFont(22);
  }

}

// Draw
//---------------------------------------------------------------------------
TCanvas* DrawRatio(TString name, TString title, TGraphAsymmErrors* GraphData, TGraphAsymmErrors * GraphSim, TGraphAsymmErrors* GraphRatio, TString dataName = "Data", TString simName = "Simulation")
{ 
  Float_t fracOfHeight = 0.3;
  Float_t rightMargin = 0.03;
  TCanvas *c = new TCanvas(name.Data(), "ratio",1200,800);

  c->Divide(1,2,0,0);
  
  c->cd(1);
  
  gPad->SetPad(0., fracOfHeight, 0.99, 0.99);
  gPad->SetTopMargin(0.03);
  gPad->SetRightMargin(rightMargin);
  
  Double_t minData = GraphData->GetYaxis()->GetXmin();
  Double_t minSim = GraphSim->GetYaxis()->GetXmin();
  Double_t maxData = GraphData->GetYaxis()->GetXmax();
  Double_t maxSim = GraphSim->GetYaxis()->GetXmax();
  Double_t minA = TMath::Min(minData,minSim);
  if (minA < 0.1) minA -= 0.1;
  Double_t maxA = TMath::Max(maxData,maxSim);
  
  GraphData->SetMinimum(minA);
  GraphData->SetMaximum(maxA);
  //GraphData->SetMinimum(minA-0.1*minA);
  //GraphData->SetMaximum(maxA+0.1*maxA);
  GraphData->SetTitle(title.Data());
  GraphData->GetYaxis()->SetLabelSize(0.051);
  GraphData->GetYaxis()->SetTitleSize(0.051);
  GraphData->GetYaxis()->SetTitleOffset(0.8);
  GraphData->GetYaxis()->CenterTitle(kTRUE);
  GraphData->GetXaxis()->SetLabelOffset(0.1);
  GraphData->SetMarkerStyle(20);
  GraphData->SetMarkerSize(0.6);
  GraphData->SetMarkerColor(4);
  GraphData->SetLineColor(4);
  GraphData->Draw("AP");
  
  
  GraphSim->SetMarkerStyle(20);
  GraphSim->SetMarkerSize(0.6);
  GraphSim->SetMarkerColor(2);
  GraphSim->SetLineColor(2);
  GraphSim->Draw("Psame");
   //------
//  GraphSim.Draw("same");
//  
  TLegend *legend = new TLegend (0.8, 0.8, 0.95, 0.95);
  legend->SetTextSize(0.06);
  legend->AddEntry(GraphSim, Form(" %s",simName.Data()), "ep");
  legend->AddEntry(GraphData, Form(" %s",dataName.Data()), "ep");
  legend->Draw("same");
//

//  GraphRatio->Draw("psamey+");
  
  c->cd(2);
  
  gPad->SetPad(0., 0., 0.99, fracOfHeight);
  gPad->SetRightMargin(rightMargin);
  gPad->SetBottomMargin(0.08/fracOfHeight);
  gPad->SetGridy();
  
  //------
//  GraphRatio.SetName(Form("%s_over_%s",dataName.Data(),simName.Data()));
////  GraphRatio.SetStats(0);
//  GraphRatio.SetTitle("");

  GraphRatio->SetLineStyle(1);
  GraphRatio->SetLineColor(1);
  GraphRatio->SetMarkerStyle(20);
  GraphRatio->SetMarkerSize(0.4);
  GraphRatio->SetMarkerColor(1);
  GraphRatio->GetXaxis()->SetLabelSize(0.11);
  GraphRatio->GetXaxis()->SetTitleSize(0.12);
//  GraphRatio.GetXaxis()->SetTitle("GeV/c   ");
//  GraphRatio->GetXaxis()->SetTitleOffset(-0.6);
  GraphRatio->GetXaxis()->CenterTitle(kTRUE);
  GraphRatio->GetXaxis()->SetLabelFont(22);
  GraphRatio->GetXaxis()->SetTitleFont(22);
  
  GraphRatio->GetYaxis()->SetLabelSize(0.07);
  GraphRatio->GetYaxis()->SetTitleSize(0.07);
  GraphRatio->GetYaxis()->CenterTitle(kTRUE);
  GraphRatio->GetYaxis()->SetTitleOffset(0.37);
//  GraphRatio->GetYaxis()->SetLabelFont(100);
  GraphRatio->GetYaxis()->SetLabelFont(22);
//  GraphRatio.SetMinimum(0.86);
//  GraphRatio.SetMaximum(1.14);
  GraphRatio->Draw("ap");
  //------
//  GraphRatio.Draw("same");
//  
//  
  TLegend *legend2 = new TLegend (0.70, 0.3, 0.95, 0.37);
  legend2->AddEntry(GraphRatio, Form(" %s / %s ",dataName.Data(), simName.Data()), "ep");
  legend2->Draw("same");
//
  //gStyle->SetFillColor(0);
  c->Update();
  //------
  
  return c;
}

//-------------------------------------------------------------------------------------------------------------------------------------------
void ComparisonDataMC(TString fileNameData, TString fileNameSim)
{
  // Open input Data files
  TFile *fileData = new TFile(fileNameData.Data(), "read");
  if (!fileData || !fileData->IsOpen()) {
    printf("cannot open file %s \n",fileNameData.Data());
    return;
  }
  // Open input Sim files
  TFile *fileSim = new TFile(fileNameSim.Data(), "read");
  if (!fileSim || !fileSim->IsOpen()) {
    printf("cannot open file %s \n",fileNameSim.Data());
    return;
  }
  
  // Get Global Data and Sim graphs
  
//  TGraphAsymmErrors *effVScentData = static_cast<TGraphAsymmErrors*>(fileData->FindObjectAny("trackingEffVscentrality"));
//  if (!effVScentData) {
//    printf("Efficiency vs centrality from data not found\n");
//    return;
//  }
//  TGraphAsymmErrors *effVScentSim = static_cast<TGraphAsymmErrors*>(fileSim->FindObjectAny("trackingEffVscentrality"));
//  if (!effVScentSim) {
//    printf("Efficiency vs centrality from sim not found\n");
//    return;
//  }
  
  TGraphAsymmErrors *effVSyData = static_cast<TGraphAsymmErrors*>(fileData->FindObjectAny("trackingEffVsy"));
  if (!effVSyData) {
    printf("Efficiency vs rapidity from data not found\n");
    return;
  }
  TGraphAsymmErrors *effVSySim = static_cast<TGraphAsymmErrors*>(fileSim->FindObjectAny("trackingEffVsy"));
  if (!effVSySim) {
    printf("Efficiency vs rapidity from sim not found\n");
    return;
  }

  
  TGraphAsymmErrors *effVSptData = static_cast<TGraphAsymmErrors*>(fileData->FindObjectAny("trackingEffVspt"));
  if (!effVSptData) {
    printf("Efficiency vs pt from data not found\n");
    return;
  }
  TGraphAsymmErrors *effVSptSim = static_cast<TGraphAsymmErrors*>(fileSim->FindObjectAny("trackingEffVspt"));
  if (!effVSptSim) {
    printf("Efficiency vs pt from sim not found\n");
    return;
  }
  
  TGraphAsymmErrors *effVSphiData = static_cast<TGraphAsymmErrors*>(fileData->FindObjectAny("trackingEffVsphi"));
  if (!effVSphiData) {
    printf("Efficiency vs phi from data not found\n");
    return;
  }
  TGraphAsymmErrors *effVSphiSim = static_cast<TGraphAsymmErrors*>(fileSim->FindObjectAny("trackingEffVsphi"));
  if (!effVSphiSim) {
    printf("Efficiency vs phi from sim not found\n");
    return;
  }
  
  // Create an array list with the global ratios
  TObjArray globalRatios;
  
  // Create an array with the global plots of the individual efficencies and the ratios
  TObjArray globalRatiosAndEff;
  
//  globalRatios.Add(CreateRatioGraph("RatioEffVsCent","data/sim tracking efficiency versus centrality",effVScentData,effVScentSim));
  
  //---- Eff vs y
  TGraphAsymmErrors* effVSyDataCopy = static_cast<TGraphAsymmErrors*>(effVSyData->Clone()); // We make clones to do not modify them
  TGraphAsymmErrors* effVSySimCopy = static_cast<TGraphAsymmErrors*>(effVSySim->Clone());
  
  TGraphAsymmErrors *ratioY = CreateRatioGraph("RatioEffVsY","data/sim tracking efficiency versus rapidity",*effVSyData,*effVSySim);
  globalRatios.Add(ratioY);
  
  TGraphAsymmErrors* ratioYCopy = static_cast<TGraphAsymmErrors*>(ratioY->Clone());
 
  globalRatiosAndEff.Add(DrawRatio("RatioEffVSyAndEff","Comparison Data&MC tracking efficiency versus rapidity", effVSyDataCopy,effVSySimCopy,ratioYCopy));
  //-----

  //-----Eff vs Pt
  TGraphAsymmErrors* effVSptDataCopy = static_cast<TGraphAsymmErrors*>(effVSptData->Clone()); // We make clones to do not modify them
  TGraphAsymmErrors* effVSptSimCopy = static_cast<TGraphAsymmErrors*>(effVSptSim->Clone());

  TGraphAsymmErrors *ratioPt = CreateRatioGraph("RatioEffVsPt","data/sim tracking efficiency versus Pt",*effVSptData,*effVSptSim);
  globalRatios.Add(ratioPt);
  
  TGraphAsymmErrors* ratioPtCopy = static_cast<TGraphAsymmErrors*>(ratioPt->Clone());
  
   globalRatiosAndEff.Add(DrawRatio("RatioEffVSptAndEff","Comparison Data&MC tracking efficiency versus Pt",effVSptDataCopy,effVSptSimCopy,ratioPtCopy));
  //-----
  
  //----Eff vs phi
  TGraphAsymmErrors* effVSphiDataCopy = static_cast<TGraphAsymmErrors*>(effVSphiData->Clone()); // We make clones to do not modify them
  TGraphAsymmErrors* effVSphiSimCopy = static_cast<TGraphAsymmErrors*>(effVSphiSim->Clone());
  
  TGraphAsymmErrors *ratioPhi = CreateRatioGraph("RatioEffVsPhi","data/sim tracking efficiency versus phi",*effVSphiData,*effVSphiSim);
  globalRatios.Add(ratioPhi);
  
  TGraphAsymmErrors* ratioPhiCopy = static_cast<TGraphAsymmErrors*>(ratioPhi->Clone());
  
  globalRatiosAndEff.Add(DrawRatio("RatioEffVSphiAndEff","Comparison Data&MC tracking efficiency versus phi",effVSphiDataCopy,effVSphiSimCopy,ratioPhiCopy));
  //-------
  
  //----Eff vs y vs phi

  TH2F *effVSyVSphiData = static_cast<TH2F*>(fileData->FindObjectAny("trackingEffVsphi-y"));
  if (!effVSyVSphiData) {
    printf("Efficiency vs rapidity vs phi from data not found\n");
    return;
  }
  TH2F *effVSyVSphiSim = static_cast<TH2F*>(fileSim->FindObjectAny("trackingEffVsphi-y"));
  if (!effVSyVSphiSim) {
    printf("Efficiency vs rapidity vs phi from sim not found\n");
    return;
  }
  Int_t nBins2dX = effVSyVSphiData->GetXaxis()->GetNbins();
  Int_t nBins2dY = effVSyVSphiData->GetYaxis()->GetNbins();
  Double_t effData2D,effSim2D,ratio2D;
  
  TH2F *effVSphiVSyRatio = new TH2F("RatioEffVSphiVSy","EffData/EffSim vs phi vs y",nBins2dX, effVSyVSphiData->GetXaxis()->GetBinLowEdge(1), effVSyVSphiData->GetXaxis()->GetBinUpEdge(nBins2dX),nBins2dY, effVSyVSphiData->GetYaxis()->GetBinLowEdge(1), effVSyVSphiData->GetYaxis()->GetBinUpEdge(nBins2dY));
  effVSphiVSyRatio->GetXaxis()->SetTitle("phi");
  effVSphiVSyRatio->GetYaxis()->SetTitle("y");
    
  for (Int_t i = 1 ; i <= nBins2dX ; i++ )
  {
    for (Int_t j = 1 ; j <= nBins2dY ; j++ )
    {
      effData2D = effVSyVSphiData->GetBinContent(i,j);
      effSim2D = effVSyVSphiSim->GetBinContent(i,j);

      if (effData2D > 0. && effSim2D > 0.)
      {
        ratio2D = effData2D/effSim2D;
//        ratio2DErrh = rat*TMath::Sqrt(effDErrh*effDErrh/effD*effD + effSErrl*effSErrl/effS*effS);
//        ratio2DErrl = rat*TMath::Sqrt(effDErrl*effDErrl/effD*effD + effSErrh*effSErrh/effS*effS);
      }
      if (effData2D == 0 && effSim2D == 0)
      {
        ratio2D = 1.;
//        ratio2DErrh = 0.;
//        ratio2DErrl = 0.;
      }
      if (effData2D == 0 && effSim2D > 0.)
      {
        ratio2D = 0.;
//        ratio2DErrh = 0.;
//        ratio2DErrl = 0.;
      }
      if (effData2D > 0. && effSim2D == 0)
      {
        ratio2D = 2.;
//        ratio2DErrh = 0.;
//        ratio2DErrl = 0.;
      }
      effVSphiVSyRatio->SetBinContent(i,j,ratio2D);
    }
  }
  
  
  TH2F *effVSphiVSyRatioRapBins = new TH2F();
  effVSphiVSyRatioRapBins->GetXaxis()->SetTitle("phi");
  effVSphiVSyRatioRapBins->GetYaxis()->SetTitle("y");
  effVSphiVSyRatioRapBins->SetName("RatioEffVSphiVSyRapBins");
  effVSphiVSyRatioRapBins->SetTitle("EffData/EffSim vs phi vs y");

  
  Int_t nxBins = effVSphiVSyRatio->GetXaxis()->GetNbins();
  Int_t nyBins = effVSphiVSyRatio->GetYaxis()->GetNbins();
  
  Double_t xBinEdge[nxBins+1];
  Double_t yBinEdge[nyBins+1];
  
  for (Int_t ybin = 0 ; ybin <= nyBins ; ybin++)
  {
    yBinEdge[ybin] = 2*TMath::ATan(TMath::Exp((effVSphiVSyRatio->GetYaxis()->GetBinLowEdge(ybin+1))));
  }
  for (Int_t xbin = 0 ; xbin <= nxBins ; xbin++)
  {
    xBinEdge[xbin] = effVSphiVSyRatio->GetXaxis()->GetBinLowEdge(xbin+1);
  }
  
  effVSphiVSyRatioRapBins->SetBins(nxBins,xBinEdge,nyBins,yBinEdge);
  
  for (Int_t xbin = 1 ; xbin <= nxBins ; xbin++)
  {
    for (Int_t ybin = 1 ; ybin <= nyBins ; ybin++)
    {
      effVSphiVSyRatioRapBins->SetBinContent(xbin,ybin,effVSphiVSyRatio->GetBinContent(xbin,ybin));
    }
  }
  globalRatiosAndEff.Add(effVSphiVSyRatio);
  globalRatiosAndEff.Add(effVSphiVSyRatioRapBins);

 
  //--------
  
  // Get Chamber and DE Data and Sim graphs
//  TObjArray *listChEffVSrunData = static_cast<TObjArray*>(fileData->FindObjectAny("ChambersEffVSrun"));
//  if (!listChEffVSrunData) {
//    printf("list of Chamber efficiencies vs run from data not found\n");
//    return;
//  }
//  TObjArray *listChEffVSrunSim = static_cast<TObjArray*>(fileSim->FindObjectAny("ChambersEffVSrun"));
//  if (!listChEffVSrunSim) {
//    printf("list of Chamber efficiencies vs run from sim not found\n");
//    return;
//  }

  TObjArray *listChEffVSDEData = static_cast<TObjArray*>(fileData->FindObjectAny("ChambersEffVSDE"));
  if (!listChEffVSDEData) {
    printf("list of Chamber efficiencies per DE from data not found\n");
    return;
  }
  TObjArray *listChEffVSDESim = static_cast<TObjArray*>(fileSim->FindObjectAny("ChambersEffVSDE"));
  if (!listChEffVSDESim) {
    printf("list of Chamber efficiencies per DE from sim not found\n");
    return;
  }
  
//  TObjArray *listDEEffVSrunData = static_cast<TObjArray*>(fileData->FindObjectAny("DEEffVSrun"));
//  if (!listDEEffVSrunData) {
//    printf("list of DE efficiencies vs run from data not found\n");
//    return;
//  }
//  TObjArray *listDEEffVSrunSim = static_cast<TObjArray*>(fileSim->FindObjectAny("DEEffVSrun"));
//  if (!listDEEffVSrunSim) {
//    printf("list of DE efficiencies vs run from sim not found\n");
//    return;
//  }
  
  // Graph for global efficiency vs run
  TGraphAsymmErrors* gData ;//= static_cast<TGraphAsymmErrors*>(listChEffVSrunData->At(0));
  TGraphAsymmErrors* gSim ;//= static_cast<TGraphAsymmErrors*>(listChEffVSrunSim->At(0));
  
  
  //----Eff vs run  
//  TGraphAsymmErrors *ratioEffvsrRun = CreateRatioGraph("RatioEffVsRun","data/sim tracking efficiency versus run",*gData,*gSim);
//  globalRatios.Add(ratioEffvsrRun);
//  
//  TGraphAsymmErrors* ratioEffvsrRunCopy = static_cast<TGraphAsymmErrors*>(ratioEffvsrRun->Clone());
//  
//  globalRatiosAndEff.Add(DrawRatio("RatioEffVsRunAndEff","Comparison Data&MC tracking efficiency versus run",gData,gSim,ratioEffvsrRunCopy));
  //-------

  //globalRatios.Add(CreateRatioGraph("RatioEffVsRun","data/sim tracking efficiency versus run",*gData,*gSim));
  
  // Create a list with the Chamber and DE ratios
//  TObjArray chamberVSrunRatios;
//  TObjArray deVSrunRatios;
  TObjArray chamberVSdeRatios;
  
//  TObjArray chamberVSrunRatiosAndEff;
//  TObjArray deVSrunRatiosAndEff;
  TObjArray chamberVSdeRatiosAndEff;
  
  // Compute the ratios for Chamber vs run
//  for (Int_t nList = 1 ; nList < listChEffVSrunData->GetEntries() ; nList++)
//  {
//    gData = static_cast<TGraphAsymmErrors*>(listChEffVSrunData->At(nList));
//    gSim = static_cast<TGraphAsymmErrors*>(listChEffVSrunSim->At(nList));
//    if (!gData || !gSim )
//    {
//      printf("Error readig from Chamber efficiency vs run list \n");
//      return;
//    }
//    //----Eff of Chs vs run
//    TString name =  Form("RatioEffCh%dVsRun",nList); TString title = Form("Chamber %d data/sim tracking efficiency versus run",nList);
//    
//    TGraphAsymmErrors *ratioEffChVsrRun = CreateRatioGraph(name.Data(),title.Data(),*gData,*gSim);
//    chamberVSrunRatios.Add(ratioEffChVsrRun);
//    
//    TGraphAsymmErrors* ratioEffChVsrRunCopy = static_cast<TGraphAsymmErrors*>(ratioEffChVsrRun->Clone());
//    
//    TString nameRatio =  Form("RatioEffCh%dVsRunAndEff",nList); TString titleRatio = Form("Comparison Data&MC Ch%d tracking efficiency versus run",nList);
//    chamberVSrunRatiosAndEff.Add(DrawRatio(nameRatio.Data(),titleRatio.Data(),gData,gSim,ratioEffChVsrRunCopy));
//    //-------
//
//    
////    chamberVSrunRatios.Add(CreateRatioGraph(,,*gData,*gSim));
//    
//  }
  
  //Load the mapping for the DE histos
  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  AliCDBManager::Instance()->SetRun(0);
  AliMUONCDB::LoadMapping();
  AliMpDEIterator deit;

  // Loop over Chambers
  for (Int_t ich = 0 ; ich < 10 ; ich++)
  {
    // Compute the ratios for DE vs run
    deit.First(ich);
  
//    while ( !deit.IsDone() )
//    {
//      TString currentDEName = Form("EffDE%dVSrun",deit.CurrentDEId());
//      gData = static_cast<TGraphAsymmErrors*>(listDEEffVSrunData->FindObject(currentDEName.Data()));
//      gSim = static_cast<TGraphAsymmErrors*>(listDEEffVSrunSim->FindObject(currentDEName.Data()));
//      
//      TString name =  Form("RatioEffDE%dVsRun",deit.CurrentDEId()); TString title = Form("DE %d data/sim tracking efficiency versus run",deit.CurrentDEId());
//      if (!gData || !gSim )
//      {
//        printf("Error readig from DE efficiency vs run list \n");
//        return;
//      }
//      //----Eff of DEs vs run
//      TGraphAsymmErrors *ratioEffDEvsRun = CreateRatioGraph(name.Data(),title.Data(),*gData,*gSim);
//      deVSrunRatios.Add(ratioEffDEvsRun);
//      
//      TGraphAsymmErrors* ratioEffDEvsRunCopy = static_cast<TGraphAsymmErrors*>(ratioEffDEvsRun->Clone());
//      
//      TString nameRatio =  Form("RatioEffDE%dVsRunAndEff",deit.CurrentDEId()); TString titleRatio = Form("Comparison Data&MC DE%d tracking efficiency versus run",deit.CurrentDEId());
//      deVSrunRatiosAndEff.Add(DrawRatio(nameRatio.Data(),titleRatio.Data(),gData,gSim,ratioEffDEvsRunCopy));
//      //-------
//
////      deVSrunRatios.Add(CreateRatioGraph(name.Data(),title.Data(),*gData,*gSim));
//      
//      deit.Next();
//    }
  
    // Compute the ratios for Ch vs DE
    gData = static_cast<TGraphAsymmErrors*>(listChEffVSDEData->FindObject(Form("Ch%dEffVSDE",ich+1)));
    gSim = static_cast<TGraphAsymmErrors*>(listChEffVSDESim->FindObject(Form("Ch%dEffVSDE",ich+1)));
    
    if (!gData || !gSim )
    {
      printf("Error reading from Chamber efficiency per DE list \n");
      return;
    }
    TString name =  Form("RatioEffCh%dVsDE",ich+1); TString title = Form("Chamber %d data/sim tracking efficiency versus DE",ich+1);
    //----Eff of CHs vs DE
    TGraphAsymmErrors *ratioEffChvsDE = CreateRatioGraph(name.Data(),title.Data(),*gData,*gSim);
    chamberVSdeRatios.Add(ratioEffChvsDE);
    
    TGraphAsymmErrors* ratioEffChvsDECopy = static_cast<TGraphAsymmErrors*>(ratioEffChvsDE->Clone());
    
    TString nameRatio =  Form("RatioEffCh%dVsDEAndEff",ich+1); TString titleRatio = Form("Comparison Data&MC Ch%d tracking efficiency versus DE",ich+1);
    chamberVSdeRatiosAndEff.Add(DrawRatio(nameRatio.Data(),titleRatio.Data(),gData,gSim,ratioEffChvsDECopy));
    //-------

    
//    chamberVSdeRatios.Add(CreateRatioGraph(name.Data(),title.Data(),*gData,*gSim));
    
  }

  //Beautify graphs
  BeautifyGraphs(globalRatios,"EffData/EffSim");
//  BeautifyGraphs(deVSrunRatios,"EffData/EffSim");
//  BeautifyGraphs(chamberVSrunRatios,"EffData/EffSim");
  BeautifyGraphs(chamberVSdeRatios,"EffData/EffSim");

//  BeautifyGraphs(globalRatiosAndEff,"EffData/EffSim");
//  BeautifyGraphs(deVSrunRatiosAndEff,"EffData/EffSim");
//  BeautifyGraphs(chamberVSrunRatiosAndEff,"EffData/EffSim");
//  BeautifyGraphs(chamberVSdeRatiosAndEff,"EffData/EffSim");

  // set bin labels
//  SetRunLabel(deVSrunRatios,irun,runs);
//  SetRunLabel(chamberVSrunRatios,irun,runs);
//  SetRunLabel(globalRatios,irun,runs,1); //Write it in such a way the number is the position on the list of the graph you want to label
//  
  // save output
  TFile* file = new TFile("EffComparison.root","update");
  
  globalRatios.Write("GlobalEffRatios", TObject::kOverwrite | TObject::kSingleKey);
//  chamberVSrunRatios.Write("ChambersEffVSrunRatios", TObject::kOverwrite | TObject::kSingleKey);
//  deVSrunRatios.Write("DEEffVSrunRatios", TObject::kOverwrite | TObject::kSingleKey);
  chamberVSdeRatios.Write("ChamberEffperDERatios", TObject::kOverwrite | TObject::kSingleKey);
  
  globalRatiosAndEff.Write("GlobalEffRatiosAndEffs", TObject::kOverwrite | TObject::kSingleKey);
//  chamberVSrunRatiosAndEff.Write("ChambersEffVSrunRatiosAndEff", TObject::kOverwrite | TObject::kSingleKey);
//  deVSrunRatiosAndEff.Write("DEEffVSrunRatiosAndEff", TObject::kOverwrite | TObject::kSingleKey);
  chamberVSdeRatiosAndEff.Write("ChamberEffperDERatiosAndEff", TObject::kOverwrite | TObject::kSingleKey);

  
   
  file->Close();
 
  fileData->Close();
  fileSim->Close();
   


}
