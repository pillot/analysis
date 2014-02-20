/*
 *  CompareResolution.C
 *  aliroot
 *
 *  Created by Philippe Pillot on 15/10/10.
 *  Copyright 2010 SUBATECH. All rights reserved.
 *
 */

#if !defined(__CINT__) || defined(__MAKECINT__)

#include <Riostream.h>
#include "TStyle.h"
#include "TFile.h"
#include "TList.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TAxis.h"

#endif

Bool_t drawLegend = kTRUE;

void CompareResolution(const char* fileName1, const char* fileName2,
		       const char* fileName3 = 0x0, const char* fileName4 = 0x0,
		       Bool_t print = kFALSE)
{
  /// Compare chamber and DE systematic shift and resolution between the 2 results
  /// fileName1 labelled "before"
  /// fileName2 labelled "after"
  
  gStyle->SetFillColor(0);
  
  Int_t nFiles = 2;
  if (fileName4) nFiles = 4;
  else if (fileName3) nFiles = 3;
  TString fileName[4] = {fileName1, fileName2, fileName3, fileName4};
  Int_t color[4] = {2, 4, 6, 8};
  if (!fileName4) color[2] = color[3];
  
  const Int_t nObjs = 4;
  TString objName[3][nObjs] = {{"gResidualXPerChMean_ClusterIn", "gCombinedResidualXPerChSigma", "gResidualYPerChMean_ClusterIn", "gCombinedResidualYPerChSigma"},
			       {"gResidualXPerHalfChMean_ClusterIn", "gCombinedResidualXPerHalfChSigma", "gResidualYPerHalfChMean_ClusterIn", "gCombinedResidualYPerHalfChSigma"},
			       {"gResidualXPerDEMean_ClusterIn", "gResidualYPerDEMean_ClusterIn", "gCombinedResidualXPerDESigma", "gCombinedResidualYPerDESigma"}};
  TString mgName[6] = {"mgResidualXMeanVsP_ClusterIn", "mgResidualYMeanVsP_ClusterIn",
		       "mgResidualXMeanVsAngle_ClusterIn", "mgResidualYMeanVsAngle_ClusterIn",
		       "mgHChResidualXMeanVsAngle_ClusterIn", "mgHChResidualYMeanVsAngle_ClusterIn"};
  
  // loop over files and get objects
  TFile **outFile = new TFile*[nFiles];
  TGraphErrors **gObj[3][nObjs];
  for (Int_t i=0; i<3; i++) for (Int_t j=0; j<nObjs; j++) gObj[i][j] = new TGraphErrors*[nFiles];
  TMultiGraph **mgObj[6] = {0x0, 0x0, 0x0, 0x0, 0x0, 0x0};
  for (Int_t i=0; i<6; i++) mgObj[i] = new TMultiGraph*[nFiles];
  for (Int_t iFile=0; iFile<nFiles; iFile++) {
    // open file
    outFile[iFile] = TFile::Open(fileName[iFile],"READ");
    if (!outFile[iFile] || !outFile[iFile]->IsOpen()) return;
    // get results
    TObjArray* list = static_cast<TObjArray*>(outFile[iFile]->FindObjectAny("ChamberRes"));
    if (!list) return;
    for (Int_t i=0; i<3; i++) for (Int_t j=0; j<nObjs; j++) {
      gObj[i][j][iFile] = static_cast<TGraphErrors*>(list->FindObject(objName[i][j].Data()));
      if (!gObj[i][j][iFile]) continue;
    }
    for (Int_t i=0; i<6; i++) {
      mgObj[i][iFile] = static_cast<TMultiGraph*>(list->FindObject(mgName[i].Data()));
      if (!mgObj[i][iFile]) continue;
    }
  }
  
  // print
  if (print) {
    
    printf("\n");
    Double_t dummy, val;
    for (Int_t iFile=0; iFile<nFiles; iFile++) {
      
      printf("------ in file %d ------\n",iFile);
      printf("\nhalf-chamber systematic shifts:\n");
      printf(" - non-bending:");
      for (Int_t i = 0; i < 20; i++) {
	if (!gObj[1][0][iFile]) continue;
	gObj[1][0][iFile]->GetPoint(i, dummy, val);
	printf((i==0)?" %5.3f":", %5.3f",val);
      }
      printf("\n -     bending:");
      for (Int_t i = 0; i < 20; i++) {
	if (!gObj[1][2][iFile]) continue;
	gObj[1][2][iFile]->GetPoint(i, dummy, val);
	printf((i==0)?" %6.4f":", %6.4f",val);
      }
      printf("\n\nhalf-chamber resolution:\n");
      printf(" - non-bending:");
      for (Int_t i = 0; i < 20; i++) {
	if (!gObj[1][1][iFile]) continue;
	gObj[1][1][iFile]->GetPoint(i, dummy, val);
	printf((i==0)?" %5.3f":", %5.3f",val);
      }
      printf("\n -     bending:");
      for (Int_t i = 0; i < 20; i++) {
	if (!gObj[1][3][iFile]) continue;
	gObj[1][3][iFile]->GetPoint(i, dummy, val);
	printf((i==0)?" %6.4f":", %6.4f",val);
      }
      printf("\n\n");
      
    }
    
  }
  
  // display
//  TString legend[3] = {"before", "after", ""};
//  TString legend[3] = {"old", "new", ""};
//  TString legend[3] = {"no match", "match", ""};
//  TString legend[3] = {"ref", "4cf", "5cf"};
//  TString legend[3] = {"ref", "new_0", "new_1"};
//  TString legend[3] = {"w/o GMS", "w/ GMS", ""};
//  TString legend[3] = {"Align BOFF", "Align BON", ""};
//  TString legend[3] = {"p-p", "Pb-Pb", ""};
//  TString legend[3] = {"Vanik", "Javier", ""};
//  TString legend[3] = {"all st", "st345", ""};
//  TString legend[3] = {"all", "mu+", "mu-"};
//  TString legend[3] = {"BOFF", "BON", "GMS"};
//  TString legend[3] = {"Javier+GMS", "Ruben+GMS", "Ruben+GMS"};
//  TString legend[3] = {"pp GMSpp", "PbPb GMSpp", "PbPb GMSPbPb"};
//  TString legend[4] = {"20GeV woMono", "20 GeV woMono2", "woMono", "woMono2"};
//  TString legend[4] = {"PbPb woMono", "PbPb woMono improve", "pp woMono", "pp woMono improve"};
//  TString legend[2] = {"PbPb GMSRuben", "PbPb GMSJavier2"};
//  TString legend[4] = {"pp GMSRuben", "PbPb GMSRuben", "PbPb GMSJavier", "PbPb GMSJavier2"};
//  TString legend[3] = {"Javier1", "Javier2", "Javier"};
//  TString legend[2] = {"w/ mono-cathodes", "w/o mono-cathodes"};
  TString legend[3] = {"no cut", "> 20 GeV/c"};
  TLegend *lRes = new TLegend(0.70,0.80,0.99,0.99);
  for (Int_t i=0; i<nFiles; i++) {
    if (!gObj[0][0][i]) continue;
    lRes->AddEntry(gObj[0][0][i],legend[i].Data(),"PL");
  }
  
  // comparison
  TString cName[3] = {"cResPerCh", "cResPerHalfCh", "cResPerDE"};
  Int_t size[2][3] = {{600, 1200, 1200}, {500, 500, 800}};
  Int_t divide[2][3] = {{2, 2, 1}, {2, 2, 4}};
  Float_t titleSize[2][3] = {{0.05, 0.05, 0.05}, {0.05, 0.05, 0.08}};
  Float_t labelSize[2][3] = {{0.06, 0.07, 0.06}, {0.06, 0.06, 0.08}};
  Float_t titleOffset[2][3] = {{1., 1., 1.}, {1., 1., 0.4}};
  Float_t range[3][nObjs][2] = {{{-0.015, 0.015}, {0.02, 0.1}, {-0.015, 0.015}, {0., 0.1}},
				{{-0.02, 0.02}, {0.02, 0.1}, {-0.015, 0.015}, {0., 0.05}},
				{{-0.1, 0.1}, {-0.1, 0.1}, {0., 0.5}, {0., 0.7}}};
  for (Int_t i=0; i<3; i++) {
    TCanvas* c = new TCanvas(cName[i].Data(), cName[i].Data(), size[0][i], size[1][i]);
    c->Divide(divide[0][i], divide[1][i]);
    for (Int_t j=0; j<nObjs; j++) {
      c->cd(j+1);
      for (Int_t iFile=0; iFile<nFiles; iFile++) {
	if (!gObj[i][j][iFile]) continue;
	gObj[i][j][iFile]->SetMarkerColor(color[iFile]);
	gObj[i][j][iFile]->SetLineColor(color[iFile]);
	if (iFile == 0) {
	  gObj[i][j][iFile]->GetXaxis()->SetTitleSize(titleSize[0][i]);
	  gObj[i][j][iFile]->GetXaxis()->SetLabelSize(labelSize[0][i]);
	  gObj[i][j][iFile]->GetXaxis()->SetTitleOffset(titleOffset[0][i]);
	  gObj[i][j][iFile]->GetYaxis()->SetTitleSize(titleSize[1][i]);
	  gObj[i][j][iFile]->GetYaxis()->SetLabelSize(labelSize[1][i]);
	  gObj[i][j][iFile]->GetYaxis()->SetTitleOffset(titleOffset[1][i]);
	  gObj[i][j][iFile]->GetYaxis()->SetRangeUser(range[i][j][0],range[i][j][1]);
	  gObj[i][j][iFile]->DrawClone("ap");
	} else gObj[i][j][iFile]->DrawClone("p");
	if (drawLegend) lRes->DrawClone();
      }
    }
  }
  
  // comparison versus P and angle per chamber
  TString gName[4] = {"gShiftXVsP_ch", "gShiftYVsP_ch", "gShiftXVsAngle_ch", "gShiftYVsAngle_ch"};
  TString cName2[2] = {"cShiftVsP", "cShiftVsAngle"};
  Int_t size2[2][2] = {{1200, 1200}, {900, 900}};
  Int_t divide2[2][2] = {{5, 5}, {4, 4}};
  Float_t titleSize2[2][2] = {{0.05, 0.05}, {0.05, 0.05}};
  Float_t labelSize2[2][2] = {{0.06, 0.06}, {0.06, 0.06}};
//  Float_t range2[2][2] = {{-0.05, -0.1}, {0.05, 0.1}};
  Float_t range2[2][2] = {{-0.025, -0.05}, {0.025, 0.05}};
  for (Int_t i=0; i<2; i++) {
    TCanvas* c = new TCanvas(cName2[i].Data(), cName2[i].Data(), size2[0][i], size2[1][i]);
    c->Divide(divide2[0][i], divide2[1][i]);
    for (Int_t j=0; j<2; j++) {
      for (Int_t k=0; k<10; k++) {
	c->cd(10*j+k+1);
	for (Int_t iFile=0; iFile<nFiles; iFile++) {
	  if (!mgObj[2*i+j][iFile]) continue;
	  TGraphErrors *g = static_cast<TGraphErrors*>(mgObj[2*i+j][iFile]->GetListOfGraphs()->FindObject(Form("%s%d",gName[2*i+j].Data(),k+1)));
	  g->SetTitle(g->GetName());
	  g->SetMarkerColor(color[iFile]);
	  g->SetLineColor(color[iFile]);
	  if (iFile == 0) {
	    g->GetXaxis()->SetTitleSize(titleSize2[0][i]);
	    g->GetXaxis()->SetLabelSize(labelSize2[0][i]);
	    g->GetYaxis()->SetTitleSize(titleSize2[1][i]);
	    g->GetYaxis()->SetLabelSize(labelSize2[1][i]);
	    g->GetYaxis()->SetRangeUser(range2[0][i],range2[1][i]);
	    g->DrawClone("ap");
	  } else g->DrawClone("p");
	  if (drawLegend) lRes->DrawClone();
	}
      }
    }
  }
  
  // comparison versus angle per half-chamber
  TString gName3[2] = {"gShiftXVsAngle_halfCh", "gShiftYVsAngle_halfCh"};
  TString cName3[2] = {"cHChShiftXVsAngle", "cHChShiftYVsAngle"};
  TString side[2] = {"I", "O"};
  Int_t size3[2][2] = {{1200, 1200}, {900, 900}};
  Int_t divide3[2][2] = {{5, 5}, {4, 4}};
  Float_t titleSize3[2][2] = {{0.05, 0.05}, {0.05, 0.05}};
  Float_t labelSize3[2][2] = {{0.06, 0.06}, {0.06, 0.06}};
//  Float_t range3[2][2] = {{-0.1, -0.1}, {0.1, 0.1}};
  Float_t range3[2][2] = {{-0.05, -0.05}, {0.05, 0.05}};
  for (Int_t i=0; i<2; i++) {
    TCanvas* c = new TCanvas(cName3[i].Data(), cName3[i].Data(), size3[0][i], size3[1][i]);
    c->Divide(divide3[0][i], divide3[1][i]);
    for (Int_t j=0; j<10; j++) {
      for (Int_t k=0; k<2; k++) {
	c->cd(2*j+k+1);
	for (Int_t iFile=0; iFile<nFiles; iFile++) {
	  if (!mgObj[4+i][iFile]) continue;
	  TGraphErrors *g = static_cast<TGraphErrors*>(mgObj[4+i][iFile]->GetListOfGraphs()->FindObject(Form("%s%d%s",gName3[i].Data(),j+1,side[k].Data())));
	  g->SetTitle(g->GetName());
	  g->SetMarkerColor(color[iFile]);
	  g->SetLineColor(color[iFile]);
	  if (iFile == 0) {
	    g->GetXaxis()->SetTitleSize(titleSize3[0][i]);
	    g->GetXaxis()->SetLabelSize(labelSize3[0][i]);
	    g->GetYaxis()->SetTitleSize(titleSize3[1][i]);
	    g->GetYaxis()->SetLabelSize(labelSize3[1][i]);
	    g->GetYaxis()->SetRangeUser(range3[0][i],range3[1][i]);
	    g->DrawClone("ap");
	  } else g->DrawClone("p");
	  if (drawLegend) lRes->DrawClone();
	}
      }
    }
  }
  
  // close file and clean memory
  delete lRes;
  for (Int_t i=0; i<3; i++) for (Int_t j=0; j<nObjs; j++) delete[] gObj[i][j];
  for (Int_t i=0; i<6; i++) delete[] mgObj[i];
  for (Int_t iFile=0; iFile<nFiles; iFile++) outFile[iFile]->Close();
  delete[] outFile;
  
}
