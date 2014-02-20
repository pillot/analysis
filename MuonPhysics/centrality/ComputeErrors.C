/*
 *  ComputeErrors.C
 *  aliroot_dev
 *
 *  Created by philippe pillot on 09/12/11.
 *  Copyright 2011 SUBATECH. All rights reserved.
 *
 */

void ComputeErrors(Bool_t weight = kTRUE)
{
  /// read the files GlauberInfo*.root produced by the other macro to extract the error on weighted NPart
  /// error is computed as sqrt(Sum(Np_x-Np_0)^2) where Np_0 is the reference (sigma64_mind4_r662_a546)
  
  const Int_t nFiles = 13;
  TString fileNames[nFiles] = {
  "GlauberInfo_sigma64_mind4_r662_a546.root",
  "GlauberInfo_sigma69_mind4_r662_a546.root",
  "GlauberInfo_sigma59_mind4_r662_a546.root",
  "GlauberInfo_sigma64_mind0_r662_a546.root",
  "GlauberInfo_sigma64_mind8_r662_a546.root",
  "GlauberInfo_sigma64_mind4_r668_a546.root",
  "GlauberInfo_sigma64_mind4_r656_a546.root",
  "GlauberInfo_sigma64_mind4_r662_a556.root",
  "GlauberInfo_sigma64_mind4_r662_a536.root",
  "GlauberInfo_sigma64_mind4_r668_a556.root",
  "GlauberInfo_sigma64_mind4_r668_a536.root",
  "GlauberInfo_sigma64_mind4_r656_a556.root",
  "GlauberInfo_sigma64_mind4_r656_a536.root"
}
  
  // get the reference values
  TFile *fileRef = TFile::Open(fileNames[0].Data(), "READ");
  TGraphErrors *gNpartRef = fileRef->FindObjectAny("gNpart");
  Int_t nValues = gNpartRef->GetN();
  Double_t *nPartRef = gNpartRef->GetY();
  TGraphErrors *gNcollRef = fileRef->FindObjectAny("gNcoll");
  Double_t *nCollRef = gNcollRef->GetY();
  
  // loop over centrality bins
  Double_t *wNpartRef = new Double_t[nValues];
  Double_t *wNpartErr = new Double_t[nValues];
  Double_t *nCollErr = new Double_t[nValues];
  for (Int_t iCent=0; iCent<nValues; iCent++) {
    
    // weighted Npart
    wNpartRef[iCent] = (weight) ? nPartRef[iCent]/nCollRef[iCent] : nPartRef[iCent];
    
    // error
    wNpartErr[iCent] = 0.;
    nCollErr[iCent] = 0.;
    
  }
  
  // loop over other inputs
  for (Int_t iFile=1; iFile<nFiles; iFile++) {
    
    // get other values
    TFile *file = TFile::Open(fileNames[iFile].Data(), "READ");
    TGraphErrors *gNpart = file->FindObjectAny("gNpart");
    Double_t *nPart = gNpart->GetY();
    TGraphErrors *gNcoll = file->FindObjectAny("gNcoll");
    Double_t *nColl = gNcoll->GetY();
    
    // loop over centrality bins
    for (Int_t iCent=0; iCent<nValues; iCent++) {
      
      // weighted Npart
      Double_t wNpart = (weight) ? nPart[iCent]/nColl[iCent] : nPart[iCent];
      
      // error
      wNpartErr[iCent] += (wNpart - wNpartRef[iCent]) * (wNpart - wNpartRef[iCent]);
      nCollErr[iCent] += (nColl[iCent] - nCollRef[iCent]) * (nColl[iCent] - nCollRef[iCent]);
      
    }
    
    file->Close();
    delete file;
  }
  
  // print results
  for (Int_t iCent=0; iCent<nValues; iCent++) {
    printf("bin %d: wNpart = %f ± %f\t\t Ncoll = %f ± %f\n", iCent,
	   wNpartRef[iCent], TMath::Sqrt(wNpartErr[iCent]),
	   nCollRef[iCent], TMath::Sqrt(nCollErr[iCent]));
  }
  
  fileRef->Close();
  delete fileRef;
  
}

