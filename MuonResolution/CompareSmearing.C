
/*
 *  CompareSmearing.C
 *
 *  Created by Philippe Pillot on 13/12/17.
 *  Copyright 2017 SUBATECH
 *  This software is made available under the terms of the GNU GPL 3.0
 *
 */

#if !defined(__CINT__) || defined(__MAKECINT__)

#include <stdarg.h>
#include <TROOT.h>
#include <TString.h>
#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TObjArray.h>

#endif

//-----------------------------------------------------------------------
void CompareSmearing(char *fileRef ...)
{
  /// compare the reconstructed distributions after different smearing
  /// parameters format is "fileName[:Legend]"
  
  TH1F *hpTRec0 = 0x0, *hetaRec0 = 0x0;
  Int_t color = 1;
  
  TCanvas *cRec = new TCanvas("cRec", "cRec");
  cRec->Divide(2, 2);
  TLegend *lRec = new TLegend(0.60,0.50,0.85,0.7);
  
  va_list ap;
  va_start(ap, fileRef);
  
  for (Int_t i = 0; ; ++i) {
    
    TString fileNameLegend = (i == 0) ? fileRef : va_arg(ap,char*);
    if (fileNameLegend.IsNull()) break;
    
    TObjArray *oNameLegend = fileNameLegend.Tokenize(":");
    TString fileName = oNameLegend->At(0)->GetName();
    TString fileLegend;
    if (oNameLegend->GetEntries() > 1) fileLegend = oNameLegend->At(1)->GetName();
    else fileLegend = TString::Format("file%d",i);
    
    TFile *file = TFile::Open(fileName.Data(),"READ");
    if (!file || !file->IsOpen()) return;
    
    TObjArray* recHistoList = static_cast<TObjArray*>(file->FindObjectAny("RecHistos"));
    if (!recHistoList) return;
    
    TH1F *hpTRec = static_cast<TH1F*>(recHistoList->FindObject("hpTRec")->Clone());
    TH1F *hetaRec = static_cast<TH1F*>(recHistoList->FindObject("hetaRec")->Clone());
    if (!hpTRec || !hetaRec) return;
    hpTRec->SetDirectory(0);
    hetaRec->SetDirectory(0);
    
    file->Close();
    
    lRec->AddEntry(hpTRec,fileLegend.Data(),"PL");
    
    gROOT->SetSelectedPad(cRec->cd(1));
    gPad->SetLogy();
    hpTRec->SetLineColor(color);
//    hpTRec->Rebin(10);
    hpTRec->Rebin(4);
    hpTRec->Draw((gPad->GetListOfPrimitives()->GetEntries() > 0) ? "sames" : "");
    gROOT->SetSelectedPad(cRec->cd(2));
    gPad->SetLogy();
    hetaRec->SetLineColor(color);
    hetaRec->Draw((gPad->GetListOfPrimitives()->GetEntries() > 0) ? "sames" : "");
//    hetaRec->Rebin(23);
    hetaRec->Rebin(5);

    if (i == 0) {
      
      hpTRec0 = hpTRec;
      hetaRec0 = hetaRec;
      
    } else {
      
      TH1F *hpTRat = static_cast<TH1F*>(hpTRec->Clone());
      hpTRat->SetNameTitle("hpTRat", "hpTRat");
      hpTRat->Divide(hpTRec0);
      hpTRat->SetStats(kFALSE);
      TH1F *hetaRat = static_cast<TH1F*>(hetaRec->Clone());
      hetaRat->SetNameTitle("hetaRat", "hetaRat");
      hetaRat->Divide(hetaRec0);
      hetaRat->SetStats(kFALSE);

      gROOT->SetSelectedPad(cRec->cd(3));
      hpTRat->SetLineColor(color);
      hpTRat->Draw((gPad->GetListOfPrimitives()->GetEntries() > 0) ? "same" : "");
      gROOT->SetSelectedPad(cRec->cd(4));
      hetaRat->SetLineColor(color);
      hetaRat->Draw((gPad->GetListOfPrimitives()->GetEntries() > 0) ? "same" : "");
      
    }
    
    ++color;
    if (color == 5 || color == 10) ++color;
    
  }
  
  va_end(ap);
  
  gROOT->SetSelectedPad(cRec->cd(1));
  lRec->Draw();
  
}

