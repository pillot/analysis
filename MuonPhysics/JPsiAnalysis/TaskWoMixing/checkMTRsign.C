#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include "TSystem.h"
#include "TROOT.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1.h"

#include "AliCFGridSparse.h"
#include "AliAnalysisMuonUtility.h"
#include "AliMergeableCollection.h"

#include "AliCDBManager.h"
#include "AliMUONTriggerDisplay.h"

#endif

void checkMTRsign ( Int_t trackSign, Int_t trigSign,
                   TString identifier = "/PhysSelPass/CPBI1MUL-B-NOPF-MUON,CPBI1MUL-B-NOPF-MUON&CPBI1MLL-B-NOPF-MUON,CPBI1MLL-B-NOPF-MUON/0_10,10_20,20_30,30_40,40_50,50_60,60_70,70_80,80_90",
//                   TString identifier = "/PhysSelPass/CPBI1MUL-B-NOPF-MUON,CPBI1MUL-B-NOPF-MUON&CPBI1MLL-B-NOPF-MUON,CPBI1MLL-B-NOPF-MUON/0_10,10_20",
//                   TString identifier = "/PhysSelPass/CPBI1MUL-B-NOPF-MUON,CPBI1MUL-B-NOPF-MUON&CPBI1MLL-B-NOPF-MUON,CPBI1MLL-B-NOPF-MUON/40_50,50_60,60_70,70_80,80_90",
//                   TString identifier = "/PhysSelPass/CPBI2_B1-B-NOPF-ALLNOTRD/0_10,10_20",
//                   TString identifier = "/PhysSelPass/CPBI2_B1-B-NOPF-ALLNOTRD/40_50,50_60,60_70,70_80,80_90",
                   TString filename="all" )
{
  Int_t varCharge = 2;
  Int_t varTrigSign = 4;
  Int_t varLoCircuit = 5;
  Int_t nVars = 6;
  
  TFile* outfile = TFile::Open("mismatch.root","RECREATE");
  
  AliMergeableCollection* mcol[2];
  if (filename == "all") {
    
    TFile* file = TFile::Open("Bplus.root");
    if ( ! file ) {
      printf("Error: cannot find Bplus.root\n");
      return;
    }
    mcol[0] = static_cast<AliMergeableCollection*>(file->FindObjectAny("MTRSignOut"));
    file = TFile::Open("Bminus.root");
    if ( ! file ) {
      printf("Error: cannot find Bminus.root\n");
      return;
    }
    mcol[1] = static_cast<AliMergeableCollection*>(file->FindObjectAny("MTRSignOut"));
    if ( ! mcol[0] || ! mcol[1] ) {
      printf("Error: cannot find mergeable collection\n");
      file->ls();
      return;
    }
    
  } else {
    
    TFile* file = TFile::Open(filename.Data());
    if ( ! file ) {
      printf("Error: cannot find %s\n", filename.Data());
      return;
    }
    mcol[0] = static_cast<AliMergeableCollection*>(file->FindObjectAny("MTRSignOut"));
    if ( ! mcol[0] ) {
      printf("Error: cannot find mergeable collection\n");
      file->ls();
      return;
    }
    mcol[1] = 0x0;
    
  }
  
  TString fullName = Form("%s/MTRSign", identifier.Data());
  fullName.ReplaceAll("//","/");
  AliCFGridSparse* gridSparse[2];
  gridSparse[0] = static_cast<AliCFGridSparse*>(mcol[0]->GetSum(fullName.Data()));
  gridSparse[1] = mcol[1] ? static_cast<AliCFGridSparse*>(mcol[1]->GetSum(fullName.Data())) : 0x0;
  if ( ! gridSparse[0] || (mcol[1] && !gridSparse[1]) ) {
    printf("Cannot find %s\n",fullName.Data());
    mcol[0]->Print("*");
    return;
  }
  
//  AliAnalysisMuonUtility::SetSparseRange(gridSparse[0],0,"",2.,14.99); // REMEMBER TO CUT
//  if (gridSparse[1]) AliAnalysisMuonUtility::SetSparseRange(gridSparse[1],0,"",2.,14.99); // REMEMBER TO CUT
  //      AliAnalysisMuonUtility::SetSparseRange(gridSparse,1,"",20.,100.); // REMEMBER TO CUT

  Int_t ishift = 0;
  Int_t shiftPx = 50;
  TCanvas* can = 0x0;

  TString canName = "canCorr";
  can = new TCanvas(canName.Data(), canName.Data(), 10+shiftPx*ishift++, 10, 900, 600);
  TH1* histo = gridSparse[0]->Project(varCharge,varTrigSign);
  histo->Sumw2();
  if (gridSparse[1]) histo->Add(gridSparse[1]->Project(varCharge,varTrigSign));
  Double_t nTracks = histo->Integral();
  histo->SetName("trackChargeVsTrigSign");
  histo->Scale(100./nTracks);
  histo->SetStats(kFALSE);
  histo->Draw("TEXTCOLZ");
  outfile->cd();
  histo->Write();
  
  Int_t trackSignBin = ( trackSign < 0 ) ? 1 : 2;
  Int_t optrackSignBin = ( trackSign < 0 ) ? 2 : 1;
  Int_t trigSignBinMin = ( trigSign < 0 ) ? 1 : 3;
  Int_t trigSignBinMax = ( trigSign < 0 ) ? 1 : 3;
  
  TObjArray* histoList[2] = {0x0,0x0};
  TString baseName = "";
  for ( Int_t itype=0; itype<2; itype++ ) {
    histoList[itype] = new TObjArray(nVars);
    if ( itype == 0 ) {
      baseName = "All";
      // Apply cuts
      AliAnalysisMuonUtility::SetSparseRange(gridSparse[0],varCharge,"",trackSignBin,trackSignBin,"USEBIN");
      if (gridSparse[1]) AliAnalysisMuonUtility::SetSparseRange(gridSparse[1],varCharge,"",optrackSignBin,optrackSignBin,"USEBIN");
    }
    else {
      baseName = "BadSign";
      // Apply cuts
      AliAnalysisMuonUtility::SetSparseRange(gridSparse[0],varCharge,"",trackSignBin,trackSignBin,"USEBIN");
      AliAnalysisMuonUtility::SetSparseRange(gridSparse[0],varTrigSign,"",trigSignBinMin,trigSignBinMax,"USEBIN");
      if (gridSparse[1]) {
        AliAnalysisMuonUtility::SetSparseRange(gridSparse[1],varCharge,"",optrackSignBin,optrackSignBin,"USEBIN");
        AliAnalysisMuonUtility::SetSparseRange(gridSparse[1],varTrigSign,"",trigSignBinMin,trigSignBinMax,"USEBIN");
      }
    }
    for ( Int_t ivar = 0; ivar<nVars; ivar++ ) {
      if ( ivar == varCharge || ivar == varTrigSign ) continue;
      histo = gridSparse[0]->Project(ivar);
      histo->Sumw2();
      if (gridSparse[1]) histo->Add(gridSparse[1]->Project(ivar));
      histoList[itype]->AddAtAndExpand(histo,ivar);
    }
  }
  
  
  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  AliCDBManager::Instance()->SetRun(0);
  AliMUONTriggerDisplay trigDispl;
  for ( Int_t ivar = 0; ivar<nVars; ivar++ ) {
    histo = static_cast<TH1*>(histoList[1]->At(ivar));
    if ( ! histo ) continue;
    TH1* auxHisto = static_cast<TH1*>(histoList[0]->At(ivar));
    histo->Divide(auxHisto);
    histo->Scale(100.);
    canName = Form("can_proj_%i",ivar);
    can = new TCanvas(canName.Data(), canName.Data(), 10+shiftPx*ishift++, 10, 900, 600);
    //if ( ivar < 2 ) can->SetLogy();
    histo->SetStats(kFALSE);
    if ( ivar == varLoCircuit ) {
      TH1* displayHisto = trigDispl.GetDisplayHistogram(histo,Form("%s_display",histo->GetName()),AliMUONTriggerDisplay::kDisplayBoards,0,11,"Fraction of tracks per board");
      displayHisto->SetStats(kFALSE);
      displayHisto->Draw("COLZ");
      displayHisto = trigDispl.GetBoardNumberHisto("boardNumbers");
      displayHisto->Draw("TEXTSAME");
      outfile->cd();
      displayHisto->Write();
    }
    else {
      histo->Draw("e");
      outfile->cd();
      histo->Write();
    }
  }
  for ( Int_t itype=0; itype<2; itype++ ) {
    delete histoList[itype];
  }
  
}
