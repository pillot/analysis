//
//  DrawOccupancy.C
//  aliroot_dev
//
//  Created by philippe pillot on 30/04/2014.
//  Copyright (c) 2014 Philippe Pillot. All rights reserved.
//

#include <stdio.h>

#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>

#include "AliCDBManager.h"

#include "AliMUONCDB.h"
#include "AliMUONVClusterStore.h" // need something from the library rec to find symbol AliMUONCDB at the compilation !?
#include "AliMUONCalibrationData.h"
#include "AliMUONRecoParam.h"
#include "AliMUONDigitCalibrator.h"
#include "AliMUONConstants.h"
#include "AliMUONVDigit.h"
#include "AliMUONVDigitStore.h"
#include "AliMUON2DMap.h"
#include "AliMUONVCalibParam.h"
#include "AliMUONCalibParamND.h"
#include "AliMUONTrackerData.h"

#include "AliMpConstants.h"
#include "AliMpCDB.h"
#include "AliMpDEIterator.h"
#include "AliMpVSegmentation.h"
#include "AliMpSegmentation.h"
#include "AliMpVPadIterator.h"
#include "AliMpPad.h"


const Int_t nDE[10] = {4, 4, 4, 4, 18, 18, 26, 26, 26, 26};

const Double_t maxOcc = 0.2; // maximum occupancy

void GetNPads(Int_t nPadsTot[10][3][3], Int_t nPadsDETot[10][26][3]);
void GetNWorkingPads(Int_t run, Int_t nPadsTot[10][3][3], Int_t nPadsDETot[10][26][3], AliMUON2DMap& occStore);
void CountFiredPads(AliMUONVDigitStore *digitStore, Int_t nPadsCurrent[10][3][3], Int_t nPadsDECurrent[10][26][3],
                    AliMUON2DMap &occStore, Long64_t nEvents);
Bool_t isOnTheRight(Int_t deId);


//------------------------------------------------------------------
void DrawOccupancy(const char *digitFileName, Int_t run, Long64_t event = -1)
{
  /// compute detector occupancy from reconstructed digits

  AliCDBManager* man = AliCDBManager::Instance();
  //  man->SetDefaultStorage("raw://");
  //  man->SetDefaultStorage("local:///Users/pillot/Work/Alice/Data/2015/CDBMirror/alice/data/2015/OCDB");
  //  man->SetDefaultStorage("local:///Users/pillot/Work/Alice/Data/2011/CDBMirror/alice/data/2011/OCDB");
  man->SetDefaultStorage("local:///dev/null");
  man->SetSnapshotMode("OCDB.root");
  man->SetRun(run);
  if (!AliMUONCDB::LoadMapping()) {
    return;
  }

  AliMUON2DMap occStore(kTRUE);

  // get the total number of pads
  Int_t nPadsTot[10][3][3];
  Int_t nPadsDETot[10][26][3];
  GetNPads(nPadsTot, nPadsDETot);

  // get the number of working pads
  Int_t nPadsWorking[10][3][3];
  Int_t nPadsDEWorking[10][26][3];
  GetNWorkingPads(run, nPadsWorking, nPadsDEWorking, occStore);

  Int_t nPadsSt12Working = 0;
  for (Int_t iCh = 0; iCh < 4; ++iCh) {
    nPadsSt12Working += nPadsWorking[iCh][2][2];
  }

  Int_t nPadsSt345Working = 0;
  for (Int_t iCh = 4; iCh < 10; ++iCh) {
    nPadsSt345Working += nPadsWorking[iCh][2][2];
  }

  for (Int_t iCh = 0; iCh < 10; ++iCh) {
    printf("ch%02d: nPadsTot = %d, nWorkingPads = %d\n", iCh + 1, nPadsTot[iCh][2][2], nPadsWorking[iCh][2][2]);
  }

  Int_t nPadsSpectroTot[3][3];
  Int_t nPadsSpectroWorking[3][3];
  for (Int_t iPlane = 0; iPlane < 3; ++iPlane) {
    for (Int_t iSide = 0; iSide < 3; ++iSide) {
      nPadsSpectroTot[iPlane][iSide] = 0;
      nPadsSpectroWorking[iPlane][iSide] = 0;
      for (Int_t iCh = 0; iCh < 10; ++iCh) {
        nPadsSpectroTot[iPlane][iSide] += nPadsTot[iCh][iPlane][iSide];
        nPadsSpectroWorking[iPlane][iSide] += nPadsWorking[iCh][iPlane][iSide];
      }
    }
  }
  printf("spectro: nPadsTot = %d, nWorkingPads = %d\n", nPadsSpectroTot[2][2], nPadsSpectroWorking[2][2]);

  // read digits
  TFile* digitFile = TFile::Open(digitFileName);
  if (!digitFile) return;
  TTree* treeD = static_cast<TTree*>(digitFile->Get("TreeD"));
  if (!treeD) return;
  AliMUONVDigitStore* digitStore = AliMUONVDigitStore::Create(*treeD);
  digitStore->Connect(*treeD);
  
  Int_t nPads[10][3][3];
  TH1F *hOcc[10][3][3];
  TH1F *hAsym[10][3];
  Int_t nPadsDE[10][26][3];
  TH1F *hOccDE[10][26][3];
  TH1F *hDEBad[10][3];
  const char* plane[3] = {"Bending", "NonBending", ""};
  const char* side[3] = {"Left", "Right", ""};
  for (Int_t iCh = 0; iCh < 10; ++iCh) {
    for (Int_t iPlane = 0; iPlane < 3; ++iPlane) {
      for (Int_t iSide = 0; iSide < 3; ++iSide) {
        hOcc[iCh][iPlane][iSide] = new TH1F(Form("hOccCh%d%s%s", iCh+1, plane[iPlane], side[iSide]),
          Form("occupancy of chamber %d %s %s;occupancy;#events", iCh+1, plane[iPlane], side[iSide]),
          5000, 0., 50.);
        hOcc[iCh][iPlane][iSide]->SetDirectory(0);
        nPads[iCh][iPlane][iSide] = 0;
      }
      hAsym[iCh][iPlane] = new TH1F(Form("hAsymCh%d%s", iCh+1, plane[iPlane]),
        Form("asymmetry of chamber %d %s;asymmetry;#events", iCh+1, plane[iPlane]),
        10000, -50., 50.);
      hAsym[iCh][iPlane]->SetDirectory(0);
      for (Int_t iDE = 0; iDE < nDE[iCh]; ++iDE) {
        hOccDE[iCh][iDE][iPlane] = new TH1F(Form("hOccDE%d%02d%s", iCh+1, iDE, plane[iPlane]),
                                            Form("occupancy of DE %d%02d %s;occupancy;#events", iCh+1, iDE, plane[iPlane]),
                                            1000, 0., 100.);
        hOccDE[iCh][iDE][iPlane]->SetDirectory(0);
        nPadsDE[iCh][iDE][iPlane] = 0;
      }
    }
    for (Int_t iSide = 0; iSide < 3; ++iSide) {
      hDEBad[iCh][iSide] = new TH1F(Form("hDEBadCh%d%s", iCh+1, side[iSide]),
                                    Form("number of DE over max occupancy in chamber %d %s", iCh+1, side[iSide]),
                                    nDE[iCh]+1, -0.5, nDE[iCh]+0.5);
      hDEBad[iCh][iSide]->SetDirectory(0);
    }
  }
  
  Int_t nPadsSpectro[3][3];
  TH1F *hOccSpectro[3][3];
  TH2F *hAsymSpectro[3];
  for (Int_t iPlane = 0; iPlane < 3; ++iPlane) {
    for (Int_t iSide = 0; iSide < 3; ++iSide) {
      hOccSpectro[iPlane][iSide] = new TH1F(Form("hOccSpectro%s%s", plane[iPlane], side[iSide]),
        Form("occupancy of the spectrometer %s %s;occupancy;#events", plane[iPlane], side[iSide]),
        5000, 0., 50.);
      hOccSpectro[iPlane][iSide]->SetDirectory(0);
      nPadsSpectro[iPlane][iSide] = 0;
    }
    hAsymSpectro[iPlane] = new TH2F(Form("hAsymSpectro%s", plane[iPlane]),
                                    Form("asymmetry of the spectrometer %s;occupancy right + left;occupancy right - left", plane[iPlane]),
                                    4000, 0., 20., 4000, -20., 20.);
    hAsymSpectro[iPlane]->SetDirectory(0);
  }
  TH2F* hOccSt12vs345 = new TH2F("hOccSt12vs345", "occupancy stations 1-2 vs stations 3-4-5;occ st345;occ st12",
                                 5000, 0., 50., 5000, 0., 50.);
  hOccSt12vs345->SetDirectory(0);
  TH2F* hOccSt12vsnBadDE = new TH2F("hOccSt12vsnBadDE", "occupancy stations 1-2 vs number of DE over max occupancy;#badDE;occ st12",
                                    51, -0.5, 50.5, 5000, 0., 50.);
  hOccSt12vsnBadDE->SetDirectory(0);
  TH2F* hOccSt345vsnBadDE = new TH2F("hOccSt345vsnBadDE", "occupancy stations 3-4-5 vs number of DE over max occupancy;#badDE;occ st345",
                                     51, -0.5, 50.5, 5000, 0., 50.);
  hOccSt345vsnBadDE->SetDirectory(0);
  TH1F* hDESpectroBad[3];
  TH1F *hChBad[3];
  for (Int_t iSide = 0; iSide < 3; ++iSide) {
    hDESpectroBad[iSide] = new TH1F(Form("hDEBad%s", side[iSide]),
                                    Form("number of DE over max occupancy %s", side[iSide]),
                                    157, -0.5, 156.5);
    hDESpectroBad[iSide]->SetDirectory(0);
    hChBad[iSide] = new TH1F(Form("hChBad%s", side[iSide]),
                             Form("number of chamber with DE over max occupancy %s", side[iSide]),
                             11, -0.5, 10.5);
    hChBad[iSide]->SetDirectory(0);
  }
  TH2F* hDESpectroBadLR = new TH2F("hDEBadLR", "number of DE over max occupancy;right;left",
                                   79, -0.5, 78.5, 79, -0.5, 78.5);
  hDESpectroBadLR->SetDirectory(0);

  // loop over events
  Int_t nPadsCurrent[10][3][3] = {0};
  Int_t nPadsDECurrent[10][26][3] = {0};
  Int_t nPadsSpectroCurrent[3][3] = {0};
  Int_t nDEBadCurrent[10][3] = {0};
  Long64_t nEvents = treeD->GetEntries();
  Long64_t firstEv = 0, lastEv = nEvents-1;
  if (event >= 0 && event < nEvents) {
    firstEv = lastEv = event;
    nEvents = 1;
  }
  for (Long64_t iEv = firstEv; iEv <= lastEv; ++iEv) {
    
    treeD->GetEntry(iEv);
    
    CountFiredPads(digitStore, nPadsCurrent, nPadsDECurrent, occStore, nEvents);
    /*
    if (!(nPadsDECurrent[0][0][0] > 0.06*nPadsDETot[0][0][0] ||
          nPadsDECurrent[0][0][1] > 0.06*nPadsDETot[0][0][1])) continue;
    else printf("high occ ch1... event = %lld\n", iEv);
    */
    // mean occupancy stations 1-2
    Int_t nPadsSt12Current = 0;
    for (Int_t iCh = 0; iCh < 4; ++iCh) {
      nPadsSt12Current += nPadsCurrent[iCh][2][2];
    }
    Float_t occSt12Current = (nPadsSt12Working > 0) ? 100. * nPadsSt12Current / nPadsSt12Working : 0.;
    // if (nPadsSt12Current == 0) {
    //   digitStore->Clear();
    //   continue;
    // }

    // mean occupancy stations 3-4-5
    Int_t nPadsSt345Current = 0;
    for (Int_t iCh = 4; iCh < 10; ++iCh) {
      nPadsSt345Current += nPadsCurrent[iCh][2][2];
    }
    Float_t occSt345Current = (nPadsSt345Current > 0) ? 100. * nPadsSt345Current / nPadsSt345Working : 0.;

    // check occupancy per DE
    Int_t nChBadCurrent[3] = {0, 0, 0};
    Int_t nDESpectroBadCurrent[3] = {0, 0, 0};
    for (Int_t iCh = 0; iCh < 10; ++iCh) {
      for (Int_t iSide = 0; iSide < 3; ++iSide) {
        nDEBadCurrent[iCh][iSide] = 0;
      }
      for (Int_t iDE = 0; iDE < nDE[iCh]; ++iDE) {
        if (100. * nPadsDECurrent[iCh][iDE][0] / nPadsDETot[iCh][iDE][0] > 100. * maxOcc ||
            100. * nPadsDECurrent[iCh][iDE][1] / nPadsDETot[iCh][iDE][1] > 100. * maxOcc) {
          // if (100. * nPadsDECurrent[iCh][iDE][2] / nPadsDETot[iCh][iDE][2] > 100. * maxOcc) {
          nDEBadCurrent[iCh][isOnTheRight(100 * (iCh + 1) + iDE) ? 1 : 0]++;
          nDEBadCurrent[iCh][2]++;
        }
      }
      for (Int_t iSide = 0; iSide < 3; ++iSide) {
        nDESpectroBadCurrent[iSide] += nDEBadCurrent[iCh][iSide];
        if (nDEBadCurrent[iCh][iSide] > 0) {
          nChBadCurrent[iSide]++;
        }
      }
    }

    // reject events with too many high occupancy DE
    // if (nDESpectroBadCurrent[2] > 4) {
    //   digitStore->Clear();
    //   continue;
    // }

    // fill histograms for high occupancy DE
    for (Int_t iSide = 0; iSide < 3; ++iSide) {
      for (Int_t iCh = 0; iCh < 10; ++iCh) {
        hDEBad[iCh][iSide]->Fill(nDEBadCurrent[iCh][iSide]);
      }
      hDESpectroBad[iSide]->Fill(nDESpectroBadCurrent[iSide]);
      hChBad[iSide]->Fill(nChBadCurrent[iSide]);
    }
    hDESpectroBadLR->Fill(nDESpectroBadCurrent[1], nDESpectroBadCurrent[0]);

    // chamber and DE occupancy and chamber asymmetry
    for (Int_t iCh = 0; iCh < 10; ++iCh) {
      for (Int_t iPlane = 0; iPlane < 3; ++iPlane) {
        for (Int_t iSide = 0; iSide < 3; ++iSide) {
          nPads[iCh][iPlane][iSide] += nPadsCurrent[iCh][iPlane][iSide];
          hOcc[iCh][iPlane][iSide]->Fill(100. * nPadsCurrent[iCh][iPlane][iSide] / nPadsTot[iCh][iPlane][iSide]);
        }
        Float_t diff = 100.*(((Float_t)nPadsCurrent[iCh][iPlane][1])/((Float_t)nPadsTot[iCh][iPlane][1]) -
                             ((Float_t)nPadsCurrent[iCh][iPlane][0])/((Float_t)nPadsTot[iCh][iPlane][0]));
        hAsym[iCh][iPlane]->Fill(diff);
        for (Int_t iDE = 0; iDE < nDE[iCh]; ++iDE) {
          nPadsDE[iCh][iDE][iPlane] += nPadsDECurrent[iCh][iDE][iPlane];
          hOccDE[iCh][iDE][iPlane]->Fill(100. * nPadsDECurrent[iCh][iDE][iPlane] / nPadsDETot[iCh][iDE][iPlane]);
        }
      }
    }

    // spectrometer occupancy and asymmetry
    for (Int_t iPlane = 0; iPlane < 3; ++iPlane) {
      for (Int_t iSide = 0; iSide < 3; ++iSide) {
        nPadsSpectroCurrent[iPlane][iSide] = 0;
        for (Int_t iCh = 0; iCh < 10; ++iCh) {
          nPadsSpectroCurrent[iPlane][iSide] += nPadsCurrent[iCh][iPlane][iSide];
        }
        nPadsSpectro[iPlane][iSide] += nPadsSpectroCurrent[iPlane][iSide];
        hOccSpectro[iPlane][iSide]->Fill(100. * nPadsSpectroCurrent[iPlane][iSide] / nPadsSpectroTot[iPlane][iSide]);
      }
      Float_t diff = 100.*(((Float_t)nPadsSpectroCurrent[iPlane][1])/((Float_t)nPadsSpectroTot[iPlane][1]) -
                           ((Float_t)nPadsSpectroCurrent[iPlane][0])/((Float_t)nPadsSpectroTot[iPlane][0]));
      hAsymSpectro[iPlane]->Fill(100.*(((Float_t)nPadsSpectroCurrent[iPlane][1])/((Float_t)nPadsSpectroTot[iPlane][1]) +
                                       ((Float_t)nPadsSpectroCurrent[iPlane][0])/((Float_t)nPadsSpectroTot[iPlane][0])), diff);
      if (iPlane == 2) {
        if (diff < -2.) printf("asym left:  %lld\n", iEv);
        else if (diff > 2.) printf("asym right: %lld\n", iEv);
      }
    }
    hOccSt12vs345->Fill(occSt345Current, occSt12Current);
    hOccSt12vsnBadDE->Fill(nDESpectroBadCurrent[2], occSt12Current);
    hOccSt345vsnBadDE->Fill(nDESpectroBadCurrent[2], occSt345Current);

    digitStore->Clear();
  }
  
  // print results
  for (Int_t iCh = 0; iCh < 10; ++iCh)
    printf("ch%02d: <nPads> = %07.2f --> occ = %04.2f%%, diff = %04.2f%%\n", iCh+1, ((Float_t)nPads[iCh][2][2])/((Float_t)nEvents),
           100.*((Float_t)nPads[iCh][2][2])/((Float_t)nEvents)/((Float_t)nPadsTot[iCh][2][2]),
           100.*(((Float_t)nPads[iCh][2][1])/((Float_t)nPadsTot[iCh][2][1]) -
            ((Float_t)nPads[iCh][2][0])/((Float_t)nPadsTot[iCh][2][0]))/((Float_t)nEvents));
  printf("spectro: <nPads> = %07.2f --> occ = %04.2f%%, diff = %04.2f%%\n", ((Float_t)nPadsSpectro[2][2])/((Float_t)nEvents),
         100.*((Float_t)nPadsSpectro[2][2])/((Float_t)nEvents)/((Float_t)nPadsSpectroTot[2][2]),
         100.*(((Float_t)nPadsSpectro[2][1])/((Float_t)nPadsSpectroTot[2][1]) -
          ((Float_t)nPadsSpectro[2][0])/((Float_t)nPadsSpectroTot[2][0]))/((Float_t)nEvents));
  
  AliMUONTrackerData occupancy("occupancy", "occupancy", 1, kFALSE);
  occupancy.SetDimensionName(0, "occ");
  occupancy.Add(occStore);
  
  // display results
  TCanvas *cOcc = new TCanvas("cOcc", "occupancy per chamber", 1200, 1200);
  cOcc->Divide(5,6);
  for (Int_t iCh = 0; iCh < 10; ++iCh) {
    for (Int_t iPlane = 0; iPlane < 3; ++iPlane) {
      cOcc->cd(10*iPlane+iCh+1);
      gPad->SetLogy();
      hOcc[iCh][iPlane][2]->SetLineColor(1);
      hOcc[iCh][iPlane][2]->Draw();
      hOcc[iCh][iPlane][0]->SetLineColor(2);
      hOcc[iCh][iPlane][0]->Draw("same");
      hOcc[iCh][iPlane][1]->SetLineColor(4);
      hOcc[iCh][iPlane][1]->Draw("same");
    }
  }
  TCanvas *cAsym = new TCanvas("cAsym", "asymmetry per chamber", 1200, 1200);
  cAsym->Divide(5,6);
  for (Int_t iCh = 0; iCh < 10; ++iCh) {
    for (Int_t iPlane = 0; iPlane < 3; ++iPlane) {
      cAsym->cd(10*iPlane+iCh+1);
      gPad->SetLogy();
      hAsym[iCh][iPlane]->Draw();
    }
  }
  TCanvas *cOccSpectro = new TCanvas("cOccSpectro", "occupancy", 1200, 800);
  cOccSpectro->Divide(3,2);
  for (Int_t iPlane = 0; iPlane < 3; ++iPlane) {
    cOccSpectro->cd(iPlane+1);
    gPad->SetLogy();
    hOccSpectro[iPlane][2]->SetLineColor(1);
    hOccSpectro[iPlane][2]->Draw();
    hOccSpectro[iPlane][0]->SetLineColor(2);
    hOccSpectro[iPlane][0]->Draw("same");
    hOccSpectro[iPlane][1]->SetLineColor(4);
    hOccSpectro[iPlane][1]->Draw("same");
  }
  cOccSpectro->cd(4);
  gPad->SetLogz();
  hOccSt12vs345->Draw("colz");
  cOccSpectro->cd(5);
  gPad->SetLogz();
  hOccSt12vsnBadDE->Draw("colz");
  cOccSpectro->cd(6);
  gPad->SetLogz();
  hOccSt345vsnBadDE->Draw("colz");
  TCanvas* cAsymSpectro = new TCanvas("cAsymSpectro", "asymmetry", 1200, 400);
  cAsymSpectro->Divide(3,1);
  for (Int_t iPlane = 0; iPlane < 3; ++iPlane) {
    cAsymSpectro->cd(iPlane+1);
    gPad->SetLogz();
    hAsymSpectro[iPlane]->Draw();
  }
  TCanvas *cOccDE[10];
  for (Int_t iCh = 0; iCh < 10; ++iCh) {
    cOccDE[iCh] = new TCanvas(Form("cOccCh%d",iCh+1), Form("occupancy per DE in chamber %d",iCh+1), 1200, 1200);
    cOccDE[iCh]->Divide(2,nDE[iCh]/2);
    Int_t iDELeft = 0, iDERight = 0;
    for (Int_t iDE = 0; iDE < nDE[iCh]; ++iDE) {
      if (isOnTheRight(100*(iCh+1)+iDE))
        cOccDE[iCh]->cd((iDE > nDE[iCh]/2) ? 2*((5*nDE[iCh]-2)/4-iDE+1) : 2*((nDE[iCh]-2)/4-iDE+1));
      else cOccDE[iCh]->cd(2*(iDE-(nDE[iCh]-2)/4)-1);
      gPad->SetLogy();
      hOccDE[iCh][iDE][2]->SetLineColor(1);
      if (iCh > 3) {
        hOccDE[iCh][iDE][2]->GetXaxis()->SetLabelSize(0.15);
        hOccDE[iCh][iDE][2]->GetYaxis()->SetLabelSize(0.15);
      }
      hOccDE[iCh][iDE][2]->Draw();
      hOccDE[iCh][iDE][0]->SetLineColor(2);
      hOccDE[iCh][iDE][0]->Draw("same");
      hOccDE[iCh][iDE][1]->SetLineColor(4);
      hOccDE[iCh][iDE][1]->Draw("same");
    }
  }
  TCanvas *cDEBad = new TCanvas("cDEBad", "number of DE over max occupancy per chamber", 1200, 400);
  cDEBad->Divide(5,2);
  for (Int_t iCh = 0; iCh < 10; ++iCh) {
    cDEBad->cd(iCh+1);
    gPad->SetLogy();
    hDEBad[iCh][2]->SetLineColor(1);
    hDEBad[iCh][2]->Draw();
    hDEBad[iCh][0]->SetLineColor(2);
    hDEBad[iCh][0]->Draw("same");
    hDEBad[iCh][1]->SetLineColor(4);
    hDEBad[iCh][1]->Draw("same");
  }
  TCanvas *cDESpectroBad = new TCanvas("cDESpectroBad", "number of DE over max occupancy", 1200, 400);
  cDESpectroBad->Divide(3,1);
  cDESpectroBad->cd(1);
  gPad->SetLogy();
  hDESpectroBad[2]->SetLineColor(1);
  hDESpectroBad[2]->Draw();
  hDESpectroBad[0]->SetLineColor(2);
  hDESpectroBad[0]->Draw("same");
  hDESpectroBad[1]->SetLineColor(4);
  hDESpectroBad[1]->Draw("same");
  cDESpectroBad->cd(2);
  gPad->SetLogy();
  hChBad[2]->SetLineColor(1);
  hChBad[2]->Draw();
  hChBad[0]->SetLineColor(2);
  hChBad[0]->Draw("same");
  hChBad[1]->SetLineColor(4);
  hChBad[1]->Draw("same");
  cDESpectroBad->cd(3);
  gPad->SetLogz();
  hDESpectroBadLR->Draw("colz");

  // save occupancy to a file
  TFile *outFile = TFile::Open("occupancy.root", "UPDATE");
  if (outFile && outFile->IsOpen()) {
    occupancy.Write(0x0, TObject::kOverwrite);
    cOcc->Write(0x0, TObject::kOverwrite);
    outFile->Close();
  }
  
  // display occupancy
//  gSystem->Exec("mchview --use occupancy.root");
  
}

//------------------------------------------------------------------
void GetNPads(Int_t nPadsTot[10][3][3], Int_t nPadsDETot[10][26][3])
{
  /// get the total number of pads

  for (Int_t iCh = 0; iCh < AliMUONConstants::NTrackingCh(); iCh++) {

    for (Int_t iPlane = 0; iPlane < 3; ++iPlane) {
      for (Int_t iSide = 0; iSide < 3; ++iSide) {
        nPadsTot[iCh][iPlane][iSide] = 0;
      }
      for (Int_t iDE = 0; iDE < nDE[iCh]; ++iDE) {
        nPadsDETot[iCh][iDE][iPlane] = 0;
      }
    }

    AliMpDEIterator deIt;
    deIt.First(iCh);
    while (!deIt.IsDone()) {

      Int_t deId = deIt.CurrentDEId();
      Int_t iSide = isOnTheRight(deId) ? 1 : 0;
      Int_t iDE = deId % 100;

      const AliMpVSegmentation* seg[2] =
        {AliMpSegmentation::Instance()->GetMpSegmentation(deId, AliMp::kCath0),
         AliMpSegmentation::Instance()->GetMpSegmentation(deId, AliMp::kCath1)};

      for (Int_t iCath = 0; iCath < 2; iCath++) {

        Int_t iPlane = seg[iCath]->PlaneType();
        Int_t nPads = seg[iCath]->NofPads();

        nPadsTot[iCh][iPlane][iSide] += nPads;
        nPadsTot[iCh][iPlane][2] += nPads;
        nPadsTot[iCh][2][iSide] += nPads;
        nPadsTot[iCh][2][2] += nPads;

        nPadsDETot[iCh][iDE][iPlane] += nPads;
        nPadsDETot[iCh][iDE][2] += nPads;
      }

      deIt.Next();
    }
  }
}

//------------------------------------------------------------------
void GetNWorkingPads(Int_t run, Int_t nPadsTot[10][3][3], Int_t nPadsDETot[10][26][3], AliMUON2DMap& occStore)
{
  /// get the number of working pads

  // prepare digit calibrator
  AliMUONRecoParam* recoParam = AliMUONCDB::LoadRecoParam();
  AliMUONCalibrationData* calibrationData = new AliMUONCalibrationData(run);
  AliMUONDigitCalibrator* digitCalibrator = new AliMUONDigitCalibrator(*calibrationData, recoParam);

  // loop over DEs
  Double_t activeArea[10][2];
  Double_t minPadSize[10][2];
  Double_t maxPadSize[10][2];
  for (Int_t iCh = 0; iCh < AliMUONConstants::NTrackingCh(); iCh++) {

    for (Int_t iPlane = 0; iPlane < 3; ++iPlane) {
      for (Int_t iSide = 0; iSide < 3; ++iSide)
        nPadsTot[iCh][iPlane][iSide] = 0;
      for (Int_t iDE = 0; iDE < nDE[iCh]; ++iDE)
        nPadsDETot[iCh][iDE][iPlane] = 0;
    }
    for (Int_t iPlane = 0; iPlane < 2; ++iPlane) {
      activeArea[iCh][iPlane] = 0.;
      minPadSize[iCh][iPlane] = 999.;
      maxPadSize[iCh][iPlane] = 0.;
    }

    AliMpDEIterator deIt;
    deIt.First(iCh);
    while (!deIt.IsDone()) {

      Int_t deId = deIt.CurrentDEId();
      Int_t iSide = isOnTheRight(deId) ? 1 : 0;
      Int_t iDE = deId % 100;

      const AliMpVSegmentation* seg[2] =
        {AliMpSegmentation::Instance()->GetMpSegmentation(deId, AliMp::kCath0),
         AliMpSegmentation::Instance()->GetMpSegmentation(deId, AliMp::kCath1)};

      // loop over cathods
      for (Int_t iCath = 0; iCath < 2; iCath++) {

        Int_t iPlane = seg[iCath]->PlaneType();

        // loop over pads
        AliMpVPadIterator* padIt = seg[iCath]->CreateIterator();
        padIt->First();
        while (!padIt->IsDone()) {

          AliMpPad pad = padIt->CurrentItem();
          Int_t manuId = pad.GetManuId();
          Int_t manuChannel = pad.GetManuChannel();

          if (digitCalibrator->IsValidDigit(deId, manuId, manuChannel)) {

            // make sure all valid pads are taken into account
            AliMUONVCalibParam* c = static_cast<AliMUONVCalibParam*>(occStore.FindObject(deId, manuId));
            if (!c) {
              c = new AliMUONCalibParamND(1, AliMpConstants::ManuNofChannels(), deId, manuId, AliMUONVCalibParam::InvalidFloatValue());
              occStore.Add(c);
            }
            c->SetValueAsDouble(manuChannel, 0, 0.);

            // count valid digits
            ++nPadsTot[iCh][iPlane][iSide];
            ++nPadsTot[iCh][iPlane][2];
            ++nPadsTot[iCh][2][iSide];
            ++nPadsTot[iCh][2][2];

            ++nPadsDETot[iCh][iDE][iPlane];
            ++nPadsDETot[iCh][iDE][2];

            // increase active area
            activeArea[iCh][iPlane] += 4. * pad.GetDimensionX() * pad.GetDimensionY();
            if (iPlane == 1) {
              if (2. * pad.GetDimensionX() + 1.e-4 < minPadSize[iCh][1])
                minPadSize[iCh][1] = 2. * pad.GetDimensionX();
              if (2. * pad.GetDimensionX() - 1.e-4 > maxPadSize[iCh][1])
                maxPadSize[iCh][1] = 2. * pad.GetDimensionX();
            } else {
              if (2. * pad.GetDimensionY() + 1.e-4 < minPadSize[iCh][0])
                minPadSize[iCh][0] = 2. * pad.GetDimensionY();
              if (2. * pad.GetDimensionY() - 1.e-4 > maxPadSize[iCh][0])
                maxPadSize[iCh][0] = 2. * pad.GetDimensionY();
            }
          }

          padIt->Next();
        }
      }

      deIt.Next();
    }

    printf("active area ch%0d: NB = %6.0f cm2; B = %6.0f cm2\n", iCh + 1, activeArea[iCh][1], activeArea[iCh][0]);
    printf("min pad size ch%0d: NB = %3.1f cm; B = %3.1f cm\n", iCh + 1, minPadSize[iCh][1], minPadSize[iCh][0]);
    printf("max pad size ch%0d: NB = %3.1f cm; B = %3.1f cm\n", iCh + 1, maxPadSize[iCh][1], maxPadSize[iCh][0]);
  }
}

//------------------------------------------------------------------
void CountFiredPads(AliMUONVDigitStore *digitStore, Int_t nPadsCurrent[10][3][3], Int_t nPadsDECurrent[10][26][3],
                    AliMUON2DMap &occStore, Long64_t nEvents)
{
  /// print digits occupancy per chamber
  
  for (Int_t iCh = 0; iCh < 10; ++iCh)
    for (Int_t iPlane = 0; iPlane < 3; ++iPlane) {
      for (Int_t iSide = 0; iSide < 3; ++iSide)
        nPadsCurrent[iCh][iPlane][iSide] = 0;
      for (Int_t iDE = 0; iDE < nDE[iCh]; ++iDE)
        nPadsDECurrent[iCh][iDE][iPlane] = 0;
    }
  
  // loop over digits
  AliMUONVDigit* digit = 0x0;
  TIter nextDigit(digitStore->CreateTrackerIterator());
  while ((digit = static_cast<AliMUONVDigit*>(nextDigit()))) {
    
    if (digit->Charge() <= 0.) continue;
    
    Int_t deId = digit->DetElemId();
    Int_t chId = deId/100-1;
    Int_t iSide = isOnTheRight(deId) ? 1 : 0;
    Int_t iDE = deId%100;
    Int_t manuId = digit->ManuId();
    Int_t iPlane = (manuId & AliMpConstants::ManuMask(AliMp::kNonBendingPlane)) ? 1 : 0;
    Int_t manuChannel = digit->ManuChannel();
    
    // count fired digits
    AliMUONVCalibParam *c = static_cast<AliMUONVCalibParam*>(occStore.FindObject(deId, manuId));
    if (!c) {
      c = new AliMUONCalibParamND(1, AliMpConstants::ManuNofChannels(), deId, manuId, AliMUONVCalibParam::InvalidFloatValue());
      occStore.Add(c);
    }
    Double_t occ = c->ValueAsDouble(manuChannel, 0);
    if (occ >= AliMUONVCalibParam::InvalidFloatValue()) {
      printf("invalid pad fired!\n");
      occ = 0.;
    }
    c->SetValueAsDouble(manuChannel, 0, occ + 1./nEvents);
    
    ++nPadsCurrent[chId][iPlane][iSide];
    ++nPadsCurrent[chId][iPlane][2];
    ++nPadsCurrent[chId][2][iSide];
    ++nPadsCurrent[chId][2][2];
    
    ++nPadsDECurrent[chId][iDE][iPlane];
    ++nPadsDECurrent[chId][iDE][2];
    
  }
  
}

//------------------------------------------------------------------
Bool_t isOnTheRight(Int_t deId)
{
  /// return kTRUE if this DE is on the right side of the spectrometer
  
  Int_t chId = deId/100-1;
  Int_t deNum = deId%100;
  
  if (chId < 4) return (deNum == 0 || deNum == 3);
  else if (chId < 6) return (deNum < 5 || deNum > 13);
  else return (deNum < 7 || deNum > 19);
  
}

