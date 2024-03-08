//
//  DrawPreClustersSize.C
//  aliroot_dev
//
//  Created by philippe pillot on 17/02/2014.
//  Copyright (c) 2014 Philippe Pillot. All rights reserved.
//

#include <stdio.h>
#include <set>
#include <functional>

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TPad.h>

#include "AliCodeTimer.h"
#include "AliCDBManager.h"

#include "AliMUONCDB.h"
#include "AliMUONConstants.h"
#include "AliMUONVDigit.h"
#include "AliMUONVDigitStore.h"
#include "AliMUONVCluster.h"
#include "AliMUONVClusterStore.h"

#include "AliMpConstants.h"
#include "AliMpDEIterator.h"
#include "AliMpVSegmentation.h"
#include "AliMpSegmentation.h"

const Int_t nDE[10] = {4, 4, 4, 4, 18, 18, 26, 26, 26, 26};

const Double_t maxOccEv = 0.2; // maximum occupancy per DE to reject the full event
const Double_t maxOccDE = 0.2; // maximum occupancy per DE to reject it

void GetNPads(Int_t nPadsDETot[10][26][2]);
void CountFiredPads(AliMUONVDigitStore* digitStore, Int_t nPadsDECurrent[10][26][2]);
Bool_t isOnTheRight(Int_t deId);

//------------------------------------------------------------------
void DrawPreClustersSize(const char* clusterFileName, const char* digitFileName)
{
  /// Draw the distribution of preclusters' size

  AliCodeTimerAutoGeneral("", 0);

  // load mapping
  AliCDBManager* man = AliCDBManager::Instance();
  //man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  man->SetDefaultStorage("local:///dev/null");
  man->SetSnapshotMode("OCDB.root");
  man->SetRun(295584);
  // man->SetDefaultStorage("local://$ALIROOT_OCDB_ROOT/OCDB");
  // man->SetRun(0);
  if (!AliMUONCDB::LoadMapping()) {
    return;
  }

  // total number of pads per DE
  Int_t nPadsDETot[10][26][2] = {0};
  GetNPads(nPadsDETot);

  // create histograms
  TH1F* hSpectro = new TH1F("hSpectro", "hSpectro", 10000, 0, 10000);
  TH1F* hCh[10];
  for (Int_t ich = 0; ich < 10; ich++) {
    hCh[ich] = new TH1F(Form("hCh%d", ich + 1), Form("hCh%d", ich + 1), 10000, 0, 10000);
  }
  TH2F* hChargesSt12 = new TH2F("hChargesSt12", "hChargesSt12;#digits;charge", 5000, 0, 5000, 700, 0, 700);
  TH2F* hChargesSt345 = new TH2F("hChargesSt345", "hChargesSt345;#digits;charge", 5000, 0, 5000, 700, 0, 700);
  TH2F* hMaxChargeSt12 = new TH2F("hMaxChargeSt12", "hMaxChargeSt12;#digits;max charge", 5000, 0, 5000, 700, 0, 700);
  TH2F* hMaxChargeSt345 = new TH2F("hMaxChargeSt345", "hMaxChargeSt345;#digits;max charge", 5000, 0, 5000, 700, 0, 700);
  TH2F* hChargeAsymm2D = new TH2F("hChargeAsymm2D", "charge asymmetry vs charge;charge;asymmetry", 200, 0., 20000., 201, -1.005, 1.005);

  // read preclusters
  TFile* clusterFile = new TFile(clusterFileName);
  if (!clusterFile)
    return;
  TTree* treeR = static_cast<TTree*>(clusterFile->Get("TreeR"));
  if (!treeR)
    return;
  AliMUONVClusterStore* clusterStore = AliMUONVClusterStore::Create(*treeR);
  clusterStore->Connect(*treeR);

  // read digits
  TFile* digitFile = TFile::Open(digitFileName);
  if (!digitFile)
    return;
  TTree* treeD = static_cast<TTree*>(digitFile->Get("TreeD"));
  if (!treeD)
    return;
  AliMUONVDigitStore* digitStore = AliMUONVDigitStore::Create(*treeD);
  digitStore->Connect(*treeD);

  // loop over events
  Int_t nDigTot = 0;
  Int_t nPreClMore1000Dig = 0;
  Int_t nPreClMore10000Dig = 0;
  Int_t nDigInPreClMore1000Dig = 0;
  Int_t nDigInPreClMore10000Dig = 0;
  Int_t nPadsDECurrent[10][26][2] = {0};
  Double_t occDECurrent[10][26][2] = {0};
  Bool_t badDECurrent[10][26] = {false};
  Int_t nRejectedDE[10][26] = {0};
  Int_t nRejectedEvents = 0;
  Int_t nRejectedDEEvents = 0;
  std::multiset<float, std::greater<float>> charges{}; // list of charge of digits in decreasing order
  Long64_t nEvents = treeR->GetEntries();
  for (Long64_t iEv = 0; iEv < nEvents; ++iEv) {

    // if (iEv != 6755) {
    //   continue;
    // }

    printf("\rprocessing... %lld%%", 100 * (iEv + 1) / nEvents);

    clusterStore->Clear();
    digitStore->Clear();
    treeR->GetEntry(iEv);
    treeD->GetEntry(iEv);

    // check correlated DE occupancy in stations 3-4-5
    CountFiredPads(digitStore, nPadsDECurrent);
    int nBadDECurrent[2] = {0, 0};
    for (Int_t iCh = 4; iCh < 10; ++iCh) {
      for (Int_t iDE = 0; iDE < nDE[iCh]; ++iDE) {
        for (Int_t iPlane = 0; iPlane < 2; ++iPlane) {
          occDECurrent[iCh][iDE][iPlane] = static_cast<Double_t>(nPadsDECurrent[iCh][iDE][iPlane]) / static_cast<Double_t>(nPadsDETot[iCh][iDE][iPlane]);
        }
        if (occDECurrent[iCh][iDE][0] > maxOccEv || occDECurrent[iCh][iDE][1] > maxOccEv) {
          ++nBadDECurrent[isOnTheRight(100 * (iCh + 1) + iDE) ? 1 : 0];
        }
      }
    }
    if (nBadDECurrent[0] + nBadDECurrent[1] > 4) {
      ++nRejectedEvents;
      continue;
    }

    // check individual DE occupancy in stations 3-4-5
    Bool_t rejectedDEEvent = false;
    for (Int_t iCh = 4; iCh < 10; ++iCh) {
      for (Int_t iDE = 0; iDE < nDE[iCh]; ++iDE) {
        badDECurrent[iCh][iDE] = (occDECurrent[iCh][iDE][0] > maxOccDE || occDECurrent[iCh][iDE][1] > maxOccDE);
        if (badDECurrent[iCh][iDE]) {
          ++nRejectedDE[iCh][iDE];
          rejectedDEEvent = true;
        }
      }
    }
    if (rejectedDEEvent) {
      ++nRejectedDEEvents;
    }

    // loop over preclusters
    AliMUONVCluster* cluster = 0x0;
    TIter nextCluster(clusterStore->CreateIterator());
    while ((cluster = static_cast<AliMUONVCluster*>(nextCluster()))) {

      Int_t deId = cluster->GetDetElemId();
      Int_t chId = deId / 100 - 1;
      Int_t iDE = deId % 100;

      // check DE occupancy in stations 3-4-5
      if (chId > 3 && badDECurrent[chId][iDE]) {
        continue;
      }

      Int_t nDigits = cluster->GetNDigits();
      nDigTot += nDigits;

      hSpectro->Fill(nDigits);
      hCh[chId]->Fill(nDigits);

      if (nDigits > 1000) {
        nPreClMore1000Dig++;
        nDigInPreClMore1000Dig += nDigits;
      }

      if (nDigits > 10000) {
        nPreClMore10000Dig++;
        nDigInPreClMore10000Dig += nDigits;
      }

      charges.clear();
      double chargeTot[2] = {0., 0.};
      for (int iDigit = 0; iDigit < nDigits; ++iDigit) {
        AliMUONVDigit* digit = digitStore->FindObject(cluster->GetDigitId(iDigit));
        if (!digit) {
          printf("missing digit\n");
          return;
        }
        charges.insert(digit->Charge());
        Int_t iPlane = (digit->ManuId() & AliMpConstants::ManuMask(AliMp::kNonBendingPlane)) ? 1 : 0;
        chargeTot[iPlane] += digit->Charge();
      }

      TH2F* h = (chId < 4) ? hChargesSt12 : hChargesSt345;
      for (const auto charge : charges) {
        h->Fill(nDigits, charge);
      }
      h = (chId < 4) ? hMaxChargeSt12 : hMaxChargeSt345;
      h->Fill(nDigits, *charges.begin());
      double chargeTotMean = 0.5 * (chargeTot[1] + chargeTot[0]);
      double chargeTotAsymm = (chargeTot[1] - chargeTot[0]) / (chargeTot[1] + chargeTot[0]);
      hChargeAsymm2D->Fill(chargeTotMean, chargeTotAsymm);
    }
  }

  printf("\rtotal number of digits = %d\n", nDigTot);
  printf("number of preclusters > 1000 digits = %d (%d digits (%4.2f%%) inside those)\n",
         nPreClMore1000Dig, nDigInPreClMore1000Dig, 100. * nDigInPreClMore1000Dig / nDigTot);
  printf("number of preclusters > 10000 digits = %d (%d digits (%4.2f%%) inside those)\n",
         nPreClMore10000Dig, nDigInPreClMore10000Dig, 100. * nDigInPreClMore10000Dig / nDigTot);
  printf("number of rejected events = %d / %lld (%4.2f%%)\n", nRejectedEvents, nEvents, 100. * nRejectedEvents / nEvents);
  printf("number of events with rejected DE = %d / %lld (%4.2f%%)\n", nRejectedDEEvents, nEvents, 100. * nRejectedDEEvents / nEvents);
  printf("DEId ");
  for (Int_t iDE = 0; iDE < 26; ++iDE) {
    printf("  %2d   ", iDE);
  }
  for (Int_t iCh = 4; iCh < 10; ++iCh) {
    printf("\nCh%02d ", iCh + 1);
    for (Int_t iDE = 0; iDE < nDE[iCh]; ++iDE) {
      printf("%5.2f%% ", 100. * nRejectedDE[iCh][iDE] / nEvents);
    }
  }
  printf("\n");

  // draw results
  TCanvas* cSpectro = new TCanvas("cSpectro", "cSpectro", 10, 10, 500, 500);
  cSpectro->cd();
  gPad->SetLogy();
  hSpectro->Draw();

  TCanvas* cCh = new TCanvas("cCh", "cCh", 10, 10, 1250, 500);
  cCh->Divide(5, 2);
  for (Int_t ich = 0; ich < 10; ich++) {
    cCh->cd(ich + 1);
    gPad->SetLogy();
    hCh[ich]->Draw();
  }

  TCanvas* cCharges = new TCanvas("cCharges", "cCharges", 10, 10, 1000, 1000);
  cCharges->Divide(2, 2);
  cCharges->cd(1);
  gPad->SetLogz();
  hChargesSt12->Draw("colz");
  cCharges->cd(2);
  gPad->SetLogz();
  hChargesSt345->Draw("colz");
  cCharges->cd(3);
  gPad->SetLogz();
  hMaxChargeSt12->Draw("colz");
  cCharges->cd(4);
  gPad->SetLogz();
  hMaxChargeSt345->Draw("colz");

  TCanvas* cChargeAsymm = new TCanvas("cChargeAsymm", "cChargeAsymm");
  gPad->SetLogz();
  hChargeAsymm2D->Draw("colz");

  AliCodeTimer::Instance()->Print();
}

//------------------------------------------------------------------
void GetNPads(Int_t nPadsDETot[10][26][2])
{
  /// get the total number of pads

  for (Int_t iCh = 0; iCh < AliMUONConstants::NTrackingCh(); ++iCh) {

    for (Int_t iDE = 0; iDE < nDE[iCh]; ++iDE) {
      for (Int_t iPlane = 0; iPlane < 2; ++iPlane) {
        nPadsDETot[iCh][iDE][iPlane] = 0;
      }
    }

    AliMpDEIterator deIt;
    deIt.First(iCh);
    while (!deIt.IsDone()) {

      Int_t deId = deIt.CurrentDEId();
      Int_t iDE = deId % 100;

      const AliMpVSegmentation* seg[2] =
        {AliMpSegmentation::Instance()->GetMpSegmentation(deId, AliMp::kCath0),
         AliMpSegmentation::Instance()->GetMpSegmentation(deId, AliMp::kCath1)};

      for (Int_t iCath = 0; iCath < 2; iCath++) {
        nPadsDETot[iCh][iDE][seg[iCath]->PlaneType()] += seg[iCath]->NofPads();
      }

      deIt.Next();
    }
  }
}

//------------------------------------------------------------------
void CountFiredPads(AliMUONVDigitStore* digitStore, Int_t nPadsDECurrent[10][26][2])
{
  /// get the total number of fired pads

  for (Int_t iCh = 0; iCh < 10; ++iCh) {
    for (Int_t iDE = 0; iDE < nDE[iCh]; ++iDE) {
      for (Int_t iPlane = 0; iPlane < 2; ++iPlane) {
        nPadsDECurrent[iCh][iDE][iPlane] = 0;
      }
    }
  }

  // loop over digits
  AliMUONVDigit* digit = 0x0;
  TIter nextDigit(digitStore->CreateTrackerIterator());
  while ((digit = static_cast<AliMUONVDigit*>(nextDigit()))) {

    if (digit->Charge() <= 0.) {
      continue;
    }

    Int_t deId = digit->DetElemId();
    Int_t chId = deId / 100 - 1;
    Int_t iDE = deId % 100;
    Int_t iPlane = (digit->ManuId() & AliMpConstants::ManuMask(AliMp::kNonBendingPlane)) ? 1 : 0;

    ++nPadsDECurrent[chId][iDE][iPlane];
  }
}

//------------------------------------------------------------------
Bool_t isOnTheRight(Int_t deId)
{
  /// return kTRUE if this DE is on the right side of the spectrometer

  Int_t chId = deId / 100 - 1;
  Int_t deNum = deId % 100;

  if (chId < 4) {
    return (deNum == 0 || deNum == 3);
  } else if (chId < 6) {
    return (deNum < 5 || deNum > 13);
  } else {
    return (deNum < 7 || deNum > 19);
  }
}
