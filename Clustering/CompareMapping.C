#include <stdexcept>

#include <TH1D.h>
#include <TCanvas.h>

#include "AliCDBManager.h"

#include "AliMpDEIterator.h"
#include "AliMpVSegmentation.h"
#include "AliMpSegmentation.h"
#include "AliMpVPadIterator.h"
#include "AliMpPad.h"

#include "AliMUONCDB.h"
#include "AliMUONPad.h"
#include "AliMUONConstants.h"

#include "MCHMappingInterface/Segmentation.h"

void CompareMapping()
{
  /// Compare the mapping between O2 and AliRoot
  /// O2 mapping need to be loaded before: gSystem->Load("libO2MCHMappingImpl3")

  // load AliRoot mapping
  AliCDBManager* man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALIROOT_OCDB_ROOT/OCDB");
  man->SetRun(0);
  if (!AliMUONCDB::LoadMapping()) {
    return;
  }

  TH1D* hXDiff = new TH1D("hXDiff", "hXDiff", 2001, -0.00010005, 0.00010005);
  TH1D* hYDiff = new TH1D("hYDiff", "hYDiff", 2001, -0.00010005, 0.00010005);
  TH1D* hXSizeDiff = new TH1D("hXSizeDiff", "hXSizeDiff", 2001, -0.0000010005, 0.0000010005);
  TH1D* hYSizeDiff = new TH1D("hYSizeDiff", "hYSizeDiff", 2001, -0.0000010005, 0.0000010005);

  // loop over every pad and compare position and size
  for (int iCh = 0; iCh < AliMUONConstants::NTrackingCh(); ++iCh) {

    AliMpDEIterator deIt;
    deIt.First(iCh);
    while (!deIt.IsDone()) {

      int deId = deIt.CurrentDEId();
      const AliMpVSegmentation* seg[2] =
        {AliMpSegmentation::Instance()->GetMpSegmentation(deId, AliMp::kCath0),
         AliMpSegmentation::Instance()->GetMpSegmentation(deId, AliMp::kCath1)};

      for (int iCath = 0; iCath < 2; ++iCath) {

        AliMpVPadIterator* padIt = seg[iCath]->CreateIterator();
        padIt->First();
        while (!padIt->IsDone()) {

          AliMpPad pad = padIt->CurrentItem();

          int manuId = pad.GetManuId();
          int manuCh = pad.GetManuChannel();
          auto& segO2 = o2::mch::mapping::segmentation(deId);
          int padId = segO2.findPadByFEE(manuId, manuCh);
          if (!segO2.isValid(padId)) {
            throw std::runtime_error("pad not found");
          }

          hXDiff->Fill(segO2.padPositionX(padId) - pad.GetPositionX());
          hYDiff->Fill(segO2.padPositionY(padId) - pad.GetPositionY());
          hXSizeDiff->Fill(segO2.padSizeX(padId) / 2. - pad.GetDimensionX());
          hYSizeDiff->Fill(segO2.padSizeY(padId) / 2. - pad.GetDimensionY());

          double x1 = round(pad.GetPositionX() * 1.e4) / 1.e4;
          double y1 = round(pad.GetPositionY() * 1.e5) / 1.e5;
          double dx1 = round(pad.GetDimensionX() * 1.e7) / 1.e7;
          double dy1 = round(pad.GetDimensionY() * 1.e6) / 1.e6;

          double x2 = round(segO2.padPositionX(padId) * 1.e4) / 1.e4;
          double y2 = round(segO2.padPositionY(padId) * 1.e5) / 1.e5;
          double dx2 = round(segO2.padSizeX(padId) / 2. * 1.e7) / 1.e7;
          double dy2 = round(segO2.padSizeY(padId) / 2. * 1.e6) / 1.e6;

          if (x1 != x2) {
            printf("DE %d: x = %.17f != %.17f\n", deId, x1, x2);
          }
          if (y1 != y2) {
            printf("DE %d: y = %.17f != %.17f\n", deId, y1, y2);
          }
          if (dx1 != dx2) {
            printf("DE %d: dx = %.17f != %.17f\n", deId, dx1, dx2);
          }
          if (dy1 != dy2) {
            printf("DE %d: dy = %.17f != %.17f\n", deId, dy1, dy2);
          }

          padIt->Next();
        }
      }
      deIt.Next();
    }
  }

  // draw differences
  TCanvas* c = new TCanvas("c", "c", 10, 10, 600, 600);
  c->Divide(2, 2);
  c->cd(1);
  hXDiff->Draw();
  c->cd(2);
  hYDiff->Draw();
  c->cd(3);
  hXSizeDiff->Draw();
  c->cd(4);
  hYSizeDiff->Draw();
}