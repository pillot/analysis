//
//  CompareNeighbours.C
//  aliroot_dev
//
//  Created by philippe pillot on 18/02/2014.
//  Copyright (c) 2014 Philippe Pillot. All rights reserved.
//

#include <list>

#include <TFile.h>
#include <TTree.h>
#include <TExMap.h>
#include "TObjArray.h"

#include "AliCDBManager.h"
#include "AliCodeTimer.h"

#include "AliMUONCDB.h"
#include "AliMUONCalibParamNI.h"
#include "AliMUONCalibrationData.h"
#include "AliMUONConstants.h"
#include "AliMUONVDigit.h"
#include "AliMUON2DMap.h"
#include "AliMUONTrackerData.h"

#include "AliMpConstants.h"
#include "AliMpDEIterator.h"
#include "AliMpVSegmentation.h"
#include "AliMpSegmentation.h"
#include "AliMpVPadIterator.h"
#include "AliMpPad.h"
#include "AliMpManuIterator.h"


struct mpPad {
  UInt_t id; // unique ID
  Int_t nNeighbours1; // number of neighbours
  Int_t neighbours1[10]; // indices of neighbours in array stored in mpDE
  Int_t nNeighbours2; // number of neighbours
  Int_t neighbours2[10]; // indices of neighbours in array stored in mpDE
  Float_t area[2][2]; // 2D area
};

struct mpDE {
  Int_t id; // unique ID
  Int_t nPads[2]; // number of pads on each plane
  mpPad *pads[2]; // array of pads on each plane
  TExMap padIndices[2]; // indices of pads from their ID
};

const Int_t nDEs = 156;
mpDE mpDEs[nDEs];


void MakeNeighbourStore(AliMUONVStore& neighbourStore);
void CreateMapping(AliMUONVStore* neighbours);
void FindNeighbours(Int_t iDE, Int_t iPlane);
Bool_t AreOverlapping(Float_t area1[2][2], Float_t area2[2][2], Float_t precision);
void DrawDifferences();


//------------------------------------------------------------------
void CompareNeighbours()
{
  /// fill neighbours from OCDB or from overlapping pad and show the difference
  
  // load mapping locally and create the local structure
  AliCDBManager* man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  man->SetRun(0);
  if (!AliMUONCDB::LoadMapping()) return;
  //AliMUONVStore* neighbours = AliMUONCalibrationData::CreateNeighbours(0);
  AliMUONVStore* neighbours = new AliMUON2DMap(true);
  MakeNeighbourStore(*neighbours);
  CreateMapping(neighbours);
  delete neighbours;
  
  DrawDifferences();
  
  AliCodeTimer::Instance()->Print();
  
}

//------------------------------------------------------------------
void MakeNeighbourStore(AliMUONVStore& neighbourStore)
{
  /// Fill the neighbours store with, for each channel, a TObjArray of its
  /// neighbouring pads (including itself)
  
  AliCodeTimerAutoGeneral("",0);
  
  if (!AliMUONCDB::CheckMapping()) return;
  
  AliInfoGeneral("AliMUONCDB", "Generating NeighbourStore. This will take a while. Please be patient.");
  
  Int_t nchannels(0);
  
  TObjArray tmp;
  
  Int_t detElemId;
  Int_t manuId;
  
  AliMpManuIterator it;
  
  while ( it.Next(detElemId,manuId) )
  {
    const AliMpVSegmentation* seg =
    AliMpSegmentation::Instance()->GetMpSegmentationByElectronics(detElemId,manuId);
    
    AliMUONVCalibParam* calibParam = static_cast<AliMUONVCalibParam*>(neighbourStore.FindObject(detElemId,manuId));
    if (!calibParam)
    {
      Int_t dimension(13);
      Int_t size(AliMpConstants::ManuNofChannels());
      Int_t defaultValue(-1);
      Int_t packingFactor(size);
      
      calibParam = new AliMUONCalibParamNI(dimension,size,detElemId,manuId,defaultValue,packingFactor);
      Bool_t ok = neighbourStore.Add(calibParam);
      if (!ok)
      {
        AliErrorGeneral("AliMUONCDB", Form("Could not set DetElemId=%d manuId=%d",detElemId,manuId));
        return;
      }
    }
    
    for ( Int_t manuChannel = 0; manuChannel < AliMpConstants::ManuNofChannels(); ++manuChannel )
    {
      AliMpPad pad = seg->PadByLocation(manuId,manuChannel,kFALSE);
      
      if (pad.IsValid())
      {
        ++nchannels;
        
        seg->GetNeighbours(pad,tmp,true,true);
        Int_t nofPadNeighbours = tmp.GetEntriesFast();
        
        for ( Int_t i = 0; i < nofPadNeighbours; ++i )
        {
          AliMpPad* p = static_cast<AliMpPad*>(tmp.UncheckedAt(i));
          Int_t x;
          //          Bool_t ok =
          calibParam->PackValues(p->GetManuId(),p->GetManuChannel(),x);
          //          if (!ok)
          //          {
          //            AliError("Could not pack value. Something is seriously wrong. Please check");
          //            StdoutToAliError(pad->Print(););
          //            return -1;
          //          }
          calibParam->SetValueAsInt(manuChannel,i,x);
        }
      }
    }
  }
  
}

//------------------------------------------------------------------
void CreateMapping(AliMUONVStore* neighbours)
{
  /// fill the mpDE and mpPad structures once for all
  
  AliCodeTimerAutoGeneral("",0);
  
  // loop over DEs
  Int_t iDE = 0;
  Int_t packValue2 = -1, manuId2 = -1, manuChannel2 = -1;
  UInt_t padId = 0;
  for (Int_t iCh = 0; iCh < AliMUONConstants::NTrackingCh(); iCh++) {
    AliMpDEIterator deIt;
    deIt.First(iCh);
    while (!deIt.IsDone()) {
      
      Int_t deId = deIt.CurrentDEId();
      mpDEs[iDE].id = deId;
      
      const AliMpVSegmentation* seg[2] =
      { AliMpSegmentation::Instance()->GetMpSegmentation(deId,AliMp::kCath0),
        AliMpSegmentation::Instance()->GetMpSegmentation(deId,AliMp::kCath1)
      };
      
      // loop over cathods
      for (Int_t iCath = 0; iCath < 2; iCath++) {
        
        Int_t iPlane = seg[iCath]->PlaneType();
        Int_t nPads = seg[iCath]->NofPads();
        mpDEs[iDE].nPads[iPlane] = nPads;
        mpDEs[iDE].pads[iPlane] = new mpPad[nPads];
        
        // 1st loop over pads to associate an index to their Id
        Int_t iPad = 0;
        AliMpVPadIterator* padIt = seg[iCath]->CreateIterator();
        padIt->First();
        while (!padIt->IsDone()) {
          
          AliMpPad pad = padIt->CurrentItem();
          
          mpDEs[iDE].pads[iPlane][iPad].nNeighbours1 = 0;
          mpDEs[iDE].pads[iPlane][iPad].nNeighbours2 = 0;
          mpDEs[iDE].pads[iPlane][iPad].area[0][0] = pad.GetPositionX() - pad.GetDimensionX();
          mpDEs[iDE].pads[iPlane][iPad].area[0][1] = pad.GetPositionX() + pad.GetDimensionX();
          mpDEs[iDE].pads[iPlane][iPad].area[1][0] = pad.GetPositionY() - pad.GetDimensionY();
          mpDEs[iDE].pads[iPlane][iPad].area[1][1] = pad.GetPositionY() + pad.GetDimensionY();
          
          padId = AliMUONVDigit::BuildUniqueID(deId, pad.GetManuId(), pad.GetManuChannel(), iCath);
          mpDEs[iDE].pads[iPlane][iPad].id = padId;
          mpDEs[iDE].padIndices[iPlane].Add(padId, iPad);
          
          iPad++;
          padIt->Next();
          
        }
        
        // 2nd loop over pads to fill the neighbours
        iPad = 0;
        padIt->First();
        while (!padIt->IsDone()) {
          
          AliMpPad pad = padIt->CurrentItem();
          AliMUONVCalibParam* calibParam = static_cast<AliMUONVCalibParam*>(neighbours->FindObject(deId, pad.GetManuId()));
          
          for (Int_t iNeighbour = 1; iNeighbour < calibParam->Dimension(); iNeighbour++) {
            
            packValue2 = calibParam->ValueAsInt(pad.GetManuChannel(),iNeighbour);
            if (packValue2 <= 0) continue;
            
            calibParam->UnpackValue(packValue2, manuId2, manuChannel2);
            padId = AliMUONVDigit::BuildUniqueID(deId, manuId2, manuChannel2, iCath);
            
            mpDEs[iDE].pads[iPlane][iPad].neighbours1[mpDEs[iDE].pads[iPlane][iPad].nNeighbours1++] = mpDEs[iDE].padIndices[iPlane].GetValue(padId);
            
          }
          
          iPad++;
          padIt->Next();
          
        }
        
        // 2nd way of finding neighbours
        FindNeighbours(iDE, iPlane);
        
      }
      
      iDE++;
      deIt.Next();
      
    }
    
  }
  
}

//------------------------------------------------------------------
void FindNeighbours(Int_t iDE, Int_t iPlane)
{
  /// fill the mpPad neighbours structures of given DE/plane
  
  AliCodeTimerAutoGeneral("",0);
  
  // loop over pads
  for (Int_t iPad1 = 0; iPad1 < mpDEs[iDE].nPads[iPlane]; ++iPad1) {
    
    // loop over next pads to find the neighbours
    for (Int_t iPad2 = iPad1+1; iPad2 < mpDEs[iDE].nPads[iPlane]; ++iPad2) {
      
      if (AreOverlapping(mpDEs[iDE].pads[iPlane][iPad1].area, mpDEs[iDE].pads[iPlane][iPad2].area, 1.e-4)) {
        
        mpDEs[iDE].pads[iPlane][iPad1].neighbours2[mpDEs[iDE].pads[iPlane][iPad1].nNeighbours2++] = iPad2;
        mpDEs[iDE].pads[iPlane][iPad2].neighbours2[mpDEs[iDE].pads[iPlane][iPad2].nNeighbours2++] = iPad1;
        
      }
      
    }
    
  }
  
}

//------------------------------------------------------------------
Bool_t AreOverlapping(Float_t area1[2][2], Float_t area2[2][2], Float_t precision)
{
  /// check if the two areas overlap
  /// precision in cm: positive = increase pad size / negative = decrease pad size
  
  return (!((area1[0][0] - area2[0][1] > precision) || (area2[0][0] - area1[0][1] > precision) ||
            (area1[1][0] - area2[1][1] > precision) || (area2[1][0] - area1[1][1] > precision)));
  
}

//------------------------------------------------------------------
void DrawDifferences()
{
  /// draw the neighbours found to be different between the 2 methods
  
  AliCodeTimerAutoGeneral("",0);
  
  AliMUON2DMap digitStore1(kTRUE);
  AliMUON2DMap digitStore2(kTRUE);
  
  // loop over DEs
  for (Int_t iDE = 1; iDE <= nDEs; iDE++) {
    
    // loop over planes
    for (Int_t iPlane = 0; iPlane < 2; iPlane++) {
      
      // loop over pads
      for (Int_t iPad = 0; iPad < mpDEs[iDE].nPads[iPlane]; ++iPad) {
        
        Bool_t sameNeighbours = (mpDEs[iDE].pads[iPlane][iPad].nNeighbours1 == mpDEs[iDE].pads[iPlane][iPad].nNeighbours2);
        
        if (sameNeighbours) {
          
          // loop over neighbours from method 1
          for (Int_t iNeighbour1 = 0; iNeighbour1 < mpDEs[iDE].pads[iPlane][iPad].nNeighbours1; iNeighbour1++) {
            
            // loop over neighbours from method 2
            Bool_t matchFound = kFALSE;
            for (Int_t iNeighbour2 = 0; iNeighbour2 < mpDEs[iDE].pads[iPlane][iPad].nNeighbours2; iNeighbour2++) {
              
              if (mpDEs[iDE].pads[iPlane][iPad].neighbours2[iNeighbour2] == mpDEs[iDE].pads[iPlane][iPad].neighbours1[iNeighbour1]) {
                
                matchFound = kTRUE;
                break;
                
              }
              
            }
            
            if (!matchFound) {
              
              sameNeighbours = kFALSE;
              break;
              
            }
            
          }
          
        }
        
        if (!sameNeighbours) {
          
          // register the pad with different neighbours
          Int_t manuId = AliMUONVDigit::ManuId(mpDEs[iDE].pads[iPlane][iPad].id);
          Int_t manuChannel = AliMUONVDigit::ManuChannel(mpDEs[iDE].pads[iPlane][iPad].id);
          
          AliMUONVCalibParam* c = static_cast<AliMUONVCalibParam*>(digitStore1.FindObject(mpDEs[iDE].id, manuId));
          if (!c) {
            c = new AliMUONCalibParamNI(1, AliMpConstants::ManuNofChannels(), mpDEs[iDE].id, manuId);
            digitStore1.Add(c);
          }
          c->SetValueAsInt(manuChannel, 0, c->ValueAsInt(manuChannel, 0) + 10);
          
          c = static_cast<AliMUONVCalibParam*>(digitStore2.FindObject(mpDEs[iDE].id, manuId));
          if (!c) {
            c = new AliMUONCalibParamNI(1, AliMpConstants::ManuNofChannels(), mpDEs[iDE].id, manuId);
            digitStore2.Add(c);
          }
          c->SetValueAsInt(manuChannel, 0, c->ValueAsInt(manuChannel, 0) + 10);
          
          // loop over neighbours from method 1
          for (Int_t iNeighbour1 = 0; iNeighbour1 < mpDEs[iDE].pads[iPlane][iPad].nNeighbours1; iNeighbour1++) {
            
            UInt_t padId = mpDEs[iDE].pads[iPlane][mpDEs[iDE].pads[iPlane][iPad].neighbours1[iNeighbour1]].id;
            manuId = AliMUONVDigit::ManuId(padId);
            manuChannel = AliMUONVDigit::ManuChannel(padId);
            
            // register the neighbour
            c = static_cast<AliMUONVCalibParam*>(digitStore1.FindObject(mpDEs[iDE].id, manuId));
            if (!c) {
              c = new AliMUONCalibParamNI(1, AliMpConstants::ManuNofChannels(), mpDEs[iDE].id, manuId);
              digitStore1.Add(c);
            }
            c->SetValueAsInt(manuChannel, 0, c->ValueAsInt(manuChannel, 0) + 1);
            
          }
          
          // loop over neighbours from method 2
          for (Int_t iNeighbour2 = 0; iNeighbour2 < mpDEs[iDE].pads[iPlane][iPad].nNeighbours2; iNeighbour2++) {
            
            UInt_t padId = mpDEs[iDE].pads[iPlane][mpDEs[iDE].pads[iPlane][iPad].neighbours2[iNeighbour2]].id;
            manuId = AliMUONVDigit::ManuId(padId);
            manuChannel = AliMUONVDigit::ManuChannel(padId);
            
            // register the neighbour
            c = static_cast<AliMUONVCalibParam*>(digitStore2.FindObject(mpDEs[iDE].id, manuId));
            if (!c) {
              c = new AliMUONCalibParamNI(1, AliMpConstants::ManuNofChannels(), mpDEs[iDE].id, manuId);
              digitStore2.Add(c);
            }
            c->SetValueAsInt(manuChannel, 0, c->ValueAsInt(manuChannel, 0) + 1);
            
          }
          
        }
        
      }
      
    }
    
  }
  
  // create the tracker data
  AliMUONTrackerData digitData1("neighbours1", "neighbours1", 1, kTRUE);
  digitData1.SetDimensionName(0, "index");
  digitData1.Add(digitStore1);
  AliMUONTrackerData digitData2("neighbours2", "neighbours2", 1, kTRUE);
  digitData2.SetDimensionName(0, "index");
  digitData2.Add(digitStore2);
  
  // save it to a file
  TFile *outFile = TFile::Open("neighbours.root", "UPDATE");
  if (outFile && outFile->IsOpen()) {
    digitData1.Write(0x0, TObject::kOverwrite);
    digitData2.Write(0x0, TObject::kOverwrite);
    outFile->Close();
  }
  
}

