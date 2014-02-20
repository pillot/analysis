/*
 *  ModifyMUONRecoParam.C
 *  aliroot_dev
 *
 *  Created by philippe pillot on 14/10/13.
 *  Copyright 2013 SUBATECH. All rights reserved.
 *
 */

void ModifyMUONRecoParam(TString ocdbIn,
			 Int_t startRun = 0,
			 Int_t endRun = AliCDBRunRange::Infinity(),
			 TString ocdbOut = "")
{
  /// Get the MUON recoParam from the default storage, modify them and save them at the given location
  
  AliCDBManager* man = AliCDBManager::Instance();
  man->SetDefaultStorage(ocdbIn.Data());
  man->SetRun(startRun);
  
  AliCDBEntry* entry = AliCDBManager::Instance()->Get("MUON/Calib/RecoParam");
  if (!entry) return;
  
  TObject* o = entry->GetObject();
  if ( o->IsA() == TObjArray::Class() ) {
    
    TObjArray* array = static_cast<TObjArray*>(o);
    for ( Int_t i = 0; i <= array->GetLast(); ++i ) {
      
      AliMUONRecoParam* p = static_cast<AliMUONRecoParam*>(array->At(i));
      
      // modify parameter(s) here
      if (p->GetEventSpecie() != 5) continue;
      //p->MakeMoreTrackCandidates(kTRUE);
      //p->SetMaxTrackCandidates(50000);
      //p->SetManuOccupancyLimits(-1.,0.015);
      //p->SetPadGoodnessMask(0x400be00);
      p->DiscardMonoCathodClusters(kTRUE, 10., 10.);
      
    }
    
  } else {
    
    AliMUONRecoParam* p = static_cast<AliMUONRecoParam*>(o);
    
    // modify parameter(s) here
    if (p->GetEventSpecie() != 5) continue;
    //p->MakeMoreTrackCandidates(kTRUE);
    //p->SetMaxTrackCandidates(50000);
    //p->SetManuOccupancyLimits(-1.,0.015);
    //p->SetPadGoodnessMask(0x400be00);
    p->DiscardMonoCathodClusters(kTRUE, 10., 10.);
    
  }
  
  if (!ocdbOut.IsNull()) man->SetDefaultStorage(ocdbOut.Data());
  AliMUONCDB::WriteToCDB(o, "MUON/Calib/RecoParam", startRun, endRun,
			 "reconstruction parameters for MUON", "L. Aphecetche and P. Pillot");
  
}
