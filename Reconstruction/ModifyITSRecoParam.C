/*
 *  ModifyITSRecoParam.C
 *  aliroot_dev
 *
 *  Created by philippe pillot on 10/10/13.
 *  Copyright 2013 SUBATECH. All rights reserved.
 *
 */

void ModifyITSRecoParam(TString ocdbIn, Int_t run, TString ocdbOut = "")
{
  /// Get the ITS recoParam from the default storage, modify them and save them at the given location
  
  AliCDBManager* man = AliCDBManager::Instance();
  man->SetDefaultStorage(ocdbIn.Data());
  man->SetRun(run);
  
  AliCDBEntry* entry = AliCDBManager::Instance()->Get("ITS/Calib/RecoParam");
  if (!entry) return;
  
  TObject* o = entry->GetObject();
  if ( o->IsA() == TObjArray::Class() ) {
    
    TObjArray* array = static_cast<TObjArray*>(o);
    for ( Int_t i = 0; i <= array->GetLast(); ++i ) {
      
      AliITSRecoParam* p = static_cast<AliITSRecoParam*>(array->At(i));
      
      // modify parameter(s) here
      p->ReconstructOnlySPD();
      
    }
    
  } else {
    
    AliITSRecoParam* p = static_cast<AliITSRecoParam*>(o);
    
    // modify parameter(s) here
    p->ReconstructOnlySPD();
    
  }
  
  if (!ocdbOut.IsNull()) man->SetDefaultStorage(ocdbOut.Data());
  AliMUONCDB::WriteToCDB(o, "ITS/Calib/RecoParam", run, run, "Reconstruction parameters ITS", "Andrea Dainese");
  
}
