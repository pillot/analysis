/*
 *  ChangeRejectList
 *
 *  Get a rejectlist from OCDB and "patch" it
 *
 *  Author: Laurent Aphecetche
 *
 */

#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliMUONRejectList.h"
#include "AliMUONCDB.h"
#include "AliCDBRunRange.h"
#include "AliMpCDB.h"

//______________________________________________________________________________
void ChangeRejectList(Int_t inputRun=266437,
                      const char* inputOCDB="alien://folder=/alice/data/2016/OCDB",
                      const char* outputOCDB="alien://folder=/alice/cern.ch/user/p/ppillot/OCDB2016Test")
{
  AliCDBManager::Instance()->SetDefaultStorage(inputOCDB);
  AliCDBManager::Instance()->SetRun(inputRun);
  
  AliMpCDB::LoadAll();

  AliCDBEntry* entry = AliCDBManager::Instance()->Get("MUON/Calib/RejectList");

  AliMUONRejectList* rl = static_cast<AliMUONRejectList*>(entry->GetObject());

  // repaired one
//  rl->SetHVProbability("MchHvLvRight/Chamber03Right/Quad4Sect1",0.0);

  // alignment ok
//  rl->SetDetectionElementProbability(806,0.0);
  
  // new hole
//  rl->SetPCBProbability(1023,3);
  rl->SetBusPatchProbability(540, 1.);
  rl->SetBusPatchProbability(528, 1.);
  rl->SetManuProbability(303,    8, 1.);
  rl->SetManuProbability(303,   34, 1.);
  rl->SetManuProbability(303,   60, 1.);
  rl->SetManuProbability(303,   85, 1.);
  rl->SetManuProbability(303,  110, 1.);
  rl->SetManuProbability(303,  136, 1.);
  rl->SetManuProbability(303,  162, 1.);
  rl->SetManuProbability(303,  187, 1.);
  rl->SetManuProbability(303,  203, 1.);
  rl->SetManuProbability(303,  216, 1.);
  rl->SetManuProbability(303,  228, 1.);
  rl->SetManuProbability(303, 1032, 1.);
  rl->SetManuProbability(303, 1058, 1.);
  rl->SetManuProbability(303, 1084, 1.);
  rl->SetManuProbability(303, 1109, 1.);
  rl->SetManuProbability(303, 1134, 1.);
  rl->SetManuProbability(303, 1160, 1.);
  rl->SetManuProbability(303, 1186, 1.);
  rl->SetManuProbability(303, 1211, 1.);
  rl->SetManuProbability(303, 1227, 1.);
  rl->SetManuProbability(303, 1240, 1.);
  rl->SetManuProbability(303, 1253, 1.);
  
  AliCDBManager::Instance()->SetDefaultStorage(outputOCDB);
  
  AliMUONCDB::WriteToCDB(rl, "MUON/Calib/RejectList", 0, AliCDBRunRange::Infinity(),
                         "reject list for MUON, updated for 2015, from comparison of cluster maps between reconstructed data and simulations", "L. Aphecetche and P. Pillot");
}

