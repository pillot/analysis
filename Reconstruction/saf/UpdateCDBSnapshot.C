#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TSystem.h>
#include <TMap.h>
#include <TObjArray.h>
#include <TString.h>
#include <TIterator.h>
#include <TMath.h>
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include <AliCDBId.h>
#endif

// will load only objects specific for this run (~55/202) provided  the mirror is up to data


void UpdateCDBSnapshot(int run, const char* mirrorBase);
const char* PreloadFromStorage(AliCDBStorage *storage,const char* mirrorBase, Bool_t isDef);

void UpdateCDBSnapshot(int run, const char* mirrorBase)
{
  AliCDBManager* man = AliCDBManager::Instance();
  if (!man->IsDefaultStorageSet()) {
    printf("Default storage is not set, setting to raw\n"); 
    man->SetDefaultStorage("raw://");
  }
  if (man->GetRun() < 0) {
    printf("Run number is not set, setting to %d\n", run); 
    man->SetRun(run);
  }
  //
  TString locDefFolder="";
  AliCDBStorage *defStorage = man->GetDefaultStorage();
  printf("Preloading objects from default storage\n");
  if (defStorage->GetType()=="alien") {
    defStorage->QueryCDB(run);
    locDefFolder = PreloadFromStorage(defStorage,mirrorBase,kTRUE);
  }
  else printf("... Not needed: default storage is not on alien.\n");
  TObjArray specArr;
  specArr.SetOwner(kTRUE);
  //
  // check if specific storages are involved
  const TMap* stMap = man->GetStorageMap();
  TIter nextSt(stMap);
  TObjString *str;
  printf("Preloading objects from specific storages\n");
  while ((str=(TObjString*)nextSt())) {
    TString calType = str->GetString();
    if (calType=="default") continue;
    AliCDBStorage* storage = man->GetSpecificStorage(calType.Data());    
    if (!storage || storage->GetType()!="alien") continue;
    storage->QueryCDB(run,calType.Data());
    const char* newPath = PreloadFromStorage(storage,mirrorBase,kFALSE);
    specArr.AddLast(new TNamed(calType.Data(),newPath));
    //
  }
  //
  // redefine defauld storage to local one
  if (!locDefFolder.IsNull()) {
    printf("Redefining default storage to %s\n",locDefFolder.Data());
    man->UnsetDefaultStorage(); 
    man->SetDefaultStorage(locDefFolder.Data());
  }
  //
  // redefine specific storages to local ones
  TIter nextSp(&specArr);
  TNamed* specPath = 0;
  printf("Redefining specific storages to local snapshot\n");
  while ((specPath=(TNamed*)nextSp())) man->SetSpecificStorage(specPath->GetName(),specPath->GetTitle());
  //
}

const char* PreloadFromStorage(AliCDBStorage *storage,const char* mirrorBase, Bool_t isDef)
{
  // check if objects to be used are already preloaded locally, if not, fetch them from alien CDB
  //
  static TString drainFolder;
  drainFolder = Form("local://%s/%s",mirrorBase, storage->GetBaseFolder().Data());
  TObjArray* arrCDBID = storage->GetQueryCDBList();
  TIter nxt(arrCDBID);
  AliCDBId* cdbID = 0;
  AliCDBManager* man = AliCDBManager::Instance();
  const TMap* stMap = man->GetStorageMap();
  man->SetDrain(drainFolder.Data());
  while ((cdbID=(AliCDBId*)nxt())) {
    if (isDef && stMap->GetValue(cdbID->GetPath())) {
      printf("Entry %s will be taken from specific storage, skipping\n",cdbID->GetPath().Data());
      continue; // defined in the specific storage
    }
    TString locPath = Form("%s/%s%s/Run%d_%d_v%d_s%d.root",mirrorBase, 
			   storage->GetBaseFolder().Data(),cdbID->GetPath().Data(),
			   cdbID->GetFirstRun(),cdbID->GetLastRun(),cdbID->GetVersion(),
			   TMath::Max(0,cdbID->GetSubVersion()));
    printf("Check for %s\n",locPath.Data());
    if (gSystem->AccessPathName(locPath.Data())) man->Get(cdbID->GetPath(),man->GetRun(),cdbID->GetVersion());
    else printf("found local copy %s\n",locPath.Data());
  }
  man->UnsetDrain();
  //
  return drainFolder.Data();
}
