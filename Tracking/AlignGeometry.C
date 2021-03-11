#include <string>

#include <TError.h>
#include <TSystem.h>
#include <TGeoManager.h>

#include "AliCDBManager.h"
#include "AliGeomManager.h"

using namespace std;

void AlignGeometry(int run, string ocdb = "OCDB.root")
{
  /// load geometry, apply alignment and save the new aligned geometry
  /// parameter "ocdb" can a snapshot (e.g. OCDB.root) or a path (e.g. raw:// or local://./OCDB)

  // set OCDB location
  AliCDBManager* man = AliCDBManager::Instance();
  if (ocdb.rfind(".root") != string::npos) {
    if (gSystem->AccessPathName(ocdb.c_str(), kFileExists) == 0) {
      man->SetDefaultStorage("local:///dev/null");
      man->SetSnapshotMode(ocdb.c_str());
    } else {
      Error("AlignGeometry", "snapshot %s not found", ocdb.c_str());
      return;
    }
  } else if (ocdb == "raw://" || ocdb.find("alien://") != string::npos || ocdb.find("local://") != string::npos) {
    man->SetDefaultStorage(ocdb.c_str());
  } else {
    Error("AlignGeometry", "%s is not a valid OCDB path", ocdb.c_str());
    return;
  }
  man->SetRun(run);

  // load the geometry and apply alignment
  AliGeomManager::LoadGeometry();
  if (!AliGeomManager::GetGeometry() || !AliGeomManager::ApplyAlignObjsFromCDB("MUON")) {
    Error("AlignGeometry", "unable to load geometry and apply alignment");
    return;
  }

  // save the new aligned geometry
  AliGeomManager::GetGeometry()->Export("AlignedGeometry.root");
}
