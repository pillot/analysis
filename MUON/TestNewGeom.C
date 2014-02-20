
TestNewGeom(){
    // get original geometry transformer
    AliCDBManager* cdbm = AliCDBManager::Instance();
	TString defaultStorage("local://$ALICE_ROOT/OCDB");
	cdbm->SetDefaultStorage(defaultStorage.Data());
	cdbm->SetRun(119841);
    AliGeomManager::LoadGeometry();
    if (!AliGeomManager::GetGeometry()) return;


	TString fOldAlignStorage("local:///Users/philippe/Work/Alice/Work/Data/pp7TeV/LHC10c/PositiveField/run119841/pass2_GMS_realign_R/OCDB");
	TString fNewAlignStorage("local:///Users/philippe/Work/Alice/Work/aliroot/SHUTTLE/TestShuttle/TestCDB");
	
    if (fOldAlignStorage != "none") {
      if (!fOldAlignStorage.IsNull()) cdbm->SetSpecificStorage("MUON/Align/Data",fOldAlignStorage.Data());
      else cdbm->SetSpecificStorage("MUON/Align/Data",defaultStorage.Data());
      AliGeomManager::ApplyAlignObjsFromCDB("MUON");
    }
    AliMUONGeometryTransformer * fOldGeoTransformer = new AliMUONGeometryTransformer();
    fOldGeoTransformer->LoadGeometryData();
    cout << "Old" << endl;
	for(int i=0;i<16;i++){
		cout<<" iCh = "<<i<<endl;
    	fOldGeoTransformer->GetModuleTransformer(i)->GetTransformation()->Print();
    }
	//return;
	
    // get new geometry transformer
    cdbm->UnloadFromCache("GRP/Geometry/Data");
   if (fOldAlignStorage != "none") cdbm->UnloadFromCache("MUON/Align/Data");
    AliGeomManager::GetGeometry()->UnlockGeometry();
    AliGeomManager::LoadGeometry();
    if (!AliGeomManager::GetGeometry()) return;
    if (!fNewAlignStorage.IsNull()) cdbm->SetSpecificStorage("MUON/Align/Data",fNewAlignStorage.Data());
    else cdbm->SetSpecificStorage("MUON/Align/Data",defaultStorage.Data());
    AliGeomManager::ApplyAlignObjsFromCDB("MUON");
    AliMUONGeometryTransformer *  fNewGeoTransformer = new AliMUONGeometryTransformer();
    fNewGeoTransformer->LoadGeometryData();
    cout << "New" << endl;
	for(int i=0;i<16;i++){
		cout<<" iCh = "<<i<<endl;
    	fNewGeoTransformer->GetModuleTransformer(i)->GetTransformation()->Print();
}
}