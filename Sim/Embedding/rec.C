// Start with specific muon settings

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliCDBManager.h"
#include "AliReconstruction.h"
#include <TGrid.h>
#include <TSystem.h>
#include "AliITSRecoParam.h"
#endif

AliITSRecoParam* GetSpecialITSRecoParam()
{
AliITSRecoParam *itsRecoParam = AliITSRecoParam::GetHighFluxParam();
itsRecoParam->SetTrackleterPhiWindowL2(0.07);
itsRecoParam->SetTrackleterZetaWindowL2(0.4);
itsRecoParam->SetTrackleterPhiWindowL1(0.10);
itsRecoParam->SetTrackleterZetaWindowL1(0.6);
itsRecoParam->SetTrackleterPhiWindow(0.06);
itsRecoParam->SetTrackleterThetaWindow(0.025);
itsRecoParam->SetTrackleterScaleDThetaBySin2T(kTRUE);
itsRecoParam->SetTrackleterRemoveClustersFromOverlaps(kTRUE);

itsRecoParam->SetVertexerZ();  // this is the present default. We plan to use 3D vertexer for pass 2. If you want to select 3D you can replace this line with
itsRecoParam->ReconstructOnlySPD();

return itsRecoParam;
}

void rec(const char *filename="raw.root")
{
  /////////////////////////////////////////////////////////////////////////////////////////
  //
  // Reconstruction script for 2011 RAW data - muon fast reco
  //
  /////////////////////////////////////////////////////////////////////////////////////////



  // Set the CDB storage location
  AliCDBManager * man = AliCDBManager::Instance();
  //man->SetDefaultStorage("alien://folder=/alice/data/2011/OCDB");
  man->SetDefaultStorage("local:///Users/pillot/Work/Alice/Data/2012/LHC12h/raw/190147/saf/CDBMirror/alice/data/2012/OCDB");
//  man->SetSpecificStorage("MUON/Align/Data","local://copy2011pass1CDB");
  man->SetSpecificStorage("ITS/Calib/RecoParam","local:///Users/pillot/Work/Alice/Data/2012/LHC12h/raw/190147/OCDB_RecoParams");
  man->SetSpecificStorage("MUON/Calib/RecoParam","local:///Users/pillot/Work/Alice/Data/2012/LHC12h/raw/190147/OCDB_newParams");

  AliReconstruction rec;

  // Generate or use the local OCDB.root file
	//  rec.SetFromCDBSnapshot("OCDB.root");

  // Set reconstruction flags (skip detectors here if neded with -<detector name>
	//  rec.SetRunReconstruction("MUON ITS VZERO ZDC T0");
	rec.SetRunReconstruction("MUON ITS");

  // QA options
//  rec.SetRunQA("Global MUON:ALL") ;
  rec.SetRunQA("MUON:ALL");
  rec.SetQAWriteExpert(AliQAv1::kMUON);
//	rec.SetRunQA(":") ;
  rec.SetRunGlobalQA(kFALSE);
  rec.SetQARefDefaultStorage("local://$ALICE_ROOT/QAref") ;

  // MUON only reco - recoparameters
  //rec.SetRecoParam("ITS",GetSpecialITSRecoParam());

  // AliReconstruction settings
  rec.SetWriteESDfriend(kFALSE);
  //rec.SetWriteAlignmentData();
//  rec.SetInput(filename);
  //rec.SetUseTrackingErrorsForAlignment("ITS");

  // switch off cleanESD
  rec.SetCleanESD(kFALSE);

  // Specific reco params for ZDC (why isn't this automatic?)
	//  rec.SetRecoParam("ZDC",AliZDCRecoParamPbPb::GetHighFluxParam(2760));

  //Ignore SetStopOnError
  rec.SetStopOnError(kFALSE);

	// Ignore SetStopOnMissingTriggerFile
	rec.SetStopOnMissingTriggerFile(kFALSE);
	
  AliLog::Flush();
  rec.Run();

}

