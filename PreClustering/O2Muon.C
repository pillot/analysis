#include "O2Muon.h"
#include "AliRawReader.h"
#include <iostream>
#include "AliRawEventHeaderBase.h"
#include "AliMUONCalibrationData.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBEntry.h"
#include "AliMUONRecoParam.h"
#include "AliMUONDigitCalibrator.h"
#include "AliMUONVDigitStore.h"
#include "AliMUONDigitMaker.h"
#include "TFile.h"
#include "TTree.h"
#include "AliMUONGeometryTransformer.h"
#include "AliMUONSimpleClusterServer.h"
#include "AliMUONVClusterFinder.h"
#include "AliMUONPreClusterFinder.h"
#include "AliMUONClusterStoreV2.h"
#include "AliMpArea.h"
#include "AliMpDDLStore.h"
#include "AliMpSegmentation.h"
#include "AliMpCDB.h"
#include "AliLog.h"
#include "AliGeomManager.h"
#include "TSystem.h"
#include "AliCodeTimer.h"
#include "AliMUONClusterFinderMLEM.h"
#include "AliMUONReconstructor.h"
#include <fstream>

O2Muon::O2Muon(const char* ocdbPath) : TObject(), mOCDBPath(ocdbPath)
{
  
}

O2Muon::~O2Muon()
{
  
}

AliMUONRecoParam* O2Muon::getRecoParam(int runNumber)
{
  AliMUONRecoParam* recoParam(0x0);
  
  AliCDBEntry* e = AliCDBManager::Instance()->Get("MUON/Calib/RecoParam",runNumber);
  if (e)
  {
    TObject* o = e->GetObject();
    if ( o->IsA() == TObjArray::Class() )
    {
      TIter next(static_cast<TObjArray*>(o));
      AliMUONRecoParam* p;
      while ( ( p = static_cast<AliMUONRecoParam*>(next()) ))
      {
        if ( p->IsDefault()) recoParam = p;
      }
    }
    else
    {
      recoParam = static_cast<AliMUONRecoParam*>(o);
    }
  }
  return recoParam;
}


int O2Muon::makeDigitFile(const char* rawDataInputFile, const char* digitOutputFile, Bool_t calibrate)
{
  /** Create a Root file with calibrated MCH digits from a raw data file.
    * Not meant to be fast, just a re-use of the existing classes to get the job done.
    *
    * @param in rawDataInputFile rawdata input file used as a source of raw data
    * @param out digitOutputFile output file with digits
    * @param calibrate : whether or not the digits are calibrated
    */
  
  AliRawReader* rawReader = AliRawReader::Create(rawDataInputFile);
  
  Bool_t ok = rawReader->NextEvent();
  
  if (!ok) return -1;
  
  int runNumber = rawReader->GetRunNumber();

  if (runNumber<0)
  {
    return -2;
  }
  
  rawReader->RewindEvents();

  std::cout << "RUN=" << runNumber << std::endl;

  prepareOCDB(runNumber);
  
  AliMUONRecoParam* recoParam = getRecoParam(runNumber);

  AliMUONDigitMaker digitMaker(kFALSE);
  
  AliMUONCalibrationData* calibrationData(0x0);
  AliMUONDigitCalibrator* digitCalibrator(0x0);
  
  if ( calibrate )
  {
    calibrationData = new AliMUONCalibrationData(runNumber);
    digitCalibrator  = new AliMUONDigitCalibrator(*calibrationData,recoParam);
  }
  
  AliMUONVDigitStore* digitStore = AliMUONVDigitStore::Create("AliMUONDigitStoreV2R");
  
  if (!digitStore)
  {
    return -3;
  }
  
  TFile* output = new TFile(digitOutputFile,"RECREATE");
  
  TTree* treeD = new TTree("TreeD","Digits");

  ok = digitStore->Connect(*treeD,true);
  
  if (!ok)
  {
    return -4;
  }
  
  int nofEvents = 0;
  int nofPhysicsEvents = 0;
  int meanNofDigits = 0;
  
  while ( rawReader->NextEvent() )
  {
    ++nofEvents;

    Int_t eventType = rawReader->GetType();
    
    if (eventType != AliRawEventHeaderBase::kPhysicsEvent )
    {
      continue;
    }

    ++nofPhysicsEvents;

    digitMaker.Raw2Digits(rawReader,digitStore);

    if ( calibrate )
    {
      digitCalibrator->Calibrate(*digitStore);
    }
    
    treeD->Fill();
    
    meanNofDigits += digitStore->GetSize();
    
    digitStore->Clear();
  }
  
  output->cd();
  treeD->Write();
  output->Close();
  delete output;
  
  if (nofPhysicsEvents>0)
  {
    meanNofDigits /= nofPhysicsEvents;
  }
  
  std::cout << "nofEvents=" << nofEvents << " nofPhysicsEvents=" << nofPhysicsEvents << " < nofDigit per event >=" << meanNofDigits << std::endl;
  
  delete digitCalibrator;
  delete calibrationData;
  
  return 0;
}

int O2Muon::makeClustering(const char* digitInputFile, const char* clusterFinderType, const char* outputLogFile, int runNumber)
{
  /** Make pre-clusters out of digits. Calibrated or not, that's not really relevant for the pre-clustering,
   *  except for the fact that in the calibrated case some pads will have a charge of zero, and so should not
   *  be used (the pads are not removed from the digitstore at calibration stage to avoid reallocating 
   *  the digitstore)
   */
  
  prepareOCDB(runNumber);
  
  TFile* digitFile = TFile::Open(digitInputFile);
  
  if (!digitFile)
  {
    return -1;
  }
  
  TTree* treeD = static_cast<TTree*>(digitFile->Get("TreeD"));
  
  if (!treeD)
  {
    return -2;
  }
  
  AliMUONVDigitStore* digitStore = AliMUONVDigitStore::Create(*treeD);
  
  digitStore->Connect(*treeD);

  AliMUONVClusterFinder* clusterFinder = AliMUONReconstructor::CreateClusterFinder(clusterFinderType);

  AliGeomManager::LoadGeometry(Form("%s/test/QA/geometry.root",gSystem->Getenv("$ALICE_ROOT")));
  
  AliMUONGeometryTransformer transformer;
  transformer.LoadGeometryData();
  
  AliMUONSimpleClusterServer clusterServer(clusterFinder,transformer);
  
  AliMUONVClusterStore* clusterStore = new AliMUONClusterStoreV2();
  AliMpArea area; // invalid area to clusterize everything
  
  AliMUONRecoParam* recoParam = getRecoParam(runNumber);

  AliCodeTimerAuto("total",0);
  Int_t nofClusters(0);
  Long64_t nofEntries = treeD->GetEntries();
  
  for ( Long64_t i = 0; i < nofEntries; ++i )
  {
    treeD->GetEntry(i);
    
    TIter next(digitStore->CreateIterator());
    clusterServer.UseDigits(next,digitStore);
    
    for ( int chamberId = 0; chamberId < 10; ++chamberId )
    {
      clusterServer.Clusterize(chamberId,*clusterStore,area,recoParam);
    }

    nofClusters += clusterStore->GetSize();
    clusterStore->Clear();
  }
  
  std::ofstream out(outputLogFile);
  std::streambuf* coutbuf = std::cout.rdbuf(); //save old buf
  std::cout.rdbuf(out.rdbuf()); //redirect std::cout to out file

  AliCodeTimer::Instance()->Print();

  std::cout << "Mean number of clusters per event : " << (nofEntries > 0 ? nofClusters*1.0/nofEntries : 0) << std::endl;
  
  std::cout.rdbuf(coutbuf); //reset to standard output
  
  delete clusterStore;
  delete digitFile;
  
  return 0;
}

void O2Muon::prepareOCDB(int runNumber)
{
  if ( mOCDBPath.size() > 0 )
  {
    AliCDBStorage* storage = AliCDBManager::Instance()->GetDefaultStorage();
    
    if ( ( storage && ( storage->GetURI() != mOCDBPath.c_str() ) ) || (!storage) )
    {
      AliCDBManager::Instance()->SetDefaultStorage(mOCDBPath.c_str());
    }
    
    AliCDBManager::Instance()->SetRun(runNumber);
    
    if ( AliMpDDLStore::Instance(false) ) {
      AliCDBManager::Instance()->UnloadFromCache("MUON/Calib/DDLStore");
      delete AliMpDDLStore::Instance();
    }
    
    if ( AliMpSegmentation::Instance(false) ) {
      AliCDBManager::Instance()->UnloadFromCache("MUON/Calib/Mapping");
      delete AliMpSegmentation::Instance();
    }
    
    // Load mapping
    if ( ! AliMpCDB::LoadDDLStore() ) {
      AliFatal("Could not access mapping from OCDB !");
    }
    
  }
  
}


