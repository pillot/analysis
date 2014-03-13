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
#include "AliMUONPolygon.h"
#include "AliMUONContourMaker.h"
#include "AliMUONContour.h"
#include "AliMpPad.h"
#include "AliMpVSegmentation.h"
#include "AliMUONContourMakerTest.h"
#include "TCanvas.h"
#include "AliTriggerConfiguration.h"
#include "AliTriggerClass.h"
#include "AliMUONTrackStoreV1.h"
#include "TGeoGlobalMagField.h"
#include "AliGRPObject.h"
#include "AliMagF.h"
#include "AliMUONTrackReconstructorK.h"
#include "AliMUONRawStreamTrackerHP.h"

//_________________________________________________________________________________________________
O2Muon::O2Muon(const char* ocdbPath) : TObject(), mOCDBPath(ocdbPath) {
  
}

//_________________________________________________________________________________________________
O2Muon::~O2Muon() {
  
}

//_________________________________________________________________________________________________
int O2Muon::getRunNumber(AliRawReader& rawReader) const {
  /** Get run number from rawreader. 
   WARNING : rawReader is rewinded during the operation !
   */
   
  Bool_t ok = rawReader.NextEvent();
  
  if (!ok) return -1;
  
  int runNumber = rawReader.GetRunNumber();
  
  if (runNumber<0)
  {
    return -2;
  }

  rawReader.RewindEvents();

  return runNumber;
}

//_________________________________________________________________________________________________
AliMUONRecoParam* O2Muon::getRecoParam(int runNumber) {
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

//_________________________________________________________________________________________________
int O2Muon::filterRaw(const char* rawDataInputFile, const char* triggerClass) {

  TString input(rawDataInputFile);
  
  if (strlen(triggerClass)>0) {
    input += "?Trigger=";
    input += triggerClass;
  }

  AliRawReader* rawReader = AliRawReader::Create(input.Data());

  if ( !rawReader ) {
    return -1;
  }
  
  Int_t runNumber = getRunNumber(*rawReader);
  
  if ( runNumber < 0 ) {
    return -2;
  }
  
  prepareOCDB(runNumber,rawReader);

  int n(0);
  
  TFile* inputFile = TFile::Open(rawDataInputFile);
  TTree* inputTree = static_cast<TTree*>(inputFile->Get("RAW"));
  
  TString outputRawFile(gSystem->BaseName(rawDataInputFile));
  
  if (strlen(triggerClass)>0) {
    outputRawFile.ReplaceAll(".root",Form(".%s.root",triggerClass));
  }

  TFile* outputFile = TFile::Open(outputRawFile.Data(),"recreate");
  TTree* outputTree = inputTree->CloneTree(0);
  
  while ( rawReader->NextEvent() ) {
    ++n;
    inputTree->CopyAddresses(outputTree);
    outputTree->Fill();
  }

  outputTree->AutoSave();
  delete outputFile;
  delete inputFile;
  delete rawReader;
  
  std::cout << "n=" << n << std::endl;

  return 0;
}

//_________________________________________________________________________________________________
int O2Muon::decodeEvents(const char* rawDataInputFile) {
  
  AliRawReader* rawReader = AliRawReader::Create(rawDataInputFile);
  Long64_t nofEvents(0);

  {
    AliCodeTimerAuto("total",0);
    
    AliMUONRawStreamTrackerHP stream(rawReader);
    
    while ( rawReader->NextEvent() )  {
      ++nofEvents;
      
      Int_t busPatch;
      UShort_t manuId, adc;
      UChar_t manuChannel;
      
      AliCodeTimerAutoGeneral("decode",1);

      stream.First();
      
      while ( stream.Next(busPatch,manuId,manuChannel,adc) )
      {
        adc *= 2;
      }
      
    }
  }
  
  delete rawReader;
  
//  std::ofstream out(outputLogFile);
//  std::streambuf* coutbuf = std::cout.rdbuf(); //save old buf
//  std::cout.rdbuf(out.rdbuf()); //redirect std::cout to out file
  
  AliCodeTimer::Instance()->Print();
  
  std::cout << "Number of events : " << nofEvents << std::endl;
//  
//  std::cout.rdbuf(coutbuf); //reset to standard output
//  

  return 0;
}

//_________________________________________________________________________________________________
int O2Muon::makeDigitFile(const char* rawDataInputFile, const char* digitOutputFile,
                          Bool_t calibrate) {
  /** Create a Root file with calibrated MCH digits from a raw data file.
    * Not meant to be fast, just a re-use of the existing classes to get the job done.
    *
    * @param in rawDataInputFile rawdata input file used as a source of raw data
    * @param out digitOutputFile output file with digits
    * @param calibrate : whether or not the digits are calibrated
    */
  
  AliRawReader* rawReader = AliRawReader::Create(rawDataInputFile);
  
  int runNumber = getRunNumber(*rawReader);
  
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
  
  TObjArray *fields = TString(rawDataInputFile).Tokenize("?");
  TString sDigitOutputFile(digitOutputFile);
  
  if (fields->GetEntries()>1) {
    sDigitOutputFile.ReplaceAll(".root",Form(".%s.root",static_cast<TObjString*>(fields->At(1))->String().Data()));
  }
  
  delete fields;
  
  TFile* output = new TFile(sDigitOutputFile.Data(),"RECREATE");
  
  TTree* treeD = new TTree("TreeD","Digits");

  Bool_t ok = digitStore->Connect(*treeD,true);
  
  if (!ok)
  {
    return -4;
  }
  
  int nofEvents = 0;
  int nofPhysicsEvents = 0;
  int meanNofDigits = 0;
  
  AliCodeTimerAuto("total",0);
  
  while ( rawReader->NextEvent() ) //&& nofEvents < 10 )
  {
    ++nofEvents;

    Int_t eventType = rawReader->GetType();
    
    if (eventType != AliRawEventHeaderBase::kPhysicsEvent )
    {
      continue;
    }

//    AliInfo(Form("Event %d classMask %lld",nofEvents,rawReader->GetClassMask()));
    
    ++nofPhysicsEvents;

    {
      AliCodeTimerAuto("Raw2Digits",1);
      digitMaker.Raw2Digits(rawReader,digitStore);
    }

    if ( calibrate )
    {
      AliCodeTimerAuto("Calibrate",1);
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
  delete digitStore;
  delete rawReader;
  
  return 0;
}

//_________________________________________________________________________________________________
int O2Muon::makeDigitFiles(const char* rawDataInputFileList, const char* triggerClass,
                           Bool_t calibrate) {
  std::ifstream in(rawDataInputFileList);
  char line[1024];
  
  while (in.getline(line,1024,'\n'))
  {
    TString input(line);
    
    
    TString filename(gSystem->BaseName(line));
    filename.ReplaceAll(".root",".digits.root");
    
    if (strlen(triggerClass)>0) {
      input += "?Trigger=";
      input += triggerClass;
    }
      
    std::cout << input << " -> " << filename.Data() << std::endl;
    
    makeDigitFile(input,filename.Data(),calibrate);
  }
  return 0;
}

//_________________________________________________________________________________________________
int O2Muon::makeClustering(const char* digitInputFile, const char* clusterOutputFile,
                           const char* clusterFinderType, const char* outputLogFile,
                           int runNumber) {
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
  
  TFile* clusterFile = TFile::Open(clusterOutputFile,"RECREATE");
  TTree* treeR = new TTree("TreeR","Cluster Container");
  clusterStore->Connect(*treeR,kTRUE);
  
  AliMpArea area; // invalid area to clusterize everything
  
  AliMUONRecoParam* recoParam = getRecoParam(runNumber);

  AliCodeTimerAuto("total",0);
  Int_t nofClusters(0);
  Long64_t nofEntries = treeD->GetEntries();
  
//  std::cout << "Ready to start. Hit a key" << std::endl;
//  Int_t key;
//  std::cin >> key;
  
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
    
    treeR->Fill();
    
    clusterStore->Clear();
    digitStore->Clear();
  }
  
  clusterFile->cd();
  treeR->Write();
  delete clusterFile;
  
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

//_________________________________________________________________________________________________
int O2Muon::setupMagneticField() {
  
  AliCDBEntry* grpEntry = AliCDBManager::Instance()->Get("GRP/GRP/Data");
  
  AliGRPObject* grpData = dynamic_cast<AliGRPObject*>(grpEntry->GetObject());

  if ( !TGeoGlobalMagField::Instance()->IsLocked() ) {
    // Construct the field map out of the information retrieved from GRP.
    Bool_t ok = kTRUE;
    // L3
    Float_t l3Current = grpData->GetL3Current((AliGRPObject::Stats)0);
    if (l3Current == AliGRPObject::GetInvalidFloat()) {
      AliError("GRP/GRP/Data entry:  missing value for the L3 current !");
      ok = kFALSE;
    }
    
    Char_t l3Polarity = grpData->GetL3Polarity();
    if (l3Polarity == AliGRPObject::GetInvalidChar()) {
      AliError("GRP/GRP/Data entry:  missing value for the L3 polarity !");
      ok = kFALSE;
    }
    
    // Dipole
    Float_t diCurrent = grpData->GetDipoleCurrent((AliGRPObject::Stats)0);
    if (diCurrent == AliGRPObject::GetInvalidFloat()) {
      AliError("GRP/GRP/Data entry:  missing value for the dipole current !");
      ok = kFALSE;
    }
    
    Char_t diPolarity = grpData->GetDipolePolarity();
    if (diPolarity == AliGRPObject::GetInvalidChar()) {
      AliError("GRP/GRP/Data entry:  missing value for the dipole polarity !");
      ok = kFALSE;
    }
    
    // read special bits for the polarity convention and map type
    Int_t  polConvention = grpData->IsPolarityConventionLHC() ? AliMagF::kConvLHC : AliMagF::kConvDCS2008;
    Bool_t uniformB = grpData->IsUniformBMap();
    
    if (ok) {
      
      TString beamType = grpData->GetBeamType();
      if (beamType==AliGRPObject::GetInvalidString()) {
        AliError("GRP/GRP/Data entry:  missing value for the beam type ! Using UNKNOWN");
        beamType = "UNKNOWN";
      }
      
      Float_t beamEnergy = grpData->GetBeamEnergy();
      if (beamEnergy==AliGRPObject::GetInvalidFloat()) {
        AliError("GRP/GRP/Data entry:  missing value for the beam energy ! Using 0");
        beamEnergy = 0;
      }

      AliMagF* fld = AliMagF::CreateFieldMap(TMath::Abs(l3Current) * (l3Polarity ? -1:1),
                                             TMath::Abs(diCurrent) * (diPolarity ? -1:1),
                                             polConvention,uniformB,beamEnergy, beamType.Data());
      if (fld) {
        TGeoGlobalMagField::Instance()->SetField( fld );
        TGeoGlobalMagField::Instance()->Lock();
        AliInfo("Running with the B field constructed out of GRP !");
      }
      else AliFatal("Failed to create a B field map !");
    }
    else AliFatal("B field is neither set nor constructed from GRP ! Exitig...");
  }
  
  return 0;
}

//_________________________________________________________________________________________________
int O2Muon::makeClusteringAndTracking(const char* digitInputFile, const char* trackOutputFile,
                                      const char* clusterFinderType, const char* outputLogFile,
                                      int runNumber) {
  /** Make clusters and tracking out of digits.
   */
  
  prepareOCDB(runNumber);
  
  setupMagneticField();
  
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
  
  AliMUONRecoParam* recoParam = getRecoParam(runNumber);
  
  AliMUONTrackReconstructorK trackReco(recoParam,&clusterServer,&transformer);

  AliMUONVClusterStore* clusterStore = new AliMUONClusterStoreV2();
  
  AliMUONVTrackStore* trackStore = new AliMUONTrackStoreV1;
  
  TFile* trackFile = TFile::Open(trackOutputFile,"RECREATE");
  TTree* treeT = new TTree("TreeT","Track Container");
  trackStore->Connect(*treeT,kTRUE);
  
  AliMpArea area; // invalid area to clusterize everything
  
  AliCodeTimerAuto("total",0);
  Int_t nofClusters(0);
  Int_t nofTracks(0);
  Long64_t nofEntries = treeD->GetEntries();
  
  //  std::cout << "Ready to start. Hit a key" << std::endl;
  //  Int_t key;
  //  std::cin >> key;
  
  for ( Long64_t i = 0; i < nofEntries; ++i )
  {
    treeD->GetEntry(i);
    
    TIter next(digitStore->CreateIterator());
    clusterServer.UseDigits(next,digitStore);
    
    {
      AliCodeTimerAuto("Clustering",1);
      for ( int chamberId = 0; chamberId < 10; ++chamberId )
      {
        clusterServer.Clusterize(chamberId,*clusterStore,area,recoParam);
      }
    }
    
    nofClusters += clusterStore->GetSize();
    
    {
      AliCodeTimerAuto("Tracking",2);
      trackReco.EventReconstruct(*clusterStore,*trackStore);
    }
    
    treeT->Fill();
    
    nofTracks += trackStore->GetSize();
    
    clusterStore->Clear();
    digitStore->Clear();
    trackStore->Clear();
  }
  
  trackFile->cd();
  treeT->Write();
  delete trackFile;
  
  std::ofstream out(outputLogFile);
  std::streambuf* coutbuf = std::cout.rdbuf(); //save old buf
  std::cout.rdbuf(out.rdbuf()); //redirect std::cout to out file
  
  AliCodeTimer::Instance()->Print();
  
  std::cout << "Mean number of clusters per event : " << (nofEntries > 0 ? nofClusters*1.0/nofEntries : 0) << std::endl;
  std::cout << "Mean number of tracks per event : " << (nofEntries > 0 ? nofTracks*1.0/nofEntries : 0) << std::endl;
  
  std::cout.rdbuf(coutbuf); //reset to standard output
  
  delete clusterStore;
  delete digitFile;
  delete digitStore;
  delete trackStore;
  
  return 0;
}

//_________________________________________________________________________________________________
void O2Muon::prepareOCDB(int runNumber, AliRawReader* rawReader) {
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
    
    if ( rawReader ) {
      
      AliCDBEntry* entryCTP = AliCDBManager::Instance()->Get("GRP/CTP/Config");
      
      // Load trigger classes for this run
      
      AliTriggerConfiguration* config = static_cast<AliTriggerConfiguration*>(entryCTP->GetObject());
      const TObjArray& classesArray = config->GetClasses();
      TIter nextTriggerClass(&classesArray);
      AliTriggerClass* trclass;
      
      while ( ( trclass = static_cast<AliTriggerClass*>(nextTriggerClass())) ) {
        if (trclass && trclass->GetMask()>0) {
          Int_t trindex = TMath::Nint(TMath::Log2(trclass->GetMask()));
//              AliInfo(Form("class %s %lld index %d",trclass->GetName(),trclass->GetMask(),trindex));
          rawReader->LoadTriggerClass(trclass->GetName(),trindex);
        }
      }
    }
    

  }
  
}

//_________________________________________________________________________________________________
int O2Muon::showClusters(const char* clusterInputFile, Int_t runNumber, Bool_t interactive) {
  prepareOCDB(runNumber);

  TFile* clusterFile = TFile::Open(clusterInputFile);
  if (!clusterFile) {
    return -2;
  }
  
  TTree* treeR = static_cast<TTree*>(clusterFile->Get("TreeR"));
  if (!treeR) {
    return -3;
  }
  
  AliMUONVClusterStore* clusterStore = AliMUONVClusterStore::Create(*treeR);
  
  clusterStore->Connect(*treeR);
  
  Long64_t nofEntries = treeR->GetEntries();
  
  //  std::cout << "Ready to start. Hit a key" << std::endl;
  //  Int_t key;
  //  std::cin >> key;
  
  AliMUONVCluster* cluster = 0x0;
  AliMUONContour* contour = 0x0;
  TObjArray contours;
  contours.SetOwner(kTRUE);
  
  AliCodeTimerAuto("",0);
  
  for ( Long64_t i = 0; i < nofEntries; ++i ) {
    
    treeR->GetEntry(i);
    TIter next(clusterStore->CreateIterator());
    while ( ( cluster = static_cast<AliMUONVCluster*>(next()) ) )
    {
      contour = createClusterContour(*cluster);
      if (contour) {
        contours.Add(contour);
        
        if ( interactive && contour->NumberOfVertices() > 10 ) {
          AliMUONContourMakerTest t;
          TObjArray tmp;
          tmp.Add(contour);
          Double_t xmin, ymin, xmax, ymax;
          
          t.GetBoundingBox(tmp,xmin,ymin,xmax,ymax,kFALSE);
          TCanvas* c = new TCanvas;
          c->Range(xmin,ymin,xmax,ymax);
          t.Plot(*contour,2);
          c->Update();
          
          int n;
          std::cin >> n;
          if (n==0) return 0;
          delete c;
        }
      }
    }
    
    clusterStore->Clear();
  }
  
  std::cout << contours.GetEntries() << std::endl;
  
  AliCodeTimer::Instance()->Print();
  
  return 0;
}

//_________________________________________________________________________________________________
AliMUONContour* O2Muon::createClusterContour(const AliMUONVCluster& cluster) {
  /// Create the contour of a given cluster
  
  
  AliCodeTimerAuto("",0);
  
  TObjArray polygons(cluster.GetNDigits()); // array of AliMUONPolygon objects
  polygons.SetOwner(kTRUE);
  
  UInt_t firstDigitId = cluster.GetDigitId(0);

  Int_t detElemId = AliMUONVDigit::DetElemId(firstDigitId);
  
  const AliMpVSegmentation* seg = AliMpSegmentation::Instance()->GetMpSegmentationByElectronics(AliMUONVDigit::DetElemId(firstDigitId),AliMUONVDigit::ManuId(firstDigitId));
  
  Int_t ndigits = cluster.GetNDigits();
  
  if (ndigits > 200 ) return 0x0;
  
  for ( Int_t i = 0; i < ndigits; ++i ) {
    
    Double_t x,y,dx,dy;
    
    UInt_t digitId = cluster.GetDigitId(i);

    if ( AliMUONVDigit::Cathode(digitId) != AliMUONVDigit::Cathode(firstDigitId) ) continue;
    
    AliMpPad pad = seg->PadByLocation(AliMUONVDigit::ManuId(digitId),AliMUONVDigit::ManuChannel(digitId),kFALSE);
    
    if (!pad.IsValid())
    {
      std::cout << Form("Got an invalid pad for DE %4d manuId %4d channel %4d",detElemId,AliMUONVDigit::ManuId(digitId),AliMUONVDigit::ManuChannel(digitId)) << std::endl;
      cluster.Print();
      continue;
    }
    x = pad.GetPositionX();
    y = pad.GetPositionY();
    dx = pad.GetDimensionX();
    dy = pad.GetDimensionY();
    
    AliMUONPolygon* pol = new AliMUONPolygon(x,y,dx,dy);
    
    polygons.Add(pol);
  }
  
  AliMUONContourMaker maker;
  
  AliMUONContour* contour = maker.CreateContour(polygons);
  
  if (!contour || !contour->IsValid() ) {
    AliError(Form("Failed to properly create contour for cluster. Contour = %p",contour));
    if ( contour ) {
      AliError(Form("nofVertices=%d area.isvalid=%d",contour->NumberOfVertices(),contour->Area().IsValid()));
      StdoutToAliError(contour->Area().Print(););
    }
    delete contour;
    return 0x0;
  }
  
  return contour;
}


