#ifndef O2MUON_H
#define O2MUON_H

#include "TObject.h"
#include <string>

class AliMUONRecoParam;
class AliMUONContour;
class AliMUONVCluster;
class AliRawReader;

class O2Muon : public TObject
{
public:
  O2Muon(const char* ocdbpath="local:///Users/laurent/Alice/OCDBcopy2011");
  virtual ~O2Muon();
  
  /** simple loop over events with decoding of the MCH data */
  int decodeEvents(const char* rawDataInputFile);

  /** Create a new raw data file with only the events matching the given trigger.
   */
  int filterRaw(const char* rawDataInputFile="/alice/data/2011/LHC11h/000169099/raw/11000169099001.28.FILTER_RAWMUON.root", const char* triggerClass="CPBI2_B1-B-NOPF-ALLNOTRD");
  
  /** Create a Root file with (calibrated) MCH digits from a raw data file.
    * Not meant to be fast, just a re-use of the existing classes to get the job done. 
    */
  int makeDigitFile(const char* rawDataInputFile="/alice/data/2011/LHC11h/000169099/raw/11000169099001.28.FILTER_RAWMUON.root", const char* digitOutputFile="digits.root", Bool_t calibrate=kTRUE);

  int makeDigitFiles(const char* rawDataInputFileList, const char* triggerClass="", Bool_t calibrate=kTRUE);
  
  /// Make pre-clusters out of digits
  int makeClustering(const char* digitInputFile="digits.root",
                     const char* clusterOutputFile="clusters.root",
                     const char* clusterFinderType="PRECLUSTER",
                     const char* outputLogFile="precluster.log",int runNumber=196474);

  /// Make a contour of the cluster
  AliMUONContour* createClusterContour(const AliMUONVCluster& cluster);

  /// Make a display of the clusters
  int showClusters(const char* clusterInputFile="clusters.root", Int_t runNumber=196474, Bool_t interactive=kTRUE);

  /// Go from digits to tracks
  int makeClusteringAndTracking(const char* digitInputFile, const char* trackOutputFile,
                                const char* clusterFinderType, const char* outputLogFile,
                                int runNumber);

private:

  int getRunNumber(AliRawReader& rawReader) const;

  AliMUONRecoParam* getRecoParam(int runNumber);

  void prepareOCDB(int runNumber, AliRawReader* rawReader=0x0);

  int setupMagneticField();
  
private:
  std::string mOCDBPath;
  
  ClassDef(O2Muon,0)
};

#endif
