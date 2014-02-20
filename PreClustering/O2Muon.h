#ifndef O2MUON_H
#define O2MUON_H

#include "TObject.h"
#include <string>

class AliMUONRecoParam;

class O2Muon : public TObject
{
public:
  O2Muon(const char* ocdbpath="local:///Users/laurent/Alice/OCDBcopy2013");
  virtual ~O2Muon();
  
  /** Create a Root file with (calibrated) MCH digits from a raw data file.
    * Not meant to be fast, just a re-use of the existing classes to get the job done. 
    */
  int makeDigitFile(const char* rawDataInputFile="/alice/data/2013/LHC13f/000196474/raw/13000196474082.15filtered.root", const char* digitOutputFile="digits.root", Bool_t calibrate=kTRUE);

  /// Make pre-clusters out of digits
  int makeClustering(const char* digitInputFile="digits.root", const char* clusterFinderType="PRECLUSTER",
                     const char* outputLogFile="precluster.log",int runNumber=196474);
  
private:
  
  AliMUONRecoParam* getRecoParam(int runNumber);

  void prepareOCDB(int runNumber);
  
private:
  std::string mOCDBPath;
  
  ClassDef(O2Muon,0)
};

#endif
