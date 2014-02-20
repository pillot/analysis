#ifndef AF_H
#define AF_H

class TFileCollection;
class TProof;
#include <vector>
#include "TString.h"
class TList;

void MakeRawCollection(const char* runlist="runList.txt", 
                       const char* collectionName="toto", 
                       const char* collectionScript="makeesdcollection.sh");

class AF : public TObject
{
public:
  AF(const char* af="saf", const char* username="laphecet");
  
  void DryRun(Bool_t flag) { fDryRun = flag; }
  
  Bool_t DryRun() const { return fDryRun; }

  void MergedOnly(Bool_t flag) { fMergedOnly = flag; }
  
  Bool_t MergedOnly() const { return fMergedOnly; }
  
  void PrivateProduction(const char* name, Bool_t simpleRunNumbers=false) { fPrivateProduction = name; fSimpleRunNumbers=simpleRunNumbers; }

  Bool_t IsSimpleRunNumbers() const { return fSimpleRunNumbers; }
  
  Bool_t IsPrivateProduction() const { return fPrivateProduction.Length()>0; }

  TString PrivateProduction() const { return fPrivateProduction; }

  TFileCollection* CreateCollectionFromRunList(const char* collectionType,
                                               const std::vector<int> runs,
                                               const char* dataType,
                                               const char* esdpass="pass1_plusplusplus",
                                               int aodPassNumber=1,
                                               const char* basename="/alice/data/2010/LHC10e",
                                               bool computeTotalSize=true,
                                               Int_t fileLimit=-1);

  TFileCollection* CreateCollectionFromRunList(const char* collectionType,
                                               const char* runlist,
                                               const char* dataType,
                                               const char* esdpass="pass1_plusplusplus",
                                               int aodPassNumber=1,
                                               const char* basename="/alice/data/2010/LHC10e",
                                               bool computeTotalSize=true,
                                               Int_t fileLimit=-1);
  
  TFileCollection* CreateCollectionFromRunList(const char* collectionType,
                                               Int_t runNumber,
                                               const char* dataType,
                                               const char* esdpass="pass1_plusplusplus",
                                               int aodPassNumber=1,
                                               const char* basename="/alice/data/2010/LHC10e",
                                               bool computeTotalSize=true,
                                               Int_t fileLimit=-1);
  
  void CreateDataSets(const std::vector<int>& runs,
                      const char* dataType = "aodmuon",
                      const char* esdpass="pass2",
                      int aodPassNumber=49,
                      const char* basename="/alice/data/2010/LHC10h");
  
  void CreateDataSets(const char* runlist = "aod049.list",
                      const char* dataType = "aodmuon",
                      const char* esdpass="pass2",
                      int aodPassNumber=49,
                      const char* basename="/alice/data/2010/LHC10h");
  
  void CreateDataSets(Int_t runNumber,
                      const char* dataType = "aodmuon",
                      const char* esdpass="pass2",
                      int aodPassNumber=49,
                      const char* basename="/alice/data/2010/LHC10h");

  void DeleteDataSets(const char* dsList, Bool_t forceDataRemoval=kFALSE);
  void DeleteDataSets(const TList& list, Bool_t forceDataRemoval=kFALSE);

  void GetOneDataSetSize(const char* dsname, Int_t& nFiles, Int_t& nCorruptedFiles, Long64_t& size, Bool_t showDetails=kFALSE);

  void GroupDatasets();

  void GroupDatasets(const TList& list);
  void GroupDatasets(const char* dslist);

  void GetDataSetsSize(const char* dslist, Bool_t showDetails=kFALSE);
  
  void MergeOneDataSet(const char* dsname);
  
  void MergeDataSets(const char* dsList);
  
  void CompareRunLists(const char* runlist1, const char* runlist2);
  
  void StopDataSets(const char* dsMask);

  void ShowDataSetContent(const char* dsset);
  void ShowDataSetContent(const TList& dsset);
  
  void ShowDataSets(const char* runlist);
  
private:
  void ReadIntegers(const char* filename, std::vector<int>& integers);
  void DeleteOneDataSet(TProof* proof, const char* dsName, Bool_t forceDataRemoval=kFALSE);
  void StopOneDataSet(TProof* proof, const char* dsname);
  
private:
  TString fConnect; // Connect string (username@afmaster)
  Bool_t fDryRun; // whether to do real things or just show what would be done
  Bool_t fMergedOnly; // pick only the merged AODs when merging is in the same directory as non-merged...
  TString fPrivateProduction; // dataset(s) basename for private productions
  Bool_t fSimpleRunNumbers; // true if run number are not with format %09d but %d instead
  
//  const char* connect = "laphecet@nansafmaster.in2p3.fr/?N",
  
  ClassDef(AF,3)
};

#endif
