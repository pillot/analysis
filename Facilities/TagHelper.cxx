#include "TagHelper.h"

#include "TChain.h"
#include "TGrid.h"
#include "AliTagAnalysis.h"
#include "TGridResult.h"
#include "TString.h"
#include "AliRunTagCuts.h"
#include "AliEventTagCuts.h"
#include "TSystem.h"
#include "Riostream.h"
#include "AliLHCTagCuts.h"
#include "AliDetectorTagCuts.h"

namespace TagHelper
{

  //______________________________________________________________________________
bool CreateXMLCollectionFromRunList(const char* collectionName, 
                                    const char* runlist,
                                    const char* type,
                                    int passNumber)
{
  /// From a list of run, create a collection with only events containing at least one muon
  /// \param collectionName : output name of the collection (without extension, which will be .xml)
  /// \param runlist : text file containing one integer per line = one run number per line
  /// \param type : files to consider, either ESD or AOD
  /// \param passNumber : 1 or 2 most probably (to distinguish which reco pass is used)
  
  if (!gGrid) TGrid::Connect("alien://");
  if (!gGrid) 
  {
    return 0x0;
  }
  
  TString stype(type);
  stype.ToUpper();

  if ( stype != "ESD" && stype != "AOD" )
  {
    cout << "Only ESD or AOD type supported" << endl;
    return false;
  }
  
  ifstream in(gSystem->ExpandPathName(runlist));
  Int_t runNumber;
  Int_t ntagfiles(0);
  
  AliTagAnalysis tagAnalysis("ESD");
  
  while ( in >> runNumber ) 
  {
    TGridResult *res = gGrid->Query("/alice/data/2009",
				    Form("%09d/%ss/pass%d/*%d*/Run%d.Event*.ESD.tag.root",
					 runNumber,stype.Data(),passNumber,runNumber,runNumber));
//    TGridResult *res = gGrid->Query(" /alice/sim/LHC10a3",
//                                    Form("%d/*/Run%d.Event*.ESD.tag.root", runNumber,runNumber));
    Int_t nFiles = res->GetEntries();
    if (!nFiles) 
    {
      continue;
    }
    
    for (Int_t i = 0; i < nFiles; ++i) 
    {
      TString filename = res->GetKey(i, "turl");
      if(filename == "") continue;
      tagAnalysis.AddTagsFile(filename.Data(),kFALSE);
      ++ntagfiles;
    }
    delete res;
  }
  
  cout << ntagfiles << " tag files added" << endl;
  
  AliRunTagCuts runCuts;
  AliEventTagCuts eventCuts;
  AliLHCTagCuts lhcCuts;
  AliDetectorTagCuts detCuts;
  
  //eventCuts.SetNMuonRange(1,99999);
  
  return tagAnalysis.CreateXMLCollection(collectionName,&runCuts,&lhcCuts,&detCuts,&eventCuts);  
  
}

//______________________________________________________________________________
TChain* CreateChainFromXMLCollection(const char* collectionName, const char* type)
{
  if (!gGrid) TGrid::Connect("alien://");
  if (!gGrid) 
  {
    return 0x0;
  }
  
  TString stype(type);
  stype.ToUpper();
  
  if ( stype != "ESD" && stype != "AOD" )
  {
    cout << "Only ESD or AOD type supported" << endl;
    return 0x0;
  }
  
  AliTagAnalysis tagAnalysis(stype.Data());
  
  return tagAnalysis.CreateChainFromCollection(collectionName,stype=="ESD" ? "esdTree":"aodTree");
}

}
