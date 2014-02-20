#ifndef TAGHELPER_H
#define TAGHELPER_H

class TChain;

namespace TagHelper
{
  bool CreateXMLCollectionFromRunList(const char* collectionName, 
                                      const char* runlist,
                                      const char* type="ESD",
                                      int passNumber=2);
  
  TChain* CreateChainFromXMLCollection(const char* collectionName, const char* type="ESD");
};

#endif
