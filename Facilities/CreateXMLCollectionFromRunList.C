/*
 *  CreateXMLCollectionFromRunList.C
 *  aliroot
 *
 *  Created by Philippe Pillot on 12/01/10.
 *  Copyright 2010 SUBATECH. All rights reserved.
 *
 */

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "/Users/philippe/Work/Alice/Work/Data/Macro/TagHelper.h"
#include "TROOT.h"
#include "TSystem.h"
#endif

void CreateXMLCollectionFromRunList(const char* collectionName, 
				    const char* runlist,
				    const char* type,
				    int passNumber)
{
  /// create xml collection from provided run list using tag, via the TagHelper:
  /// \param collectionName : output name of the collection (without extension, which will be .xml)
  /// \param runlist : text file containing one integer per line = one run number per line
  /// \param type : files to consider, either ESD or AOD
  /// \param passNumber : 1 or 2 most probably (to distinguish which reco pass is used)
  
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gROOT->LoadMacro("/Users/philippe/Work/Alice/Work/Data/Macro/TagHelper.cxx++");
  
  gROOT->ProcessLine(Form("TagHelper::CreateXMLCollectionFromRunList(\"%s\",\"%s\",\"%s\",%d);",
			  collectionName, runlist, type, passNumber));
  
}

