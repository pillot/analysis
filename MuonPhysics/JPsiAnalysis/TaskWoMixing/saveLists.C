//
//  saveLists.C
//  aliroot_dev
//
//  Created by philippe pillot on 04/07/2014.
//  Copyright (c) 2014 Philippe Pillot. All rights reserved.
//

//______________________________________________________________________________
void saveLists(TString fileName = "Output.root", TString containerName = "cOut_MULorMLL")
{
  /// save the lists of problematic events
  
  gROOT->LoadMacro("$WORK/Macros/Facilities/runTaskFacilities.C");
  LoadAlirootLocally("CORRFW", "", "AliAnalysisTaskJPsi");
  
  TFile* file = TFile::Open(fileName.Data(),"READ");
  TList* container = static_cast<TList*>(file->FindObjectAny(containerName.Data()));
  
  TString listName[12] = {"trgClassMissTrgL0Ev", "trgL0MissTrgClassEv", "trgClassMissTrgOffEv", "trgOffMissTrgClassEv", "trgOffMissTrgL0Ev", "trgL0MissTrgOffEv", "OSTrkOSTrgMLLOnlyEv", "OSTrkLSTrgMLLOnlyEv", "OSTrkLSTrgMULEv", "OSTrkOSTrgFakeMLLOnlyEv", "OSTrkLSTrgFakeMLLOnlyEv", "OSTrkLSTrgFakeMULEv"};
  
  Char_t overwrite = '\0';
  
  for (Int_t i = 0; i < 12; ++i) {
    
    TString outName = Form("%s.txt", listName[i].Data());
    if (!gSystem->AccessPathName(outName.Data())) {
      if (overwrite != 'a' && overwrite != 'k') {
        overwrite = '\0';
        while (overwrite != 'y' && overwrite != 'n' && overwrite != 'a' && overwrite != 'k') {
          cout<<Form("file %s exist in current directory. Overwrite? [y=yes, n=no, a=all, k=keep all] ",outName.Data())<<flush;
          cin>>overwrite;
        }
      }
      if (overwrite == 'y' || overwrite == 'a') {
        gSystem->Exec(Form("rm -f %s", outName.Data()));
        gSystem->Exec(Form("touch %s", outName.Data()));
      } else continue;
    } else gSystem->Exec(Form("touch %s", outName.Data()));
    
    MyList *myList = static_cast<MyList*>(cOut_MULorMLL->FindObject(listName[i].Data()));
    TIter nextBadFile(myList->GetList());
    TObjString *badFile = 0x0;
    while ((badFile = static_cast<TObjString*>(nextBadFile())))
      gSystem->Exec(Form("echo \"%s\" >> %s", badFile->GetName(), outName.Data()));
    
    gSystem->Exec(Form("sort %s > __tmp__", outName.Data()));
    gSystem->Exec(Form("mv -f __tmp__ %s", outName.Data()));
    
  }
  
  file->Close();
  
}

