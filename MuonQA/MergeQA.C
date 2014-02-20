/*
 *  MergeQA.C
 *  aliroot
 *
 *  Created by Philippe Pillot on 19/07/10.
 *  Copyright 2010 SUBATECH. All rights reserved.
 *
 */

#if !defined(__CINT__) || defined(__MAKECINT__)

#include <Riostream.h>

// ROOT includes
#include "TString.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TGrid.h"
#include "TGridResult.h"
#include "TMap.h"
#include "TMath.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TH1.h"
#include "TChain.h"
#include "TKey.h"
#include "THashList.h"
#include "TSystem.h"
#endif

Bool_t MergeRecursive(TDirectory*, TList*);

void MergeQA(const char* alienBaseDir, const char* fileName = "AnalysisResults.root", const char* runList = 0x0)
{
  /// Merge all files "fileName" found in "alienBaseDir" skipping the kParamContainer.
  /// Example: alienBaseDir = "/alice/cern.ch/user/p/ppillot/pp7TeV/LHC10d/MuonQA/pass1/results/".
  /// If runList != 0x0: only merge results of the selected runs.
  
  // Load common libraries
  gSystem->Load("libTree");
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");
  gSystem->Load("libPhysics");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCORRFW");
  gSystem->Load("libPWG3base");
  gSystem->Load("libPWG3muon");
  
  if (!TGrid::Connect("alien://")) {
    Error("MergeQA","cannot connect to grid");
    return;
  }
  
  // list runs to be merged together
  TString currRun;
  TObjArray runs;
  runs.SetOwner();
  if (runList) {
    
    // only the ones in the runList
    ifstream inFile(runList);
    if (!inFile.is_open()) {
      Error("MergeQA",Form("unable to open file %s", runList));
      return;
    }
    
    while (!inFile.eof()) {
      currRun.ReadLine(inFile, kTRUE);
      if (currRun.IsNull()) continue;
      if (!currRun.IsDigit()) {
	Error("MergeQA","invalid run number: %s", currRun.Data());
	return;
      }
      runs.AddLast(new TObjString(Form("%09d", currRun.Atoi())));
    }
    
    inFile.close();
    
  } else {
    
    // all runs
    runs.AddLast(new TObjString("*"));
    
  }
  
  // create output file
  TFile* outputFile = TFile::Open(fileName, "RECREATE");
  if (!outputFile) {
    Error("MergeQA", "cannot open the MERGER output file %s", fileName);
    return;
  }
  
  // loop over runs
  TList fileList;
  for ( Int_t irun=0; irun<runs.GetEntriesFast(); irun++ ) {
    
    TString run = ((TObjString*)runs.UncheckedAt(irun))->GetString();
    
    // get the list of files to merge
    TString command = Form("find %s/ %s/*%s", alienBaseDir, run.Data(), fileName);
    TGridResult *res = gGrid->Command(command);
    if (!res) {
      Error("MergeQA",Form("no result for the command: %s",command.Data()));
      return;
    }
    
    // Loop over 'find' results and get next LFN
    TIter nextmap(res);
    TMap *map = 0;
    while ((map=(TMap*)nextmap())) {
      
      TObjString *objs = dynamic_cast<TObjString*>(map->GetValue("turl"));
      if (!objs || !objs->GetString().Length()) {
	Error("MergeQA","turl not found for the run %s... SKIPPING", run.Data());
	continue;
      }
      
      // Add file to be merged
      fileList.Add(TFile::Open(objs->GetString()));
      
    }
    
    delete res;
    
  }
  
  // merge files
  Bool_t result = MergeRecursive(outputFile, &fileList);
  if (!result) Error("MergeQA", "error during merge of your ROOT files");
  else outputFile->Close();
  
}

//______________________________________________________________________________
Bool_t MergeRecursive(TDirectory *target, TList *sourcelist)
{
  // Merge all objects in a directory
  // NB. This function is a copy of the hadd function MergeROOTFile
  
  Bool_t fHistoOneGo = kFALSE;
  Bool_t fNoTrees = kFALSE;
  Bool_t fFastMethod = kTRUE;
  
  // Get the dir name
  TString path(target->GetPath());
  path.Remove(0, path.Last(':') + 2);
  
  //gain time, do not add the objects in the list in memory
  Bool_t addDirStat = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  
  TDirectory *first_source = (TDirectory*)sourcelist->First();
  
  Int_t nguess = sourcelist->GetSize()+1000;
  THashList allNames(nguess);
  ((THashList*)target->GetList())->Rehash(nguess);
  ((THashList*)target->GetListOfKeys())->Rehash(nguess);
  
  while (first_source) {
    TDirectory *current_sourcedir = first_source->GetDirectory(path);
    if (!current_sourcedir) {
      first_source = (TDirectory*)sourcelist->After(first_source);
      continue;
    }
    
    // loop over all keys in this directory
    TChain *globChain = 0;
    TIter nextkey( current_sourcedir->GetListOfKeys() );
    TKey *key, *oldkey=0;
    
    while ( (key = (TKey*)nextkey())) {
      if (current_sourcedir == target) break;
      
      // skip the kParamContainer
      if (!strcmp(key->GetName(),"general2")) continue;
      
      // keep only the highest cycle number for each key
      if (oldkey && !strcmp(oldkey->GetName(),key->GetName())) continue;
      if (!strcmp(key->GetClassName(),"TProcessID")) {key->ReadObj(); continue;}
      if (allNames.FindObject(key->GetName())) continue;
      TClass *cl = TClass::GetClass(key->GetClassName());
      if (!cl || !cl->InheritsFrom(TObject::Class())) {
	Info("MergeRecursive", "cannot merge object type, name: %s title: %s",
	     key->GetName(), key->GetTitle());
	continue;
      }
      allNames.Add(new TObjString(key->GetName()));
      
      // read object from first source file
      current_sourcedir->cd();
      TObject *obj = key->ReadObj();
      if (!obj) {
	Info("MergeRecursive", "could not read object for key {%s, %s}",
	     key->GetName(), key->GetTitle());
	continue;
      }
      
      if (obj->IsA()->InheritsFrom("TH1")) {
	// descendant of TH1 -> merge it
	
	TH1 *h1 = (TH1*)obj;
	TList listH;
	
	// loop over all source files and add the content of the
	// correspondant histogram to the one pointed to by "h1"
	TFile *nextsource = (TFile*)sourcelist->After( first_source );
	while ( nextsource ) {
	  // make sure we are at the correct directory level by cd'ing to path
	  TDirectory *ndir = nextsource->GetDirectory(path);
	  if (ndir) {
	    ndir->cd();
	    TKey *key2 = (TKey*)gDirectory->GetListOfKeys()->FindObject(key->GetName());
	    if (key2) {
	      TObject *hobj = key2->ReadObj();
	      hobj->ResetBit(kMustCleanup);
	      listH.Add(hobj);
	      // Run the merging now, if required
	      if (!fHistoOneGo) {
		h1->Merge(&listH);
		listH.Delete();
	      }
	    }
	  }
	  nextsource = (TFile*)sourcelist->After( nextsource );
	}
	// Merge the list, if still to be done
	if (fHistoOneGo) {
	  h1->Merge(&listH);
	  listH.Delete();
	}
      } else if ( obj->IsA()->InheritsFrom( "TTree" ) ) {
	
	// loop over all source files create a chain of Trees "globChain"
	if (!fNoTrees) {
	  TString obj_name;
	  if (path.Length()) {
	    obj_name = path + "/" + obj->GetName();
	  } else {
	    obj_name = obj->GetName();
	  }
	  globChain = new TChain(obj_name);
	  globChain->Add(first_source->GetName());
	  TFile *nextsource = (TFile*)sourcelist->After( first_source );
	  while ( nextsource ) {
	    //do not add to the list a file that does not contain this Tree
	    TFile *curf = TFile::Open(nextsource->GetName());
	    if (curf) {
	      Bool_t mustAdd = kFALSE;
	      if (curf->FindKey(obj_name)) {
		mustAdd = kTRUE;
	      } else {
		//we could be more clever here. No need to import the object
		//we are missing a function in TDirectory
		TObject *aobj = curf->Get(obj_name);
		if (aobj) { mustAdd = kTRUE; delete aobj;}
	      }
	      if (mustAdd) {
		globChain->Add(nextsource->GetName());
	      }
	    }
	    delete curf;
	    nextsource = (TFile*)sourcelist->After( nextsource );
	  }
	}
      } else if ( obj->IsA()->InheritsFrom( "TDirectory" ) ) {
	// it's a subdirectory
	
	//cout << "Found subdirectory " << obj->GetName() << endl;
	// create a new subdir of same name and title in the target file
	target->cd();
	TDirectory *newdir = target->mkdir( obj->GetName(), obj->GetTitle() );
	
	// newdir is now the starting point of another round of merging
	// newdir still knows its depth within the target file via
	// GetPath(), so we can still figure out where we are in the recursion
	MergeRecursive( newdir, sourcelist);
	
      } else if (obj->InheritsFrom(TObject::Class()) &&
		 obj->IsA()->GetMethodWithPrototype("Merge", "TCollection*") ) {
	// Object implements Merge(TCollection*)
	
	TList listH;
	TString listHargs;
	listHargs.Form("((TCollection*)0x%lx)",&listH);
	
	// Loop over all source files and merge same-name object
	TFile *nextsource = (TFile*)sourcelist->After( first_source );
	while (nextsource) {
	  // make sure we are at the correct directory level by cd'ing to path
	  TDirectory *ndir = nextsource->GetDirectory(path);
	  if (ndir) {
	    ndir->cd();
	    TKey *key2 = (TKey*)gDirectory->GetListOfKeys()->FindObject(key->GetName());
	    if (key2) {
	      TObject *hobj = key2->ReadObj();
	      // Set ownership for collections
	      if (hobj->InheritsFrom(TCollection::Class())) {
		((TCollection*)hobj)->SetOwner(); 
	      }   
	      hobj->ResetBit(kMustCleanup);
	      listH.Add(hobj);
	      Int_t error = 0;
	      obj->Execute("Merge", listHargs.Data(), &error);
	      if (error) {
		Error("MergeRecursive", "calling Merge() on '%s' with the corresponding object in '%s'",
		      obj->GetName(), nextsource->GetName());
	      }
	      listH.Delete();
	    }
	  }
	  nextsource = (TFile*)sourcelist->After( nextsource );
	}
      } else {
	// Object is of no type that we can merge 
	Warning("MergeRecursive", "cannot merge object type (n:'%s', t:'%s') - "
		"Merge(TCollection *) not implemented",
		obj->GetName(), obj->GetTitle());
	
	// Loop over all source files and write similar objects directly to the output file
	TFile *nextsource = (TFile*)sourcelist->After( first_source );
	while (nextsource) {
	  // make sure we are at the correct directory level by cd'ing to path
	  TDirectory *ndir = nextsource->GetDirectory(path);
	  if (ndir) {
	    ndir->cd();
	    TKey *key2 = (TKey*)gDirectory->GetListOfKeys()->FindObject(key->GetName());
	    if (key2) {
	      TObject *nobj = key2->ReadObj();
	      nobj->ResetBit(kMustCleanup);
	      if (target->WriteTObject(nobj, key2->GetName(), "SingleKey") <= 0) {
		Warning("MergeRecursive", "problems copying object (n:'%s', t:'%s') to output file ",
			obj->GetName(), obj->GetTitle());
	      }
	      delete nobj;
	    }
	  }
	  nextsource = (TFile*)sourcelist->After( nextsource );
	}
      }
      
      // now write the merged histogram (which is "in" obj) to the target file
      // note that this will just store obj in the current directory level,
      // which is not persistent until the complete directory itself is stored
      // by "target->Write()" below
      target->cd();
      
      //!!if the object is a tree, it is stored in globChain...
      if(obj->IsA()->InheritsFrom( "TDirectory" )) {
	//printf("cas d'une directory\n");
      } else if(obj->IsA()->InheritsFrom( "TTree" )) {
	if (!fNoTrees) {
	  if (globChain) {
	    globChain->ls();
	    if (fFastMethod) globChain->Merge(target->GetFile(),0,"keep fast");
	    else             globChain->Merge(target->GetFile(),0,"keep");
	    delete globChain;
	  }
	}
      } else if (obj->IsA()->InheritsFrom( "TCollection" )) {
	obj->Write( key->GetName(), TObject::kSingleKey );
	((TCollection*)obj)->SetOwner();
      } else {
	obj->Write( key->GetName() );
      }
      if (obj->IsA()->InheritsFrom("TCollection")) ((TCollection*)obj)->Delete();
      oldkey = key;
      delete obj;
    } // while ( ( TKey *key = (TKey*)nextkey() ) )
    first_source = (TDirectory*)sourcelist->After(first_source);
  }
  // save modifications to target file
  target->SaveSelf(kTRUE);
  TH1::AddDirectory(addDirStat);
  return kTRUE;
}

