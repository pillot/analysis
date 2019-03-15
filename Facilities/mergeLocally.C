/*
 *  merge.C
 *  aliroot_dev
 *
 *  Created by philippe pillot on 04/11/11.
 *  Copyright 2011 SUBATECH. All rights reserved.
 *
 */

#if !defined(__CINT__) || defined(__MAKECINT__)

#include <Riostream.h>

// ROOT includes
#include "TROOT.h"
#include "TSystem.h"
#include "TClassRef.h"
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
#include "TFileMerger.h"
#include "TFileMergeInfo.h"
#include "TClonesArray.h"
#include "TList.h"

#endif

TList paramContainers;
Int_t nFilesMax = 10; // maximum number of files merged in one step

TList *CreateFileList(TList &fileNames);
Bool_t MergeRecursive(TDirectory *target, TList *sourcelist, Int_t type);
void Merge(TList &files, TString &fileName, Bool_t skipParamContainers);

void mergeLocally(TString fileList = "fileList.txt", Bool_t skipParamContainers = kFALSE, Bool_t runList = kFALSE)
{
  /// merge all files in the list
  /// *** to be compiled ***
  
  paramContainers.SetOwner();
  // performance task
  paramContainers.AddLast(new TObjString("Efficiency"));
  paramContainers.AddLast(new TObjString("MomentumAtVtx"));
  paramContainers.AddLast(new TObjString("SlopeAtVtx"));
  paramContainers.AddLast(new TObjString("EtaAtVtx"));
  paramContainers.AddLast(new TObjString("PhiAtVtx"));
  paramContainers.AddLast(new TObjString("MomentumAtFirstCl"));
  paramContainers.AddLast(new TObjString("SlopeAtFirstCl"));
  paramContainers.AddLast(new TObjString("DCA"));
  paramContainers.AddLast(new TObjString("Clusters"));
  // efficiency task
  paramContainers.AddLast(new TObjString("ClustersCounters"));
  // QA task
  paramContainers.AddLast(new TObjString("general2"));
  // resolution task
  paramContainers.AddLast(new TObjString("LocalChi2"));
  paramContainers.AddLast(new TObjString("ChamberRes"));

  // load potentially needed libraries
    gROOT->LoadMacro("$HOME/Work/Alice/Macros/Facilities/runTaskFacilities.C");
		TString extraLibs="";
		TString extraPkgs="";
//  TString extraPkgs="PWGPPMUONdep:PWGPPMUONlite";
//  TString extraPkgs="PWGPPMUONdep:PWGPPMUONlite:PWGmuondep";
    gROOT->ProcessLineFast(Form("LoadAlirootLocally(\"%s\", \"\", \"\", \"%s\")",extraLibs.Data(),extraPkgs.Data()));
	//  gROOT->ProcessLineFast(Form("LoadAlirootLocally(\"%s\", \"include\", \"AliAnalysisTaskJPsi:AliAnalysisTaskMTRSign\")",extraLibs.Data()));
  
  // open the file list
  ifstream inFile(fileList.Data());
  if (!inFile.is_open()) {
    Error("merge", "unable to open file %s", fileList.Data());
    return;
  }
  
  // loop over files
  TString fileName;
  TClonesArray fileLists("TList", 100);
  TList *lastList = 0x0;
  Int_t nFiles = 0;
  while (!inFile.eof()) {
    
    // get the current file
    TString currFile;
    currFile.ReadLine(inFile, kTRUE);
    if (currFile.IsNull()) continue;
    
    // in case the file list is in fact a run list
    if (runList) {
      currFile.Prepend("./runs/");
//      currFile += "/AnalysisResults.root";
      currFile += "/QAresults.root";
//      currFile += "/chamberResolution_step2.root";
    }
    
    // pick up the fileName if not already done
    if (fileName.IsNull()) {
      TObjArray* tmp = currFile.Tokenize("/");
      fileName = ((TObjString*)tmp->Last())->String();
      delete tmp;
    }
    
    // create a new list if needed
    if (nFiles%nFilesMax == 0) {
      lastList = new (fileLists[fileLists.GetEntriesFast()]) TList();
      lastList->SetOwner(kTRUE);
    }
    
    // add the file for merging
    lastList->Add(new TObjString(currFile.Data()));
    nFiles++;
    
  }
  
  // close the run list
  inFile.close();
  
  // nothing found - skip this output
  if (fileLists.GetEntriesFast() == 0) {
    Error("merge","no file found");
    return;
  }
  
  // loop over lists of files and merge them till we get only one file
  gSystem->Exec("rm -f __merge*__*");
  TList *fileNames = static_cast<TList*>(fileLists.First());
  Int_t nMerges = 0;
  while (fileNames) {
    
    TList *nextFileNames = static_cast<TList*>(fileLists.After(fileNames));
    
    // intermediate or final file name
    TString tmpFileName(fileName);
    if (nextFileNames) tmpFileName.Prepend(Form("__merge%d__", nMerges++));
    
    // create the list of files
    TList *files = CreateFileList(*fileNames);
    
    // merge
    Merge(*files, tmpFileName, skipParamContainers);
    
    // close input files
    TIter nextFile(files);
    TFile *file;
    while ((file = static_cast<TFile*>(nextFile()))) file->Close();
    delete files;
    
    // add intermediate file to the list of files to be merged. Create a new list if needed
    if (nextFileNames) {
      if (nFiles%nFilesMax == 0) lastList = new (fileLists[fileLists.GetEntriesFast()]) TList();
      lastList->Add(new TObjString(Form("%s/%s", gSystem->pwd(), tmpFileName.Data())));
      nFiles++;
    }
    
    fileNames = nextFileNames;
    
  }
  
  gSystem->Exec("rm -f __merge*__*");
  
}

//______________________________________________________________________________
TList *CreateFileList(TList &fileNames)
{
  /// create the list of files
  
  TList *files = new TList;
  
  TIter nextFileName(&fileNames);
  TObjString *fileName;
  while ((fileName = static_cast<TObjString*>(nextFileName())))
    files->Add(TFile::Open(fileName->GetName()));
  
  return files;
  
}

//______________________________________________________________________________
void Merge(TList &files, TString &fileName, Bool_t skipParamContainers)
{
  /// merge all files in the list in a file of the given name
  
  printf("merging %d files...\n", files.GetEntries());
  files.Print();
  
  if (skipParamContainers) {
    
    // create output file
    TFile* outputFile = TFile::Open(fileName.Data(), "RECREATE");
    if (!outputFile) {
      Error("merge", "cannot open the MERGER output file %s", fileName.Data());
      return;
    }
    
    // merge
    if (!MergeRecursive(outputFile, &files, TFileMerger::kRegular|TFileMerger::kAll))
      Error("merge", "could not merge all files");
    
    // close ouput file
    outputFile->Close();
    
  } else {
    
    // create TFileMerger
    TFileMerger *fm = new TFileMerger(kFALSE);
    fm->SetFastMethod(kTRUE);
    
    // add files
    TIter nextFile(&files);
    TFile *file;
    while ((file = static_cast<TFile*>(nextFile()))) fm->AddFile(file);
    
    // open output file
    fm->OutputFile(fileName);
    
    // merge (closing of output file is done at the end of merging)
    if (!fm->Merge()) Error("merge", "could not merge all files");
    
    // clean memory
    delete fm;
    
  }
  
}

//______________________________________________________________________________
Bool_t MergeRecursive(TDirectory *target, TList *sourcelist, Int_t type)
{
  // Merge all objects in a directory
  // The type is defined by the bit values in EPartialMergeType:
  //   kRegular      : normal merge, overwritting the output file (default)
  //   kIncremental  : merge the input file with the (existing) content of the output file (if already exising)
  //   kAll          : merge all type of objects (default)
  //   kResetable    : merge only the objects with a MergeAfterReset member function.
  //   kNonResetable : merge only the objects without a MergeAfterReset member function.
  
  TString fMsgPrefix = "";
  Int_t  fPrintLevel = 0;
  Bool_t fHistoOneGo = kFALSE;
  Bool_t fNoTrees = kFALSE;
  Bool_t fFastMethod = kTRUE;
  Bool_t fCompressionChange = kTRUE;
  TClassRef R__TH1_Class("TH1");
  TClassRef R__TTree_Class("TTree");
  
  Bool_t status = kTRUE;
  if (fPrintLevel > 0) {
    Printf("%s Target path: %s",fMsgPrefix.Data(),target->GetPath());
  }
  
  // Get the dir name
  TString path(target->GetPath());
  // coverity[unchecked_value] 'target' is from a file so GetPath always returns path starting with filename: 
  path.Remove(0, path.Last(':') + 2);
  
  Int_t nguess = sourcelist->GetSize()+1000;
  THashList allNames(nguess);
  allNames.SetOwner(kTRUE);
  ((THashList*)target->GetList())->Rehash(nguess);
  ((THashList*)target->GetListOfKeys())->Rehash(nguess);
  
  TFileMergeInfo info(target);
  
  if ((fFastMethod && !fCompressionChange)) {
    info.fOptions.Append(" fast");
  }
  
  TFile      *current_file;
  TDirectory *current_sourcedir;
  if (type & TFileMerger::kIncremental) {
    current_file      = 0;
    current_sourcedir = target;
  } else {
    current_file      = (TFile*)sourcelist->First();
    current_sourcedir = current_file->GetDirectory(path);
  }
  
  while (current_file || current_sourcedir) {
    // When current_sourcedir != 0 and current_file == 0 we are going over the target
    // for an incremental merge.
    if (current_sourcedir && (current_file == 0 || current_sourcedir != target)) {
      
      // loop over all keys in this directory
      TIter nextkey( current_sourcedir->GetListOfKeys() );
      TKey *key;
      TString oldkeyname;
      
      while ( (key = (TKey*)nextkey())) {
	
        // only merge MUON outputs
        if (!path.Contains("MUON") && !strstr(key->GetName(),"MUON")) continue;

	// skip the kParamContainer
	Bool_t pcFound = kFALSE;
	TIter nextpc(&paramContainers);
	TObjString *pc = 0x0;
	while ((pc=static_cast<TObjString*>(nextpc()))) {
	  if (!strcmp(key->GetName(),pc->GetName())) {
	    pcFound = kTRUE;
	    break;
	  }
	}
	if (pcFound) continue;
	
	// Keep only the highest cycle number for each key.  They are stored in the (hash) list
	// consecutively and in decreasing order of cycles, so we can continue until the name
	// changes.
	if (oldkeyname == key->GetName()) continue;
	// Read in but do not copy directly the processIds.
	if (strcmp(key->GetClassName(),"TProcessID") == 0) { key->ReadObj(); continue;}
	// If we have already seen this object [name], we already processed
	// the whole list of files for this objects and we can just skip it
	// and any related cycles.
	if (allNames.FindObject(key->GetName())) {
	  oldkeyname = key->GetName();
	  continue;
	}
	
	TClass *cl = TClass::GetClass(key->GetClassName());
	if (!cl || !cl->InheritsFrom(TObject::Class())) {
	  Info("MergeRecursive", "cannot merge object type, name: %s title: %s",
	       key->GetName(), key->GetTitle());
	  continue;
	}
	allNames.Add(new TObjString(key->GetName()));
	
	if (fNoTrees && cl->InheritsFrom(R__TTree_Class)) {
	  // Skip the TTree objects and any related cycles.
	  oldkeyname = key->GetName();
	  continue;
	}
	
	if (!(type&TFileMerger::kResetable && type&TFileMerger::kNonResetable)) {
	  // If neither or both are requested at the same time, we merger both types.
	  if (!(type&TFileMerger::kResetable)) {
	    if (cl->GetResetAfterMerge()) {
	      // Skip the object with a reset after merge routine (TTree and other incrementally mergeable objects)
	      oldkeyname = key->GetName();
	      continue;                  
	    }
	  }
	  if (!(type&TFileMerger::kNonResetable)) {
	    if (!cl->GetResetAfterMerge()) {
	      // Skip the object without a reset after merge routine (Histograms and other non incrementally mergeable objects)
	      oldkeyname = key->GetName();
	      continue;                  
	    }
	  }
	}
	// read object from first source file
	TObject *obj;
	if (type & TFileMerger::kIncremental) {
	  obj = current_sourcedir->GetList()->FindObject(key->GetName());
	  if (!obj) {
	    obj = key->ReadObj();
	  }
	} else {
	  obj = key->ReadObj();
	}
	if (!obj) {
	  Info("MergeRecursive", "could not read object for key {%s, %s}",
	       key->GetName(), key->GetTitle());
	  continue;
	}
	
	if ( obj->IsA()->InheritsFrom( TDirectory::Class() ) ) {
	  // it's a subdirectory
	  
	  target->cd();
	  TDirectory *newdir;
	  if (type & TFileMerger::kIncremental) {
	    newdir = target->GetDirectory(obj->GetName());
	    if (!newdir) {
	      newdir = target->mkdir( obj->GetName(), obj->GetTitle() );
	    }
	  } else {
	    newdir = target->mkdir( obj->GetName(), obj->GetTitle() );
	  }
	  
	  // newdir is now the starting point of another round of merging
	  // newdir still knows its depth within the target file via
	  // GetPath(), so we can still figure out where we are in the recursion
	  status = MergeRecursive(newdir, sourcelist, type);
	  if (!status) return status;
	  
	} else if (obj->IsA()->GetMerge()) {
	  
	  TList inputs;
	  Bool_t oneGo = fHistoOneGo && obj->IsA()->InheritsFrom(R__TH1_Class);
	  
	  // Loop over all source files and merge same-name object
	  TFile *nextsource = current_file ? (TFile*)sourcelist->After( current_file ) : (TFile*)sourcelist->First();
	  if (nextsource == 0) {
	    // There is only one file in the list
	    ROOT::MergeFunc_t func = obj->IsA()->GetMerge();
	    func(obj, &inputs, &info);
	    info.fIsFirst = kFALSE;
	  } else {
	    do {
	      // make sure we are at the correct directory level by cd'ing to path
	      TDirectory *ndir = nextsource->GetDirectory(path);
	      if (ndir) {
		ndir->cd();
		TKey *key2 = (TKey*)ndir->GetListOfKeys()->FindObject(key->GetName());
		if (key2) {
		  TObject *hobj = key2->ReadObj();
		  if (!hobj) {
		    Info("MergeRecursive", "could not read object for key {%s, %s}; skipping file %s",
			 key->GetName(), key->GetTitle(), nextsource->GetName());
		    nextsource = (TFile*)sourcelist->After(nextsource);
		    continue;
		  }
		  // Set ownership for collections
		  if (hobj->InheritsFrom(TCollection::Class())) {
		    ((TCollection*)hobj)->SetOwner();
		  }
		  hobj->ResetBit(kMustCleanup);
		  inputs.Add(hobj);
		  if (!oneGo) {
		    ROOT::MergeFunc_t func = obj->IsA()->GetMerge();
		    Long64_t result = func(obj, &inputs, &info);
		    info.fIsFirst = kFALSE;
		    if (result < 0) {
		      Error("MergeRecursive", "calling Merge() on '%s' with the corresponding object in '%s'",
			    obj->GetName(), nextsource->GetName());
		    }
		    inputs.Delete();
		  }
		}
	      }
	      nextsource = (TFile*)sourcelist->After( nextsource );
	    } while (nextsource);
	    // Merge the list, if still to be done
	    if (oneGo || info.fIsFirst) {
	      ROOT::MergeFunc_t func = obj->IsA()->GetMerge();
	      func(obj, &inputs, &info);
	      info.fIsFirst = kFALSE;
	      inputs.Delete();
	    }
	  }
	} else if (obj->InheritsFrom(TObject::Class()) &&
		   obj->IsA()->GetMethodWithPrototype("Merge", "TCollection*,TFileMergeInfo*") ) {
	  // Object implements Merge(TCollection*,TFileMergeInfo*) and has a reflex dictionary ... 
	  
	  TList listH;
	  TString listHargs;
	  listHargs.Form("(TCollection*)0x%lx,(TFileMergeInfo*)0x%lx", (ULong_t)&listH,(ULong_t)&info);
	  
	  // Loop over all source files and merge same-name object
	  TFile *nextsource = current_file ? (TFile*)sourcelist->After( current_file ) : (TFile*)sourcelist->First();
	  if (nextsource == 0) {
	    // There is only one file in the list
	    Int_t error = 0;
	    obj->Execute("Merge", listHargs.Data(), &error);
	    info.fIsFirst = kFALSE;
	    if (error) {
	      Error("MergeRecursive", "calling Merge() on '%s' with the corresponding object in '%s'",
		    obj->GetName(), key->GetName());
	    }
	  } else {
	    while (nextsource) {
	      // make sure we are at the correct directory level by cd'ing to path
	      TDirectory *ndir = nextsource->GetDirectory(path);
	      if (ndir) {
		ndir->cd();
		TKey *key2 = (TKey*)ndir->GetListOfKeys()->FindObject(key->GetName());
		if (key2) {
		  TObject *hobj = key2->ReadObj();
		  if (!hobj) {
		    Info("MergeRecursive", "could not read object for key {%s, %s}; skipping file %s",
			 key->GetName(), key->GetTitle(), nextsource->GetName());
		    nextsource = (TFile*)sourcelist->After(nextsource);
		    continue;
		  }
		  // Set ownership for collections
		  if (hobj->InheritsFrom(TCollection::Class())) {
		    ((TCollection*)hobj)->SetOwner();
		  }
		  hobj->ResetBit(kMustCleanup);
		  listH.Add(hobj);
		  Int_t error = 0;
		  obj->Execute("Merge", listHargs.Data(), &error);
		  info.fIsFirst = kFALSE;
		  if (error) {
		    Error("MergeRecursive", "calling Merge() on '%s' with the corresponding object in '%s'",
			  obj->GetName(), nextsource->GetName());
		  }
		  listH.Delete();
		}
	      }
	      nextsource = (TFile*)sourcelist->After( nextsource );
	    }
	    // Merge the list, if still to be done
	    if (info.fIsFirst) {
	      Int_t error = 0;
	      obj->Execute("Merge", listHargs.Data(), &error);
	      info.fIsFirst = kFALSE;
	      listH.Delete();
	    }
	  }
	} else if (obj->InheritsFrom(TObject::Class()) &&
		   obj->IsA()->GetMethodWithPrototype("Merge", "TCollection*") ) {
	  // Object implements Merge(TCollection*) and has a reflex dictionary ...
	  
	  TList listH;
	  TString listHargs;
	  listHargs.Form("((TCollection*)0x%lx)", (ULong_t)&listH);
	  
	  // Loop over all source files and merge same-name object
	  TFile *nextsource = current_file ? (TFile*)sourcelist->After( current_file ) : (TFile*)sourcelist->First();
	  if (nextsource == 0) {
	    // There is only one file in the list
	    Int_t error = 0;
	    obj->Execute("Merge", listHargs.Data(), &error);
	    if (error) {
	      Error("MergeRecursive", "calling Merge() on '%s' with the corresponding object in '%s'",
		    obj->GetName(), key->GetName());
	    }
	  } else {
	    while (nextsource) {
	      // make sure we are at the correct directory level by cd'ing to path
	      TDirectory *ndir = nextsource->GetDirectory(path);
	      if (ndir) {
		ndir->cd();
		TKey *key2 = (TKey*)ndir->GetListOfKeys()->FindObject(key->GetName());
		if (key2) {
		  TObject *hobj = key2->ReadObj();
		  if (!hobj) {
		    Info("MergeRecursive", "could not read object for key {%s, %s}; skipping file %s",
			 key->GetName(), key->GetTitle(), nextsource->GetName());
		    nextsource = (TFile*)sourcelist->After(nextsource);
		    continue;
		  }
		  // Set ownership for collections
		  if (hobj->InheritsFrom(TCollection::Class())) {
		    ((TCollection*)hobj)->SetOwner();
		  }
		  hobj->ResetBit(kMustCleanup);
		  listH.Add(hobj);
		  Int_t error = 0;
		  obj->Execute("Merge", listHargs.Data(), &error);
		  info.fIsFirst = kFALSE;
		  if (error) {
		    Error("MergeRecursive", "calling Merge() on '%s' with the corresponding object in '%s'",
			  obj->GetName(), nextsource->GetName());
		  }
		  listH.Delete();
		}
	      }
	      nextsource = (TFile*)sourcelist->After( nextsource );
	    }
	    // Merge the list, if still to be done
	    if (info.fIsFirst) {
	      Int_t error = 0;
	      obj->Execute("Merge", listHargs.Data(), &error);
	      info.fIsFirst = kFALSE;
	      listH.Delete();
	    }
	  }
	} else {
	  // Object is of no type that we can merge
	  Bool_t warned = kFALSE;
	  
	  // Loop over all source files and write similar objects directly to the output file
	  TFile *nextsource = current_file ? (TFile*)sourcelist->After( current_file ) : (TFile*)sourcelist->First();
	  while (nextsource) {
	    // make sure we are at the correct directory level by cd'ing to path
	    TDirectory *ndir = nextsource->GetDirectory(path);
	    if (ndir) {
	      ndir->cd();
	      TKey *key2 = (TKey*)ndir->GetListOfKeys()->FindObject(key->GetName());
	      if (key2) {
		if (warned) {
		  Warning("MergeRecursive", "cannot merge object type (n:'%s', t:'%s') - "
			  "Merge(TCollection *) not implemented",
			  obj->GetName(), obj->GetTitle());
		  warned = kTRUE;
		}
		TObject *nobj = key2->ReadObj();
		if (!nobj) {
		  Info("MergeRecursive", "could not read object for key {%s, %s}; skipping file %s",
		       key->GetName(), key->GetTitle(), nextsource->GetName());
		  nextsource = (TFile*)sourcelist->After(nextsource);
		  continue;
		}
		nobj->ResetBit(kMustCleanup);
		if (target->WriteTObject(nobj, key2->GetName(), "SingleKey") <= 0) {
		  Warning("MergeRecursive", "problems copying object (n:'%s', t:'%s') to output file ",
			  obj->GetName(), obj->GetTitle());
		  status = kFALSE;
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
	// by "target->SaveSelf()" below
	target->cd();
	
	oldkeyname = key->GetName();
	//!!if the object is a tree, it is stored in globChain...
	if(obj->IsA()->InheritsFrom( TDirectory::Class() )) {
	  //printf("cas d'une directory\n");
	} else if (obj->IsA()->InheritsFrom( TCollection::Class() )) {
	  if ( obj->Write( oldkeyname, TObject::kSingleKey | TObject::kOverwrite ) <= 0 ) {
	    status = kFALSE;
	  }
	  ((TCollection*)obj)->SetOwner();
	} else {
	  if ( obj->Write( oldkeyname, TObject::kOverwrite ) <= 0) {
	    status = kFALSE;
	  }
	}
	if (obj->IsA()->InheritsFrom(TCollection::Class())) ((TCollection*)obj)->Delete();
	delete obj;
	info.Reset();
      } // while ( ( TKey *key = (TKey*)nextkey() ) )
    }
    current_file = current_file ? (TFile*)sourcelist->After(current_file) : (TFile*)sourcelist->First();
    if (current_file) {
      current_sourcedir = current_file->GetDirectory(path);
    } else {
      current_sourcedir = 0;
    }
  }
  // save modifications to the target directory.
  if (!(type&TFileMerger::kIncremental)) {
    // In case of incremental build, we will call Write on the top directory/file, so we do not need
    // to call SaveSelf explicilty.
    target->SaveSelf(kTRUE);
  }
  return status;
}

