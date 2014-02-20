//______________________________________________________________________________
Int_t GetRunNumber(TFileCollection *fc)
{
  /// extract the run number from the urls in the collection
  
  Int_t run = -1, irun = -1;
  
  TIter next(fc->GetList());
  TFileInfo *fi = 0x0;
  while ((fi = static_cast<TFileInfo*>(next()))) {
    
    TString fileName = gSystem->BaseName(fi->GetCurrentUrl()->GetFile());
    
    TString srun = fileName(5,6);
    if (srun.IsDec()) irun = srun.Atoi();
    else {
      Error("GetRunNumber", Form("cannot extract run number from file name \"%s\"", filename.Data()));
      return -1;
    }
    
    if (run >= 0 && irun != run) {
      Error("GetRunNumber", Form("different run numbers extracted from file names (%d vs %d)", irun, run));
      return -1;
    }
    
    run = irun;
    
  }
  
  return run;
  
}

//______________________________________________________________________________
TString CreateEntryList(TFileCollection *fc)
{
  
  TChain chain("RAW");
  if (!chain.AddFileInfoList((TCollection*)(fc->GetList()))) {
    Error("SetEntryList","Bad file list in collection, the chain is empty");
    return "";
  }
  
  TEntryList el;
  el.Enter(129, &chain);
  el.Enter(131, &chain);
  el.Enter(141, &chain);
  el.Enter(165, &chain);
  el.Enter(167, &chain);
  el.Enter(276, &chain);
  el.Enter(1081, &chain);
  el.Enter(1219, &chain);
  el.Enter(1256, &chain);
  el.Enter(1435, &chain);
  el.Enter(1542, &chain);
  el.Enter(2004, &chain);
  el.Enter(2284, &chain);
  el.Enter(3451, &chain);
  el.Enter(3558, &chain);
  el.Enter(3789, &chain);
  el.Enter(3903, &chain);
  el.Enter(3978, &chain);
  el.Enter(4185, &chain);
  el.Enter(4365, &chain);
  el.Enter(4401, &chain);
  el.Enter(4434, &chain);
  el.Enter(4593, &chain);
  el.Enter(4950, &chain);
  el.Enter(5014, &chain);
  el.Enter(5064, &chain);
  el.Enter(5208, &chain);
  el.Enter(5215, &chain);
  el.Enter(5263, &chain);
  el.Enter(5295, &chain);
  el.Enter(5602, &chain);
  el.Enter(5904, &chain);
  el.Enter(6382, &chain);
  el.Enter(6422, &chain);
  el.Enter(6479, &chain);
  el.Enter(6573, &chain);
  el.Enter(6899, &chain);
  el.Enter(7016, &chain);
  el.Enter(7243, &chain);
  el.Enter(7434, &chain);
  el.Enter(7532, &chain);
  el.Enter(8064, &chain);
  el.Enter(8109, &chain);
  el.Enter(8140, &chain);
  el.Enter(8642, &chain);
  el.Enter(8651, &chain);
  el.Enter(8678, &chain);
  el.Enter(8705, &chain);
  
  el.Print("all");
  if (!el.GetLists()) {
    Error("SetEntryList","must have several files selected for this to work");
    return "";
  }
  
  TFile *f = TFile::Open("collection.root", "RECREATE");
  el.Write();
  f->Close();
  
  return "collection.root";
  
}

//______________________________________________________________________________
void runLocal(TString input, TString cdbSnapshot = "CDBMirror", Int_t run = -1)
{
  
  gROOT->LoadMacro("$WORK/Macros/Facilities/runTaskFacilities.C");
  
  // copy files needed for the reconstruction
  TList pathList; pathList.SetOwner();
  pathList.Add(new TObjString("$WORK/Macros/Reconstruction/saf"));
  TList fileList; fileList.SetOwner();
  fileList.Add(new TObjString("runLocal.C"));
  fileList.Add(new TObjString("rec.C"));
  fileList.Add(new TObjString("UpdateCDBSnapshot.C"));
  CopyFileLocally(pathList, fileList);
  
  // check input
  if (input.EndsWith(".root")) {
    
    // check if the run number is provided for CDB snapshot if needed
    if (!cdbSnapshot.IsNull() && run < 0) {
      Error("runLocal", "a run number is needed to update CDB snapshot");
      return;
    }
    Info("runLocal", Form("Will update/use local CDB snapshot %s",cdbSnapshot.Data()));
    
  } else {
    
    // control if the collection is correct
    TFileCollection *fc = new TFileCollection();
    if (fc->AddFromFile(input.Data()) <= 0) {
      Error("runLocal", "input must be either a root file or a text file containing a list of root files");
      return;
    }
    
    // extract the run number for CDB snapshot if needed
    if (!cdbSnapshot.IsNull()) {
      Int_t run2 = GetRunNumber(fc);
      if (run < 0) {
	if (run2 < 0) {
	  Error("runLocal", "cannot extract a run number from file names in the collection. Needed to update CDB snapshot");
	  return;
	} else run = run2;
      } else if (run2 >= 0 && run2 != run) {
	Error("runLocal", "run number extracted from file names in the collection different from the given one");
	return;
      }
      Info("runLocal", Form("Will update/use local CDB snapshot %s",cdbSnapshot.Data()));
    }
    
    // create a collection of events
    input = CreateEntryList(fc);
    if (input.IsNull()) return;
    
    delete fc;
    
  }
  
  // Run reconstruction
  gROOT->LoadMacro("rec.C");
  gROOT->ProcessLine(Form("rec(\"%s\",\"\",\"%s\",%d);", input.Data(), cdbSnapshot.Data(), run));
  
}

