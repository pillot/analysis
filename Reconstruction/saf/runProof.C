//______________________________________________________________________________
Bool_t CreateCollection(TString input)
{
  /// create (or copy) the collection of urls to input raw data files
  
  if (gSystem->AccessPathName(input.Data())) {
    
    gSystem->Exec("rm -f __collection__");
    
    TFileCollection *fc = gProof->GetDataSet(input.Data());
    if (!fc) {
      Error("CreateCollection", "dataset not found");
      return kFALSE;
    }
    
    TIter next(fc->GetList());
    TFileInfo *fi = 0x0;
    while ((fi = static_cast<TFileInfo*>(next()))) {
      
      const char *url = (fi->GetCurrentUrl()) ? fi->GetCurrentUrl()->GetUrl() : 0x0;
      if (!url) {
	Warning("CreateCollection", "found TFileInfo with empty Url - ignoring");
	continue;
      }
      
      gSystem->Exec(Form("echo %s >> __collection__", url));
      
    }
    
    delete fc;
    
  } else {
    
    gSystem->Exec(Form("cp -f %s __collection__", input.Data()));
    
  }
  
  return kTRUE;
  
}

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
void runProof(TString aaf, TString input, TString output, TString cdbSnapshot = "CDBMirror",
	      TString rootVersion = "v5-34-05", TString alirootVersion = "v5-04-58-AN")
{
  
  // check runing mode
  if (aaf != "prooflite" && aaf != "saf") {
    printf("authorized modes are \"prooflite\" or \"saf\"\n");
    printf("authorized inputs are a dataset or a .txt file containing urls\n");
    return;
  }
  
  gROOT->LoadMacro("$WORK/Macros/Facilities/runTaskFacilities.C");
  
  // copy files needed for the reconstruction
  TList pathList; pathList.SetOwner();
  pathList.Add(new TObjString("$WORK/Macros/Reconstruction/saf"));
  TList fileList; fileList.SetOwner();
  fileList.Add(new TObjString("runProof.C"));
  fileList.Add(new TObjString("rec.C"));
  fileList.Add(new TObjString("UpdateCDBSnapshot.C"));
  CopyFileLocally(pathList, fileList);
  
  // prepare proof
  TString exec(gROOT->GetPath());
  LoadAlirootOnProof(aaf, rootVersion, alirootVersion, "", "", "", exec.BeginsWith("aliroot"), "REC");
  
  // check input
  Int_t run = -1;
  if (input.BeginsWith("raw://run")) {
    
    // extract the run number for CDB snapshot if needed
    if (!cdbSnapshot.IsNull()) {
      TString srun(input);
      srun.ReplaceAll("raw://run", "");
      if (srun.IsDec()) run = srun.Atoi();
    }
    
  } else {
    
    // create the collection of urls
    if (!CreateCollection(input)) return;
    
    // control if the collection is correct
    TFileCollection *fc = new TFileCollection();
    if (fc->AddFromFile("__collection__") <= 0) return;
    input = "__collection__";
    
    // extract the run number for CDB snapshot if needed
    if (!cdbSnapshot.IsNull()) run = GetRunNumber(fc);
    
    delete fc;
    
  }
  
  // check the run number for CDB snapshot if needed
  if (!cdbSnapshot.IsNull()) {
    if (run < 0) {
      Error("runProof", "cannot extract a run number from input. Needed to update CDB snapshot");
      return;
    }
    Info("runProof", Form("Will update/use local CDB snapshot %s",cdbSnapshot.Data()));
  }
  
  // check ouput
  TString dsout = "";
  if (!output.IsNull()) {
    dsout = gSystem->BaseName(output.Data());
    TString dsoutfull = Form("/%s/%s/%s", gProof->GetGroup(), gProof->GetUser(), dsout.Data());
    if (gProof->ExistsDataSet(dsoutfull.Data())) {
      Error("runProof", Form("dataset %s already exist", dsoutfull.Data()));
      return;
    } else {
      Info("runProof", Form("will save ouput in dataset %s", dsoutfull.Data()));
    }
  }
  if (aaf == "saf" && dsout.IsNull()) {
    Error("runProof", "a valid output dataset is needed to run on saf");
    return;
  }
  
  // change the packetizer strategy to try to use all workers
  // (with the adaptative packetizer (default), for some reason, it starts with 1000 events
  // per packets so if the number of events is < 1000*nWorkers not all the workers are used)
  Int_t numWorkers = gProof->GetParallel();
  gProof->SetParameter("PROOF_PacketizerStrategy", (Int_t)0);
  gProof->SetParameter("PROOF_PacketAsAFraction", (Int_t)(1234567890./numWorkers));
  
  // Run reconstruction
  gROOT->LoadMacro("rec.C");
  gROOT->ProcessLine(Form("rec(\"%s\",\"%s\",\"%s\",%d);", input.Data(), dsout.Data(), cdbSnapshot.Data(), run));
  
  // save logs
  TProof::Mgr(gProof->GetUrl())->GetSessionLogs()->Save("*","run.log");

  // Check the produced dataset
  if (!output.IsNull()) {
    TFileCollection *coll = gProof->GetDataSet(output.Data());
    if (coll) {
      Int_t nEvents = coll->GetTotalEntries("/esdTree");
      if (nEvents > 0) {
	cout << "===========================================================================" << endl;
	cout << nEvents << " events reconstructed and stored in the dataset " << output.Data() << endl;
	cout << "===========================================================================" << endl;
	cout << "The dataset is:" << endl;
	coll->Print();
	cout << "===========================================================================" << endl;
      }
    }
  }
  
}

