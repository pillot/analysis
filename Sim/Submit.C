/*
 *  Submit.C
 *  aliroot
 *
 *  Created by Philippe Pillot on 17/07/10.
 *  Copyright 2010 SUBATECH. All rights reserved.
 *
 */

// JPsi
Double_t nGenEvents = 4000000.;
//Double_t nGenFactor = -1.;
//Double_t nGenFactor = 0.0008; // LHC13de
//Double_t nGenFactor = 0.0012; // LHC13f
//Double_t nGenFactor = 0.5; // LHC12hi
//Double_t nGenFactor = 0.65; // LHC15n
//Double_t nGenFactor = 0.3; // LHC15o
//Double_t nGenFactor = 0.04; // LHC15o
Double_t nGenFactor = 0.22; // LHC17pq
Int_t maxEventsPerChunk = 5000;

// Beauty
//Int_t nGenEvents = 16932452;
//Int_t maxEventsPerChunk = 5000;

Bool_t FileExists(const char *lfn);
Bool_t DirectoryExists(const char *dirname);

//______________________________________________________________________________
void Submit(const char* outDir, TString jdl = "run.jdl", const char* runList = "runListProd.txt", Bool_t submit = kFALSE)
{
  /// Submit multiple production jobs with the format "submit jdl 000run#.xml 000run#".
  /// Example:
  /// - outDir = "Sim/LHC10h/JPsiPbPb276/AlignRawVtxRaw/ESDs"
  /// - runListProd.txt must contains the list of run number and the associated number of interesting events
  
  TString homeDir = "/alice/cern.ch/user/p/ppillot";
  
  if (!TGrid::Connect("alien://")) {
    Error("Submit","cannot connect to grid");
    return;
  }
  
  TString outputDir = Form("%s/%s", homeDir.Data(), outDir);
  
  if (!DirectoryExists(outputDir.Data())) {
    Error("Submit", "directory %s does not exist", outputDir.Data());
    return;
  }
  
  gGrid->Cd(outputDir.Data());
  
  if (!FileExists(jdl.Data())) {
    Error("Submit", "file %s does not exist in %s", jdl.Data(), outputDir.Data());
    return;
  }
  
   
  ifstream inFile(runList);
  if (!inFile.is_open()) {
    Error("Submit", "cannot open file %s", runList);
    return;
  }
  
  
  Double_t input[2][1000];
  Int_t nRuns = 0;
  Double_t totEvt = 0.;
  TString line;
  while (! inFile.eof() ) {
    
    line.ReadLine(inFile,kTRUE);
    if(line.IsNull()) continue;
    
    TObjArray *param = line.Tokenize(" ");
    if (param->GetEntries() != 2) {
      Error("Submit", "bad input line %s", line.Data());
      continue;
    }
     
    input[0][nRuns] = ((TObjString*)param->UncheckedAt(0))->String().Atof();
    input[1][nRuns] = ((TObjString*)param->UncheckedAt(1))->String().Atof();
    
    totEvt += input[1][nRuns];
    nRuns++;
    
    delete param;
  }
  
  inFile.close();
  
  
  Double_t ratio = (nGenFactor > 0.) ? nGenFactor : nGenEvents / totEvt;
  cout << endl;
  cout << "total number of selected events = " << totEvt << endl;
  if (nGenFactor <= 0.) cout << "required number of generated events = " << nGenEvents << endl;
  cout << "number of generated events per event = " << ratio << endl;
  cout << endl;
  cout << "run\tchunks\tevents" << endl;
  cout << "----------------------" << endl;
  
  Double_t nEvtRun = 0.;
  Int_t nChunk = 0;
  Int_t nEvtChunk = 0;
  Int_t nJobs = 0;
  Int_t nEvts = 0;
  for (int i=0; i<nRuns; i++){
    
    nEvtRun = ratio * input[1][i];
    nChunk = 1;  
    while (nEvtRun/nChunk+0.5 >= maxEventsPerChunk+1) nChunk++;
    nEvtChunk = (Int_t) (nEvtRun/nChunk + 0.5);
    nJobs += nChunk;
    nEvts += nChunk*nEvtChunk;
    cout << input[0][i] << "\t" << nChunk << "\t" << nEvtChunk << endl;
    
    TString query = Form("submit %s %d %d %d", jdl.Data(), ((Int_t) input[0][i]), nChunk, nEvtChunk);
    printf("%s ...", query.Data());
    fflush(stdout);
    
    if (! submit) {
      cout << endl;
      continue;
    }
    
    TGridResult *res = gGrid->Command(query);
    
    if (res) {
      
      TString cjobId1 = res->GetKey(0,"jobId");
      
      if (!cjobId1.Length()) {
		printf(" FAILED\n\n");
		gGrid->Stdout();
		gGrid->Stderr();
      } else printf(" DONE\n   --> the job Id is: %s \n\n", cjobId1.Data());
      
      delete res;
      
    } else printf(" FAILED\n\n");
    
  }
  
  cout << endl;
  cout << "total number of jobs = " << nJobs << endl;
  cout << "total number of generated events = " << nEvts << endl;
  cout << "difference compared to expected = " << nEvts - ratio*totEvt
       << " (" << 100.*(nEvts-ratio*totEvt)/(ratio*totEvt) << " %)" << endl;
  cout << endl;
  
}

//______________________________________________________________________________
Bool_t FileExists(const char *lfn)
{
  // Returns true if file exists.
  if (!gGrid) return kFALSE;
  TGridResult *res = gGrid->Ls(lfn);
  if (!res) return kFALSE;
  TMap *map = dynamic_cast<TMap*>(res->At(0));
  if (!map) {
    delete res;
    return kFALSE;
  }   
  TObjString *objs = dynamic_cast<TObjString*>(map->GetValue("name"));
  if (!objs || !objs->GetString().Length()) {
    delete res;
    return kFALSE;
  }
  delete res;   
  return kTRUE;
}

//______________________________________________________________________________
Bool_t DirectoryExists(const char *dirname)
{
  // Returns true if directory exists. Can be also a path.
  if (!gGrid) return kFALSE;
  // Check if dirname is a path
  TString dirstripped = dirname;
  dirstripped = dirstripped.Strip();
  dirstripped = dirstripped.Strip(TString::kTrailing, '/');
  TString dir = gSystem->BaseName(dirstripped);
  dir += "/";
  TString path = gSystem->DirName(dirstripped);
  TGridResult *res = gGrid->Ls(path, "-F");
  if (!res) return kFALSE;
  TIter next(res);
  TMap *map;
  TObject *obj;
  while ((map=dynamic_cast<TMap*>(next()))) {
    obj = map->GetValue("name");
    if (!obj) break;
    if (dir == obj->GetName()) {
      delete res;
      return kTRUE;
    }
  }
  delete res;
  return kFALSE;
}      

