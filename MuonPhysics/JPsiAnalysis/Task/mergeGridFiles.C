#if !defined(__CINT__) || defined(__MAKECINT__)

#include <Riostream.h>

// ROOT includes
#include "TString.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TFileMerger.h"
#include "TGrid.h"
#include "TList.h"
#include "TSystem.h"
#include "TGridResult.h"
#include "TMap.h"
#endif

Bool_t AddFileList(TString, TString&, Bool_t);
void ReadListFromFile(TString, TString&, Bool_t);
Int_t GetLastStage(TGridResult*);

void mergeGridFiles(TString outFilename, TString inFileList, Bool_t isLocal = kFALSE, Int_t nFilesPerStep = 100, Bool_t copyLocal = kTRUE)
{
  TString fileList = "";
  ReadListFromFile(inFileList, fileList, isLocal);

  if ( fileList.IsNull() ) {
    printf("List of merging files is null: Nothing done!\n");
    return;
  }

  TString logFilename = "toBeMerged.txt";

  if ( ! isLocal && ! gGrid )
    TGrid::Connect("alien://");

  if ( isLocal ) copyLocal = kFALSE;

  TString currList = fileList;
  TString currOutput = "", mergedFiles = "";
  TFileMerger* fileMerger = 0x0;

  for ( Int_t istep = 0; istep < 100; istep++ ) {
    TObjArray* array = currList.Tokenize(" ");
    currList = "";
    Int_t nFiles = array->GetEntries();
    Int_t subStep = -1;
    for (Int_t ifile = 0; ifile < nFiles; ifile++ ) {
      if ( ! fileMerger )
	fileMerger = new TFileMerger(copyLocal);
      fileMerger->AddFile(array->At(ifile)->GetName());
      Int_t nFilesToMerge = fileMerger->GetMergeList()->GetEntries();
      if ( nFilesToMerge % nFilesPerStep != 0 && ifile < nFiles - 1 ) 
	continue;
      // The following part is executed only at the end of each step
      currOutput = outFilename;
      if ( nFiles > nFilesPerStep ) {
	subStep++;
	currOutput.ReplaceAll(".root",Form("_%i_%i.root", istep, subStep));
	AddFileList(currOutput, currList, kTRUE);
      }
      fileMerger->OutputFile(currOutput.Data());
      fileMerger->Merge();
      printf("\nMerged in %s:\n", currOutput.Data());
      mergedFiles = "";
      for ( Int_t ientry=0; ientry<nFilesToMerge; ientry++ )
	mergedFiles += Form("%s ", fileMerger->GetMergeList()->At(ientry)->GetName());
      //printf("%s\n\n", mergedFiles.Data());    // pb if string too long. Upto 2048

      // remove merged files
      if ( istep > 0 )
	gSystem->Exec(Form("rm %s", mergedFiles.Data()));

      delete fileMerger;
      fileMerger = 0x0;

      // Write log file to keep trace of files to be merged
      ofstream logFile(logFilename.Data());
      TString logString = "";
      for ( Int_t jfile = ifile + 1; jfile < nFiles; jfile++ ) {
	logString += Form("%s ", array->At(jfile)->GetName());
      }
      logString.Append(currList.Data());
      logString.ReplaceAll(" ", "\n");
      logFile << logString.Data() << endl;;
      logFile.close();
    } // loop on files


    delete array;
    printf("Step %i completed!\n", istep);

    if ( nFiles <= nFilesPerStep ) break;
  } // loop on steps

  gSystem->Exec(Form("rm %s", logFilename.Data()));
}


//___________________________________________________
Bool_t AddFileList(TString filename, TString& fileList, Bool_t isLocal)
{
  if ( filename.IsNull() || ! filename.Contains(".root") ) return kFALSE;
  if ( ! isLocal && ! filename.Contains("alien") )
    filename.Prepend("alien://");
  
  if ( ! fileList.IsNull() )
    fileList.Append(" ");
  fileList.Append(filename.Data());

  return kTRUE;
}

//___________________________________________________
void ReadListFromFile(TString filename, TString& fileList, Bool_t isLocal)
{
  ifstream inFile(filename.Data());
  TString currLine = "";
  if ( inFile.is_open() ) {
    while ( ! inFile.eof() ) {
      currLine.ReadLine(inFile,kTRUE); // Read line
      AddFileList(currLine, fileList, isLocal);
    }
    inFile.close();
  }
}


//___________________________________________________
void completeProd(TString runListName="runList.txt", TString outTaskFilename="QAresults.root", TString baseDir="/alice/data/2010/LHC10h", TString prodDir = "")
{
  TString outFilename = "/tmp/completeFileList.txt";

    // Get run list from file
  ifstream inFile(runListName.Data());
  TObjArray runList;
  runList.SetOwner();
  TString currRun;
  Int_t minRun = 99999999;
  Int_t maxRun = -1;
  if (inFile.is_open()) {
    while (! inFile.eof() ) {
      currRun.ReadLine(inFile,kTRUE); // Read line
      if ( currRun.IsNull() || ! currRun.IsDigit() ) continue;
      Int_t currRunInt = currRun.Atoi();
      minRun = TMath::Min(currRunInt, minRun);
      maxRun = TMath::Max(currRunInt, maxRun);
      runList.Add(new TObjString(Form("%d", currRunInt)));
    }
    inFile.close();
  }

  outFilename.ReplaceAll(".txt", Form("_%i_%i.txt", minRun, maxRun));

  //TString filePattern[3] = {"", "Stage*/","*/"};
  TString filePattern[2] = {"","*/"};

  ofstream outFile(outFilename.Data());

  if ( outTaskFilename.Contains("QAresults.root") ) {
    TString loadLibs[3] = {"libTENDER.so", "libPWG1.so", "libPWG3base.so"};
    for ( Int_t ilib=0; ilib<3; ilib++ ) {
      Int_t exitVal = gSystem->Load(loadLibs[ilib].Data());
      if ( exitVal < 0 ) {
        printf("Please run with aliroot if you're merging QA objects!\n");
        return;
      }
    }
  }

  if ( ! gGrid )
    TGrid::Connect("alien://");

  baseDir.ReplaceAll("alien://","");

  TMap* map = 0x0;
  TString stageName = "";
  for ( Int_t irun=0; irun<runList.GetEntries(); irun++ ) {
    TString currRunString = ((TObjString*)runList.At(irun))->GetString();

    TString tmpFilename = Form("/tmp/mergeListRun%s.txt", currRunString.Data());
    ofstream tmpFile(tmpFilename.Data());
    TString mergeFilename = "";

    for ( Int_t ipattern=0; ipattern<2; ipattern++ ) {
      TString command = Form("find %s/*%s %s/%s%s", baseDir.Data(), currRunString.Data(), prodDir.Data(), filePattern[ipattern].Data(), outTaskFilename.Data());

      printf("%s\n", command.Data()); // REMEMBER TO CUT

      TGridResult* res = gGrid->Command(command);

      if ( ! res || res->GetEntries() == 0 ) continue;

      Int_t mergeStage = ( ipattern == 1 ) ? GetLastStage(res)  : -1;
      stageName = Form("Stage_%i", mergeStage);

      TIter nextmap(res);
      while ( ( map = (TMap*)nextmap() ) ) {
        // Loop 'find' results and get next LFN
        TObjString *objs = dynamic_cast<TObjString*>(map->GetValue("turl"));
        if (!objs || !objs->GetString().Length()) 
          continue;

        mergeFilename = objs->GetString();

        if ( mergeStage > 0 && ! mergeFilename.Contains(stageName.Data()) ) continue;

        tmpFile << mergeFilename.Data() << endl;
      } // loop on grid lfns

      delete res;

      tmpFile.close();

      if ( ipattern == 1 ) {
        mergeFilename = outTaskFilename;
        mergeFilename.ReplaceAll(".root", Form("_%s.root", currRunString.Data()));
        mergeGridFiles(mergeFilename, tmpFilename, kFALSE, 50);
      }

      if ( ! mergeFilename.Contains("alien://") )
        outFile << gSystem->pwd() << "/";
      outFile << mergeFilename.Data() << endl;
      gSystem->Exec(Form("rm %s", tmpFilename.Data()));
      break;
    } // loop on pattern
  } // loop on runs
  printf("\nOutput written in\n%s\n", outFilename.Data());
}


Int_t GetLastStage(TGridResult* res)
{
  Int_t lastStage = 0;

  TMap* map = 0x0;
  TIter nextmap(res);
  TString filename = "", currToken = "";
  while ( ( map = (TMap*)nextmap() ) ) {
    // Loop 'find' results and get next LFN
    TObjString *objs = dynamic_cast<TObjString*>(map->GetValue("turl"));
    if (!objs || !objs->GetString().Length()) 
      continue;

    filename = objs->GetString();
    
    if ( ! filename.Contains("Stage_") ) continue;

    TObjArray* array = filename.Tokenize("/");
    array->SetOwner();
    for ( Int_t ientry=0; ientry<array->GetEntries(); ientry++ ) {
      currToken = array->At(ientry)->GetName();
      if ( currToken.Contains("Stage_") ) {
        currToken.Remove(0,6);
        Int_t currStage = currToken.Atoi();
        lastStage = TMath::Max(currStage, lastStage);
        break;
      }
    }
    delete array;
  } // loop on grid lfns

  return lastStage;
}
