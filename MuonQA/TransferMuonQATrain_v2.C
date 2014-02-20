//--------------------------------------------------------------------------
// Macro for QA monitoring.
//
// In case it is not run with full aliroot, it needs the following libraries to compile:
//  - libSTplotEERBase.so
//  - libESD.so
//  - libAOD.so
//  - libANALYSIS.so
//  - libANALYSISalice.so
//  - libCORRFW.so
//  - libPWG3muon.so
//
// The macro reads results produced by the PWG1 QA train and produce monitoring plots.
//
// Author: Cynthia Hadjidakis - IPN Orsay
//--------------------------------------------------------------------------

#if !defined(__CINT__) || defined(__MAKECINT__)

#include <Riostream.h>

// ROOT includes
#include "TMath.h"
#include "TGrid.h"
#include "TFile.h"
#include "TH1.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TPad.h"

// ALIROOT includes
#include "../PWG3/base/AliCounterCollection.h"

#endif

TObjArray * GetListOfFiles(const char* baseDir, const char * trainName, const char* inFile);
TObjArray * GetListOfRuns(const char* runList, TObjArray *&listoffiles);

// .x TransferMuonQATrain_v2.C("alien:///alice/data/2010/LHC10e","QA50","/Users/cynthia/Documents/alice/data/MuonQA/LHC10e/pass2/runlist_period3_test3.txt")
//--------------------------------------------------------------------------
Bool_t TransferMuonQATrain_v2(const char* baseDir, const char * trainName, const char* runList)
{
  #if defined(__CINT__) && !defined(__MAKECINT__)
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
  #endif
  
  // Cosmetics and configuration
  gStyle->SetFillColor(10);
  gStyle->SetPadGridX(kTRUE);
  gStyle->SetPadGridY(kTRUE);
  gStyle->SetPadRightMargin(0.01);

  TString LHCPeriod = "LHC10e"; 

  TString OutFileName = "QA_";  OutFileName += LHCPeriod;
  TString OutFileNamePDF=  OutFileName.Data();  OutFileNamePDF+= ".pdf";
  TString OutFileNamePDF_open = OutFileNamePDF.Data(); OutFileNamePDF_open += "[";  
  TString OutFileNamePDF_close= OutFileNamePDF.Data(); OutFileNamePDF_close += "]";  
  TString OutFileNameROOT=  OutFileName.Data();  OutFileNameROOT+= ".root";

  if (0){ // Equivalent to the fast read option
    gEnv->SetValue("XNet.ConnectTimeout",10);
    gEnv->SetValue("XNet.RequestTimeout",10);
    gEnv->SetValue("XNet.MaxRedirectCount",2);
    gEnv->SetValue("XNet.ReconnectTimeout",10);
    gEnv->SetValue("XNet.FirstConnectMaxCnt",1);
  }

	const char* inputFile = "QAresults.root";

	TString sbaseDir = baseDir;
	if (sbaseDir.Contains("alien:") && !TGrid::Connect("alien://")) {
    Error("MuonQATrain","cannot connect to grid");
    return 0;
  }	

  //----------------------------------------------------------- //
  //          Build the list of files, the list of runs         //
	//					to be processed
  //----------------------------------------------------------- //
	
	TObjArray *listoffiles = (TObjArray*) GetListOfFiles(baseDir,trainName,inputFile);
	if(!listoffiles) return kFALSE;
  TObjArray *runs = (TObjArray*) GetListOfRuns(runList,listoffiles);
	if(!runs||!listoffiles){
		Error("TransferMuonQATrain","cannot get a list of selected runs or files");
		return kFALSE;
	}
	printf("TransferMuonQATrain: Files found to be processed = %d\n",runs->GetEntriesFast());
	if(runs->GetEntriesFast()==0)	return kFALSE;

	//-------------------------------------//
	//    Loop over the files on grid      //
	//-------------------------------------//

	TFile *file = 0;
	TIter next0(listoffiles);
	TObject *nextfile;
	TString snextfile;
	Int_t nFiles = -1;
	AliCounterCollection *trackCounters = 0;
	AliCounterCollection *eventCounters = 0;
	AliCounterCollection *mergedTrackCounters = 0;
  AliCounterCollection *mergedEventCounters = 0;	
	TObjArray*  outputList; 
	TObjArray*  outputListExpert;
	TObjArray*  outputListNorm;  
	TFile	outputHistoFile;	
	TString fileName;

	TString command = "mkdir results";
	cout<<"Shell command = "<<command<<endl;
	gSystem->Exec(command.Data());
		
	while ((nextfile=next0())) {

		nFiles++;
		snextfile = nextfile->GetName();
		//Open the file
		file = TFile::Open(snextfile.Data());
		if(!file) continue;
		
		outputList = (TObjArray *)file->Get("MUON_QA/general1");
		outputListExpert = (TObjArray *)file->Get("MUON_QA/expert");
		outputListNorm = (TObjArray *)file->Get("MUON_QA/general2");
		
		trackCounters = (AliCounterCollection *) file->Get("MUON_QA/trackCounters");
		eventCounters = (AliCounterCollection *) file->Get("MUON_QA/eventCounters");
		
		if(!outputList || !outputListExpert || !outputListNorm || !trackCounters || !eventCounters){
		  Error("TransferMonQATrain","Object not found for that file");
			continue;
		}

		//-------------------------------------//
		//    Merge the AliCounterCollection
		//-------------------------------------//
		if(nFiles==0){
			mergedTrackCounters = (AliCounterCollection*) trackCounters->Clone();
			mergedEventCounters = (AliCounterCollection*) eventCounters->Clone();
		}		
		else{
			mergedTrackCounters->Add(trackCounters);
			mergedEventCounters->Add(eventCounters);
		}
	

		//-------------------------------------//
		//     Save the MUONQA histos in AnalysisResults.root for each run number
		//-------------------------------------//
		fileName = "AnalysisResults.root";
		outputHistoFile.Open(fileName,"recreate");
		TDirectoryFile *dir = new TDirectoryFile("MUON_QA","MUON_QA");
		outputHistoFile.Cd("MUON_QA");

		TDirectory *dir0 = (TDirectory*) file->GetDirectory("MUON_QA");
		TObjArray* general1 = static_cast<TObjArray*>(dir0->FindObjectAny("general1"));
		TObjArray* expert = static_cast<TObjArray*>(dir0->FindObjectAny("expert"));
		TObjArray* general2 = static_cast<TObjArray*>(dir0->FindObjectAny("general2"));
		
		if(general1) general1->Write("general1",TObject::kSingleKey);
		if(general2) general2->Write("general2",TObject::kSingleKey);		
		if(expert) expert->Write("expert",TObject::kSingleKey);
		outputHistoFile.Close();
		
		command = "mkdir -p results/";
		command+=((TObjString*)runs->UncheckedAt(nFiles))->GetString();
		cout<<"Shell command = "<<command<<endl;
		gSystem->Exec(command.Data());
		command=" mv AnalysisResults.root results/";
		command+=((TObjString*)runs->UncheckedAt(nFiles))->GetString();
		command+= "/.";
		cout<<"Shell command = "<<command<<endl;
		gSystem->Exec(command.Data());
  } //end of loop over files
	
	//-------------------------------------//
	//      Save the AliCounterCollection in MergedAnalysisResults.root
	//-------------------------------------//

	TFile outputFile("MergedAnalysisResults.root","recreate");
	TDirectoryFile *dir = new TDirectoryFile("MUON_QA","MUON_QA");
	outputFile.Cd("MUON_QA");

	mergedTrackCounters->Write();
	mergedEventCounters->Write();
	outputFile.Close();
	
	command = "mv MergedAnalysisResults.root results/.";
	cout<<"Shell command = "<<command<<endl;
	gSystem->Exec(command.Data());
	
	return kTRUE;
}

TObjArray * GetListOfRuns(const char* runList, TObjArray *&listoffiles)
{

	TObjArray * runs = new TObjArray();
  runs->SetOwner();
	
	if(!runList){
		Error("GetListOfruns","runList is not defined... exit");
		return 0;
	}
	else {
   // only the ones in the runList
    ifstream inFile(runList);
    if (!inFile.is_open()) {
      Error("GetListOfRuns",Form("unable to open file %s", runList));
      return 0;
    }
    
    TString currRun;
    while (!inFile.eof()) {
      currRun.ReadLine(inFile, kTRUE);
      if (currRun.IsNull()) continue;
      if (!currRun.IsDigit()) {
				Error("GetListOfRuns","invalid run number: %s", currRun.Data());
				return 0;
      }
      runs->AddLast(new TObjString(Form("%09d", currRun.Atoi())));
		}
    
    inFile.close();
	}
	
	printf("GetListOfRuns: Nr of runs in the runlist = %d \n",runs->GetEntriesFast());

	if(runList && listoffiles){	
		//Filter the selected runs and modify listoffiles
		TObjArray*  runsFound = new TObjArray();
		runsFound->SetOwner();	

		//filter the selected runs	
		TIter next0(listoffiles);
		TObject *nextfile;
		TString snextfile;

		TObjArray *listoffilestmp = new TObjArray();	
		listoffilestmp->SetOwner();
		while ((nextfile=next0())) {//loop over files found on alien
			snextfile = nextfile->GetName();
			for ( Int_t irun=0; irun<runs->GetEntriesFast(); irun++ ) { //loop over selected runs
				TString run = ((TObjString*)runs->UncheckedAt(irun))->GetString();
				if(snextfile.Contains(run)){
					listoffilestmp->Add(nextfile);
					runsFound->AddLast(new TObjString(Form("%09d", run.Atoi())));
				}
			}
		}
		runs = runsFound;
		listoffiles->Clear();
		listoffiles = (TObjArray*) listoffilestmp->Clone();
		
		printf("GetListOfRuns Nr of selected runs corresponding to the list of files = %d \n",runs->GetEntriesFast());		
	}

return runs;
		
}

TObjArray* GetListOfFiles(const char* baseDir, const char * trainName, const char* inFile)
{

	TString sbaseDir = baseDir;
	TString strainName = trainName;
	TString inputFile = inFile;
	TString command;
	
	if(!sbaseDir.Contains("alien://")){
		Error("GetListOfFiles","Not implemented for files not on alien-->exit");
		return 0;
	}
	
	sbaseDir.ReplaceAll("alien://", "");

	TObjArray *listoffiles = new TObjArray();
	
	if (sbaseDir.Contains(".xml")) {
	// Read files pointed by the xml 
      TGridCollection *coll = (TGridCollection*)gROOT->ProcessLine(Form("TAlienCollection::Open(\"%s\");", sbasedir));
      if (!coll) {
         ::Error("GetListOfFiles", "Input XML collection empty.");
         return 0;
      }
      // Iterate grid collection
      while (coll->Next()) {
         TString fname = gSystem->DirName(coll->GetTURL());
         fname += "/";
         fname += inputFile;      
         listoffiles->Add(new TNamed(fname.Data(),""));
      }   
	}
	else {   
      command = Form("find %s/ *%s/%s", sbaseDir.Data(), strainName.Data(), inputFile.Data());
      printf("command: %s\n", command.Data());
      TGridResult *res = gGrid->Command(command);
      if (!res) {
         ::Error("GetListOfFiles","No result for the find command\n");
         delete listoffiles;
         return 0;
      }     
      TIter nextmap(res);
      TMap *map = 0;
      while ((map=(TMap*)nextmap())) {
         TObjString *objs = dynamic_cast<TObjString*>(map->GetValue("turl"));
         if (!objs || !objs->GetString().Length()) {
            // Nothing found - skip this output
            delete res;
            delete listoffiles;
            return 0;
         }
         listoffiles->Add(new TNamed(objs->GetName(),""));
      }
      delete res;
	}
	if (!listoffiles->GetEntries()) {
      ::Error("GetListOfFiles","No files from the find command=%s\n",command.Data());
      delete listoffiles;
      return 0;
   }     
	else printf("GetListOfFiles: Number of files found %d\n",listoffiles->GetEntries());

	return listoffiles;

}
