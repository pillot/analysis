#if !defined(__CINT__) || defined(__MAKECINT__)

#include <iostream>

#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TGrid.h>
#include <TSystem.h>
#include <TClonesArray.h>

#include "AliESDEvent.h"
#include "AliESDRun.h"
#include "AliESDHeader.h"
#include "AliESDZDC.h"
#include "AliESDVZERO.h"
#include "AliMultiplicity.h"

#endif

Bool_t CheckESD(const char* fileNameE="AliESDsTmp.root",const char* fileNameR="alien:///alice/data/2010/LHC10h/000137549/ESDs/pass1_4plus/10000137549001.100/AliESDs.root") {
	
	Int_t nev = 10;
	Int_t nskip = 0;
	const char* fileNameM = "AliESDs.root";
	const char *trgname;
	
	if (gSystem->Getenv("DC_ESDFILE")) {
    fileNameR = gSystem->Getenv("DC_ESDFILE");
  }
  if (gSystem->Getenv("DC_NEVENTS")) {
    nev = atoi(gSystem->Getenv("DC_NEVENTS"));
  }
	if (gSystem->Getenv("DC_TRGNAME")) {
    trgname = gSystem->Getenv("DC_TRGNAME");
  }
	if (gSystem->Getenv("DC_EEVENT")) {
		nskip = atoi(gSystem->Getenv("DC_EEVENT"));
	}
	
	if (!gGrid) TGrid::Connect("alien://");
  TFile *_fileE = TFile::Open(fileNameE);
  TFile *_fileR = TFile::Open(fileNameR);
  AliESDEvent* esdE = new AliESDEvent();
  AliESDEvent* esdR = new AliESDEvent();
	
	AliESDRun* runM = 0x0;
	AliESDHeader* headM = 0x0;
	AliESDZDC* zdcM = 0x0;
	AliESDVZERO* vzeroM = 0x0;
	AliMultiplicity *multM = 0x0;
	
	TTree* treeESDE;
  TTree* treeESDR;
  _fileE->GetObject("esdTree", treeESDE);
  _fileR->GetObject("esdTree", treeESDR);
  esdE->ReadFromTree(treeESDE);
  esdR->ReadFromTree(treeESDR);
	
	// Output ESDs
	TFile *_fileM = TFile::Open(fileNameM, "RECREATE");
	TTree *treeESDM = new TTree("esdTree", "Tree with ESD objects");
  AliESDEvent *esdM = new AliESDEvent();
  esdM->CreateStdContent();
	esdM->WriteToTree(treeESDM);
	
  Int_t nEventsE = treeESDE->GetEntries();
  Int_t nEventsR = treeESDR->GetEntries();
  Int_t nEvents = TMath::Min(nEventsE,nEventsR);
	nEvents = TMath::Min(nEvents,nev);
	Int_t lastSelectedEvent = -1;
	if (nskip>0) {
		printf("We will skip %d events\n",nskip);
		lastSelectedEvent = nskip-1;
	}
  for (Int_t iEv=0; iEv<nEvents; iEv++) {
    printf("Event %i\n",iEv);
		
    treeESDE->GetEvent(iEv);
		if (trgname) {
			for (Int_t jEv=lastSelectedEvent+1; jEv<nEventsR; jEv++) {
				treeESDR->GetEvent(jEv);
				
				TString firedTrigClasses = esdR->GetFiredTriggerClasses();
				if (firedTrigClasses.Contains(trgname)) {
					lastSelectedEvent = jEv;
					break;
				}
			}
		}	else {
			treeESDR->GetEvent(iEv);
			lastSelectedEvent = iEv;
		}
		printf("and Event %i\n",lastSelectedEvent);
		
    // Copy the esd event from the Merged (Signal) reconstrcution
		*esdM = *esdE;
		
		// Replace some branches with those from original esd event
		runM = dynamic_cast<AliESDRun *>(esdM->FindListObject("AliESDRun"));
		*runM = *(dynamic_cast<AliESDRun *>(esdR->FindListObject("AliESDRun")));
		headM = dynamic_cast<AliESDHeader *>(esdM->FindListObject("AliESDHeader"));
		*headM = *(dynamic_cast<AliESDHeader *>(esdR->FindListObject("AliESDHeader")));
		zdcM = dynamic_cast<AliESDZDC *>(esdM->FindListObject("AliESDZDC"));
		*zdcM = *(dynamic_cast<AliESDZDC *>(esdR->FindListObject("AliESDZDC")));
		vzeroM = dynamic_cast<AliESDVZERO *>(esdM->FindListObject("AliESDVZERO"));
		*vzeroM = *(dynamic_cast<AliESDVZERO *>(esdR->FindListObject("AliESDVZERO")));
		multM = dynamic_cast<AliMultiplicity *>(esdM->FindListObject("AliMultiplicity"));
		*multM = *(dynamic_cast<AliMultiplicity *>(esdR->FindListObject("AliMultiplicity")));
		
		treeESDM->Fill();
  }
	_fileM->cd();
	treeESDM->Write(treeESDM->GetName(),TObject::kOverwrite);
	_fileM->Close();
	
	// result of check
  Info("CheckESD", "check of ESD was successfull");
  return kTRUE;
}
