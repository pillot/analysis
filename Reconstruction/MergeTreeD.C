// #include <TFile.h>
// #include <TTree.h>
// #include <TError.h>

// #include "AliRunLoader.h"
// #include "AliLoader.h"

// #include "AliMUONVDigitStore.h"

//-----------------------------------------------------------------------
void MergeTreeD(const char *outFileName = "digits.root")
{
  /// The reconstruction produces one TreeD per event, also containing trigger info
  /// Merge them into one single TreeD and keep only the digit stores
  /// To read the digits in MUON.Digits.root one needs galice.root too

  // prepare to read the digits
  AliRunLoader *rl = AliRunLoader::Open("galice.root", "MUONLoader");
  AliLoader *muonLoader = rl->GetDetectorLoader("MUON");
  if (muonLoader->LoadDigits("READ") != 0) {
    Error("MergeTreeD", "unable to load digits");
    return;
  }
  AliMUONVDigitStore* digitStore = AliMUONVDigitStore::Create(*(muonLoader->TreeD()));

  // prepare to write the digits in the new file
  TFile *outFile = TFile::Open(outFileName,"RECREATE");
  TTree *outTreeD = new TTree("TreeD","Digit Container");

  int nEvents = rl->GetNumberOfEvents();
  for (int iEvent = 0; iEvent < nEvents; ++iEvent) {

    // load the current event
    if (!(rl->GetEvent(iEvent) == 0)) {
      Error("MergeTreeD", "unable to load event %d", iEvent);
      return;
    }

    // get the digit store
    TTree* inTreeD = muonLoader->TreeD();
    digitStore->Clear();
    digitStore->Connect(*inTreeD);
    inTreeD->GetEvent(0);

    // write it in the new tree
    digitStore->Connect(*outTreeD);
    outTreeD->Fill();
  }

  outFile->cd();
  outTreeD->Write();
  outFile->Close();
}
