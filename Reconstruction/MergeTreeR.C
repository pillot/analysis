// #include <TFile.h>
// #include <TTree.h>
// #include <TError.h>

// #include "AliRunLoader.h"
// #include "AliLoader.h"

// #include "AliMUONVClusterStore.h"

//-----------------------------------------------------------------------
void MergeTreeR(const char *outFileName = "clusters.root")
{
  /// The reconstruction produces one TreeR per event, also containing trigger info
  /// Merge them into one single TreeR and keep only the cluster stores
  /// To read the clusters in MUON.RecPoints.root one needs galice.root too

  // prepare to read the clusters
  AliRunLoader *rl = AliRunLoader::Open("galice.root", "MUONLoader");
  AliLoader *muonLoader = rl->GetDetectorLoader("MUON");
  muonLoader->LoadRecPoints("READ");

  // prepare to write the clusters in the new file
  TFile *outFile = TFile::Open(outFileName,"RECREATE");
  TTree *outTreeR = new TTree("TreeR","Cluster Container");

  int nEvents = rl->GetNumberOfEvents();
  AliMUONVClusterStore* clusterStore = 0x0;
  for (int iEvent = 0; iEvent < nEvents; ++iEvent) {
    
    // load the current event
    if (!(rl->GetEvent(iEvent) == 0)) {
      Error("MergeTreeR", "unable to load event %d", iEvent);
      return;
    }

    // get the cluster store
    TTree* inTreeR = muonLoader->TreeR();
    if (!clusterStore) {
      clusterStore = AliMUONVClusterStore::Create(*inTreeR);
    } else {
	    clusterStore->Clear();
    }
	  clusterStore->Connect(*inTreeR);
	  inTreeR->GetEvent(0);

    // write it in the new tree
    clusterStore->Connect(*outTreeR);
    outTreeR->Fill();
  }
  
  outFile->cd();
  outTreeR->Write();
  outFile->Close();
}
