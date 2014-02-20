#include "TFile.h"
#include "TTree.h"
#include "TList.h"
#include "TBranch.h"
#include "TObjArray.h"
#include "TGrid.h"

void branchSize(const char* file, const char* treename="aodTree")
{
  if (TString(file).Contains("alien://"))
  {
    TGrid::Connect("alien://");
  }
  
  Int_t entry = 1;
  TFile *f1 = TFile::Open(file);

  if (!f1) return;
  
  TTree *tree = (TTree*)f1->Get(treename);

  if (!tree) return;
  
  TObjArray *tb = tree->GetListOfBranches();

  if (!tb) return;
  
  Int_t fBytes = (Int_t)tree->GetEntry(entry);
  printf("\n ###################################\n");
  printf("%20s ReadBytes %10d  ",tree->GetName(),fBytes);
  
  ULong64_t totBytes = tree->GetTotBytes();
  ULong64_t zipBytes = tree->GetZipBytes();
  
  printf("TotalBytes %4.0f KB AfterCompression %4.0f KB Nentries %lld\n",
         totBytes/1024.0,zipBytes/1024.0,tree->GetEntries());
  printf("###################################\n");
  TBranch* b;
  TIter next(tb);
  while ( ( b = static_cast<TBranch*>(next()) ) )
  { 
    Int_t bytes = b->GetEntry(entry);
    printf("%20s ReadBytes %10d  ",b->GetName(),bytes);
    printf("TotalBytes %10d (%4.0f %%) AfterCompression %10d (%4.0f %%)\n",
           (Int_t)b->GetTotBytes("*"),
           (totBytes>0) ? b->GetTotBytes("*")*100.0/totBytes : 0,
           (Int_t)b->GetZipBytes("*"),
           (zipBytes>0) ? b->GetZipBytes("*")*100.0/zipBytes : 0
           );
  }

  TBranch* br = tree->BranchRef();
  
  if (br)
  {

    printf("%20s ",br->GetName());
    printf("AfterCompression %10d (%5.2f %% of file size)\n",
           (Int_t)br->GetZipBytes("*"),
           br->GetZipBytes("*")*100.0/f1->GetSize()
           );
    
  }
  
  
  delete f1;
}
