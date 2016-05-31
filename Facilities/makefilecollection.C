#include "Riostream.h"
#include "TProof.h"
#include "TFileCollection.h"
#include "TFile.h"

void makefilecollection(const char* dslist)
{
  if (!gProof)
  {
    std::cerr << "Need a Proof connection first !" << std::endl;
    return;
  }
  std::ifstream in(dslist);
  
  std::string line;
  
  TFileCollection* fc = 0x0;
  while (std::getline(in,line))
  {
    TString name(line);
    
    //		if (name.Contains("245343")) continue;
    
    name.ReplaceAll("/","_");
    name.ReplaceAll("=","_");
    name.ReplaceAll(";","_");
    name.ReplaceAll("*","X");
    name += ".root";
    
    //		line += ";ForceUpdate";
    
    if (fc) fc->Add(gProof->GetDataSet(line.c_str()));
    else fc = gProof->GetDataSet(line.c_str());
    
  }
  
  TFile* f = TFile::Open("ds.root","recreate");
  fc->Write("dataset");
  delete f;
  
}