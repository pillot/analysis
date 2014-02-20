#include "Riostream.h"
#include "TUrl.h"
#include "TSystem.h"
#include "TString.h"
#include "TFile.h"
#include "TGrid.h"

void CopyFromRemote(const char* txtfile="saf.aods.txt")
{
  char line[1024];
  ifstream in(gSystem->ExpandPathName(txtfile));
  
  while ( in.getline(line,1024,'\n') )
  {
  	TUrl url(line);
  	
    if ( TString(url.GetProtocol()) == "alien" ) 
    {
      if (!gGrid) 
      {
        TGrid::Connect("alien://");
        if (!gGrid) 
        {
          cout << "Cannot get alien connection" << endl;
          return;          
        }
      }
      
    }
  	TString file(url.GetFile());
  	
  	TString dir(gSystem->DirName(file));
  	
  	gSystem->mkdir(dir.Data(),kTRUE);

    if ( gSystem->AccessPathName(file.Data())==kFALSE)
    {
      cout << "Skipping copy of " << file.Data() << " as it already exists" << endl;
    }
    else
    {
  	  TFile::Cp(line,file.Data());
      if ( TString(line).Contains("root_archive.zip") )
      {
        gSystem->Exec(Form("unzip %s -d %s",file.Data(),gSystem->DirName(file.Data())));
        gSystem->Exec(Form("rm %s",file.Data()));
      }
  	}
  }
}