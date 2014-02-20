//______________________________________________________________________________
void LoadLocalLibs()
{
    gSystem->Load("libVMC");
    gSystem->Load("libMinuit");
    gSystem->Load("libTree");
    gSystem->Load("libProofPlayer");
    gSystem->Load("libXMLParser");
    gSystem->Load("libSTEERBase");
    gSystem->Load("libESD");
    gSystem->Load("libAOD");
    gSystem->Load("libANALYSIS");
    gSystem->Load("libANALYSISalice");
    gSystem->Load("libPWG3base");

//gSystem->Load("$HOME/Code/Analysis/la_mumu/libla_mumu");

	gROOT->LoadMacro("AliHistogramCollection.cxx+g");
  	gROOT->LoadMacro("AliAnalysisMuMu.cxx+g");
  	gROOT->LoadMacro("AliAnalysisMuMuFromAOD.cxx+g");
}

//______________________________________________________________________________
TChain* CreateLocalChain(const char* filelist)
{
	TChain* c = new TChain("aodTree");
	
	char line[1024];
	
	ifstream in(filelist);
	while ( in.getline(line,1024,'\n') )
	{
		c->Add(line);
	}
	return c;
}


//______________________________________________________________________________
/*
void aod(const char* dataset="/MUON/laphecet/DATA_LHC10h_AOD026_000137844", 
	     const char* where="laphecet@alice-caf.cern.ch",
	     const char* workers="workers=2x")
*/
 void aod(const char* dataset="ds.aod049.once.txt", 
	     const char* where="laphecet@nansafmaster.in2p3.fr",
	     const char* workers="workers=8x")
/*
//DATA_LHC10h_AOD033_000137844
//MUON/laphecet/DATA_LHC10h_AOD033_137XXX

//void aod(const char* dataset="/default/shabetai/LHC11a_000146806_AllAODs",
//         const char* where="laphecet@nansafmaster.in2p3.fr:1093/?N",
//         const char* workers="workers=8x")

/*
void aod(const char* dataset="/MUON/laphecet/DATA_LHC10h_AOD033_137XXX",
         const char* where="laphecet@alice-caf.cern.ch:1093/?N",
         const char* workers="workers=1x")
*/

 {
  // Create the analysis manager

   TProof* p(0x0);
//   TString alirootMode("PAR");
   TString alirootMode("");
   
   Bool_t prooflite = (strlen(where)==0) || TString(where).Contains("workers");
   
   if (prooflite)
   {
     cout << "Will work in LITE mode" << endl;
   }
   
   if ( dataset )
   {
     p = TProof::Open(where,workers);
     
     if (!p)
     {
       cout << "Cannot connect to Proof : " << where << endl;
       return;
     }
     
     p->GetManager()->SetROOTVersion("VO_ALICE@ROOT::v5-28-00a");
     
     alirootMode.ToUpper();
     
     if ( alirootMode == "PAR" ) 
     {
       cout << "Will work with PAR files" << endl;
       
       std::vector<std::string> pars;
       
       pars.push_back("STEERBase");
       pars.push_back("ESD");
       pars.push_back("AOD");
       pars.push_back("ANALYSIS");
       pars.push_back("ANALYSISalice");
       pars.push_back("CORRFW");
       pars.push_back("PWG3base");
       //  pars.push_back("PWG3muon");
       
       Bool_t ok(kTRUE);
       
       for ( std::vector<std::string>::size_type i = 0; i < pars.size(); ++i )
       {
         std::string package = pars[i];
         
         if ( gProof->UploadPackage(package.c_str()) ) 
         {
           ok = kFALSE;
         }
         
         if ( gProof->EnablePackage(package.c_str(),"",kTRUE) ) 
         {
           ok = kFALSE;
         }
         
         if (!ok) 
         {
           cout << "Problem with PAR " << package.c_str() << endl;
           return;           
         }
       }
     }
     else 
     {
       TList* list = new TList();
       
       list->Add(new TNamed("ALIROOT_EXTRA_LIBS", "PWG3base"));
       list->Add(new TNamed("ALIROOT_EXTRA_INCLUDES", "PWG3/base"));
       
       //    list->Add(new TNamed("ALIROOT_ENABLE_ALIEN", "1"));
       
       if (!alirootMode.IsNull())
       {
         list->Add(new TNamed("ALIROOT_MODE", alirootMode.Data()));  
       }
       
       if (!prooflite)
       {
         p->EnablePackage("VO_ALICE@AliRoot::v4-21-23-AN", list, kTRUE);
       }
       else
       {
         //      list->Add(new TNamed("ALIROOT_LOCAL_PATH",gSystem->Getenv("ALICE_ROOT")));       
         p->UploadPackage("AliRootProofLite");
         if (p->EnablePackage("AliRootProofLite",list)) return;
       }
     }
     
     // compile task on workers
     // compile task on workers
     p->Load("AliHistogramCollection.cxx++g");
     p->Load("AliAnalysisMuMu.cxx++g");
     p->Load("AliAnalysisMuMuFromAOD.cxx++g");
     
   }
   else 
   {
     LoadLocalLibs();
   }
   
  AliAnalysisManager *mgr = new AliAnalysisManager("MuMu");
  
  AliInputEventHandler* input = new AliAODInputHandler;
//  input->SetInactiveBranches("v0s cascades jets Cells Clusters vertices");

  mgr->SetInputEventHandler(input);
  
  TString sds(dataset);

  TList triggers;
  triggers.SetOwner(kTRUE);

   triggers.Add(new TObjString("CINT7-B-NOPF-ALLNOTRD"));
   triggers.Add(new TObjString("CMUSH7-B-NOPF-MUON"));	
   triggers.Add(new TObjString("CMUL7-B-NOPF-MUON"));
   triggers.Add(new TObjString("CMUU7-B-NOPF-MUON"));

//  triggers.Add(new TObjString("CMBAC-B-NOPF-ALL"));
//   triggers.Add(new TObjString("CMBACS2-B-NOPF-ALL"));
//   triggers.Add(new TObjString("CMBACS2-B-NOPF-ALLNOTRD"));

//   triggers.Add(new TObjString("MBBG1"));
//  triggers.Add(new TObjString("CINT1B-ABCE-NOPF-ALL"));

  TString outputname("test.MuMu.AOD.1.root");
  
  if ( dataset ) 
  {
  	TString af("local");
  	
  	if ( gProof )
  	{
  	  af="unknown";
  	  TString master(gProof->GetSessionTag());
	  if (master.Contains("lx")) af = "caf";
	  if (master.Contains("nansaf")) af = "saf";
//	  if (master.Contains("skaf")) af = "skaf";
  	}
  	outputname = Form("%s.%s.root",gSystem->BaseName(dataset),af.Data());
  	cout << outputname << endl;
  }
  
  AddTaskMuMuFromAOD(outputname.Data(),&triggers);
  
//  mgr->SetDebugLevel(10); 
  
  if (!mgr->InitAnalysis()) 
  {
  	cout << "Could not InitAnalysis" << endl;
    return;
  }
  mgr->PrintStatus();
  if ( dataset ) 
  {
    mgr->StartAnalysis("proof",dataset);
  }
  else
  {
//    TChain* c = new TChain("aodTree");
//c->Add("/Users/laurent/Alice/Data/AODs/LHC10h/000137161/ESDs/pass1_plusplusplus/AOD026/0001/AliAOD.Muons.root");
//	TChain* c = CreateLocalChain("aods.local.txt");
	TChain* c = CreateLocalChain("list.txt");
//   	mgr->SetNSysInfo(1);
    TStopwatch timer;
    mgr->StartAnalysis("local",c);
    timer.Print();
//	mgr->ProfileTask("AliAnalysisMuMu-fromAOD");
  }
   
  AliCodeTimer::Instance()->Print();
}

