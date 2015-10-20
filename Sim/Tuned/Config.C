//
// Configuration for the first physics production 2008
//

// One can use the configuration macro in compiled mode by
// root [0] gSystem->Load("libgeant321");
// root [0] gSystem->SetIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/include\
//                   -I$ALICE_ROOT -I$ALICE/geant3/TGeant3");
// root [0] .x grun.C(1,"Config.C++")

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TRandom.h>
#include <TDatime.h>
#include <TSystem.h>
#include <TVirtualMC.h>
#include <TVirtualMCDecayer.h>
#include <TGeant3TGeo.h>
#include "AliRunLoader.h"
#include "AliRun.h"
#include "AliConfig.h"
#include "AliDecayerPythia.h"
#include "AliMagF.h"
#include "AliBODY.h"
#include "AliMAG.h"
#include "AliABSOv3.h"
#include "AliDIPOv3.h"
#include "AliHALLv3.h"
#include "AliFRAMEv2.h"
#include "AliSHILv3.h"
#include "AliPIPEv3.h"
#include "AliITSv11.h"
#include "AliZDCv4.h"
#include "AliFMDv1.h"
#include "AliMUONv1.h"
#include "AliT0v1.h"
#include "AliVZEROv7.h"
#endif 


enum PDC06Proc_t 
{
  kMuon, kRunMax
};

const char * pprRunName[] = {
  "kMuon"
};

//--- Functions ---
void ProcessEnvironmentVars();

// Generator, beam energy
static PDC06Proc_t   proc     = kMuon;
static Float_t       energy   = 2760; // energy in CMS
//========================//
// Set Random Number seed //
//========================//
TDatime dt;
static UInt_t seed    = dt.Get();

void Config()
{
  // Get settings from environment variables
  ProcessEnvironmentVars();

  gRandom->SetSeed(seed);
  cerr<<"Config.C: Seed for random number generation= "<<seed<<endl; 

  // Libraries required by geant321
#if defined(__CINT__)
  gSystem->AddIncludePath("-I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");
  gSystem->Load("libgeant321");
#endif

  new TGeant3TGeo("C++ Interface to Geant3");

  //=======================================================================
  //  Create the output file

   
  AliRunLoader* rl=0x0;

  cout<<"Config.C: Creating Run Loader ..."<<endl;
  rl = AliRunLoader::Open("galice.root",
			  AliConfig::GetDefaultEventFolderName(),
			  "recreate");
  if (rl == 0x0)
    {
      gAlice->Fatal("Config.C","Can not instatiate the Run Loader");
      return;
    }
  rl->SetCompressionLevel(2);
  rl->SetNumberOfEventsPerFile(20000);
  gAlice->SetRunLoader(rl);
  // gAlice->SetGeometryFromFile("geometry.root");
  // gAlice->SetGeometryFromCDB();
  
  //
  //=======================================================================
  // ************* STEERING parameters FOR ALICE SIMULATION **************
  // --- Specify event type to be tracked through the ALICE setup
  // --- All positions are in cm, angles in degrees, and P and E in GeV


    gMC->SetProcess("DCAY",1);
    gMC->SetProcess("PAIR",1);
    gMC->SetProcess("COMP",1);
    gMC->SetProcess("PHOT",1);
    gMC->SetProcess("PFIS",0);
    gMC->SetProcess("DRAY",0);
    gMC->SetProcess("ANNI",1);
    gMC->SetProcess("BREM",1);
    gMC->SetProcess("MUNU",1);
    gMC->SetProcess("CKOV",1);
    gMC->SetProcess("HADR",1);
    gMC->SetProcess("LOSS",2);
    gMC->SetProcess("MULS",1);
    gMC->SetProcess("RAYL",1);

    Float_t cut = 1.e-3;        // 1MeV cut by default
    Float_t tofmax = 1.e10;

    gMC->SetCut("CUTGAM", cut);
    gMC->SetCut("CUTELE", cut);
    gMC->SetCut("CUTNEU", cut);
    gMC->SetCut("CUTHAD", cut);
    gMC->SetCut("CUTMUO", cut);
    gMC->SetCut("BCUTE",  cut); 
    gMC->SetCut("BCUTM",  cut); 
    gMC->SetCut("DCUTE",  cut); 
    gMC->SetCut("DCUTM",  cut); 
    gMC->SetCut("PPCUTM", cut);
    gMC->SetCut("TOFMAX", tofmax); 

  //======================//
  // Set External decayer //
  //======================//
  TVirtualMCDecayer* decayer = new AliDecayerPythia();
  decayer->SetForceDecay(kAll);
  decayer->Init();
  gMC->SetExternalDecayer(decayer);

  //=========================//
  // Generator Configuration //
  //=========================//
  AliGenerator* gener = 0x0;
  
  if (proc == kMuon) {
    gSystem->AddIncludePath("-I$ALICE_ROOT/include");
    gSystem->AddIncludePath("-I$ALICE_ROOT/EVGEN");
    gROOT->LoadMacro("MuonGenerator.C+");
    gener = MuonGenerator();
  }  
  
  if (!gener) {
    cout<<"Generator is not set !"<<endl;
    return;
  }
  
  // Size of the interaction diamond
  /*
  // Longitudinal
  Float_t sigmaz  = 5.4 / TMath::Sqrt(2.); // [cm]
  
  // Transverse
  Float_t betast  = 3.5;                      // beta* [m]
  Float_t eps     = 3.75e-6;                   // emittance [m]
  Float_t gamma   = energy / 2.0 / 0.938272;  // relativistic gamma [1]
  Float_t sigmaxy = TMath::Sqrt(eps * betast / gamma) / TMath::Sqrt(2.) * 100.;  // [cm]
  
  printf("\n \n Diamond size x-y: %10.3e z: %10.3e\n \n", sigmaxy, sigmaz);
  
  gener->SetSigma(sigmaxy, sigmaxy, sigmaz);      // Sigma in (X,Y,Z) (cm) on IP position
  gener->SetVertexSmear(kPerEvent);
  */
  gener->SetOrigin(0., 0., 0.); //vertex position
  gener->SetSigma(0., 0., 0.); //Sigma in (X,Y,Z) (cm) on IP position
  
  gener->Init();
  
  rl->CdGAFile();
  
  
  //=========================//
  //     magnetic field      //
  //=========================//
//  TGeoGlobalMagField::Instance()->SetField(new AliMagF("Maps","Maps", -1., -1, AliMagF::k5kG));
  
  
  //=========================//
  //        material         //
  //=========================//
  Int_t iABSO  = 1;
  Int_t iDIPO  = 1;
  Int_t iFMD   = 1;
  Int_t iFRAME = 1;
  Int_t iHALL  = 1;
  Int_t iITS   = 0;
  Int_t iMAG   = 1;
  Int_t iMUON  = 1;
  Int_t iPIPE  = 1;
  Int_t iSHIL  = 1;
  Int_t iT0    = 1;
  Int_t iVZERO = 1;
  Int_t iZDC   = 0;

    //=================== Alice BODY parameters =============================
    AliBODY *BODY = new AliBODY("BODY", "Alice envelop");

    if (iMAG)
    {
        //=================== MAG parameters ============================
        // --- Start with Magnet since detector layouts may be depending ---
        // --- on the selected Magnet dimensions ---
        AliMAG *MAG = new AliMAG("MAG", "Magnet");
    }

    if (iABSO)
    {
        //=================== ABSO parameters ============================
        AliABSO *ABSO = new AliABSOv3("ABSO", "Muon Absorber");
    }

    if (iDIPO)
    {
        //=================== DIPO parameters ============================

        AliDIPO *DIPO = new AliDIPOv3("DIPO", "Dipole version 3");
    }

    if (iHALL)
    {
        //=================== HALL parameters ============================

        AliHALL *HALL = new AliHALLv3("HALL", "Alice Hall");
    }

    if (iFRAME)
    {
        //=================== FRAME parameters ============================

        AliFRAMEv2 *FRAME = new AliFRAMEv2("FRAME", "Space Frame");
	FRAME->SetHoles(1);
    }

    if (iSHIL)
    {
        //=================== SHIL parameters ============================

        AliSHIL *SHIL = new AliSHILv3("SHIL", "Shielding Version 3");
    }

    if (iPIPE)
    {
        //=================== PIPE parameters ============================

        AliPIPE *PIPE = new AliPIPEv3("PIPE", "Beam Pipe");
    }
 
    if (iITS)
    {
        //=================== ITS parameters ============================

        AliITS *ITS  = new AliITSv11("ITS","ITS v11");    
    }

    if (iZDC)
    {
        //=================== ZDC parameters ============================
	
        AliZDC *ZDC = new AliZDCv4("ZDC", "normal ZDC");
        //Collimators aperture
        ZDC->SetVCollSideCAperture(0.85);
        ZDC->SetVCollSideCCentre(0.);
        ZDC->SetVCollSideAAperture(0.75);
        ZDC->SetVCollSideACentre(0.);
        //Detector position
        ZDC->SetYZNC(1.6);
        ZDC->SetYZNA(1.6);
        ZDC->SetYZPC(1.6);
        ZDC->SetYZPA(1.6);
    }

    if (iFMD)
    {
        //=================== FMD parameters ============================

	AliFMD *FMD = new AliFMDv1("FMD", "normal FMD");
    }

    if (iMUON)
    {
        //=================== MUON parameters ===========================
        // New MUONv1 version (geometry defined via builders)
	AliMUON *MUON = new AliMUONv1("MUON", "default");
	// activate trigger efficiency by cells
	MUON->SetTriggerEffCells(1);
        MUON->SetTriggerResponseV1(2);
    }

    if (iT0)
    {
        //=================== T0 parameters ============================
        AliT0 *T0 = new AliT0v1("T0", "T0 Detector");
    }

     if (iVZERO)
    {
        //=================== ACORDE parameters ============================

        AliVZERO *VZERO = new AliVZEROv7("VZERO", "normal VZERO");
    }
}


void ProcessEnvironmentVars()
{
    // Run type
    if (gSystem->Getenv("CONFIG_RUN_TYPE")) {
      for (Int_t iRun = 0; iRun < kRunMax; iRun++) {
	if (strcmp(gSystem->Getenv("CONFIG_RUN_TYPE"), pprRunName[iRun])==0) {
	  proc = (PDC06Proc_t)iRun;
	  cout<<"Run type set to "<<pprRunName[iRun]<<endl;
	}
      }
    }

    // Energy
    if (gSystem->Getenv("CONFIG_ENERGY")) {
      energy = atoi(gSystem->Getenv("CONFIG_ENERGY"));
      cout<<"Energy set to "<<energy<<" GeV"<<endl;
    }

    // Random Number seed
    if (gSystem->Getenv("CONFIG_SEED")) {
      seed = atoi(gSystem->Getenv("CONFIG_SEED"));
    }
}
