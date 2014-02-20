// FOR SINGLE MUON EFFICIENCY
// Remember to define the directory and option
// gAlice->SetConfigFunction("Config('$HOME','box');");
// april 3rd: added L3 magnet 
#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TRandom.h>
#include <TDatime.h>
#include <TSystem.h>
#include <TVirtualMC.h>
#include <TGeant3TGeo.h>
#include "STEER/AliRunLoader.h"
#include "STEER/AliRun.h"
#include "STEER/AliConfig.h"
#include "PYTHIA6/AliDecayerPythia.h"
#include "PYTHIA6/AliGenPythia.h"
#include "TDPMjet/AliGenDPMjet.h"
#include "STEER/AliMagFCheb.h"
#include "STRUCT/AliBODY.h"
#include "STRUCT/AliMAG.h"
#include "STRUCT/AliABSOv3.h"
#include "STRUCT/AliDIPOv3.h"
#include "STRUCT/AliHALLv3.h"
#include "STRUCT/AliFRAMEv2.h"
#include "STRUCT/AliSHILv3.h"
#include "STRUCT/AliPIPEv3.h"
#include "ITS/AliITSv11Hybrid.h"
#include "TPC/AliTPCv2.h"
#include "TOF/AliTOFv6T0.h"
#include "HMPID/AliHMPIDv3.h"
#include "ZDC/AliZDCv3.h"
#include "TRD/AliTRDv1.h"
#include "TRD/AliTRDgeometry.h"
#include "FMD/AliFMDv1.h"
#include "MUON/AliMUONv1.h"
#include "PHOS/AliPHOSv1.h"
#include "PHOS/AliPHOSSimParam.h"
#include "PMD/AliPMDv1.h"
#include "T0/AliT0v1.h"
#include "EMCAL/AliEMCALv2.h"
#include "ACORDE/AliACORDEv1.h"
#include "VZERO/AliVZEROv7.h"
#endif

//--- Functions ---
class AliGenPythia;
void ProcessEnvironmentVars();

TDatime dt;
static UInt_t seed = dt.Get();

// Comment line
static TString comment;

void Config()
{
 
  // Get settings from environment variables
    ProcessEnvironmentVars();

    gRandom->SetSeed(seed);
    cerr<<"Seed for random number generation= "<<seed<<endl; 
 
    //===============================================================
  //  Libraries required by geant321
#if defined(__CINT__)
    gSystem->Load("liblhapdf");      // Parton density functions
    gSystem->Load("libEGPythia6");   // TGenerator interface
    gSystem->Load("libpythia6");     // Pythia
    gSystem->Load("libAliPythia6");  // ALICE specific implementations
    gSystem->Load("libgeant321");  
#endif  

    new TGeant3TGeo("C++ Interface to Geant3");

  //  Create the output file    
    AliRunLoader* rl=0x0;
    cout <<"Config.C: Creating Run Loader ..."<< endl;
    rl = AliRunLoader::Open(
	"galice.root", AliConfig::GetDefaultEventFolderName(), "recreate");
    if (rl == 0x0) {
    gAlice->Fatal("Config.C","Can not instatiate the Run Loader");
    return;
    }
    rl->SetCompressionLevel(2);
    rl->SetNumberOfEventsPerFile(20000);
    gAlice->SetRunLoader(rl);

  //=======================================================================
  // Set External decayer
    TVirtualMCDecayer *decayer = new AliDecayerPythia();
    decayer->SetForceDecay(kAll);
    decayer->Init();
    gMC->SetExternalDecayer(decayer);

  //=======================================================================
  // ******* GEANT STEERING parameters FOR ALICE SIMULATION *******
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

// BeautyMNR -PbPb 3.94 TeV --> 4 for Energy
// BeautyMNR -pp 7 TeV --> 7 for Energy

    //AliGenCorrHF *gener = new AliGenCorrHF(1, 5, 4); //beauty
    AliGenCorrHF *gener = new AliGenCorrHF(1, 4, 7); //charm
    gener->SetMomentumRange(0,9999);
    gener->SetCutOnChild(1);     
    gener->SetChildThetaRange(160.0,180.0);
    gener->SetOrigin(0,0,0);       
    gener->SetSigma(0,0,0);    
    gener->SetForceDecay(kSemiMuonic);
    gener->SetTrackingFlag(1);
    gener->Init();



/*// charmATLAS
    AliGenCorrHF *gener = new AliGenCorrHF("CharmPP7PythiaMB6.4Atlas.root",1, 4, 7);
    gener->SetMomentumRange(0,9999);
    gener->SetCutOnChild(1);     
    gener->SetChildThetaRange(160.0,180.0);
    gener->SetOrigin(0,0,0);       
    gener->SetSigma(0,0,0);    
    gener->SetForceDecay(kSemiMuonic);
    gener->SetTrackingFlag(1);
    gener->Init();

// beautyMNR
    AliGenCorrHF *gener = new AliGenCorrHF(1, 5, 7);
    gener->SetMomentumRange(0,9999);
    gener->SetCutOnChild(1);     
    gener->SetChildThetaRange(160.0,180.0);
    gener->SetOrigin(0,0,0);       
    gener->SetSigma(0,0,0);    
    gener->SetForceDecay(kSemiMuonic);
    gener->SetTrackingFlag(1);
    gener->Init();

// beautyATLAS
    AliGenCorrHF *gener = new AliGenCorrHF("BeautyPP7PythiaMB6.4Atlas.root",1, 5, 7);
    gener->SetMomentumRange(0,9999);
    gener->SetCutOnChild(1);     
    gener->SetChildThetaRange(160.0,180.0);
    gener->SetOrigin(0,0,0);       
    gener->SetSigma(0,0,0);    
    gener->SetForceDecay(kSemiMuonic);
    gener->SetTrackingFlag(1);
    gener->Init();   */

  //============================================================= 
  // Field (L3 0.5 T) outside dimuon spectrometer 
  //AliMagF *field = new AliMagF("Maps","Maps", 2, 1., 1., 10., AliMagF::k5kG);
    //   AliMagF *field = new AliMagF("Maps","Maps", -1., -1., AliMagF::k5kG,AliMagF::kBeamTypepp, 7000/2.0);
    //TGeoGlobalMagField::Instance()->SetField(field);

    rl->CdGAFile();
    
    Int_t iABSO  = 1;
    Int_t iACORDE= 0;
    Int_t iDIPO  = 1;
    Int_t iEMCAL = 0;
    Int_t iFMD   = 1;
    Int_t iFRAME = 1;
    Int_t iHALL  = 1;
    Int_t iITS   = 0;
    Int_t iMAG   = 1;
    Int_t iMUON  = 1;
    Int_t iPHOS  = 0;
    Int_t iPIPE  = 1;
    Int_t iPMD   = 0;
    Int_t iHMPID = 0;
    Int_t iSHIL  = 1;
    Int_t iT0    = 0;
    Int_t iTOF   = 0;
    Int_t iTPC   = 0;
    Int_t iTRD   = 0;
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

    if (iTPC)
    {
      //============================ TPC parameters =====================

        AliTPC *TPC = new AliTPCv2("TPC", "Default");
    }


    if (iTOF) {
        //=================== TOF parameters ============================

	AliTOF *TOF = new AliTOFv6T0("TOF", "normal TOF");
    }


    if (iHMPID)
    {
        //=================== HMPID parameters ===========================

        AliHMPID *HMPID = new AliHMPIDv3("HMPID", "normal HMPID");

    }


    if (iZDC)
    {
        //=================== ZDC parameters ============================

        AliZDC *ZDC = new AliZDCv4("ZDC", "normal ZDC");
	ZDC->SetSpectatorsTrack();	
        ZDC->SetLumiLength(0.);
    }

    if (iTRD)
    {
        //=================== TRD parameters ============================

        AliTRD *TRD = new AliTRDv1("TRD", "TRD slow simulator");
        AliTRDgeometry *geoTRD = TRD->GetGeometry();
	// Partial geometry: modules at 0,1,7,8,9,16,17
	// starting at 3h in positive direction
	geoTRD->SetSMstatus(2,0);
	geoTRD->SetSMstatus(3,0);
	geoTRD->SetSMstatus(4,0);
        geoTRD->SetSMstatus(5,0);
	geoTRD->SetSMstatus(6,0);
        geoTRD->SetSMstatus(11,0);
        geoTRD->SetSMstatus(12,0);
        geoTRD->SetSMstatus(13,0);
        geoTRD->SetSMstatus(14,0);
        geoTRD->SetSMstatus(15,0);
        geoTRD->SetSMstatus(16,0);
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
	MUON->SetTriggerEffCells(1);  // not needed if raw masks
        MUON->SetTriggerResponseV1(2);
    }

    if (iPHOS)
    {
        //=================== PHOS parameters ===========================

     AliPHOS *PHOS = new AliPHOSv1("PHOS", "noCPV_Modules123");

    }


    if (iPMD)
    {
        //=================== PMD parameters ============================

        AliPMD *PMD = new AliPMDv1("PMD", "normal PMD");
    }

    if (iT0)
    {
        //=================== T0 parameters ============================
        AliT0 *T0 = new AliT0v1("T0", "T0 Detector");
    }

    if (iEMCAL)
    {
        //=================== EMCAL parameters ============================

        AliEMCAL *EMCAL = new AliEMCALv2("EMCAL", "EMCAL_COMPLETEV1");
    }

     if (iACORDE)
    {
        //=================== ACORDE parameters ============================

        AliACORDE *ACORDE = new AliACORDEv1("ACORDE", "normal ACORDE");
    }

     if (iVZERO)
    {
        //=================== ACORDE parameters ============================

        AliVZERO *VZERO = new AliVZEROv7("VZERO", "normal VZERO");
    }
}

Float_t EtaToTheta(Float_t arg){
  return (180./TMath::Pi())*2.*atan(exp(-arg));
}

void ProcessEnvironmentVars()
{

    // Random Number seed
    if (gSystem->Getenv("CONFIG_SEED")) {
	seed = atoi(gSystem->Getenv("CONFIG_SEED"));
    }
}
