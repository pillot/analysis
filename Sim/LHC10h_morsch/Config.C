//
// Configuration for PWG3 barrel open charm and beauty PbPb 2010
//
// 1 HIJING event 2.76 TeV with b<30fm and central dNch/deta=2000
// + 
// N (20) PYTHIA pp 2.76 TeV Perugia0
// 20% ccbar pair per event
//     at least one in |y|<1.5
//     D mesons decay hadronically
// 20% bbbar pair per event
//     at least one in |y|<1.5
//     D mesons decay hadronically
// 20% ccbar pair per event
//     decays not forced
//     at least one electron from charm in |y|<1.2
// 20% bbbar pair per event
//     decays not forced
//     at least one electron from charm or beauty in |y|<1.2
//  10% J/psi(|y|<1.0)->e+e-
//  10% B(|y|<2.0)->J/psi(|y|<2.0)->e+e-
//
// One can use the configuration macro in compiled mode by
// root [0] gSystem->Load("libgeant321");
// root [0] gSystem->SetIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/include\
//                   -I$ALICE_ROOT -I$ALICE/geant3/TGeant3");
// root [0] .x grun.C(1,"Config.C++")
//
//
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


enum PDC06Proc_t 
{
  kPythia6, kPythia6D6T, kPythia6ATLAS, kPythia6ATLAS_Flat, kPythiaPerugia0, kPhojet, kPythiaPerugia0chadr, kPythiaPerugia0bchadr, kPythiaPerugia0cele, kPythiaPerugia0bele, kPythiaPerugia0Jpsi2e, kPythiaPerugia0BtoJpsi2e, kHijing, kHijing2000, kHijing2000HF, kHydjet, kDpmjet, kAmptHF, kAmpt, kRunMax
};

const char * pprRunName[] = {
  "kPythia6", "kPythia6D6T", "kPythia6ATLAS", "kPythia6ATLAS_Flat", "kPythiaPerugia0", "kPhojet",  "kPythiaPerugia0chadr", "kPythiaPerugia0bchadr", "kPythiaPerugia0cele", "kPythiaPerugia0bele", "kPythiaPerugia0Jpsi2e", "kPythiaPerugia0BtoJpsi2e", "kHijing", "kHijing2000", "kHijing2000HF", "kHydjet", "kDpmjet", "kAmptHF", "kAmpt"
};

enum Mag_t
{
  kNoField, k5kG, kFieldMax
};

const char * pprField[] = {
  "kNoField", "k5kG"
};

enum PprTrigConf_t
{
    kDefaultPPTrig, kDefaultPbPbTrig
};

const char * pprTrigConfName[] = {
    "p-p","Pb-Pb"
};

//--- Functions ---
class AliGenPythia;
AliGenerator *MbPythia();
AliGenerator *MbPythiaTuneD6T();
AliGenerator *MbPhojet();
AliGenerator *Hijing();
AliGenerator *Hijing2000();
AliGenerator *Hijing2000HF(Int_t typeHF);
AliGenerator *Hydjet();
AliGenerator *Dpmjet();
AliGenerator *Ampt();

void ProcessEnvironmentVars();

// Geterator, field, beam energy
static PDC06Proc_t   proc     = kHijing2000HF;
static Mag_t         mag      = k5kG;
static Float_t       energy   = 2760.; // energy in CMS
static Float_t       bMin     = 0.;
static Float_t       bMax =   = 30.;
static PprTrigConf_t strig = kDefaultPbPbTrig; // default pp trigger configuration
static Double_t      JpsiPol  = 0; // Jpsi polarisation
static Bool_t        JpsiHarderPt = kFALSE; // Jpsi harder pt spectrum (8.8 TeV)
//========================//
// Set Random Number seed //
//========================//
TDatime dt;
static UInt_t seed    = dt.Get();

// Comment line
static TString comment;

void Config()
{
  // Get settings from environment variables
  ProcessEnvironmentVars();

  gRandom->SetSeed(seed);
  cerr<<"Seed for random number generation= "<<seed<<endl; 

  // Libraries required by geant321
#if defined(__CINT__)
  gSystem->Load("liblhapdf");      // Parton density functions
  gSystem->Load("libEGPythia6");   // TGenerator interface
  if (proc == kPythia6 || proc == kPhojet || proc == kDpmjet) {
    gSystem->Load("libpythia6");        // Pythia 6.2
    gSystem->Load("libAliPythia6");     // ALICE specific implementations
  } else if (proc != kHydjet) {
    gSystem->Load("libpythia6.4.21");   // Pythia 6.4
    gSystem->Load("libAliPythia6");     // ALICE specific implementations	
  }

  if (proc == kHijing || proc == kHijing2000 || proc == kHijing2000HF) {
	  gSystem->Load("libhijing");	
  	  gSystem->Load("libTHijing");
  } else if (proc == kHydjet)  {
	  gSystem->Load("libTUHKMgen");
  } else if (proc == kDpmjet) {
	  gSystem->Load("libdpmjet");
          gSystem->Load("libTDPMjet");
  } else if (proc == kAmptHF || proc == kAmpt) {
	  gSystem->Load("libampt");
       	  gSystem->Load("libTAmpt");
  }

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
  rl->SetNumberOfEventsPerFile(1000);
  gAlice->SetRunLoader(rl);
  // gAlice->SetGeometryFromFile("geometry.root");
  // gAlice->SetGeometryFromCDB();
  
    // Set the trigger configuration
    AliSimulation::Instance()->SetTriggerConfig(pprTrigConfName[strig]);
    cout<<"Trigger configuration is set to  "<<pprTrigConfName[strig]<<endl;

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



    // RANDOM SELECTION OF ONE OF THE SIX GENERATION TYPES
    //
    Int_t typeHF = -1;
    Float_t randHF = gRandom->Rndm();
    if(randHF<0.2) {
      typeHF=0;
    } else if (randHF>=0.2 && randHF<0.4) {
      typeHF=1;
    } else if (randHF>=0.4 && randHF<0.6) {
      typeHF=2;
    } else if (randHF>=0.6 && randHF<0.8) {
      typeHF=3;
    } else if (randHF>=0.8 && randHF<0.9) {
      typeHF=4;
    } else {
      typeHF=5;
    }

    //======================//
    // Set External decayer //
    //======================//
    if (proc != kHydjet) {
      TVirtualMCDecayer* decayer = new AliDecayerPythia();
      if(proc == kHijing2000HF && (typeHF==0 || typeHF==1)) {
	decayer->SetForceDecay(kHadronicD);
      } else {
	decayer->SetForceDecay(kAll);
      }
      decayer->Init();
      gMC->SetExternalDecayer(decayer);  
    }

  //=========================//
  // Generator Configuration //
  //=========================//
  AliGenerator* gener = 0x0;
  
  if (proc == kPythia6) {
      gener = MbPythia();
  } else if (proc == kPythia6D6T) {
      gener = MbPythiaTuneD6T();
  } else if (proc == kPythia6ATLAS) {
      gener = MbPythiaTuneATLAS();
  } else if (proc == kPythiaPerugia0) {
      gener = MbPythiaTunePerugia0();
  } else if (proc == kPythia6ATLAS_Flat) {
      gener = MbPythiaTuneATLAS_Flat();
  } else if (proc == kPhojet) {
      gener = MbPhojet();
  } else if (proc == kHijing) {
      gener = Hijing();	
  } else if (proc == kHijing2000) {
      gener = Hijing2000();	
  } else if (proc == kHijing2000HF || proc == kAmptHF) {
      gener = Hijing2000HF(typeHF);	
  } else if (proc == kHydjet) {
      gener = Hydjet();	
  } else if (proc == kDpmjet) {
      gener = Dpmjet();	
  } else if (proc == kAmpt) {
      gener = Ampt();	 
  }
  
  
  //
  //
  // Size of the interaction diamond
  // Longitudinal
  Float_t sigmaz  = 5.4 / TMath::Sqrt(2.); // [cm]
  
  //
  // Transverse
  Float_t betast  = 3.5;                      // beta* [m]
  Float_t eps     = 3.75e-6;                   // emittance [m]
  Float_t gamma   = energy / 2.0 / 0.938272;  // relativistic gamma [1]
  Float_t sigmaxy = TMath::Sqrt(eps * betast / gamma) / TMath::Sqrt(2.) * 100.;  // [cm]

  printf("\n \n Diamond size x-y: %10.3e z: %10.3e\n \n", sigmaxy, sigmaz);
    
  gener->SetSigma(sigmaxy, sigmaxy, sigmaz);      // Sigma in (X,Y,Z) (cm) on IP position
  gener->SetVertexSmear(kPerEvent);
  gener->Init();

  printf("\n \n Comment: %s \n \n", comment.Data());

   //	
   // FIELD
   //

  TGeoGlobalMagField::Instance()->SetField(new AliMagF("Maps","Maps", -1., -1., AliMagF::k5kG,
     	   	AliMagF::kBeamTypeAA, 1380.));


  rl->CdGAFile();
  
  Int_t iABSO  = 1;
  Int_t iACORDE= 0;
  Int_t iDIPO  = 1;
  Int_t iEMCAL = 1;
  Int_t iFMD   = 1;
  Int_t iFRAME = 1;
  Int_t iHALL  = 1;
  Int_t iITS   = 1;
  Int_t iMAG   = 1;
  Int_t iMUON  = 1;
  Int_t iPHOS  = 1;
  Int_t iPIPE  = 1;
  Int_t iPMD   = 1;
  Int_t iHMPID = 1;
  Int_t iSHIL  = 1;
  Int_t iT0    = 1;
  Int_t iTOF   = 1;
  Int_t iTPC   = 1;
  Int_t iTRD   = 1;
  Int_t iVZERO = 1;
  Int_t iZDC   = 1;
  

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

	AliITS *ITS  = new AliITSv11Hybrid("ITS","ITS v11Hybrid");
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

        AliZDC *ZDC = new AliZDCv3("ZDC", "normal ZDC");
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
	MUON->SetTriggerEffCells(1); // not needed if raw masks 
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

        AliEMCAL *EMCAL = new AliEMCALv2("EMCAL", "EMCAL_FIRSTYEARV1");
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
//
//           PYTHIA
//

AliGenerator* MbPythia()
{
      comment = comment.Append(" pp: Pythia low-pt");
//
//    Pythia
      AliGenPythia* pythia = new AliGenPythia(-1); 
      pythia->SetMomentumRange(0, 999999.);
      pythia->SetThetaRange(0., 180.);
      pythia->SetYRange(-12.,12.);
      pythia->SetPtRange(0,1000.);
      pythia->SetProcess(kPyMb);
      pythia->SetEnergyCMS(energy);
      
      return pythia;
}

AliGenerator* MbPythiaTuneD6T()
{
      comment = comment.Append(" pp: Pythia low-pt");
//
//    Pythia
      AliGenPythia* pythia = new AliGenPythia(-1); 
      pythia->SetMomentumRange(0, 999999.);
      pythia->SetThetaRange(0., 180.);
      pythia->SetYRange(-12.,12.);
      pythia->SetPtRange(0,1000.);
      pythia->SetProcess(kPyMb);
      pythia->SetEnergyCMS(energy);
//    Tune
//    109     D6T : Rick Field's CDF Tune D6T (NB: needs CTEQ6L pdfs externally)
      pythia->SetTune(109); // F I X 
      pythia->SetStrucFunc(kCTEQ6l);
//
      return pythia;
}

AliGenerator* MbPythiaTunePerugia0()
{
      comment = comment.Append(" pp: Pythia low-pt (Perugia0)");
//
//    Pythia
      AliGenPythia* pythia = new AliGenPythia(-1); 
      pythia->SetMomentumRange(0, 999999.);
      pythia->SetThetaRange(0., 180.);
      pythia->SetYRange(-12.,12.);
      pythia->SetPtRange(0,1000.);
      pythia->SetProcess(kPyMb);
      pythia->SetEnergyCMS(energy);
//    Tune
//    320     Perugia 0
      pythia->SetTune(320); 
      pythia->UseNewMultipleInteractionsScenario();
//
      return pythia;
}


AliGenerator* MbPythiaTuneATLAS()
{
      comment = comment.Append(" pp: Pythia low-pt");
//
//    Pythia
      AliGenPythia* pythia = new AliGenPythia(-1); 
      pythia->SetMomentumRange(0, 999999.);
      pythia->SetThetaRange(0., 180.);
      pythia->SetYRange(-12.,12.);
      pythia->SetPtRange(0,1000.);
      pythia->SetProcess(kPyMb);
      pythia->SetEnergyCMS(energy);
//    Tune
//    C   306 ATLAS-CSC: Arthur Moraes' (new) ATLAS tune (needs CTEQ6L externally)
      pythia->SetTune(306);
      pythia->SetStrucFunc(kCTEQ6l);
//
      return pythia;
}

AliGenerator* PythiaJets()
{
      comment = comment.Append(" pp: Pythia low-pt");
//
//    Pythia
      AliGenPythia* pythia = new AliGenPythia(-1); 
      pythia->SetMomentumRange(0, 999999.);
      pythia->SetThetaRange(0., 180.);
      pythia->SetYRange(-12.,12.);
      pythia->SetPtRange(0,1000.);
      pythia->SetProcess(kPyJets);
      pythia->SetEnergyCMS(energy);
      pythia->SetStrucFunc(kCTEQ6l);
      pythia->SetPtHard(50, 1000);
//
      return pythia;
}

AliGenerator* MbPythiaTuneATLAS_Flat()
{
      AliGenPythia* pythia = MbPythiaTuneATLAS();
      
      comment = comment.Append("; flat multiplicity distribution");
      
      // set high multiplicity trigger
      // this weight achieves a flat multiplicity distribution
      TH1 *weight = new TH1D("weight","weight",201,-0.5,200.5);
      weight->SetBinContent(1,5.49443);
      weight->SetBinContent(2,8.770816);
      weight->SetBinContent(6,0.4568624);
      weight->SetBinContent(7,0.2919915);
      weight->SetBinContent(8,0.6674189);
      weight->SetBinContent(9,0.364737);
      weight->SetBinContent(10,0.8818444);
      weight->SetBinContent(11,0.531885);
      weight->SetBinContent(12,1.035197);
      weight->SetBinContent(13,0.9394057);
      weight->SetBinContent(14,0.9643193);
      weight->SetBinContent(15,0.94543);
      weight->SetBinContent(16,0.9426507);
      weight->SetBinContent(17,0.9423649);
      weight->SetBinContent(18,0.789456);
      weight->SetBinContent(19,1.149026);
      weight->SetBinContent(20,1.100491);
      weight->SetBinContent(21,0.6350525);
      weight->SetBinContent(22,1.351941);
      weight->SetBinContent(23,0.03233504);
      weight->SetBinContent(24,0.9574557);
      weight->SetBinContent(25,0.868133);
      weight->SetBinContent(26,1.030998);
      weight->SetBinContent(27,1.08897);
      weight->SetBinContent(28,1.251382);
      weight->SetBinContent(29,0.1391099);
      weight->SetBinContent(30,1.192876);
      weight->SetBinContent(31,0.448944);
      weight->SetBinContent(32,1);
      weight->SetBinContent(33,1);
      weight->SetBinContent(34,1);
      weight->SetBinContent(35,1);
      weight->SetBinContent(36,0.9999997);
      weight->SetBinContent(37,0.9999997);
      weight->SetBinContent(38,0.9999996);
      weight->SetBinContent(39,0.9999996);
      weight->SetBinContent(40,0.9999995);
      weight->SetBinContent(41,0.9999993);
      weight->SetBinContent(42,1);
      weight->SetBinContent(43,1);
      weight->SetBinContent(44,1);
      weight->SetBinContent(45,1);
      weight->SetBinContent(46,1);
      weight->SetBinContent(47,0.9999999);
      weight->SetBinContent(48,0.9999998);
      weight->SetBinContent(49,0.9999998);
      weight->SetBinContent(50,0.9999999);
      weight->SetBinContent(51,0.9999999);
      weight->SetBinContent(52,0.9999999);
      weight->SetBinContent(53,0.9999999);
      weight->SetBinContent(54,0.9999998);
      weight->SetBinContent(55,0.9999998);
      weight->SetBinContent(56,0.9999998);
      weight->SetBinContent(57,0.9999997);
      weight->SetBinContent(58,0.9999996);
      weight->SetBinContent(59,0.9999995);
      weight->SetBinContent(60,1);
      weight->SetBinContent(61,1);
      weight->SetBinContent(62,1);
      weight->SetBinContent(63,1);
      weight->SetBinContent(64,1);
      weight->SetBinContent(65,0.9999999);
      weight->SetBinContent(66,0.9999998);
      weight->SetBinContent(67,0.9999998);
      weight->SetBinContent(68,0.9999999);
      weight->SetBinContent(69,1);
      weight->SetBinContent(70,1);
      weight->SetBinContent(71,0.9999997);
      weight->SetBinContent(72,0.9999995);
      weight->SetBinContent(73,0.9999994);
      weight->SetBinContent(74,1);
      weight->SetBinContent(75,1);
      weight->SetBinContent(76,1);
      weight->SetBinContent(77,1);
      weight->SetBinContent(78,0.9999999);
      weight->SetBinContent(79,1);
      weight->SetBinContent(80,1);
      weight->SetEntries(526);
        
      Int_t limit = weight->GetRandom();
      pythia->SetTriggerChargedMultiplicity(limit, 1.4);
      
      comment = comment.Append(Form("; multiplicity threshold set to %d in |eta| < 1.4", limit));

      return pythia;
}

AliGenerator* MbPhojet()
{
      comment = comment.Append(" pp: Pythia low-pt");
//
//    DPMJET
#if defined(__CINT__)
  gSystem->Load("libdpmjet");      // Parton density functions
  gSystem->Load("libTDPMjet");      // Parton density functions
#endif
      AliGenDPMjet* dpmjet = new AliGenDPMjet(-1); 
      dpmjet->SetMomentumRange(0, 999999.);
      dpmjet->SetThetaRange(0., 180.);
      dpmjet->SetYRange(-12.,12.);
      dpmjet->SetPtRange(0,1000.);
      dpmjet->SetProcess(kDpmMb);
      dpmjet->SetEnergyCMS(energy);

      return dpmjet;
}


AliGenerator* MbPythiaTunePerugia0chadr()
{
      comment = comment.Append(" pp: Pythia (Perugia0) chadr (1 ccbar per event, 1 c-hadron in |y|<1.5, chadrons decay to hadrons");
//
//    Pythia
      AliGenPythia* pythia = new AliGenPythia(-1);
      pythia->SetMomentumRange(0, 999999.);
      pythia->SetThetaRange(0., 180.);
      pythia->SetYRange(-1.,1.);
      pythia->SetPtRange(0,1000.);
      pythia->SetProcess(kPyCharmppMNRwmi);
      pythia->SetEnergyCMS(energy);
//    Tune
//    320     Perugia 0
      pythia->SetTune(320);
      pythia->UseNewMultipleInteractionsScenario();
//
//    decays
      pythia->SetForceDecay(kHadronicD);

      return pythia;
}

AliGenerator* MbPythiaTunePerugia0bchadr()
{
      comment = comment.Append(" pp: Pythia (Perugia0) bchadr (1 bbbar per event, 1 c-hadron in |y|<1.5, chadrons decay to hadrons");
//
//    Pythia
      AliGenPythia* pythia = new AliGenPythia(-1);
      pythia->SetMomentumRange(0, 999999.);
      pythia->SetThetaRange(0., 180.);
      pythia->SetYRange(-1.5,1.5);
      pythia->SetPtRange(0,1000.);
      pythia->SetProcess(kPyBeautyppMNRwmi);
      pythia->SetEnergyCMS(energy);
//    Tune
//    320     Perugia 0
      pythia->SetTune(320);
      pythia->UseNewMultipleInteractionsScenario();
//
//    decays
      pythia->SetForceDecay(kHadronicD);

      return pythia;
}

AliGenerator* MbPythiaTunePerugia0cele()
{
      comment = comment.Append(" pp: Pythia (Perugia0) cele (1 ccbar per event, 1 electron in |y|<1.2");
//
//    Pythia
      AliGenPythia* pythia = new AliGenPythia(-1);
      pythia->SetMomentumRange(0, 999999.);
      pythia->SetThetaRange(0., 180.);
      //pythia->SetYRange(-2.,2.);
      pythia->SetPtRange(0,1000.);
      pythia->SetProcess(kPyCharmppMNRwmi);
      pythia->SetEnergyCMS(energy);
//    Tune
//    320     Perugia 0
      pythia->SetTune(320);
      pythia->UseNewMultipleInteractionsScenario();
//
//    decays
      pythia->SetCutOnChild(1);
      pythia->SetPdgCodeParticleforAcceptanceCut(11);
      pythia->SetChildYRange(-1.2,1.2);
      pythia->SetChildPtRange(0,10000.);

      return pythia;
}

AliGenerator* MbPythiaTunePerugia0bele()
{
      comment = comment.Append(" pp: Pythia (Perugia0) bele (1 bbbar per event, 1 electron in |y|<1.2");
//
//    Pythia
      AliGenPythia* pythia = new AliGenPythia(-1);
      pythia->SetMomentumRange(0, 999999.);
      pythia->SetThetaRange(0., 180.);
      //pythia->SetYRange(-2.,2.);
      pythia->SetPtRange(0,1000.);
      pythia->SetProcess(kPyBeautyppMNRwmi);
      pythia->SetEnergyCMS(energy);
//    Tune
//    320     Perugia 0
      pythia->SetTune(320);
      pythia->UseNewMultipleInteractionsScenario();
//
//    decays
      pythia->SetCutOnChild(1);
      pythia->SetPdgCodeParticleforAcceptanceCut(11);
      pythia->SetChildYRange(-1.2,1.2);
      pythia->SetChildPtRange(0,10000.);

      return pythia;
}

AliGenerator* MbPythiaTunePerugia0Jpsi2e()
{
  comment = comment.Append("Jpsi forced to dielectrons");
  AliGenParam *jpsi=0x0;
  if(JpsiHarderPt) jpsi = new AliGenParam(1, AliGenMUONlib::kJpsi, "CDF pp 8.8", "Jpsi");  // 8.8 TeV
  else jpsi = new AliGenParam(1, AliGenMUONlib::kJpsi, "CDF pp 7", "Jpsi");  // 7 TeV
  jpsi->SetPtRange(0.,999.);
  jpsi->SetYRange(-1.0, 1.0);
  jpsi->SetPhiRange(0.,360.);
  jpsi->SetForceDecay(kDiElectron);
  return jpsi;
}

AliGenerator* MbPythiaTunePerugia0BtoJpsi2e()
{
      comment = comment.Append(" pp: Pythia (Perugia0) BtoJpsi (1 bbbar per event, 1 b-hadron in |y|<2, 1 J/psi in |y|<2");
//
//    Pythia
      AliGenPythia* pythia = new AliGenPythia(-1);
      pythia->SetMomentumRange(0, 999999.);
      pythia->SetThetaRange(0., 180.);
      pythia->SetYRange(-2.,2.);
      pythia->SetPtRange(0,1000.);
      pythia->SetProcess(kPyBeautyppMNRwmi);
      pythia->SetEnergyCMS(energy);
//    Tune
//    320     Perugia 0
      pythia->SetTune(320);
      pythia->UseNewMultipleInteractionsScenario();
//
//    decays
      pythia->SetCutOnChild(1);
      pythia->SetPdgCodeParticleforAcceptanceCut(443);
      pythia->SetChildYRange(-2,2);
      pythia->SetChildPtRange(0,10000.);
      //
//    decays
      pythia->SetForceDecay(kBJpsiDiElectron);

      return pythia;
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

    // Field
    if (gSystem->Getenv("CONFIG_FIELD")) {
      for (Int_t iField = 0; iField < kFieldMax; iField++) {
	if (strcmp(gSystem->Getenv("CONFIG_FIELD"), pprField[iField])==0) {
	  mag = (Mag_t)iField;
	  cout<<"Field set to "<<pprField[iField]<<endl;
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

  // Impact param
    if (gSystem->Getenv("CONFIG_BMIN")) {
      bMin = atof(gSystem->Getenv("CONFIG_BMIN"));
    }

    if (gSystem->Getenv("CONFIG_BMAX")) {
      bMax = atof(gSystem->Getenv("CONFIG_BMAX"));
    }
    cout<<"Impact parameter in ["<<bMin<<","<<bMax<<"]"<<endl;
}

AliGenerator* Hijing()
{
    AliGenHijing *gener = new AliGenHijing(-1);
// centre of mass energy 
    gener->SetEnergyCMS(2760.);
    gener->SetImpactParameterRange(bMin, bMax);	
// reference frame
    gener->SetReferenceFrame("CMS");
// projectile
     gener->SetProjectile("A", 208, 82);
     gener->SetTarget    ("A", 208, 82);
// tell hijing to keep the full parent child chain
     gener->KeepFullEvent();
// enable jet quenching
     gener->SetJetQuenching(1);
// enable shadowing
     gener->SetShadowing(1);
// Don't track spectators
     gener->SetSpectators(0);
// kinematic selection
     gener->SetSelectAll(0);
     return gener;
}

AliGenerator* Hijing2000()
{
    AliGenHijing *gener = (AliGenHijing*) Hijing();
    gener->SetJetQuenching(0);	
    gener->SetPtHardMin (2.3);
    return gener;
}

AliGenerator* Hijing2000HF(Int_t typeHF)
{
  comment = comment.Append(" PbPb: Hjing2000 + pythia events for HF signals");

  AliGenCocktail *cocktail = new AliGenCocktail();
  cocktail->SetProjectile("A", 208, 82);
  cocktail->SetTarget    ("A", 208, 82);
  //
  // 1 Hijing event  
  // provides underlying event and collision geometry 
  if  (proc == kHijing2000HF) { 
  	AliGenHijing *hijing = Hijing2000();
  	cocktail->AddGenerator(hijing,"hijing",1);
  }
  if  (proc == kAmptHF) { 
  	AliGenAmpt *ampt = Ampt();
  	cocktail->AddGenerator(ampt,"ampt",1);
  }
  //
  // N Pythia Heavy Flavor events
  // N is determined from impact parameter according to the following formula 
  TFormula* formula = new TFormula("Signals", 
				   "20. * (x < 5.) + 80./3.*(1. - x/20.)*(x>5.)");
  //
  AliGenerator* pythiaHF = 0x0;      
  switch(typeHF) {
  case 0:
    pythiaHF = MbPythiaTunePerugia0chadr();
    break;
  case 1:
    pythiaHF = MbPythiaTunePerugia0bchadr();
    break;
  case 2:
    pythiaHF = MbPythiaTunePerugia0cele();
    break;
  case 3:
    pythiaHF = MbPythiaTunePerugia0bele();
    break;
  case 4:
    pythiaHF = MbPythiaTunePerugia0Jpsi2e();
    break;
  case 5:
    pythiaHF = MbPythiaTunePerugia0BtoJpsi2e();
    break;
  default:
    pythiaHF = MbPythiaTunePerugia0chadr();
    break;
  }
  cocktail->AddGenerator(pythiaHF, "pythiaHF", 1, formula); 
  //
  // Hard scattering pt > 50 GeV
  //
  AliGenerator* pythiaJets = PythiaJets();
  cocktail->AddGenerator(pythiaJets,"pythiaJets",1); 
  //
  // Pi0
  // Flat pt spectrum in range 0..20
  // Set pseudorapidity range from -1.2 to 1.2
  // 
  Float_t thmin          = (180./TMath::Pi())*2.*atan(exp(-1.2));  
  Float_t thmax          = (180./TMath::Pi())*2.*atan(exp( 1.2));  
  AliGenBox *pi0 = new AliGenBox(20);
  pi0->SetPart(111);
  pi0->SetPtRange(0.,20.);
  pi0->SetPhiRange(0,360);
  pi0->SetThetaRange(thmin,thmax);
  cocktail->AddGenerator(pi0,"pi0", 1); 
  //
  // Jpsi->mu+ mu-
  //
  jpsi2m = new AliGenParam(1, AliGenMUONlib::kJpsi, "CDF pp 3.94", "Jpsi");  // 7 TeV
  jpsi2m->SetPtRange(0.,999.);
  jpsi2m->SetYRange(-4.2, -2.3);
  jpsi2m->SetPhiRange(0.,360.);
  jpsi2m->SetForceDecay(kDiMuon);
  jpsi2m->SetCutOnChild(1);
  jpsi2m->SetChildPhiRange(0.,360.);
  jpsi2m->SetChildThetaRange(168.5,178.5);
  cocktail->AddGenerator(jpsi2m, "Jpsi2M", 1); 
  //
  // Chi_c -> J/Psi + gamma, J/Psi -> e+e-
  //
  AliGenParam* genChic = new AliGenParam(1, AliGenMUONlib::kChic,"default","Chic");
  genChic->SetMomentumRange(0, 999.);        // Wide cut on the momentum
  genChic->SetPtRange(0, 100.);              // Wide cut on Pt
  genChic->SetCutOnChild(1);                 // Enable cuts on decay products
  genChic->SetChildPhiRange(0., 360.);
  genChic->SetChildThetaRange(thmin, thmax); // In the acceptance of the Central Barrel
  genChic->SetForceDecay(kChiToJpsiGammaToElectronElectron); // Chi_c -> J/Psi + gamma, J/Psi -> e+e-
  cocktail->AddGenerator(genChic, "Chi_c", 1); 
  //	
  return cocktail;
}

AliGenerator* Hydjet()
{
  AliGenUHKM *genHi = new AliGenUHKM(-1);
  genHi->SetAllParametersLHC();
  genHi->SetProjectile("A", 208, 82);
  genHi->SetTarget    ("A", 208, 82);
  genHi->SetEcms(2760);
  genHi->SetEnergyCMS(2760.);
  genHi->SetBmin(bMin);
  genHi->SetBmax(bMax);
  genHi->SetPyquenPtmin(9);
  return genHi;
}

AliGenerator* Dpmjet()
{
  AliGenDPMjet* dpmjet = new AliGenDPMjet(-1); 
  dpmjet->SetEnergyCMS(energy);
  dpmjet->SetProjectile("A", 208, 82);
  dpmjet->SetTarget    ("A", 208, 82);
  dpmjet->SetImpactParameterRange(bMin, bMax);
  dpmjet->SetPi0Decay(0);
  return dpmjet;
}

AliGenerator* Ampt()
{

  AliGenAmpt *genHi = new AliGenAmpt(-1);
  genHi->SetEnergyCMS(2760);
  genHi->SetReferenceFrame("CMS");
  genHi->SetProjectile("A", 208, 82);
  genHi->SetTarget    ("A", 208, 82);
  genHi->SetPtHardMin (2);
  genHi->SetImpactParameterRange(bMin,bMax);
  genHi->SetJetQuenching(0); // enable jet quenching
  genHi->SetShadowing(1);    // enable shadowing
  genHi->SetDecaysOff(1);    // neutral pion and heavy particle decays switched off
  genHi->SetSpectators(0);   // track spectators 
  genHi->KeepFullEvent();
  genHi->SetSelectAll(0);
  return genHi;
}
