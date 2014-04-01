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

#include <TGeant3TGeo.h>

#include "STEER/AliRunLoader.h"

#include "STEER/AliRun.h"

#include "STEER/AliConfig.h"

#include "PYTHIA6/AliDecayerPythia.h"

#include "PYTHIA6/AliGenPythia.h"

#include "TDPMjet/AliGenDPMjet.h"

#include "STEER/AliMagF.h"

#include "STRUCT/AliBODY.h"

#include "STRUCT/AliMAG.h"

#include "STRUCT/AliABSOv3.h"

#include "STRUCT/AliDIPOv3.h"

#include "STRUCT/AliHALLv3.h"

#include "STRUCT/AliFRAMEv2.h"

#include "STRUCT/AliSHILv3.h"

#include "STRUCT/AliPIPEv3.h"

//#include "ITS/AliITSv11Hybrid.h"

#include "ITS/AliITSv11.h"

#include "TPC/AliTPCv2.h"

#include "TOF/AliTOFv6T0.h"

#include "HMPID/AliHMPIDv3.h"

#include "ZDC/AliZDCv3.h"

#include "TRD/AliTRDv1.h"

#include "TRD/AliTRDgeometry.h"

#include "FMD/AliFMDv1.h"

#include "MUON/AliMUONv1.h"

#include "PHOS/AliPHOSv1.h"


#include "PHOS/AliPHOSv1.h"

#include "PMD/AliPMDv1.h"

#include "T0/AliT0v1.h"

#include "EMCAL/AliEMCALv2.h"

#include "ACORDE/AliACORDEv1.h"

#include "VZERO/AliVZEROv7.h"

#endif



enum PDC06Proc_t
{
  kPythia6, kPhojet, kJPsiPbPb, kJPsiPbPb2760, kJPsiHptPbPb2760, kPsiPPbPb2760, kBSignalPbPb2760, kUpsiPbPb2760, kLMRPbPb2760, kRunMax
};

const char * pprRunName[] = {
  "kPythia6", "kPhojet", "kJPsiPbPb", "kJPsiPbPb2760", "kJPsiHptPbPb2760", "kPsiPPbPb2760", "kBSignalPbPb2760", "kUpsiPbPb2760", "kLMRPbPb2760"
};

//--- Magnetic Field ---

enum Mag_t
{
  kNoField, k5kG, kFieldMax
};

const char * pprField[] = {
  "kNoField", "k5kG"
};

//--- Embedding run ---

enum EmbedRun_t
{
  kBackground, kMerged, kSignal, kEmbRunMax
};

const char *embedRun[] = {
  "kBackground", "kMerged", "kSignal"
};

//--- Functions ---

class AliGenPythia;
AliGenerator *MbPythia();
AliGenerator *MbPhojet();
AliGenerator *JpsiPbPb();
AliGenerator* JPsiPbPb2760();
AliGenerator* JPsiHptPbPb2760();
AliGenerator* PsiPPbPb2760();
AliGenerator* BSignalPbPb2760();
AliGenerator* UpsiPPbPb2760();
AliGenerator *LMRPbPb2760();

void ProcessEnvironmentVars();

// Geterator, field, beam energy

static PDC06Proc_t   proc     = kPsiPPbPb2760;
static Mag_t         mag      = k5kG;
static Float_t       energy   = 2760; // energy in CMS

static EmbedRun_t    embedrun = kSignal;
//========================//

// Set Random Number seed //

//========================//

TDatime dt;
static UInt_t seed    = dt.Get();

// Comment line

static TString comment;

Float_t EtaToTheta(Float_t arg);

void Config()
{
	
	
  // Get settings from environment variables

  ProcessEnvironmentVars();
	
  gRandom->SetSeed(seed);
  cerr<<"Seed for random number generation= "<<seed<<endl;
	
  // Libraries required by geant321

#if defined(__CINT__)

  gSystem->Load("liblhapdf");      // Parton density functions

  gSystem->Load("libEG");
  gSystem->Load("libEGPythia6");   // TGenerator interface

  gSystem->Load("libpythia6");     // Pythia

  gSystem->Load("libAliPythia6");  // ALICE specific implementations

  gSystem->Load("libgeant321");
  gSystem->Load("libTTherminator");
#endif

	
  new TGeant3TGeo("C++ Interface to Geant3");
	
  //=======================================================================

  //  Create the output file

	
	
  AliRunLoader* rl=0;;
	
  cout<<"Config.C: Creating Run Loader ..."<<endl;
	//  rl = AliRunLoader::Open("galice.root");

  rl = AliRunLoader::Open("galice.root",
													AliConfig::GetDefaultEventFolderName(),
													"recreate");
	if (rl == 0x0)
	{
		gAlice->Fatal("Config.C","Can not instatiate the Run Loader");
		return;
	}
  rl->SetCompressionLevel(2);
  rl->SetNumberOfEventsPerFile(10000);
  gAlice->SetRunLoader(rl);
  
  
  // Set the trigger configuration

  if ((embedrun == kBackground) || (embedrun == kMerged)) {
		//    AliSimulation::Instance()->SetTriggerConfig("Pb-Pb");

		AliSimulation::Instance()->SetTriggerConfig("");
    cout<<"Trigger configuration is set to  Pb-Pb"<<endl;
  }
  else {
    // Set the trigger configuration: proton-proton

    //    AliSimulation::Instance()->SetTriggerConfig("p-p");

    AliSimulation::Instance()->SetTriggerConfig("Pb-Pb");
    cout<<"Trigger configuration is set to  Pb-Pb"<<endl;
  }
	
  //

  // Set External decayer

  TVirtualMCDecayer *decayer = new AliDecayerPythia();
  
  decayer->SetForceDecay(kAll);
  decayer->Init();
  gMC->SetExternalDecayer(decayer);
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
  
  if ((embedrun == kMerged) || (embedrun == kSignal) || (embedrun == kBackground)) {
    //=========================//

    // Generator Configuration //

    //=========================//

    AliGenerator* gener = 0x0;
    if (proc == kPythia6) {
      gener = MbPythia();
    } else if (proc == kPhojet) {
      gener = MbPhojet();
    } else if (proc == kJPsiPbPb) {
      gener = JPsiPbPb();
    } else if (proc == kJPsiPbPb2760) {
      gener = JPsiPbPb2760();
    } else if (proc == kJPsiHptPbPb2760) {
      gener = JPsiHptPbPb2760();
    } else if (proc == kPsiPPbPb2760) {
      gener = PsiPPbPb2760();
    } else if (proc == kBSignalPbPb2760) {
      gener = BSignalPbPb2760();
    } else if (proc == kUpsiPbPb2760) {
      gener = UpsiPbPb2760();
    } else if (proc == kLMRPbPb2760) {
      gener = LMRPbPb2760();
    }
    
  }
  else {
    AliGenCocktail *gener = new AliGenCocktail();
    gener->SetPhiRange(0, 360);
    // Set pseudorapidity range from -8 to 8.

    Float_t thmin = EtaToTheta(1);   // theta min. <---> eta max

    Float_t thmax = EtaToTheta(-1);  // theta max. <---> eta min

    gener->SetThetaRange(thmin,thmax);
    gener->SetProjectile("A",208,82);
    gener->SetTarget("A",208,82);
		
    AliGenTherminator *genther = new AliGenTherminator();
    genther->SetFileName("event.out");
    genther->SetEventNumberInFile(1);
    genther->SetTemperature(0.145);
    genther->SetMiuI(-0.0009);
    genther->SetMiuS(0.000);
    genther->SetMiuB(0.0008);
    genther->SetAlfaRange(8.0);
    genther->SetRapRange(4.0);
    genther->SetRhoMax(7.74);
    genther->SetTau(9.74);
    genther->SetModel("Lhyquid3D");
    genther->SetLhyquidSet("LHC500C2030");
		
    gener->AddGenerator(genther, "THERMINATOR LHYQUID3D", 1);
  }
  
  
	/*

	 // PRIMARY VERTEX

	 //

	 gener->SetOrigin(0., 0., 0.);    // vertex position

	 //

	 //

	 // Size of the interaction diamond

	 // Longitudinal

	 Float_t sigmaz;

	 

	 if (embedrun == kBackground) {

	 sigmaz  = 7.55 / TMath::Sqrt(2.); // [cm]

	 }

	 else {

	 Float_t sigmaz  = 5.4 / TMath::Sqrt(2.); // [cm]

	 if (energy == 900)

	 sigmaz  = 10.5 / TMath::Sqrt(2.); // [cm]

	 }

	 

	 //

	 // Transverse

	 Float_t betast  = 10;                 // beta* [m]

	 Float_t eps     = 3.75e-6;            // emittance [m]

	 Float_t gamma   = energy / 2.0 / 0.938272;  // relativistic gamma [1]

	 Float_t sigmaxy = TMath::Sqrt(eps * betast / gamma) / TMath::Sqrt(2.) * 100.;  // [cm]

	 printf("\n \n Diamond size x-y: %10.3e z: %10.3e\n \n", sigmaxy, sigmaz);

	 

	 gener->SetSigma(sigmaxy, sigmaxy, sigmaz);      // Sigma in (X,Y,Z) (cm) on IP position

	 gener->SetCutVertexZ(3.);        // Truncate at 3 sigma

	 gener->SetVertexSmear(kPerEvent);

	 */
	
  // PRIMARY VERTEX (will be overwritten by ip from real event)

  //

  gener->SetOrigin(0., 0., 0.);    // vertex position

  gener->SetSigma(0., 0., 0.); //Sigma in (X,Y,Z) (cm) on IP position

  gener->Init();
	
	/*

	 // Taken from GRP instead

	 // FIELD

	 //

	 // Field

	 

	 //  AliMagF* field = 0x0;

	 if (mag == kNoField) {

	 comment = comment.Append(" | L3 field 0.0 T");

	 TGeoGlobalMagField::Instance()->SetField(new AliMagF("Maps","Maps", 0., 0., AliMagF::k5kGUniform));

	 } else if (mag == k5kG) {

	 comment = comment.Append(" | L3 field 0.5 T");

	 TGeoGlobalMagField::Instance()->SetField(new AliMagF("Maps","Maps", -1., -1., AliMagF::k5kG));

	 }

	 printf("\n \n Comment: %s \n \n", comment.Data());

	 //  TGeoGlobalMagField::Instance()->SetField(field);

	 */
  rl->CdGAFile();
  
  Int_t iABSO  = 1;
  Int_t iACORDE= 0;
  Int_t iDIPO  = 1;
  Int_t iEMCAL = 0;
  Int_t iFMD   = 1;
  Int_t iFRAME = 1;
  Int_t iHALL  = 1;
  Int_t iITS   = 1;
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

		
		//			AliITS *ITS  = new AliITSv11Hybrid("ITS","ITS v11Hybrid");

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

		
		AliZDC *ZDC = new AliZDCv3("ZDC", "normal ZDC");
	}
	
	if (iTRD)
	{
		//=================== TRD parameters ============================

		
		AliTRD *TRD = new AliTRDv1("TRD", "TRD slow simulator");
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

		const char* digitstore="AliMUONDigitStoreV2S";
		AliMUON *MUON = new AliMUONv1("MUON", "default");
		// activate trigger efficiency by cells

		MUON->SetTriggerEffCells(1);
		// activate MTR cluster size

		MUON->SetTriggerResponseV1(2);
		if (embedrun == kBackground) {
			digitstore="AliMUONDigitStoreV2R";
			MUON->SetConvertTrigger(true);
		}
		MUON->SetDigitStoreClassName(digitstore);
		cout << "MUON DigitStore is " << MUON->DigitStoreClassName().Data() << endl;
		
		// Noise-only digits in tracker/trigger (0=no noise, 1=default (noise in tracker), 2=noise in tracker and trigger):

		cout << "****** DISABLING NOISE GENERATION AS WE DO EMBEDDING ******" << endl;
		MUON->SetDigitizerWithNoise(0);
				
	}
	
	if (iPHOS)
	{
		//=================== PHOS parameters ===========================

		AliPHOS *PHOS = new AliPHOSv1("PHOS", "IHEP");
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

		
		AliEMCAL *EMCAL = new AliEMCALv2("EMCAL", "EMCAL_COMPLETE");
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
	comment = comment.Append(" pp at 14 TeV: Pythia low-pt");
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

AliGenerator* MbPhojet()
{
	comment = comment.Append(" pp at 14 TeV: Pythia low-pt");
	//

	//    DPMJET

#if defined(__CINT__)

#endif

	gSystem->Load("libdpmjet");      // Parton density functions

	gSystem->Load("libTDPMjet");      // Parton density functions

	
	AliGenDPMjet* dpmjet = new AliGenDPMjet(-1);
	dpmjet->SetMomentumRange(0, 999999.);
	dpmjet->SetThetaRange(0., 180.);
	dpmjet->SetYRange(-12.,12.);
	dpmjet->SetPtRange(0,1000.);
	dpmjet->SetProcess(kDpmMb);
	dpmjet->SetEnergyCMS(energy);
	
	return dpmjet;
}

AliGenerator* JPsiPbPb()
{
	AliGenParam *jpsiPbPb = new AliGenParam(1, AliGenMUONlib::kJpsi,"CDF PbPb 3.94");
	jpsiPbPb->SetMomentumRange(0,999);
	jpsiPbPb->SetPtRange(0,50.);
	jpsiPbPb->SetYRange(-4.5,-2.);
	jpsiPbPb->SetPhiRange(0., 360.);
	jpsiPbPb->SetCutOnChild(1);
	jpsiPbPb->SetChildPhiRange(0.,360.);
	jpsiPbPb->SetChildThetaRange(0.,180.);
	jpsiPbPb->SetForceDecay(kDiMuon);
	jpsiPbPb->SetTrackingFlag(1);
	
	return jpsiPbPb;
}

AliGenerator* JPsiPbPb2760()
{
	AliGenParam *jpsiPbPb2760 = new AliGenParam(1, AliGenMUONlib::kJpsi,"PbPb 2.76");
	jpsiPbPb2760->SetMomentumRange(0,999);
	jpsiPbPb2760->SetPtRange(0.,999.);
	jpsiPbPb2760->SetYRange(-4.2,-2.3);
	jpsiPbPb2760->SetPhiRange(0., 360.);
	jpsiPbPb2760->SetCutOnChild(1);
	jpsiPbPb2760->SetChildPhiRange(0.,360.);
	jpsiPbPb2760->SetChildThetaRange(0.,180.);
	jpsiPbPb2760->SetForceDecay(kDiMuon);
	jpsiPbPb2760->SetTrackingFlag(1);
	
	return jpsiPbPb2760;
}

AliGenerator* JPsiHptPbPb2760()
{
	AliGenParam *jpsiHptPbPb2760 = new AliGenParam(1, AliGenMUONlib::kJpsi,"PbPb 2.76");
	jpsiHptPbPb2760->SetMomentumRange(0,999);
	jpsiHptPbPb2760->SetPtRange(4.,999.);
	jpsiHptPbPb2760->SetYRange(-4.2,-2.3);
	jpsiHptPbPb2760->SetPhiRange(0., 360.);
	jpsiHptPbPb2760->SetCutOnChild(1);
	jpsiHptPbPb2760->SetChildPhiRange(0.,360.);
	jpsiHptPbPb2760->SetChildThetaRange(0.,180.);
	jpsiHptPbPb2760->SetForceDecay(kDiMuon);
	jpsiHptPbPb2760->SetTrackingFlag(1);
	
	return jpsiHptPbPb2760;
}

AliGenerator* PsiPPbPb2760()
{
	AliGenParam *psipPbPb2760 = new AliGenParam(1, AliGenMUONlib::kPsiP,"PbPb 2.76");
	psipPbPb2760->SetMomentumRange(0,999);
	psipPbPb2760->SetPtRange(0.,999.);
	psipPbPb2760->SetYRange(-4.2,-2.3);
	psipPbPb2760->SetPhiRange(0., 360.);
	psipPbPb2760->SetCutOnChild(1);
	psipPbPb2760->SetChildPhiRange(0.,360.);
	psipPbPb2760->SetChildThetaRange(0.,180.);
	psipPbPb2760->SetForceDecay(kDiMuon);
	psipPbPb2760->SetTrackingFlag(1);
	
	return psipPbPb2760;
}

AliGenerator* BSignalPbPb2760()
{
	AliGenCorrHF *bsignalPbPb2760 = new AliGenCorrHF(1,5,3);
	bsignalPbPb2760->SetMomentumRange(0,999);
	bsignalPbPb2760->SetPtRange(0.,999.);
	bsignalPbPb2760->SetYRange(-8,0.);
	bsignalPbPb2760->SetPhiRange(0., 360.);
	bsignalPbPb2760->SetCutOnChild(1);
	bsignalPbPb2760->SetChildPhiRange(0.,360.);
	bsignalPbPb2760->SetChildThetaRange(0.,180.);
	bsignalPbPb2760->SetForceDecay(kSemiMuonic);
	bsignalPbPb2760->SetTrackingFlag(1);
	
	return bsignalPbPb2760;
}

AliGenerator* UpsiPbPb2760()
{
	AliGenParam *upsiPbPb2760 = new AliGenParam(1, AliGenMUONlib::kUpsilon,"PbPb 2.76");
	upsiPbPb2760->SetMomentumRange(0,999);
	upsiPbPb2760->SetPtRange(0.,999.);
	upsiPbPb2760->SetYRange(-4.2,-2.3);
	upsiPbPb2760->SetPhiRange(0., 360.);
	upsiPbPb2760->SetCutOnChild(1);
	upsiPbPb2760->SetChildPhiRange(0.,360.);
	upsiPbPb2760->SetChildThetaRange(0.,180.);
	upsiPbPb2760->SetForceDecay(kDiMuon);
	upsiPbPb2760->SetTrackingFlag(1);
	
	return upsiPbPb2760;
}

AliGenerator* LMRPbPb2760()
{
  AliGenMUONLMR *lmrPbPb2760 = new AliGenMUONLMR();
  lmrPbPb2760->SetCMSEnergy(AliGenMUONLMR::kCMS2760GeV);
  lmrPbPb2760->SetMomentumRange(0,999);
  lmrPbPb2760->SetPhiRange(0.,360);
  lmrPbPb2760->SetThetaRange(0.,180);
  lmrPbPb2760->SetYRange(-4.2,-2.3);
  lmrPbPb2760->SetPtRange(0.,100);
  lmrPbPb2760->SetCutOnChild(1);
  lmrPbPb2760->SetChildThetaRange(0.,180.);
  lmrPbPb2760->SetChildPhiRange(0.,360.);
  lmrPbPb2760->SetChildPtRange(0.,100);
  lmrPbPb2760->SetScaleMultiplicity(0,0); // Turn OFF Pion simulation

  lmrPbPb2760->SetScaleMultiplicity(1,0); // Turn OFF Kaon simulation

  lmrPbPb2760->SetPtParams(2,1,0.641,2.62); // setting eta params

  lmrPbPb2760->SetPtParams(3,1,1.3551,3.16); // setting rho params

  lmrPbPb2760->SetPtParams(4,1,1.3551,3.16); // setting omega params

  lmrPbPb2760->SetPtParams(5,1,1.0811,2.74); // setting phi params

  lmrPbPb2760->SetPtParams(6,1,0.72,2.5); // setting eta prime params

  
  return lmrPbPb2760;
	
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
	
	// Embedding job type

	if (gSystem->Getenv("CONFIG_EMBEDDING")) {
		for (Int_t iEmb = 0; iEmb < kEmbRunMax; iEmb++) {
			if (strcmp(gSystem->Getenv("CONFIG_EMBEDDING"), embedRun[iEmb])==0) {
				embedrun = (EmbedRun_t)iEmb;
				cout<<"Embedding run set to "<<embedRun[embedrun]<<endl;
			}
		}
		
	}
	
}

Float_t EtaToTheta(Float_t arg){
  return (180./TMath::Pi())*2.*atan(exp(-arg));
}
