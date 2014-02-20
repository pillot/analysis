
#include <iostream>

// ROOT includes
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TString.h"
#include "TObjArray.h"
#include "TMath.h"
#include "TBrowser.h"
#include "TGraph.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TGeoManager.h"

// STEER includes
#include "AliLog.h"
#include "AliGeomManager.h"
#include "AliCDBManager.h"

#include "AliESDEvent.h"
#include "AliESDMuonTrack.h"
#include "AliESDMuonCluster.h"
#include "AliESDInputHandler.h"

// ANALYSIS includes
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisDataSlot.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskESDCheck.h"

// MUON includes
#include "AliMUONCDB.h"
#include "AliMUONConstants.h"
#include "AliMUONESDInterface.h"
#include "AliMUONRecoParam.h"
#include "AliMUONTrack.h"
#include "AliMUONTrackParam.h"
#include "AliMUONTrackExtrap.h"
#include "AliMUONVCluster.h"
#include "AliMpDEIterator.h"

#include "AliCounterCollection.h"

ClassImp(AliAnalysisTaskESDCheck)

const Int_t AliAnalysisTaskESDCheck::fgkNTriggerClass = 10;

const char* AliAnalysisTaskESDCheck::fgkTriggerClass[10] =
{
"CBEAMB-ABCE-NOPF-ALL",
"CSMBB-ABCE-NOPF-ALL",
"CINT1A-ABCE-NOPF-ALL",
"CINT1B-ABCE-NOPF-ALL",
"CINT1C-ABCE-NOPF-ALL",
"CINT1-E-NOPF-ALL",
"CMUS1A-ABCE-NOPF-MUON",
"CMUS1B-ABCE-NOPF-MUON",
"CMUS1C-ABCE-NOPF-MUON",
"CMUS1-E-NOPF-MUON"
};

const char* AliAnalysisTaskESDCheck::fgkTriggerShortName[11] =
{
"CBEAMB",
"CSMBB",
"CINT1A",
"CINT1B",
"CINT1C",
"CINT1-E",
"CMUS1A",
"CMUS1B",
"CMUS1C",
"CMUS1-E",
"Other"
};

//________________________________________________________________________
AliAnalysisTaskESDCheck::AliAnalysisTaskESDCheck(const char *name) :
  AliAnalysisTaskSE(name), 
  fList(0x0),
  fListExpert(0x0),
  fTrackCounters(0x0),
  fEventCounters(0x0),
  fStSelect(kAll),
  fFullPrintout(kFALSE),
  fSelectCharge(0),
  fOCDBLoaded(kFALSE),
  fcurrentNEvnt(0),
  fTotalNEvnt(0),
  fcurrentNTotEvnt(0),
  fTotalNTotEvnt(0),
  fPreviousRun(-1)
//  fTest(0x0)
{
  //
  /// Constructor.
  //
  for (Int_t i=0; i<3; i++) { // loop over track types
    fcurrentNTracks[i] = 0;
    fTotalNTracks[i] = 0;
  }
  
  for (Int_t i=0; i<11; i++) { // loop over trigger classes
    fCurrentNTrig[i] = 0;
    fTotalNTrig[i] = 0;
    fCurrentNTotTrig[i] = 0;
    fTotalNTotTrig[i] = 0;
    for (Int_t j=0; j<3; j++) { // loop over track types
      fCurrentStat[i][j] = 0;
      fTotalStat[i][j] = 0;
    }
  }
  
  // Output slot #1 writes into a TObjArray container
  DefineOutput(1,TObjArray::Class());
  // Output slot #2 writes into a TObjArray container
  DefineOutput(2,TObjArray::Class());
  // Output slot #3 writes track counters
  DefineOutput(3,AliCounterCollection::Class());
  // Output slot #4 writes event counters
  DefineOutput(4,AliCounterCollection::Class());
  // Output slot #5 test
//  DefineOutput(5,TH1F::Class());
}

//________________________________________________________________________
AliAnalysisTaskESDCheck::~AliAnalysisTaskESDCheck()
{
  /// Destructor.
  delete fList;
  delete fListExpert;
  delete fTrackCounters;
  delete fEventCounters;
//  delete fTest;
}
/*
//___________________________________________________________________________
void AliAnalysisTaskESDCheck::LocalInit()
{
  /// Initialize objects locally
  Printf("   LocalInit of task %s\n", GetName());
  fTest = new TH1F("test", "number of tracks;n_{tracks}", 20, 0., 20.);
  PostData(5, fTest);
}
*/
//___________________________________________________________________________
void AliAnalysisTaskESDCheck::UserCreateOutputObjects() {
  //
  /// Create histograms
  /// Called once
  //
  Printf("   UserCreateOutputObjects of task %s\n", GetName());

  fList = new TObjArray(2000);
  fList->SetOwner();
  fListExpert = new TObjArray(2000);
  fListExpert->SetOwner();

  Int_t nCh = AliMUONConstants::NTrackingCh();
  Int_t nDE = 1100;
  
  // track info
  TH1F* hESDnTracks = new TH1F("hESDnTracks", "number of tracks;n_{tracks}", 20, 0., 20.);
  fList->AddAtAndExpand(hESDnTracks, kESDnTracks);
  
  TH1F* hESDMatchTrig = new TH1F("hESDMatchTrig", "number of tracks matched with trigger;n_{tracks}", 20, 0., 20.);
  fList->AddAtAndExpand(hESDMatchTrig, kESDMatchTrig);
  
  TH1F* hESDSign = new TH1F("hESDSign", "track sign;sign", 3, -1.5, 1.5);
  fList->AddAtAndExpand(hESDSign, kESDSign);
  
  TH1F* hESDDCA = new TH1F("hESDDCA", "DCA distribution;DCA (cm)", 500, 0., 500.);
  fList->AddAtAndExpand(hESDDCA, kESDDCA);
  
  TH1F* hESDMomentum = new TH1F("hESDMomentum", "momentum distribution;p (GeV/c)", 300, 0., 300.);
  fList->AddAtAndExpand(hESDMomentum, kESDMomentum);
  
  TH1F* hESDMomentumUncorrected = new TH1F("hESDMomentumUncorrected", "uncorrected momentum distribution;p (GeV/c)", 300, 0., 300.);
  fList->AddAtAndExpand(hESDMomentumUncorrected, kESDMomentumUncorrected);
  
  TGraphAsymmErrors* gESDMomentumUncorrected = new TGraphAsymmErrors(30);
  fList->AddAtAndExpand(gESDMomentumUncorrected, kESDMomentumUncorrectedG);
  gESDMomentumUncorrected->SetName("gESDMomentumUncorrected");
  gESDMomentumUncorrected->SetTitle("uncorrected muon momentum distribution;p (GeV/c);entries / 10 GeV/c");
  
  TGraph* gESDMomentumError = new TGraph(30);
  fList->AddAtAndExpand(gESDMomentumError, kESDMomentumErrorG);
  gESDMomentumError->SetName("gESDMomentumError");
  gESDMomentumError->SetTitle("uncorrected muon momentum resolution;p (GeV/c); #sigma_{p} (GeV/c)");
  
  TGraph* gESDMomentumRelativeError = new TGraph(30);
  fList->AddAtAndExpand(gESDMomentumRelativeError, kESDMomentumRelativeErrorG);
  gESDMomentumRelativeError->SetName("gESDMomentumRelativeError");
  gESDMomentumRelativeError->SetTitle("uncorrected muon momentum relative resolution;p (GeV/c); #sigma_{p}/p (%)");
  
  TGraph* gESDMomentumRelativeErrorMin = new TGraph(30);
  fList->AddAtAndExpand(gESDMomentumRelativeErrorMin, kESDMomentumRelativeErrorMinG);
  gESDMomentumRelativeErrorMin->SetName("gESDMomentumRelativeErrorMin");
  gESDMomentumRelativeErrorMin->SetTitle("uncorrected muon momentum relative resolution - low edge;p (GeV/c); #sigma_{p}/p (%)");
  
  TGraph* gESDMomentumRelativeErrorMax = new TGraph(30);
  fList->AddAtAndExpand(gESDMomentumRelativeErrorMax, kESDMomentumRelativeErrorMaxG);
  gESDMomentumRelativeErrorMax->SetName("gESDMomentumRelativeErrorMax");
  gESDMomentumRelativeErrorMax->SetTitle("uncorrected muon momentum relative resolution - high edge;p (GeV/c); #sigma_{p}/p (%)");
  
  TH2F* hESDUncorrMomentumRecoError = new TH2F("hESDUncorrMomentumRecoError", "uncorrected muon momentum reconstructed resolution vs p;p (GeV/c); #sigma_{p}/p (%)",
					       300, 0., 300., 1000, 0., 10.);
  fList->AddAtAndExpand(hESDUncorrMomentumRecoError, kESDUncorrMomentumRecoError);
  
  TH2F* hESDMomentumRecoError = new TH2F("hESDMomentumRecoError", "muon momentum reconstructed resolution at vertex vs p;p (GeV/c); #sigma_{p}/p (%)",
					 300, 0., 300., 1000, 0., 10.);
  fList->AddAtAndExpand(hESDMomentumRecoError, kESDMomentumRecoError);
  
  TH2F* hESDUncorrSlopeXRecoError = new TH2F("hESDUncorrSlopeXRecoError", "uncorrected muon slope_{X} reconstructed resolution vs p;p (GeV/c); #sigma_{slope_{X}}",
					       300, 0., 300., 1000, 0., 0.003);
  fList->AddAtAndExpand(hESDUncorrSlopeXRecoError, kESDUncorrSlopeXRecoError);
  
  TH2F* hESDSlopeXRecoError = new TH2F("hESDSlopeXRecoError", "muon slope_{X} reconstructed resolution at vertex vs p;p (GeV/c); #sigma_{slope_{X}}",
					     300, 0., 300., 1000, 0., 0.02);
  fList->AddAtAndExpand(hESDSlopeXRecoError, kESDSlopeXRecoError);
  
  TH2F* hESDUncorrSlopeYRecoError = new TH2F("hESDUncorrSlopeYRecoError", "uncorrected muon slope_{Y} reconstructed resolution vs p;p (GeV/c); #sigma_{slope_{Y}}",
					     300, 0., 300., 1000, 0., 0.003);
  fList->AddAtAndExpand(hESDUncorrSlopeYRecoError, kESDUncorrSlopeYRecoError);
  
  TH2F* hESDSlopeYRecoError = new TH2F("hESDSlopeYRecoError", "muon slope_{Y} reconstructed resolution at vertex vs p;p (GeV/c); #sigma_{slope_{Y}}",
				       300, 0., 300., 1000, 0., 0.02);
  fList->AddAtAndExpand(hESDSlopeYRecoError, kESDSlopeYRecoError);
  
  TH1F* hESDPt = new TH1F("hESDPt", "transverse momentum distribution;p_{t} (GeV/c)", 300, 0., 30);
  fList->AddAtAndExpand(hESDPt, kESDPt);
  
  TH2F* hESDUncorrPtRecoError = new TH2F("hESDUncorrPtRecoError", "uncorrected muon transverse momentum reconstructed resolution at vertex vs p_{t};p_{t} (GeV/c); #sigma_{p_{t}}/p_{t} (%)",
					 300, 0., 30., 300, 0., 30.);
  fList->AddAtAndExpand(hESDUncorrPtRecoError, kESDUncorrPtRecoError);
  
  TH2F* hESDPtRecoError = new TH2F("hESDPtRecoError", "muon transverse momentum reconstructed resolution at vertex vs p_{t};p_{t} (GeV/c); #sigma_{p_{t}}/p_{t} (%)",
				   300, 0., 30., 300, 0., 30.);
  fList->AddAtAndExpand(hESDPtRecoError, kESDPtRecoError);
  
  TH2F* hESDPDCARecoError = new TH2F("hESDPDCARecoError", "muon p #times DCA reconstructed resolution vs p;p (GeV/c); #sigma_{p #times DCA} (GeV #times cm)",
				   300, 0., 300., 1000, 0., 200.);
  fList->AddAtAndExpand(hESDPDCARecoError, kESDPDCARecoError);
  
  TH1F* hESDRapidity = new TH1F("hESDRapidity", "rapidity distribution;rapidity", 200, -4.5, -2.);
  fList->AddAtAndExpand(hESDRapidity, kESDRapidity);
  
  TH1F* hESDChi2 = new TH1F("hESDChi2", "normalized #chi^{2} distribution;#chi^{2} / ndf", 500, 0., 50.);
  fList->AddAtAndExpand(hESDChi2, kESDChi2);
  
  TH1F* hESDProbChi2 = new TH1F("hESDProbChi2", "distribution of probability of #chi^{2};prob(#chi^{2})", 100, 0., 1.);
  fList->AddAtAndExpand(hESDProbChi2, kESDProbChi2);
  
  TH1F* hESDThetaX = new TH1F("hESDThetaX", "#theta_{X} distribution;#theta_{X} (degree)", 360, -180., 180);
  fList->AddAtAndExpand(hESDThetaX, kESDThetaX);
  
  TH1F* hESDThetaY = new TH1F("hESDThetaY", "#theta_{Y} distribution;#theta_{Y} (degree)", 360, -180., 180);
  fList->AddAtAndExpand(hESDThetaY, kESDThetaY);
  
  // cluster info
  for (Int_t i = 0; i < nCh; i++) {
    Float_t rMax = AliMUONConstants::Rmax(i/2);
    TH2F* hESDClusterHitMap = new TH2F(Form("hESDClusterHitMap%d",i+1), Form("cluster position distribution in chamber %d;X (cm);Y (cm)",i+1),
				       100, -rMax, rMax, 100, -rMax, rMax);
    fListExpert->AddAtAndExpand(hESDClusterHitMap, kESDClusterHitMap+i);
  }
  
  TH1F* hESDnClustersPerTrack = new TH1F("hESDnClustersPerTrack", "number of associated clusters per track;n_{clusters}", 20, 0., 20.);
  fList->AddAtAndExpand(hESDnClustersPerTrack, kESDnClustersPerTrack);
  
  TH1F* hESDnChamberHitPerTrack = new TH1F("hESDnChamberHitPerTrack", "number of chambers hit per track;n_{chamber hit}", 15, 0., 15.);
  fList->AddAtAndExpand(hESDnChamberHitPerTrack, kESDnChamberHitPerTrack);
  
  TH1F* hESDnClustersPerCh = new TH1F("hESDnClustersPerCh", "averaged number of clusters per chamber per track;chamber ID;<n_{clusters}>", nCh, -0.5, nCh-0.5);
  hESDnClustersPerCh->Sumw2();
  hESDnClustersPerCh->SetOption("P");
  hESDnClustersPerCh->SetMarkerStyle(kFullDotMedium);
  hESDnClustersPerCh->SetMarkerColor(kRed);
  fList->AddAtAndExpand(hESDnClustersPerCh, kESDnClustersPerCh);
  
  TH1F* hESDnClustersPerDE = new TH1F("hESDnClustersPerDE", "averaged number of clusters per DE per track;DetElem ID;<n_{clusters}>", nDE+1, -0.5, nDE+0.5);
  hESDnClustersPerDE->Sumw2();
  hESDnClustersPerDE->SetOption("P");
  hESDnClustersPerDE->SetMarkerStyle(kFullDotMedium);
  hESDnClustersPerDE->SetMarkerColor(kRed);
  fList->AddAtAndExpand(hESDnClustersPerDE, kESDnClustersPerDE);
  
  for (Int_t i = 0; i < nCh; i++) {
    TH1F* hESDClusterChargeInCh = new TH1F(Form("hESDClusterChargeInCh%d",i+1), Form("cluster charge distribution in chamber %d;charge (ADC counts)",i+1), 500, 0., 5000.);
    fListExpert->AddAtAndExpand(hESDClusterChargeInCh, kESDClusterChargeInCh+i);
  }
  
  TH1F* hESDClusterChargePerChMean = new TH1F("hESDClusterChargePerChMean", "cluster mean charge per chamber;chamber ID;<charge> (ADC counts)", nCh, -0.5, nCh-0.5);
  hESDClusterChargePerChMean->SetOption("P");
  hESDClusterChargePerChMean->SetMarkerStyle(kFullDotMedium);
  hESDClusterChargePerChMean->SetMarkerColor(kRed);
  fList->AddAtAndExpand(hESDClusterChargePerChMean, kESDClusterChargePerChMean);
  
  TH1F* hESDClusterChargePerChSigma = new TH1F("hESDClusterChargePerChSigma", "cluster charge dispersion per chamber;chamber ID;#sigma_{charge} (ADC counts)", nCh, -0.5, nCh-0.5);
  hESDClusterChargePerChSigma->SetOption("P");
  hESDClusterChargePerChSigma->SetMarkerStyle(kFullDotMedium);
  hESDClusterChargePerChSigma->SetMarkerColor(kRed);
  fList->AddAtAndExpand(hESDClusterChargePerChSigma, kESDClusterChargePerChSigma);
  
  TH1F* hESDClusterChargePerDE = new TH1F("hESDClusterChargePerDE", "cluster mean charge per DE;DetElem ID;<charge> (ADC counts)", nDE+1, -0.5, nDE+0.5);
  hESDClusterChargePerDE->SetOption("P");
  hESDClusterChargePerDE->SetMarkerStyle(kFullDotMedium);
  hESDClusterChargePerDE->SetMarkerColor(kRed);
  fList->AddAtAndExpand(hESDClusterChargePerDE, kESDClusterChargePerDE);
  
  for (Int_t i = 0; i < nCh; i++) {
    TH1F* hESDClusterSizeInCh = new TH1F(Form("hESDClusterSizeInCh%d",i+1), Form("cluster size distribution in chamber %d;size (n_{pads})",i+1), 200, 0., 200.);
    fListExpert->AddAtAndExpand(hESDClusterSizeInCh, kESDClusterSizeInCh+i);
  }
  
  TH1F* hESDClusterSizePerChMean = new TH1F("hESDClusterSizePerChMean", "cluster mean size per chamber;chamber ID;<size> (n_{pads})", nCh, -0.5, nCh-0.5);
  hESDClusterSizePerChMean->SetOption("P");
  hESDClusterSizePerChMean->SetMarkerStyle(kFullDotMedium);
  hESDClusterSizePerChMean->SetMarkerColor(kRed);
  fList->AddAtAndExpand(hESDClusterSizePerChMean, kESDClusterSizePerChMean);
  
  TH1F* hESDClusterSizePerChSigma = new TH1F("hESDClusterSizePerChSigma", "cluster size dispersion per chamber;chamber ID;#sigma_{size} (n_{pads})", nCh, -0.5, nCh-0.5);
  hESDClusterSizePerChSigma->SetOption("P");
  hESDClusterSizePerChSigma->SetMarkerStyle(kFullDotMedium);
  hESDClusterSizePerChSigma->SetMarkerColor(kRed);
  fList->AddAtAndExpand(hESDClusterSizePerChSigma, kESDClusterSizePerChSigma);
  
  TH1F* hESDClusterSizePerDE = new TH1F("hESDClusterSizePerDE", "cluster mean size per DE;DetElem ID;<size> (n_{pads})", nDE+1, -0.5, nDE+0.5);
  hESDClusterSizePerDE->SetOption("P");
  hESDClusterSizePerDE->SetMarkerStyle(kFullDotMedium);
  hESDClusterSizePerDE->SetMarkerColor(kRed);
  fList->AddAtAndExpand(hESDClusterSizePerDE, kESDClusterSizePerDE);
  
  // cluster - track info
  for (Int_t i = 0; i < nCh; i++) {
    TH1F* hESDResidualXInCh = new TH1F(Form("hESDResidualXInCh%d",i+1), Form("cluster-track residual-X distribution in chamber %d;#Delta_{X} (cm)",i+1), 1000, -5., 5.);
    fListExpert->AddAtAndExpand(hESDResidualXInCh, kESDResidualXInCh+i);
    
    TH1F* hESDResidualYInCh = new TH1F(Form("hESDResidualYInCh%d",i+1), Form("cluster-track residual-Y distribution in chamber %d;#Delta_{Y} (cm)",i+1), 1000, -5., 5.);
    fListExpert->AddAtAndExpand(hESDResidualYInCh, kESDResidualYInCh+i);
    
    TH1F* hESDLocalChi2XInCh = new TH1F(Form("hESDLocalChi2XInCh%d",i+1), Form("local chi2-X distribution in chamber %d;local #chi^{2}_{X}",i+1), 1000, 0., 25);
    fListExpert->AddAtAndExpand(hESDLocalChi2XInCh, kESDLocalChi2XInCh+i);
    
    TH1F* hESDLocalChi2YInCh = new TH1F(Form("hESDLocalChi2YInCh%d",i+1), Form("local chi2-Y distribution in chamber %d;local #chi^{2}_{Y}",i+1), 1000, 0., 25);
    fListExpert->AddAtAndExpand(hESDLocalChi2YInCh, kESDLocalChi2YInCh+i);
    
    TH1F* hESDLocalChi2InCh = new TH1F(Form("hESDLocalChi2InCh%d",i+1), Form("local chi2 (~0.5*(#chi^{2}_{X}+#chi^{2}_{Y})) distribution in chamber %d;local #chi^{2}",i+1), 1000, 0., 25);
    fListExpert->AddAtAndExpand(hESDLocalChi2InCh, kESDLocalChi2InCh+i);
  }
  
  TH1F* hESDResidualXPerChMean = new TH1F("hESDResidualXPerChMean", "cluster-track residual-X per Ch: mean;chamber ID;<#Delta_{X}> (cm)", nCh, -0.5, nCh-0.5);
  hESDResidualXPerChMean->SetOption("P");
  hESDResidualXPerChMean->SetMarkerStyle(kFullDotMedium);
  hESDResidualXPerChMean->SetMarkerColor(kRed);
  fList->AddAtAndExpand(hESDResidualXPerChMean, kESDResidualXPerChMean);
  
  TH1F* hESDResidualYPerChMean = new TH1F("hESDResidualYPerChMean", "cluster-track residual-Y per Ch: mean;chamber ID;<#Delta_{Y}> (cm)", nCh, -0.5, nCh-0.5);
  hESDResidualYPerChMean->SetOption("P");
  hESDResidualYPerChMean->SetMarkerStyle(kFullDotMedium);
  hESDResidualYPerChMean->SetMarkerColor(kRed);
  fList->AddAtAndExpand(hESDResidualYPerChMean, kESDResidualYPerChMean);
  
  TH1F* hESDResidualXPerChSigma = new TH1F("hESDResidualXPerChSigma", "cluster-track residual-X per Ch: sigma;chamber ID;#sigma_{X} (cm)", nCh, -0.5, nCh-0.5);
  hESDResidualXPerChSigma->SetOption("P");
  hESDResidualXPerChSigma->SetMarkerStyle(kFullDotMedium);
  hESDResidualXPerChSigma->SetMarkerColor(kRed);
  fList->AddAtAndExpand(hESDResidualXPerChSigma, kESDResidualXPerChSigma);
  
  TH1F* hESDResidualYPerChSigma = new TH1F("hESDResidualYPerChSigma", "cluster-track residual-Y per Ch: sigma;chamber ID;#sigma_{Y} (cm)", nCh, -0.5, nCh-0.5);
  hESDResidualYPerChSigma->SetOption("P");
  hESDResidualYPerChSigma->SetMarkerStyle(kFullDotMedium);
  hESDResidualYPerChSigma->SetMarkerColor(kRed);
  fList->AddAtAndExpand(hESDResidualYPerChSigma, kESDResidualYPerChSigma);
  
  TH1F* hESDLocalChi2XPerChMean = new TH1F("hESDLocalChi2XPerChMean", "local chi2-X per Ch: mean;chamber ID;<local #chi^{2}_{X}>", nCh, -0.5, nCh-0.5);
  hESDLocalChi2XPerChMean->SetOption("P");
  hESDLocalChi2XPerChMean->SetMarkerStyle(kFullDotMedium);
  hESDLocalChi2XPerChMean->SetMarkerColor(kRed);
  fList->AddAtAndExpand(hESDLocalChi2XPerChMean, kESDLocalChi2XPerChMean);
  
  TH1F* hESDLocalChi2YPerChMean = new TH1F("hESDLocalChi2YPerChMean", "local chi2-Y per Ch: mean;chamber ID;<local #chi^{2}_{Y}>", nCh, -0.5, nCh-0.5);
  hESDLocalChi2YPerChMean->SetOption("P");
  hESDLocalChi2YPerChMean->SetMarkerStyle(kFullDotMedium);
  hESDLocalChi2YPerChMean->SetMarkerColor(kRed);
  fList->AddAtAndExpand(hESDLocalChi2YPerChMean, kESDLocalChi2YPerChMean);
  
  TH1F* hESDLocalChi2PerChMean = new TH1F("hESDLocalChi2PerChMean", "local chi2 (~0.5*(#chi^{2}_{X}+#chi^{2}_{Y})) per Ch: mean;chamber ID;<local #chi^{2}>", nCh, -0.5, nCh-0.5);
  hESDLocalChi2PerChMean->SetOption("P");
  hESDLocalChi2PerChMean->SetMarkerStyle(kFullDotMedium);
  hESDLocalChi2PerChMean->SetMarkerColor(kRed);
  fList->AddAtAndExpand(hESDLocalChi2PerChMean, kESDLocalChi2PerChMean);
  
  TH1F* hESDResidualXPerDEMean = new TH1F("hESDResidualXPerDEMean", "cluster-track residual-X per DE: mean;DetElem ID;<#Delta_{X}> (cm)", nDE+1, -0.5, nDE+0.5);
  hESDResidualXPerDEMean->SetOption("P");
  hESDResidualXPerDEMean->SetMarkerStyle(kFullDotMedium);
  hESDResidualXPerDEMean->SetMarkerColor(kRed);
  fList->AddAtAndExpand(hESDResidualXPerDEMean, kESDResidualXPerDEMean);
  
  TH1F* hESDResidualYPerDEMean = new TH1F("hESDResidualYPerDEMean", "cluster-track residual-Y per DE: mean;DetElem ID;<#Delta_{Y}> (cm)", nDE+1, -0.5, nDE+0.5);
  hESDResidualYPerDEMean->SetOption("P");
  hESDResidualYPerDEMean->SetMarkerStyle(kFullDotMedium);
  hESDResidualYPerDEMean->SetMarkerColor(kRed);
  fList->AddAtAndExpand(hESDResidualYPerDEMean, kESDResidualYPerDEMean);
  
  TH1F* hESDResidualXPerDESigma = new TH1F("hESDResidualXPerDESigma", "cluster-track residual-X per DE: sigma;DetElem ID;#sigma_{X} (cm)", nDE+1, -0.5, nDE+0.5);
  hESDResidualXPerDESigma->SetOption("P");
  hESDResidualXPerDESigma->SetMarkerStyle(kFullDotMedium);
  hESDResidualXPerDESigma->SetMarkerColor(kRed);
  fList->AddAtAndExpand(hESDResidualXPerDESigma, kESDResidualXPerDESigma);
  
  TH1F* hESDResidualYPerDESigma = new TH1F("hESDResidualYPerDESigma", "cluster-track residual-Y per DE: sigma;DetElem ID;#sigma_{Y} (cm)", nDE+1, -0.5, nDE+0.5);
  hESDResidualYPerDESigma->SetOption("P");
  hESDResidualYPerDESigma->SetMarkerStyle(kFullDotMedium);
  hESDResidualYPerDESigma->SetMarkerColor(kRed);
  fList->AddAtAndExpand(hESDResidualYPerDESigma, kESDResidualYPerDESigma);
  
  TH1F* hESDLocalChi2XPerDEMean = new TH1F("hESDLocalChi2XPerDEMean", "local chi2-X per DE: mean;DetElem ID;<local #chi^{2}_{X}>", nDE+1, -0.5, nDE+0.5);
  hESDLocalChi2XPerDEMean->SetOption("P");
  hESDLocalChi2XPerDEMean->SetMarkerStyle(kFullDotMedium);
  hESDLocalChi2XPerDEMean->SetMarkerColor(kRed);
  fList->AddAtAndExpand(hESDLocalChi2XPerDEMean, kESDLocalChi2XPerDEMean);
  
  TH1F* hESDLocalChi2YPerDEMean = new TH1F("hESDLocalChi2YPerDEMean", "local chi2-Y per DE: mean;DetElem ID;<local #chi^{2}_{Y}>", nDE+1, -0.5, nDE+0.5);
  hESDLocalChi2YPerDEMean->SetOption("P");
  hESDLocalChi2YPerDEMean->SetMarkerStyle(kFullDotMedium);
  hESDLocalChi2YPerDEMean->SetMarkerColor(kRed);
  fList->AddAtAndExpand(hESDLocalChi2YPerDEMean, kESDLocalChi2YPerDEMean);

  TH1F* hESDLocalChi2PerDEMean = new TH1F("hESDLocalChi2PerDEMean", "local chi2 (~0.5*(#chi^{2}_{X}+#chi^{2}_{Y})) per DE: mean;DetElem ID;<local #chi^{2}>", nDE+1, -0.5, nDE+0.5);
  hESDLocalChi2PerDEMean->SetOption("P");
  hESDLocalChi2PerDEMean->SetMarkerStyle(kFullDotMedium);
  hESDLocalChi2PerDEMean->SetMarkerColor(kRed);
  fList->AddAtAndExpand(hESDLocalChi2PerDEMean, kESDLocalChi2PerDEMean);
  
  // intermediate histograms
  TH1F* hESDnTotClustersPerCh = new TH1F("hESDnTotClustersPerCh", "total number of associated clusters per chamber;chamber ID;#Sigma(n_{clusters})", nCh, -0.5, nCh-0.5);
  fListExpert->AddAtAndExpand(hESDnTotClustersPerCh, kESDnTotClustersPerCh);
  TH1F* hESDnTotClustersPerDE = new TH1F("hESDnTotClustersPerDE", "total number of associated clusters per DE;DetElem ID;#Sigma(n_{clusters})", nDE+1, -0.5, nDE+0.5);
  fListExpert->AddAtAndExpand(hESDnTotClustersPerDE, kESDnTotClustersPerDE);
  TH1F* hESDnTotFullClustersPerDE = new TH1F("hESDnTotFullClustersPerDE", "total number of associated clusters containing pad info per DE;DetElem ID;#Sigma(n_{full clusters})", nDE+1, -0.5, nDE+0.5);
  fListExpert->AddAtAndExpand(hESDnTotFullClustersPerDE, kESDnTotFullClustersPerDE);
  TH1F* hESDSumClusterChargePerDE = new TH1F("hESDSumClusterChargePerDE", "sum of cluster charge per DE;DetElem ID;#Sigma(charge) (ADC counts)", nDE+1, -0.5, nDE+0.5);
  fListExpert->AddAtAndExpand(hESDSumClusterChargePerDE, kESDSumClusterChargePerDE);
  TH1F* hESDSumClusterSizePerDE = new TH1F("hESDSumClusterSizePerDE", "sum of cluster size per DE;DetElem ID;#Sigma(size) (n_{pads})", nDE+1, -0.5, nDE+0.5);
  fListExpert->AddAtAndExpand(hESDSumClusterSizePerDE, kESDSumClusterSizePerDE);
  TH1F* hESDSumResidualXPerDE = new TH1F("hESDSumResidualXPerDE", "sum of cluster-track residual-X per DE;DetElem ID;#Sigma(#Delta_{X}) (cm)", nDE+1, -0.5, nDE+0.5);
  fListExpert->AddAtAndExpand(hESDSumResidualXPerDE, kESDSumResidualXPerDE);
  TH1F* hESDSumResidualYPerDE = new TH1F("hESDSumResidualYPerDE", "sum of cluster-track residual-Y per DE;DetElem ID;#Sigma(#Delta_{Y}) (cm)", nDE+1, -0.5, nDE+0.5);
  fListExpert->AddAtAndExpand(hESDSumResidualYPerDE, kESDSumResidualYPerDE);
  TH1F* hESDSumResidualX2PerDE = new TH1F("hESDSumResidualX2PerDE", "sum of cluster-track residual-X**2 per DE;DetElem ID;#Sigma(#Delta_{X}^{2}) (cm^{2})", nDE+1, -0.5, nDE+0.5);
  fListExpert->AddAtAndExpand(hESDSumResidualX2PerDE, kESDSumResidualX2PerDE);
  TH1F* hESDSumResidualY2PerDE = new TH1F("hESDSumResidualY2PerDE", "sum of cluster-track residual-Y**2 per DE;DetElem ID;#Sigma(#Delta_{Y}^{2}) (cm^{2})", nDE+1, -0.5, nDE+0.5);
  fListExpert->AddAtAndExpand(hESDSumResidualY2PerDE, kESDSumResidualY2PerDE);
  TH1F* hESDSumLocalChi2XPerDE = new TH1F("hESDSumLocalChi2XPerDE", "sum of local chi2-X per DE;DetElem ID;#Sigma(local #chi^{2}_{X})", nDE+1, -0.5, nDE+0.5);
  fListExpert->AddAtAndExpand(hESDSumLocalChi2XPerDE, kESDSumLocalChi2XPerDE);
  TH1F* hESDSumLocalChi2YPerDE = new TH1F("hESDSumLocalChi2YPerDE", "sum of local chi2-Y per DE;DetElem ID;#Sigma(local #chi^{2}_{Y})", nDE+1, -0.5, nDE+0.5);
  fListExpert->AddAtAndExpand(hESDSumLocalChi2YPerDE, kESDSumLocalChi2YPerDE);
  TH1F* hESDSumLocalChi2PerDE = new TH1F("hESDSumLocalChi2PerDE", "sum of local chi2 (~0.5*(#chi^{2}_{X}+#chi^{2}_{Y})) per DE;DetElem ID;#Sigma(local #chi^{2})", nDE+1, -0.5, nDE+0.5);
  fListExpert->AddAtAndExpand(hESDSumLocalChi2PerDE, kESDSumLocalChi2PerDE);
  
  // initialize track counters
  fTrackCounters = new AliCounterCollection("trackCounters");
  fTrackCounters->AddRubric("track", "tracker/trigger/matched/any");
  TString triggerClassNames = "/";
  for (Int_t i=0; i<=AliAnalysisTaskESDCheck::fgkNTriggerClass; i++)
    triggerClassNames += Form("%s/",AliAnalysisTaskESDCheck::fgkTriggerShortName[i]);
  triggerClassNames += "any/";
  fTrackCounters->AddRubric("trigger", triggerClassNames.Data());
  fTrackCounters->AddRubric("run", 1000000);
  fTrackCounters->AddRubric("selected", "yes/no");
  fTrackCounters->Init();
  
  // initialize event counters
  fEventCounters = new AliCounterCollection("eventCounters");
  fEventCounters->AddRubric("event", "muon/any");
  fEventCounters->AddRubric("trigger", triggerClassNames.Data());
  fEventCounters->AddRubric("run", 1000000);
  fEventCounters->AddRubric("selected", "yes/no");
  fEventCounters->Init();
  
  // Post data at least once per task to ensure data synchronisation (required for merging)
  PostData(1, fList);
  PostData(2, fListExpert);
  PostData(3, fTrackCounters);
  PostData(4, fEventCounters);
}

//________________________________________________________________________
void AliAnalysisTaskESDCheck::UserExec(Option_t *) {
  //
  /// Main loop
  /// Called for each event
  //
  
  // check if OCDB properly loaded
  if (!fOCDBLoaded) return;
  
  // check physics selection
  TString selected = (fInputHandler && fInputHandler->IsEventSelected()) ? "selected:yes" : "selected:no";
  
  AliESDEvent* fESD = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!fESD) {
    Printf("ERROR: fESD not available");
    return;
  }
  
  Int_t nTracks = (Int_t) fESD->GetNumberOfMuonTracks(); 
  Int_t nTrackerTracks = 0;
  Int_t nTriggerTracks = 0;
  Int_t nTrackMatchTrig = 0;
  AliMUONTrack track;
  
  // fill event statistics
  fcurrentNTotEvnt++;
  fTotalNTotEvnt++;
  fEventCounters->Count(Form("event:any/trigger:any/run:%d/%s", fCurrentRunNumber, selected.Data()));
  
  Bool_t triggerFired = kFALSE;
  for (Int_t i=0; i<10; i++) {
    if (fESD->IsTriggerClassFired(AliAnalysisTaskESDCheck::fgkTriggerClass[i])) {
      fCurrentNTotTrig[i]++;
      fTotalNTotTrig[i]++;
      fEventCounters->Count(Form("event:any/trigger:%s/run:%d/%s", AliAnalysisTaskESDCheck::fgkTriggerShortName[i], fCurrentRunNumber, selected.Data()));
      triggerFired = kTRUE;
    }
  }
  if (!triggerFired) {
    fCurrentNTotTrig[10]++;
    fTotalNTotTrig[10]++;
    fEventCounters->Count(Form("event:any/trigger:other/run:%d/%s", fCurrentRunNumber, selected.Data()));
  }
  
  // loop over tracks
  for (Int_t iTrack = 0; iTrack < nTracks; ++iTrack) {
    
    // get the ESD track and skip "ghosts"
    AliESDMuonTrack* esdTrack = fESD->GetMuonTrack(iTrack);
    if (!esdTrack->ContainTrackerData()) {
      nTriggerTracks++;
      continue;
    }
    
    nTrackerTracks++;
    
    if (esdTrack->ContainTriggerData()) {
      nTriggerTracks++;
      nTrackMatchTrig++;
    }
    
    // select on track charge
    if (fSelectCharge*esdTrack->Charge() < 0) continue;
    
    //if (!esdTrack->ContainTriggerData()) cout<<"tracker/trigger matching failed for event: "<<
    // AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()->GetTree()->GetTree()->GetReadEntry()<<endl;
    
    ((TH1F*)fList->UncheckedAt(kESDMomentum))->Fill(esdTrack->P());
    ((TH1F*)fList->UncheckedAt(kESDMomentumUncorrected))->Fill(esdTrack->PUncorrected());
    ((TH1F*)fList->UncheckedAt(kESDPt))->Fill(esdTrack->Pt());
    ((TH1F*)fList->UncheckedAt(kESDRapidity))->Fill(esdTrack->Y());
    Int_t ndf = 2 * esdTrack->GetNHit() - 5;
    ((TH1F*)fList->UncheckedAt(kESDChi2))->Fill(esdTrack->GetChi2()/ndf);
    ((TH1F*)fList->UncheckedAt(kESDProbChi2))->Fill(TMath::Prob(esdTrack->GetChi2(),ndf));
    ((TH1F*)fList->UncheckedAt(kESDThetaX))->Fill(ChangeThetaRange(esdTrack->GetThetaXUncorrected()));
    ((TH1F*)fList->UncheckedAt(kESDThetaY))->Fill(ChangeThetaRange(esdTrack->GetThetaYUncorrected()));
    ((TH1F*)fList->UncheckedAt(kESDnClustersPerTrack))->Fill(esdTrack->GetNHit());
    ((TH1F*)fList->UncheckedAt(kESDSign))->Fill(esdTrack->Charge());
    ((TH1F*)fList->UncheckedAt(kESDDCA))->Fill(esdTrack->GetDCA());
    
    Int_t nChamberHit = 0;
    for (Int_t ich=0; ich<10; ich++) if (esdTrack->IsInMuonClusterMap(ich)) nChamberHit++;
    ((TH1F*)fList->UncheckedAt(kESDnChamberHitPerTrack))->Fill(nChamberHit);
    
    // get corresponding MUON track
    AliMUONESDInterface::ESDToMUON(*esdTrack,track);
    //if (!SelectTrack(track)) continue;
    
    // momentum error at first cluster
    //AliMUONTrackParam* trackParam0 = static_cast<AliMUONTrackParam*>(track.GetTrackParamAtCluster()->First());
    //AliMUONESDInterface::SetParamAtFirstCluster(*trackParam0,*esdTrack);
    //AliMUONESDInterface::SetParamCov(*trackParam0,*esdTrack);
    Double_t cov[21];
    esdTrack->GetCovarianceXYZPxPyPz(cov);
    Double_t pUX = esdTrack->PxUncorrected();
    Double_t pUY = esdTrack->PyUncorrected();
    Double_t pUZ = esdTrack->PzUncorrected();
    Double_t pU  = esdTrack->PUncorrected();
    Double_t sigmaPU = TMath::Sqrt(pUX * (pUX*cov[9]  + pUY*cov[13] + pUZ*cov[18]) +
				   pUY * (pUX*cov[13] + pUY*cov[14] + pUZ*cov[19]) +
				   pUZ * (pUX*cov[18] + pUY*cov[19] + pUZ*cov[20])) / pU;
    ((TH2F*)fList->UncheckedAt(kESDUncorrMomentumRecoError))->Fill(pU,100.*sigmaPU/pU);
    
    // transverse momentum error at first cluster
    Double_t pTU  = TMath::Sqrt(pUX*pUX + pUY*pUY);
    Double_t sigmaPTU = TMath::Sqrt(pUX * (pUX*cov[9]  + pUY*cov[13]) + pUY * (pUX*cov[13] + pUY*cov[14])) / pTU;
    ((TH2F*)fList->UncheckedAt(kESDUncorrPtRecoError))->Fill(pTU,100.*sigmaPTU/pTU);
    
    // slopeX-Y error at first cluster
    AliMUONTrackParam trackParamAtFirstCluster;
    AliMUONESDInterface::GetParamAtFirstCluster(*esdTrack,trackParamAtFirstCluster);
    AliMUONESDInterface::GetParamCov(*esdTrack,trackParamAtFirstCluster);
    Int_t firstCh = 0;
    while (firstCh < 10 && !esdTrack->IsInMuonClusterMap(firstCh)) firstCh++;
    AliMUONTrackExtrap::AddMCSEffect(&trackParamAtFirstCluster, AliMUONConstants::ChamberThicknessInX0(firstCh)/2., -1.);
    const TMatrixD& paramCovU = trackParamAtFirstCluster.GetCovariances();
    ((TH2F*)fList->UncheckedAt(kESDUncorrSlopeXRecoError))->Fill(pU,TMath::Sqrt(paramCovU(1,1)));
    ((TH2F*)fList->UncheckedAt(kESDUncorrSlopeYRecoError))->Fill(pU,TMath::Sqrt(paramCovU(3,3)));
    
    // momentum error at vertex
    AliMUONTrackParam trackParamAtVtx(trackParamAtFirstCluster);
    AliMUONTrackExtrap::ExtrapToVertex(&trackParamAtVtx, esdTrack->GetNonBendingCoor(), esdTrack->GetBendingCoor(), esdTrack->GetZ(), 0., 0.);
    TMatrixD paramCovP(trackParamAtVtx.GetCovariances());
    Cov2CovP(trackParamAtVtx.GetParameters(),paramCovP);
    Double_t p = esdTrack->P();
    Double_t sigmaP = TMath::Sqrt(paramCovP(4,4));
    ((TH2F*)fList->UncheckedAt(kESDMomentumRecoError))->Fill(p,100.*sigmaP/p);
    
    // transverse momentum error at vertex
    AliESDMuonTrack tmpTrack;
    AliMUONESDInterface::SetParamAtFirstCluster(trackParamAtVtx, tmpTrack);
    AliMUONESDInterface::SetParamCov(trackParamAtVtx, tmpTrack);
    Double_t covAtVtx[21];
    tmpTrack.GetCovarianceXYZPxPyPz(covAtVtx);
    Double_t px = esdTrack->Px();
    Double_t py = esdTrack->Py();
    Double_t pt = esdTrack->Pt();
    Double_t sigmaPt = TMath::Sqrt(px * (px*covAtVtx[9]  + py*covAtVtx[13]) + py * (px*covAtVtx[13] + py*covAtVtx[14])) / pt;
    ((TH2F*)fList->UncheckedAt(kESDPtRecoError))->Fill(pt,100.*sigmaPt/pt);
    
    // slopeX-Y error at vertex
    ((TH2F*)fList->UncheckedAt(kESDSlopeXRecoError))->Fill(p,TMath::Sqrt(paramCovP(1,1)));
    ((TH2F*)fList->UncheckedAt(kESDSlopeYRecoError))->Fill(p,TMath::Sqrt(paramCovP(3,3)));
    
    // p*DCA
    AliMUONTrackParam trackParamAtDCA(trackParamAtFirstCluster);
    AliMUONTrackExtrap::ExtrapToVertexWithoutBranson(&trackParamAtDCA, esdTrack->GetZ());
    const TMatrixD& covDCA = trackParamAtDCA.GetCovariances();
    Double_t xDCA = esdTrack->GetNonBendingCoorAtDCA();
    Double_t yDCA = esdTrack->GetBendingCoorAtDCA();
    Double_t pDCA = esdTrack->PAtDCA();
    Double_t sigmaDCA = TMath::Sqrt(xDCA*xDCA*covDCA(0,0) + yDCA*yDCA*covDCA(2,2) + 2.*xDCA*yDCA*covDCA(0,2)) / esdTrack->GetDCA();
    ((TH1F*)fList->UncheckedAt(kESDPDCARecoError))->Fill(p,0.5*(pDCA+pU)*sigmaDCA);
    
    // loop over clusters
    AliMUONTrackParam* trackParam = static_cast<AliMUONTrackParam*>(track.GetTrackParamAtCluster()->First());
    while (trackParam) {
      
      AliMUONVCluster* cluster = trackParam->GetClusterPtr();
      Int_t chId = cluster->GetChamberId();
      Int_t deID = cluster->GetDetElemId();
      Double_t residualX = cluster->GetX() - trackParam->GetNonBendingCoor();
      Double_t residualY = cluster->GetY() - trackParam->GetBendingCoor();
      Double_t sigmaResidualX2 = cluster->GetErrX2() - trackParam->GetCovariances()(0,0);
      Double_t sigmaResidualY2 = cluster->GetErrY2() - trackParam->GetCovariances()(2,2);
      Double_t localChi2X = (sigmaResidualX2 > 0.) ? residualX*residualX/sigmaResidualX2 : 0.;
      Double_t localChi2Y = (sigmaResidualY2 > 0.) ? residualY*residualY/sigmaResidualY2 : 0.;
      Double_t localChi2 = 0.5 * trackParam->GetLocalChi2();
      
      ((TH1F*)fListExpert->UncheckedAt(kESDClusterHitMap+chId))->Fill(cluster->GetX(), cluster->GetY());
      
      ((TH1F*)fListExpert->UncheckedAt(kESDnTotClustersPerCh))->Fill(chId);
      ((TH1F*)fListExpert->UncheckedAt(kESDnTotClustersPerDE))->Fill(deID);
      
      ((TH1F*)fListExpert->UncheckedAt(kESDClusterChargeInCh+chId))->Fill(cluster->GetCharge());
      ((TH1F*)fListExpert->UncheckedAt(kESDSumClusterChargePerDE))->Fill(deID, cluster->GetCharge());
      
      if (cluster->GetNDigits() > 0) { // discard clusters with pad not stored in ESD
	((TH1F*)fListExpert->UncheckedAt(kESDnTotFullClustersPerDE))->Fill(deID);
        ((TH1F*)fListExpert->UncheckedAt(kESDClusterSizeInCh+chId))->Fill(cluster->GetNDigits());
	((TH1F*)fListExpert->UncheckedAt(kESDSumClusterSizePerDE))->Fill(deID, cluster->GetNDigits());
      }
      
      ((TH1F*)fListExpert->UncheckedAt(kESDResidualXInCh+chId))->Fill(residualX);
      ((TH1F*)fListExpert->UncheckedAt(kESDResidualYInCh+chId))->Fill(residualY);
      ((TH1F*)fListExpert->UncheckedAt(kESDSumResidualXPerDE))->Fill(deID, residualX);
      ((TH1F*)fListExpert->UncheckedAt(kESDSumResidualYPerDE))->Fill(deID, residualY);
      ((TH1F*)fListExpert->UncheckedAt(kESDSumResidualX2PerDE))->Fill(deID, residualX*residualX);
      ((TH1F*)fListExpert->UncheckedAt(kESDSumResidualY2PerDE))->Fill(deID, residualY*residualY);
      
      ((TH1F*)fListExpert->UncheckedAt(kESDLocalChi2XInCh+chId))->Fill(localChi2X);
      ((TH1F*)fListExpert->UncheckedAt(kESDLocalChi2YInCh+chId))->Fill(localChi2Y);
      ((TH1F*)fListExpert->UncheckedAt(kESDLocalChi2InCh+chId))->Fill(localChi2);
      ((TH1F*)fListExpert->UncheckedAt(kESDSumLocalChi2XPerDE))->Fill(deID, localChi2X);
      ((TH1F*)fListExpert->UncheckedAt(kESDSumLocalChi2YPerDE))->Fill(deID, localChi2Y);
      ((TH1F*)fListExpert->UncheckedAt(kESDSumLocalChi2PerDE))->Fill(deID, localChi2);
      
      trackParam = static_cast<AliMUONTrackParam*>(track.GetTrackParamAtCluster()->After(trackParam));
    }
    
  }
  
  ((TH1F*)fList->UncheckedAt(kESDnTracks))->Fill(nTrackerTracks);
  ((TH1F*)fList->UncheckedAt(kESDMatchTrig))->Fill(nTrackMatchTrig);
  
  if (nTracks > 0) {
    
    if (nTriggerTracks < 10) {
      fcurrentNEvnt++;
      fTotalNEvnt++;
      
      // fill event counters
      fEventCounters->Count(Form("event:muon/trigger:any/run:%d/%s", fCurrentRunNumber, selected.Data()));
      
      if (fFullPrintout) cout<<Form(" %d\t%4d\t\t%4d\t\t%4d\t", fESD->GetEventNumberInFile(), nTrackerTracks, nTriggerTracks, nTrackMatchTrig);
      
      fcurrentNTracks[0] += nTrackerTracks;
      fcurrentNTracks[1] += nTriggerTracks;
      fcurrentNTracks[2] += nTrackMatchTrig;
      fTotalNTracks[0] += nTrackerTracks;
      fTotalNTracks[1] += nTriggerTracks;
      fTotalNTracks[2] += nTrackMatchTrig;
      
      // fill track counters
      fTrackCounters->Count(Form("track:tracker/trigger:any/run:%d/%s", fCurrentRunNumber, selected.Data()), nTrackerTracks);
      fTrackCounters->Count(Form("track:trigger/trigger:any/run:%d/%s", fCurrentRunNumber, selected.Data()), nTriggerTracks);
      fTrackCounters->Count(Form("track:matched/trigger:any/run:%d/%s", fCurrentRunNumber, selected.Data()), nTrackMatchTrig);
      fTrackCounters->Count(Form("track:any/trigger:any/run:%d/%s", fCurrentRunNumber, selected.Data()), nTrackerTracks+nTriggerTracks);
      
      Bool_t triggerFiredForTrack = kFALSE;
      for (Int_t i=0; i<AliAnalysisTaskESDCheck::fgkNTriggerClass; i++) {
	
	if (fESD->IsTriggerClassFired(AliAnalysisTaskESDCheck::fgkTriggerClass[i])) {
	  
	  fCurrentNTrig[i]++;
	  fTotalNTrig[i]++;
	  
	  // fill event counters
	  fEventCounters->Count(Form("event:muon/trigger:%s/run:%d/%s", AliAnalysisTaskESDCheck::fgkTriggerShortName[i], fCurrentRunNumber, selected.Data()));
	  
	  fCurrentStat[i][0] += nTrackerTracks;
	  fCurrentStat[i][1] += nTriggerTracks;
	  fCurrentStat[i][2] += nTrackMatchTrig;
	  fTotalStat[i][0] += nTrackerTracks;
	  fTotalStat[i][1] += nTriggerTracks;
	  fTotalStat[i][2] += nTrackMatchTrig;
	  
	  // fill track counters
	  fTrackCounters->Count(Form("track:tracker/trigger:%s/run:%d/%s", AliAnalysisTaskESDCheck::fgkTriggerShortName[i], fCurrentRunNumber, selected.Data()), nTrackerTracks);
	  fTrackCounters->Count(Form("track:trigger/trigger:%s/run:%d/%s", AliAnalysisTaskESDCheck::fgkTriggerShortName[i], fCurrentRunNumber, selected.Data()), nTriggerTracks);
	  fTrackCounters->Count(Form("track:matched/trigger:%s/run:%d/%s", AliAnalysisTaskESDCheck::fgkTriggerShortName[i], fCurrentRunNumber, selected.Data()), nTrackMatchTrig);
	  fTrackCounters->Count(Form("track:any/trigger:%s/run:%d/%s", AliAnalysisTaskESDCheck::fgkTriggerShortName[i], fCurrentRunNumber, selected.Data()), nTrackerTracks+nTriggerTracks);
	  
	  triggerFiredForTrack = kTRUE;
	  
	  if (fFullPrintout) cout<<"\t  x";
	} else if (fFullPrintout) cout<<"\t  .";
	
      }
      
      if (!triggerFiredForTrack) {
	
	fCurrentNTrig[10]++;
	fTotalNTrig[10]++;
	
	// fill event counters
	fEventCounters->Count(Form("event:muon/trigger:other/run:%d/%s", fCurrentRunNumber, selected.Data()));
	
	fCurrentStat[10][0] += nTrackerTracks;
	fCurrentStat[10][1] += nTriggerTracks;
	fCurrentStat[10][2] += nTrackMatchTrig;
	fTotalStat[10][0] += nTrackerTracks;
	fTotalStat[10][1] += nTriggerTracks;
	fTotalStat[10][2] += nTrackMatchTrig;
	
	// fill counters
	fTrackCounters->Count(Form("track:tracker/trigger:Other/run:%d/%s", fCurrentRunNumber, selected.Data()), nTrackerTracks);
	fTrackCounters->Count(Form("track:trigger/trigger:Other/run:%d/%s", fCurrentRunNumber, selected.Data()), nTriggerTracks);
	fTrackCounters->Count(Form("track:matched/trigger:Other/run:%d/%s", fCurrentRunNumber, selected.Data()), nTrackMatchTrig);
	fTrackCounters->Count(Form("track:any/trigger:Other/run:%d/%s", fCurrentRunNumber, selected.Data()), nTrackerTracks+nTriggerTracks);
	
      }
      
      if (nTrackerTracks>1 && fFullPrintout) cout<<"\t\tDIMUON!!!";
      
      if (fFullPrintout) cout<<endl;
    } else if (fFullPrintout) cout<<Form(" %d\ttrigger readout problem (%d tracks)\n", fESD->GetEventNumberInFile(), nTriggerTracks);
    
  }
  
  // double loop over tracks to look for connection
  AliMUONTrack track1, track2;
  for (Int_t iTrack1 = 0; iTrack1 < nTracks; ++iTrack1) {
    AliESDMuonTrack* esdTrack1 = fESD->GetMuonTrack(iTrack1);
    
    if (!esdTrack1->ContainTrackerData()) continue;
    
    AliMUONESDInterface::ESDToMUON(*esdTrack1,track1,kFALSE);
    
    for (Int_t iTrack2 = iTrack1+1; iTrack2 < nTracks; ++iTrack2) {
      AliESDMuonTrack* esdTrack2 = fESD->GetMuonTrack(iTrack2);
      
      if (!esdTrack2->ContainTrackerData()) continue;
      
      AliMUONESDInterface::ESDToMUON(*esdTrack2,track2,kFALSE);
      
      Bool_t connected = kFALSE;
      for (Int_t i=0; i<track1.GetNClusters(); i++) {
	UInt_t clId1 = ((AliMUONTrackParam*)track1.GetTrackParamAtCluster()->UncheckedAt(i))->GetClusterPtr()->GetUniqueID();
	for (Int_t j=0; j<track2.GetNClusters(); j++) {
	  UInt_t clId2 = ((AliMUONTrackParam*)track2.GetTrackParamAtCluster()->UncheckedAt(j))->GetClusterPtr()->GetUniqueID();
	  if (clId1 == clId2) {
	    cout<<Form(" %d\tWarning: connected tracks!\n", fESD->GetEventNumberInFile());
	    connected = kTRUE;
	    break;
	  }
	}
	if (connected) break;
      }
      
    }
    
  }
    
  // Post final data. It will be written to a file with option "RECREATE"
  PostData(1, fList);
  PostData(2, fListExpert);
  PostData(3, fTrackCounters);
  PostData(4, fEventCounters);
}

//________________________________________________________________________
void AliAnalysisTaskESDCheck::NotifyRun() {
  //
  /// Called each time the run change
  /// print global statistics of the previous run and reset counters
  //
  
  // load necessary data from OCDB
  //fOCDBLoaded = kFALSE;
  if (!fOCDBLoaded) {
    AliCDBManager::Instance()->SetRun(fCurrentRunNumber);
    if (!AliMUONCDB::LoadMapping()) return;
    if (!AliMUONCDB::LoadField()) return;
    AliMUONRecoParam* recoParam = AliMUONCDB::LoadRecoParam();
    if (!recoParam) return;
    AliMUONESDInterface::ResetTracker(recoParam);
    AliGeomManager::LoadGeometry();
    if (!AliGeomManager::GetGeometry()) return;
    fOCDBLoaded = kTRUE;
  }
  
  // print statistic of previous run
  if (fcurrentNTotEvnt > 0) {
    PrintCurrentStat();
    cout<<"statistics without selection:"<<endl;
    fEventCounters->Print("trigger/event",Form("run:%d", fPreviousRun));
    fTrackCounters->Print("trigger/track",Form("run:%d", fPreviousRun));
    cout<<"statistics of selected events:"<<endl;
    fEventCounters->Print("trigger/event",Form("run:%d/selected:yes", fPreviousRun));
    fTrackCounters->Print("trigger/track",Form("run:%d/selected:yes", fPreviousRun));
  }
  
  fPreviousRun = fCurrentRunNumber;
  
  // reset statistics
  for (Int_t i=0; i<3; i++) fcurrentNTracks[i] = 0;
  for (Int_t i=0; i<11; i++) {
    fCurrentNTrig[i] = 0;
    fCurrentNTotTrig[i] = 0;
    for (Int_t j=0; j<3; j++) fCurrentStat[i][j] = 0;
  }
  fcurrentNEvnt = 0;
  fcurrentNTotEvnt = 0;
  
  if (fFullPrintout) {
    cout<<endl<<"run "<<fCurrentRunNumber<<":"<<endl;
    cout<<"event\tnTracker\tnTrigger\tnMatched";
    for (Int_t i=0; i<10; i++) cout<<Form("\t%s",AliAnalysisTaskESDCheck::fgkTriggerShortName[i]);
    cout<<endl;
  }
  
}

//________________________________________________________________________
void AliAnalysisTaskESDCheck::Terminate(Option_t *) {
  //
  /// Normalize histograms
  /// Draw result to the screen
  /// Called once at the end of the query.
  //
  
  // recover output objects
  fList = static_cast<TObjArray*> (GetOutputData(1));
  fListExpert = static_cast<TObjArray*> (GetOutputData(2));
  if (!fList || !fListExpert) return;
  fTrackCounters = static_cast<AliCounterCollection*> (GetOutputData(3));
  fEventCounters = static_cast<AliCounterCollection*> (GetOutputData(4));
  /*
  fTest = static_cast<TH1F*> (GetOutputData(5));
  fTest->Add((TH1F*)fList->UncheckedAt(kESDnTracks));
  */
  // global statistic
  if (fTrackCounters && fEventCounters) {
    PrintCurrentStat();
    cout<<"statistics without selection:"<<endl;
    fEventCounters->Print("trigger/event",Form("run:%d", fPreviousRun));
    fTrackCounters->Print("trigger/track",Form("run:%d", fPreviousRun));
    cout<<"statistics of selected events:"<<endl;
    fEventCounters->Print("trigger/event",Form("run:%d/selected:yes", fPreviousRun));
    fTrackCounters->Print("trigger/track",Form("run:%d/selected:yes", fPreviousRun));
    PrintTotalStat();
    cout<<"whole statistics without selection:"<<endl;
    fEventCounters->Print("trigger/event");
    fTrackCounters->Print("trigger/track");
    cout<<"whole statistics of selected events:"<<endl;
    fEventCounters->Print("trigger/event","selected:yes");
    fTrackCounters->Print("trigger/track","selected:yes");
    
    if (!gROOT->IsBatch()) {
      new TCanvas();
      fEventCounters->Draw("event","trigger","");
      new TCanvas();
      fTrackCounters->Draw("track","trigger","");
      new TCanvas();
      fEventCounters->Draw("event","trigger","selected:yes");
      new TCanvas();
      fTrackCounters->Draw("track","trigger","selected:yes");
    }
  }
  
  Double_t nTracks = ((TH1F*)fList->UncheckedAt(kESDnClustersPerTrack))->GetEntries();
  if (nTracks <= 0) {
    AliInfo("no track found");
    return;
  }
  
  // compute expected momentum resolution from sagitta and fill graphs
  static const Double_t k = 750.; // p = k/sagitta (100 GeV/c muon have ~ 7.5 mm sagitta)
  static const Double_t sigmaS = 4.; // 4 mm sagitta resolution
  // need geometry for energy loss calculation
  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  AliCDBManager::Instance()->SetRun(0);
  if (!AliMUONCDB::LoadMapping()) return;
  AliGeomManager::LoadGeometry();
  if (!AliGeomManager::GetGeometry()) return;
  // fill graphs
  TH1* hESDMomentumUncorrected = ((TH1*)fList->UncheckedAt(kESDMomentumUncorrected));
  TGraphAsymmErrors* gESDMomentumUncorrected = ((TGraphAsymmErrors*)fList->UncheckedAt(kESDMomentumUncorrectedG));
  TGraph* gESDMomentumError = ((TGraph*)fList->UncheckedAt(kESDMomentumErrorG));
  TGraph* gESDMomentumRelativeError = ((TGraph*)fList->UncheckedAt(kESDMomentumRelativeErrorG));
  TGraph* gESDMomentumRelativeErrorMin = ((TGraph*)fList->UncheckedAt(kESDMomentumRelativeErrorMinG));
  TGraph* gESDMomentumRelativeErrorMax = ((TGraph*)fList->UncheckedAt(kESDMomentumRelativeErrorMaxG));
  AliMUONTrackParam param;
  param.SetNonBendingCoor(50.);
  param.SetBendingCoor(50.);
  param.SetZ(-530.);
  for(Int_t ib = 1; ib<=hESDMomentumUncorrected->GetNbinsX(); ib++) {
    Double_t pU = hESDMomentumUncorrected->GetBinCenter(ib);
    Double_t pUErr = 0.45 * sigmaS / k * pU * pU; // 0.45 because overestimated by a factor ~ 2
    Double_t pUErrLow = 0.45 * (pU - k / (k/pU + sigmaS));
    Double_t pUErrHigh = 0.45 * (k / TMath::Max((k/pU - sigmaS),2.) - pU);
    gESDMomentumUncorrected->SetPoint(ib-1,pU,hESDMomentumUncorrected->GetBinContent(ib));
    gESDMomentumUncorrected->SetPointError(ib-1,TMath::Max(pUErrLow,5.),TMath::Max(pUErrHigh,5.),hESDMomentumUncorrected->GetBinError(ib),hESDMomentumUncorrected->GetBinError(ib));
    Double_t p = pU;
    param.SetInverseBendingMomentum(1./pU);
    p += AliMUONTrackExtrap::TotalMomentumEnergyLoss(&param, 0., 0., 0.);
    gESDMomentumError->SetPoint(ib-1,p,pUErrLow);
    gESDMomentumRelativeError->SetPoint(ib-1,p,100.*pUErr/p);
    gESDMomentumRelativeErrorMin->SetPoint(ib-1,p,100.*pUErrLow/p);
    gESDMomentumRelativeErrorMax->SetPoint(ib-1,p,100.*pUErrHigh/p);
  }
  gESDMomentumUncorrected->SetFillColor(5);
  
  // normalize histograms and fill summary plots
  TH1* hESDnClustersPerCh = ((TH1F*)fList->UncheckedAt(kESDnClustersPerCh));
  TH1* hESDnTotClustersPerCh = ((TH1F*)fListExpert->UncheckedAt(kESDnTotClustersPerCh));
  TH1* hESDnClustersPerDE = ((TH1F*)fList->UncheckedAt(kESDnClustersPerDE));
  TH1* hESDnTotClustersPerDE = ((TH1F*)fListExpert->UncheckedAt(kESDnTotClustersPerDE));
  TH1* hESDnTotFullClustersPerDE = ((TH1F*)fListExpert->UncheckedAt(kESDnTotFullClustersPerDE));
  TH1* hESDClusterChargePerChMean = ((TH1F*)fList->UncheckedAt(kESDClusterChargePerChMean));
  TH1* hESDClusterChargePerChSigma = ((TH1F*)fList->UncheckedAt(kESDClusterChargePerChSigma));
  TH1* hESDClusterSizePerChMean = ((TH1F*)fList->UncheckedAt(kESDClusterSizePerChMean));
  TH1* hESDClusterSizePerChSigma = ((TH1F*)fList->UncheckedAt(kESDClusterSizePerChSigma));
  TH1* hESDResidualXPerChMean = ((TH1F*)fList->UncheckedAt(kESDResidualXPerChMean));
  TH1* hESDResidualXPerChSigma = ((TH1F*)fList->UncheckedAt(kESDResidualXPerChSigma));
  TH1* hESDResidualYPerChMean = ((TH1F*)fList->UncheckedAt(kESDResidualYPerChMean));
  TH1* hESDResidualYPerChSigma = ((TH1F*)fList->UncheckedAt(kESDResidualYPerChSigma));
  TH1* hESDLocalChi2XPerChMean = ((TH1F*)fList->UncheckedAt(kESDLocalChi2XPerChMean));
  TH1* hESDLocalChi2YPerChMean = ((TH1F*)fList->UncheckedAt(kESDLocalChi2YPerChMean));
  TH1* hESDLocalChi2PerChMean = ((TH1F*)fList->UncheckedAt(kESDLocalChi2PerChMean));
  TH1* hESDSumClusterChargePerDE = ((TH1F*)fListExpert->UncheckedAt(kESDSumClusterChargePerDE));
  TH1* hESDClusterChargePerDE = ((TH1F*)fList->UncheckedAt(kESDClusterChargePerDE));
  TH1* hESDSumClusterSizePerDE = ((TH1F*)fListExpert->UncheckedAt(kESDSumClusterSizePerDE));
  TH1* hESDClusterSizePerDE = ((TH1F*)fList->UncheckedAt(kESDClusterSizePerDE));
  TH1* hESDSumResidualXPerDE = ((TH1F*)fListExpert->UncheckedAt(kESDSumResidualXPerDE));
  TH1* hESDSumResidualX2PerDE = ((TH1F*)fListExpert->UncheckedAt(kESDSumResidualX2PerDE));
  TH1* hESDResidualXPerDEMean = ((TH1F*)fList->UncheckedAt(kESDResidualXPerDEMean));
  TH1* hESDResidualXPerDESigma = ((TH1F*)fList->UncheckedAt(kESDResidualXPerDESigma));
  TH1* hESDSumResidualYPerDE = ((TH1F*)fListExpert->UncheckedAt(kESDSumResidualYPerDE));
  TH1* hESDSumResidualY2PerDE = ((TH1F*)fListExpert->UncheckedAt(kESDSumResidualY2PerDE));
  TH1* hESDResidualYPerDEMean = ((TH1F*)fList->UncheckedAt(kESDResidualYPerDEMean));
  TH1* hESDResidualYPerDESigma = ((TH1F*)fList->UncheckedAt(kESDResidualYPerDESigma));
  TH1* hESDSumLocalChi2XPerDE = ((TH1F*)fListExpert->UncheckedAt(kESDSumLocalChi2XPerDE));
  TH1* hESDLocalChi2XPerDEMean = ((TH1F*)fList->UncheckedAt(kESDLocalChi2XPerDEMean));
  TH1* hESDSumLocalChi2YPerDE = ((TH1F*)fListExpert->UncheckedAt(kESDSumLocalChi2YPerDE));
  TH1* hESDLocalChi2YPerDEMean = ((TH1F*)fList->UncheckedAt(kESDLocalChi2YPerDEMean));
  TH1* hESDSumLocalChi2PerDE = ((TH1F*)fListExpert->UncheckedAt(kESDSumLocalChi2PerDE));
  TH1* hESDLocalChi2PerDEMean = ((TH1F*)fList->UncheckedAt(kESDLocalChi2PerDEMean));
  
  hESDnClustersPerCh->Add(hESDnTotClustersPerCh, 1./nTracks);
  hESDnClustersPerDE->Add(hESDnTotClustersPerDE, 1./nTracks);
  
  // loop over chambers
  for (Int_t iCh = 0; iCh < AliMUONConstants::NTrackingCh(); iCh++) {
    
    TH1* hESDClusterChargeInCh = ((TH1F*)fListExpert->UncheckedAt(kESDClusterChargeInCh+iCh));
    Double_t sigmaCharge = hESDClusterChargeInCh->GetRMS();
    hESDClusterChargePerChMean->SetBinContent(iCh+1, hESDClusterChargeInCh->GetMean());
    hESDClusterChargePerChMean->SetBinError(iCh+1, hESDClusterChargeInCh->GetMeanError());
    hESDClusterChargePerChSigma->SetBinContent(iCh+1, sigmaCharge);
    hESDClusterChargePerChSigma->SetBinError(iCh+1, hESDClusterChargeInCh->GetRMSError());
    
    TH1* hESDClusterSizeInCh = ((TH1F*)fListExpert->UncheckedAt(kESDClusterSizeInCh+iCh));
    Double_t sigmaSize = hESDClusterSizeInCh->GetRMS();
    hESDClusterSizePerChMean->SetBinContent(iCh+1, hESDClusterSizeInCh->GetMean());
    hESDClusterSizePerChMean->SetBinError(iCh+1, hESDClusterSizeInCh->GetMeanError());
    hESDClusterSizePerChSigma->SetBinContent(iCh+1, sigmaSize);
    hESDClusterSizePerChSigma->SetBinError(iCh+1, hESDClusterSizeInCh->GetRMSError());
    
    TH1* hESDResidualXInCh = ((TH1F*)fListExpert->UncheckedAt(kESDResidualXInCh+iCh));
    hESDResidualXInCh->GetXaxis()->SetRangeUser(-3.*hESDResidualXInCh->GetRMS(), 3.*hESDResidualXInCh->GetRMS());
    Double_t sigmaResidualX = hESDResidualXInCh->GetRMS();
    hESDResidualXPerChMean->SetBinContent(iCh+1, hESDResidualXInCh->GetMean());
    hESDResidualXPerChMean->SetBinError(iCh+1, hESDResidualXInCh->GetMeanError());
    hESDResidualXPerChSigma->SetBinContent(iCh+1, sigmaResidualX);
    hESDResidualXPerChSigma->SetBinError(iCh+1, hESDResidualXInCh->GetRMSError());
    hESDResidualXInCh->GetXaxis()->SetRange(0,0);
    
    TH1* hESDResidualYInCh = ((TH1F*)fListExpert->UncheckedAt(kESDResidualYInCh+iCh));
    hESDResidualYInCh->GetXaxis()->SetRangeUser(-3.*hESDResidualYInCh->GetRMS(), 3.*hESDResidualYInCh->GetRMS());
    Double_t sigmaResidualY = hESDResidualYInCh->GetRMS();
    hESDResidualYPerChMean->SetBinContent(iCh+1, hESDResidualYInCh->GetMean());
    hESDResidualYPerChMean->SetBinError(iCh+1, hESDResidualYInCh->GetMeanError());
    hESDResidualYPerChSigma->SetBinContent(iCh+1, sigmaResidualY);
    hESDResidualYPerChSigma->SetBinError(iCh+1, hESDResidualYInCh->GetRMSError());
    hESDResidualYInCh->GetXaxis()->SetRange(0,0);
    
    TH1* hESDLocalChi2XInCh = ((TH1F*)fListExpert->UncheckedAt(kESDLocalChi2XInCh+iCh));
    Double_t sigmaLocalChi2X = hESDLocalChi2XInCh->GetRMS();
    hESDLocalChi2XPerChMean->SetBinContent(iCh+1, hESDLocalChi2XInCh->GetMean());
    hESDLocalChi2XPerChMean->SetBinError(iCh+1, hESDLocalChi2XInCh->GetMeanError());
    
    TH1* hESDLocalChi2YInCh = ((TH1F*)fListExpert->UncheckedAt(kESDLocalChi2YInCh+iCh));
    Double_t sigmaLocalChi2Y = hESDLocalChi2YInCh->GetRMS();
    hESDLocalChi2YPerChMean->SetBinContent(iCh+1, hESDLocalChi2YInCh->GetMean());
    hESDLocalChi2YPerChMean->SetBinError(iCh+1, hESDLocalChi2YInCh->GetMeanError());
    
    TH1* hESDLocalChi2InCh = ((TH1F*)fListExpert->UncheckedAt(kESDLocalChi2InCh+iCh));
    Double_t sigmaLocalChi2 = hESDLocalChi2InCh->GetRMS();
    hESDLocalChi2PerChMean->SetBinContent(iCh+1, hESDLocalChi2InCh->GetMean());
    hESDLocalChi2PerChMean->SetBinError(iCh+1, hESDLocalChi2InCh->GetMeanError());
    
    // loop over DE into chamber iCh
    AliMpDEIterator it;
    it.First(iCh);
    while ( !it.IsDone()) {
      
      Int_t iDE = it.CurrentDEId();
      
      Double_t nClusters = hESDnTotClustersPerDE->GetBinContent(iDE+1);
      if (nClusters > 1) {
	
	hESDClusterChargePerDE->SetBinContent(iDE+1, hESDSumClusterChargePerDE->GetBinContent(iDE+1)/nClusters);
	hESDClusterChargePerDE->SetBinError(iDE+1, sigmaCharge/TMath::Sqrt(nClusters));
	
	Double_t meanResX = hESDSumResidualXPerDE->GetBinContent(iDE+1)/nClusters;
	hESDResidualXPerDEMean->SetBinContent(iDE+1, meanResX);
	hESDResidualXPerDEMean->SetBinError(iDE+1, sigmaResidualX/TMath::Sqrt(nClusters));
	hESDResidualXPerDESigma->SetBinContent(iDE+1, TMath::Sqrt(hESDSumResidualX2PerDE->GetBinContent(iDE+1)/nClusters - meanResX*meanResX));
	hESDResidualXPerDESigma->SetBinError(iDE+1, sigmaResidualX/TMath::Sqrt(2.*nClusters));
	
	Double_t meanResY = hESDSumResidualYPerDE->GetBinContent(iDE+1)/nClusters;
	hESDResidualYPerDEMean->SetBinContent(iDE+1, meanResY);
	hESDResidualYPerDEMean->SetBinError(iDE+1, sigmaResidualY/TMath::Sqrt(nClusters));
	hESDResidualYPerDESigma->SetBinContent(iDE+1, TMath::Sqrt(hESDSumResidualY2PerDE->GetBinContent(iDE+1)/nClusters - meanResY*meanResY));
	hESDResidualYPerDESigma->SetBinError(iDE+1, sigmaResidualY/TMath::Sqrt(2.*nClusters));
	
	hESDLocalChi2XPerDEMean->SetBinContent(iDE+1, hESDSumLocalChi2XPerDE->GetBinContent(iDE+1)/nClusters);
	hESDLocalChi2XPerDEMean->SetBinError(iDE+1, sigmaLocalChi2X/TMath::Sqrt(nClusters));
	
	hESDLocalChi2YPerDEMean->SetBinContent(iDE+1, hESDSumLocalChi2YPerDE->GetBinContent(iDE+1)/nClusters);
	hESDLocalChi2YPerDEMean->SetBinError(iDE+1, sigmaLocalChi2Y/TMath::Sqrt(nClusters));
	
	hESDLocalChi2PerDEMean->SetBinContent(iDE+1, hESDSumLocalChi2PerDE->GetBinContent(iDE+1)/nClusters);
	hESDLocalChi2PerDEMean->SetBinError(iDE+1, sigmaLocalChi2/TMath::Sqrt(nClusters));
	
      } else {
	
	hESDClusterChargePerDE->SetBinContent(iDE+1, hESDSumClusterChargePerDE->GetBinContent(iDE+1));
	hESDClusterChargePerDE->SetBinError(iDE+1, hESDClusterChargeInCh->GetXaxis()->GetXmax());
	
	hESDResidualXPerDEMean->SetBinContent(iDE+1, hESDSumResidualXPerDE->GetBinContent(iDE+1));
	hESDResidualXPerDEMean->SetBinError(iDE+1, hESDResidualXInCh->GetXaxis()->GetXmax());
	hESDResidualXPerDESigma->SetBinContent(iDE+1, 0.);
	hESDResidualXPerDESigma->SetBinError(iDE+1, hESDResidualXInCh->GetXaxis()->GetXmax());
	
	hESDResidualYPerDEMean->SetBinContent(iDE+1, hESDSumResidualYPerDE->GetBinContent(iDE+1));
	hESDResidualYPerDEMean->SetBinError(iDE+1, hESDResidualYInCh->GetXaxis()->GetXmax());
	hESDResidualYPerDESigma->SetBinContent(iDE+1, 0.);
	hESDResidualYPerDESigma->SetBinError(iDE+1, hESDResidualYInCh->GetXaxis()->GetXmax());
	
	hESDLocalChi2XPerDEMean->SetBinContent(iDE+1, hESDSumLocalChi2XPerDE->GetBinContent(iDE+1));
	hESDLocalChi2XPerDEMean->SetBinError(iDE+1, hESDLocalChi2XInCh->GetXaxis()->GetXmax());
	
	hESDLocalChi2YPerDEMean->SetBinContent(iDE+1, hESDSumLocalChi2YPerDE->GetBinContent(iDE+1));
	hESDLocalChi2YPerDEMean->SetBinError(iDE+1, hESDLocalChi2YInCh->GetXaxis()->GetXmax());
	
	hESDLocalChi2PerDEMean->SetBinContent(iDE+1, hESDSumLocalChi2PerDE->GetBinContent(iDE+1));
	hESDLocalChi2PerDEMean->SetBinError(iDE+1, hESDLocalChi2InCh->GetXaxis()->GetXmax());
	
      }
      
      Double_t nFullClusters = hESDnTotFullClustersPerDE->GetBinContent(iDE+1);
      if (nFullClusters > 1) {
	
	hESDClusterSizePerDE->SetBinContent(iDE+1, hESDSumClusterSizePerDE->GetBinContent(iDE+1)/nFullClusters);
	hESDClusterSizePerDE->SetBinError(iDE+1, sigmaSize/TMath::Sqrt(nFullClusters));
	
      } else {
	
	hESDClusterSizePerDE->SetBinContent(iDE+1, hESDSumClusterSizePerDE->GetBinContent(iDE+1));
	hESDClusterSizePerDE->SetBinError(iDE+1, hESDClusterSizeInCh->GetXaxis()->GetXmax());
	
      }
      
      it.Next();
    }
    
  }
  
  TFile *histoFile = new TFile("histo.root", "RECREATE");
  histoFile->mkdir("general","general");
  histoFile->cd("general");
  fList->Write();
  histoFile->mkdir("expert","expert");
  histoFile->cd("expert");
  fListExpert->Write();
  histoFile->Close();
  
  if (!gROOT->IsBatch()) new TBrowser();
  
}

//________________________________________________________________________
Double_t AliAnalysisTaskESDCheck::ChangeThetaRange(Double_t theta)
{
  if(theta < -2.5) return (theta / TMath::Pi() + 1.) * 180.;
  else if(theta > 2.5) return (theta / TMath::Pi() - 1.) * 180.;
  else return theta / TMath::Pi() * 180.;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskESDCheck::SelectTrack(AliMUONTrack &track)
{
  /// Select track according to cluster location
  
  if (fStSelect == kAll) return kTRUE;
  
  Int_t currentCh, previousCh = -1, nChHitInSt45 = 0;
  Bool_t clusterInSt[5];
  for (Int_t iSt = 0; iSt < 5; iSt++) clusterInSt[iSt] = kFALSE;
  
  AliMUONTrackParam* trackParam = static_cast<AliMUONTrackParam*>(track.GetTrackParamAtCluster()->First());
  while (trackParam) {
    
    currentCh = trackParam->GetClusterPtr()->GetChamberId();
    
    clusterInSt[currentCh/2] = kTRUE;
    
    if (currentCh > 5 && currentCh != previousCh) {
      nChHitInSt45++;
      previousCh = currentCh;
    }
    
    trackParam = static_cast<AliMUONTrackParam*>(track.GetTrackParamAtCluster()->After(trackParam));
  }
  
  if (fStSelect == k1245 && clusterInSt[0] && clusterInSt[1] && nChHitInSt45 >= 2) return kTRUE;
  
  else if (fStSelect == k245 && clusterInSt[1]/* && nChHitInSt45 >= 3*/) return kTRUE;
  
  return kFALSE;
  
}

//__________________________________________________________________________
void AliAnalysisTaskESDCheck::Cov2CovP(const TMatrixD &param, TMatrixD &cov)
{
  /// change coordinate system: (X, SlopeX, Y, SlopeY, q/Pyz) -> (X, SlopeX, Y, SlopeY, q*PTot)
  /// parameters (param) are given in the (X, SlopeX, Y, SlopeY, q/Pyz) coordinate system
  
  // charge * total momentum
  Double_t qPTot = TMath::Sqrt(1. + param(1,0)*param(1,0) + param(3,0)*param(3,0)) /
  TMath::Sqrt(1. + param(3,0)*param(3,0)) / param(4,0);
  
  // Jacobian of the opposite transformation
  TMatrixD jacob(5,5);
  jacob.UnitMatrix();
  jacob(4,1) = qPTot * param(1,0) / (1. + param(1,0)*param(1,0) + param(3,0)*param(3,0));
  jacob(4,3) = - qPTot * param(1,0) * param(1,0) * param(3,0) /
  (1. + param(3,0)*param(3,0)) / (1. + param(1,0)*param(1,0) + param(3,0)*param(3,0));
  jacob(4,4) = - qPTot / param(4,0);
  
  // compute covariances in new coordinate system
  TMatrixD tmp(cov,TMatrixD::kMultTranspose,jacob);
  cov.Mult(jacob,tmp);
}

//__________________________________________________________________________
void AliAnalysisTaskESDCheck::PrintCurrentStat()
{
  /// print the statistics of the current run if not empty
  if (fFullPrintout) {
    cout<<Form("summary\t%4d\t\t%4d\t\t%4d\t", fcurrentNTracks[0], fcurrentNTracks[1], fcurrentNTracks[2]);
    for (Int_t i=0; i<10; i++) cout<<Form("\t%3d", fCurrentNTrig[i]);
    cout<<endl;
  }
  cout<<endl;
  cout<<Form("trigger\tnEvTot\tnEvInMu\t nTrk\t nTrg\tnMatch")<<endl;
  for (Int_t i=0; i<11; i++) {
    cout<<Form("%s\t%6d\t%5d\t%4d\t%4d\t%4d", AliAnalysisTaskESDCheck::fgkTriggerShortName[i],
	       fCurrentNTotTrig[i], fCurrentNTrig[i], fCurrentStat[i][0], fCurrentStat[i][1], fCurrentStat[i][2])<<endl;
  }
  cout<<Form("Any\t%6d\t%5d\t%4d\t%4d\t%4d",
	     fcurrentNTotEvnt, fcurrentNEvnt, fcurrentNTracks[0], fcurrentNTracks[1], fcurrentNTracks[2])<<endl;
  cout<<endl;
}

//__________________________________________________________________________
void AliAnalysisTaskESDCheck::PrintTotalStat()
{
  /// print the overall statistics
  if (fFullPrintout) {
    cout<<Form("total\t%4d\t\t%4d\t\t%4d\t", fTotalNTracks[0], fTotalNTracks[1], fTotalNTracks[2]);
    for (Int_t i=0; i<10; i++) cout<<Form("\t%3d", fTotalNTrig[i]);
    cout<<endl;
  }
  cout<<endl;
  cout<<Form("trigger\tnEvTot\tnEvInMu\t nTrk\t nTrg\tnMatch")<<endl;
  for (Int_t i=0; i<11; i++) {
    cout<<Form("%s\t%6d\t%5d\t%4d\t%4d\t%4d", AliAnalysisTaskESDCheck::fgkTriggerShortName[i],
	       fTotalNTotTrig[i], fTotalNTrig[i], fTotalStat[i][0], fTotalStat[i][1], fTotalStat[i][2])<<endl;
  }
  cout<<Form("Any\t%6d\t%5d\t%4d\t%4d\t%4d",
	     fTotalNTotEvnt, fTotalNEvnt, fTotalNTracks[0], fTotalNTracks[1], fTotalNTracks[2])<<endl;
  cout<<endl;
}

