/*
 *  MUONResolution.C
 *  aliroot_dev
 *
 *  Created by philippe pillot on 11/07/13.
 *  Copyright 2013 SUBATECH. All rights reserved.
 *
 */

// ROOT includes
#include <Riostream.h>
#include "TMath.h"
#include "TLorentzVector.h"
#include "TMatrixD.h"
#include "TObjArray.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TF1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGeoManager.h"

// STEER includes
#include "AliCDBManager.h"
#include "AliGeomManager.h"
#include "AliLog.h"
#include "AliESDEvent.h"
#include "AliESDMuonTrack.h"

// MUON includes
#include "AliMUONCDB.h"
#include "AliMUONConstants.h"
#include "AliMUONTrack.h"
#include "AliMUONRecoCheck.h"
#include "AliMUONTrackParam.h"
#include "AliMUONRecoParam.h"
#include "AliMUONVTrackStore.h"
#include "AliMUONVCluster.h"
#include "AliMUONTrackExtrap.h"
#include "AliMUONESDInterface.h"
#include "AliMUONVTriggerTrackStore.h"
#include "AliMUONTriggerTrack.h"
#include "AliMpDEIterator.h"

TF1* fGaus = 0x0;
TF1* fPGaus = 0x0;

void FitGausResVsP(TH2* h, Int_t nBins, const Double_t mean0, const Double_t sigma0, const char* fitting, TGraphAsymmErrors* gMean, TGraphAsymmErrors* gSigma);
void FitPGausResVsP(TH2* h, Int_t nBins, const Double_t mean0, const Double_t sigma0, const char* fitting, TGraphAsymmErrors* gMean, TGraphAsymmErrors* gSigma);
Double_t CrystalBallExtended(Double_t *x, Double_t *par);


//------------------------------------------------------------------------------------
void MUONResolution()
{
  /// fit the JPsi mass peak
  /// get the track x,y resolution at vertex (DCA)
  /// get the track x,y resolution at the begining of the absorber
  
//  Double_t zAbsIn = AliMUONConstants::AbsZBeg();
  Double_t zAbsIn = -76.;
  AliMUONRecoCheck rc("AliESDs.root", "./generated/");
  
  // load necessary data from OCDB
  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  AliCDBManager::Instance()->SetSpecificStorage("GRP/GRP/Data", Form("local://%s",gSystem->pwd()));
  AliCDBManager::Instance()->SetRun(rc.GetRunNumber());
  if (!AliMUONCDB::LoadField()) return;
  AliMUONTrackExtrap::SetField();
  AliGeomManager::LoadGeometry("geometry.root");
  if (!AliGeomManager::GetGeometry()) return;
  AliMUONRecoParam* recoParam = AliMUONCDB::LoadRecoParam();
  if (!recoParam) return;
  AliMUONESDInterface::ResetTracker(recoParam);
  
  // get sigma cut from recoParam to associate clusters with TrackRefs in case the label are not used
  Double_t sigmaCut = (recoParam->ImproveTracks()) ? recoParam->GetSigmaCutForImprovement() : recoParam->GetSigmaCutForTracking();
  // compute the mask of requested stations from recoParam
  UInt_t requestedStationMask = 0;
  for (Int_t i = 0; i < 5; i++) if (recoParam->RequestStation(i)) requestedStationMask |= ( 1 << i );
  // get from recoParam whether a track need 2 chambers hit in the same station (4 or 5) or not to be reconstructible
  Bool_t request2ChInSameSt45 = !recoParam->MakeMoreTrackCandidates();
  
  // create output histograms
  TFile *histoFile = new TFile("MUONResolution.root", "RECREATE");
  
  TH1F* hMass = new TH1F("hMass", "#mu^{+}#mu^{-} invariant mass distribution;mass (GeV/c^{2})", 200, 2., 4.);
  hMass->SetDirectory(0);
  
  const Int_t pNBins = 30;
  const Double_t pEdges[2] = {0., 300.};
  
  // at DCA
  histoFile->mkdir("pDCA","pDCA");
  histoFile->cd("pDCA");
  
  const Int_t pDCANBins = 500;
  const Double_t pDCAEdges[2] = {0., 1000.};
  const Double_t pResDCAEdges[2] = {0., 200.};
  
  TH1D *hPDCAX_2_3_Deg = new TH1D("hPDCAX_2_3_Deg","p #times DCAX for tracks within [2,3[ degrees at absorber end;p #times DCAX (GeV #times cm)",2*pDCANBins, -pDCAEdges[1], pDCAEdges[1]);
  TH1D *hPDCAX_3_10_Deg = new TH1D("hPDCAX_3_10_Deg","p #times DCAX for tracks within [3,10[ degrees at absorber end;p #times DCAX (GeV #times cm)",2*pDCANBins, -pDCAEdges[1], pDCAEdges[1]);
  TH2D *hPDCAXVsP_2_3_Deg = new TH2D("hPDCAXVsP_2_3_Deg","p #times DCAX versus p for tracks within [2,3[ degrees at absorber end;p (GeV/c);p #times DCAX (GeV #times cm)",pNBins,pEdges[0],pEdges[1], 2*pDCANBins, -pDCAEdges[1], pDCAEdges[1]);
  TH2D *hPDCAXVsP_3_10_Deg = new TH2D("hPDCAXVsP_3_10_Deg","p #times DCAX versus p for tracks within [3,10[ degrees at absorber end;p (GeV/c);p #times DCAX (GeV #times cm)",pNBins,pEdges[0],pEdges[1], 2*pDCANBins, -pDCAEdges[1], pDCAEdges[1]);
  TGraphAsymmErrors* gMeanPDCAXVsP_2_3_Deg = new TGraphAsymmErrors(pNBins);
  gMeanPDCAXVsP_2_3_Deg->SetName("gMeanPDCAXVsP_2_3_Deg");
  gMeanPDCAXVsP_2_3_Deg->SetTitle("<p #times DCAX> versus p for tracks within [2,3[ degrees at absorber end;p (GeV/c);<p #times DCAX> (GeV #times cm)");
  TGraphAsymmErrors* gSigmaPDCAXVsP_2_3_Deg = new TGraphAsymmErrors(pNBins);
  gSigmaPDCAXVsP_2_3_Deg->SetName("gSigmaPDCAXVsP_2_3_Deg");
  gSigmaPDCAXVsP_2_3_Deg->SetTitle("#sigma_{p #times DCAX} versus p for tracks within [2,3[ degrees at absorber end;p (GeV/c);#sigma_{p #times DCAX} (GeV #times cm)");
  TGraphAsymmErrors* gMeanPDCAXVsP_3_10_Deg = new TGraphAsymmErrors(pNBins);
  gMeanPDCAXVsP_3_10_Deg->SetName("gMeanPDCAXVsP_3_10_Deg");
  gMeanPDCAXVsP_3_10_Deg->SetTitle("<p #times DCAX> versus p for tracks within [3,10[ degrees at absorber end;p (GeV/c);<p #times DCAX> (GeV #times cm)");
  TGraphAsymmErrors* gSigmaPDCAXVsP_3_10_Deg = new TGraphAsymmErrors(pNBins);
  gSigmaPDCAXVsP_3_10_Deg->SetName("gSigmaPDCAXVsP_3_10_Deg");
  gSigmaPDCAXVsP_3_10_Deg->SetTitle("#sigma_{p #times DCAX} versus p for tracks within [3,10[ degrees at absorber end;p (GeV/c);#sigma_{p #times DCAX} (GeV #times cm)");
  TH2D *hPResDCAXVsP_2_3_Deg = new TH2D("hPResDCAXVsP_2_3_Deg","p #times #sigma_{DCAX} versus p for tracks within [2,3[ degrees at absorber end;p (GeV/c);p #times #sigma_{DCAX} (GeV #times cm)",pNBins,pEdges[0],pEdges[1], pDCANBins, pResDCAEdges[0], pResDCAEdges[1]);
  TH2D *hPResDCAXVsP_3_10_Deg = new TH2D("hPResDCAXVsP_3_10_Deg","p #times #sigma_{DCAX} versus p for tracks within [3,10[ degrees at absorber end;p (GeV/c);p #times #sigma_{DCAX} (GeV #times cm)",pNBins,pEdges[0],pEdges[1], pDCANBins, pResDCAEdges[0], pResDCAEdges[1]);
  
  TH1D *hPDCAY_2_3_Deg = new TH1D("hPDCAY_2_3_Deg","p #times DCAY for tracks within [2,3[ degrees at absorber end;p #times DCAY (GeV #times cm)",2*pDCANBins, -pDCAEdges[1], pDCAEdges[1]);
  TH1D *hPDCAY_3_10_Deg = new TH1D("hPDCAY_3_10_Deg","p #times DCAY for tracks within [3,10[ degrees at absorber end;p #times DCAY (GeV #times cm)",2*pDCANBins, -pDCAEdges[1], pDCAEdges[1]);
  TH2D *hPDCAYVsP_2_3_Deg = new TH2D("hPDCAYVsP_2_3_Deg","p #times DCAY versus p for tracks within [2,3[ degrees at absorber end;p (GeV/c);p #times DCAY (GeV #times cm)",pNBins,pEdges[0],pEdges[1], 2*pDCANBins, -pDCAEdges[1], pDCAEdges[1]);
  TH2D *hPDCAYVsP_3_10_Deg = new TH2D("hPDCAYVsP_3_10_Deg","p #times DCAY versus p for tracks within [3,10[ degrees at absorber end;p (GeV/c);p #times DCAY (GeV #times cm)",pNBins,pEdges[0],pEdges[1], 2*pDCANBins, -pDCAEdges[1], pDCAEdges[1]);
  TGraphAsymmErrors* gMeanPDCAYVsP_2_3_Deg = new TGraphAsymmErrors(pNBins);
  gMeanPDCAYVsP_2_3_Deg->SetName("gMeanPDCAYVsP_2_3_Deg");
  gMeanPDCAYVsP_2_3_Deg->SetTitle("<p #times DCAY> versus p for tracks within [2,3[ degrees at absorber end;p (GeV/c);<p #times DCAY> (GeV #times cm)");
  TGraphAsymmErrors* gSigmaPDCAYVsP_2_3_Deg = new TGraphAsymmErrors(pNBins);
  gSigmaPDCAYVsP_2_3_Deg->SetName("gSigmaPDCAYVsP_2_3_Deg");
  gSigmaPDCAYVsP_2_3_Deg->SetTitle("#sigma_{p #times DCAY} versus p for tracks within [2,3[ degrees at absorber end;p (GeV/c);#sigma_{p #times DCAY} (GeV #times cm)");
  TGraphAsymmErrors* gMeanPDCAYVsP_3_10_Deg = new TGraphAsymmErrors(pNBins);
  gMeanPDCAYVsP_3_10_Deg->SetName("gMeanPDCAYVsP_3_10_Deg");
  gMeanPDCAYVsP_3_10_Deg->SetTitle("<p #times DCAY> versus p for tracks within [3,10[ degrees at absorber end;p (GeV/c);<p #times DCAY> (GeV #times cm)");
  TGraphAsymmErrors* gSigmaPDCAYVsP_3_10_Deg = new TGraphAsymmErrors(pNBins);
  gSigmaPDCAYVsP_3_10_Deg->SetName("gSigmaPDCAYVsP_3_10_Deg");
  gSigmaPDCAYVsP_3_10_Deg->SetTitle("#sigma_{p #times DCAY} versus p for tracks within [3,10[ degrees at absorber end;p (GeV/c);#sigma_{p #times DCAY} (GeV #times cm)");
  TH2D *hPResDCAYVsP_2_3_Deg = new TH2D("hPResDCAYVsP_2_3_Deg","p #times #sigma_{DCAY} versus p for tracks within [2,3[ degrees at absorber end;p (GeV/c);p #times #sigma_{DCAY} (GeV #times cm)",pNBins,pEdges[0],pEdges[1], pDCANBins, pResDCAEdges[0], pResDCAEdges[1]);
  TH2D *hPResDCAYVsP_3_10_Deg = new TH2D("hPResDCAYVsP_3_10_Deg","p #times #sigma_{DCAY} versus p for tracks within [3,10[ degrees at absorber end;p (GeV/c);p #times #sigma_{DCAY} (GeV #times cm)",pNBins,pEdges[0],pEdges[1], pDCANBins, pResDCAEdges[0], pResDCAEdges[1]);
  
  TH1D *hPDCA_2_3_Deg = new TH1D("hPDCA_2_3_Deg","p #times DCA for tracks within [2,3[ degrees at absorber end;p #times DCA (GeV #times cm)",pDCANBins, pDCAEdges[0], pDCAEdges[1]);
  TH1D *hPDCA_3_10_Deg = new TH1D("hPDCA_3_10_Deg","p #times DCA for tracks within [3,10[ degrees at absorber end;p #times DCA (GeV #times cm)",pDCANBins, pDCAEdges[0], pDCAEdges[1]);
  TH2D *hPDCAVsP_2_3_Deg = new TH2D("hPDCAVsP_2_3_Deg","p #times DCA versus p for tracks within [2,3[ degrees at absorber end;p (GeV/c);p #times DCA (GeV #times cm)",pNBins,pEdges[0],pEdges[1], pDCANBins, pDCAEdges[0], pDCAEdges[1]);
  TH2D *hPDCAVsP_3_10_Deg = new TH2D("hPDCAVsP_3_10_Deg","p #times DCA versus p for tracks within [3,10[ degrees at absorber end;p (GeV/c);p #times DCA (GeV #times cm)",pNBins,pEdges[0],pEdges[1], pDCANBins, pDCAEdges[0], pDCAEdges[1]);
  TGraphAsymmErrors* gMeanPDCAVsP_2_3_Deg = new TGraphAsymmErrors(pNBins);
  gMeanPDCAVsP_2_3_Deg->SetName("gMeanPDCAVsP_2_3_Deg");
  gMeanPDCAVsP_2_3_Deg->SetTitle("<p #times DCA> versus p for tracks within [2,3[ degrees at absorber end;p (GeV/c);<p #times DCA> (GeV #times cm)");
  TGraphAsymmErrors* gSigmaPDCAVsP_2_3_Deg = new TGraphAsymmErrors(pNBins);
  gSigmaPDCAVsP_2_3_Deg->SetName("gSigmaPDCAVsP_2_3_Deg");
  gSigmaPDCAVsP_2_3_Deg->SetTitle("#sigma_{p #times DCA} versus p for tracks within [2,3[ degrees at absorber end;p (GeV/c);#sigma_{p #times DCA} (GeV #times cm)");
  TGraphAsymmErrors* gMeanPDCAVsP_3_10_Deg = new TGraphAsymmErrors(pNBins);
  gMeanPDCAVsP_3_10_Deg->SetName("gMeanPDCAVsP_3_10_Deg");
  gMeanPDCAVsP_3_10_Deg->SetTitle("<p #times DCA> versus p for tracks within [3,10[ degrees at absorber end;p (GeV/c);<p #times DCA> (GeV #times cm)");
  TGraphAsymmErrors* gSigmaPDCAVsP_3_10_Deg = new TGraphAsymmErrors(pNBins);
  gSigmaPDCAVsP_3_10_Deg->SetName("gSigmaPDCAVsP_3_10_Deg");
  gSigmaPDCAVsP_3_10_Deg->SetTitle("#sigma_{p #times DCA} versus p for tracks within [3,10[ degrees at absorber end;p (GeV/c);#sigma_{p #times DCA} (GeV #times cm)");
  TH2D *hPResDCAVsP_2_3_Deg = new TH2D("hPResDCAVsP_2_3_Deg","p #times #sigma_{DCA} versus p for tracks within [2,3[ degrees at absorber end;p (GeV/c);p #times #sigma_{DCA} (GeV #times cm)",pNBins,pEdges[0],pEdges[1], pDCANBins, pResDCAEdges[0], pResDCAEdges[1]);
  TH2D *hPResDCAVsP_3_10_Deg = new TH2D("hPResDCAVsP_3_10_Deg","p #times #sigma_{DCA} versus p for tracks within [3,10[ degrees at absorber end;p (GeV/c);p #times #sigma_{DCA} (GeV #times cm)",pNBins,pEdges[0],pEdges[1], pDCANBins, pResDCAEdges[0], pResDCAEdges[1]);
  
  // at the begining of the absorber without branson correction
  histoFile->mkdir("pDXYAtAbsIn","pDXYAtAbsIn");
  histoFile->cd("pDXYAtAbsIn");
  
  TH1D *hPDXAbsIn_2_3_Deg = new TH1D("hPDXAbsIn_2_3_Deg","p #times DX at absorber begin for tracks within [2,3[ degrees at absorber end;p #times DX (GeV #times cm)",2*pDCANBins, -pDCAEdges[1], pDCAEdges[1]);
  TH1D *hPDXAbsIn_3_10_Deg = new TH1D("hPDXAbsIn_3_10_Deg","p #times DX at absorber begin for tracks within [3,10[ degrees at absorber end;p #times DX (GeV #times cm)",2*pDCANBins, -pDCAEdges[1], pDCAEdges[1]);
  TH2D *hPDXAbsInVsP_2_3_Deg = new TH2D("hPDXAbsInVsP_2_3_Deg","p #times DX at absorber begin versus p for tracks within [2,3[ degrees at absorber end;p (GeV/c);p #times DX (GeV #times cm)",pNBins,pEdges[0],pEdges[1], 2*pDCANBins, -pDCAEdges[1], pDCAEdges[1]);
  TH2D *hPDXAbsInVsP_3_10_Deg = new TH2D("hPDXAbsInVsP_3_10_Deg","p #times DX at absorber begin versus p for tracks within [3,10[ degrees at absorber end;p (GeV/c);p #times DX (GeV #times cm)",pNBins,pEdges[0],pEdges[1], 2*pDCANBins, -pDCAEdges[1], pDCAEdges[1]);
  TGraphAsymmErrors* gMeanPDXAbsInVsP_2_3_Deg = new TGraphAsymmErrors(pNBins);
  gMeanPDXAbsInVsP_2_3_Deg->SetName("gMeanPDXAbsInVsP_2_3_Deg");
  gMeanPDXAbsInVsP_2_3_Deg->SetTitle("<p #times DX at absorber begin> versus p for tracks within [2,3[ degrees at absorber end;p (GeV/c);<p #times DX> (GeV #times cm)");
  TGraphAsymmErrors* gSigmaPDXAbsInVsP_2_3_Deg = new TGraphAsymmErrors(pNBins);
  gSigmaPDXAbsInVsP_2_3_Deg->SetName("gSigmaPDXAbsInVsP_2_3_Deg");
  gSigmaPDXAbsInVsP_2_3_Deg->SetTitle("#sigma_{p #times DX at absorber begin} versus p for tracks within [2,3[ degrees at absorber end;p (GeV/c);#sigma_{p #times DX} (GeV #times cm)");
  TGraphAsymmErrors* gMeanPDXAbsInVsP_3_10_Deg = new TGraphAsymmErrors(pNBins);
  gMeanPDXAbsInVsP_3_10_Deg->SetName("gMeanPDXAbsInVsP_3_10_Deg");
  gMeanPDXAbsInVsP_3_10_Deg->SetTitle("<p #times DX at absorber begin> versus p for tracks within [3,10[ degrees at absorber end;p (GeV/c);<p #times DX> (GeV #times cm)");
  TGraphAsymmErrors* gSigmaPDXAbsInVsP_3_10_Deg = new TGraphAsymmErrors(pNBins);
  gSigmaPDXAbsInVsP_3_10_Deg->SetName("gSigmaPDXAbsInVsP_3_10_Deg");
  gSigmaPDXAbsInVsP_3_10_Deg->SetTitle("#sigma_{p #times DX at absorber begin} versus p for tracks within [3,10[ degrees at absorber end;p (GeV/c);#sigma_{p #times DX} (GeV #times cm)");
  TH2D *hPResDXAbsInVsP_2_3_Deg = new TH2D("hPResDXAbsInVsP_2_3_Deg","p #times #sigma_{DX} at absorber begin versus p for tracks within [2,3[ degrees at absorber end;p (GeV/c);p #times #sigma_{DX} (GeV #times cm)",pNBins,pEdges[0],pEdges[1], pDCANBins, pResDCAEdges[0], pResDCAEdges[1]);
  TH2D *hPResDXAbsInVsP_3_10_Deg = new TH2D("hPResDXAbsInVsP_3_10_Deg","p #times #sigma_{DX} at absorber begin versus p for tracks within [3,10[ degrees at absorber end;p (GeV/c);p #times #sigma_{DX} (GeV #times cm)",pNBins,pEdges[0],pEdges[1], pDCANBins, pResDCAEdges[0], pResDCAEdges[1]);
  
  TH1D *hPDYAbsIn_2_3_Deg = new TH1D("hPDYAbsIn_2_3_Deg","p #times DY at absorber begin for tracks within [2,3[ degrees at absorber end;p #times DY (GeV #times cm)",2*pDCANBins, -pDCAEdges[1], pDCAEdges[1]);
  TH1D *hPDYAbsIn_3_10_Deg = new TH1D("hPDYAbsIn_3_10_Deg","p #times DY at absorber begin for tracks within [3,10[ degrees at absorber end;p #times DY (GeV #times cm)",2*pDCANBins, -pDCAEdges[1], pDCAEdges[1]);
  TH2D *hPDYAbsInVsP_2_3_Deg = new TH2D("hPDYAbsInVsP_2_3_Deg","p #times DY at absorber begin versus p for tracks within [2,3[ degrees at absorber end;p (GeV/c);p #times DY (GeV #times cm)",pNBins,pEdges[0],pEdges[1], 2*pDCANBins, -pDCAEdges[1], pDCAEdges[1]);
  TH2D *hPDYAbsInVsP_3_10_Deg = new TH2D("hPDYAbsInVsP_3_10_Deg","p #times DY at absorber begin versus p for tracks within [3,10[ degrees at absorber end;p (GeV/c);p #times DY (GeV #times cm)",pNBins,pEdges[0],pEdges[1], 2*pDCANBins, -pDCAEdges[1], pDCAEdges[1]);
  TGraphAsymmErrors* gMeanPDYAbsInVsP_2_3_Deg = new TGraphAsymmErrors(pNBins);
  gMeanPDYAbsInVsP_2_3_Deg->SetName("gMeanPDYAbsInVsP_2_3_Deg");
  gMeanPDYAbsInVsP_2_3_Deg->SetTitle("<p #times DY at absorber begin> versus p for tracks within [2,3[ degrees at absorber end;p (GeV/c);<p #times DY> (GeV #times cm)");
  TGraphAsymmErrors* gSigmaPDYAbsInVsP_2_3_Deg = new TGraphAsymmErrors(pNBins);
  gSigmaPDYAbsInVsP_2_3_Deg->SetName("gSigmaPDYAbsInVsP_2_3_Deg");
  gSigmaPDYAbsInVsP_2_3_Deg->SetTitle("#sigma_{p #times DY at absorber begin} versus p for tracks within [2,3[ degrees at absorber end;p (GeV/c);#sigma_{p #times DY} (GeV #times cm)");
  TGraphAsymmErrors* gMeanPDYAbsInVsP_3_10_Deg = new TGraphAsymmErrors(pNBins);
  gMeanPDYAbsInVsP_3_10_Deg->SetName("gMeanPDYAbsInVsP_3_10_Deg");
  gMeanPDYAbsInVsP_3_10_Deg->SetTitle("<p #times DY at absorber begin> versus p for tracks within [3,10[ degrees at absorber end;p (GeV/c);<p #times DY> (GeV #times cm)");
  TGraphAsymmErrors* gSigmaPDYAbsInVsP_3_10_Deg = new TGraphAsymmErrors(pNBins);
  gSigmaPDYAbsInVsP_3_10_Deg->SetName("gSigmaPDYAbsInVsP_3_10_Deg");
  gSigmaPDYAbsInVsP_3_10_Deg->SetTitle("#sigma_{p #times DY at absorber begin} versus p for tracks within [3,10[ degrees at absorber end;p (GeV/c);#sigma_{p #times DY} (GeV #times cm)");
  TH2D *hPResDYAbsInVsP_2_3_Deg = new TH2D("hPResDYAbsInVsP_2_3_Deg","p #times #sigma_{DY} at absorber begin versus p for tracks within [2,3[ degrees at absorber end;p (GeV/c);p #times #sigma_{DY} (GeV #times cm)",pNBins,pEdges[0],pEdges[1], pDCANBins, pResDCAEdges[0], pResDCAEdges[1]);
  TH2D *hPResDYAbsInVsP_3_10_Deg = new TH2D("hPResDYAbsInVsP_3_10_Deg","p #times #sigma_{DY} at absorber begin versus p for tracks within [3,10[ degrees at absorber end;p (GeV/c);p #times #sigma_{DY} (GeV #times cm)",pNBins,pEdges[0],pEdges[1], pDCANBins, pResDCAEdges[0], pResDCAEdges[1]);
  
  TH1D *hPDAbsIn_2_3_Deg = new TH1D("hPDAbsIn_2_3_Deg","p #times D at absorber begin for tracks within [2,3[ degrees at absorber end;p #times D (GeV #times cm)",pDCANBins, pDCAEdges[0], pDCAEdges[1]);
  TH1D *hPDAbsIn_3_10_Deg = new TH1D("hPDAbsIn_3_10_Deg","p #times D at absorber begin for tracks within [3,10[ degrees at absorber end;p #times D (GeV #times cm)",pDCANBins, pDCAEdges[0], pDCAEdges[1]);
  TH2D *hPDAbsInVsP_2_3_Deg = new TH2D("hPDAbsInVsP_2_3_Deg","p #times D at absorber begin versus p for tracks within [2,3[ degrees at absorber end;p (GeV/c);p #times D (GeV #times cm)",pNBins,pEdges[0],pEdges[1], pDCANBins, pDCAEdges[0], pDCAEdges[1]);
  TH2D *hPDAbsInVsP_3_10_Deg = new TH2D("hPDAbsInVsP_3_10_Deg","p #times D at absorber begin versus p for tracks within [3,10[ degrees at absorber end;p (GeV/c);p #times D (GeV #times cm)",pNBins,pEdges[0],pEdges[1], pDCANBins, pDCAEdges[0], pDCAEdges[1]);
  TGraphAsymmErrors* gMeanPDAbsInVsP_2_3_Deg = new TGraphAsymmErrors(pNBins);
  gMeanPDAbsInVsP_2_3_Deg->SetName("gMeanPDAbsInVsP_2_3_Deg");
  gMeanPDAbsInVsP_2_3_Deg->SetTitle("<p #times D at absorber begin> versus p for tracks within [2,3[ degrees at absorber end;p (GeV/c);<p #times D> (GeV #times cm)");
  TGraphAsymmErrors* gSigmaPDAbsInVsP_2_3_Deg = new TGraphAsymmErrors(pNBins);
  gSigmaPDAbsInVsP_2_3_Deg->SetName("gSigmaPDAbsInVsP_2_3_Deg");
  gSigmaPDAbsInVsP_2_3_Deg->SetTitle("#sigma_{p #times D at absorber begin} versus p for tracks within [2,3[ degrees at absorber end;p (GeV/c);#sigma_{p #times D} (GeV #times cm)");
  TGraphAsymmErrors* gMeanPDAbsInVsP_3_10_Deg = new TGraphAsymmErrors(pNBins);
  gMeanPDAbsInVsP_3_10_Deg->SetName("gMeanPDAbsInVsP_3_10_Deg");
  gMeanPDAbsInVsP_3_10_Deg->SetTitle("<p #times D at absorber begin> versus p for tracks within [3,10[ degrees at absorber end;p (GeV/c);<p #times D> (GeV #times cm)");
  TGraphAsymmErrors* gSigmaPDAbsInVsP_3_10_Deg = new TGraphAsymmErrors(pNBins);
  gSigmaPDAbsInVsP_3_10_Deg->SetName("gSigmaPDAbsInVsP_3_10_Deg");
  gSigmaPDAbsInVsP_3_10_Deg->SetTitle("#sigma_{p #times D at absorber begin} versus p for tracks within [3,10[ degrees at absorber end;p (GeV/c);#sigma_{p #times D} (GeV #times cm)");
  TH2D *hPResDAbsInVsP_2_3_Deg = new TH2D("hPResDAbsInVsP_2_3_Deg","p #times #sigma_{D} at absorber begin versus p for tracks within [2,3[ degrees at absorber end;p (GeV/c);p #times #sigma_{D} (GeV #times cm)",pNBins,pEdges[0],pEdges[1], pDCANBins, pResDCAEdges[0], pResDCAEdges[1]);
  TH2D *hPResDAbsInVsP_3_10_Deg = new TH2D("hPResDAbsInVsP_3_10_Deg","p #times #sigma_{D} at absorber begin versus p for tracks within [3,10[ degrees at absorber end;p (GeV/c);p #times #sigma_{D} (GeV #times cm)",pNBins,pEdges[0],pEdges[1], pDCANBins, pResDCAEdges[0], pResDCAEdges[1]);
  
  // at the begining of the absorber with branson correction
  histoFile->mkdir("pDXYAtAbsInWithB","pDXYAtAbsInWithB");
  histoFile->cd("pDXYAtAbsInWithB");
  
  const Double_t pDEdges[2] = {0., 50.};
  const Double_t pResDEdges[2] = {0., 10.};
  
  TH1D *hPDXAbsInWithB_2_3_Deg = new TH1D("hPDXAbsInWithB_2_3_Deg","p #times DX at absorber begin with Branson correction for tracks within [2,3[ degrees at absorber end;p #times DX (GeV #times cm)",2*pDCANBins, -pDEdges[1], pDEdges[1]);
  TH1D *hPDXAbsInWithB_3_10_Deg = new TH1D("hPDXAbsInWithB_3_10_Deg","p #times DX at absorber begin with Branson correction for tracks within [3,10[ degrees at absorber end;p #times DX (GeV #times cm)",2*pDCANBins, -pDEdges[1], pDEdges[1]);
  TH2D *hPDXAbsInWithBVsP_2_3_Deg = new TH2D("hPDXAbsInWithBVsP_2_3_Deg","p #times DX at absorber begin with Branson correction versus p for tracks within [2,3[ degrees at absorber end;p (GeV/c);p #times DX (GeV #times cm)",pNBins,pEdges[0],pEdges[1], 2*pDCANBins, -pDEdges[1], pDEdges[1]);
  TH2D *hPDXAbsInWithBVsP_3_10_Deg = new TH2D("hPDXAbsInWithBVsP_3_10_Deg","p #times DX at absorber begin with Branson correction versus p for tracks within [3,10[ degrees at absorber end;p (GeV/c);p #times DX (GeV #times cm)",pNBins,pEdges[0],pEdges[1], 2*pDCANBins, -pDEdges[1], pDEdges[1]);
  TGraphAsymmErrors* gMeanPDXAbsInWithBVsP_2_3_Deg = new TGraphAsymmErrors(pNBins);
  gMeanPDXAbsInWithBVsP_2_3_Deg->SetName("gMeanPDXAbsInWithBVsP_2_3_Deg");
  gMeanPDXAbsInWithBVsP_2_3_Deg->SetTitle("<p #times DX at absorber begin with Branson correction> versus p for tracks within [2,3[ degrees at absorber end;p (GeV/c);<p #times DX> (GeV #times cm)");
  TGraphAsymmErrors* gSigmaPDXAbsInWithBVsP_2_3_Deg = new TGraphAsymmErrors(pNBins);
  gSigmaPDXAbsInWithBVsP_2_3_Deg->SetName("gSigmaPDXAbsInWithBVsP_2_3_Deg");
  gSigmaPDXAbsInWithBVsP_2_3_Deg->SetTitle("#sigma_{p #times DX at absorber begin with Branson correction} versus p for tracks within [2,3[ degrees at absorber end;p (GeV/c);#sigma_{p #times DX} (GeV #times cm)");
  TGraphAsymmErrors* gMeanPDXAbsInWithBVsP_3_10_Deg = new TGraphAsymmErrors(pNBins);
  gMeanPDXAbsInWithBVsP_3_10_Deg->SetName("gMeanPDXAbsInWithBVsP_3_10_Deg");
  gMeanPDXAbsInWithBVsP_3_10_Deg->SetTitle("<p #times DX at absorber begin with Branson correction> versus p for tracks within [3,10[ degrees at absorber end;p (GeV/c);<p #times DX> (GeV #times cm)");
  TGraphAsymmErrors* gSigmaPDXAbsInWithBVsP_3_10_Deg = new TGraphAsymmErrors(pNBins);
  gSigmaPDXAbsInWithBVsP_3_10_Deg->SetName("gSigmaPDXAbsInWithBVsP_3_10_Deg");
  gSigmaPDXAbsInWithBVsP_3_10_Deg->SetTitle("#sigma_{p #times DX at absorber begin with Branson correction} versus p for tracks within [3,10[ degrees at absorber end;p (GeV/c);#sigma_{p #times DX} (GeV #times cm)");
  TH2D *hPResDXAbsInWithBVsP_2_3_Deg = new TH2D("hPResDXAbsInWithBVsP_2_3_Deg","p #times #sigma_{DX} at absorber begin with Branson correction versus p for tracks within [2,3[ degrees at absorber end;p (GeV/c);p #times #sigma_{DX} (GeV #times cm)",pNBins,pEdges[0],pEdges[1], pDCANBins, pResDEdges[0], pResDEdges[1]);
  TH2D *hPResDXAbsInWithBVsP_3_10_Deg = new TH2D("hPResDXAbsInWithBVsP_3_10_Deg","p #times #sigma_{DX} at absorber begin with Branson correction versus p for tracks within [3,10[ degrees at absorber end;p (GeV/c);p #times #sigma_{DX} (GeV #times cm)",pNBins,pEdges[0],pEdges[1], pDCANBins, pResDEdges[0], pResDEdges[1]);
  
  TH1D *hPDYAbsInWithB_2_3_Deg = new TH1D("hPDYAbsInWithB_2_3_Deg","p #times DY at absorber begin with Branson correction for tracks within [2,3[ degrees at absorber end;p #times DY (GeV #times cm)",2*pDCANBins, -pDEdges[1], pDEdges[1]);
  TH1D *hPDYAbsInWithB_3_10_Deg = new TH1D("hPDYAbsInWithB_3_10_Deg","p #times DY at absorber begin with Branson correction for tracks within [3,10[ degrees at absorber end;p #times DY (GeV #times cm)",2*pDCANBins, -pDEdges[1], pDEdges[1]);
  TH2D *hPDYAbsInWithBVsP_2_3_Deg = new TH2D("hPDYAbsInWithBVsP_2_3_Deg","p #times DY at absorber begin with Branson correction versus p for tracks within [2,3[ degrees at absorber end;p (GeV/c);p #times DY (GeV #times cm)",pNBins,pEdges[0],pEdges[1], 2*pDCANBins, -pDEdges[1], pDEdges[1]);
  TH2D *hPDYAbsInWithBVsP_3_10_Deg = new TH2D("hPDYAbsInWithBVsP_3_10_Deg","p #times DY at absorber begin with Branson correction versus p for tracks within [3,10[ degrees at absorber end;p (GeV/c);p #times DY (GeV #times cm)",pNBins,pEdges[0],pEdges[1], 2*pDCANBins, -pDEdges[1], pDEdges[1]);
  TGraphAsymmErrors* gMeanPDYAbsInWithBVsP_2_3_Deg = new TGraphAsymmErrors(pNBins);
  gMeanPDYAbsInWithBVsP_2_3_Deg->SetName("gMeanPDYAbsInWithBVsP_2_3_Deg");
  gMeanPDYAbsInWithBVsP_2_3_Deg->SetTitle("<p #times DY at absorber begin with Branson correction> versus p for tracks within [2,3[ degrees at absorber end;p (GeV/c);<p #times DY> (GeV #times cm)");
  TGraphAsymmErrors* gSigmaPDYAbsInWithBVsP_2_3_Deg = new TGraphAsymmErrors(pNBins);
  gSigmaPDYAbsInWithBVsP_2_3_Deg->SetName("gSigmaPDYAbsInWithBVsP_2_3_Deg");
  gSigmaPDYAbsInWithBVsP_2_3_Deg->SetTitle("#sigma_{p #times DY at absorber begin with Branson correction} versus p for tracks within [2,3[ degrees at absorber end;p (GeV/c);#sigma_{p #times DY} (GeV #times cm)");
  TGraphAsymmErrors* gMeanPDYAbsInWithBVsP_3_10_Deg = new TGraphAsymmErrors(pNBins);
  gMeanPDYAbsInWithBVsP_3_10_Deg->SetName("gMeanPDYAbsInWithBVsP_3_10_Deg");
  gMeanPDYAbsInWithBVsP_3_10_Deg->SetTitle("<p #times DY at absorber begin with Branson correction> versus p for tracks within [3,10[ degrees at absorber end;p (GeV/c);<p #times DY> (GeV #times cm)");
  TGraphAsymmErrors* gSigmaPDYAbsInWithBVsP_3_10_Deg = new TGraphAsymmErrors(pNBins);
  gSigmaPDYAbsInWithBVsP_3_10_Deg->SetName("gSigmaPDYAbsInWithBVsP_3_10_Deg");
  gSigmaPDYAbsInWithBVsP_3_10_Deg->SetTitle("#sigma_{p #times DY at absorber begin with Branson correction} versus p for tracks within [3,10[ degrees at absorber end;p (GeV/c);#sigma_{p #times DY} (GeV #times cm)");
  TH2D *hPResDYAbsInWithBVsP_2_3_Deg = new TH2D("hPResDYAbsInWithBVsP_2_3_Deg","p #times #sigma_{DY} at absorber begin with Branson correction versus p for tracks within [2,3[ degrees at absorber end;p (GeV/c);p #times #sigma_{DY} (GeV #times cm)",pNBins,pEdges[0],pEdges[1], pDCANBins, pResDEdges[0], pResDEdges[1]);
  TH2D *hPResDYAbsInWithBVsP_3_10_Deg = new TH2D("hPResDYAbsInWithBVsP_3_10_Deg","p #times #sigma_{DY} at absorber begin with Branson correction versus p for tracks within [3,10[ degrees at absorber end;p (GeV/c);p #times #sigma_{DY} (GeV #times cm)",pNBins,pEdges[0],pEdges[1], pDCANBins, pResDEdges[0], pResDEdges[1]);
  
  TH1D *hPDAbsInWithB_2_3_Deg = new TH1D("hPDAbsInWithB_2_3_Deg","p #times D at absorber begin with Branson correction for tracks within [2,3[ degrees at absorber end;p #times D (GeV #times cm)",pDCANBins, pDEdges[0], pDEdges[1]);
  TH1D *hPDAbsInWithB_3_10_Deg = new TH1D("hPDAbsInWithB_3_10_Deg","p #times D at absorber begin with Branson correction for tracks within [3,10[ degrees at absorber end;p #times D (GeV #times cm)",pDCANBins, pDEdges[0], pDEdges[1]);
  TH2D *hPDAbsInWithBVsP_2_3_Deg = new TH2D("hPDAbsInWithBVsP_2_3_Deg","p #times D at absorber begin with Branson correction versus p for tracks within [2,3[ degrees at absorber end;p (GeV/c);p #times D (GeV #times cm)",pNBins,pEdges[0],pEdges[1], pDCANBins, pDEdges[0], pDEdges[1]);
  TH2D *hPDAbsInWithBVsP_3_10_Deg = new TH2D("hPDAbsInWithBVsP_3_10_Deg","p #times D at absorber begin with Branson correction versus p for tracks within [3,10[ degrees at absorber end;p (GeV/c);p #times D (GeV #times cm)",pNBins,pEdges[0],pEdges[1], pDCANBins, pDEdges[0], pDEdges[1]);
  TGraphAsymmErrors* gMeanPDAbsInWithBVsP_2_3_Deg = new TGraphAsymmErrors(pNBins);
  gMeanPDAbsInWithBVsP_2_3_Deg->SetName("gMeanPDAbsInWithBVsP_2_3_Deg");
  gMeanPDAbsInWithBVsP_2_3_Deg->SetTitle("<p #times D at absorber begin with Branson correction> versus p for tracks within [2,3[ degrees at absorber end;p (GeV/c);<p #times D> (GeV #times cm)");
  TGraphAsymmErrors* gSigmaPDAbsInWithBVsP_2_3_Deg = new TGraphAsymmErrors(pNBins);
  gSigmaPDAbsInWithBVsP_2_3_Deg->SetName("gSigmaPDAbsInWithBVsP_2_3_Deg");
  gSigmaPDAbsInWithBVsP_2_3_Deg->SetTitle("#sigma_{p #times D at absorber begin with Branson correction} versus p for tracks within [2,3[ degrees at absorber end;p (GeV/c);#sigma_{p #times D} (GeV #times cm)");
  TGraphAsymmErrors* gMeanPDAbsInWithBVsP_3_10_Deg = new TGraphAsymmErrors(pNBins);
  gMeanPDAbsInWithBVsP_3_10_Deg->SetName("gMeanPDAbsInWithBVsP_3_10_Deg");
  gMeanPDAbsInWithBVsP_3_10_Deg->SetTitle("<p #times D at absorber begin with Branson correction> versus p for tracks within [3,10[ degrees at absorber end;p (GeV/c);<p #times D> (GeV #times cm)");
  TGraphAsymmErrors* gSigmaPDAbsInWithBVsP_3_10_Deg = new TGraphAsymmErrors(pNBins);
  gSigmaPDAbsInWithBVsP_3_10_Deg->SetName("gSigmaPDAbsInWithBVsP_3_10_Deg");
  gSigmaPDAbsInWithBVsP_3_10_Deg->SetTitle("#sigma_{p #times D at absorber begin with Branson correction} versus p for tracks within [3,10[ degrees at absorber end;p (GeV/c);#sigma_{p #times D} (GeV #times cm)");
  TH2D *hPResDAbsInWithBVsP_2_3_Deg = new TH2D("hPResDAbsInWithBVsP_2_3_Deg","p #times #sigma_{D} at absorber begin with Branson correction versus p for tracks within [2,3[ degrees at absorber end;p (GeV/c);p #times #sigma_{D} (GeV #times cm)",pNBins,pEdges[0],pEdges[1], pDCANBins, pResDEdges[0], pResDEdges[1]);
  TH2D *hPResDAbsInWithBVsP_3_10_Deg = new TH2D("hPResDAbsInWithBVsP_3_10_Deg","p #times #sigma_{D} at absorber begin with Branson correction versus p for tracks within [3,10[ degrees at absorber end;p (GeV/c);p #times #sigma_{D} (GeV #times cm)",pNBins,pEdges[0],pEdges[1], pDCANBins, pResDEdges[0], pResDEdges[1]);
  
  // loop over events
  Int_t nevents = rc.NumberOfEvents();
  for (Int_t ievent = 0; ievent < nevents; ievent++)
  {
    if ((ievent+1)%100 == 0) cout<<"\rEvent processing... "<<ievent+1<<flush;
    
    AliMUONVTrackStore* trackStore = rc.ReconstructedTracks(ievent, kFALSE);
    AliMUONVTrackStore* trackRefStore = rc.ReconstructibleTracks(ievent, requestedStationMask, request2ChInSameSt45);
    
    // loop over trackRef
    TIter next(trackRefStore->CreateIterator());
    AliMUONTrack* trackRef;
    while ( ( trackRef = static_cast<AliMUONTrack*>(next()) ) )
    {
      
      // loop over trackReco and look for compatible track
      AliMUONTrack* trackMatched = 0x0;
      Int_t nMatchClusters = 0;
      TIter next2(trackStore->CreateIterator());
      AliMUONTrack* trackReco;
      while ( ( trackReco = static_cast<AliMUONTrack*>(next2()) ) ) {
	if (trackReco->Match(*trackRef, sigmaCut, nMatchClusters)) {
	  trackMatched = trackReco;
	  break;
	}
      }
      if (!trackMatched) continue;
      
      // get momentum at vertex
      AliMUONTrackParam *trackParam = trackMatched->GetTrackParamAtVertex();
      Double_t p = trackParam->P();
      
      // get track position at the end of the absorber
      AliMUONTrackParam trackParamAtAbsEnd(*((AliMUONTrackParam*)trackMatched->GetTrackParamAtCluster()->First()));
      AliMUONTrackExtrap::ExtrapToZ(&trackParamAtAbsEnd, AliMUONConstants::AbsZEnd());
      Double_t pU = trackParamAtAbsEnd.P();
      Double_t xAbs = trackParamAtAbsEnd.GetNonBendingCoor();
      Double_t yAbs = trackParamAtAbsEnd.GetBendingCoor();
      Double_t dAbs = TMath::Sqrt(xAbs*xAbs + yAbs*yAbs);
      Double_t aAbs = TMath::ATan(-dAbs/AliMUONConstants::AbsZEnd()) * TMath::RadToDeg();
      
      // get track parameters at DCA
      AliMUONTrackParam trackParamAtDCA(*((AliMUONTrackParam*) trackMatched->GetTrackParamAtCluster()->First()));
      AliMUONTrackExtrap::ExtrapToVertexWithoutBranson(&trackParamAtDCA, 0.);
      Double_t pMean = 0.5*(p+pU);
      Double_t xDCA = trackParamAtDCA.GetNonBendingCoor();
      Double_t yDCA = trackParamAtDCA.GetBendingCoor();
      Double_t dca = TMath::Sqrt(xDCA*xDCA + yDCA*yDCA);
      Double_t pDCAX = pMean*xDCA;
      Double_t pDCAY = pMean*yDCA;
      Double_t pDCA = pMean*dca;
      const TMatrixD &covDCA = trackParamAtDCA.GetCovariances();
      Double_t pSigmaDCAX = pMean*TMath::Sqrt(covDCA(0,0));
      Double_t pSigmaDCAY = pMean*TMath::Sqrt(covDCA(2,2));
      Double_t sigmaDCA = TMath::Sqrt(xDCA*xDCA*covDCA(0,0) + yDCA*yDCA*covDCA(2,2)) / dca;
      Double_t pSigmaDCA = pMean*sigmaDCA;
      
      // get trackRef parameters at the begining of the absorber
      AliMUONTrackParam trackRefParamAtAbsIn(*(trackRef->GetTrackParamAtVertex()));
      AliMUONTrackExtrap::ExtrapToZ(&trackRefParamAtAbsIn, zAbsIn);
      Double_t xRefAbsIn = trackRefParamAtAbsIn.GetNonBendingCoor();
      Double_t yRefAbsIn = trackRefParamAtAbsIn.GetBendingCoor();
      
      // get track parameters at the begining of the absorber without branson correction
      AliMUONTrackParam trackParamAtAbsIn(trackParamAtDCA);
      AliMUONTrackExtrap::ExtrapToZCov(&trackParamAtAbsIn, zAbsIn);
      Double_t dxAbsIn = trackParamAtAbsIn.GetNonBendingCoor() - xRefAbsIn;
      Double_t dyAbsIn = trackParamAtAbsIn.GetBendingCoor() - yRefAbsIn;
      Double_t dAbsIn = TMath::Sqrt(dxAbsIn*dxAbsIn + dyAbsIn*dyAbsIn);
      Double_t pDAbsInX = pMean*dxAbsIn;
      Double_t pDAbsInY = pMean*dyAbsIn;
      Double_t pDAbsIn = pMean*dAbsIn;
      const TMatrixD &covDAbsIn = trackParamAtAbsIn.GetCovariances();
      Double_t pSigmaAbsInX = pMean*TMath::Sqrt(covDAbsIn(0,0));
      Double_t pSigmaAbsInY = pMean*TMath::Sqrt(covDAbsIn(2,2));
      Double_t sigmaAbsIn = TMath::Sqrt(dxAbsIn*dxAbsIn*covDAbsIn(0,0) + dyAbsIn*dyAbsIn*covDAbsIn(2,2)) / dAbsIn;
      Double_t pSigmaAbsIn = pMean*sigmaAbsIn;
      
      // get track parameters at the begining of the absorber with branson correction
      AliMUONTrackParam trackParamAtAbsInWithB(*((AliMUONTrackParam*) trackMatched->GetTrackParamAtCluster()->First()));
      AliMUONTrackExtrap::ExtrapToVertex(&trackParamAtAbsInWithB, 0., 0., 0., 0., 0.);
      AliMUONTrackExtrap::ExtrapToZCov(&trackParamAtAbsInWithB, zAbsIn);
      Double_t dxAbsInWithB = trackParamAtAbsInWithB.GetNonBendingCoor() - xRefAbsIn;
      Double_t dyAbsInWithB = trackParamAtAbsInWithB.GetBendingCoor() - yRefAbsIn;
      Double_t dAbsInWithB = TMath::Sqrt(dxAbsInWithB*dxAbsInWithB + dyAbsInWithB*dyAbsInWithB);
      Double_t pDAbsInWithBX = pMean*dxAbsInWithB;
      Double_t pDAbsInWithBY = pMean*dyAbsInWithB;
      Double_t pDAbsInWithB = pMean*dAbsInWithB;
      const TMatrixD &covDAbsInWithB = trackParamAtAbsInWithB.GetCovariances();
      Double_t pSigmaAbsInWithBX = pMean*TMath::Sqrt(covDAbsInWithB(0,0));
      Double_t pSigmaAbsInWithBY = pMean*TMath::Sqrt(covDAbsInWithB(2,2));
      Double_t sigmaAbsInWithB = TMath::Sqrt(dxAbsInWithB*dxAbsInWithB*covDAbsInWithB(0,0) + dyAbsInWithB*dyAbsInWithB*covDAbsInWithB(2,2)) / dAbsInWithB;
      Double_t pSigmaAbsInWithB = pMean*sigmaAbsInWithB;
      
      // plot pDCA vs p
      if (aAbs > 2. && aAbs < 3.) {
	
	hPDCAX_2_3_Deg->Fill(pDCAX);
	hPDCAY_2_3_Deg->Fill(pDCAY);
	hPDCA_2_3_Deg->Fill(pDCA);
	hPDCAXVsP_2_3_Deg->Fill(pMean, pDCAX);
	hPDCAYVsP_2_3_Deg->Fill(pMean, pDCAY);
	hPDCAVsP_2_3_Deg->Fill(pMean, pDCA);
	hPResDCAXVsP_2_3_Deg->Fill(pMean, pSigmaDCAX);
	hPResDCAYVsP_2_3_Deg->Fill(pMean, pSigmaDCAY);
	hPResDCAVsP_2_3_Deg->Fill(pMean, pSigmaDCA);
	
	hPDXAbsIn_2_3_Deg->Fill(pDAbsInX);
	hPDYAbsIn_2_3_Deg->Fill(pDAbsInY);
	hPDAbsIn_2_3_Deg->Fill(pDAbsIn);
	hPDXAbsInVsP_2_3_Deg->Fill(pMean, pDAbsInX);
	hPDYAbsInVsP_2_3_Deg->Fill(pMean, pDAbsInY);
	hPDAbsInVsP_2_3_Deg->Fill(pMean, pDAbsIn);
	hPResDXAbsInVsP_2_3_Deg->Fill(pMean, pSigmaAbsInX);
	hPResDYAbsInVsP_2_3_Deg->Fill(pMean, pSigmaAbsInY);
	hPResDAbsInVsP_2_3_Deg->Fill(pMean, pSigmaAbsIn);
	
	hPDXAbsInWithB_2_3_Deg->Fill(pDAbsInWithBX);
	hPDYAbsInWithB_2_3_Deg->Fill(pDAbsInWithBY);
	hPDAbsInWithB_2_3_Deg->Fill(pDAbsInWithB);
	hPDXAbsInWithBVsP_2_3_Deg->Fill(pMean, pDAbsInWithBX);
	hPDYAbsInWithBVsP_2_3_Deg->Fill(pMean, pDAbsInWithBY);
	hPDAbsInWithBVsP_2_3_Deg->Fill(pMean, pDAbsInWithB);
	hPResDXAbsInWithBVsP_2_3_Deg->Fill(pMean, pSigmaAbsInWithBX);
	hPResDYAbsInWithBVsP_2_3_Deg->Fill(pMean, pSigmaAbsInWithBY);
	hPResDAbsInWithBVsP_2_3_Deg->Fill(pMean, pSigmaAbsInWithB);
	
      }
      else if (aAbs >= 3. && aAbs < 10.) {
	
	hPDCAX_3_10_Deg->Fill(pDCAX);
	hPDCAY_3_10_Deg->Fill(pDCAY);
	hPDCA_3_10_Deg->Fill(pDCA);
	hPDCAXVsP_3_10_Deg->Fill(pMean, pDCAX);
	hPDCAYVsP_3_10_Deg->Fill(pMean, pDCAY);
	hPDCAVsP_3_10_Deg->Fill(pMean, pDCA);
	hPResDCAXVsP_3_10_Deg->Fill(pMean, pSigmaDCAX);
	hPResDCAYVsP_3_10_Deg->Fill(pMean, pSigmaDCAY);
	hPResDCAVsP_3_10_Deg->Fill(pMean, pSigmaDCA);
	
	hPDXAbsIn_3_10_Deg->Fill(pDAbsInX);
	hPDYAbsIn_3_10_Deg->Fill(pDAbsInY);
	hPDAbsIn_3_10_Deg->Fill(pDAbsIn);
	hPDXAbsInVsP_3_10_Deg->Fill(pMean, pDAbsInX);
	hPDYAbsInVsP_3_10_Deg->Fill(pMean, pDAbsInY);
	hPDAbsInVsP_3_10_Deg->Fill(pMean, pDAbsIn);
	hPResDXAbsInVsP_3_10_Deg->Fill(pMean, pSigmaAbsInX);
	hPResDYAbsInVsP_3_10_Deg->Fill(pMean, pSigmaAbsInY);
	hPResDAbsInVsP_3_10_Deg->Fill(pMean, pSigmaAbsIn);
	
	hPDXAbsInWithB_3_10_Deg->Fill(pDAbsInWithBX);
	hPDYAbsInWithB_3_10_Deg->Fill(pDAbsInWithBY);
	hPDAbsInWithB_3_10_Deg->Fill(pDAbsInWithB);
	hPDXAbsInWithBVsP_3_10_Deg->Fill(pMean, pDAbsInWithBX);
	hPDYAbsInWithBVsP_3_10_Deg->Fill(pMean, pDAbsInWithBY);
	hPDAbsInWithBVsP_3_10_Deg->Fill(pMean, pDAbsInWithB);
	hPResDXAbsInWithBVsP_3_10_Deg->Fill(pMean, pSigmaAbsInWithBX);
	hPResDYAbsInWithBVsP_3_10_Deg->Fill(pMean, pSigmaAbsInWithBY);
	hPResDAbsInWithBVsP_3_10_Deg->Fill(pMean, pSigmaAbsInWithB);
	
      }
      
    } // end loop track ref.
    
    TLorentzVector v1, v2, v;
    AliESDEvent *esd = const_cast<AliESDEvent*>(rc.GetESDEvent());
    Int_t nTracks = (Int_t)esd->GetNumberOfMuonTracks(); 
    
    // first loop over reconstructed tracks
    for (Int_t iTrack1 = 0; iTrack1 < nTracks; iTrack1++) {
      
      AliESDMuonTrack* muonTrack1 = esd->GetMuonTrack(iTrack1);
      
      // skip ghosts
      if (!muonTrack1->ContainTrackerData()) continue;
      
      Short_t charge1 = muonTrack1->Charge();
      muonTrack1->LorentzP(v1);
      
      // second loop over reconstructed tracks
      for (Int_t iTrack2 = iTrack1 + 1; iTrack2 < nTracks; iTrack2++) {
	
	AliESDMuonTrack* muonTrack2 = esd->GetMuonTrack(iTrack2);
	
	// skip ghosts
	if (!muonTrack2->ContainTrackerData()) continue;
	
	// skip like sign pairs
	if (charge1 * muonTrack2->Charge() > 0) continue;
	
	muonTrack2->LorentzP(v2);
	
	// invariant mass
	v = v1 + v2;
	hMass->Fill(v.M());
	
      }
      
    }
    
  } // end loop on event  
  cout<<"\rEvent processing... "<<nevents<<" done"<<endl;
  
  TF1 *fCB = new TF1("fCB", CrystalBallExtended, 0, 6, 7);
  fGaus = new TF1("fGaus","gaus");
  fPGaus = new TF1("fPGaus","x*gaus");
  
  // compute p*DCA(XY) resolution in the region [2,3] deg at absorber end
  fGaus->SetParameters(1., 0., 90.);
  hPDCAX_2_3_Deg->Fit("fGaus");
  Double_t sDCAX23 = fGaus->GetParameter(2);
  fGaus->SetParameters(1., 0., 90.);
  hPDCAY_2_3_Deg->Fit("fGaus");
  Double_t sDCAY23 = fGaus->GetParameter(2);
  fPGaus->SetParameters(1., 0., 90.);
  hPDCA_2_3_Deg->Fit("fPGaus");
  FitGausResVsP(hPDCAXVsP_2_3_Deg, pNBins, 0., 90., "p*DCAX (tracks in [2,3] deg.)", gMeanPDCAXVsP_2_3_Deg, gSigmaPDCAXVsP_2_3_Deg);
  FitGausResVsP(hPDCAYVsP_2_3_Deg, pNBins, 0., 90., "p*DCAY (tracks in [2,3] deg.)", gMeanPDCAYVsP_2_3_Deg, gSigmaPDCAYVsP_2_3_Deg);
  FitPGausResVsP(hPDCAVsP_2_3_Deg, pNBins, 0., 90., "p*DCA (tracks in [2,3] deg.)", gMeanPDCAVsP_2_3_Deg, gSigmaPDCAVsP_2_3_Deg);
  
  // compute p*DCA(XY) resolution in the region [3,10] deg at absorber end
  fGaus->SetParameters(1., 0., 55.);
  hPDCAX_3_10_Deg->Fit("fGaus");
  Double_t sDCAX310 = fGaus->GetParameter(2);
  fGaus->SetParameters(1., 0., 55.);
  hPDCAY_3_10_Deg->Fit("fGaus");
  Double_t sDCAY310 = fGaus->GetParameter(2);
  fPGaus->SetParameters(1., 0., 55.);
  hPDCA_3_10_Deg->Fit("fPGaus");
  FitGausResVsP(hPDCAXVsP_3_10_Deg, pNBins, 0., 55., "p*DCAX (tracks in [3,10] deg.)", gMeanPDCAXVsP_3_10_Deg, gSigmaPDCAXVsP_3_10_Deg);
  FitGausResVsP(hPDCAYVsP_3_10_Deg, pNBins, 0., 55., "p*DCAY (tracks in [3,10] deg.)", gMeanPDCAYVsP_3_10_Deg, gSigmaPDCAYVsP_3_10_Deg);
  FitPGausResVsP(hPDCAVsP_3_10_Deg, pNBins, 0., 55., "p*DCA (tracks in [3,10] deg.)", gMeanPDCAVsP_3_10_Deg, gSigmaPDCAVsP_3_10_Deg);
  
  printf("sigma pDCAXY ~ %f GeV.cm (2-3°) / %f GeV.cm (3-10°)\n", 0.5*(sDCAX23+sDCAY23), 0.5*(sDCAX310+sDCAY310));
  
  // compute p*D(XY) resolution in the region [2,3] deg at absorber end without Branson correction
  fGaus->SetParameters(1., 0., 75.);
  hPDXAbsIn_2_3_Deg->Fit("fGaus");
  Double_t sDX23 = fGaus->GetParameter(2);
  fGaus->SetParameters(1., 0., 75.);
  hPDYAbsIn_2_3_Deg->Fit("fGaus");
  Double_t sDY23 = fGaus->GetParameter(2);
  fPGaus->SetParameters(1., 0., 75.);
  hPDAbsIn_2_3_Deg->Fit("fPGaus");
  FitGausResVsP(hPDXAbsInVsP_2_3_Deg, pNBins, 0., 75., "p*DX (tracks in [2,3] deg.)", gMeanPDXAbsInVsP_2_3_Deg, gSigmaPDXAbsInVsP_2_3_Deg);
  FitGausResVsP(hPDYAbsInVsP_2_3_Deg, pNBins, 0., 75., "p*DY (tracks in [2,3] deg.)", gMeanPDYAbsInVsP_2_3_Deg, gSigmaPDYAbsInVsP_2_3_Deg);
  FitPGausResVsP(hPDAbsInVsP_2_3_Deg, pNBins, 0., 75., "p*D (tracks in [2,3] deg.)", gMeanPDAbsInVsP_2_3_Deg, gSigmaPDAbsInVsP_2_3_Deg);
  
  // compute p*D(XY) resolution in the region [3,10] deg at absorber end without Branson correction
  fGaus->SetParameters(1., 0., 45.);
  hPDXAbsIn_3_10_Deg->Fit("fGaus");
  Double_t sDX310 = fGaus->GetParameter(2);
  fGaus->SetParameters(1., 0., 45.);
  hPDYAbsIn_3_10_Deg->Fit("fGaus");
  Double_t sDY310 = fGaus->GetParameter(2);
  fPGaus->SetParameters(1., 0., 45.);
  hPDAbsIn_3_10_Deg->Fit("fPGaus");
  FitGausResVsP(hPDXAbsInVsP_3_10_Deg, pNBins, 0., 45., "p*DX (tracks in [3,10] deg.)", gMeanPDXAbsInVsP_3_10_Deg, gSigmaPDXAbsInVsP_3_10_Deg);
  FitGausResVsP(hPDYAbsInVsP_3_10_Deg, pNBins, 0., 45., "p*DY (tracks in [3,10] deg.)", gMeanPDYAbsInVsP_3_10_Deg, gSigmaPDYAbsInVsP_3_10_Deg);
  FitPGausResVsP(hPDAbsInVsP_3_10_Deg, pNBins, 0., 45., "p*D (tracks in [3,10] deg.)", gMeanPDAbsInVsP_3_10_Deg, gSigmaPDAbsInVsP_3_10_Deg);
  
  printf("sigma pDXY ~ %f GeV.cm (2-3°) / %f GeV.cm (3-10°)\n", 0.5*(sDX23+sDY23), 0.5*(sDX310+sDY310));
  
  // compute p*D(XY) resolution in the region [2,3] deg at absorber end with Branson correction
  fGaus->SetParameters(1., 0., 3.);
  hPDXAbsInWithB_2_3_Deg->Fit("fGaus");
  Double_t sDXwB23 = fGaus->GetParameter(2);
  fGaus->SetParameters(1., 0., 3.);
  hPDYAbsInWithB_2_3_Deg->Fit("fGaus");
  Double_t sDYwB23 = fGaus->GetParameter(2);
  fPGaus->SetParameters(1., 0., 3.);
  hPDAbsInWithB_2_3_Deg->Fit("fPGaus");
  FitGausResVsP(hPDXAbsInWithBVsP_2_3_Deg, pNBins, 0., 3., "p*DX (tracks in [2,3] deg.) with Branson", gMeanPDXAbsInWithBVsP_2_3_Deg, gSigmaPDXAbsInWithBVsP_2_3_Deg);
  FitGausResVsP(hPDYAbsInWithBVsP_2_3_Deg, pNBins, 0., 3., "p*DY (tracks in [2,3] deg.) with Branson", gMeanPDYAbsInWithBVsP_2_3_Deg, gSigmaPDYAbsInWithBVsP_2_3_Deg);
  FitPGausResVsP(hPDAbsInWithBVsP_2_3_Deg, pNBins, 0., 3., "p*D (tracks in [2,3] deg.) with Branson", gMeanPDAbsInWithBVsP_2_3_Deg, gSigmaPDAbsInWithBVsP_2_3_Deg);
  
  // compute p*D(XY) resolution in the region [3,10] deg at absorber end with Branson correction
  fGaus->SetParameters(1., 0., 2.5);
  hPDXAbsInWithB_3_10_Deg->Fit("fGaus");
  Double_t sDXwB310 = fGaus->GetParameter(2);
  fGaus->SetParameters(1., 0., 2.5);
  hPDYAbsInWithB_3_10_Deg->Fit("fGaus");
  Double_t sDYwB310 = fGaus->GetParameter(2);
  fPGaus->SetParameters(1., 0., 2.5);
  hPDAbsInWithB_3_10_Deg->Fit("fPGaus");
  FitGausResVsP(hPDXAbsInWithBVsP_3_10_Deg, pNBins, 0., 2.5, "p*DX (tracks in [3,10] deg.) with Branson", gMeanPDXAbsInWithBVsP_3_10_Deg, gSigmaPDXAbsInWithBVsP_3_10_Deg);
  FitGausResVsP(hPDYAbsInWithBVsP_3_10_Deg, pNBins, 0., 2.5, "p*DY (tracks in [3,10] deg.) with Branson", gMeanPDYAbsInWithBVsP_3_10_Deg, gSigmaPDYAbsInWithBVsP_3_10_Deg);
  FitPGausResVsP(hPDAbsInWithBVsP_3_10_Deg, pNBins, 0., 2.5, "p*D (tracks in [3,10] deg.) with Branson", gMeanPDAbsInWithBVsP_3_10_Deg, gSigmaPDAbsInWithBVsP_3_10_Deg);
  
  printf("sigma pDXYwB ~ %f GeV.cm (2-3°) / %f GeV.cm (3-10°)\n", 0.5*(sDXwB23+sDYwB23), 0.5*(sDXwB310+sDYwB310));
  
  // fit JPsi
  TCanvas *cMass = new TCanvas("cMass", "cMass");
  cMass->cd(1);
  gPad->SetLogy();
  hMass->Draw("e0");
  Double_t nEntries = hMass->GetEntries();
  fCB->SetParameters(nEntries, 3., 0.08, 1., 5., 1., 5.);
  hMass->Fit(fCB, "RIM+", "e0");
  Double_t nJPsi = fCB->Integral(0., 6.);
  Double_t binWidth = hMass->GetBinWidth(1);
  printf("JPsi fit: chi2/ndf = %4.2f, sigma = %5.3f ± %5.3f, NJPsi = %d ± %d, Entries = %d (%+3.1f%%)\n",
	 fCB->GetChisquare()/fCB->GetNDF(), fCB->GetParameter(2), fCB->GetParError(2),
	 (Int_t)(nJPsi/binWidth), (Int_t)(fCB->IntegralError(0., 6.)/binWidth),
	 (Int_t)nEntries, 100.*(nJPsi/binWidth-nEntries)/nEntries);
  
  // save output
  histoFile->cd();
  histoFile->Write();
  hMass->Write();
  histoFile->cd("pDCA");
  gMeanPDCAXVsP_2_3_Deg->Write();
  gSigmaPDCAXVsP_2_3_Deg->Write();
  gMeanPDCAYVsP_2_3_Deg->Write();
  gSigmaPDCAYVsP_2_3_Deg->Write();
  gMeanPDCAVsP_2_3_Deg->Write();
  gSigmaPDCAVsP_2_3_Deg->Write();
  gMeanPDCAXVsP_3_10_Deg->Write();
  gSigmaPDCAXVsP_3_10_Deg->Write();
  gMeanPDCAYVsP_3_10_Deg->Write();
  gSigmaPDCAYVsP_3_10_Deg->Write();
  gMeanPDCAVsP_3_10_Deg->Write();
  gSigmaPDCAVsP_3_10_Deg->Write();
  histoFile->cd("pDXYAtAbsIn");
  gMeanPDXAbsInVsP_2_3_Deg->Write();
  gSigmaPDXAbsInVsP_2_3_Deg->Write();
  gMeanPDYAbsInVsP_2_3_Deg->Write();
  gSigmaPDYAbsInVsP_2_3_Deg->Write();
  gMeanPDAbsInVsP_2_3_Deg->Write();
  gSigmaPDAbsInVsP_2_3_Deg->Write();
  gMeanPDXAbsInVsP_3_10_Deg->Write();
  gSigmaPDXAbsInVsP_3_10_Deg->Write();
  gMeanPDYAbsInVsP_3_10_Deg->Write();
  gSigmaPDYAbsInVsP_3_10_Deg->Write();
  gMeanPDAbsInVsP_3_10_Deg->Write();
  gSigmaPDAbsInVsP_3_10_Deg->Write();
  histoFile->cd("pDXYAtAbsInWithB");
  gMeanPDXAbsInWithBVsP_2_3_Deg->Write();
  gSigmaPDXAbsInWithBVsP_2_3_Deg->Write();
  gMeanPDYAbsInWithBVsP_2_3_Deg->Write();
  gSigmaPDYAbsInWithBVsP_2_3_Deg->Write();
  gMeanPDAbsInWithBVsP_2_3_Deg->Write();
  gSigmaPDAbsInWithBVsP_2_3_Deg->Write();
  gMeanPDXAbsInWithBVsP_3_10_Deg->Write();
  gSigmaPDXAbsInWithBVsP_3_10_Deg->Write();
  gMeanPDYAbsInWithBVsP_3_10_Deg->Write();
  gSigmaPDYAbsInWithBVsP_3_10_Deg->Write();
  gMeanPDAbsInWithBVsP_3_10_Deg->Write();
  gSigmaPDAbsInWithBVsP_3_10_Deg->Write();
  histoFile->Close();
  
}


//------------------------------------------------------------------------------------
void FitGausResVsP(TH2* h, Int_t nBins, const Double_t mean0, const Double_t sigma0,
		     const char* fitting, TGraphAsymmErrors* gMean, TGraphAsymmErrors* gSigma)
{
  /// generic function to fit residuals versus momentum with a gaussian
  
  Int_t rebinFactorX = TMath::Max(h->GetNbinsX()/nBins, 1);
  for (Int_t i = rebinFactorX; i <= h->GetNbinsX(); i+=rebinFactorX) {
    cout<<Form("\rFitting %s... %d/%d",fitting,i/rebinFactorX,nBins)<<flush;
    TH1D *tmp = h->ProjectionY("tmp",i-rebinFactorX+1,i,"e");
    fGaus->SetParameters(1., mean0, sigma0);
/*    tmp->Fit("fGaus","WWNQ");
    Int_t rebin = static_cast<Int_t>(TMath::Min(0.1*tmp->GetNbinsX(),TMath::Max(0.5*fGaus->GetParameter(2)/tmp->GetBinWidth(1),1.)));
*/    Int_t rebin = static_cast<Int_t>(0.5*sigma0*tmp->GetNbinsX()/(tmp->GetBinLowEdge(tmp->GetNbinsX()+1)-tmp->GetBinLowEdge(1)));
    while (tmp->GetNbinsX()%rebin!=0) rebin--;
    tmp->Rebin(rebin);
    tmp->Fit("fGaus","NQ");
    h->GetXaxis()->SetRange(i-rebinFactorX+1,i);
    Double_t p = (tmp->GetEntries() > 0) ? h->GetMean() : 0.5 * (h->GetBinLowEdge(i-rebinFactorX+1) + h->GetBinLowEdge(i+1));
    h->GetXaxis()->SetRange();
    Double_t pErr[2] = {p-h->GetBinLowEdge(i-rebinFactorX+1), h->GetBinLowEdge(i+1)-p};
    gMean->SetPoint(i/rebinFactorX-1, p, fGaus->GetParameter(1));
    gMean->SetPointError(i/rebinFactorX-1, pErr[0], pErr[1], fGaus->GetParError(1), fGaus->GetParError(1));
    gSigma->SetPoint(i/rebinFactorX-1, p, fGaus->GetParameter(2));
    gSigma->SetPointError(i/rebinFactorX-1, pErr[0], pErr[1], fGaus->GetParError(2), fGaus->GetParError(2));
    delete tmp;
  }
  cout<<Form("\rFitting %s... %d/%d",fitting,nBins,nBins)<<endl;
  
}

//------------------------------------------------------------------------------------
void FitPGausResVsP(TH2* h, Int_t nBins, const Double_t mean0, const Double_t sigma0,
		    const char* fitting, TGraphAsymmErrors* gMean, TGraphAsymmErrors* gSigma)
{
  /// generic function to fit p*residuals versus momentum with p*gaussian
  
  Int_t rebinFactorX = TMath::Max(h->GetNbinsX()/nBins, 1);
  for (Int_t i = rebinFactorX; i <= h->GetNbinsX(); i+=rebinFactorX) {
    cout<<Form("\rFitting %s... %d/%d",fitting,i/rebinFactorX,nBins)<<flush;
    TH1D *tmp = h->ProjectionY("tmp",i-rebinFactorX+1,i,"e");
    fPGaus->SetParameters(1., mean0, sigma0);
    Int_t rebin = static_cast<Int_t>(0.5*sigma0*tmp->GetNbinsX()/(tmp->GetBinLowEdge(tmp->GetNbinsX()+1)-tmp->GetBinLowEdge(1)));
    while (tmp->GetNbinsX()%rebin!=0) rebin--;
    tmp->Rebin(rebin);
    tmp->Fit("fPGaus","NQ");
    h->GetXaxis()->SetRange(i-rebinFactorX+1,i);
    Double_t p = (tmp->GetEntries() > 0) ? h->GetMean() : 0.5 * (h->GetBinLowEdge(i-rebinFactorX+1) + h->GetBinLowEdge(i+1));
    h->GetXaxis()->SetRange();
    Double_t pErr[2] = {p-h->GetBinLowEdge(i-rebinFactorX+1), h->GetBinLowEdge(i+1)-p};
    gMean->SetPoint(i/rebinFactorX-1, p, fPGaus->GetParameter(1));
    gMean->SetPointError(i/rebinFactorX-1, pErr[0], pErr[1], fPGaus->GetParError(1), fPGaus->GetParError(1));
    gSigma->SetPoint(i/rebinFactorX-1, p, fPGaus->GetParameter(2));
    gSigma->SetPointError(i/rebinFactorX-1, pErr[0], pErr[1], fPGaus->GetParError(2), fPGaus->GetParError(2));
    delete tmp;
  }
  cout<<Form("\rFitting %s... %d/%d",fitting,nBins,nBins)<<endl;
  
}

//------------------------------------------------------------------------------------
Double_t CrystalBallExtended(Double_t *x, Double_t *par)
{
  //par[0] = Normalization
  //par[1] = mean
  //par[2] = sigma
  //par[3] = alpha
  //par[4] = n
  //par[5] = alpha'
  //par[6] = n'  
  
  Double_t t = (x[0]-par[1])/par[2];
  if (par[3] < 0) t = -t;
  
  Double_t absAlpha = fabs((Double_t)par[3]);
  Double_t absAlpha2 = fabs((Double_t)par[5]);
  
  if (t >= -absAlpha && t < absAlpha2) // gaussian core
  {
    return par[0]*(exp(-0.5*t*t));
  }
  
  if (t < -absAlpha) //left tail
  {
    Double_t a =  TMath::Power(par[4]/absAlpha,par[4])*exp(-0.5*absAlpha*absAlpha);
    Double_t b = par[4]/absAlpha - absAlpha;
    return par[0]*(a/TMath::Power(b - t, par[4]));
  }
  
  if (t >= absAlpha2) //right tail
  {
    
    Double_t c =  TMath::Power(par[6]/absAlpha2,par[6])*exp(-0.5*absAlpha2*absAlpha2);
    Double_t d = par[6]/absAlpha2 - absAlpha2;
    return par[0]*(c/TMath::Power(d + t, par[6]));
  }
  
  return 0. ; 
} 
