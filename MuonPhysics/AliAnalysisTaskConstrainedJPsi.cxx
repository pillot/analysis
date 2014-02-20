/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// ROOT includes
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TObjArray.h"
#include "TList.h"
#include "TString.h"
#include "TMath.h"
#include "TGeoGlobalMagField.h"
#include "TGeoManager.h"
#include "TDatabasePDG.h"
#include "TParticle.h"

// STEER includes
#include "AliCDBManager.h"
#include "AliGeomManager.h"
#include "AliESDEvent.h"
#include "AliESDMuonTrack.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"

// ANALYSIS includes
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskConstrainedJPsi.h"

// MUON includes
#include "AliMUONCDB.h"
#include "AliMUONTrackParam.h"
#include "AliMUONTrackExtrap.h"
#include "AliMUONESDInterface.h"

ClassImp(AliAnalysisTaskConstrainedJPsi)

//________________________________________________________________________
AliAnalysisTaskConstrainedJPsi::AliAnalysisTaskConstrainedJPsi() :
AliAnalysisTaskSE(), 
fMMu(TDatabasePDG::Instance()->GetParticle("mu-")->Mass()),
fMJPsi(TDatabasePDG::Instance()->GetParticle("J/psi")->Mass()),
fData(0x0),
fConstData(0x0),
fResK(0x0),
fConstResK(0x0),
fRes(0x0),
fConstRes(0x0),
fDefaultStorage(""),
fMuonPairCuts(0x0),
fVtxZRes(1.e-4)
{
  /// Dummy constructor
  for (Int_t i = 0; i < 5; i++) fChResCorrFactor[i] = 1.;
}

//________________________________________________________________________
AliAnalysisTaskConstrainedJPsi::AliAnalysisTaskConstrainedJPsi(const char *name) :
AliAnalysisTaskSE(name), 
fMMu(TDatabasePDG::Instance()->GetParticle("mu-")->Mass()),
fMJPsi(TDatabasePDG::Instance()->GetParticle("J/psi")->Mass()),
fData(0x0),
fConstData(0x0),
fResK(0x0),
fConstResK(0x0),
fRes(0x0),
fConstRes(0x0),
fDefaultStorage("raw://"),
fMuonPairCuts(0x0),
fVtxZRes(1.e-4)
{
  /// Constructor
  for (Int_t i = 0; i < 5; i++) fChResCorrFactor[i] = 1.;
  
  // Output slot #1 writes into a TObjArray container
  DefineOutput(1,TObjArray::Class());
  // Output slot #2 writes into a TObjArray container
  DefineOutput(2,TObjArray::Class());
  // Output slot #3 writes into a TObjArray container
  DefineOutput(3,TObjArray::Class());
  // Output slot #4 writes into a TObjArray container
  DefineOutput(4,TObjArray::Class());
  // Output slot #5 writes into a TObjArray container
  DefineOutput(5,TObjArray::Class());
  // Output slot #6 writes into a TObjArray container
  DefineOutput(6,TObjArray::Class());
}

//________________________________________________________________________
AliAnalysisTaskConstrainedJPsi::~AliAnalysisTaskConstrainedJPsi()
{
  /// Destructor
  if (!AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
    delete fData;
    delete fConstData;
    delete fResK;
    delete fConstResK;
    delete fRes;
    delete fConstRes;
  }
  delete fMuonPairCuts;
}

//___________________________________________________________________________
void AliAnalysisTaskConstrainedJPsi::UserCreateOutputObjects()
{
  /// Create histograms
  
  fData = new TObjArray(20);
  fData->SetOwner();
  fConstData = new TObjArray(20);
  fConstData->SetOwner();
  fResK = new TObjArray(20);
  fResK->SetOwner();
  fConstResK = new TObjArray(20);
  fConstResK->SetOwner();
  fRes = new TObjArray(20);
  fRes->SetOwner();
  fConstRes = new TObjArray(20);
  fConstRes->SetOwner();
  
  TH1F *h1 = 0x0;
  TH2F *h2 = 0x0;
  TObjArray *dataObj[2] = {fData, fConstData};
  TObjArray *resKObj[2] = {fResK, fConstResK};
  TObjArray *resObj[2] = {fRes, fConstRes};
  TString dataName[2] = {"", "C"};
  TString dataTitle[2] = {"", "constrained "};
  
  for (Int_t i = 0; i < 2; i++) {
    
    // kinematics
    h1 = new TH1F(Form("h%sMass",dataName[i].Data()),
		  Form("%sinvariant mass distribution;mass (GeV/c^{2})",dataTitle[i].Data()),
		  600, 0., 6.);
    dataObj[i]->AddAtAndExpand(h1, kMass);
    
    h1 = new TH1F(Form("h%sP",dataName[i].Data()),
		  Form("%smomentum distribution;p (GeV/c)",dataTitle[i].Data()),
		  300, 0., 300.);
    dataObj[i]->AddAtAndExpand(h1, kP);
    
    h1 = new TH1F(Form("h%sPt",dataName[i].Data()),
		  Form("%stransverse momentum distribution;p_{t} (GeV/c)",dataTitle[i].Data()),
		  300, 0., 30.);
    dataObj[i]->AddAtAndExpand(h1, kPt);
    
    h1 = new TH1F(Form("h%sY",dataName[i].Data()),
		  Form("%srapidity distribution;y",dataTitle[i].Data()),
		  150, -4.5, -2.);
    dataObj[i]->AddAtAndExpand(h1, kY);
    
    h1 = new TH1F(Form("h%sZvtx",dataName[i].Data()),
		  Form("%sz-vertex distribution;z (cm)",dataTitle[i].Data()),
		  400, -20., 20.);
    dataObj[i]->AddAtAndExpand(h1, kZvtx);
    
    // resolutions from Kalman
    h1 = new TH1F(Form("h%sMassResK",dataName[i].Data()),
		  Form("%sinvariant mass resolution;#sigma_{mass} (GeV/c^{2})",dataTitle[i].Data()),
		  100, 0., 1.);
    resKObj[i]->AddAtAndExpand(h1, kMassRes);
    
    h2 = new TH2F(Form("h%sPResKVsP",dataName[i].Data()),
		  Form("%smomentum resolution vs p;p (GeV/c); #sigma_{p}/p (%%)",dataTitle[i].Data()),
		  300, 0., 300., 250, 0., 10.);
    resKObj[i]->AddAtAndExpand(h2, kPResVsP);
    
    h2 = new TH2F(Form("h%sPtResKVsPt",dataName[i].Data()),
		  Form("%stransverse momentum resolution vs p_{t};p_{t} (GeV/c); #sigma_{p_{t}}/p_{t} (%%)",dataTitle[i].Data()),
		  300, 0., 30., 300, 0., 30.);
    resKObj[i]->AddAtAndExpand(h2, kPtResVsPt);
    
    h1 = new TH1F(Form("h%sYResK",dataName[i].Data()),
		  Form("%srapidity resolution;#sigma_{y}",dataTitle[i].Data()),
		  100, 0., 0.1);
    resKObj[i]->AddAtAndExpand(h1, kYRes);
    
    h1 = new TH1F(Form("h%sZvtxResK",dataName[i].Data()),
		  Form("%sz-vertex resolution;#sigma_z (cm)",dataTitle[i].Data()),
		  100, 0., 10.);
    resKObj[i]->AddAtAndExpand(h1, kZvtxRes);
    
    h1 = new TH1F(Form("h%sChi2",dataName[i].Data()),
		  Form("%snormalized #chi^{2} distribution;#chi^{2} / ndf",dataTitle[i].Data()),
		  500, 0., 50.);
    resKObj[i]->AddAtAndExpand(h1, kChi2);
    
    // resolutions compared to MC
    h1 = new TH1F(Form("h%sMassRes",dataName[i].Data()),
		  Form("%sinvariant mass resolution;#Delta_{mass} (GeV/c^{2})",dataTitle[i].Data()),
		  200, -2., 2.);
    resObj[i]->AddAtAndExpand(h1, kMassRes);
    
    h2 = new TH2F(Form("h%sPResVsP",dataName[i].Data()),
		  Form("%smomentum resolution vs p;p (GeV/c); #Delta_{p}/p (%%)",dataTitle[i].Data()),
		  300, 0., 300., 250, -10., 10.);
    resObj[i]->AddAtAndExpand(h2, kPResVsP);
    
    h2 = new TH2F(Form("h%sPtResVsPt",dataName[i].Data()),
		  Form("%stransverse momentum resolution vs p_{t};p_{t} (GeV/c); #Delta_{p_{t}}/p_{t} (%%)",dataTitle[i].Data()),
		  300, 0., 30., 300, -30., 30.);
    resObj[i]->AddAtAndExpand(h2, kPtResVsPt);
    
    h1 = new TH1F(Form("h%sYRes",dataName[i].Data()),
		  Form("%srapidity resolution;#Delta_{y}",dataTitle[i].Data()),
		  200, -0.2, 0.2);
    resObj[i]->AddAtAndExpand(h1, kYRes);
    
    h1 = new TH1F(Form("h%sZvtxRes",dataName[i].Data()),
		  Form("%sz-vertex resolution;#Delta_z (cm)",dataTitle[i].Data()),
		  200, -20., 20.);
    resObj[i]->AddAtAndExpand(h1, kZvtxRes);
    
  }
  
  // Post data at least once per task to ensure data synchronisation (required for merging)
  PostData(1, fData);
  PostData(2, fConstData);
  PostData(3, fResK);
  PostData(4, fConstResK);
  PostData(5, fRes);
  PostData(6, fConstRes);
}

//________________________________________________________________________
void AliAnalysisTaskConstrainedJPsi::NotifyRun()
{
  /// load OCDB objects
  
  // set OCDB location
  AliCDBManager* cdbm = AliCDBManager::Instance();
  if (cdbm->IsDefaultStorageSet()) printf("ConstrainedJPsiTask: CDB default storage already set!\n");
  else cdbm->SetDefaultStorage(fDefaultStorage.Data());
  if (cdbm->GetRun() > -1) printf("ConstrainedJPsiTask: run number already set!\n");
  else cdbm->SetRun(fCurrentRunNumber);
  
  // load geometry for track extrapolation
  if (!AliGeomManager::GetGeometry()) {
    AliGeomManager::LoadGeometry();
    if (!AliGeomManager::GetGeometry()) return;  
    if (!AliGeomManager::ApplyAlignObjsFromCDB("MUON")) return;
  }
  
  // load magnetic field for track extrapolation
  if (!TGeoGlobalMagField::Instance()->GetField()) {
    if (!AliMUONCDB::LoadField()) return;
  }
  
  // set the magnetic field for track extrapolations
  AliMUONTrackExtrap::SetField();
  
  // get the pairCuts for this run
  if (fMuonPairCuts) fMuonPairCuts->SetRun(fCurrentRunNumber);
}

//________________________________________________________________________
void AliAnalysisTaskConstrainedJPsi::UserExec(Option_t *)
{
  /// Called for each event
  
  AliESDEvent* esd = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!esd) return;
  
  // Load MC event if any
  AliMCEventHandler* mcH = static_cast<AliMCEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  TParticle *jPsi = mcH ? mcH->MCEvent()->Stack()->Particle(0) : 0x0;
  Double_t pMC = 0., pTMC = 0., yMC = 0., mMC = 0., vZMC = 0.;
  if (jPsi) {
    pMC = jPsi->P();
    pTMC = jPsi->Pt();
    yMC = jPsi->Y();
    mMC = jPsi->GetMass();
    vZMC = jPsi->Vz();
  }
  
  // first loop over tracks
  TList muons;
  muons.SetOwner(kFALSE);
  AliMUONTrackParam param1AtVtx, param2AtVtx, param1AtVtxdZ, param2AtVtxdZ;
  Int_t previousTrack1 = -1;
  Int_t nTracks = (Int_t)esd->GetNumberOfMuonTracks();
  for (Int_t iTrack1 = 0; iTrack1 < nTracks; iTrack1++) {
    
    // get the ESD track seleted for physics
    AliESDMuonTrack* esdTrack1 = esd->GetMuonTrack(iTrack1);
    if (!esdTrack1->ContainTrackerData()) continue; // make sure to skip ghosts
    muons.Clear();
    muons.AddLast(esdTrack1);
    
    // second loop over tracks
    for (Int_t iTrack2 = iTrack1 + 1; iTrack2 < nTracks; iTrack2++) {
      
      // get the ESD track seleted for physics
      AliESDMuonTrack* esdTrack2 = esd->GetMuonTrack(iTrack2);
      if (!esdTrack2->ContainTrackerData()) continue; // make sure to skip ghosts
      if (muons.GetSize() == 2) muons.RemoveLast();
      muons.AddLast(esdTrack2);
      
      // apply standard pair cuts if any
      if (fMuonPairCuts && !fMuonPairCuts->IsSelected(&muons)) continue;
      
      // get track parameters at first cluster and extrapolate them to the vertex and to the vertex - sigma_z
      Double_t vertex[3] = {esdTrack1->GetNonBendingCoor(), esdTrack1->GetBendingCoor(), esdTrack1->GetZ()};
      Double_t vertexdZ[3] = {vertex[0], vertex[1], vertex[2]-fVtxZRes};
      if (iTrack1 != previousTrack1) { // do it only once per track
	GetTrackParamAtVtx(*esdTrack1, vertex, param1AtVtx);
	GetTrackParamAtVtx(*esdTrack1, vertexdZ, param1AtVtxdZ);
	previousTrack1 = iTrack1;
      }
      GetTrackParamAtVtx(*esdTrack2, vertex, param2AtVtx);
      GetTrackParamAtVtx(*esdTrack2, vertexdZ, param2AtVtxdZ);
      
      // compute J/Psi parameters
      TMatrixD paramMu(7,1);
      TMatrixD paramJPsi(7,1);
      GetJPsiParam(param1AtVtx, param2AtVtx, vertex[2], paramMu, paramJPsi);
      Double_t p2 = paramJPsi(0,0)*paramJPsi(0,0) + paramJPsi(1,0)*paramJPsi(1,0) + paramJPsi(2,0)*paramJPsi(2,0);
      Double_t e  = TMath::Sqrt(p2 + paramJPsi(3,0)*paramJPsi(3,0));
      Double_t p  = TMath::Sqrt(p2);
      Double_t pT = TMath::Sqrt(paramJPsi(0,0)*paramJPsi(0,0) + paramJPsi(1,0)*paramJPsi(1,0));
      Double_t y  = 0.5 * TMath::Log((e+paramJPsi(2,0)) / (e-paramJPsi(2,0)));
      
      // fill histos
      ((TH1F*)fData->UncheckedAt(kP))->Fill(p);
      ((TH1F*)fData->UncheckedAt(kPt))->Fill(pT);
      ((TH1F*)fData->UncheckedAt(kY))->Fill(y);
      ((TH1F*)fData->UncheckedAt(kMass))->Fill(paramJPsi(3,0));
      ((TH1F*)fData->UncheckedAt(kZvtx))->Fill(paramJPsi(6,0));
      if (jPsi) {
	((TH2F*)fRes->UncheckedAt(kPResVsP))->Fill(pMC,100.*(p-pMC)/pMC);
	((TH2F*)fRes->UncheckedAt(kPtResVsPt))->Fill(pTMC,100.*(pT-pTMC)/pTMC);
	((TH1F*)fRes->UncheckedAt(kYRes))->Fill(y-yMC);
	((TH1F*)fRes->UncheckedAt(kMassRes))->Fill(paramJPsi(3,0)-mMC);
	((TH1F*)fRes->UncheckedAt(kZvtxRes))->Fill(paramJPsi(6,0)-vZMC);
      }
      
      // get muon covariance matrix
      TMatrixD covMu(7,7);
      covMu.Zero();
      TMatrixD covMu1(5,5);
      Cov2CovP(param1AtVtx, covMu1);
      covMu.SetSub(0,0,covMu1.GetSub(2,4,2,4));
      TMatrixD covMu2(5,5);
      Cov2CovP(param2AtVtx, covMu2);
      covMu.SetSub(3,3,covMu2.GetSub(2,4,2,4));
      covMu(6,6) = fVtxZRes*fVtxZRes;
      
      // compute J/Psi parameters at vertex + sigma_z
      TMatrixD dummy(7,1);
      TMatrixD paramJPsidZ(7,1);
      GetJPsiParam(param1AtVtxdZ, param2AtVtxdZ, vertexdZ[2], dummy, paramJPsidZ);
      
      // compute J/Psi covariance matrix
      TMatrixD covJPsi(7,7);
      CovMu2CovJPsi(paramMu, covMu, paramJPsi, paramJPsidZ, vertexdZ[2]-vertex[2], covJPsi);
      
      // fill histos
      ((TH2F*)fResK->UncheckedAt(kPResVsP))->Fill(p, 100.*SigmaP(paramJPsi, covJPsi)/p);
      ((TH2F*)fResK->UncheckedAt(kPtResVsPt))->Fill(pT, 100.*SigmaPt(paramJPsi, covJPsi)/pT);
      ((TH1F*)fResK->UncheckedAt(kYRes))->Fill(SigmaY(paramJPsi, covJPsi));
      ((TH1F*)fResK->UncheckedAt(kMassRes))->Fill(TMath::Sqrt(covJPsi(3,3)));
      ((TH1F*)fResK->UncheckedAt(kZvtxRes))->Fill(TMath::Sqrt(covJPsi(6,6)));
      ((TH1F*)fResK->UncheckedAt(kChi2))->Fill(0.5*(esdTrack1->GetNormalizedChi2()+esdTrack2->GetNormalizedChi2()));
      
      // compute new J/Psi parameters and covariances adding J/Psi mass constraint
      TMatrixD cParamJPsi(7,1);
      TMatrixD cCovJPsi(7,7);
      Double_t addChi2 = AddMassConstraint(paramJPsi, covJPsi, cParamJPsi, cCovJPsi);
      Double_t cPT = TMath::Sqrt(cParamJPsi(0,0)*cParamJPsi(0,0) + cParamJPsi(1,0)*cParamJPsi(1,0));
      Double_t cP  = TMath::Sqrt(cPT*cPT + cParamJPsi(2,0)*cParamJPsi(2,0));
      Double_t cE  = TMath::Sqrt(cP*cP + cParamJPsi(3,0)*cParamJPsi(3,0));
      Double_t cY  = 0.5 * TMath::Log((cE+cParamJPsi(2,0)) / (cE-cParamJPsi(2,0)));
      
      // fill histos
      ((TH1F*)fConstData->UncheckedAt(kP))->Fill(cP);
      ((TH1F*)fConstData->UncheckedAt(kPt))->Fill(cPT);
      ((TH1F*)fConstData->UncheckedAt(kY))->Fill(cY);
      ((TH1F*)fConstData->UncheckedAt(kMass))->Fill(cParamJPsi(3,0));
      ((TH1F*)fConstData->UncheckedAt(kZvtx))->Fill(cParamJPsi(6,0));
      ((TH2F*)fConstResK->UncheckedAt(kPResVsP))->Fill(cP, 100.*SigmaP(cParamJPsi, cCovJPsi)/cP);
      ((TH2F*)fConstResK->UncheckedAt(kPtResVsPt))->Fill(cPT, 100.*SigmaPt(cParamJPsi, cCovJPsi)/cPT);
      ((TH1F*)fConstResK->UncheckedAt(kYRes))->Fill(SigmaY(cParamJPsi, cCovJPsi));
      ((TH1F*)fConstResK->UncheckedAt(kMassRes))->Fill(TMath::Sqrt(cCovJPsi(3,3)));
      ((TH1F*)fConstResK->UncheckedAt(kZvtxRes))->Fill(TMath::Sqrt(cCovJPsi(6,6)));
      ((TH1F*)fConstResK->UncheckedAt(kChi2))->Fill(0.5*(esdTrack1->GetNormalizedChi2()+esdTrack2->GetNormalizedChi2())+addChi2);
      if (jPsi) {
	((TH2F*)fConstRes->UncheckedAt(kPResVsP))->Fill(pMC,100.*(cP-pMC)/pMC);
	((TH2F*)fConstRes->UncheckedAt(kPtResVsPt))->Fill(pTMC,100.*(cPT-pTMC)/pTMC);
	((TH1F*)fConstRes->UncheckedAt(kYRes))->Fill(cY-yMC);
	((TH1F*)fConstRes->UncheckedAt(kMassRes))->Fill(cParamJPsi(3,0)-mMC);
	((TH1F*)fConstRes->UncheckedAt(kZvtxRes))->Fill(cParamJPsi(6,0)-vZMC);
      }
      
    }
    
  }
  
  // Post final data. It will be written to a file with option "RECREATE"
  PostData(1, fData);
  PostData(2, fConstData);
  PostData(3, fResK);
  PostData(4, fConstResK);
  PostData(5, fRes);
  PostData(6, fConstRes);
}

//________________________________________________________________________
void AliAnalysisTaskConstrainedJPsi::Terminate(Option_t *)
{
  /// display results
  
  fRes = dynamic_cast<TObjArray*>(GetOutputData(5));
  fConstRes = dynamic_cast<TObjArray*>(GetOutputData(6));
  if (!fRes || !fConstRes) return;
  
  TH1 *hPRes = ((TH2F*)fRes->UncheckedAt(kPResVsP))->ProjectionY("hPRes",1,300,"e");
  TH1 *hCPRes = ((TH2F*)fConstRes->UncheckedAt(kPResVsP))->ProjectionY("hCPRes",1,300,"e");
  
  TCanvas *cP = new TCanvas("cP","momentum resolution");
  cP->cd();
  hPRes->SetLineColor(2);
  hPRes->Draw();
  hCPRes->Draw("sames");
  
  TH1 *hPtRes = ((TH2F*)fRes->UncheckedAt(kPtResVsPt))->ProjectionY("hPtRes",1,300,"e");
  TH1 *hCPtRes = ((TH2F*)fConstRes->UncheckedAt(kPtResVsPt))->ProjectionY("hCPtRes",1,300,"e");
  
  TCanvas *cPt = new TCanvas("cPt","pT resolution");
  cPt->cd();
  hPtRes->SetLineColor(2);
  hPtRes->Draw();
  hCPtRes->Draw("sames");
  
}

//__________________________________________________________________________
void AliAnalysisTaskConstrainedJPsi::CorrectChRes(TMatrixD &cov)
{
  /// correct the covariance matrix in the coordinate system (X, SlopeX, Y, SlopeY, q/Pyz)
  /// to account for the requested modifications of the chamber resolution
  for (Int_t i = 0; i < 5; i++) {
    for (Int_t j = 0; j < 5; j++) {
      cov(i,j) *= fChResCorrFactor[i]*fChResCorrFactor[j];
    }
  }
}

//__________________________________________________________________________
void AliAnalysisTaskConstrainedJPsi::GetTrackParamAtVtx(const AliESDMuonTrack &esdTrack, const Double_t vtx[3],
							AliMUONTrackParam &param)
{
  /// return track parameters and covariances extrapolated to the given vertex
  
  // set track parameters and covariances at first cluster
  AliMUONESDInterface::GetParamAtFirstCluster(esdTrack, param);
  AliMUONESDInterface::GetParamCov(esdTrack, param);
  
  // account for the requested modifications of the chamber resolution
  TMatrixD cov(param.GetCovariances());
  CorrectChRes(cov);
  param.SetCovariances(cov);
  
  // extrapolate to vertex
  AliMUONTrackExtrap::ExtrapToVertex(&param, vtx[0], vtx[1], vtx[2], 0., 0.);
}

//__________________________________________________________________________
void AliAnalysisTaskConstrainedJPsi::GetJPsiParam(const AliMUONTrackParam &param1, const AliMUONTrackParam &param2,
						  Double_t vtxZ, TMatrixD &paramMu, TMatrixD &paramJPsi)
{
  /// return muons and J/Psi parameters
  
  // get muon parameters
  paramMu(0,0) = param1.Px();
  paramMu(1,0) = param1.Py();
  paramMu(2,0) = param1.Pz();
  paramMu(3,0) = param2.Px();
  paramMu(4,0) = param2.Py();
  paramMu(5,0) = param2.Pz();
  paramMu(6,0) = vtxZ;
  
  // compute J/Psi parameters
  Double_t pX = paramMu(0,0) + paramMu(3,0);
  Double_t pY = paramMu(1,0) + paramMu(4,0);
  Double_t pZ = paramMu(2,0) + paramMu(5,0);
  Double_t p1 = param1.P();
  Double_t p2 = param2.P();
  Double_t e  = TMath::Sqrt(p1*p1 + fMMu*fMMu) + TMath::Sqrt(p2*p2 + fMMu*fMMu);
  Double_t pp = pX*pX + pY*pY + pZ*pZ;
  paramJPsi(0,0) = pX;
  paramJPsi(1,0) = pY;
  paramJPsi(2,0) = pZ;
  paramJPsi(3,0) = TMath::Sqrt(e*e - pp);
  paramJPsi(4,0) = paramMu(0,0) - paramMu(3,0);
  paramJPsi(5,0) = paramMu(1,0) - paramMu(4,0);
  paramJPsi(6,0) = paramMu(6,0);
}

//__________________________________________________________________________
void AliAnalysisTaskConstrainedJPsi::Cov2CovP(const AliMUONTrackParam &param, TMatrixD &covP)
{
  /// change coordinate system: (X, SlopeX, Y, SlopeY, q/Pyz) -> (X, Y, pX, pY, pZ)
  /// parameters (param) are given in the (X, SlopeX, Y, SlopeY, q/Pyz) coordinate system
  
  // Get useful parameters
  Double_t slopeX = param.GetNonBendingSlope();
  Double_t slopeY = param.GetBendingSlope();
  Double_t qOverPYZ = param.GetInverseBendingMomentum();
  Double_t pZ = param.Pz();
  
  // compute Jacobian to change the coordinate system from (X,SlopeX,Y,SlopeY,c/pYZ) to (X,Y,pX,pY,pZ)
  Double_t dpZdSlopeY = - qOverPYZ * qOverPYZ * pZ * pZ * pZ * slopeY;
  Double_t dpZdQOverPYZ = (qOverPYZ != 0.) ? - pZ / qOverPYZ : - FLT_MAX;
  TMatrixD jacob(5,5);
  jacob.Zero();
  jacob(0,0) = 1.;
  jacob(1,2) = 1.;
  jacob(2,1) = pZ;
  jacob(2,3) = slopeX * dpZdSlopeY;
  jacob(2,4) = slopeX * dpZdQOverPYZ;
  jacob(3,3) = pZ + slopeY * dpZdSlopeY;
  jacob(3,4) = slopeY * dpZdQOverPYZ;
  jacob(4,3) = dpZdSlopeY;
  jacob(4,4) = dpZdQOverPYZ;
  
  // compute covariances in new coordinate system
  TMatrixD tmp(param.GetCovariances(),TMatrixD::kMultTranspose,jacob);
  covP.Mult(jacob,tmp);
}

//__________________________________________________________________________
void AliAnalysisTaskConstrainedJPsi::CovMu2CovJPsi(const TMatrixD &paramMu, const TMatrixD &covMu,
						   const TMatrixD &paramJPsi, const TMatrixD &paramJPsidZ,
						   Double_t dZ, TMatrixD &covJPsi)
{
  /// change coordinate system: (pX1, pY1, pZ1, pX2, pY2, pZ2, z) -> (pX, pY, pZ, M, pX1-pX2, pY1-pY2, z)
  
  // compute Jacobian to change the coordinate system
  Double_t p1 = TMath::Sqrt(paramMu(0,0)*paramMu(0,0) + paramMu(1,0)*paramMu(1,0) + paramMu(2,0)*paramMu(2,0));
  Double_t p2 = TMath::Sqrt(paramMu(3,0)*paramMu(3,0) + paramMu(4,0)*paramMu(4,0) + paramMu(5,0)*paramMu(5,0));
  Double_t e1 = TMath::Sqrt(p1*p1 + fMMu*fMMu);
  Double_t e2 = TMath::Sqrt(p2*p2 + fMMu*fMMu);
  Double_t e  = e1 + e2;
  TMatrixD jacob(7,7);
  jacob.Zero();
  jacob(0,0) = 1.;
  jacob(0,3) = 1.;
  jacob(0,6) = (paramJPsidZ(0,0) - paramJPsi(0,0)) / dZ;
  jacob(1,1) = 1.;
  jacob(1,4) = 1.;
  jacob(1,6) = (paramJPsidZ(1,0) - paramJPsi(1,0)) / dZ;
  jacob(2,2) = 1.;
  jacob(2,5) = 1.;
  jacob(2,6) = (paramJPsidZ(2,0) - paramJPsi(2,0)) / dZ;
  jacob(3,0) = (paramMu(0,0)*e/e1 - paramJPsi(0,0)) / paramJPsi(3,0);
  jacob(3,1) = (paramMu(1,0)*e/e1 - paramJPsi(1,0)) / paramJPsi(3,0);
  jacob(3,2) = (paramMu(2,0)*e/e1 - paramJPsi(2,0)) / paramJPsi(3,0);
  jacob(3,3) = (paramMu(3,0)*e/e2 - paramJPsi(0,0)) / paramJPsi(3,0);
  jacob(3,4) = (paramMu(4,0)*e/e2 - paramJPsi(1,0)) / paramJPsi(3,0);
  jacob(3,5) = (paramMu(5,0)*e/e2 - paramJPsi(2,0)) / paramJPsi(3,0);
  jacob(3,6) = (paramJPsidZ(3,0) - paramJPsi(3,0)) / dZ;
  jacob(4,0) = 1.;
  jacob(4,3) = -1.;
  jacob(4,6) = (paramJPsidZ(4,0) - paramJPsi(4,0)) / dZ;
  jacob(5,1) = 1.;
  jacob(5,4) = -1.;
  jacob(5,6) = (paramJPsidZ(5,0) - paramJPsi(5,0)) / dZ;
  jacob(6,6) = 1.;
  
  // compute covariances in new coordinate system
  TMatrixD tmp(covMu,TMatrixD::kMultTranspose,jacob);
  covJPsi.Mult(jacob,tmp);
}

//__________________________________________________________________________
Double_t AliAnalysisTaskConstrainedJPsi::SigmaP(const TMatrixD &paramJPsi, const TMatrixD &covJPsi)
{
  /// compute the momentum resolution
  
  // compute Jacobian to change the coordinate system
  Double_t p = TMath::Sqrt(paramJPsi(0,0)*paramJPsi(0,0) + paramJPsi(1,0)*paramJPsi(1,0) + paramJPsi(2,0)*paramJPsi(2,0));
  TMatrixD jacob(1,7);
  jacob.Zero();
  jacob(0,0) = paramJPsi(0,0) / p;
  jacob(0,1) = paramJPsi(1,0) / p;
  jacob(0,2) = paramJPsi(2,0) / p;
  
  // compute the resolution
  TMatrixD tmp(covJPsi,TMatrixD::kMultTranspose,jacob);
  TMatrixD sigma(jacob,TMatrixD::kMult,tmp);
  
  return TMath::Sqrt(sigma(0,0));
}

//__________________________________________________________________________
Double_t AliAnalysisTaskConstrainedJPsi::SigmaPt(const TMatrixD &paramJPsi, const TMatrixD &covJPsi)
{
  /// compute the trasverse momentum resolution
  
  // compute Jacobian to change the coordinate system
  Double_t pt = TMath::Sqrt(paramJPsi(0,0)*paramJPsi(0,0) + paramJPsi(1,0)*paramJPsi(1,0));
  TMatrixD jacob(1,7);
  jacob.Zero();
  jacob(0,0) = paramJPsi(0,0) / pt;
  jacob(0,1) = paramJPsi(1,0) / pt;
  
  // compute the resolution
  TMatrixD tmp(covJPsi,TMatrixD::kMultTranspose,jacob);
  TMatrixD sigma(jacob,TMatrixD::kMult,tmp);
  
  return TMath::Sqrt(sigma(0,0));
}

//__________________________________________________________________________
Double_t AliAnalysisTaskConstrainedJPsi::SigmaY(const TMatrixD &paramJPsi, const TMatrixD &covJPsi)
{
  /// compute the rapidity resolution
  
  // compute Jacobian to change the coordinate system
  Double_t e2 = paramJPsi(0,0)*paramJPsi(0,0) + paramJPsi(1,0)*paramJPsi(1,0) +
                paramJPsi(2,0)*paramJPsi(2,0) + paramJPsi(3,0)*paramJPsi(3,0);
  Double_t e = TMath::Sqrt(e2);
  Double_t dYdPFactor = paramJPsi(2,0) / e / (e2 - paramJPsi(2,0)*paramJPsi(2,0));
  TMatrixD jacob(1,7);
  jacob.Zero();
  jacob(0,0) = -paramJPsi(0,0) * dYdPFactor;
  jacob(0,1) = -paramJPsi(1,0) * dYdPFactor;
  jacob(0,2) = 1 / e;
  jacob(0,3) = -paramJPsi(3,0) * dYdPFactor;
  
  // compute the resolution
  TMatrixD tmp(covJPsi,TMatrixD::kMultTranspose,jacob);
  TMatrixD sigma(jacob,TMatrixD::kMult,tmp);
  
  return TMath::Sqrt(sigma(0,0));
}

//__________________________________________________________________________
Double_t AliAnalysisTaskConstrainedJPsi::AddMassConstraint(const TMatrixD &paramJPsi, const TMatrixD &covJPsi,
							   TMatrixD &cParamJPsi, TMatrixD &cCovJPsi)
{
  /// compute new track parameters and their covariances including the J/Psi mass
  /// as a new measurement with infinite resolution using kalman filter.
  /// Return the additional track chi2
  
  // Set mass measurement parameters (m)
  TMatrixD massParam(7,1);
  massParam.Zero();
  massParam(3,0) = fMJPsi;
  
  // Compute the actual parameter weight (W)
  TMatrixD paramWeight(covJPsi);
  if (paramWeight.Determinant() != 0) {
    paramWeight.Invert();
  } else {
    AliWarning(" Determinant = 0");
    printf("paramWeight:\n");
    paramWeight.Print();
    return 1.e10;
  }
  
  // Compute the new measurement weight (U)
  TMatrixD massWeight(7,7);
  massWeight.Zero();
  massWeight(3,3) = 1.e10;
  
  // Compute the new parameters covariance matrix ( (W+U)^-1 )
  cCovJPsi.Plus(paramWeight,massWeight);
  if (cCovJPsi.Determinant() != 0) {
    cCovJPsi.Invert();
  } else {
    AliWarning(" Determinant = 0");
    printf("cCovJPsi:\n");
    cCovJPsi.Print();
    return 1.e10;
  }
  
  // Compute the new parameters (p' = ((W+U)^-1)U(m-p) + p)
  TMatrixD tmp(massParam,TMatrixD::kMinus,paramJPsi); // m-p
  TMatrixD tmp2(massWeight,TMatrixD::kMult,tmp); // U(m-p)
  cParamJPsi.Mult(cCovJPsi,tmp2); // ((W+U)^-1)U(m-p)
  cParamJPsi += paramJPsi; // ((W+U)^-1)U(m-p) + p
  
  // Compute the additional chi2 (= ((p'-p)^-1)W(p'-p) + ((p'-m)^-1)U(p'-m))
  tmp = cParamJPsi; // p'
  tmp -= paramJPsi; // (p'-p)
  TMatrixD tmp3(paramWeight,TMatrixD::kMult,tmp); // W(p'-p)
  TMatrixD addChi2(tmp,TMatrixD::kTransposeMult,tmp3); // ((p'-p)^-1)W(p'-p)
  tmp = cParamJPsi; // p'
  tmp -= massParam; // (p'-m)
  TMatrixD tmp4(massWeight,TMatrixD::kMult,tmp); // U(p'-m)
  addChi2 += TMatrixD(tmp,TMatrixD::kTransposeMult,tmp4); // ((p'-p)^-1)W(p'-p) + ((p'-m)^-1)U(p'-m)
  
  return addChi2(0,0);
  
}

