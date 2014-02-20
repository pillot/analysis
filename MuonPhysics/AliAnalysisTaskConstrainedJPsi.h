#ifndef ALIANALYSISTASKCONSTRAINEDJPSI_H
#define ALIANALYSISTASKCONSTRAINEDJPSI_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/// \ingroup muondep
/// \class AliAnalysisTaskConstrainedJPsi
/// \brief task to compute J/Psi kinematics constrained by the J/Psi PDG mass
//Author: Philippe Pillot - SUBATECH Nantes

#include "TMatrixD.h"
#include "AliAnalysisTaskSE.h"
#include "AliMuonPairCuts.h"

class AliMUONTrackParam;

class AliAnalysisTaskConstrainedJPsi : public AliAnalysisTaskSE {
public:
  
  AliAnalysisTaskConstrainedJPsi();
  AliAnalysisTaskConstrainedJPsi(const char *name);
  virtual ~AliAnalysisTaskConstrainedJPsi();
  
  // Interface methods
  virtual void UserCreateOutputObjects();
  virtual void NotifyRun();
  virtual void UserExec(Option_t *);
  virtual void Terminate(Option_t *);
  
  /// Set location of the default OCDB storage (if not set use "raw://")
  void SetDefaultStorage(const char* ocdbPath) { fDefaultStorage = ocdbPath; }
  
  // set standard cuts to select pairs to be considered
  void SetMuonPairCuts(AliMuonPairCuts &pairCuts);
  
  // set the chamber resolution correction factors
  void SetChResCorrFactor(Double_t coorX, Double_t corrY);
  
  // set the vertex z resolution
  void SetVtxZRes(Double_t val);
  
private:
  
  /// Not implemented
  AliAnalysisTaskConstrainedJPsi(const AliAnalysisTaskConstrainedJPsi& rhs);
  /// Not implemented
  AliAnalysisTaskConstrainedJPsi& operator = (const AliAnalysisTaskConstrainedJPsi& rhs);
  
  // correct the covariance matrix in the coordinate system (X, SlopeX, Y, SlopeY, q/Pyz)
  void CorrectChRes(TMatrixD &cov);
  
  // return track parameters and covariances extrapolated to the given vertex
  void GetTrackParamAtVtx(const AliESDMuonTrack &esdTrack, const Double_t vtx[3], AliMUONTrackParam &param);
  
  // return muons and J/Psi parameters
  void GetJPsiParam(const AliMUONTrackParam &param1, const AliMUONTrackParam &param2,
		    Double_t vtxZ, TMatrixD &paramMu, TMatrixD &paramJPsi);
  
  // change coordinate system: (X, SlopeX, Y, SlopeY, q/Pyz) -> (X, Y, pX, pY, pZ)
  void Cov2CovP(const AliMUONTrackParam &param, TMatrixD &covP);
  
  // change coordinate system: (pX1, pY1, pZ1, pX2, pY2, pZ2, z) -> (pX, pY, pZ, M, pX1-pX2, pY1-pY2, z)
  void CovMu2CovJPsi(const TMatrixD &paramMu, const TMatrixD &covMu, const TMatrixD &paramJPsi,
		     const TMatrixD &paramJPsidZ, Double_t dZ, TMatrixD &covJPsi);
  
  // compute the momentum resolution
  Double_t SigmaP(const TMatrixD &paramJPsi, const TMatrixD &covJPsi);
  
  // compute the trasverse momentum resolution
  Double_t SigmaPt(const TMatrixD &paramJPsi, const TMatrixD &covJPsi);
  
  // compute the rapidity resolution
  Double_t SigmaY(const TMatrixD &paramJPsi, const TMatrixD &covJPsi);
  
  // compute new track parameters and covariances constrained by the J/Psi mass
  Double_t AddMassConstraint(const TMatrixD &paramJPsi, const TMatrixD &covJPsi,
			     TMatrixD &cParamJPsi, TMatrixD &cCovJPsi);
  
private:
  
  // J/Psi kinematics distributions
  enum eData {
    kMass = 0, ///< invariant mass
    kP    = 1, ///< P
    kPt   = 2, ///< Pt
    kY    = 3, ///< rapidity
    kZvtx = 4  ///< z-vertex
  };
  
  // JPsi kinematics resolutions given by the Kalman filter
  enum eRes {
    kMassRes   = 0, ///< invariant mass
    kPResVsP   = 1, ///< P versus P
    kPtResVsPt = 2, ///< Pt versus Pt
    kYRes      = 3, ///< rapidity
    kZvtxRes   = 4, ///< z-vertex
    kChi2      = 5  ///< normalized chi2
  };
  
  Double_t fMMu;   //!< muon mass
  Double_t fMJPsi; //!< J/Psi mass
  
  TObjArray*  fData;      //!< List of reference data (not constrained)
  TObjArray*  fConstData; //!< List of constrained data
  TObjArray*  fResK;      //!< List of reference resolution from Kalman (not constrained)
  TObjArray*  fConstResK; //!< List of constrained resolution from Kalman
  TObjArray*  fRes;       //!< List of reference resolution compared to MC (not constrained)
  TObjArray*  fConstRes;  //!< List of constrained resolution compared to MC
  
  TString  fDefaultStorage; ///< location of the default OCDB storage
  
  AliMuonPairCuts* fMuonPairCuts; ///< cuts to select pairs to be considered
  
  Double_t fChResCorrFactor[5]; ///< chamber resolution correction factors
  Double_t fVtxZRes;            ///< vertex z resolution
  
  ClassDef(AliAnalysisTaskConstrainedJPsi, 1);
};


//________________________________________________________________________
inline void AliAnalysisTaskConstrainedJPsi::SetMuonPairCuts(AliMuonPairCuts &pairCuts)
{
  /// set standard cuts to select pairs to be considered
  delete fMuonPairCuts;
  fMuonPairCuts = new AliMuonPairCuts(pairCuts);
}


//________________________________________________________________________
inline void AliAnalysisTaskConstrainedJPsi::SetChResCorrFactor(Double_t coorX, Double_t corrY)
{
  /// set the chamber resolution correction factors
  /// in the coordinate system (X, SlopeX, Y, SlopeY, q/Pyz)
  fChResCorrFactor[0] = coorX;
  fChResCorrFactor[1] = coorX;
  fChResCorrFactor[2] = corrY;
  fChResCorrFactor[3] = corrY;
  fChResCorrFactor[4] = corrY;
}


//________________________________________________________________________
inline void AliAnalysisTaskConstrainedJPsi::SetVtxZRes(Double_t val)
{
  /// set the vertex z resolution.
  /// The value cannot be too low otherwise it is too different from other
  /// resolutions and the covariance matrix sometimes cannot be inverted (Det=0)
  if (val < 1.e-4) fVtxZRes = 1.e-4;
  else fVtxZRes = val;
}

#endif

