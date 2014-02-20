/*
 *  MuonGenerator.C
 *  aliroot_dev
 *
 *  Created by philippe pillot on 05/03/13.
 *  Copyright 2013 SUBATECH. All rights reserved.
 *
 */


#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TRandom.h"
#include "AliGenerator.h"
#include "AliGenParam.h"
#endif 


static Double_t PtMuon( const Double_t *x, const Double_t */*dummy*/ );
static Double_t PtMuonLizardo( const Double_t *x, const Double_t */*dummy*/ );
static Double_t YMuon( const Double_t *x, const Double_t */*dummy*/ );
static Double_t YMuonLizardo( const Double_t *x, const Double_t */*dummy*/ );
static Double_t V2Muon( const Double_t */*dummy*/, const Double_t */*dummy*/ );
static Int_t IpMuon( TRandom *ran );


//-------------------------------------------------------------------------
AliGenerator* MuonGenerator()
{
  printf("using MuonGenerator()...\n");
  //  AliGenParam *singleMu = new AliGenParam(10, AliGenMUONlib::kMuon,"","");
  //AliGenParam *singleMu = new AliGenParam(10,-1,PtMuon,YMuon,V2Muon,IpMuon);
  AliGenParam *singleMu = new AliGenParam(3,-1,PtMuonLizardo,YMuonLizardo,V2Muon,IpMuon);
  singleMu->SetMomentumRange(0., 1.e6);
  singleMu->SetPtRange(0.8, 999.);
//  singleMu->SetYRange(-4.3, -2.3);
  singleMu->SetYRange(-4.2, -2.3);
  singleMu->SetPhiRange(0., 360.);
  //for test only
  //singleMu->SetPtRange(1,999.);
  //singleMu->SetYRange(-4., -2.5);
  singleMu->SetForceDecay(kNoDecay);
  singleMu->SetTrackingFlag(1);
  
  return singleMu;
}

//-------------------------------------------------------------------------
Double_t PtMuon( const Double_t *x, const Double_t */*dummy*/ )
{
  // muon pT
  Double_t pT = *x;
//  Double_t p[6] = {371.909, 0.84614, 0.560486, 9.34831, 0.000474983, -0.853963}; // LHC13de tune0
  Double_t p[6] = {371.665, 0.845642, 0.56192, 9.34859, 0.000474519, -0.851091}; // LHC13de tune1
//  Double_t p[6] = {455.614, 0.942071, 0.706755, 8.69856, 0.000168775, -0.925487}; // LHC13f tune1
  return p[0] * (1. / TMath::Power(p[1] + TMath::Power(pT,p[2]), p[3]) + p[4] * TMath::Exp(p[5]*pT));
}

//-------------------------------------------------------------------------
Double_t PtMuonLizardo( const Double_t *x, const Double_t */*dummy*/ )
{
  Double_t pT=*x;
  Double_t p[4] = {2.24269, 1., 4.80675, 0.943831};
  return p[0] / TMath::Power( p[1] + TMath::Power(pT,p[2]), p[3] );
}

//-------------------------------------------------------------------------
Double_t YMuon( const Double_t *x, const Double_t */*dummy*/ )
{
  // muon y
  Double_t y = *x;
//  Double_t p[8] = {0.539134, 1, 0, 0.0499378, 0, -0.00450342, 0, 2}; // LHC13de tune0
  Double_t p[8] = {0.777922, 1, 0, -0.0184202, 0, -0.00107081, 0, 2}; // LHC13de tune1
//  Double_t p[8] = {1.29511, 1., 0., -0.0767846, 0., 0.00176313, 0., 2.}; // LHC13f tune1
  Double_t arg = y/p[7];
  return p[0] * (p[1] * (1. + p[2]*y + p[3]*y*y + p[4]*y*y*y + p[5]*y*y*y*y) + p[6]*TMath::Exp(-0.5*arg*arg));
}

//-------------------------------------------------------------------------
Double_t YMuonLizardo( const Double_t *x, const Double_t */*dummy*/ )
{
  Double_t y = *x;
  //pol3 only valid in y= -4;-2.5
  Double_t p[4] = {-0.0659365, -0.173503, -0.0711377, -0.00823051};
  return p[0] + p[1]*y + p[2]*y*y + p[3]*y*y*y;
}

//-------------------------------------------------------------------------
Double_t V2Muon( const Double_t */*dummy*/, const Double_t */*dummy*/ )
{
  //muon v2
  return 0.;
}

//-------------------------------------------------------------------------
Int_t IpMuon( TRandom *ran )
{
  //muon composition
  if (ran->Rndm() < 0.5 ) {
    return 13;
  }
  else {
    return -13;
  }
}

