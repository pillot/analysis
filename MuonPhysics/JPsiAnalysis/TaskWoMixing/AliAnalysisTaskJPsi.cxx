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

#include <TH1.h>
#include <TList.h>
#include <TArrayF.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TObjString.h>

#include "AliLog.h"
#include "AliCounterCollection.h"
#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliVParticle.h"
#include "AliESDMuonTrack.h"

#include "AliMuonTrackCuts.h"

#include "AliAnalysisManager.h"
#include "AliAnalysisMuonUtility.h"
#include "AliAnalysisTaskJPsi.h"

ClassImp(MyList)
ClassImp(AliAnalysisTaskJPsi)

//________________________________________________________________________
MyList::MyList(const char *name, const char *title)
: TNamed(name,title),
fList(new TList)
{
  // constructor
  fList->SetOwner();
}

//________________________________________________________________________
MyList::~MyList()
{
  // destructor
  delete fList;
}

//________________________________________________________________________
void MyList::Add(TObject* obj)
{
  /// add a new object to this list
  if (!obj) return;
  fList->Add(obj);
}

//________________________________________________________________________
Long64_t MyList::Merge(TCollection* list)
{
  /// add the objects in the collection of lists to this one
  
  if (!list || list->IsEmpty()) return (Long64_t)fList->GetSize();
  
  TIter next(list);
  TObject* obj = 0x0;
  while ((obj = next())) {
    
    // check that "obj" is an object of the class MyList
    MyList *list2 = dynamic_cast<MyList*>(obj);
    if (!list2) {
      AliFatal(Form("object named \"%s\" is a %s instead of a MyList!", obj->GetName(), obj->ClassName()));
      continue;
    }
    
    // add objects to this one
    TIter next2(list2->fList);
    TObject* obj2 = 0x0;
    while ((obj2 = next2())) fList->Add(obj2);
    
    // the other list should no longer be owner of the objects
    list2->fList->SetOwner(kFALSE);
    
  }
  
  return (Long64_t)fList->GetSize();
  
}

//________________________________________________________________________
void MyList::Print(Option_t*) const
{
  /// print the internal list
  printf("internal content of %s list:\n",GetName());
  fList->Print();
}


//________________________________________________________________________
AliAnalysisTaskJPsi::AliAnalysisTaskJPsi()
: AliAnalysisTaskSE(),
fOutput(0),
fHistCent(0),
fEventCounters(0),
fMuonTrackCuts(0),
fTrackCutsSet(kFALSE),
fSelectTrgSign(kFALSE),
fSelectSameTrgSignFake(kFALSE),
fDeltaDev(0),
fUseMCLabel(kFALSE),
fnCent(0),
fRecordEvWithTrgIssues(kFALSE),
fTrgClassMissTrgL0Ev(0x0),
fTrgL0MissTrgClassEv(0x0),
fTrgClassMissTrgOffEv(0x0),
fTrgOffMissTrgClassEv(0x0),
fTrgOffMissTrgL0Ev(0x0),
fTrgL0MissTrgOffEv(0x0),
fOSTrkOSTrgMLLOnlyEv(0x0),
fOSTrkLSTrgMLLOnlyEv(0x0),
fOSTrkLSTrgMULEv(0x0),
fOSTrkOSTrgFakeMLLOnlyEv(0x0),
fOSTrkLSTrgFakeMLLOnlyEv(0x0),
fOSTrkLSTrgFakeMULEv(0x0)
{
  // Dummy constructor ALWAYS needed for I/O.
}

//________________________________________________________________________
AliAnalysisTaskJPsi::AliAnalysisTaskJPsi(const char *name)
: AliAnalysisTaskSE(name),
fOutput(0),
fHistCent(0),
fEventCounters(0),
fMuonTrackCuts(0),
fTrackCutsSet(kFALSE),
fSelectTrgSign(kFALSE),
fSelectSameTrgSignFake(kFALSE),
fDeltaDev(0),
fUseMCLabel(kFALSE),
//fnCent(10),
fnCent(1),
fRecordEvWithTrgIssues(kFALSE),
fTrgClassMissTrgL0Ev(0x0),
fTrgL0MissTrgClassEv(0x0),
fTrgClassMissTrgOffEv(0x0),
fTrgOffMissTrgClassEv(0x0),
fTrgOffMissTrgL0Ev(0x0),
fTrgL0MissTrgOffEv(0x0),
fOSTrkOSTrgMLLOnlyEv(0x0),
fOSTrkLSTrgMLLOnlyEv(0x0),
fOSTrkLSTrgMULEv(0x0),
fOSTrkOSTrgFakeMLLOnlyEv(0x0),
fOSTrkLSTrgFakeMLLOnlyEv(0x0),
fOSTrkLSTrgFakeMULEv(0x0)
{
  Float_t ptLowEdge[8] = {0., 1., 2., 3., 4., 5., 6., 7.};
  Float_t ptUpEdge[8] = {1., 2., 3., 4., 5., 6., 7., 8.};
  fPtBin[0].Set(8,ptLowEdge);
  fPtBin[1].Set(8,ptUpEdge);

  Float_t yLowEdge[3] = {-2.5, -3., -3.5};
  Float_t yUpEdge[3] = {-3., -3.5, -4.};
  fYBin[0].Set(3,yLowEdge);
  fYBin[1].Set(3,yUpEdge);
  
  fCentBinRange[0][0] = -999.;
  fCentBinRange[0][1] = 999.;
  fCentBinName[0] = "any";
  for (Int_t iCent = 1; iCent < 10; ++iCent) {
    fCentBinRange[iCent][0] = 10.*(iCent-1);
    fCentBinRange[iCent][1] = 10.*iCent;
    fCentBinName[iCent] = Form("%d%d", 10*(iCent-1), 10*iCent);
  }
  
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
AliAnalysisTaskJPsi::~AliAnalysisTaskJPsi()
{
  if (!AliAnalysisManager::GetAnalysisManager()->IsProofMode()) delete fOutput;
  delete fMuonTrackCuts;
}

//________________________________________________________________________
void AliAnalysisTaskJPsi::NotifyRun()
{
  // Set run number for cuts
  if (fTrackCutsSet) return;
  if (!fMuonTrackCuts) AliFatal("No AliMuonTrackCuts");
  fMuonTrackCuts->SetRun(fInputHandler);
  fTrackCutsSet = kTRUE;
}

//________________________________________________________________________
void AliAnalysisTaskJPsi::UserCreateOutputObjects()
{
  // Create histograms
  fOutput = new TList();
  fOutput->SetOwner();
  
  Int_t nbins = 560;  // -> 25 Mev/c2
  TString hName;
  const Int_t npt = fPtBin[0].GetSize();
  const Int_t ny = fYBin[0].GetSize();
  
  //########################################
  //########################################
  // histo w/o sharp pt muon cut from tracker
  //########################################
  //########################################
  
  //########################################
  // histo vs Centrality  & pt bins,  2.5<y<4
  //########################################
  
  // classic dimuon histo
  TH1F *hDimuPM_pt_cent[npt][fnCent];
  TH1F *hDimuPP_pt_cent[npt][fnCent];
  TH1F *hDimuMM_pt_cent[npt][fnCent];
/*  // Mix dimuon histo
  TH1F *hDimuPM_Mix_pt_cent[npt][fnCent];
  TH1F *hDimuPP_Mix_pt_cent[npt][fnCent];
  TH1F *hDimuMM_Mix_pt_cent[npt][fnCent];
*/
  //########################################
  // histo vs Centrality & y bins,  0<pt<8
  //########################################
  
  // classic dimuon histo
  TH1F *hDimuPM_y_cent[ny][fnCent];
  TH1F *hDimuPP_y_cent[ny][fnCent];
  TH1F *hDimuMM_y_cent[ny][fnCent];
/*  // Mix dimuon histo
  TH1F *hDimuPM_Mix_y_cent[ny][fnCent];
  TH1F *hDimuPP_Mix_y_cent[ny][fnCent];
  TH1F *hDimuMM_Mix_y_cent[ny][fnCent];
*/
  //########################################
  // histo vs Centrality & y bins,  0.3<pt<8
  //########################################
  
  // classic dimuon histo
  TH1F *hDimuPM_y_ptphoto_cent[ny][fnCent];
  TH1F *hDimuPP_y_ptphoto_cent[ny][fnCent];
  TH1F *hDimuMM_y_ptphoto_cent[ny][fnCent];
/*  // Mix dimuon histo
  TH1F *hDimuPM_Mix_y_ptphoto_cent[ny][fnCent];
  TH1F *hDimuPP_Mix_y_ptphoto_cent[ny][fnCent];
  TH1F *hDimuMM_Mix_y_ptphoto_cent[ny][fnCent];
*/
  for (Int_t i = 0; i < fnCent; i++)
  {
    //########################################
    // histo vs Centrality  & pt bins,  2.5<y<4
    //########################################
    
    for (Int_t j = 0; j < npt; j++)
    {
      // Classic dimuon analysis (ref)
      hName = Form("hDimuPM_pt_%i_%s", j, fCentBinName[i].Data());
      hDimuPM_pt_cent[j][i] = new TH1F(hName.Data(),hName.Data(),nbins,0.,14.);
      fOutput->Add(hDimuPM_pt_cent[j][i]);
      hName = Form("hDimuPP_pt_%i_%s", j, fCentBinName[i].Data());
      hDimuPP_pt_cent[j][i] = new TH1F(hName.Data(),hName.Data(),nbins,0.,14.);
      fOutput->Add(hDimuPP_pt_cent[j][i]);
      hName = Form("hDimuMM_pt_%i_%s", j, fCentBinName[i].Data());
      hDimuMM_pt_cent[j][i] = new TH1F(hName.Data(),hName.Data(),nbins,0.,14.);
      fOutput->Add(hDimuMM_pt_cent[j][i]);
/*      // Mix dimuon analysis
      hName = Form("hDimuPM_pt_%i_%s_Mix", j,fCentBinName[i].Data());
      hDimuPM_Mix_pt_cent[j][i] = new TH1F(hName.Data(),hName.Data(),nbins,0.,14.);
      fOutput->Add(hDimuPM_Mix_pt_cent[j][i]);
      hName = Form("hDimuPP_pt_%i_%s_Mix", j,fCentBinName[i].Data());
      hDimuPP_Mix_pt_cent[j][i] = new TH1F(hName.Data(),hName.Data(),nbins,0.,14.);
      fOutput->Add(hDimuPP_Mix_pt_cent[j][i]);
      hName = Form("hDimuMM_pt_%i_%s_Mix", j,fCentBinName[i].Data());
      hDimuMM_Mix_pt_cent[j][i] = new TH1F(hName.Data(),hName.Data(),nbins,0.,14.);
      fOutput->Add(hDimuMM_Mix_pt_cent[j][i]);
*/
      if (j < ny) {
        //########################################
        // histo vs Centrality & y bins,  0<pt<8
        //########################################
        
        // Classic dimuon analysis (ref)
        hName = Form("hDimuPM_y_%i_%s", j, fCentBinName[i].Data());
        hDimuPM_y_cent[j][i] = new TH1F(hName.Data(),hName.Data(),nbins,0.,14.);
        fOutput->Add(hDimuPM_y_cent[j][i]);
        hName = Form("hDimuPP_y_%i_%s", j, fCentBinName[i].Data());
        hDimuPP_y_cent[j][i] = new TH1F(hName.Data(),hName.Data(),nbins,0.,14.);
        fOutput->Add(hDimuPP_y_cent[j][i]);
        hName = Form("hDimuMM_y_%i_%s", j, fCentBinName[i].Data());
        hDimuMM_y_cent[j][i] = new TH1F(hName.Data(),hName.Data(),nbins,0.,14.);
        fOutput->Add(hDimuMM_y_cent[j][i]);
/*        // Mix dimuon analysis
        hName = Form("hDimuPM_y_%i_%s_Mix", j,fCentBinName[i].Data());
        hDimuPM_Mix_y_cent[j][i] = new TH1F(hName.Data(),hName.Data(),nbins,0.,14.);
        fOutput->Add(hDimuPM_Mix_y_cent[j][i]);
        hName = Form("hDimuPP_y_%i_%s_Mix", j,fCentBinName[i].Data());
        hDimuPP_Mix_y_cent[j][i] = new TH1F(hName.Data(),hName.Data(),nbins,0.,14.);
        fOutput->Add(hDimuPP_Mix_y_cent[j][i]);
        hName = Form("hDimuMM_y_%i_%s_Mix", j,fCentBinName[i].Data());
        hDimuMM_Mix_y_cent[j][i] = new TH1F(hName.Data(),hName.Data(),nbins,0.,14.);
        fOutput->Add(hDimuMM_Mix_y_cent[j][i]);
*/
        //########################################
        // histo vs Centrality & y bins,  0.3<pt<8
        //########################################
        
        // Classic dimuon analysis (ref)
        hName = Form("hDimuPM_ypt_%i_%s", j, fCentBinName[i].Data());
        hDimuPM_y_ptphoto_cent[j][i] = new TH1F(hName.Data(),hName.Data(),nbins,0.,14.);
        fOutput->Add(hDimuPM_y_ptphoto_cent[j][i]);
        hName = Form("hDimuPP_ypt_%i_%s", j, fCentBinName[i].Data());
        hDimuPP_y_ptphoto_cent[j][i] = new TH1F(hName.Data(),hName.Data(),nbins,0.,14.);
        fOutput->Add(hDimuPP_y_ptphoto_cent[j][i]);
        hName = Form("hDimuMM_ypt_%i_%s", j, fCentBinName[i].Data());
        hDimuMM_y_ptphoto_cent[j][i] = new TH1F(hName.Data(),hName.Data(),nbins,0.,14.);
        fOutput->Add(hDimuMM_y_ptphoto_cent[j][i]);
/*        // Mix dimuon analysis
        hName = Form("hDimuPM_ypt_%i_%s_Mix", j,fCentBinName[i].Data());
        hDimuPM_Mix_y_ptphoto_cent[j][i] = new TH1F(hName.Data(),hName.Data(),nbins,0.,14.);
        fOutput->Add(hDimuPM_Mix_y_ptphoto_cent[j][i]);
        hName = Form("hDimuPP_ypt_%i_%s_Mix", j,fCentBinName[i].Data());
        hDimuPP_Mix_y_ptphoto_cent[j][i] = new TH1F(hName.Data(),hName.Data(),nbins,0.,14.);
        fOutput->Add(hDimuPP_Mix_y_ptphoto_cent[j][i]);
        hName = Form("hDimuMM_ypt_%i_%s_Mix", j,fCentBinName[i].Data());
        hDimuMM_Mix_y_ptphoto_cent[j][i] = new TH1F(hName.Data(),hName.Data(),nbins,0.,14.);
        fOutput->Add(hDimuMM_Mix_y_ptphoto_cent[j][i]);
*/
      }
      
    }
    
  }
  
  //########################################
  //########################################
  // histo w/ sharp pt muon cut from tracker
  //########################################
  //########################################
  
  //########################################
  // histo vs Centrality  & pt bins,  2.5<y<4
  //########################################
  
  // classic dimuon histo
  TH1F *hDimuPM_ptMuonCut_cent[npt][fnCent];
  TH1F *hDimuPP_ptMuonCut_cent[npt][fnCent];
  TH1F *hDimuMM_ptMuonCut_cent[npt][fnCent];
/*  // Mix dimuon histo
  TH1F *hDimuPM_Mix_ptMuonCut_cent[npt][fnCent];
  TH1F *hDimuPP_Mix_ptMuonCut_cent[npt][fnCent];
  TH1F *hDimuMM_Mix_ptMuonCut_cent[npt][fnCent];
*/
  //########################################
  // histo vs Centrality & y bins,  0<pt<8
  //########################################
  
  // classic dimuon histo
  TH1F *hDimuPM_yMuonCut_cent[ny][fnCent];
  TH1F *hDimuPP_yMuonCut_cent[ny][fnCent];
  TH1F *hDimuMM_yMuonCut_cent[ny][fnCent];
/*  // Mix dimuon histo
  TH1F *hDimuPM_Mix_yMuonCut_cent[ny][fnCent];
  TH1F *hDimuPP_Mix_yMuonCut_cent[ny][fnCent];
  TH1F *hDimuMM_Mix_yMuonCut_cent[ny][fnCent];
*/
  //########################################
  // histo vs Centrality & y bins,  0.3<pt<8
  //########################################
  
  // classic dimuon histo
  TH1F *hDimuPM_yMuonCut_ptphoto_cent[ny][fnCent];
  TH1F *hDimuPP_yMuonCut_ptphoto_cent[ny][fnCent];
  TH1F *hDimuMM_yMuonCut_ptphoto_cent[ny][fnCent];
/*  // Mix dimuon histo
  TH1F *hDimuPM_Mix_yMuonCut_ptphoto_cent[ny][fnCent];
  TH1F *hDimuPP_Mix_yMuonCut_ptphoto_cent[ny][fnCent];
  TH1F *hDimuMM_Mix_yMuonCut_ptphoto_cent[ny][fnCent];
*/
  for (Int_t i = 0; i < fnCent; i++)
  {
    //########################################
    // histo vs Centrality  & pt bins,  2.5<y<4
    //########################################
    
    for (Int_t j = 0; j < npt; j++)
    {
      // Classic dimuon analysis (ref)
      hName = Form("hDimuPM_ptMuonCut_%i_%s", j, fCentBinName[i].Data());
      hDimuPM_ptMuonCut_cent[j][i] = new TH1F(hName.Data(),hName.Data(),nbins,0.,14.);
      fOutput->Add(hDimuPM_ptMuonCut_cent[j][i]);
      hName = Form("hDimuPP_ptMuonCut_%i_%s", j, fCentBinName[i].Data());
      hDimuPP_ptMuonCut_cent[j][i] = new TH1F(hName.Data(),hName.Data(),nbins,0.,14.);
      fOutput->Add(hDimuPP_ptMuonCut_cent[j][i]);
      hName = Form("hDimuMM_ptMuonCut_%i_%s", j, fCentBinName[i].Data());
      hDimuMM_ptMuonCut_cent[j][i] = new TH1F(hName.Data(),hName.Data(),nbins,0.,14.);
      fOutput->Add(hDimuMM_ptMuonCut_cent[j][i]);
/*      // Mix dimuon analysis
      hName = Form("hDimuPM_ptMuonCut_%i_%s_Mix", j,fCentBinName[i].Data());
      hDimuPM_Mix_ptMuonCut_cent[j][i] = new TH1F(hName.Data(),hName.Data(),nbins,0.,14.);
      fOutput->Add(hDimuPM_Mix_ptMuonCut_cent[j][i]);
      hName = Form("hDimuPP_ptMuonCut_%i_%s_Mix", j,fCentBinName[i].Data());
      hDimuPP_Mix_ptMuonCut_cent[j][i] = new TH1F(hName.Data(),hName.Data(),nbins,0.,14.);
      fOutput->Add(hDimuPP_Mix_ptMuonCut_cent[j][i]);
      hName = Form("hDimuMM_ptMuonCut_%i_%s_Mix", j,fCentBinName[i].Data());
      hDimuMM_Mix_ptMuonCut_cent[j][i] = new TH1F(hName.Data(),hName.Data(),nbins,0.,14.);
      fOutput->Add(hDimuMM_Mix_ptMuonCut_cent[j][i]);
*/
      if (j < ny) {
        //########################################
        // histo vs Centrality & y bins,  0<pt<8
        //########################################
        
        // Classic dimuon analysis (ref)
        hName = Form("hDimuPM_yMuonCut_%i_%s", j, fCentBinName[i].Data());
        hDimuPM_yMuonCut_cent[j][i] = new TH1F(hName.Data(),hName.Data(),nbins,0.,14.);
        fOutput->Add(hDimuPM_yMuonCut_cent[j][i]);
        hName = Form("hDimuPP_yMuonCut_%i_%s", j, fCentBinName[i].Data());
        hDimuPP_yMuonCut_cent[j][i] = new TH1F(hName.Data(),hName.Data(),nbins,0.,14.);
        fOutput->Add(hDimuPP_yMuonCut_cent[j][i]);
        hName = Form("hDimuMM_yMuonCut_%i_%s", j, fCentBinName[i].Data());
        hDimuMM_yMuonCut_cent[j][i] = new TH1F(hName.Data(),hName.Data(),nbins,0.,14.);
        fOutput->Add(hDimuMM_yMuonCut_cent[j][i]);
/*        // Mix dimuon analysis
        hName = Form("hDimuPM_yMuonCut_%i_%s_Mix", j,fCentBinName[i].Data());
        hDimuPM_Mix_yMuonCut_cent[j][i] = new TH1F(hName.Data(),hName.Data(),nbins,0.,14.);
        fOutput->Add(hDimuPM_Mix_yMuonCut_cent[j][i]);
        hName = Form("hDimuPP_yMuonCut_%i_%s_Mix", j,fCentBinName[i].Data());
        hDimuPP_Mix_yMuonCut_cent[j][i] = new TH1F(hName.Data(),hName.Data(),nbins,0.,14.);
        fOutput->Add(hDimuPP_Mix_yMuonCut_cent[j][i]);
        hName = Form("hDimuMM_yMuonCut_%i_%s_Mix", j,fCentBinName[i].Data());
        hDimuMM_Mix_yMuonCut_cent[j][i] = new TH1F(hName.Data(),hName.Data(),nbins,0.,14.);
        fOutput->Add(hDimuMM_Mix_yMuonCut_cent[j][i]);
*/
        //########################################
        // histo vs Centrality & y bins,  0.3<pt<8
        //########################################
        
        // Classic dimuon analysis (ref)
        hName = Form("hDimuPM_yptMuonCut_%i_%s", j, fCentBinName[i].Data());
        hDimuPM_yMuonCut_ptphoto_cent[j][i] = new TH1F(hName.Data(),hName.Data(),nbins,0.,14.);
        fOutput->Add(hDimuPM_yMuonCut_ptphoto_cent[j][i]);
        hName = Form("hDimuPP_yptMuonCut_%i_%s", j, fCentBinName[i].Data());
        hDimuPP_yMuonCut_ptphoto_cent[j][i] = new TH1F(hName.Data(),hName.Data(),nbins,0.,14.);
        fOutput->Add(hDimuPP_yMuonCut_ptphoto_cent[j][i]);
        hName = Form("hDimuMM_yptMuonCut_%i_%s", j, fCentBinName[i].Data());
        hDimuMM_yMuonCut_ptphoto_cent[j][i] = new TH1F(hName.Data(),hName.Data(),nbins,0.,14.);
        fOutput->Add(hDimuMM_yMuonCut_ptphoto_cent[j][i]);
/*        // Mix dimuon analysis
        hName = Form("hDimuPM_yptMuonCut_%i_%s_Mix", j,fCentBinName[i].Data());
        hDimuPM_Mix_yMuonCut_ptphoto_cent[j][i] = new TH1F(hName.Data(),hName.Data(),nbins,0.,14.);
        fOutput->Add(hDimuPM_Mix_yMuonCut_ptphoto_cent[j][i]);
        hName = Form("hDimuPP_yptMuonCut_%i_%s_Mix", j,fCentBinName[i].Data());
        hDimuPP_Mix_yMuonCut_ptphoto_cent[j][i] = new TH1F(hName.Data(),hName.Data(),nbins,0.,14.);
        fOutput->Add(hDimuPP_Mix_yMuonCut_ptphoto_cent[j][i]);
        hName = Form("hDimuMM_yptMuonCut_%i_%s_Mix", j,fCentBinName[i].Data());
        hDimuMM_Mix_yMuonCut_ptphoto_cent[j][i] = new TH1F(hName.Data(),hName.Data(),nbins,0.,14.);
        fOutput->Add(hDimuMM_Mix_yMuonCut_ptphoto_cent[j][i]);
*/
      }
      
    }
    
  }
  
  // event counters
  fEventCounters = new AliCounterCollection("eventCounters");
  fEventCounters->AddRubric("run", 1000);
  fEventCounters->AddRubric("trigger", "MUL/MLL/MUL&MLL/other");
  fEventCounters->AddRubric("input", "0MUL/0MLL/0MUL&0MLL/other");
  fEventCounters->AddRubric("offline", "OS/LS/OS&LS/other");
  fEventCounters->AddRubric("offline2", "OS/LS/OS&LS/other");
  TString centKeys = "any";
  for (Int_t iCent = 1; iCent < fnCent; ++iCent) centKeys += Form("/%s", fCentBinName[iCent].Data());
  fEventCounters->AddRubric("centrality", centKeys.Data());
  fEventCounters->Init();
  fOutput->Add(fEventCounters);
  
  // Create histograms for check
  fHistCent = new TH1F("fHistCent", "Cent of events mixed", 100, 0, 100);
  fOutput->Add(fHistCent);
  
  if (fRecordEvWithTrgIssues) {
    
    // List of ESD files & events with l0 trigger input missing according to trigger class
    fTrgClassMissTrgL0Ev = new MyList("trgClassMissTrgL0Ev");
    fOutput->Add(fTrgClassMissTrgL0Ev);
    
    // List of ESD files & events with trigger class missing according to l0 trigger input
    fTrgL0MissTrgClassEv = new MyList("trgL0MissTrgClassEv");
    fOutput->Add(fTrgL0MissTrgClassEv);
    
    // List of ESD files & events with trigger track pair deviation sign missing according to trigger class
    fTrgClassMissTrgOffEv = new MyList("trgClassMissTrgOffEv");
    fOutput->Add(fTrgClassMissTrgOffEv);
    
    // List of ESD files & events with trigger class missing according to trigger track pair deviation sign
    fTrgOffMissTrgClassEv = new MyList("trgOffMissTrgClassEv");
    fOutput->Add(fTrgOffMissTrgClassEv);
    
    // List of ESD files & events with l0 trigger input missing according to trigger track pair deviation sign
    fTrgOffMissTrgL0Ev = new MyList("trgOffMissTrgL0Ev");
    fOutput->Add(fTrgOffMissTrgL0Ev);
    
    // List of ESD files & events with trigger track pair deviation sign missing according to l0 trigger input
    fTrgL0MissTrgOffEv = new MyList("trgL0MissTrgOffEv");
    fOutput->Add(fTrgL0MissTrgOffEv);
    
    // List of ESD files & events with selected OS muon pair(s) matching OS trigger track pair in MLL&!MUL class
    fOSTrkOSTrgMLLOnlyEv = new MyList("OSTrkOSTrgMLLOnlyEv");
    fOutput->Add(fOSTrkOSTrgMLLOnlyEv);
    
    // List of ESD files & events with selected OS muon pair(s) matching LS trigger track pair in MLL&!MUL class
    fOSTrkLSTrgMLLOnlyEv = new MyList("OSTrkLSTrgMLLOnlyEv");
    fOutput->Add(fOSTrkLSTrgMLLOnlyEv);
    
    // List of ESD files & events with selected OS muon pair(s) matching LS trigger track pair in MUL class
    fOSTrkLSTrgMULEv = new MyList("OSTrkLSTrgMULEv");
    fOutput->Add(fOSTrkLSTrgMULEv);
    
    // List of ESD files & events with selected OS muon pair(s) matching the same trigger track (->?S) in MLL&!MUL class
    fOSTrkOSTrgFakeMLLOnlyEv = new MyList("OSTrkOSTrgFakeMLLOnlyEv");
    fOutput->Add(fOSTrkOSTrgFakeMLLOnlyEv);
    
    // List of ESD files & events with selected OS muon pair(s) matching the same trigger track (->LS) in MLL&!MUL class
    fOSTrkLSTrgFakeMLLOnlyEv = new MyList("OSTrkLSTrgFakeMLLOnlyEv");
    fOutput->Add(fOSTrkLSTrgFakeMLLOnlyEv);
    
    // List of ESD files & events with selected OS muon pair(s) matching the same trigger track (->LS) in MUL class
    fOSTrkLSTrgFakeMULEv = new MyList("OSTrkLSTrgFakeMULEv");
    fOutput->Add(fOSTrkLSTrgFakeMULEv);
    
  }
  
  PostData(1, fOutput); // Post data for ALL output slots >0 here, to get at least an empty histogram
}

//________________________________________________________________________
void AliAnalysisTaskJPsi::UserExec(Option_t *)
{
  AliVEvent *evt = InputEvent();
  Bool_t isESD = kFALSE;
  if (dynamic_cast<AliESDEvent*>(evt)) isESD = kTRUE;
  else if (!dynamic_cast<AliAODEvent*>(evt)) return;
  
  // select physics events
  if (evt->GetEventType() != 7) return;
  /*
  if (!isESD) {
    if (!static_cast<AliAODEvent*>(evt)->GetHeader()->GetESDFileName().Contains("11000169099005.154")) return;
    //printf("file %s\n",static_cast<AliAODEvent*>(evt)->GetHeader()->GetESDFileName().Data());
    //printf("event %d\n",static_cast<AliAODEvent*>(evt)->GetHeader()->GetEventNumberESDFile());
  }
  */
  // Init
  Bool_t ptMuonCutBoth = kFALSE;
  const Int_t npt = fPtBin[0].GetSize();
  const Int_t ny = fYBin[0].GetSize();
  
  TString esdFile = isESD ? CurrentFileName() : static_cast<AliAODHeader*>(static_cast<AliAODEvent*>(evt)->GetHeader())->GetESDFileName().Data();
  Int_t esdEvent = isESD ? static_cast<AliESDEvent*>(evt)->GetEventNumberInFile() : static_cast<AliAODHeader*>(static_cast<AliAODEvent*>(evt)->GetHeader())->GetEventNumberESDFile();
  
  // Get Centrality
  Float_t cent = evt->GetCentrality()->GetCentralityPercentile("V0M");
  fHistCent->Fill(cent);
  
  // Get trigger class
  TString trigger = evt->GetFiredTriggerClasses();
//  Bool_t trigMUL = trigger.Contains("CPBI1MUL-B-NOPF-MUON");
//  Bool_t trigMLL = trigger.Contains("CPBI1MLL-B-NOPF-MUON");
  Bool_t trigMUL = trigger.Contains("CMUU7-B-NOPF-");
  Bool_t trigMLL = trigger.Contains("CMUL7-B-NOPF-");
  TString trigKey = "trigger:";
  if (trigMUL && trigMLL) trigKey += "MUL&MLL";
  else if (trigMUL) trigKey += "MUL";
  else if (trigMLL) trigKey += "MLL";
  else trigKey += "other";
  
  // Get trigger input
  UInt_t l0input = AliAnalysisMuonUtility::GetL0TriggerInputs(evt);
  Bool_t inputMUL = TESTBIT(l0input,5);
  Bool_t inputMLL = TESTBIT(l0input,7);
  TString inputKey = "input:";
  if (inputMUL && inputMLL) inputKey += "0MUL&0MLL";
  else if (inputMUL) inputKey += "0MUL";
  else if (inputMLL) inputKey += "0MLL";
  else inputKey += "other";
  
  // trigger class from trigger tracks
  Bool_t offlineMUL = kFALSE;
  Bool_t offlineMLL = kFALSE;
  Bool_t offlineMUL2 = kFALSE;
  Bool_t offlineMLL2 = kFALSE;
  Int_t ntracks = AliAnalysisMuonUtility::GetNTracks(evt);
  Int_t *loCircuit = new Int_t[ntracks];
  Int_t nTrgTracks = 0;
  if (isESD) for (Int_t i = 0; i < ntracks; ++i) {
    
    AliVParticle *tracki = AliAnalysisMuonUtility::GetTrack(i, evt);
    if (!AliAnalysisMuonUtility::MatchLpt(tracki)) continue;
    
    // make sure this track is not already accounted for
    Bool_t trackExist = kFALSE;
    Int_t loCircuiti = AliAnalysisMuonUtility::GetLoCircuit(tracki);
    for (Int_t k = 0; k < nTrgTracks; ++k) {
      if (loCircuiti == loCircuit[k]) {
        trackExist = kTRUE;
        break;
      }
    }
    if (trackExist) continue;
    loCircuit[nTrgTracks++] = loCircuiti;
    
    Int_t trgDevSigni = TriggerDevSign(tracki);
    Int_t trgDevSigni2 = TriggerDevSign(tracki, fDeltaDev);
    
    for (Int_t j = i+1; j < ntracks; ++j) {
      
      AliVParticle *trackj = AliAnalysisMuonUtility::GetTrack(j, evt);
      if (!AliAnalysisMuonUtility::MatchLpt(trackj)) continue;
      
      // make sure this track is not already accounted for
      trackExist = kFALSE;
      Int_t loCircuitj = AliAnalysisMuonUtility::GetLoCircuit(trackj);
      for (Int_t k = 0; k < nTrgTracks; ++k) {
        if (loCircuitj == loCircuit[k]) {
          trackExist = kTRUE;
          break;
        }
      }
      if (trackExist) continue;
      
      Int_t trgDevSignj = TriggerDevSign(trackj);
      Int_t trgDevSignj2 = TriggerDevSign(trackj, fDeltaDev);
      
      Int_t trgSign = trgDevSigni * trgDevSignj;
      if (trgSign <= 0) offlineMUL = kTRUE;
      if (trgSign >= 0) offlineMLL = kTRUE;
      
      Int_t trgSign2 = trgDevSigni2 * trgDevSignj2;
      if (trgSign2 <= 0) offlineMUL2 = kTRUE;
      if (trgSign2 >= 0) offlineMLL2 = kTRUE;
      
    }
  }
  delete[] loCircuit;
  TString offlineKey = "offline:";
  if (offlineMUL && offlineMLL) offlineKey += "OS&LS";
  else if (offlineMUL) offlineKey += "OS";
  else if (offlineMLL) offlineKey += "LS";
  else offlineKey += "other";
  TString offlineKey2 = "offline2:";
  if (offlineMUL2 && offlineMLL2) offlineKey2 += "OS&LS";
  else if (offlineMUL2) offlineKey2 += "OS";
  else if (offlineMLL2) offlineKey2 += "LS";
  else offlineKey2 += "other";
  
  // Count
  for (Int_t iCent = 0; iCent < fnCent; iCent++) {
    if (cent > fCentBinRange[iCent][0] && cent <= fCentBinRange[iCent][1]) {
      fEventCounters->Count(Form("run:%d/%s/%s/%s/%s/centrality:%s", fCurrentRunNumber, trigKey.Data(), inputKey.Data(), offlineKey.Data(), offlineKey2.Data(), fCentBinName[iCent].Data()));
      if (iCent > 0) break;
    }
  }
  
  if (fRecordEvWithTrgIssues) {
    
    // find l0 trigger input missing according to trigger class
    if ((trigMUL && !inputMUL) || (trigMLL && !inputMLL)) {
      
      // record faulty esd file/event
      TObjString *badESDFile = FindFile(esdFile, fTrgClassMissTrgL0Ev);
      if (badESDFile) badESDFile->SetString(Form("%s,%d",badESDFile->GetName(),esdEvent));
      else fTrgClassMissTrgL0Ev->Add(new TObjString(Form("%s@%d",esdFile.Data(),esdEvent)));
      /*
       // print event details
       if (!isESD) printf("ESD file: %s | event: %d\n", esdFile.Data(), esdEvent);
       printf("event %d: type = %d | trigger class = %s | trigger input l0 = ", fEntry, evt->GetEventType(), trigger.Data());
       for (Int_t i = 31; i >= 0; --i) printf("%d",TESTBIT(l0input,i) ? 1 : 0);
       UInt_t l1input = AliAnalysisMuonUtility::GetL1TriggerInputs(evt);
       printf(" | l1 = ");
       for (Int_t i = 31; i >= 0; --i) printf("%d",TESTBIT(l1input,i) ? 1 : 0);
       UInt_t l2input = AliAnalysisMuonUtility::GetL2TriggerInputs(evt);
       printf(" | l2 = ");
       for (Int_t i = 31; i >= 0; --i) printf("%d",TESTBIT(l2input,i) ? 1 : 0);
       printf("\n");
       if (isESD) printf("trigger input = %s\n",static_cast<AliESDEvent*>(evt)->GetHeader()->GetFiredTriggerInputs().Data());
       */
    }
    
    // find trigger class missing according to l0 trigger input
    if ((trigMUL || trigMLL) && ((inputMUL && !trigMUL) || (inputMLL && !trigMLL))) {
      TObjString *badESDFile = FindFile(esdFile, fTrgL0MissTrgClassEv);
      if (badESDFile) badESDFile->SetString(Form("%s,%d",badESDFile->GetName(),esdEvent));
      else fTrgL0MissTrgClassEv->Add(new TObjString(Form("%s@%d",esdFile.Data(),esdEvent)));
    }
    
    // find trigger track pair deviation sign missing according to trigger class
    if (isESD && ((trigMUL && !offlineMUL) || (trigMLL && !offlineMLL))) {
      
      // record faulty esd file/event
      TObjString *badESDFile = FindFile(esdFile, fTrgClassMissTrgOffEv);
      if (badESDFile) badESDFile->SetString(Form("%s,%d",badESDFile->GetName(),esdEvent));
      else fTrgClassMissTrgOffEv->Add(new TObjString(Form("%s@%d",esdFile.Data(),esdEvent)));
      /*
       // print event details
       printf("event %d: type = %d | trigger class = %s | trigger sign = ", fEntry, evt->GetEventType(), trigger.Data());
       if (offlineMUL && offlineMLL) printf("Opposite Sign & Like Sign\n");
       else if (offlineMUL) printf("Opposite Sign\n");
       else if (offlineMLL) printf("Like Sign\n");
       else printf("none\n");
       */
    }
    
    // find trigger class missing according to trigger track pair deviation sign
    if (isESD && (trigMUL || trigMLL) && ((offlineMUL && !trigMUL) || (offlineMLL && !trigMLL))) {
      TObjString *badESDFile = FindFile(esdFile, fTrgOffMissTrgClassEv);
      if (badESDFile) badESDFile->SetString(Form("%s,%d",badESDFile->GetName(),esdEvent));
      else fTrgOffMissTrgClassEv->Add(new TObjString(Form("%s@%d",esdFile.Data(),esdEvent)));
    }
    
    // find l0 trigger input missing according to trigger track pair deviation sign
    if (isESD && ((offlineMUL && !inputMUL) || (offlineMLL && !inputMLL))) {
      TObjString *badESDFile = FindFile(esdFile, fTrgOffMissTrgL0Ev);
      if (badESDFile) badESDFile->SetString(Form("%s,%d",badESDFile->GetName(),esdEvent));
      else fTrgOffMissTrgL0Ev->Add(new TObjString(Form("%s@%d",esdFile.Data(),esdEvent)));
    }
    
    // find trigger track pair deviation sign missing according to l0 trigger input
    if (isESD && (trigMUL || trigMLL || offlineMUL || offlineMLL) && ((inputMUL && !offlineMUL) || (inputMLL && !offlineMLL))) {
      TObjString *badESDFile = FindFile(esdFile, fTrgL0MissTrgOffEv);
      if (badESDFile) badESDFile->SetString(Form("%s,%d",badESDFile->GetName(),esdEvent));
      else fTrgL0MissTrgOffEv->Add(new TObjString(Form("%s@%d",esdFile.Data(),esdEvent)));
    }
    
  }
  
  // Track i loop
  for (Int_t i = 0; i < ntracks; ++i)
  {
    AliVParticle *tracki = AliAnalysisMuonUtility::GetTrack(i, evt);
    if (!fMuonTrackCuts->IsSelected(tracki)) continue;
    if (fUseMCLabel && tracki->GetLabel() < 0) continue;
    
    TLorentzVector p(tracki->Px(),tracki->Py(),tracki->Pz(),TMath::Sqrt(MuonMass2()+tracki->P()*tracki->P()));
    Int_t loCircuiti = AliAnalysisMuonUtility::GetLoCircuit(tracki);
    Int_t trgDevSigni = TriggerDevSign(tracki, fDeltaDev);
    
    // Track j loop
    for (Int_t j = i+1; j < ntracks; ++j)
    {
      AliVParticle *trackj = AliAnalysisMuonUtility::GetTrack(j, evt);
      if (!fMuonTrackCuts->IsSelected(trackj)) continue;
      if (fUseMCLabel && trackj->GetLabel() < 0) continue;
      
      TLorentzVector pj(trackj->Px(),trackj->Py(),trackj->Pz(),TMath::Sqrt(MuonMass2()+trackj->P()*trackj->P()));
      Int_t loCircuitj = AliAnalysisMuonUtility::GetLoCircuit(trackj);
      Int_t trgDevSignj = TriggerDevSign(trackj, fDeltaDev);
      
      // Tracker sharp pt muon cut both
      if (p.Pt() > 1. && pj.Pt() > 1.) ptMuonCutBoth = kTRUE;
      else ptMuonCutBoth = kFALSE;
      
      // Create dimuon kinematics
      pj += p;
      if (!PairRapidityCut(pj)) continue;
      
      // pt max due to pp reference
      //if (pj.Pt()>8.) continue;
      
      Float_t massDimuon = pj.M();
      Float_t ptDimuon = pj.Pt();
      Float_t yDimuon = pj.Rapidity();
      
      // Charge from trigger (OS = -1, LS = 1 and unknown = 0)
      Int_t trgSign = trgDevSigni * trgDevSignj;
      
      // Charge from tracker
      TString sign;
      if ( tracki->Charge() > 0. && trackj->Charge() > 0. ) sign = "PP";
      else if ( tracki->Charge() < 0. && trackj->Charge() < 0. ) sign = "MM";
      else sign = "PM";
      
//      if (fDeltaDev == 0 && fRecordEvWithTrgIssues && cent > 0. && cent <= 90. && sign == "PM" && massDimuon > 2.8 && massDimuon < 3.3) {
      if (fDeltaDev == 0 && fRecordEvWithTrgIssues && sign == "PM" && massDimuon > 2.8 && massDimuon < 3.3) {
        
        if (loCircuiti != loCircuitj) {
          
          // record selected OS muon pair(s) matching OS trigger track pair in MLL&!MUL class
          if (trgSign <= 0 && trigMLL && !trigMUL) {
            TObjString *badESDFile = FindFile(esdFile, fOSTrkOSTrgMLLOnlyEv);
            if (badESDFile) badESDFile->SetString(Form("%s,%d",badESDFile->GetName(),esdEvent));
            else fOSTrkOSTrgMLLOnlyEv->Add(new TObjString(Form("%s@%d",esdFile.Data(),esdEvent)));
          }
          
          // record selected OS muon pair(s) matching LS trigger track pair in MLL&!MUL class
          if (trgSign > 0 && trigMLL && !trigMUL) {
            TObjString *badESDFile = FindFile(esdFile, fOSTrkLSTrgMLLOnlyEv);
            if (badESDFile) badESDFile->SetString(Form("%s,%d",badESDFile->GetName(),esdEvent));
            else fOSTrkLSTrgMLLOnlyEv->Add(new TObjString(Form("%s@%d",esdFile.Data(),esdEvent)));
          }
          
          // record selected OS muon pair(s) matching LS trigger track pair in MUL class
          if (trgSign > 0 && trigMUL) {
            TObjString *badESDFile = FindFile(esdFile, fOSTrkLSTrgMULEv);
            if (badESDFile) badESDFile->SetString(Form("%s,%d",badESDFile->GetName(),esdEvent));
            else fOSTrkLSTrgMULEv->Add(new TObjString(Form("%s@%d",esdFile.Data(),esdEvent)));
          }
          
        } else {
          
          // record selected OS muon pair(s) matching the same trigger track (->?S) in MLL&!MUL class
          if (trgSign <= 0 && trigMLL && !trigMUL) {
            TObjString *badESDFile = FindFile(esdFile, fOSTrkOSTrgFakeMLLOnlyEv);
            if (badESDFile) badESDFile->SetString(Form("%s,%d",badESDFile->GetName(),esdEvent));
            else fOSTrkOSTrgFakeMLLOnlyEv->Add(new TObjString(Form("%s@%d",esdFile.Data(),esdEvent)));
          }
          
          // record selected OS muon pair(s) matching the same trigger track (->LS) in MLL&!MUL class
          if (trgSign > 0 && trigMLL && !trigMUL) {
            TObjString *badESDFile = FindFile(esdFile, fOSTrkLSTrgFakeMLLOnlyEv);
            if (badESDFile) badESDFile->SetString(Form("%s,%d",badESDFile->GetName(),esdEvent));
            else fOSTrkLSTrgFakeMLLOnlyEv->Add(new TObjString(Form("%s@%d",esdFile.Data(),esdEvent)));
          }
          
          // record selected OS muon pair(s) matching the same trigger track (->LS) in MUL class
          if (trgSign > 0 && trigMUL) {
            TObjString *badESDFile = FindFile(esdFile, fOSTrkLSTrgFakeMULEv);
            if (badESDFile) badESDFile->SetString(Form("%s,%d",badESDFile->GetName(),esdEvent));
            else fOSTrkLSTrgFakeMULEv->Add(new TObjString(Form("%s@%d",esdFile.Data(),esdEvent)));
          }
          
        }
        
      }
      
      // select pairs whose sign is correctly identified by the trigger
      if (fSelectTrgSign && ((sign == "PM" && trgSign > 0) || (sign != "PM" && trgSign < 0))) continue;
      
      // select pairs matching the same trigger track and producing LS & !OS
      if (fSelectSameTrgSignFake && (loCircuiti != loCircuitj || trgSign <= 0)) continue;
      
      // Loop on centrality bins
      for (Int_t iCent = 0; iCent < fnCent; iCent++) {
        if (cent > fCentBinRange[iCent][0] && cent <= fCentBinRange[iCent][1]) {
          
          // Loop on pt bins
          for (Int_t ipt = 0; ipt < npt; ipt++) {
            if (ptDimuon > fPtBin[0].At(ipt) && ptDimuon <= fPtBin[1].At(ipt)) {
              ( (TH1F*) fOutput->FindObject( Form("hDimu%s_pt_%i_%s", sign.Data(), ipt, fCentBinName[iCent].Data()) ) )->Fill(massDimuon);
              if (ptMuonCutBoth) ( (TH1F*) fOutput->FindObject( Form("hDimu%s_ptMuonCut_%i_%s", sign.Data(), ipt, fCentBinName[iCent].Data()) ) )->Fill(massDimuon);
            }
          }
          
          // Loop on y bins
          for (Int_t iy = 0; iy < ny; iy++) {
            if (yDimuon < fYBin[0].At(iy) && yDimuon >= fYBin[1].At(iy)) {
              ( (TH1F*) fOutput->FindObject( Form("hDimu%s_y_%i_%s", sign.Data(), iy, fCentBinName[iCent].Data()) ) )->Fill(massDimuon);
              if (ptMuonCutBoth) ( (TH1F*) fOutput->FindObject( Form("hDimu%s_yMuonCut_%i_%s", sign.Data(), iy, fCentBinName[iCent].Data()) ) )->Fill(massDimuon);
              if (ptDimuon > 0.3 && ptDimuon <= 8.) {
                ( (TH1F*) fOutput->FindObject( Form("hDimu%s_ypt_%i_%s", sign.Data(), iy, fCentBinName[iCent].Data()) ) )->Fill(massDimuon);
                if (ptMuonCutBoth) ( (TH1F*) fOutput->FindObject( Form("hDimu%s_yptMuonCut_%i_%s", sign.Data(), iy, fCentBinName[iCent].Data()) ) )->Fill(massDimuon);
                
              }
            }
          }
          
          if (iCent > 0) break;
          
        }
      } // end loop centrality bins
      
    } // end loop track j
  } // end loop track i
  
  
  PostData(1, fOutput);
  
}

//_____________________________________________________________________________
void AliAnalysisTaskJPsi::Terminate(Option_t *)
{
  fOutput = dynamic_cast<TList*>(GetOutputData(1));
  if (!fOutput) return;
  
  TString extension = GetName();
  extension.Remove(0, extension.Index("_"));
  
  TH1F *hDimu = new TH1F(Form("hDimuPM_cent0-90%s",extension.Data()),"hDimuPM_cent0-90",560,0.,14.);
  TH1F *hDimu58 = new TH1F(Form("hDimuPM_pT5-8_cent50-90%s",extension.Data()),"hDimuPM_pT5-8_cent50-90",560,0.,14.);
  
  // Loop on centrality bins
  for (Int_t iCent = 1; iCent < fnCent; iCent++) {
    
    hDimu->Add((TH1F*)fOutput->FindObject(Form("hDimuPM_y_0_%s", fCentBinName[iCent].Data())));
    
    if (iCent > 5) hDimu58->Add((TH1F*)fOutput->FindObject(Form("hDimuPM_pt_3_%s", fCentBinName[iCent].Data())));
    
  }
  
  TCanvas *c = new TCanvas(Form("cDimu%s",extension.Data()),Form("cDimu%s",extension.Data()),800,400);
  c->Divide(2,1);
  c->cd(1);
  gPad->SetLogy();
  hDimu->Draw();
  c->cd(2);
  gPad->SetLogy();
  hDimu58->Draw();
  
}

//_____________________________________________________________________________
Double_t AliAnalysisTaskJPsi::MuonMass2() const
{
  static Double_t m2 = 1.11636129640000012e-02;
  return m2;
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskJPsi::PairRapidityCut(const TLorentzVector& pair) const
{
  Double_t y = pair.Rapidity();
  return ( y < -2.5 && y > -4.0 );
}

//_____________________________________________________________________________
Int_t AliAnalysisTaskJPsi::TriggerDevSign(AliVParticle *track, UInt_t deltaDev) const
{
  /// get the sign (±1) of track deviation in the trigger (0 = unknown)
  
  AliESDMuonTrack *esdTrack = dynamic_cast<AliESDMuonTrack*>(track);
  if (!esdTrack) return 0;
  
  Int_t deviation = esdTrack->LoDev();
  if (deviation > 15+((Int_t)deltaDev)) return 1;
  else if (deviation < 15-((Int_t)deltaDev)) return -1;
  else return 0;
  
}

//_____________________________________________________________________________
TObjString* AliAnalysisTaskJPsi::FindFile(TString &file, MyList *badFiles) const
{
  /// look for the given file in the list of problematic ones
  
  TIter nextBadFile(badFiles->GetList());
  TObjString *badFile = 0x0;
  while ((badFile = static_cast<TObjString*>(nextBadFile())))
    if (badFile->String().BeginsWith(file)) return badFile;
  return 0x0;
  
}

