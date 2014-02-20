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

#include <TList.h>
#include <TH1.h>
#include <TArrayF.h>

#include "AliLog.h"
#include "AliAnalysisManager.h"
#include "AliAODTrack.h"
#include "AliAODEvent.h"
#include "AliCounterCollection.h"
#include "AliMuonTrackCuts.h"

#include "AliAnalysisTaskEx02.h"
#include <AliMultiInputEventHandler.h>
#include <AliMixInputEventHandler.h>

ClassImp(AliAnalysisTaskEx02)


//________________________________________________________________________
AliAnalysisTaskEx02::AliAnalysisTaskEx02() 
   : AliAnalysisTaskSE(),
     fOutput(0),
     fHistCent(0),
     fEventCounters(0),
     fMuonTrackCuts(0)
{
   // Dummy constructor ALWAYS needed for I/O.
}

//________________________________________________________________________
AliAnalysisTaskEx02::AliAnalysisTaskEx02(const char *name) 
   : AliAnalysisTaskSE(name),
     fOutput(0),
     fHistCent(0),
     fEventCounters(0),
     fMuonTrackCuts(0)
{
  Float_t ptLowEdge[8] = {0., 1., 2., 3., 4., 5., 6., 7.};
  Float_t ptUpEdge[8] = {1., 2., 3., 4., 5., 6., 7., 8.};
  fPtBin[0].Set(8,ptLowEdge);
  fPtBin[1].Set(8,ptUpEdge);
  
  Float_t yLowEdge[3] = {-2.5, -3., -3.5};
  Float_t yUpEdge[3] = {-3., -3.5, -4.};
  fYBin[0].Set(3,yLowEdge);
  fYBin[1].Set(3,yUpEdge);
  
  DefineOutput(1, TList::Class());     
}

//________________________________________________________________________
AliAnalysisTaskEx02::~AliAnalysisTaskEx02()
{
  if (fOutput && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
	delete fOutput;
  }
  delete fMuonTrackCuts;
}

//________________________________________________________________________
void AliAnalysisTaskEx02::NotifyRun()
{
  // Set run number for cuts
  if (fMuonTrackCuts==0x0) AliFatal("No AliMuonTrackCuts");
  fMuonTrackCuts->SetRun(fCurrentRunNumber);
}

//________________________________________________________________________
void AliAnalysisTaskEx02::UserCreateOutputObjects()
{
  // Create histograms  
  fOutput = new TList();
  fOutput->SetOwner(); 
  
  Int_t nbins = 560;  // -> 25 Mev/c2
  TString hName;
  const Int_t nCent=9;
  TString centBinName[nCent] = {"010", "1020", "2030", "3040", "4050", "5060", "6070", "7080", "8090"};
  const Int_t npt = fPtBin[0].GetSize();   
  const Int_t ny = fYBin[0].GetSize(); 
  
  //########################################
  // histo vs Centrality  & pt bins,  2.5<y<4
  //########################################
  
  // classic dimuon histo
  TH1F *hDimuPM_pt_cent[npt][nCent];
  TH1F *hDimuPP_pt_cent[npt][nCent];
  TH1F *hDimuMM_pt_cent[npt][nCent];
  // Mix dimuon histo
  TH1F *hDimuPM_Mix_pt_cent[npt][nCent];
  TH1F *hDimuPP_Mix_pt_cent[npt][nCent];
  TH1F *hDimuMM_Mix_pt_cent[npt][nCent];
 
  //########################################
  // histo vs Centrality & y bins,  0<pt<8  
  //########################################
  
  // classic dimuon histo
  TH1F *hDimuPM_y_cent[ny][nCent];
  TH1F *hDimuPP_y_cent[ny][nCent];
  TH1F *hDimuMM_y_cent[ny][nCent];
  // Mix dimuon histo
  TH1F *hDimuPM_Mix_y_cent[ny][nCent];
  TH1F *hDimuPP_Mix_y_cent[ny][nCent];
  TH1F *hDimuMM_Mix_y_cent[ny][nCent];
  
 
  for (Int_t i = 0; i < nCent; i++) 
  {
	//########################################
	// histo vs Centrality  & pt bins,  2.5<y<4
	//########################################
	
	for (Int_t j = 0; j < npt; j++) 
	{
	  // Classic dimuon analysis (ref)
	  hName = Form("hDimuPM_pt_%i_%s", j, centBinName[i].Data());
	  hDimuPM_pt_cent[j][i] = new TH1F(hName.Data(),hName.Data(),nbins,0.,14.);
	  fOutput->Add(hDimuPM_pt_cent[j][i]);
	  hName = Form("hDimuPP_pt_%i_%s", j, centBinName[i].Data());
	  hDimuPP_pt_cent[j][i] = new TH1F(hName.Data(),hName.Data(),nbins,0.,14.);
	  fOutput->Add(hDimuPP_pt_cent[j][i]);
	  hName = Form("hDimuMM_pt_%i_%s", j, centBinName[i].Data());
	  hDimuMM_pt_cent[j][i] = new TH1F(hName.Data(),hName.Data(),nbins,0.,14.);
	  fOutput->Add(hDimuMM_pt_cent[j][i]);
	  // Mix dimuon analysis
	  hName = Form("hDimuPM_pt_%i_%s_Mix", j,centBinName[i].Data());
	  hDimuPM_Mix_pt_cent[j][i] = new TH1F(hName.Data(),hName.Data(),nbins,0.,14.);
	  fOutput->Add(hDimuPM_Mix_pt_cent[j][i]);
	  hName = Form("hDimuPP_pt_%i_%s_Mix", j,centBinName[i].Data());
	  hDimuPP_Mix_pt_cent[j][i] = new TH1F(hName.Data(),hName.Data(),nbins,0.,14.);
	  fOutput->Add(hDimuPP_Mix_pt_cent[j][i]);
	  hName = Form("hDimuMM_pt_%i_%s_Mix", j,centBinName[i].Data());
	  hDimuMM_Mix_pt_cent[j][i] = new TH1F(hName.Data(),hName.Data(),nbins,0.,14.);
	  fOutput->Add(hDimuMM_Mix_pt_cent[j][i]);
	  
	
	  //########################################
	  // histo vs Centrality & y bins,  0<pt<8  
	  //########################################
	
	  if (j < ny) {
		// Classic dimuon analysis (ref)
		hName = Form("hDimuPM_y_%i_%s", j, centBinName[i].Data());
		hDimuPM_y_cent[j][i] = new TH1F(hName.Data(),hName.Data(),nbins,0.,14.);
		fOutput->Add(hDimuPM_y_cent[j][i]);
		hName = Form("hDimuPP_y_%i_%s", j, centBinName[i].Data());
		hDimuPP_y_cent[j][i] = new TH1F(hName.Data(),hName.Data(),nbins,0.,14.);
		fOutput->Add(hDimuPP_y_cent[j][i]);
		hName = Form("hDimuMM_y_%i_%s", j, centBinName[i].Data());
		hDimuMM_y_cent[j][i] = new TH1F(hName.Data(),hName.Data(),nbins,0.,14.);
		fOutput->Add(hDimuMM_y_cent[j][i]);
		// Mix dimuon analysis
		hName = Form("hDimuPM_y_%i_%s_Mix", j,centBinName[i].Data());
		hDimuPM_Mix_y_cent[j][i] = new TH1F(hName.Data(),hName.Data(),nbins,0.,14.);
		fOutput->Add(hDimuPM_Mix_y_cent[j][i]);
		hName = Form("hDimuPP_y_%i_%s_Mix", j,centBinName[i].Data());
		hDimuPP_Mix_y_cent[j][i] = new TH1F(hName.Data(),hName.Data(),nbins,0.,14.);
		fOutput->Add(hDimuPP_Mix_y_cent[j][i]);
		hName = Form("hDimuMM_y_%i_%s_Mix", j,centBinName[i].Data());
		hDimuMM_Mix_y_cent[j][i] = new TH1F(hName.Data(),hName.Data(),nbins,0.,14.);
		fOutput->Add(hDimuMM_Mix_y_cent[j][i]);
	  }
	}
  }
  
  // event counters
  fEventCounters = new AliCounterCollection("eventCounters");
  fEventCounters->AddRubric("centrality", "010/1020/2030/3040/4050/5060/6070/7080/8090");
  fEventCounters->Init();
  fOutput->Add(fEventCounters);
  
  // Create histograms for check
  fHistCent = new TH1F("fHistCent", "Cent of events mixed", 100, 0, 100);
  fOutput->Add(fHistCent);  
  
  PostData(1, fOutput); // Post data for ALL output slots >0 here, to get at least an empty histogram
}

//________________________________________________________________________
void AliAnalysisTaskEx02::UserExec(Option_t *)
{

   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   AliMultiInputEventHandler *inEvHMainMulti = dynamic_cast<AliMultiInputEventHandler *>(mgr->GetInputEventHandler());
   if (inEvHMainMulti) {
	 AliInputEventHandler *inEvMain = dynamic_cast<AliInputEventHandler *>(inEvHMainMulti->GetFirstInputEventHandler());

     // Get AOD Event
	 AliAODEvent* aodEvent = dynamic_cast<AliAODEvent*>(inEvMain->GetEvent());
	 if (aodEvent) {
	   
	   // Init
	   const Int_t nCent=9;
	   Double_t centBinRange[nCent][2] = {{0., 10.}, {10., 20.}, {20., 30.}, {30., 40.}, {40., 50.}, {50., 60.}, {60., 70.}, {70., 80.}, {80.,90.}};
	   TString centBinName[nCent] = {"010", "1020", "2030", "3040", "4050", "5060", "6070", "7080", "8090"};
	  
	   const Int_t npt = fPtBin[0].GetSize();   
	   const Int_t ny = fYBin[0].GetSize();
	   
	   // Get Centrality
	   AliCentrality* acent = aodEvent->GetCentrality();
	   Float_t cent=0.;
	   if (acent) {
		 cent = acent->GetCentralityPercentile("V0M");
		 fHistCent->Fill(cent);	   
	   }
	   
	   // Count
	   TString trigger = aodEvent->GetFiredTriggerClasses();
	   Bool_t isTrigger = trigger.Contains("CPBI1MUL-B-NOPF-MUON");
	   if(isTrigger) {
		 for (Int_t iCent = 0; iCent < nCent; iCent++) {		   
		   if (cent > centBinRange[iCent][0] && cent <= centBinRange[iCent][1]) {
			 fEventCounters->Count(Form("centrality:%s",centBinName[iCent].Data()));
		   }
		 }
	   }
	   
		  
	   // Track i loop
	   Int_t ntracks = aodEvent->GetNumberOfTracks();
	   for (Int_t i = 0; i < ntracks; ++i)
	   {
		 AliAODTrack* tracki = aodEvent->GetTrack(i);
		 if (!fMuonTrackCuts->IsSelected(tracki)) continue;
		 TLorentzVector p(tracki->Px(),tracki->Py(),tracki->Pz(),TMath::Sqrt(MuonMass2()+tracki->P()*tracki->P()));
		 
		 // Track j loop
		 for (Int_t j = i+1; j < ntracks; ++j)
		 {
		   AliAODTrack* trackj = aodEvent->GetTrack(j);
		   if (!fMuonTrackCuts->IsSelected(trackj)) continue;
		   TLorentzVector pj(trackj->Px(),trackj->Py(),trackj->Pz(),TMath::Sqrt(MuonMass2()+trackj->P()*trackj->P()));
		   
		   // Create dimuon kinematics 
		   pj += p;
		   if (!PairRapidityCut(pj)) continue;
		   
		   // pt max due to pp reference
		   if (pj.Pt()>8.) continue;
		   
		   Float_t massDimuon = pj.M();
		   Float_t ptDimuon = pj.Pt();
		   Float_t yDimuon = pj.Rapidity();
		   
		   // Charge
		   TString sign;
		   if ( tracki->Charge() > 0. && trackj->Charge() > 0. ) sign = "PP";
		   else if ( tracki->Charge() < 0. && trackj->Charge() < 0. ) sign = "MM";
		   else sign = "PM";
		   
		   // Loop on centrality bins
		   for (Int_t iCent = 0; iCent < nCent; iCent++) {		   
			 if (cent > centBinRange[iCent][0] && cent <= centBinRange[iCent][1]) {
			   
			   // Loop on pt bins
			   for (Int_t ipt = 0; ipt < npt; ipt++) {		   
				 if (ptDimuon > fPtBin[0].At(ipt) && ptDimuon <= fPtBin[1].At(ipt)) {
				   ( (TH1F*) fOutput->FindObject( Form("hDimu%s_pt_%i_%s", sign.Data(), ipt, centBinName[iCent].Data()) ) )->Fill(massDimuon);
				 }
			   }
			   
			   // Loop on y bins
			   for (Int_t iy = 0; iy < ny; iy++) {		   
				 if (yDimuon < fYBin[0].At(iy) && yDimuon >= fYBin[1].At(iy)) {
				   ( (TH1F*) fOutput->FindObject( Form("hDimu%s_y_%i_%s", sign.Data(), iy, centBinName[iCent].Data()) ) )->Fill(massDimuon);
				 }
			   }
			   
			   break;
			 
			 }
		   } // end loop centrality bins
		   			   
		 } // end loop track j
	   } // end loop track i
	 }
   }
  
  PostData(1, fOutput);
}

//________________________________________________________________________
void AliAnalysisTaskEx02::UserExecMix(Option_t *)
{
  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliMultiInputEventHandler *inEvHMain = dynamic_cast<AliMultiInputEventHandler *>(mgr->GetInputEventHandler());
  if (inEvHMain) {
	AliMixInputEventHandler *mixEH = dynamic_cast<AliMixInputEventHandler *>(inEvHMain->GetFirstMultiInputHandler());
	if (!mixEH) return;
	if (mixEH->CurrentBinIndex() < 0) {
	  AliDebug(AliLog::kDebug + 1, "Current event mixEH->CurrentEntry() == -1");
	  return ;
	}
	AliDebug(AliLog::kDebug, Form("Mixing %lld %d [%lld,%lld] %d", mixEH->CurrentEntry(), mixEH->CurrentBinIndex(), mixEH->CurrentEntryMain(), mixEH->CurrentEntryMix(), mixEH->NumberMixed()));
	AliInputEventHandler *ihMainCurrent = inEvHMain->GetFirstInputEventHandler();
	
	AliMultiInputEventHandler *inEvHMixedCurrent = mixEH->GetFirstMultiInputHandler();
	AliInputEventHandler *ihMixedCurrent = inEvHMixedCurrent->GetFirstInputEventHandler();
	
	// Get AOD Events
	AliAODEvent *aodEvent = dynamic_cast<AliAODEvent *>(ihMainCurrent->GetEvent());
	if (aodEvent) {
	  AliAODEvent *aodEventMix = dynamic_cast<AliAODEvent *>(ihMixedCurrent->GetEvent());
	  AliDebug(AliLog::kDebug, Form("Multi=%d MultiMix=%d", aodEvent->GetNumberOfTracks(), aodEventMix->GetNumberOfTracks()));
	  
	  // Init
	  const Int_t nCent=9;
	  Double_t centBinRange[nCent][2] = {{0., 10.}, {10., 20.}, {20., 30.}, {30., 40.}, {40., 50.}, {50., 60.}, {60., 70.}, {70., 80.}, {80.,90.}};
	  TString centBinName[nCent] = {"010", "1020", "2030", "3040", "4050", "5060", "6070", "7080", "8090"};
	  
	  const Int_t npt = fPtBin[0].GetSize();   
	  const Int_t ny = fYBin[0].GetSize();
	  
	  // Get Centrality
	  AliCentrality* acent = aodEvent->GetCentrality();
	  Float_t cent=0.;
	  if (acent) cent = acent->GetCentralityPercentile("V0M");
	  
	  
	  // Track i loop of current event
	  Int_t ntracks = aodEvent->GetNumberOfTracks();
	  Int_t ntracksMix = aodEventMix->GetNumberOfTracks();
	  for (Int_t i = 0; i < ntracks; ++i) {
		AliAODTrack* tracki = aodEvent->GetTrack(i);
		if (!fMuonTrackCuts->IsSelected(tracki)) continue;
		TLorentzVector p(tracki->Px(),tracki->Py(),tracki->Pz(),TMath::Sqrt(MuonMass2()+tracki->P()*tracki->P()));
		
		// Track j loop of Mix event
		for (Int_t jMix = 0; jMix < ntracksMix; jMix++)
		{
		  AliAODTrack* trackjMix = aodEventMix->GetTrack(jMix);
		  if (!fMuonTrackCuts->IsSelected(trackjMix)) continue;
		  TLorentzVector pjMix(trackjMix->Px(),trackjMix->Py(),trackjMix->Pz(),TMath::Sqrt(MuonMass2()+trackjMix->P()*trackjMix->P()));
		  
		  // Create dimuon kinematics 
		  pjMix += p;
		  if (!PairRapidityCut(pjMix)) continue;
		  
		  // pt max due to pp reference
		  if (pjMix.Pt()>8.) continue;
		  
		  Float_t massDimuon = pjMix.M();
		  Float_t ptDimuon = pjMix.Pt();
		  Float_t yDimuon = pjMix.Rapidity();
	
		  // Charge
		  TString sign;
		  if ( tracki->Charge() > 0. && trackjMix->Charge() > 0. ) sign = "PP";
		  else if ( tracki->Charge() < 0. && trackjMix->Charge() < 0. ) sign = "MM";
		  else sign = "PM";
		  
		  // Loop on centrality bins
		  for (Int_t iCent = 0; iCent < nCent; iCent++) {		   
			if (cent > centBinRange[iCent][0] && cent <= centBinRange[iCent][1]) {
			  
			  // Loop on pt bins
			  for (Int_t ipt = 0; ipt < npt; ipt++) {		   
				if (ptDimuon > fPtBin[0].At(ipt) && ptDimuon <= fPtBin[1].At(ipt)) {
				  ( (TH1F*) fOutput->FindObject( Form("hDimu%s_pt_%i_%s_Mix", sign.Data(), ipt, centBinName[iCent].Data()) ) )->Fill(massDimuon);
				}
			  }
			  
			  // Loop on y bins
			  for (Int_t iy = 0; iy < ny; iy++) {		   
				if (yDimuon < fYBin[0].At(iy) && yDimuon >= fYBin[1].At(iy)) {
				  ( (TH1F*) fOutput->FindObject( Form("hDimu%s_y_%i_%s_Mix", sign.Data(), iy, centBinName[iCent].Data()) ) )->Fill(massDimuon);
				}
			  }
			  
			  break;
			  
			} 
		  } // end loop centrality bins

		} // end loop track j of Mix event
	  } // end loop track i of Current event 
	}
  }

PostData(1, fOutput);
}

//_____________________________________________________________________________
Double_t AliAnalysisTaskEx02::MuonMass2() const
{
  static Double_t m2 = 1.11636129640000012e-02;
  return m2;
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskEx02::TrackPtCut(const AliAODTrack& track) const
{
  return track.Pt() > 1.;
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskEx02::PairRapidityCut(const TLorentzVector& pair) const
{
  Double_t y = pair.Rapidity();
  Bool_t ok = ( y < -2.5 && y > -4.0 );
  
  return ok;
}

