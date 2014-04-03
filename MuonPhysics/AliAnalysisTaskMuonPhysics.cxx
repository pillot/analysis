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

#include <Riostream.h>

// ROOT includes
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TFormula.h"
#include "TCanvas.h"
#include "TList.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TMath.h"

// STEER includes
#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliVParticle.h"
#include "AliESDMuonTrack.h"
#include "AliAODTrack.h"
#include "AliCentrality.h"
#include "AliESDVertex.h"
#include "AliVVertex.h"
#include "AliMultiplicity.h"
#include "AliAODTracklets.h"

// ANALYSIS includes
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisDataSlot.h"
#include "AliAnalysisDataContainer.h"
#include "AliAnalysisManager.h"
#include "AliCounterCollection.h"
#include "AliAnalysisTaskMuonPhysics.h"
#include "AliAnalysisMuonUtility.h"
/*
// MUON includes
#include "AliMUONTrack.h"
#include "AliMUONESDInterface.h"
*/
ClassImp(AliAnalysisTaskMuonPhysics)

//________________________________________________________________________
AliAnalysisTaskMuonPhysics::AliAnalysisTaskMuonPhysics() :
AliAnalysisTaskSE(), 
fList(0x0),
fTrackCounters(0x0),
fEventCounters(0x0),
fTrigCounters(0x0),
fMuonTrackCuts(0x0),
fMuonTrackCuts2(0x0),
fUseMCLabel(kFALSE),
fSelectBadTracks(kFALSE),
fCentMin(-FLT_MAX),
fCentMax(FLT_MAX),
fVsRun(kFALSE)
{
  // Dummy constructor
}

//________________________________________________________________________
AliAnalysisTaskMuonPhysics::AliAnalysisTaskMuonPhysics(const char *name) :
AliAnalysisTaskSE(name), 
fList(0x0),
fTrackCounters(0x0),
fEventCounters(0x0),
fTrigCounters(0x0),
fMuonTrackCuts(0x0),
fMuonTrackCuts2(0x0),
fUseMCLabel(kFALSE),
fSelectBadTracks(kFALSE),
fCentMin(-FLT_MAX),
fCentMax(FLT_MAX),
fVsRun(kFALSE)
{
  /// Constructor
  
  fBranchNames = "ESD:AliESDRun.,AliESDHeader.,MuonTracks,MuonClusters,MuonPads,SPDVertex.,AliESDTZERO.";
  
  // Output slot #1 writes into a TObjArray container
  DefineOutput(1,TObjArray::Class());
  // Output slot #2 writes track counters
  DefineOutput(2,AliCounterCollection::Class());
  // Output slot #3 writes event counters
  DefineOutput(3,AliCounterCollection::Class());
  // Output slot #4 writes trigger track counters
  DefineOutput(4,AliCounterCollection::Class());
}

//________________________________________________________________________
AliAnalysisTaskMuonPhysics::~AliAnalysisTaskMuonPhysics()
{
  /// Destructor
  if (!AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
    delete fList;
    delete fTrackCounters;
    delete fEventCounters;
    delete fTrigCounters;
  }
  delete fMuonTrackCuts;
  delete fMuonTrackCuts2;
}

//___________________________________________________________________________
void AliAnalysisTaskMuonPhysics::UserCreateOutputObjects()
{
  /// Create histograms and counters
  
  fList = new TObjArray(2000);
  fList->SetOwner();
  
  TH1F* hPt = new TH1F("hPt", "transverse momentum distribution;p_{t} (GeV/c)", 300, 0., 30.);
  fList->AddAtAndExpand(hPt, kPt);
  
  TH1F* hRapidity = new TH1F("hRapidity", "rapidity distribution;rapidity", 200, -4.5, -2.);
  fList->AddAtAndExpand(hRapidity, kRapidity);
  
  TH1F* hDCA = new TH1F("hDCA", "DCA distribution;DCA (cm)", 500, 0., 500.);
  fList->AddAtAndExpand(hDCA, kDCA);
  
  TH1F* hChi2 = new TH1F("hChi2", "normalized #chi^{2} distribution;#chi^{2} / ndf", 500, 0., 50.);
  fList->AddAtAndExpand(hChi2, kChi2);
  
  TH1F* hNClustersPerTrack = new TH1F("hNClustersPerTrack", "number of associated clusters per track;n_{clusters}", 20, 0., 20.);
  fList->AddAtAndExpand(hNClustersPerTrack, kNClustersPerTrack);
  
  TH1F* hNChamberHitPerTrack = new TH1F("hNChamberHitPerTrack", "number of chambers hit per track;n_{chamber hit}", 15, 0., 15.);
  fList->AddAtAndExpand(hNChamberHitPerTrack, kNChamberHitPerTrack);
  
  TH1F* hMass = new TH1F("hMass", "#mu^{+}#mu^{-} invariant mass distribution;mass (GeV/c^{2})", 2000, 0., 20.);
  fList->AddAtAndExpand(hMass, kMass);
  
  TH1F* hPUncorrected = new TH1F("hPUncorrected", "uncorrected momentum distribution;p_{uncorr} (GeV/c)", 300, 0., 300.);
  fList->AddAtAndExpand(hPUncorrected, kPUncorrected);
  
  TH1F* hRAbs = new TH1F("hRAbs", "track position at the end of the absorber;R_{abs} (cm)", 1000, 0., 100.);
  fList->AddAtAndExpand(hRAbs, kRAbs);
  
  TH1F* hDCAX = new TH1F("hDCAX", "DCA distribution in x direction;DCA-X (cm)", 2000, -500., 500.);
  fList->AddAtAndExpand(hDCAX, kDCAX);
  
  TH1F* hDCAY = new TH1F("hDCAY", "DCA distribution in y direction;DCA-Y (cm)", 2000, -500., 500.);
  fList->AddAtAndExpand(hDCAY, kDCAY);
  
  TH1F* hNTracks = new TH1F("hNTracks", "number of tracker tracks;n_{tracks}", 50, 0., 50.);
  fList->AddAtAndExpand(hNTracks, kNTracks);
  
  TH1F* hNTrigAll = new TH1F("hNTrigAll", "number of trigger tracks;n_{tracks}", 50, 0., 50.);
  fList->AddAtAndExpand(hNTrigAll, kNTrigAll);
  
  TH1F* hNTrigLow = new TH1F("hNTrigLow", "number of low-pt trigger tracks;n_{tracks}", 50, 0., 50.);
  fList->AddAtAndExpand(hNTrigLow, kNTrigLow);
  
  TH1F* hNTrigHigh = new TH1F("hNTrigHigh", "number of high-pt trigger tracks;n_{tracks}", 50, 0., 50.);
  fList->AddAtAndExpand(hNTrigHigh, kNTrigHigh);
  
  TH1F* hNMatchTracks = new TH1F("hNMatchTracks", "number of tracks matched with trigger;n_{tracks}", 50, 0., 50.);
  fList->AddAtAndExpand(hNMatchTracks, kNMatchTracks);
  
  TH1F* hChi2Trig = new TH1F("hChi2Trig", "matching #chi^{2} distribution;#chi^{2} / ndf", 500, 0., 50.);
  fList->AddAtAndExpand(hChi2Trig, kChi2Trig);
  
  TH1F* hPDCA23 = new TH1F("hPDCA23", "p #times DCA distribution in 2#circ < #theta_{abs} < 3#circ;p #times DCA (GeV.cm/c)", 2500, 0., 5000.);
  fList->AddAtAndExpand(hPDCA23, kPDCA23);
  
  TH1F* hPDCA310 = new TH1F("hPDCA310", "p #times DCA distribution in 3#circ < #theta_{abs} < 10#circ;p #times DCA (GeV.cm/c)", 2500, 0., 5000.);
  fList->AddAtAndExpand(hPDCA310, kPDCA310);
  
  TH1F* hPtMuPlus = new TH1F("hPtMuPlus", "#mu^{+} transverse momentum distribution;p_{t} (GeV/c)", 500, 0., 100.);
  fList->AddAtAndExpand(hPtMuPlus, kPtMuPlus);
  
  TH1F* hPtMuMinus = new TH1F("hPtMuMinus", "#mu^{-} transverse momentum distribution;p_{t} (GeV/c)", 500, 0., 100.);
  fList->AddAtAndExpand(hPtMuMinus, kPtMuMinus);
  
  TH1F* hSPDzVtx = new TH1F("hSPDzVtx", "SPD z vertex;z (cm)", 500, -25., 25.);
  fList->AddAtAndExpand(hSPDzVtx, kSPDzVtx);
  
  TH1F* hT0zVtx = new TH1F("hT0zVtx", "T0 z vertex;z (cm)", 500, -25., 25.);
  fList->AddAtAndExpand(hT0zVtx, kT0zVtx);
  
  TH1F* hT0SPDDeltaZVtx = new TH1F("hT0SPDDeltaZVtx", "T0 - SPD #Delta_z vertex;#Delta_z (cm)", 500, -25., 25.);
  fList->AddAtAndExpand(hT0SPDDeltaZVtx, kT0SPDDeltaZVtx);
  
  TH2F* hT0VsSPDZVtx = new TH2F("hT0VsSPDZVtx", "T0 versus SPD z vertex;SPD z (cm);T0 z (cm)", 500, -25., 25., 500, -25., 25.);
  fList->AddAtAndExpand(hT0VsSPDZVtx, kT0VsSPDZVtx);
  
  TH1F* hPtDimu = new TH1F("hPtDimu", "#dimuon transverse momentum distribution;p_{t} (GeV/c)", 200, 0., 20.);
  fList->AddAtAndExpand(hPtDimu, kPtDimu);
  
  TH1F* hInvBendingP = new TH1F("hInvBendingP", "#mu inverse bending momentum distribution;1/p_{yz} (c/GeV)", 2000, -0.1, 0.1);
  fList->AddAtAndExpand(hInvBendingP, kInvBendingP);
  
  TH2F* hDCA23VsP = new TH2F("hDCA23VsP", "DCA distribution versus p in 2#circ < #theta_{abs} < 3#circ;p(GeV/c);DCA (cm)", 160, 0., 800., 400, 0., 200.);
  fList->AddAtAndExpand(hDCA23VsP, kDCA23VsP);
  
  TH2F* hDCA310VsP = new TH2F("hDCA310VsP", "DCA distribution versus p in 3#circ < #theta_{abs} < 10#circ;p(GeV/c);DCA (cm)", 160, 0., 800., 400, 0., 200.);
  fList->AddAtAndExpand(hDCA310VsP, kDCA310VsP);
  
  TH1F* hCent = new TH1F("hCent", "centrality percentile;centrality (%)", 10000, 0., 100.);
  fList->AddAtAndExpand(hCent, kCent);
  
  TH1F* hMult = new TH1F("hMult", "tracklet multiplicity; n tracklets", 10000, 0., 10000.);
  fList->AddAtAndExpand(hMult, kMult);
  
  TH1F* hMultSelect = new TH1F("hMultSelect", "tracklet multiplicity for events with selected tracks; n tracklets", 10000, 0., 10000.);
  fList->AddAtAndExpand(hMultSelect, kMultSelect);
  
  // centrality binning
  TString centbins = "any";
  for (Int_t i=0; i<100; i++) centbins += Form("/%d",i);
  
  // initialize track counters
  fTrackCounters = new AliCounterCollection(GetOutputSlot(2)->GetContainer()->GetName());
  fTrackCounters->AddRubric("run", 100000);
  fTrackCounters->AddRubric("track", "trackeronly/matched");
  fTrackCounters->AddRubric("trig", "any/low/high");
  fTrackCounters->AddRubric("rabs", "yes/no");
  fTrackCounters->AddRubric("eta", "yes/no");
  fTrackCounters->AddRubric("pdca", "yes/no");
  fTrackCounters->AddRubric("chi2", "yes/no");
  fTrackCounters->AddRubric("pt", "any/0.5GeV/1GeV/2GeV/4GeV");
  fTrackCounters->AddRubric("cent", centbins.Data());
  fTrackCounters->Init();
  
  // initialize event counters
  fEventCounters = new AliCounterCollection(GetOutputSlot(3)->GetContainer()->GetName());
  fEventCounters->AddRubric("event", "all/trigger/selected");
  fEventCounters->AddRubric("cent", centbins.Data());
  fEventCounters->AddRubric("spdvtx", "yes/no");
  fEventCounters->AddRubric("t0vtx", "yes/no");
  fEventCounters->Init();
  
  // initialize trigger track counters
  fTrigCounters = new AliCounterCollection(GetOutputSlot(4)->GetContainer()->GetName());
  fTrigCounters->AddRubric("run", 100000);
  fTrigCounters->AddRubric("board", 242);
  fTrigCounters->AddRubric("match", 4);
  fTrigCounters->AddRubric("ntrig", 1000);
  fTrigCounters->AddRubric("cent", centbins.Data());
  fTrigCounters->Init();
  
  // Post data at least once per task to ensure data synchronisation (required for merging)
  PostData(1, fList);
  PostData(2, fTrackCounters);
  PostData(3, fEventCounters);
  PostData(4, fTrigCounters);
}

//________________________________________________________________________
void AliAnalysisTaskMuonPhysics::NotifyRun()
{
  /// Prepare processing of new run: load corresponding OADB objects...
  
  // get the trackCuts for this run
  if (!fMuonTrackCuts) AliFatal("You must specify the requested selections (AliMuonTrackCut obj is missing)");
  if (fMuonTrackCuts->GetPassName().IsNull())
    fMuonTrackCuts->SetCustomParamFromRun(fCurrentRunNumber, AliAnalysisMuonUtility::GetPassName(fInputHandler));
  else
    fMuonTrackCuts->SetCustomParamFromRun(fCurrentRunNumber, fMuonTrackCuts->GetPassName());
  fMuonTrackCuts->CustomParam()->SetChi2NormCut(3.5);
  
  // get the trackCuts2 for this run
  if (fMuonTrackCuts2) fMuonTrackCuts2->SetRun(fInputHandler);
  
}

//________________________________________________________________________
void AliAnalysisTaskMuonPhysics::UserExec(Option_t *)
{
  /// Called for each event
  
  AliVEvent *evt = InputEvent();
  Bool_t isESD = kFALSE;
  if (dynamic_cast<AliESDEvent*>(evt)) isESD = kTRUE;
  else if (!dynamic_cast<AliAODEvent*>(evt)) return;
  
  if (isESD) LoadBranches();
  /*
  Int_t ievent = static_cast<AliESDEvent*>(evt)->GetEventNumberInFile();
  if (ievent == 129) return;
  if (ievent == 131) return;
  if (ievent == 141) return;
  //if (ievent == 165) return;
  if (ievent == 167) return;
  if (ievent == 276) return;
  if (ievent == 1081) return;
  if (ievent == 1219) return;
  if (ievent == 1256) return;
  if (ievent == 1435) return;
  if (ievent == 1542) return;
  if (ievent == 2004) return;
  if (ievent == 2284) return;
  if (ievent == 3451) return;
  if (ievent == 3558) return;
  if (ievent == 3789) return;
  if (ievent == 3903) return;
  if (ievent == 3978) return;
  if (ievent == 4185) return;
  if (ievent == 4365) return;
  if (ievent == 4401) return;
  if (ievent == 4434) return;
  if (ievent == 4593) return;
  if (ievent == 4950) return;
  if (ievent == 5014) return;
  if (ievent == 5064) return;
  //if (ievent == 5208) return;
  if (ievent == 5215) return;
  if (ievent == 5263) return;
  if (ievent == 5295) return;
  if (ievent == 5602) return;
  if (ievent == 5904) return;
  if (ievent == 6382) return;
  if (ievent == 6422) return;
  if (ievent == 6479) return;
  if (ievent == 6573) return;
  if (ievent == 6899) return;
  if (ievent == 7016) return;
  if (ievent == 7243) return;
  if (ievent == 7434) return;
  if (ievent == 7532) return;
  if (ievent == 8064) return;
  if (ievent == 8109) return;
  if (ievent == 8140) return;
  if (ievent == 8642) return;
  //if (ievent == 8651) return;
  if (ievent == 8678) return;
  if (ievent == 8705) return;
  */
  // get the centrality percentile
  Float_t centrality = evt->GetCentrality()->GetCentralityPercentile("V0M");
//  Float_t centrality = evt->GetCentrality()->GetCentralityPercentile("V0A");
  if (centrality <= fCentMin || centrality > fCentMax) return;
  TString centKey = (centrality >0 && centrality < 100) ? Form("%d", (Int_t)centrality) : "";
  
  // get the tracklet multiplicity
  Int_t ntracklets = 0;
  if (isESD) {
    const AliMultiplicity *spdmult = static_cast<AliESDEvent*>(evt)->GetMultiplicity();
    if (spdmult) for (Int_t i = 0; i < spdmult->GetNumberOfTracklets(); i++)
      if (TMath::Abs(spdmult->GetEta(i)) < 0.8) ntracklets++;
  } else {
    const AliAODTracklets *spdmult = static_cast<AliAODEvent*>(evt)->GetTracklets();
    if (spdmult) for (Int_t i = 0; i < spdmult->GetNumberOfTracklets(); i++)
      if (TMath::Abs(-TMath::Log(TMath::Tan(spdmult->GetTheta(i)/2.))) < 0.8) ntracklets++;
  }
  
  // get the SPD vertex
  TString vtxSPD = "spdvtx:";
  Double_t vertex[3] = {0., 0., 0.};
  const AliVVertex* vert = AliAnalysisMuonUtility::GetVertexSPD(evt);
  if (vert && GetVtxStatus(*vert)) {
    vert->GetXYZ(vertex);
    vtxSPD += "yes";
  } else vtxSPD += "no";
  
  // get the T0 vertex
  TString vtxT0 = "t0vtx:yes";
  Double_t zVtxT0 = isESD ? static_cast<AliESDEvent*>(evt)->GetT0zVertex() : static_cast<AliAODEvent*>(evt)->GetTZEROData()->GetT0VertexRaw();
  if (TMath::Abs(zVtxT0) > 98.) vtxT0 = "t0vtx:no";
  
  // total number of events
  fEventCounters->Count(Form("event:all/cent:any/%s/%s",vtxSPD.Data(),vtxT0.Data()));
  if (!centKey.IsNull()) fEventCounters->Count(Form("event:all/cent:%s/%s/%s",centKey.Data(),vtxSPD.Data(),vtxT0.Data()));
  
  // first loop over tracks
  Int_t nTracks = AliAnalysisMuonUtility::GetNTracks(evt);
  /*
  // cut events according to the set of trigger tracks
  Int_t n = 0;
  Bool_t adjacentTracks = kFALSE;
  for (Int_t i = 0; i < nTracks; i++) {
    
    AliVParticle *track1 = AliAnalysisMuonUtility::GetTrack(i, evt);
    
    Int_t loCircuit1 = AliAnalysisMuonUtility::GetLoCircuit(track1);
    if (loCircuit1 == 0) continue;
    n++;
    
    for (Int_t j = 0; j < i; j++) {
      
      AliVParticle *track2 = AliAnalysisMuonUtility::GetTrack(j, evt);
      
      Int_t loCircuit2 = AliAnalysisMuonUtility::GetLoCircuit(track2);
      if (loCircuit2 == 0) continue;
      
      if (TMath::Abs(loCircuit2-loCircuit1) <= 1) adjacentTracks = kTRUE;
      
    }
    
  }
  //  if (n > 5) return;
  if (adjacentTracks) return;
  */
  Int_t nTrkTracks = 0, nTrgTracks = 0, nTrgTracksLow = 0, nTrgTracksHigh = 0, nMatchTracks = 0;
  Int_t (*loCircuit)[2] = new Int_t[nTracks][2];
  Bool_t containTriggerTrack = kFALSE;
  Bool_t containSelectedTrack = kFALSE;
  for (Int_t iTrack1 = 0; iTrack1 < nTracks; iTrack1++) {
    
    // get the track
    AliVParticle *track1 = AliAnalysisMuonUtility::GetTrack(iTrack1, evt);
    if (!isESD && !static_cast<AliAODTrack*>(track1)->IsMuonTrack()) continue;
    
    Bool_t useTrk = isESD ? ((fSelectBadTracks && track1->TestBit(BIT(20))) ||
			     (!fSelectBadTracks && !track1->TestBit(BIT(20)))) : kTRUE;
    
    Bool_t useTrg = isESD ? ((fSelectBadTracks && track1->TestBit(BIT(21))) ||
			     (!fSelectBadTracks && !track1->TestBit(BIT(21)))) : kTRUE;
    
    Int_t loCircuit1 = useTrg ? AliAnalysisMuonUtility::GetLoCircuit(track1) : 0;
    Int_t matchTrig1 = useTrg ? AliAnalysisMuonUtility::GetMatchTrigger(track1) : 0;
    
    // count the number of tracker tracks, matched or not
    if (useTrk && AliAnalysisMuonUtility::IsMuonTrack(track1)) {
      nTrkTracks++;
      if (matchTrig1 > 0) nMatchTracks++;
    }
    
    // count the number of trigger tracks
    if (loCircuit1) {
      
      // make sure this track is not already accounted for
      Bool_t trackExist = kFALSE;
      for (Int_t i=0; i<nTrgTracks; i++) {
	if (loCircuit1 == loCircuit[i][0]) {
	  trackExist = kTRUE;
	  break;
	}
      }
      
      // count according trigger level
      if (!trackExist) {
	containTriggerTrack = kTRUE;
	loCircuit[nTrgTracks][0] = loCircuit1;
	loCircuit[nTrgTracks][1] = matchTrig1;
	nTrgTracks++;
	if (matchTrig1 > 1) nTrgTracksLow++;
	if (matchTrig1 > 2) nTrgTracksHigh++;
      }
      
    }
    
    // skip track passing the second set of cuts
    Bool_t isRejected1 = (fMuonTrackCuts2 && fMuonTrackCuts2->IsSelected(track1));
    
    // get the ESD track seleted for physics
    Bool_t isSelected1 = IsSelected(*track1, isESD, !isRejected1, centKey, fSelectBadTracks);
    
    if (isSelected1 && !isRejected1) {
      
      containSelectedTrack = kTRUE;
      
      // fill single muon histograms
      Double_t pT = track1->Pt();
      ((TH1F*)fList->UncheckedAt(kPt))->Fill(pT);
      if (track1->Charge() > 0) ((TH1F*)fList->UncheckedAt(kPtMuPlus))->Fill(pT);
      else ((TH1F*)fList->UncheckedAt(kPtMuMinus))->Fill(pT);
      ((TH1F*)fList->UncheckedAt(kRapidity))->Fill(track1->Y());
      Double_t dcaX = AliAnalysisMuonUtility::GetXatDCA(track1) - AliAnalysisMuonUtility::GetXatVertex(track1);
      Double_t dcaY = AliAnalysisMuonUtility::GetYatDCA(track1) - AliAnalysisMuonUtility::GetYatVertex(track1);
      ((TH1F*)fList->UncheckedAt(kDCA))->Fill(TMath::Sqrt(dcaX*dcaX + dcaY*dcaY));
      ((TH1F*)fList->UncheckedAt(kChi2))->Fill(AliAnalysisMuonUtility::GetChi2perNDFtracker(track1));
      ((TH1F*)fList->UncheckedAt(kNClustersPerTrack))->Fill(isESD ? static_cast<AliESDMuonTrack*>(track1)->GetNHit() : 0.);
      Int_t nChamberHit = 0;
//      for (Int_t ich=0; ich<10; ich++) if (AliAnalysisMuonUtility::IsTrkChamberHit(ich, track1)) nChamberHit++;
      for (Int_t ich=0; ich<10; ich++)
	if ((isESD && static_cast<AliESDMuonTrack*>(track1)->IsInMuonClusterMap(ich)) ||
	    (!isESD && static_cast<AliAODTrack*>(track1)->HitsMuonChamber(ich))) nChamberHit++;
      ((TH1F*)fList->UncheckedAt(kNChamberHitPerTrack))->Fill(nChamberHit);
      ((TH1F*)fList->UncheckedAt(kPUncorrected))->Fill(isESD ? static_cast<AliESDMuonTrack*>(track1)->PUncorrected() : 0.);
      ((TH1F*)fList->UncheckedAt(kRAbs))->Fill(AliAnalysisMuonUtility::GetRabs(track1));
      ((TH1F*)fList->UncheckedAt(kDCAX))->Fill(dcaX);
      ((TH1F*)fList->UncheckedAt(kDCAY))->Fill(dcaY);
//      ((TH1F*)fList->UncheckedAt(kChi2Trig))->Fill(AliAnalysisMuonUtility::GetChi2MatchTrigger(track1));
      ((TH1F*)fList->UncheckedAt(kChi2Trig))->Fill(isESD ? static_cast<AliESDMuonTrack*>(track1)->GetChi2MatchTrigger() :
						   static_cast<AliAODTrack*>(track1)->GetChi2MatchTrigger());
      Double_t averageP = fMuonTrackCuts->GetAverageMomentum(track1);
      Double_t corrDCA = fMuonTrackCuts->GetCorrectedDCA(track1).Mag();
      Double_t pdca = averageP*corrDCA;
      Double_t thetaAbs = AliAnalysisMuonUtility::GetThetaAbsDeg(track1);
      if (thetaAbs > 2 && thetaAbs < 3) {
	((TH1F*)fList->UncheckedAt(kPDCA23))->Fill(pdca);
	((TH2F*)fList->UncheckedAt(kDCA23VsP))->Fill(averageP,corrDCA);
      }
      else if (thetaAbs >= 3 && thetaAbs < 10) {
	((TH1F*)fList->UncheckedAt(kPDCA310))->Fill(pdca);
	((TH2F*)fList->UncheckedAt(kDCA310VsP))->Fill(averageP,corrDCA);
      }
      if (isESD) ((TH1F*)fList->UncheckedAt(kInvBendingP))->Fill(static_cast<AliESDMuonTrack*>(track1)->GetInverseBendingMomentumUncorrected());
      
    } else if (!isSelected1 && !(isESD && fSelectBadTracks && IsSelected(*track1, isESD, kFALSE, centKey, kFALSE))) continue;
    
    TLorentzVector muV1(track1->Px(), track1->Py(), track1->Pz(), track1->E());
    Short_t charge1 = track1->Charge();
    
    // second loop over tracks
    for (Int_t iTrack2 = iTrack1 + 1; iTrack2 < nTracks; iTrack2++) {
      
      // get the track seleted for physics
      AliVParticle *track2 = AliAnalysisMuonUtility::GetTrack(iTrack2, evt);
      Bool_t isSelected2 = IsSelected(*track2, isESD, kFALSE, centKey, fSelectBadTracks);
      
      if (isSelected2 || (isESD && fSelectBadTracks && isSelected1 && IsSelected(*track2, isESD, kFALSE, centKey, kFALSE))) {
	
	// skip dimuons with both tracks passing the second set of cuts
	Bool_t isRejected2 = (fMuonTrackCuts2 && fMuonTrackCuts2->IsSelected(track2));
	if (isRejected1 && isRejected2) continue;
	
	// select opposite charge dimuons
	if (charge1 * track2->Charge() > 0) continue;
	
	TLorentzVector muV2(track2->Px(), track2->Py(), track2->Pz(), track2->E());
	TLorentzVector dimuV = muV1 + muV2;
	
	// fill dimuon histograms
	((TH1F*)fList->UncheckedAt(kMass))->Fill(dimuV.M());
	((TH1F*)fList->UncheckedAt(kPtDimu))->Fill(dimuV.Pt());
	
      }
      
    }
    
  }
  
  // fill trigger track counters
  TString runKey = fVsRun ? Form("run:%d", fCurrentRunNumber) : "run:any";
  for (Int_t i = 0; i < nTrgTracks; i++) {
    fTrigCounters->Count(Form("%s/board:%d/match:%d/ntrig:%d/cent:any",runKey.Data(),loCircuit[i][0],loCircuit[i][1],nTrgTracks));
    if (!centKey.IsNull())
      fTrigCounters->Count(Form("%s/board:%d/match:%d/ntrig:%d/cent:%s",runKey.Data(),loCircuit[i][0],loCircuit[i][1],nTrgTracks,centKey.Data()));
  }
  
  // clean memory
  delete[] loCircuit;
  
  // fill histo per event
  ((TH1F*)fList->UncheckedAt(kNTracks))->Fill(nTrkTracks);
  ((TH1F*)fList->UncheckedAt(kNTrigAll))->Fill(nTrgTracks);
  ((TH1F*)fList->UncheckedAt(kNTrigLow))->Fill(nTrgTracksLow);
  ((TH1F*)fList->UncheckedAt(kNTrigHigh))->Fill(nTrgTracksHigh);
  ((TH1F*)fList->UncheckedAt(kNMatchTracks))->Fill(nMatchTracks);
  if (vtxSPD.Contains("yes")) ((TH1F*)fList->UncheckedAt(kSPDzVtx))->Fill(vertex[2]);
  if (vtxT0.Contains("yes")) ((TH1F*)fList->UncheckedAt(kT0zVtx))->Fill(zVtxT0);
  if (vtxSPD.Contains("yes") && vtxT0.Contains("yes")
      && zVtxT0 > -15. && zVtxT0 < 15. && vertex[2] > -10. && vertex[2] < 10.
      ) ((TH1F*)fList->UncheckedAt(kT0SPDDeltaZVtx))->Fill(zVtxT0 - vertex[2]);
  if (vtxSPD.Contains("yes") && vtxT0.Contains("yes"))
    ((TH2F*)fList->UncheckedAt(kT0VsSPDZVtx))->Fill(vertex[2], zVtxT0);
  ((TH1F*)fList->UncheckedAt(kCent))->Fill(centrality);
  if (vtxSPD.Contains("yes") && vertex[2] > -10. && vertex[2] < 10.) {
    ((TH1F*)fList->UncheckedAt(kMult))->Fill(ntracklets);
    if (containSelectedTrack) ((TH1F*)fList->UncheckedAt(kMultSelect))->Fill(ntracklets);
  }
  
  // number of events with at least one trigger track
  if (containTriggerTrack) {
    fEventCounters->Count(Form("event:trigger/cent:any/%s/%s",vtxSPD.Data(),vtxT0.Data()));
    if (!centKey.IsNull()) fEventCounters->Count(Form("event:trigger/cent:%s/%s/%s",centKey.Data(),vtxSPD.Data(),vtxT0.Data()));
  }
  
  // number of events with at least one selected track
  if (containSelectedTrack) {
    fEventCounters->Count(Form("event:selected/cent:any/%s/%s",vtxSPD.Data(),vtxT0.Data()));
    if (!centKey.IsNull()) fEventCounters->Count(Form("event:selected/cent:%s/%s/%s",centKey.Data(),vtxSPD.Data(),vtxT0.Data()));
  }
  
  // Post final data. It will be written to a file with option "RECREATE"
  PostData(1, fList);
  PostData(2, fTrackCounters);
  PostData(3, fEventCounters);
  PostData(4, fTrigCounters);
}

//________________________________________________________________________
void AliAnalysisTaskMuonPhysics::Terminate(Option_t *)
{
  fList = dynamic_cast<TObjArray*>(GetOutputData(1));
  if (!fList) return;
  
  // reset track cuts
  fMuonTrackCuts->SetCustomParamFromRun(169859, "pass2_muon");
  fMuonTrackCuts->Print();
  
  // pDCA cut functions
  Double_t nSigmaCut = fMuonTrackCuts->GetMuonTrackCutsParam().GetNSigmaPdca();
  Double_t sigmaMeasCut[2] = { fMuonTrackCuts->GetMuonTrackCutsParam().GetSigmaPdca23(), fMuonTrackCuts->GetMuonTrackCutsParam().GetSigmaPdca310()};
  Double_t relPResolution = fMuonTrackCuts->GetMuonTrackCutsParam().GetRelPResolution();
  Double_t angleResolution = 535.*fMuonTrackCuts->GetMuonTrackCutsParam().GetSlopeResolution();
  TString cutFormula = "[1]*TMath::Sqrt( ( [0] / ( 1. - [1]*[2]*x / ( 1.+[1]*[2]*x ) ) ) * ( [0] / ( 1. - [1]*[2]*x / ( 1.+[1]*[2]*x ) ) ) + [3]*[3]*x*x) / x";
  TF1* cutFunction23 = new TF1("fpDCACut23",cutFormula.Data(), 0.1, 800.);
  cutFunction23->SetParameters(sigmaMeasCut[0], nSigmaCut, relPResolution, angleResolution);
  TF1* cutFunction310 = new TF1("fpDCACut310",cutFormula.Data(), 0.1, 800.);
  cutFunction310->SetParameters(sigmaMeasCut[1], nSigmaCut, relPResolution, angleResolution);
  
  // draw DCA vs p with cut functions
  TCanvas *c = new TCanvas("cDCAvsP","cDCAvsP",800,400);
  c->Divide(2,1);
  c->cd(1);
  gPad->SetLogz();
  ((TH2F*)fList->UncheckedAt(kDCA23VsP))->Draw("COLZ");
  cutFunction23->SetLineWidth(2);
  cutFunction23->SetLineColor(2);
  cutFunction23->DrawClone("same");
  cutFunction23->SetLineWidth(1);
  cutFunction23->SetParameters(80., nSigmaCut, 0.0004, 535.*0.0004);
  cutFunction23->DrawClone("same");
  cutFunction23->SetParameters(80., nSigmaCut, 0.0005, 535.*0.0005);
  cutFunction23->DrawClone("same");
  c->cd(2);
  gPad->SetLogz();
  ((TH2F*)fList->UncheckedAt(kDCA310VsP))->Draw("COLZ");
  cutFunction310->SetLineWidth(2);
  cutFunction310->SetLineColor(2);
  cutFunction310->DrawClone("same");
  cutFunction310->SetLineWidth(1);
  cutFunction310->SetParameters(54., nSigmaCut, 0.0004, 535.*0.0004);
  cutFunction310->DrawClone("same");
  cutFunction310->SetParameters(54., nSigmaCut, 0.0005, 535.*0.0005);
  cutFunction310->DrawClone("same");
  
}

//________________________________________________________________________
Bool_t AliAnalysisTaskMuonPhysics::IsSelected(AliVParticle& track, Bool_t isESD, Bool_t fillCounters,
					      TString &centKey, Bool_t selectBadTracks)
{
  /// return kTRUE if the track passes all the cuts for physics
  
  // skip ghosts
  if (!AliAnalysisMuonUtility::IsMuonTrack(&track)) return kFALSE;
  
  // skip tracks without MC label
  if (fUseMCLabel && track.GetLabel() < 0) return kFALSE;
  
  // skip or keep only bad tracks
  if (isESD && !selectBadTracks && track.TestBit(BIT(20))) return kFALSE;
  if (isESD && selectBadTracks && !track.TestBit(BIT(20)) && !track.TestBit(BIT(21))) return kFALSE;
  /*
  // skip tracks with not enough clusters to be still reconstructible
  AliMUONTrack muonTrack;
  AliMUONESDInterface::ESDToMUON(esdTrack, muonTrack, kFALSE);
  if (!muonTrack.IsValid(0x1F, kTRUE)) return kFALSE;
  *//*
  // skip tracks with less than 3 chambers hit in station 4 and 5
  Int_t nFiredChInSt45 = 0;
  for (Int_t iCh = 6; iCh < 10; iCh++) if (AliAnalysisMuonUtility::IsTrkChamberHit(iCh, &track)) nFiredChInSt45++;
  if (nFiredChInSt45 < 3) return kFALSE;
  */
  // skip tracks with pT < 2 GeV/c
  //if (track.Pt() < 60.) return kFALSE;
  
  // skip tracks with p < 200 GeV/c
  //if (track.P() < 100.) return kFALSE;
  
  // skip tracks with DCA > 20cm || < 5cm
  //Double_t corrDCA = fMuonTrackCuts->GetCorrectedDCA(&track).Mag();
  //if (corrDCA > 20.) return kFALSE;
  
  // skip decays
  //if (isESD && track.TestBit(BIT(22))) return kFALSE;
  
  //select tracks on left/right side of the spectrometer
  //if (isESD && track.Px() < 0. && static_cast<AliESDMuonTrack&>(track).GetNonBendingCoorUncorrected() < 0.) return kFALSE;
  
  UInt_t selectionMask = fMuonTrackCuts->GetSelectionMask(&track);
  
  Bool_t useTrg = isESD ? ((selectBadTracks && track.TestBit(BIT(21))) ||
			   (!selectBadTracks && !track.TestBit(BIT(21)))) : kTRUE;
  
  // cut on trigger matching
  Bool_t selectTrig = kFALSE;
  TString trackKey = "track:";
  if (useTrg && AliAnalysisMuonUtility::MatchApt(&track)) {
    selectTrig = kTRUE;
    trackKey += "matched";
  } else trackKey += "trackeronly";
  
  // fill counters according to selecting cuts if requested
  if (fillCounters) {
    
    // run
    TString runKey = fVsRun ? Form("run:%d", fCurrentRunNumber) : "run:any";
    
    // list of trigger types
    TList listOfTrigKeys;
    listOfTrigKeys.SetOwner();
    listOfTrigKeys.AddLast(new TObjString("trig:any"));
    if (selectionMask & AliMuonTrackCuts::kMuMatchLpt) listOfTrigKeys.AddLast(new TObjString("trig:low"));
    if (selectionMask & AliMuonTrackCuts::kMuMatchHpt) listOfTrigKeys.AddLast(new TObjString("trig:high"));
    
    // cut on ThetaAbs
    TString rabsKey = "rabs:";
    if (selectionMask & AliMuonTrackCuts::kMuThetaAbs) rabsKey += "yes";
    else rabsKey += "no";
    
    // cut on eta
    TString etaKey = "eta:";
    if (selectionMask & AliMuonTrackCuts::kMuEta) etaKey += "yes";
    else etaKey += "no";
    
    // cut on pDCA
    TString pdcaKey = "pdca:";
    if (selectionMask & AliMuonTrackCuts::kMuPdca) pdcaKey += "yes";
    else pdcaKey += "no";
    
    // cut on chi2
    TString chi2Key = "chi2:";
    if (selectionMask & AliMuonTrackCuts::kMuTrackChiSquare) chi2Key += "yes";
    else chi2Key += "no";
    
    TIter nextTrigKey(&listOfTrigKeys);
    TObjString *trigKey = 0x0;
    while ((trigKey = static_cast<TObjString*>(nextTrigKey()))) {
      fTrackCounters->Count(Form("%s/%s/%s/%s/%s/%s/%s/pt:any/cent:any", runKey.Data(), trackKey.Data(), trigKey->GetName(), rabsKey.Data(), etaKey.Data(), pdcaKey.Data(), chi2Key.Data()));
      if (!centKey.IsNull()) fTrackCounters->Count(Form("%s/%s/%s/%s/%s/%s/%s/pt:any/cent:%s", runKey.Data(), trackKey.Data(), trigKey->GetName(), rabsKey.Data(), etaKey.Data(), pdcaKey.Data(), chi2Key.Data(), centKey.Data()));
      Double_t pt = track.Pt();
      if (pt > 0.5) {
	fTrackCounters->Count(Form("%s/%s/%s/%s/%s/%s/%s/pt:0.5GeV/cent:any", runKey.Data(), trackKey.Data(), trigKey->GetName(), rabsKey.Data(), etaKey.Data(), pdcaKey.Data(), chi2Key.Data()));
	if (!centKey.IsNull()) fTrackCounters->Count(Form("%s/%s/%s/%s/%s/%s/%s/pt:0.5GeV/cent:%s", runKey.Data(), trackKey.Data(), trigKey->GetName(), rabsKey.Data(), etaKey.Data(), pdcaKey.Data(), chi2Key.Data(), centKey.Data()));
      }
      if (pt > 1.) {
	fTrackCounters->Count(Form("%s/%s/%s/%s/%s/%s/%s/pt:1GeV/cent:any", runKey.Data(), trackKey.Data(), trigKey->GetName(), rabsKey.Data(), etaKey.Data(), pdcaKey.Data(), chi2Key.Data()));
	if (!centKey.IsNull()) fTrackCounters->Count(Form("%s/%s/%s/%s/%s/%s/%s/pt:1GeV/cent:%s", runKey.Data(), trackKey.Data(), trigKey->GetName(), rabsKey.Data(), etaKey.Data(), pdcaKey.Data(), chi2Key.Data(), centKey.Data()));
      }
      if (pt > 2.) {
	fTrackCounters->Count(Form("%s/%s/%s/%s/%s/%s/%s/pt:2GeV/cent:any", runKey.Data(), trackKey.Data(), trigKey->GetName(), rabsKey.Data(), etaKey.Data(), pdcaKey.Data(), chi2Key.Data()));
	if (!centKey.IsNull()) fTrackCounters->Count(Form("%s/%s/%s/%s/%s/%s/%s/pt:2GeV/cent:%s", runKey.Data(), trackKey.Data(), trigKey->GetName(), rabsKey.Data(), etaKey.Data(), pdcaKey.Data(), chi2Key.Data(), centKey.Data()));
      }
      if (pt > 4.) {
	fTrackCounters->Count(Form("%s/%s/%s/%s/%s/%s/%s/pt:4GeV/cent:any", runKey.Data(), trackKey.Data(), trigKey->GetName(), rabsKey.Data(), etaKey.Data(), pdcaKey.Data(), chi2Key.Data()));
	if (!centKey.IsNull()) fTrackCounters->Count(Form("%s/%s/%s/%s/%s/%s/%s/pt:4GeV/cent:%s", runKey.Data(), trackKey.Data(), trigKey->GetName(), rabsKey.Data(), etaKey.Data(), pdcaKey.Data(), chi2Key.Data(), centKey.Data()));
      }
    }
  }
  
  UInt_t filterMask = fMuonTrackCuts->GetFilterMask();
  Bool_t filterTrig = (filterMask & (AliMuonTrackCuts::kMuMatchApt | AliMuonTrackCuts::kMuMatchLpt | AliMuonTrackCuts::kMuMatchHpt));
  if (!filterTrig && isESD && selectBadTracks && !track.TestBit(BIT(20))) return kFALSE;
  return ((!filterTrig || selectTrig) && ((selectionMask & filterMask) == filterMask));
}

//________________________________________________________________________
Bool_t AliAnalysisTaskMuonPhysics::GetVtxStatus(const AliVVertex &vtx) const
{
  /// return the status of the vertex (copy of AliVertex::GetStatus())
  TString title = vtx.GetTitle();
  if(vtx.GetNContributors()>0 || (title.Contains("cosmics") && !title.Contains("failed"))) return 1;
  if(title.Contains("smearMC")) return 1;
  return 0;
}

