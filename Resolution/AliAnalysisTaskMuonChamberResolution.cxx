
// ROOT includes
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TMultiGraph.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <Riostream.h>
#include <TString.h>
#include <TGeoManager.h>

// STEER includes
#include "AliESDEvent.h"
#include "AliESDMuonTrack.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliGeomManager.h"

// ANALYSIS includes
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisDataSlot.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskMuonChamberResolution.h"

// MUON includes
#include "AliMUONCDB.h"
#include "AliMUONRecoParam.h"
#include "AliMUONESDInterface.h"
#include "AliMUONVTrackReconstructor.h"
#include "AliMUONTrack.h"
#include "AliMUONTrackParam.h"
#include "AliMUONTrackExtrap.h"
#include "AliMUONVCluster.h"
#include "AliMUONGeometryTransformer.h"
#include "AliMUONGeometryModuleTransformer.h"
#include "AliMUONGeometryDetElement.h"
#include "AliMpDEIterator.h"

#ifndef SafeDelete
#define SafeDelete(x) if (x != NULL) { delete x; x = NULL; }
#endif

ClassImp(AliAnalysisTaskMuonChamberResolution)

const Int_t AliAnalysisTaskMuonChamberResolution::fgkMinEntries = 10;

//________________________________________________________________________
AliAnalysisTaskMuonChamberResolution::AliAnalysisTaskMuonChamberResolution() :
  AliAnalysisTaskSE(), 
  fResiduals(NULL),
  fResidualsVsP(NULL),
  fSummary(NULL),
  fCanvases(NULL),
  fNEvents(0),
  fMinMomentum(0.),
  fMatchTrig(kFALSE),
  fExtrapMode(1),
  fCorrectForSystematics(kTRUE),
  fOCDBLoaded(kFALSE),
  fNDE(0),
  fReAlign(kFALSE),
  fOldAlignStorage(""),
  fNewAlignStorage(""),
  fOldGeoTransformer(NULL),
  fNewGeoTransformer(NULL)
{
  /// Default constructor
}

//________________________________________________________________________
AliAnalysisTaskMuonChamberResolution::AliAnalysisTaskMuonChamberResolution(const char *name) :
AliAnalysisTaskSE(name), 
fResiduals(NULL),
fResidualsVsP(NULL),
fSummary(NULL),
fCanvases(NULL),
fNEvents(0),
fMinMomentum(0.),
fMatchTrig(kFALSE),
fExtrapMode(1),
fCorrectForSystematics(kTRUE),
fOCDBLoaded(kFALSE),
fNDE(0),
fReAlign(kFALSE),
fOldAlignStorage(""),
fNewAlignStorage(""),
fOldGeoTransformer(NULL),
fNewGeoTransformer(NULL)
{
  /// Constructor
  
  for (Int_t i = 0; i < AliMUONConstants::NTrackingCh(); i++) SetStartingResolution(i, -1., -1.);
  
  // Output slot #1 writes into a TObjArray container
  DefineOutput(1,TObjArray::Class());
  // Output slot #2 writes into a TObjArray container
  DefineOutput(2,TObjArray::Class());
  // Output slot #3 writes into a TObjArray container
  DefineOutput(3,TObjArray::Class());
}

//________________________________________________________________________
AliAnalysisTaskMuonChamberResolution::~AliAnalysisTaskMuonChamberResolution()
{
  /// Destructor
  SafeDelete(fResiduals);
  SafeDelete(fResidualsVsP);
  SafeDelete(fSummary);
  SafeDelete(fCanvases);
  SafeDelete(fOldGeoTransformer);
  SafeDelete(fNewGeoTransformer);
}

//___________________________________________________________________________
void AliAnalysisTaskMuonChamberResolution::UserCreateOutputObjects()
{
  /// Create histograms
  
  // do it once the OCDB has been loaded (i.e. from NotifyRun())
  if (!fOCDBLoaded) return;
  
  fResiduals = new TObjArray(1000);
  fResiduals->SetOwner();
  fResidualsVsP = new TObjArray(1000);
  fResidualsVsP->SetOwner();
  TH2F* h2;
  
  // find the highest chamber resolution and set histogram bins
  const AliMUONRecoParam* recoParam = AliMUONESDInterface::GetTracker()->GetRecoParam();
  Double_t maxSigma[2] = {-1., -1.};
  for (Int_t i = 0; i < 10; i++) {
    if (recoParam->GetDefaultNonBendingReso(i) > maxSigma[0]) maxSigma[0] = recoParam->GetDefaultNonBendingReso(i);
    if (recoParam->GetDefaultBendingReso(i) > maxSigma[1]) maxSigma[1] = recoParam->GetDefaultBendingReso(i);
  }
  const char* axes[2] = {"X", "Y"};
  const Int_t nBins = 2000;
  const Int_t nSigma = 10;
  const Int_t pNBins = 20;
  const Double_t pEdges[2] = {0., 50.};
  
  for (Int_t ia = 0; ia < 2; ia++) {
    
    Double_t maxRes = nSigma*maxSigma[ia];
    
    h2 = new TH2F(Form("hResidual%sPerCh_ClusterIn",axes[ia]), Form("cluster-track residual-%s distribution per chamber (cluster attached to the track);chamber ID;#Delta_{%s} (cm)",axes[ia],axes[ia]), 10, 0.5, 10.5, nBins, -maxRes, maxRes);
    fResiduals->AddAtAndExpand(h2, kResidualPerCh_ClusterIn+ia);
    h2 = new TH2F(Form("hResidual%sPerCh_ClusterOut",axes[ia]), Form("cluster-track residual-%s distribution per chamber (cluster not attached to the track);chamber ID;#Delta_{%s} (cm)",axes[ia],axes[ia]), 10, 0.5, 10.5, nBins, -2.*maxRes, 2.*maxRes);
    fResiduals->AddAtAndExpand(h2, kResidualPerCh_ClusterOut+ia);
    
    h2 = new TH2F(Form("hResidual%sPerHalfCh_ClusterIn",axes[ia]), Form("cluster-track residual-%s distribution per half chamber (cluster attached to the track);half chamber ID;#Delta_{%s} (cm)",axes[ia],axes[ia]), 20, 0.5, 20.5, nBins, -maxRes, maxRes);
    fResiduals->AddAtAndExpand(h2, kResidualPerHalfCh_ClusterIn+ia);
    h2 = new TH2F(Form("hResidual%sPerHalfCh_ClusterOut",axes[ia]), Form("cluster-track residual-%s distribution per half chamber (cluster not attached to the track);half chamber ID;#Delta_{%s} (cm)",axes[ia],axes[ia]), 20, 0.5, 20.5, nBins, -2.*maxRes, 2.*maxRes);
    fResiduals->AddAtAndExpand(h2, kResidualPerHalfCh_ClusterOut+ia);
    
    h2 = new TH2F(Form("hResidual%sPerDE_ClusterIn",axes[ia]), Form("cluster-track residual-%s distribution per DE (cluster not attached to the track);DE ID;#Delta_{%s} (cm)",axes[ia],axes[ia]), fNDE, 0.5, fNDE+0.5, nBins, -maxRes, maxRes);
    for (Int_t i = 1; i <= fNDE; i++) h2->GetXaxis()->SetBinLabel(i, Form("%d",fDEIds[i]));
    fResiduals->AddAtAndExpand(h2, kResidualPerDE_ClusterIn+ia);
    h2 = new TH2F(Form("hResidual%sPerDE_ClusterOut",axes[ia]), Form("cluster-track residual-%s distribution per DE (cluster not attached to the track);DE ID;#Delta_{%s} (cm)",axes[ia],axes[ia]), fNDE, 0.5, fNDE+0.5, nBins, -2.*maxRes, 2.*maxRes);
    for (Int_t i = 1; i <= fNDE; i++) h2->GetXaxis()->SetBinLabel(i, Form("%d",fDEIds[i]));
    fResiduals->AddAtAndExpand(h2, kResidualPerDE_ClusterOut+ia);
    
    h2 = new TH2F(Form("hTrackRes%sPerCh",axes[ia]), Form("track #sigma_{%s} per Ch;chamber ID;#sigma_{%s} (cm)",axes[ia],axes[ia]), 10, 0.5, 10.5, nBins, 0., maxRes);
    fResiduals->AddAtAndExpand(h2, kTrackResPerCh+ia);
    h2 = new TH2F(Form("hTrackRes%sPerHalfCh",axes[ia]), Form("track #sigma_{%s} per half Ch;half chamber ID;#sigma_{%s} (cm)",axes[ia],axes[ia]), 20, 0.5, 20.5, nBins, 0., maxRes);
    fResiduals->AddAtAndExpand(h2, kTrackResPerHalfCh+ia);
    h2 = new TH2F(Form("hTrackRes%sPerDE",axes[ia]), Form("track #sigma_{%s} per DE;DE ID;#sigma_{%s} (cm)",axes[ia],axes[ia]), fNDE, 0.5, fNDE+0.5, nBins, 0., maxRes);
    for (Int_t i = 1; i <= fNDE; i++) h2->GetXaxis()->SetBinLabel(i, Form("%d",fDEIds[i]));
    fResiduals->AddAtAndExpand(h2, kTrackResPerDE+ia);
    
    h2 = new TH2F(Form("hMCS%sPerCh",axes[ia]), Form("MCS %s-dispersion of extrapolated track per Ch;chamber ID;#sigma_{%s} (cm)",axes[ia],axes[ia]), 10, 0.5, 10.5, nBins, 0., 0.2);
    fResiduals->AddAtAndExpand(h2, kMCSPerCh+ia);
    h2 = new TH2F(Form("hMCS%sPerHalfCh",axes[ia]), Form("MCS %s-dispersion of extrapolated track per half Ch;half chamber ID;#sigma_{%s} (cm)",axes[ia],axes[ia]), 20, 0.5, 20.5, nBins, 0., 0.2);
    fResiduals->AddAtAndExpand(h2, kMCSPerHalfCh+ia);
    h2 = new TH2F(Form("hMCS%sPerDE",axes[ia]), Form("MCS %s-dispersion of extrapolated track per DE;DE ID;#sigma_{%s} (cm)",axes[ia],axes[ia]), fNDE, 0.5, fNDE+0.5, nBins, 0., 0.2);
    for (Int_t i = 1; i <= fNDE; i++) h2->GetXaxis()->SetBinLabel(i, Form("%d",fDEIds[i]));
    fResiduals->AddAtAndExpand(h2, kMCSPerDE+ia);
    
    h2 = new TH2F(Form("hClusterRes2%sPerCh",axes[ia]), Form("cluster #sigma_{%s}^{2} per Ch;chamber ID;#sigma_{%s}^{2} (cm^{2})",axes[ia],axes[ia]), 10, 0.5, 10.5, nSigma*nBins, -0.1*maxRes*maxRes, maxRes*maxRes);
    fResiduals->AddAtAndExpand(h2, kClusterRes2PerCh+ia);
    
    for (Int_t i = 0; i < AliMUONConstants::NTrackingCh(); i++) {
      h2 = new TH2F(Form("hResidual%sInCh%dVsP_ClusterIn",axes[ia],i+1), Form("cluster-track residual-%s distribution in chamber %d versus momentum (cluster attached to the track);p (GeV/c^{2});#Delta_{%s} (cm)",axes[ia],i+1,axes[ia]), pNBins, pEdges[0], pEdges[1], nBins, -maxRes, maxRes);
      fResidualsVsP->AddAtAndExpand(h2, kResidualInChVsP_ClusterIn+10*ia+i);
      h2 = new TH2F(Form("hResidual%sInCh%dVsP_ClusterOut",axes[ia],i+1), Form("cluster-track residual-%s distribution in chamber %d versus momentum (cluster not attached to the track);p (GeV/c^{2});#Delta_{%s} (cm)",axes[ia],i+1,axes[ia]), pNBins, pEdges[0], pEdges[1], nBins, -2.*maxRes, 2.*maxRes);
      fResidualsVsP->AddAtAndExpand(h2, kResidualInChVsP_ClusterOut+10*ia+i);
    }
    
  }
  
  // Post data at least once per task to ensure data synchronisation (required for merging)
  PostData(1, fResiduals);
  PostData(2, fResidualsVsP);
}

//________________________________________________________________________
void AliAnalysisTaskMuonChamberResolution::UserExec(Option_t *)
{
  /// Main event loop
  
  // check if OCDB properly loaded
  if (!fOCDBLoaded) return;
  
  AliESDEvent* esd = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!esd) return;
  
  if ((++fNEvents)%100 == 0) cout<<"\rEvent processing... "<<fNEvents<<"\r"<<flush;
  
  // get tracker to refit
  AliMUONVTrackReconstructor* tracker = AliMUONESDInterface::GetTracker();
  
  // loop over tracks
  Int_t nTracks = (Int_t) esd->GetNumberOfMuonTracks(); 
  for (Int_t iTrack = 0; iTrack < nTracks; ++iTrack) {
    
    // get the ESD track
    AliESDMuonTrack* esdTrack = esd->GetMuonTrack(iTrack);
    
    // skip ghost tracks
    if (!esdTrack->ContainTrackerData()) continue;
    
    // skip tracks not matched with trigger if required
    if (fMatchTrig && !esdTrack->ContainTriggerData()) continue;
    
    // skip low momentum tracks
    if (esdTrack->PUncorrected() < fMinMomentum) continue;
    
    // get the corresponding MUON track
    AliMUONTrack track;
    AliMUONESDInterface::ESDToMUON(*esdTrack, track, kFALSE);
    
    // change the cluster resolution (and position)
    ModifyClusters(track);
    
    // refit the track
    if (!tracker->RefitTrack(track, kFALSE)) break;
    
    // save track unchanged
    AliMUONTrack referenceTrack(track);
    Double_t p = static_cast<AliMUONTrackParam*>(track.GetTrackParamAtCluster()->First())->P();
    
    // loop over clusters
    Int_t nClusters = track.GetNClusters();
    for (Int_t iCluster=0; iCluster<nClusters; iCluster++) {
      
      // Get current, previous and next trackParams
      AliMUONTrackParam* trackParam = static_cast<AliMUONTrackParam*>(track.GetTrackParamAtCluster()->UncheckedAt(iCluster));
      AliMUONTrackParam* previousTrackParam = static_cast<AliMUONTrackParam*>(track.GetTrackParamAtCluster()->Before(trackParam));
      AliMUONTrackParam* nextTrackParam = static_cast<AliMUONTrackParam*>(track.GetTrackParamAtCluster()->After(trackParam));
      
      // save current trackParam and remove it from the track
      AliMUONTrackParam currentTrackParam(*trackParam);
      track.RemoveTrackParamAtCluster(trackParam);
      
      // get cluster info
      AliMUONVCluster* cluster = currentTrackParam.GetClusterPtr();
      Int_t chId = cluster->GetChamberId();
      Int_t halfChId = (cluster->GetX() > 0) ? 2*chId : 2*chId+1;
      Int_t deId = cluster->GetDetElemId();
      
      // make sure the track has another cluster in the same station and can still be refitted
      Bool_t refit = track.IsValid( 1 << (chId/2) );
      if (refit) {
	
	// refit the track and proceed if everything goes fine
	if (tracker->RefitTrack(track, kFALSE)) {
	  
	  // compute residuals
	  AliMUONTrackParam* referenceTrackParam = static_cast<AliMUONTrackParam*>(referenceTrack.GetTrackParamAtCluster()->UncheckedAt(iCluster));
	  Double_t deltaX = cluster->GetX() - referenceTrackParam->GetNonBendingCoor();
	  Double_t deltaY = cluster->GetY() - referenceTrackParam->GetBendingCoor();
	  
	  // fill histograms of residuals with cluster still attached to the track
	  ((TH2F*)fResiduals->UncheckedAt(kResidualPerCh_ClusterIn))->Fill(chId+1, deltaX);
	  ((TH2F*)fResiduals->UncheckedAt(kResidualPerCh_ClusterIn+1))->Fill(chId+1, deltaY);
	  ((TH2F*)fResiduals->UncheckedAt(kResidualPerHalfCh_ClusterIn))->Fill(halfChId+1, deltaX);
	  ((TH2F*)fResiduals->UncheckedAt(kResidualPerHalfCh_ClusterIn+1))->Fill(halfChId+1, deltaY);
	  ((TH2F*)fResiduals->UncheckedAt(kResidualPerDE_ClusterIn))->Fill(fDEIndices[deId], deltaX);
	  ((TH2F*)fResiduals->UncheckedAt(kResidualPerDE_ClusterIn+1))->Fill(fDEIndices[deId], deltaY);
	  ((TH2F*)fResidualsVsP->UncheckedAt(kResidualInChVsP_ClusterIn+chId))->Fill(p, deltaX);
	  ((TH2F*)fResidualsVsP->UncheckedAt(kResidualInChVsP_ClusterIn+10+chId))->Fill(p, deltaY);
	  
	  // find the track parameters closest to the current cluster position
	  Double_t dZWithPrevious = (previousTrackParam) ? TMath::Abs(previousTrackParam->GetClusterPtr()->GetZ() - cluster->GetZ()) : FLT_MAX;
	  Int_t previousChId = (previousTrackParam) ? previousTrackParam->GetClusterPtr()->GetChamberId() : -1;
	  Double_t dZWithNext = (nextTrackParam) ? TMath::Abs(nextTrackParam->GetClusterPtr()->GetZ() - cluster->GetZ()) : FLT_MAX;
	  AliMUONTrackParam* startingTrackParam = 0x0;
	  if ((fExtrapMode == 0 && dZWithPrevious < dZWithNext) ||
	      (fExtrapMode == 1 && previousTrackParam && !(chId/2 == 2 && previousChId/2 == 1) &&
	       !(chId/2 == 3 && previousChId/2 == 2))) startingTrackParam = previousTrackParam;
	  else startingTrackParam = nextTrackParam;
	  
	  // reset current parameters
	  currentTrackParam.SetParameters(startingTrackParam->GetParameters());
	  currentTrackParam.SetZ(startingTrackParam->GetZ());
	  currentTrackParam.SetCovariances(startingTrackParam->GetCovariances());
	  currentTrackParam.ResetPropagator();
	  
	  // extrapolate to the current cluster position and fill histograms of residuals if everything goes fine
	  if (AliMUONTrackExtrap::ExtrapToZCov(&currentTrackParam, currentTrackParam.GetClusterPtr()->GetZ(), kTRUE)) {
	    
	    // compute MCS dispersion on the first chamber
	    TMatrixD mcsCov(5,5);
	    if (startingTrackParam == nextTrackParam && chId == 0) {
	      AliMUONTrackParam trackParamForMCS;
	      trackParamForMCS.SetParameters(nextTrackParam->GetParameters());
	      AliMUONTrackExtrap::AddMCSEffect(&trackParamForMCS,AliMUONConstants::ChamberThicknessInX0(nextTrackParam->GetClusterPtr()->GetChamberId()),-1.);
	      const TMatrixD &propagator = currentTrackParam.GetPropagator();
	      TMatrixD tmp(trackParamForMCS.GetCovariances(),TMatrixD::kMultTranspose,propagator);
	      mcsCov.Mult(propagator,tmp);
	    } else mcsCov.Zero();
	    
	    // compute residuals
	    Double_t trackResX2 = currentTrackParam.GetCovariances()(0,0) + mcsCov(0,0);
	    Double_t trackResY2 = currentTrackParam.GetCovariances()(2,2) + mcsCov(2,2);
	    deltaX = cluster->GetX() - currentTrackParam.GetNonBendingCoor();
	    deltaY = cluster->GetY() - currentTrackParam.GetBendingCoor();
	    
	    // fill histograms with cluster not attached to the track
	    ((TH2F*)fResiduals->UncheckedAt(kResidualPerCh_ClusterOut))->Fill(chId+1, deltaX);
	    ((TH2F*)fResiduals->UncheckedAt(kResidualPerCh_ClusterOut+1))->Fill(chId+1, deltaY);
	    ((TH2F*)fResiduals->UncheckedAt(kResidualPerHalfCh_ClusterOut))->Fill(halfChId+1, deltaX);
	    ((TH2F*)fResiduals->UncheckedAt(kResidualPerHalfCh_ClusterOut+1))->Fill(halfChId+1, deltaY);
	    ((TH2F*)fResiduals->UncheckedAt(kResidualPerDE_ClusterOut))->Fill(fDEIndices[deId], deltaX);
	    ((TH2F*)fResiduals->UncheckedAt(kResidualPerDE_ClusterOut+1))->Fill(fDEIndices[deId], deltaY);
	    ((TH2F*)fResidualsVsP->UncheckedAt(kResidualInChVsP_ClusterOut+chId))->Fill(p, deltaX);
	    ((TH2F*)fResidualsVsP->UncheckedAt(kResidualInChVsP_ClusterOut+10+chId))->Fill(p, deltaY);
	    ((TH2F*)fResiduals->UncheckedAt(kTrackResPerCh))->Fill(chId+1, TMath::Sqrt(trackResX2));
	    ((TH2F*)fResiduals->UncheckedAt(kTrackResPerCh+1))->Fill(chId+1, TMath::Sqrt(trackResY2));
	    ((TH2F*)fResiduals->UncheckedAt(kTrackResPerHalfCh))->Fill(halfChId+1, TMath::Sqrt(trackResX2));
	    ((TH2F*)fResiduals->UncheckedAt(kTrackResPerHalfCh+1))->Fill(halfChId+1, TMath::Sqrt(trackResY2));
	    ((TH2F*)fResiduals->UncheckedAt(kTrackResPerDE))->Fill(fDEIndices[deId], TMath::Sqrt(trackResX2));
	    ((TH2F*)fResiduals->UncheckedAt(kTrackResPerDE+1))->Fill(fDEIndices[deId], TMath::Sqrt(trackResY2));
	    ((TH2F*)fResiduals->UncheckedAt(kMCSPerCh))->Fill(chId+1, TMath::Sqrt(mcsCov(0,0)));
	    ((TH2F*)fResiduals->UncheckedAt(kMCSPerCh+1))->Fill(chId+1, TMath::Sqrt(mcsCov(2,2)));
	    ((TH2F*)fResiduals->UncheckedAt(kMCSPerHalfCh))->Fill(halfChId+1, TMath::Sqrt(mcsCov(0,0)));
	    ((TH2F*)fResiduals->UncheckedAt(kMCSPerHalfCh+1))->Fill(halfChId+1, TMath::Sqrt(mcsCov(2,2)));
	    ((TH2F*)fResiduals->UncheckedAt(kMCSPerDE))->Fill(fDEIndices[deId], TMath::Sqrt(mcsCov(0,0)));
	    ((TH2F*)fResiduals->UncheckedAt(kMCSPerDE+1))->Fill(fDEIndices[deId], TMath::Sqrt(mcsCov(2,2)));
	    ((TH2F*)fResiduals->UncheckedAt(kClusterRes2PerCh))->Fill(chId+1, deltaX*deltaX - trackResX2);
	    ((TH2F*)fResiduals->UncheckedAt(kClusterRes2PerCh+1))->Fill(chId+1, deltaY*deltaY - trackResY2);
	  }
	  
	}
	
      }
      
      // restore the track
      track.AddTrackParamAtCluster(currentTrackParam, *(currentTrackParam.GetClusterPtr()), kTRUE);
      
    }
    
  }
  
  // Post final data. It will be written to a file with option "RECREATE"
  PostData(1, fResiduals);
  PostData(2, fResidualsVsP);
}

//________________________________________________________________________
void AliAnalysisTaskMuonChamberResolution::NotifyRun()
{
  /// load necessary data from OCDB corresponding to the first run number and initialize analysis
  
  if (fOCDBLoaded) return;
  
  AliCDBManager* cdbm = AliCDBManager::Instance();
  cdbm->SetRun(fCurrentRunNumber);
  
  if (!AliMUONCDB::LoadField()) return;
  
  if (!AliMUONCDB::LoadMapping()) return;
  
  AliMUONRecoParam* recoParam = AliMUONCDB::LoadRecoParam();
  if (!recoParam) return;
  
  fOCDBLoaded = kTRUE;
  
  AliMUONESDInterface::ResetTracker(recoParam);
  
  for (Int_t i = 0; i < AliMUONConstants::NTrackingCh(); i++) {
    
    // set the cluster resolution to default if not already set and create output objects
    if (fClusterResNB[i] < 0.) fClusterResNB[i] = recoParam->GetDefaultNonBendingReso(i);
    if (fClusterResB[i] < 0.) fClusterResB[i] = recoParam->GetDefaultBendingReso(i);
    
    // fill correspondence between DEId and indices for histo (starting from 1)
    AliMpDEIterator it;
    it.First(i);
    while (!it.IsDone()) {
      fNDE++;
      fDEIndices[it.CurrentDEId()] = fNDE;
      fDEIds[fNDE] = it.CurrentDEId();
      it.Next();
    }
    
  }
  
  if (fReAlign) {
    
    // recover default storage name
    TString defaultStorage(cdbm->GetDefaultStorage()->GetType());
    if (defaultStorage == "alien") defaultStorage += Form("://folder=%s", cdbm->GetDefaultStorage()->GetBaseFolder().Data());
    else defaultStorage += Form("://%s", cdbm->GetDefaultStorage()->GetBaseFolder().Data());
    
    // reset existing geometry/alignment if any
    if (cdbm->GetEntryCache()->Contains("GRP/Geometry/Data")) cdbm->UnloadFromCache("GRP/Geometry/Data");
    if (cdbm->GetEntryCache()->Contains("MUON/Align/Data")) cdbm->UnloadFromCache("MUON/Align/Data");
    if (AliGeomManager::GetGeometry()) AliGeomManager::GetGeometry()->UnlockGeometry();
    
    // get original geometry transformer
    AliGeomManager::LoadGeometry();
    if (!AliGeomManager::GetGeometry()) return;
    if (fOldAlignStorage != "none") {
      if (!fOldAlignStorage.IsNull()) cdbm->SetSpecificStorage("MUON/Align/Data",fOldAlignStorage.Data());
      else cdbm->SetSpecificStorage("MUON/Align/Data",defaultStorage.Data());
      AliGeomManager::ApplyAlignObjsFromCDB("MUON");
    }
    fOldGeoTransformer = new AliMUONGeometryTransformer();
    fOldGeoTransformer->LoadGeometryData();
    cout << "Old" << endl;
    fOldGeoTransformer->GetModuleTransformer(0)->GetTransformation()->Print();
    
    // get new geometry transformer
    cdbm->UnloadFromCache("GRP/Geometry/Data");
    if (fOldAlignStorage != "none") cdbm->UnloadFromCache("MUON/Align/Data");
    AliGeomManager::GetGeometry()->UnlockGeometry();
    AliGeomManager::LoadGeometry();
    if (!AliGeomManager::GetGeometry()) return;
    if (!fNewAlignStorage.IsNull()) cdbm->SetSpecificStorage("MUON/Align/Data",fNewAlignStorage.Data());
    else cdbm->SetSpecificStorage("MUON/Align/Data",defaultStorage.Data());
    AliGeomManager::ApplyAlignObjsFromCDB("MUON");
    fNewGeoTransformer = new AliMUONGeometryTransformer();
    fNewGeoTransformer->LoadGeometryData();
    cout << "New" << endl;
    fNewGeoTransformer->GetModuleTransformer(0)->GetTransformation()->Print();
    
  }
  
  // print starting chamber resolution
  printf("\nstarting chamber resolution:\n");
  printf(" - non-bending:");
  for (Int_t i = 0; i < AliMUONConstants::NTrackingCh(); i++) printf((i==0)?" %5.3f":", %5.3f",fClusterResNB[i]);
  printf("\n -     bending:");
  for (Int_t i = 0; i < AliMUONConstants::NTrackingCh(); i++) printf((i==0)?" %6.4f":", %6.4f",fClusterResB[i]);
  printf("\n\n");
  
  UserCreateOutputObjects();
  
}

//________________________________________________________________________
void AliAnalysisTaskMuonChamberResolution::Terminate(Option_t *)
{
  /// compute final results
  
  // recover output objects
  fResiduals = static_cast<TObjArray*> (GetOutputData(1));
  fResidualsVsP = static_cast<TObjArray*> (GetOutputData(2));
  if (!fResiduals || !fResidualsVsP) return;
  
  // summary graphs
  fSummary = new TObjArray(1000);
  fSummary->SetOwner();
  TGraphErrors* g;
  TMultiGraph* mg;
  
  const char* axes[2] = {"X", "Y"};
  Double_t newClusterRes[2][10], newClusterResErr[2][10];
  fNDE = ((TH2F*)fResiduals->UncheckedAt(kResidualPerDE_ClusterIn))->GetXaxis()->GetNbins();
  
  for (Int_t ia = 0; ia < 2; ia++) {
    
    g = new TGraphErrors(AliMUONConstants::NTrackingCh());
    g->SetName(Form("gResidual%sPerChMean_ClusterIn",axes[ia]));
    g->SetTitle(Form("cluster-track residual-%s per Ch: mean (cluster in);chamber ID;<#Delta_{%s}> (cm)",axes[ia],axes[ia]));
    g->SetMarkerStyle(kFullDotLarge);
    fSummary->AddAtAndExpand(g, kResidualPerChMean_ClusterIn+ia);
    g = new TGraphErrors(AliMUONConstants::NTrackingCh());
    g->SetName(Form("gResidual%sPerChMean_ClusterOut",axes[ia]));
    g->SetTitle(Form("cluster-track residual-%s per Ch: mean (cluster out);chamber ID;<#Delta_{%s}> (cm)",axes[ia],axes[ia]));
    g->SetMarkerStyle(kFullDotLarge);
    fSummary->AddAtAndExpand(g, kResidualPerChMean_ClusterOut+ia);
    
    g = new TGraphErrors(2*AliMUONConstants::NTrackingCh());
    g->SetName(Form("gResidual%sPerHalfChMean_ClusterIn",axes[ia]));
    g->SetTitle(Form("cluster-track residual-%s per half Ch: mean (cluster in);half chamber ID;<#Delta_{%s}> (cm)",axes[ia],axes[ia]));
    g->SetMarkerStyle(kFullDotLarge);
    fSummary->AddAtAndExpand(g, kResidualPerHalfChMean_ClusterIn+ia);
    g = new TGraphErrors(2*AliMUONConstants::NTrackingCh());
    g->SetName(Form("gResidual%sPerHalfChMean_ClusterOut",axes[ia]));
    g->SetTitle(Form("cluster-track residual-%s per half Ch: mean (cluster out);half chamber ID;<#Delta_{%s}> (cm)",axes[ia],axes[ia]));
    g->SetMarkerStyle(kFullDotLarge);
    fSummary->AddAtAndExpand(g, kResidualPerHalfChMean_ClusterOut+ia);
    
    g = new TGraphErrors(fNDE);
    g->SetName(Form("gResidual%sPerDEMean_ClusterIn",axes[ia]));
    g->SetTitle(Form("cluster-track residual-%s per DE: mean (cluster in);DE ID;<#Delta_{%s}> (cm)",axes[ia],axes[ia]));
    g->SetMarkerStyle(kFullDotLarge);
    fSummary->AddAtAndExpand(g, kResidualPerDEMean_ClusterIn+ia);
    g = new TGraphErrors(fNDE);
    g->SetName(Form("gResidual%sPerDEMean_ClusterOut",axes[ia]));
    g->SetTitle(Form("cluster-track residual-%s per DE: mean (cluster out);DE ID;<#Delta_{%s}> (cm)",axes[ia],axes[ia]));
    g->SetMarkerStyle(kFullDotLarge);
    fSummary->AddAtAndExpand(g, kResidualPerDEMean_ClusterOut+ia);
    
    g = new TGraphErrors(AliMUONConstants::NTrackingCh());
    g->SetName(Form("gResidual%sPerChSigma_ClusterIn",axes[ia]));
    g->SetTitle(Form("cluster-track residual-%s per Ch: sigma (cluster in);chamber ID;#sigma_{%s} (cm)",axes[ia],axes[ia]));
    g->SetMarkerStyle(kFullDotLarge);
    fSummary->AddAtAndExpand(g, kResidualPerChSigma_ClusterIn+ia);
    g = new TGraphErrors(AliMUONConstants::NTrackingCh());
    g->SetName(Form("gResidual%sPerChSigma_ClusterOut",axes[ia]));
    g->SetTitle(Form("cluster-track residual-%s per Ch: sigma (cluster out);chamber ID;#sigma_{%s} (cm)",axes[ia],axes[ia]));
    g->SetMarkerStyle(kFullDotLarge);
    fSummary->AddAtAndExpand(g, kResidualPerChSigma_ClusterOut+ia);
    
    g = new TGraphErrors(AliMUONConstants::NTrackingCh());
    g->SetName(Form("gResidual%sPerChDispersion_ClusterOut",axes[ia]));
    g->SetTitle(Form("cluster-track residual-%s per Ch: dispersion (cluster out);chamber ID;#sigma_{%s} (cm)",axes[ia],axes[ia]));
    g->SetMarkerStyle(kFullDotLarge);
    fSummary->AddAtAndExpand(g, kResidualPerChDispersion_ClusterOut+ia);
    
    g = new TGraphErrors(AliMUONConstants::NTrackingCh());
    g->SetName(Form("gCombinedResidual%sPerChSigma",axes[ia]));
    g->SetTitle(Form("combined cluster-track residual-%s per Ch: sigma;chamber ID;#sigma_{%s} (cm)",axes[ia],axes[ia]));
    g->SetMarkerStyle(kFullDotLarge);
    fSummary->AddAtAndExpand(g, kCombinedResidualPerChSigma+ia);
    
    g = new TGraphErrors(2*AliMUONConstants::NTrackingCh());
    g->SetName(Form("gCombinedResidual%sPerHalfChSigma",axes[ia]));
    g->SetTitle(Form("combined cluster-track residual-%s per half Ch: sigma;half chamber ID;#sigma_{%s} (cm)",axes[ia],axes[ia]));
    g->SetMarkerStyle(kFullDotLarge);
    fSummary->AddAtAndExpand(g, kCombinedResidualPerHalfChSigma+ia);
    
    g = new TGraphErrors(fNDE);
    g->SetName(Form("gCombinedResidual%sPerDESigma",axes[ia]));
    g->SetTitle(Form("combined cluster-track residual-%s per DE: sigma;DE ID;#sigma_{%s} (cm)",axes[ia],axes[ia]));
    g->SetMarkerStyle(kFullDotLarge);
    fSummary->AddAtAndExpand(g, kCombinedResidualPerDESigma+ia);
    
    mg = new TMultiGraph(Form("mgCombinedResidual%sSigmaVsP",axes[ia]),Form("cluster %s-resolution per chamber versus momentum;p (GeV/c^{2});#sigma_{%s} (cm)",axes[ia],axes[ia]));
    for (Int_t i = 0; i < AliMUONConstants::NTrackingCh(); i++) {
      g = new TGraphErrors(((TH2F*)fResidualsVsP->UncheckedAt(kResidualInChVsP_ClusterIn+10*ia+i))->GetNbinsX());
      g->SetName(Form("gRes%sVsP_ch%d",axes[ia],i+1));
      g->SetMarkerStyle(kFullDotMedium);
      g->SetMarkerColor(i+1+i/9);
      mg->Add(g,"p");
    }
    fSummary->AddAtAndExpand(mg, kCombinedResidualSigmaVsP+ia);
    
    g = new TGraphErrors(AliMUONConstants::NTrackingCh());
    g->SetName(Form("gTrackRes%sPerCh",axes[ia]));
    g->SetTitle(Form("track <#sigma_{%s}> per Ch;chamber ID;#sigma_{%s} (cm)",axes[ia],axes[ia]));
    g->SetMarkerStyle(kFullDotLarge);
    fSummary->AddAtAndExpand(g, kTrackResPerChMean+ia);
    
    g = new TGraphErrors(AliMUONConstants::NTrackingCh());
    g->SetName(Form("gMCS%sPerCh",axes[ia]));
    g->SetTitle(Form("MCS %s-dispersion of extrapolated track per Ch;chamber ID;#sigma_{%s} (cm)",axes[ia],axes[ia]));
    g->SetMarkerStyle(kFullDotLarge);
    fSummary->AddAtAndExpand(g, kMCSPerChMean+ia);
    
    g = new TGraphErrors(AliMUONConstants::NTrackingCh());
    g->SetName(Form("gClusterRes%sPerCh",axes[ia]));
    g->SetTitle(Form("cluster <#sigma_{%s}> per Ch;chamber ID;#sigma_{%s} (cm)",axes[ia],axes[ia]));
    g->SetMarkerStyle(kFullDotLarge);
    fSummary->AddAtAndExpand(g, kClusterResPerCh+ia);
    
    g = new TGraphErrors(2*AliMUONConstants::NTrackingCh());
    g->SetName(Form("gClusterRes%sPerHalfCh",axes[ia]));
    g->SetTitle(Form("cluster <#sigma_{%s}> per half Ch;half chamber ID;#sigma_{%s} (cm)",axes[ia],axes[ia]));
    g->SetMarkerStyle(kFullDotLarge);
    fSummary->AddAtAndExpand(g, kClusterResPerHalfCh+ia);
    
    g = new TGraphErrors(fNDE);
    g->SetName(Form("gClusterRes%sPerDE",axes[ia]));
    g->SetTitle(Form("cluster <#sigma_{%s}> per DE;DE ID;#sigma_{%s} (cm)",axes[ia],axes[ia]));
    g->SetMarkerStyle(kFullDotLarge);
    fSummary->AddAtAndExpand(g, kClusterResPerDE+ia);
    
    g = new TGraphErrors(AliMUONConstants::NTrackingCh());
    g->SetName(Form("gCalcClusterRes%sPerCh",axes[ia]));
    g->SetTitle(Form("calculated cluster <#sigma_{%s}> per Ch;chamber ID;#sigma_{%s} (cm)",axes[ia],axes[ia]));
    g->SetMarkerStyle(kFullDotLarge);
    fSummary->AddAtAndExpand(g, kCalcClusterResPerCh+ia);
    
    // compute residual mean and dispersion per chamber and half chamber
    Double_t meanIn, meanInErr, meanOut, meanOutErr, sigmaIn, sigmaInErr, sigmaOut, sigmaOutErr;
    Double_t sigmaTrack, sigmaTrackErr, sigmaMCS, sigmaMCSErr, clusterRes, clusterResErr, sigmaCluster, sigmaClusterErr;
    for (Int_t i = 0; i < AliMUONConstants::NTrackingCh(); i++) {
      
      // method 1
      TH1D *tmp = ((TH2F*)fResiduals->UncheckedAt(kResidualPerCh_ClusterIn+ia))->ProjectionY("tmp",i+1,i+1,"e");
      GetMean(tmp, meanIn, meanInErr, (TGraphErrors*)fSummary->UncheckedAt(kResidualPerChMean_ClusterIn+ia), i, i+1);
      GetRMS(tmp, sigmaIn, sigmaInErr, (TGraphErrors*)fSummary->UncheckedAt(kResidualPerChSigma_ClusterIn+ia), i, i+1);
      delete tmp;
      
      tmp = ((TH2F*)fResiduals->UncheckedAt(kResidualPerCh_ClusterOut+ia))->ProjectionY("tmp",i+1,i+1,"e");
      GetMean(tmp, meanOut, meanOutErr, (TGraphErrors*)fSummary->UncheckedAt(kResidualPerChMean_ClusterOut+ia), i, i+1);
      GetRMS(tmp, sigmaOut, sigmaOutErr, (TGraphErrors*)fSummary->UncheckedAt(kResidualPerChSigma_ClusterOut+ia), i, i+1);
      delete tmp;
      
      if (fCorrectForSystematics) {
	sigmaIn = TMath::Sqrt(sigmaIn*sigmaIn + meanIn*meanIn);
	sigmaInErr = (sigmaIn>0) ? TMath::Sqrt(sigmaIn*sigmaIn*sigmaInErr*sigmaInErr + meanIn*meanIn*meanInErr*meanInErr) / sigmaIn : 0.;
	sigmaOut = TMath::Sqrt(sigmaOut*sigmaOut + meanOut*meanOut);
	sigmaOutErr = (sigmaOut>0) ? TMath::Sqrt(sigmaOut*sigmaOut*sigmaOutErr*sigmaOutErr + meanOut*meanOut*meanOutErr*meanOutErr) / sigmaOut : 0.;
      }
      ((TGraphErrors*)fSummary->UncheckedAt(kResidualPerChDispersion_ClusterOut+ia))->SetPoint(i, i+1, sigmaOut);
      ((TGraphErrors*)fSummary->UncheckedAt(kResidualPerChDispersion_ClusterOut+ia))->SetPointError(i, 0., sigmaOutErr);
      
      clusterRes = TMath::Sqrt(sigmaIn*sigmaOut);
//      clusterResErr = (clusterRes > 0.) ? 0.5 * TMath::Sqrt(sigmaInErr*sigmaInErr*sigmaOut*sigmaOut + sigmaIn*sigmaIn*sigmaOutErr*sigmaOutErr) / clusterRes : 0.;
      clusterResErr = TMath::Sqrt(sigmaInErr*sigmaOutErr);
      ((TGraphErrors*)fSummary->UncheckedAt(kCombinedResidualPerChSigma+ia))->SetPoint(i, i+1, clusterRes);
      ((TGraphErrors*)fSummary->UncheckedAt(kCombinedResidualPerChSigma+ia))->SetPointError(i, 0., clusterResErr);
      newClusterRes[ia][i] = clusterRes;
      newClusterResErr[ia][i] = clusterResErr;
      
      // method 2
      tmp = ((TH2F*)fResiduals->UncheckedAt(kTrackResPerCh+ia))->ProjectionY("tmp",i+1,i+1,"e");
      GetMean(tmp, sigmaTrack, sigmaTrackErr, (TGraphErrors*)fSummary->UncheckedAt(kTrackResPerChMean+ia), i, i+1, kFALSE);
      delete tmp;
      
      tmp = ((TH2F*)fResiduals->UncheckedAt(kMCSPerCh+ia))->ProjectionY("tmp",i+1,i+1,"e");
      GetMean(tmp, sigmaMCS, sigmaMCSErr, (TGraphErrors*)fSummary->UncheckedAt(kMCSPerChMean+ia), i, i+1, kFALSE);
      delete tmp;
      
      sigmaCluster = sigmaOut*sigmaOut - sigmaTrack*sigmaTrack;
      if (sigmaCluster > 0.) {
	sigmaCluster = TMath::Sqrt(sigmaCluster);
	sigmaClusterErr = TMath::Sqrt(sigmaOut*sigmaOut*sigmaOutErr*sigmaOutErr + sigmaTrack*sigmaTrack*sigmaTrackErr*sigmaTrackErr) / sigmaCluster;
      } else {
	sigmaCluster = 0.;
	sigmaClusterErr = 0.;
      }
      ((TGraphErrors*)fSummary->UncheckedAt(kClusterResPerCh+ia))->SetPoint(i, i+1, sigmaCluster);
      ((TGraphErrors*)fSummary->UncheckedAt(kClusterResPerCh+ia))->SetPointError(i, 0., sigmaClusterErr);
      
      // method 3
      tmp = ((TH2F*)fResiduals->UncheckedAt(kClusterRes2PerCh+ia))->ProjectionY("tmp",i+1,i+1,"e");
      ZoomRight(tmp);
      clusterRes = tmp->GetMean();
      if (clusterRes > 0.) {
	((TGraphErrors*)fSummary->UncheckedAt(kCalcClusterResPerCh+ia))->SetPoint(i, i+1, TMath::Sqrt(clusterRes));
	((TGraphErrors*)fSummary->UncheckedAt(kCalcClusterResPerCh+ia))->SetPointError(i, 0., 0.5 * tmp->GetMeanError() / TMath::Sqrt(clusterRes));
      } else {
	((TGraphErrors*)fSummary->UncheckedAt(kCalcClusterResPerCh+ia))->SetPoint(i, i+1, 0.);
	((TGraphErrors*)fSummary->UncheckedAt(kCalcClusterResPerCh+ia))->SetPointError(i, 0., 0.);
      }
      delete tmp;
      
      // method 1 versus p
      FillSigmaClusterVsP((TH2F*)fResidualsVsP->UncheckedAt(kResidualInChVsP_ClusterIn+10*ia+i),
			  (TH2F*)fResidualsVsP->UncheckedAt(kResidualInChVsP_ClusterOut+10*ia+i),
			  (TGraphErrors*)((TMultiGraph*)fSummary->UncheckedAt(kCombinedResidualSigmaVsP+ia))->GetListOfGraphs()->FindObject(Form("gRes%sVsP_ch%d",axes[ia],i+1)));
      
      // compute residual mean and dispersion per half chamber
      for (Int_t j = 0; j < 2; j++) {
	Int_t k = 2*i+j;
	
	// method 1
	tmp = ((TH2F*)fResiduals->UncheckedAt(kResidualPerHalfCh_ClusterIn+ia))->ProjectionY("tmp",k+1,k+1,"e");
	GetMean(tmp, meanIn, meanInErr, (TGraphErrors*)fSummary->UncheckedAt(kResidualPerHalfChMean_ClusterIn+ia), k, k+1);
	GetRMS(tmp, sigmaIn, sigmaInErr);
	delete tmp;
	
	tmp = ((TH2F*)fResiduals->UncheckedAt(kResidualPerHalfCh_ClusterOut+ia))->ProjectionY("tmp",k+1,k+1,"e");
	GetMean(tmp, meanOut, meanOutErr, (TGraphErrors*)fSummary->UncheckedAt(kResidualPerHalfChMean_ClusterOut+ia), k, k+1);
	GetRMS(tmp, sigmaOut, sigmaOutErr);
	delete tmp;
	
	if (fCorrectForSystematics) {
	  sigmaIn = TMath::Sqrt(sigmaIn*sigmaIn + meanIn*meanIn);
	  sigmaInErr = (sigmaIn>0) ? TMath::Sqrt(sigmaIn*sigmaIn*sigmaInErr*sigmaInErr + meanIn*meanIn*meanInErr*meanInErr) / sigmaIn : 0.;
	  sigmaOut = TMath::Sqrt(sigmaOut*sigmaOut + meanOut*meanOut);
	  sigmaOutErr = (sigmaOut>0) ? TMath::Sqrt(sigmaOut*sigmaOut*sigmaOutErr*sigmaOutErr + meanOut*meanOut*meanOutErr*meanOutErr) / sigmaOut : 0.;
	}
	
	clusterRes = TMath::Sqrt(sigmaIn*sigmaOut);
//	clusterResErr = (clusterRes > 0.) ? 0.5 * TMath::Sqrt(sigmaInErr*sigmaInErr*sigmaOut*sigmaOut + sigmaIn*sigmaIn*sigmaOutErr*sigmaOutErr) / clusterRes : 0.;
	clusterResErr = TMath::Sqrt(sigmaInErr*sigmaOutErr);
	((TGraphErrors*)fSummary->UncheckedAt(kCombinedResidualPerHalfChSigma+ia))->SetPoint(k, k+1, clusterRes);
	((TGraphErrors*)fSummary->UncheckedAt(kCombinedResidualPerHalfChSigma+ia))->SetPointError(k, 0., clusterResErr);
	
	// method 2
	tmp = ((TH2F*)fResiduals->UncheckedAt(kTrackResPerHalfCh+ia))->ProjectionY("tmp",k+1,k+1,"e");
	GetMean(tmp, sigmaTrack, sigmaTrackErr, 0x0, 0, 0, kFALSE);
	delete tmp;
	
	tmp = ((TH2F*)fResiduals->UncheckedAt(kMCSPerHalfCh+ia))->ProjectionY("tmp",k+1,k+1,"e");
	GetMean(tmp, sigmaMCS, sigmaMCSErr, 0x0, 0, 0, kFALSE);
	delete tmp;
	
	sigmaCluster = sigmaOut*sigmaOut - sigmaTrack*sigmaTrack;
	if (sigmaCluster > 0.) {
	  sigmaCluster = TMath::Sqrt(sigmaCluster);
	  sigmaClusterErr = TMath::Sqrt(sigmaOut*sigmaOut*sigmaOutErr*sigmaOutErr + sigmaTrack*sigmaTrack*sigmaTrackErr*sigmaTrackErr) / sigmaCluster;
	} else {
	  sigmaCluster = 0.;
	  sigmaClusterErr = 0.;
	}
	((TGraphErrors*)fSummary->UncheckedAt(kClusterResPerHalfCh+ia))->SetPoint(k, k+1, sigmaCluster);
	((TGraphErrors*)fSummary->UncheckedAt(kClusterResPerHalfCh+ia))->SetPointError(k, 0., sigmaClusterErr);
	
      }
      
    }
    
    // compute residual mean and dispersion per DE
    for (Int_t i = 0; i < fNDE; i++) {
      
      // method 1
      TH1D *tmp = ((TH2F*)fResiduals->UncheckedAt(kResidualPerDE_ClusterIn+ia))->ProjectionY("tmp",i+1,i+1,"e");
      GetMean(tmp, meanIn, meanInErr, (TGraphErrors*)fSummary->UncheckedAt(kResidualPerDEMean_ClusterIn+ia), i, i+1);
      GetRMS(tmp, sigmaIn, sigmaInErr);
      delete tmp;
      
      tmp = ((TH2F*)fResiduals->UncheckedAt(kResidualPerDE_ClusterOut+ia))->ProjectionY("tmp",i+1,i+1,"e");
      GetMean(tmp, meanOut, meanOutErr, (TGraphErrors*)fSummary->UncheckedAt(kResidualPerDEMean_ClusterOut+ia), i, i+1);
      GetRMS(tmp, sigmaOut, sigmaOutErr);
      delete tmp;
      
      if (fCorrectForSystematics) {
	sigmaIn = TMath::Sqrt(sigmaIn*sigmaIn + meanIn*meanIn);
	sigmaInErr = (sigmaIn>0) ? TMath::Sqrt(sigmaIn*sigmaIn*sigmaInErr*sigmaInErr + meanIn*meanIn*meanInErr*meanInErr) / sigmaIn : 0.;
	sigmaOut = TMath::Sqrt(sigmaOut*sigmaOut + meanOut*meanOut);
	sigmaOutErr = (sigmaOut>0) ? TMath::Sqrt(sigmaOut*sigmaOut*sigmaOutErr*sigmaOutErr + meanOut*meanOut*meanOutErr*meanOutErr) / sigmaOut : 0.;
      }
      
      clusterRes = TMath::Sqrt(sigmaIn*sigmaOut);
//      clusterResErr = (clusterRes > 0.) ? 0.5 * TMath::Sqrt(sigmaInErr*sigmaInErr*sigmaOut*sigmaOut + sigmaIn*sigmaIn*sigmaOutErr*sigmaOutErr) / clusterRes : 0.;
      clusterResErr = TMath::Sqrt(sigmaInErr*sigmaOutErr);
      ((TGraphErrors*)fSummary->UncheckedAt(kCombinedResidualPerDESigma+ia))->SetPoint(i, i+1, clusterRes);
      ((TGraphErrors*)fSummary->UncheckedAt(kCombinedResidualPerDESigma+ia))->SetPointError(i, 0., clusterResErr);
      
      // method 2
      tmp = ((TH2F*)fResiduals->UncheckedAt(kTrackResPerDE+ia))->ProjectionY("tmp",i+1,i+1,"e");
      GetMean(tmp, sigmaTrack, sigmaTrackErr, 0x0, 0, 0, kFALSE);
      delete tmp;
      
      tmp = ((TH2F*)fResiduals->UncheckedAt(kMCSPerDE+ia))->ProjectionY("tmp",i+1,i+1,"e");
      GetMean(tmp, sigmaMCS, sigmaMCSErr, 0x0, 0, 0, kFALSE);
      delete tmp;
      
      sigmaCluster = sigmaOut*sigmaOut - sigmaTrack*sigmaTrack;
      if (sigmaCluster > 0.) {
	sigmaCluster = TMath::Sqrt(sigmaCluster);
	sigmaClusterErr = TMath::Sqrt(sigmaOut*sigmaOut*sigmaOutErr*sigmaOutErr + sigmaTrack*sigmaTrack*sigmaTrackErr*sigmaTrackErr) / sigmaCluster;
      } else {
	sigmaCluster = 0.;
	sigmaClusterErr = 0.;
      }
      ((TGraphErrors*)fSummary->UncheckedAt(kClusterResPerDE+ia))->SetPoint(i, i+1, sigmaCluster);
      ((TGraphErrors*)fSummary->UncheckedAt(kClusterResPerDE+ia))->SetPointError(i, 0., sigmaClusterErr);
      
    }
    
    // set graph labels
    TAxis* xAxis = ((TH2F*)fResiduals->UncheckedAt(kResidualPerDE_ClusterOut+ia))->GetXaxis();
    ((TGraphErrors*)fSummary->UncheckedAt(kResidualPerDEMean_ClusterIn+ia))->GetXaxis()->Set(fNDE, 0.5, fNDE+0.5);
    ((TGraphErrors*)fSummary->UncheckedAt(kResidualPerDEMean_ClusterOut+ia))->GetXaxis()->Set(fNDE, 0.5, fNDE+0.5);
    ((TGraphErrors*)fSummary->UncheckedAt(kCombinedResidualPerDESigma+ia))->GetXaxis()->Set(fNDE, 0.5, fNDE+0.5);
    ((TGraphErrors*)fSummary->UncheckedAt(kClusterResPerDE+ia))->GetXaxis()->Set(fNDE, 0.5, fNDE+0.5);
    for (Int_t i = 1; i <= fNDE; i++) {
      const char* label = xAxis->GetBinLabel(i);
      ((TGraphErrors*)fSummary->UncheckedAt(kResidualPerDEMean_ClusterIn+ia))->GetXaxis()->SetBinLabel(i, label);
      ((TGraphErrors*)fSummary->UncheckedAt(kResidualPerDEMean_ClusterOut+ia))->GetXaxis()->SetBinLabel(i, label);
      ((TGraphErrors*)fSummary->UncheckedAt(kCombinedResidualPerDESigma+ia))->GetXaxis()->SetBinLabel(i, label);
      ((TGraphErrors*)fSummary->UncheckedAt(kClusterResPerDE+ia))->GetXaxis()->SetBinLabel(i, label);
    }
    
  }
  
  // display
  fCanvases = new TObjArray(1000);
  fCanvases->SetOwner();
  TCanvas* cResPerCh = new TCanvas("cResPerCh","cResPerCh",1200,500);
  cResPerCh->Divide(4,2);
  for (Int_t ia = 0; ia < 2; ia++) {
    cResPerCh->cd(1+4*ia);
    g = (TGraphErrors*)fSummary->UncheckedAt(kResidualPerChMean_ClusterOut+ia);
    g->Draw("ap");
    g->SetMarkerColor(2);
    g->SetLineColor(2);
    g = (TGraphErrors*)fSummary->UncheckedAt(kResidualPerChMean_ClusterIn+ia);
    g->Draw("p");
    g->SetMarkerColor(4);
    g->SetLineColor(4);
    cResPerCh->cd(2+4*ia);
    g = (TGraphErrors*)fSummary->UncheckedAt(kResidualPerChSigma_ClusterOut+ia);
    g->Draw("ap");
    g->SetMinimum(0.);
    g->SetMarkerColor(2);
    g->SetLineColor(2);
    g = (TGraphErrors*)fSummary->UncheckedAt(kResidualPerChSigma_ClusterIn+ia);
    g->Draw("p");
    g->SetMarkerColor(4);
    g->SetLineColor(4);
    g = (TGraphErrors*)fSummary->UncheckedAt(kMCSPerChMean+ia);
    g->Draw("p");
    g->SetMarkerColor(5);
    g->SetLineColor(5);
    g = (TGraphErrors*)fSummary->UncheckedAt(kCombinedResidualPerChSigma+ia);
    g->Draw("p");
    g->SetMarkerColor(3);
    g->SetLineColor(3);
    cResPerCh->cd(3+4*ia);
    g = (TGraphErrors*)fSummary->UncheckedAt(kResidualPerChDispersion_ClusterOut+ia);
    g->Draw("ap");
    g->SetMinimum(0.);
    g->SetMarkerColor(2);
    g->SetLineColor(2);
    g = (TGraphErrors*)fSummary->UncheckedAt(kMCSPerChMean+ia);
    g->Draw("p");
    g = (TGraphErrors*)fSummary->UncheckedAt(kTrackResPerChMean+ia);
    g->Draw("p");
    g->SetMarkerColor(4);
    g->SetLineColor(4);
    g = (TGraphErrors*)fSummary->UncheckedAt(kClusterResPerCh+ia);
    g->Draw("p");
    cResPerCh->cd(4+4*ia);
    g = (TGraphErrors*)fSummary->UncheckedAt(kCombinedResidualPerChSigma+ia);
    g->Draw("ap");
    g->SetMinimum(0.);
    g = (TGraphErrors*)fSummary->UncheckedAt(kClusterResPerCh+ia);
    g->Draw("p");
    g = (TGraphErrors*)fSummary->UncheckedAt(kCalcClusterResPerCh+ia);
    g->Draw("p");
    g->SetMarkerColor(6);
    g->SetLineColor(6);
  }
  fCanvases->AddAtAndExpand(cResPerCh, kResPerCh);
  
  TCanvas* cResPerHalfCh = new TCanvas("cResPerHalfCh","cResPerHalfCh",1200,500);
  cResPerHalfCh->Divide(2,2);
  for (Int_t ia = 0; ia < 2; ia++) {
    cResPerHalfCh->cd(1+2*ia);
    g = (TGraphErrors*)fSummary->UncheckedAt(kResidualPerHalfChMean_ClusterOut+ia);
    g->Draw("ap");
    g->SetMarkerColor(2);
    g->SetLineColor(2);
    g = (TGraphErrors*)fSummary->UncheckedAt(kResidualPerHalfChMean_ClusterIn+ia);
    g->Draw("p");
    g->SetMarkerColor(4);
    g->SetLineColor(4);
    cResPerHalfCh->cd(2+2*ia);
    g = (TGraphErrors*)fSummary->UncheckedAt(kCombinedResidualPerHalfChSigma+ia);
    g->Draw("ap");
    g->SetMinimum(0.);
    g->SetMarkerColor(3);
    g->SetLineColor(3);
    g = (TGraphErrors*)fSummary->UncheckedAt(kClusterResPerHalfCh+ia);
    g->Draw("p");
  }
  fCanvases->AddAtAndExpand(cResPerHalfCh, kResPerHalfCh);
  
  TCanvas* cResPerDE = new TCanvas("cResPerDE","cResPerDE",1200,800);
  cResPerDE->Divide(1,4);
  for (Int_t ia = 0; ia < 2; ia++) {
    cResPerDE->cd(1+ia);
    g = (TGraphErrors*)fSummary->UncheckedAt(kResidualPerDEMean_ClusterOut+ia);
    g->Draw("ap");
    g->SetMarkerColor(2);
    g->SetLineColor(2);
    g = (TGraphErrors*)fSummary->UncheckedAt(kResidualPerDEMean_ClusterIn+ia);
    g->Draw("p");
    g->SetMarkerColor(4);
    g->SetLineColor(4);
    cResPerDE->cd(3+ia);
    g = (TGraphErrors*)fSummary->UncheckedAt(kCombinedResidualPerDESigma+ia);
    g->Draw("ap");
    g->SetMinimum(0.);
    g->SetMarkerColor(3);
    g->SetLineColor(3);
    g = (TGraphErrors*)fSummary->UncheckedAt(kClusterResPerDE+ia);
    g->Draw("p");
  }
  fCanvases->AddAtAndExpand(cResPerDE, kResPerDE);
  
  TCanvas* cResPerChVsP = new TCanvas("cResPerChVsP","cResPerChVsP");
  cResPerChVsP->Divide(1,2);
  for (Int_t ia = 0; ia < 2; ia++) {
    cResPerChVsP->cd(1+ia);
    mg = (TMultiGraph*)fSummary->UncheckedAt(kCombinedResidualSigmaVsP+ia);
    mg->Draw("ap");
  }
  fCanvases->AddAtAndExpand(cResPerChVsP, kResPerChVsP);
  
  // print results
  printf("\nchamber resolution:\n");
  printf(" - non-bending:");
  for (Int_t i = 0; i < AliMUONConstants::NTrackingCh(); i++) printf((i==0)?" %5.3f":", %5.3f",newClusterRes[0][i]);
  printf("\n -     bending:");
  for (Int_t i = 0; i < AliMUONConstants::NTrackingCh(); i++) printf((i==0)?" %6.4f":", %6.4f",newClusterRes[1][i]);
  printf("\n\n");
  
  Double_t x, y;
  printf("\nhalf chamber systematic shifts:\n");
  printf(" - non-bending:");
  for (Int_t i = 0; i < 20; i++) {
    ((TGraphErrors*)fSummary->UncheckedAt(kResidualPerHalfChMean_ClusterIn))->GetPoint(i,x,y);
    printf((i==0)?" %5.3f":", %5.3f",y);
  }
  printf("\n -     bending:");
  for (Int_t i = 0; i < 20; i++) {
    ((TGraphErrors*)fSummary->UncheckedAt(kResidualPerHalfChMean_ClusterIn+1))->GetPoint(i,x,y);
    printf((i==0)?" %6.4f":", %6.4f",y);
  }
  printf("\n\n");
  
  // Post final data.
  PostData(3, fSummary);
}

//________________________________________________________________________
void AliAnalysisTaskMuonChamberResolution::ModifyClusters(AliMUONTrack& track)
{
  /// Reset the clusters resolution from the ones given to the task and change
  /// the cluster position according to the new alignment parameters if required
  
  Double_t gX,gY,gZ,lX,lY,lZ;
  
  // loop over clusters
  Int_t nClusters = track.GetNClusters();
  for (Int_t iCluster=0; iCluster<nClusters; iCluster++) {
    
    AliMUONVCluster* cl = static_cast<AliMUONTrackParam*>(track.GetTrackParamAtCluster()->UncheckedAt(iCluster))->GetClusterPtr();
    
    // change their resolution
    cl->SetErrXY(fClusterResNB[cl->GetChamberId()], fClusterResB[cl->GetChamberId()]);
    
    // change their position
    if (fReAlign) {
      gX = cl->GetX();
      gY = cl->GetY();
      gZ = cl->GetZ();
      //fOldGeoTransformer->GetModuleTransformerByDEId(cl->GetDetElemId())->GetTransformation()->Print();
      //fOldGeoTransformer->GetModuleTransformerByDEId(cl->GetDetElemId())->GetDetElement(cl->GetDetElemId())->GetGlobalTransformation()->Print();
      //cout<<"chamber "<<cl->GetChamberId()<<" before: "<<gX<<", "<<gY<<", "<<gZ<<endl;
      fOldGeoTransformer->Global2Local(cl->GetDetElemId(),gX,gY,gZ,lX,lY,lZ);
      fNewGeoTransformer->Local2Global(cl->GetDetElemId(),lX,lY,lZ,gX,gY,gZ);
      cl->SetXYZ(gX,gY,gZ);
      //fNewGeoTransformer->GetModuleTransformerByDEId(cl->GetDetElemId())->GetTransformation()->Print();
      //fNewGeoTransformer->GetModuleTransformerByDEId(cl->GetDetElemId())->GetDetElement(cl->GetDetElemId())->GetGlobalTransformation()->Print();
      //cout<<"          after: "<<gX<<", "<<gY<<", "<<gZ<<endl;
    }
    
  }
  
}

//________________________________________________________________________
void AliAnalysisTaskMuonChamberResolution::Zoom(TH1* h, Double_t fractionCut)
{
  /// Reduce the range of the histogram by removing a given fration of the statistic at each edge
  ZoomLeft(h, fractionCut);
  ZoomRight(h, fractionCut);
}

//________________________________________________________________________
void AliAnalysisTaskMuonChamberResolution::ZoomLeft(TH1* h, Double_t fractionCut)
{
  /// Reduce the range of the histogram by removing a given fration of the statistic on the left side
  Int_t maxEventsCut = (Int_t) (fractionCut * h->GetEntries());
  Int_t nBins = h->GetNbinsX();
  
  // set low edge  
  Int_t minBin;
  Int_t eventsCut = 0;
  for (minBin = 1; minBin <= nBins; minBin++) {
    eventsCut += (Int_t) h->GetBinContent(minBin);
    if (eventsCut > maxEventsCut) break;
  }
  
  // set new axis range
  h->GetXaxis()->SetRange(--minBin, h->GetXaxis()->GetLast());
}

//________________________________________________________________________
void AliAnalysisTaskMuonChamberResolution::ZoomRight(TH1* h, Double_t fractionCut)
{
  /// Reduce the range of the histogram by removing a given fration of the statistic on the right side
  Int_t maxEventsCut = (Int_t) (fractionCut * h->GetEntries());
  Int_t nBins = h->GetNbinsX();
  
  // set high edge
  Int_t maxBin;
  Int_t eventsCut = 0;
  for (maxBin = nBins; maxBin >= 1; maxBin--) {
    eventsCut += (Int_t) h->GetBinContent(maxBin);
    if (eventsCut > maxEventsCut) break;
  }
  
  // set new axis range
  h->GetXaxis()->SetRange(h->GetXaxis()->GetFirst(), ++maxBin);
}

//________________________________________________________________________
void AliAnalysisTaskMuonChamberResolution::GetMean(TH1* h, Double_t& mean, Double_t& meanErr, TGraphErrors* g, Int_t i, Double_t x, Bool_t zoom)
{
  /// Fill graph with the mean value of the histogram and the corresponding error (zooming if required)
  Int_t firstBin = h->GetXaxis()->GetFirst();
  Int_t lastBin = h->GetXaxis()->GetLast();
  if (zoom) Zoom(h);
  mean = (h->GetEntries() > fgkMinEntries) ? h->GetMean() : 0.;
  meanErr = (h->GetEntries() > fgkMinEntries) ? h->GetMeanError() : 0.;
  if (g) {
    g->SetPoint(i, x, mean);
    g->SetPointError(i, 0., meanErr);
  }
  if (zoom) h->GetXaxis()->SetRange(firstBin,lastBin);
}

//________________________________________________________________________
void AliAnalysisTaskMuonChamberResolution::GetRMS(TH1* h, Double_t& rms, Double_t& rmsErr, TGraphErrors* g, Int_t i, Double_t x, Bool_t zoom)
{
  /// Return the RMS of the histogram and the corresponding error (zooming if required) and fill graph if !=0x0
  Int_t firstBin = h->GetXaxis()->GetFirst();
  Int_t lastBin = h->GetXaxis()->GetLast();
  if (zoom) Zoom(h);
  rms = (h->GetEntries() > fgkMinEntries) ? h->GetRMS() : 0.;
  rmsErr = (h->GetEntries() > fgkMinEntries) ? h->GetRMSError() : 0.;
  if (g) {
    g->SetPoint(i, x, rms);
    g->SetPointError(i, 0., rmsErr);
  }
  if (zoom) h->GetXaxis()->SetRange(firstBin,lastBin);
}

//________________________________________________________________________
void AliAnalysisTaskMuonChamberResolution::FillSigmaClusterVsP(TH2* hIn, TH2* hOut, TGraphErrors* g, Bool_t zoom)
{
  /// Fill graph with cluster resolution from combined residuals with cluster in/out (zooming if required)
  Double_t sigmaIn, sigmaInErr, sigmaOut, sigmaOutErr, clusterRes, clusterResErr;
  for (Int_t j = 1; j <= hIn->GetNbinsX(); j++) {
    TH1D* tmp = hIn->ProjectionY("tmp",j,j,"e");
    GetRMS(tmp, sigmaIn, sigmaInErr, 0x0, 0, 0., zoom);
    delete tmp;
    tmp = hOut->ProjectionY("tmp",j,j,"e");
    GetRMS(tmp, sigmaOut, sigmaOutErr, 0x0, 0, 0., zoom);
    delete tmp;
    Double_t p = 0.5 * (hIn->GetBinLowEdge(j) + hIn->GetBinLowEdge(j+1));
    Double_t pErr = p - hIn->GetBinLowEdge(j);
    clusterRes = TMath::Sqrt(sigmaIn*sigmaOut);
    clusterResErr = (clusterRes > 0.) ? 0.5 * TMath::Sqrt(sigmaInErr*sigmaInErr*sigmaOut*sigmaOut + sigmaIn*sigmaIn*sigmaOutErr*sigmaOutErr) / clusterRes : 0.;
    g->SetPoint(j, p, clusterRes);
    g->SetPointError(j, pErr, clusterResErr);
  }
}

