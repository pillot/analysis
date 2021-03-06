// $Id: MUON_displaySimu.C 37552 2009-12-03 14:53:09Z ivana $
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

/// \ingroup evemacros
/// \file MUON_displaySimu.C
///
/// \author B. Vulpescu, LPC

class AliEveMUONData;
class AliEveEventManager;

AliEveMUONData     *g_muon_data       = 0;

Int_t  g_currentEvent = -1;
Bool_t g_fromRaw      = kFALSE;


void MUON_displaySimu(Bool_t fromRaw = kFALSE, Bool_t showTracks = kTRUE, Bool_t clustersFromESD = kTRUE)
{
  //
  // display from simulated digits (or produced raw data) 
  // tracks: ESD, Refs, MC
  // 

  if (!AliMpSegmentation::Instance()) AliMpCDB::LoadMpSegmentation();
  if (!AliMpDDLStore::Instance())     AliMpCDB::LoadDDLStore();

  // set the magnetic field for track extrapolations
  AliEveEventManager::AssertMagField();
  
  // prepare to recover track parameters at each cluster
  AliMUONESDInterface::ResetTracker(AliMUONCDB::LoadRecoParam());

  TTree* dt = 0;
  TTree* ct = 0;
  TTree* ht = 0;

  if (AliEveEventManager::GetMaster() == 0) {
    printf("No alieve event: use alieve_init(...) \n");
    return;
  }

  if (g_currentEvent == AliEveEventManager::GetMaster()->GetEventId()) {
    if (g_fromRaw == fromRaw) {
      printf("Same event... \n");
      return;
    } else {
      if (g_fromRaw) {
	printf("Same event with digits.\n");
	AliEveEventManager::GetMaster()->GotoEvent(g_currentEvent);
      } else {
	printf("Same event with raw.\n");
	AliEveEventManager::GetMaster()->GotoEvent(g_currentEvent);
      }
    }
  }

  g_fromRaw = fromRaw;

  TString dataPath = TString(AliEveEventManager::GetMaster()->GetTitle());
  dataPath.Append("/raw.root");

  AliRunLoader* rl =  AliEveEventManager::AssertRunLoader();
  g_muon_data = new AliEveMUONData;

  if (!fromRaw) {
    rl->LoadDigits("MUON");
    dt = rl->GetTreeD("MUON", false);
    if (dt == 0) {
      cout << "No digits produced!" << endl;
    } else {
      cout << "With aliroot digits!" << endl;
      g_muon_data->LoadDigits(dt);
    }
  } else {
    if (gSystem->AccessPathName(dataPath.Data(),kFileExists)) {
      cout << "No raw data produced!" << endl;
    } else {
      cout << "With raw digits!" << endl;
      g_muon_data->LoadRaw(dataPath.Data());
    }
  }

  TString esdDataPath = TString(AliEveEventManager::GetMaster()->GetTitle());
  esdDataPath.Append("/AliESDs.root");
  if (clustersFromESD) {
    g_muon_data->LoadRecPointsFromESD(esdDataPath.Data());
  } else {
    rl->LoadRecPoints("MUON");
    ct = rl->GetTreeR("MUON", false);
    g_muon_data->LoadRecPoints(ct);
  }
/*  
  rl->LoadHits("MUON");
  ht = rl->GetTreeH("MUON", false);
  g_muon_data->LoadHits(ht);
*/  
  g_currentEvent = AliEveEventManager::GetMaster()->GetEventId();

  gStyle->SetPalette(1, 0);

  gEve->DisableRedraw();

  TEveElementList* l = new TEveElementList("MUONChambers");
  l->SetTitle("MUON chambers");
  l->SetMainColor(2);
  gEve->AddElement(l);

  for (Int_t ic = 0; ic < 14; ic++)
  {
    AliEveMUONChamber* mucha = new AliEveMUONChamber(ic);

    mucha->SetFrameColor(2);
    mucha->SetChamberID(ic);

    mucha->SetDataSource(g_muon_data);

    gEve->AddElement(mucha, l);
  }

  if (showTracks) {
    MUON_ESD_tracks();
    MUON_Ref_tracks();
    MUON_MC_tracks();
  }

  gEve->Redraw3D(kTRUE);
  gEve->EnableRedraw();
}

//______________________________________________________________________________
void MUON_ESD_tracks()
{
  AliESDEvent* esd = AliEveEventManager::AssertESD();

  TEveTrackList* lt = new TEveTrackList("ESD-Tracks");
  lt->SetMainColor(6);
  //lt->SetMUON();

  gEve->AddElement(lt);

  AliESDMuonTrack *mt;
  TEveRecTrack rt;
  Int_t nMuonTracks = esd->GetNumberOfMuonTracks();
  Int_t nTrack = 0;
  for (Int_t n = 0; n < nMuonTracks; n++)
  {
    mt = esd->GetMuonTrack(n);

    nTrack++;

    rt.fLabel = n;

    AliEveMUONTrack* track = new AliEveMUONTrack(&rt, lt->GetPropagator());

    track->MakeESDTrack(mt);

    gEve->AddElement(track, lt);
  }

}

//______________________________________________________________________________
void MUON_Ref_tracks()
{
  TString dataPathESD = TString(AliEveEventManager::GetMaster()->GetTitle());
  dataPathESD.Append("/AliESDs.root");
  TString dataPathSIM = TString(AliEveEventManager::GetMaster()->GetTitle());
  dataPathSIM.Append("/");

  AliMUONRecoCheck recoCheck(dataPathESD.Data(),dataPathSIM.Data());
  AliMUONVTrackStore* trackRefStore = recoCheck.TrackRefs(AliEveEventManager::GetMaster()->GetEventId());
  TIter next(trackRefStore->CreateIterator());
  AliMUONTrack* trackRef;

  TEveTrackList* lt = new TEveTrackList("Ref-Tracks");
  lt->SetMainColor(6);

  gEve->AddElement(lt);

  TEveRecTrack rt;
  Int_t i = 0;
  while ( ( trackRef = static_cast<AliMUONTrack*>(next()) ) )
  {
    rt.fLabel = i++;

    AliEveMUONTrack* track = new AliEveMUONTrack(&rt, lt->GetPropagator());

    track->MakeRefTrack(trackRef);

    gEve->AddElement(track, lt);
  }

}

//______________________________________________________________________________
void MUON_MC_tracks()
{
  Double_t RADDEG = 180.0/TMath::Pi();

  AliRunLoader* rl =  AliEveEventManager::AssertRunLoader();
  rl->LoadKinematics();
  AliStack* stack = rl->Stack();

  Int_t nPrimary = stack->GetNprimary();
  Int_t nTracks  = stack->GetNtrack();

  TEveTrackList* lt = new TEveTrackList("MC-Tracks");
  lt->SetMainColor(6);
  //lt->SetMUON();

  gEve->AddElement(lt);

  Int_t pdgCode;
  TParticle *part;
  TEveRecTrack rt;

  Int_t nHitTracks = g_muon_data->GetNTrackList();
  Int_t index;
  for (Int_t i = 0; i < nHitTracks; i++)
  {
    index = g_muon_data->GetTrack(i);
    if (index >= nTracks) {
      cout << "TEveHit track index larger than number in stack!" << endl;
      continue;
    }

    part = stack->Particle(index);
    if (part->P() < 0.001) continue;  // skip momenta < 1.0 MeV/c
    rt.fLabel = i;

    AliEveMUONTrack* track = new AliEveMUONTrack(&rt, lt->GetPropagator());

    track->MakeMCTrack(part);

    gEve->AddElement(track, lt);
  }

  rl->UnloadKinematics();
}

