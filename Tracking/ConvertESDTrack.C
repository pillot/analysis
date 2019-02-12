#include <fstream>
#include <iostream>

#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TGeoGlobalMagField.h>

#include "AliCDBManager.h"
#include "AliGeomManager.h"

#include "AliESDEvent.h"
#include "AliESDMuonTrack.h"
#include "AliESDMuonCluster.h"

#include "AliMUONCDB.h"
#include "AliMUONTrack.h"
#include "AliMUONTrackParam.h"
#include "AliMUONRecoParam.h"
#include "AliMUONESDInterface.h"

#include "Field/MagneticField.h"

#include "MCHBase/ClusterBlock.h"
#include "MCHBase/TrackBlock.h"

using namespace std;
using namespace o2::mch;

void ConvertESDTrack(TString esdFileName, TString outFileName = "AliESDs.in", Bool_t refit = kFALSE, Int_t LastEvent = -1)
{
  /// convert ESD tracks+clusters into O2 structures
  /// saved in a binary file with the following format:
  ///
  /// #tracks in event 1
  /// TrackParamStruct of 1st track
  /// #clusters in track 1
  /// ClusterStruct of 1st cluster
  /// ClusterStruct of 2nd cluster
  /// ...
  /// ClusterStruct of nth cluster
  /// TrackParamStruct of 2nd track
  /// #clusters in track 2
  /// ...
  /// #tracks in event 2
  /// ...

  // open the ESD file
  TFile* esdFile = TFile::Open(esdFileName);
  if (!esdFile || !esdFile->IsOpen()) {
    Error("ConvertESDTrack", "opening ESD file %s failed", esdFileName.Data());
    return;
  }
  AliESDEvent* esd = new AliESDEvent();
  TTree* tree = (TTree*) esdFile->Get("esdTree");
  if (!tree) {
    Error("ConvertESDTrack", "no ESD tree found");
    return;
  }
  esd->ReadFromTree(tree);
  
  // output file
  ofstream out(outFileName.Data(),ios::out | ios::binary);
  if (!out.is_open()) {
    return;
  }
  
  Int_t nevents = (LastEvent >= 0) ? TMath::Min(LastEvent+1, (Int_t)tree->GetEntries()) : (Int_t)tree->GetEntries();

  TrackParamStruct sTrackParam;
  ClusterStruct sCluster;
  AliMUONTrack muonTrack;

  for (Int_t iEvent = 0; iEvent < nevents; iEvent++) {
    
    // get the ESD event
    if (tree->GetEvent(iEvent) <= 0) {
      Error("ConvertESDTrack", "no ESD object found for event %d", iEvent);
      return;
    }
    out.write((char*)&iEvent,sizeof(Int_t));
    
    // need OCDB access for refitting
    if (iEvent == 0 && refit) {
      AliCDBManager::Instance()->SetDefaultStorage("local://./OCDB");
      AliCDBManager::Instance()->SetRun(esd->GetRunNumber());
      auto field =
      o2::field::MagneticField::createFieldMap(-30000., -5999.95, o2::field::MagneticField::kConvLHC, false, 3500., "A-A",
                                               "$(O2_ROOT)/share/Common/maps/mfchebKGI_sym.root");
      TGeoGlobalMagField::Instance()->SetField(field);
      TGeoGlobalMagField::Instance()->Lock();
      AliGeomManager::LoadGeometry();
      if (!AliGeomManager::GetGeometry()) return;
      if (!AliGeomManager::ApplyAlignObjsFromCDB("MUON")) return;
      AliMUONRecoParam* recoParam = AliMUONCDB::LoadRecoParam();
      if (!recoParam) return;
      //recoParam->UseSmoother(false);
      AliMUONESDInterface::ResetTracker(recoParam);
    }

    // count the number of tracks to be stored, excluding ghost,
    // as well as the total number of bytes requested to store
    // the event, excluding event number and total size
    Int_t nTracks = (Int_t)esd->GetNumberOfMuonTracks();
    Int_t nTrkTracks = 0;
    Int_t size = 0;
    for (Int_t iTrack = 0; iTrack < nTracks; iTrack++) {
      AliESDMuonTrack *track = esd->GetMuonTrack(iTrack);
      if (!track->ContainTrackerData()) continue;
      ++nTrkTracks;
      size += sizeof(TrackParamStruct) + sizeof(Int_t) + track->GetNClusters()*sizeof(ClusterStruct);
      
    }
    if (nTrkTracks > 0) size += sizeof(Int_t);
    out.write((char*)&size,sizeof(Int_t));
    
    if (nTrkTracks < 1) continue;
    
    out.write((char*)&nTrkTracks,sizeof(Int_t));
    
    for (Int_t iTrack = 0; iTrack < nTracks; iTrack++) {
      
      AliESDMuonTrack *track = esd->GetMuonTrack(iTrack);
      if (!track->ContainTrackerData()) continue;

      // refit the track if requested
      AliMUONESDInterface::ESDToMUON(*track, muonTrack, refit);
      AliMUONTrackParam *param = static_cast<AliMUONTrackParam*>(muonTrack.GetTrackParamAtCluster()->First());
      /*
      // store track parameters at vertex
      sTrackParam.x = track->GetNonBendingCoor();
      sTrackParam.y = track->GetBendingCoor();
      sTrackParam.z = track->GetZ();
      sTrackParam.px = track->Px();
      sTrackParam.py = track->Py();
      sTrackParam.pz = track->Pz();
      sTrackParam.sign = track->Charge();
      out.write((char*)&sTrackParam,sizeof(TrackParamStruct));
      */
      // store track parameters at first cluster
      sTrackParam.x = param->GetNonBendingCoor();
      sTrackParam.y = param->GetBendingCoor();
      sTrackParam.z = param->GetZ();
      sTrackParam.px = param->Px();
      sTrackParam.py = param->Py();
      sTrackParam.pz = param->Pz();
      sTrackParam.sign = param->GetCharge();
      out.write((char*)&sTrackParam,sizeof(TrackParamStruct));

      Int_t nClusters = track->GetNClusters();
      out.write((char*)&nClusters,sizeof(Int_t));
      
      for (Int_t iCl = 0; iCl < nClusters; iCl++) {
        
        AliESDMuonCluster *cluster = esd->FindMuonCluster(track->GetClusterId(iCl));
        
        // store cluster information
        sCluster.x = cluster->GetX();
        sCluster.y = cluster->GetY();
        sCluster.z = cluster->GetZ();
        sCluster.ex = cluster->GetErrX();
        sCluster.ey = cluster->GetErrY();
        sCluster.uid = cluster->GetUniqueID();
        out.write((char*)&sCluster,sizeof(ClusterStruct));
        
      }
      
    }
    
  }
  
  out.close();
  
}

