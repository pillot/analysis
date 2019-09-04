#include <fstream>
#include <iostream>
#include <chrono>
#include <list>

#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TIterator.h>
#include <TGeoGlobalMagField.h>

#include "AliCDBManager.h"
#include "AliGRPManager.h"
#include "AliGRPObject.h"
#include "AliGeomManager.h"

#include "AliESDEvent.h"
#include "AliESDMuonCluster.h"

#include "AliMpCDB.h"

#include "AliMUONCDB.h"
#include "AliMUONRecoParam.h"
#include "AliMUONESDInterface.h"
#include "AliMUONVCluster.h"
#include "AliMUONVClusterStore.h"
#include "AliMUONTrack.h"
#include "AliMUONTrackParam.h"
#include "AliMUONVTrackStore.h"
#include "AliMUONTracker.h"
#include "AliMUONVTrackReconstructor.h"
#include "AliMUONGeometryTransformer.h"

#include "Field/MagneticField.h"

#include "MCHBase/ClusterBlock.h"
#include "MCHBase/TrackBlock.h"

using namespace std;
using namespace o2::mch;

bool SetMagField();
AliMUONVTrackReconstructor* CreateTrackReconstructor(int runNumber);
void MUONToO2(AliMUONVCluster& cluster, ClusterStruct& o2Cluster);
void WriteClusters(AliMUONVClusterStore& clusterStore, ofstream& outFile);
void MUONToO2(AliMUONTrackParam& trackParam, TrackParamStruct& o2TrackParam);
void WriteTracks(AliMUONVTrackStore& trackStore, ofstream& outFile);

//------------------------------------------------------------------
void ConvertMUONClusters(int runNumber, TString inFileName, TString outFileName = "clusters.in",
                         bool findTracks = kFALSE, int event = -1)
{
  /// convert MUON clusters from ESD or RecPoints into O2 structures
  /// saved in a binary file with the following format:
  ///
  /// event number
  /// number of clusters in event 1
  /// ClusterStruct of 1st cluster
  /// ClusterStruct of 2nd cluster
  /// ...
  /// ClusterStruct of nth cluster
  /// event number
  /// number of clusters in event 2
  /// ...
  ///
  /// if findTracks = kTRUE: reconstruct tracks from clusters with AliRoot track finder
  /// and save the result in a binary file with the following format:
  ///
  /// event number
  /// event size (number of bytes)
  /// number of tracks in event 1
  /// TrackParamStruct of 1st track
  /// number of clusters in track 1
  /// ClusterStruct of 1st cluster
  /// ClusterStruct of 2nd cluster
  /// ...
  /// ClusterStruct of nth cluster
  /// TrackParamStruct of 2nd track
  /// number of clusters in track 2
  /// ...
  /// event number
  /// event size (number of bytes)
  /// number of tracks in event 2
  /// ...

  // open the input file
  TFile* inFile = TFile::Open(inFileName);
  if (!inFile || !inFile->IsOpen()) {
    Error("ConvertMUONCluster", "opening file %s failed", inFileName.Data());
    return;
  }

  // get the input data
  AliESDEvent* esd = nullptr;
  AliMUONVClusterStore* clusterStore = nullptr;
  TTree *tree = static_cast<TTree*>(inFile->Get("esdTree"));
  if (tree) {
    esd = new AliESDEvent();
    esd->ReadFromTree(tree);
    clusterStore = AliMUONESDInterface::NewClusterStore();
  } else {
    tree = static_cast<TTree*>(inFile->Get("TreeR"));
    if (!tree) {
      Error("ConvertMUONCluster", "neither esdTree nor TreeR found in the input file");
      return;
    }
    clusterStore = AliMUONVClusterStore::Create(*tree);
    clusterStore->Connect(*tree);
  }

  // open output files
  ofstream outClusterFile(outFileName.Data(), ios::out | ios::binary);
  if (!outClusterFile.is_open()) {
    return;
  }
  ofstream* outTrackFile(nullptr);
  if (findTracks) {
    outTrackFile = new ofstream("AliRootTracks.out", ios::out | ios::binary);
    if (!outTrackFile->is_open()) {
      return;
    }
  }

  AliMUONVCluster* cluster = clusterStore->CreateCluster(0, 0, 0);
  AliMUONVTrackStore* trackStore = findTracks ? AliMUONESDInterface::NewTrackStore() : nullptr;
  AliMUONVTrackReconstructor* trackReconstructor = findTracks ? CreateTrackReconstructor(runNumber) : nullptr;
  if (findTracks && !trackReconstructor) {
    return;
  }
  std::chrono::duration<double> trackingTime{};

  int nevents = (event >= 0) ? TMath::Min(event + 1, (int)tree->GetEntries()) : (int)tree->GetEntries();
  int iEvent = (event >= 0) ? event : 0;
  for (; iEvent < nevents; iEvent++) {

    // read the event
    clusterStore->Clear();
    if (tree->GetEntry(iEvent) <= 0) {
      Error("ConvertMUONCluster", "no entry object found for event %d", iEvent);
      return;
    }
    outClusterFile.write((char*)&iEvent, sizeof(int));

    if (esd) {
      // add the clusters to the cluster store manually in case they are taken from ESD
      int nClusters = esd->GetNumberOfMuonClusters();
      for (int iCluster = 0; iCluster < nClusters; ++iCluster) {
        AliESDMuonCluster* esdCluster = esd->GetMuonCluster(iCluster);
        AliMUONESDInterface::ESDToMUON(*esdCluster, *cluster);
        // reset the resolution of the mono-cathod clusters before running the reconstruction
        // make sure to use float precision so that the values are identical both in O2 and AliRoot format
        if (trackReconstructor) {
          cluster->SetErrXY(trackReconstructor->GetRecoParam()->GetDefaultNonBendingReso(cluster->GetChamberId()),
                            trackReconstructor->GetRecoParam()->GetDefaultBendingReso(cluster->GetChamberId()));
        } else {
          cluster->SetErrXY(0.2f, 0.2f);
        }
        clusterStore->Add(*cluster);
      }
    } else {
      // cluster resolution data members are ex and ey in O2 while they are ex^2 and ey^2 in AliRoot
      // by resetting it to a float the resolution used both in O2 and AliRoot is the same and so are the tracks
      TIter nextCl(clusterStore->CreateIterator());
      AliMUONVCluster* cl(nullptr);
      int iCl(0);
      std::list<uint32_t> clustersToRemove{};
      while ((cl = static_cast<AliMUONVCluster*>(nextCl()))) {
        cl->SetErrXY((float)cl->GetErrX(), (float)cl->GetErrY());
        // find duplicate clusters (same position but different Id)
        TIter nextCl2(clusterStore->CreateIterator());
        AliMUONVCluster* cl2(nullptr);
        int iCl2(-1);
        while ((cl2 = static_cast<AliMUONVCluster*>(nextCl2()))) {
          if (++iCl2 > iCl && cl2->GetX() == cl->GetX() && cl2->GetY() == cl->GetY() && cl2->GetZ() == cl->GetZ()) {
            clustersToRemove.push_back(cl->GetUniqueID());
            break;
          }
        }
        ++iCl;
      }
      // remove duplicate clusters
      for (const auto& id : clustersToRemove) {
        clusterStore->Remove(*(clusterStore->FindObject(id)));
      }
    }

    // get the number of clusters effectively stored (without duplicates)
    int nClusters = clusterStore->GetSize();
    outClusterFile.write((char*)&nClusters, sizeof(int));

    // write the clusters in the binary file
    WriteClusters(*clusterStore, outClusterFile);

    if (findTracks) {

      // reconstruct tracks from clusters
      trackStore->Clear();
      auto tStart = std::chrono::high_resolution_clock::now();
      trackReconstructor->EventReconstruct(*clusterStore, *trackStore);
      auto tEnd = std::chrono::high_resolution_clock::now();
      trackingTime += tEnd - tStart;
      if (clusterStore->GetSize() > 0 && trackStore->GetSize() == 0) {
        cout << "no track has been retreived in event " << iEvent << endl;
      }
/*
      // refit the tracks to retreive the exact same parameters as the refitted ESD tracks
      AliMUONTrack* track(nullptr);
      TIter next(trackStore->CreateIterator());
      while ((track = static_cast<AliMUONTrack*>(next()))) {
        trackReconstructor->RefitTrack(*track, kFALSE);
      }
*/
      // write the tracks in the second binary file
      outTrackFile->write((char*)&iEvent, sizeof(int));
      WriteTracks(*trackStore, *outTrackFile);
    }
  }

  cout << "tracking duration = " << trackingTime.count() << " s" << endl;

  // close files and cleanup
  outClusterFile.close();
  if (findTracks) {
    outTrackFile->close();
  }
  delete cluster;
  delete clusterStore;
  delete trackStore;
  delete trackReconstructor;
}

//------------------------------------------------------------------
bool SetMagField()
{
  /// Set the magnetic field using O2 maps and GRP info
  
  AliGRPManager grpMan;
  if (!grpMan.ReadGRPEntry()) {
    Error("SetMagField", "failed to load GRP Data from OCDB");
    return false;
  }
  
  const AliGRPObject *grpData = grpMan.GetGRPData();
  if (!grpData) {
    Error("SetMagField", "GRP Data is not loaded");
    return false;
  }
  
  float l3Current = grpData->GetL3Current((AliGRPObject::Stats)0);
  if (l3Current == AliGRPObject::GetInvalidFloat()) {
    Error("SetMagField", "GRP/GRP/Data entry:  missing value for the L3 current !");
    return false;
  }

  char l3Polarity = grpData->GetL3Polarity();
  if (l3Polarity == AliGRPObject::GetInvalidChar()) {
    Error("SetMagField", "GRP/GRP/Data entry:  missing value for the L3 polarity !");
    return false;
  }
  
  float diCurrent = grpData->GetDipoleCurrent((AliGRPObject::Stats)0);
  if (diCurrent == AliGRPObject::GetInvalidFloat()) {
    Error("SetMagField", "GRP/GRP/Data entry:  missing value for the dipole current !");
    return false;
  }

  char diPolarity = grpData->GetDipolePolarity();
  if (diPolarity == AliGRPObject::GetInvalidChar()) {
    Error("SetMagField", "GRP/GRP/Data entry:  missing value for the dipole polarity !");
    return false;
  }
  
  float beamEnergy = grpData->GetBeamEnergy();
  if (beamEnergy==AliGRPObject::GetInvalidFloat()) {
    Error("SetMagField", "GRP/GRP/Data entry:  missing value for the beam energy !");
    return false;
  }

  TString beamType = grpData->GetBeamType();
  if (beamType==AliGRPObject::GetInvalidString()) {
    Error("SetMagField", "GRP/GRP/Data entry:  missing value for the beam type !");
    return false;
  }
  
  Info("SetMagField", "l3Current = %f, diCurrent = %f", TMath::Abs(l3Current) * (l3Polarity ? -1:1),
                                                        TMath::Abs(diCurrent) * (diPolarity ? -1:1));
  
  auto field = o2::field::MagneticField::createFieldMap(TMath::Abs(l3Current) * (l3Polarity ? -1:1),
                                                        TMath::Abs(diCurrent) * (diPolarity ? -1:1),
                                                        o2::field::MagneticField::kConvLHC, false, beamEnergy, beamType.Data(),
                                                        "$(O2_ROOT)/share/Common/maps/mfchebKGI_sym.root");
  TGeoGlobalMagField::Instance()->SetField(field);
  TGeoGlobalMagField::Instance()->Lock();

  return true;
}

//------------------------------------------------------------------
AliMUONVTrackReconstructor* CreateTrackReconstructor(int runNumber)
{
  /// access OCDB and prepare track finding

  AliCDBManager* man = AliCDBManager::Instance();
  if (gSystem->AccessPathName("OCDB.root", kFileExists)==0) {
    man->SetDefaultStorage("local:///dev/null");
    man->SetSnapshotMode("OCDB.root");
  } else {
    man->SetDefaultStorage("local://./OCDB");
  }
  man->SetRun(runNumber);

  if (!SetMagField()) {
//  if (!AliMUONCDB::LoadField()) {
    return nullptr;
  }

  AliGeomManager::LoadGeometry();
  if (!AliGeomManager::GetGeometry() || !AliGeomManager::ApplyAlignObjsFromCDB("MUON")) {
    return nullptr;
  }

  AliMUONRecoParam* recoParam = AliMUONCDB::LoadRecoParam();
  if (!recoParam) {
    return nullptr;
  }

  bool discardMonoCathodClusters = false;
  recoParam->DiscardMonoCathodClusters(discardMonoCathodClusters);
  //recoParam->MakeMoreTrackCandidates(true);
  //recoParam->RequestStation(0, false);
  //recoParam->RequestStation(1, false);
  //recoParam->RequestStation(2, false);
  //recoParam->RequestStation(3, false);
  //recoParam->RequestStation(4, false);
  //recoParam->RecoverTracks(false);
  //recoParam->ComplementTracks(false);
  //recoParam->ImproveTracks(false);
  recoParam->SetMaxTrackCandidates(1000000);

  AliMUONGeometryTransformer* transformer = nullptr;
  if (discardMonoCathodClusters) {
    if (!AliMpCDB::LoadDDLStore()) {
      return nullptr;
    }
    transformer = new AliMUONGeometryTransformer;
    transformer->LoadGeometryData();
  }

  return AliMUONTracker::CreateTrackReconstructor(recoParam, nullptr, transformer);
}

//------------------------------------------------------------------
void MUONToO2(AliMUONVCluster& cluster, ClusterStruct& o2Cluster)
{
  /// copy the cluster information from the MUON cluster into the O2 cluster
  o2Cluster.x = cluster.GetX();
  o2Cluster.y = cluster.GetY();
  o2Cluster.z = cluster.GetZ();
  o2Cluster.ex = cluster.GetErrX();
  o2Cluster.ey = cluster.GetErrY();
  o2Cluster.uid = cluster.GetUniqueID();
}

//------------------------------------------------------------------
void WriteClusters(AliMUONVClusterStore& clusterStore, ofstream& outFile)
{
  /// loop over clusters per chamber (that way they are ordered the same way during the AliRoot track finding)
  /// convert them to O2 format and write them in the output file
  ClusterStruct o2Cluster{};
  for (int iChamber = 0; iChamber < 10; ++iChamber) {
    TIter nextInCh(clusterStore.CreateChamberIterator(iChamber, iChamber));
    AliMUONVCluster* cluster(nullptr);
    while ((cluster = static_cast<AliMUONVCluster*>(nextInCh()))) {
      MUONToO2(*cluster, o2Cluster);
      outFile.write((char*)&o2Cluster, sizeof(ClusterStruct));
    }
  }
}

//------------------------------------------------------------------
void MUONToO2(AliMUONTrackParam& trackParam, TrackParamStruct& o2TrackParam)
{
  /// copy the track parameters from the MUON trackParam into the O2 trackParam
  o2TrackParam.x = trackParam.GetNonBendingCoor();
  o2TrackParam.y = trackParam.GetBendingCoor();
  o2TrackParam.z = trackParam.GetZ();
  o2TrackParam.px = trackParam.Px();
  o2TrackParam.py = trackParam.Py();
  o2TrackParam.pz = trackParam.Pz();
  o2TrackParam.sign = trackParam.GetCharge();
}

//------------------------------------------------------------------
void WriteTracks(AliMUONVTrackStore& trackStore, ofstream& outFile)
{
  /// write the tracks in the output file

  int size(0); // total number of bytes requested to store the event, excluding event number and total size
  int nTracks(trackStore.GetSize());
  if (nTracks == 0) {
    outFile.write((char*)&size, sizeof(int));
    return;
  } else {
    size = sizeof(int);
  }
  AliMUONTrack* track(nullptr);
  TIter next(trackStore.CreateIterator());
  while ((track = static_cast<AliMUONTrack*>(next()))) {
    size += sizeof(TrackParamStruct) + sizeof(double) + sizeof(int) + track->GetNClusters() * sizeof(ClusterStruct);
  }
  outFile.write((char*)&size, sizeof(int));

  outFile.write((char*)&nTracks, sizeof(int));

  next.Reset();
  TrackParamStruct o2TrackParam{};
  ClusterStruct o2Cluster{};
  while ((track = static_cast<AliMUONTrack*>(next()))) {
    AliMUONTrackParam* param(static_cast<AliMUONTrackParam*>(track->GetTrackParamAtCluster()->First()));
    MUONToO2(*param, o2TrackParam);
    outFile.write((char*)&o2TrackParam, sizeof(TrackParamStruct));

    double chi2 = track->GetGlobalChi2();
    outFile.write((char*)&chi2, sizeof(double));

    int nClusters(track->GetNClusters());
    outFile.write((char*)&nClusters, sizeof(int));

    for (int iCl = 0; iCl < nClusters; ++iCl) {
      AliMUONVCluster* cluster(static_cast<AliMUONTrackParam*>(track->GetTrackParamAtCluster()->UncheckedAt(iCl))->GetClusterPtr());
      MUONToO2(*cluster, o2Cluster);
      outFile.write((char*)&o2Cluster, sizeof(ClusterStruct));
    }
  }
}
