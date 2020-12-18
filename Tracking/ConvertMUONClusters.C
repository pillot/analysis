#include <fstream>
#include <iostream>
#include <chrono>
#include <list>
#include <vector>

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
#include "AliMpDDLStore.h"
#include "AliMpDetElement.h"
#include "AliMpPad.h"
#include "AliMpCathodType.h"
#include "AliMpPlaneType.h"
#include "AliMpVSegmentation.h"
#include "AliMpSegmentation.h"

#include "AliMUONCDB.h"
#include "AliMUONRecoParam.h"
#include "AliMUONESDInterface.h"
#include "AliMUONVDigit.h"
#include "AliMUONVCluster.h"
#include "AliMUONVClusterStore.h"
#include "AliMUONTrack.h"
#include "AliMUONTrackParam.h"
#include "AliMUONVTrackStore.h"
#include "AliMUONTracker.h"
#include "AliMUONVTrackReconstructor.h"
#include "AliMUONGeometryTransformer.h"

#include "Field/MagneticField.h"

#include "DataFormatsMCH/TrackMCH.h"
#include "MCHBase/ClusterBlock.h"

using namespace std;
using namespace o2::mch;

AliMUONGeometryTransformer* transformer(nullptr);

bool SetMagField();
AliMUONGeometryTransformer* LoadGeometry(int runNumber);
AliMUONVTrackReconstructor* CreateTrackReconstructor();
void ChangeMonoCathodClusterRes(AliMUONVClusterStore& clusterStore);
void MUONToO2(AliMUONVCluster& cluster, ClusterStruct& o2Cluster);
void WriteClusters(AliMUONVClusterStore& clusterStore, ofstream& outFile);
void MUONToO2(AliMUONTrack& track, TrackMCH& o2Track, std::vector<ClusterStruct>& o2Clusters);
void WriteTracks(AliMUONVTrackStore& trackStore, ofstream& outFile);

//------------------------------------------------------------------
void ConvertMUONClusters(int runNumber, TString inFileName, TString outFileName = "clusters.in",
                         bool findTracks = kFALSE, int event = -1)
{
  /// convert MUON clusters from ESD or RecPoints into O2 structures
  /// saved in a binary file with the following format:
  ///
  /// number of clusters in event 1
  /// number of associated digits (= 0)
  /// ClusterStruct of 1st cluster
  /// ClusterStruct of 2nd cluster
  /// ...
  /// ClusterStruct of nth cluster
  /// number of clusters in event 2
  /// ...
  ///
  /// if findTracks = kTRUE: reconstruct tracks from clusters with AliRoot track finder
  /// and save the result in a binary file with the following format:
  ///
  /// number of tracks at vertex in event 1 (= 0)
  /// number of MCH tracks in event 1
  /// number of associated clusters in event 1
  /// list of TrackMCH
  /// list of ClusterStruct
  /// number of tracks at vertex in event 2 (= 0)
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

  // load geometry and mapping, to find mono-cathod clusters
  transformer = LoadGeometry(runNumber);
  if (!transformer) {
    Error("ConvertMUONCluster", "cannot get the geometry transformer");
    return;
  }

  AliMUONVCluster* cluster = clusterStore->CreateCluster(0, 0, 0);
  AliMUONVTrackStore* trackStore = findTracks ? AliMUONESDInterface::NewTrackStore() : nullptr;
  AliMUONVTrackReconstructor* trackReconstructor = findTracks ? CreateTrackReconstructor() : nullptr;
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

    if (esd) {
      // add the clusters to the cluster store manually in case they are taken from ESD
      int nClusters = esd->GetNumberOfMuonClusters();
      for (int iCluster = 0; iCluster < nClusters; ++iCluster) {
        AliESDMuonCluster* esdCluster = esd->GetMuonCluster(iCluster);
        AliMUONESDInterface::ESDToMUON(*esdCluster, *cluster);
        // reset the resolution of the mono-cathod clusters and others in double precision before running the reconstruction
        cluster->SetErrXY(0.2, 0.2);
        clusterStore->Add(*cluster);
      }
    } else {
      TIter nextCl(clusterStore->CreateIterator());
      AliMUONVCluster* cl(nullptr);
      int iCl(0);
      std::list<uint32_t> clustersToRemove{};
      while ((cl = static_cast<AliMUONVCluster*>(nextCl()))) {
        // reset the cluster resolution in double precision before running the reconstruction
        cl->SetErrXY(0.2, 0.2);
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
    }

    // change the resolution of all mono-cathod clusters
    ChangeMonoCathodClusterRes(*clusterStore);

    // write the clusters in the binary file
    WriteClusters(*clusterStore, outClusterFile);

    if (findTracks) {

      // refit the tracks to retreive the exact same parameters as in O2
      AliMUONTrack* track(nullptr);
      TIter next(trackStore->CreateIterator());
      while ((track = static_cast<AliMUONTrack*>(next()))) {
        for (int iCl = 0; iCl < track->GetNClusters(); ++iCl) {
          AliMUONVCluster* cl = static_cast<AliMUONTrackParam*>(track->GetTrackParamAtCluster()->UncheckedAt(iCl))->GetClusterPtr();
          // reset the cluster resolution in float, as in O2
          cl->SetErrXY((float)cl->GetErrX(), (float)cl->GetErrY());
        }
        trackReconstructor->RefitTrack(*track, kFALSE);
      }

      // write the tracks in the second binary file
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
AliMUONGeometryTransformer* LoadGeometry(int runNumber)
{
  /// load the geometry and mapping from the OCDB

  // set OCDB location
  AliCDBManager* man = AliCDBManager::Instance();
  if (gSystem->AccessPathName("OCDB.root", kFileExists) == 0) {
    man->SetDefaultStorage("local:///dev/null");
    man->SetSnapshotMode("OCDB.root");
  } else {
    man->SetDefaultStorage("local://./OCDB");
  }
  man->SetRun(runNumber);

  // load the geometry
  AliGeomManager::LoadGeometry();
  if (!AliGeomManager::GetGeometry() || !AliGeomManager::ApplyAlignObjsFromCDB("MUON")) {
    return nullptr;
  }

  // load the mapping
  if (!AliMpCDB::LoadDDLStore()) {
    return nullptr;
  }

  // create the geometry transformer
  AliMUONGeometryTransformer* transformer = new AliMUONGeometryTransformer();
  transformer->LoadGeometryData();

  return transformer;
}

//------------------------------------------------------------------
AliMUONVTrackReconstructor* CreateTrackReconstructor()
{
  /// load the magnetic field and recoParam from the OCDB and prepare track finding

  if (!SetMagField()) {
//  if (!AliMUONCDB::LoadField()) {
    return nullptr;
  }

  AliMUONRecoParam* recoParam = AliMUONCDB::LoadRecoParam();
  if (!recoParam) {
    return nullptr;
  }

  recoParam->DiscardMonoCathodClusters(true, 10., 10.);
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

  return AliMUONTracker::CreateTrackReconstructor(recoParam, nullptr, transformer);
}

//------------------------------------------------------------------
void ChangeMonoCathodClusterRes(AliMUONVClusterStore& clusterStore)
{
  /// assign a different resolution to the mono-cathod clusters in the direction of the missing plane

  // loop over clusters in stations 3, 4 and 5
  TIter next(clusterStore.CreateChamberIterator(4, 9));
  AliMUONVCluster* cluster(nullptr);
  while ((cluster = static_cast<AliMUONVCluster*>(next()))) {

    // get the cathod corresponding to the bending/non-bending plane
    int deId = cluster->GetDetElemId();
    AliMpDetElement* de = AliMpDDLStore::Instance()->GetDetElement(deId, false);
    AliMp::CathodType cath1 = de->GetCathodType(AliMp::kBendingPlane);
    AliMp::CathodType cath2 = de->GetCathodType(AliMp::kNonBendingPlane);

    // get the corresponding segmentation
    const AliMpVSegmentation* seg1 = AliMpSegmentation::Instance()->GetMpSegmentation(deId, cath1);
    const AliMpVSegmentation* seg2 = AliMpSegmentation::Instance()->GetMpSegmentation(deId, cath2);

    // get local coordinate of the cluster
    double lX(0.), lY(0.), lZ(0.);
    transformer->Global2Local(deId, cluster->GetX(), cluster->GetY(), cluster->GetZ(), lX, lY, lZ);

    // find pads below the cluster
    AliMpPad pad1 = seg1->PadByPosition(lX, lY, false);
    AliMpPad pad2 = seg2->PadByPosition(lX, lY, false);

    // build their ID if pads are valid
    uint32_t padId1 = (pad1.IsValid()) ? AliMUONVDigit::BuildUniqueID(deId, pad1.GetManuId(), pad1.GetManuChannel(), cath1) : 0;
    uint32_t padId2 = (pad2.IsValid()) ? AliMUONVDigit::BuildUniqueID(deId, pad2.GetManuId(), pad2.GetManuChannel(), cath2) : 0;

    // check if the cluster contains these pads
    bool hasNonBending(false);
    bool hasBending(false);
    for (int iDigit = 0; iDigit < cluster->GetNDigits(); ++iDigit) {
      if (cluster->GetDigitId(iDigit) == padId1) {
        hasBending = true;
        if (hasNonBending) break;
      } else if (cluster->GetDigitId(iDigit) == padId2) {
        hasNonBending = true;
        if (hasBending) break;
      }
    }

    // modify the cluster resolution if needed
    if (!hasNonBending) cluster->SetErrXY(10., cluster->GetErrY());
    if (!hasBending) cluster->SetErrXY(cluster->GetErrX(), 10.);
  }
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

  // get the number of clusters effectively stored (without duplicates)
  int nClusters = clusterStore.GetSize();
  outFile.write((char*)&nClusters, sizeof(int));

  // write the number of associated digits
  int nDigits(0);
  outFile.write((char*)&nDigits, sizeof(int));

  // write the clusters
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
void MUONToO2(AliMUONTrack& track, TrackMCH& o2Track, std::vector<ClusterStruct>& o2Clusters)
{
  /// copy the track info from the MUON track into the O2 TrackMCH

  AliMUONTrackParam* param = static_cast<AliMUONTrackParam*>(track.GetTrackParamAtCluster()->First());
  o2Track.setZ(param->GetZ());
  o2Track.setParameters(param->GetParameters());
  o2Track.setCovariances(param->GetCovariances());
  o2Track.setChi2(track.GetGlobalChi2());
  o2Track.setClusterRef(o2Clusters.size(), track.GetNClusters());

  for (int iCl = 0; iCl < track.GetNClusters(); ++iCl) {
    o2Clusters.emplace_back();
    MUONToO2(*static_cast<AliMUONTrackParam*>(track.GetTrackParamAtCluster()->UncheckedAt(iCl))->GetClusterPtr(), o2Clusters.back());
  }
}

//------------------------------------------------------------------
void WriteTracks(AliMUONVTrackStore& trackStore, ofstream& outFile)
{
  /// write the tracks in O2 format in the output file

  // count the number of clusters attached to the tracks
  int nClusters(0);
  AliMUONTrack* track(nullptr);
  TIter next(trackStore.CreateIterator());
  while ((track = static_cast<AliMUONTrack*>(next()))) {
    nClusters += track->GetNClusters();
  }

  // write the number of tracks at vertex (= 0), MCH tracks and attached clusters
  int nTracksAtVtx = 0;
  outFile.write(reinterpret_cast<char*>(&nTracksAtVtx), sizeof(int));
  int nTracks = trackStore.GetSize();
  outFile.write(reinterpret_cast<char*>(&nTracks), sizeof(int));
  outFile.write(reinterpret_cast<char*>(&nClusters), sizeof(int));

  if (nTracks == 0) {
    return;
  }

  // write the MCH tracks and store the attached clusters
  std::vector<ClusterStruct> o2Clusters{};
  o2Clusters.reserve(nClusters);
  TrackMCH o2Track{};
  next.Reset();
  while ((track = static_cast<AliMUONTrack*>(next()))) {
    MUONToO2(*track, o2Track, o2Clusters);
    outFile.write((char*)&o2Track, sizeof(TrackMCH));
  }

  // write the attached clusters
  outFile.write(reinterpret_cast<char*>(o2Clusters.data()), o2Clusters.size() * sizeof(ClusterStruct));
}
