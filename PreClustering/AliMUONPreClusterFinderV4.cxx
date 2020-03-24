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

// $Id$

#include "AliMUONPreClusterFinderV4.h"

#include "AliMUONCluster.h"
#include "AliMUONPad.h"

#include "AliMpArea.h"

#include "AliLog.h"

#include <Riostream.h>
#include <TObjArray.h>
#include <TVector2.h>

//-----------------------------------------------------------------------------
/// \class AliMUONPreClusterFinderV4
///
/// Implementation of AliMUONVClusterFinder
///
/// This class simply find adjacent pads to form preclusters
///
/// \author Philippe Pillot
//-----------------------------------------------------------------------------

ClassImp(AliMUONPreClusterFinderV4)

//_____________________________________________________________________________
AliMUONPreClusterFinderV4::AliMUONPreClusterFinderV4()
: AliMUONVClusterFinder(),
  fClusters("AliMUONCluster"),
  fPads(0x0),
  fDetElemId(0),
  fArea()
{
  /// ctor
}

//_____________________________________________________________________________
AliMUONPreClusterFinderV4::~AliMUONPreClusterFinderV4()
{
  /// dtor : note we're owner of the preclusters, but not of the pads
}

//_____________________________________________________________________________
Bool_t AliMUONPreClusterFinderV4::Prepare(Int_t detElemId, TObjArray* pads[2], const AliMpArea& area)
{
  /// Prepare for preclustering, by giving access to digit lists

  fClusters.Delete();
  
  fPads = pads;
  fDetElemId = detElemId;
  fArea = area;
  
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliMUONPreClusterFinderV4::UsePad(const AliMUONPad& pad)
{
  /// Add a pad to the list of pads to be considered

  if ( pad.DetElemId() != fDetElemId ) {
    AliError(Form("Cannot add pad from DE %d to this cluster finder which is "
                  "currently dealing with DE %d",pad.DetElemId(),fDetElemId));
    return kFALSE;
  }
  
  fPads[pad.Cathode()]->Add(new AliMUONPad(pad)); 

  return kTRUE;
}

//_____________________________________________________________________________
AliMUONCluster* AliMUONPreClusterFinderV4::NextCluster()
{
  /// Builds the next precluster, and returns it
  
  // Get a pad to start a new precluster, if any
  AliMUONPad* pad = GetNextPad(0);
  if (!pad) pad = GetNextPad(1);
  if (!pad) return 0x0;

  // Build the precluster recursively
  AliMUONCluster* cluster = NewCluster();
  AddPad(*cluster, pad);
/*
  // Remove this precluster if made of only 1 pad and proceed to the next one
  if (cluster->Multiplicity() <= 1) {
    fClusters.RemoveAt(fClusters.GetLast());
    return NextCluster();
  }
*/
  return cluster;
}

//_____________________________________________________________________________
AliMUONPad* AliMUONPreClusterFinderV4::GetNextPad(Int_t cathode) const
{
/// Return the next unused pad of given cathode, which is within fArea

  TIter next(fPads[cathode]);

  if (!fArea.IsValid()) {
    return static_cast<AliMUONPad*>(next());
  } else {
    AliMUONPad* pad;
    while ((pad = static_cast<AliMUONPad*>(next()))) {
      AliMpArea padArea(pad->X(), pad->Y(), pad->DX(), pad->DY());
      if (fArea.Overlap(padArea)) return pad;
    }
  }

  return 0x0;
}

//_____________________________________________________________________________
AliMUONCluster* AliMUONPreClusterFinderV4::NewCluster()
{
  /// Create a new (empty) cluster

  Int_t id = fClusters.GetLast()+1;
  AliMUONCluster* cluster = new (fClusters[id]) AliMUONCluster;
  cluster->SetUniqueID(id);

  return cluster;
}

//_____________________________________________________________________________
void AliMUONPreClusterFinderV4::AddPad(AliMUONCluster& cluster, AliMUONPad* pad)
{
  /// Add a pad and its neighbours to a precluster (recursive method)

  static Double_t precision = 1E-4; // cm
  static TVector2 precisionAdjustment(precision, precision);    

  // add the current pad to the cluster and remove it from the list
  AliMUONPad* addedPad = cluster.AddPad(*pad);
  Int_t cathode = pad->Cathode();
  TObject* o = fPads[cathode]->Remove(pad);
  delete o;
  
  // add neighbouring pads on the same cathod
  TIter next1(fPads[cathode]);
  AliMUONPad* testPad(0x0);
  while ((testPad = static_cast<AliMUONPad*>(next1()))) {
    if (AliMUONPad::AreNeighbours(*testPad, *addedPad)) {
      AddPad(cluster, testPad);
    }
  }

  // add neighbouring pads on the other cathod
  TIter next2(fPads[1-cathode]);
  while ((testPad = static_cast<AliMUONPad*>(next2()))) {
    if (AliMUONPad::AreOverlapping(*testPad, *addedPad, precisionAdjustment)) {
      AddPad(cluster, testPad);
    }
  }
}
