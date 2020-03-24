#ifndef ALIMUONPRECLUSTERFINDERV4_H
#define ALIMUONPRECLUSTERFINDERV4_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup rec
/// \class AliMUONPreClusterFinderV4
/// \brief A basic pre-cluster finder
/// 
// Author Philippe Pillot, Subatech

#ifndef AliMUONVCLUSTERFINDER_H
#  include "AliMUONVClusterFinder.h"
#endif
#ifndef ALI_MP_AREA_H
#  include "AliMpArea.h"
#endif
#ifndef ROOT_TClonesArray
#  include <TClonesArray.h>
#endif

class TStopwatch;
class AliMUONPad;
class TObjArray;

class AliMUONPreClusterFinderV4 : public AliMUONVClusterFinder
{
public:
  AliMUONPreClusterFinderV4();
  virtual ~AliMUONPreClusterFinderV4();

  using AliMUONVClusterFinder::Prepare;

  virtual Bool_t Prepare(Int_t detElemId, TObjArray* pads[2], const AliMpArea& area);

  virtual Bool_t UsePad(const AliMUONPad& pad);

  virtual AliMUONCluster* NextCluster();

private:
  /// Not implemented
  AliMUONPreClusterFinderV4(const AliMUONPreClusterFinderV4& rhs);
  /// Not implemented
  AliMUONPreClusterFinderV4& operator=(const AliMUONPreClusterFinderV4& rhs);

  AliMUONPad* GetNextPad(Int_t cathode) const;

  AliMUONCluster* NewCluster();

  void AddPad(AliMUONCluster& cluster, AliMUONPad* pad);

private:
  TClonesArray fClusters; //!<! the clusters we've found (owner)
  TObjArray** fPads; //!<! the pads corresponding to the digits (not owner)
  Int_t fDetElemId; //!<! which DE we're considering
  AliMpArea fArea; //!<! area into which to consider pads to *start* a cluster

  ClassDef(AliMUONPreClusterFinderV4,1) // A basic pre-cluster finder
};

#endif
