//
//  PreClusterizerV2.h
//  aliroot_dev
//
//  Created by philippe pillot on 13/03/2014.
//  Copyright (c) 2014 Philippe Pillot. All rights reserved.
//

#ifndef PRECLUSTERIZERV2_H
#define PRECLUSTERIZERV2_H

#include <vector>
#include <TExMap.h>

struct mpPad {
  UChar_t nNeighbours; // number of neighbours
  UShort_t neighbours[10]; // indices of neighbours in array stored in mpDE
  Float_t area[2][2]; // 2D area
  AliMUONVDigit *digit; // pointer to the associated digit
  Bool_t useMe; // kFALSE if no digit attached or already visited
};

struct mpDE {
  Int_t id; // unique ID
  UChar_t iPlanes[2]; // plane type corresponding to both cathods
  UShort_t nPads[2]; // number of pads on each plane
  mpPad *pads; // array of pads on both planes
  TExMap padIndices[2]; // indices of pads from their ID
  UShort_t nFiredPads[2]; // number of fired pads on each plane
  std::vector<UShort_t> firedPads[2]; // indices of fired pads on each plane
  UShort_t nOrderedPads[2]; // current number of fired pads in the following arrays
  std::vector<UShort_t> orderedPads[2]; // indices of fired pads ordered after preclustering and merging
};

struct preCluster {
  UShort_t firstPad; // index of first associated pad in the orderedPads array
  UShort_t lastPad; // index of last associated pad in the orderedPads array
  Float_t area[2][2]; // 2D area containing the precluster
  Bool_t useMe; // kFALSE if precluster already merged to another one
  Bool_t storeMe; // kTRUE if precluster to be saved (merging result)
};

class AliMUONVStore;
class AliMUONVDigitStore;

class PreClusterizerV2 : public TObject
{
public:
  
  PreClusterizerV2(const char* ocdbpath="local://$ALICE_ROOT/OCDB", Int_t run = 0);
  virtual ~PreClusterizerV2();
  
  void PreClusterize(const char *digitFileName, const char *clusterFileName, const char *method);
  
private:
  
  static const UChar_t nDEs = 156;
  
private:
  
  void CreateMapping(AliMUONVStore* neighbours);
  void FindNeighbours(mpDE &de, UChar_t iPlane);
  
  Int_t LoadDigits(AliMUONVDigitStore* digitStore);
  void ResetPads();
  
  Int_t PreClusterizeFIFO(std::vector<preCluster*> preClusters[nDEs][2], UShort_t nPreClusters[nDEs][2]);
  
  Int_t PreClusterizeRecursive(std::vector<preCluster*> preClusters[nDEs][2], UShort_t nPreClusters[nDEs][2]);
  void AddPad(mpDE &de, UShort_t iPad, preCluster &cl);
  
  Int_t MergePreClusters(std::vector<preCluster*> preClusters[nDEs][2], UShort_t nPreClusters[nDEs][2]);
  void MergePreClusters(preCluster &cl, std::vector<preCluster*> preClusters[2], UShort_t nPreClusters[2],
                        mpDE &de, UChar_t iPlane, preCluster *&mergedCl);
  preCluster* UsePreClusters(preCluster *cl, mpDE &de);
  void MergePreClusters(preCluster &cl1, preCluster &cl2, mpDE &de);
  
  Bool_t AreOverlapping(preCluster &cl1, preCluster &cl2, mpDE &de, Float_t precision);
  Bool_t AreOverlapping(Float_t area1[2][2], Float_t area2[2][2], Float_t precision);
  
  void StorePreClusters(std::vector<preCluster*> preClusters[nDEs][2], UShort_t nPreClusters[nDEs][2],
                        AliMUONVClusterStore *clusterStore);
  
private:
  
  mpDE mpDEs[nDEs];
  TExMap deIndices;
  
  ClassDef(PreClusterizerV2,0)
};

#endif
