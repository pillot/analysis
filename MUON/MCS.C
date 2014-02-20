/*
 *  MCS.C
 *  aliroot
 *
 *  Created by Philippe Pillot on 01/10/09.
 *  Copyright 2009 SUBATECH. All rights reserved.
 *
 */

#include <Riostream.h>

#include "TString.h"
#include "TMath.h"
#include "TGeoManager.h"
#include "TGeoNode.h"
#include "TGeoMaterial.h"
#include "TGeoShape.h"

void MCS()
{
  /// check the amplitude of the Multiple Coulomb Scettering in tracking chambers
  
  Double_t OneOverX0MeanCh[10] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  Double_t totalLengthCh[10] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  Double_t OneOverX0Mean = 0.;
  Double_t totalLength = 0.;
  
  // Import TGeo geometry
  if (!gGeoManager) {
    TGeoManager::Import("geometry.root");
    if (!gGeoManager) {
      cout<<"getting geometry from file geometry.root failed"<<endl;
      return;
    }
  }
  
  // Initialize starting point and direction
  Double_t trackXYZIn[3] = {35., -45., -510.};
  Double_t trackXYZOut[3] = {100., -138., -1500.};
  Double_t pathLength = TMath::Sqrt((trackXYZOut[0] - trackXYZIn[0])*(trackXYZOut[0] - trackXYZIn[0])+
				    (trackXYZOut[1] - trackXYZIn[1])*(trackXYZOut[1] - trackXYZIn[1])+
				    (trackXYZOut[2] - trackXYZIn[2])*(trackXYZOut[2] - trackXYZIn[2]));
  Double_t b[3];
  b[0] = (trackXYZOut[0] - trackXYZIn[0]) / pathLength;
  b[1] = (trackXYZOut[1] - trackXYZIn[1]) / pathLength;
  b[2] = (trackXYZOut[2] - trackXYZIn[2]) / pathLength;
  TGeoNode *currentnode = gGeoManager->InitTrack(trackXYZIn, b);
  if (!currentnode) {
    cout<<"start point out of geometry"<<endl;
    return;
  }
  
  Double_t x0 = 0.;  // radiation-length (cm-1)
  Double_t localPathLength = 0.;
  Double_t remainingPathLength = pathLength;
  Int_t chId = 0;
  Bool_t enterCh = kFALSE;
  do {
    // Get material properties
    TGeoMaterial *material = currentnode->GetVolume()->GetMedium()->GetMaterial();
    TString mName(material->GetName());
    mName.ToUpper();
    x0 = material->GetRadLen();
    
    // Get path length within this material
    gGeoManager->FindNextBoundary(remainingPathLength);
    localPathLength = gGeoManager->GetStep() + 1.e-6;
    // Check if boundary within remaining path length. If so, make sure to cross the boundary to prepare the next step
    if (localPathLength >= remainingPathLength) localPathLength = remainingPathLength;
    else {
      currentnode = gGeoManager->Step();
      if (!currentnode) {
        cout<<"navigation failed"<<endl;
	return;
      }
      if (!gGeoManager->IsEntering()) {
        // make another small step to try to enter in new absorber slice
        gGeoManager->SetStep(0.001);
	currentnode = gGeoManager->Step();
	if (!gGeoManager->IsEntering() || !currentnode) {
          cout<<"navigation failed"<<endl;
	  return;
	}
        localPathLength += 0.001;
      }
    }
    
    // check if entering a new chamber
    if (localPathLength > 8.) {
      if(!mName.Contains("AIR")) {
	material->Print();
	cout<<"changing chamber while crossing material"<<endl;
	return;
      }
      if (enterCh) chId++;
      if (chId > 9) break;
      enterCh = kFALSE;
      remainingPathLength -= localPathLength;
      continue;
    } else {
      if (!enterCh) cout<<"-------------------------------------------------------------- enter chamber "<<chId+1<<endl;
      enterCh = kTRUE;
    }
    
    // print parameters
    if(!mName.Contains("AIR")) {
      if(mName.Contains("DIPO")) {
	cout<<"the track cross the dipole"<<endl;
	return;
      }
      //material->Print();
      //cout<<"localPathLength = "<<localPathLength<<endl;
      OneOverX0Mean += localPathLength / x0;
      totalLength += localPathLength;
      OneOverX0MeanCh[chId] += localPathLength / x0;
      totalLengthCh[chId] += localPathLength;
    }
    
    // prepare next step
    remainingPathLength -= localPathLength;
  } while (remainingPathLength > TGeoShape::Tolerance());
  
  // results
  cout<<endl<<endl;
  cout<<"chamber   thickness (cm)    x/X0 (%)"<<endl;
  cout<<"------------------------------------"<<endl;
  for (Int_t i=0; i<10; i++) printf("  %2d          %4.2f            %4.2f\n",i+1,totalLengthCh[i],100.*OneOverX0MeanCh[i]);
  cout<<"------------------------------------"<<endl;
  printf("  tot         %4.1f            %4.1f\n",totalLength,100.*OneOverX0Mean);
  cout<<endl;
  
}

