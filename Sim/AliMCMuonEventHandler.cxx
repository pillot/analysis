/************************************************************************** 
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved  *
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

/* $Id: AliMCMuonEventHandler.cxx 48111 2011-03-06 06:50:22Z morsch $ */
//---------------------------------------------------------------------------------
//                          Class AliMCMuonEventHandler
// This class gives access to MC truth during the analysis.
// Monte Carlo truth is containe in the kinematics tree (produced particles) and 
// the tree of reference hits.
//      
// Origin: Andreas Morsch, CERN, andreas.morsch@cern.ch 
//---------------------------------------------------------------------------------



#include "AliMCMuonEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliPDG.h"
#include "AliTrackReference.h"
#include "AliHeader.h"
#include "AliStack.h"
#include "AliLog.h"

#include <TTree.h>
#include <TFile.h>
#include <TList.h>
#include <TParticle.h>
#include <TString.h>
#include <TClonesArray.h>
#include <TDirectoryFile.h>

ClassImp(AliMCMuonEventHandler)

AliMCMuonEventHandler::AliMCMuonEventHandler() :
AliMCEventHandler()
{
  //
  // Default constructor
  //
  // Be sure to add all particles to the PDG database
  AliPDG::AddParticlesToPdgDataBase();
}

AliMCMuonEventHandler::AliMCMuonEventHandler(const char* name, const char* title) :
AliMCEventHandler(name, title)
{
  //
  // Constructor
  //
}
AliMCMuonEventHandler::~AliMCMuonEventHandler()
{ 
  // Destructor
}

Bool_t AliMCMuonEventHandler::Notify(const char *path)
{
  // Notify about directory change
  // The directory is taken from the 'path' argument
  // Reconnect trees
  TString fileName(path);
  if(fileName.Contains("AliESDsSignal.root")){
    fileName.ReplaceAll("AliESDsSignal.root", "");
  }
  else if(fileName.Contains("AliAOD.root")){
    fileName.ReplaceAll("AliAOD.root", "");
  }
  else if(fileName.Contains("galice.root")){
    // for running with galice and kinematics alone...
    fileName.ReplaceAll("galice.root", "");
  }
  else if (fileName.BeginsWith("root:")) {
    fileName.Append("?ZIP=");
  }
  
  *fPathName = fileName;
  AliInfo(Form("Path: -%s-\n", fPathName->Data()));
  
  ResetIO();
  InitIO("");
  
  // Handle subsidiary handlers
  if (fSubsidiaryHandlers) {
    TIter next(fSubsidiaryHandlers);
    AliMCMuonEventHandler *handler;
    while((handler = (AliMCMuonEventHandler*) next())) {
      TString* spath = handler->GetInputPath();
      if (spath->Contains("merged")) {
	if (! fPathName->IsNull()) {
	  handler->Notify(Form("%s/../.", fPathName->Data()));
	} else {
	  handler->Notify("../");
	}
      }
    }
  }
  
  return kTRUE;
}

