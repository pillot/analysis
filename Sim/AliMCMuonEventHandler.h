// -*- mode: C++ -*-
#ifndef ALIMCMUONEVENTHANDLER_H
#define ALIMCMUONEVENTHANDLER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/* $Id: AliMCMuonEventHandler.h 48111 2011-03-06 06:50:22Z morsch $ */

//-------------------------------------------------------------------------
//                          Class AliMCEvent
// This class gives access to MC truth during the analysis.
// Monte Carlo truth is contained in the kinematics tree (produced particles) and 
// the tree of reference hits.
//      
// Origin: Andreas Morsch, CERN, andreas.morsch@cern.ch 
//-------------------------------------------------------------------------
#include "AliMCEventHandler.h"
#include "AliHeader.h"
#include <TExMap.h>

class TFile;
class TTree;
class TList;

class TParticle;
class TString;
class TClonesArray;
class TDirectoryFile;

class AliMCEvent;



class AliMCMuonEventHandler : public AliMCEventHandler
{
public:
  
  
  AliMCMuonEventHandler();
  AliMCMuonEventHandler(const char* name, const char* title);
  virtual ~AliMCMuonEventHandler();
  
  virtual Bool_t       Notify() { return AliVEventHandler::Notify(); };
  virtual Bool_t       Notify(const char* path);
  
  
private:
  AliMCMuonEventHandler(const AliMCMuonEventHandler& handler);             
  AliMCMuonEventHandler& operator=(const AliMCMuonEventHandler& handler);  
private:
  
  ClassDef(AliMCMuonEventHandler,1)  //MC Truth EventHandler class
};
#endif 

