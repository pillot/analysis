#ifndef ALIVERTEXFROMFILE_H
#define ALIVERTEXFROMFILE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

// Generator for vertices taken from a binary file

#include "AliVertexGenerator.h"

class AliVertexFromFile: public AliVertexGenerator {
 public:
  AliVertexFromFile();
  AliVertexFromFile(const char* fileName, Int_t eventsPerEntry = 1);
  virtual ~AliVertexFromFile();

  virtual TVector3 GetVertex(Bool_t& isGood);
  virtual Float_t GetLastVertexTime() const {return 0.;}
  
 private:
  AliVertexFromFile(const AliVertexFromFile &vgf);
  AliVertexFromFile & operator=(const AliVertexFromFile &);

  void ReadNextVertex();

  std::ifstream*   fFile;           //! file with vertices
  TVector3         fVertex;         //! position of the current vertex
  Int_t            fEvent;          //! current event number
  Int_t            fEventsPerEntry; // number of events with same vertex

  ClassDef(AliVertexFromFile, 1)     // generator for vertices taken from a file
};

#endif














