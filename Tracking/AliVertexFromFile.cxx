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

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Generator for vertices taken from a binary file                           //
//                                                                           //
// The file name is passed as argument to the constructor.                   //
// If a second argument is given, this determines the number                 //
// of events for which the same vertex is used.                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <fstream>

#include "AliVertexFromFile.h"
#include "AliLog.h"


ClassImp(AliVertexFromFile)


//_____________________________________________________________________________
AliVertexFromFile::AliVertexFromFile() :
  fFile(NULL),
  fVertex(0,0,0),
  fEvent(0),
  fEventsPerEntry(0)
{
// default constructor: initialize data members
}

//_____________________________________________________________________________
AliVertexFromFile::AliVertexFromFile(const char* fileName, Int_t eventsPerEntry) :
  fFile(NULL),
  fVertex(0,0,0),
  fEvent(0),
  fEventsPerEntry(eventsPerEntry)
{
// main constructor:
// fileName is the name of the binary file containing the vertices
// eventsPerEntry is the number of events for which the same vertex is used
  fFile = new std::ifstream(fileName,ios::binary);
  if (!fFile || !fFile->is_open()) {
    AliFatal(Form("could not open file %s", fileName));
    return;
  }
}

//_____________________________________________________________________________
AliVertexFromFile::~AliVertexFromFile()
{
// clean up
  if (fFile) fFile->close();
  delete fFile;
}

//_____________________________________________________________________________
TVector3 AliVertexFromFile::GetVertex()
{
// get the vertex from the binary file
  if (fEventsPerEntry < 2 || fEvent%fEventsPerEntry == 0) {
    ReadNextVertex();
  }
  ++fEvent;
  return fVertex;
}

//_____________________________________________________________________________
void AliVertexFromFile::ReadNextVertex()
{
// read the next vertex in the input file
// set it to (0,0,0) when reaching end of file
  int event(-1);
  if (!fFile || !fFile->read(reinterpret_cast<char*>(&event), sizeof(int))) {
    AliError(Form("no more vertex to read"));
    fVertex.SetXYZ(0., 0., 0.);
  } else {
    fFile->read(reinterpret_cast<char*>(&fVertex[0]), sizeof(double));
    fFile->read(reinterpret_cast<char*>(&fVertex[1]), sizeof(double));
    fFile->read(reinterpret_cast<char*>(&fVertex[2]), sizeof(double));
  }
}
