#include <fstream>

#include <TFile.h>
#include <TTree.h>
#include <TError.h>
#include <TString.h>
#include <TObjString.h>
#include <TGrid.h>
#include <TGridResult.h>
#include <TMap.h>

#include "AliESDEvent.h"
#include "AliESDVertex.h"

void AddVertices(TFile& esdFile, ofstream& outFile);

//---------------------------------------------------------------------------------
void ExtractSPDVertices(TString fileName = "AliESDs.root", TString gridLocation = "")
{
  /// save the SPD vertex of every events from the input ESD file(s) into a binary file with the format below
  /// if a grid location is provided, vertices are taken from all ESD files found in gridLocation/*
  ///
  /// event number
  /// x
  /// y
  /// z
  /// event number
  /// ...

  // open the output file
  ofstream outFile("vertices.in", ios::out | ios::binary);
  if (!outFile.is_open()) {
    Error("ExtractSPDVertices", "opening file vertices.in failed");
    return;
  }

  if (gridLocation.IsNull()) {

    // add vertices from the local ESD file
    TFile* esdFile = TFile::Open(fileName);
    if (!esdFile || !esdFile->IsOpen()) {
      Error("ExtractSPDVertices", "opening file %s failed", fileName.Data());
      return;
    }
    AddVertices(*esdFile, outFile);
    esdFile->Close();

  } else {

    // connect to alien
    if (!TGrid::Connect("alien://")) {
      Error("ExtractSPDVertices","cannot connect to grid");
      return;
    }

    // find all files "fileName" found in gridLocation/*
    TGridResult *res = gGrid->Command(Form("find %s */%s", gridLocation.Data(), fileName.Data()));
    if (!res || res->GetSize() == 0) {
      Error("ExtractSPDVertices", "no file found in %s", gridLocation.Data());
      delete res;
      return;
    }

    TIter resIt(res);
    TMap *map = nullptr;
    while ((map = static_cast<TMap*>(resIt.Next()))) {

      // get the turl of the current ESD file
      TObjString *objs = dynamic_cast<TObjString*>(map->GetValue("turl"));
      if (!objs || !objs->GetString().Length()) {
        Error("ExtractSPDVertices","turl not found in %s... skipping", gridLocation.Data());
        continue;
      }

      // add vertices from the current ESD file
      TFile* esdFile = TFile::Open(objs->GetName());
      if (!esdFile || !esdFile->IsOpen()) {
        Error("ExtractSPDVertices", "opening file %s failed", fileName.Data());
        return;
      }
      AddVertices(*esdFile, outFile);
      esdFile->Close();
    }

    delete res;
  }

  outFile.close();
}

//---------------------------------------------------------------------------------
void AddVertices(TFile& esdFile, ofstream& outFile)
{
  /// add the SPD vertex of every events from esdFile into outFile

  static int event = -1;

  TTree *tree = static_cast<TTree*>(esdFile.Get("esdTree"));
  if (!tree) {
    Error("AddVertices", "esdTree not found in the input file");
    return;
  }

  AliESDEvent* esd = new AliESDEvent();
  esd->ReadFromTree(tree);

  for (int64_t iEvent = 0; iEvent < tree->GetEntries(); ++iEvent) {

    if (tree->GetEntry(iEvent) <= 0) {
      Error("AddVertices", "no entry object found for event %lld", iEvent);
      return;
    }

    ++event;

    double vertex[3] = {0., 0., 0.};
    const AliESDVertex* esdVert = esd->GetVertex(); 
    if (esdVert->GetNContributors() > 0 || !strcmp(esdVert->GetTitle(),"vertexer: smearMC")) {
      esdVert->GetXYZ(vertex);
    }

    outFile.write((char*)&event, sizeof(int));
    outFile.write((char*)&vertex[0], sizeof(double));
    outFile.write((char*)&vertex[1], sizeof(double));
    outFile.write((char*)&vertex[2], sizeof(double));
  }

  delete esd;
}
