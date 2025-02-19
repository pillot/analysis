#include <fstream>

#include "TGrid.h"
#include "TString.h"
#include "TFileMerger.h"

void merge(const char* fileList, const char* outFileName)
{
  /// merge files in the input fileList into a file named outFileName

  TFileMerger fm(false);

  // add the files to merge
  std::ifstream in;
  in.open(fileList);
  TString line;
  bool connectedToAliEn = false;
  while (in.good()) {
    in >> line;
    if (line.Length() == 0) {
      continue;
    }
    if (line.BeginsWith("alien:") && !connectedToAliEn) {
      printf("Connecting to AliEn...");
      TGrid::Connect("alien:");
      connectedToAliEn = true; // Only try once
    }
    fm.AddFile(line.Data());
  }

  // nothing found - skip this output
  if (!fm.GetMergeList() || !fm.GetMergeList()->GetSize()) {
    printf("no file found\n");
    exit(1);
  }

  // open output file
  fm.OutputFile(outFileName);

  // merge
  printf("merging %d files...\n", fm.GetMergeList()->GetSize());
  if (!fm.Merge()) {
    printf("could not merge all files\n");
    exit(1);
  }
}
