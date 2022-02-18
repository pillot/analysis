#include <initializer_list>
#include <string>
#include <tuple>
#include <utility>

#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>

std::tuple<TFile*, TTreeReader*> LoadData(const char* fileName, const char* treeName);

//_________________________________________________________________________________________________
void check(std::initializer_list<std::pair<std::string, std::string>> filetreenames)
{
  /// check consistency between files

  int nTF0 = -1;

  for (auto names : filetreenames) {

    if (names.first.empty() || names.second.empty()) {
      continue;
    }

    auto [f, r] = LoadData(names.first.c_str(), names.second.c_str());
    int nTF = r->GetEntries(false);

    if (nTF <= 0) {
      printf("* check failed: #TF in %s = %d\n", names.first.c_str(), nTF);
      exit(1);
    }

    if (nTF0 < 0) {
      nTF0 = nTF;
    } else if (nTF != nTF0) {
      printf("* check failed: #TF in %s = %d != %d \n", names.first.c_str(), nTF, nTF0);
      exit(1);
    }

    f->Close();
  }
}

//_________________________________________________________________________________________________
std::tuple<TFile*, TTreeReader*> LoadData(const char* fileName, const char* treeName)
{
  /// open the input file and get the intput tree

  TFile* f = TFile::Open(fileName, "READ");
  if (!f || f->IsZombie()) {
    printf("* check failed: cannot open file %s\n", fileName);
    exit(1);
  }

  TTreeReader* r = new TTreeReader(treeName, f);
  if (r->IsZombie()) {
    printf("* check failed: tree %s not found\n", treeName);
    exit(1);
  }

  return std::make_tuple(f, r);
}
