#include <iostream>
#include <fstream>
#include <vector>

#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TIterator.h>

#include "AliMUONVDigit.h"
#include "AliMUONVDigitStore.h"

#include "MCHBase/Digit.h"

using namespace o2::mch;
using namespace std;

//------------------------------------------------------------------
void ConvertAliRootDigitsV2(TString inFileName, TString outFileName = "digits.v2.in")
{
  /// Convert AliRoot digits to O2 digits
  /// saved in a binary file with the following format:
  ///
  /// Number of digits in event 1
  /// Digit 1
  /// Digit 2
  /// ...
  /// Digit n
  /// Number of digits in event 2
  /// ...

  // read digits
  TFile* digitFile = TFile::Open(inFileName);
  if (!digitFile) return;
  TTree* treeD = static_cast<TTree*>(digitFile->Get("TreeD"));
  if (!treeD) return;
  AliMUONVDigitStore* digitStore = AliMUONVDigitStore::Create(*treeD);
  if (!digitStore->Connect(*treeD)) return;
  
  // output file
  ofstream out(outFileName.Data(),ios::out | ios::binary);
  if (!out.is_open()) return;

  // loop over events
  Long64_t nEvents = treeD->GetEntries();
  std::vector<Digit> digits;
  int sizeofDigit = sizeof(Digit);
  for (Long64_t iEv = 0; iEv < nEvents; ++iEv) {
    
    treeD->GetEntry(iEv);

    // fill the vector of digits
    auto digitIt = digitStore->CreateIterator();
    AliMUONVDigit *digit(nullptr);
    while ((digit = static_cast<AliMUONVDigit*>(digitIt->Next()))) {
      if (digit->Charge() <= 0) continue;
//      if (digit->Charge() > 610) cout << digit->Charge() << endl;
      unsigned long charge = (digit->Charge() / 1024) * ULONG_MAX;
      Digit::Time time{static_cast<uint64_t>(digit->IsSaturated() ? 1 : 0)};
      digits.push_back({digit->DetElemId(), static_cast<int>(digit->GetUniqueID()), charge, time});
    }

    // write the number of digits
    int nDigits = digits.size();
    out.write((char*)&nDigits, sizeof(int));

    // write digits if any
    if (nDigits > 0) {
      out.write((char*)digits.data(), nDigits * sizeofDigit);
    }
    
    digits.clear();
    digitStore->Clear();
  }
}
