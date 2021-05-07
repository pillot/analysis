#include <cstring>
#include <iostream>
#include <fstream>
#include <limits>
#include <vector>

#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TIterator.h>

#include "AliMUONVDigit.h"
#include "AliMUONVDigitStore.h"

#include "DataFormatsMCH/Digit.h"

using namespace o2::mch;
using namespace std;

//------------------------------------------------------------------
void ConvertAliRootDigitsV2(TString inFileName, TString outFileName = "digits.v3.in")
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

  // write the file format
  uint64_t fileFormat = 2305844383603244847;
  out.write((char*)&fileFormat, sizeof(uint64_t));

  // loop over events
  Long64_t nEvents = treeD->GetEntries();
  std::vector<Digit> digits;
  int sizeofDigit = sizeof(Digit);
  for (Long64_t iEv = 0; iEv < nEvents; ++iEv) {
    
    treeD->GetEntry(iEv);

    // fill the vector of digits
//    auto digitIt = digitStore->CreateIterator();        // ! this changes the order of the output digits !
    auto digitIt = digitStore->CreateTrackerIterator(); // ! this changes the order of the output digits !
    AliMUONVDigit *digit(nullptr);
    while ((digit = static_cast<AliMUONVDigit*>(digitIt->Next()))) {
      if (digit->Charge() <= 0) continue;
//      if (digit->Charge() > 610) cout << digit->Charge() << endl;
//      unsigned long charge = (static_cast<double>(digit->Charge()) / 1024) * static_cast<double>(std::numeric_limits<unsigned long>::max());
//      unsigned long charge = (static_cast<double>(digit->Charge()) / 0x8000000000) * static_cast<double>(std::numeric_limits<unsigned long>::max());
      float charge = digit->Charge();
      uint32_t adc(0);
      std::memcpy(&adc, &charge, sizeof(charge));
      digits.push_back({digit->DetElemId(), static_cast<int>(digit->GetUniqueID()), adc, 0, 0, digit->IsSaturated()});
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
