#include <cstring>
#include <iostream>
#include <fstream>
#include <limits>
#include <vector>

#include <TROOT.h>
#include <TSystem.h>
#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TIterator.h>

#include "AliMUONVDigit.h"
#include "AliMUONVDigitStore.h"

using namespace std;

struct DigitD0 {
  int32_t tfTime{0};      /// time since the beginning of the time frame, in bunch crossing units
  uint16_t nofSamples{0}; /// number of samples in the signal + saturated bit
  int detID{0};           /// ID of the Detection Element to which the digit corresponds to
  int padID{0};           /// PadIndex to which the digit corresponds to
  uint32_t adc{0};        /// Amplitude of signal

  void setNofSamples(uint16_t n) { nofSamples = (nofSamples & 0x8000) + (n & 0x7FFF); }
  uint16_t getNofSamples() const { return (nofSamples & 0x7FFF); }

  void setSaturated(bool sat) { nofSamples = sat ? nofSamples | 0x8000 : nofSamples & 0x7FFF; }
  bool isSaturated() const { return ((nofSamples & 0x8000) > 0); }
};

int DigitId2PadId(UInt_t digitId, bool impl4);

//------------------------------------------------------------------
void ConvertAliRootDigitsV2(TString inFileName, TString outFileName = "digits.v3.in", bool impl4 = true)
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

  // load the digitId converter linked with the requested mapping implementation
  gSystem->Load(impl4 ? "libO2MCHMappingImpl4" : "libO2MCHMappingImpl3");
  gROOT->LoadMacro("/Users/PILLOT/Work/Alice/Macros/PreClustering/ConvertDigitId.C++");

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
  std::vector<DigitD0> digits;
  int sizeofDigit = sizeof(DigitD0);
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
      int padId = DigitId2PadId(digit->GetUniqueID(), impl4);
      if (padId < 0) {
        continue;
      }
      digits.push_back({0, 0, digit->DetElemId(), padId, adc});
      digits.back().setSaturated(digit->IsSaturated());
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
