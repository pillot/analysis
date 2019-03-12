#include <fstream>
#include <vector>

#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TIterator.h>

#include "AliMUONVDigit.h"
#include "AliMUONVDigitStore.h"

#include "MCHBase/DigitBlock.h"

using namespace o2::mch;

//------------------------------------------------------------------
void ConvertAliRootDigits(TString inFileName, TString outFileName = "digits.in")
{
  /// Convert AliRoot digits to O2 digits
  /// saved in a binary file with the following format:
  ///
  /// DigitBlock 1
  /// DigitStruct 1
  /// DigitStruct 2
  /// ...
  /// DigitStruct n
  /// DigitBlock 2
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

  // digit block header
  DigitBlock digitBlock;
  digitBlock.header.fType = 2003;
  digitBlock.header.fRecordWidth = sizeof(DigitStruct);
  digitBlock.header.fNrecords = 0;

  // loop over events
  Long64_t nEvents = treeD->GetEntries();
  std::vector<DigitStruct> digits;
  for (Long64_t iEv = 0; iEv < nEvents; ++iEv) {
    
    treeD->GetEntry(iEv);

    // fill the vector of digits
    auto digitIt = digitStore->CreateIterator();
    AliMUONVDigit *digit(nullptr);
    while ((digit = static_cast<AliMUONVDigit*>(digitIt->Next()))) {
      if (digit->Charge() <= 0) continue;
      digits.push_back({digit->GetUniqueID(), 0, static_cast<uint16_t>(digit->ADC())});
    }

    // write digit block header
    digitBlock.header.fNrecords = digits.size();
    out.write((char*)&digitBlock, sizeof(DigitBlock));

    // write digits if any
    if (digitBlock.header.fNrecords > 0) {
      out.write((char*)digits.data(), digitBlock.header.fNrecords * digitBlock.header.fRecordWidth);
    }
    
    digits.clear();
    digitStore->Clear();
  }
}
