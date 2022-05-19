#include <algorithm>
#include <iterator>
#include <utility>
#include <vector>

#include <TFile.h>
#include <TTree.h>
#include <TIterator.h>

#include "AliMUONVDigit.h"
#include "AliMUONVDigitStore.h"

#include "Framework/Logger.h"
#include "DataFormatsMID/ColumnData.h"
#include "DataFormatsMID/ROFRecord.h"

using namespace o2;


void digitsToColumnData(const AliMUONVDigitStore& digitStore, std::vector<mid::ColumnData>& midDigits);
int convertFromLegacyDeId(int detElemId);
std::pair<int, int> boardToPattern(int boardId, int cathode);

//_________________________________________________________________________________________________
void ConvertAliRootDigits(const char* inFileName, const char* outFileName = "middigits.root", int nROFsPerTF = 1000)
{
  /// Convert AliRoot trigger digits to O2 MID digits

  // read AliRoot digits
  TFile* inFile = TFile::Open(inFileName);
  if (!inFile || inFile->IsZombie()) {
    LOG(error) << "opening file " << inFileName << " failed";
    exit(-1);
  }
  TTree* inTree = static_cast<TTree*>(inFile->Get("TreeD"));
  if (!inTree) {
    LOG(error) << "tree TreeD not found";
    exit(-1);
  }
  AliMUONVDigitStore* digitStore = AliMUONVDigitStore::Create(*inTree);
  if (!digitStore->Connect(*inTree)) {
    exit(-1);
  }

  // prepare O2 output
  TFile* outFile = TFile::Open(outFileName, "RECREATE");
  TTree* outTree = new TTree("o2sim", "o2sim");
  std::vector<mid::ColumnData> midDigits{};
  outTree->Branch("MIDDigit", &midDigits);
  std::vector<mid::ROFRecord> midROFs{};
  outTree->Branch("MIDROFRecords", &midROFs);

  int nROFsInCurrentTF = 0;
  for (int iEv = 0; iEv < inTree->GetEntries(); ++iEv) {

    digitStore->Clear();

    inTree->GetEntry(iEv);

    auto offset = midDigits.size();
    digitsToColumnData(*digitStore, midDigits);

    midROFs.emplace_back(InteractionRecord{0, static_cast<uint32_t>(iEv)}, mid::EventType::Standard,
                         offset, midDigits.size() - offset);

    if (++nROFsInCurrentTF == nROFsPerTF) {
      outTree->Fill();
      midDigits.clear();
      midROFs.clear();
      nROFsInCurrentTF = 0;
    }
  }

  if (nROFsInCurrentTF > 0) {
    outTree->Fill();
  }

  outFile->cd();
  outTree->SetEntries();
  outTree->Write();
  outFile->Close();
  delete outFile;
  inFile->Close();
  delete inFile;
}

//_________________________________________________________________________________________________
void digitsToColumnData(const AliMUONVDigitStore& digitStore, std::vector<mid::ColumnData>& midDigits)
{
  /// converts digits in the old Run2 format to ColumnData

  TIter next(digitStore.CreateTriggerIterator());
  AliMUONVDigit* digit(nullptr);
  auto offset = midDigits.size();
  while ((digit = static_cast<AliMUONVDigit*>(next()))) {

    // decodes the digit given as a uniqueId in the Run2 format
    auto digitId = digit->GetUniqueID();
    int deId = convertFromLegacyDeId(digitId & 0xFFF);
    int boardId = (digitId & 0xFFF000) >> 12;
    int cathode = (digitId & 0x40000000) >> 30;
    auto [icolumn, iline] = boardToPattern(boardId, cathode);
    int channel = (digitId & 0x3F000000) >> 24;

    // get the column data if already there, or create it
    auto itFirst = std::next(midDigits.begin(), offset);
    auto itDigit = std::find_if(itFirst, midDigits.end(), [deId, cId = icolumn](const mid::ColumnData& d) {
      return d.deId == deId && d.columnId == cId;
    });
    if (itDigit == midDigits.end()) {
      itDigit = midDigits.emplace(midDigits.end(), mid::ColumnData{(uint8_t) deId, (uint8_t)icolumn});
    }

    // update the pattern
    itDigit->patterns[(cathode == 1) ? 4 : iline] |= (1 << channel);
  }
}

//_________________________________________________________________________________________________
int convertFromLegacyDeId(int detElemId)
{
  /// converts the detection element ID (Run2 format) into the new ID (Run3 format)
  int ich = (detElemId / 100 - 11);
  int irpc = (detElemId % 100 + 4) % 18;
  if (irpc >= 9) {
    irpc = 17 - irpc;
    ich += 4;
  }
  return ich * 9 + irpc;
}

//_________________________________________________________________________________________________
std::pair<int, int> boardToPattern(int boardId, int cathode)
{
  /// converts old Run2 local board Id into the new format

  static const int endBoard[7] = {16, 38, 60, 76, 92, 108, 117};
  static const std::vector<int> lines[3] = {{3, 19, 41, 63, 79, 95, 5, 21, 43, 65, 81, 97, 7, 23, 45,
                                             67, 83, 99, 27, 49, 69, 85, 101, 9, 31, 53, 71, 87, 103,
                                             13, 35, 57, 73, 89, 105, 15, 37, 59, 75, 91, 107},
                                            {8, 24, 46, 28, 50, 10, 32, 54},
                                            {25, 47, 29, 51, 11, 33, 55}};

  int halfBoardId = (boardId - 1) % 117 + 1;
  int icolumn = 0;
  for (; icolumn < 7; ++icolumn) {
    if (halfBoardId <= endBoard[icolumn]) {
      break;
    }
  }

  if (cathode == 1) {
    return std::make_pair(icolumn, 0);
  }

  for (int iline = 0; iline < 3; ++iline) {
    if (std::count(lines[iline].begin(), lines[iline].end(), halfBoardId) > 0) {
      return std::make_pair(icolumn, iline + 1);
    }
  }

  return std::make_pair(icolumn, 0);
}
