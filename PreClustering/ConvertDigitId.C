#include <array>

#include "AliMUONVDigit.h"

#include "MCHMappingInterface/Segmentation.h"

//_________________________________________________________________________________________________
int DigitId2PadId(UInt_t digitId, bool impl4)
{
  /// convert the run2 digit ID to run3 pad ID (using the requested implementation)

  // conversion matrix between the original channels numbering of the RUN2 readout electronics and the final version of the RUN3 DualSAMPA-based readout
  static const std::array<int, 64> refManu2ds_st12 = {
    36, 35, 34, 33, 32, 37, 38, 43, 45, 47, 49, 50, 53, 41, 39, 40,
    63, 62, 61, 60, 59, 58, 56, 57, 54, 55, 52, 51, 48, 46, 44, 42,
    31, 30, 29, 28, 27, 26, 25, 24, 22, 23, 20, 18, 17, 15, 13, 11,
    4, 3, 2, 1, 0, 5, 6, 10, 12, 14, 16, 19, 21, 8, 7, 9};
  static const std::array<int, 64> refManu2ds_st345_v5 = {
    63, 62, 61, 60, 59, 57, 56, 53, 51, 50, 47, 45, 44, 41, 38, 35,
    36, 33, 34, 37, 32, 39, 40, 42, 43, 46, 48, 49, 52, 54, 55, 58,
    7, 8, 5, 2, 6, 1, 3, 0, 4, 9, 10, 15, 17, 18, 22, 25,
    31, 30, 29, 28, 27, 26, 24, 23, 20, 21, 16, 19, 12, 14, 11, 13};

  int deId = AliMUONVDigit::DetElemId(digitId);
  int sampaId = AliMUONVDigit::ManuId(digitId);
  int sampaCh = AliMUONVDigit::ManuChannel(digitId);
  if (impl4) {
    sampaCh = (deId < 500) ? refManu2ds_st12[sampaCh] : refManu2ds_st345_v5[sampaCh];
  }

  auto& segmentation = o2::mch::mapping::segmentation(deId);
  int padId = segmentation.findPadByFEE(sampaId, sampaCh);

  if (!segmentation.isValid(padId)) {
    printf("invalid pad: padID = %d, DE = %d, sampa = %d, channel = %d\n", padId, deId, sampaId, sampaCh);
    return -1;
  }

  return padId;
}

//_________________________________________________________________________________________________
uint32_t PadId2DigitId(int deId, int padId, bool impl4)
{
  /// convert the run3 pad ID to run2 digit ID (using the loaded implementation)

  // cathode number of the bending plane for each DE
  static const std::array<std::vector<int>, 10> bendingCathodes{
    {{0, 1, 0, 1},
     {0, 1, 0, 1},
     {0, 1, 0, 1},
     {0, 1, 0, 1},
     {0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1},
     {0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1},
     {0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1},
     {0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1},
     {0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1},
     {0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1}}};

  // inverse conversion matrix between the original channels numbering of the RUN2 readout electronics and the final version of the RUN3 DualSAMPA-based readout
  static const std::array<int, 64> refDs2manu_st12 = {
    52, 51, 50, 49, 48, 53, 54, 62, 61, 63, 55, 47, 56, 46, 57, 45,
    58, 44, 43, 59, 42, 60, 40, 41, 39, 38, 37, 36, 35, 34, 33, 32,
    4, 3, 2, 1, 0, 5, 6, 14, 15, 13, 31, 7, 30, 8, 29, 9,
    28, 10, 11, 27, 26, 12, 24, 25, 22, 23, 21, 20, 19, 18, 17, 16};
  static const std::array<int, 64> refDs2manu_st345_v5 = {
    39, 37, 35, 38, 40, 34, 36, 32, 33, 41, 42, 62, 60, 63, 61, 43,
    58, 44, 45, 59, 56, 57, 46, 55, 54, 47, 53, 52, 51, 50, 49, 48,
    20, 17, 18, 15, 16, 19, 14, 21, 22, 13, 23, 24, 12, 11, 25, 10,
    26, 27, 9, 8, 28, 7, 29, 30, 6, 5, 31, 4, 3, 2, 1, 0};

  auto& segmentation = o2::mch::mapping::segmentation(deId);
  if (!segmentation.isValid(padId)) {
    printf("invalid pad: padID = %d, DE = %d\n", padId, deId);
    return 0;
  }

  int cathode = bendingCathodes[deId / 100 - 1][deId % 100];
  if (!segmentation.isBendingPad(padId)) {
    cathode = 1 - cathode;
  }
  int manuID = segmentation.padDualSampaId(padId);
  int manuCh = segmentation.padDualSampaChannel(padId);
  if (impl4) {
    manuCh = (deId < 500) ? refDs2manu_st12[manuCh] : refDs2manu_st345_v5[manuCh];
  }

  return (deId) | (manuID << 12) | (manuCh << 24) | (cathode << 30);
}
