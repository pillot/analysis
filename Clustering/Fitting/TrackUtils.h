#include "MCHBase/TrackBlock.h"

using o2::mch::TrackParamStruct;

//_________________________________________________________________________________________________
struct TrackLite {
  TrackParamStruct param;
  double rAbs;
  double dca;
  double pUncorr;
  int nClusters;
  double chi2;
  double matchChi2;
};
