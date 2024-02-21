#include <map>

#include <TGeoManager.h>

#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/Logger.h"
#include "MCHTracking/TrackExtrap.h"

using o2::mch::TrackExtrap;

//_________________________________________________________________________________________________
void InitFromCCDB(int run, bool local, bool loadGeom, bool loadField)
{
  /// load necessary objects from CCDB and prepare track extrapolation if needed

  auto& ccdb = o2::ccdb::BasicCCDBManager::instance();

  if (local) {

    // setup CCDB from local snapshot
    ccdb.setURL("file://ccdb");

    // in this case we cannot retrieve the timestamp of the run so we must know it
    static const std::map<int, long> timeStamps{{529691, 1669611720310}};
    if (auto ts = timeStamps.find(run); ts != timeStamps.end()) {
      ccdb.setTimestamp(ts->second);
    } else {
      LOG(error) << "unknown time stamp for run " << run;
      exit(-1);
    }

  } else {

    // retrieve the timestamp of the run from the ALICE CCDB server
    auto [tStart, tEnd] = ccdb.getRunDuration(run);
    ccdb.setTimestamp(tEnd);
  }

  // load geometry if requested
  if (loadGeom) {
    auto geom = ccdb.get<TGeoManager>("GLO/Config/GeometryAligned");
  }

  // load magnetic field and prepare track extrapolation if requested
  if (loadField) {
    auto grp = ccdb.get<o2::parameters::GRPMagField>("GLO/Config/GRPMagField");
    o2::base::Propagator::initFieldFromGRP(grp);
    TrackExtrap::setField();
  }
}
