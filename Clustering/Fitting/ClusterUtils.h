#include <cmath>

#include <TGeoManager.h>

#include "DataFormatsMCH/Cluster.h"
#include "MathUtils/Cartesian.h"
#include "MCHGeometryTransformer/Transformations.h"

using o2::mch::Cluster;

//_________________________________________________________________________________________________
Cluster MakeCluster(uint32_t uid, float x, float y, float ex = 0.2, float ey = 0.2)
{
  /// create a cluster in the global coordinate system from a local position on a DE
  /// the DE ID is extracted from the uid, which must be defined as in Cluster::buildUniqueId(...)
  /// the geometry must be loaded (e.g. from the CCDB) prior to use this function

  static auto transformation = o2::mch::geo::transformationFromTGeoManager(*gGeoManager);

  auto de = Cluster::getDEId(uid);
  o2::math_utils::Point3D<float> local{x, y, 0.};
  auto global = transformation(de)(local);

  return Cluster{global.x(), global.y(), global.z(), ex, ey, uid, 0, 0};
}
