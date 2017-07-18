//
//

#include "TSDFVolume.h"

TSDFVolume::TSDFVolume(const Eigen::Vector3i& r)
{
  resolution = r;

  volume.create(resolution.y() * resolution.z(), resolution.x(), sizeof(int));

  const Eigen::Vector3f defaultVolumeSize = Eigen::Vector3f::Constant (3.f); //meters
  const float defaultTruncDist = 0.03f; //meters

  setSize(defaultVolumeSize);
  setTsdfTruncDist(defaultTruncDist);

  reset();
}

void TSDFVolume::setSize(const Eigen::Vector3f& vs)
{
  size = vs;
  setTsdfTruncDist(truncDist);
}

void TSDFVolume::setTsdfTruncDist(float distance)
{
  float cx = size(0) / resolution(0);
  float cy = size(1) / resolution(1);
  float cz = size(2) / resolution(2);

  truncDist = std::max (distance, 2.1f * std::max (cx, std::max (cy, cz)));
}

void TSDFVolume::reset()
{

}
