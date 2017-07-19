//
//

#include "CyclicalBuffer.h"


/** \brief Constructor for a cubic CyclicalBuffer.
  * \param[in] distance_threshold distance between cube center and target point at which we decide to shift.
  * \param[in] cube_size physical size (in meters) of the volume (here, a cube) represented by the TSDF buffer.
  * \param[in] nb_voxels_per_axis number of voxels per axis of the volume represented by the TSDF buffer.
  */
CyclicalBuffer::CyclicalBuffer(const double threshold, const double cubeSize, const int nbVoxelsPerAxis)
{
  distanceThreshold = threshold;
  buffer.volume_size.x() = cubeSize;
  buffer.volume_size.y() = cubeSize;
  buffer.volume_size.z() = cubeSize;
  buffer.voxels_size.x() = nbVoxelsPerAxis;
  buffer.voxels_size.y() = nbVoxelsPerAxis;
  buffer.voxels_size.z() = nbVoxelsPerAxis;
}


/** \brief Set the physical size represented by the default TSDF volume.
 * \param[in] size_x size of the volume on X axis, in meters.
 * \param[in] size_y size of the volume on Y axis, in meters.
 * \param[in] size_z size of the volume on Z axis, in meters.
 */
void CyclicalBuffer::setVolumeSize(const double size_x, const double size_y, const double size_z)
{
  buffer.volume_size.x() = size_x;
  buffer.volume_size.y() = size_y;
  buffer.volume_size.z() = size_z;
}

/** \brief Set the physical size represented by the default TSDF volume.
 * \param[in] size size of the volume on all axis, in meters.
 */
void CyclicalBuffer::setVolumeSize(const double size)
{
  buffer.volume_size.x() = size;
  buffer.volume_size.y() = size;
  buffer.volume_size.z() = size;
}

/** \brief Sets the distance threshold between cube's center and target point that triggers a shift.
  * \param[in] threshold the distance in meters at which to trigger shift.
  */
void CyclicalBuffer::setDistanceThreshold(const double threshold)
{ 
  distanceThreshold = threshold;
}


const tsdf_buffer* CyclicalBuffer::getBuffer()
{
  return &buffer;
}



