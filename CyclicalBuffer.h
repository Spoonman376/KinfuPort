//
//

#ifndef CyclicalBuffer_h
#define CyclicalBuffer_h

#include <stdio.h>


#include "TSDFBuffer.h"

class CyclicalBuffer
{
protected:

  /** \brief structure that contains all TSDF buffer's addresses */
  tsdf_buffer buffer;

  double distanceThreshold;

public:

  /** \brief Constructor for a cubic CyclicalBuffer.
    * \param[in] distance_threshold distance between cube center and target point at which we decide to shift.
    * \param[in] cube_size physical size (in meters) of the volume (here, a cube) represented by the TSDF buffer.
    * \param[in] nb_voxels_per_axis number of voxels per axis of the volume represented by the TSDF buffer.
    */
  CyclicalBuffer(const double threshold, const double cubeSize = 3.f, const int nbVoxelsPerAxis = 512);

  /** \brief Set the physical size represented by the default TSDF volume.
   * \param[in] size_x size of the volume on X axis, in meters.
   * \param[in] size_y size of the volume on Y axis, in meters.
   * \param[in] size_z size of the volume on Z axis, in meters.
   */ 
  void setVolumeSize(const double size_x, const double size_y, const double size_z);

  /** \brief Set the physical size represented by the default TSDF volume.
   * \param[in] size size of the volume on all axis, in meters.
   */
  void setVolumeSize(const double size);

  /** \brief Sets the distance threshold between cube's center and target point that triggers a shift.
    * \param[in] threshold the distance in meters at which to trigger shift.
    */
  void setDistanceThreshold(const double threshold);

  /** \brief get a pointer to the tsdf_buffer structure.
    * \return a pointer to the tsdf_buffer used by cyclical buffer object.
    */
  const tsdf_buffer* getBuffer();




};

#endif /* CyclicalBuffer_h */
