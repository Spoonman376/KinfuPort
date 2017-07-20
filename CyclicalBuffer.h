/*
 * Software License Agreement (BSD License)
 *
 *  Point Cloud Library (PCL) - www.pointclouds.org
 *  Copyright (c) 2011, Willow Garage, Inc.
 *
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 *   * Neither the name of Willow Garage, Inc. nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 *  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 *  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 *  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 */

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
