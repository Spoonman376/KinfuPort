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

#ifndef TSDFVolume_h
#define TSDFVolume_h

#include <stdio.h>
#include <vector>

#include <pcl/point_types.h>
#include <pcl/point_cloud.h>

//#include "Eigen.h"

#include <opencv2/opencv.hpp>

#include "TSDFBuffer.h"
#include "PortedGPUFunctions.h"


class TSDFVolume
{

protected:
  /** \brief tsdf volume size in meters */
  Vector3f size;
  
  /** \brief tsdf volume resolution */
  Vector3i resolution;

  /** \brief tsdf volume data container */
  cv::Mat volume;

  /** \brief tsdf truncation distance */
  float truncDist;

public:

  TSDFVolume(const Vector3i& resolution);

  cv::Mat data();

  /** \brief Resets tsdf volume data to uninitialized state */
  void reset();

  /** \brief Sets Tsdf volume size for each dimention
    * \param[in] size size of tsdf volume in meters
    */
  void setSize(const Vector3f& size);

  Vector3f getSize();
  
  /** \brief Sets Tsdf truncation distance. Must be greater than 2 * volume_voxel_size
    * \param[in] distance TSDF truncation distance 
    */
  void setTsdfTruncDist(float distance);

  float getTsdfTruncDist();
};


void integrateTSDFVolume(const cv::Mat& depth, const Intr& intr, const Vector3f& volume_size,
                         const Matrix3frm& Rcurr_inv, const Vector3f& tcurr,
                         float tranc_dist, cv::Mat volume, const tsdf_buffer* buffer,
                         cv::Mat& depthScaled
                         );

void shiftTsdfPointer(short2 ** value, const tsdf_buffer* buffer);

void unpackTsdf(short2 value, float& tsdf, int& weight);

float unpackTsdf (short2 value);

void packTsdf(float tsdf, int weight, short2& value);

#endif /* TSDFVolume_h */
