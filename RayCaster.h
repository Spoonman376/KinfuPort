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

#ifndef RayCaster_h
#define RayCaster_h

#include <stdio.h>

#include "Eigen.h"

#include <opencv2/opencv.hpp>

#include "TSDFVolume.h"


class RayCaster
{
  enum { CTA_SIZE_X = 32, CTA_SIZE_Y = 8 };

  Matrix3frm Rcurr;
  Vector3f tcurr;

  float time_step;
  Vector3f volume_size;

  Vector3f cell_size;
  int cols, rows;

  cv::Mat volume; // <short2>

  Intr intr;

  cv::Mat *nmap; // <float>
  cv::Mat *vmap; // <float>

  Vector3f getRayNext(int x, int y) const;

  bool checkInds(const Vector3i& g) const;

  float readTsdf(int x, int y, int z, const tsdf_buffer buffer) const;

  Vector3i getVoxel(Vector3f point) const;

  float interpolateTrilineary(const Vector3f& origin, const Vector3f& dir, float time, const tsdf_buffer buffer) const;
  
  float interpolateTrilineary(const Vector3f& point, const tsdf_buffer buffer) const;

  void operator () (const tsdf_buffer buffer) const;
};

float getMinTime(const Vector3f& volume_max, const Vector3f& origin, const Vector3f& dir);

float getMaxTime(const Vector3f& volume_max, const Vector3f& origin, const Vector3f& dir);





#endif /* RayCaster_h */
