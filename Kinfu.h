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


#ifndef Kinfu_h
#define Kinfu_h

#include <stdio.h>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <opencv2/opencv.hpp>
#include <openni2/OpenNI.h>

//#include <pcl/common/transforms.h>

#include "FramedTransformation.h"
#include "PortedGPUFunctions.h"


typedef Eigen::Transform<float, 3, Eigen::Affine> Affine3f;
typedef Eigen::Matrix<float, 3, 3, Eigen::RowMajor> Matrix3frm;
typedef Eigen::Vector3f Vector3f;


class KinfuTracker
{
protected:

  int globalTime;

  // 
  int LEVELS = 3;

  const float sigmaColor = 30;      //in mm
  const float sigmaSpace = 4.5;     // in pixels


  /** Vector of camera rotation matrices for each moment of time. */
  std::vector<Matrix3frm> rotationMats;
  
  /** Vector camera translations for each moment of time. */
  std::vector<Vector3f> translationVecs;

  /** Depth pyramid */
  std::vector<cv::Mat> depthsCur;

  /** Vertex maps pyramid for current frame in current coordinate space */
  std::vector<cv::Mat> vMapsCur;

  /** Vertex maps pyramid for previous frame in global coordinate space. */
  std::vector<cv::Mat> vMapsGPrev;

  /** Vertex maps pyramid for current frame in global coordinate space. */
  std::vector<cv::Mat> vMapGCur;


  /** Normal maps pyramid for current frame in current coordinate space. */
  std::vector<cv::Mat> nMapsCur;

  /** Normal maps pyramid for previous frame in global coordinate space. */
  std::vector<cv::Mat> nMapsGPrev;

  /** Normal maps pyramid for current frame in global coordinate space. */
  std::vector<cv::Mat> nMapsGCur;


  /** \brief Intrinsic parameters of depth camera. */
  float fx, fy, cx, cy;
  
  float maxIntegrationDistance, maxICPDistance; // meters

  /** \br-ief Height and width of input depth image. */
  int rows, cols;

  void reset();
  void setDepthIntrinsics(float fx_, float fy_, float cx_ = -1, float cy_ = -1);;


public:
  KinfuTracker();

  int getGlobalTime();

  bool rgbdodometry(const cv::Mat& image0, const cv::Mat& _depth0, const cv::Mat& validMask0,
                    const cv::Mat& image1, const cv::Mat& _depth1, const cv::Mat& validMask1,
                    const cv::Mat& cameraMatrix, float minDepth, float maxDepth, float maxDepthDiff,
                    const std::vector<int>& iterCounts, const std::vector<float>& minGradientMagnitudes,
                    const openni::VideoFrameRef& depth, const openni::VideoFrameRef& pcolor, FramedTransformation* frame_ptr
                    );

};



#endif /* Kinfu_h */
