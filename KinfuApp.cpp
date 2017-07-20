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

#include "KinfuApp.h"


KinfuApp::KinfuApp()
{
  initalizeOpenNI();
}

KinfuApp::~KinfuApp()
{
  closeOpenNI();
}

void KinfuApp::mainLoop()
{
  device.open(oniFileName.c_str());
  cout << OpenNI::getExtendedError() << endl;

  depthStream.create(device, SENSOR_DEPTH);
  colourStream.create(device, SENSOR_COLOR);

  device.getPlaybackControl()->setSpeed(-1);

  int numFrames = device.getPlaybackControl()->getNumberOfFrames(depthStream);
  frameCtr = 0;

  depthStream.start();
  colourStream.start();

  for (int i = 0; i < numFrames; ++i) {
    //get a frame
    VideoFrameRef depthFrame, colourFrame;

    depthStream.readFrame(&depthFrame);
    colourStream.readFrame(&colourFrame);

    frameCtr++;

    cout << "Processing frame " << frameCtr << " out of " << numFrames << endl;

    // execute
    //execute(depthFrame, colourFrame);
  }
}


void KinfuApp::execute(VideoFrameRef depth, VideoFrameRef colour)
{

  processFramedTransformation(depth.getFrameIndex());

  // normalize color of grayImage1_
  cv::Scalar scale0 = cv::mean( grayImage0 );
  cv::Scalar scale1 = cv::mean( grayImage1 );
  grayImage1.convertTo( grayImage1, -1, scale0.val[ 0 ] / scale1.val[ 0 ] );

  vector<int> iterCounts(4);
  iterCounts[0] = 7;
  iterCounts[1] = 7;
  iterCounts[2] = 7;
  iterCounts[3] = 10;

  vector<float> minGradMagnitudes(4);
  minGradMagnitudes[0] = 12;
  minGradMagnitudes[1] = 5;
  minGradMagnitudes[2] = 3;
  minGradMagnitudes[3] = 1;
  const float minDepth = 0.f; //in meters
  const float maxDepth = 7.5f; //in meters
  const float maxDepthDiff = 0.07f; //in meters

  float vals[] = {camera.fx, 0.0,        camera.cx,
                  0.0,        camera.fy, camera.cy,
                  0.0,       0.0,        1.0       };

  const cv::Mat cameraMatrix = cv::Mat(3,3,CV_32FC1,vals);
  const cv::Mat distCoeff(1,5,CV_32FC1,cv::Scalar(0));

  //run kinfu algorithm
  bool hasImage = kinfu.rgbdodometry(grayImage0, depthFilter0, cv::Mat(),
                                     grayImage1, depthFilter1, cv::Mat(),
                                     cameraMatrix, minDepth, maxDepth, maxDepthDiff,
                                     iterCounts, minGradMagnitudes,
                                     depth, colour, &framedTransformation
                                     );


  // Add Odometry
//  if ( kinfu_->getGlobalTime() > 0 ) {
//    // global_time_ == 0 only when lost and reset, in this case, we lose one frame
//    kinfuTraj.data_.push_back(FramedTransformation(kinfuTraj.data.size(), kinfu_->getGlobalTime() - 1, frame_id_, kinfu_->getCameraPose().matrix()));
//  }

}


void KinfuApp::processFramedTransformation(int frameId)
{
  if(frameId > 0 && frameId % fragmentRate == fragmentStart + 1)
    framedTransformation.flag = FramedTransformation::ResetFlag;

  else if(frameId > 0 && frameId % fragmentRate == fragmentStart)
    framedTransformation.flag = FramedTransformation::SavePointCloudFlag;

  else
    framedTransformation.flag = 0;

}




bool KinfuApp::initalizeOpenNI()
{
  if (OpenNI::initialize() != STATUS_OK) {
    cout << "Failed to Initialize OpenNI" << endl;
    cout << OpenNI::getExtendedError() << endl;
    return false;
  }
  return true;
}

void KinfuApp::closeOpenNI()
{
  OpenNI::shutdown();
}
