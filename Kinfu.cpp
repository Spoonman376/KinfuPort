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


#include "Kinfu.h"

KinfuTracker::KinfuTracker()
{
	rotationMats.reserve(30000);
	translationVecs.reserve(30000);

  maxICPDistance = 2.5;

  reset();
}


int KinfuTracker::getGlobalTime()
{
  return globalTime;
}


void KinfuTracker::reset()
{
  globalTime = 0;
}

void KinfuTracker::setDepthIntrinsics(float fx_, float fy_, float cx_, float cy_)
{
	fx = fx_;
	fy = fy_;
	cx = (cx_ == -1) ? cols/2-0.5f : cx_;
	cy = (cy_ == -1) ? rows/2-0.5f : cy_;
}



bool KinfuTracker::rgbdodometry(const cv::Mat& image0, const cv::Mat& depth0, const cv::Mat& validMask0,
                                const cv::Mat& image1, const cv::Mat& depth1, const cv::Mat& validMask1,
                                const cv::Mat& cameraMatrix, float minDepth, float maxDepth, float maxDepthDiff,
                                const std::vector<int>& iterCounts, const std::vector<float>& minGradientMagnitudes,
                                const openni::VideoFrameRef& depthFrame, const openni::VideoFrameRef& colorFrame,
                                FramedTransformation* frameTransform
                                )
{
	if (frameTransform != NULL && (frameTransform->flag & frameTransform->ResetFlag)) {
		reset();

    std::string s;

    if (frameTransform->type == frameTransform->DirectApply) {
			Affine3f aff_rgbd(frameTransform->transformation);
			rotationMats[0] = aff_rgbd.linear();
			translationVecs[0] = aff_rgbd.translation();
		}
	}

  const cv::Mat depth = cv::Mat();

  Intr intr(fx, fy, cx, cy, maxIntegrationDistance);

	if (!(frameTransform != NULL && ( frameTransform->flag & frameTransform->IgnoreRegistrationFlag))) {
		//depth_raw.copyTo(depths_curr_[0]);
		bilateralKernel(depth, depthsCur[0], 0.5f / sigmaSpace, 0.5f / sigmaColor);

		if (maxICPDistance > 0)
			truncateDepthKernel(depthsCur[0], maxICPDistance);

		for (int i = 1; i < LEVELS; ++i)
			pyrDownKernel(depthsCur[i-1], depthsCur[i], sigmaColor);

		for (int i = 0; i < LEVELS; ++i)
		{
			createVMap(intr(i), depthsCur[i], vMapsCur[i]);

			computeNormalsEigen(vMapsCur[i], nMapsCur[i]);
		}
	}

	//can't perform more on first frame
	if (globalTime == 0)
	{
		Matrix3frm initialCamRot = rotationMats[0]; //[Ri|ti] - pos of camera, i.e.
		Matrix3frm initialCamRotInv = initialCamRot.inverse();
		Vector3f   initialCamTrans = translationVecs[0]; //transform from camera to global coo space for (i-1)th camera pose

    Eigen33::Mat33 deviceInitalCamRot;
    memcpy(deviceInitalCamRot.data, deviceInitalCamRot.data, sizeof(float) * 9);

    Eigen33::Mat33 deviceInitalCamRotInv;
    memcpy(deviceInitalCamRot.data, deviceInitalCamRotInv.data, sizeof(float) * 9);

    float3 deviceInitialCamTrans(initialCamTrans.x(), initialCamTrans.y(), initialCamTrans.z());

//		float3 device_volume_size = device_cast<const float3>(tsdf_volume_->getSize());
//
//		device::integrateTsdfVolume(depth_raw, intr, device_volume_size, device_initial_cam_rot_inv, device_initial_cam_trans, tsdf_volume_->getTsdfTruncDist(), tsdf_volume_->data(), getCyclicalBufferStructure (), depthRawScaled_);

		for (int i = 0; i < LEVELS; ++i)
			transformMaps(vMapsCur[i], nMapsCur[i], deviceInitalCamRot, deviceInitialCamTrans, vMapsGPrev[i], nMapsGPrev[i]);

		++globalTime;
		return (false);
	}


	// GET PREVIOUS GLOBAL TRANSFORM
	// Previous global rotation
	const int sobelSize = 3;
	const double sobelScale = 1./8;

	cv::Mat depth0Clone = depth0.clone(),
          depth1Clone = depth1.clone();

	// check RGB-D input data
	CV_Assert( !image0.empty() );
	CV_Assert( image0.type() == CV_8UC1 );
	CV_Assert( depth0.type() == CV_32FC1 && depth0.size() == image0.size() );

	CV_Assert( image1.size() == image0.size() );
	CV_Assert( image1.type() == CV_8UC1 );
	CV_Assert( depth1.type() == CV_32FC1 && depth1.size() == image0.size() );

	// check masks
	CV_Assert( validMask0.empty() || (validMask0.type() == CV_8UC1 && validMask0.size() == image0.size()) );
	CV_Assert( validMask1.empty() || (validMask1.type() == CV_8UC1 && validMask1.size() == image0.size()) );

	// check camera params
	CV_Assert( cameraMatrix.type() == CV_32FC1 && cameraMatrix.size() == cv::Size(3,3) );

	// other checks
	CV_Assert( iterCounts.empty() || minGradientMagnitudes.empty() ||
		minGradientMagnitudes.size() == iterCounts.size() );

	vector<int> defaultIterCounts;
	vector<float> defaultMinGradMagnitudes;
	vector<int> const* iterCountsPtr = &iterCounts;
	vector<float> const* minGradientMagnitudesPtr = &minGradientMagnitudes;
	if( iterCounts.empty() || minGradientMagnitudes.empty() )
	{
		defaultIterCounts.resize(4);
		defaultIterCounts[0] = 7;
		defaultIterCounts[1] = 7;
		defaultIterCounts[2] = 7;
		defaultIterCounts[3] = 10;

		defaultMinGradMagnitudes.resize(4);
		defaultMinGradMagnitudes[0] = 12;
		defaultMinGradMagnitudes[1] = 5;
		defaultMinGradMagnitudes[2] = 3;
		defaultMinGradMagnitudes[3] = 1;

		iterCountsPtr = &defaultIterCounts;
		minGradientMagnitudesPtr = &defaultMinGradMagnitudes;
	}




	++globalTime;
	return true;
}
