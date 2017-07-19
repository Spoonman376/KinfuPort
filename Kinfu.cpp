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
	const Vector3i volumeResolution (VOLUME_X, VOLUME_Y, VOLUME_Z);
  tsdfVolume = new TSDFVolume(volumeResolution);
  cyclical = new CyclicalBuffer(DISTANCE_THRESHOLD, VOLUME_SIZE, VOLUME_X);

	rotationMats.reserve(30000);
	translationVecs.reserve(30000);

  maxICPDistance = 2.5;

	depthsCur.resize(LEVELS);
	vMapsGCur.resize(LEVELS);
	nMapsGCur.resize(LEVELS);

	vMapsGPrev.resize(LEVELS);
	nMapsGPrev.resize(LEVELS);

	vMapsCur.resize(LEVELS);
	nMapsCur.resize(LEVELS);

  reset();
}

KinfuTracker::~KinfuTracker()
{
  if (tsdfVolume != nullptr)
    delete tsdfVolume;
}


int KinfuTracker::getGlobalTime()
{
  return globalTime;
}

Eigen::Affine3f KinfuTracker::getCameraPose(int time)
{
	if (time > (int)rotationMats.size() || time < 0)
		time = rotationMats.size() - 1;

	Eigen::Affine3f aff;
	aff.linear() = rotationMats[time];
	aff.translation() = translationVecs[time];
	return aff;
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

void KinfuTracker::preprocessDepth(cv::Mat depth0, cv::Mat depth1, const cv::Mat& validMask0,
                     const cv::Mat& validMask1, float minDepth, float maxDepth
                     )
{
	CV_DbgAssert(depth0.size() == depth1.size());

	for( int y = 0; y < depth0.rows; y++ )
	{
		for( int x = 0; x < depth0.cols; x++ )
		{
			float& d0 = depth0.at<float>(y,x);
			if( !cvIsNaN(d0) && (d0 > maxDepth || d0 < minDepth || d0 <= 0 || (!validMask0.empty() && !validMask0.at<uchar>(y,x))) )
				d0 = std::numeric_limits<float>::quiet_NaN();

			float& d1 = depth1.at<float>(y,x);
			if( !cvIsNaN(d1) && (d1 > maxDepth || d1 < minDepth || d1 <= 0 || (!validMask1.empty() && !validMask1.at<uchar>(y,x))) )
				d1 = std::numeric_limits<float>::quiet_NaN();
		}
	}
}



void KinfuTracker::buildPyramids(const cv::Mat& image0, const cv::Mat& image1, const cv::Mat& depth0,
                                 const cv::Mat& depth1, const cv::Mat& cameraMatrix, int sobelSize, double sobelScale,
                                 const vector<float>& minGradMagnitudes, vector<cv::Mat>& pyramidImage0,
                                 vector<cv::Mat>& pyramidDepth0, vector<cv::Mat>& pyramidImage1,
                                 vector<cv::Mat>& pyramidDepth1, vector<cv::Mat>& pyramid_dI_dx1,
                                 vector<cv::Mat>& pyramid_dI_dy1, vector<cv::Mat>& pyramidTexturedMask1,
                                 vector<cv::Mat>& pyramidCameraMatrix
                                 )
{
	const int pyramidMaxLevel = (int)minGradMagnitudes.size() - 1;

	buildPyramid(image0, pyramidImage0, pyramidMaxLevel);
	buildPyramid(image1, pyramidImage1, pyramidMaxLevel);

	pyramid_dI_dx1.resize(pyramidImage1.size());
	pyramid_dI_dy1.resize(pyramidImage1.size());
	pyramidTexturedMask1.resize(pyramidImage1.size());

	pyramidCameraMatrix.reserve(pyramidImage1.size());

	cv::Mat cameraMatrix_dbl;
	cameraMatrix.convertTo(cameraMatrix_dbl, CV_64FC1);

	for( size_t i = 0; i < pyramidImage1.size(); i++ )
	{
		Sobel(pyramidImage1[i], pyramid_dI_dx1[i], CV_16S, 1, 0, sobelSize);
		Sobel(pyramidImage1[i], pyramid_dI_dy1[i], CV_16S, 0, 1, sobelSize);

		const cv::Mat& dx = pyramid_dI_dx1[i];
		const cv::Mat& dy = pyramid_dI_dy1[i];

		cv::Mat texturedMask(dx.size(), CV_8UC1, cvScalar(0));
		const float minScalesGradMagnitude2 = (float)((minGradMagnitudes[i] * minGradMagnitudes[i]) / (sobelScale * sobelScale));
		for(int y = 0; y < dx.rows; y++)
		{
			for(int x = 0; x < dx.cols; x++)
			{
				float m2 = (float)(dx.at<short>(y,x)*dx.at<short>(y,x) + dy.at<short>(y,x)*dy.at<short>(y,x));
				if( m2 >= minScalesGradMagnitude2 )
					texturedMask.at<uchar>(y,x) = 255;
			}
		}
		pyramidTexturedMask1[i] = texturedMask;
		cv::Mat levelCameraMatrix = i == 0 ? cameraMatrix_dbl : 0.5f * pyramidCameraMatrix[i-1];
		levelCameraMatrix.at<double>(2,2) = 1.;
		pyramidCameraMatrix.push_back( levelCameraMatrix);
	}

	buildPyramid(depth0, pyramidDepth0, pyramidMaxLevel);
	buildPyramid(depth1, pyramidDepth1, pyramidMaxLevel);
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

      cv::Affine3f aff;
		}
	}

  const cv::Mat depth = cv::Mat(depthFrame.getWidth(), depthFrame.getHeight(), depthFrame.getDataSize(), (void*)depthFrame.getData());

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

    integrateTSDFVolume(depth, intr, tsdfVolume->getSize(), initialCamRotInv, initialCamTrans, tsdfVolume->getTsdfTruncDist(), tsdfVolume->data(), cyclical->getBuffer(), depthRawScaled);

		for (int i = 0; i < LEVELS; ++i)
			transformMaps(vMapsCur[i], nMapsCur[i], initialCamRot, initialCamTrans, vMapsGPrev[i], nMapsGPrev[i]);

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

	preprocessDepth(depth0, depth1, validMask0, validMask1, minDepth, maxDepth);

	vector<cv::Mat> pyramidImage0, pyramidDepth0, pyramidImage1, pyramidDepth1, pyramid_dI_dx1, pyramid_dI_dy1, pyramidTexturedMask1,	pyramidCameraMatrix;
    
	buildPyramids(image0, image1, depth0, depth1, cameraMatrix, sobelSize, sobelScale, *minGradientMagnitudesPtr,	pyramidImage0, pyramidDepth0, pyramidImage1, pyramidDepth1, pyramid_dI_dx1, pyramid_dI_dy1, pyramidTexturedMask1, pyramidCameraMatrix);

  cv::Mat resultRt = cv::Mat::eye(4,4,CV_64FC1);
	cv::Mat currRt, ksi;

	Matrix3frm camRotGlobalPrev = rotationMats[globalTime - 1];   // [Ri|ti] - pos of camera, i.e.

  // Previous global translation
	Vector3f camTransGlobalPrev = translationVecs[globalTime - 1];   // transform from camera to global coo space for (i-1)th camera pose

	if (frameTransform != NULL && (frameTransform->type == frameTransform->InitializeOnly)) {
		Eigen::Affine3f aff_rgbd(frameTransform->transformation);
		camRotGlobalPrev = aff_rgbd.linear();
		camTransGlobalPrev = aff_rgbd.translation();
	}
  else if (frameTransform != NULL && (frameTransform->type == frameTransform->IncrementalOnly)) {
		Eigen::Affine3f aff_rgbd(frameTransform->transformation * getCameraPose().matrix());
		camRotGlobalPrev = aff_rgbd.linear();
		camTransGlobalPrev = aff_rgbd.translation();
	}


	// Previous global inverse rotation
	Matrix3frm camRotGlobalPrevInv = camRotGlobalPrev.inverse ();  // Rprev.t();

	// GET CURRENT GLOBAL TRANSFORM
	Matrix3frm camRotGlobalCur = camRotGlobalPrev;                 // transform to global coo for ith camera pose
	Vector3f camTransGlobalCur = camTransGlobalPrev;

  Vector3f camTransLocalPrev = camTransGlobalPrev - cyclical->getBuffer()->origin_metric;



	//LOCAL PREVIOUS TRANSFORM
//	Mat33&  device_cam_rot_local_prev_inv = device_cast<Mat33> (cam_rot_global_prev_inv);
//	Mat33&  device_cam_rot_local_prev = device_cast<Mat33> (cam_rot_global_prev); 
//
//	float3& device_cam_trans_local_prev_tmp = device_cast<float3> (cam_trans_global_prev);
//	float3 device_cam_trans_local_prev;
//	device_cam_trans_local_prev.x = device_cam_trans_local_prev_tmp.x - (getCyclicalBufferStructure ())->origin_metric.x;
//	device_cam_trans_local_prev.y = device_cam_trans_local_prev_tmp.y - (getCyclicalBufferStructure ())->origin_metric.y;
//	device_cam_trans_local_prev.z = device_cam_trans_local_prev_tmp.z - (getCyclicalBufferStructure ())->origin_metric.z;
//	float3 device_volume_size = device_cast<const float3> (tsdf_volume_->getSize());
//
//  raycast (intr, device_cam_rot_local_prev, device_cam_trans_local_prev, tsdf_volume_->getTsdfTruncDist (), device_volume_size, tsdf_volume_->data (), getCyclicalBufferStructure (), vmaps_g_prev_[0], nmaps_g_prev_[0]);


	++globalTime;
	return true;
}























