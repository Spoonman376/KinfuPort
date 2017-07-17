//
//

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





	++globalTime;
	return true;
}
