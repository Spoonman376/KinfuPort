//
//

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
