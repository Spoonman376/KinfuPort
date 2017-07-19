//
//

#ifndef PortedGPUFunctions_h
#define PortedGPUFunctions_h

#include <stdio.h>
#include <limits>

#include <opencv2/opencv.hpp>

#include "Eigen.h"


//#include "Utils.h"


using namespace std;

struct Intr;


void bilateralKernel(const cv::Mat src, cv::Mat& dst, float sigma_space2_inv_half, float sigma_color2_inv_half);

void pyrDownKernel(const cv::Mat src, cv::Mat& dst, float sigma_color);

void truncateDepthKernel(cv::Mat depth, int max_distance_mm);

void createVMap(const Intr& intr, const cv::Mat& depth, cv::Mat& vMap);

void computeNormalsEigen(const cv::Mat& vMap, cv::Mat& nMap);
 
void transformMaps(const cv::Mat& vMapSrc, const cv::Mat& nMapSrc, const Matrix3frm rMat, const Vector3f tVec, cv::Mat vMapDst, cv::Mat nMapDst);

struct Intr
{
  float fx, fy, cx, cy, trunc_dist;
  Intr () {};
  Intr (float fx_, float fy_, float cx_, float cy_, float trunc_dist_)
    : fx(fx_), fy(fy_), cx(cx_), cy(cy_), trunc_dist(trunc_dist_) {};

  Intr operator()(int level_index) const
  { 
    int div = 1 << level_index; 
    return (Intr (fx / div, fy / div, cx / div, cy / div, trunc_dist));
  }
};


#endif /* PortedGPUFunctions_h */
