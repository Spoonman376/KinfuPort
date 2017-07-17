//
//

#ifndef PortedGPUFunctions_h
#define PortedGPUFunctions_h

#include <stdio.h>
#include <limits>

#include <opencv2/opencv.hpp>

#include <Eigen/Core>


using namespace std;

struct Intr;

void bilateralKernel(const cv::Mat src, cv::Mat& dst, float sigma_space2_inv_half, float sigma_color2_inv_half);

void pyrDownKernel(const cv::Mat src, cv::Mat& dst, float sigma_color);

void truncateDepthKernel(cv::Mat depth, int max_distance_mm);

void createVMap(const Intr& intr, const cv::Mat& depth, cv::Mat& vMap);

void computeNormalsEigen(const cv::Mat& vMap, cv::Mat& nMap);


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

//struct float3
//{
//  float x, y, z;
//
//  float3()
//  {
//    x = y = z = 0.0;
//  }
//
//  float3(float a, float b, float c)
//  {
//    x = a;
//    y = b;
//    z = c;
//  }
//
//  float3 operator -(float3 s)
//  {
//    return float3(x - s.x, y - s.y, z - s.z);
//  }
//
//  float3 operator +(float3 s)
//  {
//    return float3(x + s.x, y + s.y, z + s.z);
//  }
//
//  void operator *=(float s)
//  {
//    x *= s;
//    y *= s;
//    z *= s;
//  }
//
//
//};



#endif /* PortedGPUFunctions_h */
