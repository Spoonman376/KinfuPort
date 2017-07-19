//
//

#ifndef TSDFVolume_h
#define TSDFVolume_h

#include <stdio.h>
#include <vector>

#include <pcl/point_types.h>
#include <pcl/point_cloud.h>

#include "Eigen.h"

#include <opencv2/opencv.hpp>

#include "TSDFBuffer.h"
#include "PortedGPUFunctions.h"

using namespace Eigen;


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

  Eigen::Vector3f getSize();
  
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

void packTsdf(float tsdf, int weight, short2& value);


#endif /* TSDFVolume_h */
