//
//

#ifndef TSDFVolume_h
#define TSDFVolume_h

#include <stdio.h>
#include <vector>

#include <pcl/point_types.h>
#include <pcl/point_cloud.h>

#include <Eigen/Core>

#include <opencv2/opencv.hpp>

using namespace Eigen;


class TSDFVolume
{

protected:
  /** \brief tsdf volume size in meters */
  Eigen::Vector3f size;
  
  /** \brief tsdf volume resolution */
  Eigen::Vector3i resolution;

  /** \brief tsdf volume data container */
  cv::Mat volume;

  /** \brief tsdf truncation distance */
  float truncDist;

public:

  TSDFVolume(const Eigen::Vector3i& resolution);


  /** \brief Resets tsdf volume data to uninitialized state */
  void reset();

  /** \brief Sets Tsdf volume size for each dimention
    * \param[in] size size of tsdf volume in meters
    */
  void setSize(const Eigen::Vector3f& size);
  
  /** \brief Sets Tsdf truncation distance. Must be greater than 2 * volume_voxel_size
    * \param[in] distance TSDF truncation distance 
    */
  void setTsdfTruncDist (float distance);


};


#endif /* TSDFVolume_h */
