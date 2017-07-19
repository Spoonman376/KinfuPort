//
//

#include "TSDFVolume.h"

const int DIVISOR = 32767;     // SHRT_MAX;

TSDFVolume::TSDFVolume(const Vector3i& r)
{
  resolution = r;

  volume.create(resolution.y() * resolution.z(), resolution.x(), sizeof(int));

  const Vector3f defaultVolumeSize = Vector3f::Constant (3.f); //meters
  const float defaultTruncDist = 0.03f; //meters

  setSize(defaultVolumeSize);
  setTsdfTruncDist(defaultTruncDist);

  reset();
}

cv::Mat TSDFVolume::data()
{
  return volume;
}

void TSDFVolume::setSize(const Vector3f& vs)
{
  size = vs;
  setTsdfTruncDist(truncDist);
}

Vector3f TSDFVolume::getSize()
{
  return size;
}


void TSDFVolume::setTsdfTruncDist(float distance)
{
  float cx = size(0) / resolution(0);
  float cy = size(1) / resolution(1);
  float cz = size(2) / resolution(2);

  truncDist = std::max (distance, 2.1f * std::max (cx, std::max (cy, cz)));
}

float TSDFVolume::getTsdfTruncDist()
{
  return truncDist;
}


void TSDFVolume::reset()
{

}





/** Other Functions
  *
  */

void integrateTSDFVolume(const cv::Mat& depth, const Intr& intr, const Vector3f& volume_size,
                         const Matrix3frm &Rcurr_inv, const Vector3f& tcurr,
                         float tranc_dist, cv::Mat volume, const tsdf_buffer* buffer,
                         cv::Mat& depthScaled
                         )
{
  depthScaled.create(depth.rows, depth.cols, depth.type());

  //scales depth along ray and converts mm -> meters. 
  for (int x = 0; x < depth.cols; ++x) {
    for (int y = 0; y < depth.rows; ++y) {
      int Dp = depth.ptr(y)[x];

      float xl = (x - intr.cx) / intr.fx;
      float yl = (y - intr.cy) / intr.fy;
      float lambda = sqrtf (xl * xl + yl * yl + 1);

      float res = Dp * lambda/1000.f; //meters
      if ( intr.trunc_dist > 0 && res > intr.trunc_dist )
        depthScaled.ptr (y)[x] = 0;
      else
      depthScaled.ptr(y)[x] = res;
    }
  }

  Vector3f cell_size;
  cell_size.x() = volume_size.x() / buffer->voxels_size.x();
  cell_size.y() = volume_size.y() / buffer->voxels_size.y();
  cell_size.z() = volume_size.z() / buffer->voxels_size.z();

  Matrix3frm mat;

  for (int x = 0; x < buffer->voxels_size.x(); ++x) {
    for (int y = 0; y < buffer->voxels_size.y(); ++y) {

      float v_g_x = (x + 0.5f) * cell_size.x() - tcurr.x();
      float v_g_y = (y + 0.5f) * cell_size.y() - tcurr.y();
      float v_g_z = (0 + 0.5f) * cell_size.z() - tcurr.z();

      float v_g_part_norm = v_g_x * v_g_x + v_g_y * v_g_y;

      // these indexes might be wrong based on if Eigen is row / column major
      float v_x = (Rcurr_inv(0,0) * v_g_x + Rcurr_inv(0,1) * v_g_y + Rcurr_inv(0,2) * v_g_z) * intr.fx;
      float v_y = (Rcurr_inv(1,0) * v_g_x + Rcurr_inv(1,1) * v_g_y + Rcurr_inv(1,2) * v_g_z) * intr.fy;
      float v_z = (Rcurr_inv(2,0) * v_g_x + Rcurr_inv(2,1) * v_g_y + Rcurr_inv(2,2) * v_g_z);

      float z_scaled = 0;

      float Rcurr_inv_0_z_scaled = Rcurr_inv(0,2) * cell_size.z() * intr.fx;
      float Rcurr_inv_1_z_scaled = Rcurr_inv(1,2) * cell_size.z() * intr.fy;

      float tranc_dist_inv = 1.0f / tranc_dist;

      short2* pos = (short2*)volume.ptr(y) + x;

      // shift the pointer to relative indices
      shiftTsdfPointer(&pos, buffer);

      int elem_step = (int)(volume.step * buffer->voxels_size.y() / sizeof(short2));

//#pragma unroll
      for (int z = 0; z < buffer->voxels_size.z();
           ++z,
           v_g_z += cell_size.z(),
           z_scaled += cell_size.z(),
           v_x += Rcurr_inv_0_z_scaled,
           v_y += Rcurr_inv_1_z_scaled,
           pos += elem_step)
      {
        
        // As the pointer is incremented in the for loop, we have to make sure that the pointer is never outside the memory
        if(pos > buffer->tsdf_memory_end)
          pos -= (buffer->tsdf_memory_end - buffer->tsdf_memory_start + 1);
        
        float inv_z = 1.0f / (v_z + Rcurr_inv(2,2) * z_scaled);
        if (inv_z < 0)
            continue;

        // project to current cam
		// old code
        Eigen::Vector2i coo = Eigen::Vector2i(std::round((v_x * inv_z + intr.cx) / 2) * 2,
                                              std::round((v_y * inv_z + intr.cy) / 2) * 2);

        if (coo.x() >= 0 && coo.y() >= 0 && coo.x() < depthScaled.cols && coo.y() < depthScaled.rows) //6
        {
          float Dp_scaled = depthScaled.ptr (coo.y())[coo.x()]; //meters

          float sdf = Dp_scaled - sqrtf (v_g_z * v_g_z + v_g_part_norm);

          if (Dp_scaled != 0 && sdf >= -tranc_dist) //meters
          {
            float tsdf = fmin (1.0f, sdf * tranc_dist_inv);

            //read and unpack
            float tsdf_prev;
            int weight_prev;
            unpackTsdf (*pos, tsdf_prev, weight_prev);

            const int Wrk = 1;

            float tsdf_new = (tsdf_prev * weight_prev + Wrk * tsdf) / (weight_prev + Wrk);
            int weight_new = std::min(weight_prev + Wrk, 1 << 7);

            packTsdf(tsdf_new, weight_new, *pos);
          }
        }
      }       // for(int z = 0; z < VOLUME_Z; ++z)
    }
  }
}


void shiftTsdfPointer(short2** value, const tsdf_buffer* buffer)
{
  ///Shift the pointer by (@origin - @start)
  *value += (buffer->tsdf_rolling_buff_origin - buffer->tsdf_memory_start);
  
  ///If we land outside of the memory, make sure to "modulo" the new value
  if(*value > buffer->tsdf_memory_end)
  {
    *value -= (buffer->tsdf_memory_end - buffer->tsdf_memory_start + 1);
  }       
}

void unpackTsdf (short2 value, float& tsdf, int& weight)
{
  weight = value.y;
  tsdf = (float)(value.x / 2) * 2 / DIVISOR;
}

void packTsdf (float tsdf, int weight, short2& value)
{
  int fixedp = std::max(-DIVISOR, std::min(DIVISOR, (int)(tsdf * DIVISOR)));
  value = short2(fixedp, weight);
}






