//
//

#include "PortedGPUFunctions.h"


const float qnan = std::numeric_limits<float>::quiet_NaN();


void bilateralKernel (const cv::Mat src, cv::Mat & dst, float sigma_space2_inv_half, float sigma_color2_inv_half)
{
  const int R = 6;       //static_cast<int>(sigma_space * 1.5);
  const int D = R * 2 + 1;
  const int maxInt = numeric_limits<int>::max();

  for (int x = 0; x < src.cols; ++x) {
    for (int y = 0; y < src.rows; ++y) {

      int value = src.ptr(y)[x];

      int tx = min (x - D / 2 + D, src.cols - 1);
      int ty = min (y - D / 2 + D, src.rows - 1);

      float sum1 = 0;
      float sum2 = 0;

      for (int cy = max (y - D / 2, 0); cy < ty; ++cy)
      {
        for (int cx = max (x - D / 2, 0); cx < tx; ++cx)
        {
          int tmp = src.ptr(cy)[cx];

          float space2 = (x - cx) * (x - cx) + (y - cy) * (y - cy);
          float color2 = (value - tmp) * (value - tmp);

          float weight = exp(-(space2 * sigma_space2_inv_half + color2 * sigma_color2_inv_half));

          sum1 += tmp * weight;
          sum2 += weight;
        }
      }
      int res = round((sum1 / sum2) / 2) * 2;
      dst.ptr(y)[x] = max(0, min(res, maxInt));
    }
  }
}

void pyrDownKernel (const cv::Mat src, cv::Mat& dst, float sigma_color)
{
  const int D = 5;

  for (int x = 0; x < src.cols; ++x) {
    for (int y = 0; y < src.rows; ++y) {

      int center = src.ptr (2 * y)[2 * x];

      int tx = min(2 * x - D / 2 + D, src.cols - 1);
      int ty = min(2 * y - D / 2 + D, src.rows - 1);

      int sum = 0;
      int count = 0;

      for (int cy = max(0, 2 * y - D / 2); cy < ty; ++cy) {
        for (int cx = max(0, 2 * x - D / 2); cx < tx; ++cx) {
          int val = src.ptr(cy)[cx];
          if (abs(val - center) < 3 * sigma_color) {
            sum += val;
            ++count;
          }
        }
      }
      dst.ptr(y)[x] = sum / count;
    }
  }
}

void truncateDepthKernel(cv::Mat depth, int maxDistance)
{
  for (int x = 0; x < depth.cols; ++x)
    for (int y = 0; y < depth.rows; ++y)
			if(depth.ptr(y)[x] > maxDistance * 1000.f)
				depth.ptr(y)[x] = 0;
}



void createVMap (const Intr& intr, const cv::Mat& depth, cv::Mat& vMap)
{
  vMap.create(depth.rows * 3, depth.cols, CV_32FC1);

  for (int x = 0; x < depth.cols; ++x) {
    for (int y = 0; y < depth.rows; ++y) {

      float z = depth.ptr (y)[x] / 1000.f; // load and convert: mm -> meters

      if (z != 0)
      {
        float vx = z * (x - intr.cx) / intr.fx;
        float vy = z * (y - intr.cy) / intr.fy;
        float vz = z;

        vMap.ptr (y                 )[x] = vx;
        vMap.ptr (y + depth.rows    )[x] = vy;
        vMap.ptr (y + depth.rows * 2)[x] = vz;
      }
      else
        vMap.ptr (y)[x] = numeric_limits<float>::quiet_NaN();
    }
  }

}


void computeNormalsEigen(const cv::Mat& vMap, cv::Mat& nMap)
{
  int localRadius = 1;

  int rows = vMap.rows / 3;
  int cols = vMap.cols;

  nMap = cv::Mat(vMap.rows, vMap.cols, vMap.type());

  for (int x = 0; x < rows; ++x) {
    for (int y = 0; y < cols; ++y) {

      nMap.ptr(y)[x] = qnan;

      if (isnan(vMap.ptr(y)[x]))
        return;

      int ty = std::min(y + localRadius, rows - 1);
      int tx = std::min(x + localRadius, cols - 1);

      float3 centroid(0,0,0);
      int counter = 0;
      for (int cy = std::max(y - localRadius, 0); cy < ty; ++cy) {
        for (int cx = std::max(x - localRadius, 0); cx < tx; ++cx) {
          float v_x = vMap.ptr(cy)[cx];
          if (!isnan (v_x))
          {
            centroid.x += v_x;
            centroid.y += vMap.ptr(cy + rows)[cx];
            centroid.z += vMap.ptr(cy + 2 * rows)[cx];
            ++counter;
          }
        }
      }

      if (counter < localRadius * localRadius)
        return;

      centroid *= (1.f / counter);

      float cov[] = {0, 0, 0, 0, 0, 0};

      for (int cy = max(y - localRadius, 0); cy < ty; ++cy) {
        for (int cx = max(x - localRadius, 0); cx < tx; ++cx) {
          float3 v;
          v.x = vMap.ptr(cy)[cx];
          if (isnan (v.x))
            continue;

          v.y = vMap.ptr(cy + rows)[cx];
          v.z = vMap.ptr(cy + 2 * rows)[cx];

          float3 d = v - centroid;

          cov[0] += d.x * d.x;               //cov (0, 0)
          cov[1] += d.x * d.y;               //cov (0, 1)
          cov[2] += d.x * d.z;               //cov (0, 2)
          cov[3] += d.y * d.y;               //cov (1, 1)
          cov[4] += d.y * d.z;               //cov (1, 2)
          cov[5] += d.z * d.z;               //cov (2, 2)
        }
      }

      Eigen33 eigen33 (cov);
      typedef Eigen33::Mat33 Mat33;

      Mat33 tmp;
      Mat33 vec_tmp;
      Mat33 evecs;
      float3 evals;
      eigen33.compute (tmp, vec_tmp, evecs, evals);

      float3 n = normalized (evecs[0]);

      nMap.ptr (y           )[x] = n.x;
      nMap.ptr (y + rows    )[x] = n.y;
      nMap.ptr (y + rows * 2)[x] = n.z;
    }
  }
}


void transformMaps(const cv::Mat& vMapSrc, const cv::Mat& nMapSrc, const Eigen33::Mat33 rMat, const float3 tVec, cv::Mat vMapDst, cv::Mat nMapDst)
{
  int cols = vMapSrc.cols;
  int rows = vMapSrc.rows / 3;

  vMapDst = cv::Mat(rows * 3, cols, vMapSrc.type());
  nMapDst = cv::Mat(rows * 3, cols, nMapSrc.type());

  for (int x = 0; x < cols; ++x) {
    for (int y =0; y < rows; ++y) {

      //vetexes
      float3 vSrc, vDst = float3(qnan, qnan, qnan);
      vSrc.x = vMapSrc.ptr(y)[x];

      if (!isnan(vSrc.x))
      {
        vSrc.y = vMapSrc.ptr(y + rows)[x];
        vSrc.z = vMapSrc.ptr(y + 2 * rows)[x];

        vDst = rMat * vSrc + tVec;

        vMapDst.ptr(y + rows)[x] = vDst.y;
        vMapDst.ptr(y + 2 * rows)[x] = vDst.z;
      }

      vMapDst.ptr (y)[x] = vDst.x;

      //normals
      float3 nSrc, nDst = float3(qnan, qnan, qnan);
      nSrc.x = nMapSrc.ptr(y)[x];

      if (!isnan(nSrc.x))
      {
        nSrc.y = nMapSrc.ptr(y + rows)[x];
        nSrc.z = nMapSrc.ptr(y + 2 * rows)[x];

        nDst = rMat * nSrc;

        nMapDst.ptr (y + rows)[x] = nDst.y;
        nMapDst.ptr (y + 2 * rows)[x] = nDst.z;
      }

      nMapDst.ptr (y)[x] = nDst.x;

    }
  }
}


/*
__global__ void
tsdf23 (const PtrStepSz<float> depthScaled, PtrStep<short2> volume,
        const float tranc_dist, const Mat33 Rcurr_inv, const float3 tcurr, const Intr intr, const float3 cell_size, const pcl::gpu::tsdf_buffer buffer)
{
  int x = threadIdx.x + blockIdx.x * blockDim.x;
  int y = threadIdx.y + blockIdx.y * blockDim.y;

  if (x >= buffer.voxels_size.x || y >= buffer.voxels_size.y)
    return;

  float v_g_x = (x + 0.5f) * cell_size.x - tcurr.x;
  float v_g_y = (y + 0.5f) * cell_size.y - tcurr.y;
  float v_g_z = (0 + 0.5f) * cell_size.z - tcurr.z;

  float v_g_part_norm = v_g_x * v_g_x + v_g_y * v_g_y;

  float v_x = (Rcurr_inv.data[0].x * v_g_x + Rcurr_inv.data[0].y * v_g_y + Rcurr_inv.data[0].z * v_g_z) * intr.fx;
  float v_y = (Rcurr_inv.data[1].x * v_g_x + Rcurr_inv.data[1].y * v_g_y + Rcurr_inv.data[1].z * v_g_z) * intr.fy;
  float v_z = (Rcurr_inv.data[2].x * v_g_x + Rcurr_inv.data[2].y * v_g_y + Rcurr_inv.data[2].z * v_g_z);

  float z_scaled = 0;

  float Rcurr_inv_0_z_scaled = Rcurr_inv.data[0].z * cell_size.z * intr.fx;
  float Rcurr_inv_1_z_scaled = Rcurr_inv.data[1].z * cell_size.z * intr.fy;

  float tranc_dist_inv = 1.0f / tranc_dist;

  short2* pos = volume.ptr (y) + x;
  
  // shift the pointer to relative indices
  shift_tsdf_pointer(&pos, buffer);
  
  int elem_step = volume.step * buffer.voxels_size.y / sizeof(short2);

//#pragma unroll
  for (int z = 0; z < buffer.voxels_size.z;
       ++z,
       v_g_z += cell_size.z,
       z_scaled += cell_size.z,
       v_x += Rcurr_inv_0_z_scaled,
       v_y += Rcurr_inv_1_z_scaled,
       pos += elem_step)
  {
    
    // As the pointer is incremented in the for loop, we have to make sure that the pointer is never outside the memory
    if(pos > buffer.tsdf_memory_end)
      pos -= (buffer.tsdf_memory_end - buffer.tsdf_memory_start + 1);
    
    float inv_z = 1.0f / (v_z + Rcurr_inv.data[2].z * z_scaled);
    if (inv_z < 0)
        continue;

    // project to current cam
// old code
    int2 coo =
    {
      __float2int_rn (v_x * inv_z + intr.cx),
      __float2int_rn (v_y * inv_z + intr.cy)
    };

    if (coo.x >= 0 && coo.y >= 0 && coo.x < depthScaled.cols && coo.y < depthScaled.rows)         //6
    {
      float Dp_scaled = depthScaled.ptr (coo.y)[coo.x]; //meters

      float sdf = Dp_scaled - sqrtf (v_g_z * v_g_z + v_g_part_norm);

      if (Dp_scaled != 0 && sdf >= -tranc_dist) //meters
      {
        float tsdf = fmin (1.0f, sdf * tranc_dist_inv);

        //read and unpack
        float tsdf_prev;
        int weight_prev;
        unpack_tsdf (*pos, tsdf_prev, weight_prev);

        const int Wrk = 1;

        float tsdf_new = (tsdf_prev * weight_prev + Wrk * tsdf) / (weight_prev + Wrk);
        int weight_new = min (weight_prev + Wrk, Tsdf::MAX_WEIGHT);

        pack_tsdf (tsdf_new, weight_new, *pos);
      }
    }
  }
}


*/













