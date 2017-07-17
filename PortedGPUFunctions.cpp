//
//

#include "PortedGPUFunctions.h"
#include "Utils.h"



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

  for (int x = 0; x < rows; ++x) {
    for (int y = 0; y < cols; ++y) {

      nMap.ptr(y)[x] = numeric_limits<float>::quiet_NaN();

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

      typedef Eigen33::Mat33 Mat33;
      Eigen33 eigen33 (cov);

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


/*

__global__ void
computeNmapKernelEigen (int rows, int cols, const PtrStep<float> vmap, PtrStep<float> nmap)
{
  int u = threadIdx.x + blockIdx.x * blockDim.x;
  int v = threadIdx.y + blockIdx.y * blockDim.y;

  if (u >= cols || v >= rows)
    return;

  nmap.ptr (v)[u] = numeric_limits<float>::quiet_NaN ();

  if (isnan (vmap.ptr (v)[u]))
    return;

  int ty = min (v - ky / 2 + ky, rows - 1);
  int tx = min (u - kx / 2 + kx, cols - 1);

  float3 centroid = make_float3 (0.f, 0.f, 0.f);
  int counter = 0;
  for (int cy = max (v - ky / 2, 0); cy < ty; cy += STEP)
    for (int cx = max (u - kx / 2, 0); cx < tx; cx += STEP)
    {
      float v_x = vmap.ptr (cy)[cx];
      if (!isnan (v_x))
      {
        centroid.x += v_x;
        centroid.y += vmap.ptr (cy + rows)[cx];
        centroid.z += vmap.ptr (cy + 2 * rows)[cx];
        ++counter;
      }
    }

  if (counter < kx * ky / 2)
    return;

  centroid *= 1.f / counter;

  float cov[] = {0, 0, 0, 0, 0, 0};

  for (int cy = max (v - ky / 2, 0); cy < ty; cy += STEP)
    for (int cx = max (u - kx / 2, 0); cx < tx; cx += STEP)
    {
      float3 v;
      v.x = vmap.ptr (cy)[cx];
      if (isnan (v.x))
        continue;

      v.y = vmap.ptr (cy + rows)[cx];
      v.z = vmap.ptr (cy + 2 * rows)[cx];

      float3 d = v - centroid;

      cov[0] += d.x * d.x;               //cov (0, 0)
      cov[1] += d.x * d.y;               //cov (0, 1)
      cov[2] += d.x * d.z;               //cov (0, 2)
      cov[3] += d.y * d.y;               //cov (1, 1)
      cov[4] += d.y * d.z;               //cov (1, 2)
      cov[5] += d.z * d.z;               //cov (2, 2)
    }

  typedef Eigen33::Mat33 Mat33;
  Eigen33 eigen33 (cov);

  Mat33 tmp;
  Mat33 vec_tmp;
  Mat33 evecs;
  float3 evals;
  eigen33.compute (tmp, vec_tmp, evecs, evals);

  float3 n = normalized (evecs[0]);

  u = threadIdx.x + blockIdx.x * blockDim.x;
  v = threadIdx.y + blockIdx.y * blockDim.y;

  nmap.ptr (v       )[u] = n.x;
  nmap.ptr (v + rows)[u] = n.y;
  nmap.ptr (v + 2 * rows)[u] = n.z;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
pcl::device::computeNormalsEigen (const MapArr& vmap, MapArr& nmap)
{
  int cols = vmap.cols ();
  int rows = vmap.rows () / 3;

  nmap.create (vmap.rows (), vmap.cols ());

  dim3 block (32, 8);
  dim3 grid (1, 1, 1);
  grid.x = divUp (cols, block.x);
  grid.y = divUp (rows, block.y);

  computeNmapKernelEigen<<<grid, block>>>(rows, cols, vmap, nmap);
  cudaSafeCall (cudaGetLastError ());
  cudaSafeCall (cudaDeviceSynchronize ());
}


*/













