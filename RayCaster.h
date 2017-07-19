//
//

#ifndef RayCaster_h
#define RayCaster_h

#include <stdio.h>

#include "Eigen.h"

#include <opencv2/opencv.hpp>

#include "PortedGPUFunctions.h"

struct RayCaster
{
  enum { CTA_SIZE_X = 32, CTA_SIZE_Y = 8 };

  Matrix3frm Rcurr;
  Vector3f tcurr;

  float time_step;
  Vector3f volume_size;

  Vector3f cell_size;
  int cols, rows;

  cv::Mat volume; // <short2>

  Intr intr;

  cv::Mat nmap; // <float>
  cv::Mat vmap; // <float>

  Vector3f getRayNext (int x, int y) const
  {
    Vector3f rayNext;
    rayNext.x() = (x - intr.cx) / intr.fx;
    rayNext.y() = (y - intr.cy) / intr.fy;
    rayNext.z() = 1;
    return rayNext;
  }

  bool checkInds (const Vector3i& g) const
  {
    return (g.x() >= 0 && g.y() >= 0 && g.z() >= 0 && g.x() < VOLUME_X && g.y() < VOLUME_Y && g.z() < VOLUME_Z);
  }

  __device__ __forceinline__ float
  readTsdf (int x, int y, int z, pcl::gpu::tsdf_buffer buffer) const
  {
    const short2* tmp_pos = &(volume.ptr (buffer.voxels_size.y * z + y)[x]);
    short2* pos = const_cast<short2*> (tmp_pos);
    shift_tsdf_pointer(&pos, buffer);
    return unpack_tsdf (*pos);
  }

  __device__ __forceinline__ int3
  getVoxel (float3 point) const
  {
    int vx = __float2int_rd (point.x / cell_size.x);        // round to negative infinity
    int vy = __float2int_rd (point.y / cell_size.y);
    int vz = __float2int_rd (point.z / cell_size.z);

    return make_int3 (vx, vy, vz);
  }

  __device__ __forceinline__ float
  interpolateTrilineary (const float3& origin, const float3& dir, float time, pcl::gpu::tsdf_buffer buffer) const
  {
    return interpolateTrilineary (origin + dir * time, buffer);
  }

  __device__ __forceinline__ float
  interpolateTrilineary (const float3& point, pcl::gpu::tsdf_buffer buffer) const
  {
    int3 g = getVoxel (point);

    if (g.x <= 0 || g.x >= buffer.voxels_size.x - 1)
      return numeric_limits<float>::quiet_NaN ();

    if (g.y <= 0 || g.y >= buffer.voxels_size.y - 1)
      return numeric_limits<float>::quiet_NaN ();

    if (g.z <= 0 || g.z >= buffer.voxels_size.z - 1)
      return numeric_limits<float>::quiet_NaN ();

/*      //OLD CODE
    float vx = (g.x + 0.5f) * cell_size.x;
    float vy = (g.y + 0.5f) * cell_size.y;
    float vz = (g.z + 0.5f) * cell_size.z;

    g.x = (point.x < vx) ? (g.x - 1) : g.x;
    g.y = (point.y < vy) ? (g.y - 1) : g.y;
    g.z = (point.z < vz) ? (g.z - 1) : g.z;

    float a = (point.x - (g.x + 0.5f) * cell_size.x) / cell_size.x;
    float b = (point.y - (g.y + 0.5f) * cell_size.y) / cell_size.y;
    float c = (point.z - (g.z + 0.5f) * cell_size.z) / cell_size.z;

    float res = readTsdf (g.x + 0, g.y + 0, g.z + 0, buffer) * (1 - a) * (1 - b) * (1 - c) +
                readTsdf (g.x + 0, g.y + 0, g.z + 1, buffer) * (1 - a) * (1 - b) * c +
                readTsdf (g.x + 0, g.y + 1, g.z + 0, buffer) * (1 - a) * b * (1 - c) +
                readTsdf (g.x + 0, g.y + 1, g.z + 1, buffer) * (1 - a) * b * c +
                readTsdf (g.x + 1, g.y + 0, g.z + 0, buffer) * a * (1 - b) * (1 - c) +
                readTsdf (g.x + 1, g.y + 0, g.z + 1, buffer) * a * (1 - b) * c +
                readTsdf (g.x + 1, g.y + 1, g.z + 0, buffer) * a * b * (1 - c) +
                readTsdf (g.x + 1, g.y + 1, g.z + 1, buffer) * a * b * c;
*/
    //NEW CODE
float a = point.x/ cell_size.x - (g.x + 0.5f); if (a<0) { g.x--; a+=1.0f; };
    float b = point.y/ cell_size.y - (g.y + 0.5f); if (b<0) { g.y--; b+=1.0f; };
    float c = point.z/ cell_size.z - (g.z + 0.5f); if (c<0) { g.z--; c+=1.0f; };

    float res = (1 - a) * (
        (1 - b) * (
          readTsdf (g.x + 0, g.y + 0, g.z + 0, buffer) * (1 - c) +
          readTsdf (g.x + 0, g.y + 0, g.z + 1, buffer) *      c 
          )
        + b * (
          readTsdf (g.x + 0, g.y + 1, g.z + 0, buffer) * (1 - c) +
          readTsdf (g.x + 0, g.y + 1, g.z + 1, buffer) *      c  
          )
        )
      + a * (
        (1 - b) * (
          readTsdf (g.x + 1, g.y + 0, g.z + 0, buffer) * (1 - c) +
          readTsdf (g.x + 1, g.y + 0, g.z + 1, buffer) *      c 
          )
        + b * (
          readTsdf (g.x + 1, g.y + 1, g.z + 0, buffer) * (1 - c) +
          readTsdf (g.x + 1, g.y + 1, g.z + 1, buffer) *      c 
          )
        )
      ;
    return res;
  }


  __device__ __forceinline__ void
  operator () (pcl::gpu::tsdf_buffer buffer) const
  {
    int x = threadIdx.x + blockIdx.x * CTA_SIZE_X;
    int y = threadIdx.y + blockIdx.y * CTA_SIZE_Y;

    if (x >= cols || y >= rows)
      return;

    vmap.ptr (y)[x] = numeric_limits<float>::quiet_NaN ();
    nmap.ptr (y)[x] = numeric_limits<float>::quiet_NaN ();

    float3 ray_start = tcurr;
    float3 ray_next = Rcurr * get_ray_next (x, y) + tcurr;

    float3 ray_dir = normalized (ray_next - ray_start);

    //ensure that it isn't a degenerate case
    ray_dir.x = (ray_dir.x == 0.f) ? 1e-15 : ray_dir.x;
    ray_dir.y = (ray_dir.y == 0.f) ? 1e-15 : ray_dir.y;
    ray_dir.z = (ray_dir.z == 0.f) ? 1e-15 : ray_dir.z;

    // computer time when entry and exit volume
    float time_start_volume = getMinTime (volume_size, ray_start, ray_dir);
    float time_exit_volume = getMaxTime (volume_size, ray_start, ray_dir);

    const float min_dist = 0.f;         //in meters
    time_start_volume = fmax (time_start_volume, min_dist);
    if (time_start_volume >= time_exit_volume)
      return;

    float time_curr = time_start_volume;
    int3 g = getVoxel (ray_start + ray_dir * time_curr);
    g.x = max (0, min (g.x, buffer.voxels_size.x - 1));
    g.y = max (0, min (g.y, buffer.voxels_size.y - 1));
    g.z = max (0, min (g.z, buffer.voxels_size.z - 1));

    float tsdf = readTsdf (g.x, g.y, g.z, buffer);

    //infinite loop guard
    const float max_time = 3 * (volume_size.x + volume_size.y + volume_size.z);

    for (; time_curr < max_time; time_curr += time_step)
    {
      float tsdf_prev = tsdf;

      int3 g = getVoxel (  ray_start + ray_dir * (time_curr + time_step)  );
      if (!checkInds (g))
        break;

      tsdf = readTsdf (g.x, g.y, g.z, buffer);

      if (tsdf_prev < 0.f && tsdf >= 0.f)
        break;

      if (tsdf_prev >= 0.f && tsdf < 0.f)           //zero crossing
      {
        float Ftdt = interpolateTrilineary (ray_start, ray_dir, time_curr + time_step, buffer);
        if (isnan (Ftdt))
          break;

        float Ft = interpolateTrilineary (ray_start, ray_dir, time_curr, buffer);
        if (isnan (Ft))
          break;

        //float Ts = time_curr - time_step * Ft/(Ftdt - Ft);
        float Ts = time_curr - time_step * Ft / (Ftdt - Ft);

        float3 vetex_found = ray_start + ray_dir * Ts;

        vmap.ptr (y       )[x] = vetex_found.x;
        vmap.ptr (y + rows)[x] = vetex_found.y;
        vmap.ptr (y + 2 * rows)[x] = vetex_found.z;

        int3 g = getVoxel ( ray_start + ray_dir * time_curr );
        if (g.x > 1 && g.y > 1 && g.z > 1 && g.x < buffer.voxels_size.x - 2 && g.y < buffer.voxels_size.y - 2 && g.z < buffer.voxels_size.z - 2)
        {
          float3 t;
          float3 n;

          t = vetex_found;
          t.x += cell_size.x;
          float Fx1 = interpolateTrilineary (t, buffer);

          t = vetex_found;
          t.x -= cell_size.x;
          float Fx2 = interpolateTrilineary (t, buffer);

          n.x = (Fx1 - Fx2);

          t = vetex_found;
          t.y += cell_size.y;
          float Fy1 = interpolateTrilineary (t, buffer);

          t = vetex_found;
          t.y -= cell_size.y;
          float Fy2 = interpolateTrilineary (t, buffer);

          n.y = (Fy1 - Fy2);

          t = vetex_found;
          t.z += cell_size.z;
          float Fz1 = interpolateTrilineary (t, buffer);

          t = vetex_found;
          t.z -= cell_size.z;
          float Fz2 = interpolateTrilineary (t, buffer);

          n.z = (Fz1 - Fz2);

          n = normalized (n);

          nmap.ptr (y       )[x] = n.x;
          nmap.ptr (y + rows)[x] = n.y;
          nmap.ptr (y + 2 * rows)[x] = n.z;
        }
        break;
      }

    }          /* for(;;)  */
  }
};




#endif /* RayCaster_h */
