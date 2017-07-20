/*
 * Software License Agreement (BSD License)
 *
 *  Point Cloud Library (PCL) - www.pointclouds.org
 *  Copyright (c) 2011, Willow Garage, Inc.
 *
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 *   * Neither the name of Willow Garage, Inc. nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 *  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 *  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 *  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 */

#include "RayCaster.h"

Vector3f RayCaster::getRayNext(int x, int y) const
{
  Vector3f rayNext;
  rayNext.x() = (x - intr.cx) / intr.fx;
  rayNext.y() = (y - intr.cy) / intr.fy;
  rayNext.z() = 1;
  return rayNext;
}

bool RayCaster::checkInds(const Vector3i& g) const
{
  return (g.x() >= 0 && g.y() >= 0 && g.z() >= 0 && g.x() < VOLUME_X && g.y() < VOLUME_Y && g.z() < VOLUME_Z);
}

float RayCaster::readTsdf(int x, int y, int z, const tsdf_buffer buffer) const
{
  const short2* tmp_pos = (short2*)&(volume.ptr(buffer.voxels_size.y() * z + y)[x]);
  short2* pos = const_cast<short2*> (tmp_pos);
  shiftTsdfPointer(&pos, &buffer);
  return unpackTsdf(*pos);
}

Vector3i RayCaster::getVoxel(Vector3f point) const
{
  int vx = std::floor(point.x() / cell_size.x());        // round to negative infinity
  int vy = std::floor(point.y() / cell_size.y());
  int vz = std::floor(point.z() / cell_size.z());

  return Vector3i(vx, vy, vz);
}

float RayCaster::interpolateTrilineary(const Vector3f& origin, const Vector3f& dir, float time, const tsdf_buffer buffer) const
{
  return interpolateTrilineary (origin + dir * time, buffer);
}

float RayCaster::interpolateTrilineary(const Vector3f& point, const tsdf_buffer buffer) const
{
  Vector3i g = getVoxel(point);

  if (g.x() <= 0 || g.x() >= buffer.voxels_size.x() - 1)
    return std::numeric_limits<float>::quiet_NaN();

  if (g.y() <= 0 || g.y() >= buffer.voxels_size.y() - 1)
    return std::numeric_limits<float>::quiet_NaN();

  if (g.z() <= 0 || g.z() >= buffer.voxels_size.z() - 1)
    return std::numeric_limits<float>::quiet_NaN();

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
  float a = point.x() / cell_size.x() - (g.x() + 0.5f);
  if (a<0)
    {g.x()--; a += 1.0f;}

  float b = point.y() / cell_size.y() - (g.y() + 0.5f);
  if (b<0)
    {g.y()--; b += 1.0f;};

  float c = point.z() / cell_size.z() - (g.z ()+ 0.5f);
  if (c<0)
    {g.z()--; c += 1.0f;};

  return (1 - a) * ((1 - b) * (readTsdf (g.x() + 0, g.y() + 0, g.z() + 0, buffer) * (1 - c) +
                               readTsdf (g.x() + 0, g.y() + 0, g.z() + 1, buffer) *      c)
                       + b  * (readTsdf (g.x() + 0, g.y() + 1, g.z() + 0, buffer) * (1 - c) +
                               readTsdf (g.x() + 0, g.y() + 1, g.z() + 1, buffer) *      c)
                    )
            + a  * ((1 - b) * (readTsdf (g.x() + 1, g.y() + 0, g.z() + 0, buffer) * (1 - c) +
                               readTsdf (g.x() + 1, g.y() + 0, g.z() + 1, buffer) *      c)
                       + b  * (readTsdf (g.x() + 1, g.y() + 1, g.z() + 0, buffer) * (1 - c) +
                               readTsdf (g.x() + 1, g.y() + 1, g.z() + 1, buffer) *      c)
                    );
}



void RayCaster::operator () (const tsdf_buffer buffer) const
{
  for (int x = 0; x < cols; ++x) {
    for (int y = 0; y < rows; ++y) {

      vmap->ptr(y)[x] = std::numeric_limits<float>::quiet_NaN();
      nmap->ptr(y)[x] = std::numeric_limits<float>::quiet_NaN();

      Vector3f ray_start = tcurr;
      Vector3f ray_next = Rcurr * getRayNext(x, y) + tcurr;

      Vector3f ray_dir = (ray_next - ray_start).normalized();

      //ensure that it isn't a degenerate case
      ray_dir.x() = (ray_dir.x() == 0.f) ? 1e-15 : ray_dir.x();
      ray_dir.y() = (ray_dir.y() == 0.f) ? 1e-15 : ray_dir.y();
      ray_dir.z() = (ray_dir.z() == 0.f) ? 1e-15 : ray_dir.z();

      // computer time when entry and exit volume
      float time_start_volume = getMinTime(volume_size, ray_start, ray_dir);
      float time_exit_volume = getMaxTime(volume_size, ray_start, ray_dir);

      const float min_dist = 0.f;         //in meters
      time_start_volume = std::fmax(time_start_volume, min_dist);
      if (time_start_volume >= time_exit_volume)
        return;

      float time_curr = time_start_volume;
      Vector3i g = getVoxel(ray_start + ray_dir * time_curr);
      g.x() = std::max (0, std::min (g.x(), buffer.voxels_size.x() - 1));
      g.y() = std::max (0, std::min (g.y(), buffer.voxels_size.y() - 1));
      g.z() = std::max (0, std::min (g.z(), buffer.voxels_size.z() - 1));

      float tsdf = readTsdf (g.x(), g.y(), g.z(), buffer);

      //infinite loop guard
      const float max_time = 3 * (volume_size.x() + volume_size.y() + volume_size.z());

      for (; time_curr < max_time; time_curr += time_step)
      {
        float tsdf_prev = tsdf;

        Vector3i g = getVoxel(ray_start + ray_dir * (time_curr + time_step));
        if (!checkInds(g))
          break;

        tsdf = readTsdf(g.x(), g.y(), g.z(), buffer);

        if (tsdf_prev < 0.f && tsdf >= 0.f)
          break;

        if (tsdf_prev >= 0.f && tsdf < 0.f)           //zero crossing
        {
          float Ftdt = interpolateTrilineary(ray_start, ray_dir, time_curr + time_step, buffer);
          if (isnan(Ftdt))
            break;

          float Ft = interpolateTrilineary(ray_start, ray_dir, time_curr, buffer);
          if (isnan(Ft))
            break;

          //float Ts = time_curr - time_step * Ft/(Ftdt - Ft);
          float Ts = time_curr - time_step * Ft / (Ftdt - Ft);

          Vector3f vetex_found = ray_start + ray_dir * Ts;

          vmap->ptr(y       )[x] = vetex_found.x();
          vmap->ptr(y + rows)[x] = vetex_found.y();
          vmap->ptr(y + 2 * rows)[x] = vetex_found.z();

          Vector3i g = getVoxel ( ray_start + ray_dir * time_curr );
          if (g.x() > 1 && g.y() > 1 && g.z() > 1 && g.x() < buffer.voxels_size.x() - 2 && g.y() < buffer.voxels_size.y() - 2 && g.z() < buffer.voxels_size.z() - 2)
          {
            Vector3f t;
            Vector3f n;

            t = vetex_found;
            t.x() += cell_size.x();
            float Fx1 = interpolateTrilineary(t, buffer);

            t = vetex_found;
            t.x() -= cell_size.x();
            float Fx2 = interpolateTrilineary(t, buffer);

            n.x() = (Fx1 - Fx2);

            t = vetex_found;
            t.y() += cell_size.y();
            float Fy1 = interpolateTrilineary(t, buffer);

            t = vetex_found;
            t.y() -= cell_size.y();
            float Fy2 = interpolateTrilineary(t, buffer);

            n.y() = (Fy1 - Fy2);

            t = vetex_found;
            t.z() += cell_size.z();
            float Fz1 = interpolateTrilineary(t, buffer);

            t = vetex_found;
            t.z() -= cell_size.z();
            float Fz2 = interpolateTrilineary(t, buffer);

            n.z() = (Fz1 - Fz2);

            n.normalize();

            nmap->ptr(y       )[x] = n.x();
            nmap->ptr(y + rows)[x] = n.y();
            nmap->ptr(y + 2 * rows)[x] = n.z();
          }
          break;
        }
      }          /* for(;;)  */
    }
  }
}




float getMinTime(const Vector3f& volume_max, const Vector3f& origin, const Vector3f& dir)
{
  float txmin = ((dir.x() > 0 ? 0.f : volume_max.x()) - origin.x()) / dir.x();
  float tymin = ((dir.y() > 0 ? 0.f : volume_max.y()) - origin.y()) / dir.y();
  float tzmin = ((dir.z() > 0 ? 0.f : volume_max.z()) - origin.z()) / dir.z();

  return std::fmax(std::fmax(txmin, tymin), tzmin);
}

float getMaxTime(const Vector3f& volume_max, const Vector3f& origin, const Vector3f& dir)
{
  float txmax = ((dir.x() > 0 ? volume_max.x() : 0.f) - origin.x()) / dir.x();
  float tymax = ((dir.y() > 0 ? volume_max.y() : 0.f) - origin.y()) / dir.y();
  float tzmax = ((dir.z() > 0 ? volume_max.z() : 0.f) - origin.z()) / dir.z();

  return std::fmin(std::fmin(txmax, tymax), tzmax);
}
