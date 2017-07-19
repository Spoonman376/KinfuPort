//
//

#ifndef TSDFBuffer_h
#define TSDFBuffer_h

#include <stdio.h>

#include "Eigen.h"


struct short2
{
  short x, y;

  short2(short a, short b) {x = a; y = b;}
};

struct uchar4
{
  unsigned char x, y, z, w;

  uchar4(unsigned char a, unsigned char b, unsigned char c, unsigned char d)
  {
    x = a; y = b; z = c; w = d;
  }
};

struct tsdf_buffer
{
  /** \brief Address of the first element of the TSDF volume in memory*/  
  short2* tsdf_memory_start;
  /** \brief Address of the last element of the TSDF volume in memory*/          
  short2* tsdf_memory_end;
  /** \brief Memory address of the origin of the rolling buffer. MUST BE UPDATED AFTER EACH SHIFT.*/
  short2* tsdf_rolling_buff_origin;

  uchar4* color_memory_start;
  uchar4* color_memory_end;
  uchar4* color_rolling_buff_origin;

  /** \brief Internal cube origin for rollign buffer.*/
  Vector3i origin_GRID;
  /** \brief Cube origin in world coordinates.*/
  Vector3f origin_GRID_global;
  /** \brief Current metric origin of the cube, in world coordinates.*/ 
  Vector3f origin_metric;
  /** \brief Size of the volume, in meters.*/
  Vector3f volume_size; //3.0
  /** \brief Number of voxels in the volume, per axis*/
  Vector3i voxels_size; //512

  /** \brief Default constructor*/
  tsdf_buffer () 
  {
    tsdf_memory_start = 0;  tsdf_memory_end = 0; tsdf_rolling_buff_origin = 0; 
    color_memory_start = 0; color_memory_end = 0; color_rolling_buff_origin = 0;
    origin_GRID.x() = 0; origin_GRID.y() = 0; origin_GRID.z() = 0;
    origin_GRID_global.x() = 0.f; origin_GRID_global.y() = 0.f; origin_GRID_global.z() = 0.f;
    origin_metric.x() = 0.f; origin_metric.y() = 0.f; origin_metric.z() = 0.f;
    volume_size.x() = 3.f; volume_size.y() = 3.f; volume_size.z() = 3.f;
    voxels_size.x() = 512; voxels_size.y() = 512; voxels_size.z() = 512;
  }          

};


#endif /* TSDFBuffer_h */
