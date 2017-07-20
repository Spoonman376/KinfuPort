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

#ifndef Eigen_h
#define Eigen_h

#include <Eigen/Core>

#include "Utils.h"

typedef Eigen::Transform<float, 3, Eigen::Affine> Affine3f;
typedef Eigen::Matrix<float, 3, 3, Eigen::RowMajor> Matrix3frm;

typedef Eigen::Vector3f Vector3f;
typedef Eigen::Vector3i Vector3i;


template<> template<>
Eigen::MatrixBase<Eigen::Block<Matrix3frm, 1, 3, true> >::cross_product_return_type<Eigen::Block<Matrix3frm, 1, 3, true> >::type Eigen::MatrixBase<Eigen::Block<Matrix3frm, 1, 3, true> >::cross<Eigen::Block<Matrix3frm, 1, 3, true> >(Eigen::MatrixBase<Eigen::Block<Matrix3frm, 1, 3, true> > const&) const;




struct Eigen33
{
public:
  Eigen33(volatile float* mat_pkg_arg) : mat_pkg(mat_pkg_arg) {}

  static Vector3f unitOrthogonal (const Vector3f &src);

  void computeRoots2(const float& b, const float& c, Vector3f& roots);

  void computeRoots3(float c0, float c1, float c2, Vector3f& roots);

  void compute(Matrix3frm& tmp, Matrix3frm& vec_tmp, Matrix3frm& evecs, Vector3f& evals);


private:
  volatile float* mat_pkg;

  float m00() const { return mat_pkg[0]; }
  float m01() const { return mat_pkg[1]; }
  float m02() const { return mat_pkg[2]; }
  float m10() const { return mat_pkg[1]; }
  float m11() const { return mat_pkg[3]; }
  float m12() const { return mat_pkg[4]; }
  float m20() const { return mat_pkg[2]; }
  float m21() const { return mat_pkg[4]; }
  float m22() const { return mat_pkg[5]; }

  Vector3f row0() const { return Vector3f( m00(), m01(), m02() ); }
  Vector3f row1() const { return Vector3f( m10(), m11(), m12() ); }
  Vector3f row2() const { return Vector3f( m20(), m21(), m22() ); }

  static bool isMuchSmallerThan (float x, float y);

};


#endif /* Eigen_h */
