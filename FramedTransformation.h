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

#ifndef FramedTransformation_h
#define FramedTransformation_h

#include <stdio.h>
#include "Eigen.h"

using namespace std;

struct FramedTransformation
{
  enum RegistrationType { Kinfu = 0, DirectApply = 1, InitializeOnly = 2, IncrementalOnly = 3 };
  enum ActionFlag {
    ResetFlag = 0x1,                // if reset at the very beginning
    IgnoreRegistrationFlag = 0x2,		// if discard the registration
    IgnoreIntegrationFlag = 0x4,		// if discard integration
    PushMatrixHashFlag = 0x8,       // if push the transformation matrix into the hash table
    SavePointCloudFlag = 0x10,      // if save point cloud after execution
    SaveAbsoluteMatrix = 0x20,      // if save absolute matrix, work with IgnoreIntegrationFlag
    ExtractSLACMatrix = 0x40,       // if extract SLAC matrix
  };

  int id1;
  int id2;
  int frame;
  RegistrationType type;
  int flag;
  Eigen::Matrix4f transformation;

  FramedTransformation() : type( Kinfu ), flag( 0 ) {}
  FramedTransformation( int id1, int id2, int f, Eigen::Matrix4f t ) : id1( id1 ), id2( id2 ), frame( f ), transformation( t ), type( DirectApply ), flag( 0 ) {}
  FramedTransformation( int id1, int id2, int f, Eigen::Matrix4f t, RegistrationType tp, int flg )
    : id1( id1 ), id2( id2 ), frame( f ), transformation( t ), type( tp ), flag( flg ) {}
};


#endif /* FramedTransformation_h */
