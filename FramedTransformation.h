//
//

#ifndef FramedTransformation_h
#define FramedTransformation_h

#include <stdio.h>
#include <Eigen/Core>

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
