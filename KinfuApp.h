//
//

#ifndef KinfuApp_h
#define KinfuApp_h

#include <stdio.h>
#include <iostream>
#include <openni2/OpenNI.h>
#include <pcl/pcl_base.h>
#include <opencv2/opencv.hpp>


#include "CameraParams.h"
//#include "FramedTransformation.h"
#include "Kinfu.h"


using namespace std;
using namespace openni;
using namespace cv;

class KinfuApp
{
protected:
  Device device;
  VideoStream depthStream, colourStream;

  int frameCtr, fragmentRate, fragmentStart;
  Mat grayImage0, grayImage1;
  Mat depthFilter0, depthFilter1;

  CameraParams camera;
  FramedTransformation framedTransformation;

  KinfuTracker kinfu;

public:
  string oniFileName, cameraParamFileName;

  KinfuApp();
  ~KinfuApp();
  void mainLoop();




protected:
  void execute(VideoFrameRef depth, VideoFrameRef colour);

  void processFramedTransformation(int frameId);

  bool initalizeOpenNI();
  void closeOpenNI();





};




#endif /* KinfuApp_h */
