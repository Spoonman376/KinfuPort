//
//

#ifndef CameraParams_h
#define CameraParams_h

#include <stdio.h>
#include <iostream>


struct CameraParams
{
public:
	float fx, fy, cx, cy, ICPTrunc, integrationTrunc;

	CameraParams();

	void loadFromFile(std::string filename);
};

#endif /* CameraParams_h */
