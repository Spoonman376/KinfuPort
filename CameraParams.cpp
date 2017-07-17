//
//

#include "CameraParams.h"


CameraParams::CameraParams()
{
  //Just some initial values
  fx = 525.0;
  fy = 525.0;
  cx = 319.5;
  cy = 239.5;
  ICPTrunc = 2.5;
  integrationTrunc = 2.5;
}

void CameraParams::loadFromFile(std::string filename) {
  FILE * f = fopen( filename.c_str(), "r" );
  if ( f != NULL ) {
    char buffer[1024];
    while ( fgets( buffer, 1024, f ) != NULL ) {
      if ( strlen( buffer ) > 0 && buffer[ 0 ] != '#' ) {
        sscanf( buffer, "%lf", &fx);
        fgets( buffer, 1024, f );
        sscanf( buffer, "%lf", &fy);
        fgets( buffer, 1024, f );
        sscanf( buffer, "%lf", &cx);
        fgets( buffer, 1024, f );
        sscanf( buffer, "%lf", &cy);
        fgets( buffer, 1024, f );
        sscanf( buffer, "%lf", &ICPTrunc);
        fgets( buffer, 1024, f );
        sscanf( buffer, "%lf", &integrationTrunc);
      }
    }
    fclose ( f );
  }
}
