/////////////////////////
//
// Originally submitted to the Elog by Jonathan Mulholland May 6th 2013
//
// Butchered into class form by Kyle Grammer
//
// g++ compilation line: g++ -o coils.exe coils.cxx -lgsl -lgslcblas -lm
// no dependency files necessary
// must have installed the gnu scientific libraries
#include <vector>
#include <gsl/gsl_integration.h>


#ifndef MAG_FIELD_MAP_H
#define MAG_FIELD_MAP_H

 
class MagFieldMap
{
private:
  double angleL[11];                        // tilt of the coil (degrees) (0 degree aligns with the z-axis)
  double xL[11],yL[11],zL[11];              // the middle of the centeral axis of the coil (meters)
  double lengthL[11],thickL[11],inRL[11];   // length of the coil, thickness of the winding, inner radius of the coil (meters)
  double densityL[11];                      // current density: amps per m^2
  //  For breaking each coil into loops for the BiotSavart treatment
  int divL[11];   // how many segments to break the length of the coil into 
  int divR[11];   // how many segments to break the thickness of the winding into
  double current;                           // magnet current (Amps)
  double mu0;               // permeability of free space
  double d2r;              // convert degrees to radians
  //double position[] = {0.00,0.00,0.00,0.0}; // (x,y,z,r) position of field point and loop radius passed to integrand
 
 
public:

  MagFieldMap(); // private default constructor
  ~MagFieldMap(); // private default destructor

  Vector3 FieldEval(Vector3 origin) const;

  static double bx (double theta, void * params)
  {
    double x = *((double *) params);
    double y = *((double *) params+1);
    double z = *((double *) params+2);
    double R = *((double *) params+3);
    double bx = -R/4.0/3.141593*z*cos(theta)/pow(gsl_hypot3(x-R*cos(theta),y-R*sin(theta),z),3);
    return bx;
  }

  static double by (double theta, void * params)
  {
    double x = *((double *) params);
    double y = *((double *) params+1);
    double z = *((double *) params+2);
    double R = *((double *) params+3);
    double by = -R/4.0/3.141593*z*sin(theta)/pow(gsl_hypot3(x-R*cos(theta),y-R*sin(theta),z),3);
    return by;
  }

  static double bz (double theta, void * params)
  {
    double x = *((double *) params);
   double y = *((double *) params+1);
   double z = *((double *) params+2);
   double R = *((double *) params+3);
   double bz = R/4.0/3.141593 * (-y*sin(theta) - x*cos(theta) + R )/pow(gsl_hypot3(x-R*cos(theta),y-R*sin(theta),z),3);
   return bz;
 }

};

#endif
