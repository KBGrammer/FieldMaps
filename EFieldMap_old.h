
#include "EFieldPoint.h"
#include <vector>


#ifndef E_FIELD_MAP_H
#define E_FIELD_MAP_H

 
class EFieldMap
{
private:
  std::vector < EFieldPoint > map;

  std::vector < int > bins;
  std::vector < double > keys;
  double bin_spacing;
  double stepsize;
 
 
public:
  //EFieldMap();
  EFieldMap(double tstep) {stepsize = tstep; } // private default constructor
 
  void ReadEFieldMap();

  void PrintMap(int count);

  void PrintKeysBins();

  void ListPointsinRad(Vector3 vect, double rad);
  
  Vector3 Interpolate_Field(Vector3 vect) const;
  Vector3 Interpolate_Field_Gradient(Vector3 vect, Vector3 &gradient);
  Vector3 Interpolate_Gradient(Vector3 vect, Vector3 &gradient);

};

#endif
