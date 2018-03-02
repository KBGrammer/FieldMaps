#include "Vector3.h"

#ifndef E_FIELD_POINT_H
#define E_FIELD_POINT_H
 
class EFieldPoint
{
private:
  Vector3 m_coords;
  Vector3 m_efield;
 
 
public:
    EFieldPoint() { } // private default constructor
    EFieldPoint(Vector3 coords, Vector3 efield);
 
    void SetEFieldPoint(Vector3 coords,	Vector3 efield);
 
    Vector3 GetPoint() const;
    Vector3 GetField() const;

};
 
#endif
