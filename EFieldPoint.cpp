#include "EFieldPoint.h"
#include "Vector3.h"
 
// EFieldPoint constructor
EFieldPoint::EFieldPoint(Vector3 coords, Vector3 efield)
{
  SetEFieldPoint(coords, efield);
}
 
// EFieldPoint member function
void EFieldPoint::SetEFieldPoint(Vector3 coords, Vector3 efield)
{
  m_coords = coords;
  m_efield = efield;
}

Vector3 EFieldPoint::GetPoint() const
{
  return m_coords;
}

Vector3 EFieldPoint::GetField() const
{
  return m_efield;
}
