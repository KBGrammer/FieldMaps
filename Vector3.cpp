#include "Vector3.h"
#include <cmath>
#include <iostream>
#include "mymath.h"

Vector3::Vector3()
{
  
}

Vector3::Vector3(const Vector3 & vect)
{
  set_coords(vect.xc(), vect.yc(), vect.zc());
  mag = sqrt(x_coord * x_coord + y_coord * y_coord + z_coord * z_coord);
}

Vector3::Vector3(double x1, double y1, double z1)
{
  x_coord = x1; y_coord = y1; z_coord = z1;
  mag = sqrt(x1*x1+y1*y1+z1*z1);
}

Vector3::~Vector3()
{

}

Vector3& Vector3::operator=(const Vector3& rhs)
{
  if (this == &rhs) 
    return *this;

  double tx, ty, tz;
  rhs.return_coords(tx, ty, tz);
  set_coords(tx, ty, tz);
  
  return *this;
}

void Vector3::set_coords(double x1, double y1, double z1)
{
  x_coord = x1; y_coord = y1; z_coord = z1;
  mag = sqrt(x1*x1+y1*y1+z1*z1);
}

void Vector3::return_coords(double & x1, double & y1, double & z1) const
{
  x1 = x_coord;
  y1 = y_coord;
  z1 = z_coord;
}

void Vector3::scalefactor(double sca)
{
  set_coords(x_coord*sca, y_coord*sca, z_coord*sca);
  mag = mag * sca;
}
void Vector3::scalecomponents(double s1, double s2, double s3)
{
  x_coord *= s1;
  y_coord *= s2;
  z_coord *= s3;
}

double Vector3::xc() const
{
  return x_coord;
}

double Vector3::yc() const
{
  return y_coord;
}

double Vector3::zc() const
{
  return z_coord;
}

double Vector3::distance(const Vector3 & vec)
{
  return sqrt(squ(x_coord - vec.xc()) +
	      squ(y_coord - vec.yc()) +
	      squ(z_coord - vec.zc()));
}

double Vector3::distance(const Vector3 & vec) const
{
  return sqrt(squ(x_coord - vec.xc()) +
	      squ(y_coord - vec.yc()) +
	      squ(z_coord - vec.zc()));
}

double Vector3::distancesq(const Vector3 & vec)
{
  return squ(x_coord - vec.xc()) +
    squ(y_coord - vec.yc()) +
    squ(z_coord - vec.zc());
}

double Vector3::dot(const Vector3 & vec)
{
  double tx, ty, tz;
  vec.return_coords(tx, ty, tz);
  return (x_coord*tx+
	  y_coord*ty+
	  z_coord*tz);
}

Vector3 Vector3::cross(const Vector3 & vec)
{
  double tx, ty, tz;
  vec.return_coords(tx, ty, tz);

  std::vector <Vector3> vectors;

  // vectors.push_back(Vector3(y_coord*tz,z_coord*tx,x_coord*ty));
  // vectors.push_back(Vector3(-z_coord*ty,-x_coord*tz,-y_coord*tx));
  
  // Vector3 crossprod = kahansum(vectors);
  
  Vector3 crossprod((y_coord*tz-z_coord*ty),
  		    (z_coord*tx-x_coord*tz),
  		    (x_coord*ty-y_coord*tx));
  return crossprod;
}

Vector3 Vector3::projection(const Vector3 & vec)
{
  return (((this->dot(vec)) / (vec.magnitude()*vec.magnitude())) * vec);
}

double Vector3::magnitude()
{
  return mag;
}

double Vector3::magnitude() const
{
  return mag;
}

// Overloaded friends!!

std::ostream& operator<< (std::ostream &out, Vector3 &tVect)
{
    // Since operator<< is a friend of the Point class, we can access
    // Point's members directly.
  out << tVect.xc() << " " << tVect.yc() << " " << tVect.zc();
  return out;
}

std::istream& operator>> (std::istream &in, Vector3 &tVect)
{
    // Since operator<< is a friend of the Point class, we can access
    // Point's members directly.
  double tx, ty, tz;
  in >> tx >> ty >> tz;
  tVect.set_coords(tx, ty, tz);
  return in;
}

Vector3 operator+(const Vector3 &v1, const Vector3 &v2)
{
  return Vector3(v1.xc()+v2.xc(), v1.yc()+v2.yc(), v1.zc()+v2.zc());
}

Vector3 operator+(const Vector3 &v1, const double & scalar)
{
  return Vector3(v1.xc()+scalar, v1.yc()+scalar, v1.zc()+scalar);
}

Vector3 operator-(const Vector3 &v1, const Vector3 &v2)
{
  return Vector3(v1.xc()-v2.xc(), v1.yc()-v2.yc(), v1.zc()-v2.zc());
}

Vector3 operator-(const Vector3 &v1, const double & scalar)
{
  return Vector3(v1.xc()-scalar, v1.yc()-scalar, v1.zc()-scalar);
}

Vector3 operator*(const Vector3 &v1, const double &scalar)
{
  return Vector3(v1.xc()*scalar, v1.yc()*scalar, v1.zc()*scalar);
}

Vector3 operator*(const double &scalar, const Vector3 &v1)
{
  return Vector3(v1.xc()*scalar, v1.yc()*scalar, v1.zc()*scalar);
}

Vector3 operator/(const Vector3 &v1, const double &scalar)
{
  return Vector3(v1.xc()/scalar, v1.yc()/scalar, v1.zc()/scalar);
}

