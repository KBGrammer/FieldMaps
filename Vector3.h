#include <iostream>

#ifndef VECTOR3_H
#define VECTOR3_H

class Vector3 {
 private:
  double coords[3];
  double x_coord;
  double y_coord;
  double z_coord;
  double mag;

 public:
  Vector3();
  Vector3(const Vector3 & vect);
  Vector3(double tx, double ty, double tz);
  Vector3& operator=(const Vector3& rhs);
  ~Vector3();

  void set_coords(double tx, double ty, double tz);
  void return_coords(double & tx, double & ty, double & tz) const;
  double xc() const;
  double yc() const;
  double zc() const;
  double distance(const Vector3 & vec);
  double distance(const Vector3 & vec) const;
  double distancesq(const Vector3 & vec);
  void scalefactor(double sca);
  void scalecomponents(double s1, double s2, double s3);
  // void make_unitvector();
  double dot(const Vector3 & vec);
  double magnitude();
  double magnitude() const;
  Vector3 cross(const Vector3 & vec);
  Vector3 projection(const Vector3 & vec);



  friend std::ostream& operator<< (std::ostream &out, Vector3 &tVect);
  friend std::istream& operator>> (std::istream &in, Vector3 &tVect);
  friend Vector3 operator+(const Vector3 &v1, const Vector3 &v2);
  friend Vector3 operator+(const Vector3 &v1, const double & scalar);
  friend Vector3 operator-(const Vector3 &v1, const Vector3 &v2);
  friend Vector3 operator-(const Vector3 &v1, const double & scalar);
  friend Vector3 operator*(const Vector3 &v1, const double & scalar);
  friend Vector3 operator*(const double & scalar, const Vector3 &v1);
  friend Vector3 operator/(const Vector3 &v1, const double & scalar);

};


#endif
