#include "Vector3.h"
#include <cmath>
#include <iostream>
#include "mymath.h"
#include "Matrix.h"


Matrix::Matrix()
{
  
}

Matrix::Matrix(std::vector <double> els)
{
  number_elements = els.size();
  for(int i = 0; i < number_elements; i++) {
    elements.push_back(els[i]);
  }
}

Matrix::~Matrix()
{

}

Vector3 Matrix::Multiply(Vector3 vec)
{
  double tx = elements[0]*vec.xc() +
    elements[1] * vec.yc() +
    elements[2] * vec.zc();
  double ty = elements[3]*vec.xc() +
    elements[4] * vec.yc() +
    elements[5] * vec.zc();
  double tz = elements[6]*vec.xc() +
    elements[6] * vec.yc() +
    elements[7] * vec.zc();
  
  Vector3 result(tx,ty,tz);
  return result;
}
