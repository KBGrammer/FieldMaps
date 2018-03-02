#include <iostream>
#include <vector>
#include "Vector3.h"

#ifndef MATRIX_H
#define MATRIX_H

class Matrix {
 private:

  std::vector <double> elements;
  int number_elements;

 public:

  Matrix();
  Matrix(std::vector <double> els);
  ~Matrix();

  Vector3 Multiply(Vector3 vec);

};

#endif
