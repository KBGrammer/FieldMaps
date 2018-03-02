#include <vector>
#include "Vector3.h"

#ifndef H_MYMATH
#define H_MYMATH

double squ(double a);
double inv(double a);
double variance(double sum_x, double sum_xsq);
double intpow(double a, int p);
double gammln(double x);
unsigned long long fact(unsigned long long num);

/* double kahansum(std::vector <double> & params); */
Vector3 kahansum(std::vector <Vector3> & params);
/* template <typename T> */
/* T kahansum(std::vector <T> & params) */
/* { */
/*   T sum = params[0]; */
/*   T c = 0.0;                  // A running compensation for lost low-order bits. */
/*   for(int i = 1; i < params.size(); i++) { */
/*     T y = params[i] - c; */
/*     T t = sum + y; */
/*     c = (t - sum) - y; */
/*     sum = t; */
/*   } */
/*   return sum+c; */
/* } */

#endif
