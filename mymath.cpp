#include "mymath.h"
#include <math.h>
#include <iostream>
#include <vector>

using namespace::std;

double squ(double a)
{
  return a*a;
}

double inv(double a)
{
  return 1/a;
}

double variance(double sum_x, double sum_xsq)
{
  return sum_xsq - squ(sum_x);
}

double intpow(double a, int p)
{
  double temp = a;
  if(p == 1) {
    return a;
  }
  if(p == 0) {
    return 1;
  }
  if(p > 1) {
    for(int i = 2; i <= p; i++) {
      temp = temp * a;
    }
    return temp;
  }
  else {
    temp = 1;
    for(int i = -1; i >= p; i--) {
      temp = temp / a;
    }
    return temp;
  }
}

unsigned long long fact(unsigned long long num)
{
  unsigned long long temp = 1;
  if(num != 0 && num != 1) {
    for(unsigned long long i = 2; i < num; i++) {
      temp = temp * i;
    }
  }
  return temp;
}

double gammln(const double xx)
{
  double x,tmp,y,ser;
  static const double cof[14]={57.1562356658629235,-59.5979603554754912,
                               14.1360979747417471,-0.491913816097620199,
			       .339946499848118887e-4,
                               .465236289270485756e-4,-.983744753048795646e-4,
			       .158088703224912494e-3,
                               -.210264441724104883e-3,.217439618115212643e-3,
			       -.164318106536763890e-3,
                               .844182239838527433e-4,-.261908384015814087e-4,
			       .368991826595316234e-5};
  if (xx <= 0) throw("bad arg in gammln");
  y=x=xx;
  tmp = x+5.24218750000000000;
  tmp = (x+0.5)*log(tmp)-tmp;
  ser = 0.999999999999997092;
  for (int j=0;j<14;j++) ser += cof[j]/++y;
  return tmp+log(2.5066282746310005*ser/x);
}

// double kahansum(double a, double b)
// {
//   double sum = a;
//   double c = 0.0;                  // A running compensation for lost low-order bits.
//   double y = b - c;
//   double t = sum + y;
//   c = (t - sum) - y;
//   sum = t + c;
//   return sum;
// }

// double kahansum(std::vector <double> & params)
// {
//   double sum = params[0];
//   double c = 0.0;                  // A running compensation for lost low-order bits.
//   for(int i = 1; i < params.size(); i++) {
//     double y = params[i];// - c;
//     double t = sum + y;
//     c = (t - sum) - y;
//     sum = t;
//   }
//   //std::cerr << sum << ' ' << c << std::endl;
//   return sum+c;
// }

Vector3 kahansum(std::vector <Vector3> & params)
{
  Vector3 sum = params[0];
  //Vector3 sum1 = params[0];
  // A running compensation for lost low-order bits.
  Vector3 c(0.0,0.0,0.0);

  for(int i = 1; i < params.size(); i++) {
    Vector3 y = params[i] - c;
    Vector3 t = sum + y;
    c = (t - sum) - y;
    sum = t;
    //std::cerr << c << std::endl;
    //sum1 = sum1 + params[i];
  }

  //std::cerr << sum << ' ' << sum1 << std::endl;

  return sum;
}

