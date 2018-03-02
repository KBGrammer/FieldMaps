#include "rng.h"
#include <math.h>
// From NR book.
// Period = 3*10^57

Rng::Rng(unsigned long long j) :v(4101842887655102017ULL), w(1)
{
  exc_div = 1/(pow(2,64)-1);
  inc_div = 1/(pow(2,64)-2);
  u = j ^ v;
  rint64();
  v = u;
  rint64();
  w = v;
  rint64();
}

Rng::~Rng()
{

}

// Returns random 64 bit integer on [0,2^64-2]
unsigned long long Rng::rint64()
{
  u = u * 2862933555777941757ULL + 7046029254386353087ULL;
  v ^= v >> 17;
  v ^= v << 31;
  v ^= v >> 8;
  unsigned long long x = u ^ (u << 21);
  x ^= x >> 35;
  x ^= x << 4;
  return ((x + v) ^  w) - 1; // -1 to allow zero
}

// Random double on [0,1] inclusive
double Rng::rdouble()
{
  return inc_div * rint64();
}

// Random float on [0,1] inclusive
float Rng::rfloat()
{
  return float(inc_div*rint64());
}

// Random double on [0,1)
double Rng::rdouble_exc()
{
  return exc_div * rint64();
}

// Random integer on [0,2^32-1]
unsigned int Rng::rint32()
{
  return (unsigned int)rint64();
}


GRng::GRng(unsigned long long l) :Rng(l)
{
  mu = 0;
  sigma = 0.1;
}

GRng::GRng(double m, double s, unsigned long long l) :Rng(l)
{
  mu = m;
  sigma = s;
}

GRng::~GRng()
{

}

void GRng::set_gauss(double m, double s)
{
  mu = m;
  sigma = s;
}

/*
  Normal deviates by ratio of uniforms
  Numerical recipes 3rd edition page 369
 */
double GRng::rgauss()
{
  double a,b,c,d,e;

  do {
    a = rdouble();
    b = 1.7156*(rdouble()-0.5);
    c = a - 0.449871;
    d = fabs(b)+0.386595;
    e = c*c + d*(0.19600*d-0.25472*c);
  }while(e > 0.27597 && (e > 0.27846 || b*b > -4.*log(a)*a*a));
  return mu+sigma*b/a;
}

MTRng::MTRng(unsigned long int j)
{
  initialize(j);
}

void MTRng::initialize(unsigned long int seed)
{
  index = 0;
  MT[0] = seed;
  for(int i = 1; i < 624; i++) {
    MT[i] = 0xffffffff & (1812433253 * MT[i-1] ^ (MT[i-1] << 30) + i);
  }
}

void MTRng::generate_numbers()
{
  for(int i = 0; i < 624; i++) {
    unsigned long int y = ((MT[i] & 0x80000000) >> 31) 
      + ( MT[(i+1) % 624] & 0x8fffffff);

    if(y&1) 
      MT[i] = MT[(i+397) % 624] ^ (y >> 1);
    else
      MT[i] = MT[(i+397) % 624] ^ (y >> 1) ^ 2567483615ul; // 0x9908b0df; ul indicates unsigned long in C
  }
  index = 0;
}

unsigned long int MTRng::rint32()
{
  if(index >= 624)
    generate_numbers();

  unsigned long int y = MT[index];
  y = y ^ (y >> 11);
  y = y ^ ((y << 7) & 2636928640ul); // 0x9d2c5680
  y = y ^ ((y << 15) & 4022730752ul); // 0xefc60000
  y = y ^ (y >> 18);
  index++; // increment index
  return y;
}

unsigned long long MTRng::rint64()
{
  unsigned long long y1 = rint32();
  unsigned long long y2 = rint32();
  
  return (y1 << 32) + y2;
}

double MTRng::rdouble()
{
  return rint32() / (pow(2,32)-2);
}
