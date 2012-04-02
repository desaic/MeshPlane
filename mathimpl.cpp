#include "math.hpp"
#include "randomc.h"

Vec3 operator*(real_t c, const Vec3 & v)
{
  return v*c;
}

real_t logsum(std::vector<real_t > & x)
{
  real_t m = x[0];
  for(size_t ii=0; ii<x.size(); ii++) {
    if(x[ii]>m) {
      m=x[ii];
    }
  }
  real_t sum=0;
  for(size_t ii=0; ii<x.size(); ii++) {
    sum+=std::exp(x[ii]-m);
  }
  return std::log(sum)+m;
}

int bsrch(std::vector<real_t> & val , real_t x)
{
  int bottom = 0;
  int top = val.size()- 1;
  int middle;

  while (bottom <= top) {

    middle = (bottom + top) / 2;

    if (val[middle] < x) {
      bottom = middle + 1;
    } else {
      top = middle - 1;
    }
  }
  return top;
}

int sample_cdf(std::vector<real_t > &cdf)
{
  double r = rndg.Random();
  r*=cdf[cdf.size()-1];
  int label = bsrch(cdf, r);
  return label;
}
