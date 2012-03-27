#include "bp.hpp"
#include "randomc.h"
#include <iostream>

real_t BP::bpdistance(Trig & t , Plane & p)
{
  real_t cost = (t.n-p.n).norm();
  cost += distw*((t.c-p.c).norm());
  cost /= (1+distw);
  return cost;
}

void BP::update_distance(std::vector<real_t > & dist,
                         std::vector<real_t > & cdf,
                      std::vector<Plane> & plane,
                      Mesh & m, int ll)
{
  real_t sum=0;
  for(size_t ii=0;ii<m.t.size();ii++){
    real_t d = bpdistance(m.t[ii], plane[ll]);
    if(d<dist[ii] ){
      dist[ii]=d;
      m.t[ii].label=ll;
    }
    sum+=dist[ii]*dist[ii];
    cdf[ii]=sum;
  }
}



BP::BP(Mesh & m):distw(1.0f)
{
  plane.resize(m.nLabel);
  m.get_normal_center();
  std::vector<real_t > dist(m.t.size(), 0.0f) ;
  std::vector<real_t > cdf(m.t.size(),0.0f);
//1. pick a center uniform at random
  std::cout<<"kmeans++ init\n";
  int r = rndg.IRandom(0, m.t.size()-1 );
  plane[0].n = m.t[r].n;
  plane[0].c = m.t[r].c;
  real_t sum=0;
  for(size_t ii=0;ii<m.t.size();ii++){
    dist[ii]=bpdistance(m.t[ii],plane[0]);
    m.t[ii].label=0;
    sum+=dist[ii]*dist[ii];
    cdf[ii]=sum;
  }

  for(int kk=1;kk<m.nLabel;kk++){
    int center = sample_cdf(cdf);
    plane[kk].n=m.t[center].n;
    plane[kk].c=m.t[center].c;
    //std::cout<<kk<<"\n";
    update_distance(dist, cdf ,plane,m,kk);
  }

}
