#include "kmeans.hpp"
#include "randomc.h"
#include <list>
#include <vector>
#include <queue>
#include <iostream>
extern float distw;
extern float distCenterW;
real_t tpdistance(Trig & t , Plane & p);

void update_distance(std::vector<real_t > & dist,
                      std::vector<real_t > & cdf,
                      std::vector<Plane> & plane,
                      Mesh & m, int ll)
{
  real_t sum=0;
  for(size_t ii=0;ii<m.t.size();ii++){
    real_t d = tpdistance(m.t[ii], plane[ll]);
    if(d<dist[ii] ){
      dist[ii]=d;
      m.t[ii].label=ll;
    }
    sum+=dist[ii]*dist[ii];
    cdf[ii]=sum;
  }
}

real_t tpdistance(Trig & t , Plane & p)
{
  real_t cost = (t.n-p.n).norm();
  float nlen = (t.c - p.c).dot(p.n) ;
  cost += distw*std::abs( nlen);
  cost += distCenterW * (t.c - p.c - nlen * p.n).norm();
  cost /= (1+distw + distCenterW);
  return cost;
}

struct TrigDist{
  int idx;
  int label;
  real_t dist;
  TrigDist(int _idx=0, int _label=0, real_t _dist=0):
  idx(_idx),label(_label),dist(_dist){}
  bool operator<(const TrigDist & b) const {
    if(idx<b.idx){
      return true;
    }
    if(idx==b.idx){
      return label<b.label;
    }
    return false;
  }
};

class CmpDist
{
  public:
  CmpDist(){}
  bool operator() (const TrigDist& lhs, const TrigDist&rhs) const
  {
    return lhs.dist>rhs.dist;
  }
};

void findCentroids(Mesh & m, std::vector<Plane> & plane, std::vector<int> &centroid)
{
  m.get_normal_center();
  get_plane(m,plane);
  std::vector<real_t> minDist ;
  centroid.resize(plane.size());
  centroid.assign(centroid.size(),-1);
  minDist.resize(centroid.size(),100);
  for(size_t ii=0;ii<m.t.size();ii++){
    int ll = m.t[ii].label;
    real_t dist = tpdistance(m.t[ii],plane[ll]);
    if(centroid[ll]<0
       ||dist < minDist[ll]){
      centroid[ll]=ii;
      minDist[ll]=dist;
      continue;
    }
  }
}

void flood(Mesh & m, const std::vector<int> & centroids,
           std::vector<Plane > & plane)
{
  CmpDist cmpdist;
  std::priority_queue< TrigDist, std::vector<TrigDist>, CmpDist> pq;
  std::vector<bool> labeled(m.t.size());
  for(size_t ii=0;ii<centroids.size();ii++){
    int cidx=centroids[ii];
    if(cidx<0){
      std::cout<<"wtf\n";
    }
    real_t dist = tpdistance(m.t[cidx], plane[ii]);
    pq.push(TrigDist(cidx, m.t[cidx].label, dist) );
  }
  std::map<TrigDist , bool> processed;
  while(!pq.empty()){
    TrigDist trigdist= pq.top();
    pq.pop();
    int tidx=trigdist.idx;
    int label=trigdist.label;
    if(!labeled[tidx]){
      m.t[tidx].label=label;
      labeled[tidx]=true;
    }else{
      continue;
    }

    for(size_t ii=0;ii<m.adjMat[tidx].size();ii++){
      int nbrIdx=m.adjMat[tidx][ii];
      real_t nbrDist=tpdistance(m.t[nbrIdx],plane[label]);
      TrigDist nbrTd(nbrIdx,label,nbrDist );
      if(processed[nbrTd]){
        continue;
      }
      processed[nbrTd]=true;
      pq.push(nbrTd);
    }
  }
}

void initKmeans(Mesh & m)
{
  std::vector<Plane> plane;
  std::vector<int> centroids;
  m.adjlist();
  plane.resize(m.nLabel);
  centroids.resize(m.nLabel);
  m.get_normal_center();
  std::vector<real_t > dist(m.t.size(), 0.0f) ;
  std::vector<real_t > cdf(m.t.size(),0.0f);
//1. pick a center uniform at random
  std::cout<<"kmeans++ init\n";
  int r = rndg.IRandom(0, m.t.size()-1 );
  plane[0].n = m.t[r].n;
  plane[0].c = m.t[r].c;
  centroids[0]=r;
  real_t sum=0;
  for(size_t ii=0;ii<m.t.size();ii++){
    dist[ii]=tpdistance(m.t[ii],plane[0]);
    m.t[ii].label=0;
    sum+=dist[ii]*dist[ii];
    cdf[ii]=sum;
  }

  for(int kk=1;kk<m.nLabel;kk++){
    //kmeans++
	  //int center = sample_cdf(cdf);
	//furthest sampling
	  float maxDist = 0;
	  int center = 0;
	  for (int j = 0; j < dist.size(); j++) {
		  if (dist[j] > maxDist) {
			  maxDist = dist[j];
			  center = j;
		  }
	  }
	  
    centroids[kk]=center;
    plane[kk].n=m.t[center].n;
    plane[kk].c=m.t[center].c;
 //   std::cout<<kk<<": "<< center <<"\n";
	//std::cout << plane[kk].c[0] << " " << plane[kk].c[1] << " " << plane[kk].c[2] << "\n";
    update_distance(dist, cdf ,plane,m,kk);
  }
  m.planes = plane;
  flood(m,centroids, plane);
}

void runKmeans(Mesh & m)
{
  std::vector<Plane> plane;
  std::vector<int>centroid;
  findCentroids(m,plane,centroid);
  flood(m,centroid,plane);
}
