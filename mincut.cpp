#include "mincut.hpp"
#include <vector>
#include <stdlib.h>
#include <map>
#include <utility>
#include "cgd.hpp"
#define SWAP(a,b)  (a)^=(b);(b)^=(a);(a)^=(b);
int minc_nlabel=50;
float dataCostW=100;
float smoothW=10;
float distw = 10;
void data_cost(Mesh & m, int nLabel, std::vector<Plane>&plane,
	       std::vector<std::vector< float > > & datac)
{
  //recalculate center and normals for each face since the mesh is morphed

  m.get_normal_center();
  randcenter(m,plane,nLabel);
  datac.resize(m.t.size());

  float mn = -1;
  float mx = -1;
  for (unsigned int ii=0;ii<m.t.size();ii++){
    datac[ii].resize(nLabel);
    for (int jj =0 ;jj<nLabel;jj++){
      Vec3 plane_d = plane[jj].n.dot(plane[jj].c);
      Vec3 d = m.t[ii].c.dot(plane[jj].c);
      float cost = distw * (plane_d-d).L1n();
      cost += (m.t[ii].n- plane[jj].n).L1n();
      cost /= (1+distw);
      datac[ii][jj]=cost;
      if(mn<0 || mn>cost){
	mn=cost;
      }
      if(mx<0 || mx<cost){
	mx=cost;
      }
    }
  }
  float scale = mx-mn;
  for (unsigned int ii=0;ii<m.t.size();ii++){
    for (int jj =0 ;jj<nLabel;jj++){
      datac[ii][jj] -= mn;
      datac[ii][jj] /= scale;
    }
  }
}

float distance( Trig & a, Trig &b)
{
  float d = (a.n-b.n).L1n();
  return d;
}

/**@param smoothc smoothc[i][j] is the smoothness cost for
   vertices i and j. i < j.
 */
void smooth_cost(Mesh& m,
		 std::map<std::pair< int , int> , float > & smoothc)
{
  float mn = -1;
  float mx = -1;

  for(unsigned int ii=0;ii< m.t.size();ii++){
    for(unsigned int nbr=0;nbr<m.adjMat[nbr].size();nbr++ ){
      unsigned int nbrIdx = m.adjMat[ii][nbr];
      if(nbrIdx<ii){
	continue;
      }
      float cost = distance( m.t[ii], m.t[nbrIdx]);
      cost = 1/(1+cost);
      smoothc[std::make_pair(ii,nbrIdx)]=cost;
      if(mn<0 || mn>cost){
	mn=cost;
      }
      if(mx<0 || mx<cost){
	mx=cost;
      }
    }
  }

  float scale = mx-mn;

  std::map<std::pair< int , int> , float >::iterator it;
  for(it = smoothc.begin();it!=smoothc.end();it++){
    float cost = it->second;
    cost = (cost - mn)/scale;
    it->second = cost;
  }
}

struct DataCostVector
  :public GCoptimization::DataCostFunctor{
public:
  DataCostVector(std::vector<std::vector<float> > * _v):v(_v){}
  inline GCoptimization::EnergyTermType
  compute(GCoptimization::SiteID s,
	  GCoptimization::LabelID l){return (int)((*v)[s][l]);}
  std::vector<std::vector<float> > * v;
};

struct SmoothCostMap
  :public GCoptimization::SmoothCostFunctor{
public:
  SmoothCostMap(std::map<std::pair<int, int>,float > *_smoothc):smoothc(_smoothc){}
  inline GCoptimization::EnergyTermType compute
  (GCoptimization::SiteID s1,
   GCoptimization::SiteID s2,
   GCoptimization::LabelID l1,
   GCoptimization::LabelID l2){
    if(l1==l2){return 0;}
    int is1=s1,is2=(int)s2;
    if(is1>is2){
      SWAP(is1,is2);
    }
    return (int)((*smoothc)[std::make_pair(is1,is2)]);
  }
  std::map<std::pair<int, int>,float > *smoothc;
};

void scale(std::vector<std::vector<float> >&a, float scale)
{
  for(unsigned int ii=0;ii<a.size();ii++){
    for(unsigned int jj=0;jj<a[ii].size();jj++){
      a[ii][jj]*=scale;
    }
  }
}

void scale(std::map<std::pair<int,int>,float >&a, float scale)
{
  std::map<std::pair<int,int>,float >::iterator it;
  for(it = a.begin();it!=a.end();it++){
    it->second *= scale;
  }
}

#include <list>
void bfs(Mesh & m)
{
  std::vector<bool> processed (m.t.size(),0);
  std::vector<int> newlabel(m.t.size(),0);
  size_t cnt=0;
  int curidx=0;
  int label =0 ;
  std::list<int>q;
  while(cnt<m.t.size()){
    while(processed[curidx]){
      curidx++;
    }
    q.push_back(curidx);
    processed[curidx]=1;
    cnt++;
    while(!q.empty()){
      int idx = q.front();
      q.pop_front();
      newlabel[idx]=label;
      for(size_t nbr=0; nbr<m.adjMat[idx].size();nbr++){
	int nbrIdx = m.adjMat[idx][nbr];
	if(processed[nbrIdx]){
	  continue;
	}
	if(m.t[idx].label == m.t[nbrIdx].label){
	  q.push_back(nbrIdx);
	  processed[nbrIdx]=1;
	  cnt++;
	}
      }
    }
    label++;
  }

  m.nLabel=label;
  m.assign_color();
  for(size_t ii=0;ii<newlabel.size();ii++){
    m.t[ii].label=newlabel[ii];
  }
}

void runMincut(Mesh &mesh)
{
  std::vector<std::vector<float> >datac;
  std::map<std::pair<int, int> , float >smoothc;
  std::vector<Plane> plane;
  int NITER = 5;

  GCoptimizationGeneralGraph * gc = new GCoptimizationGeneralGraph(mesh.t.size(),minc_nlabel);

  DataCostVector datacFn(&datac);
  SmoothCostMap smoothFn(&smoothc);
  gc->setDataCostFunctor(&datacFn);
  gc->setSmoothCostFunctor(&smoothFn);
  mesh.adjlist();
  for(size_t ii=0;ii<mesh.t.size();ii++){
    for(size_t nbr=0;nbr<mesh.adjMat[ii].size();nbr++){
      gc->setNeighbors(ii,mesh.adjMat[ii][nbr]);
    }
  }
  for (int iter = 0;iter<NITER;iter++){
    printf("%d\n",iter);
    data_cost(mesh, minc_nlabel, plane, datac);
    smooth_cost(mesh,smoothc);
    scale(datac, dataCostW);
    scale(smoothc, smoothW);
    gc->expansion();
    for(size_t ii =0 ;ii<mesh.t.size();ii++){
      mesh.t[ii].label= gc->whatLabel(ii);
    }
  }
  bfs(mesh);
  delete gc;
}

