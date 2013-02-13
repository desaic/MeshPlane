#include "mincut.hpp"
#include <vector>
#include <stdlib.h>
#include "mesh.hpp"
#include <map>
#define SWAP(a,b)  (a)^=(b);(b)^=(a);(a)^=(b);
int minc_nlabel=50;
float dataCostW=100;
float smoothW=10;
float distw = 10;
real_t saliency_weight=1;
float L1n(Vec3f v)
{
  return std::abs(v[0]+v[1]+v[2]);
}
real_t mcdistance( Plane & p, Trig &t)
{

  float plane_d = dot(p.n,p.c);
  float d = dot(t.c,p.n);
  real_t  cost = distw * (plane_d-d);
  cost += L1n(t.n- p.n);
      cost /= (1+distw);
      cost*=t.A;
  return cost;
}
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
      real_t cost = mcdistance(plane[jj],m.t[ii]);
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

float mcdistance( Trig & a, Trig &b)
{
  float d = L1n(a.n-b.n);
  return d;
}
/**@param smoothc smoothc[i][j] is the smoothness cost for
   vertices i and j. i < j.
 */
void smooth_cost(Mesh& m,
		 std::map<EdgeId , float > & smoothc)
{
  float mn = -1;
  float mx = -1;

  for(unsigned int ii=0;ii< m.t.size();ii++){
    for(unsigned int nbr=0;nbr<m.adjMat[ii].size();nbr++ ){
      unsigned int nbrIdx = m.adjMat[ii][nbr];
      if(nbrIdx<ii){
	      continue;
      }
      EdgeId eid(ii,nbrIdx);
      float cost = L1n( m.t[ii].n - m.t[nbrIdx].n);//(area1+area2)*
      cost = 4-cost;//1/(1+cost);

      if(m.saliency.find(eid) != m.saliency.end()){
       real_t salw=(m.saliency[eid]);
        cost+=salw*saliency_weight;
      }
      if(m.usr_weit.find(eid) != m.usr_weit.end()){
        real_t usrw=(m.usr_weit[eid]);
        cost+=usrw*saliency_weight;
      }

      cost *=(m.t[ii].A+m.t[nbrIdx].A)/2;

      smoothc[eid]=cost;
      if(mn<0 || mn>cost){
	      mn=cost;
      }
      if(mx<0 || mx<cost){
	      mx=cost;
      }
    }
  }

  float scale = mx-mn;

  std::map<EdgeId , float >::iterator it;
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
  SmoothCostMap(std::map<EdgeId,float > *_smoothc):smoothc(_smoothc){}
  inline GCoptimization::EnergyTermType compute
  (GCoptimization::SiteID s1,
   GCoptimization::SiteID s2,
   GCoptimization::LabelID l1,
   GCoptimization::LabelID l2){
    if(l1==l2){return 0;}
    if(smoothc->find(EdgeId(s1,s2))==smoothc->end()){
      printf("not an edge\n");
      return 0;
    }
    return (int)((*smoothc)[EdgeId(s1,s2)]);
  }
  std::map<EdgeId,float > *smoothc;
};

void scale(std::vector<std::vector<float> >&a, float scale)
{
  for(unsigned int ii=0;ii<a.size();ii++){
    for(unsigned int jj=0;jj<a[ii].size();jj++){
      a[ii][jj]*=scale;
    }
  }
}

void scale(std::map<EdgeId,float >&a, float scale)
{
  std::map<EdgeId,float >::iterator it;
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
int MC_ITER=5;

void runMincut(Mesh &mesh)
{
  std::vector<std::vector<float> >datac;
  std::map<EdgeId , float >smoothc;
  std::vector<Plane> plane;

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
  for (int iter = 0;iter<MC_ITER;iter++){
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

