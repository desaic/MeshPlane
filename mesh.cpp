#include "mesh.hpp"
#include <fstream>
#include <iostream>
#include <algorithm>
#ifdef _WIN32
#include <windows.h>
#endif
#include <GL/gl.h>
#include <cstdlib>
#include <GL/glu.h>
#include <utility>
#include "mesh_query.h"
#include <map>
#include <sstream>
struct Edge {
  std::pair<int,int >e;
  Edge(int a=0, int b=0 ) {
    if(a>b) {
      int tmp =a;
      a=b;
      b=tmp;
    }
    e.first=a;
    e.second=b;
  }
  bool operator ==(const Edge & _e) {
    return e==_e.e;
  }
};
class EdgeCmp {
public:
  bool operator()(const Edge & a, const Edge & b) {
    return a.e<b.e;
  }
};

void Mesh::save(const char * filename)
{
  std::ofstream out;
  out.open(filename);
  out<<v.size()<<"\n"<<t.size()<<"\n";
  for(size_t ii=0; ii<v.size(); ii++) {
    for(int jj=0; jj<3; jj++) {
      out<<v[ii][jj]<<" ";
    }
    out<<"\n";
  }
  for(size_t ii=0; ii<t.size(); ii++) {
    out<<"3";
    for(int jj=0; jj<3; jj++) {
      out<<" "<<t[ii][jj];
    }
    out<<" "<<t[ii].label;
    out<<"\n";
  }
  out.close();
}

bool Mesh::self_intersect()
{ double * vertex = new double[3*v.size()];
  int * trig = new int[t.size()*3];
  int idx=0;
  for(size_t ii=0; ii<v.size(); ii++) {
    for(int jj=0; jj<3; jj++) {
      vertex[idx]=v[ii][jj];
      idx++;
    }
  }
  idx=0;
  for(size_t ii=0; ii<t.size(); ii++) {
    for(int jj=0; jj<3; jj++) {
      trig[idx]=t[ii][jj];
      idx++;
    }
  }
  MeshObject * mo = construct_mesh_object( v.size() , vertex , t.size(), trig);
  std::map<Edge, bool, EdgeCmp> m;
  double point0[3];
  double point1[3];
  int triangle_index;
  bool ret = false;
  double ss, tt, aa, bb, cc;
  for(size_t ii=0; ii<t.size(); ii++) {
    for(int jj=0; jj<3; jj++) {
      int a = t[ii][jj];
      int b = t[ii][(jj+1)%3];
      for(int kk=0; kk<3; kk++) {
        point0[kk]=v[a][kk];
        point1[kk]=v[b][kk];
      }
      Edge e(a,b);
      if(m[e]) {
        continue;
      }
      if(segment_intersects_mesh(point0,point1, mo, &triangle_index,&ss,&tt,&aa,&bb,&cc)) {
        printf("%d %d int trig %lu intersects trig %d\n",a,b,ii, triangle_index);
        bad[ii]=true;
        bad[triangle_index]=true;
        ret=true;
      }
    }
  }
  destroy_mesh_object(mo);
  return ret;
}

void Mesh::get_normal_center()
{
  for (unsigned int ii=0; ii<t.size(); ii++) {
    Trig & tt= t[ii];
    Vec3 b=v[tt[2]]-v[tt[0]];
    tt.n=(v[tt[1]]-v[tt[0]]).cross(b);
    tt.n/=tt.n.norm();
    for (int ii=0; ii<3; ii++) {
      tt.c+=v[tt[ii]];
    }
    tt.c/=3;
  }
}

void findEdge(int va, int vb, std::vector<Trig>& t,
              std::vector<std::vector< int > > &vtlist, int * wings)
{
  int cnt=0;
  for(size_t ii=0; ii<vtlist[va].size(); ii++) {
    int tIdx = vtlist[va][ii];
    for(int jj=0; jj<3; jj++) {
      if(t[tIdx][jj]==vb) {
        wings[cnt]=tIdx;
        cnt++;
        if(cnt==2) {
          return;
        }
      }
    }
  }
}

#include <set>
void Mesh::save_plane(const char * filename)
{
  std::vector<bool>processed(t.size());
  std::vector<std::vector< int > > vtlist(v.size());
  std::vector<std::set<int> > vlabel(v.size());
// std::vector< std::vector<std::vector<int> > > lines();
  lines.resize(nLabel);
//  std::vector<Plane> planes;
  for(size_t ii=0; ii<t.size(); ii++) {
    for(int jj=0; jj<3; jj++) {
      int vidx=t[ii][jj];
      vtlist[vidx].push_back(ii);
      vlabel[vidx].insert(t[ii].label);
    }
  }
  adjlist();
  get_normal_center();
  randcenter(*this, planes, nLabel);
  int wings[2];
  for(size_t ii=0; ii<t.size(); ii++) {
    if(processed[ii]) {
      continue;
    }
    processed[ii]=true;
    int boundary= false;
    int prev=t[ii][0], cur;
    int label = t[ii].label;

    for(size_t jj=0; jj<adjMat[ii].size();jj++) {
      int nbrIdx = adjMat [ii][jj];
      if(t[nbrIdx].label!=label) {
        boundary= true;
        for(int kk=0;kk<3;kk++){
          for(int ll=0;ll<3;ll++){
            if(t[ii][kk]==t[nbrIdx][ll]){
              prev = t[ii][kk];
              cur=prev;
              break;
            }
          }

        }
        break;
      }
    }
    if(!boundary) {
      continue;
    }

    std::vector<int> vertlist;
    int next = cur;
    int v0=cur;

    while(1) {
      if(vlabel[cur].size()>2) {
        vertlist.push_back(cur);
      }
      bool foundBd = false;
      for(size_t jj=0; jj<vtlist[cur].size(); jj++) {
        int nbrTidx = vtlist[cur][jj];
        if(t[nbrTidx].label!=label) {
          continue;
        }
        for(int kk=0; kk<3; kk++) {
          int nbrVidx = t[nbrTidx][kk];
          if(nbrVidx==cur || nbrVidx == prev
              || vlabel[nbrVidx].size()<2) {
            continue;
          }

          findEdge(cur,nbrVidx,t,vtlist,wings);
          if(t[wings[0]].label==t[wings[1]].label) {
            continue;
          }
          foundBd=true;
          next = nbrVidx;
          processed[nbrTidx]=true;
          break;
        }
        if(foundBd) {
          break;
        }
      }
      if( (!foundBd) || next==v0) {
        break;
      }
      prev = cur;
      cur = next;
    }
    if(vertlist.size()>2) {
      lines[label].push_back(vertlist);
    }
  }

  std::ofstream out;
  out.open(filename);
  out<<lines.size()<<"\n";
  for(size_t ii=0; ii<lines.size(); ii++) {
    out<<lines[ii].size()<<"\n";
    if(lines[ii].size()<1){
      out<<"\n";
      continue;
    }
    Vec3 v0=v[lines[ii][0][0]];

    Vec3 ax, ay;
    Vec3 n=planes[ii].n;
    Vec3 arbit(1,0,0);
    if(std::abs(n[0])>0.9) {
      arbit=Vec3(0,1,0);
    }
    ax = arbit.cross(n);
    ax/=ax.norm();
    ay=n.cross(ax);
    ay/=ay.norm();

    for(size_t jj=0; jj<lines[ii].size(); jj++) {
      out<<lines[ii][jj].size()<<"\n";
      for(size_t kk=0; kk<lines[ii][jj].size(); kk++) {
        Vec3 d = v[lines[ii][jj][kk]]-v0;
        float x = (d.dot(ax));
        float y = (d.dot(ay));
        out<<x<<","<<y<<" ";
      }
      out<<"\n";
    }
    out<<"\n";
  }
  out.close();
}

Mesh::Mesh(const char * filename, int _nLabel)
  :nLabel(_nLabel),highlight(0)
{
  std::ifstream f ;
  f.open(filename);
  if(!f.is_open()) {
    std::cout<<"cannot open "<<filename<<"\n";
    return;
  }
  int nvert, ntrig;
  f>>nvert;
  f>>ntrig;
  v.resize(nvert);
  t.resize(ntrig);
  for (int ii =0; ii<nvert; ii++) {
    for (int jj=0; jj<3; jj++) {
      f>>v[ii][jj];
    }
  }
  std::string line;
  std::getline(f,line);
  for (int ii =0; ii<ntrig; ii++) {
    int nidx;
    std::getline(f,line);
    std::stringstream ss(line);
    ss>>nidx;
    for (int jj=0; jj<3; jj++) {
      ss>>t[ii][jj];
    }
    if(ss.eof()) {
      t[ii].label=0;
    } else {
      ss>>t[ii].label;
    }
    nLabel=std::max(t[ii].label,nLabel);
  }
  nLabel++;

  real_t mn[3]= {1,1,1};
  real_t mx[3]= {-1,-1,-1};

  //scale and translate to [0 , 1]
  for (unsigned int dim = 0; dim<3; dim++) {
    for( int ii=0; ii<nvert; ii++) {
      mn [dim]= std::min(v[ii][dim],mn[dim]);
      mx[dim] = std::max(v[ii][dim],mx[dim]);
    }
    real_t translate = -0.5*(mx[dim]+mn[dim]);
    //  translate = -mn[dim];
    for(int ii=0; ii<nvert; ii++) {
      v[ii][dim]=(v[ii][dim]+translate);
    }
  }

  real_t scale = 1/(mx[0]-mn[0]);
  for(unsigned int dim=1; dim<3; dim++) {
    scale=std::min(1/(mx[dim]-mn[dim]),scale);
  }

  for(int ii=0; ii<nvert; ii++) {
    for (unsigned int dim = 0; dim<3; dim++) {
      v[ii][dim]=v[ii][dim]*scale;
    }
  }
  assign_color();
  compute_norm();
  v0=v;
}

void Mesh::assign_color()
{
  int nc0=color.size();
  color.resize(nLabel);
  for (int ii=nc0; ii<nLabel; ii++) {
    for (int jj=0; jj<3; jj++) {
      color[ii][jj]=rand()/(float)(RAND_MAX);
    }
  }
}

void Mesh::compute_norm()
{
  n.resize(v.size());
  for(unsigned int ii=0; ii<t.size(); ii++) {
    Vec3 a = v[t[ii][0]] - v[t[ii][1]];
    Vec3 b = v[t[ii][2]] - v[t[ii][0]];
    b=a.cross(b);
    for(int jj=0; jj<3; jj++) {
      n[t[ii][jj]]+=b;
    }
  }
  for(unsigned int ii=0; ii<v.size(); ii++) {
    n[ii]/= n[ii].norm();
  }
}
void Mesh::draw(std::vector<Vec3>&v)
{
  //	glDisable(GL_LIGHTING);

  glBegin(GL_TRIANGLES);
  GLfloat specular[4]= {0.51f,0.51f,0.51f,1.0f};
  GLfloat ambient[4]= {0.1f,0.1f,0.1f,1.0f};

  glMaterialfv(GL_FRONT,GL_SPECULAR,specular);
  glMaterialfv(GL_FRONT,GL_AMBIENT,ambient);
  GLfloat s=10;
  glMaterialfv(GL_FRONT_AND_BACK,GL_SHININESS,&s);
  for(unsigned int ii=0; ii<t.size(); ii++) {
    int l = t[ii].label;
    GLfloat diffuse[4]= {color[l][0],color[l][1],color[l][2],1.0f};

    glColor3f(color[l][0],color[l][1],color[l][2]);

    Vec3 a = v[t[ii][0]] - v[t[ii][1]];
    Vec3 b = v[t[ii][2]] - v[t[ii][0]];
    b=a.cross(b);
    b= b/b.norm();

    glMaterialfv(GL_FRONT,GL_DIFFUSE,diffuse);
    //glNormal3f(n[t[ii][0]][0],n[t[ii][0]][1],n[t[ii][0]][2]);
    glNormal3f(b.x[0],b.x[1],b.x[2]);


    glVertex3f(v[t[ii][0]][0],v[t[ii][0]][1],v[t[ii][0]][2]);
    //glNormal3f(n[t[ii][1]][0],n[t[ii][1]][1],n[t[ii][1]][2]);
    glVertex3f(v[t[ii][1]][0],v[t[ii][1]][1],v[t[ii][1]][2]);
    //glNormal3f(n[t[ii][2]][0],n[t[ii][2]][1],n[t[ii][2]][2]);
    glVertex3f(v[t[ii][2]][0],v[t[ii][2]][1],v[t[ii][2]][2]);

  }
  glEnd();

  /*  std::map<int, bool>::iterator it;
    for(it = bad.begin();it!=bad.end();it++){
      int ii = it->first;
      glVertex3f(2*v[t[ii][0]][0],2*v[t[ii][0]][1],2*v[t[ii][0]][2]);
      glVertex3f(2*v[t[ii][1]][0],2*v[t[ii][1]][1],2*v[t[ii][1]][2]);
      glVertex3f(2*v[t[ii][2]][0],2*v[t[ii][2]][1],2*v[t[ii][2]][2]);

    }
    //glEnable(GL_LIGHTING);
    glBegin(GL_LINES);
    for(it = bad.begin();it!=bad.end();it++){
      int ii = it->first;
      glVertex3f(v[t[ii][0]][0],v[t[ii][0]][1],v[t[ii][0]][2]);
      glVertex3f(2*v[t[ii][0]][0],2*v[t[ii][0]][1],2*v[t[ii][0]][2]);

      glVertex3f(v[t[ii][1]][0],v[t[ii][1]][1],v[t[ii][1]][2]);
      glVertex3f(2*v[t[ii][1]][0],2*v[t[ii][1]][1],2*v[t[ii][1]][2]);

      glVertex3f(v[t[ii][2]][0],v[t[ii][2]][1],v[t[ii][2]][2]);
      glVertex3f(2*v[t[ii][2]][0],2*v[t[ii][2]][1],2*v[t[ii][2]][2]);

    }
    glEnd();
  */
}
void Mesh::drawLines()
{
  glDisable(GL_LIGHTING);
  glBegin(GL_LINES);
  glColor3f(1.0,0.5,0.5);
  for(size_t ii=0; ii<lines.size(); ii++) {

    for(size_t jj=0; jj<lines[ii].size(); jj++) {
      for(size_t kk=0; kk<lines[ii][jj].size(); kk++) {
        size_t kk1 = kk+1;
        if(kk1==lines[ii][jj].size()) {
          kk1=0;
        }
        int v0=lines[ii][jj][kk];
        int v1=lines[ii][jj][kk1];
        glVertex3f(v[v0][0],v[v0][1],v[v0][2]);
        glVertex3f(v[v1][0],v[v1][1],v[v1][2]);

        if((int)ii==highlight && planes.size()>ii){
        Vec3 v0hl=v[v0]+planes[ii].n;
        Vec3 v1hl=v[v1]+planes[ii].n;
        glVertex3f(v[v0][0],v[v0][1],v[v0][2]);
        glVertex3f(v0hl[0],v0hl[1],v0hl[2]);
        glVertex3f(v0hl[0],v0hl[1],v0hl[2]);
        glVertex3f(v1hl[0],v1hl[1],v1hl[2]);

        }
      }
    }
  }

  glEnd();
  glEnable(GL_LIGHTING);
}
bool is_nbr(Trig & a, Trig&b, int vert)
{
  for (int ii=0; ii<3; ii++) {

    int va=a[ii];
    if(va<=vert) {
      continue;
    }

    for (unsigned int jj=0; jj<3; jj++) {
      int vb=b[jj];
      if(vb<=vert) {
        continue;
      }
      if(va==vb) {
        return true;
      }
    }
  }
  return false;
}


void Mesh::adjlist()
{
  if(adjMat.size()==t.size()) {
    return;
  }
  std::vector<std::vector<int> >trigList;
  trigList.resize(v.size());
  for (unsigned int ii=0; ii<t.size(); ii++) {
    for (unsigned int jj=0; jj<3; jj++) {
      int vidx=t[ii][jj];
      trigList[vidx].push_back(ii);
    }
  }
  adjMat.resize(t.size());
  for (unsigned int ii=0; ii<v.size(); ii++) {
    int n_nbr=trigList[ii].size();
    for (int jj=0; jj<n_nbr; jj++) {
      int tj=trigList[ii][jj];
      for (int kk=(jj+1); kk<n_nbr; kk++) {
        int tk=trigList[ii][kk];
        if(is_nbr(t[tj],t[tk],ii)) {
          adjMat[tj].push_back(tk);
          adjMat[tk].push_back(tj);
        }

      }
    }
  }
}
#include "randomc.h"
//pick random triangles as centers of clusters
void randcenter(Mesh & m,std::vector<Plane>&plane, int nLabel)
{
  std::vector<int>count;
  plane.resize(nLabel);
  count.resize(m.t.size());
  for (unsigned int ii=0; ii<m.t.size(); ii++) {
    size_t ll = m.t[ii].label;
    if(ll>=plane.size()) {
      continue;
    }
    plane[ll].n += m.t[ii].n;
    plane[ll].c += m.t[ii].c;
    count[ll]++;
  }

  for(int ii=0; ii<nLabel; ii++) {
    if(count[ii]>0) {
      plane[ii].n/=plane[ii].n.norm();
      plane[ii].c/=count[ii];
    } else {
      int r = rndg.IRandom(0, m.v.size()-1 );
      plane[ii].n = m.t[r].n;
      plane[ii].c = m.t[r].c;
    }
  }
}
