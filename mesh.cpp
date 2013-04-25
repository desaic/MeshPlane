#include "mesh.hpp"
#include "voxel.hpp"

#include <fstream>
#include <iostream>
#include <algorithm>
#ifdef _WIN32
#include <windows.h>
#endif
#include <GL/gl.h>
#include <GL/glu.h>
#include <stdio.h>
#include <cstdlib>
#include <utility>
#include "mesh_query.h"
#include <map>
#include <sstream>
#include <string.h>
#include "imageio.h"
#include "util.h"

bool contains(int * a ,  int x, size_t size);

bool Mesh::isNbr(int ta, int tb)
{
  for(int ii = 0; ii<3;ii++){
    if(contains(t[ta].x,t[tb][ii],3)){
      return true;
    }
  }
  return false;
}

void BBoxUnion(const std::vector<Vec3f >& v, Vec3f & mn, Vec3f & mx)
{
  for(unsigned int ii = 0;ii<v.size();ii++){
    for(int dim = 0 ; dim<3;dim++){
      if(v[ii][dim]<mn[dim]){
        mn[dim] = v[ii][dim];
      }
      if(v[ii][dim]>mx[dim]){
        mx[dim] = v[ii][dim];
      }
    }
  }
}

void BBox(const std::vector<Vec3f >& v, Vec3f & mn, Vec3f & mx)
{
  mn = v[0];
  mx = v[0];
  for(unsigned int ii = 1 ;ii<v.size();ii++){
    for(int dim = 0 ; dim<3;dim++){
      if(v[ii][dim]<mn[dim]){
        mn[dim] = v[ii][dim];
      }
      if(v[ii][dim]>mx[dim]){
        mx[dim] = v[ii][dim];
      }
    }
  }
}

void BBox(const Mesh & m, Vec3f & mn, Vec3f & mx)
{
  BBox(m.v,mn,mx);
}

bool contains(int * a ,  int x, size_t size)
{
  for(size_t ii=0; ii<size; ii++) {
    if(a[ii]==x) {
      return true;
    }
  }
  return false;
}
struct LabeledEdge {
  std::pair<int,int >e;
  int label;
  int tidx;
  LabeledEdge(int a=0, int b=0 , int _label=0, int _tidx=0) :
    label(_label), tidx(_tidx) {
    if(a>b) {
      int tmp =a;
      a=b;
      b=tmp;
    }
    e.first=a;
    e.second=b;
  }
  bool operator<(const LabeledEdge& _e)const {
    if(e < _e.e) {
      return true;
    }
    if(e>_e.e) {
      return false;
    }
    if(label<_e.label) {
      return true;
    }
    return false;
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
  std::map<EdgeId, bool> m;
  double point0[3];
  double point1[3];
  int triangle_index;
  bool ret = false;
  double ss, tt, aa, bb, cc;
  bad.clear();
  for(size_t ii=0; ii<t.size(); ii++) {
    for(int jj=0; jj<3; jj++) {
      int a = t[ii][jj];
      int b = t[ii][(jj+1)%3];
      for(int kk=0; kk<3; kk++) {
        point0[kk]=v[a][kk];
        point1[kk]=v[b][kk];
      }
      EdgeId e(a,b);
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
  real_t max_area=0;
  for (unsigned int ii=0; ii<t.size(); ii++) {
    Trig & tt= t[ii];
    Vec3f b=v[tt[2]]-v[tt[0]];
    tt.n=cross(v[tt[1]]-v[tt[0]],b);
    tt.A=mag(tt.n);
	//if(tt.A<0.000001){
//std::cout<<"bad trig\n";
	//}
    if(tt.A>max_area){
      max_area=tt.A;
    }
    tt.n/=mag(tt.n);
    tt.c=Vec3f(0,0,0);
    for (int ii=0; ii<3; ii++) {
      tt.c+=v[tt[ii]];
    }
    tt.c/=3;
  }
  //scale area so that maxmum area is 1
  for(size_t ii=0;ii<t.size();ii++){
    t[ii].A/=max_area;
  }
}

int findEdge(int va, int vb, std::vector<Trig>& t,
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
          return cnt ;
        }
      }
    }
  }
  return cnt;
}
#include <pthread.h>
pthread_mutex_t meshm = PTHREAD_MUTEX_INITIALIZER;
#include <set>

bool isOnEdge(int vidx, std::vector<std::vector<int > > & adjMat,
              std::vector<std::vector<int > > & vtlist,
              std::vector<Trig > & t)
{
  int prevT=vtlist[vidx][0];
  for(size_t ii=1; ii<vtlist[vidx].size(); ii++) {
    int tidx = vtlist[vidx][ii];
    bool found = false;
    int nbrTidx;
    for(size_t jj=0; jj<adjMat[tidx].size(); jj++) {
      nbrTidx=adjMat[tidx][jj];
      if(nbrTidx==prevT) {
        continue;
      }
      if(contains(t[nbrTidx].x, vidx,3)) {
        found = true;
        break;
      }
    }
    if(!found) {
      return true;
    }
    prevT=nbrTidx;
  }
  return false;
}

bool isOnEdge(const LabeledEdge & le, std::vector<std::vector<int > > & adjMat,
              std::vector<Trig > & t)
{
  int jj=le.e.first;
  int jj0=le.e.second;
  bool foundNbr = false;
  int tidx=le.tidx;

//find neighboring triangle that share the same edge
  int nbrIdx=0;
  for(size_t kk=0; kk<adjMat[tidx].size(); kk++) {
    nbrIdx=adjMat[tidx][kk];
    if(contains(t[nbrIdx].x, jj,3)
        &&contains(t[nbrIdx].x, jj0,3)) {
      foundNbr = true;
      break;
    }
  }
  return !foundNbr;
}


void Mesh::fix_inner_cluster()
{
  std::vector<bool>processed(t.size());
  //all clusters inside another single cluster is merged into the outer cluster
  for(size_t ii=0; ii<t.size(); ii++) {
    if(processed[ii]) {
      continue;
    }
    std::map<int,bool>nbrSet;
    std::vector<int> que;
    size_t front = 0;
    que.push_back(ii);
    int label=t[ii].label;
    while(front<que.size()) {
      int idx=que[front];
      front++;
      processed[idx]=true;
      for(size_t jj=0; jj<adjMat[idx].size(); jj++) {
        int nbrIdx = adjMat[idx][jj];
        if(t[nbrIdx].label!=label) {
          nbrSet[t[nbrIdx].label]=true;
          continue;
        }
        if(processed[nbrIdx]) {
          continue;
        }
        processed[nbrIdx]=true;
        que.push_back(nbrIdx);
      }
    }
    if(nbrSet.size()==1) {
      label = nbrSet.begin()->first;
      for(size_t jj=0; jj<que.size(); jj++) {
        t[que[jj]].label=label;
      }
    }
  }
}

void Mesh::compute_plane()
{
  std::vector<std::vector< int > > vtlist(v.size());
  std::vector<std::set<int> > vlabel(v.size());
  std::set<LabeledEdge> edgeVisited;
  std::set<LabeledEdge> edges;
  std::vector< std::vector<std::vector<int> > > local_lines(nLabel);
  for(int ii=0; ii<nLabel; ii++) {
    local_lines[ii].resize(0);
  }

  adjlist();
  fix_inner_cluster();
  for(size_t ii=0; ii<t.size(); ii++) {
    int jj0=2;
    for(int jj=0; jj<3; jj++) {
      int vidx=t[ii][jj];
      vtlist[vidx].push_back(ii);
      vlabel[vidx].insert(t[ii].label);
      //find neighboring triangle that share the same edge
      bool foundNbr = false;
      bool foundDifNbr=false;
      int nbrIdx=0;
      for(size_t kk=0; kk<adjMat[ii].size(); kk++) {
        nbrIdx=adjMat[ii][kk];
        if(contains(t[nbrIdx].x, t[ii][jj],3)
            &&contains(t[nbrIdx].x, t[ii][jj0],3)) {
          foundNbr = true;
          if(t[nbrIdx].label!=t[ii].label) {
            foundDifNbr=true;
            break;
          }
        }
      }
      if(!foundNbr
          ||foundDifNbr) {
        edges.insert(LabeledEdge(t[ii][jj],t[ii][jj0],t[ii].label,ii));
      }
      jj0=jj;
    }
  }
  get_normal_center();
  get_plane(*this, planes);
  std::vector<bool>processed(t.size());
  for(size_t ii=0; ii<t.size(); ii++) {
    if(processed[ii]){
      continue;
    }
    processed[ii]=true;
    int jj0=2;

    for(int jj=0; jj<3; jj++) {
      LabeledEdge le(t[ii][jj],t[ii][jj0],t[ii].label,ii);
      jj0=jj;
      if(edges.find(le)==edges.end()) {
        continue;
      }
      std::vector<int> vertlist;
      int label=le.label;
      int cur=le.e.first;
      int prev=le.e.second;
      if(vlabel[prev].size()>2 || isOnEdge(le,adjMat, t) ) {
          vertlist.push_back(prev);
      }
      edgeVisited.insert(le);
      int next = cur;
      int tIdx=0;
      while(1) {
        if(vlabel[cur].size()>2
           &&(vertlist.size()==0||vertlist[0]!=cur)){
          vertlist.push_back(cur);
        }
        if(isOnEdge(le,adjMat,t)){
          if(vertlist.size()==0||vertlist[vertlist.size()-1]!=prev){
            vertlist.push_back(prev);
          }
          if(vertlist.size()==0||vertlist[0]!=cur){
            vertlist.push_back(cur);
          }
        }
        bool foundBd = false;
        for(size_t jj=0; jj<vtlist[cur].size(); jj++) {
          int nbrTidx = vtlist[cur][jj];
          if(t[nbrTidx].label!=label) {
            continue;
          }
          for(int kk=0; kk<3; kk++) {
            int nbrVidx = t[nbrTidx][kk];
            if(nbrVidx==cur || nbrVidx == prev) {
              continue;
            }
            le =LabeledEdge(cur, nbrVidx, label,nbrTidx);
            if(edgeVisited.find(le)!=edgeVisited.end()) {
              continue;
            }
            if(edges.find(le)==edges.end()){
              continue;
            }
            foundBd=true;
            next = nbrVidx;
            processed[nbrTidx]=true;
            edgeVisited.insert(le);
            tIdx=nbrTidx;
            break;
          }
          if(foundBd) {
            break;
          }
        }
        if(!foundBd) {
          break;
        }

        prev = cur;
        cur = next;
        le =LabeledEdge(prev, cur, label,tIdx);
      }
      if(vertlist.size()>2) {
        local_lines[label].push_back(vertlist);
      }
    }
  }

  pthread_mutex_lock(&meshm);
  lines = local_lines;
  pthread_mutex_unlock(&meshm);
}

void Mesh::load_tex(const char * filename) {
  int width,height;
  unsigned char * buf = imageio_load_image(filename, &width,&height);
  if(!buf) {
    return;
  }
  tex_buf=buf;
  glEnable(GL_TEXTURE_2D);
  glGenTextures(1,&texture);
  glBindTexture( GL_TEXTURE_2D, texture );
  //glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_LUMINANCE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER,
                  GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,
                  GL_NEAREST);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA,
               width,  height,  0, GL_RGBA, GL_UNSIGNED_BYTE, buf);
  printf("gl error %d\n",glGetError());
  printf("max size %d\n",GL_MAX_TEXTURE_SIZE);
}
void Mesh::save_plane(const char * filename)
{
  std::ofstream out;
  out.open(filename);
  out<<lines.size()<<"\n";

  for(size_t ii=0; ii<lines.size(); ii++) {
    out<<lines[ii].size()<<"\n";
    if(lines[ii].size()<1) {
      out<<"\n";
      continue;
    }
    if(lines[ii][0].size()<3){
      continue;
    }
    int jj0=lines[ii][0].size()-1;
    Vec3f n=Vec3f(0,0,0);
    for(size_t jj=0;jj<lines[ii][0].size();jj++){
      for(int kk=0;kk<3;kk++){
        int kk1 = (kk+1)%3;
        int kk2 = (kk+2)%3;
        n[kk] += (  v[lines[ii][0][jj0]] [kk1]
                  - v[lines[ii][0][jj]] [kk1])*
                  ( v[lines[ii][0][jj0]][kk2]
                   +v[lines[ii][0][jj]] [kk2]);
      }
      jj0=jj;
    }
    n/=mag(n);
    if(dot(n,planes[ii].n)<0){
      n=-n;
    }
    //n=planes[ii].n;
    Vec3f v0=v[lines[ii][0][0]];
//    Vec3f v1=v[lines[ii][0][1]];
//    Vec3f v2=v[lines[ii][0][2]];

    Vec3f ax, ay;
    Vec3f arbit(1,0,0);
    if(std::abs(n[0])>0.9) {
      arbit=Vec3f(0,1,0);
    }
    ax = cross(arbit,n);
    ax/=mag(ax);
    ay=cross(n,ax);
    ay/=mag(ay);
    out<<n.v[0]<<" "<<n.v[1]<<" "<<n.v[2]<<"\n";
    //save transformation from world coordinates to plane coordinates
    //normal is z axis
    out<<ax.v[0]<<" "<<ax.v[1]<<" "<<ax.v[2]<<"\n";
    out<<ay.v[0]<<" "<<ay.v[1]<<" "<<ay.v[2]<<"\n";
    out<<v0.v[0]<<" "<<v0.v[1]<<" "<<v0.v[2]<<"\n";
    for(size_t jj=0; jj<lines[ii].size(); jj++) {
      out<<lines[ii][jj].size()<<"\n";
      for(size_t kk=0; kk<lines[ii][jj].size(); kk++) {
        Vec3f d = v[lines[ii][jj][kk]]-v0;
        float x = dot(d,(ax));
        float y = dot(d,(ay));
        out<<x<<","<<y<<" ";
      }
      out<<"\n";
      for(size_t kk=0; kk<lines[ii][jj].size(); kk++) {
        out<<lines[ii][jj][kk]<<" ";
      }
      out<<"\n";
    }
    out<<"\n";
  }
  out.close();
}

void Mesh::save_off(const char * filename)
{
  std::ofstream out;
  out.open(filename);
  out<<"OFF\n";
  out<<v.size()<<" ";
  int npoly=0;
  for(size_t ii=0; ii<lines.size(); ii++) {
    for(size_t jj=0; jj<lines[ii].size(); jj++) {
      if(lines[ii][jj].size()>2){
        npoly++;
      }  
    }
  }
  out<<npoly<<" 0\n\n";

  for(size_t ii=0;ii<v.size();ii++){
    out<< std::fixed <<v[ii][0]<<" "<<v[ii][1]<<" "<<v[ii][2]<<"\n";
  }

  for(size_t ii=0; ii<lines.size(); ii++) {
    for(size_t jj=0; jj<lines[ii].size(); jj++) {
      out<<lines[ii][jj].size()<<" ";
      if(lines[ii][jj].size()<3){
        continue;
      }
      for(size_t kk=0; kk<lines[ii][jj].size(); kk++) {
        out<<" "<<lines[ii][jj][kk];
      }
      out<<"\n";
    }
  }
  out<<"\n";
  out.close();
}

void Mesh::read_obj(std::ifstream & f)
{
  std::string line;
  std::string vTok("v");
  std::string fTok("f");
  std::string texTok("vt");
  char bslash='/',space=' ';
  std::string tok;
  while(1) {
    std::getline(f,line);
    if(f.eof()) {
      break;
    }
    if(line.size()<3) {
      continue;
    }
    if(line.at(0)=='#') {
      continue;
    }
    std::stringstream ss(line);
    ss>>tok;
    if(tok==vTok) {
      Vec3f vec;
      ss>>vec[0]>>vec[1]>>vec[2];
      v.push_back(vec);
    } else if(tok==fTok) {
      if(line.find(bslash)!=std::string::npos) {
        std::replace(line.begin(),line.end(),bslash,space);
        std::stringstream facess(line);
        Trig trig;
        facess>>tok;
        for(int ii=0; ii<3; ii++) {
          facess>>trig[ii]>>trig.texId[ii];
          trig[ii]--;
          trig.texId[ii]--;
        }
        t.push_back(trig);
      } else {
        Trig trig;
        for(int ii=0; ii<3; ii++) {
          ss>>trig[ii];
          trig[ii]--;
          trig.texId[ii]=0;
        }
        t.push_back(trig);
      }
    } else if(tok==texTok) {
      Vec3f texcoord;
      ss>>texcoord[0];
      ss>>texcoord[1];
      tex.push_back(texcoord);
    }
  }
}

void Mesh::read_ply(std::ifstream & f)
{
  std::string line;
  std::string vertLine("element vertex");
  std::string faceLine("element face");
  std::string texLine("property float s");
  std::string endHeaderLine("end_header");
  while(true) {
    std::getline(f,line);
    if(std::string::npos!=line.find(vertLine)) {
      break;
    }
  }
  std::string token;
  std::stringstream ss(line);
  ss>>token>>token;
  int nvert;
  ss>>nvert;
  bool hasTex=false;
  while(true) {
    std::getline(f,line);
    if(std::string::npos!=line.find(faceLine)) {
      break;
    }
    if(std::string::npos!=line.find(texLine)) {
      hasTex=true;
    }
  }
  std::stringstream ss1(line);
  ss1>>token>>token;
  int nface;
  ss1>>nface;
  while(true) {
    std::getline(f,line);
    if(std::string::npos!=line.find(endHeaderLine)) {
      break;
    }
  }

  v.resize(nvert);
  t.resize(nface);
  if(hasTex) {
    tex.resize(nvert);
  }
  for (int ii =0; ii<nvert; ii++) {
    for (int jj=0; jj<3; jj++) {
      f>>v[ii][jj];
    }
    if(hasTex) {
      for (int jj=0; jj<2; jj++) {
        f>>tex[ii][jj];
      }
      tex[ii][1]=1-tex[ii][1];;
    }
  }
  for (int ii =0; ii<nface; ii++) {
    int nidx;
    f>>nidx;
    t[ii].label=0;
    for (int jj=0; jj<3; jj++) {
      f>>t[ii][jj];
    }
  }
}

void Mesh::read_ply2(std::ifstream&f)
{
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
      if(ss.fail()) {
        t[ii].label=0;
      }
    }
    nLabel=std::max(t[ii].label,nLabel);
  }
  nLabel++;
}

void int2b(unsigned int x, GLubyte * b)
{
  for(int ii=0;ii<4;ii++){
    b[ii]=x&0xff;
    x=x>>8;
  }
}

unsigned int b2int(GLubyte * b)
{
  unsigned int x=0;
  for(int ii=3;ii>=0;ii--){
    x=x<<8;
    x+=b[ii];
  }
  return x;
}

void Mesh::drawCol()
{
  glDisable(GL_LIGHTING);
  glBegin(GL_TRIANGLES);
  for(unsigned int ii=0; ii<t.size(); ii++) {
    unsigned int l = t[ii].label;
    GLubyte b[4];
    int2b(l,b);
    glColor4ub(b[0],b[1],b[2],b[3]);
    glVertex3f(v[t[ii][0]][0],v[t[ii][0]][1],v[t[ii][0]][2]);
    glVertex3f(v[t[ii][1]][0],v[t[ii][1]][1],v[t[ii][1]][2]);
    glVertex3f(v[t[ii][2]][0],v[t[ii][2]][1],v[t[ii][2]][2]);

  }
  glEnd();
}

void Mesh::save_obj(const char * filename)
{
  std::ofstream out(filename);
  if(!out.good()){
    std::cout<<"cannot open output file"<<filename<<"\n";
    return;
  }
  std::string vTok("v");
  std::string fTok("f");
  std::string texTok("vt");
  char bslash='/';
  std::string tok;
  for(size_t ii=0;ii<v.size();ii++){
    out<<vTok<<" "<<v[ii][0]<<" "<<v[ii][1]<<" "<<v[ii][2]<<"\n";
  }
  if(tex.size()>0){
    for(size_t ii=0;ii<tex.size();ii++){
      out<<texTok<<" "<<tex[ii][0]<<" "<<tex[ii][1]<<"\n";
    }
    for(size_t ii=0;ii<t.size();ii++){
      out<<fTok<<" "<<t[ii][0]+1<<bslash<<t[ii].texId[0]+1<<" "
      <<t[ii][1]+1<<bslash<<t[ii].texId[1]+1<<" "
      <<t[ii][2]+1<<bslash<<t[ii].texId[2]+1<<"\n";
    }
  }else{
    for(size_t ii=0;ii<t.size();ii++){
      out<<fTok<<" "<<t[ii][0]+1<<" "<<t[ii][1]+1<<" "<<t[ii][2]+1<<"\n";
    }
  }

  out.close();
}

Mesh::Mesh(const char * filename, int _nLabel , bool _auto)
  :nLabel(_nLabel),highlight(1000),remap_tex(0),fbo(0),checkIntersect(0),
  tex_buf(0),autoscale(_auto)

{
  std::ifstream f ;
  f.open(filename);
  if(!f.is_open()) {
    std::cout<<"cannot open "<<filename<<"\n";
    return;
  }

  switch(filename[strlen(filename)-1]) {
  case 'y':
    read_ply(f);
    break;
  case 'j':
    read_obj(f);
    break;
  default:
    read_ply2(f);
  }
  assign_color();
  float mn[3]= {1,1,1};
  float mx[3]= {-1,-1,-1};

  //scale and translate to [0 , 1]
  if(autoscale){
  for (unsigned int dim = 0; dim<3; dim++) {
    for( size_t ii=0; ii<v.size(); ii++) {
      mn[dim]= std::min(v[ii][dim],mn[dim]);
      mx[dim] = std::max(v[ii][dim],mx[dim]);
    }
    float translate = -0.5*(mx[dim]+mn[dim]);
    //  translate = -mn[dim];
    for(size_t ii=0; ii<v.size(); ii++) {
      v[ii][dim]=(v[ii][dim]+translate);
    }
  }

  float scale = 1/(mx[0]-mn[0]);
  for(unsigned int dim=1; dim<3; dim++) {
    scale=std::min(1/(mx[dim]-mn[dim]),scale);
  }

  for(size_t ii=0; ii<v.size(); ii++) {
    for (unsigned int dim = 0; dim<3; dim++) {
      v[ii][dim]=v[ii][dim]*scale;
    }
  }
  }
  v0=v;
  compute_norm();

  f.close();
}

void Mesh::assign_color()
{
  color.push_back(Vec3f(0.3,0.3,0.3));
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
  for(unsigned int ii=0; ii<v.size(); ii++) {
    n[ii] = Vec3f(0,0,0);
  }
  for(unsigned int ii=0; ii<t.size(); ii++) {
    Vec3f a = v[t[ii][0]] - v[t[ii][1]];
    Vec3f b = v[t[ii][2]] - v[t[ii][0]];
    b=cross(a,b);
    for(int jj=0; jj<3; jj++) {
      if(t[ii][jj]>=(int)n.size()){
        std::cout<<"Input obj file buggy. Vertex index larger than number of vertices\n";
        std::cout<<"triangle: "<<ii<<"\n";
      }
      n[t[ii][jj]]+=b;
    }
  }
  for(unsigned int ii=0; ii<v.size(); ii++) {
    n[ii]/= mag(n[ii]);
  }
}

void Mesh::drawPlane(int k)
{
  if(k>(int)planes.size()-1) {
    k=planes.size()-1;
  }
  Vec3f nn = planes[k].n;
  Vec3f v0=planes[k].c;
  Vec3f ax(0), ay(0);
  Vec3f arbit(1,0,0);
  if(std::abs(nn[0])>0.9) {
    arbit=Vec3f(0,1,0);
  }
  ax = cross(arbit,nn);
  normalize(ax);
  ax=cross(nn,ax);
  normalize(ay);

  glDisable(GL_LIGHTING);
  glBegin(GL_TRIANGLES);
  if(tex_buf) {
    glBindTexture(GL_TEXTURE_2D,texture);
  }

  for(unsigned int ii=0; ii<t.size(); ii++) {
    int l = t[ii].label;
    if(l!=k) {
      continue;
    }
    if(tex_buf && tex.size()>0) {
      for(int jj=0; jj<3; jj++) {
        Vec3f vert=v[t[ii][jj]];
        vert-=v0;
        vert=Vec3f(dot(vert,ax),dot(vert,ay),dot(vert,nn));
        glNormal3f(0,0,1);
        glTexCoord2f(tex[t[ii].texId[jj]][0],tex[t[ii].texId[jj]][1]);
        glVertex3f(vert[0],vert[1],-1);
      }
    } else {
      glNormal3f(n[t[ii][0]][0],n[t[ii][0]][1],n[t[ii][0]][2]);
      glVertex3f(v[t[ii][0]][0],v[t[ii][0]][1],v[t[ii][0]][2]);
      glNormal3f(n[t[ii][1]][0],n[t[ii][1]][1],n[t[ii][1]][2]);
      glVertex3f(v[t[ii][1]][0],v[t[ii][1]][1],v[t[ii][1]][2]);
      glNormal3f(n[t[ii][2]][0],n[t[ii][2]][1],n[t[ii][2]][2]);
      glVertex3f(v[t[ii][2]][0],v[t[ii][2]][1],v[t[ii][2]][2]);

    }

  }
  glEnd();
}

void Mesh::draw(std::vector<Vec3f>&v)
{
//  glDisable(GL_LIGHTING);
//  glDisable(GL_TEXTURE_2D);

  glBegin(GL_TRIANGLES);
  GLfloat specular[4]= {0.51f,0.51f,0.51f,1.0f};
  GLfloat ambient[4]= {0.1f,0.1f,0.1f,1.0f};

  glMaterialfv(GL_FRONT,GL_SPECULAR,specular);
  glMaterialfv(GL_FRONT,GL_AMBIENT,ambient);
  GLfloat s=10;
  glMaterialfv(GL_FRONT,GL_SHININESS,&s);
  if(tex_buf) {
    glBindTexture(GL_TEXTURE_2D,texture);
  }

  for(unsigned int ii=0; ii<t.size(); ii++) {
   unsigned int  l = t[ii].label;
  /*  real_t sal=1;
    for(size_t jj=0;jj<adjMat[ii].size();jj++){
      EdgeId eid(ii,adjMat[ii][jj]);
      if(usr_weit.find(eid)!=usr_weit.end()){
        sal+=usr_weit[eid];
      }
    }
    sal/=3;
    if(sal<0.2){
      sal=0.2;
    }
    sal=1-sal;
    glColor3f(sal,sal,sal);

    GLfloat diffuse[4]= {sal,sal,sal,1.0f};
    glMaterialfv(GL_FRONT,GL_DIFFUSE,diffuse);
*/
    Vec3f a = v[t[ii][1]] - v[t[ii][0]];
    Vec3f b = v[t[ii][2]] - v[t[ii][0]];
    b=-cross(a,b);
    normalize(b);
    if(tex_buf && tex.size()>0) {
      glNormal3f(n[t[ii][0]][0],n[t[ii][0]][1],n[t[ii][0]][2]);
      glTexCoord2f(tex[t[ii].texId[0]][0],tex[t[ii].texId[0]][1]);
      glVertex3f(v[t[ii][0]][0],v[t[ii][0]][1],v[t[ii][0]][2]);

      glTexCoord2f(tex[t[ii].texId[1]][0],tex[t[ii].texId[1]][1]);
      glNormal3f(n[t[ii][1]][0],n[t[ii][1]][1],n[t[ii][1]][2]);
      glVertex3f(v[t[ii][1]][0],v[t[ii][1]][1],v[t[ii][1]][2]);

      glTexCoord2f(tex[t[ii].texId[2]][0],tex[t[ii].texId[2]][1]);
      glNormal3f(n[t[ii][2]][0],n[t[ii][2]][1],n[t[ii][2]][2]);
      glVertex3f(v[t[ii][2]][0],v[t[ii][2]][1],v[t[ii][2]][2]);
    } else {
      if(l<color.size()){
        GLfloat clr[4]= {(float)color[l][0],(float)color[l][1],
                         (float)color[l][2],1.0f};
        glColor3f(color[l][0],color[l][1],color[l][2]);
        glMaterialfv(GL_FRONT,GL_DIFFUSE,clr);
      }
      glNormal3f(b[0],b[1],b[2]);
      //glNormal3f(n[t[ii][0]][0],n[t[ii][0]][1],n[t[ii][0]][2]);
      glVertex3f(v[t[ii][0]][0],v[t[ii][0]][1],v[t[ii][0]][2]);
     // glNormal3f(n[t[ii][1]][0],n[t[ii][1]][1],n[t[ii][1]][2]);
      glVertex3f(v[t[ii][1]][0],v[t[ii][1]][1],v[t[ii][1]][2]);
     // glNormal3f(n[t[ii][2]][0],n[t[ii][2]][1],n[t[ii][2]][2]);
      glVertex3f(v[t[ii][2]][0],v[t[ii][2]][1],v[t[ii][2]][2]);

    }
    if(t[ii].label==highlight) {
      glVertex3f(v[t[ii][0]][0]+b.v[0],v[t[ii][0]][1]+b.v[1],v[t[ii][0]][2]+b.v[2]);

      glVertex3f(v[t[ii][1]][0]+b.v[0],v[t[ii][1]][1]+b.v[1],v[t[ii][1]][2]+b.v[2]);

      glVertex3f(v[t[ii][2]][0]+b.v[0],v[t[ii][2]][1]+b.v[1],v[t[ii][2]][2]+b.v[2]);
    }
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

void Mesh::draw_tex()
{
  if(!tex_buf || tex.size()==0) {
    return;
  }

  glDisable(GL_LIGHTING);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glTranslatef(-0.5,-0.5,0);
  glBegin(GL_TRIANGLES);
  if(tex_buf) {
    glBindTexture(GL_TEXTURE_2D,texture);
  }

  for(unsigned int ii=0; ii<t.size(); ii++) {
    if(remap_tex==0){
      glTexCoord2f(tex[t[ii].texId[0]][0],tex[t[ii].texId[0]][1]);
      glVertex3f(tex[t[ii].texId[0]][0],tex[t[ii].texId[0]][1],-1);

      glTexCoord2f(tex[t[ii].texId[1]][0],tex[t[ii].texId[1]][1]);
      glVertex3f(tex[t[ii].texId[1]][0],tex[t[ii].texId[1]][1],-1);

      glTexCoord2f(tex[t[ii].texId[2]][0],tex[t[ii].texId[2]][1]);
      glVertex3f(tex[t[ii].texId[2]][0],tex[t[ii].texId[2]][1],-1);
    }else{
      glTexCoord2f(tex[t[ii].texId[0]][0],tex[t[ii].texId[0]][1]);
      glVertex3f(remap_tex->tex[remap_tex->t[ii].texId[0]][0],
                 remap_tex->tex[remap_tex->t[ii].texId[0]][1],-1);

      glTexCoord2f(tex[t[ii].texId[1]][0],tex[t[ii].texId[1]][1]);
      glVertex3f(remap_tex->tex[remap_tex->t[ii].texId[1]][0],
                 remap_tex->tex[remap_tex->t[ii].texId[1]][1],-1);

      glTexCoord2f(tex[t[ii].texId[2]][0],tex[t[ii].texId[2]][1]);
      glVertex3f(remap_tex->tex[remap_tex->t[ii].texId[2]][0],
                 remap_tex->tex[remap_tex->t[ii].texId[2]][1],-1);
    }
  }
  glEnd();
}

void Mesh::drawLines()
{
  glDisable(GL_LIGHTING);
  glDisable(GL_TEXTURE_2D);
  glBegin(GL_LINES);
  glColor3f(0.19,0.19,0.29);
  pthread_mutex_lock(&meshm);

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

        if((int)ii==highlight && planes.size()>ii) {
          Vec3f v0hl=v[v0]+planes[ii].n;
          Vec3f v1hl=v[v1]+planes[ii].n;
          glVertex3f(v[v0][0],v[v0][1],v[v0][2]);
          glVertex3f(v0hl[0],v0hl[1],v0hl[2]);
          glVertex3f(v0hl[0],v0hl[1],v0hl[2]);
          glVertex3f(v1hl[0],v1hl[1],v1hl[2]);

        }
      }
    }
  }
  pthread_mutex_unlock(&meshm);
  glEnd();
  glEnable(GL_TEXTURE_2D);
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

void update_distance(std::vector<real_t > & dist,
                      std::vector<real_t > & cdf,
                      std::vector<Plane> & plane,
                      Mesh & m, int ll);

#include "randomc.h"
//pick random triangles as centers of clusters
void randcenter(Mesh & m,std::vector<Plane>&plane, int nLabel)
{
  std::vector<int>count;
  plane.resize(nLabel);
  count.resize(m.t.size());
  std::vector<real_t > dist(m.t.size(), 0.0f) ;
  std::vector<real_t > cdf(m.t.size(),0.0f);
  real_t sum=0;
  real_t totalA=0;
  for (unsigned int ii=0; ii<m.t.size(); ii++) {
    size_t ll = m.t[ii].label;
    if(ll>=plane.size()) {
      continue;
    }
    plane[ll].n += m.t[ii].n;
    plane[ll].c += m.t[ii].c;
    plane[ll].A+=m.t[ii].A;
    totalA+=m.t[ii].A;
    dist[ii]=mcdistance(plane[ll],m.t[ii]);
    sum+=dist[ii]*dist[ii];
    cdf[ii]=sum;
    count[ll]++;
  }

  for(int ii=0; ii<nLabel; ii++) {
    if(count[ii]>0) {
      normalize(plane[ii].n);
      plane[ii].c/=count[ii];

    }
  }
  for(int ii=0; ii<nLabel; ii++) {
    //if(count[ii]!=0 && plane[ii].A<totalA/300.0){
    //  std::cout<<"tinyplane\n";
    //}
    if(count[ii]==0 || plane[ii].A<totalA/300.0){

      int r=rndg.IRandom(0,m.t.size()-1);//sample_cdf(cdf);
      m.t[r].label=ii;
      plane[ii].n = m.t[r].n;
      plane[ii].c = m.t[r].c;
     // std::cout<<r<<" "<<ii<<" "<<plane[ii].c[0]<<" "<<plane[ii].n[0]<<"\n";
      update_distance(dist, cdf, plane,m,ii);
    }
  }
}
/**@param m assume triangle norms and centers are already computed
 */
void get_plane(Mesh & m , std::vector<Plane> & plane)
{
  plane.resize(m.nLabel);
  std::vector<float > cnt (m.nLabel,0);
  for(size_t ii=0; ii<m.t.size(); ii++) {
    if(m.t[ii].A<0.000001){
      continue;
    }
    Vec3f a = m.v[m.t[ii][1]] -  m.v[m.t[ii][0]];
    Vec3f b = m.v[m.t[ii][2]] -  m.v[m.t[ii][0]];
    Vec3f n = cross(a,b);
    float area = mag(n);
    int label = m.t[ii].label;
    cnt[label]+=area;
    plane[label].n += n;
    a=(m.t[ii].c * area);
    plane[label].c += a;
  }

  for(int ii=0; ii<m.nLabel; ii++) {
    if(cnt[ii]>0){
      plane[ii].c/=cnt[ii];
      // plane[ii].n/=cnt[ii];
      normalize(plane[ii].n);
    }
  }
}

GLcharARB * read_entire_file(const char * filename, int * len )
{
  FILE * file = fopen(filename, "r");
  if(!file){
    printf("cannot open shader %s\n", filename);
    return 0;
  }
  GLcharARB * buf=0;
  fseek(file, 0, SEEK_END);
	size_t length = ftell(file);
	fseek(file, 0, SEEK_SET);
	buf = new GLcharARB[length+1];
  length = fread( buf, 1, length, file);
  buf[length]=0;
  *len=length;
  return buf;
}

void Mesh::init_select(const char * shaderfile)
{
  glewInit();
  /*if(!glCreateShaderObjectARB){
    select_shader=0;
    return;
  }
  select_shader=glCreateShaderObjectARB(GL_FRAGMENT_SHADER);
  int len =0;
  GLcharARB * buf = read_entire_file(shaderfile, &len);
  glShaderSourceARB(select_shader, 1,  (const GLcharARB**) (&buf), 0);
  glCompileShaderARB(select_shader);

  int  ret_code = 0;
  glGetObjectParameterivARB(select_shader, GL_OBJECT_COMPILE_STATUS_ARB, &ret_code);
    if (ret_code == 0)
    {
      char error_string[512];
      glGetInfoLogARB(select_shader, sizeof(error_string), 0, error_string);
      printf("%s\n",error_string);
      select_shader=0;
    }

  delete [] buf;
  */
  int tex_wid=1280, tex_hig=720;

  glGenTextures(1, &fbot);
  glBindTexture(GL_TEXTURE_2D, fbot);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, tex_wid, tex_hig, 0,
             GL_RGBA, GL_UNSIGNED_BYTE, 0);
  glBindTexture(GL_TEXTURE_2D, 0);

  GLuint rboId;
  glGenRenderbuffers(1, &rboId);
  glBindRenderbuffer(GL_RENDERBUFFER, rboId);
  glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT,
                      tex_wid, tex_hig);
  glBindRenderbuffer(GL_RENDERBUFFER, 0);

  glGenFramebuffers(1, &fbo);
  glBindFramebuffer(GL_FRAMEBUFFER, fbo);

  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0,
                       GL_TEXTURE_2D, fbot, 0);
  // attach the renderbuffer to depth attachment point
  glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT,
                          GL_RENDERBUFFER, rboId);

  glBindFramebuffer(GL_FRAMEBUFFER, 0);
}
