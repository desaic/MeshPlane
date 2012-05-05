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
#include <string.h>
#include "imageio.h"
bool contains(int * a ,  int x, size_t size)
{
  for(size_t ii=0;ii<size;ii++){
    if(a[ii]==x){
      return true;
    }
  }
  return false;
}

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

struct LabeledEdge{
  std::pair<int,int >e;
  int label;
  LabeledEdge(int a=0, int b=0 , int _label=0) :
  label(_label){
    if(a>b) {
      int tmp =a;
      a=b;
      b=tmp;
    }
    e.first=a;
    e.second=b;
  }
  bool operator<(const LabeledEdge& _e)const{
    if(e < _e.e){
      return true;
    }
    if(e>_e.e){
      return false;
    }
    if(label<_e.label){
      return true;
    }
    return false;
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
#include <pthread.h>
pthread_mutex_t meshm = PTHREAD_MUTEX_INITIALIZER;
#include <set>

bool isOnEdge(int vidx, std::vector<std::vector<int > > & adjMat,
              std::vector<std::vector<int > > & vtlist,
              std::vector<Trig > & t)
{
  int prevT=vtlist[vidx][0];
  for(size_t ii=1;ii<vtlist[vidx].size();ii++){
    int tidx = vtlist[vidx][ii];
    bool found = false;
    int nbrTidx;
    for(size_t jj=0;jj<adjMat[tidx].size();jj++){
      nbrTidx=adjMat[tidx][jj];
      if(nbrTidx==prevT){
        continue;
      }
      if(contains(t[nbrTidx].x, vidx,3)){
        found = true;
        break;
      }
    }
    if(!found){
      return true;
    }
    prevT=nbrTidx;
  }
  return false;
}

void Mesh::fix_inner_cluster()
{
  std::vector<bool>processed(t.size());
  printf("0\n");
  //all clusters inside another single cluster is merged into the outer cluster
  for(size_t ii=0;ii<t.size();ii++){
    if(processed[ii]){
      continue;
    }
    std::map<int,bool>nbrSet;
    std::vector<int> que;
    size_t front = 0;
    que.push_back(ii);
    int label=t[ii].label;
    while(front<que.size()){
      int idx=que[front];
      front++;
      processed[idx]=true;
      for(size_t jj=0;jj<adjMat[idx].size();jj++){
        int nbrIdx = adjMat[idx][jj];
        if(t[nbrIdx].label!=label){
          nbrSet[t[nbrIdx].label]=true;
          continue;
        }
        if(processed[nbrIdx]){
          continue;
        }
        processed[nbrIdx]=true;
        que.push_back(nbrIdx);
      }
    }
    if(nbrSet.size()==1){
      label = nbrSet.begin()->first;
      for(size_t jj=0;jj<que.size();jj++){
        t[que[jj]].label=label;
      }
    }
  }
  printf("1\n");
}

void Mesh::compute_plane()
{
  std::vector<bool>processed(t.size());
  std::vector<std::vector< int > > vtlist(v.size());
  std::vector<std::set<int> > vlabel(v.size());
  std::set<LabeledEdge> edgeVisited;
  std::vector< std::vector<std::vector<int> > > local_lines(nLabel);
  for(int ii=0;ii<nLabel;ii++){
    local_lines[ii].resize(0);
  }

  adjlist();
  fix_inner_cluster();
  for(size_t ii=0; ii<t.size(); ii++) {
    for(int jj=0; jj<3; jj++) {
      int vidx=t[ii][jj];
      vtlist[vidx].push_back(ii);
      vlabel[vidx].insert(t[ii].label);
    }
  }
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
/*
    if(adjMat[ii].size()<3){
      boundary=true;
      for(int jj=0;jj<3;jj++){
        if(isOnEdge(t[ii][jj],adjMat,vtlist, t)){
          cur = t[ii][jj];
          break;
        }
      }
    }
*/
    if(!boundary) {
      continue;
    }

    std::vector<int> vertlist;
    int next = cur;
    int v0=cur;
    if(label==0){
//      std::cout<<"debug\n";
    }
    if(ii==4682){
      std::cout<<"debug\n";
    }
//  int prevTidx=ii;
    while(1) {
      if(vlabel[cur].size()>2 || isOnEdge(cur,adjMat, vtlist,t) ){
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
          LabeledEdge le (cur, nbrVidx, label);
          if(edgeVisited.find(le)!=edgeVisited.end()){
            continue;
          }
          findEdge(cur,nbrVidx,t,vtlist,wings);
          if(t[wings[0]].label==t[wings[1]].label) {
            continue;
          }
          foundBd=true;
          next = nbrVidx;
          processed[nbrTidx]=true;
          edgeVisited.insert(le);
          //prevTidx=nbrTidx;
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
      local_lines[label].push_back(vertlist);
    }
  }

  pthread_mutex_lock(&meshm);
  lines = local_lines;
  pthread_mutex_unlock(&meshm);
}
real_t Mesh::area(const Trig & trig)
{
  Vec3 n=(v[trig[1]]-v[trig[0]]).cross(v[trig[2]]-v[trig[0]]);
  return n.norm()/2;
}
void Mesh::load_tex(const char * filename){
  int width,height;
  unsigned char * buf = imageio_load_image(filename, &width,&height);
  if(!buf){
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
}
void Mesh::load_ptex(const char * filename)
{
  Ptex::String errMsg("cannot open ptx\n");
  ptx=PtexTexture::open(filename,errMsg);
}

void Mesh::save_plane(const char * filename)
{
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
    out<<n.x[0]<<" "<<n.x[1]<<" "<<n.x[2]<<"\n";
    //save transformation from world coordinates to plane coordinates
    //normal is z axis
    out<<ax.x[0]<<" "<<ax.x[1]<<" "<<ax.x[2]<<"\n";
    out<<ay.x[0]<<" "<<ay.x[1]<<" "<<ay.x[2]<<"\n";
    out<<v0.x[0]<<" "<<v0.x[1]<<" "<<v0.x[2]<<"\n";
    for(size_t jj=0; jj<lines[ii].size(); jj++) {
      out<<lines[ii][jj].size()<<"\n";
      for(size_t kk=0; kk<lines[ii][jj].size(); kk++) {
        Vec3 d = v[lines[ii][jj][kk]]-v0;
        float x = (d.dot(ax));
        float y = (d.dot(ay));
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

void Mesh::read_obj(std::ifstream & f)
{
  std::string line;
  std::string vTok("v");
  std::string fTok("f");
  std::string texTok("vt");
  char bslash='/',space=' ';
  std::string tok;
  while(1){
    std::getline(f,line);
    if(f.eof()){
      break;
    }
    if(line.size()<3){
      continue;
    }
    if(line.at(0)=='#'){
      continue;
    }
    std::stringstream ss(line);
    ss>>tok;
    if(tok==vTok){
      Vec3 vec;
      ss>>vec[0]>>vec[1]>>vec[2];
      v.push_back(vec);
    }else if(tok==fTok){
      if(line.find(bslash)!=std::string::npos){
        std::replace(line.begin(),line.end(),bslash,space);
        std::stringstream facess(line);
        Trig trig;
        facess>>tok;
        for(int ii=0;ii<3;ii++){
          facess>>trig[ii]>>trig.texId[ii];
          trig[ii]--;
          trig.texId[ii]--;
        }
        t.push_back(trig);
      }else{
        Trig trig;
        for(int ii=0;ii<3;ii++){
          ss>>trig[ii];
          trig[ii]--;
          trig.texId[ii]=0;
        }
        t.push_back(trig);
      }
    }else if(tok==texTok){
      Vec3 texcoord;
      ss>>texcoord[0];
      ss>>texcoord[1];
      tex.push_back(texcoord);
    }
  }

  color.push_back(Vec3(1,1,1));
  nLabel=1;
}

void Mesh::read_ply(std::ifstream & f)
{
  std::string line;
  std::string vertLine("element vertex");
  std::string faceLine("element face");
  std::string texLine("property float s");
  std::string endHeaderLine("end_header");
  while(true){
    std::getline(f,line);
    if(std::string::npos!=line.find(vertLine)){
      break;
    }
  }
  std::string token;
  std::stringstream ss(line);
  ss>>token>>token;
  int nvert;
  ss>>nvert;
  bool hasTex=false;
  while(true){
    std::getline(f,line);
    if(std::string::npos!=line.find(faceLine)){
      break;
    }
    if(std::string::npos!=line.find(texLine)){
      hasTex=true;
    }
  }
  std::stringstream ss1(line);
  ss1>>token>>token;
  int nface;
  ss1>>nface;
  while(true){
    std::getline(f,line);
    if(std::string::npos!=line.find(endHeaderLine)){
      break;
    }
  }

  v.resize(nvert);
  t.resize(nface);
  if(hasTex){
    tex.resize(nvert);
  }
  for (int ii =0; ii<nvert; ii++) {
    for (int jj=0; jj<3; jj++) {
      f>>v[ii][jj];
    }
    if(hasTex){
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
  color.push_back(Vec3(1,1,1));
  nLabel=1;
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
      if(ss.fail()){
        t[ii].label=0;
      }
    }
    nLabel=std::max(t[ii].label,nLabel);
  }
  nLabel++;

}
Mesh::Mesh(const char * filename, int _nLabel)
  :nLabel(_nLabel),highlight(100),ptx(0),tex_buf(0)
{
  std::ifstream f ;
  f.open(filename);
  if(!f.is_open()) {
    std::cout<<"cannot open "<<filename<<"\n";
    return;
  }

  switch(filename[strlen(filename)-1]){
  case 'y':
    read_ply(f);
    break;
  case 'j':
    read_obj(f);
    break;
  default:
    read_ply2(f);
    assign_color();
  }

  real_t mn[3]= {1,1,1};
  real_t mx[3]= {-1,-1,-1};

  //scale and translate to [0 , 1]
  for (unsigned int dim = 0; dim<3; dim++) {
    for( size_t ii=0; ii<v.size(); ii++) {
      mn [dim]= std::min(v[ii][dim],mn[dim]);
      mx[dim] = std::max(v[ii][dim],mx[dim]);
    }
    real_t translate = -0.5*(mx[dim]+mn[dim]);
    //  translate = -mn[dim];
    for(size_t ii=0; ii<v.size(); ii++) {
      v[ii][dim]=(v[ii][dim]+translate);
    }
  }

  real_t scale = 1/(mx[0]-mn[0]);
  for(unsigned int dim=1; dim<3; dim++) {
    scale=std::min(1/(mx[dim]-mn[dim]),scale);
  }

  for(size_t ii=0; ii<v.size(); ii++) {
    for (unsigned int dim = 0; dim<3; dim++) {
      v[ii][dim]=v[ii][dim]*scale;
    }
  }
  v0=v;
  compute_norm();

  f.close();
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

void Mesh::drawPlane(int k)
{
  if(k>(int)planes.size()-1){
    k=planes.size()-1;
  }
  Vec3 nn = planes[k].n;
  Vec3 v0=planes[k].c;
  Vec3 ax, ay;
  Vec3 arbit(1,0,0);
  if(std::abs(nn[0])>0.9) {
    arbit=Vec3(0,1,0);
  }
  ax = arbit.cross(nn);
  ax/=ax.norm();
  ay=nn.cross(ax);
  ay/=ay.norm();

  glDisable(GL_LIGHTING);
  glBegin(GL_TRIANGLES);
  if(tex_buf){
    glBindTexture(GL_TEXTURE_2D,texture);
  }

  for(unsigned int ii=0; ii<t.size(); ii++) {
    int l = t[ii].label;
    if(l!=k){
      continue;
    }
    if(tex_buf && tex.size()>0){
      for(int jj=0;jj<3;jj++){
        Vec3 vert=v[t[ii][jj]];
        vert-=v0;
        vert=Vec3(vert.dot(ax),vert.dot(ay),vert.dot(nn));
        glNormal3f(0,0,1);
        glTexCoord2f(tex[t[ii].texId[jj]][0],tex[t[ii].texId[jj]][1]);
        glVertex3f(vert[0],vert[1],-1);
      }
    }else{
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

void Mesh::draw(std::vector<Vec3>&v)
{
 // glDisable(GL_LIGHTING);
  glBegin(GL_TRIANGLES);
  GLfloat specular[4]= {0.51f,0.51f,0.51f,1.0f};
  GLfloat ambient[4]= {0.1f,0.1f,0.1f,1.0f};

  glMaterialfv(GL_FRONT,GL_SPECULAR,specular);
  glMaterialfv(GL_FRONT,GL_AMBIENT,ambient);
  GLfloat s=10;
  glMaterialfv(GL_FRONT_AND_BACK,GL_SHININESS,&s);
  if(tex_buf){
    glBindTexture(GL_TEXTURE_2D,texture);
  }

  for(unsigned int ii=0; ii<t.size(); ii++) {
    int l = t[ii].label;
/*    if(ptx && (int)ii< ptx->numFaces()){
      float ptxColor[4];
      ptx->getPixel(ii,0,0,ptxColor,0,4);
      glColor3f(ptxColor[0],ptxColor[1],ptxColor[2]);
      GLfloat diffuse[4]= {ptxColor[0],ptxColor[1],ptxColor[2],1.0f};
      glMaterialfv(GL_FRONT,GL_DIFFUSE,diffuse);
    }
    else{
  */
    GLfloat diffuse[4]= {color[l][0],color[l][1],color[l][2],1.0f};
    glColor3f(color[l][0],color[l][1],color[l][2]);
   // glMaterialfv(GL_FRONT,GL_DIFFUSE,diffuse);

    Vec3 a = v[t[ii][1]] - v[t[ii][0]];
    Vec3 b = v[t[ii][2]] - v[t[ii][0]];
    b=a.cross(b);
    b= b/b.norm();
    if(tex_buf && tex.size()>0){
      glNormal3f(n[t[ii][0]][0],n[t[ii][0]][1],n[t[ii][0]][2]);
      glTexCoord2f(tex[t[ii].texId[0]][0],tex[t[ii].texId[0]][1]);
      glVertex3f(v[t[ii][0]][0],v[t[ii][0]][1],v[t[ii][0]][2]);

      glTexCoord2f(tex[t[ii].texId[1]][0],tex[t[ii].texId[1]][1]);
      glNormal3f(n[t[ii][1]][0],n[t[ii][1]][1],n[t[ii][1]][2]);
      glVertex3f(v[t[ii][1]][0],v[t[ii][1]][1],v[t[ii][1]][2]);

      glTexCoord2f(tex[t[ii].texId[2]][0],tex[t[ii].texId[2]][1]);
      glNormal3f(n[t[ii][2]][0],n[t[ii][2]][1],n[t[ii][2]][2]);
      glVertex3f(v[t[ii][2]][0],v[t[ii][2]][1],v[t[ii][2]][2]);
    }else{
      glNormal3f(n[t[ii][0]][0],n[t[ii][0]][1],n[t[ii][0]][2]);
      glVertex3f(v[t[ii][0]][0],v[t[ii][0]][1],v[t[ii][0]][2]);
      glNormal3f(n[t[ii][1]][0],n[t[ii][1]][1],n[t[ii][1]][2]);
      glVertex3f(v[t[ii][1]][0],v[t[ii][1]][1],v[t[ii][1]][2]);
      glNormal3f(n[t[ii][2]][0],n[t[ii][2]][1],n[t[ii][2]][2]);
      glVertex3f(v[t[ii][2]][0],v[t[ii][2]][1],v[t[ii][2]][2]);

    }
    if(t[ii].label==highlight){
      glVertex3f(v[t[ii][0]][0]+b.x[0],v[t[ii][0]][1]+b.x[1],v[t[ii][0]][2]+b.x[2]);

      glVertex3f(v[t[ii][1]][0]+b.x[0],v[t[ii][1]][1]+b.x[1],v[t[ii][1]][2]+b.x[2]);

      glVertex3f(v[t[ii][2]][0]+b.x[0],v[t[ii][2]][1]+b.x[1],v[t[ii][2]][2]+b.x[2]);
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
void Mesh::drawLines()
{
  glDisable(GL_LIGHTING);
  glBegin(GL_LINES);
  glColor3f(0.9,0.9,0.4);
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
  pthread_mutex_unlock(&meshm);
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
/**@param m assume triangle norms and centers are already computed
 */
void get_plane(Mesh & m , std::vector<Plane> & plane)
{
  plane.resize(m.nLabel);
  std::vector<float > cnt (m.nLabel,0);
  for(size_t ii=0; ii<m.t.size(); ii++) {
    Vec3 a = m.v[m.t[ii][1]] -  m.v[m.t[ii][0]];
    Vec3 b = m.v[m.t[ii][2]] -  m.v[m.t[ii][0]];
    Vec3 n = a.cross(b);
    float area = n.norm();
    int label = m.t[ii].label;
    cnt[label]+=area;
    plane[label].n += n;
    a=(m.t[ii].c * area);
    plane[label].c += a;
  }

  for(int ii=0; ii<m.nLabel; ii++) {
    plane[ii].c/=cnt[ii];
   // plane[ii].n/=cnt[ii];
    plane[ii].n/=plane[ii].n.norm();
  }
}
