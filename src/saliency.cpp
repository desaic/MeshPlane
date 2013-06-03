#include "mesh.hpp"
#include "saliency.hpp"
#include "imageio.h"
#include <cstdlib>  //for abs
#define IX(ii,jj)  (4* ((ii)*width+(jj))  )
#define IXR(ii,jj) (4* ((ii)*width+(jj))  )
#define IXG(ii,jj) (4* ((ii)*width+(jj)) +1)
#define IXB(ii,jj) (4* ((ii)*width+(jj)) +2)
#define IXA(ii,jj) (4* ((ii)*width+(jj)) +3)

#include "util.h"
class Mapper {
public:
  unsigned char * buf;
  int width, height;
  Mapper(unsigned char* _buf=0, int _width=0, int _height=0):
    buf(_buf),width(_width),height(_height) {}
  virtual void map(int x, int y) {};
};

class MaxMapper:public Mapper {
public:
  unsigned int max;
  MaxMapper(unsigned char* _buf=0, int _width=0, int _height=0):Mapper(_buf,_width,_height),
    max(0){}
  void map(int x, int y ) {
    max=std::max((unsigned int)buf[IX(y,x)],max);
  }
};

class SumMapper:public Mapper {
public:
  unsigned int sum;
  unsigned int cnt;
  SumMapper(unsigned char* _buf=0, int _width=0, int _height=0):Mapper(_buf,_width,_height),
    sum(0) ,cnt(0){}
  void map(int x, int y ) {
    sum+=(unsigned int)buf[IX(y,x)];
    cnt++;
  }
};

void line(int x0,int  y0, int x1, int y1, Mapper & mapper );

void trace_edge(Mesh & m, int tidx, int kk0,int kk1, Mapper & mapper)
{
  Vec3 v0,v1;
  v0=m.tex[ m.t[tidx].texId[kk0] ];
  v1=m.tex[ m.t[tidx].texId[kk1] ];
  int x0,y0,x1,y1;
  x0=(int)(v0[0]*mapper.width);
  y0=(int)(v0[1]*mapper.height);
  x1=(int)(v1[0]*mapper.width);
  y1=(int)(v1[1]*mapper.height);
  line(x0,y0,x1,y1,mapper);
}

void usr_map(const char * filename, Mesh & m)
{
  int width, height;
  unsigned char * buf = imageio_load_image(filename, &width, &height);
  m.adjlist();

  for(size_t tidx=0; tidx<m.t.size(); tidx++) {
    for(size_t jj=0; jj<m.adjMat[tidx].size(); jj++) {
      int nbr=m.adjMat[tidx][jj];
      int kk=0;
      bool foundExclude=false;
      for(kk=0; kk<3; kk++) {
        int vidx=m.t[tidx][kk];
        if(!contains(m.t[nbr].x, vidx,3)) {
          foundExclude=true;
          break;
        }
      }
      if(!foundExclude) {
        //shouldn't really happen
        continue;
      }
      SumMapper mapper(buf,width,height);
      trace_edge(m,tidx,(kk+1)%3,(kk+2)%3,mapper);
      trace_edge(m,nbr,(kk+1)%3,(kk+2)%3,mapper);
      m.usr_weit[EdgeId(tidx,nbr)]=(mapper.sum/255.0)/mapper.cnt;
    }
  }
  delete []buf;
}

void saliency_map(const char * filename, Mesh & m)
{
  int width, height;
  unsigned char * buf = imageio_load_image(filename, &width, &height);
  m.adjlist();
  MaxMapper mapper(buf,width,height);
  for(size_t tidx=0; tidx<m.t.size(); tidx++) {
    mapper.max=0;
    for(size_t jj=0; jj<m.adjMat[tidx].size(); jj++) {
      int nbr=m.adjMat[tidx][jj];
      int kk=0;
      bool foundExclude=false;
      for(kk=0; kk<3; kk++) {
        int vidx=m.t[tidx][kk];
        if(!contains(m.t[nbr].x, vidx,3)) {
          foundExclude=true;
          break;
        }
      }
      if(!foundExclude) {
        //shouldn't really happen
        continue;
      }
      trace_edge(m,tidx,(kk+1)%3,(kk+2)%3,mapper);
      trace_edge(m,nbr,(kk+1)%3,(kk+2)%3,mapper);
      m.saliency[EdgeId(tidx,nbr)]=mapper.max/255.0;
    }
  }
  delete []buf;
}


void line(int x0,int  y0, int x1, int y1, Mapper & mapper ) {
  int dx = std::abs(x1-x0);
  int dy = std::abs(y1-y0);
  int sx,sy,err;
  if (x0 < x1 ) {
    sx = 1;
  } else {
    sx = -1;
  }
  if (y0 < y1) {
    sy = 1;
  } else {
    sy = -1;
  }
  err = dx-dy;

  while(1) {
    mapper.map(x0,y0);
    if (x0 == x1 && y0 == y1) {
      break;
    }
    int e2 = 2*err;
    if (e2 > -dy) {
      err = err - dy;
      x0 = x0 + sx;
    }
    if (e2 <  dx ) {
      err = err + dx;
      y0 = y0 + sy;
    }
  }
}
