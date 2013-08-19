#include "voxel.hpp"
#include "trigAABBIntersect.hpp"
#include "triTriIntersect.hpp"
#include <set>
#include <GL/gl.h>
#include <GL/glut.h>
Mesh Voxel::cube;

bool triTriIntersect(int ia, int ib, Mesh & m)
{
  int res = NoDivTriTriIsect(m.v[m.t[ia][0]].v,m.v[m.t[ia][1]].v,m.v[m.t[ia][2]].v,
                             m.v[m.t[ib][0]].v,m.v[m.t[ib][1]].v,m.v[m.t[ib][2]].v);
  return res;
}

int Voxel::intersect(int tidx, Mesh & m)
{
  //bounding box of a triangle
  int tmn[3], tmx[3];
  bbox(tidx,m,tmn,tmx);
  IntSet::const_iterator it;
  for(int ix=tmn[0];ix<=(tmx[0]);ix++){
    for(int iy=tmn[1];iy<=(tmx[1]);iy++){
      for(int iz=tmn[2];iz<=(tmx[2]);iz++){
        GridIdx gi(ix,iy,iz);
        if(!trigCubeIntersect(tidx,m,gi)){
          continue;
        }
        const IntSet & trigSet = (*grid)(ix,iy,iz);
        for(it = trigSet.cbegin(); it!=trigSet.cend(); it++ ){
          int ib = (*it);
          if(m.isNbr(tidx,ib)){
            continue;
          }
          if(triTriIntersect(tidx, ib,m)){
            return ib;
          }
        }
      }
    }
  }
  return -1;
}

void Voxel::remove(int tidx,Mesh & m )
{
  //bounding box of a triangle
  int tmn[3], tmx[3];
  bbox(tidx,m,tmn,tmx);
  for(int ix=tmn[0];ix<=(tmx[0]);ix++){
    for(int iy=tmn[1];iy<=(tmx[1]);iy++){
      for(int iz=tmn[2];iz<=(tmx[2]);iz++){
        GridIdx gi(ix,iy,iz);
        if(trigCubeIntersect(tidx,m,gi)){
          (*grid)(ix,iy,iz).erase(tidx);
        }
      }
    }
  }
}

void Voxel::insert(int tidx,Mesh & m)
{
  rasterize(tidx,m,*grid);
}

Voxel::~Voxel()
{
  if(grid!= 0){
    delete grid;
  }
}

Voxel::Voxel(Vec3f mn, Vec3f mx, int _res)
:grid(0),res(_res)
{
  if(cube.v.size()==0){
    std::ifstream infile("cube.ply2");
    Voxel::cube.read_ply2(infile);
    infile.close();
  }
  grid = new Grid(res);
  gridlen = (1.0f/res)*(mx-mn);
  orig = mn;
}

void Voxel::draw()
{
  glEnable(GL_LIGHTING);
  for(int ii =0 ;ii<res;ii++){
    for(int jj =0 ;jj<res;jj++){
      for(int kk =0 ;kk<res;kk++){
        if( (*grid)(ii,jj,kk).size()==0){
          continue;
        }
        glPushMatrix();
        glTranslatef(ii*gridlen[0],jj*gridlen[1],kk*gridlen[2]);
        glutSolidCube(gridlen[0]);
        glPopMatrix();
      }
    }
  }
}

bool Voxel::trigCubeIntersect(int tidx, const Mesh & m ,
                       GridIdx & cube)
{
  float boxcenter[3]={(float)((0.5+cube[0])*gridlen[0]+orig[0]),
                      (float)((0.5+cube[1])*gridlen[1]+orig[1]),
                      (float)((0.5+cube[2])*gridlen[2]+orig[2])};
  float boxhalfsize[3]={(float)gridlen[0]/(1.99f),
                      (float)gridlen[1]/(1.99f),
                      (float)gridlen[2]/(1.99f)};
  float triverts[3][3];
  for(int ii=0;ii<3;ii++){
    for(int jj=0;jj<3;jj++){
      triverts[ii][jj]=m.v[m.t[tidx][ii]][jj];
    }
  }
  return triBoxOverlap(boxcenter,boxhalfsize, triverts);
}

int Voxel::clampIdx(int idx)
{
  idx = std::max(0,idx);
  idx = std::min(idx, res-1);
  return idx;
}

void Voxel::vec2grid(Vec3f & v,  GridIdx & grid)
{
  for(int ii=0; ii<VOXEL_DIM; ii++) {
    grid[ii]= (int)((v[ii]-orig[ii])/gridlen[ii]);
    grid[ii] = clampIdx(grid[ii]);
  }
}

void Voxel::rasterize(int tidx, Mesh & m , Grid & grid)
{
  //bounding box of a triangle
  int tmn[3], tmx[3];
  bbox(tidx,m,tmn,tmx);
  for(int ix=tmn[0];ix<=(tmx[0]);ix++){
    for(int iy=tmn[1];iy<=(tmx[1]);iy++){
      for(int iz=tmn[2];iz<=(tmx[2]);iz++){
        GridIdx gi(ix,iy,iz);
        if(trigCubeIntersect(tidx,m,gi)){
          //if(grid(ix,iy,iz).size()==0){
          //grid.set(ix,iy,iz,IntSet());
          //}
          grid(ix,iy,iz).insert(tidx);
        }
      }
    }
  }
}

void Voxel::bbox(int ii, Mesh & m , int * tmn, int * tmx)
{
  GridIdx vidx;
  vec2grid(m.v[m.t[ii][0]],vidx);
  for(int jj=0; jj<3; jj++) {
    tmn[jj]=vidx[jj]-1;
    tmx[jj]=vidx[jj];
  }
  for(int jj=1; jj<3; jj++) {
    vec2grid(m.v[m.t[ii][jj]],vidx);
    for(int kk=0; kk<VOXEL_DIM; kk++) {
      if(vidx[kk]-1<tmn[kk]) {
        tmn[kk]=vidx[kk]-1;
      }
      if(vidx[kk]>tmx[kk]) {
        tmx[kk]=vidx[kk];
      }
    }
  }
  for(int ii = 0;ii<VOXEL_DIM;ii++){
    tmn[ii] = std::max(0,tmn[ii]);
    tmx[ii] = std::min(res-1,tmx[ii]); 
  }
}
