#include "voxel.hpp"
#include "trigAABBIntersect.hpp"
#include <set>
#include <GL/gl.h>

void Voxel::draw()
{
}

bool Voxel::trigCubeIntersect(int tidx, const Mesh & m ,
                       GridIdx & cube)
{
  float boxcenter[3]={(float)((0.5+cube[0])*gridlen),
                      (float)((0.5+cube[1])*gridlen),
                      (float)((0.5+cube[2])*gridlen)};
  float boxhalfsize[3]={(float)gridlen/(1.99f),
                        (float)gridlen/(1.99f),
                        (float)gridlen/(1.99f)};
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

void Voxel::vec2grid(Vec3 & v,  GridIdx & grid)
{
  for(int ii=0; ii<VOXEL_DIM; ii++) {
    grid[ii]= (int)(v[ii]/gridlen);
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
          if(grid(ix,iy,iz).size()==0){
            grid.set(ix,iy,iz,IntSet());
          }
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
}

///@param m should be scaled to [0.2 0.8]
Mesh Voxel::cube;
Voxel::Voxel(Mesh & m):
  res(128),
  grid(res,IntSet(0))
{
  gridlen=1.0f/res;

  for(size_t ii=0; ii<m.t.size(); ii++) {
    rasterize(ii,m,grid);
  }

  if(cube.v.size()==0){
    std::ifstream infile("cube.ply2");
    Voxel::cube.read_ply2(infile);
    infile.close();
  }
  for(size_t ii=0;ii<cube.v.size();ii++){
    cube.v[ii]/=(float)res;
  }
}

