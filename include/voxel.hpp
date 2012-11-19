#ifndef VOXEL_HPP
#define VOXEL_HPP
#include "vec.h"
#include "mesh.hpp"
#include <vector>
#include "mesh.hpp"
#include "octree.h"
#include <unordered_set>
#define VOXEL_DIM 3

struct GridIdx{
  int i[3];
  GridIdx(int a=0,int b=0,int c=0){
    i[0]=a;
    i[1]=b;
    i[2]=c;
  }
  int & operator[](int idx){return i[idx];}
  int operator[](int idx)const{return i[idx];}

  bool operator<(const GridIdx & grid)const{
    for(int ii=0;ii<VOXEL_DIM;ii++){
      if(i[ii]==grid.i[ii]){
        continue;
      }
      if(i[ii]<grid.i[ii]){
        return true;
      }
      return false;
    }
    return false;
  }
};

typedef std::unordered_set<int> IntSet;
typedef Octree<IntSet> Grid;

class Voxel{
public:
  Voxel(Mesh & m);
  int res;
  real_t gridlen;
  Grid grid;
  void draw();
  ///@brief insert a triangle
  ///@param tidx triangle index
  void insert(Mesh & m , int tidx);
  void remove(Mesh & m , int tidx);
private:
  static Mesh cube;
  int clampIdx(int idx);
  void bbox(int ii, Mesh & m , int * tmn, int * tmx);
  void rasterize(int tidx, Mesh & m , Grid & grid);
  void vec2grid(Vec3 & v, GridIdx & grid);
  bool trigCubeIntersect(int tidx, const Mesh & m ,
                         GridIdx & cube);
};
#endif
