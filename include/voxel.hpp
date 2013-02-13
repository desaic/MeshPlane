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
  ///@param _res must be power of 2
  Voxel(Vec3f mn, Vec3f mx,int _res = 128);
  ~Voxel();
  Grid * grid;
  Vec3f gridlen;
  Vec3f orig;
  int res;
  void draw();
  ///@brief insert a triangle
  ///@param tidx triangle index
  void insert(int tidx, Mesh & m);
  void remove(int tidx, Mesh & m);
  int intersect(int tidx, Mesh & m);
private:
  static Mesh cube;
  int clampIdx(int idx);
  void bbox(int ii, Mesh & m , int * tmn, int * tmx);
  void rasterize(int tidx, Mesh & m , Grid & grid);
  ///@brief computes grid index given coordinate
  void vec2grid(Vec3f & v, GridIdx & grid);
  bool trigCubeIntersect(int tidx, const Mesh & m ,
                         GridIdx & cube);
};
#endif
