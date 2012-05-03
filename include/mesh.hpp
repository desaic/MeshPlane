#ifndef MESH_H
#define MESH_H
#include <vector>
#include <fstream>
#include "math.hpp"
#include "Ptexture.h"
#include <GL/gl.h>
//by default counterclockwise winding is front face
struct Trig{
  Trig(){x[0]=0;x[1]=0;x[2]=0;x[3]=0;}
  int & operator[](const int i) {return x[i];}
  int operator[](const int i)const {return x[i];}
  int x[4];
  int label;
  Vec3 n,c,a;
};


struct Plane{
  Vec3 c;
  Vec3 n;
  Plane(){};
  ~Plane(){};
};
#include <map>
class Mesh{
public:
  int nLabel;
  std::map<int,bool>bad;
  bool self_intersect();
  std::vector<Vec3>v;
  std::vector<Vec3>tex;
  std::vector<Vec3>v0;
  std::vector<Trig>t;
  std::vector<Vec3>n;
  std::vector<Vec3>v4;
  std::vector<std::vector<int> >  adjMat;
  Mesh():v(0),t(0){}
  Mesh(const char * filename,int _nLabel=50);
  void adjlist();
  void draw(std::vector<Vec3>&v);
  void drawLines();
  void assign_color();
  void save(const char * filename);
  void save_plane(const char * filename);
  bool has_v4;
  real_t area(const Trig & trig);
  void get_normal_center();
  std::vector<Plane>planes;
  int highlight;
  void compute_plane();
  void read_ply(std::ifstream & f);
  void read_ply2(std::ifstream & f);
  void load_tex(const char * filename);
  void load_ptex(const char * filename);
private:
  void compute_norm();
  void fix_inner_cluster();
  PtexTexture *ptx;
  std::vector<Vec3>color;
  std::vector< std::vector<std::vector<int> > > lines;
  GLuint texture;
  unsigned char * tex_buf;
};
void randcenter(Mesh & m,std::vector<Plane>&plane, int nLabel);
/**@param m assume triangle norms and centers are already computed
 */
void get_plane(Mesh & m , std::vector<Plane> & plane);

#endif
