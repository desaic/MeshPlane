#ifndef MESH_H
#define MESH_H
#include <vector>
#include <fstream>

#ifdef _WIN32
#define NOMINMAX
#include <Windows.h>
#endif 

#include <GL/glew.h>
#include "math.hpp"

//#include <GL/gl.h>
#include <map>
//by default counterclockwise winding is front face
struct Trig{
  Trig():label(0){x[0]=0;x[1]=0;x[2]=0;x[3]=0;
    texId[0]=0;texId[1]=0;texId[2]=0;}
  int & operator[](const int i) {return x[i];}
  int operator[](const int i)const {return x[i];}
  int x[4];
  int label;
  Vec3 n,c;
  real_t A;
  int texId[3];
};


struct Plane{
  Vec3 c;
  Vec3 n;
  real_t A;
  Plane():A(0){};
  ~Plane(){};
};

struct EdgeId{
  int v[2];
  int & operator[](int idx){
    return v[idx];
  }
  int operator[](int idx)const {
    return v[idx];
  }
  EdgeId(int v0=0,int v1=0){
    if(v0>v1){
      v[0]=v1;
      v[1]=v0;
    }else{
      v[0]=v0;
      v[1]=v1;
    }
  }
  bool operator<(const EdgeId & ea)const{
    if(v[0]<ea.v[0]){
      return true;
    }
    if(v[0]==ea.v[0]){
      return v[1]<ea.v[1];
    }
    return false;
  }
};


class Mesh{
public:
  int nLabel;
  int targetLabel;
  std::map<int,bool>bad;
  bool self_intersect();
  std::vector<Vec3>v;
  std::vector<Vec3>tex;
  std::vector<Vec3>v0;
  std::vector<Trig>t;
  std::vector<Vec3>n;
  std::vector<Vec3>v4;
  std::vector<std::vector<int> >  adjMat;
  Mesh():v(0),t(0),remap_tex(0){}
  Mesh(const char * filename,int _nLabel=50, bool _auto=true);
  void adjlist();
  void draw(std::vector<Vec3>&v);
  void drawCol();
  void drawLines();
  void drawPlane(int k);
  void assign_color();
  void save(const char * filename);
  void save_plane(const char * filename);
  void save_off(const char * filename);
  bool has_v4;
  void get_normal_center();
  std::vector<Plane>planes;
  int highlight;
  void compute_plane();
  void read_ply(std::ifstream & f);
  void read_ply2(std::ifstream & f);
  void read_obj(std::ifstream & f);
  void save_obj(const char * filename);
  void load_tex(const char * filename);
  void load_ptex(const char * filename);
  void draw_tex();
  void init_select(const char * shaderfile);
  Mesh * remap_tex;
  std::map<EdgeId, real_t > saliency;
  std::map<EdgeId, real_t > usr_weit;
  GLuint fbo;
  bool checkIntersect;
private:
  void compute_norm();
  void fix_inner_cluster();
  std::vector<Vec3>color;
  std::vector< std::vector<std::vector<int> > > lines;
  GLuint texture;
  unsigned char * tex_buf;

  GLuint fbot;
  GLhandleARB select_shader;

  bool autoscale;
};
void randcenter(Mesh & m,std::vector<Plane>&plane, int nLabel);
/**@param m assume triangle norms and centers are already computed
 */
void get_plane(Mesh & m , std::vector<Plane> & plane);
real_t mcdistance( Plane & p, Trig &t);
unsigned int b2int(GLubyte * b);
#endif
