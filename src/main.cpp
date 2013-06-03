#ifdef _WIN32
#include <windows.h>
#endif

#include "cgd.hpp"
#include "kmeans.hpp"
#include "mesh.hpp"
#include "mesh_query.h"
#include "mincut.hpp"
#include <imageio.h>
#include "quat.h"
#include "saliency.hpp"
#include "voxel.hpp"

#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/freeglut.h>
#include <cmath>
#include <pthread.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <sstream>
static Mesh * m;

static Quat rot;
static Voxel * voxel;
static int planeId=0;
static bool draw_tex=false;
static bool draw_uv=false;
static int imgnum=0;
bool running=false;
int ldown;
int oldx,oldy;
extern int minc_nlabel;

struct Cam{
  Cam():rotx(0),roty(0){
    for (int ii=0;ii<3;ii++){
      eye[ii]=0.0;
      at[ii]=0.0;
    }
    eye[2]=2;
  }
  GLdouble eye[3];
  GLdouble at[3];
  float rotx,roty;
};
static Cam* cam;
void init(void)
{
  glClearColor (1, 1, 1, 0.0);
  // glShadeModel (GL_SMOOTH);
  //glShadeModel (GL_FLAT );
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_LIGHT1);
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_ALPHA);
  // glEnable(GL_NORMALIZE);

  GLfloat white[]={1.0,1.0,1.0,1.0};
  glLightfv (GL_LIGHT1, GL_DIFFUSE, white);
  glLightfv (GL_LIGHT1, GL_SPECULAR, white);

  if(draw_tex||draw_uv){
    glMatrixMode (GL_PROJECTION);
    glLoadIdentity();
    glOrtho(-0.5,0.5,-0.5,0.5,10,0.1);
  }
}

void display(void)
{
  glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  if(draw_tex){
 //   m->drawPlane(planeId);
    glFlush ();
    if(planeId<(int)m->planes.size()){
      std::stringstream ss;
      ss<<planeId<<".png";
      imageio_save_screenshot(ss.str().c_str());
      planeId++;
    }
    return;
  }
  else if(draw_uv){
    m->draw_tex();
    glFlush();
    if(planeId==0){
      imageio_save_screenshot("remap_tex.png");
      planeId++;
    }
    return;
  }
  if(running){
    wPt+=1;
    if(imgnum>230){
      wPt+=200;
    }
    if(imgnum%10==0){
    //  runMincut(*m);
      //runKmeans(*m);
    }
    m->compute_plane();
    cgd(*m);
  }
  gluLookAt(cam->eye[0], cam->eye[1],cam->eye[2],
	    cam->at[0],cam->at[1], cam->at[2],
	    0.0, 1.0, 0.0);

  GLfloat position[] = { 1.0, -1, 1, 1.0 };
  GLfloat position1[] = { -1.0, -1, -1, 1.0 };

  glPushMatrix();
  glDisable(GL_LIGHTING);
  glColor3f(1,1,1);
  glTranslatef(position[0],position[1],position[2]);
 // glutWireCube(0.2);
  glPopMatrix();

  glPushMatrix();
  glTranslatef(position1[0],position1[1],position1[2]);
 // glutWireCube(0.2);
  glEnable(GL_LIGHTING);
  glPopMatrix();

  glLightfv (GL_LIGHT0, GL_POSITION, position);
  glLightfv (GL_LIGHT1, GL_POSITION, position1);

  Vec3 axis ;
  real_t angle;
  glPushMatrix();
  //glTranslatef(-1,0,0);
  rot.to_angle_axis(axis,&angle);
  angle=angle*180/3.14159;

  glRotatef(angle,axis[0],axis[1],axis[2]);
  if(m!=0){
    m->draw(m->v);
//    m->drawPlane(1);
//    m->drawPlane(2);
 //   m->drawPlane(3);
 //   m->drawPlane(4);
 //   m->drawPlane(5);
 //   m->drawPlane(6);
  }
  if(voxel!=0){
    voxel->draw();
  }
  glBindFramebuffer(GL_FRAMEBUFFER, m->fbo);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  m->drawCol();
  glBindFramebuffer(GL_FRAMEBUFFER, 0);
  glPopMatrix();

  //glPushMatrix();
  //glTranslatef(1,0,0);
  //glRotatef(angle,axis[0],axis[1],axis[2]);

//    m->draw(m->v0);
  //m->drawLines();
  glPopMatrix();
  glFlush ();

  if(running){

    char buf[32]="anim/";
    sprintf(buf+5,"%04d",imgnum);
    strcat(buf,".png");
    imageio_save_screenshot(buf);
    printf("%d\n",imgnum);
    imgnum++;
  }

}

void reshape (int w, int h)
{
  glViewport (0, 0, (GLsizei) w, (GLsizei) h);
  glMatrixMode (GL_PROJECTION);
  glLoadIdentity();
  if(draw_tex||draw_uv){
    glOrtho(-0.5,0.5,-0.5,0.5,10,0.1);
  }else{
    gluPerspective(40.0, (GLfloat) w/(GLfloat) h, 0.1, 20.0);
  }
}

GLdouble  norm(GLdouble * v)
{
  GLdouble sum=0;
  for(int ii=0;ii<3;ii++){
    GLdouble d=cam->eye[ii]-cam->at[ii];
    sum+=d*d;
  }
  sum=sqrt(sum);
  return sum;
}

void add(GLdouble * a, GLdouble * b,GLdouble * ans,GLdouble c)
{
  for(int ii=0;ii<3;ii++){
    ans[ii]=a[ii]+c*b[ii];
  }
}
void mul(GLdouble * a,GLdouble c)
{
  for(int ii=0;ii<3;ii++){
    a[ii]*=c;
  }
}
void keyboard(unsigned char key,int x, int y)
{

  Vec3 axis ;
  real_t angle;
  switch(key){
  case 'w':
    cam->eye[2]-=0.1;
    break;
  case 's':
    cam->eye[2]+=0.1;
    break;
  case 'a':
    cam->eye[0]-=0.1;
    break;
  case 'd':
    cam->eye[0]+=0.1;
    break;
  case 'z':
    rot.to_angle_axis(axis,&angle);
    angle=angle*180/3.14159;
    std::cout<<axis[0]<<" "<<axis[1]<<" "<<axis[2]<<"\n";
    std::cout<<angle<<"\n";
    break;
  case 'x':
    imageio_save_screenshot("screenshot.png");
    break;
  case 'o':
    m->save_obj("mesh.obj");
    break;
  case 'f':
    m->save_off("off.off");
    break;
  }
  glutPostRedisplay();
}


void mouse(int button, int state, int x, int y)
{
  switch (button) {
  case GLUT_LEFT_BUTTON:
    switch (state){
    case GLUT_DOWN:
      ldown=1;
      oldx=x;
      oldy=y;
      break;
    case GLUT_UP:
      ldown=0;
      glBindFramebuffer(GL_FRAMEBUFFER,m->fbo);
      glFlush();
      unsigned char buf[4];
      glReadPixels(x, 720-y, 1, 1, GL_RGBA, GL_UNSIGNED_BYTE, buf);
      unsigned int l = b2int((GLubyte*)buf);
      std::cout<<l<<"\n";
      glBindFramebuffer(GL_FRAMEBUFFER,0);
      break;
    }
    break;
  default:
    break;
  }
}
void motion (int x, int y)
{

  if(ldown){
    Quat q(Vec3(0,1,0),(x-oldx)/60.0);
    rot =q*rot;
    rot=rot.normalize();
    q=Quat(Vec3(1,0,0),(y-oldy)/60.0);
    rot =q*rot;
    rot=rot.normalize();
    oldx=x;
    oldy=y;
  }

}

void animate(int t)
{
  glutTimerFunc(60, animate, 0);

  glutPostRedisplay();

}

void* iterate(void* arg){
  int ITER=100;
  MC_ITER=1;
  Mesh * m=(Mesh*)arg;
  wS=0.1;
  wI=1;
  wV0=50;
  wPt=0.5;
  vW=1;
  dataCostW=2000;
//  smoothW=400;
  saliency_weight=5;
  distw=15;

 // BP bp(*m);
  initKmeans(*m);
  for(int ii=0;ii<ITER;ii++){
    wPt+=3;
    printf("iter %d\n",ii);
    runMincut(*m);
    //runKmeans(*m);
    printf("cut\n");
    m->compute_plane();
    cgd(*m);
//    weighted_avg(*m);
  }
  wPt=1000;
  cgd(*m);
  wPt=2000;
  cgd(*m);
  
  m->save("planar_output.ply2");
  m->compute_plane();
  m->save_plane("plane.txt");
  return 0;
}
void *scanHighlight(void * arg)
{
  Mesh * m = (Mesh*)arg;
  while(1){
    std::string s;
    std::getline(std::cin,s);
    std::stringstream ss(s);
    ss>>m->highlight;
  }
  return 0;
}

void testSelfInter()
{
  m=new Mesh ("intersect.obj",1, true);
  Vec3f mn,mx;
  
  BBox(*m,mn,mx);
  voxel=new Voxel(mn,mx);
  std::cout<<"Use Grid\n";
  for(unsigned int ii = 0;ii<m->t.size();ii++){
    int inter = voxel->intersect(ii,*m);
    if(inter>=0){
      std::cout<<inter<<" and "<<ii<<"\n";
    }
    voxel->insert(ii,*m);
  }
  
  bool inter = m->self_intersect();
  std::cout<<"expect intersection\n"<<inter<<"\n";
  delete m;
  m=new Mesh ("nointer.obj",1, true);
  inter = m->self_intersect();
  std::cout<<"expect no intersection\n"<<inter<<"\n";
  delete m;  
}

int main(int argc, char** argv)
{
  if(argc<2){
    printf("%s filename\n",argv[0]);
    exit(0);
  }

  if(strcmp(argv[1],"-test")==0){
    testSelfInter();
 //   return 0;
  }
  
  glutInit(&argc, argv);
  int nLabel=50;
  bool run=false;
  bool copy_uv=false;
  const char * tex_file="";
  const char * label_file="";
  const char * uvfile="";
  const char * salfile="";
  const char * usrfile="";
  bool autoscale = true;
  bool checkIntersect=false;
  for(int ii=0;ii<argc;ii++){
    if(strcmp(argv[ii], "-k")==0){
      ii++;
      nLabel=atoi(argv[ii]);
    }
    if(strcmp(argv[ii],"-r")==0){
      run=true;
    }
    if(strcmp(argv[ii],"-t")==0){
      ii++;
      tex_file=argv[ii];
    }
    if(strcmp(argv[ii],"-l")==0){
      ii++;
      label_file=argv[ii];
    }
    if(strcmp(argv[ii],"-uv")==0){
      draw_uv=true;
    }
    if(strcmp(argv[ii],"-tex")==0){
      draw_tex=true;
    }
    if(strcmp(argv[ii],"-uvfile")==0){
      ii++;
      uvfile=argv[ii];
    }
    if(strcmp(argv[ii],"-copyuv")==0){
      copy_uv=true;
    }
    if(strcmp(argv[ii],"-sal")==0){
      ii++;
      salfile=argv[ii];
    }
    if(strcmp(argv[ii],"-noscale")==0){
      autoscale=false;
    }
    if(strcmp(argv[ii],"-inter")==0){
      checkIntersect=true;
    }
    if(strcmp(argv[ii], "-smoothcost")==0){
      ii++;
      smoothW=atoi(argv[ii]);
    }
    if(strcmp(argv[ii],"-usr")==0){
      ii++;
      usrfile=argv[ii];
    }
  }

  glutInitDisplayMode (GLUT_SINGLE | GLUT_RGB | GLUT_DEPTH|GLUT_ALPHA );
  if(draw_tex||draw_uv){
    glutInitWindowSize (2400, 2400);
  }else{
    glutInitWindowSize (1280 , 720);
  }
  //glutInitWindowPosition (100, 100);
  glutCreateWindow (argv[0]);
  init ();
  glutDisplayFunc(display);
  glutReshapeFunc(reshape);
  glutMouseFunc(mouse);
  glutMotionFunc(motion);
  glutKeyboardFunc(keyboard);
  glutTimerFunc(0.1, animate, 0);

  srand(123456);
  m=new Mesh (argv[1],nLabel, autoscale);
  m->checkIntersect = checkIntersect;
  m->init_select("shader/select.glsl");
  if(tex_file[0]){
    m->load_tex(tex_file);
  }
  minc_nlabel=nLabel;
  if(label_file[0]){
    Mesh mtemp(label_file,1);
    m->v=mtemp.v;
    m->nLabel=mtemp.nLabel;
    for(size_t ii=0;ii<m->t.size();ii++){
      m->t[ii].label=mtemp.t[ii].label;
    }
    m->assign_color();
  }

  if(uvfile[0]){
    Mesh * uvmesh=new Mesh(uvfile,1);
    m->remap_tex=uvmesh;

    if(copy_uv){
      m->tex=uvmesh->tex;
      for(size_t ii=0;ii<m->t.size();ii++){
        for(int jj=0;jj<3;jj++){
          m->t[ii].texId[jj]=uvmesh->t[ii].texId[jj];
        }
      }
    }
    m->save_obj("remap_tex.obj");
  }

  if(salfile[0]){
    saliency_map(salfile,*m);
  }
  if(usrfile[0]){
    usr_map(usrfile,*m);
  }
  m->compute_plane();
  if(draw_tex||draw_uv){
    glMatrixMode (GL_PROJECTION);
    glLoadIdentity();
    glOrtho(-0.5,0.5,-0.5,0.5,10,0.1);
  }

  m->save_plane("plane.txt");

  //draw_tex=true;
  if(run){

    wS=1;
    wI=1;
    wV0=20;
    wPt=1;
    vW=1;
    dataCostW=1000;
    //smoothW=2500;
    saliency_weight=10;
    distw=3;
    MC_ITER=1;
    m->compute_plane();
    //  running=true;
    pthread_t thread;
    pthread_create(&thread, 0, iterate,(void*)m);
    pthread_detach(thread);
  }
  pthread_t hlthread;
  pthread_create(&hlthread,0,scanHighlight,(void*)m);
  pthread_detach(hlthread);
  cam=new Cam();
  rot=Quat(Vec3(1,0,0),0);
  ldown=0;
  glutMainLoop();
  return 0;
}

