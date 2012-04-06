#ifdef _WIN32
#include <windows.h>
#endif
#include <GL/gl.h>
#include <GL/glu.h>
#include <cmath>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <GL/freeglut.h>
#include "mesh.hpp"
#include "poly.hpp"
#include "mincut.hpp"
#include <string.h>
#include "mesh_query.h"
#include "bp.hpp"
static Mesh * m;
//static Poly * p;
struct Cam{
  Cam():rotx(0),roty(0){
    for (int ii=0;ii<3;ii++){
      eye[ii]=0.5;
      at[ii]=0.0;
    }
    eye[2]=3;
  }
  GLdouble eye[3];
  GLdouble at[3];
  float rotx,roty;
};
static Cam* cam;
void init(void)
{
  glClearColor (0.2, 0.4, 0.9, 0.0);
  // glShadeModel (GL_SMOOTH);
  glShadeModel (GL_FLAT );
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_LIGHT1);
  glEnable(GL_DEPTH_TEST);
  // glEnable(GL_NORMALIZE);


  GLfloat white[]={1.0,1.0,1.0,1.0};
  glLightfv (GL_LIGHT1, GL_DIFFUSE, white);
  glLightfv (GL_LIGHT1, GL_SPECULAR, white);
}

/*  Here is where the light position is reset after the modeling
 *  transformation (glRotated) is called.  This places the
 *  light at a new position in world coordinates.  The cube
 *  represents the position of the light.
 */
void display(void)
{
  glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();



  gluLookAt(cam->eye[0], cam->eye[1],cam->eye[2],
	    cam->at[0],cam->at[1], cam->at[2],
	    0.0, 1.0, 0.0);

	    GLfloat position[] = { 0.0, -.5, 2, 1.0 };
  GLfloat position1[] = { 0.0, .5, -2, 1.0 };

  glPushMatrix();
  glDisable(GL_LIGHTING);
  glColor3f(1,1,1);
  glTranslatef(position[0],position[1],position[2]);
  glutWireCube(0.2);
  glPopMatrix();

  glPushMatrix();
  glTranslatef(position1[0],position1[1],position1[2]);
  glutWireCube(0.2);
  glEnable(GL_LIGHTING);
  glPopMatrix();

  glLightfv (GL_LIGHT0, GL_POSITION, position);
  glLightfv (GL_LIGHT1, GL_POSITION, position1);

  glPushMatrix();
  glTranslatef(-1,0,0);
  glRotatef(180*cam->rotx/3.14,1,0,0);
  glRotatef(180*cam->roty/3.14,0,1,0);
  m->draw(m->v);
  glPopMatrix();

  glPushMatrix();
  glTranslatef(1,0,0);
  glRotatef(180*cam->rotx/3.14,1,0,0);
  glRotatef(180*cam->roty/3.14,0,1,0);

//    m->draw(m->v0);
  m->drawLines();
  glPopMatrix();
  glFlush ();

}

void reshape (int w, int h)
{
  glViewport (0, 0, (GLsizei) w, (GLsizei) h);
  glMatrixMode (GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(40.0, (GLfloat) w/(GLfloat) h, 1.0, 20.0);
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

  case 'x':
    cam->rotx+=0.1;
    if(cam->rotx>6.28){
      cam->rotx-=6.28;
    }
    break;
  case 'z':
    cam->rotx-=0.1;
    if(cam->rotx<-6.28){
      cam->rotx+=6.28;
    }
    break;
  case 'c':
    cam->roty+=0.1;
    if(cam->roty>6.28){
      cam->roty-=6.28;
    }
    break;
  case 'v':
    cam->roty-=0.1;
    if(cam->roty<-6.28){
      cam->roty+=6.28;
    }
    break;

  }
  glutPostRedisplay();
}

void mouse(int button, int state, int x, int y)
{
  switch (button) {
  case GLUT_LEFT_BUTTON:

    break;
  default:
    break;
  }
}

void animate(int t)
{
  glutTimerFunc(60, animate, 0);
  glutPostRedisplay();
}

#include <pthread.h>
#include "cgd.hpp"
extern int minc_nlabel;
void* iterate(void* arg){
  int ITER=50;
  MC_ITER=1;
  Mesh * m=(Mesh*)arg;
  wS=30;
  wI=1;
  wV0=200;
  wPt=0.5;

  vW=30;
  dataCostW=300;
  smoothW=100;
	distw=1;
	BP bp(*m);
  m->compute_plane();

  for(int ii=0;ii<ITER;ii++){
    wPt+=1;
    printf("iter %d\n",ii);
    runMincut(*m);
    m->compute_plane();
    cgd(*m);
//    weighted_avg(*m);
  }
  wPt=2000;
  cgd(*m);
  m->save("planar_output.ply2");

  m->save_plane("plane.txt");
  return 0;
}
#include <iostream>
#include <string>
#include <sstream>
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
int main(int argc, char** argv)
{
  if(argc<2){
    printf("%s filename\n",argv[0]);
    exit(0);
  }
  glutInit(&argc, argv);
  glutInitDisplayMode (GLUT_SINGLE | GLUT_RGB | GLUT_DEPTH);
  glutInitWindowSize (500, 500);
  glutInitWindowPosition (100, 100);
  glutCreateWindow (argv[0]);
  init ();
  glutDisplayFunc(display);
  glutReshapeFunc(reshape);
  glutMouseFunc(mouse);
  glutKeyboardFunc(keyboard);
  glutTimerFunc(0.1, animate, 0);
  int nLabel=50;
  bool run=false;
  for(int ii=0;ii<argc;ii++){
    if(strcmp(argv[ii], "-k")==0){
      ii++;
      nLabel=atoi(argv[ii]);
    }
    if(strcmp(argv[ii],"-r")==0){
      ii++;
      run=true;
    }
  }
  m=new Mesh (argv[1],nLabel);
  minc_nlabel=nLabel;
  m->compute_plane();
  //m->save_plane("plane.txt");
  if(run){
    pthread_t thread;
    pthread_create(&thread, 0, iterate,(void*)m);
    pthread_detach(thread);
  }
  pthread_t hlthread;
  pthread_create(&hlthread,0,scanHighlight,(void*)m);
  pthread_detach(hlthread);

  cam=new Cam();
  glutMainLoop();
  return 0;
//p=new Poly(m);

}
