#include "cgd.hpp"
#include "mat.hpp"
#include "sparseccs.h"
#include "sparseMatrix.h"
#include "CGSolver.h"
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
    plane[ii].n/=cnt[ii];
    plane[ii].n/=plane[ii].n.norm();
  }
}
#include <fstream>
void add_v4(Mesh & m )
{
  m.v4.resize(m.t.size());
  for (size_t ii=0; ii<m.t.size(); ii++) {
    Vec3 a = m.v[m.t[ii][1]] -  m.v[m.t[ii][0]];
    Vec3 b = m.v[m.t[ii][2]] -  m.v[m.t[ii][0]];
    Vec3 n = a.cross(b);
    double area = n.norm();
    n/=sqrt(area);
    m.v4[ii]=m.v[m.t[ii][0]] + n;
    m.t[ii].x[3]=m.v.size()+ii;
  }
}

void vmat(Mesh & m, std::vector<Mat3>&A)
{
  A.resize(m.t.size(),Mat3::Zero);
  for (size_t ii=0; ii<m.t.size(); ii++) {
    Mat3 mat(m.v[m.t[ii][1]] - m.v[m.t[ii][0]],
             m.v[m.t[ii][2]] - m.v[m.t[ii][0]],
             m.v4[ii]        - m.v[m.t[ii][0]]);
    mat.gauss_pivot(A[ii]);
  }
}
/**@brief compressed column storage*/
struct CCS {
public:
  std::vector<real>vals;
  std::vector<int>rowInds;
  std::vector<int>colPtr;
  CCS() {
    colPtr.push_back(0);
  }
};

#include <stdio.h>
/**@brief rowVal although it is a map,
it is used as a binary heap here*/
void addrow(std::map<int, real>&rowVal, CCS&ccs, int row_offset=0)
{
  std::map<int,real>::iterator it;
  int nCols=ccs.colPtr.size()-1;
  ccs.colPtr.push_back(ccs.colPtr[nCols]+rowVal.size());

  for(it = rowVal.begin(); it!=rowVal.end(); it++) {
    ccs.rowInds.push_back(it->first + row_offset);
    ccs.vals.push_back(it->second);
  }
}

/**
  Creates a column in the matrix only for the first axis.
  The calling function would then create two more copies
  for the other two axis
  @param _ii row number for m
*/
void add_coef(Trig& t,int _ii, Mat3 & m,
              std::map<int,real> &val)
{
  float sum=0;
  for(int i=0; i<3; i++) {
    sum -= m._m[i][_ii];
  }
  val[t[0]] += sum;
  for(int i=0; i<3; i++) {
    val[t[i+1]] += m._m[i][_ii];
  }
}

void scale(std::map<int,real>&a, float c )
{
  std::map<int,real>::iterator it;
  for(it = a.begin(); it!=a.end(); it++) {
    it->second*=c;
  }
}
/**@brief result is added to a*/
void add_times(std::map<int,real>&a, std::map<int,real>&b,float c)
{
  std::map<int,real>::iterator it;
  for(it = b.begin(); it!=b.end(); it++) {
    a[it->first]+=c * it->second;
  }
}
int nvar;
/**@brief neighboring transformations
should be similar
@param tIdx index of the triangle being proecssed
@param mat vector of transformations for each triangle

assumes m already contains adjacency matrix
*/
void smooth_mat(Mesh & m,
                std::vector<Mat3>&mat,
                float wS, CCS&ccs)
{
  //each triangle has 9 constraints for
  //each number in the transformation matrix
for(int axis=0; axis<3; axis++) {

    for(size_t tIdx=0; tIdx<m.t.size(); tIdx++) {

      for(int ii=0; ii<3; ii++) {
        std::map<int,real>rowVal;
        std::map<int,real>neiborRowVal;
        add_coef(m.t[tIdx],ii,mat[tIdx],rowVal);
        for(size_t nbr =0; nbr<m.adjMat[tIdx].size(); nbr++) {
          int nbrIdx = m.adjMat[tIdx][nbr];
          add_coef(m.t[nbrIdx],ii,mat[nbrIdx],neiborRowVal);
        }
        add_times(rowVal, neiborRowVal,-(1.0/m.adjMat[tIdx].size()));
        scale(rowVal, wS);

        addrow(rowVal, ccs, nvar*axis);
      }
    }
  }
}

#define VAR_IDX(i) ((axis)*nvar+(i))

void planar_mat(Mesh & m, std::vector<Plane>&plane, float wPt, CCS& ccs)
{
//  int nvar = m.t.size()+ m.v.size();
for(size_t tIdx=0; tIdx<m.t.size(); tIdx++) {

  Plane & p = plane[m.t[tIdx].label];
  for(int row = 0; row<2; row++) {
    std::map<int, real> val;
    for(int axis=0; axis<3; axis++) {
      int idx= VAR_IDX(m.t[tIdx][row]);
      val[idx]=p.n.x[axis]   * wPt;

      idx= VAR_IDX(m.t[tIdx][row+1]);
      val[idx]=-p.n.x[axis]  * wPt;
    }
    addrow(val,ccs);
  }
  }
}
/**@brief
  @param b will be filled with 1 0 0 0 1 0 0 0 1
*/
void identity_mat(Mesh & m , std::vector<Mat3> &mat, float wI, CCS&ccs, double * b)
{
  int idx=0;
 // int nvar = m.v.size()+m.t.size();
  for(int axis=0; axis<3; axis++) {

  for(size_t ii=0; ii<m.t.size(); ii++) {
    for(int jj=0; jj<3; jj++) {

        if(jj==axis) {
          b[idx]=wI;
        } else {
          b[idx]=0;
        }
        idx++;
      }
    }
  }
  for(int axis=0; axis<3; axis++) {

  for(size_t ii=0; ii<m.t.size(); ii++) {
    for(int jj=0; jj<3; jj++) {
      std::map<int,real>val;
      add_coef(m.t[ii],jj,mat[ii],val);
      scale(val,wI);

      addrow(val, ccs, nvar*axis);
    }
    }
  }
}

int find(int * a,int b,  size_t size)
{
  for(size_t ii=0; ii<size; ii++) {
    if(a[ii]==b) {
      return ii;
    }
  }
  return -1;
}

/**@brief transformed vertices should stay close to
original vertices*/
#include <set>
float vW=30;
void vertex_mat(Mesh &m , float vW0, CCS& ccs , double * bb)
{
  //weight for vertices at the intersection of at least three clusters
  //weight for regular vertices is 1

  //maximum number of clusters adjacent to a vertex
  size_t maxLabel=3;

  std::vector<bool>processed (m.v.size());
  std::vector<bool>important(m.v.size());

  for(size_t ii=0; ii<m.t.size(); ii++) {
    for(int vidx = 0; vidx<3; vidx++) {
      int idx = m.t[ii][vidx];
      if(processed[idx]) {
        continue;
      }
      std::set<int>labels;
      processed[idx]=1;
      labels.insert(m.t[ii].label);
      bool found = 1;
      size_t prev = ii;
      size_t cur = ii;
      while(found) {
        found =0 ;
        for (size_t nbr = 0; nbr<m.adjMat[ii].size(); nbr++) {
          size_t nbrIdx = m.adjMat[cur][nbr];
          if(nbrIdx==prev || nbrIdx == ii) {
            continue;
          }
          if(find(m.t[nbrIdx].x,idx,3)>=0) {
            found=1;
            prev=cur;
            cur=nbrIdx;
            break;
          }
        }
        if(found) {
          int lab = m.t[cur].label;
          labels.insert(lab);
          if(labels.size()>=maxLabel) {
            break;
          }
        }
      }

      if(labels.size()>=maxLabel) {
        important[idx]=true;
      }

    }
  }
  int bidx=0;
  for(int axis=0; axis<3; axis++) {

  for(size_t idx=0;idx<m.v.size();idx++){
    real f;
    std::map<int,real>val;
    if(important[idx]){
      f=vW*vW0;
    } else {
       f=vW0;
    }
    val[idx]=f;
    addrow(val,ccs,axis*nvar);
    bb[bidx]=f*m.v0[idx][axis];
    bidx++;
  }
  }
}

void vertex2arr(const Mesh & m, double * x)
{
  int idx=0;
  for(int axis=0; axis<3; axis++) {
    for(size_t ii=0; ii<m.v.size(); ii++) {
      x[idx]=m.v[ii].get(axis);
      idx++;
    }
    for(size_t ii=0; ii<m.v4.size(); ii++) {
      x[idx]=m.v4[ii].get(axis);
      idx++;
    }
  }
}

void array2vertex(const double * x, Mesh &m)
{
  int idx=0;
  std::vector<Vec3> v0 = m.v;
  for(int axis=0; axis<3; axis++) {
    for(size_t ii=0; ii<m.v.size(); ii++) {
      m.v[ii][axis] = x[idx];
      idx++;
    }

    for(size_t ii=0; ii<m.v4.size(); ii++) {
      m.v4[ii][axis]=x[idx];
      idx++;
    }
  }
  m.self_intersect();
  std::map<int,bool>::iterator it;
  for(it = m.bad.begin();it!=m.bad.end();it++){
    for(int ii=0;ii<3;ii++){
      int vidx=m.t[it->first][ii];
      m.v[vidx]=v0[vidx];
    }
  }
}

void printAB(CCS&ccs, std::vector<double > & b, const double *x)
{
  std::ofstream out;
  out.open("AB.txt");
  //out<<ccs.vals.size()<<" "<<(ccs.colPtr.size()-1)<<"\n";
  for(size_t ii=0; ii<ccs.colPtr.size()-1; ii++) {
    for(int jj=ccs.colPtr[ii]; jj<ccs.colPtr[ii+1]; jj++) {
      int row = ccs.rowInds[jj];
      real val = ccs.vals[jj];
      out<<(row+1)<<" "<<(ii+1)<<" "<<val<<"\n";
    }
  }
  out<<"B\n";//<<b.size()<<"\n";
  for(size_t ii=0; ii<b.size(); ii++) {
    out<<b[ii]<<"\n";
  }
  out.close();
}

  float wS=0.1;
  float wI=0.1;
  float wV0=1;
  float wPt=2;

//static int iter=0;
void cgd(Mesh & m)
{
  int height=20*m.t.size()+3*m.v.size();
//  int height=2*m.t.size()+3*m.v.size();
//  int width = 3*(m.v.size());
  int width = 3*(m.v.size()+m.t.size());
  //nvar=m.v.size();
  nvar=m.v.size()+m.t.size();
  std::vector<Plane>plane;
  get_plane(m,plane);
  add_v4(m);
  std::vector<Mat3>mat;
  vmat(m,mat);
  m.adjlist();
  CCS ccs;
  std::vector<double >b(height, 0.0);
  //b for smoothness matrix is 0
  smooth_mat(m,mat,wS,ccs);
  identity_mat(m, mat, wI, ccs , &b[9*m.t.size()] );
  planar_mat(m,plane,wPt,ccs);
  vertex_mat(m,wV0,ccs,&b[20*m.t.size()]);
//  vertex_mat(m,wV0,ccs,&b[2*m.t.size()]);

  SparseCCS* s = new SparseCCS(width, height, &ccs.vals[0], &ccs.rowInds[0], &ccs.colPtr[0]);
  SparseCCS* st = s->transposed();
  SparseCCS* sst=s->multiply_LM(*st);
  SparseMatrixOutline *outline = new SparseMatrixOutline(width);
printAB(ccs,b,0);
  double * ATb= s->operator* (&b[0]);

//  int rowcnt=sst->rows_count();
  int colcnt=sst->columns_count();
  const real * vals=sst->values();
  const int * row_indices=sst->row_indices();
  const int * col_pointers=sst->column_pointers();
  for(int ii=0; ii<colcnt; ii++) {
    for(int jj=col_pointers[ii]; jj<col_pointers[ii+1]; jj++) {
      outline->AddEntry(row_indices[jj],ii,vals[jj]);
    }
  }
  // delete s;
  delete st;
  delete sst;

  SparseMatrix A(outline);
  delete outline;
  CGSolver solver(&A);
  double * x=new double [width];
  vertex2arr(m,x);

  double eps = 1E-9;
  int maxIter = 500;
  int verbose = 0;
  int ret = solver.SolveLinearSystemWithJacobiPreconditioner(x, ATb, eps, maxIter, verbose);//
  if(ret<0) {
    printf("optimization error\n");
  }
  A.CheckLinearSystemSolution(x,ATb);
  printf("\n");
  array2vertex(x,m);

  delete []x;
  delete []ATb;
}

void weighted_avg(Mesh& m)
{
  std::vector<int>nbrCnt(m.v.size());
  std::vector<Vec3>nbrPos(m.v.size());
  for(size_t ii=0;ii<m.v.size();ii++){
    Vec3 sum;
    for(int jj=0;jj<3;jj++){
      sum += m.v[m.t[ii][jj]];
    }
    for(int jj=0;jj<3;jj++){
      nbrPos[m.t[ii][jj]] += (sum-m.v[m.t[ii][jj]]);
      nbrCnt[m.t[ii][jj]]+=2;
    }
  }
  std::vector<Vec3> v0 = m.v;
  for(size_t ii=0;ii<m.v.size();ii++){
    nbrPos[ii]/=nbrCnt[ii];
  }
}
