#ifndef CGD_HPP
#define CGD_HPP
#include "mesh.hpp"
#include "sparseMatrix.h"
#include "CGSolver.h"
extern  float wS;
extern  float wI;
extern  float wV0;
extern  float wPt,wP,wN,wV,w0;
extern float vW;
void cgd(Mesh & m);
void weighted_avg(Mesh& m);
#endif
