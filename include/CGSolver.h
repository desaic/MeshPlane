/*************************************************************************
 *                                                                       *
 * "sparseMatrix" library , Copyright (C) 2007 CMU, 2009 MIT             *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code author: Jernej Barbic                                            *
 * http://www.jernejbarbic.com/code                                      *
 * Research: Jernej Barbic, Doug L. James, Jovan Popovic                 *
 * Funding: NSF, Link Foundation, Singapore-MIT GAMBIT Game Lab          *
 * Version 2.0                                                           *
 *                                                                       *
 * This library is free software; you can redistribute it and/or         *
 * modify it under the terms of the BSD-style license that is            *
 * included with this library in the file LICENSE.txt                    *
 *                                                                       *
 * This library is distributed in the hope that it will be useful,       *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the file     *
 * LICENSE.TXT for more details.                                         *
 *                                                                       *
 *************************************************************************/

/*
  A conjugate gradient solver built on top of the sparse matrix class.
  There are two solver versions: without preconditioning, and with 
  Jacobi preconditioning.

  You can either provide a sparse matrix, or a callback function to
  multiply x |--> A * x .

  The sparse matrix must be symmetric and positive-definite.

  The CG solvers were implemented by following Jonathan Shewchuk's 
  An Introduction to the Conjugate Gradient Method Without the Agonizing Pain:
  http://www.cs.cmu.edu/~jrs/jrspapers.html#cg

  See also sparseMatrix.h .

  See example.cpp for example usage.
*/

#ifndef _CGSOLVER_H_
#define _CGSOLVER_H_

#include "sparseMatrix.h"

class CGSolver
{
public:

  // standard constructor
  CGSolver(SparseMatrix * A);

  // this constructor version makes it possible to only provide a
  // "black-box" mtx-vec multiplication routine (no need to explicitly give the matrix):
  // given x, the routine must compute A * x, and store it into Ax
  // "data" should not be used/touched by custom user routines
  typedef void (*blackBoxProductType)(const void * data, const double * x, double * Ax);
  CGSolver(int n, blackBoxProductType callBackFunction, void * data);

  ~CGSolver();

  // solves A * x = b
  // x should (on input) contain the initial guess
  // on output, x will be modified to the new value (the solution)
  // 0 < eps < 1 is the error tolerance: solver converges when the L2 residual errors is less than eps the initial L2 residual error
  // maxIter is the max number of CG iterations
  // return value is the number of iterations performed 
  // if solver did not converge, the return value will be given a negative sign
  int SolveLinearSystem(double * x, const double * b, double eps=1e-6, int maxIter=1000, int verbose=0);

  // same as above, except it uses Jacobi preconditioning
  // the employed error metric is M^{-1}-weighted L2 residual error (see Shewchuk)
  int SolveLinearSystemWithJacobiPreconditioner(double * x, const double * b, double eps, int maxIter, int verbose=0);

  double DotProduct(double * x, double * y); // length of vectors x,y equals n (dimension of A)

protected:
  int n;
  blackBoxProductType multiplicator;
  void * multiplicatorData;
  SparseMatrix * A; 
  double * r, * d, * q;
  double * invDiagonal;

  double TriVectorProduct(double * x, double * y, double * z);

  static void DefaultMultiplicator(const void * data, const double * x, double * Ax);

  void InitBuffers();
};

#endif

