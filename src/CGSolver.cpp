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

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "CGSolver.h"

CGSolver::CGSolver(SparseMatrix * A_): A(A_) 
{
  n = A->Getn();
  InitBuffers();
  multiplicator = CGSolver::DefaultMultiplicator;
  multiplicatorData = (void*)A;
}

CGSolver::CGSolver(int n_, blackBoxProductType callBackFunction_, void * data_): n(n_), multiplicator(callBackFunction_), multiplicatorData(data_), A(NULL)
{
  InitBuffers();
}

CGSolver::~CGSolver()
{
  free(r);
  free(d);
  free(q);
  free(invDiagonal);
}

void CGSolver::DefaultMultiplicator(const void * data, const double * x, double * Ax)
{
  SparseMatrix * A = (SparseMatrix*)data;
  A->MultiplyVector(x, Ax);
}

void CGSolver::InitBuffers()
{
  r = (double*) malloc (sizeof(double) * n);
  d = (double*) malloc (sizeof(double) * n);
  q = (double*) malloc (sizeof(double) * n);
  invDiagonal = (double*) malloc (sizeof(double) * n);
}

double CGSolver::DotProduct(double * x, double * y)
{
  double result = 0;
  for(int i=0; i<n; i++)
    result += x[i] * y[i];

  return result;
}

double CGSolver::TriVectorProduct(double * x, double * y, double * z)
{
  double result = 0;
  for(int i=0; i<n; i++)
    result += x[i] * y[i] * z[i];

  return result;
}

int CGSolver::SolveLinearSystem(double * x, const double * b, double eps, int maxIter, int verbose)
{
  int iter=1;
  multiplicator(multiplicatorData, x, r); //A->MultiplyVector(x,r);
  for (int i=0; i<n; i++)
  {
    r[i] = b[i] - r[i];
    d[i] = r[i];
  }

  double deltaNew = DotProduct(r,r);
  double delta0 = deltaNew;

  while ((iter <= maxIter) && (deltaNew > eps * eps * delta0))
  {
    if (verbose)
    {
      printf("CG iteration %d: current L2 error vs initial error=%G\n", iter, sqrt(deltaNew/delta0));
    }

    multiplicator(multiplicatorData, d, q); //A->MultiplyVector(d,q); // q = A * d
    double dtq = DotProduct(d,q);
    double alpha = deltaNew / dtq;
    //printf("deltaNew=%G dtq=%G alpha=%G\n", deltaNew, dtq, alpha);

    for(int i=0; i<n; i++)
      x[i] += alpha * d[i];

    if (iter % 50 == 0)
    {
      // periodically compute the exact residual
      multiplicator(multiplicatorData, x, r); //A->MultiplyVector(x,r);
      for (int i=0; i<n; i++)
        r[i] = b[i] - r[i];
    }
    else
    {
      for (int i=0; i<n; i++)
        r[i] = r[i] - alpha * q[i];
    }

    double deltaOld = deltaNew;
    deltaNew = DotProduct(r,r);
    double beta = deltaNew / deltaOld;

    for (int i=0; i<n; i++)
      d[i] = r[i] + beta * d[i];

    iter++;
  }

  return (iter-1) * ((deltaNew > eps * eps * delta0) ? -1 : 1);
}

int CGSolver::SolveLinearSystemWithJacobiPreconditioner(double * x, const double * b, double eps, int maxIter, int verbose)
{
  // extract diagonal entries
  A->BuildDiagonalIndices(); // note: if indices are already built, this call will do nothing (you can therefore also call BuildDiagonalIndices() once and for all before calling SolveLinearSystemWithJacobiPreconditioner); in any case, BuildDiagonalIndices() is fast (a single linear traversal of all matrix elements)

  A->GetDiagonal(invDiagonal);
  for(int i=0; i<n; i++)
    invDiagonal[i] = 1.0 / invDiagonal[i];

  int iter=1;
  multiplicator(multiplicatorData, x, r); //A->MultiplyVector(x,r);
  for (int i=0; i<n; i++)
  {
    r[i] = b[i] - r[i];
    d[i] = invDiagonal[i] * r[i];
  }

  double deltaNew = TriVectorProduct(r,r,invDiagonal);
  double delta0 = deltaNew;

  while ((iter <= maxIter) && (deltaNew > eps * eps * delta0))
  {
    if (verbose)
    {
      printf("CG iteration %d: current M^{-1}-L2 error vs initial error=%G\n", iter, sqrt(deltaNew/delta0));
    }

    multiplicator(multiplicatorData, d, q); //A->MultiplyVector(d,q); // q = A * d
    double dtq = DotProduct(d,q);
    double alpha = deltaNew / dtq;

    for(int i=0; i<n; i++)
      x[i] += alpha * d[i];

    if (iter % 50 == 0)
    {
      // periodically compute the exact residual
      multiplicator(multiplicatorData, x, r); //A->MultiplyVector(x,r);
      for (int i=0; i<n; i++)
        r[i] = b[i] - r[i];
    }
    else
    {
      for (int i=0; i<n; i++)
        r[i] = r[i] - alpha * q[i];
    }

    double deltaOld = deltaNew;
    deltaNew = TriVectorProduct(r,r,invDiagonal);
    double beta = deltaNew / deltaOld;

    for (int i=0; i<n; i++)
      d[i] = invDiagonal[i] * r[i] + beta * d[i];

    iter++;
  }

  if (deltaNew < 0)
  {
    printf("Warning: deltaNew=%G is negative. Input matrix might not be SPD. Solution could be incorrect.\n", deltaNew);
  }

  return (iter-1) * ((deltaNew > eps * eps * delta0) ? -1 : 1);
}

