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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <set>
#include <algorithm>
#include "sparseMatrix.h"
using namespace std;

SparseMatrixOutline::SparseMatrixOutline(int n): n_(n)
{
  MakeEmpty();
}

SparseMatrixOutline::~SparseMatrixOutline()
{
  // deallocate individual row maps
  for(int i=0; i<n_; i++)
    entries_[i].clear();

  // deallocate vector entries_
  entries_.clear();
}

SparseMatrixOutline::SparseMatrixOutline(int n, double diagonal): n_(n)
{
  MakeEmpty();

  pair<int,double> entry;

  for(int i=0; i<n_; i++)
  {
    entry.first = i;
    entry.second = diagonal;
    entries_[i].insert(entry);
  }
}

SparseMatrixOutline::SparseMatrixOutline(int n, double * diagonal): n_(n)
{
  MakeEmpty();

  pair<int,double> entry;

  for(int i=0; i<n_; i++)
  {
    entry.first = i;
    entry.second = diagonal[i];
    entries_[i].insert(entry);
  }
}

SparseMatrixOutline::SparseMatrixOutline(char * filename, int expand)
{
  if (expand <= 0)
  {
    printf("Error: invalid expand factor %d in SparseMatrixOutline constructor.\n", expand);
    throw 1;
  }
  
  FILE * inputMatrix = fopen(filename,"ra");
  if (!inputMatrix)
  {
    printf("Error: couldn't open input sparse matrix file %s.\n",filename);
    throw 2;
  }

  // read input size 
  int m,n;
  if(fscanf(inputMatrix,"%d\n%d\n",&m,&n)==0){
    printf("error reading file\n");
    
  }

  n_ = expand * m;

  printf("Loading matrix from %s... Size is %d x %d .\n",filename,n_,expand*n);fflush(NULL);

  MakeEmpty();

  char s[4096];
  while (fgets(s,4096,inputMatrix) != NULL)
  {
    int i1,j1;
    double x;
    sscanf(s, "%d %d %lf\n", &i1, &j1, &x);
    for(int e=0; e<expand; e++)
      AddEntry(expand * i1 + e, expand * j1 + e, x);
  }

  fclose(inputMatrix);
}

void SparseMatrixOutline::MakeEmpty()
{
  // allocate empty lists, one for very row
  entries_.clear();
  map<int,double> emptyMap;
  for (int i=0; i<n_; i++)
    entries_.push_back(emptyMap);
}

void SparseMatrixOutline::AddEntry(int i, int j, double value)
{
  map<int,double>::iterator pos = entries_[i].find(j);
  if (pos != entries_[i].end())
    pos->second += value;
  else
  {
    pair<int,double> entry(j,value);
    entries_[i].insert(entry);
  }
}

void SparseMatrixOutline::MultiplyRow(int row, double scalar)
{
  for(map<int,double>::iterator iter = entries_[row].begin(); iter != entries_[row].end(); iter++)
    iter->second *= scalar;
}

void SparseMatrixOutline::AddBlock3x3Entry(int i, int j, double * matrix3x3)
{
  for(int k=0; k<3; k++)
    for(int l=0; l<3; l++)
      AddEntry(3*i+k,3*j+l,matrix3x3[3*k+l]);
}

double SparseMatrixOutline::GetEntry(int i, int j) const
{
  map<int,double>::const_iterator pos = entries_[i].find(j);
  if (pos != entries_[i].end())
    return (pos->second);
  else
    return 0;
}

int SparseMatrixOutline::GetNumColumns() const
{
  int numColumns = -1;
  for(int i=0; i< n_; i++)
  {
    map<int,double>::const_iterator j1;
    // over all row entries 
    for(j1=entries_[i].begin(); j1 != entries_[i].end(); j1++) 
      if (j1->first > numColumns)
        numColumns = j1->first;
  }
  return numColumns + 1;
}

int SparseMatrixOutline::Save(char * filename, int oneIndexed) const
{
  FILE * fout = fopen(filename,"wa");
  if (!fout)
    return 1;

  fprintf(fout,"%d\n%d\n",n_,GetNumColumns());
  for(int i=0; i< n_; i++)
  {
    map<int,double>::const_iterator j1;
    // over all row entries 
    for(j1=entries_[i].begin(); j1 != entries_[i].end(); ++j1) 
      fprintf(fout,"%d %d %.15f\n",i,j1->first + oneIndexed, j1->second + oneIndexed);

  }
  fclose(fout);

  return 0;
}

void SparseMatrixOutline::Print() const
{
  for (int i=0; i<n_; i++)
  {
    for (int j=0; j<n_; j++)
      printf("%f ",GetEntry(i,j));
    printf("\n");
  }
}

int SparseMatrixOutline::GetNumEntries() const
{
  int num = 0;
  for(int i=0; i<n_; i++)
    num += entries_[i].size();
  return num;
}

SparseMatrix::SparseMatrix(char * filename)
{
  SparseMatrixOutline sparseMatrixOutline(filename);
  InitFromOutline(&sparseMatrixOutline);
}

SparseMatrix::SparseMatrix(SparseMatrixOutline * sparseMatrixOutline)
{
  InitFromOutline(sparseMatrixOutline);
}

// construct matrix from the outline
void SparseMatrix::InitFromOutline(SparseMatrixOutline * sparseMatrixOutline)
{
  n_ = sparseMatrixOutline->Getn();
  Allocate();

  for(int i=0; i<n_; i++)
  {
    rowLength_[i] = sparseMatrixOutline->entries_[i].size();
    cindices_[i] = (int*) malloc (sizeof(int) * rowLength_[i]);
    centries_[i] = (double*) malloc (sizeof(double) * rowLength_[i]);

    map<int,double>::iterator pos;
    int j = 0;
    int prev = -1;
    for(pos = sparseMatrixOutline->entries_[i].begin(); pos != sparseMatrixOutline->entries_[i].end(); pos++)
    {
      cindices_[i][j] = pos->first;
      if (cindices_[i][j] <= prev)
        printf("Warning: entries not sorted in a row in a sparse matrix.\n");
      prev = cindices_[i][j];
      centries_[i][j] = pos->second;
      j++;
    }
  }
}

// allocator
void SparseMatrix::Allocate()
{
  // compact representation for fast matrix multiplications
  rowLength_ = (int*) malloc(sizeof(int) * n_);
  cindices_ = (int**) malloc(sizeof(int*) * n_);
  centries_ = (double**) malloc(sizeof(double*) * n_);
  diagonalIndices_ = NULL;
  subMatrixIndices_ = NULL;
  subMatrixIndexLengths_ = NULL;
  superMatrixIndices_ = NULL;
  superRows_ = NULL;
  transposeIndex_ = NULL;
}

// destructor
SparseMatrix::~SparseMatrix()
{
  for(int i=0; i<n_; i++)
  {
    free(cindices_[i]);
    free(centries_[i]);
  }

  if (subMatrixIndices_ != NULL)
  {
    for(int i=0; i<n_; i++)
      free(subMatrixIndices_[i]);
    free(subMatrixIndices_);
    free(subMatrixIndexLengths_);
  }

  if (superRows_ != NULL)
  {
    for(int i=0; i<n_; i++)
      free(superMatrixIndices_[i]);
    free(superMatrixIndices_);
    free(superRows_);
  }

  free(rowLength_);
  free(cindices_);
  free(centries_);
  free(diagonalIndices_);
  FreeTranspositionIndices();
}

SparseMatrix::SparseMatrix(const SparseMatrix & source) // copy constructor
{
  //printf("Copy constructor:\n");fflush(NULL);
  n_ = source.Getn();

  // compact representation for fast matrix-vtx multiplications
  rowLength_ = (int*) malloc(sizeof(int) * n_);
  cindices_ = (int**) malloc(sizeof(int*) * n_);
  centries_ = (double**) malloc(sizeof(double*) * n_);


  for(int i=0; i<n_; i++)
  {
    rowLength_[i] = source.rowLength_[i];
    cindices_[i] = (int*) malloc (sizeof(int) * rowLength_[i]);
    centries_[i] = (double*) malloc (sizeof(double) * rowLength_[i]);

    for(int j=0; j < rowLength_[i]; j++)
    {
      cindices_[i][j] = source.cindices_[i][j];
      centries_[i][j] = source.centries_[i][j];
    }
  }

  subMatrixIndices_ = NULL; 
  subMatrixIndexLengths_ = NULL;
  if (source.subMatrixIndices_ != NULL)
  {
    subMatrixIndices_ = (int**) malloc(sizeof(int*) * n_);
    subMatrixIndexLengths_ = (int*) malloc(sizeof(int) * n_);

    for(int i=0; i<n_; i++)
    {
      subMatrixIndexLengths_[i] = source.subMatrixIndexLengths_[i];
      subMatrixIndices_[i] = (int*) malloc(sizeof(int) * subMatrixIndexLengths_[i]);
      for(int j=0; j < subMatrixIndexLengths_[i]; j++)
      {
        subMatrixIndices_[i][j] = source.subMatrixIndices_[i][j];
      }
    }
  }

  superRows_ = NULL;
  superMatrixIndices_ = NULL;
  if (source.superRows_ != NULL)
  {
    superRows_ = (int*) malloc(sizeof(int) * n_);
    superMatrixIndices_ = (int**) malloc(sizeof(int*) * n_);
    for(int i=0; i<n_; i++)
    {
      superRows_[i] = source.superRows_[i];
      superMatrixIndices_[i] = (int*) malloc(sizeof(int) * rowLength_[i]);
      for(int j=0; j < rowLength_[i]; j++)
      {
        superMatrixIndices_[i][j] = source.superMatrixIndices_[i][j];
      }
    }
  }

  diagonalIndices_ = NULL; 
  if (source.diagonalIndices_ != NULL)
  {
    diagonalIndices_ = (int*) malloc (sizeof(int) * n_);
    memcpy(diagonalIndices_, source.diagonalIndices_, sizeof(int) * n_);
  }

  transposeIndex_ = NULL;
  if (source.transposeIndex_ != NULL)
  {
    transposeIndex_ = (int**) malloc (sizeof(int*) * n_);
    for(int i=0; i<n_; i++)
    {
      transposeIndex_[i] = (int*) malloc (sizeof(int) * rowLength_[i]);
      for(int j=0; j<rowLength_[i]; j++)
        transposeIndex_[i][j] = source.transposeIndex_[i][j];
    }
  }
}

void SparseMatrix::MultiplyVector(const double * input, double * result) const
{
  for(int i=0; i<n_; i++)
  {
    result[i] = 0;
    for(int j=0; j < rowLength_[i]; j++)
    {
      result[i] += input[cindices_[i][j]] * centries_[i][j];
    }
  }
}

void SparseMatrix::MultiplyVectorAdd(const double * input, double * result) const
{
  for(int i=0; i<n_; i++)
  {
    for(int j=0; j < rowLength_[i]; j++)
    {
      result[i] += input[cindices_[i][j]] * centries_[i][j];
    }
  }
}

void SparseMatrix::TransposeMultiplyVector(const double * input, double * result, int resultLength) const
{
  for(int i=0; i<resultLength; i++)
    result[i] = 0;

  for(int i=0; i<n_; i++)
  {
    for(int j=0; j < rowLength_[i]; j++)
    {
      result[cindices_[i][j]] += input[i] * centries_[i][j];
    }
  }
}

void SparseMatrix::TransposeMultiplyVectorAdd(const double * input, double * result) const
{
  for(int i=0; i<n_; i++)
  {
    for(int j=0; j < rowLength_[i]; j++)
    {
      result[cindices_[i][j]] += input[i] * centries_[i][j];
    }
  }
}

void SparseMatrix::MultiplyMatrix(int numRows, int numColumns, const double * input, double * result) const
{
  for(int column=0; column<numColumns; column++)
    MultiplyVector(&input[numRows * column], &result[n_ * column]);
}

void SparseMatrix::MultiplyMatrixAdd(const double * input, int numColumns, double * result) const
{
  for(int column=0; column<numColumns; column++)
    MultiplyVectorAdd(&input[n_ * column], &result[n_ * column]);
}

double SparseMatrix::QuadraticForm(const double * input) const
{
  double result = 0;
 
  for(int i=0; i<n_; i++)
  {
    for(int j=0; j < rowLength_[i]; j++)
    {
      int index = cindices_[i][j];
      if (index < i)
        continue;
      if (index == i)
        result += centries_[i][j] * input[i] * input[index];
      else
        result += 2.0 * centries_[i][j] * input[i] * input[index];
    }
  }
  
  return result;
} 

void SparseMatrix::NormalizeVector(double * vector) const
{
  double norm = sqrt(QuadraticForm(vector));
  for(int i=0; i<n_; i++)
    vector[i] /= norm;
}

SparseMatrix SparseMatrix::operator+(const SparseMatrix & mat2) const
{
  SparseMatrix result(*this);
  for(int i=0; i<n_; i++)
    for(int j=0; j < rowLength_[i]; j++)
      result.centries_[i][j] += mat2.centries_[i][j];
  return result;
}

SparseMatrix SparseMatrix::operator-(const SparseMatrix & mat2) const
{
  SparseMatrix result(*this);
  for(int i=0; i<n_; i++)
    for(int j=0; j < rowLength_[i]; j++)
      result.centries_[i][j] -= mat2.centries_[i][j];
  return result;
}

SparseMatrix operator* (const double alpha, const SparseMatrix & mat2)
{
  SparseMatrix result(mat2);
  for(int i=0; i<result.n_; i++)
    for(int j=0; j < result.rowLength_[i]; j++)
      result.centries_[i][j] *= alpha;
  return result;
}

SparseMatrix & SparseMatrix::operator*=(const double alpha)
{
  for(int i=0; i<n_; i++)
    for(int j=0; j < rowLength_[i]; j++)
      centries_[i][j] *= alpha;
  return *this;
}

SparseMatrix & SparseMatrix::operator+=(const SparseMatrix & mat2)
{   
  for(int i=0; i<n_; i++)
    for(int j=0; j < rowLength_[i]; j++)
      centries_[i][j] += mat2.centries_[i][j];
  return *this;
}

SparseMatrix & SparseMatrix::operator-=(const SparseMatrix & mat2)
{  
  for(int i=0; i<n_; i++)
    for(int j=0; j < rowLength_[i]; j++)
      centries_[i][j] -= mat2.centries_[i][j]; 
  return *this;
}

SparseMatrix & SparseMatrix::operator=(const SparseMatrix & source)
{
  for(int i=0; i<n_; i++)
  {
    for(int j=0; j < rowLength_[i]; j++)
      centries_[i][j] = source.centries_[i][j];
  }

  return *this;
}

void SparseMatrix::ScalarMultiply(const double alpha, SparseMatrix * dest)
{
  if (dest == NULL)
    dest = this;

  for(int i=0; i<n_; i++)
    for(int j=0; j < rowLength_[i]; j++)
      dest->centries_[i][j] = centries_[i][j] * alpha;
}

void SparseMatrix::ResetToZero()
{
  for(int i=0; i<n_; i++)
    memset(centries_[i], 0, sizeof(double) * rowLength_[i]);
}

void SparseMatrix::Print(int sparsePrint) const
{
  if (sparsePrint)
  {
    for (int i=0; i<n_; i++)
      for(int j=0; j< rowLength_[i]; j++)
        printf("%d %d %G\n", i, cindices_[i][j], centries_[i][j]);
  }
  else
  {
    for (int i=0; i<n_; i++)
    {
      int index = 0;
      for(int j=0; j< rowLength_[i]; j++)
      {
        while (index < cindices_[i][j])
        {
          index++;
          printf("%f,",0.0);
        }
        printf("%f,",centries_[i][j]);
        index++;
      }

      while (index < n_)
      {
        index++;
        printf("%f,",0.0);
      }
    
      printf("\n");
    } 
  }
}

// finds the position in row i of element with column index jAbsolute
int SparseMatrix::GetInverseIndex(int i, int jAbsolute) const
{
  for(int j=0; j < rowLength_[i]; j++)
    if (cindices_[i][j] == jAbsolute)
      return j;

  return -1;
}

void SparseMatrix::BuildDiagonalIndices()
{
  if (diagonalIndices_ != NULL)
    return;

  diagonalIndices_ = (int*) malloc (sizeof(int) * n_);
  for(int i=0; i<n_; i++)
    diagonalIndices_[i] = GetInverseIndex(i,i);
}

void SparseMatrix::FreeDiagonalIndices()
{
  free(diagonalIndices_);
}

void SparseMatrix::GetDiagonal(double * diagonal)
{
  if (diagonalIndices_ != NULL)
  {
    for(int i=0; i<n_; i++)
      diagonal[i] = centries_[i][diagonalIndices_[i]];
  }
  else
  {
    for(int i=0; i<n_; i++)
      for(int j=0; j<GetRowLength(i); j++)
      {
        if (GetColumnIndex(i, j) == i)
          diagonal[i] = centries_[i][j];
      }
  }
}

void SparseMatrix::AddDiagonalMatrix(double * diagonalMatrix)
{
  if (diagonalIndices_ != NULL)
  {
    for(int i=0; i<n_; i++)
      centries_[i][diagonalIndices_[i]] += diagonalMatrix[i];
  }
  else
  {
    for(int i=0; i<n_; i++)
      for(int j=0; j<GetRowLength(i); j++)
      {
        if (GetColumnIndex(i, j) == i)
          centries_[i][j] += diagonalMatrix[i];
      }
  }
}

int SparseMatrix::GetNumEntries() const
{
  int num = 0;
  for(int i=0; i<n_; i++)
    num += rowLength_[i];

  return num;
}

double SparseMatrix::SumEntries() const
{
  double sum=0;
  for(int i=0; i<n_; i++)
    for(int j=0; j<rowLength_[i]; j++)
      sum += centries_[i][j];

  return sum;
}

void SparseMatrix::SumRowEntries(double * rowSums) const
{
  for(int i=0; i<n_; i++)
  {
    double sum=0;
    for(int j=0; j<rowLength_[i]; j++)
      sum += centries_[i][j];
    rowSums[i] = sum;
  }
}

void SparseMatrix::MakeLinearDataArray(double * data) const
{
  int count=0;
  for(int i=0; i<n_; i++)
  {
    for(int j=0; j<rowLength_[i]; j++)
    {
      data[count] = centries_[i][j];
      count++;
    }
  }   
}

void SparseMatrix::MakeLinearRowIndexArray(double * indices) const
{
  int count=0;
  for(int i=0; i<n_; i++)
  {
    for(int j=0; j<rowLength_[i]; j++)
    {
      indices[count] = i;
      count++;
    }
  }   
}

void SparseMatrix::MakeLinearColumnIndexArray(double * indices) const
{
  int count=0;
  for(int i=0; i<n_; i++)
  {
    for(int j=0; j<rowLength_[i]; j++)
    {
      indices[count] = cindices_[i][j];
      count++;
    }
  }   
}

void SparseMatrix::MakeLinearRowIndexArray(int * indices) const
{
  int count=0;
  for(int i=0; i<n_; i++)
  {
    for(int j=0; j<rowLength_[i]; j++)
    {
      indices[count] = i;
      count++;
    }
  }   
}

void SparseMatrix::MakeLinearColumnIndexArray(int * indices) const
{
  int count=0;
  for(int i=0; i<n_; i++)
  {
    for(int j=0; j<rowLength_[i]; j++)
    {
      indices[count] = cindices_[i][j];
      count++;
    }
  }   
}

void SparseMatrix::FreeTranspositionIndices()
{
  if (transposeIndex_ == NULL)
    return;

  for(int i=0; i<n_; i++)
    free(transposeIndex_[i]);
  free(transposeIndex_);

  transposeIndex_ = NULL;
}

void SparseMatrix::BuildTranspositionIndices()
{
  if (transposeIndex_ != NULL)
    return;

  transposeIndex_ = (int**) malloc (sizeof(int*) * n_);

  int * buffer = (int*) calloc (GetNumColumns(), sizeof(int));

  for(int i=0; i<n_; i++)
  {
    transposeIndex_[i] = (int*) malloc (sizeof(int) * rowLength_[i]);
    for(int j=0; j<rowLength_[i]; j++)
    {
      transposeIndex_[i][j] = buffer[cindices_[i][j]];
      buffer[cindices_[i][j]]++;
    }   
  }  

  free(buffer);
}

double SparseMatrix::SkewSymmetricCheck()
{
  double maxEntry = 0;  

  BuildTranspositionIndices();  

  for(int i=0; i<n_; i++)  
  {    
    for(int j=0; j<GetRowLength(i); j++)    
    {      
      double entry1 = GetEntry(i, j);      
      int tindex = TransposedIndex(i, j);      
      double entry2 = GetEntry(GetColumnIndex(i,j), tindex);      

      // entry1 + entry2 should be zero          
      if (fabs(entry1 + entry2) > maxEntry)
        maxEntry = fabs(entry1 + entry2);
    }  
  }  

  FreeTranspositionIndices();

  return maxEntry;
}

void SparseMatrix::SymmetrizeMatrix()
{
  for(int i=0; i<n_; i++)
  {
    for(int j=0; j<rowLength_[i]; j++)
    {
      int jAbs = cindices_[i][j];

      if (jAbs >= i)
        break; 

      // copy elt (jAbs,i) into position (i,jAbs)
      centries_[i][j] = centries_[jAbs][TransposedIndex(i,j)];
    }
  }
}

double SparseMatrix::GetMaxAbsEntry() const
{
  double max = 0;
  for(int i=0; i<n_; i++)
  {
    for(int j=0; j<rowLength_[i]; j++)
    {
      double el = fabs(GetEntry(i,j));
      if (el > max)
        max = el;
    }
  }
  
  return max;
}

// solve M * x = b
// ASSUMES the sparse matrix is diagonal !!!
void SparseMatrix::DiagonalSolve(double * rhs) const
{
  for(int i=0; i< n_; i++)
    rhs[i] /= centries_[i][0]; // the diagonal element
}

void SparseMatrix::BuildRenumberingVector(int nConstrained, int nSuper, int numFixedDOFs, int * fixedDOFs, int ** superDOFs, int oneIndexed)
{
  // superRows_[i] is the row index in the super matrix corresponsing to row i of constrained matrix
  (*superDOFs) = (int*) malloc (sizeof(int) * nConstrained);
  int constrainedDOF = 0;
  int superDOF = 0;
  for(int i=0; i<numFixedDOFs; i++)
  {
    int nextSuperDOF = fixedDOFs[i];
    nextSuperDOF -= oneIndexed;
    if ( (nextSuperDOF >= nSuper) || (nextSuperDOF < 0) )
    {
      printf("Error: invalid fixed super DOF %d specified.\n", nextSuperDOF);
      exit(1);
    }

    while (superDOF < nextSuperDOF)
    {
      if (constrainedDOF >= nConstrained)
      {
        printf("Error: too many DOFs specified.\n");
        exit(1);
      }
      (*superDOFs)[constrainedDOF] = superDOF; 
      constrainedDOF++;
      superDOF++;
    }

    superDOF++; // skip the deselected DOF
  }
  while (superDOF < nSuper)
  {
    if (constrainedDOF >= nConstrained)
    {
      printf("Error: too many DOFs specified.\n");
      exit(1);
    }
    (*superDOFs)[constrainedDOF] = superDOF; 

    constrainedDOF++;
    superDOF++;
  }
}

void SparseMatrix::BuildSuperMatrixIndices(int numFixedRowColumns, int * fixedRowColumns, SparseMatrix * superMatrix, int oneIndexed)
{
  BuildSuperMatrixIndices(numFixedRowColumns, fixedRowColumns, numFixedRowColumns, fixedRowColumns, superMatrix, oneIndexed); 
}

void SparseMatrix::BuildSuperMatrixIndices(int numFixedRows, int * fixedRows, int numFixedColumns, int * fixedColumns, SparseMatrix * superMatrix, int oneIndexed)
{
  int numColumns = GetNumColumns();
  int numSuperColumns = superMatrix->GetNumColumns();
 
  if ((n_ + numFixedRows != superMatrix->n_) || (numColumns + numFixedColumns != numSuperColumns) )
  {
    printf("Error in BuildSuperMatrixIndices: number of constrained DOFs does not match the size of the two matrices.\n");
    printf("my num rows: %d num fixed rows: %d super matrix num rows: %d\n", n_, numFixedRows, superMatrix->n_);
    printf("my num columns: %d num fixed columns: %d super matrix num columns: %d\n", numColumns, numFixedColumns, numSuperColumns);
    exit(1);
  }

  // build row renumbering function:
  BuildRenumberingVector(n_, superMatrix->n_, numFixedRows, fixedRows, &superRows_, oneIndexed);
  // build column renumbering function:
  int * superColumns_;
  BuildRenumberingVector(numColumns, numSuperColumns, numFixedColumns, fixedColumns, &superColumns_, oneIndexed);

  // superRows_[i] is the row index in the super matrix corresponsing to row i of constrained matrix
  // superColumns_[i] is the absolute column index in the super matrix corresponsing to absolute column i of constrained matrix

  // build column indices
  superMatrixIndices_ = (int**) malloc (sizeof(int*) * n_);
  for(int i=0; i < n_; i++)
  {
    superMatrixIndices_[i] = (int*) malloc (sizeof(int) *  rowLength_[i]);
    for(int j=0; j < rowLength_[i]; j++)
    {
      int iConstrained = i;
      int jConstrainedAbsolute = cindices_[iConstrained][j];
      int iSuper = superRows_[iConstrained];
      int jSuperAbsolute = superColumns_[jConstrainedAbsolute];
      int jSuper = superMatrix->GetInverseIndex(iSuper, jSuperAbsolute);
      if (jSuper < 0)
      {
        printf("Error in BuildSuperMatrixIndices: failed to compute inverse index.\n");
        printf("i=%d j=%d iConstrained=%d jConstrainedAbsolute=%d iSuper=%d jSuperAbsolute=%d jSuper=%d\n", i, j, iConstrained, jConstrainedAbsolute, iSuper, jSuperAbsolute, jSuper);
        fflush(NULL);
        exit(1);
      }
      superMatrixIndices_[i][j] = jSuper;
    }
  } 

  free(superColumns_);
}

void SparseMatrix::AssignSuperMatrix(SparseMatrix * superMatrix)
{
  for(int i=0; i<n_; i++)
  {
    double * row = superMatrix->centries_[superRows_[i]];
    int * indices = superMatrixIndices_[i];
    for(int j=0; j < rowLength_[i]; j++)
      centries_[i][j] = row[indices[j]];
  }
}

void SparseMatrix::BuildSubMatrixIndices(SparseMatrix & mat2)
{
  subMatrixIndices_ = (int**) malloc (sizeof(int*) * n_);
  subMatrixIndexLengths_ = (int*) malloc (sizeof(int*) * n_);

  for(int i=0; i<n_; i++)
  {
    subMatrixIndices_[i] = (int*) malloc (sizeof(int) * mat2.rowLength_[i]);
    subMatrixIndexLengths_[i] = mat2.rowLength_[i];
    int * indices = mat2.cindices_[i];
    for(int j=0; j < mat2.rowLength_[i]; j++)
    {
      // finds the position in row i of element with column index jAbsolute
      // int inverseIndex(int i, int jAbsolute);
      subMatrixIndices_[i][j] = GetInverseIndex(i,indices[j]);
    }
  }
}

SparseMatrix & SparseMatrix::AddSubMatrix(double factor, SparseMatrix & mat2)
{
  for(int i=0; i<n_; i++)
  {
    int * indices = subMatrixIndices_[i];
    for(int j=0; j < mat2.rowLength_[i]; j++)
      centries_[i][indices[j]] += factor * mat2.centries_[i][j];
  }

  return *this;
}

int SparseMatrix::GetNumLowerTriangleEntries() const
{
  int num = 0;
  for(int i=0; i<n_; i++)
  {
    for(int j=0; j < rowLength_[i]; j++)
    {
      if (cindices_[i][j] <= i)
        num++;
    }
  }
  return num;
}

int SparseMatrix::GetNumUpperTriangleEntries() const
{
  int num = 0;
  for(int i=0; i<n_; i++)
  {
    for(int j=0; j < rowLength_[i]; j++)
    {
      if (cindices_[i][j] >= i)
        num++;
    }
  }
  return num;
}

int SparseMatrix::GenerateNAGFormat(double * a, int * irow, int * icol, int * istr) const
{
  int num = 0;
  for(int i=0; i<n_; i++)
  {
    istr[i] = num; // starting address of row i
    for(int j=0; j < rowLength_[i]; j++)
    {
      if (cindices_[i][j] <= i) // over lower triangle
      {
        a[num] = centries_[i][j];
        irow[num] = i+1; // NAG is 1-indexed
        icol[num] = cindices_[i][j]+1; // NAG is 1-indexed
        num++;
      }
    }
  }
  
  istr[n_] = num;

  return num;
} 

void SparseMatrix::GenerateCompressedRowMajorFormat(double * a, int * ia, int * ja, int upperTriangleOnly, int oneIndexed) const
{
  int count = 0;
  for(int row=0; row<n_; row++)
  {
    if (ia != NULL)
      ia[row] = count + oneIndexed;

    int rowLength = GetRowLength(row);
    for(int j=0; j< rowLength; j++)
    {
      if ((!upperTriangleOnly) || (cindices_[row][j] >= row))
      {
        if (a != NULL)
          a[count] = centries_[row][j];         
        if (ja != NULL)
          ja[count] = cindices_[row][j] + oneIndexed; 
        count++;
      }
    }
  }

  if (ia != NULL)
    ia[n_] = count + oneIndexed;
}

void SparseMatrix::GenerateCompressedRowMajorFormat_four_array(double * values, int * columns, int * pointerB, int * pointerE, int upperTriangleOnly, int oneIndexed) const
{
  int count = 0;
  for(int row=0; row<n_; row++)
  {
    if (pointerB != NULL)
      pointerB[row] = count + oneIndexed;

    int rowLength = GetRowLength(row);
    for(int j=0; j< rowLength; j++)
    {
      if ((!upperTriangleOnly) || (cindices_[row][j] >= row))
      {
        if (values != NULL)
          values[count] = centries_[row][j];         
        if (columns != NULL)
          columns[count] = cindices_[row][j] + oneIndexed; 
        count++;
      }
    }

    if (pointerE != NULL)
      pointerE[row] = count + oneIndexed;
  }
}

int SparseMatrix::Save(char * filename, int oneIndexed) const
{
  FILE * fout = fopen(filename,"wa");
  if (!fout)
    return 1;

  fprintf(fout,"%d\n%d\n",n_,GetNumColumns());
  for(int i=0; i< n_; i++)
  {
    for(int j=0; j < rowLength_[i]; j++)
    {
      int index = cindices_[i][j]; 
      double entry = centries_[i][j];
      fprintf(fout,"%d %d %.15G\n",i + oneIndexed, index + oneIndexed, entry);
    }
  }
  fclose(fout);

  return 0;
}

int SparseMatrix::SaveToMatlabFormat(char * filename) const
{
  FILE * fout = fopen(filename,"wa");
  if (!fout)
    return 1;

  for(int i=0; i< n_; i++)
  {
    for(int j=0; j < rowLength_[i]; j++)
    {
      int index = cindices_[i][j]; 
      double entry = centries_[i][j];
      fprintf(fout,"%d %d %.15G\n",i + 1, index + 1, entry);
    }
  }
  fclose(fout);

  return 0;
}

void SparseMatrix::RemoveColumn(int index)
{
  // remove column 'index'
  for(int i=0; i<n_; i++)
  {
    // over all rows
    for(int j=0; j<rowLength_[i]; j++)
    {
      // seek for entry 'index'
      if (cindices_[i][j] == index) // found
      {
        // shift all elements ahead one step back
        for(int k=j; k<rowLength_[i]-1; k++)
        {
          cindices_[i][k] = cindices_[i][k+1];
          centries_[i][k] = centries_[i][k+1];
        } 
        rowLength_[i]--;
      }
    }

    // decrease indices for DOFs above index
    for(int j=0; j<rowLength_[i]; j++)
    {
      if(cindices_[i][j] > index)     
      {
        // decrease index
        cindices_[i][j]--;
      }
    }   
  }
}

void SparseMatrix::RemoveRowColumn(int index)
{
  // remove row 'index'
  free(centries_[index]);
  free(cindices_[index]);

  for(int i=index; i<n_-1; i++)
  {
    centries_[i] = centries_[i+1];
    cindices_[i] = cindices_[i+1];
    rowLength_[i] = rowLength_[i+1];
  }

  // remove column 'index'
  for(int i=0; i<n_-1; i++)
  {
    // over all rows
    for(int j=0; j<rowLength_[i]; j++)
    {
      // seek for entry 'index'
      if (cindices_[i][j] == index) // found
      {
        // shift all elements ahead one step back
        for(int k=j; k<rowLength_[i]-1; k++)
        {
          cindices_[i][k] = cindices_[i][k+1];
          centries_[i][k] = centries_[i][k+1];
        } 
        rowLength_[i]--;
      }
    }

    // decrease indices for DOFs above index
    for(int j=0; j<rowLength_[i]; j++)
    {
      if(cindices_[i][j] > index)     
      {
        // decrease index
        cindices_[i][j]--;
      }
    }   
  }

  n_--;
}

void SparseMatrix::RemoveRowsColumnsSlow(int numRemovedRowsColumns, int * removedRowsColumns, int oneIndexed)
{
  for(int i=0; i<numRemovedRowsColumns; i++)
    RemoveRowColumn(removedRowsColumns[i]-i-oneIndexed);
}

void SparseMatrix::RemoveColumns(int numColumns, int * columns, int oneIndexed)
{
  for(int i=0; i<numColumns; i++)
    RemoveColumn(columns[i]-i-oneIndexed);
}

void SparseMatrix::RemoveRowsColumns(int numRemovedRowsColumns, int * removedRowsColumns, int oneIndexed)
{
  // the removed dofs must be pre-sorted
  // build a map from old dofs to new ones
  vector<int> oldToNew(n_);
  int dof = 0;
  int dofCount = 0;
  for(int i=0; i<numRemovedRowsColumns; i++)
  {
    while (dof < removedRowsColumns[i] - oneIndexed)
    {
      oldToNew[dof] = dofCount;
      dofCount++;
      dof++;
    }
    oldToNew[dof] = -1;
    dof++;
  }
  while (dof < n_)
  {
    oldToNew[dof] = dofCount;
    dofCount++;
    dof++;
  }

  // now, traverse all rows and renumber the entries
  int targetRow = 0;
  for(int sourceRow = 0; sourceRow < n_; sourceRow++)
  {
    if (oldToNew[sourceRow] == -1)
    {
      free(cindices_[sourceRow]);    
      free(centries_[sourceRow]);    
      continue;
    }

    int targetIndex = 0;
    for(int sourceIndex=0; sourceIndex<rowLength_[sourceRow]; sourceIndex++)
    {
      int oldIndex = cindices_[sourceRow][sourceIndex];
      int newIndex = oldToNew[oldIndex];
      if (newIndex == -1)
        continue;
      cindices_[sourceRow][targetIndex] = newIndex;
      centries_[sourceRow][targetIndex] = centries_[sourceRow][sourceIndex];
      targetIndex++;
    }

    cindices_[sourceRow] = (int*) realloc(cindices_[sourceRow], sizeof(int) * targetIndex);
    centries_[sourceRow] = (double*) realloc(centries_[sourceRow], sizeof(double) * targetIndex);

    cindices_[targetRow] = cindices_[sourceRow];
    centries_[targetRow] = centries_[sourceRow];
    rowLength_[targetRow] = targetIndex;
    targetRow++;
  }

  n_ -= numRemovedRowsColumns;
  centries_ = (double**) realloc(centries_, sizeof(double*) * n_);
  cindices_ = (int**) realloc(cindices_, sizeof(double*) * n_);
  rowLength_ = (int*) realloc(rowLength_, sizeof(int) * n_);
}

double SparseMatrix::GetInfinityNorm() const
{
  double norm = 0.0;

  for(int i=0; i< n_; i++)
  {
    double absRowSum = 0;

    for(int j=0; j<rowLength_[i]; j++)
    {
      absRowSum += fabs(centries_[i][j]);
    }

    if (absRowSum > norm)
      norm = absRowSum;
  }

  return norm;
}

void SparseMatrix::DoOneGaussSeidelIteration(double * x, const double * b) const
{
  for(int i=0; i<n_; i++)
  {
    double buffer = b[i];
    int diagIndex = -1;
    for(int j=0; j<rowLength_[i]; j++)
    {
      int column = cindices_[i][j];
      if (column != i)
        buffer -= centries_[i][j] * x[column];
      else
        diagIndex = j;
    }
    x[i] = buffer / centries_[i][diagIndex];
  }
}

void SparseMatrix::ComputeResidual(const double * x, const double * b, double * residual) const
{
  MultiplyVector(x,residual);
  for(int i=0; i<n_; i++)
    residual[i] -= b[i];
}

double SparseMatrix::CheckLinearSystemSolution(const double * x, const double * b, int verbose, double * buffer) const
{
  double * bufferv = NULL;

  if (buffer == NULL)
  {
    bufferv = (double*) malloc (sizeof(double) * n_);
    buffer = bufferv;
  }

  MultiplyVector(x,buffer);

  double inftyNorm = 0;
  double inftyNorm_b = 0;
  for(int i=0; i<n_; i++)
  {
    if (fabs(buffer[i] - b[i]) > inftyNorm)
      inftyNorm = fabs(buffer[i] - b[i]);

    if (fabs(b[i]) > inftyNorm_b)
      inftyNorm_b = fabs(b[i]);
  }

  if (verbose)
  {
    printf("Infinity residual norm ||Ax-b|| is %G. ||b|| is %G.\n", inftyNorm, inftyNorm_b);
    printf("Relative infinity residual norm ||Ax-b||/||b|| is %G.\n", inftyNorm / inftyNorm_b);
  }

  free(bufferv);
  
  return inftyNorm / inftyNorm_b;
}

void SparseMatrix::MakeDenseMatrix(double * denseMatrix) const
{
  memset(denseMatrix, 0, GetNumRows() * GetNumColumns());
  for(int i=0; i< n_; i++)
    for(int j=0; j<rowLength_[i]; j++)
      denseMatrix[n_ * cindices_[i][j] + i] = centries_[i][j];
}

void SparseMatrix::MultiplyRow(int row, double scalar) // multiplies all elements in row 'row' with scalar 'scalar'
{
  for(int j=0; j<rowLength_[row]; j++)
    centries_[row][j] *= scalar;
      
}

int SparseMatrix::GetNumColumns() const 
{
  int numColumns = -1;
  for(int i=0; i< n_; i++)
  {
    for(int j=0; j<rowLength_[i]; j++)
    {
      if (cindices_[i][j] > numColumns)
        numColumns = cindices_[i][j];
    }
  }
  return numColumns + 1;
}

void SparseMatrix::IncreaseNumRows(int newNumRows)
{
  int newn_ = n_ + newNumRows;

  rowLength_ = (int*) realloc (rowLength_, sizeof(int) * newn_);
  for(int i=0; i<newNumRows; i++)
    rowLength_[n_ + i] = 0;

  cindices_ = (int**) realloc (cindices_, sizeof(int*) * newn_);
  for(int i=0; i<newNumRows; i++)
    cindices_[n_ + i] = NULL;

  centries_ = (double**) realloc (centries_, sizeof(double*) * newn_);
  for(int i=0; i<newNumRows; i++)
    centries_[n_ + i] = NULL;

  n_ = newn_;
}

SparseMatrix SparseMatrix::ConjugateMatrix(SparseMatrix & S)
{
  SparseMatrixOutline outline(S.GetNumColumns());

  for(int i=0; i< n_; i++)
  {
    for(int j=0; j<rowLength_[i]; j++)
    {
      int I = i;
      int J = cindices_[i][j];
      double scalar = centries_[i][j];

      // compute tensor product of rows I and J of S
      for(int k=0; k<S.rowLength_[I]; k++)
        for(int l=0; l<S.rowLength_[J]; l++)
        {
          int K = S.cindices_[I][k];
          int L = S.cindices_[J][l];
          outline.AddEntry(K, L, scalar * S.centries_[I][k] * S.centries_[J][l]);
        }
    }
  }
 
  return SparseMatrix(&outline);;
}

SparseMatrix SparseMatrix::Transpose()
{
  SparseMatrixOutline outline(GetNumColumns());

  for(int i=0; i< n_; i++)
    for(int j=0; j<rowLength_[i]; j++)
      outline.AddEntry(cindices_[i][j], i, centries_[i][j]);
 
  return SparseMatrix(&outline);;
}

void SparseMatrix::SetRows(SparseMatrix * source, int startRow, int startColumn) 
{
  for(int i=0; i<source->GetNumRows(); i++)
  {
    int row = startRow + i;
    if (row >= n_)
      return;

    rowLength_[row] = source->GetRowLength(i);
    cindices_[row] = (int*) realloc (cindices_[row], sizeof(int) * rowLength_[row]);
    centries_[row] = (double*) realloc (centries_[row], sizeof(double) * rowLength_[row]);
    for(int j=0; j<rowLength_[row]; j++)
    {
      cindices_[row][j] = startColumn + source->cindices_[i][j];
      centries_[row][j] = source->centries_[i][j];
    }
  }
}

