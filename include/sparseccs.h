/* ***********************************************
 *
 * Copyright (c) 2010 Boris Petrov

 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:

 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.

 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *
 * ***********************************************
 */


#ifndef SPARSE_CCS_H
#define SPARSE_CCS_H

#define THRESHOLD 0.0000000001

#include <math.h>
#include <stdio.h>
#include <cstdlib>

typedef double real;

class SparseCCS
{
public:

	SparseCCS(int rows, int columns, real* vals, int* rowInds, int* colPtr);
	~SparseCCS();

	// sets the value at that row and column - this modifies the matrix and could be O(number_of_nonzeros) in the worst case!
	void set_cell_value(int row, int col, real value);

	// Sparse matrix transposition.
	SparseCCS* transposed() const;
	// Sparse matrix inversion - the matrix should be symmetric, lower triangular, positive definite. Cholesky decomposition is used.
	SparseCCS* symmetric_lower_inverse() const;

	// Multiplies the matrix with a column vector.
	real* operator *(const real* right) const;
//	double* operator *(const double* right) const;
	// Returns the element at this row and column.
	real operator ()(int row, int column) const;

	// Equality and inequality.
	bool operator ==(const SparseCCS&) const;
	bool operator !=(const SparseCCS&) const;

	// Multiplies all elements in the matrix by a number IN PLACE!
	void multiply_by_number(real);

	SparseCCS* cholesky_decompose_lower_triangular_returning_lower_triangular() const;
	SparseCCS* ldlt_decompose_lower_triangular_returning_lower_triangular(real*& D) const;
	SparseCCS* invert_diagonal_matrix() const;
	SparseCCS* invert_lower_triangular_cholesky_decomposed() const;
	SparseCCS* invert_lower_triangular_cholesky_decomposed_returning_lower_triangular() const;
	SparseCCS* invert_lower_triangular_ldlt_decomposed(const real* D) const;
	SparseCCS* invert_lower_triangular_ldlt_decomposed_returning_lower_triangular(const real* D) const;
	real* solve_eqn(const real* b) const;
	real* solve_ldlt(const real* D, const real* b) const;
	void solve_ldlt_in_place(const SparseCCS* LT, const real* D, real* b, int up_to = 0) const;

	// Adds two sparse matrices.
	SparseCCS* add(const SparseCCS& second) const;

	SparseCCS* multiply_F(const SparseCCS& second) const;
	SparseCCS* multiply_LM(const SparseCCS& second) const;
	SparseCCS* multiply_diagonal_dense(const real* second, int secondCols) const;
	SparseCCS* multiply_diagonal_sparse(const real* second, int secondCols) const;
	SparseCCS* multiply_returning_unordered_F(const SparseCCS& second) const;
	real* multiply_returning_diagonal(const SparseCCS& second) const;
	SparseCCS* multiply_returning_lower_triangular_F(const SparseCCS& second) const;
	SparseCCS* multiply_returning_lower_triangular_LM(const SparseCCS& second) const;

	SparseCCS* multiply_three_F(const SparseCCS& second, const SparseCCS& third) const;
	SparseCCS* multiply_three_LM(const SparseCCS& second, const SparseCCS& third) const;
	real* multiply_three_returning_diagonal(const SparseCCS& second, const SparseCCS& third) const;
	SparseCCS* multiply_three_returning_lower_triangular_F(const SparseCCS& second, const SparseCCS& third) const;
	SparseCCS* multiply_three_returning_lower_triangular_LM(const SparseCCS& second, const SparseCCS& third) const;

	inline int columns_count() const { return cols; }
	inline int rows_count() const { return rows; }

	inline real* values() const { return vals; }
	inline int* row_indices() const { return rowind; }
	inline int* column_pointers() const { return colptr; }

	void write_matrix_file(const char *) const;
	static SparseCCS* read_matrix_file(const char *);

	static SparseCCS* deep_copy(const SparseCCS* matrix);

	static SparseCCS* generate_random(int rows, int cols, real XMin, real XMax, double fillPercent);
	static SparseCCS* generate_random_lower(int size, real XMin, real XMax, double fillPercent);

private:

	SparseCCS(int rows, int columns, int nnz)
	{
		this->rows = rows;
		cols = columns;

		vals = new real[nnz];
		colptr = new int[cols + 1];
		rowind = new int[nnz];
	}

	SparseCCS(int rows, int columns)
	{
		this->rows = rows;
		cols = columns;

		vals = NULL;
		rowind = NULL;
		colptr = new int[cols + 1];
	}

	inline void setNNZ(int nnz)
	{
		vals = new real[nnz];
		rowind = new int[nnz];
	}

	int cols;
	int rows;

	int* colptr;	// length cols + 1
	int* rowind;	// length colptr[cols]
	real* vals;		// length colptr[cols]
private:
	template <typename T>
	class Vector
	{
	public:
		inline Vector(int initialSize)
		{
			allocated = initialSize;
			count = 0;
			values = new T[allocated];
		}

		inline ~Vector()
		{
			delete[] values;
		}

		void add(const T& elem)
		{
			if (count == allocated)
			{
				T *temp = new T[allocated];
				for (int i = 0; i < allocated; i++)
				{
					temp[i] = values[i];
				}
				delete[] values;
				values = new T[allocated * 2];
				for (int i = 0; i < allocated; i++)
				{
					values[i] = temp[i];
				}
				delete[] temp;
				allocated *= 2;
			}
			values[count++] = elem;
		}

		inline const T& operator[] (int index) const
		{
			return values[index];
		}

		inline T& operator[] (int index)
		{
			return values[index];
		}

		inline int size() const
		{
			return count;
		}

	private:
		T *values;
		int count;
		int allocated;
	};
};

#endif // SPARSE_CCS_H
