/* Library for reading, writing and manipulating matrices in C.

  Copyright (c) 2024 Debajyoti Debnath

  Licensed under the Apache License, Version 2.0 (the "License");
  you may not use this file except in compliance with the License.
  You may obtain a copy of the License at

      http://www.apache.org/licenses/LICENSE-2.0

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.
*/

#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <stdbool.h>
#include <stdio.h>

// Matrix for double precision data
typedef struct {
    unsigned int nrows, ncols;
    double *data;
} Matrix;

// Matrix for integer data
typedef struct {
    unsigned int nrows, ncols;
    int *data;
} IntMatrix;

// Operation status codes
typedef enum {
    MATRIX_SUCCESS = 0,
    MATRIX_ERR_NULL_PTR,
    MATRIX_ERR_MEMORY_ALLOCATION,
    MATRIX_ERR_DIMENSION_MISMATCH,
    MATRIX_ERR_INVALID_DIMENSION,
    MATRIX_ERR_RANGE_INVALID,
    MATRIX_ERR_RANGE_INVALID_STEP,
    MATRIX_ERR_INDEX_OUT_OF_BOUNDS,
    MATRIX_ERR_TOO_MANY_INTS_TO_GENERATE,
    MATRIX_ERROR_ZERO_STD_DEV
} MatrixStatusCode;

// Functions for integer matrices
MatrixStatusCode intmat_create(IntMatrix *matrix, int nrow, int ncol);
MatrixStatusCode intmat_copy(IntMatrix *copy, const IntMatrix *mat);
MatrixStatusCode intmat_range(IntMatrix *mat, int low, int high, int step,
                              unsigned int dimension);
MatrixStatusCode intmat_print(const IntMatrix *matrix);
MatrixStatusCode intmat_fill(IntMatrix *matrix, const int value);
MatrixStatusCode intmat_fill_random(IntMatrix *matrix, int low, int high,
                                    bool replace);
MatrixStatusCode intmat_scale(IntMatrix *mat, const int fac);
MatrixStatusCode intmat_add_scalar(IntMatrix *mat, const int scalar);
MatrixStatusCode intmat_add(IntMatrix *mat_a, const IntMatrix *mat_b);
MatrixStatusCode intmat_sub(IntMatrix *mat_a, const IntMatrix *mat_b);
MatrixStatusCode intmat_mul(const IntMatrix *mat_a, bool transpose_a,
                            const IntMatrix *mat_b, bool transpose_b,
                            IntMatrix *result);
MatrixStatusCode intmat_repeat(IntMatrix *repeated, const IntMatrix *vec,
                               unsigned int dimension, unsigned int repeats);
MatrixStatusCode intmat_vec_add(IntMatrix *mat, const IntMatrix *vec);
MatrixStatusCode intmat_vec_sub(IntMatrix *mat, const IntMatrix *vec);
MatrixStatusCode intmat_gather(const IntMatrix *from, IntMatrix *to,
                               const IntMatrix *indices,
                               unsigned int dimension);
MatrixStatusCode intmat_destroy(IntMatrix *matrix);

// Functions for double matrices
MatrixStatusCode mat_create(Matrix *matrix, int nrow, int ncol);
MatrixStatusCode mat_copy(Matrix *copy, const Matrix *mat);
MatrixStatusCode mat_print(const Matrix *matrix);
MatrixStatusCode mat_range(Matrix *mat, double low, double high, double step,
                           unsigned int dimension);
MatrixStatusCode mat_fill(Matrix *matrix, const double value);
MatrixStatusCode mat_fill_random(Matrix *matrix);
MatrixStatusCode mat_fill_random_gaussian(Matrix *matrix, Matrix *means,
                                          Matrix *stds);
MatrixStatusCode mat_scale(Matrix *mat, const double fac);
MatrixStatusCode mat_abs_sum(double* sum, Matrix *mat);
MatrixStatusCode mat_norm(double* norm, Matrix *mat);
MatrixStatusCode mat_add_scalar(Matrix *mat, const double scalar);
MatrixStatusCode mat_add(Matrix *mat_a, const Matrix *mat_b);
MatrixStatusCode mat_sub(Matrix *mat_a, const Matrix *mat_b);
MatrixStatusCode mat_mul(const Matrix *mat_a, bool transpose_a,
                         const Matrix *mat_b, bool transpose_b, Matrix *result);
MatrixStatusCode mat_repeat(Matrix *repeated, const Matrix *vec,
                            unsigned int dimension, unsigned int repeats);
MatrixStatusCode mat_vec_add(Matrix *mat, const Matrix *vec);
MatrixStatusCode mat_vec_sub(Matrix *mat, const Matrix *vec);
MatrixStatusCode mat_gather(const Matrix *from, Matrix *to,
                            const IntMatrix *indices, unsigned int dimension);
MatrixStatusCode mat_destroy(Matrix *matrix);

#endif // _MATRIX_H_
