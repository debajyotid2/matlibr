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

#ifndef MATRIX_H
#define MATRIX_H

#ifdef __cplusplus
extern "C" {
#endif

#include "blas_helpers.h"
#include "matrix_types.h"
#include <stdbool.h>

// Functions for int32 matrices
MatrixStatusCode mati32_create(MatrixI32 *matrix, int nrow, int ncol);
MatrixStatusCode mati32_assign(MatrixI32 *matrix, int nrow, int ncol, int32_t value);
MatrixStatusCode mati32_copy(MatrixI32 *copy, const MatrixI32 *mat);
MatrixStatusCode mati32_range(MatrixI32 *mat, int32_t low, int32_t high, int32_t step,
                              unsigned int dimension);
MatrixStatusCode mati32_print(const MatrixI32 *matrix);
MatrixStatusCode mati32_fill(MatrixI32 *matrix, const int32_t value);
MatrixStatusCode mati32_set_data(MatrixI32 *matrix, const int32_t *src, size_t count);
MatrixStatusCode mati32_scale(MatrixI32 *mat, const int32_t fac);
MatrixStatusCode mati32_add_scalar(MatrixI32 *mat, const int32_t scalar);
MatrixStatusCode mati32_add(MatrixI32 *mat_a, const MatrixI32 *mat_b);
MatrixStatusCode mati32_sub(MatrixI32 *mat_a, const MatrixI32 *mat_b);
MatrixStatusCode mati32_mul(const MatrixI32 *mat_a, bool transpose_a, const MatrixI32 *mat_b,
                            bool transpose_b, MatrixI32 *result);
MatrixStatusCode mati32_repeat(MatrixI32 *repeated, const MatrixI32 *vec, unsigned int dimension,
                               unsigned int repeats);
MatrixStatusCode mati32_vec_add(MatrixI32 *mat, const MatrixI32 *vec);
MatrixStatusCode mati32_vec_sub(MatrixI32 *mat, const MatrixI32 *vec);
MatrixStatusCode mati32_gather(const MatrixI32 *from, MatrixI32 *to, const MatrixI32 *indices,
                               unsigned int dimension);
MatrixStatusCode mati32_destroy(MatrixI32 *matrix);

// Functions for int64 matrices
MatrixStatusCode mati64_create(MatrixI64 *matrix, int nrow, int ncol);
MatrixStatusCode mati64_assign(MatrixI64 *matrix, int nrow, int ncol, int64_t value);
MatrixStatusCode mati64_copy(MatrixI64 *copy, const MatrixI64 *mat);
MatrixStatusCode mati64_range(MatrixI64 *mat, int64_t low, int64_t high, int64_t step,
                              unsigned int dimension);
MatrixStatusCode mati64_print(const MatrixI64 *matrix);
MatrixStatusCode mati64_fill(MatrixI64 *matrix, const int64_t value);
MatrixStatusCode mati64_set_data(MatrixI64 *matrix, const int64_t *src, size_t count);
MatrixStatusCode mati64_scale(MatrixI64 *mat, const int64_t fac);
MatrixStatusCode mati64_add_scalar(MatrixI64 *mat, const int64_t scalar);
MatrixStatusCode mati64_add(MatrixI64 *mat_a, const MatrixI64 *mat_b);
MatrixStatusCode mati64_sub(MatrixI64 *mat_a, const MatrixI64 *mat_b);
MatrixStatusCode mati64_mul(const MatrixI64 *mat_a, bool transpose_a, const MatrixI64 *mat_b,
                            bool transpose_b, MatrixI64 *result);
MatrixStatusCode mati64_repeat(MatrixI64 *repeated, const MatrixI64 *vec, unsigned int dimension,
                               unsigned int repeats);
MatrixStatusCode mati64_vec_add(MatrixI64 *mat, const MatrixI64 *vec);
MatrixStatusCode mati64_vec_sub(MatrixI64 *mat, const MatrixI64 *vec);
MatrixStatusCode mati64_gather(const MatrixI64 *from, MatrixI64 *to, const MatrixI32 *indices,
                               unsigned int dimension);
MatrixStatusCode mati64_destroy(MatrixI64 *matrix);

// Functions for float32 (float) matrices
MatrixStatusCode matf32_create(MatrixF32 *matrix, int nrow, int ncol);
MatrixStatusCode matf32_assign(MatrixF32 *matrix, int nrow, int ncol, float value);
MatrixStatusCode matf32_copy(MatrixF32 *copy, const MatrixF32 *mat);
MatrixStatusCode matf32_range(MatrixF32 *mat, float low, float high, float step,
                              unsigned int dimension);
MatrixStatusCode matf32_print(const MatrixF32 *matrix);
MatrixStatusCode matf32_fill(MatrixF32 *matrix, const float value);
MatrixStatusCode matf32_set_data(MatrixF32 *matrix, const float *src, size_t count);
MatrixStatusCode matf32_scale(MatrixF32 *mat, const float fac);
MatrixStatusCode matf32_add_scalar(MatrixF32 *mat, const float scalar);
MatrixStatusCode matf32_add(MatrixF32 *mat_a, const MatrixF32 *mat_b);
MatrixStatusCode matf32_sub(MatrixF32 *mat_a, const MatrixF32 *mat_b);
MatrixStatusCode matf32_mul(const MatrixF32 *mat_a, bool transpose_a, const MatrixF32 *mat_b,
                            bool transpose_b, MatrixF32 *result);
MatrixStatusCode matf32_repeat(MatrixF32 *repeated, const MatrixF32 *vec, unsigned int dimension,
                               unsigned int repeats);
MatrixStatusCode matf32_vec_add(MatrixF32 *mat, const MatrixF32 *vec);
MatrixStatusCode matf32_vec_sub(MatrixF32 *mat, const MatrixF32 *vec);
MatrixStatusCode matf32_gather(const MatrixF32 *from, MatrixF32 *to, const MatrixI32 *indices,
                               unsigned int dimension);
MatrixStatusCode matf32_destroy(MatrixF32 *matrix);

// Functions for float64 (double) matrices
MatrixStatusCode matf64_create(MatrixF64 *matrix, int nrow, int ncol);
MatrixStatusCode matf64_assign(MatrixF64 *matrix, int nrow, int ncol, double value);
MatrixStatusCode matf64_copy(MatrixF64 *copy, const MatrixF64 *mat);
MatrixStatusCode matf64_range(MatrixF64 *mat, double low, double high, double step,
                              unsigned int dimension);
MatrixStatusCode matf64_print(const MatrixF64 *matrix);
MatrixStatusCode matf64_fill(MatrixF64 *matrix, const double value);
MatrixStatusCode matf64_set_data(MatrixF64 *matrix, const double *src, size_t count);
MatrixStatusCode matf64_scale(MatrixF64 *mat, const double fac);
MatrixStatusCode matf64_add_scalar(MatrixF64 *mat, const double scalar);
MatrixStatusCode matf64_add(MatrixF64 *mat_a, const MatrixF64 *mat_b);
MatrixStatusCode matf64_sub(MatrixF64 *mat_a, const MatrixF64 *mat_b);
MatrixStatusCode matf64_mul(const MatrixF64 *mat_a, bool transpose_a, const MatrixF64 *mat_b,
                            bool transpose_b, MatrixF64 *result);
MatrixStatusCode matf64_repeat(MatrixF64 *repeated, const MatrixF64 *vec, unsigned int dimension,
                               unsigned int repeats);
MatrixStatusCode matf64_vec_add(MatrixF64 *mat, const MatrixF64 *vec);
MatrixStatusCode matf64_vec_sub(MatrixF64 *mat, const MatrixF64 *vec);
MatrixStatusCode matf64_gather(const MatrixF64 *from, MatrixF64 *to, const MatrixI32 *indices,
                               unsigned int dimension);
MatrixStatusCode matf64_destroy(MatrixF64 *matrix);

// Functions for `IntMat` matrices (for backwards compatibility)
#define intmat_create(matrix, nrow, ncol) mati32_create((matrix), (nrow), (ncol))
#define intmat_assign(matrix, nrow, ncol, value) mati32_assign((matrix), (nrow), (ncol), (value))
#define intmat_copy(copy, mat) mati32_copy((copy), (mat))
#define intmat_range(mat, low, high, step, dimension)       \
    mati32_range((mat), (low), (high), (step), (dimension))
#define intmat_print(matrix) mati32_print((matrix))
#define intmat_fill(matrix, value) mati32_fill((matrix), (value))
#define intmat_set_data(matrix, src, count) mati32_set_data((matrix), (src), (count))
#define intmat_scale(mat, fac) mati32_scale((mat), (fac))
#define intmat_add_scalar(mat, scalar) mati32_add_scalar((mat), (scalar))
#define intmat_add(mat_a, mat_b) mati32_add((mat_a), (mat_b))
#define intmat_sub(mat_a, mat_b) mati32_sub((mat_a), (mat_b))
#define intmat_mul(mat_a, transpose_a, mat_b, transpose_b, result)       \
    mati32_mul((mat_a), (transpose_a), (mat_b), (transpose_b), (result))
#define intmat_repeat(repeated, vec, dimension, repeats)     \
    mati32_repeat((repeated), (vec), (dimension), (repeats))
#define intmat_vec_add(mat, vec) mati32_vec_add((mat), (vec))
#define intmat_vec_sub(mat, vec) mati32_vec_sub((mat), (vec))
#define intmat_gather(from, to, indices, dimension)     \
    mati32_gather((from), (to), (indices), (dimension))
#define intmat_destroy(matrix) mati32_destroy((matrix))

MatrixStatusCode intmat_fill_random(IntMatrix *mat, int low, int high, bool replace);

// Functions for `Matrix` matrices (for backwards compatibility)
#define mat_create(matrix, nrow, ncol) matf64_create((matrix), (nrow), (ncol))
#define mat_assign(matrix, nrow, ncol, value) matf64_assign((matrix), (nrow), (ncol), (value))
#define mat_copy(copy, mat) matf64_copy((copy), (mat))
#define mat_range(mat, low, high, step, dimension)          \
    matf64_range((mat), (low), (high), (step), (dimension))
#define mat_print(matrix) matf64_print((matrix))
#define mat_fill(matrix, value) matf64_fill((matrix), (value))
#define mat_set_data(matrix, src, count) matf64_set_data((matrix), (src), (count))
#define mat_scale(mat, fac) matf64_scale((mat), (fac))
#define mat_add_scalar(mat, scalar) matf64_add_scalar((mat), (scalar))
#define mat_add(mat_a, mat_b) matf64_add((mat_a), (mat_b))
#define mat_sub(mat_a, mat_b) matf64_sub((mat_a), (mat_b))
#define mat_mul(mat_a, transpose_a, mat_b, transpose_b, result)          \
    matf64_mul((mat_a), (transpose_a), (mat_b), (transpose_b), (result))
#define mat_repeat(repeated, vec, dimension, repeats)        \
    matf64_repeat((repeated), (vec), (dimension), (repeats))
#define mat_vec_add(mat, vec) matf64_vec_add((mat), (vec))
#define mat_vec_sub(mat, vec) matf64_vec_sub((mat), (vec))
#define mat_gather(from, to, indices, dimension) matf64_gather((from), (to), (indices), (dimension))
#define mat_destroy(matrix) matf64_destroy((matrix))

MatrixStatusCode mat_fill_random(Matrix *mat);
MatrixStatusCode mat_fill_random_gaussian(Matrix *matrix, Matrix *means, Matrix *stds);
MatrixStatusCode mat_abs_sum(double *sum, Matrix *mat);
MatrixStatusCode mat_norm(double *norm, Matrix *mat);

#ifdef __cplusplus
}
#endif

#endif // MATRIX_H
