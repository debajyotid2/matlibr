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

#include "matrix_types.h"
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

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

// Functions for int32 (integer) matrices (for backwards compatibility)
MatrixStatusCode intmat_create(IntMatrix *matrix, int nrow, int ncol);
MatrixStatusCode intmat_assign(IntMatrix *matrix, int nrow, int ncol, int value);
MatrixStatusCode intmat_copy(IntMatrix *copy, const IntMatrix *mat);
MatrixStatusCode intmat_range(IntMatrix *mat, int low, int high, int step, unsigned int dimension);
MatrixStatusCode intmat_print(const IntMatrix *matrix);
MatrixStatusCode intmat_fill(IntMatrix *matrix, const int value);
MatrixStatusCode intmat_fill_random(IntMatrix *matrix, int low, int high, bool replace);
MatrixStatusCode intmat_set_data(IntMatrix *matrix, const int *src, size_t count);
MatrixStatusCode intmat_scale(IntMatrix *mat, const int fac);
MatrixStatusCode intmat_add_scalar(IntMatrix *mat, const int scalar);
MatrixStatusCode intmat_add(IntMatrix *mat_a, const IntMatrix *mat_b);
MatrixStatusCode intmat_sub(IntMatrix *mat_a, const IntMatrix *mat_b);
MatrixStatusCode intmat_mul(const IntMatrix *mat_a, bool transpose_a, const IntMatrix *mat_b,
                            bool transpose_b, IntMatrix *result);
MatrixStatusCode intmat_repeat(IntMatrix *repeated, const IntMatrix *vec, unsigned int dimension,
                               unsigned int repeats);
MatrixStatusCode intmat_vec_add(IntMatrix *mat, const IntMatrix *vec);
MatrixStatusCode intmat_vec_sub(IntMatrix *mat, const IntMatrix *vec);
MatrixStatusCode intmat_gather(const IntMatrix *from, IntMatrix *to, const IntMatrix *indices,
                               unsigned int dimension);
MatrixStatusCode intmat_destroy(IntMatrix *matrix);

// Functions for float64 (double) matrices (for backwards compatibility)
MatrixStatusCode mat_create(Matrix *matrix, int nrow, int ncol);
MatrixStatusCode mat_assign(Matrix *matrix, int nrow, int ncol, double value);
MatrixStatusCode mat_copy(Matrix *copy, const Matrix *mat);
MatrixStatusCode mat_range(Matrix *mat, double low, double high, double step,
                           unsigned int dimension);
MatrixStatusCode mat_print(const Matrix *matrix);
MatrixStatusCode mat_fill(Matrix *matrix, const double value);
MatrixStatusCode mat_fill_random(Matrix *matrix);
MatrixStatusCode mat_set_data(Matrix *matrix, const double *src, size_t count);
MatrixStatusCode mat_scale(Matrix *mat, const double fac);
MatrixStatusCode mat_add_scalar(Matrix *mat, const double scalar);
MatrixStatusCode mat_add(Matrix *mat_a, const Matrix *mat_b);
MatrixStatusCode mat_sub(Matrix *mat_a, const Matrix *mat_b);
MatrixStatusCode mat_mul(const Matrix *mat_a, bool transpose_a, const Matrix *mat_b,
                         bool transpose_b, Matrix *result);
MatrixStatusCode mat_repeat(Matrix *repeated, const Matrix *vec, unsigned int dimension,
                            unsigned int repeats);
MatrixStatusCode mat_vec_add(Matrix *mat, const Matrix *vec);
MatrixStatusCode mat_vec_sub(Matrix *mat, const Matrix *vec);
MatrixStatusCode mat_gather(const Matrix *from, Matrix *to, const IntMatrix *indices,
                            unsigned int dimension);
MatrixStatusCode mat_destroy(Matrix *matrix);

MatrixStatusCode mat_fill_random_gaussian(Matrix *matrix, Matrix *means, Matrix *stds);
MatrixStatusCode mat_abs_sum(double *sum, Matrix *mat);
MatrixStatusCode mat_norm(double *norm, Matrix *mat);

#ifdef __cplusplus
}
#endif

#endif // MATRIX_H
