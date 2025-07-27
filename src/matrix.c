/* Library for reading, writing and manipulating matrix data.

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

#include <cblas.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "matrix.h"

// Function for displaying and handling errors
static void error_handler(const char* function, int line, 
                          MatrixStatusCode status_code, const char* msg) {
    fprintf(stderr, "MATRIX LIB ERROR: %s,%d: - error code %d: %s\n", 
            function, line, status_code, msg);
}

// clang-format off
// Macro for error handling
#define HANDLE_ERROR(status_code, msg)                  \
error_handler(__func__, __LINE__, status_code, msg)     

// Create a matrix
#define DEFINE_MATRIX_CREATE(FUNC_NAME, MATRIX_TYPE, DATA_TYPE)                                 \
MatrixStatusCode FUNC_NAME(MATRIX_TYPE* matrix, int nrows, int ncols) {                         \
    if (matrix==NULL) {                                                                         \
        HANDLE_ERROR(MATRIX_ERR_NULL_PTR, "Null pointer received");                             \
        return MATRIX_ERR_NULL_PTR;                                                             \
    }                                                                                           \
    if (nrows <= 0 || ncols <= 0) {                                                             \
        HANDLE_ERROR(MATRIX_ERR_INVALID_DIMENSION, "Number of rows/columns cannot be <= 0.");   \
        return MATRIX_ERR_INVALID_DIMENSION;                                                    \
    }                                                                                           \
    matrix->nrows = (unsigned int)nrows;                                                        \
    matrix->ncols = (unsigned int)ncols;                                                        \
    matrix->data = (DATA_TYPE *)calloc(matrix->nrows * matrix->ncols, sizeof(DATA_TYPE));       \
    return MATRIX_SUCCESS;                                                                      \
}

// Copy a matrix
#define DEFINE_MATRIX_COPY(FUNC_NAME, MATRIX_TYPE, COPY_FUNC, CREATE_FUNC)                      \
MatrixStatusCode FUNC_NAME(MATRIX_TYPE *copy, const MATRIX_TYPE *mat) {                         \
    if (mat==NULL || copy==NULL) {                                                              \
        HANDLE_ERROR(MATRIX_ERR_NULL_PTR, "Null pointer received.");                            \
        return MATRIX_ERR_NULL_PTR;                                                             \
    }                                                                                           \
    MatrixStatusCode status = CREATE_FUNC(copy, mat->nrows, mat->ncols);                        \
    if (status != MATRIX_SUCCESS) {                                                             \
        return status;                                                                          \
    }                                                                                           \
    COPY_FUNC(mat->nrows * mat->ncols, mat->data, 1, copy->data, 1);                            \
    return MATRIX_SUCCESS;                                                                      \
}

// Fill a matrix with a single value
#define DEFINE_MATRIX_FILL(FUNC_NAME, MATRIX_TYPE, DATA_TYPE)                                   \
MatrixStatusCode FUNC_NAME(MATRIX_TYPE* mat, const DATA_TYPE value) {                           \
    if (mat==NULL || mat->data==NULL) {                                                         \
        HANDLE_ERROR(MATRIX_ERR_NULL_PTR, "Null pointer received.");                            \
        return MATRIX_ERR_NULL_PTR;                                                             \
    }                                                                                           \
    size_t total_elements = mat->nrows * mat->ncols;                                            \
    for (size_t i=0; i<total_elements; ++i) {                                                   \
        mat->data[i] = value;                                                                   \
    }                                                                                           \
    return MATRIX_SUCCESS;                                                                      \
}

// Scale a matrix by a scalar
#define DEFINE_MATRIX_SCALE(FUNC_NAME, MATRIX_TYPE, DATA_TYPE, SCALING_FUNC)                    \
MatrixStatusCode FUNC_NAME(MATRIX_TYPE *mat, const DATA_TYPE fac) {                             \
    if (mat==NULL || mat->data==NULL) {                                                         \
        HANDLE_ERROR(MATRIX_ERR_NULL_PTR, "Null pointer received.");                            \
        return MATRIX_ERR_NULL_PTR;                                                             \
    }                                                                                           \
    SCALING_FUNC(mat->nrows * mat->ncols, fac, mat->data, 1);                                   \
    return MATRIX_SUCCESS;                                                                      \
}

// Print a matrix
#define DEFINE_MATRIX_PRINT(FUNC_NAME, MATRIX_TYPE, FORMAT_SPECIFIER)                           \
MatrixStatusCode FUNC_NAME(const MATRIX_TYPE* mat) {                                            \
    if (mat==NULL || mat->data==NULL) {                                                         \
        HANDLE_ERROR(MATRIX_ERR_NULL_PTR, "Null pointer received.");                            \
        return MATRIX_ERR_NULL_PTR;                                                             \
    }                                                                                           \
    size_t total_elements = mat->nrows*mat->ncols;                                              \
    for (size_t i = 0; i < total_elements; i++) {                                               \
        printf(FORMAT_SPECIFIER " ", mat->data[i]);                                             \
        if (i%mat->ncols == 0) {                                                                \
            printf("\n");                                                                       \
        }                                                                                       \
    }                                                                                           \
    return MATRIX_SUCCESS;                                                                      \
}

// Create a row (1)/column (0) vector (according to specified
// dimension) with elements from "low" to "high"
// (excluded) in "step" steps.
#define DEFINE_MATRIX_RANGE(FUNC_NAME, MATRIX_TYPE, DATA_TYPE, CREATE_FUNC)                         \
MatrixStatusCode FUNC_NAME(MATRIX_TYPE* out, DATA_TYPE low, DATA_TYPE high,                         \
                           DATA_TYPE step, unsigned int dimension) {                                \
    if (out==NULL) {                                                                                \
        HANDLE_ERROR(MATRIX_ERR_NULL_PTR, "Null pointer received.");                                \
        return MATRIX_ERR_NULL_PTR;                                                                 \
    }                                                                                               \
    unsigned int n_elem;                                                                            \
    MatrixStatusCode status;                                                                        \
    if (high <= low) {                                                                              \
        HANDLE_ERROR(MATRIX_ERR_RANGE_INVALID, "Upper limit cannot be <= lower limit.");            \
        return MATRIX_ERR_RANGE_INVALID;                                                            \
    }                                                                                               \
    if (step <= 0) {                                                                                \
        HANDLE_ERROR(MATRIX_ERR_RANGE_INVALID_STEP, "Step cannot be zero or negative.");            \
        return MATRIX_ERR_RANGE_INVALID_STEP;                                                       \
    }                                                                                               \
    n_elem = (int)((high - low) / step);                                                            \
    if (n_elem==0) {                                                                                \
        return MATRIX_SUCCESS;     /* Return immediately if no elements can be created. */          \
    }                                                                                               \
    switch (dimension) {                                                                            \
    case 0:                                                                                         \
        status = CREATE_FUNC(out, n_elem, 1);                                                       \
        if (status != MATRIX_SUCCESS) {                                                             \
            return status;                                                                          \
        }                                                                                           \
        break;                                                                                      \
    case 1:                                                                                         \
        status = CREATE_FUNC(out, 1, n_elem);                                                       \
        if (status != MATRIX_SUCCESS) {                                                             \
            return status;                                                                          \
        }                                                                                           \
        break;                                                                                      \
    default:                                                                                        \
        HANDLE_ERROR(MATRIX_ERR_INVALID_DIMENSION, "Dimension must be either 0 or 1.");             \
        return MATRIX_ERR_INVALID_DIMENSION;                                                        \
    }                                                                                               \
    for (int i = 0; i < n_elem; ++i) {                                                              \
        out->data[i] = low + i * step;                                                              \
    }                                                                                               \
    return MATRIX_SUCCESS;                                                                          \
}

// Repeat a vector along a given dimension
#define DEFINE_MATRIX_REPEAT(FUNC_NAME, MATRIX_TYPE, DATA_TYPE, CREATE_FUNC, COPY_FUNC)                 \
MatrixStatusCode FUNC_NAME(MATRIX_TYPE* repeated, const MATRIX_TYPE *vec,                               \
                            unsigned int dimension, unsigned int repeats) {                             \
    if (vec==NULL || repeated==NULL) {                                                                  \
        HANDLE_ERROR(MATRIX_ERR_NULL_PTR, "Null pointer received.");                                    \
        return MATRIX_ERR_NULL_PTR;                                                                     \
    }                                                                                                   \
    unsigned int nidx, idx_fac, inc;                                                                    \
    MatrixStatusCode status;                                                                            \
    /* Check if vec is a vector */                                                                      \
    if (vec->nrows > 1 && vec->ncols > 1) {                                                             \
        HANDLE_ERROR(MATRIX_ERR_INVALID_DIMENSION, "Only a one dimensional vector can be repeated.");   \
        return MATRIX_ERR_INVALID_DIMENSION;                                                            \
    }                                                                                                   \
    switch (dimension) {                                                                                \
    case 0:                                                                                             \
        if (vec->nrows != 1) {                                                                          \
            HANDLE_ERROR(MATRIX_ERR_INVALID_DIMENSION, "Can only repeat a row vector along rows.");     \
            return MATRIX_ERR_INVALID_DIMENSION;                                                        \
        }                                                                                               \
        nidx = vec->ncols;                                                                              \
        idx_fac = vec->ncols;                                                                           \
        inc = 1;                                                                                        \
        status = CREATE_FUNC(repeated, repeats, vec->ncols);                                            \
        if (status != MATRIX_SUCCESS) {                                                                 \
            return status;                                                                              \
        }                                                                                               \
        break;                                                                                          \
    case 1:                                                                                             \
        if (vec->ncols != 1) {                                                                          \
            HANDLE_ERROR(MATRIX_ERR_INVALID_DIMENSION, "Can only repeat a column vector along rows.");  \
            return MATRIX_ERR_INVALID_DIMENSION;                                                        \
        }                                                                                               \
        nidx = vec->nrows;                                                                              \
        idx_fac = 1;                                                                                    \
        inc = repeats;                                                                                  \
        status = CREATE_FUNC(repeated, vec->nrows, repeats);                                            \
        if (status != MATRIX_SUCCESS) {                                                                 \
            return status;                                                                              \
        }                                                                                               \
        break;                                                                                          \
    default:                                                                                            \
        HANDLE_ERROR(MATRIX_ERR_INVALID_DIMENSION, "Dimension must be either 0 or 1.\n");               \
        return MATRIX_ERR_INVALID_DIMENSION;                                                            \
    }                                                                                                   \
    for (size_t i = 0; i < repeats; i++) {                                                              \
        COPY_FUNC(nidx, vec->data, 1, &(repeated->data[idx_fac * i]), inc);                             \
    }                                                                                                   \
    return MATRIX_SUCCESS;                                                                              \
}

// Add a scalar to a matrix
#define DEFINE_MATRIX_ADD_SCALAR(FUNC_NAME, MATRIX_TYPE, DATA_TYPE)                                     \
MatrixStatusCode FUNC_NAME(MATRIX_TYPE *mat, const DATA_TYPE scalar) {                                  \
    if (mat==NULL || mat->data==NULL) {                                                                 \
        HANDLE_ERROR(MATRIX_ERR_NULL_PTR, "Null pointer received.");                                    \
        return MATRIX_ERR_NULL_PTR;                                                                     \
    }                                                                                                   \
    if (scalar == 0) {                                                                                  \
        return MATRIX_SUCCESS;                                                                          \
    }                                                                                                   \
    int nelem = mat->nrows * mat->ncols;                                                                \
    for (int i=0; i<nelem; ++i) {                                                                       \
        mat->data[i] += scalar;                                                                         \
    }                                                                                                   \
    return MATRIX_SUCCESS;                                                                              \
}

// Add two matrices
// Addition is performed as A := A+B
#define DEFINE_MATRIX_ADD(FUNC_NAME, MATRIX_TYPE, AXPY_FUNC)                                            \
MatrixStatusCode FUNC_NAME(MATRIX_TYPE *mat_a, const MATRIX_TYPE *mat_b) {                              \
    if (mat_a==NULL || mat_b==NULL) {                                                                   \
        HANDLE_ERROR(MATRIX_ERR_NULL_PTR, "Null pointer received.");                                    \
        return MATRIX_ERR_NULL_PTR;                                                                     \
    }                                                                                                   \
    /* Ensure that both matrices are of same shape */                                                   \
    if (mat_a->nrows != mat_b->nrows ||                                                                 \
        mat_a->ncols != mat_b->ncols) {                                                                 \
        HANDLE_ERROR(MATRIX_ERR_DIMENSION_MISMATCH, "Matrices A and B must be of same dimension.");     \
        return MATRIX_ERR_DIMENSION_MISMATCH;                                                           \
    }                                                                                                   \
    AXPY_FUNC(mat_b->nrows * mat_b->ncols, 1, mat_b->data, 1, mat_a->data, 1);                          \
    return MATRIX_SUCCESS;                                                                              \
}

// Subtract two matrices
// Subtraction is performed as A := A - B
#define DEFINE_MATRIX_SUB(FUNC_NAME, MATRIX_TYPE, AXPY_FUNC)                                            \
MatrixStatusCode FUNC_NAME(MATRIX_TYPE *mat_a, const MATRIX_TYPE *mat_b) {                              \
    if (mat_a==NULL || mat_b==NULL) {                                                                   \
        HANDLE_ERROR(MATRIX_ERR_NULL_PTR, "Null pointer received.");                                    \
        return MATRIX_ERR_NULL_PTR;                                                                     \
    }                                                                                                   \
    /* Ensure that both matrices are of same shape */                                                   \
    if (mat_a->nrows != mat_b->nrows ||                                                                 \
        mat_a->ncols != mat_b->ncols) {                                                                 \
        HANDLE_ERROR(MATRIX_ERR_DIMENSION_MISMATCH, "Matrices A and B must be of same dimension.");     \
        return MATRIX_ERR_DIMENSION_MISMATCH;                                                           \
    }                                                                                                   \
    AXPY_FUNC(mat_b->nrows * mat_b->ncols, -1, mat_b->data, 1, mat_a->data, 1);                         \
    return MATRIX_SUCCESS;                                                                              \
}

// Add a vector to a matrix
// Addition is done as: A := A + B
// where vector B is repeated along the number
// of dimensions as required to match A's dimensions.
#define DEFINE_MATRIX_VEC_ADD(FUNC_NAME, MATRIX_TYPE, AXPY_FUNC)                                        \
MatrixStatusCode FUNC_NAME(MATRIX_TYPE *mat, const MATRIX_TYPE *vec) {                                  \
    if (mat==NULL || vec==NULL) {                                                                       \
        HANDLE_ERROR(MATRIX_ERR_NULL_PTR, "Got null pointer for matrix or vector.");                    \
        return MATRIX_ERR_NULL_PTR;                                                                     \
    }                                                                                                   \
    if (vec->nrows > 1 && vec->ncols > 1) {                                                             \
        HANDLE_ERROR(MATRIX_ERR_INVALID_DIMENSION, "Second argument must be a row or column vector.");  \
        return MATRIX_ERR_INVALID_DIMENSION;                                                            \
    }                                                                                                   \
    /* Row vector addition */                                                                           \
    if (vec->nrows == 1) {                                                                              \
        if (vec->ncols != mat->ncols) {                                                                 \
            HANDLE_ERROR(MATRIX_ERR_DIMENSION_MISMATCH,                                                 \
                         "Column dimension mismatch between matrix and row vector.");                   \
            return MATRIX_ERR_DIMENSION_MISMATCH;                                                       \
        }                                                                                               \
        for (size_t i=0; i < mat->nrows; ++i) {                                                         \
            AXPY_FUNC(mat->ncols, 1, vec->data, 1, &mat->data[i*mat->ncols], 1);                        \
        }                                                                                               \
    } else if (vec->ncols == 1) {                                                                       \
        if (vec->nrows != mat->nrows) {                                                                 \
            HANDLE_ERROR(MATRIX_ERR_DIMENSION_MISMATCH,                                                 \
                         "Row dimension mismatch between matrix and column vector.");                   \
            return MATRIX_ERR_DIMENSION_MISMATCH;                                                       \
        }                                                                                               \
        for (size_t i=0; i<mat->nrows; ++i) {                                                           \
            for (size_t j=0; j < mat->ncols; ++j) {                                                     \
                mat->data[i*mat->ncols + j] += vec->data[i];                                            \
            }                                                                                           \
        }                                                                                               \
    }                                                                                                   \
    return MATRIX_SUCCESS;                                                                              \
}

// Subtract a vector from a matrix
// Subtraction is done as: A := A - B
// where vector B is repeated along the number
// of dimensions as required to match A's dimensions.
#define DEFINE_MATRIX_VEC_SUB(FUNC_NAME, MATRIX_TYPE, AXPY_FUNC)                                        \
MatrixStatusCode FUNC_NAME(MATRIX_TYPE *mat, const MATRIX_TYPE *vec) {                                  \
    if (mat==NULL || vec==NULL) {                                                                       \
        HANDLE_ERROR(MATRIX_ERR_NULL_PTR, "Got null pointer for matrix or vector.");                    \
        return MATRIX_ERR_NULL_PTR;                                                                     \
    }                                                                                                   \
    if (vec->nrows > 1 && vec->ncols > 1) {                                                             \
        HANDLE_ERROR(MATRIX_ERR_INVALID_DIMENSION, "Second argument must be a row or column vector.");  \
        return MATRIX_ERR_INVALID_DIMENSION;                                                            \
    }                                                                                                   \
    /* Row vector addition */                                                                           \
    if (vec->nrows == 1) {                                                                              \
        if (vec->ncols != mat->ncols) {                                                                 \
            HANDLE_ERROR(MATRIX_ERR_DIMENSION_MISMATCH,                                                 \
                         "Column dimension mismatch between matrix and row vector.");                   \
            return MATRIX_ERR_DIMENSION_MISMATCH;                                                       \
        }                                                                                               \
        for (size_t i=0; i < mat->nrows; ++i) {                                                         \
            AXPY_FUNC(mat->ncols, -1, vec->data, 1, &mat->data[i*mat->ncols], 1);                       \
        }                                                                                               \
    } else if (vec->ncols == 1) {                                                                       \
        if (vec->nrows != mat->nrows) {                                                                 \
            HANDLE_ERROR(MATRIX_ERR_DIMENSION_MISMATCH,                                                 \
                         "Row dimension mismatch between matrix and column vector.");                   \
            return MATRIX_ERR_DIMENSION_MISMATCH;                                                       \
        }                                                                                               \
        for (size_t i=0; i<mat->nrows; ++i) {                                                           \
            for (size_t j=0; j < mat->ncols; ++j) {                                                     \
                mat->data[i*mat->ncols + j] -= vec->data[i];                                            \
            }                                                                                           \
        }                                                                                               \
    }                                                                                                   \
    return MATRIX_SUCCESS;                                                                              \
}

// Multiply two matrices A and B. Matrices are multiplied
// after transforming them. Matrix dimensions must be such that
//     dim(transform(A)) = m x k
//     dim(transform(B)) = k' x n
#define DEFINE_MATRIX_MUL(FUNC_NAME, MATRIX_TYPE, DATA_TYPE, GEMM_FUNC, FILL_FUNC)                      \
MatrixStatusCode FUNC_NAME(const MATRIX_TYPE *mat_a, bool transpose_a,                                  \
                           const MATRIX_TYPE *mat_b, bool transpose_b,                                  \
                           MATRIX_TYPE *result) {                                                       \
    if (mat_a==NULL || mat_b==NULL || result==NULL) {                                                   \
        HANDLE_ERROR(MATRIX_ERR_NULL_PTR, "Got null pointer.");                                         \
        return MATRIX_ERR_NULL_PTR;                                                                     \
    }                                                                                                   \
    unsigned int m, n, k, k_prime;                                                                      \
    unsigned int lda, ldb;                                                                              \
    DATA_TYPE alpha = 1, beta = 0;                                                                      \
    CBLAS_TRANSPOSE trans_a = transpose_a ? CblasTrans : CblasNoTrans;                                  \
    CBLAS_TRANSPOSE trans_b = transpose_b ? CblasTrans : CblasNoTrans;                                  \
    m = transpose_a ? mat_a->ncols : mat_a->nrows;                                                      \
    k = transpose_a ? mat_a->nrows : mat_a->ncols;                                                      \
    k_prime = transpose_b ? mat_b->ncols : mat_b->nrows;                                                \
    n = transpose_b ? mat_b->nrows : mat_b->ncols;                                                      \
    /* Ensure correct dimensions for matrix multiplication */                                           \
    if (k != k_prime) {                                                                                 \
        HANDLE_ERROR(MATRIX_ERR_DIMENSION_MISMATCH,                                                     \
            "IntMatrix dimensions must satisfy k = k' for multiplying m "                               \
            "x k and k' x n matrices.");                                                                \
        return MATRIX_ERR_DIMENSION_MISMATCH;                                                           \
    }                                                                                                   \
    if (result->nrows != m || result->ncols != n) {                                                     \
        HANDLE_ERROR(MATRIX_ERR_DIMENSION_MISMATCH, "Incorrect dimensions of result matrix.");          \
        return MATRIX_ERR_DIMENSION_MISMATCH;                                                           \
    }                                                                                                   \
    FILL_FUNC(result, (DATA_TYPE)0);                                                                    \
    lda = transpose_a ? m : k;                                                                          \
    ldb = transpose_b ? k : n;                                                                          \
    GEMM_FUNC(trans_a, trans_b, m, n, k, alpha, mat_a->data, lda, mat_b->data,                          \
          ldb, beta, result->data, n);                                                                  \
    return MATRIX_SUCCESS;                                                                              \
}

// Gather rows/columns from "from" and store in
// "to" according to specified indices.
#define DEFINE_MATRIX_GATHER(FUNC_NAME, MATRIX_TYPE, INT_MATRIX_TYPE, COPY_FUNC)                        \
MatrixStatusCode FUNC_NAME(const MATRIX_TYPE *from, MATRIX_TYPE *to, const INT_MATRIX_TYPE *indices,    \
                           unsigned int dimension) {                                                    \
    if (from==NULL || to==NULL || indices==NULL) {                                                      \
        HANDLE_ERROR(MATRIX_ERR_NULL_PTR, "Got null pointer for matrices or indices.");                 \
        return MATRIX_ERR_NULL_PTR;                                                                     \
    }                                                                                                   \
    if (indices->ncols != 1) {                                                                          \
        HANDLE_ERROR(MATRIX_ERR_INVALID_DIMENSION, "'indices' must be a row vector.");                  \
        return MATRIX_ERR_INVALID_DIMENSION;                                                            \
    }                                                                                                   \
    switch (dimension) {                                                                                \
    case 0:  /* Gather rows */                                                                          \
        if (to->ncols != from->ncols) {                                                                 \
            HANDLE_ERROR(MATRIX_ERR_DIMENSION_MISMATCH,                                                 \
                         "'to' and 'from' matrices must have same number of columns.");                 \
            return MATRIX_ERR_DIMENSION_MISMATCH;                                                       \
        }                                                                                               \
        if (to->nrows != indices->nrows) {                                                              \
            HANDLE_ERROR(MATRIX_ERR_DIMENSION_MISMATCH,                                                 \
                         "'to' must have the same number of rows as number of indices in 'indices'.");  \
            return MATRIX_ERR_DIMENSION_MISMATCH;                                                       \
        }                                                                                               \
        for (size_t i=0; i<indices->nrows; ++i) {                                                       \
            if (indices->data[i] >= from->nrows) {                                                      \
                HANDLE_ERROR(MATRIX_ERR_INDEX_OUT_OF_BOUNDS, "Row index out of bounds.");               \
                return MATRIX_ERR_INDEX_OUT_OF_BOUNDS;                                                  \
            }                                                                                           \
            COPY_FUNC(from->ncols, &from->data[indices->data[i] * from->ncols], 1,                      \
                      &to->data[i*to->ncols], 1);                                                       \
        }                                                                                               \
        break;                                                                                          \
    case 1: /* Gather columns */                                                                        \
        if (to->nrows != from->nrows) {                                                                 \
            HANDLE_ERROR(MATRIX_ERR_DIMENSION_MISMATCH,                                                 \
                         "'to' and 'from' matrices must have same number of rows.");                    \
            return MATRIX_ERR_DIMENSION_MISMATCH;                                                       \
        }                                                                                               \
        if (to->ncols != indices->nrows) {                                                              \
            HANDLE_ERROR(MATRIX_ERR_DIMENSION_MISMATCH,                                                 \
                     "'to' must have the same number of columns as number of indices in 'indices'.");   \
            return MATRIX_ERR_DIMENSION_MISMATCH;                                                       \
        }                                                                                               \
        for (size_t i=0; i<indices->nrows; ++i) {                                                       \
            if (indices->data[i] >= from->ncols) {                                                      \
                HANDLE_ERROR(MATRIX_ERR_INDEX_OUT_OF_BOUNDS, "Column index out of bounds.");            \
                return MATRIX_ERR_INDEX_OUT_OF_BOUNDS;                                                  \
            }                                                                                           \
            for (size_t j=0; j<from->nrows; ++j) {                                                      \
                to->data[j*to->ncols+i] = from->data[j*from->ncols+indices->data[i]];                   \
            }                                                                                           \
        }                                                                                               \
        break;                                                                                          \
    default:                                                                                            \
        HANDLE_ERROR(MATRIX_ERR_INVALID_DIMENSION, "Dimension must be either rows(0) or columns(1).");  \
        return MATRIX_ERR_INVALID_DIMENSION;                                                            \
    }                                                                                                   \
    return MATRIX_SUCCESS;                                                                              \
}

// Destroy a matrix
#define DEFINE_MATRIX_DESTROY(FUNC_NAME, MATRIX_TYPE)                                               \
MatrixStatusCode FUNC_NAME(MATRIX_TYPE *matrix) {                                                   \
    if (matrix == NULL) {                                                                           \
        HANDLE_ERROR(MATRIX_ERR_NULL_PTR, "Null pointer received.");                                \
        return MATRIX_ERR_NULL_PTR;                                                                 \
    }                                                                                               \
    if ((matrix->data) != NULL) {                                                                   \
        free(matrix->data);                                                                         \
    }                                                                                               \
    matrix->data = NULL;                                                                            \
    return MATRIX_SUCCESS;                                                                          \
}

/************************************************************/
/*******Basic C implementations of BLAS functions************/
/*******for which integer or double implementations**********/
/*******are not available in OpenBLAS.***********************/
/************************************************************/

/************************************************************/
/***********Integer matrix operations************************/
/************************************************************/

// Scales a vector x := alpha * x. (?scal)
static void iscal(const unsigned int num_elem, const int alpha, int *x,
           const unsigned int incx) {
    if (x == NULL || incx == 0) {
        return;
    }
    for (size_t i = 0; i < num_elem; i++) {
        x[i * incx] *= alpha;
    }
}

// Copies a vector y := x. (?copy)
static void icopy(const unsigned int num_elem, const int *x, const unsigned int incx,
           int *y, const unsigned int incy) {
    if (x == NULL || y == NULL) {
        return;
    }

    if (incx == 0 || incy == 0) {
        return;
    }

    if (incx == 1 && incy == 1) {
        memcpy(y, x, num_elem * sizeof(int));
        return;
    }

    for (size_t i = 0; i < num_elem; i++) {
        y[i * incy] = x[i * incx];
    }
}

// Scales a vector x and adds it to another vector y.(?axpy)
// y := alpha * x + y
static void iaxpy(const unsigned int num_elem, const int alpha, const int *x,
           const unsigned int incx, int *y, const unsigned int incy) {
    if (x == NULL || y == NULL) {
        return;
    }

    if (incx == 0 || incy == 0) {
        return;
    }

    for (size_t i = 0; i < num_elem; i++) {
        y[i * incy] += alpha * x[i * incx];
    }
}

// General dense matrix multiplication for integer matrices. (?gemm)
// C := alpha * op(A) op(B) + beta * C
// Conventions for ?gemm in the BLAS standard are used.
// Please refer to the documentation for details.
// (e.g. "https://www.intel.com/content/www/us/en
// /docs/onemkl/developer-reference-c/2023-1/cblas-gemm
// -001.html#GUID-97718E5C-6E0A-44F0-B2B1-A551F0F164B2")
static void igemm(const CBLAS_TRANSPOSE transa, const CBLAS_TRANSPOSE transb,
           const unsigned int m, const unsigned int n, const unsigned int k,
           const int alpha, const int *a, const unsigned int lda, const int *b,
           const unsigned int ldb, const int beta, int *c,
           const unsigned int ldc) {
    int *a_trans = NULL;
    int *b_trans = NULL;
    unsigned int min_lda, min_ldb, min_ldc;
    bool a_is_transposed = (transa == CblasTrans);
    bool b_is_transposed = (transb == CblasTrans);

    // Validation
    if (a == NULL || b == NULL || c == NULL) {
        perror("ERROR! One of the matrices is a null pointer.");
        return;
    }
    if ((transa == CblasConjTrans || transa == CblasConjNoTrans) ||
        (transb == CblasConjTrans || transb == CblasConjNoTrans)) {
        perror("ERROR! Conjugate transpose is not supported for integer GEMM.");
        return;
    }
    if (a_is_transposed) {
        min_lda = m > 1 ? m : 1;
    } else {
        min_lda = k > 1 ? k : 1;
    }
    if (lda < min_lda) {
        perror("ERROR! Invalid LDA (leading dimension of A).");
        return;
    }
    if (b_is_transposed) {
        min_ldb = k > 1 ? k : 1;
    } else {
        min_ldb = n > 1 ? n : 1;
    }
    if (ldb < min_ldb) {
        perror("ERROR! Invalid LDB (leading dimension of B).");
        return;
    }
    min_ldc = n > 1 ? n : 1;
    if (ldc < min_ldc) {
        perror("ERROR! Invalid LDC (leading dimension of C).");
        return;
    }

    // Quick return
    if (m == 0 || n == 0 || ((alpha == 0 || k == 0) && beta == 1)) {
        return;
    }

    // Scale C by beta
    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j < n; ++j) {
            c[i * ldc + j] = (beta == 0) ? 0: c[i * ldc + j] * beta;
        }
    }

    if (alpha == 0 || k == 0) {
        return;
    }

    // Do C += alpha * op(A) * op(B)
    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j < n; ++j) {
            int dot_prod = 0;
            for (size_t l = 0; l < k; ++l) {
                int a_val = a_is_transposed ? a[l * lda + i] : a[i * lda + l];
                int b_val = b_is_transposed ? b[j * ldb + l] : b[l * ldb + j];
                dot_prod += a_val * b_val;
            }
            c[i * ldc + j] += alpha * dot_prod;
        }
    }
}

// Sparse BLAS-like function for gathering elements from a
// sparse storage vector to a dense storage vector according
// to supplied indices.
static void iusga(const unsigned int num_elem, const int *y, const unsigned int incy,
           int *x, const unsigned int *idxs) {
    if (x == NULL || y == NULL || idxs == NULL) {
        perror("ERROR! Null pointer in array argument(s).");
        return;
    }
    for (size_t i = 0; i < num_elem; i++) {
        x[i] = y[idxs[i * incy]];
    }
}

// Sparse BLAS-like function for gathering elements from a
// sparse storage vector to a dense storage vector according
// to supplied indices. Double precision arrays supported.
static void dusga(const unsigned int num_elem, const double *y,
           const unsigned int incy, double *x, const unsigned int *idxs) {
    if (x == NULL || y == NULL || idxs == NULL) {
        perror("ERROR! Null pointer in array argument(s).");
        return;
    }
    for (size_t i = 0; i < num_elem; i++) {
        x[i] = y[idxs[i * incy]];
    }
}

// Wrapper function for cblas_dgemm
void double_gemm_wrapper(const CBLAS_TRANSPOSE transa, const CBLAS_TRANSPOSE transb,
           const unsigned int m, const unsigned int n, const unsigned int k,
           const double alpha, const double *a, const unsigned int lda, const double *b,
           const unsigned int ldb, const double beta, double *c,
           const unsigned int ldc) {
    cblas_dgemm(CblasRowMajor, transa, transb,
                m, n, k,
                alpha, a, lda,
                b, ldb, beta,
                c, ldc);
}

/************************************************************************/
/***************Functions for IntMatrix (integer data)******************/
/************************************************************************/

DEFINE_MATRIX_CREATE(intmat_create, IntMatrix, int)
DEFINE_MATRIX_COPY(intmat_copy, IntMatrix, icopy, intmat_create)
DEFINE_MATRIX_FILL(intmat_fill, IntMatrix, int)
DEFINE_MATRIX_SCALE(intmat_scale, IntMatrix, int, iscal)
DEFINE_MATRIX_REPEAT(intmat_repeat, IntMatrix, int, intmat_create, icopy)
DEFINE_MATRIX_PRINT(intmat_print, IntMatrix, "%d")
DEFINE_MATRIX_RANGE(intmat_range, IntMatrix, int, intmat_create)
DEFINE_MATRIX_ADD_SCALAR(intmat_add_scalar, IntMatrix, int)
DEFINE_MATRIX_ADD(intmat_add, IntMatrix, iaxpy)
DEFINE_MATRIX_SUB(intmat_sub, IntMatrix, iaxpy)
DEFINE_MATRIX_VEC_ADD(intmat_vec_add, IntMatrix, iaxpy)
DEFINE_MATRIX_VEC_SUB(intmat_vec_sub, IntMatrix, iaxpy)
DEFINE_MATRIX_MUL(intmat_mul, IntMatrix, int, igemm, intmat_fill)
DEFINE_MATRIX_GATHER(intmat_gather, IntMatrix, IntMatrix, icopy)
DEFINE_MATRIX_DESTROY(intmat_destroy, IntMatrix)

// Fill a matrix with random integers between low and high (exclusive)
// with or without replacement.
MatrixStatusCode intmat_fill_random(IntMatrix *mat, int low, int high, bool replace,
                                    unsigned int seed) {
    int *temp_ints = NULL;

    if (low >= high) {
        HANDLE_ERROR(MATRIX_ERR_RANGE_INVALID, "low >= high.");
        return MATRIX_ERR_RANGE_INVALID;
    }
    
    srand(seed);

    // Random numbers with replacement
    if (replace) {
        for (size_t i = 0; i < mat->nrows; i++) {
            for (size_t j = 0; j < mat->ncols; j++) {
                mat->data[i * mat->ncols + j] = low + rand() % (high - low);
            }
        }
        return MATRIX_SUCCESS;
    }

    if (mat->nrows * mat->ncols > (size_t)(high - low)) {
        HANDLE_ERROR(MATRIX_ERR_TOO_MANY_INTS_TO_GENERATE, 
                     "Too many numbers to generate without replacement.");
        return MATRIX_ERR_TOO_MANY_INTS_TO_GENERATE;
    }

    // Random numbers without replacement
    // 1. Generate numbers from low to high (exclusive) and store in arr.
    // 2. Shuffle arr (e.g. using Fisher-Yates, 
    //              https://en.wikipedia.org/wiki/Fisher%E2%80%93Yates_shuffle).
    // 3. Select first nrows * ncols numbers.
    temp_ints = (int *)calloc(high - low, sizeof(int));
    for (int i = low; i < high; i++) {
        temp_ints[i - low] = i;
    }
    for (size_t i = high - low - 1; i > 0; i--) {
        int idx = rand() % i;
        int temp = temp_ints[idx];
        temp_ints[idx] = temp_ints[i];
        temp_ints[i] = temp;
    }
    for (size_t i = 0; i < mat->nrows; i++) {
        for (size_t j = 0; j < mat->ncols; j++) {
            mat->data[i * mat->ncols + j] = temp_ints[i * mat->ncols + j];
        }
    }
    free(temp_ints);
    return MATRIX_SUCCESS;
}

/************************************************************************/
/***************Functions for Matrix (double precision data)*************/
/************************************************************************/

DEFINE_MATRIX_CREATE(mat_create, Matrix, double)
DEFINE_MATRIX_COPY(mat_copy, Matrix, cblas_dcopy, mat_create)
DEFINE_MATRIX_FILL(mat_fill, Matrix, double)
DEFINE_MATRIX_SCALE(mat_scale, Matrix, double, cblas_dscal)
DEFINE_MATRIX_REPEAT(mat_repeat, Matrix, double, mat_create, cblas_dcopy)
DEFINE_MATRIX_PRINT(mat_print, Matrix, "%g")
DEFINE_MATRIX_RANGE(mat_range, Matrix, double, mat_create)
DEFINE_MATRIX_ADD_SCALAR(mat_add_scalar, Matrix, double)
DEFINE_MATRIX_ADD(mat_add, Matrix, cblas_daxpy)
DEFINE_MATRIX_SUB(mat_sub, Matrix, cblas_daxpy)
DEFINE_MATRIX_VEC_ADD(mat_vec_add, Matrix, cblas_daxpy)
DEFINE_MATRIX_VEC_SUB(mat_vec_sub, Matrix, cblas_daxpy)
DEFINE_MATRIX_MUL(mat_mul, Matrix, double, double_gemm_wrapper, mat_fill)
DEFINE_MATRIX_GATHER(mat_gather, Matrix, IntMatrix, cblas_dcopy)
DEFINE_MATRIX_DESTROY(mat_destroy, Matrix)


// Fill a matrix with random numbers between 0.0 and 1.0 (half-open)
void mat_fill_random(Matrix *mat, unsigned int seed) {
    srand(seed);

    for (size_t i = 0; i < mat->nrows; i++)
        for (size_t j = 0; j < mat->ncols; j++)
            mat->data[i * mat->ncols + j] = (double)rand() / (double)(RAND_MAX);
}

// Fill a matrix with random numbers from a Gaussian distribution
// with mean "mean" and standard deviation "std", using the Box-Muller
// transform (https://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform)
void mat_fill_random_gaussian(Matrix *mat, Matrix *means, Matrix *stds,
                              unsigned int seed) {
    srand(seed);

    if (mat == NULL || means == NULL || stds == NULL) {
        perror("ERROR: Null values in argument matrices.");
        return;
    }

    if (means->nrows != mat->ncols || stds->nrows != mat->ncols ||
        means->ncols != 1 || stds->ncols != 1) {
        perror("ERROR: Means/Stddevs matrix is of improper dimensions.");
        return;
    }

    const double two_pi = 2.0 * M_PI;
    double u_1, u_2, mag;

    for (size_t i = 0; i < mat->nrows; i++) {
        for (size_t j = 0; j < mat->ncols; j++) {
            double mean = means->data[j];
            double std = stds->data[j];

            if (std == 0.0) {
                perror("ERROR: Standard deviation cannot be zero.");
                return;
            }

            // Use Box-Muller transform to generate the random number
            // from a standard normal distribution
            do {
                u_1 = (double)rand() / (double)(RAND_MAX);
            } while (u_1 == 0.0);

            u_2 = (double)rand() / (double)(RAND_MAX);

            mag = std * sqrt(-2.0 * log(u_1));

            mat->data[i * mat->ncols + j] =
                rand() % 2 ? mag * cos(two_pi * u_2) + mean
                           : mag * sin(two_pi * u_2) + mean;
        }
    }
}

// Sum of absolute values of matrix elements
double mat_abs_sum(Matrix *mat) {
    if (mat == NULL)
        return 0.0;

    double sum = 0.0;
    sum = cblas_dasum(mat->nrows * mat->ncols, mat->data, 1);
    return sum;
}

// Euclidean norm of a matrix across all rows and
// columns.
double mat_norm(Matrix *mat) {
    if (mat == NULL) {
        perror("ERROR: Null pointer in argument matrix.");
        return 0.0;
    }
    return cblas_dnrm2(mat->nrows * mat->ncols, mat->data, 1);
}
