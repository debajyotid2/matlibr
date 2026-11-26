#include "blas_helpers.h"

// Scales a vector x := alpha * x. (?scal)
#define DEFINE_SCAL(FUNC_NAME, DATA_TYPE)                                      \
    void FUNC_NAME(const size_t num_elem, const DATA_TYPE alpha, DATA_TYPE *x, \
                   const unsigned int incx) {                                  \
        if (x == NULL || incx == 0) {                                          \
            return;                                                            \
        }                                                                      \
        for (size_t i = 0; i < num_elem; i++) {                                \
            x[i * incx] *= alpha;                                              \
        }                                                                      \
    }

// Copies a vector y := x. (?copy)
#define DEFINE_COPY(FUNC_NAME, DATA_TYPE)                                     \
    void FUNC_NAME(const size_t num_elem, const DATA_TYPE *MATRIX_RESTRICT x, \
                   const unsigned int incx, DATA_TYPE *MATRIX_RESTRICT y,     \
                   const unsigned int incy) {                                 \
        if (x == NULL || y == NULL) {                                         \
            return;                                                           \
        }                                                                     \
                                                                              \
        if (incx == 0 || incy == 0) {                                         \
            return;                                                           \
        }                                                                     \
                                                                              \
        if (incx == 1 && incy == 1) {                                         \
            memcpy(y, x, num_elem * sizeof(DATA_TYPE));                       \
            return;                                                           \
        }                                                                     \
                                                                              \
        for (size_t i = 0; i < num_elem; i++) {                               \
            y[i * incy] = x[i * incx];                                        \
        }                                                                     \
    }

// Scales a vector x and adds it to another vector y.(?axpy)
// y := alpha * x + y
#define DEFINE_AXPY(FUNC_NAME, DATA_TYPE)                                       \
    void FUNC_NAME(const size_t num_elem, const DATA_TYPE alpha,                \
                   const DATA_TYPE *MATRIX_RESTRICT x, const unsigned int incx, \
                   DATA_TYPE *MATRIX_RESTRICT y, const unsigned int incy) {     \
        if (x == NULL || y == NULL) {                                           \
            return;                                                             \
        }                                                                       \
                                                                                \
        if (incx == 0 || incy == 0) {                                           \
            return;                                                             \
        }                                                                       \
                                                                                \
        for (size_t i = 0; i < num_elem; i++) {                                 \
            y[i * incy] += alpha * x[i * incx];                                 \
        }                                                                       \
    }

// General dense matrix multiplication for integer matrices. (?gemm)
// C := alpha * op(A) op(B) + beta * C
// Conventions for ?gemm in the BLAS standard are used.
// Please refer to the documentation for details.
// (e.g. "https://www.intel.com/content/www/us/en
// /docs/onemkl/developer-reference-c/2023-1/cblas-gemm
// -001.html#GUID-97718E5C-6E0A-44F0-B2B1-A551F0F164B2")
#define DEFINE_GEMM(FUNC_NAME, DATA_TYPE)                                                      \
    void FUNC_NAME(const CBLAS_ORDER order, const CBLAS_TRANSPOSE transa,                      \
                   const CBLAS_TRANSPOSE transb,                                               \
                   const unsigned int m, const unsigned int n, const unsigned int k,           \
                   const DATA_TYPE alpha, const DATA_TYPE *MATRIX_RESTRICT a,                  \
                   const unsigned int lda, const DATA_TYPE *MATRIX_RESTRICT b,                 \
                   const unsigned int ldb, const DATA_TYPE beta, DATA_TYPE *MATRIX_RESTRICT c, \
                   const unsigned int ldc) {                                                   \
        /* 'order' is unused and included only for compatibility. Only row-major data */       \
        /* buffers are supported.                                                     */       \
        unsigned int min_lda, min_ldb, min_ldc;                                                \
        bool         a_is_transposed = (transa == CblasTrans);                                 \
        bool         b_is_transposed = (transb == CblasTrans);                                 \
                                                                                               \
        if (a == NULL || b == NULL || c == NULL) {                                             \
            fprintf(stderr, "ERROR! One of the matrices is a null pointer.");                  \
            return;                                                                            \
        }                                                                                      \
        if ((transa == CblasConjTrans || transa == CblasConjNoTrans) ||                        \
            (transb == CblasConjTrans || transb == CblasConjNoTrans)) {                        \
            fprintf(stderr, "ERROR! Conjugate transpose is not supported for integer GEMM.");  \
            return;                                                                            \
        }                                                                                      \
        if (a_is_transposed) {                                                                 \
            min_lda = m > 1 ? m : 1;                                                           \
        } else {                                                                               \
            min_lda = k > 1 ? k : 1;                                                           \
        }                                                                                      \
        if (lda < min_lda) {                                                                   \
            fprintf(stderr, "ERROR! Invalid LDA (leading dimension of A).");                   \
            return;                                                                            \
        }                                                                                      \
        if (b_is_transposed) {                                                                 \
            min_ldb = k > 1 ? k : 1;                                                           \
        } else {                                                                               \
            min_ldb = n > 1 ? n : 1;                                                           \
        }                                                                                      \
        if (ldb < min_ldb) {                                                                   \
            fprintf(stderr, "ERROR! Invalid LDB (leading dimension of B).");                   \
            return;                                                                            \
        }                                                                                      \
        min_ldc = n > 1 ? n : 1;                                                               \
        if (ldc < min_ldc) {                                                                   \
            fprintf(stderr, "ERROR! Invalid LDC (leading dimension of C).");                   \
            return;                                                                            \
        }                                                                                      \
                                                                                               \
        if (m == 0 || n == 0 || ((alpha == 0 || k == 0) && beta == 1)) {                       \
            return;                                                                            \
        }                                                                                      \
                                                                                               \
        for (size_t i = 0; i < m; ++i) {                                                       \
            for (size_t j = 0; j < n; ++j) {                                                   \
                c[i * ldc + j] = (beta == 0) ? 0 : c[i * ldc + j] * beta;                      \
            }                                                                                  \
        }                                                                                      \
                                                                                               \
        if (alpha == 0 || k == 0) {                                                            \
            return;                                                                            \
        }                                                                                      \
                                                                                               \
        for (size_t i = 0; i < m; ++i) {                                                       \
            for (size_t l = 0; l < k; ++l) {                                                   \
                DATA_TYPE a_val = a_is_transposed ? a[l * lda + i] : a[i * lda + l];           \
                for (size_t j = 0; j < n; ++j) {                                               \
                    DATA_TYPE b_val = b_is_transposed ? b[j * ldb + l] : b[l * ldb + j];       \
                    c[i * ldc + j] += alpha * a_val * b_val;                                   \
                }                                                                              \
            }                                                                                  \
        }                                                                                      \
    }

// Sparse BLAS-like function for gathering elements from a
// sparse storage vector to a dense storage vector according
// to supplied indices.
void iusga(const size_t num_elem, const int *y, const unsigned int incy, int *x,
           const unsigned int *idxs) {
    if (x == NULL || y == NULL || idxs == NULL) {
        fprintf(stderr, "ERROR! Null pointer in array argument(s).");
        return;
    }
    for (size_t i = 0; i < num_elem; i++) {
        x[i] = y[idxs[i * incy]];
    }
}

// Sparse BLAS-like function for gathering elements from a
// sparse storage vector to a dense storage vector according
// to supplied indices. Double precision arrays supported.
void dusga(const size_t num_elem, const double *y, const unsigned int incy, double *x,
           const unsigned int *idxs) {
    if (x == NULL || y == NULL || idxs == NULL) {
        fprintf(stderr, "ERROR! Null pointer in array argument(s).");
        return;
    }
    for (size_t i = 0; i < num_elem; i++) {
        x[i] = y[idxs[i * incy]];
    }
}

DEFINE_AXPY(i32_axpy, int32_t)
DEFINE_COPY(i32_copy, int32_t)
DEFINE_SCAL(i32_scal, int32_t)
DEFINE_GEMM(i32_gemm, int32_t)

DEFINE_AXPY(i64_axpy, int64_t)
DEFINE_COPY(i64_copy, int64_t)
DEFINE_SCAL(i64_scal, int64_t)
DEFINE_GEMM(i64_gemm, int64_t)
