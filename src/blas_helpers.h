#ifndef BLAS_HELPERS_H
#define BLAS_HELPERS_H

#ifdef __cplusplus
extern "C" {
#endif

#include <cblas.h>
#include <stdint.h>

#ifdef __cplusplus
#    define MATRIX_RESTRICT __restrict
#else
#    define MATRIX_RESTRICT restrict
#endif

void i32_scal(const size_t num_elem, const int32_t alpha, int32_t *x, const unsigned int incx);
void i32_gemm(const CBLAS_ORDER order, const CBLAS_TRANSPOSE transa, const CBLAS_TRANSPOSE transb, const unsigned int m,
              const unsigned int n, const unsigned int k, const int32_t alpha,
              const int32_t *MATRIX_RESTRICT a, const unsigned int lda,
              const int32_t *MATRIX_RESTRICT b, const unsigned int ldb, const int32_t beta,
              int32_t *MATRIX_RESTRICT c, const unsigned int ldc);
void i32_axpy(const size_t num_elem, const int32_t alpha, const int32_t *MATRIX_RESTRICT x,
              const unsigned int incx, int32_t *MATRIX_RESTRICT y, const unsigned int incy);
void i32_copy(const size_t num_elem, const int32_t *MATRIX_RESTRICT x, const unsigned int incx,
              int32_t *MATRIX_RESTRICT y, const unsigned int incy);
void i64_scal(const size_t num_elem, const int64_t alpha, int64_t *x, const unsigned int incx);
void i64_gemm(const CBLAS_ORDER order, const CBLAS_TRANSPOSE transa, const CBLAS_TRANSPOSE transb, const unsigned int m,
              const unsigned int n, const unsigned int k, const int64_t alpha,
              const int64_t *MATRIX_RESTRICT a, const unsigned int lda,
              const int64_t *MATRIX_RESTRICT b, const unsigned int ldb, const int64_t beta,
              int64_t *MATRIX_RESTRICT c, const unsigned int ldc);
void i64_axpy(const size_t num_elem, const int64_t alpha, const int64_t *MATRIX_RESTRICT x,
              const unsigned int incx, int64_t *MATRIX_RESTRICT y, const unsigned int incy);
void i64_copy(const size_t num_elem, const int64_t *MATRIX_RESTRICT x, const unsigned int incx,
              int64_t *MATRIX_RESTRICT y, const unsigned int incy);

void iusga(const size_t num_elem, const int *y, const unsigned int incy, int *x,
           const unsigned int *idxs);
void dusga(const size_t num_elem, const double *y, const unsigned int incy, double *x,
           const unsigned int *idxs);

#ifdef __cplusplus
}
#endif

#endif
