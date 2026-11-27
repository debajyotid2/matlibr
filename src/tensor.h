#ifndef TENSOR_H
#define TENSOR_H

#ifdef __cplusplus
extern "C" {
#endif

#ifdef __cplusplus
#    define MATRIX_RESTRICT __restrict
#else
#    define MATRIX_RESTRICT restrict
#endif

#include "matrix_types.h"
#include "tensor_types.h"

#include <stdbool.h>
#include <stddef.h>

TensorStatusCode tensorf64_create(TensorF64 *t, int ndim, const unsigned int *shape);
TensorStatusCode tensorf64_copy(TensorF64 *dest, const TensorF64 *src);
TensorStatusCode tensorf64_reshape(TensorF64 *t, int new_ndim, const unsigned int *new_shape);
TensorStatusCode tensorf64_set_data(TensorF64 *t, const double *src, size_t count);
TensorStatusCode tensorf64_add(TensorF64 *a, const TensorF64 *b);
TensorStatusCode tensorf64_sub(TensorF64 *a, const TensorF64 *b);
TensorStatusCode tensorf64_mul(const TensorF64 *a, bool trans_a, const TensorF64 *b, bool trans_b,
                               TensorF64 *out, const double alpha, const double beta);
TensorStatusCode tensorf64_permute(const TensorF64 *src, TensorF64 *dest, const unsigned int *dims);
TensorStatusCode tensorf64_gather(const TensorF64 *from, TensorF64 *to, const MatrixI32 *idxs,
                                  unsigned int dim);
TensorStatusCode tensorf64_destroy(TensorF64 *t);

#ifdef __cplusplus
}
#endif
#endif // TENSOR_H
