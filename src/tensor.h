#ifndef TENSOR_H
#define TENSOR_H

#ifdef __cplusplus
extern "C" {
#endif

#include "matrix_types.h"
#include "tensor_types.h"

#include <stddef.h>
#include <stdbool.h>

TensorStatusCode tensorf64_create(TensorF64 *t, int ndim, const unsigned int *shape);
TensorStatusCode tensorf64_copy(TensorF64 *dest, const TensorF64 *src);
TensorStatusCode tensorf64_reshape(TensorF64 *t, int new_ndim, const unsigned int *new_shape);
TensorStatusCode tensorf64_set_data(TensorF64 *t, const double *src, size_t count);
TensorStatusCode tensorf64_add(TensorF64 *a, const TensorF64 *b);
TensorStatusCode tensorf64_sub(TensorF64 *a, const TensorF64 *b);
TensorStatusCode tensorf64_mul(const TensorF64 *a, bool trans_a, const TensorF64 *b, bool trans_b,
                               TensorF64 *out);
TensorStatusCode tensorf64_destroy(TensorF64 *t);

#ifdef __cplusplus
}
#endif
#endif // TENSOR_H
