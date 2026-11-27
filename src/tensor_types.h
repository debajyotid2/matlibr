#ifndef TENSOR_TYPES_H
#define TENSOR_TYPES_H

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "matrix_types.h"

#ifdef __cplusplus
extern "C" {
#endif

#define TENSOR_MAX_DIMS 8
#define TENSOR_INIT {0}

#define DEFINE_TENSOR_STRUCT(STRUCT_NAME, DATA_TYPE) \
    typedef struct {                                 \
        uint64_t     magic;                          \
        size_t       ndim;                           \
        size_t       size;                           \
        size_t       capacity;                       \
        unsigned int shape[TENSOR_MAX_DIMS];         \
        size_t       strides[TENSOR_MAX_DIMS];       \
        bool         is_view;                        \
        DATA_TYPE   *data;                           \
    } STRUCT_NAME;

DEFINE_TENSOR_STRUCT(TensorI32, int32_t)
DEFINE_TENSOR_STRUCT(TensorI64, int64_t)
DEFINE_TENSOR_STRUCT(TensorF32, float)
DEFINE_TENSOR_STRUCT(TensorF64, double)

#undef DEFINE_TENSOR_STRUCT

// Operation status codes
typedef enum {
    TENSOR_SUCCESS = 0,
    TENSOR_ERR_NULL_PTR,
    TENSOR_ERR_MEMORY_ALLOCATION,
    TENSOR_ERR_DIMENSION_MISMATCH,
    TENSOR_ERR_INVALID_DIMENSION,
    TENSOR_ERR_RANGE_INVALID,
    TENSOR_ERR_RANGE_INVALID_STEP,
    TENSOR_ERR_INDEX_OUT_OF_BOUNDS,
    TENSOR_ERR_TOO_MANY_INTS_TO_GENERATE,
    TENSOR_ERR_ZERO_STD_DEV,
    TENSOR_ERR_INPLACE_MATMUL,
    TENSOR_ERR_UNINITIALIZED,
    TENSOR_ERR_VIEW_TENSOR_OVERWRITE,
    TENSOR_ERR_NON_CONTIGUOUS,
    TENSOR_ERR_INPLACE_GATHER,
    TENSOR_ERR_UNKNOWN
} TensorStatusCode;

static inline size_t _tensor_offset(const void *data, size_t ndim, unsigned int *shape,
                                    size_t *strides, size_t nargs, const size_t *coords,
                                    const char *func, int line, ...) {
    if (data == NULL || shape == NULL || strides == NULL) {
        fprintf(stderr, "TENSOR LIB ERROR: %s,%d: %s\n", func, line, "Null pointer received.");
        abort();
    }
    if (ndim > TENSOR_MAX_DIMS || ndim == 0) {
        fprintf(stderr, "TENSOR LIB ERROR: %s,%d: %s\n", func, line,
                "Number of dimensions out of bounds.");
        abort();
    }
    if (nargs != ndim) {
        fprintf(stderr, "TENSOR LIB ERROR: %s,%d:  Expected %zu indices, got %zu\n", func, line,
                ndim, nargs);
        abort();
    }
    size_t offset = 0;
    for (int i = 0; (size_t)i < ndim; ++i) {
        int64_t dim = coords[i];
        if (dim < 0 || (size_t)dim >= shape[i]) {
            fprintf(stderr, "TENSOR LIB ERROR: %s,%d: Index %ld out of bounds.\n", func, line, dim);
            abort();
        }
        offset += (size_t)dim * strides[i];
    }
    return offset;
}

// Tensor elementwise access
// NOTE: Do not call this with side-effects in place of `tensor`.
#define TENSOR_AT(tensor, ...)                                                                  \
    (((tensor) != NULL ? (tensor)->data : NULL)[_tensor_offset(                                 \
        (tensor) != NULL ? (tensor)->data : NULL, (tensor) != NULL ? (tensor)->ndim : 0,        \
        (tensor) != NULL ? (tensor)->shape : NULL, (tensor) != NULL ? (tensor)->strides : NULL, \
        sizeof((const size_t[]){__VA_ARGS__}) / sizeof(size_t), (const size_t[]){__VA_ARGS__},  \
        __func__, __LINE__)])

#ifdef __cplusplus
}
#endif
#endif // TENSOR_TYPES_H
