#include "tensor.h"
#include "matrix.h"
#include "matrix_types.h"
#include "tensor_types.h"

#include <stdlib.h>
#include <string.h>

#define MATRIX_MAGIC 0xFEAB1365FEDC2937ULL // Keep this value updated with `matrix.c`
#define TENSOR_MAGIC 0xF633DEF2AB795EC1ULL

// Function for displaying and handling errors
static void error_handler(const char *function, int line, TensorStatusCode status_code,
                          const char *msg) {
    fprintf(stderr, "TENSOR LIB ERROR: %s,%d: - error code %d: %s\n", function, line, status_code,
            msg);
}

// clang-format off
// Macro for error handling
#define HANDLE_ERROR(status_code, msg)                  \
error_handler(__func__, __LINE__, status_code, msg)

// Macro to return on error
#define RETURN_ON_ERROR(status_code, msg)               \
    do {                                                \
        HANDLE_ERROR(status_code, msg);                 \
        return (status_code);                           \
    } while (0)

static inline TensorStatusCode tensorstatuscode_from_matrixstatuscode(MatrixStatusCode code) {
    switch (code) { 
        case MATRIX_SUCCESS: return TENSOR_SUCCESS;
        case MATRIX_ERR_NULL_PTR: return TENSOR_ERR_NULL_PTR;
        case MATRIX_ERR_MEMORY_ALLOCATION: return TENSOR_ERR_MEMORY_ALLOCATION;
        case MATRIX_ERR_DIMENSION_MISMATCH: return TENSOR_ERR_DIMENSION_MISMATCH;
        case MATRIX_ERR_INVALID_DIMENSION: return TENSOR_ERR_INVALID_DIMENSION;
        case MATRIX_ERR_RANGE_INVALID: return TENSOR_ERR_RANGE_INVALID;
        case MATRIX_ERR_RANGE_INVALID_STEP: return TENSOR_ERR_RANGE_INVALID_STEP;
        case MATRIX_ERR_INDEX_OUT_OF_BOUNDS: return TENSOR_ERR_INDEX_OUT_OF_BOUNDS;
        case MATRIX_ERR_TOO_MANY_INTS_TO_GENERATE: return TENSOR_ERR_TOO_MANY_INTS_TO_GENERATE;
        case MATRIX_ERR_ZERO_STD_DEV: return TENSOR_ERR_ZERO_STD_DEV;
        case MATRIX_ERR_INPLACE_MATMUL: return TENSOR_ERR_INPLACE_MATMUL;
        case MATRIX_ERR_UNINITIALIZED: return TENSOR_ERR_UNINITIALIZED;
        default: return TENSOR_ERR_UNKNOWN;
    }
}

static inline bool tensor_shapes_equal(
    const unsigned int* shape_a, const unsigned int* shape_b, 
    const size_t ndim_a, const size_t ndim_b) {
    if (!shape_a || !shape_b || ndim_a == 0 || ndim_b == 0 || ndim_a != ndim_b) 
        return false;
    for (int i = 0; (size_t)i < ndim_a; ++i)
        if (shape_a[i] != shape_b[i]) return false;
    return true;
}

// This is strictly for internal use only
#define DEFINE_TENSOR_VIEW_AS_MATRIX(FUNC_NAME, TENSOR_TYPE, MATRIX_TYPE)\
static inline MATRIX_TYPE FUNC_NAME(const TENSOR_TYPE* t, size_t offset, \
                unsigned int nrows, unsigned int ncols) {                \
    if (!t || !t->ndim) return (MATRIX_TYPE)MATRIX_INIT;                 \
    return (MATRIX_TYPE){                                                \
            .magic = MATRIX_MAGIC,                                       \
            .nrows = (size_t)nrows,                                      \
            .ncols = (size_t)ncols,                                      \
            .data = t->data + offset,                                    \
            .capacity = (size_t)nrows * (size_t)ncols,                   \
    };                                                                   \
}

DEFINE_TENSOR_VIEW_AS_MATRIX(tensorf64_view_as_matrix, TensorF64, MatrixF64)
#undef DEFINE_TENSOR_VIEW_AS_MATRIX

// Internal use function: check if a tensor is already initialized with memory.
// If not, allocate memory and initialize it.
#define DEFINE_TENSOR_PREPARE(FUNC_NAME, TENSOR_TYPE, DATA_TYPE)                              \
static TensorStatusCode FUNC_NAME(TENSOR_TYPE* t, size_t ndim, const unsigned int* shape,     \
            bool allow_resize) {                                                              \
    if (!t || !shape) {                                                                       \
        RETURN_ON_ERROR(TENSOR_ERR_NULL_PTR, "Null pointer received.");                       \
    }                                                                                         \
    if (ndim <= 0 || ndim > TENSOR_MAX_DIMS) {                                                \
        RETURN_ON_ERROR(TENSOR_ERR_INVALID_DIMENSION, "Invalid 'ndim'");                      \
    }                                                                                         \
    size_t size = 1;                                                                          \
    for (int i = 0; (size_t)i < ndim; ++i) {                                                  \
        if (shape[i] == 0) {                                                                  \
            RETURN_ON_ERROR(TENSOR_ERR_INVALID_DIMENSION, "Invalid dim.");                    \
        }                                                                                     \
        size *= (size_t)shape[i];                                                             \
    }                                                                                         \
    if (t->magic != TENSOR_MAGIC) {                                                           \
        t->data = (DATA_TYPE*)malloc(size * sizeof(DATA_TYPE));                               \
        if (!t->data) {                                                                       \
            memset(t->shape, (unsigned int)0, ndim * sizeof(unsigned int));                   \
            t->capacity = 0;                                                                  \
            t->size = 0;                                                                      \
            t->ndim = 0;                                                                      \
            RETURN_ON_ERROR(TENSOR_ERR_MEMORY_ALLOCATION, "Allocation failed.");              \
        }                                                                                     \
        t->ndim = (size_t)ndim;                                                               \
        memcpy(t->shape, shape, ndim * sizeof(unsigned int));                                 \
        t->size = size;                                                                       \
        t->capacity = size;                                                                   \
        t->magic = TENSOR_MAGIC;                                                              \
        return TENSOR_SUCCESS;                                                                \
    }                                                                                         \
    if (t->ndim == ndim && memcmp(t->shape, shape, ndim * sizeof(unsigned int)) == 0) {       \
        return TENSOR_SUCCESS;                                                                \
    }                                                                                         \
    if (!allow_resize) {                                                                      \
        RETURN_ON_ERROR(TENSOR_ERR_DIMENSION_MISMATCH, "Wrong dimensions for output tensor.");\
    }                                                                                         \
    if (size > t->capacity) {                                                                 \
        DATA_TYPE* new_data = (DATA_TYPE*)realloc(t->data, size * sizeof(DATA_TYPE));         \
        if (!new_data) {                                                                      \
            RETURN_ON_ERROR(TENSOR_ERR_MEMORY_ALLOCATION, "Reallocation failed.");            \
        }                                                                                     \
        t->data = new_data;                                                                   \
        t->capacity = size;                                                                   \
    }                                                                                         \
    memset(t->shape, 0, TENSOR_MAX_DIMS * sizeof(unsigned int));                              \
    memcpy(t->shape, shape, ndim * sizeof(unsigned int));                                     \
    t->ndim  = (size_t)ndim;                                                                  \
    t->size = size;                                                                           \
    return TENSOR_SUCCESS;                                                                    \
}

#define DEFINE_TENSOR_CREATE(FUNC_NAME, TENSOR_TYPE, DATA_TYPE, PREPARE_FUNC)       \
TensorStatusCode FUNC_NAME(TENSOR_TYPE *t, int ndim, const unsigned int *shape) {   \
    TensorStatusCode ret = PREPARE_FUNC(t, ndim, shape, true);                      \
    if (ret != TENSOR_SUCCESS) {                                                    \
        return ret;                                                                 \
    }                                                                               \
    memset(t->data, 0, t->size * sizeof(DATA_TYPE));                                \
    return TENSOR_SUCCESS;                                                          \
};

#define DEFINE_TENSOR_COPY(FUNC_NAME, TENSOR_TYPE, MATRIX_TYPE, PREPARE_FUNC,           \
                            MATRIX_COPY_FUNC, TENSOR_VIEW_AS_MATRIX_FUNC)               \
TensorStatusCode FUNC_NAME(TENSOR_TYPE *dest, const TENSOR_TYPE *src) {                 \
    if (dest==NULL || !src || !src->data) {                                             \
        RETURN_ON_ERROR(TENSOR_ERR_NULL_PTR, "Null pointer received.");                 \
    }                                                                                   \
    if (dest == src) {                                                                  \
        return TENSOR_SUCCESS;                                                          \
    }                                                                                   \
    TensorStatusCode ret = PREPARE_FUNC(dest, src->ndim, src->shape, true);             \
    if (ret != TENSOR_SUCCESS) {                                                        \
        return ret;                                                                     \
    }                                                                                   \
    const MATRIX_TYPE src_mat = TENSOR_VIEW_AS_MATRIX_FUNC(src, 0, 1, src->size);       \
    MATRIX_TYPE dest_mat = TENSOR_VIEW_AS_MATRIX_FUNC(dest, 0, 1, src->size);           \
    MatrixStatusCode mat_ret = MATRIX_COPY_FUNC(&dest_mat, &src_mat);                   \
    if (mat_ret != MATRIX_SUCCESS) {                                                    \
        return tensorstatuscode_from_matrixstatuscode(mat_ret);                         \
    }                                                                                   \
    return TENSOR_SUCCESS;                                                              \
};

#define DEFINE_TENSOR_RESHAPE(FUNC_NAME, TENSOR_TYPE)                                           \
TensorStatusCode FUNC_NAME(TENSOR_TYPE *t, int new_ndim,                                        \
                            const unsigned int *new_shape) {                                    \
    if (!t || !t->data || !new_shape) {                                                         \
        RETURN_ON_ERROR(TENSOR_ERR_NULL_PTR, "Null pointer received.");                         \
    }                                                                                           \
    if (new_ndim <= 0 || new_ndim > TENSOR_MAX_DIMS) {                                          \
        RETURN_ON_ERROR(TENSOR_ERR_INVALID_DIMENSION, "Invalid 'new_ndim'.");                   \
    }                                                                                           \
    size_t new_size = 1;                                                                        \
    for (int i = 0; i < new_ndim; ++i) {                                                        \
        if (new_shape[i] == 0) {                                                                \
            RETURN_ON_ERROR(TENSOR_ERR_INVALID_DIMENSION, "Invalid dimension in new shape.");   \
        }                                                                                       \
        new_size *= (size_t)new_shape[i];                                                       \
    }                                                                                           \
    if (new_size != t->size) {                                                                  \
        RETURN_ON_ERROR(TENSOR_ERR_DIMENSION_MISMATCH,                                          \
                        "New shape must have same number of elements as old.");                 \
    }                                                                                           \
    t->ndim = (size_t)new_ndim;                                                                 \
    memset(t->shape, 0, TENSOR_MAX_DIMS * sizeof(unsigned int));                                \
    memcpy(t->shape, new_shape, new_ndim * sizeof(unsigned int));                               \
    return TENSOR_SUCCESS;                                                                      \
};

#define DEFINE_TENSOR_SET_DATA(FUNC_NAME, TENSOR_TYPE, DATA_TYPE)                           \
TensorStatusCode FUNC_NAME(TENSOR_TYPE *t, const DATA_TYPE *src, size_t count) {            \
    if (!t || !t->data || !src) {                                                           \
        RETURN_ON_ERROR(TENSOR_ERR_NULL_PTR, "Null pointer received");                      \
    }                                                                                       \
    if (t->size != count) {                                                                 \
        RETURN_ON_ERROR(TENSOR_ERR_DIMENSION_MISMATCH,                                      \
                        "Source must have the same number of elements as destination.");    \
    }                                                                                       \
    memcpy(t->data, src, count * sizeof(DATA_TYPE));                                        \
    return TENSOR_SUCCESS;                                                                  \
};

#define DEFINE_TENSOR_ADD(FUNC_NAME, TENSOR_TYPE, MATRIX_TYPE, MATRIX_ADD_FUNC,     \
                                        TENSOR_VIEW_AS_MATRIX_FUNC)                 \
TensorStatusCode FUNC_NAME(TENSOR_TYPE *a, const TENSOR_TYPE *b) {                  \
    if (!a || !b || !a->data || !b->data) {                                         \
        RETURN_ON_ERROR(TENSOR_ERR_NULL_PTR, "Null pointer received");              \
    }                                                                               \
    if (!tensor_shapes_equal(a->shape, b->shape, a->ndim, b->ndim)) {               \
        RETURN_ON_ERROR(TENSOR_ERR_DIMENSION_MISMATCH,                              \
                        "Both tensors must be of same shape.");                     \
    }                                                                               \
    MATRIX_TYPE a_mat = TENSOR_VIEW_AS_MATRIX_FUNC(a, 0, 1, a->size);               \
    const MATRIX_TYPE b_mat = TENSOR_VIEW_AS_MATRIX_FUNC(b, 0, 1, b->size);         \
    MatrixStatusCode ret = MATRIX_ADD_FUNC(&a_mat, &b_mat);                         \
    if (ret != MATRIX_SUCCESS) {                                                    \
        return tensorstatuscode_from_matrixstatuscode(ret);                         \
    }                                                                               \
    return TENSOR_SUCCESS;                                                          \
};

#define DEFINE_TENSOR_SUB(FUNC_NAME, TENSOR_TYPE, MATRIX_TYPE, MATRIX_SUB_FUNC,     \
                                        TENSOR_VIEW_AS_MATRIX_FUNC)                 \
TensorStatusCode FUNC_NAME(TENSOR_TYPE *a, const TENSOR_TYPE *b) {                  \
    if (!a || !b || !a->data || !b->data) {                                         \
        RETURN_ON_ERROR(TENSOR_ERR_NULL_PTR, "Null pointer received");              \
    }                                                                               \
    if (!tensor_shapes_equal(a->shape, b->shape, a->ndim, b->ndim)) {               \
        RETURN_ON_ERROR(TENSOR_ERR_DIMENSION_MISMATCH,                              \
                        "Both tensors must be of same shape.");                     \
    }                                                                               \
    MATRIX_TYPE a_mat = TENSOR_VIEW_AS_MATRIX_FUNC(a, 0, 1, a->size);               \
    const MATRIX_TYPE b_mat = TENSOR_VIEW_AS_MATRIX_FUNC(b, 0, 1, b->size);         \
    MatrixStatusCode ret = MATRIX_SUB_FUNC(&a_mat, &b_mat);                         \
    if (ret != MATRIX_SUCCESS) {                                                    \
        return tensorstatuscode_from_matrixstatuscode(ret);                         \
    }                                                                               \
    return TENSOR_SUCCESS;                                                          \
};

#define DEFINE_TENSOR_MUL(FUNC_NAME, TENSOR_TYPE, MATRIX_TYPE, PREPARE_FUNC, MATMUL_FUNC, TENSOR_VIEW_AS_MATRIX_FUNC)      \
TensorStatusCode FUNC_NAME(const TENSOR_TYPE *a, bool trans_a,                                                             \
                           const TENSOR_TYPE *b, bool trans_b,                                                             \
                               TENSOR_TYPE *out) {                                                                         \
    if (a==NULL || b==NULL || out==NULL || !a->data || !b->data) {                                                         \
        RETURN_ON_ERROR(TENSOR_ERR_NULL_PTR, "Got null pointer.");                                                         \
    }                                                                                                                      \
    if (out == a || out == b) {                                                                                            \
        RETURN_ON_ERROR(TENSOR_ERR_INPLACE_MATMUL,                                                                         \
                "Inplace tensor multiplication not supported.");                                                           \
    }                                                                                                                      \
    if (a->ndim != b->ndim || a->ndim < 2 || b->ndim < 2) {                                                                \
        RETURN_ON_ERROR(TENSOR_ERR_DIMENSION_MISMATCH,                                                                     \
        "Shapes of A and B do not have correct number of dimensions.");                                                    \
    }                                                                                                                      \
    size_t num_matmuls = 1;                                                                                                \
    for (int i = 0; (size_t)i < a->ndim - 2; ++i) {                                                                        \
        if (a->shape[i] != b->shape[i]) {                                                                                  \
            RETURN_ON_ERROR(TENSOR_ERR_DIMENSION_MISMATCH,                                                                 \
                "Shapes of A and B must have same 'ndim - 2' dimensions.");                                                \
        }                                                                                                                  \
        num_matmuls *= a->shape[i];                                                                                        \
    }                                                                                                                      \
                                                                                                                           \
    size_t a_nrows = a->shape[a->ndim - 2];                                                                                \
    size_t a_ncols = a->shape[a->ndim - 1];                                                                                \
    size_t b_nrows = b->shape[a->ndim - 2];                                                                                \
    size_t b_ncols = b->shape[a->ndim - 1];                                                                                \
                                                                                                                           \
    size_t m = trans_a ? a_ncols: a_nrows;                                                                                 \
    size_t k = trans_a ? a_nrows: a_ncols;                                                                                 \
    size_t k_prime = trans_b ? b_ncols: b_nrows;                                                                           \
    size_t n = trans_b ? b_nrows: b_ncols;                                                                                 \
                                                                                                                           \
    if (k != k_prime) {                                                                                                    \
        RETURN_ON_ERROR(TENSOR_ERR_DIMENSION_MISMATCH,                                                                     \
            "Tensor dimensions must satisfy k = k' for multiplying "                                                       \
            "(..., m, k) and (..., k', n) matrices.");                                                                     \
    }                                                                                                                      \
    unsigned int out_shape[TENSOR_MAX_DIMS] = {0};                                                                         \
    memcpy(&out_shape, a->shape, (a->ndim - 2) * sizeof(unsigned int));                                                    \
    out_shape[a->ndim-2] = (unsigned int)m;                                                                                \
    out_shape[a->ndim-1] = (unsigned int)n;                                                                                \
                                                                                                                           \
    TensorStatusCode res = PREPARE_FUNC(out, a->ndim, out_shape, false);                                                   \
    if (res != TENSOR_SUCCESS) {                                                                                           \
        return res;                                                                                                        \
    }                                                                                                                      \
                                                                                                                           \
    MatrixStatusCode mat_res;                                                                                              \
    for (size_t matmul_idx = 0; matmul_idx < num_matmuls; ++matmul_idx) {                                                  \
        MATRIX_TYPE temp_a = TENSOR_VIEW_AS_MATRIX_FUNC(a, matmul_idx * a_nrows*a_ncols, a_nrows, a_ncols);                \
        MATRIX_TYPE temp_b = TENSOR_VIEW_AS_MATRIX_FUNC(b, matmul_idx * b_nrows*b_ncols, b_nrows, b_ncols);                \
        MATRIX_TYPE temp_out = TENSOR_VIEW_AS_MATRIX_FUNC(out, matmul_idx * m * n, m, n);                                  \
                                                                                                                           \
        mat_res = MATMUL_FUNC(&temp_a, trans_a, &temp_b, trans_b, &temp_out);                                              \
        if (mat_res != MATRIX_SUCCESS) return tensorstatuscode_from_matrixstatuscode(mat_res);                             \
    }                                                                                                                      \
    return TENSOR_SUCCESS;                                                                                                 \
}

#define DEFINE_TENSOR_DESTROY(FUNC_NAME, TENSOR_TYPE)                                               \
TensorStatusCode FUNC_NAME(TENSOR_TYPE *t) {                                                        \
    if (!t || !t->data) {                                                                           \
        RETURN_ON_ERROR(TENSOR_ERR_NULL_PTR, "Null pointer received.");                             \
    }                                                                                               \
    if (t->magic != TENSOR_MAGIC) {                                                                 \
        RETURN_ON_ERROR(TENSOR_ERR_UNINITIALIZED, "Tensor not initialized or already destroyed");   \
    }                                                                                               \
    if (t->data != NULL) {                                                                          \
        free(t->data);                                                                              \
    }                                                                                               \
    t->data = NULL;                                                                                 \
    t->capacity = 0;                                                                                \
    t->magic = 0;                                                                                   \
    memset(t->shape, (unsigned int)0, t->ndim * sizeof(unsigned int));                              \
    t->ndim = 0;                                                                                    \
    t->size = 0;                                                                                    \
    return TENSOR_SUCCESS;                                                                          \
}

/**************** TensorF64 operations ***********************/
DEFINE_TENSOR_PREPARE(tensorf64_prepare, TensorF64, double)
DEFINE_TENSOR_CREATE(tensorf64_create, TensorF64, double, tensorf64_prepare)
DEFINE_TENSOR_COPY(tensorf64_copy, TensorF64, MatrixF64, tensorf64_prepare, 
                   matf64_copy, tensorf64_view_as_matrix)
DEFINE_TENSOR_RESHAPE(tensorf64_reshape, TensorF64)
DEFINE_TENSOR_SET_DATA(tensorf64_set_data, TensorF64, double)
DEFINE_TENSOR_ADD(tensorf64_add, TensorF64, MatrixF64, matf64_add, 
                    tensorf64_view_as_matrix)
DEFINE_TENSOR_SUB(tensorf64_sub, TensorF64, MatrixF64, matf64_sub, 
                    tensorf64_view_as_matrix)
DEFINE_TENSOR_MUL(tensorf64_mul, TensorF64, MatrixF64, tensorf64_prepare,
                  matf64_mul, tensorf64_view_as_matrix)
DEFINE_TENSOR_DESTROY(tensorf64_destroy, TensorF64)
