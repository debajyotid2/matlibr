#ifndef MATRIX_TYPES_H
#define MATRIX_TYPES_H

#include <stdio.h>
#include <stdlib.h>

// Matrices must be initialized to MATRIX_INIT, otherwise the initialization may cause segmentation
// faults.
#define MATRIX_INIT {0}

// Matrix for double precision data
typedef struct {
    size_t  nrows, ncols;
    double *data;
    size_t  capacity;
} Matrix;

// Matrix for integer data
typedef struct {
    size_t nrows, ncols;
    int   *data;
    size_t capacity;
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
    MATRIX_ERR_ZERO_STD_DEV,
    MATRIX_ERR_INPLACE_MATMUL,
} MatrixStatusCode;

static inline size_t _matrix_offset(const void *matrix, size_t nrows, size_t ncols, int row,
                                    int col, const char *func, int line) {
    if (matrix == NULL) {
        fprintf(stderr, "MATRIX LIB ERROR: %s,%d: %s\n", func, line, "Null pointer received.");
        abort();
    }
    if ((int)row >= (int)nrows || (int)row < 0) {
        fprintf(stderr, "MATRIX LIB ERROR: %s,%d: %s\n", func, line, "Row index out of bounds.");
        abort();
    } else if ((int)col < 0 || (int)col >= (int)ncols) {
        fprintf(stderr, "MATRIX LIB ERROR: %s,%d: %s\n", func, line, "Column index out of bounds.");
        abort();
    }
    return ncols * (size_t)row + (size_t)col;
}

// Matrix elementwise access
// NOTE: Do not call this with an expression in place of `matrix`, otherwise the result might be
// undefined as the `matrix` argument is evaluated four times in the macro expansion.
#define MATRIX_AT(matrix, row, col)                                                               \
    (((matrix) != NULL ? (matrix)->data : NULL)[_matrix_offset(                                   \
        (matrix), (matrix) != NULL ? (matrix)->nrows : 0, (matrix) != NULL ? (matrix)->ncols : 0, \
        (row), (col), __func__, __LINE__)])

#endif // MATRIX_TYPES_H
