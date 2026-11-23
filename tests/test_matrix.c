/* Tests for module matrix.h

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

#include "matrix.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

// Helper Functions for Testing

// Checks if two IntMatrix instances are identical.
bool are_int_matrices_equal(const IntMatrix *a, const IntMatrix *b) {
    if (a == NULL || b == NULL)
        return false;
    if (a->nrows != b->nrows || a->ncols != b->ncols)
        return false;
    return memcmp(a->data, b->data, a->nrows * a->ncols * sizeof(int)) == 0;
}

// Checks if two double Matrix instances are equal within a tolerance.
bool are_double_matrices_equal(const Matrix *a, const Matrix *b, double epsilon) {
    if (a == NULL || b == NULL)
        return false;
    if (a->nrows != b->nrows || a->ncols != b->ncols)
        return false;
    for (size_t i = 0; i < a->nrows * a->ncols; ++i) {
        if (fabs(a->data[i] - b->data[i]) > epsilon)
            return false;
    }
    return true;
}

// Comparison function for qsort with integers.
int compare_ints(const void *a, const void *b) {
    int arg1 = *(const int *)a;
    int arg2 = *(const int *)b;
    if (arg1 < arg2)
        return -1;
    if (arg1 > arg2)
        return 1;
    return 0;
}

// Checks if two MatrixI32 instances are identical.
bool are_i32_matrices_equal(const MatrixI32 *a, const MatrixI32 *b) {
    if (a == NULL || b == NULL)
        return false;
    if (a->nrows != b->nrows || a->ncols != b->ncols)
        return false;
    return memcmp(a->data, b->data, a->nrows * a->ncols * sizeof(int32_t)) == 0;
}

// Checks if two MatrixI64 instances are identical.
bool are_i64_matrices_equal(const MatrixI64 *a, const MatrixI64 *b) {
    if (a == NULL || b == NULL)
        return false;
    if (a->nrows != b->nrows || a->ncols != b->ncols)
        return false;
    return memcmp(a->data, b->data, a->nrows * a->ncols * sizeof(int64_t)) == 0;
}

// Checks if two float MatrixF32 instances are equal within a tolerance.
bool are_f32_matrices_equal(const MatrixF32 *a, const MatrixF32 *b, float epsilon) {
    if (a == NULL || b == NULL)
        return false;
    if (a->nrows != b->nrows || a->ncols != b->ncols)
        return false;
    for (size_t i = 0; i < a->nrows * a->ncols; ++i) {
        if (fabsf(a->data[i] - b->data[i]) > epsilon)
            return false;
    }
    return true;
}

// Checks if two double MatrixF64 instances are equal within a tolerance.
bool are_f64_matrices_equal(const MatrixF64 *a, const MatrixF64 *b, double epsilon) {
    if (a == NULL || b == NULL)
        return false;
    if (a->nrows != b->nrows || a->ncols != b->ncols)
        return false;
    for (size_t i = 0; i < a->nrows * a->ncols; ++i) {
        if (fabs(a->data[i] - b->data[i]) > epsilon)
            return false;
    }
    return true;
}

// Comparison function for qsort with int32_t.
int compare_i32s(const void *a, const void *b) {
    int32_t arg1 = *(const int32_t *)a;
    int32_t arg2 = *(const int32_t *)b;
    if (arg1 < arg2)
        return -1;
    if (arg1 > arg2)
        return 1;
    return 0;
}

// Comparison function for qsort with int64_t.
int compare_i64s(const void *a, const void *b) {
    int64_t arg1 = *(const int64_t *)a;
    int64_t arg2 = *(const int64_t *)b;
    if (arg1 < arg2)
        return -1;
    if (arg1 > arg2)
        return 1;
    return 0;
}

TEST_CASE("MatrixI32: Creation and Destruction", "[matrixi32]") {
    MatrixI32 mat = MATRIX_INIT;

    SECTION("Creation with valid dimensions") {
        REQUIRE(mati32_create(&mat, 3, 4) == MATRIX_SUCCESS);
        REQUIRE(mat.nrows == 3);
        REQUIRE(mat.ncols == 4);
        REQUIRE(mat.data != NULL);
        for (size_t i = 0; i < 12; ++i)
            REQUIRE(mat.data[i] == 0);
        mati32_destroy(&mat);
    }

    SECTION("Creation with invalid dimensions") {
        REQUIRE(mati32_create(&mat, 0, 5) == MATRIX_ERR_INVALID_DIMENSION);
        REQUIRE(mati32_create(&mat, 5, 0) == MATRIX_ERR_INVALID_DIMENSION);
        REQUIRE(mati32_create(&mat, -1, 5) == MATRIX_ERR_INVALID_DIMENSION);
    }

    SECTION("Creation with null pointer") {
        REQUIRE(mati32_create(NULL, 3, 3) == MATRIX_ERR_NULL_PTR);
    }

    SECTION("Creation without MATRIX_INIT, destruction and reuse") {
        MatrixI32 mat_b;
        REQUIRE(mati32_create(&mat_b, 3, 3) == MATRIX_SUCCESS);
        REQUIRE(mati32_destroy(&mat_b) == MATRIX_SUCCESS);
        REQUIRE(mati32_create(&mat_b, 2, 4) == MATRIX_SUCCESS);
        REQUIRE(mati32_destroy(&mat_b) == MATRIX_SUCCESS);
    }

    SECTION("Destruction") {
        mati32_create(&mat, 2, 2);
        REQUIRE(mati32_destroy(&mat) == MATRIX_SUCCESS);
        REQUIRE(mat.data == NULL);
        REQUIRE(mati32_destroy(NULL) == MATRIX_ERR_NULL_PTR);
    }
}

TEST_CASE("MatrixI32: Core Operations", "[matrixi32]") {
    MatrixI32 mat_a = MATRIX_INIT;
    MatrixI32 result = MATRIX_INIT;
    mati32_create(&mat_a, 2, 2);
    mati32_create(&result, 2, 2);
    mat_a.data[0] = 1;
    mat_a.data[1] = 2;
    mat_a.data[2] = 3;
    mat_a.data[3] = 4;

    SECTION("Copy") {
        REQUIRE(mati32_copy(&result, &mat_a) == MATRIX_SUCCESS);
        REQUIRE(are_i32_matrices_equal(&mat_a, &result));
        REQUIRE(mat_a.data != result.data);
    }

    SECTION("Assign") {
        int     row = 0, col = 1;
        int32_t value = 42;
        REQUIRE(mati32_assign(&mat_a, row, col, value) == MATRIX_SUCCESS);
        REQUIRE(mat_a.data[row * mat_a.ncols + col] == value);
    }

    SECTION("Fill") {
        REQUIRE(mati32_fill(&mat_a, 7) == MATRIX_SUCCESS);
        for (int i = 0; i < 4; ++i)
            REQUIRE(mat_a.data[i] == 7);
    }

    SECTION("Scale") {
        REQUIRE(mati32_scale(&mat_a, 2) == MATRIX_SUCCESS);
        int32_t expected_data[] = {2, 4, 6, 8};
        for (int i = 0; i < 4; ++i)
            REQUIRE(mat_a.data[i] == expected_data[i]);
    }

    SECTION("Add Scalar") {
        REQUIRE(mati32_add_scalar(&mat_a, 10) == MATRIX_SUCCESS);
        int32_t expected_data[] = {11, 12, 13, 14};
        for (int i = 0; i < 4; ++i)
            REQUIRE(mat_a.data[i] == expected_data[i]);
    }

    mati32_destroy(&result);
    mati32_destroy(&mat_a);
}

TEST_CASE("MatrixI32: Set Data", "[matrixi32]") {
    MatrixI32 mat = MATRIX_INIT;
    mati32_create(&mat, 2, 3);

    SECTION("Assign flat array in row-major order") {
        const int32_t src[6] = {10, 20, 30, 40, 50, 60};
        REQUIRE(mati32_set_data(&mat, src, 6) == MATRIX_SUCCESS);

        for (size_t i = 0; i < 6; ++i) {
            REQUIRE(mat.data[i] == src[i]);
        }
        // Verify 2D row-major mapping (2 rows x 3 cols)
        REQUIRE(MATRIX_AT(&mat, 0, 2) == 30);
        REQUIRE(MATRIX_AT(&mat, 1, 0) == 40);
        REQUIRE(MATRIX_AT(&mat, 1, 2) == 60);
    }

    SECTION("Overwrite existing matrix data") {
        mati32_fill(&mat, 99);
        const int32_t src[6] = {-1, -2, -3, -4, -5, -6};

        REQUIRE(mati32_set_data(&mat, src, 6) == MATRIX_SUCCESS);
        for (size_t i = 0; i < 6; ++i) {
            REQUIRE(mat.data[i] == src[i]);
        }
    }

    mati32_destroy(&mat);
}

TEST_CASE("MatrixI32: Set Data Error Cases", "[matrixi32]") {
    MatrixI32 mat = MATRIX_INIT;
    mati32_create(&mat, 2, 3);
    const int32_t src[6] = {1, 2, 3, 4, 5, 6};

    SECTION("Null pointer handling") {
        MatrixI32 mat_b = MATRIX_INIT;

        REQUIRE(mati32_set_data(NULL, src, 6) == MATRIX_ERR_NULL_PTR);
        REQUIRE(mati32_set_data(&mat_b, src, 6) == MATRIX_ERR_NULL_PTR);
        REQUIRE(mati32_set_data(&mat, NULL, 6) == MATRIX_ERR_NULL_PTR);
    }

    SECTION("Element count mismatch") {
        mati32_fill(&mat, 42);

        REQUIRE(mati32_set_data(&mat, src, 5) == MATRIX_ERR_DIMENSION_MISMATCH);
        REQUIRE(mati32_set_data(&mat, src, 7) == MATRIX_ERR_DIMENSION_MISMATCH);
        REQUIRE(mati32_set_data(&mat, src, 0) == MATRIX_ERR_DIMENSION_MISMATCH);

        for (size_t i = 0; i < 6; ++i) {
            REQUIRE(mat.data[i] == 42);
        }
    }

    mati32_destroy(&mat);
}

TEST_CASE("MatrixI32: Arithmetic Operations", "[matrixi32]") {
    MatrixI32 mat_a = MATRIX_INIT;
    MatrixI32 mat_b = MATRIX_INIT;
    MatrixI32 result = MATRIX_INIT;
    MatrixI32 expected_mat = MATRIX_INIT;
    mati32_create(&mat_a, 2, 2);
    mati32_create(&mat_b, 2, 2);
    mat_a.data[0] = 1;
    mat_a.data[1] = 2;
    mat_a.data[2] = 3;
    mat_a.data[3] = 4;
    mat_b.data[0] = 5;
    mat_b.data[1] = 6;
    mat_b.data[2] = 7;
    mat_b.data[3] = 8;

    SECTION("Addition") {
        REQUIRE(mati32_add(&mat_a, &mat_b) == MATRIX_SUCCESS);
        mati32_create(&expected_mat, 2, 2);
        int32_t expected_data[] = {6, 8, 10, 12};
        memcpy(expected_mat.data, expected_data, 4 * sizeof(int32_t));
        REQUIRE(are_i32_matrices_equal(&mat_a, &expected_mat));
        mati32_destroy(&expected_mat);
    }

    SECTION("Subtraction") {
        REQUIRE(mati32_sub(&mat_a, &mat_b) == MATRIX_SUCCESS);
        mati32_create(&expected_mat, 2, 2);
        int32_t expected_data[] = {-4, -4, -4, -4};
        memcpy(expected_mat.data, expected_data, 4 * sizeof(int32_t));
        REQUIRE(are_i32_matrices_equal(&mat_a, &expected_mat));
        mati32_destroy(&expected_mat);
    }

    SECTION("Multiplication") {
        mati32_create(&result, 2, 2);
        REQUIRE(mati32_mul(&mat_a, false, &mat_b, false, &result) == MATRIX_SUCCESS);
        mati32_create(&expected_mat, 2, 2);
        int32_t expected_data[] = {19, 22, 43, 50};
        memcpy(expected_mat.data, expected_data, 4 * sizeof(int32_t));
        REQUIRE(are_i32_matrices_equal(&result, &expected_mat));
        mati32_destroy(&result);
        mati32_destroy(&expected_mat);
    }

    mati32_destroy(&mat_a);
    mati32_destroy(&mat_b);
}

TEST_CASE("MatrixI32: Vector Operations", "[matrixi32]") {
    MatrixI32 mat = MATRIX_INIT;
    MatrixI32 vec = MATRIX_INIT;
    MatrixI32 result = MATRIX_INIT;
    MatrixI32 expected_mat = MATRIX_INIT;

    mati32_create(&mat, 2, 3);
    for (int i = 0; i < 6; ++i)
        mat.data[i] = 1;

    SECTION("Add row vector") {
        mati32_create(&vec, 1, 3);
        vec.data[0] = 1;
        vec.data[1] = 2;
        vec.data[2] = 3;
        REQUIRE(mati32_vec_add(&mat, &vec) == MATRIX_SUCCESS);
        mati32_create(&expected_mat, 2, 3);
        int32_t expected_data[] = {2, 3, 4, 2, 3, 4};
        memcpy(expected_mat.data, expected_data, 6 * sizeof(int32_t));
        REQUIRE(are_i32_matrices_equal(&mat, &expected_mat));
        mati32_destroy(&vec);
        mati32_destroy(&expected_mat);
    }

    SECTION("Range") {
        REQUIRE(mati32_range(&result, 0, 10, 2, 1) == MATRIX_SUCCESS);
        mati32_create(&expected_mat, 1, 5);
        int32_t expected_data[] = {0, 2, 4, 6, 8};
        memcpy(expected_mat.data, expected_data, 5 * sizeof(int32_t));
        REQUIRE(are_i32_matrices_equal(&result, &expected_mat));
        mati32_destroy(&result);
        mati32_destroy(&expected_mat);
    }

    SECTION("Repeat") {
        mati32_create(&vec, 1, 2);
        vec.data[0] = 5;
        vec.data[1] = 10;
        REQUIRE(mati32_repeat(&result, &vec, 0, 3) == MATRIX_SUCCESS);
        mati32_create(&expected_mat, 3, 2);
        int32_t expected_data[] = {5, 10, 5, 10, 5, 10};
        memcpy(expected_mat.data, expected_data, 6 * sizeof(int32_t));
        REQUIRE(are_i32_matrices_equal(&result, &expected_mat));
        mati32_destroy(&vec);
        mati32_destroy(&result);
        mati32_destroy(&expected_mat);
    }

    mati32_destroy(&mat);
}

TEST_CASE("MatrixI64: Creation and Destruction", "[matrixi64]") {
    MatrixI64 mat = MATRIX_INIT;

    SECTION("Creation with valid dimensions") {
        REQUIRE(mati64_create(&mat, 3, 4) == MATRIX_SUCCESS);
        REQUIRE(mat.nrows == 3);
        REQUIRE(mat.ncols == 4);
        REQUIRE(mat.data != NULL);
        for (size_t i = 0; i < 12; ++i)
            REQUIRE(mat.data[i] == 0);
        mati64_destroy(&mat);
    }

    SECTION("Creation with invalid dimensions") {
        REQUIRE(mati64_create(&mat, 0, 5) == MATRIX_ERR_INVALID_DIMENSION);
        REQUIRE(mati64_create(&mat, 5, 0) == MATRIX_ERR_INVALID_DIMENSION);
        REQUIRE(mati64_create(&mat, -1, 5) == MATRIX_ERR_INVALID_DIMENSION);
    }

    SECTION("Creation with null pointer") {
        REQUIRE(mati64_create(NULL, 3, 3) == MATRIX_ERR_NULL_PTR);
    }

    SECTION("Destruction") {
        mati64_create(&mat, 2, 2);
        REQUIRE(mati64_destroy(&mat) == MATRIX_SUCCESS);
        REQUIRE(mat.data == NULL);
        REQUIRE(mati64_destroy(NULL) == MATRIX_ERR_NULL_PTR);
    }
}

TEST_CASE("MatrixI64: Core Operations", "[matrixi64]") {
    MatrixI64 mat_a = MATRIX_INIT;
    MatrixI64 result = MATRIX_INIT;
    mati64_create(&mat_a, 2, 2);
    mati64_create(&result, 2, 2);
    mat_a.data[0] = 1;
    mat_a.data[1] = 2;
    mat_a.data[2] = 3;
    mat_a.data[3] = 4;

    SECTION("Copy") {
        REQUIRE(mati64_copy(&result, &mat_a) == MATRIX_SUCCESS);
        REQUIRE(are_i64_matrices_equal(&mat_a, &result));
        REQUIRE(mat_a.data != result.data);
    }

    SECTION("Assign") {
        int     row = 0, col = 1;
        int64_t value = 42;
        REQUIRE(mati64_assign(&mat_a, row, col, value) == MATRIX_SUCCESS);
        REQUIRE(mat_a.data[row * mat_a.ncols + col] == value);
    }

    SECTION("Fill") {
        REQUIRE(mati64_fill(&mat_a, 7) == MATRIX_SUCCESS);
        for (int i = 0; i < 4; ++i)
            REQUIRE(mat_a.data[i] == 7);
    }

    SECTION("Scale") {
        REQUIRE(mati64_scale(&mat_a, 2) == MATRIX_SUCCESS);
        int64_t expected_data[] = {2, 4, 6, 8};
        for (int i = 0; i < 4; ++i)
            REQUIRE(mat_a.data[i] == expected_data[i]);
    }

    SECTION("Add Scalar") {
        REQUIRE(mati64_add_scalar(&mat_a, 10) == MATRIX_SUCCESS);
        int64_t expected_data[] = {11, 12, 13, 14};
        for (int i = 0; i < 4; ++i)
            REQUIRE(mat_a.data[i] == expected_data[i]);
    }

    mati64_destroy(&result);
    mati64_destroy(&mat_a);
}

TEST_CASE("MatrixI64: Set Data", "[matrixi64]") {
    MatrixI64 mat = MATRIX_INIT;
    mati64_create(&mat, 2, 3);

    SECTION("Assign flat array in row-major order") {
        const int64_t src[6] = {10, 20, 30, 40, 50, 60};
        REQUIRE(mati64_set_data(&mat, src, 6) == MATRIX_SUCCESS);

        for (size_t i = 0; i < 6; ++i) {
            REQUIRE(mat.data[i] == src[i]);
        }
        // Verify 2D row-major mapping (2 rows x 3 cols)
        REQUIRE(MATRIX_AT(&mat, 0, 2) == 30);
        REQUIRE(MATRIX_AT(&mat, 1, 0) == 40);
        REQUIRE(MATRIX_AT(&mat, 1, 2) == 60);
    }

    SECTION("Overwrite existing matrix data") {
        mati64_fill(&mat, 99);
        const int64_t src[6] = {-1, -2, -3, -4, -5, -6};

        REQUIRE(mati64_set_data(&mat, src, 6) == MATRIX_SUCCESS);
        for (size_t i = 0; i < 6; ++i) {
            REQUIRE(mat.data[i] == src[i]);
        }
    }

    mati64_destroy(&mat);
}

TEST_CASE("MatrixI64: Set Data Error Cases", "[matrixi64]") {
    MatrixI64 mat = MATRIX_INIT;
    mati64_create(&mat, 2, 3);
    const int64_t src[6] = {1, 2, 3, 4, 5, 6};

    SECTION("Null pointer handling") {
        MatrixI64 mat_b = MATRIX_INIT;

        REQUIRE(mati64_set_data(NULL, src, 6) == MATRIX_ERR_NULL_PTR);
        REQUIRE(mati64_set_data(&mat_b, src, 6) == MATRIX_ERR_NULL_PTR);
        REQUIRE(mati64_set_data(&mat, NULL, 6) == MATRIX_ERR_NULL_PTR);
    }

    SECTION("Element count mismatch") {
        mati64_fill(&mat, 42);

        REQUIRE(mati64_set_data(&mat, src, 5) == MATRIX_ERR_DIMENSION_MISMATCH);
        REQUIRE(mati64_set_data(&mat, src, 7) == MATRIX_ERR_DIMENSION_MISMATCH);
        REQUIRE(mati64_set_data(&mat, src, 0) == MATRIX_ERR_DIMENSION_MISMATCH);

        for (size_t i = 0; i < 6; ++i) {
            REQUIRE(mat.data[i] == 42);
        }
    }

    mati64_destroy(&mat);
}

TEST_CASE("MatrixI64: Arithmetic Operations", "[matrixi64]") {
    MatrixI64 mat_a = MATRIX_INIT;
    MatrixI64 mat_b = MATRIX_INIT;
    MatrixI64 result = MATRIX_INIT;
    MatrixI64 expected_mat = MATRIX_INIT;
    mati64_create(&mat_a, 2, 2);
    mati64_create(&mat_b, 2, 2);
    mat_a.data[0] = 1;
    mat_a.data[1] = 2;
    mat_a.data[2] = 3;
    mat_a.data[3] = 4;
    mat_b.data[0] = 5;
    mat_b.data[1] = 6;
    mat_b.data[2] = 7;
    mat_b.data[3] = 8;

    SECTION("Addition") {
        REQUIRE(mati64_add(&mat_a, &mat_b) == MATRIX_SUCCESS);
        mati64_create(&expected_mat, 2, 2);
        int64_t expected_data[] = {6, 8, 10, 12};
        memcpy(expected_mat.data, expected_data, 4 * sizeof(int64_t));
        REQUIRE(are_i64_matrices_equal(&mat_a, &expected_mat));
        mati64_destroy(&expected_mat);
    }

    SECTION("Subtraction") {
        REQUIRE(mati64_sub(&mat_a, &mat_b) == MATRIX_SUCCESS);
        mati64_create(&expected_mat, 2, 2);
        int64_t expected_data[] = {-4, -4, -4, -4};
        memcpy(expected_mat.data, expected_data, 4 * sizeof(int64_t));
        REQUIRE(are_i64_matrices_equal(&mat_a, &expected_mat));
        mati64_destroy(&expected_mat);
    }

    SECTION("Multiplication") {
        mati64_create(&result, 2, 2);
        REQUIRE(mati64_mul(&mat_a, false, &mat_b, false, &result) == MATRIX_SUCCESS);
        mati64_create(&expected_mat, 2, 2);
        int64_t expected_data[] = {19, 22, 43, 50};
        memcpy(expected_mat.data, expected_data, 4 * sizeof(int64_t));
        REQUIRE(are_i64_matrices_equal(&result, &expected_mat));
        mati64_destroy(&result);
        mati64_destroy(&expected_mat);
    }

    mati64_destroy(&mat_a);
    mati64_destroy(&mat_b);
}

TEST_CASE("MatrixI64: Vector Operations", "[matrixi64]") {
    MatrixI64 mat = MATRIX_INIT;
    MatrixI64 vec = MATRIX_INIT;
    MatrixI64 result = MATRIX_INIT;
    MatrixI64 expected_mat = MATRIX_INIT;

    mati64_create(&mat, 2, 3);
    for (int i = 0; i < 6; ++i)
        mat.data[i] = 1;

    SECTION("Add row vector") {
        mati64_create(&vec, 1, 3);
        vec.data[0] = 1;
        vec.data[1] = 2;
        vec.data[2] = 3;
        REQUIRE(mati64_vec_add(&mat, &vec) == MATRIX_SUCCESS);
        mati64_create(&expected_mat, 2, 3);
        int64_t expected_data[] = {2, 3, 4, 2, 3, 4};
        memcpy(expected_mat.data, expected_data, 6 * sizeof(int64_t));
        REQUIRE(are_i64_matrices_equal(&mat, &expected_mat));
        mati64_destroy(&vec);
        mati64_destroy(&expected_mat);
    }

    SECTION("Range") {
        REQUIRE(mati64_range(&result, 0, 10, 2, 1) == MATRIX_SUCCESS);
        mati64_create(&expected_mat, 1, 5);
        int64_t expected_data[] = {0, 2, 4, 6, 8};
        memcpy(expected_mat.data, expected_data, 5 * sizeof(int64_t));
        REQUIRE(are_i64_matrices_equal(&result, &expected_mat));
        mati64_destroy(&result);
        mati64_destroy(&expected_mat);
    }

    SECTION("Repeat") {
        mati64_create(&vec, 1, 2);
        vec.data[0] = 5;
        vec.data[1] = 10;
        REQUIRE(mati64_repeat(&result, &vec, 0, 3) == MATRIX_SUCCESS);
        mati64_create(&expected_mat, 3, 2);
        int64_t expected_data[] = {5, 10, 5, 10, 5, 10};
        memcpy(expected_mat.data, expected_data, 6 * sizeof(int64_t));
        REQUIRE(are_i64_matrices_equal(&result, &expected_mat));
        mati64_destroy(&vec);
        mati64_destroy(&result);
        mati64_destroy(&expected_mat);
    }

    mati64_destroy(&mat);
}

TEST_CASE("MatrixF32: Creation and Core Ops", "[matrixf32]") {
    MatrixF32 mat = MATRIX_INIT;
    matf32_create(&mat, 2, 2);
    REQUIRE(matf32_fill(&mat, 3.5f) == MATRIX_SUCCESS);
    for (int i = 0; i < 4; ++i) {
        REQUIRE_THAT(mat.data[i], Catch::Matchers::WithinAbs(3.5f, 1e-5));
    }

    SECTION("Assign") {
        int   row = 0, col = 1;
        float value = 42.0f;
        REQUIRE(matf32_assign(&mat, row, col, value) == MATRIX_SUCCESS);
        REQUIRE_THAT(mat.data[row * mat.ncols + col], Catch::Matchers::WithinAbs(value, 1e-5));
    }

    matf32_destroy(&mat);
}

TEST_CASE("MatrixF32: Scalar and Arithmetic Operations", "[matrixf32]") {
    MatrixF32 mat = MATRIX_INIT;
    matf32_create(&mat, 2, 3);
    for (int i = 0; i < 6; ++i)
        mat.data[i] = 1.5f;

    SECTION("Add Scalar") {
        REQUIRE(matf32_add_scalar(&mat, 10.5f) == MATRIX_SUCCESS);
        for (int i = 0; i < 6; ++i) {
            REQUIRE_THAT(mat.data[i], Catch::Matchers::WithinAbs(12.0f, 1e-5));
        }
    }

    SECTION("Addition") {
        MatrixF32 mat_b = MATRIX_INIT;
        matf32_create(&mat_b, 2, 3);
        matf32_fill(&mat_b, 0.5f);
        REQUIRE(matf32_add(&mat, &mat_b) == MATRIX_SUCCESS);
        for (int i = 0; i < 6; ++i) {
            REQUIRE_THAT(mat.data[i], Catch::Matchers::WithinAbs(2.0f, 1e-5));
        }
        matf32_destroy(&mat_b);
    }

    matf32_destroy(&mat);
}

TEST_CASE("MatrixF32: Vector Operations", "[matrixf32]") {
    MatrixF32 mat = MATRIX_INIT;
    MatrixF32 vec_row = MATRIX_INIT;
    MatrixF32 vec_col = MATRIX_INIT;
    MatrixF32 result = MATRIX_INIT;
    MatrixF32 expected_mat = MATRIX_INIT;
    MatrixI32 indices = MATRIX_INIT;
    float     tol = 1.0e-5f;

    matf32_create(&mat, 2, 3);
    for (int i = 0; i < 6; ++i)
        mat.data[i] = 1.5f;

    matf32_create(&vec_row, 1, 3);
    vec_row.data[0] = 1.0f;
    vec_row.data[1] = 2.0f;
    vec_row.data[2] = 3.0f;

    matf32_create(&vec_col, 2, 1);
    vec_col.data[0] = 10.0f;
    vec_col.data[1] = 20.0f;

    SECTION("Add row vector") {
        REQUIRE(matf32_vec_add(&mat, &vec_row) == MATRIX_SUCCESS);
        matf32_create(&expected_mat, 2, 3);
        float expected_data[] = {2.5f, 3.5f, 4.5f, 2.5f, 3.5f, 4.5f};
        memcpy(expected_mat.data, expected_data, 6 * sizeof(float));
        REQUIRE(are_f32_matrices_equal(&mat, &expected_mat, tol));
        matf32_destroy(&expected_mat);
    }

    SECTION("Subtract column vector") {
        REQUIRE(matf32_vec_sub(&mat, &vec_col) == MATRIX_SUCCESS);
        matf32_create(&expected_mat, 2, 3);
        float expected_data[] = {-8.5f, -8.5f, -8.5f, -18.5f, -18.5f, -18.5f};
        memcpy(expected_mat.data, expected_data, 6 * sizeof(float));
        REQUIRE(are_f32_matrices_equal(&mat, &expected_mat, tol));
        matf32_destroy(&expected_mat);
    }

    SECTION("Range") {
        REQUIRE(matf32_range(&result, 0.0f, 1.0f, 0.25f, 1) == MATRIX_SUCCESS);
        matf32_create(&expected_mat, 1, 4);
        float expected_data[] = {0.0f, 0.25f, 0.5f, 0.75f};
        memcpy(expected_mat.data, expected_data, 4 * sizeof(float));
        REQUIRE(are_f32_matrices_equal(&result, &expected_mat, tol));
        matf32_destroy(&result);
        matf32_destroy(&expected_mat);
    }

    SECTION("Repeat") {
        REQUIRE(matf32_repeat(&result, &vec_row, 0, 2) == MATRIX_SUCCESS);
        matf32_create(&expected_mat, 2, 3);
        float expected_data[] = {1.0f, 2.0f, 3.0f, 1.0f, 2.0f, 3.0f};
        memcpy(expected_mat.data, expected_data, 6 * sizeof(float));
        REQUIRE(are_f32_matrices_equal(&result, &expected_mat, tol));
        matf32_destroy(&result);
        matf32_destroy(&expected_mat);
    }

    SECTION("Gather") {
        mati32_create(&indices, 2, 1);
        indices.data[0] = 2;
        indices.data[1] = 0;
        matf32_create(&result, 2, 2);
        REQUIRE(matf32_gather(&mat, &result, &indices, 1) == MATRIX_SUCCESS);
        matf32_create(&expected_mat, 2, 2);
        float expected_data[] = {1.5f, 1.5f, 1.5f, 1.5f};
        memcpy(expected_mat.data, expected_data, 4 * sizeof(float));
        REQUIRE(are_f32_matrices_equal(&result, &expected_mat, tol));

        indices.data[0] = 3;
        REQUIRE(matf32_gather(&mat, &result, &indices, 1) == MATRIX_ERR_INDEX_OUT_OF_BOUNDS);
        mati32_destroy(&indices);
        matf32_destroy(&result);
        matf32_destroy(&expected_mat);
    }

    matf32_destroy(&mat);
    matf32_destroy(&vec_row);
    matf32_destroy(&vec_col);
}

TEST_CASE("MatrixF32: Set Data", "[matrixf32]") {
    MatrixF32 mat = MATRIX_INIT;
    matf32_create(&mat, 3, 2);

    SECTION("Fill with flat array") {
        const float src[6] = {1.5f, -2.25f, 3.0f, 4.75f, -5.5f, 6.125f};
        REQUIRE(matf32_set_data(&mat, src, 6) == MATRIX_SUCCESS);

        for (size_t i = 0; i < 6; ++i) {
            REQUIRE(mat.data[i] == src[i]);
        }
        REQUIRE(MATRIX_AT(&mat, 0, 1) == -2.25f);
        REQUIRE(MATRIX_AT(&mat, 1, 0) == 3.0f);
        REQUIRE(MATRIX_AT(&mat, 2, 1) == 6.125f);
    }

    SECTION("Overwrite existing matrix data") {
        matf32_fill(&mat, 7.5f);
        const float src[6] = {0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f};

        REQUIRE(matf32_set_data(&mat, src, 6) == MATRIX_SUCCESS);
        for (size_t i = 0; i < 6; ++i) {
            REQUIRE(mat.data[i] == src[i]);
        }
    }

    matf32_destroy(&mat);
}

TEST_CASE("MatrixF32: Set Data Error Cases", "[matrixf32]") {
    MatrixF32 mat = MATRIX_INIT;
    matf32_create(&mat, 3, 2);
    const float src[6] = {1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f};

    SECTION("Null pointer handling") {
        MatrixF32 mat_b = MATRIX_INIT;

        REQUIRE(matf32_set_data(NULL, src, 6) == MATRIX_ERR_NULL_PTR);
        REQUIRE(matf32_set_data(&mat_b, src, 6) == MATRIX_ERR_NULL_PTR);
        REQUIRE(matf32_set_data(&mat, NULL, 6) == MATRIX_ERR_NULL_PTR);
    }

    SECTION("Element count mismatch") {
        matf32_fill(&mat, 3.14f);

        REQUIRE(matf32_set_data(&mat, src, 4) == MATRIX_ERR_DIMENSION_MISMATCH);
        REQUIRE(matf32_set_data(&mat, src, 10) == MATRIX_ERR_DIMENSION_MISMATCH);
        REQUIRE(matf32_set_data(&mat, src, 0) == MATRIX_ERR_DIMENSION_MISMATCH);

        for (size_t i = 0; i < 6; ++i) {
            REQUIRE(mat.data[i] == 3.14f);
        }
    }

    matf32_destroy(&mat);
}

TEST_CASE("MatrixF64: Creation and Core Ops", "[matrixf64]") {
    MatrixF64 mat = MATRIX_INIT;
    matf64_create(&mat, 2, 2);
    REQUIRE(matf64_fill(&mat, 3.5) == MATRIX_SUCCESS);
    for (int i = 0; i < 4; ++i) {
        REQUIRE_THAT(mat.data[i], Catch::Matchers::WithinAbs(3.5, 1e-9));
    }

    SECTION("Assign") {
        int    row = 0, col = 1;
        double value = 42.0;
        REQUIRE(matf64_assign(&mat, row, col, value) == MATRIX_SUCCESS);
        REQUIRE_THAT(mat.data[row * mat.ncols + col], Catch::Matchers::WithinAbs(value, 1e-9));
    }

    matf64_destroy(&mat);
}

TEST_CASE("MatrixF64: Scalar and Arithmetic Operations", "[matrixf64]") {
    MatrixF64 mat = MATRIX_INIT;
    matf64_create(&mat, 2, 3);
    for (int i = 0; i < 6; ++i)
        mat.data[i] = 1.5;

    SECTION("Add Scalar") {
        REQUIRE(matf64_add_scalar(&mat, 10.5) == MATRIX_SUCCESS);
        for (int i = 0; i < 6; ++i) {
            REQUIRE_THAT(mat.data[i], Catch::Matchers::WithinAbs(12.0, 1e-9));
        }
    }

    SECTION("Addition") {
        MatrixF64 mat_b = MATRIX_INIT;
        matf64_create(&mat_b, 2, 3);
        matf64_fill(&mat_b, 0.5);
        REQUIRE(matf64_add(&mat, &mat_b) == MATRIX_SUCCESS);
        for (int i = 0; i < 6; ++i) {
            REQUIRE_THAT(mat.data[i], Catch::Matchers::WithinAbs(2.0, 1e-9));
        }
        matf64_destroy(&mat_b);
    }

    matf64_destroy(&mat);
}

TEST_CASE("MatrixF64: Vector Operations", "[matrixf64]") {
    MatrixF64 mat = MATRIX_INIT;
    MatrixF64 vec_row = MATRIX_INIT;
    MatrixF64 vec_col = MATRIX_INIT;
    MatrixF64 result = MATRIX_INIT;
    MatrixF64 expected_mat = MATRIX_INIT;
    MatrixI32 indices = MATRIX_INIT;
    double    tol = 1.0e-9;

    matf64_create(&mat, 2, 3);
    for (int i = 0; i < 6; ++i)
        mat.data[i] = 1.5;

    matf64_create(&vec_row, 1, 3);
    vec_row.data[0] = 1.0;
    vec_row.data[1] = 2.0;
    vec_row.data[2] = 3.0;

    matf64_create(&vec_col, 2, 1);
    vec_col.data[0] = 10.0;
    vec_col.data[1] = 20.0;

    SECTION("Add row vector") {
        REQUIRE(matf64_vec_add(&mat, &vec_row) == MATRIX_SUCCESS);
        matf64_create(&expected_mat, 2, 3);
        double expected_data[] = {2.5, 3.5, 4.5, 2.5, 3.5, 4.5};
        memcpy(expected_mat.data, expected_data, 6 * sizeof(double));
        REQUIRE(are_f64_matrices_equal(&mat, &expected_mat, tol));
        matf64_destroy(&expected_mat);
    }

    SECTION("Subtract column vector") {
        REQUIRE(matf64_vec_sub(&mat, &vec_col) == MATRIX_SUCCESS);
        matf64_create(&expected_mat, 2, 3);
        double expected_data[] = {-8.5, -8.5, -8.5, -18.5, -18.5, -18.5};
        memcpy(expected_mat.data, expected_data, 6 * sizeof(double));
        REQUIRE(are_f64_matrices_equal(&mat, &expected_mat, tol));
        matf64_destroy(&expected_mat);
    }

    SECTION("Range") {
        REQUIRE(matf64_range(&result, 0.0, 1.0, 0.25, 1) == MATRIX_SUCCESS);
        matf64_create(&expected_mat, 1, 4);
        double expected_data[] = {0.0, 0.25, 0.5, 0.75};
        memcpy(expected_mat.data, expected_data, 4 * sizeof(double));
        REQUIRE(are_f64_matrices_equal(&result, &expected_mat, tol));
        matf64_destroy(&result);
        matf64_destroy(&expected_mat);
    }

    SECTION("Repeat") {
        REQUIRE(matf64_repeat(&result, &vec_row, 0, 2) == MATRIX_SUCCESS);
        matf64_create(&expected_mat, 2, 3);
        double expected_data[] = {1.0, 2.0, 3.0, 1.0, 2.0, 3.0};
        memcpy(expected_mat.data, expected_data, 6 * sizeof(double));
        REQUIRE(are_f64_matrices_equal(&result, &expected_mat, tol));
        matf64_destroy(&result);
        matf64_destroy(&expected_mat);
    }

    SECTION("Gather") {
        mati32_create(&indices, 2, 1);
        indices.data[0] = 2;
        indices.data[1] = 0;
        matf64_create(&result, 2, 2);
        REQUIRE(matf64_gather(&mat, &result, &indices, 1) == MATRIX_SUCCESS);
        matf64_create(&expected_mat, 2, 2);
        double expected_data[] = {1.5, 1.5, 1.5, 1.5};
        memcpy(expected_mat.data, expected_data, 4 * sizeof(double));
        REQUIRE(are_f64_matrices_equal(&result, &expected_mat, tol));

        indices.data[0] = 3;
        REQUIRE(matf64_gather(&mat, &result, &indices, 1) == MATRIX_ERR_INDEX_OUT_OF_BOUNDS);
        mati32_destroy(&indices);
        matf64_destroy(&result);
        matf64_destroy(&expected_mat);
    }

    matf64_destroy(&mat);
    matf64_destroy(&vec_row);
    matf64_destroy(&vec_col);
}

TEST_CASE("MatrixF64: Set Data", "[matrixf64]") {
    MatrixF64 mat = MATRIX_INIT;
    matf64_create(&mat, 3, 2);

    SECTION("Fill with flat array") {
        const double src[6] = {1.5, -2.25, 3.0, 4.75, -5.5, 6.125};
        REQUIRE(matf64_set_data(&mat, src, 6) == MATRIX_SUCCESS);

        for (size_t i = 0; i < 6; ++i) {
            REQUIRE(mat.data[i] == src[i]);
        }
        REQUIRE(MATRIX_AT(&mat, 0, 1) == -2.25);
        REQUIRE(MATRIX_AT(&mat, 1, 0) == 3.0);
        REQUIRE(MATRIX_AT(&mat, 2, 1) == 6.125);
    }

    SECTION("Overwrite existing matrix data") {
        matf64_fill(&mat, 7.5);
        const double src[6] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6};

        REQUIRE(matf64_set_data(&mat, src, 6) == MATRIX_SUCCESS);
        for (size_t i = 0; i < 6; ++i) {
            REQUIRE(mat.data[i] == src[i]);
        }
    }

    matf64_destroy(&mat);
}

TEST_CASE("MatrixF64: Set Data Error Cases", "[matrixf64]") {
    MatrixF64 mat = MATRIX_INIT;
    matf64_create(&mat, 3, 2);
    const double src[6] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};

    SECTION("Null pointer handling") {
        MatrixF64 mat_b = MATRIX_INIT;

        REQUIRE(matf64_set_data(NULL, src, 6) == MATRIX_ERR_NULL_PTR);
        REQUIRE(matf64_set_data(&mat_b, src, 6) == MATRIX_ERR_NULL_PTR);
        REQUIRE(matf64_set_data(&mat, NULL, 6) == MATRIX_ERR_NULL_PTR);
    }

    SECTION("Element count mismatch") {
        matf64_fill(&mat, 3.14);

        REQUIRE(matf64_set_data(&mat, src, 4) == MATRIX_ERR_DIMENSION_MISMATCH);
        REQUIRE(matf64_set_data(&mat, src, 10) == MATRIX_ERR_DIMENSION_MISMATCH);
        REQUIRE(matf64_set_data(&mat, src, 0) == MATRIX_ERR_DIMENSION_MISMATCH);

        for (size_t i = 0; i < 6; ++i) {
            REQUIRE(mat.data[i] == 3.14);
        }
    }

    matf64_destroy(&mat);
}

TEST_CASE("IntMatrix: Creation and Destruction", "[intmatrix]") {
    IntMatrix mat = MATRIX_INIT;

    SECTION("Creation with valid dimensions") {
        REQUIRE(intmat_create(&mat, 3, 4) == MATRIX_SUCCESS);
        REQUIRE(mat.nrows == 3);
        REQUIRE(mat.ncols == 4);
        REQUIRE(mat.data != NULL);
        for (size_t i = 0; i < 12; ++i)
            REQUIRE(mat.data[i] == 0);
        intmat_destroy(&mat);
    }

    SECTION("Creation with invalid dimensions") {
        REQUIRE(intmat_create(&mat, 0, 5) == MATRIX_ERR_INVALID_DIMENSION);
        REQUIRE(intmat_create(&mat, 5, 0) == MATRIX_ERR_INVALID_DIMENSION);
        REQUIRE(intmat_create(&mat, -1, 5) == MATRIX_ERR_INVALID_DIMENSION);
    }

    SECTION("Creation with null pointer") {
        REQUIRE(intmat_create(NULL, 3, 3) == MATRIX_ERR_NULL_PTR);
    }

    SECTION("Destruction") {
        intmat_create(&mat, 2, 2);
        REQUIRE(intmat_destroy(&mat) == MATRIX_SUCCESS);
        REQUIRE(mat.data == NULL);
        REQUIRE(intmat_destroy(NULL) == MATRIX_ERR_NULL_PTR);
    }
}

TEST_CASE("IntMatrix: Core Operations", "[intmatrix]") {
    IntMatrix mat_a = MATRIX_INIT;
    IntMatrix result = MATRIX_INIT;
    intmat_create(&mat_a, 2, 2);
    intmat_create(&result, 2, 2);
    mat_a.data[0] = 1;
    mat_a.data[1] = 2;
    mat_a.data[2] = 3;
    mat_a.data[3] = 4;

    SECTION("Copy") {
        REQUIRE(intmat_copy(&result, &mat_a) == MATRIX_SUCCESS);
        REQUIRE(are_int_matrices_equal(&mat_a, &result));
        REQUIRE(mat_a.data != result.data);
    }

    SECTION("Assign") {
        int row = 0, col = 1;
        int value = 42;
        REQUIRE(intmat_assign(&mat_a, row, col, value) == MATRIX_SUCCESS);
        REQUIRE(mat_a.data[row * mat_a.ncols + col] == value);
    }

    SECTION("Fill") {
        REQUIRE(intmat_fill(&mat_a, 7) == MATRIX_SUCCESS);
        for (int i = 0; i < 4; ++i)
            REQUIRE(mat_a.data[i] == 7);
    }

    SECTION("Scale") {
        REQUIRE(intmat_scale(&mat_a, 2) == MATRIX_SUCCESS);
        int expected_data[] = {2, 4, 6, 8};
        for (int i = 0; i < 4; ++i)
            REQUIRE(mat_a.data[i] == expected_data[i]);
    }

    SECTION("Add Scalar") {
        REQUIRE(intmat_add_scalar(&mat_a, 10) == MATRIX_SUCCESS);
        int expected_data[] = {11, 12, 13, 14};
        for (int i = 0; i < 4; ++i)
            REQUIRE(mat_a.data[i] == expected_data[i]);
    }

    intmat_destroy(&result);
    intmat_destroy(&mat_a);
}

TEST_CASE("IntMatrix: Random Fill", "[intmatrix]") {
    IntMatrix mat = MATRIX_INIT;
    intmat_create(&mat, 4, 5);

    SECTION("Fill random with replacement") {
        REQUIRE(intmat_fill_random(&mat, 1, 10, true) == MATRIX_SUCCESS);
        for (size_t i = 0; i < 20; ++i) {
            REQUIRE(mat.data[i] >= 1);
            REQUIRE(mat.data[i] < 10);
        }
    }

    SECTION("Fill random without replacement") {
        REQUIRE(intmat_fill_random(&mat, 0, 20, false) == MATRIX_SUCCESS);

        size_t num_elements = mat.nrows * mat.ncols;
        int   *temp_array = (int *)malloc(num_elements * sizeof(int));
        REQUIRE(temp_array != NULL);
        memcpy(temp_array, mat.data, num_elements * sizeof(int));

        qsort(temp_array, num_elements, sizeof(int), compare_ints);

        bool duplicate_found = false;
        for (size_t i = 0; i < num_elements - 1; ++i) {
            if (temp_array[i] == temp_array[i + 1]) {
                duplicate_found = true;
                break;
            }
        }
        free(temp_array);
        REQUIRE_FALSE(duplicate_found);
    }

    SECTION("Error cases for random fill") {
        REQUIRE(intmat_fill_random(&mat, 10, 5, false) == MATRIX_ERR_RANGE_INVALID);
        REQUIRE(intmat_fill_random(&mat, 0, 19, false) == MATRIX_ERR_TOO_MANY_INTS_TO_GENERATE);
    }

    intmat_destroy(&mat);
}

TEST_CASE("IntMatrix: Set Data", "[intmatrix]") {
    IntMatrix mat = MATRIX_INIT;
    intmat_create(&mat, 2, 3);

    SECTION("Assign flat array in row-major order") {
        const int src[6] = {10, 20, 30, 40, 50, 60};
        REQUIRE(intmat_set_data(&mat, src, 6) == MATRIX_SUCCESS);

        for (size_t i = 0; i < 6; ++i) {
            REQUIRE(mat.data[i] == src[i]);
        }
        // Verify 2D row-major mapping (2 rows x 3 cols)
        REQUIRE(MATRIX_AT(&mat, 0, 2) == 30);
        REQUIRE(MATRIX_AT(&mat, 1, 0) == 40);
        REQUIRE(MATRIX_AT(&mat, 1, 2) == 60);
    }

    SECTION("Overwrite existing matrix data") {
        intmat_fill(&mat, 99);
        const int src[6] = {-1, -2, -3, -4, -5, -6};

        REQUIRE(intmat_set_data(&mat, src, 6) == MATRIX_SUCCESS);
        for (size_t i = 0; i < 6; ++i) {
            REQUIRE(mat.data[i] == src[i]);
        }
    }

    intmat_destroy(&mat);
}

TEST_CASE("IntMatrix: Set Data Error Cases", "[intmatrix]") {
    IntMatrix mat = MATRIX_INIT;
    intmat_create(&mat, 2, 3);
    const int src[6] = {1, 2, 3, 4, 5, 6};

    SECTION("Null pointer handling") {
        IntMatrix mat_b = MATRIX_INIT;

        REQUIRE(intmat_set_data(NULL, src, 6) == MATRIX_ERR_NULL_PTR);
        REQUIRE(intmat_set_data(&mat_b, src, 6) == MATRIX_ERR_NULL_PTR);
        REQUIRE(intmat_set_data(&mat, NULL, 6) == MATRIX_ERR_NULL_PTR);
    }

    SECTION("Element count mismatch") {
        intmat_fill(&mat, 42);

        REQUIRE(intmat_set_data(&mat, src, 5) == MATRIX_ERR_DIMENSION_MISMATCH);
        REQUIRE(intmat_set_data(&mat, src, 7) == MATRIX_ERR_DIMENSION_MISMATCH);
        REQUIRE(intmat_set_data(&mat, src, 0) == MATRIX_ERR_DIMENSION_MISMATCH);

        for (size_t i = 0; i < 6; ++i) {
            REQUIRE(mat.data[i] == 42);
        }
    }

    intmat_destroy(&mat);
}

TEST_CASE("IntMatrix: Arithmetic Operations", "[intmatrix]") {
    IntMatrix mat_a = MATRIX_INIT;
    IntMatrix mat_b = MATRIX_INIT;
    IntMatrix result = MATRIX_INIT;
    IntMatrix expected_mat = MATRIX_INIT;
    intmat_create(&mat_a, 2, 2);
    intmat_create(&mat_b, 2, 2);
    mat_a.data[0] = 1;
    mat_a.data[1] = 2;
    mat_a.data[2] = 3;
    mat_a.data[3] = 4;
    mat_b.data[0] = 5;
    mat_b.data[1] = 6;
    mat_b.data[2] = 7;
    mat_b.data[3] = 8;

    SECTION("Addition") {
        REQUIRE(intmat_add(&mat_a, &mat_b) == MATRIX_SUCCESS);
        intmat_create(&expected_mat, 2, 2);
        int expected_data[] = {6, 8, 10, 12};
        memcpy(expected_mat.data, expected_data, 4 * sizeof(int));
        REQUIRE(are_int_matrices_equal(&mat_a, &expected_mat));
        intmat_destroy(&expected_mat);
    }

    SECTION("Subtraction") {
        REQUIRE(intmat_sub(&mat_a, &mat_b) == MATRIX_SUCCESS);
        intmat_create(&expected_mat, 2, 2);
        int expected_data[] = {-4, -4, -4, -4};
        memcpy(expected_mat.data, expected_data, 4 * sizeof(int));
        REQUIRE(are_int_matrices_equal(&mat_a, &expected_mat));
        intmat_destroy(&expected_mat);
    }

    SECTION("Multiplication") {
        intmat_create(&result, 2, 2);
        REQUIRE(intmat_mul(&mat_a, false, &mat_b, false, &result) == MATRIX_SUCCESS);
        intmat_create(&expected_mat, 2, 2);
        int expected_data[] = {19, 22, 43, 50};
        memcpy(expected_mat.data, expected_data, 4 * sizeof(int));
        REQUIRE(are_int_matrices_equal(&result, &expected_mat));
        intmat_destroy(&result);
        intmat_destroy(&expected_mat);
    }

    intmat_destroy(&mat_a);
    intmat_destroy(&mat_b);
}

TEST_CASE("IntMatrix: Vector Operations", "[intmatrix]") {
    IntMatrix mat = MATRIX_INIT;
    IntMatrix vec = MATRIX_INIT;
    IntMatrix result = MATRIX_INIT;
    IntMatrix expected_mat = MATRIX_INIT;

    intmat_create(&mat, 2, 3);
    for (int i = 0; i < 6; ++i)
        mat.data[i] = 1;

    SECTION("Add row vector") {
        intmat_create(&vec, 1, 3);
        vec.data[0] = 1;
        vec.data[1] = 2;
        vec.data[2] = 3;
        REQUIRE(intmat_vec_add(&mat, &vec) == MATRIX_SUCCESS);
        intmat_create(&expected_mat, 2, 3);
        int expected_data[] = {2, 3, 4, 2, 3, 4};
        memcpy(expected_mat.data, expected_data, 6 * sizeof(int));
        REQUIRE(are_int_matrices_equal(&mat, &expected_mat));
        intmat_destroy(&vec);
        intmat_destroy(&expected_mat);
    }

    SECTION("Range") {
        REQUIRE(intmat_range(&result, 0, 10, 2, 1) == MATRIX_SUCCESS);
        intmat_create(&expected_mat, 1, 5);
        int expected_data[] = {0, 2, 4, 6, 8};
        memcpy(expected_mat.data, expected_data, 5 * sizeof(int));
        REQUIRE(are_int_matrices_equal(&result, &expected_mat));
        intmat_destroy(&result);
        intmat_destroy(&expected_mat);
    }

    SECTION("Repeat") {
        intmat_create(&vec, 1, 2);
        vec.data[0] = 5;
        vec.data[1] = 10;
        REQUIRE(intmat_repeat(&result, &vec, 0, 3) == MATRIX_SUCCESS);
        intmat_create(&expected_mat, 3, 2);
        int expected_data[] = {5, 10, 5, 10, 5, 10};
        memcpy(expected_mat.data, expected_data, 6 * sizeof(int));
        REQUIRE(are_int_matrices_equal(&result, &expected_mat));
        intmat_destroy(&vec);
        intmat_destroy(&result);
        intmat_destroy(&expected_mat);
    }

    intmat_destroy(&mat);
}

TEST_CASE("Matrix: Creation and Core Ops", "[matrix]") {
    Matrix mat = MATRIX_INIT;
    mat_create(&mat, 2, 2);
    REQUIRE(mat_fill(&mat, 3.5) == MATRIX_SUCCESS);
    for (int i = 0; i < 4; ++i) {
        REQUIRE_THAT(mat.data[i], Catch::Matchers::WithinAbs(3.5, 1e-9));
    }

    SECTION("Assign") {
        int    row = 0, col = 1;
        double value = 42.0;
        REQUIRE(mat_assign(&mat, row, col, value) == MATRIX_SUCCESS);
        REQUIRE_THAT(mat.data[row * mat.ncols + col], Catch::Matchers::WithinAbs(value, 1e-9));
    }

    mat_destroy(&mat);
}

TEST_CASE("Matrix: Scalar and Arithmetic Operations", "[matrix]") {
    Matrix mat = MATRIX_INIT;
    mat_create(&mat, 2, 3);
    for (int i = 0; i < 6; ++i)
        mat.data[i] = 1.5;

    SECTION("Add Scalar") {
        REQUIRE(mat_add_scalar(&mat, 10.5) == MATRIX_SUCCESS);
        for (int i = 0; i < 6; ++i) {
            REQUIRE_THAT(mat.data[i], Catch::Matchers::WithinAbs(12.0, 1e-9));
        }
    }

    SECTION("Addition") {
        Matrix mat_b = MATRIX_INIT;
        mat_create(&mat_b, 2, 3);
        mat_fill(&mat_b, 0.5);
        REQUIRE(mat_add(&mat, &mat_b) == MATRIX_SUCCESS);
        for (int i = 0; i < 6; ++i) {
            REQUIRE_THAT(mat.data[i], Catch::Matchers::WithinAbs(2.0, 1e-9));
        }
        mat_destroy(&mat_b);
    }

    mat_destroy(&mat);
}

TEST_CASE("Matrix: Vector Operations", "[matrix]") {
    Matrix    mat = MATRIX_INIT;
    Matrix    vec_row = MATRIX_INIT;
    Matrix    vec_col = MATRIX_INIT;
    Matrix    result = MATRIX_INIT;
    Matrix    expected_mat = MATRIX_INIT;
    IntMatrix indices = MATRIX_INIT;
    double    tol = 1.0e-9;

    mat_create(&mat, 2, 3);
    for (int i = 0; i < 6; ++i)
        mat.data[i] = 1.5;

    mat_create(&vec_row, 1, 3);
    vec_row.data[0] = 1.0;
    vec_row.data[1] = 2.0;
    vec_row.data[2] = 3.0;

    mat_create(&vec_col, 2, 1);
    vec_col.data[0] = 10.0;
    vec_col.data[1] = 20.0;

    SECTION("Add row vector") {
        REQUIRE(mat_vec_add(&mat, &vec_row) == MATRIX_SUCCESS);
        mat_create(&expected_mat, 2, 3);
        double expected_data[] = {2.5, 3.5, 4.5, 2.5, 3.5, 4.5};
        memcpy(expected_mat.data, expected_data, 6 * sizeof(double));
        REQUIRE(are_double_matrices_equal(&mat, &expected_mat, tol));
        mat_destroy(&expected_mat);
    }

    SECTION("Subtract column vector") {
        REQUIRE(mat_vec_sub(&mat, &vec_col) == MATRIX_SUCCESS);
        mat_create(&expected_mat, 2, 3);
        double expected_data[] = {-8.5, -8.5, -8.5, -18.5, -18.5, -18.5};
        memcpy(expected_mat.data, expected_data, 6 * sizeof(double));
        REQUIRE(are_double_matrices_equal(&mat, &expected_mat, tol));
        mat_destroy(&expected_mat);
    }

    SECTION("Range") {
        REQUIRE(mat_range(&result, 0.0, 1.0, 0.25, 1) == MATRIX_SUCCESS);
        mat_create(&expected_mat, 1, 4);
        double expected_data[] = {0.0, 0.25, 0.5, 0.75};
        memcpy(expected_mat.data, expected_data, 4 * sizeof(double));
        REQUIRE(are_double_matrices_equal(&result, &expected_mat, tol));
        mat_destroy(&result);
        mat_destroy(&expected_mat);
    }

    SECTION("Repeat") {
        REQUIRE(mat_repeat(&result, &vec_row, 0, 2) == MATRIX_SUCCESS);
        mat_create(&expected_mat, 2, 3);
        double expected_data[] = {1.0, 2.0, 3.0, 1.0, 2.0, 3.0};
        memcpy(expected_mat.data, expected_data, 6 * sizeof(double));
        REQUIRE(are_double_matrices_equal(&result, &expected_mat, tol));
        mat_destroy(&result);
        mat_destroy(&expected_mat);
    }

    SECTION("Gather") {
        intmat_create(&indices, 2, 1);
        indices.data[0] = 2;
        indices.data[1] = 0;
        mat_create(&result, 2, 2);
        REQUIRE(mat_gather(&mat, &result, &indices, 1) == MATRIX_SUCCESS);
        mat_create(&expected_mat, 2, 2);
        double expected_data[] = {1.5, 1.5, 1.5, 1.5};
        memcpy(expected_mat.data, expected_data, 4 * sizeof(double));
        REQUIRE(are_double_matrices_equal(&result, &expected_mat, tol));

        indices.data[0] = 3;
        REQUIRE(mat_gather(&mat, &result, &indices, 1) == MATRIX_ERR_INDEX_OUT_OF_BOUNDS);
        intmat_destroy(&indices);
        mat_destroy(&result);
        mat_destroy(&expected_mat);
    }

    mat_destroy(&mat);
    mat_destroy(&vec_row);
    mat_destroy(&vec_col);
}

TEST_CASE("Matrix: BLAS-based Functions", "[matrix]") {
    Matrix mat = MATRIX_INIT;
    mat_create(&mat, 2, 3);
    mat.data[0] = -1.0;
    mat.data[1] = 2.0;
    mat.data[2] = -3.0;
    mat.data[3] = 4.0;
    mat.data[4] = -5.0;
    mat.data[5] = 6.0;

    SECTION("Absolute Sum (dasum)") {
        double sum = 0.0;
        REQUIRE(mat_abs_sum(&sum, &mat) == MATRIX_SUCCESS);
        REQUIRE_THAT(sum, Catch::Matchers::WithinAbs(21.0, 1e-9));
    }

    SECTION("Euclidean Norm (dnrm2)") {
        double norm = 0.0;
        REQUIRE(mat_norm(&norm, &mat) == MATRIX_SUCCESS);
        double expected_norm_sq = 1.0 + 4.0 + 9.0 + 16.0 + 25.0 + 36.0; // 91.0
        REQUIRE_THAT(norm, Catch::Matchers::WithinAbs(sqrt(expected_norm_sq), 1e-9));
    }

    mat_destroy(&mat);
}

TEST_CASE("Matrix: Set Data", "[matrix]") {
    Matrix mat = MATRIX_INIT;
    mat_create(&mat, 3, 2);

    SECTION("Fill with flat array") {
        const double src[6] = {1.5, -2.25, 3.0, 4.75, -5.5, 6.125};
        REQUIRE(mat_set_data(&mat, src, 6) == MATRIX_SUCCESS);

        for (size_t i = 0; i < 6; ++i) {
            REQUIRE(mat.data[i] == src[i]);
        }
        REQUIRE(MATRIX_AT(&mat, 0, 1) == -2.25);
        REQUIRE(MATRIX_AT(&mat, 1, 0) == 3.0);
        REQUIRE(MATRIX_AT(&mat, 2, 1) == 6.125);
    }

    SECTION("Overwrite existing matrix data") {
        mat_fill(&mat, 7.5);
        const double src[6] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6};

        REQUIRE(mat_set_data(&mat, src, 6) == MATRIX_SUCCESS);
        for (size_t i = 0; i < 6; ++i) {
            REQUIRE(mat.data[i] == src[i]);
        }
    }

    mat_destroy(&mat);
}

TEST_CASE("Matrix: Set Data Error Cases", "[matrix]") {
    Matrix mat = MATRIX_INIT;
    mat_create(&mat, 3, 2);
    const double src[6] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};

    SECTION("Null pointer handling") {
        Matrix mat_b = MATRIX_INIT;

        REQUIRE(mat_set_data(NULL, src, 6) == MATRIX_ERR_NULL_PTR);
        REQUIRE(mat_set_data(&mat_b, src, 6) == MATRIX_ERR_NULL_PTR);
        REQUIRE(mat_set_data(&mat, NULL, 6) == MATRIX_ERR_NULL_PTR);
    }

    SECTION("Element count mismatch") {
        mat_fill(&mat, 3.14);

        REQUIRE(mat_set_data(&mat, src, 4) == MATRIX_ERR_DIMENSION_MISMATCH);
        REQUIRE(mat_set_data(&mat, src, 10) == MATRIX_ERR_DIMENSION_MISMATCH);
        REQUIRE(mat_set_data(&mat, src, 0) == MATRIX_ERR_DIMENSION_MISMATCH);

        for (size_t i = 0; i < 6; ++i) {
            REQUIRE(mat.data[i] == 3.14);
        }
    }

    mat_destroy(&mat);
}

TEST_CASE("Matrix: Random Fill", "[matrix]") {
    Matrix mat = MATRIX_INIT;
    Matrix means = MATRIX_INIT;
    Matrix stds = MATRIX_INIT;
    mat_create(&mat, 10, 10);

    SECTION("Uniform random") {
        REQUIRE(mat_fill_random(&mat) == MATRIX_SUCCESS);
        for (size_t i = 0; i < 100; ++i) {
            REQUIRE(mat.data[i] >= 0.0);
            REQUIRE(mat.data[i] <= 1.0);
        }
    }

    SECTION("Gaussian random") {
        mat_create(&means, 10, 1);
        mat_create(&stds, 10, 1);
        mat_fill(&means, 5.0);
        mat_fill(&stds, 2.0);

        REQUIRE(mat_fill_random_gaussian(&mat, &means, &stds) == MATRIX_SUCCESS);
        double sum = 0;
        for (size_t i = 0; i < 100; ++i)
            sum += mat.data[i];
        double average = sum / 100.0;
        REQUIRE_THAT(average, Catch::Matchers::WithinAbs(5.0, 1.5));

        mat_destroy(&means);
        mat_destroy(&stds);
    }

    mat_destroy(&mat);
}
