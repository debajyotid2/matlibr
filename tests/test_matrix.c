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
bool are_double_matrices_equal(const Matrix *a, const Matrix *b,
                               double epsilon) {
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

TEST_CASE("IntMatrix: Creation and Destruction", "[intmatrix]") {
    IntMatrix mat;

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
    IntMatrix mat_a, result;
    intmat_create(&mat_a, 2, 2);
    mat_a.data[0] = 1;
    mat_a.data[1] = 2;
    mat_a.data[2] = 3;
    mat_a.data[3] = 4;

    SECTION("Copy") {
        REQUIRE(intmat_copy(&result, &mat_a) == MATRIX_SUCCESS);
        REQUIRE(are_int_matrices_equal(&mat_a, &result));
        REQUIRE(mat_a.data != result.data);
        intmat_destroy(&result);
    }

    SECTION("Assign") {
        int row = 0, col = 1;
        int value = 42;
        REQUIRE(intmat_assign(&mat_a, row, col, value) == MATRIX_SUCCESS);
        REQUIRE(mat_a.data[row * mat_a.ncols + col] == value);
    }
    
    SECTION("At") {
        int row = 0, col = 1;
        int value;
        REQUIRE(intmat_at(&mat_a, row, col, &value) == MATRIX_SUCCESS);
        REQUIRE(value == 2);
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

    intmat_destroy(&mat_a);
}

TEST_CASE("IntMatrix: Random Fill", "[intmatrix]") {
    IntMatrix mat;
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
        REQUIRE(intmat_fill_random(&mat, 10, 5, false) ==
                MATRIX_ERR_RANGE_INVALID);
        REQUIRE(intmat_fill_random(&mat, 0, 19, false) ==
                MATRIX_ERR_TOO_MANY_INTS_TO_GENERATE);
    }

    intmat_destroy(&mat);
}

TEST_CASE("IntMatrix: Arithmetic Operations", "[intmatrix]") {
    IntMatrix mat_a, mat_b, result, expected_mat;
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
        REQUIRE(intmat_mul(&mat_a, false, &mat_b, false, &result) ==
                MATRIX_SUCCESS);
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
    IntMatrix mat, vec, result, expected_mat;
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
    Matrix mat;
    mat_create(&mat, 2, 2);
    REQUIRE(mat_fill(&mat, 3.5) == MATRIX_SUCCESS);
    for (int i = 0; i < 4; ++i) {
        REQUIRE_THAT(mat.data[i], Catch::Matchers::WithinAbs(3.5, 1e-9));
    }

    SECTION("Assign") {
        int    row = 0, col = 1;
        double value = 42.0;
        REQUIRE(mat_assign(&mat, row, col, value) == MATRIX_SUCCESS);
        REQUIRE_THAT(mat.data[row * mat.ncols + col],
                     Catch::Matchers::WithinAbs(value, 1e-9));
    }

    SECTION("At") {
        int row = 0, col = 1;
        double value;
        REQUIRE(mat_at(&mat, row, col, &value) == MATRIX_SUCCESS);
        REQUIRE_THAT(value, Catch::Matchers::WithinAbs(3.5, 1e-9));
    }

    mat_destroy(&mat);
}

TEST_CASE("Matrix: Scalar and Arithmetic Operations", "[matrix]") {
    Matrix mat;
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
        Matrix mat_b;
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
    Matrix    mat, vec_row, vec_col, result, expected_mat;
    IntMatrix indices;
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
        REQUIRE(mat_gather(&mat, &result, &indices, 1) ==
                MATRIX_ERR_INDEX_OUT_OF_BOUNDS);
        intmat_destroy(&indices);
        mat_destroy(&result);
        mat_destroy(&expected_mat);
    }

    mat_destroy(&mat);
    mat_destroy(&vec_row);
    mat_destroy(&vec_col);
}

TEST_CASE("Matrix: BLAS-based Functions", "[matrix]") {
    Matrix mat;
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
        REQUIRE_THAT(norm,
                     Catch::Matchers::WithinAbs(sqrt(expected_norm_sq), 1e-9));
    }

    mat_destroy(&mat);
}

TEST_CASE("Matrix: Random Fill", "[matrix]") {
    Matrix mat, means, stds;
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

        REQUIRE(mat_fill_random_gaussian(&mat, &means, &stds) ==
                MATRIX_SUCCESS);
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
