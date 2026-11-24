#include "tensor.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

// Helper Functions for Testing

// Checks if two TensorF64 instances are identical.
bool are_f64_tensors_equal(const TensorF64 *a, const TensorF64 *b, double epsilon) {
    if (a == NULL || b == NULL)
        return false;
    if (a->ndim != b->ndim)
        return false;
    if (memcmp(a->shape, b->shape, a->ndim * sizeof(unsigned int)) != 0)
        return false;
    for (size_t i = 0; i < a->size; ++i) {
        if (fabs(a->data[i] - b->data[i]) > epsilon)
            return false;
    }
    return true;
}

TEST_CASE("TensorF64: Creation and Destruction", "[tensorf64]") {
    TensorF64 t = TENSOR_INIT;

    SECTION("Creation with valid dimensions") {
        const unsigned int shape[] = {2, 3, 4};
        REQUIRE(tensorf64_create(&t, 3, shape) == TENSOR_SUCCESS);
        REQUIRE(t.ndim == 3);
        REQUIRE(t.size == 24);
        REQUIRE(t.shape[0] == 2);
        REQUIRE(t.shape[1] == 3);
        REQUIRE(t.shape[2] == 4);
        REQUIRE(t.data != NULL);
        for (size_t i = 0; i < 24; ++i)
            REQUIRE(t.data[i] == 0.0);
        tensorf64_destroy(&t);
    }

    SECTION("Creation with invalid dimensions") {
        const unsigned int zero_shape[] = {2, 0, 4};
        const unsigned int valid_shape[] = {2, 3};
        REQUIRE(tensorf64_create(&t, 0, valid_shape) == TENSOR_ERR_INVALID_DIMENSION);
        REQUIRE(tensorf64_create(&t, TENSOR_MAX_DIMS + 1, valid_shape) == TENSOR_ERR_INVALID_DIMENSION);
        REQUIRE(tensorf64_create(&t, 3, zero_shape) == TENSOR_ERR_INVALID_DIMENSION);
    }

    SECTION("Creation with null pointer and resize") {
        const unsigned int shape_a[] = {2, 3};
        const unsigned int shape_b[] = {3, 4};
        REQUIRE(tensorf64_create(NULL, 2, shape_a) == TENSOR_ERR_NULL_PTR);
        REQUIRE(tensorf64_create(&t, 2, NULL) == TENSOR_ERR_NULL_PTR);
        REQUIRE(tensorf64_create(&t, 2, shape_a) == TENSOR_SUCCESS);
        REQUIRE(tensorf64_create(&t, 2, shape_b) == TENSOR_SUCCESS);
        REQUIRE(t.size == 12);
        tensorf64_destroy(&t);
    }

    SECTION("Destruction") {
        const unsigned int shape[] = {2, 2};
        tensorf64_create(&t, 2, shape);
        REQUIRE(tensorf64_destroy(&t) == TENSOR_SUCCESS);
        REQUIRE(t.data == NULL);
        REQUIRE(tensorf64_destroy(NULL) == TENSOR_ERR_NULL_PTR);
        REQUIRE(tensorf64_destroy(&t) != TENSOR_SUCCESS);
    }
}

TEST_CASE("TensorF64: Copy and Reshape", "[tensorf64]") {
    TensorF64          t = TENSOR_INIT;
    TensorF64          copy = TENSOR_INIT;
    const unsigned int shape[] = {2, 3};
    const double       src[6] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    double             tol = 1.0e-9;

    tensorf64_create(&t, 2, shape);
    tensorf64_set_data(&t, src, 6);

    SECTION("Copy valid and self") {
        REQUIRE(tensorf64_copy(&copy, &t) == TENSOR_SUCCESS);
        REQUIRE(are_f64_tensors_equal(&t, &copy, tol));
        REQUIRE(t.data != copy.data);
        REQUIRE(tensorf64_copy(&t, &t) == TENSOR_SUCCESS);
        tensorf64_destroy(&copy);
    }

    SECTION("Copy null pointer") {
        REQUIRE(tensorf64_copy(NULL, &t) == TENSOR_ERR_NULL_PTR);
        REQUIRE(tensorf64_copy(&copy, NULL) == TENSOR_ERR_NULL_PTR);
    }

    SECTION("Reshape valid") {
        const unsigned int new_shape[] = {3, 2};
        REQUIRE(tensorf64_reshape(&t, 2, new_shape) == TENSOR_SUCCESS);
        REQUIRE(t.ndim == 2);
        REQUIRE(t.shape[0] == 3);
        REQUIRE(t.shape[1] == 2);
        REQUIRE_THAT(TENSOR_AT(&t, 1, 0), Catch::Matchers::WithinAbs(3.0, tol));
        REQUIRE_THAT(TENSOR_AT(&t, 2, 1), Catch::Matchers::WithinAbs(6.0, tol));
    }

    SECTION("Reshape error cases") {
        const unsigned int mismatch_shape[] = {4, 2};
        const unsigned int zero_shape[] = {0, 6};
        REQUIRE(tensorf64_reshape(&t, 2, mismatch_shape) == TENSOR_ERR_DIMENSION_MISMATCH);
        REQUIRE(tensorf64_reshape(&t, 0, mismatch_shape) == TENSOR_ERR_INVALID_DIMENSION);
        REQUIRE(tensorf64_reshape(&t, 2, zero_shape) == TENSOR_ERR_INVALID_DIMENSION);
        REQUIRE(tensorf64_reshape(NULL, 2, mismatch_shape) == TENSOR_ERR_NULL_PTR);
    }

    tensorf64_destroy(&t);
}

TEST_CASE("TensorF64: Set Data", "[tensorf64]") {
    TensorF64          t = TENSOR_INIT;
    const unsigned int shape[] = {2, 2, 3};
    const double       src[12] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0};
    tensorf64_create(&t, 3, shape);

    SECTION("Fill with flat array in row-major order") {
        REQUIRE(tensorf64_set_data(&t, src, 12) == TENSOR_SUCCESS);
        for (size_t i = 0; i < 12; ++i)
            REQUIRE(t.data[i] == src[i]);
        REQUIRE(TENSOR_AT(&t, 0, 1, 2) == 6.0);
        REQUIRE(TENSOR_AT(&t, 1, 0, 1) == 8.0);
    }

    SECTION("Overwrite existing tensor data") {
        const double src2[12] = {-1.0, -2.0, -3.0, -4.0, -5.0, -6.0,
                                 -7.0, -8.0, -9.0, -10.0, -11.0, -12.0};
        tensorf64_set_data(&t, src, 12);
        REQUIRE(tensorf64_set_data(&t, src2, 12) == TENSOR_SUCCESS);
        for (size_t i = 0; i < 12; ++i)
            REQUIRE(t.data[i] == src2[i]);
    }

    SECTION("Null pointer handling") {
        TensorF64 uninit = TENSOR_INIT;
        REQUIRE(tensorf64_set_data(NULL, src, 12) == TENSOR_ERR_NULL_PTR);
        REQUIRE(tensorf64_set_data(&uninit, src, 12) == TENSOR_ERR_NULL_PTR);
        REQUIRE(tensorf64_set_data(&t, NULL, 12) == TENSOR_ERR_NULL_PTR);
    }

    SECTION("Element count mismatch") {
        tensorf64_set_data(&t, src, 12);
        REQUIRE(tensorf64_set_data(&t, src, 11) == TENSOR_ERR_DIMENSION_MISMATCH);
        REQUIRE(tensorf64_set_data(&t, src, 0) == TENSOR_ERR_DIMENSION_MISMATCH);
        for (size_t i = 0; i < 12; ++i)
            REQUIRE(t.data[i] == src[i]);
    }

    tensorf64_destroy(&t);
}

TEST_CASE("TensorF64: Elementwise Add and Sub", "[tensorf64]") {
    TensorF64          a = TENSOR_INIT;
    TensorF64          b = TENSOR_INIT;
    TensorF64          expected = TENSOR_INIT;
    const unsigned int shape[] = {2, 2, 2};
    const double       data_a[8] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
    const double       data_b[8] = {8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0};
    double             tol = 1.0e-9;

    tensorf64_create(&a, 3, shape);
    tensorf64_create(&b, 3, shape);
    tensorf64_create(&expected, 3, shape);
    tensorf64_set_data(&a, data_a, 8);
    tensorf64_set_data(&b, data_b, 8);

    SECTION("Addition") {
        const double exp_data[8] = {9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0};
        tensorf64_set_data(&expected, exp_data, 8);
        REQUIRE(tensorf64_add(&a, &b) == TENSOR_SUCCESS);
        REQUIRE(are_f64_tensors_equal(&a, &expected, tol));
    }

    SECTION("Subtraction") {
        const double exp_data[8] = {-7.0, -5.0, -3.0, -1.0, 1.0, 3.0, 5.0, 7.0};
        tensorf64_set_data(&expected, exp_data, 8);
        REQUIRE(tensorf64_sub(&a, &b) == TENSOR_SUCCESS);
        REQUIRE(are_f64_tensors_equal(&a, &expected, tol));
    }

    SECTION("Add and Sub dimension mismatch") {
        TensorF64          wrong = TENSOR_INIT;
        const unsigned int wrong_shape[] = {2, 4};
        tensorf64_create(&wrong, 2, wrong_shape);
        REQUIRE(tensorf64_add(&a, &wrong) == TENSOR_ERR_DIMENSION_MISMATCH);
        REQUIRE(tensorf64_sub(&a, &wrong) == TENSOR_ERR_DIMENSION_MISMATCH);
        tensorf64_destroy(&wrong);
    }

    SECTION("Add and Sub null pointers") {
        REQUIRE(tensorf64_add(NULL, &b) == TENSOR_ERR_NULL_PTR);
        REQUIRE(tensorf64_add(&a, NULL) == TENSOR_ERR_NULL_PTR);
        REQUIRE(tensorf64_sub(NULL, &b) == TENSOR_ERR_NULL_PTR);
        REQUIRE(tensorf64_sub(&a, NULL) == TENSOR_ERR_NULL_PTR);
    }

    tensorf64_destroy(&a);
    tensorf64_destroy(&b);
    tensorf64_destroy(&expected);
}

TEST_CASE("TensorF64: Batched Multiplication", "[tensorf64]") {
    TensorF64          b = TENSOR_INIT;
    TensorF64          out = TENSOR_INIT;
    TensorF64          expected = TENSOR_INIT;
    const unsigned int shape_b[] = {2, 3, 2};
    const unsigned int shape_out[] = {2, 2, 2};
    const double       data_b[12] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 2.0, 0.0, 0.0, 2.0, 1.0, 1.0};
    const double       exp_data[8] = {22.0, 28.0, 49.0, 64.0, 23.0, 25.0, 32.0, 34.0};
    double             tol = 1.0e-9;

    tensorf64_create(&b, 3, shape_b);
    tensorf64_set_data(&b, data_b, 12);
    tensorf64_create(&expected, 3, shape_out);
    tensorf64_set_data(&expected, exp_data, 8);

    SECTION("Batched multiplication without transpose") {
        TensorF64          a = TENSOR_INIT;
        const unsigned int shape_a[] = {2, 2, 3};
        const double       data_a[12] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0};
        tensorf64_create(&a, 3, shape_a);
        tensorf64_set_data(&a, data_a, 12);

        REQUIRE(tensorf64_mul(&a, false, &b, false, &out) == TENSOR_SUCCESS);
        REQUIRE(are_f64_tensors_equal(&out, &expected, tol));
        tensorf64_destroy(&a);
    }

    SECTION("Batched multiplication with transpose") {
        TensorF64          a_t = TENSOR_INIT;
        const double       data_at[12] = {1.0, 4.0, 2.0, 5.0, 3.0, 6.0, 7.0, 10.0, 8.0, 11.0, 9.0, 12.0};
        tensorf64_create(&a_t, 3, shape_b);
        tensorf64_set_data(&a_t, data_at, 12);

        REQUIRE(tensorf64_mul(&a_t, true, &b, false, &out) == TENSOR_SUCCESS);
        REQUIRE(are_f64_tensors_equal(&out, &expected, tol));
        tensorf64_destroy(&a_t);
    }

    SECTION("Dimension and batch mismatch") {
        TensorF64          bad_batch = TENSOR_INIT;
        const unsigned int bad_shape[] = {3, 2, 3};
        tensorf64_create(&bad_batch, 3, bad_shape);
        REQUIRE(tensorf64_mul(&bad_batch, false, &b, false, &out) == TENSOR_ERR_DIMENSION_MISMATCH);
        REQUIRE(tensorf64_mul(&b, false, &b, false, &out) == TENSOR_ERR_DIMENSION_MISMATCH);
        tensorf64_destroy(&bad_batch);
    }

    SECTION("In-place and null pointer errors") {
        REQUIRE(tensorf64_mul(&b, true, &b, false, &b) == TENSOR_ERR_INPLACE_MATMUL);
        REQUIRE(tensorf64_mul(NULL, false, &b, false, &out) == TENSOR_ERR_NULL_PTR);
        REQUIRE(tensorf64_mul(&b, false, &b, false, NULL) == TENSOR_ERR_NULL_PTR);
    }

    tensorf64_destroy(&b);
    tensorf64_destroy(&out);
    tensorf64_destroy(&expected);
}
