#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/catch_test_macros.hpp>
#include <stdlib.h>

#include "tensor.h"

// Benchmarks for TensorF64 (Double Precision Operations with BLAS)
TEST_CASE("TensorF64 Performance", "[tensorf64][benchmark]") {
    // Small Tensor Benchmarks (4x32x32)
    const unsigned int small_shape[] = {4, 32, 32};
    TensorF64          sm_a = TENSOR_INIT;
    TensorF64          sm_b = TENSOR_INIT;
    TensorF64          sm_r = TENSOR_INIT;
    tensorf64_create(&sm_a, 3, small_shape);
    tensorf64_create(&sm_b, 3, small_shape);
    tensorf64_create(&sm_r, 3, small_shape);
    for (size_t i = 0; i < sm_a.size; ++i) {
        sm_a.data[i] = 1.0;
        sm_b.data[i] = 2.0;
    }

    BENCHMARK("TensorF64 (BLAS) Add [4x32x32]") { return tensorf64_add(&sm_a, &sm_b); };

    BENCHMARK("TensorF64 (BLAS) Sub [4x32x32]") { return tensorf64_sub(&sm_a, &sm_b); };

    BENCHMARK("TensorF64 (BLAS) Mul [4x32x32]") {
        return tensorf64_mul(&sm_a, false, &sm_b, false, &sm_r);
    };

    tensorf64_destroy(&sm_a);
    tensorf64_destroy(&sm_b);
    tensorf64_destroy(&sm_r);

    // Medium Tensor Benchmarks (4x256x256)
    const unsigned int med_shape[] = {4, 256, 256};
    TensorF64          mm_a = TENSOR_INIT;
    TensorF64          mm_b = TENSOR_INIT;
    TensorF64          mm_r = TENSOR_INIT;
    tensorf64_create(&mm_a, 3, med_shape);
    tensorf64_create(&mm_b, 3, med_shape);
    tensorf64_create(&mm_r, 3, med_shape);
    for (size_t i = 0; i < mm_a.size; ++i) {
        mm_a.data[i] = 1.0;
        mm_b.data[i] = 2.0;
    }

    BENCHMARK("TensorF64 (BLAS) Add [4x256x256]") { return tensorf64_add(&mm_a, &mm_b); };

    BENCHMARK("TensorF64 (BLAS) Sub [4x256x256]") { return tensorf64_sub(&mm_a, &mm_b); };

    BENCHMARK("TensorF64 (BLAS) Mul [4x256x256]") {
        return tensorf64_mul(&mm_a, false, &mm_b, false, &mm_r);
    };

    tensorf64_destroy(&mm_a);
    tensorf64_destroy(&mm_b);
    tensorf64_destroy(&mm_r);

    // Large Tensor Benchmarks (4x1024x1024)
    const unsigned int large_shape[] = {4, 1024, 1024};
    TensorF64          lm_a = TENSOR_INIT;
    TensorF64          lm_b = TENSOR_INIT;
    TensorF64          lm_r = TENSOR_INIT;
    tensorf64_create(&lm_a, 3, large_shape);
    tensorf64_create(&lm_b, 3, large_shape);
    tensorf64_create(&lm_r, 3, large_shape);
    for (size_t i = 0; i < lm_a.size; ++i) {
        lm_a.data[i] = 1.0;
        lm_b.data[i] = 2.0;
    }

    BENCHMARK("TensorF64 (BLAS) Add [4x1024x1024]") { return tensorf64_add(&lm_a, &lm_b); };

    BENCHMARK("TensorF64 (BLAS) Sub [4x1024x1024]") { return tensorf64_sub(&lm_a, &lm_b); };

    BENCHMARK("TensorF64 (BLAS) Mul [4x1024x1024]") {
        return tensorf64_mul(&lm_a, false, &lm_b, false, &lm_r);
    };

    tensorf64_destroy(&lm_a);
    tensorf64_destroy(&lm_b);
    tensorf64_destroy(&lm_r);
}

TEST_CASE("TensorF64 Creation, Copy, and Metadata Performance", "[tensorf64][benchmark]") {
    const unsigned int shape[] = {4, 512, 512};
    const unsigned int reshaped[] = {16, 256, 256};

    BENCHMARK("TensorF64 Creation [4x512x512]") {
        TensorF64 t = TENSOR_INIT;
        tensorf64_create(&t, 3, shape);
        tensorf64_destroy(&t);
        return 0;
    };

    TensorF64 f64_src = TENSOR_INIT;
    TensorF64 f64_dst = TENSOR_INIT;
    tensorf64_create(&f64_src, 3, shape);

    BENCHMARK("TensorF64 Copy [4x512x512]") {
        tensorf64_copy(&f64_dst, &f64_src);
        tensorf64_destroy(&f64_dst);
        return 0;
    };

    double *raw_src = (double *)malloc(f64_src.size * sizeof(double));
    BENCHMARK("TensorF64 Set Data [4x512x512]") {
        return tensorf64_set_data(&f64_src, raw_src, f64_src.size);
    };
    free(raw_src);

    BENCHMARK("TensorF64 Reshape [4x512x512 <-> 16x256x256]") {
        tensorf64_reshape(&f64_src, 3, reshaped);
        return tensorf64_reshape(&f64_src, 3, shape);
    };

    tensorf64_destroy(&f64_src);
}
