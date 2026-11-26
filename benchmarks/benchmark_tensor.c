#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/catch_test_macros.hpp>
#include <stdlib.h>

#include "matrix.h"
#include "tensor.h"

// Benchmarks for TensorF64 (Contiguous Arithmetic & Batched Matmul)
TEST_CASE("TensorF64 Performance", "[tensorf64][benchmark]") {
    const unsigned int shapes[3][3] = {
        {4, 32, 32},
        {4, 256, 256},
        {4, 1024, 1024}
    };
    const char *labels[] = {"4x32x32", "4x256x256", "4x1024x1024"};

    for (size_t idx = 0; idx < 3; ++idx) {
        TensorF64 a = TENSOR_INIT;
        TensorF64 b = TENSOR_INIT;
        TensorF64 r = TENSOR_INIT;
        tensorf64_create(&a, 3, shapes[idx]);
        tensorf64_create(&b, 3, shapes[idx]);
        tensorf64_create(&r, 3, shapes[idx]);
        for (size_t i = 0; i < a.size; ++i) {
            a.data[i] = 1.0;
            b.data[i] = 2.0;
        }

        BENCHMARK(std::string("TensorF64 (BLAS) Add [") + labels[idx] + "]") {
            return tensorf64_add(&a, &b);
        };

        BENCHMARK(std::string("TensorF64 (BLAS) Sub [") + labels[idx] + "]") {
            return tensorf64_sub(&a, &b);
        };

        BENCHMARK(std::string("TensorF64 (BLAS) Mul [") + labels[idx] + "]") {
            return tensorf64_mul(&a, false, &b, false, &r);
        };

        tensorf64_destroy(&a);
        tensorf64_destroy(&b);
        tensorf64_destroy(&r);
    }
}

// Benchmarks for Permuted / Strided Views (Coalescing + Odometer + Zero-Copy GEMM)
TEST_CASE("TensorF64 Strided View Performance", "[tensorf64][strided][benchmark]") {
    const unsigned int shape[] = {8, 256, 256};
    const unsigned int perm_batch_shape[] = {256, 8, 256};
    const unsigned int perm_batch[] = {1, 0, 2};
    const unsigned int perm_trans[] = {0, 2, 1};

    TensorF64 a = TENSOR_INIT;
    TensorF64 b = TENSOR_INIT;
    TensorF64 out = TENSOR_INIT;
    TensorF64 perm_a = TENSOR_INIT;
    TensorF64 perm_b = TENSOR_INIT;
    TensorF64 trans_a = TENSOR_INIT;
    TensorF64 copy_dst = TENSOR_INIT;

    tensorf64_create(&a, 3, shape);
    tensorf64_create(&b, 3, shape);
    tensorf64_create(&out, 3, shape);
    tensorf64_create(&copy_dst, 3, perm_batch_shape);
    for (size_t i = 0; i < a.size; ++i) {
        a.data[i] = 1.0;
        b.data[i] = 2.0;
    }

    tensorf64_permute(&a, &perm_a, perm_batch);
    tensorf64_permute(&b, &perm_b, perm_batch);
    tensorf64_permute(&a, &trans_a, perm_trans);

    BENCHMARK("TensorF64 Permute View [8x256x256]") {
        return tensorf64_permute(&a, &perm_a, perm_batch);
    };

    BENCHMARK("TensorF64 Strided Copy (Materialize) [256x8x256]") {
        return tensorf64_copy(&copy_dst, &perm_a);
    };

    BENCHMARK("TensorF64 Strided Add [256x8x256]") {
        return tensorf64_add(&perm_a, &perm_b);
    };

    BENCHMARK("TensorF64 Strided Batched Mul (Transposed View) [8x256x256]") {
        return tensorf64_mul(&trans_a, false, &b, false, &out);
    };

    tensorf64_destroy(&perm_a);
    tensorf64_destroy(&perm_b);
    tensorf64_destroy(&trans_a);
    tensorf64_destroy(&copy_dst);
    tensorf64_destroy(&a);
    tensorf64_destroy(&b);
    tensorf64_destroy(&out);
}

// Benchmarks for Creation, Copy, Reshape, Set Data, and Gather
TEST_CASE("TensorF64 Memory, Metadata, and Gather Performance", "[tensorf64][benchmark]") {
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
    tensorf64_create(&f64_dst, 3, shape);

    BENCHMARK("TensorF64 Copy [4x512x512]") {
        return tensorf64_copy(&f64_dst, &f64_src);
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

    MatrixI32 idxs = MATRIX_INIT;
    TensorF64 gathered = TENSOR_INIT;
    mati32_create(&idxs, 256, 1);
    for (int i = 0; i < 256; ++i) {
        idxs.data[i] = (i * 2) % 512;
    }

    BENCHMARK("TensorF64 Gather Axis 1 [4x512x512 -> 4x256x512]") {
        return tensorf64_gather(&f64_src, &gathered, &idxs, 1);
    };

    TensorF64 perm_src = TENSOR_INIT;
    TensorF64 gathered_strided = TENSOR_INIT;
    const unsigned int perm_dims[] = {1, 0, 2};
    tensorf64_permute(&f64_src, &perm_src, perm_dims);

    BENCHMARK("TensorF64 Strided Gather Axis 0 [512x4x512 -> 256x4x512]") {
        return tensorf64_gather(&perm_src, &gathered_strided, &idxs, 0);
    };

    tensorf64_destroy(&perm_src);
    tensorf64_destroy(&gathered);
    tensorf64_destroy(&gathered_strided);
    mati32_destroy(&idxs);
    tensorf64_destroy(&f64_src);
    tensorf64_destroy(&f64_dst);
}
