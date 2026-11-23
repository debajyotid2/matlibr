/* Benchmarks for module matrix.h

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

#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/catch_test_macros.hpp>

#include "matrix.h"

// Benchmarks for IntMatrix (Integer Operations)
TEST_CASE("IntMatrix Performance", "[intmatrix][benchmark]") {
    // Small Matrix Benchmarks (32x32)
    const int small_dim = 32;
    IntMatrix sm_a = MATRIX_INIT;
    IntMatrix sm_b = MATRIX_INIT;
    IntMatrix sm_r = MATRIX_INIT;
    intmat_create(&sm_a, small_dim, small_dim);
    intmat_create(&sm_b, small_dim, small_dim);
    intmat_create(&sm_r, small_dim, small_dim);
    intmat_fill(&sm_a, 1);
    intmat_fill(&sm_b, 2);

    BENCHMARK("IntMatrix Add [32x32]") { return intmat_add(&sm_a, &sm_b); };

    BENCHMARK("IntMatrix Mul [32x32]") { return intmat_mul(&sm_a, false, &sm_b, false, &sm_r); };

    intmat_destroy(&sm_a);
    intmat_destroy(&sm_b);
    intmat_destroy(&sm_r);

    // Medium Matrix Benchmarks (256x256)
    const int med_dim = 256;
    IntMatrix mm_a = MATRIX_INIT;
    IntMatrix mm_b = MATRIX_INIT;
    IntMatrix mm_r = MATRIX_INIT;
    intmat_create(&mm_a, med_dim, med_dim);
    intmat_create(&mm_b, med_dim, med_dim);
    intmat_create(&mm_r, med_dim, med_dim);
    intmat_fill(&mm_a, 1);
    intmat_fill(&mm_b, 2);

    BENCHMARK("IntMatrix Add [256x256]") { return intmat_add(&mm_a, &mm_b); };

    BENCHMARK("IntMatrix Mul [256x256]") { return intmat_mul(&mm_a, false, &mm_b, false, &mm_r); };

    intmat_destroy(&mm_a);
    intmat_destroy(&mm_b);
    intmat_destroy(&mm_r);

    // Large Matrix Benchmarks (1024x1024)
    const int large_dim = 1024;
    IntMatrix ll_a = MATRIX_INIT;
    IntMatrix ll_b = MATRIX_INIT;
    IntMatrix ll_r = MATRIX_INIT;
    intmat_create(&ll_a, large_dim, large_dim);
    intmat_create(&ll_b, large_dim, large_dim);
    intmat_create(&ll_r, large_dim, large_dim);
    intmat_fill(&ll_a, 1);
    intmat_fill(&ll_b, 2);

    BENCHMARK("IntMatrix Add [1024x1024]") { return intmat_add(&ll_a, &ll_b); };

    BENCHMARK("IntMatrix Mul [1024x1024]") {
        return intmat_mul(&ll_a, false, &ll_b, false, &ll_r);
    };

    intmat_destroy(&ll_a);
    intmat_destroy(&ll_b);
    intmat_destroy(&ll_r);
}

// Benchmarks for Matrix (Double Precision Operations with BLAS)
TEST_CASE("Matrix Performance", "[matrix][benchmark]") {
    // Small Matrix Benchmarks (32x32)
    const int small_dim = 32;
    Matrix    sm_a = MATRIX_INIT;
    Matrix    sm_b = MATRIX_INIT;
    Matrix    sm_r = MATRIX_INIT;
    mat_create(&sm_a, small_dim, small_dim);
    mat_create(&sm_b, small_dim, small_dim);
    mat_create(&sm_r, small_dim, small_dim);
    mat_fill(&sm_a, 1.0);
    mat_fill(&sm_b, 2.0);

    BENCHMARK("Matrix (BLAS) Add [32x32]") { return mat_add(&sm_a, &sm_b); };

    BENCHMARK("Matrix (BLAS) Mul [32x32]") { return mat_mul(&sm_a, false, &sm_b, false, &sm_r); };

    mat_destroy(&sm_a);
    mat_destroy(&sm_b);
    mat_destroy(&sm_r);

    // Medium Matrix Benchmarks (256x256)
    const int med_dim = 256;
    Matrix    mm_a = MATRIX_INIT;
    Matrix    mm_b = MATRIX_INIT;
    Matrix    mm_r = MATRIX_INIT;
    mat_create(&mm_a, med_dim, med_dim);
    mat_create(&mm_b, med_dim, med_dim);
    mat_create(&mm_r, med_dim, med_dim);
    mat_fill(&mm_a, 1.0);
    mat_fill(&mm_b, 2.0);

    BENCHMARK("Matrix (BLAS) Add [256x256]") { return mat_add(&mm_a, &mm_b); };

    BENCHMARK("Matrix (BLAS) Mul [256x256]") { return mat_mul(&mm_a, false, &mm_b, false, &mm_r); };

    mat_destroy(&mm_a);
    mat_destroy(&mm_b);
    mat_destroy(&mm_r);

    // Large Matrix Benchmarks (1024x1024)
    const int large_dim = 1024;
    Matrix    lm_a = MATRIX_INIT;
    Matrix    lm_b = MATRIX_INIT;
    Matrix    lm_r = MATRIX_INIT;
    mat_create(&lm_a, large_dim, large_dim);
    mat_create(&lm_b, large_dim, large_dim);
    mat_create(&lm_r, large_dim, large_dim);
    mat_fill(&lm_a, 1.0);
    mat_fill(&lm_b, 2.0);

    BENCHMARK("Matrix (BLAS) Add [1024x1024]") { return mat_add(&lm_a, &lm_b); };

    BENCHMARK("Matrix (BLAS) Mul [1024x1024]") {
        return mat_mul(&lm_a, false, &lm_b, false, &lm_r);
    };

    mat_destroy(&lm_a);
    mat_destroy(&lm_b);
    mat_destroy(&lm_r);
}

TEST_CASE("Matrix Creation and Copy Performance", "[matrix][intmatrix][benchmark]") {
    const int dim = 512;

    BENCHMARK("IntMatrix Creation [512x512]") {
        IntMatrix mat = MATRIX_INIT;
        intmat_create(&mat, dim, dim);
        intmat_destroy(&mat);
        return 0;
    };

    BENCHMARK("Matrix Creation [512x512]") {
        Matrix mat = MATRIX_INIT;
        mat_create(&mat, dim, dim);
        mat_destroy(&mat);
        return 0;
    };

    IntMatrix int_src = MATRIX_INIT;
    IntMatrix int_dst = MATRIX_INIT;
    intmat_create(&int_src, dim, dim);
    BENCHMARK("IntMatrix Copy [512x512]") {
        intmat_copy(&int_dst, &int_src);
        intmat_destroy(&int_dst);
        return 0;
    };
    intmat_destroy(&int_src);

    Matrix double_src = MATRIX_INIT;
    Matrix double_dst = MATRIX_INIT;
    mat_create(&double_src, dim, dim);
    BENCHMARK("Matrix Copy [512x512]") {
        mat_copy(&double_dst, &double_src);
        mat_destroy(&double_dst);
        return 0;
    };
    mat_destroy(&double_src);
}

// Benchmarks for MatrixI32 (32-bit Integer Operations)
TEST_CASE("MatrixI32 Performance", "[matrixi32][benchmark]") {
    // Small Matrix Benchmarks (32x32)
    const int small_dim = 32;
    MatrixI32 sm_a = MATRIX_INIT;
    MatrixI32 sm_b = MATRIX_INIT;
    MatrixI32 sm_r = MATRIX_INIT;
    mati32_create(&sm_a, small_dim, small_dim);
    mati32_create(&sm_b, small_dim, small_dim);
    mati32_create(&sm_r, small_dim, small_dim);
    mati32_fill(&sm_a, 1);
    mati32_fill(&sm_b, 2);

    BENCHMARK("MatrixI32 Add [32x32]") { return mati32_add(&sm_a, &sm_b); };

    BENCHMARK("MatrixI32 Mul [32x32]") { return mati32_mul(&sm_a, false, &sm_b, false, &sm_r); };

    mati32_destroy(&sm_a);
    mati32_destroy(&sm_b);
    mati32_destroy(&sm_r);

    // Medium Matrix Benchmarks (256x256)
    const int med_dim = 256;
    MatrixI32 mm_a = MATRIX_INIT;
    MatrixI32 mm_b = MATRIX_INIT;
    MatrixI32 mm_r = MATRIX_INIT;
    mati32_create(&mm_a, med_dim, med_dim);
    mati32_create(&mm_b, med_dim, med_dim);
    mati32_create(&mm_r, med_dim, med_dim);
    mati32_fill(&mm_a, 1);
    mati32_fill(&mm_b, 2);

    BENCHMARK("MatrixI32 Add [256x256]") { return mati32_add(&mm_a, &mm_b); };

    BENCHMARK("MatrixI32 Mul [256x256]") { return mati32_mul(&mm_a, false, &mm_b, false, &mm_r); };

    mati32_destroy(&mm_a);
    mati32_destroy(&mm_b);
    mati32_destroy(&mm_r);

    // Large Matrix Benchmarks (1024x1024)
    const int large_dim = 1024;
    MatrixI32 lm_a = MATRIX_INIT;
    MatrixI32 lm_b = MATRIX_INIT;
    MatrixI32 lm_r = MATRIX_INIT;
    mati32_create(&lm_a, large_dim, large_dim);
    mati32_create(&lm_b, large_dim, large_dim);
    mati32_create(&lm_r, large_dim, large_dim);
    mati32_fill(&lm_a, 1);
    mati32_fill(&lm_b, 2);

    BENCHMARK("MatrixI32 Add [1024x1024]") { return mati32_add(&lm_a, &lm_b); };

    BENCHMARK("MatrixI32 Mul [1024x1024]") {
        return mati32_mul(&lm_a, false, &lm_b, false, &lm_r);
    };

    mati32_destroy(&lm_a);
    mati32_destroy(&lm_b);
    mati32_destroy(&lm_r);
}

// Benchmarks for MatrixI64 (64-bit Integer Operations)
TEST_CASE("MatrixI64 Performance", "[matrixi64][benchmark]") {
    // Small Matrix Benchmarks (32x32)
    const int small_dim = 32;
    MatrixI64 sm_a = MATRIX_INIT;
    MatrixI64 sm_b = MATRIX_INIT;
    MatrixI64 sm_r = MATRIX_INIT;
    mati64_create(&sm_a, small_dim, small_dim);
    mati64_create(&sm_b, small_dim, small_dim);
    mati64_create(&sm_r, small_dim, small_dim);
    mati64_fill(&sm_a, 1);
    mati64_fill(&sm_b, 2);

    BENCHMARK("MatrixI64 Add [32x32]") { return mati64_add(&sm_a, &sm_b); };

    BENCHMARK("MatrixI64 Mul [32x32]") { return mati64_mul(&sm_a, false, &sm_b, false, &sm_r); };

    mati64_destroy(&sm_a);
    mati64_destroy(&sm_b);
    mati64_destroy(&sm_r);

    // Medium Matrix Benchmarks (256x256)
    const int med_dim = 256;
    MatrixI64 mm_a = MATRIX_INIT;
    MatrixI64 mm_b = MATRIX_INIT;
    MatrixI64 mm_r = MATRIX_INIT;
    mati64_create(&mm_a, med_dim, med_dim);
    mati64_create(&mm_b, med_dim, med_dim);
    mati64_create(&mm_r, med_dim, med_dim);
    mati64_fill(&mm_a, 1);
    mati64_fill(&mm_b, 2);

    BENCHMARK("MatrixI64 Add [256x256]") { return mati64_add(&mm_a, &mm_b); };

    BENCHMARK("MatrixI64 Mul [256x256]") { return mati64_mul(&mm_a, false, &mm_b, false, &mm_r); };

    mati64_destroy(&mm_a);
    mati64_destroy(&mm_b);
    mati64_destroy(&mm_r);

    // Large Matrix Benchmarks (1024x1024)
    const int large_dim = 1024;
    MatrixI64 lm_a = MATRIX_INIT;
    MatrixI64 lm_b = MATRIX_INIT;
    MatrixI64 lm_r = MATRIX_INIT;
    mati64_create(&lm_a, large_dim, large_dim);
    mati64_create(&lm_b, large_dim, large_dim);
    mati64_create(&lm_r, large_dim, large_dim);
    mati64_fill(&lm_a, 1);
    mati64_fill(&lm_b, 2);

    BENCHMARK("MatrixI64 Add [1024x1024]") { return mati64_add(&lm_a, &lm_b); };

    BENCHMARK("MatrixI64 Mul [1024x1024]") {
        return mati64_mul(&lm_a, false, &lm_b, false, &lm_r);
    };

    mati64_destroy(&lm_a);
    mati64_destroy(&lm_b);
    mati64_destroy(&lm_r);
}

// Benchmarks for MatrixF32 (Single Precision Operations with BLAS)
TEST_CASE("MatrixF32 Performance", "[matrixf32][benchmark]") {
    // Small Matrix Benchmarks (32x32)
    const int small_dim = 32;
    MatrixF32 sm_a = MATRIX_INIT;
    MatrixF32 sm_b = MATRIX_INIT;
    MatrixF32 sm_r = MATRIX_INIT;
    matf32_create(&sm_a, small_dim, small_dim);
    matf32_create(&sm_b, small_dim, small_dim);
    matf32_create(&sm_r, small_dim, small_dim);
    matf32_fill(&sm_a, 1.0f);
    matf32_fill(&sm_b, 2.0f);

    BENCHMARK("MatrixF32 (BLAS) Add [32x32]") { return matf32_add(&sm_a, &sm_b); };

    BENCHMARK("MatrixF32 (BLAS) Mul [32x32]") {
        return matf32_mul(&sm_a, false, &sm_b, false, &sm_r);
    };

    matf32_destroy(&sm_a);
    matf32_destroy(&sm_b);
    matf32_destroy(&sm_r);

    // Medium Matrix Benchmarks (256x256)
    const int med_dim = 256;
    MatrixF32 mm_a = MATRIX_INIT;
    MatrixF32 mm_b = MATRIX_INIT;
    MatrixF32 mm_r = MATRIX_INIT;
    matf32_create(&mm_a, med_dim, med_dim);
    matf32_create(&mm_b, med_dim, med_dim);
    matf32_create(&mm_r, med_dim, med_dim);
    matf32_fill(&mm_a, 1.0f);
    matf32_fill(&mm_b, 2.0f);

    BENCHMARK("MatrixF32 (BLAS) Add [256x256]") { return matf32_add(&mm_a, &mm_b); };

    BENCHMARK("MatrixF32 (BLAS) Mul [256x256]") {
        return matf32_mul(&mm_a, false, &mm_b, false, &mm_r);
    };

    matf32_destroy(&mm_a);
    matf32_destroy(&mm_b);
    matf32_destroy(&mm_r);

    // Large Matrix Benchmarks (1024x1024)
    const int large_dim = 1024;
    MatrixF32 lm_a = MATRIX_INIT;
    MatrixF32 lm_b = MATRIX_INIT;
    MatrixF32 lm_r = MATRIX_INIT;
    matf32_create(&lm_a, large_dim, large_dim);
    matf32_create(&lm_b, large_dim, large_dim);
    matf32_create(&lm_r, large_dim, large_dim);
    matf32_fill(&lm_a, 1.0f);
    matf32_fill(&lm_b, 2.0f);

    BENCHMARK("MatrixF32 (BLAS) Add [1024x1024]") { return matf32_add(&lm_a, &lm_b); };

    BENCHMARK("MatrixF32 (BLAS) Mul [1024x1024]") {
        return matf32_mul(&lm_a, false, &lm_b, false, &lm_r);
    };

    matf32_destroy(&lm_a);
    matf32_destroy(&lm_b);
    matf32_destroy(&lm_r);
}

// Benchmarks for MatrixF64 (Double Precision Operations with BLAS)
TEST_CASE("MatrixF64 Performance", "[matrixf64][benchmark]") {
    // Small Matrix Benchmarks (32x32)
    const int small_dim = 32;
    MatrixF64 sm_a = MATRIX_INIT;
    MatrixF64 sm_b = MATRIX_INIT;
    MatrixF64 sm_r = MATRIX_INIT;
    matf64_create(&sm_a, small_dim, small_dim);
    matf64_create(&sm_b, small_dim, small_dim);
    matf64_create(&sm_r, small_dim, small_dim);
    matf64_fill(&sm_a, 1.0);
    matf64_fill(&sm_b, 2.0);

    BENCHMARK("MatrixF64 (BLAS) Add [32x32]") { return matf64_add(&sm_a, &sm_b); };

    BENCHMARK("MatrixF64 (BLAS) Mul [32x32]") {
        return matf64_mul(&sm_a, false, &sm_b, false, &sm_r);
    };

    matf64_destroy(&sm_a);
    matf64_destroy(&sm_b);
    matf64_destroy(&sm_r);

    // Medium Matrix Benchmarks (256x256)
    const int med_dim = 256;
    MatrixF64 mm_a = MATRIX_INIT;
    MatrixF64 mm_b = MATRIX_INIT;
    MatrixF64 mm_r = MATRIX_INIT;
    matf64_create(&mm_a, med_dim, med_dim);
    matf64_create(&mm_b, med_dim, med_dim);
    matf64_create(&mm_r, med_dim, med_dim);
    matf64_fill(&mm_a, 1.0);
    matf64_fill(&mm_b, 2.0);

    BENCHMARK("MatrixF64 (BLAS) Add [256x256]") { return matf64_add(&mm_a, &mm_b); };

    BENCHMARK("MatrixF64 (BLAS) Mul [256x256]") {
        return matf64_mul(&mm_a, false, &mm_b, false, &mm_r);
    };

    matf64_destroy(&mm_a);
    matf64_destroy(&mm_b);
    matf64_destroy(&mm_r);

    // Large Matrix Benchmarks (1024x1024)
    const int large_dim = 1024;
    MatrixF64 lm_a = MATRIX_INIT;
    MatrixF64 lm_b = MATRIX_INIT;
    MatrixF64 lm_r = MATRIX_INIT;
    matf64_create(&lm_a, large_dim, large_dim);
    matf64_create(&lm_b, large_dim, large_dim);
    matf64_create(&lm_r, large_dim, large_dim);
    matf64_fill(&lm_a, 1.0);
    matf64_fill(&lm_b, 2.0);

    BENCHMARK("MatrixF64 (BLAS) Add [1024x1024]") { return matf64_add(&lm_a, &lm_b); };

    BENCHMARK("MatrixF64 (BLAS) Mul [1024x1024]") {
        return matf64_mul(&lm_a, false, &lm_b, false, &lm_r);
    };

    matf64_destroy(&lm_a);
    matf64_destroy(&lm_b);
    matf64_destroy(&lm_r);
}

TEST_CASE("Matrix Creation and Copy Performance", "[matrixi32][matrixi64][benchmark]") {
    const int dim = 512;

    BENCHMARK("MatrixI32 Creation [512x512]") {
        MatrixI32 mat = MATRIX_INIT;
        mati32_create(&mat, dim, dim);
        mati32_destroy(&mat);
        return 0;
    };

    BENCHMARK("MatrixI64 Creation [512x512]") {
        MatrixI64 mat = MATRIX_INIT;
        mati64_create(&mat, dim, dim);
        mati64_destroy(&mat);
        return 0;
    };

    MatrixI32 i32_src = MATRIX_INIT;
    MatrixI32 i32_dst = MATRIX_INIT;
    mati32_create(&i32_src, dim, dim);
    BENCHMARK("MatrixI32 Copy [512x512]") {
        mati32_copy(&i32_dst, &i32_src);
        mati32_destroy(&i32_dst);
        return 0;
    };
    mati32_destroy(&i32_src);

    MatrixI64 i64_src = MATRIX_INIT;
    MatrixI64 i64_dst = MATRIX_INIT;
    mati64_create(&i64_src, dim, dim);
    BENCHMARK("MatrixI64 Copy [512x512]") {
        mati64_copy(&i64_dst, &i64_src);
        mati64_destroy(&i64_dst);
        return 0;
    };
    mati64_destroy(&i64_src);
}

TEST_CASE("Matrix Creation and Copy Performance", "[matrixf32][matrixf64][benchmark]") {
    const int dim = 512;

    BENCHMARK("MatrixF32 Creation [512x512]") {
        MatrixF32 mat = MATRIX_INIT;
        matf32_create(&mat, dim, dim);
        matf32_destroy(&mat);
        return 0;
    };

    BENCHMARK("MatrixF64 Creation [512x512]") {
        MatrixF64 mat = MATRIX_INIT;
        matf64_create(&mat, dim, dim);
        matf64_destroy(&mat);
        return 0;
    };

    MatrixF32 f32_src = MATRIX_INIT;
    MatrixF32 f32_dst = MATRIX_INIT;
    matf32_create(&f32_src, dim, dim);
    BENCHMARK("MatrixF32 Copy [512x512]") {
        matf32_copy(&f32_dst, &f32_src);
        matf32_destroy(&f32_dst);
        return 0;
    };
    matf32_destroy(&f32_src);

    MatrixF64 f64_src = MATRIX_INIT;
    MatrixF64 f64_dst = MATRIX_INIT;
    matf64_create(&f64_src, dim, dim);
    BENCHMARK("MatrixF64 Copy [512x512]") {
        matf64_copy(&f64_dst, &f64_src);
        matf64_destroy(&f64_dst);
        return 0;
    };
    matf64_destroy(&f64_src);
}
