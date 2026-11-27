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
#include <stdlib.h>

#include "matrix.h"

// Benchmarks for MatrixI32 (32-bit Integer Operations)
TEST_CASE("MatrixI32 Performance", "[matrixi32][benchmark]") {
    const int dims[] = {32, 256, 1024};
    const char *labels[] = {"32x32", "256x256", "1024x1024"};

    for (size_t idx = 0; idx < 3; ++idx) {
        const int dim = dims[idx];
        MatrixI32 a = MATRIX_INIT;
        MatrixI32 b = MATRIX_INIT;
        MatrixI32 r = MATRIX_INIT;
        mati32_create(&a, dim, dim);
        mati32_create(&b, dim, dim);
        mati32_create(&r, dim, dim);
        mati32_fill(&a, 1);
        mati32_fill(&b, 2);

        BENCHMARK(std::string("MatrixI32 Add [") + labels[idx] + "]") {
            return mati32_add(&a, &b);
        };

        BENCHMARK(std::string("MatrixI32 Mul [") + labels[idx] + "]") {
            return mati32_mul(&a, false, &b, false, &r);
        };

        mati32_destroy(&a);
        mati32_destroy(&b);
        mati32_destroy(&r);
    }
}

// Benchmarks for MatrixI64 (64-bit Integer Operations)
TEST_CASE("MatrixI64 Performance", "[matrixi64][benchmark]") {
    const int dims[] = {32, 256, 1024};
    const char *labels[] = {"32x32", "256x256", "1024x1024"};

    for (size_t idx = 0; idx < 3; ++idx) {
        const int dim = dims[idx];
        MatrixI64 a = MATRIX_INIT;
        MatrixI64 b = MATRIX_INIT;
        MatrixI64 r = MATRIX_INIT;
        mati64_create(&a, dim, dim);
        mati64_create(&b, dim, dim);
        mati64_create(&r, dim, dim);
        mati64_fill(&a, 1);
        mati64_fill(&b, 2);

        BENCHMARK(std::string("MatrixI64 Add [") + labels[idx] + "]") {
            return mati64_add(&a, &b);
        };

        BENCHMARK(std::string("MatrixI64 Mul [") + labels[idx] + "]") {
            return mati64_mul(&a, false, &b, false, &r);
        };

        mati64_destroy(&a);
        mati64_destroy(&b);
        mati64_destroy(&r);
    }
}

// Benchmarks for MatrixF32 (Single Precision Operations with BLAS)
TEST_CASE("MatrixF32 Performance", "[matrixf32][benchmark]") {
    const int dims[] = {32, 256, 1024};
    const char *labels[] = {"32x32", "256x256", "1024x1024"};

    for (size_t idx = 0; idx < 3; ++idx) {
        const int dim = dims[idx];
        MatrixF32 a = MATRIX_INIT;
        MatrixF32 b = MATRIX_INIT;
        MatrixF32 r = MATRIX_INIT;
        matf32_create(&a, dim, dim);
        matf32_create(&b, dim, dim);
        matf32_create(&r, dim, dim);
        matf32_fill(&a, 1.0f);
        matf32_fill(&b, 2.0f);

        BENCHMARK(std::string("MatrixF32 (BLAS) Add [") + labels[idx] + "]") {
            return matf32_add(&a, &b);
        };

        BENCHMARK(std::string("MatrixF32 (BLAS) Mul [") + labels[idx] + "]") {
            return matf32_mul(&a, false, &b, false, &r);
        };

        matf32_destroy(&a);
        matf32_destroy(&b);
        matf32_destroy(&r);
    }
}

// Benchmarks for MatrixF64 (Double Precision Operations with BLAS)
TEST_CASE("MatrixF64 Performance", "[matrixf64][benchmark]") {
    const int dims[] = {32, 256, 1024};
    const char *labels[] = {"32x32", "256x256", "1024x1024"};

    for (size_t idx = 0; idx < 3; ++idx) {
        const int dim = dims[idx];
        MatrixF64 a = MATRIX_INIT;
        MatrixF64 b = MATRIX_INIT;
        MatrixF64 r = MATRIX_INIT;
        matf64_create(&a, dim, dim);
        matf64_create(&b, dim, dim);
        matf64_create(&r, dim, dim);
        matf64_fill(&a, 1.0);
        matf64_fill(&b, 2.0);

        BENCHMARK(std::string("MatrixF64 (BLAS) Add [") + labels[idx] + "]") {
            return matf64_add(&a, &b);
        };

        BENCHMARK(std::string("MatrixF64 (BLAS) Mul [") + labels[idx] + "]") {
            return matf64_mul(&a, false, &b, false, &r);
        };

        matf64_destroy(&a);
        matf64_destroy(&b);
        matf64_destroy(&r);
    }
}

TEST_CASE("Matrix Creation and Copy Performance (Integer Types)", "[matrixi32][matrixi64][benchmark]") {
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

    MatrixI32 i32_src = MATRIX_INIT, i32_dst = MATRIX_INIT;
    mati32_create(&i32_src, dim, dim);
    mati32_create(&i32_dst, dim, dim);
    BENCHMARK("MatrixI32 Copy [512x512]") {
        return mati32_copy(&i32_dst, &i32_src);
    };
    mati32_destroy(&i32_src);
    mati32_destroy(&i32_dst);

    MatrixI64 i64_src = MATRIX_INIT, i64_dst = MATRIX_INIT;
    mati64_create(&i64_src, dim, dim);
    mati64_create(&i64_dst, dim, dim);
    BENCHMARK("MatrixI64 Copy [512x512]") {
        return mati64_copy(&i64_dst, &i64_src);
    };
    mati64_destroy(&i64_src);
    mati64_destroy(&i64_dst);
}

TEST_CASE("Matrix Creation and Copy Performance (Floating-Point Types)", "[matrixf32][matrixf64][benchmark]") {
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

    MatrixF32 f32_src = MATRIX_INIT, f32_dst = MATRIX_INIT;
    matf32_create(&f32_src, dim, dim);
    matf32_create(&f32_dst, dim, dim);
    BENCHMARK("MatrixF32 Copy [512x512]") {
        return matf32_copy(&f32_dst, &f32_src);
    };
    matf32_destroy(&f32_src);
    matf32_destroy(&f32_dst);

    MatrixF64 f64_src = MATRIX_INIT, f64_dst = MATRIX_INIT;
    matf64_create(&f64_src, dim, dim);
    matf64_create(&f64_dst, dim, dim);
    BENCHMARK("MatrixF64 Copy [512x512]") {
        return matf64_copy(&f64_dst, &f64_src);
    };
    matf64_destroy(&f64_src);
    matf64_destroy(&f64_dst);
}
