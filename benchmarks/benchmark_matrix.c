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

#include <catch2/catch_test_macros.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>

#include "matrix.h"

// Benchmarks for IntMatrix (Integer Operations)
TEST_CASE("IntMatrix Performance", "[intmatrix][benchmark]") {
    // Small Matrix Benchmarks (32x32)
    const int small_dim = 32;
    IntMatrix sm_a, sm_b, sm_r;
    intmat_create(&sm_a, small_dim, small_dim);
    intmat_create(&sm_b, small_dim, small_dim);
    intmat_create(&sm_r, small_dim, small_dim);
    intmat_fill(&sm_a, 1);
    intmat_fill(&sm_b, 2);

    BENCHMARK("IntMatrix Add [32x32]") {
        return intmat_add(&sm_a, &sm_b);
    };

    BENCHMARK("IntMatrix Mul [32x32]") {
        return intmat_mul(&sm_a, false, &sm_b, false, &sm_r);
    };

    intmat_destroy(&sm_a);
    intmat_destroy(&sm_b);
    intmat_destroy(&sm_r);

    // Medium Matrix Benchmarks (256x256)
    const int med_dim = 256;
    IntMatrix mm_a, mm_b, mm_r;
    intmat_create(&mm_a, med_dim, med_dim);
    intmat_create(&mm_b, med_dim, med_dim);
    intmat_create(&mm_r, med_dim, med_dim);
    intmat_fill(&mm_a, 1);
    intmat_fill(&mm_b, 2);

    BENCHMARK("IntMatrix Add [256x256]") {
        return intmat_add(&mm_a, &mm_b);
    };

    BENCHMARK("IntMatrix Mul [256x256]") {
        return intmat_mul(&mm_a, false, &mm_b, false, &mm_r);
    };

    intmat_destroy(&mm_a);
    intmat_destroy(&mm_b);
    intmat_destroy(&mm_r);
}


// Benchmarks for Matrix (Double Precision Operations with BLAS)
TEST_CASE("Matrix Performance", "[matrix][benchmark]") {
    // Small Matrix Benchmarks (32x32)
    const int small_dim = 32;
    Matrix sm_a, sm_b, sm_r;
    mat_create(&sm_a, small_dim, small_dim);
    mat_create(&sm_b, small_dim, small_dim);
    mat_create(&sm_r, small_dim, small_dim);
    mat_fill(&sm_a, 1.0);
    mat_fill(&sm_b, 2.0);

    BENCHMARK("Matrix (BLAS) Add [32x32]") {
        return mat_add(&sm_a, &sm_b);
    };

    BENCHMARK("Matrix (BLAS) Mul [32x32]") {
        return mat_mul(&sm_a, false, &sm_b, false, &sm_r);
    };

    mat_destroy(&sm_a);
    mat_destroy(&sm_b);
    mat_destroy(&sm_r);

    // Medium Matrix Benchmarks (256x256)
    const int med_dim = 256;
    Matrix mm_a, mm_b, mm_r;
    mat_create(&mm_a, med_dim, med_dim);
    mat_create(&mm_b, med_dim, med_dim);
    mat_create(&mm_r, med_dim, med_dim);
    mat_fill(&mm_a, 1.0);
    mat_fill(&mm_b, 2.0);

    BENCHMARK("Matrix (BLAS) Add [256x256]") {
        return mat_add(&mm_a, &mm_b);
    };

    BENCHMARK("Matrix (BLAS) Mul [256x256]") {
        return mat_mul(&mm_a, false, &mm_b, false, &mm_r);
    };

    mat_destroy(&mm_a);
    mat_destroy(&mm_b);
    mat_destroy(&mm_r);

    // Large Matrix Benchmarks (1024x1024)
    const int large_dim = 1024;
    Matrix lm_a, lm_b, lm_r;
    mat_create(&lm_a, large_dim, large_dim);
    mat_create(&lm_b, large_dim, large_dim);
    mat_create(&lm_r, large_dim, large_dim);
    mat_fill(&lm_a, 1.0);
    mat_fill(&lm_b, 2.0);

    BENCHMARK("Matrix (BLAS) Add [1024x1024]") {
        return mat_add(&lm_a, &lm_b);
    };

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
        IntMatrix mat;
        intmat_create(&mat, dim, dim);
        intmat_destroy(&mat);
        return 0;
    };

    BENCHMARK("Matrix Creation [512x512]") {
        Matrix mat;
        mat_create(&mat, dim, dim);
        mat_destroy(&mat);
        return 0;
    };

    IntMatrix int_src, int_dst;
    intmat_create(&int_src, dim, dim);
    BENCHMARK("IntMatrix Copy [512x512]") {
        intmat_copy(&int_dst, &int_src);
        intmat_destroy(&int_dst);
        return 0;
    };
    intmat_destroy(&int_src);


    Matrix double_src, double_dst;
    mat_create(&double_src, dim, dim);
    BENCHMARK("Matrix Copy [512x512]") {
        mat_copy(&double_dst, &double_src);
        mat_destroy(&double_dst);
        return 0;
    };
    mat_destroy(&double_src);
}
