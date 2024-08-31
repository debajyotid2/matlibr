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

#include <stdbool.h>
#include <catch2/catch_test_macros.hpp>
#include "../src/matrix.h"

TEST_CASE("Double matrix operations.", "[matrix]")
{
    unsigned int seed = 3224;
    Matrix mymat = mat_create(10, 20);
    Matrix mymat2 = mat_create(20, 2);
    Matrix mymat3 = mat_create(10, 20);
    
    SECTION("Creating a matrix.")
    {
        REQUIRE(mymat.nrows==10);
        REQUIRE(mymat.ncols==20);
    }

    SECTION("Destroying a matrix.")
    {
        mat_destroy(&mymat);

        REQUIRE(mymat.data==NULL);
    }

    SECTION("Copying a matrix.")
    {
        mat_fill_random(&mymat, seed);
        Matrix copy = mat_copy(&mymat);
        
        for (size_t i=0; i<mymat.nrows; i++)
            for (size_t j=0; j<mymat.ncols; j++)
                REQUIRE(copy.data[i*mymat.ncols+j]==\
                        mymat.data[i*mymat.ncols+j]);

        mat_destroy(&copy);
    }

    SECTION("Copying a matrix inplace.")
    {
        mat_fill_random(&mymat, seed);
        Matrix copy = mat_create(mymat.nrows, mymat.ncols);
        mat_copy_inplace(&mymat, &copy);
        
        for (size_t i=0; i<mymat.nrows; i++)
            for (size_t j=0; j<mymat.ncols; j++)
                REQUIRE(copy.data[i*mymat.ncols+j]==\
                        mymat.data[i*mymat.ncols+j]);
        
        mat_destroy(&copy);
    }

    SECTION("Filling a matrix with a single value.")
    {
        mat_fill(&mymat, 23.0);
        
        for (size_t i=0; i<mymat.nrows; i++)
            for (size_t j=0; j<mymat.ncols; j++)
                REQUIRE(mymat.data[i*mymat.ncols+j]==23.0);
    }

    SECTION("Scaling a matrix.")
    {
        mat_fill(&mymat, -2.5);
        mat_scale(&mymat, 17.0);

        for (size_t i=0; i<mymat.nrows; i++)
            for (size_t j=0; j<mymat.ncols; j++)
                REQUIRE(mymat.data[i*mymat.ncols+j]==-42.5);
    }

    SECTION("Adding a scalar to a matrix.")
    {
        mat_fill(&mymat, -2.5);
        mat_add_scalar(&mymat, 17.0);

        for (size_t i=0; i<mymat.nrows; i++)
            for (size_t j=0; j<mymat.ncols; j++)
                REQUIRE(mymat.data[i*mymat.ncols+j]==14.5);
    }

    SECTION("Adding two matrices")
    {
        mat_fill(&mymat, -3.14);
        mat_fill(&mymat3, 6.28);

        mat_add(&mymat3, &mymat);

        for (size_t i=0; i<mymat3.nrows; i++)
            for (size_t j=0; j<mymat3.ncols; j++)
                REQUIRE(mymat3.data[i*mymat3.ncols+j]==3.14);
    }

    SECTION("Subtracting two matrices")
    {
        mat_fill(&mymat, -3.14);
        mat_fill(&mymat3, 6.28);

        mat_sub(&mymat, &mymat3);

        for (size_t i=0; i<mymat.nrows; i++)
            for (size_t j=0; j<mymat.ncols; j++)
                REQUIRE(mymat.data[i*mymat.ncols+j]==-9.42);
    }

    SECTION("Multiplying two matrices.")
    {
        mat_fill(&mymat, 20.0);
        mat_fill(&mymat2, 1.0);
        
        Matrix result = mat_mul(&mymat, false, &mymat2, false);

        for (size_t i=0; i<result.nrows; i++)
            for (size_t j=0; j<result.ncols; j++)
                REQUIRE(result.data[i*result.ncols+j]==400.0);

        mat_destroy(&result);
    }

    SECTION("Multiplying two matrices and storing result in place.")
    {
        mat_fill(&mymat, 20.0);
        mat_fill(&mymat2, 1.0);
        
        Matrix result = mat_create(mymat.nrows, mymat2.ncols);
        
        mat_mul_inplace(&mymat, false, &mymat2, false, &result);

        for (size_t i=0; i<result.nrows; i++)
            for (size_t j=0; j<result.ncols; j++)
                REQUIRE(result.data[i*result.ncols+j]==400.0);

        mat_destroy(&result);
    }

    SECTION("Adding a matrix and a vector.")
    {
        Matrix rowvec = mat_create(1, 20);
        Matrix mymat_copy;

        mat_fill_random(&rowvec, seed);
        mat_fill_random(&mymat, seed);

        mymat_copy = mat_copy(&mymat);

        mat_vec_add(&mymat, &rowvec);

        for (size_t i=0; i<mymat.nrows; i++)
            for (size_t j=0; j<mymat.ncols; j++)
                REQUIRE(mymat.data[i*mymat.ncols+j]==\
                        mymat_copy.data[i*mymat.ncols+j] + rowvec.data[j]);

        mat_destroy(&mymat_copy);
        mat_destroy(&rowvec);
    }

    SECTION("Subtracting a vector from a matrix.")
    {
        Matrix rowvec = mat_create(1, 20);
        Matrix mymat_copy;
        
        mat_fill_random(&rowvec, seed);
        mat_fill_random(&mymat, seed);

        mymat_copy = mat_copy(&mymat);
        
        mat_vec_sub(&mymat, &rowvec);

        for (size_t i=0; i<mymat.nrows; i++)
            for (size_t j=0; j<mymat.ncols; j++)
                REQUIRE(mymat.data[i*mymat.ncols+j]==\
                        mymat_copy.data[i*mymat.ncols+j] - rowvec.data[j]);

        mat_destroy(&mymat_copy);
        mat_destroy(&rowvec);
    }


    SECTION("Repeating a row vector along rows.")
    {
        Matrix rowvec = mat_create(1, 20);
        Matrix repeated;

        mat_fill_random(&rowvec, seed);
        repeated = mat_repeat(&rowvec, 0, 5);

        for (size_t i=0; i<repeated.nrows; i++)
            for (size_t j=0; j<repeated.ncols; j++)
                REQUIRE(repeated.data[i*repeated.ncols+j]==rowvec.data[j]);

        mat_destroy(&rowvec);
        mat_destroy(&repeated);
    }

    SECTION("Repeating a column vector along columns.")
    {
        Matrix colvec = mat_create(10, 1);
        Matrix repeated;

        mat_fill_random(&colvec, seed);
        repeated = mat_repeat(&colvec, 1, 5);

        for (size_t i=0; i<repeated.nrows; i++)
            for (size_t j=0; j<repeated.ncols; j++)
                REQUIRE(repeated.data[i*repeated.ncols+j]==colvec.data[i]);

        mat_destroy(&colvec);
        mat_destroy(&repeated);
    }

    mat_destroy(&mymat);
    mat_destroy(&mymat2);
    mat_destroy(&mymat3);
}

TEST_CASE("Integer matrix operations.", "[matrix]")
{
    unsigned int seed = 3224;
    IntMatrix mymat = intmat_create(10, 20);
    IntMatrix mymat2 = intmat_create(20, 2);
    IntMatrix mymat3 = intmat_create(10, 20);
    
    SECTION("Creating a matrix.")
    {
        REQUIRE(mymat.nrows==10);
        REQUIRE(mymat.ncols==20);
    }

    SECTION("Destroying a matrix.")
    {
        intmat_destroy(&mymat);

        REQUIRE(mymat.data==NULL);
    }

     SECTION("Copying a matrix.")
     {
         intmat_fill_random(&mymat, 0, 10, true, seed);
         IntMatrix copy = intmat_copy(&mymat);
         
         for (size_t i=0; i<mymat.nrows; i++)
             for (size_t j=0; j<mymat.ncols; j++)
                 REQUIRE(copy.data[i*mymat.ncols+j]==\
                         mymat.data[i*mymat.ncols+j]);
 
         intmat_destroy(&copy);
     }
 
     SECTION("Copying a matrix inplace.")
     {
         intmat_fill_random(&mymat, 0, 10, true, seed);
         IntMatrix copy = intmat_create(mymat.nrows, mymat.ncols);
         intmat_copy_inplace(&mymat, &copy);
         
         for (size_t i=0; i<mymat.nrows; i++)
             for (size_t j=0; j<mymat.ncols; j++)
                 REQUIRE(copy.data[i*mymat.ncols+j]==\
                         mymat.data[i*mymat.ncols+j]);
         
         intmat_destroy(&copy);
     }

     SECTION("Filling a matrix with a single value.")
     {
         intmat_fill(&mymat, 23);
         
         for (size_t i=0; i<mymat.nrows; i++)
             for (size_t j=0; j<mymat.ncols; j++)
                 REQUIRE(mymat.data[i*mymat.ncols+j]==23);
     }
 
     SECTION("Scaling a matrix.")
     {
         intmat_fill(&mymat, -2);
         intmat_scale(&mymat, 17);
 
         for (size_t i=0; i<mymat.nrows; i++)
             for (size_t j=0; j<mymat.ncols; j++)
                 REQUIRE(mymat.data[i*mymat.ncols+j]==-34);
     }
 
     SECTION("Adding a scalar to a matrix.")
     {
         intmat_fill(&mymat, -2);
         intmat_add_scalar(&mymat, 17);
 
         for (size_t i=0; i<mymat.nrows; i++)
             for (size_t j=0; j<mymat.ncols; j++)
                 REQUIRE(mymat.data[i*mymat.ncols+j]==15);
     }

    SECTION("Adding two matrices")
    {
        intmat_fill(&mymat, -3);
        intmat_fill(&mymat3, 6);

        intmat_add(&mymat3, &mymat);

        for (size_t i=0; i<mymat3.nrows; i++)
            for (size_t j=0; j<mymat3.ncols; j++)
                REQUIRE(mymat3.data[i*mymat3.ncols+j]==3);
    }

    SECTION("Subtracting two matrices")
    {
        intmat_fill(&mymat, -3);
        intmat_fill(&mymat3, 6);

        intmat_sub(&mymat, &mymat3);

        for (size_t i=0; i<mymat.nrows; i++)
            for (size_t j=0; j<mymat.ncols; j++)
                REQUIRE(mymat.data[i*mymat.ncols+j]==-9);
    }

     SECTION("Multiplying two matrices.")
     {
         intmat_fill(&mymat, 20);
         intmat_fill(&mymat2, 1);
         
         IntMatrix result = intmat_mul(&mymat, false, &mymat2, false);
 
         for (size_t i=0; i<result.nrows; i++)
             for (size_t j=0; j<result.ncols; j++)
                 REQUIRE(result.data[i*result.ncols+j]==400);
 
         intmat_destroy(&result);
     }

     SECTION("Multiplying two matrices and storing result in place.")
     {
         intmat_fill(&mymat, 20);
         intmat_fill(&mymat2, 1);
         
         IntMatrix result = intmat_create(mymat.nrows, mymat2.ncols);
         
         intmat_mul_inplace(&mymat, false, &mymat2, false, &result);
 
         for (size_t i=0; i<result.nrows; i++)
             for (size_t j=0; j<result.ncols; j++)
                 REQUIRE(result.data[i*result.ncols+j]==400);
 
         intmat_destroy(&result);
     }

     SECTION("Adding a matrix and a vector.")
     {
         IntMatrix rowvec = intmat_create(1, 20);
         IntMatrix myintmat_copy;
 
         intmat_fill_random(&rowvec, 0, 20, true, seed);
         intmat_fill_random(&mymat, 0, 30, true, seed);
 
         myintmat_copy = intmat_copy(&mymat);
 
         intmat_vec_add(&mymat, &rowvec);
 
         for (size_t i=0; i<mymat.nrows; i++)
             for (size_t j=0; j<mymat.ncols; j++)
                 REQUIRE(mymat.data[i*mymat.ncols+j]==\
                         myintmat_copy.data[i*mymat.ncols+j] + rowvec.data[j]);
 
         intmat_destroy(&myintmat_copy);
         intmat_destroy(&rowvec);
     }
 
     SECTION("Subtracting a vector from a matrix.")
     {
         IntMatrix rowvec = intmat_create(1, 20);
         IntMatrix myintmat_copy;
         
         intmat_fill_random(&rowvec, -1, 23, true, seed);
         intmat_fill_random(&mymat, 4, 67, true, seed);
 
         myintmat_copy = intmat_copy(&mymat);
         
         intmat_vec_sub(&mymat, &rowvec);
 
         for (size_t i=0; i<mymat.nrows; i++)
             for (size_t j=0; j<mymat.ncols; j++)
                 REQUIRE(mymat.data[i*mymat.ncols+j]==\
                         myintmat_copy.data[i*mymat.ncols+j] - rowvec.data[j]);
 
         intmat_destroy(&myintmat_copy);
         intmat_destroy(&rowvec);
     }
 
 
     SECTION("Repeating a row vector along rows.")
     {
         IntMatrix rowvec = intmat_create(1, 20);
         IntMatrix repeated;
 
         intmat_fill_random(&rowvec, -6, 36, true, seed);
         repeated = intmat_repeat(&rowvec, 0, 5);
 
         for (size_t i=0; i<repeated.nrows; i++)
             for (size_t j=0; j<repeated.ncols; j++)
                 REQUIRE(repeated.data[i*repeated.ncols+j]==rowvec.data[j]);
 
         intmat_destroy(&rowvec);
         intmat_destroy(&repeated);
     }
 
     SECTION("Repeating a column vector along columns.")
     {
         IntMatrix colvec = intmat_create(10, 1);
         IntMatrix repeated;
 
         intmat_fill_random(&colvec, -20, 20, true, seed);
         repeated = intmat_repeat(&colvec, 1, 5);
 
         for (size_t i=0; i<repeated.nrows; i++)
             for (size_t j=0; j<repeated.ncols; j++)
                 REQUIRE(repeated.data[i*repeated.ncols+j]==colvec.data[i]);
 
         intmat_destroy(&colvec);
         intmat_destroy(&repeated);
     }

    intmat_destroy(&mymat);
    intmat_destroy(&mymat2);
    intmat_destroy(&mymat3);
}
