/* Library for reading, writing and manipulating matrix data.
 
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

#include <stdio.h>
#include <stdlib.h>
#include <cblas.h>
#include <time.h>
#include <stdbool.h>
#include <math.h>
#include "matrix.h"


/************************************************************/
/*******Basic C implementations of BLAS functions************/
/*******for which integer or double implementations**********/
/*******are not available in OpenBLAS.***********************/
/************************************************************/

/************************************************************/
/***********Integer matrix operations************************/
/************************************************************/

// Scales a vector x := alpha * x. (?scal)
void iscal(const unsigned int num_elem, const int alpha,
           int *x, const unsigned int incx)
{
    if (x==NULL || incx==0)
        return;
    for (size_t i=0; i<num_elem; i++)
        x[i*incx] *= alpha;
}

// Copies a vector y := x. (?copy)
void icopy(const unsigned int num_elem, const int* x,
           const unsigned int incx, int* y,
           const unsigned int incy)
{
    if (x==NULL || y==NULL)
        return;

    if (incx==0 || incy==0)
        return;
    
    for (size_t i=0; i<num_elem; i++)
        y[i*incy] = x[i*incx];
}

// Scales a vector x and adds it to another vector y.(?axpy)
// y := alpha * x + y
void iaxpy(const unsigned int num_elem, const int alpha,
           const int* x, const unsigned int incx,
           int* y, const unsigned int incy)
{
    if (x==NULL || y==NULL)
        return;

    if (incx==0 || incy==0)
        return;

    for (size_t i=0; i<num_elem; i++)
        y[i*incy] += alpha * x[i*incx];
}

// General dense matrix multiplication for integer matrices. (?gemm)
// C := alpha * op(A) op(B) + beta * C
// Conventions for ?gemm in the BLAS standard are used.
// Please refer to the documentation for details.
// (e.g. "https://www.intel.com/content/www/us/en
// /docs/onemkl/developer-reference-c/2023-1/cblas-gemm
// -001.html#GUID-97718E5C-6E0A-44F0-B2B1-A551F0F164B2")
void igemm(const CBLAS_TRANSPOSE transa, const CBLAS_TRANSPOSE transb, 
        const unsigned int m, const unsigned int n, const unsigned int k, 
        const int alpha, const int* a, const unsigned int lda, 
        const int* b, const unsigned int ldb,
        const int beta, int* c, const unsigned int ldc)
{
    int* a_trans = NULL;
    int* b_trans = NULL;
    unsigned int min_lda, min_ldb, min_ldc;
    bool a_transpose_p = transa==CblasTrans;
    bool b_transpose_p = transb==CblasTrans;

    // Validation
    if (a==NULL || b==NULL || c==NULL)
    {
        perror("ERROR! One of the matrices is a null pointer.");
        return;
    }
    if ((transa==CblasConjTrans || transa==CblasConjNoTrans) ||
        (transb==CblasConjTrans || transb==CblasConjNoTrans))
    {
        perror("ERROR! Conjugate matrices not supported.");
        return;
    }
    if (a_transpose_p)
        min_lda = m > 1? m: 1;
    else
        min_lda = k > 1? k: 1;
    if (lda!=min_lda)
    {
        perror("ERROR! Invalid LDA.");
        return;
    }
    if (b_transpose_p)
        min_ldb = k > 1? k: 1;
    else
        min_ldb = n > 1? n: 1;
    if (ldb!=min_ldb)
    {
        perror("ERROR! Invalid LDB.");
        return;
    }
    min_ldc = n > 1? n: 1;
    if (ldc!=min_ldc)
    {
        perror("ERROR! Invalid LDC.");
        return;
    }

    // Quick return
    if (alpha==0 && beta==1)
        return;
    if (alpha==0)
    {
        for (size_t i=0; i<ldc*m; i++)
            c[i] *= beta;
        return;
    }

    // Transpositions
    if (a_transpose_p)
    {
        a_trans = (int *)calloc(m*k, sizeof(int));
        for (size_t i=0; i<m; i++)
            for (size_t j=0; j<k; j++)
                a_trans[i*k+j] = a[j*lda+i];
    }
    if (b_transpose_p)
    {
        b_trans = (int *)calloc(n*k, sizeof(int));
        for (size_t i=0; i<k; i++)
            for (size_t j=0; j<n; j++)
                b_trans[i*n+j] = b[j*ldb+i];
    }
    
    // Matrix multiply
    if (a_transpose_p && !b_transpose_p)
    {
        for (size_t i=0; i<m; ++i)
            for (size_t j=0; j<n; ++j)
                for (size_t l=0; l<k; ++l)
                    c[i*ldc+j] += alpha * a_trans[i*k+l] * b[l*ldb+j] +\
                                            beta * c[i*ldc+j];
    }
    else if (!a_transpose_p && b_transpose_p)
    {
        for (size_t i=0; i<m; ++i)
            for (size_t j=0; j<n; ++j)
                for (size_t l=0; l<k; ++l)
                    c[i*ldc+j] += alpha * a[i*lda+l] * b_trans[l*n+j] +\
                                            beta * c[i*ldc+j];
    }
    else if (a_transpose_p && b_transpose_p)
    {
        for (size_t i=0; i<m; ++i)
            for (size_t j=0; j<n; ++j)
                for (size_t l=0; l<k; ++l)
                    c[i*ldc+j] += alpha * a_trans[i*k+l] * b_trans[l*n+j] +\
                                            beta * c[i*ldc+j];
    }
    else
    {
        for (size_t i=0; i<m; ++i)
            for (size_t j=0; j<n; ++j)
                for (size_t l=0; l<k; ++l)
                {
                    c[i*ldc+j] += alpha * a[i*lda+l] * b[l*ldb+j] +\
                                            beta * c[i*ldc+j];
                }
    }

    // Free transposed matrices
    if (a_transpose_p)
        free(a_trans);
    if (b_transpose_p)
        free(b_trans);
}

// Sparse BLAS-like function for gathering elements from a
// sparse storage vector to a dense storage vector according
// to supplied indices.
void iusga(const unsigned int num_elem,
           const int* y, const unsigned int incy,
           int* x, const unsigned int* idxs)
{
    if (x==NULL || y==NULL || idxs==NULL)
    {
        perror("ERROR! Null pointer in array argument(s).");
        return;
    }
    for (size_t i=0; i<num_elem; i++)
        x[i] = y[idxs[i*incy]];
}

// Sparse BLAS-like function for gathering elements from a
// sparse storage vector to a dense storage vector according
// to supplied indices. Double precision arrays supported.
void dusga(const unsigned int num_elem,
           const double* y, const unsigned int incy,
           double* x, const unsigned int* idxs)
{
    if (x==NULL || y==NULL || idxs==NULL)
    {
        perror("ERROR! Null pointer in array argument(s).");
        return;
    }
    for (size_t i=0; i<num_elem; i++)
        x[i] = y[idxs[i*incy]];
}

/************************************************************************/
/***************Functions for IntMatrix (integer data)******************/
/************************************************************************/

// Create a matrix
IntMatrix intmat_create(int nrows, int ncols)
{
    IntMatrix matrix;
    matrix.nrows = (unsigned int)nrows; 
    matrix.ncols = (unsigned int)ncols;

    if (nrows<=0 || ncols<=0)
    {
        perror("ERROR: Number of rows/columns cannot be <= 0.");
        matrix.data = NULL;
        return matrix;
    }
    
    matrix.data = (int *)calloc(matrix.nrows * matrix.ncols, 
                                   sizeof(int));

    return matrix;
}

// Print a matrix
void intmat_print(const IntMatrix* mat)
{
    if (mat==NULL)
        return;
    if (mat->data==NULL)
        return;

    for (size_t i=0; i<mat->nrows; i++)
    {
        for (size_t j=0; j<mat->ncols; j++)
            printf("%d ", mat->data[i*mat->ncols+j]);
        printf("\n");
    }
}

// Fill a matrix with a single value
void intmat_fill(IntMatrix* mat, int value)
{
    if (mat==NULL)
        return;
    if (mat->data==NULL)
        return;
    for (size_t i=0; i<mat->nrows; i++)
        for (size_t j=0; j<mat->ncols; j++)
            mat->data[i*mat->ncols+j] = value;
}

// Fill a matrix with random integers between low and high (exclusive)
// with or without replacement.
void intmat_fill_random(IntMatrix* mat, 
                    int low, int high, bool replace, unsigned int seed)
{
    int* temp_ints = NULL;

    srand(seed);

    if (low>=high)
    {
        perror("ERROR: low >= high.");
        intmat_destroy(mat);
        return;
    }
    
    // Random numbers with replacement
    if (replace)
    {
        for (size_t i=0; i<mat->nrows; i++)
            for (size_t j=0; j<mat->ncols; j++)
                mat->data[i*mat->ncols+j] = low + rand()%(high - low);
        return;
    }

    if (mat->nrows*mat->ncols > (size_t)(high-low))
    {
        perror("ERROR: Too many numbers to generate without replacement.");
        intmat_destroy(mat);
        return;
    }
    
    // Random numbers without replacement
    // 1. Generate numbers from low to high (exclusive) and store in arr.
    // 2. Shuffle arr (e.g. using Fisher-Yates).
    // 3. Select first nrows * ncols numbers.
    temp_ints = (int *)calloc(high-low, sizeof(int));
    for (int i=low; i<high; i++)
        temp_ints[i-low] = i;
    for (size_t i=high-low-1; i>0; i--)
    {
        int idx = rand()%i;
        int temp = temp_ints[idx];
        temp_ints[idx] = temp_ints[i];
        temp_ints[i] = temp;
    }
    for (size_t i=0; i<mat->nrows; i++)
        for (size_t j=0; j<mat->ncols; j++)
            mat->data[i*mat->ncols+j] = temp_ints[i*mat->ncols+j];
    free(temp_ints);
}

// Copy a matrix
IntMatrix intmat_copy(IntMatrix* mat)
{
    IntMatrix copy = intmat_create(mat->nrows, mat->ncols);
    icopy(mat->nrows*mat->ncols, mat->data,
                1, copy.data, 1);
    return copy;
}

// Copy a matrix inplace
void intmat_copy_inplace(IntMatrix* mat, IntMatrix* copy)
{
    // Check dimensions
    if (mat->nrows!=copy->nrows || mat->ncols!=copy->ncols)
    {
        perror("ERROR: copy and original matrix must have same dimensions.");
        intmat_destroy(mat);
        return;
    }
    icopy(mat->nrows*mat->ncols, mat->data,
                1, copy->data, 1);
}

// Create a row (1)/column (0) vector (according to specified 
// dimension) with elements from "low" to "high" 
// (excluded) in "step" steps.
IntMatrix intmat_range(int low,
                  int high,
                  unsigned int step,
                  unsigned int dimension)
{
    IntMatrix out;
    bool err = false;
    unsigned int n_elem;

    if (high<=low)
    {
        perror("Upper limit cannot be <= lower limit.");
        err = true;
    }
    if (step==0)
    {
        perror("Step cannot be zero.");
        err = true;
    }

    n_elem = (high - low)/step;
    
    switch (dimension)
    {
        case 0:
            out = intmat_create(n_elem, 1);
            break;
        case 1:
            out = intmat_create(1, n_elem);
            break;
        default:
            perror("Dimension must be either 0 or 1.");
            err = true;
            break;
    }

    if (err)
    {
        intmat_destroy(&out);
        return out;
    }

    for (int i=low; i<high; i+=step)
        out.data[(i-low)/step] = i;

    return out;
}

// Scale a matrix by a scalar
void intmat_scale(IntMatrix* mat, int fac)
{
    iscal(mat->nrows*mat->ncols, fac,
          mat->data, 1);
}

// Add two matrices
// Addition is performed as A := A+B
void intmat_add(IntMatrix* intmat_a, IntMatrix* intmat_b)
{
    // Ensure that both vectors are of same length
    if (intmat_a->nrows!=intmat_b->nrows || intmat_a->ncols!=intmat_b->ncols)
    {
        perror("ERROR: matrices A and B must be of same dimension.");
        intmat_destroy(intmat_a);
        return;
    }

    iaxpy(intmat_b->nrows * intmat_b->ncols, 1, 
          intmat_b->data, 1, intmat_a->data, 1);
}

// Subtract two matrices
// Subtraction is performed as A := A - B
void intmat_sub(IntMatrix* intmat_a, IntMatrix* intmat_b)
{
    IntMatrix copy = intmat_copy(intmat_b);
    intmat_scale(&copy, -1);
    intmat_add(intmat_a, &copy);
    intmat_destroy(&copy);
}

// Add a scalar to a matrix
void intmat_add_scalar(IntMatrix* mat, int scalar)
{
    IntMatrix temp = intmat_create(mat->nrows, mat->ncols);
    intmat_fill(&temp, scalar);
    intmat_add(mat, &temp);
    intmat_destroy(&temp);
}

// Multiply two matrices A and B. Matrices are multiplied
// after transforming them. IntMatrix dimensions must be such that
//     dim(transform(A)) = m x k
//     dim(transform(B)) = k' x n
IntMatrix intmat_mul(IntMatrix* intmat_a, 
               bool transpose_a, 
               IntMatrix* intmat_b,
               bool transpose_b)
{
    unsigned int lda, ldb;
    unsigned int m, n, k, k_prime;

    IntMatrix result;
    int alpha = 1, beta = 0;
    
    CBLAS_TRANSPOSE trans_a = transpose_a? CblasTrans : CblasNoTrans;
    CBLAS_TRANSPOSE trans_b = transpose_b? CblasTrans : CblasNoTrans;
    
    m = transpose_a? intmat_a->ncols: intmat_a->nrows;
    k = transpose_a? intmat_a->nrows: intmat_a->ncols;
    k_prime = transpose_b? intmat_b->ncols: intmat_b->nrows;
    n = transpose_b? intmat_b->nrows: intmat_b->ncols;
    
    // Ensure correct dimensions for matrix multiplication
    if (k!=k_prime)
    {
        result.data = NULL;
        perror("ERROR: IntMatrix dimensions must satisfy k = k' for multiplying m x k and k' x n matrices.");
        return result;
    }
    
    result = intmat_create(m, n);
    intmat_fill(&result, 0);

    lda = transpose_a? m: k;
    ldb = transpose_b? k: n;
    
    igemm(trans_a, trans_b,
          m, n, k,
          alpha, intmat_a->data, lda,
          intmat_b->data, ldb, beta,
          result.data, n);
    return result;
}

// Multiply two matrices and store result in place
void intmat_mul_inplace(IntMatrix* intmat_a, 
                 bool transpose_a, 
                 IntMatrix* intmat_b,
                 bool transpose_b,
                 IntMatrix* result)
{
    unsigned int m, n, k, k_prime;
    unsigned int lda, ldb;
    int alpha = 1, beta = 0;

    CBLAS_TRANSPOSE trans_a = transpose_a? CblasTrans : CblasNoTrans;
    CBLAS_TRANSPOSE trans_b = transpose_b? CblasTrans : CblasNoTrans;

    m = transpose_a? intmat_a->ncols: intmat_a->nrows;
    k = transpose_a? intmat_a->nrows: intmat_a->ncols;
    k_prime = transpose_b? intmat_b->ncols: intmat_b->nrows;
    n = transpose_b? intmat_b->nrows: intmat_b->ncols;
    
    // Ensure correct dimensions for matrix multiplication
    if (k!=k_prime)
    {
        perror("ERROR: IntMatrix dimensions must satisfy k = k' for multiplying m x k and k' x n matrices.");
        intmat_destroy(result);
        return;
    }
    if (result->nrows!=m || result->ncols!=n)
    {
        perror("ERROR: Incorrect dimensions of result matrix.");
        intmat_destroy(result);
        return;
    }
    
    intmat_fill(result, 0);

    lda = transpose_a? m: k;
    ldb = transpose_b? k: n;
    
    igemm(trans_a, trans_b,
          m, n, k,
          alpha, intmat_a->data, lda,
          intmat_b->data, ldb, beta,
          result->data, n);
}

// Repeat a vector along a given dimension
IntMatrix intmat_repeat(IntMatrix* vec, 
        unsigned int dimension, unsigned int repeats)
{
    IntMatrix repeated;
    unsigned int nidx, idx_fac, inc;
    bool err = false;

    // Check if vec is a vector
    if (vec->nrows>1 && vec->ncols>1)
    {
        perror("ERROR: Only a one dimensional vector can be repeated.");
        intmat_destroy(&repeated);
        return repeated;
    }
    
    switch (dimension)
    {
        case 0:
            if (vec->nrows!=1)
            {
                perror("ERROR: Can only repeat a row vector along rows.");
                err = true;
            }
            nidx = vec->ncols;
            idx_fac = vec->ncols;
            inc = 1;
            repeated = intmat_create(repeats, vec->ncols);
            break;
        case 1:
            if (vec->ncols!=1)
            {
                perror("ERROR: Can only repeat a column vector along rows.");
                err = true;
            }
            nidx = vec->nrows;
            idx_fac = 1;
            inc = repeats;
            repeated = intmat_create(vec->nrows, repeats);
            break;
        default:
            perror("ERROR: Dimension must be either 0 or 1.\n");
            err = true;
            break;
    }
    if (err)
    {
        intmat_destroy(&repeated);
        return repeated;
    }
    
    for (size_t i=0; i<repeats; i++)
        icopy(nidx, vec->data, 1, &(repeated.data[idx_fac*i]), inc);

    return repeated;
}

// Add a vector to a matrix
// Addition is done as: A := A + B
// where vector B is repeated along the number
// of dimensions as required to match A's dimensions.
void intmat_vec_add(IntMatrix* mat, IntMatrix* vec)
{
    IntMatrix repeated;
    if (vec->nrows>1 && vec->ncols>1)
    {
        perror("ERROR: Second argument must be a vector.");
        intmat_destroy(mat);
        return;
    }
    if (vec->nrows==1)
        repeated = intmat_repeat(vec, 0, mat->nrows);
    else
        repeated = intmat_repeat(vec, 1, mat->ncols);
    intmat_add(mat, &repeated);
    intmat_destroy(&repeated);
}

// Subtract a vector from a matrix
// Subtraction is done as: A := A - B
// where vector B is repeated along the number
// of dimensions as required to match A's dimensions.
void intmat_vec_sub(IntMatrix* mat, IntMatrix* vec)
{
    IntMatrix repeated;
    if (vec->nrows>1 && vec->ncols>1)
    {
        perror("ERROR: Second argument must be a vector.");
        intmat_destroy(mat);
        return;
    }
    if (vec->nrows==1)
        repeated = intmat_repeat(vec, 0, mat->nrows);
    else
        repeated = intmat_repeat(vec, 1, mat->ncols);
    intmat_sub(mat, &repeated);
    intmat_destroy(&repeated);
}

// Gather rows/columns from "from" and store in
// "to" according to specified indices.
void intmat_gather(IntMatrix* from,
                IntMatrix* to,
                IntMatrix* indices,
                unsigned int dimension)
{
    IntMatrix idxs, idxs_repeated, ind_arg_cpy;
    switch (dimension)
    {
        case 0:
            idxs = intmat_range(0, from->ncols, 1, 1);
            idxs_repeated = intmat_repeat(&idxs, 0, 
                                    indices->nrows);
            ind_arg_cpy = intmat_copy(indices);
            intmat_scale(&ind_arg_cpy, from->ncols); 
            intmat_vec_add(&idxs_repeated, &ind_arg_cpy);
            break;
        case 1:
            idxs = intmat_range(0, from->nrows, 1, 0);
            intmat_scale(&idxs, from->ncols);
            idxs_repeated = intmat_repeat(&idxs, 1,
                                indices->ncols);
            intmat_vec_add(&idxs_repeated, indices);
            break;
        default:
            perror("Dimension must be either rows(0) or columns(1).");
            intmat_destroy(to);
            return;
    }
    iusga(idxs_repeated.nrows*idxs_repeated.ncols, from->data,
             1, to->data, (unsigned int*)idxs_repeated.data);
    intmat_destroy(&ind_arg_cpy);
    intmat_destroy(&idxs);
    intmat_destroy(&idxs_repeated);
}

// Destroy a matrix
void intmat_destroy(IntMatrix* matrix)
{
    if (matrix==NULL)
        return;
    if ((matrix->data)!=NULL)
        free(matrix->data);
    matrix->data = NULL;
}

/************************************************************************/
/***************Functions for Matrix (double precision data)*************/
/************************************************************************/

// Create a matrix
Matrix mat_create(int nrows, int ncols)
{
    Matrix matrix;
    matrix.nrows = (unsigned int)nrows; 
    matrix.ncols = (unsigned int)ncols;

    if (nrows<=0 || ncols<=0)
    {
        perror("ERROR: Number of rows/columns cannot be <= 0.");
        matrix.data = NULL;
        return matrix;
    }
    
    matrix.data = (double *)calloc(matrix.nrows * matrix.ncols, 
                                   sizeof(double));

    return matrix;
}

// Print a matrix
void mat_print(const Matrix* mat)
{
    if (mat==NULL)
        return;
    if (mat->data==NULL)
        return;

    for (size_t i=0; i<mat->nrows; i++)
    {
        for (size_t j=0; j<mat->ncols; j++)
            printf("%.4f ", mat->data[i*mat->ncols+j]);
        printf("\n");
    }
}

// Fill a matrix with a single value
void mat_fill(Matrix* mat, double value)
{
    if (mat==NULL)
        return;
    if (mat->data==NULL)
        return;
    for (size_t i=0; i<mat->nrows; i++)
        for (size_t j=0; j<mat->ncols; j++)
            mat->data[i*mat->ncols+j] = value;
}

// Fill a matrix with random numbers between 0.0 and 1.0 (half-open)
void mat_fill_random(Matrix* mat, unsigned int seed)
{
    srand(seed);

    for (size_t i=0; i<mat->nrows; i++)
        for (size_t j=0; j<mat->ncols; j++)
            mat->data[i*mat->ncols+j] = (double)rand()/(double)(RAND_MAX);
}

// Fill a matrix with random numbers from a Gaussian distribution 
// with mean "mean" and standard deviation "std"
void mat_fill_random_gaussian(Matrix* mat,
                              Matrix* means,
                              Matrix* stds,
                              unsigned int seed)
{
    double x;

    srand(seed);

    if (mat==NULL || means==NULL || stds==NULL)
    {
        perror("ERROR: Null values in argument matrices.");
        return;
    }

    if (means->nrows!=mat->ncols || stds->nrows!=mat->ncols 
            || means->ncols!=1 || stds->ncols!=1)
    {
        perror("ERROR: Means/Stddevs matrix is of improper dimensions.");
        mat_destroy(mat);
        return;
    }

    for (size_t i=0; i<mat->nrows; i++)
        for (size_t j=0; j<mat->ncols; j++)
        {
            double mean = means->data[j];
            double std = stds->data[j];

            if (std==0.0)
            {
                perror("ERROR: Standard deviation cannot be zero.");
                mat_destroy(mat);
                return;
            }

            x = (double)rand()/(double)(RAND_MAX);
            mat->data[i*mat->ncols+j] = exp(-0.5*((x-mean)/std)*((x-mean)/std))\
                                        / (std * sqrt(2*M_PI));
        }
}

// Copy a matrix
Matrix mat_copy(Matrix* mat)
{
    Matrix copy = mat_create(mat->nrows, mat->ncols);
    cblas_dcopy(mat->nrows*mat->ncols, mat->data,
                1, copy.data, 1);
    return copy;
}

// Copy a matrix inplace
void mat_copy_inplace(Matrix* mat, Matrix* copy)
{
    // Check dimensions
    if (mat->nrows!=copy->nrows || mat->ncols!=copy->ncols)
    {
        perror("ERROR: copy and original matrix must have same dimensions.");
        mat_destroy(mat);
        return;
    }
    cblas_dcopy(mat->nrows*mat->ncols, mat->data,
                1, copy->data, 1);
}

// Sum of absolute values of matrix elements
double mat_abs_sum(Matrix* mat)
{
    if (mat==NULL)
        return 0.0;

    double sum = 0.0;
    sum = cblas_dasum(mat->nrows*mat->ncols, mat->data, 1);
    return sum;
}

// Euclidean norm of a matrix across all rows and 
// columns. 
double mat_norm(Matrix* mat)
{
    if (mat==NULL)
    {
        perror("ERROR: Null pointer in argument matrix.");
        return 0.0;
    }
    return cblas_dnrm2(mat->nrows*mat->ncols, 
                        mat->data, 1);
}

// Create a row (1)/column (0) vector (according to specified 
// dimension) with elements from "low" to "high" 
// (excluded) in "step" steps.
Matrix mat_range(double low,
                  double high,
                  double step,
                  unsigned int dimension)
{
    Matrix out;
    bool err = false;
    unsigned int n_elem;

    if (high<=low)
    {
        perror("Upper limit cannot be <= lower limit.");
        err = true;
    }
    if (step<=0)
    {
        perror("Step cannot be zero or negative.");
        err = true;
    }

    n_elem = (int)((high - low)/step);
    
    switch (dimension)
    {
        case 0:
            out = mat_create(n_elem, 1);
            break;
        case 1:
            out = mat_create(1, n_elem);
            break;
        default:
            perror("Dimension must be either 0 or 1.");
            err = true;
            break;
    }

    if (err)
    {
        mat_destroy(&out);
        return out;
    }

    for (int i=low; i<high; i+=step)
        out.data[(int)((i-low)/step)] = (double)i;

    return out;
}

// Scale a matrix by a scalar
void mat_scale(Matrix* mat, double fac)
{
    cblas_dscal(mat->nrows*mat->ncols, fac,
                mat->data, 1);
}

// Add a scalar to a matrix
void mat_add_scalar(Matrix* mat, double scalar)
{
    Matrix temp = mat_create(mat->nrows, mat->ncols);
    mat_fill(&temp, scalar);
    mat_add(mat, &temp);
    mat_destroy(&temp);
}

// Multiply two matrices A and B. Matrices are multiplied
// after transforming them. Matrix dimensions must be such that
//     dim(transform(A)) = m x k
//     dim(transform(B)) = k' x n
Matrix mat_mul(Matrix* mat_a, 
               bool transpose_a, 
               Matrix* mat_b,
               bool transpose_b)
{
    unsigned int lda, ldb;
    unsigned int m, n, k, k_prime;

    Matrix result;
    
    CBLAS_TRANSPOSE trans_a = transpose_a? CblasTrans : CblasNoTrans;
    CBLAS_TRANSPOSE trans_b = transpose_b? CblasTrans : CblasNoTrans;
    
    m = transpose_a? mat_a->ncols: mat_a->nrows;
    k = transpose_a? mat_a->nrows: mat_a->ncols;
    k_prime = transpose_b? mat_b->ncols: mat_b->nrows;
    n = transpose_b? mat_b->nrows: mat_b->ncols;
    
    // Ensure correct dimensions for matrix multiplication
    if (k!=k_prime)
    {
        result.data = NULL;
        perror("ERROR: Matrix dimensions must satisfy k = k' for multiplying m x k and k' x n matrices.");
        return result;
    }
    
    result = mat_create(m, n);
    mat_fill(&result, 0.0);

    lda = transpose_a? m: k;
    ldb = transpose_b? k: n;
    
    cblas_dgemm(CblasRowMajor, trans_a, trans_b,
                m, n, k,
                1.0, mat_a->data, lda,
                mat_b->data, ldb, 0.0,
                result.data, n);
    return result;
}

// Multiply two matrices and store result in place
void mat_mul_inplace(Matrix* mat_a, 
                 bool transpose_a, 
                 Matrix* mat_b,
                 bool transpose_b,
                 Matrix* result)
{
    unsigned int m, n, k, k_prime;
    unsigned int lda, ldb;

    CBLAS_TRANSPOSE trans_a = transpose_a? CblasTrans : CblasNoTrans;
    CBLAS_TRANSPOSE trans_b = transpose_b? CblasTrans : CblasNoTrans;

    m = transpose_a? mat_a->ncols: mat_a->nrows;
    k = transpose_a? mat_a->nrows: mat_a->ncols;
    k_prime = transpose_b? mat_b->ncols: mat_b->nrows;
    n = transpose_b? mat_b->nrows: mat_b->ncols;
    
    // Ensure correct dimensions for matrix multiplication
    if (k!=k_prime)
    {
        perror("ERROR: Matrix dimensions must satisfy k = k' for multiplying m x k and k' x n matrices.");
        mat_destroy(result);
        return;
    }
    if (result->nrows!=m || result->ncols!=n)
    {
        perror("ERROR: Incorrect dimensions of result matrix.");
        mat_destroy(result);
        return;
    }
    
    mat_fill(result, 0.0);

    lda = transpose_a? m: k;
    ldb = transpose_b? k: n;
    
    cblas_dgemm(CblasRowMajor, trans_a, trans_b,
                m, n, k,
                1.0, mat_a->data, lda,
                mat_b->data, ldb, 0.0,
                result->data, n);
}

// Add two matrices
// Addition is performed as A := A+B
void mat_add(Matrix* mat_a, Matrix* mat_b)
{
    // Ensure that both vectors are of same length
    if (mat_a->nrows!=mat_b->nrows || mat_a->ncols!=mat_b->ncols)
    {
        perror("ERROR: matrices A and B must be of same dimension.");
        mat_destroy(mat_a);
        return;
    }

    cblas_daxpy(mat_b->nrows * mat_b->ncols, 1.0, 
                mat_b->data, 1, mat_a->data, 1);
}

// Subtract two matrices
// Subtraction is performed as A := A - B
void mat_sub(Matrix* mat_a, Matrix* mat_b)
{
    Matrix copy = mat_copy(mat_b);
    mat_scale(&copy, -1.0);
    mat_add(mat_a, &copy);
    mat_destroy(&copy);
}

// Repeat a vector along a given dimension
Matrix mat_repeat(Matrix* vec, unsigned int dimension, unsigned int repeats)
{
    Matrix repeated;
    unsigned int nidx, idx_fac, inc;
    bool err = false;

    // Check if vec is a vector
    if (vec->nrows>1 && vec->ncols>1)
    {
        perror("ERROR: Only a one dimensional vector can be repeated.");
        mat_destroy(&repeated);
        return repeated;
    }
    
    switch (dimension)
    {
        case 0:
            if (vec->nrows!=1)
            {
                perror("ERROR: Can only repeat a row vector along rows.");
                err = true;
            }
            nidx = vec->ncols;
            idx_fac = vec->ncols;
            inc = 1;
            repeated = mat_create(repeats, vec->ncols);
            break;
        case 1:
            if (vec->ncols!=1)
            {
                perror("ERROR: Can only repeat a column vector along rows.");
                err = true;
            }
            nidx = vec->nrows;
            idx_fac = 1;
            inc = repeats;
            repeated = mat_create(vec->nrows, repeats);
            break;
        default:
            perror("ERROR: Dimension must be either 0 or 1.\n");
            err = true;
            break;
    }
    if (err)
    {
        mat_destroy(&repeated);
        return repeated;
    }
    
    for (size_t i=0; i<repeats; i++)
        cblas_dcopy(nidx, vec->data, 1, &(repeated.data[idx_fac*i]), inc);

    return repeated;
}

// Add a vector to a matrix
// Addition is done as: A := A + B
// where vector B is repeated along the number
// of dimensions as required to match A's dimensions.
void mat_vec_add(Matrix* mat, Matrix* vec)
{
    Matrix repeated;
    if (vec->nrows>1 && vec->ncols>1)
    {
        perror("ERROR: Second argument must be a vector.");
        mat_destroy(mat);
        return;
    }
    if (vec->nrows==1)
        repeated = mat_repeat(vec, 0, mat->nrows);
    else
        repeated = mat_repeat(vec, 1, mat->ncols);
    mat_add(mat, &repeated);
    mat_destroy(&repeated);
}

// Subtract a vector from a matrix
// Subtraction is done as: A := A - B
// where vector B is repeated along the number
// of dimensions as required to match A's dimensions.
void mat_vec_sub(Matrix* mat, Matrix* vec)
{
    Matrix repeated;
    if (vec->nrows>1 && vec->ncols>1)
    {
        perror("ERROR: Second argument must be a vector.");
        mat_destroy(mat);
        return;
    }
    if (vec->nrows==1)
        repeated = mat_repeat(vec, 0, mat->nrows);
    else
        repeated = mat_repeat(vec, 1, mat->ncols);
    mat_sub(mat, &repeated);
    mat_destroy(&repeated);
}

// Gather rows/columns from "from" and store in
// "to" according to specified indices.
void mat_gather(Matrix* from,
                Matrix* to,
                IntMatrix* indices,
                unsigned int dimension)
{
    IntMatrix idxs, idxs_repeated, ind_arg_cpy;
    switch (dimension)
    {
        case 0:
            idxs = intmat_range(0, from->ncols, 1, 1);
            idxs_repeated = intmat_repeat(&idxs, 0, 
                                    indices->nrows);
            ind_arg_cpy = intmat_copy(indices);
            intmat_scale(&ind_arg_cpy, from->ncols); 
            intmat_vec_add(&idxs_repeated, &ind_arg_cpy);
            break;
        case 1:
            idxs = intmat_range(0, from->nrows, 1, 0);
            intmat_scale(&idxs, from->ncols);
            idxs_repeated = intmat_repeat(&idxs, 1,
                                indices->ncols);
            intmat_vec_add(&idxs_repeated, indices);
            break;
        default:
            perror("Dimension must be either rows(0) or columns(1).");
            mat_destroy(to);
            return;
    }
    dusga(idxs_repeated.nrows*idxs_repeated.ncols, from->data,
             1, to->data, (const unsigned int*)idxs_repeated.data);

    intmat_destroy(&ind_arg_cpy);
    intmat_destroy(&idxs);
    intmat_destroy(&idxs_repeated);
}

// Destroy a matrix
void mat_destroy(Matrix* matrix)
{
    if (matrix==NULL)
        return;
    if ((matrix->data)!=NULL)
        free(matrix->data);
    matrix->data = NULL;
}
