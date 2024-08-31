# matlibr

`matlibr` is a barebones implementation of a 2D matrix datatype in C. It also includes some basic arithmetic operations that are relevant to matrices.

## Features

`matlibr` supports
- Double precision and integer data types.
- Basic matrix arithmetic operations like addition, subtraction and multiplication, accelerated by BLAS libraries.
- Matrix-scalar operations like scaling and scalar addition.
- Matrix-vector operations like matrix-vector addition and subtraction.
- Gather operations along rows or columns for index-based access to matrix data.
- Norms of the matrix including L1 and L2 norms.

## Example usage

The following example taken from [`matrix_addition.c`](./examples/matrix_addition.c) displays the creation and destruction of matrices.

```c
#include <stdio.h>
#include <matrix.h>

int main() {
    size_t seed = 42;
    int nrow = 10, ncol = 5;
    Matrix a = mat_create(nrow, ncol);
    
    // Fill matrix a with random values between 0.0 and 1.0
    mat_fill_random(&a, seed);
    
    Matrix b = mat_copy(&a);
    mat_scale(&b, 2.34);
    
    printf("\n************Matrix addition***************\n");
    printf("A: \n");
    mat_print(&a);
    printf(" +\n B: \n");
    mat_print(&b);
    printf(" = \n");
    
    // Addition
    mat_add(&a, &b);
    mat_print(&a);
    printf("\n******************************************\n");
    
    mat_destroy(&a);
    mat_destroy(&b);
}
```

## How to build

The current build only supports Linux-based operating systems.

`g++` is the default compiler used to build sources. If it is not installed, please check [here](https://gcc.gnu.org/install/) for instructions to install `gcc`, or please use your favorite package manager (`dnf` for Fedora/RHEL, `apt` for Ubuntu, `pacman` for Arch Linux, etc.).

After [cloning the repository](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository#cloning-a-repository), we need to set up the dependencies (OpenBLAS). [OpenBLAS](https://github.com/OpenMathLib/OpenBLAS/) can be built by running

```
cd scripts
source setup_dependencies.sh <NUMBER OF THREADS>
```
It is recommended to provide the maximum number of threads on your machine for the second command for the fastest installation. 

**NB 1.** If you wish, you can install OpenBLAS directly to your system by running `source install_openblas.sh <NUMBER OF THREADS>` instead.

### Building

Once the set up is done, to build `matlibr`, please run

```
cd scripts
source build.sh <NUMBER OF THREADS>
```

This builds a static library `libmatlibr.a` which can then be linked to your projects. 

### Installation

For installation, 

```
source install.sh <NUMBER OF THREADS>
```

This installs `matlibr` to the system.

### Testing

```
cd build
ctest -V
```
This will run the unit tests defined in [`tests`](./tests).

## How to include in your project

Since `matlibr` has OpenBLAS as a dependency, it is required that OpenBLAS also be downloaded, built and installed alongside `matlibr`.

Assuming that you have followed the above steps and installed OpenBLAS and `matlibr` to your system, you can simply link the static libraries and includes at compile time.

```
g++ -I<PATH TO OPENBLAS INCLUDE DIR> -I<PATH TO MATLIBR INCLUDE DIR> -L<PATH TO OPENBLAS LIBRARY) -L(PATH TO MATLIBR LIBRARY> foo.c -o foo -lmatlibr -lopenblas
```

## License

[Apache 2.0](https://www.apache.org/licenses/LICENSE-2.0)
