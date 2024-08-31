# Examples for `matlibr`

These examples are meant to demonstrate the usage of the functions in `matlibr` in various contexts. `matlibr` supports a number of elementary matrix operations, including scaling, addition, subtraction, multiplication and indexing.

## How to build

Make sure to download dependencies of `matlibr` (ie. OpenBLAS) using the `setup_dependencies.sh` script, and then build `matlibr` first.

From the `examples` folder, to build all examples, please run

```
make build
```

### How to run

To run all examples at once, run

```
make
```

The binaries need the path to `libopenblas.so.1` at runtime. So, to run individual examples, from the binary directory (`bin`), run

```
LD_LIBRARY_PATH=<PATH_TO_OPENBLAS_LIBRARY> ./<EXAMPLE_BINARY>
```
To clean up all binaries, run

```
make clean
```
