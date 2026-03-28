#!/bin/sh

set -xe

# Set up variables
DEP_DIR="../third_party"
OPENBLAS_URL="https://github.com/OpenMathLib/OpenBLAS/releases/download/v0.3.28/OpenBLAS-0.3.28.zip"
OPENBLAS_DIR="OpenBLAS-0.3.28"
CFLAGS="-mfma"
OTHER_FLAGS="NO_SHARED=1 STATIC=1"

if [ "$1" = "" ]; then
    NUM_THREADS=1
else
    NUM_THREADS=$1
fi

# Create directory for dependencies
mkdir $DEP_DIR
cd $DEP_DIR

# Download and install OpenBLAS
wget $OPENBLAS_URL
unzip -o "${OPENBLAS_DIR}.zip"
rm "${OPENBLAS_DIR}.zip" 
cd $OPENBLAS_DIR
if [ ! -f libopenblas.a ]; then
    make libs CFLAGS="${CFLAGS}" "${OTHER_FLAGS}" TARGET=HASWELL -j "${NUM_THREADS}"
else
    echo "libopenblas.a found! Skipping build."
fi
cd ..
