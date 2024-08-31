#!/bin/sh

# Set up variables
DEP_DIR="../third_party"
OPENBLAS_URL="https://github.com/OpenMathLib/OpenBLAS/releases/download/v0.3.28/OpenBLAS-0.3.28.zip"
OPENBLAS_DIR="OpenBLAS-0.3.28"

if [ "$1" = "" ]; then
    NUM_THREADS=1
else
    NUM_THREADS=$1
fi

# Create directory for dependencies
mkdir $DEP_DIR
cd $DEP_DIR || { echo "cd $DEP_DIR failed."; exit 1; }

# Download and install OpenBLAS
wget $OPENBLAS_URL || { echo "wget $OPENBLAS_URL failed."; exit 2; }
unzip "${OPENBLAS_DIR}.zip"
rm "${OPENBLAS_DIR}.zip" 
cd $OPENBLAS_DIR || { echo "cd $OPENBLAS_DIR failed."; exit 3; }
if [ ! -f libopenblas.a ]; then
    make -j "$NUM_THREADS" || { echo "OpenBLAS build failed."; exit 4; }
else
    echo "libopenblas.a found! Skipping build."
fi
cd ..
