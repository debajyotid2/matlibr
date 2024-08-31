CC=g++
CFLAGS=-g -Wl,--whole-archive

OPENBLAS_DIR=${CURDIR}/../third_party/OpenBLAS-0.3.28
LDFLAGS=-lmatlibr -lopenblas -Wl,--no-whole-archive
LIB=-L../build/src -L$(OPENBLAS_DIR)
INCLUDE=-I../src
SRC=$(wildcard *.c)
BIN_DIR=bin
BIN=$(patsubst %.c, $(BIN_DIR)/%, $(SRC))

all: $(BIN_DIR) $(BIN) run

$(BIN_DIR):
	mkdir $(BIN_DIR)

$(BIN_DIR)/%: %.c
	$(CC) $(CFLAGS) $(INCLUDE) $(LIB) $^ -o $@ $(LDFLAGS)

run:
	for bin in $(BIN_DIR)/*; do LD_LIBRARY_PATH=$(OPENBLAS_DIR) ./$$bin; done

clean:
	rm -rf $(BIN_DIR)/*
