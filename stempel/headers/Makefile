GCC = gcc
LDFLAGS=
# gnu
OMP_FLAG = -std=c99 -O3 -march=native -mtune=native -ffast-math -funroll-loops -fargument-noalias -fopenmp -D_POSIX_C_SOURCE=200112L
LIKWID_FLAGS = -pthread
LIBS = -llikwid
LIKWID_I =$(LIKWID_INC)
LIKWID_L=$(LIKWID_LIB)

SRC = example.c
PROJECT_NAME = project

all: module

module: $(SRC)
	$(GCC) $(OMP_FLAG) ${LIKWID_FLAGS} $+ timing.c dummy.c ${LDFLAGS} -o ../bin/$(PROJECT_NAME) -DLIKWID_PERFMON -I$(LIKWID_I) -L$(LIKWID_L) $(LIBS)

.PHONY: clean

clean:
	rm -rf ../bin/$(PROJECT_NAME) *~
