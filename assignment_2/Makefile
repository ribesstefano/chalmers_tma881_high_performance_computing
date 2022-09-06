CXX := gcc
CFLAGS := -std=c17 -fopenmp
LDFLAGS = -lm

SOURCES = $(wildcard *.c)
EXECS = $(patsubst %.c, %.exe, $(SOURCES))

.PHONY: all

all: cell_distances

cell_distances: cell_distances.c
	$(CXX) $(CFLAGS) -O3 -march=x86-64 $< -o cell_distances $(LDFLAGS)

clean:
	rm -rf *.exe *.o *.s ./cell_distances
