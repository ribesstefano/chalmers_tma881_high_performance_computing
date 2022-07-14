CXX := gcc
CFLAGS := -std=gnu99 -fopenmp -m

SOURCES = $(wildcard *.c)
EXECS = $(patsubst %.c, %.exe, $(SOURCES))

.PHONY: all

all: $(EXECS)

cell_distances.exe: cell_distances.c
	$(CXX) $(CFLAGS) -O2 -march=x86-64 $< -o cell_distances.exe

clean:
	rm -rf *.exe *.o *.s
