CXX := gcc
CFLAGS := -std=c17 -fopenmp
LDFLAGS = -lm

SOURCES = $(wildcard *.c)
EXECS = $(patsubst %.c, %.exe, $(SOURCES))

.PHONY: all

all: $(EXECS)

cell_distances.exe: cell_distances.c
	$(CXX) $(CFLAGS) -O3 -march=x86-64 $< -o cell_distances.exe $(LDFLAGS)

hash_table.exe: hash_table.c
	$(CXX) $(CFLAGS) -O2 -march=x86-64 $< -o hash_table.exe $(LDFLAGS)

clean:
	rm -rf *.exe *.o *.s
