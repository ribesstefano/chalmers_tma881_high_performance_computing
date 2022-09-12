CXX := gcc
CFLAGS := -std=c17 -fopenmp
LDFLAGS = -lm

SOURCES = $(wildcard *.c)
EXECS = $(patsubst %.c, %.exe, $(SOURCES))

.PHONY: all compress

all: cell_distances

cell_distances: cell_distances.c
	$(CXX) $(CFLAGS) -O3 -march=x86-64 $< -o cell_distances $(LDFLAGS)

compress:
	tar -czvf submission.tar.gz Makefile cell_distances.c

test: compress
	rm -rf distances extracted
	python3 /home/hpc2022/cell_distances/check_submission.py submission.tar.gz
	# python3 ../reference_code/cell_distances/check_submission.py submission.tar.gz
	rm -rf distances extracted

clean:
	rm -rf *.exe *.o *.s ./cell_distances submission* distances extracted
