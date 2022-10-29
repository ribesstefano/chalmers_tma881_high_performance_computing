CXX := gcc
CFLAGS := -std=c17 -O3 -march=native -fopenmp -g
LDFLAGS = -lm

SOURCES = $(wildcard *.c)
EXECS = $(patsubst %.c, %.exe, $(SOURCES))

BLOCK_DIMS = -DBLOCK_SIZE_X=2048 -DBLOCK_SIZE_Y=2048
CFLAGS += $(BLOCK_DIMS)

.PHONY: all compress

all: cell_distances

cell_distances: cell_distances.c
	$(CXX) $(CFLAGS) $< -o cell_distances $(LDFLAGS)

link_cell_e4:
	rm cells
	ln -s /home/hpc2022/cell_distances/test_data/cell_e4 cells

link_cell_e5:
	rm cells
	ln -s /home/hpc2022/cell_distances/test_data/cell_e5 cells

compress: Makefile cell_distances.c
	tar -czvf submission.tar.gz Makefile cell_distances.c

debug: cell_distances link_cell_e4
	hyperfine --time-unit millisecond --warmup 1 --max-runs 10 "systemd-run --user --scope --property=MemoryMax=5M ./cell_distances -t1"
	perf record ./cell_distances
	perf report

perf_e4: cell_distances link_cell_e4
	hyperfine --time-unit millisecond --warmup 1 --max-runs 10 "systemd-run --user --scope --property=MemoryMax=5M ./cell_distances -t1"

perf_e5: cell_distances link_cell_e5
	hyperfine  --time-unit millisecond --warmup 1 --max-runs 10 "systemd-run --user --scope --property=MemoryMax=5M ./cell_distances -t20"

check: compress
	rm -rf distances extracted
	python3 /home/hpc2022/cell_distances/check_submission.py submission.tar.gz

clean:
	rm -rf *.exe *.o *.s ./cell_distances distances extracted submission* vgcore.* perf.data*
