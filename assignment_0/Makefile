CXX := gcc
CFLAGS := -std=gnu99

SOURCES = $(wildcard *.c)
EXECS = $(patsubst %.c, %.exe, $(SOURCES))

.PHONY: all

all: $(EXECS)

stack_alloc.exe: stack_alloc.c
	$(CXX) $(CFLAGS) $< -o $@

heap_alloc.exe: heap_alloc.c
	$(CXX) $(CFLAGS) $< -o $@

mem_fragmentation.exe: mem_fragmentation.c
	$(CXX) $(CFLAGS) $< -o $@

no_mem_fragmentation.exe: no_mem_fragmentation.c
	$(CXX) $(CFLAGS) $< -o $@

write_file.exe: write_file.c
	$(CXX) $(CFLAGS) $< -o $@

read_file.exe: read_file.c
	$(CXX) $(CFLAGS) $< -o $@

parse_args.exe: parse_args.c
	$(CXX) $(CFLAGS) $< -o $@

naive_benchmark.exe: naive_benchmark.c
	$(CXX) $(CFLAGS) -O0 $< -o naive_benchmark_O0.exe
	$(CXX) $(CFLAGS) -O1 $< -o naive_benchmark_O1.exe
	$(CXX) $(CFLAGS) -O2 $< -o naive_benchmark_O2.exe
	$(CXX) $(CFLAGS) -O3 $< -o naive_benchmark_O3.exe
	$(CXX) $(CFLAGS) -Os $< -o naive_benchmark_Os.exe
	$(CXX) $(CFLAGS) -Og $< -o naive_benchmark_Og.exe

run: stack_alloc.exe
	@time ./stack_alloc.exe

clean:
	rm -rf *.exe *.o
