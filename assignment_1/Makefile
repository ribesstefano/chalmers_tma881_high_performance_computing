CXX := gcc
CFLAGS := -std=c17

SOURCES = $(wildcard *.c)
EXECS = $(patsubst %.c, %.exe, $(SOURCES))

.PHONY: all

all: $(EXECS)

stack_alloc.exe: stack_alloc.c
	$(CXX) $(CFLAGS) $< -o $@

naive_benchmark.exe: naive_benchmark.c
	$(CXX) -S $(CFLAGS) -O0 -fverbose-asm $< -o naive_benchmark_O0.s
	$(CXX) -S $(CFLAGS) -O1 -fverbose-asm $< -o naive_benchmark_O1.s
	$(CXX) -S $(CFLAGS) -O2 -fverbose-asm $< -o naive_benchmark_O2.s
	$(CXX) -S $(CFLAGS) -O3 -fverbose-asm $< -o naive_benchmark_O3.s
	$(CXX) -S $(CFLAGS) -Os -fverbose-asm $< -o naive_benchmark_Os.s
	$(CXX) -S $(CFLAGS) -Og -fverbose-asm $< -o naive_benchmark_Og.s

valgrind_test.exe: valgrind_test.c
	$(CXX) $(CFLAGS) -O0 -g $< -o $@

valgrind: valgrind_test.exe
	valgrind --leak-check=full \
	         --show-leak-kinds=all \
	         --track-origins=yes \
	         --verbose \
	         --log-file=valgrind-out.txt \
	         ./valgrind_test.exe

gdb_test.exe: gdb_test.c
	$(CXX) $(CFLAGS) -O0 -g $< -o $@

same_file.exe: same_file.c
	$(CXX) $(CFLAGS) -O3 $< -o $@

different_file_mul.exe: different_file_mul.c
	# Do nothing for "different_file_mul.c"

different_file.exe: different_file.c
	$(CXX) $(CFLAGS) -O3 different_file_mul.c $< -o $@

inlined_manually.exe: inlined_manually.c
	$(CXX) $(CFLAGS) -O3 $< -o $@

locality.exe: locality.c
	$(CXX) $(CFLAGS) -O0 $< -o locality_O0.exe
	$(CXX) $(CFLAGS) -O2 $< -o locality_O2.exe
	$(CXX) $(CFLAGS) -O2 -march=x86-64 $< -o locality_O2arch.exe

data_dependency.exe: data_dependency.c
	$(CXX) $(CFLAGS) -O0 $< -o data_dependency_O0.exe
	$(CXX) $(CFLAGS) -O2 $< -o data_dependency_O2.exe
	$(CXX) $(CFLAGS) -O2 -march=x86-64 $< -o data_dependency_O2arch.exe

indirect_addressing.exe: indirect_addressing.c
	$(CXX) $(CFLAGS) -O0 $< -o indirect_addressing_O0.exe
	$(CXX) $(CFLAGS) -O2 $< -o indirect_addressing_O2.exe
	$(CXX) $(CFLAGS) -O2 -march=x86-64 $< -o indirect_addressing_O2arch.exe

clean:
	rm -rf *.exe *.o *.s valgrind-out.txt*
