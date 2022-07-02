CXX := gcc
CFLAGS := -std=c99

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

run: stack_alloc.exe
	@time ./stack_alloc.exe

clean:
	rm -rf *.exe *.o
