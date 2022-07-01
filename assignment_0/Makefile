CXX := gcc
CFLAGS := -std=c99

.PHONY: all

all: stack_alloc.exe heap_alloc.exe

stack_alloc.exe: stack_alloc.c
	$(CXX) $(CFLAGS) $< -o $@

heap_alloc.exe: heap_alloc.c
	$(CXX) $(CFLAGS) $< -o $@

run: stack_alloc.exe
	@time ./stack_alloc.exe

clean:
	rm -rf *.exe *.o
