CXX := gcc
CFLAGS := -std=c11
LDFLAGS = -lm -lpthread

SOURCES = $(wildcard *.c)
EXECS = $(patsubst %.c, %.exe, $(SOURCES))

.PHONY: all

all: newton

newton: newton.c
	$(CXX) $(CFLAGS) -O3 -march=x86-64 $< -o newton $(LDFLAGS)

clean:
	rm -rf *.exe *.o *.s ./newton
