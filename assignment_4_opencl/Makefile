CXX := gcc
CFLAGS := -std=c11 -O3 -march=native
LDFLAGS = -lm -lOpenCL

SOURCES = $(wildcard *.c)
EXECS = $(patsubst %.c, %.exe, $(SOURCES))

.PHONY: all

all: diffusion

diffusion: diffusion.c
	$(CXX) $(CFLAGS) $< -o diffusion $(LDFLAGS)

test: diffusion
	# ./diffusion -n2 -d0.02
	./diffusion -d0.01 -n100000

compress: Makefile diffusion.c kernel.cl
	tar -czvf submission.tar.gz Makefile diffusion.c kernel.cl

check: compress clean
	# python3 ../../diffusion_opencl/check_submission.py submission.tar.gz
	python3 /home/hpc2022/diffusion_opencl/check_submission.py submission.tar.gz

clean:
	rm -rf *.exe *.o *.s ./diffusion extracted
